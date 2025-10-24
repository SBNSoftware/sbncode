////////////////////////////////////////////////////////////////////////
/// \file  CORSIKAGenSBN.h
/// \brief Class definition for cosmic-ray generator.
///
/// \author  gputnam@fnal.gov, taken from module by: Matthew.Bass@physics.ox.ac.uk
////////////////////////////////////////////////////////////////////////
//
// ROOT includes
#include "TDatabasePDG.h"
#include "TString.h"
#include "TSystem.h" //need BaseName and DirName

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "ifdh.h" //to handle flux files
#include <sqlite3.h>

namespace evgen {

  /**
   * @brief LArSoft interface to CORSIKA event generator.
   *
   * @note A presentation on this module by the original author is archived at:
   *       https://indico.fnal.gov/event/10893/contribution/3/material/slides
   *
   * In CORSIKA jargon, a "shower" is the cascade of particles resulting from
   * a primary cosmic ray interaction.
   * This module creates a single `simb::MCTruth` object (stored as data product
   * into a `std::vector<simb::MCTruth>` with a single entry) containing all
   * the particles from cosmic ray showers crossing a _surface_ above the
   * detector.
   *
   * The generation procedure consists of selecting showers from a database of
   * *pregenerated* events, and then to adapt them to the parameters requested
   * in the module configuration. Pregenerated showers are "observed" at a
   * altitude set in CORSIKA configuration.
   *
   * Databases need to be stored as files in SQLite3 format. Multiple file
   * sources can be specified (`ShowerInputFiles` configuration parameter).
   * From each source, one database file is selected and copied locally via
   * IFDH.
   * From each source, showers are extracted proportionally to the relative flux
   * specified in the configuration (specified in `ShowerFluxConstants`,
   * see @ref CORSIKAGenSBN_Normalization "normalization" below).
   * The actual number of showers per event and per source is extracted
   * according to a Poisson distribution around the predicted average number of
   * primary cosmic rays for that source.
   *
   *
   * Flux normalization
   * -------------------
   *
   * @anchor CORSIKAGenSBN_Normalization
   *
   * CORSIKA generates showers from each specific cosmic ray type @f$ A @f$
   * (e.g. iron, proton, etc.) according to a power law distribution
   * @f$ \Phi_{A}(E) \propto E^{-\gamma_{A}} @f$ of the primary
   * particle energy @f$ E @f$ [GeV]. When sampling pregenerated events, we
   * bypass the normalization imposed by CORSIKA and gain complete
   * control on it.
   *
   * Within CORSIKAGenSBN, for each source (usually each on a different primary
   * cosmic ray type, e.g. iron, proton, etc.), the average number of generated
   * showers is
   * @f$ n_{A} = \pi S T k_{A} \int E^{-\gamma_{A}} dE @f$ with @f$ S @f$ the
   * area of the surface the flux passes across, @f$ T @f$ the exposure time,
   * the integral defined over the full energy range of the pregenerated showers
   * in the source, and @f$ k_{A} @f$ a factor specified in the configuration
   * (`ShowerFluxConstants` parameters).
   * This is the flux of primary cosmic rays, not of the observed particles
   * from their showers. Note that it depends on an area and a time interval,
   * but it is uniform with respect to translations and constant in time.
   *
   * As explained @ref CORSIKAGenSBN_Coverage "below", we consider only the
   * secondary particles that cross an "observation" surface.
   * After cosmic ray primary particles cross the flux surface (@f$ S_{\Phi} @f$
   * above), they develop into showers of particles that spread across large
   * areas. Limiting ourself to the observation of particles on a small surface
   * @f$ S_{o} @f$ has two effects. We lose the part of the showers that misses
   * that surface @f$ S_{o} @f$. Also, considering a span of time with
   * multiple showers, we miss particles from other hypothetical showers
   * whose primaries crossed outside the generation surface @f$ S_{\Phi} @f$
   * whose shower development would leak into @f$ S_{o} @f$: we did not simulate
   * those showers at all!
   * In terms of total flux of observed particles, under the assumption that the
   * flux is uniform in space, choosing @f$ S_{\Phi} @f$ the same size as
   * @f$ S_{o} @f$ makes the two effects get the same size: just as many
   * particles from primaries from @f$ S_{\Phi} @f$ leak out of @f$ S_{o} @f$,
   * as many particles from primaries from outside @f$ S_{\Phi} @f$ sneak in
   * @f$ S_{o} @f$.
   * In that case, counting _all_ the particles from the primaries crossing a
   * surface @f$ S_{\Phi} @f$ of area _S_ regardless of where they land
   * yields the right amount of cosmic ray secondary particles across an
   * observation surface @f$ S_{o} @f$ also of area _S_. Practically,
   * the particles landing outside @f$ S_{o} @f$ need to be recovered to
   * preserve the correct normalization, which is described in the
   * @ref CORSIKAGenSBN_Coverage "next section".
   *
   *
   * Surface coverage, position and timing
   * --------------------------------------
   *
   * @anchor CORSIKAGenSBN_Coverage
   *
   * The surface we detect the particles through (let's call it @f$ S_{d} @f$)
   * is defined by the smallest _rectangle_ including all cryostats in the
   * detector, and located at the height of the ceiling of the tallest cryostat.
   * This surface can be increased by specifying a positive value for
   * `ShowerAreaExtension` configuration parameter, in which case each side of
   * the rectangle will be extended by that amount.
   *
   * Showers are extracted one by one from the pregenerated samples and treated
   * independently.
   * Ideally, the detection surface @f$ S_{d} @f$ would be at the same exact
   * altitude as the observation surface set in CORSIKA (called @f$ S_{o} @f$
   * above). In practice, we go the other way around, with the assumption that
   * the shower observed at @f$ S_{d} @f$ would be very close to the one we
   * actually generated at @f$ S_{o} @f$, and teleport the generated particles
   * on @f$ S_{d} @f$. Since the cryostats may be just meters from the earth
   * surface @f$ S_{o} @f$ lies on, this is an acceptable approximation.
   *
   * All the particles of one shower are compelled within surface @f$ S_{d} @f$
   * as a first step. As explained when describing the
   * @anchor CORSIKAGenSBN_Normalization "normalization", we need to keep all the
   * shower particles, one way or the other.
   * So, particles of the shower that fell out of @f$ S_{d} @f$ are repackaged
   * into other showers and translated back in. This is equivalent to assume the
   * primary cosmic rays originating such shower would happen outside
   * our generation volume (@f$ S_{\Phi} @f$), and their shower would then spill
   * into @f$ S_{d} @f$. Since these repackaged showers are in principle
   * independent and uncorrelated, they are assigned a random time different
   * than the main shower, leveraging the assumption of constantness of the
   * flux.
   *
   * As for the azimuth, this module uses an approximation by setting north
   * direction to match the _z_ axis of LArSoft geometry (typically assumed
   * to be the direction of the beam particle).
   *
   * The particles so manipulated are then back-propagated from the observation
   * surface to an absolute height defined by `ProjectToHeight` (although for
   * particular combination of position and direction, the particles might be
   * propagated back to the edge of the world, or even outside it).
   *
   * As a final filter, only the particles whose straight projections cross any
   * of the cryostats (with some buffer volume around, defined by `BufferBox`)
   * are stored, while the other ones are discarded. Note that the actual
   * interactions that particles from the generation surface undergo may deviate
   * them enough to miss the cryostats anyway, and conversely particles that
   * have been filtered out because shooting off the cryostats might have been
   * subsequently deviated to actually cross them. This effect is not corrected
   * for at this time.
   *
   * The time of the showers is uniformly distributed within the configured
   * time interval, defined by `SampleTime` starting from `TimeOffset`.
   *
   *
   * Configuration parameters
   * =========================
   *
   * * `ShowerInputFiles` (list of paths; mandatory): a list of file paths to
   *     pregenerated CORSIKA shower files. Each entry can be a single file or
   *     use wildcards (`*`) to specify a set of files to choose among.
   *     Paths and wildcards are processed by IFDH.
   * * `ShowerFluxConstants` (list of real numbers; mandatory): for each entry
   *     @f$ A @f$ in `ShowerInputFiles`, specify the normalization factor
   *     @f$ K_{A} @f$ of their distribution [@f$ m^{-2}s^{-1} @f$]
   * * `ProjectToHeight` (real, default: `0`): the generated particles will
   *     appear to come from this height [cm]
   * * `TimeOffset` (real; default: `0`): start time of the exposure window [s],
   *     relative to the
   *     @ref DetectorClocksSimulationTime "simulation time start"
   * * `SampleTime` (real; mandatory): duration of the simulated exposure to
   *     cosmic rays [s]
   * * `ShowerAreaExtension` (real; default: `0`):
   *     extend the size of the observation surface of shower particles (_S_)
   *     by this much [cm]; e.g. 1000 will extend 10 m on each side
   * * `RandomXZShift` (real; default: `0`): the original position of each
   *     shower is randomly shifted within a square with this length as side
   *     [cm]
   * * `BufferBox` (list of six lengths, all `0` by default):
   *     extension to the volume of each cryostat for the purpose of filtering
   *     out the particles which do not cross the detector; each cryostat volume
   *     is independently extended by the same amount, specified here as
   *     *shifts* to lower _x_, higher _x_, lower _y_, higher _y_, lower _z_
   *     and higher _z_, in that order [cm] (note that to extend e.g. the
   *     negative _x_ side by 5 meters the parameter value should be -500)
   * * `SeedGenerator` (integer): force random number generator for event
   *     generation to the specified value
   * * `SeedPoisson` (integer): force random number generator for number of
   *     showers to the specified value
   * * `Seed`: alias for `SeedGenerator`
   *
   *
   * Random engines
   * ---------------
   *
   * Currently two random engines are used:
   * * a generator engine (driven by `SeedGenerator`), of general use
   * * a "Poisson" engine (driven by `SeedPoisson`), only used to determine the
   *     number of showers to be selected on each event
   *
   */
  class CORSIKAGenSBN : public art::EDProducer {
  public:
    explicit CORSIKAGenSBN(fhicl::ParameterSet const& pset);
    virtual ~CORSIKAGenSBN();

    void produce(art::Event& evt);
    void beginRun(art::Run& run);

  protected:
    void openDBs(std::string const& module_label);
    void populateNShowers();
    void populateTOffset();
    void GetSample(simb::MCTruth&);
    double wrapvar(const double var, const double low, const double high);
    double wrapvarBoxNo(const double var, const double low, const double high, int& boxno);
    /**
     * @brief Propagates a point back to the surface of a box.
     * @param xyz coordinates of the point to be propagated (`{ x, y, z }`)
     * @param dxyz direction of the point (`{ dx, dy, dz }`)
     * @param xlo lower _x_ coordinate of the target box
     * @param xhi upper _x_ coordinate of the target box
     * @param ylo lower _y_ coordinate of the target box
     * @param yhi upper _y_ coordinate of the target box
     * @param zlo lower _z_ coordinate of the target box
     * @param zhi upper _z_ coordinate of the target box
     * @param xyzout _(output, room for at least 3 numbers)_ propagated point
     *
     * The point `xyz`, assumed to be inside the box, is propagated at the level
     * of _the closest among the sides of the box_. Note that this means the
     * propagated point might still be not on the surface of the box, even if it
     * would later reach it.
     * It is anyway guaranteed that `xyzout` is not inside the target box, and
     * that at least one of its three coordinates is on the box surface.
     */
    void ProjectToBoxEdge(const double xyz[],
                          const double dxyz[],
                          const double xlo,
                          const double xhi,
                          const double ylo,
                          const double yhi,
                          const double zlo,
                          const double zhi,
                          double xyzout[]);

    int fShowerInputs = 0; ///< Number of shower inputs to process from
    std::vector<double>
      fNShowersPerEvent; ///< Number of showers to put in each event of duration fSampleTime; one per showerinput
    std::vector<int> fMaxShowers; //< Max number of showers to query, one per showerinput
    double fShowerBounds[6] = {
      0.,
      0.,
      0.,
      0.,
      0.,
      0.}; ///< Boundaries of area over which showers are to be distributed (x(min), x(max), _unused_, y, z(min), z(max) )
    double fToffset_corsika =
      0.; ///< Timing offset to account for propagation time through atmosphere, populated from db
    ifdh_ns::ifdh* fIFDH = 0; ///< (optional) flux file handling

    //fcl parameters
    double fProjectToHeight = 0.; ///< Height to which particles will be projected [cm]
    std::vector<std::string> fShowerInputFiles; ///< Set of CORSIKA shower data files to use
    std::string fShowerCopyType; // ifdh file selection and copying (default) or direct file path
    std::vector<double>
      fShowerFluxConstants;  ///< Set of flux constants to be associated with each shower data file
    double fSampleTime = 0.; ///< Duration of sample [s]
    double fToffset = 0.;    ///< Time offset of sample, defaults to zero (no offset) [s]
    std::vector<double>
      fBuffBox; ///< Buffer box extensions to cryostat in each direction (6 of them: x_lo,x_hi,y_lo,y_hi,z_lo,z_hi) [cm]
    double fShowerAreaExtension =
      0.; ///< Extend distribution of corsika particles in x,z by this much (e.g. 1000 will extend 10 m in -x, +x, -z, and +z) [cm]
    sqlite3* fdb[5]; ///< Pointers to sqlite3 database object, max of 5
    double fRandomXZShift =
      0.; ///< Each shower will be shifted by a random amount in xz so that showers won't repeatedly sample the same space [cm]
    CLHEP::HepRandomEngine& fGenEngine;
    CLHEP::HepRandomEngine& fPoisEngine;
  };
}

