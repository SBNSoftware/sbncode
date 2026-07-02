#ifndef sbncode_WireMod_Utility_WireModUtility_hh
#define sbncode_WireMod_Utility_WireModUtility_hh

//std includes
#include <vector>

//ROOT includes
#include "TFile.h"
#include "TSpline.h"
#include "TGraph2D.h"
#include "TNtuple.h"

//Framework includes
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "sbnobj/ICARUS/TPC/ChannelROI.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <limits>
#include <memory>

namespace sys {
  class WireModUtility{
    // for now make everything public, though it's probably a good idea to think about what doesn't need to be
    public:
      // detector constants, should be set by geometry service
      // the notes at the end refer to their old names in the MircoBooNE code which preceded this
      // TODO: how best to initialize the splines/graphs?
      const geo::GeometryCore* geometry;                  // save the TPC geometry
      const geo::WireReadoutGeom* wireReadout;            // new for LarSoft v10
      const detinfo::DetectorPropertiesData& detPropData; // save the detector property data
      bool   applyChannelScale;                           // do we scale with channel?
      bool   applyXScale;                                 // do we scale with X?
      bool   applyYZScale;                                // do we scale with YZ?
      bool   applyXZAngleScale;                           // do we scale with XZ angle?
      bool   applyYZAngleScale;                           // do we scale with YZ angle?
      bool   applydEdXScale;                              // do we scale with dEdx?
      bool   applyXXZAngleScale;                          // do we scale with X vs XZ angle?
      bool   applyXdQdXScale;                             // do we scale with X vs dQ/dX?
      bool   applyXZAngledQdXScale;                       // do we scale with XZ angle vs dQ/dX?
      bool   applyXXWScale;                               // do we scale with X vs ThXW?
      bool   additiveModification;                        // additive (true) vs multiplicative (false) ROI modification
      double readoutWindowTicks;                          // how many ticks are in the readout window?
      double tickOffset;                                  // do we want an offset in the ticks?

      TSpline3* spline_Charge_Channel;                    // the spline for the charge correction by channel
      TSpline3* spline_Sigma_Channel;                     // the spline for the width correction by channel
      std::vector<TSpline3*> splines_Charge_X;            // the splines for the charge correction in X
      std::vector<TSpline3*> splines_Sigma_X;             // the splines for the width correction in X
      std::vector<TSpline3*> splines_Charge_XZAngle;      // the splines for the charge correction in XZ angle
      std::vector<TSpline3*> splines_Sigma_XZAngle;       // the splines for the width correction in XZ angle
      std::vector<TSpline3*> splines_Charge_YZAngle;      // the splines for the charge correction in YZ angle
      std::vector<TSpline3*> splines_Sigma_YZAngle;       // the splines for the width correction in YZ angle
      std::vector<TSpline3*> splines_Charge_dEdX;         // the splines for the charge correction in dEdX
      std::vector<TSpline3*> splines_Sigma_dEdX;          // the splines for the width correction in dEdX
      std::vector<TGraph2D*> graph2Ds_Charge_YZ;          // the graphs for the charge correction in YZ
      std::vector<TGraph2D*> graph2Ds_Sigma_YZ;           // the graphs for the width correction in YZ

      std::vector<TGraph2D*> graph2Ds_Charge_XXZAngle;    // the graphs for the charge correction in X vs XZ angle
      std::vector<TGraph2D*> graph2Ds_Sigma_XXZAngle;     // the graphs for the width correction in X vs XZ angle
      std::vector<TGraph2D*> graph2Ds_Charge_XdQdX;       // the graphs for charge correction in X vs dQ/dX
      std::vector<TGraph2D*> graph2Ds_Sigma_XdQdX;        // the graphs for width correction in X vs dQ/dX
      std::vector<TGraph2D*> graph2Ds_Charge_XZAngledQdX; // the graphs for charge correction in XZ angle vs dQ/dX
      std::vector<TGraph2D*> graph2Ds_Sigma_XZAngledQdX;  // the graphs for width correction in XZ angle vs dQ/dX
      std::vector<TGraph2D*> graph2Ds_Charge_XXW;         // the graphs for the charge correction in XXW
      std::vector<TGraph2D*> graph2Ds_Sigma_XXW;          // the graphs for the width correction in XXW

      // lets try making a constructor here
      // assume we can get a geometry service, a detector clcok, and a detector properties
      // pass the CryoStat and TPC IDs because it's IDs all the way down
      // set some optional args fpr the booleans, the readout window, and the offset
      WireModUtility(const geo::GeometryCore* geom, 
                     const geo::WireReadoutGeom* wireRead,
                     const detinfo::DetectorPropertiesData& detProp,
                     const bool& arg_ApplyChannelScale = false, 
                     const bool& arg_ApplyXScale = true,
                     const bool& arg_ApplyYZScale = true,
                     const bool& arg_ApplyXZAngleScale = true,
                     const bool& arg_ApplyYZAngleScale = true,
                     const bool& arg_ApplydEdXScale = true,
                     const bool& arg_ApplyXXZAngleScale = false,
                     const bool& arg_ApplyXdQdXScale = false,
                     const bool& arg_ApplyXZAngledQdXScale = false,
                     const bool& arg_ApplyXXWScale = true,
                     const double& arg_TickOffset = 0)
      : geometry(geom),
        wireReadout(wireRead),
        detPropData(detProp),
        applyChannelScale(arg_ApplyChannelScale),
        applyXScale(arg_ApplyXScale),
        applyYZScale(arg_ApplyYZScale),
        applyXZAngleScale(arg_ApplyXZAngleScale),
        applyYZAngleScale(arg_ApplyYZAngleScale),
        applydEdXScale(arg_ApplydEdXScale),
        applyXXZAngleScale(arg_ApplyXXZAngleScale),
        applyXdQdXScale(arg_ApplyXdQdXScale),
        applyXZAngledQdXScale(arg_ApplyXZAngledQdXScale),
        applyXXWScale(arg_ApplyXXWScale),
        additiveModification(false),
        readoutWindowTicks(detProp.ReadOutWindowSize()),                                               // the default A2795 (ICARUS TPC readout board) readout window is 4096 samples
        tickOffset(arg_TickOffset)                                                                     // tick offset is for MC truth, default to zero and set only as necessary
      {
      }

      // typedefs
      typedef std::pair<unsigned int, unsigned int>  ROI_Key_t;
      typedef std::pair<ROI_Key_t, unsigned int> SubROI_Key_t;

      typedef struct ROIProperties
      {
        ROI_Key_t key;
        raw::ChannelID_t channel;
        geo::View_t view;
        float begin;
        float end;
        float total_q;
        float center;   //charge weighted center of ROI
        float sigma;    //charge weighted RMS of ROI
      } ROIProperties_t;

      typedef struct SubROIProperties
      {
        SubROI_Key_t key;
        raw::ChannelID_t channel;
        geo::View_t view;
        float total_q;
        float center;
        float sigma;
      } SubROIProperties_t;

      typedef struct ScaleValues
      {
        double r_Q;
        double r_sigma;
      } ScaleValues_t;

      typedef struct TruthProperties
      {
        float x;
        float x_rms;
        float x_rms_noWeight;
        float tick;
        float tick_rms;
        float tick_rms_noWeight;
        float total_energy;
        float x_min;
        float x_max;
        float tick_min; // On Plane0!!
        float tick_max; // On Plane0!!
        float y;
        float z;
        double dxdr;
        double dydr;
        double dzdr;
        double dqdr;
        double dedr;
        ScaleValues_t scales_avg[3];
      } TruthProperties_t;

      // A single SimChannel ionization deposit (sim::IDE) tagged with its readout
      // coordinates. This is the SimChannel analogue of a sim::SimEnergyDeposit used by
      // the truth-matching path; unlike an edep, a sim::IDE carries only a point position
      // (x,y,z) -- direction/step come from the matched simb::MCParticle trajectory.
      typedef struct MatchedIDE
      {
        raw::ChannelID_t channel; // readout channel the charge landed on
        geo::WireID      wire;     // wire closest to the deposit
        double           tick;     // readout tick (TDC converted to tick)
        const sim::IDE*  ide;      // trackID, numElectrons, energy, x, y, z
      } MatchedIDE_t;

      // internal containers
      std::map< ROI_Key_t, std::vector<size_t> > ROIMatchedEdepMap;
      std::map< ROI_Key_t, std::vector<size_t> > ROIMatchedHitMap;

      // SimChannel/IDE truth-matching containers (parallel to the Edep versions).
      // fIDEVec is the flat owner; ROIMatchedIDEMap stores indices into it, keyed by ROI.
      std::vector<MatchedIDE_t>                  fIDEVec;
      std::map< ROI_Key_t, std::vector<size_t> > ROIMatchedIDEMap;

      // Space-charge provider, set by the module (lar::providerFrom<SpaceChargeService>()).
      // Used to map IDE (at-the-wire) positions to/from the true MCParticle trajectory
      // frame. When null or SCE disabled, the mappings are the identity.
      const spacecharge::SpaceCharge* fSCE = nullptr;

      // some useful functions
      // geometries
      // TODO is this the most efficient for new v10 iterators?
      double planeXToTick(double xPos, const geo::PlaneGeo& plane, const geo::TPCGeo& tpcGeom, double offset = 0) {
          return detPropData.ConvertXToTicks(xPos, plane.ID()) + offset;
      }

      bool planeXInWindow(double xPos, const geo::PlaneGeo& plane, const geo::TPCGeo& tpcGeom, double offset = 0)
      {
        double tick = planeXToTick(xPos, plane, tpcGeom, offset);
        return (tick > 0 && tick <= detPropData.ReadOutWindowSize());
      }
      
      // for this function: in the future if we want to use non-gaussian functions make this take a vector of parameters
      // the another wiremod utility could overwrite the ``fitFunc'' with some non-standard function
      // would require a fair bit of remodling (ie q and sigma would need to be replace with, eg, funcVar[0] and funcVar[1] and probs a bunch of loops)
      // so lets worry about that later
      double gausFunc(double t, double mean, double sigma, double a = 1.0)
      {
        return (a / (sigma * std::sqrt(2 * util::pi()))) * std::exp(-0.5 * std::pow((t - mean)/sigma, 2));
      }

      double FoldAngle(double theta)
      {
        return (std::abs(theta) > 0.5 * util::pi()) ? util::pi() - std::abs(theta) : std::abs(theta);
      }

      double ThetaXZ_PlaneRel(double dxdr, double dydr, double dzdr, double planeAngle)
      {
        double sinPlaneAngle = std::sin(planeAngle);
        double cosPlaneAngle = std::cos(planeAngle);

        //double dydrPlaneRel = dydr * cosPlaneAngle - dzdr * sinPlaneAngle; // don't need to rotate Y for this angle
        double dzdrPlaneRel = dzdr * cosPlaneAngle + dydr * sinPlaneAngle;

        double theta = std::atan2(dxdr, dzdrPlaneRel);
        return FoldAngle(theta);
      }

      double ThetaYZ_PlaneRel(double dxdr, double dydr, double dzdr, double planeAngle)
      {
        double sinPlaneAngle = std::sin(planeAngle);
        double cosPlaneAngle = std::cos(planeAngle);

        double dydrPlaneRel = dydr * cosPlaneAngle - dzdr * sinPlaneAngle;
        double dzdrPlaneRel = dzdr * cosPlaneAngle + dydr * sinPlaneAngle;

        double theta = std::atan2(dydrPlaneRel, dzdrPlaneRel);
        return FoldAngle(theta);
      }
      // theste are set in the .cc file

      double ThetaXY_PlaneRel(double dxdr, double dydr, double dzdr, double planeAngle)
      {
        double sinPlaneAngle = std::sin(planeAngle);
        double cosPlaneAngle = std::cos(planeAngle);

        double dydrPlaneRel = dydr * cosPlaneAngle - dzdr * sinPlaneAngle;
        //double dzdrPlaneRel = dzdr * cosPlaneAngle + dydr * sinPlaneAngle; // don't need to rotate Z for this angle

        double theta = std::atan2(dxdr, dydrPlaneRel);
        return FoldAngle(theta);
      }

      double ThetaXW(double dxdr, double dydr, double dzdr, double planeAngle)
      {
        // planeAngle is the wire angle from +z (PlaneGeo::ThetaZ, numerically equal to
        // WireReadoutGeom::WireAngleToVertical). The pitch "gamma" angle is measured from the
        // wire-normal, i.e. planeAngle - pi/2 -- the same convention as the canonical LArSoft
        // pitch (Calorimetry / TrackCaloSkimmer use WireAngleToVertical - pi/2). Subtracting a
        // full pi projected the direction onto the wire itself (a cos<->sin swap), yielding the
        // wrong angle; subtract pi/2 so cosG is the projection onto the wire-normal (pitch).
        double sin = std::sin(planeAngle - 0.5 * ::util::pi<>());
        double cos = std::cos(planeAngle - 0.5 * ::util::pi<>());

        double cosG = std::abs(dydr * sin + dzdr * cos);
        double theta = std::atan(dxdr / cosG);
        return std::abs(theta);
      }
      
      // these are set in the .cc file
      ROIProperties_t CalcROIProperties(recob::Wire const&, size_t const&);
      ROIProperties_t CalcROIProperties(recob::ChannelROI const&, size_t const&);

      std::vector<std::pair<unsigned int, unsigned int>> GetTargetROIs(sim::SimEnergyDeposit const&, double offset);
      std::vector<std::pair<unsigned int, unsigned int>> GetHitTargetROIs(recob::Hit const&);

      void FillROIMatchedEdepMap(std::vector<sim::SimEnergyDeposit> const&, std::vector<recob::Wire> const&, double offset);
      void FillROIMatchedHitMap(std::vector<recob::Hit> const&, std::vector<recob::Wire> const&);

      // SimChannel/IDE versions of FillROIMatchedEdepMap. The DetectorClocksData is used to
      // convert each IDE's TDC to a readout tick (no X->tick projection / tickOffset needed,
      // since the SimChannel already carries the readout coordinate).
      void FillROIMatchedIDEMap(std::vector<sim::SimChannel> const& simchVec, std::vector<recob::Wire> const& wireVec, detinfo::DetectorClocksData const& clockData, double offset = 0);
      void FillROIMatchedIDEMap(std::vector<sim::SimChannel> const& simchVec, std::vector<recob::ChannelROI> const& chanROIVec, detinfo::DetectorClocksData const& clockData, double offset = 0);

      void FillROIMatchedEdepMap(std::vector<sim::SimEnergyDeposit> const&, std::vector<recob::ChannelROI> const&, double offset);
      void FillROIMatchedHitMap(std::vector<recob::Hit> const&, std::vector<recob::ChannelROI> const&);

      std::vector<SubROIProperties_t> CalcSubROIProperties(ROIProperties_t const&, std::vector<const recob::Hit*> const&);

      std::map<SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> MatchEdepsToSubROIs(std::vector<SubROIProperties_t> const&, std::vector<const sim::SimEnergyDeposit*> const&, double offset);

      // SimChannel/IDE analogue of MatchEdepsToSubROIs. Keyed on the IDE's native readout
      // tick (no per-plane X->tick projection), otherwise the same center+/-sigma / closest
      // sub-ROI assignment as the edep version.
      std::map<SubROI_Key_t, std::vector<const MatchedIDE_t*>> MatchIDEsToSubROIs(std::vector<SubROIProperties_t> const&, std::vector<const MatchedIDE_t*> const&);

      TruthProperties_t CalcPropertiesFromEdeps(std::vector<const sim::SimEnergyDeposit*> const&, double offset);

      // SimChannel/IDE analogue of CalcPropertiesFromEdeps. Fills TruthProperties_t from the
      // dominant track's sim::IDEs (charge-weighted position/energy) plus the matched
      // simb::MCParticle trajectory (direction / pitch / dedr / dqdr). particleMap maps
      // abs(trackID) -> MCParticle so the dominant track's trajectory can be looked up.
      TruthProperties_t CalcPropertiesFromIDEs(std::vector<const MatchedIDE_t*> const&, std::map<int, const simb::MCParticle*> const& particleMap);

      // Space-charge coordinate maps (mirroring TrackCaloSkimmer). Identity if fSCE is null
      // or SCE disabled.
      geo::Point_t WireToTrajectoryPosition(const geo::Point_t& loc, const geo::TPCID& tpc) const;          // at-the-wire -> true trajectory frame
      geo::Point_t TrajectoryToWirePosition(const geo::Point_t& loc, const geo::Vector_t& driftdir) const;  // true -> at-the-wire frame

      ScaleValues_t GetScaleValues(TruthProperties_t const&, ROIProperties_t const&);
      ScaleValues_t GetChannelScaleValues(TruthProperties_t const&, raw::ChannelID_t const&);
      ScaleValues_t GetViewScaleValues(TruthProperties_t const&, geo::View_t const&);

      void ModifyROI(std::vector<float> &,
                     ROIProperties_t const &,
                     std::vector<SubROIProperties_t> const&,
                     std::map<SubROI_Key_t, ScaleValues_t> const&,
                     double sigmaWindow = std::numeric_limits<double>::infinity());
  }; // end class
} // end namespace

#endif // sbncode_WireMod_Utility_WireModUtility_hh
