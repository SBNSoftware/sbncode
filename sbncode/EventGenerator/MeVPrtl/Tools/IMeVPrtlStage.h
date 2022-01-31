/**
 *  @file   IMeVPrtlStage.h
 *
 *  @brief  This provides an art tool interface definition for tools which can create
 *          fake particles to overlay onto input daq fragments during decoding
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IMeVPrtlStage_h
#define IMeVPrtlStage_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

// Algorithm includes
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "TVector3.h"

#include <utility>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {
/**
 *  @brief  IMeVPrtlStage interface class definiton. General interface behind each
 *  stage. Provides random number generation.
 */
class IMeVPrtlStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IMeVPrtlStage() noexcept {
      if (fEngine) delete fEngine;
    }

    IMeVPrtlStage(const char * name) {
      // setup the random number engine
      art::ServiceHandle<rndm::NuRandomService> seedSvc;
      fEngine = new CLHEP::HepJamesRandom;
      seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(fEngine), name);
      fName = name;
    }

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    virtual double MaxWeight() = 0;

    // useful helper function
    TVector3 RandomUnitVector() {
      // In order to pick a random point on a sphere -- pick a random value of _costh_, __not__ theta
      // b.c. d\Omega = d\phi dcos\theta, i.e. d\Omega != d\phi d\theta
      double costheta = CLHEP::RandFlat::shoot(fEngine, -1, 1);
      double sintheta = sqrt(1. - costheta * costheta);
      double phi = CLHEP::RandFlat::shoot(fEngine, 0, 2*M_PI);
      return TVector3(sintheta * cos(phi), sintheta * sin(phi), costheta);
    }
    double GetRandom() {
      return CLHEP::RandFlat::shoot(fEngine);
    }

    const char *Name() { return fName; }

protected:
    CLHEP::HepRandomEngine* fEngine;
    const char *fName;
};

} // namespace ldm
} // namespace evgen
#endif

