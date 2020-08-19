/**
 *  @file   IHiggsStage.h
 *
 *  @brief  This provides an art tool interface definition for tools which can create
 *          fake particles to overlay onto input daq fragments during decoding
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IHiggsStage_h
#define IHiggsStage_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

// Algorithm includes
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/JamesRandom.h"

#include <utility>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {
/**
 *  @brief  IHiggsStage interface class definiton
 */
class IHiggsStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IHiggsStage() noexcept = default;

    IHiggsStage(const char * name) {
      // setup the random number engine
      art::ServiceHandle<rndm::NuRandomService> seedSvc;
      fEngine = new CLHEP::HepJamesRandom;
      seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(fEngine), name);
    }

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    virtual float MaxWeight() = 0;
    virtual float ConstantWeight() = 0;

protected:
    CLHEP::HepRandomEngine* fEngine;
};

} // namespace ldm
} // namespace evgen
#endif

