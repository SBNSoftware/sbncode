#ifndef _SBN_WEIGHTMANAGER_H_
#define _SBN_WEIGHTMANAGER_H_

/**
 * \file WeightManager.h
 * \brief Allows to interface to EventWeight calculators
 *
 * Adapted from LArSoft's larsim EventWeight by A. Mastbaum
 * Original author: Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 */

#include "art/Framework/Principal/fwd.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "lardataobj/Simulation/sim.h"
#include "fhiclcpp/ParameterSet.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"
#include "WeightCalc.h"
#include "WeightCalcFactory.h"

namespace sbn {
  namespace evwgh {

class WeightManager {
public:
  WeightManager() {}
  ~WeightManager() {}

  /**
   * CONFIGURE FUNCTION
   *
   * 0) Looks at the weight_functions fcl parameter to get the name of the calculators   \n
   * 1) Creates the Calculators requested in step 0, and assigne a different random seed to each one \n
   * 3) The future call WeightManager::Run will run the calculators           \n
   *
   * @param cfg the input parameters for settings
   * @param the enging creator for the random seed (usually passed with *this)
   */
  template <typename EngineCreator>
  size_t Configure(fhicl::ParameterSet const& cfg, EngineCreator createEngine);

  /**
   * CORE FUNCTION: executes algorithms to assign a weight to the event as requested users. \n
   * WeightManager::Configure needs to be called first \n
   *
   * 0) Loos over all the previously emplaced calculators \n
   * 1) For each of them calculates the weights (more weight can be requested per calculator) \n
   * 3) Returns a map from "calculator name" to vector of weights calculated which is available inside EventWeightMap
   *
   * @param e the art event
   * @param inu the index of the simulated neutrino in the event
   */
  EventWeightMap Run(art::Event &e, const int inu);

  /**
   * Returns the map between calculator name and WeightCalcs
   */
  std::map<std::string, WeightCalc*> GetWeightCalcMap() { return fWeightCalcMap; }

private:
  std::map<std::string, WeightCalc*> fWeightCalcMap;  ///< A set of custom weight calculators
};


template <typename EngineCreator>
size_t WeightManager::Configure(fhicl::ParameterSet const& p, EngineCreator createEngine) {
  ::art::ServiceHandle<rndm::NuRandomService> seedservice;

  // Get list of weight functions
  auto const rw_func = p.get<std::vector<std::string> >("weight_functions");
  auto const module_label = p.get<std::string>("module_label");

  // Loop over all the functions and register them
  for (auto const& func : rw_func) {
    auto const ps_func = p.get<fhicl::ParameterSet>(func);
    std::string func_type = ps_func.get<std::string>("type");

    WeightCalc* wcalc = WeightCalcFactory::Create(func_type + "WeightCalc");
    if (wcalc == nullptr)
      throw cet::exception(__FUNCTION__) << "Function " << func << " requested in fcl file has not been registered!" << std::endl;

    if (fWeightCalcMap.find(func) != fWeightCalcMap.end())
      throw cet::exception(__FUNCTION__) << "Function " << func << " has been requested multiple times in fcl file!" << std::endl;

    // Create random engine for each rw function (name=func) (and seed it with random_seed set in the fcl)
    CLHEP::HepRandomEngine& engine = seedservice->registerAndSeedEngine(createEngine("HepJamesRandom", func),
                                                                        "HepJamesRandom", func, ps_func, "random_seed");

    wcalc->SetName(func);
    wcalc->SetType(func_type);
    wcalc->Configure(p, engine);

    fWeightCalcMap.emplace(func, wcalc);
  }

  return fWeightCalcMap.size();
}

  }  // namespace evwgh
}  // namespace sbn

#endif  // _SBN_WEIGHTMANAGER_H_
