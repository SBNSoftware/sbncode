////////////////////////////////////////////////////////////////////////
// Class:       EventWeight
// Module Type: producer
// File:        EventWeight_module.cc
//
// Generated at Fri Mar 20 09:36:11 2015 by Zarko Pavlovic using artmod
// from cetpkgsupport v1_08_04.
//
// Ported from uboonecode; see LICENSE for details.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "lardataobj/Simulation/sim.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Persistency/Common/Assns.h" 
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <iomanip>

#include "larsim/EventWeight/Base/Weight_t.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larsim/EventWeight/Base/WeightCalc.h"
#include "larsim/EventWeight/Base/WeightCalcFactory.h"
#include "larsim/EventWeight/Base/WeightManager.h"

#include "nusimdata/SimulationBase/MCTruth.h"

namespace sbncode {
  namespace evwgh {

class EventWeight : public art::EDProducer {
public:
  explicit EventWeight(fhicl::ParameterSet const & p);

  EventWeight(EventWeight const &) = delete;
  EventWeight(EventWeight &&) = delete;
  EventWeight & operator = (EventWeight const &) = delete;
  EventWeight & operator = (EventWeight &&) = delete;

  void produce(art::Event & e) override;

  void endJob() override;

private:
  ::evwgh::WeightManager _wgt_manager;
  std::string fGenieModuleLabel;
};


EventWeight::EventWeight(fhicl::ParameterSet const & p)  {
  size_t n_func = _wgt_manager.Configure(p, *this);
  fGenieModuleLabel = p.get<std::string>("genie_module_label", "generator");
  produces<std::vector<evwgh::MCEventWeight> >();         
}


void EventWeight::produce(art::Event & e) {
  std::unique_ptr<std::vector<MCEventWeight> > mcwghvec(new std::vector<MCEventWeight>);

  // Get the MC generator information out of the event       
  art::Handle<std::vector<simb::MCTruth> > mcTruthHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;

  // Actually go and get the stuff
  e.getByLabel(fGenieModuleLabel, mcTruthHandle);
  if (!mcTruthHandle.isValid()) {
    throw cet::exception(__FUNCTION__)
      << "Can't find GENIE module label with name "
      << fGenieModuleLabel << std::endl;
  }
  art::fill_ptr_vector(mclist, mcTruthHandle);
    
  // Loop over all neutrinos in this event
  for (unsigned int inu=0; inu<mclist.size(); inu++) {   
    evwgh::MCEventWeight mcwgh = _wgt_manager.Run(e, inu);
    (*mcwghvec).push_back(mcwgh);
  }

  e.put(std::move(mcwghvec));
}


void EventWeight::endJob() {
  std::map<std::string, Weight_t*> weightCalcMap = \
    _wgt_manager.GetWeightCalcMap();

  std::stringstream job_summary;
  job_summary<<std::setprecision(2);
  for (int i=1;i<=110;i++) job_summary<<"=";
  job_summary<<std::endl;
  job_summary<<std::setw(20)<<"WeightCalc"
             <<std::setw(15)<<"Type"
             <<std::setw(15)<<"#RW neutrinos"
             <<std::setw(15)<<"#Multisims"
             <<std::setw(15)<<"Min"
             <<std::setw(15)<<"Max"
             <<std::setw(15)<<"Avg"
             <<std::endl;
  for (int i=1;i<=110;i++) job_summary<<"=";
  job_summary<<std::endl;
  for (auto it=fWeightCalcMap.begin();it!=weightCalcMap.end();it++) {
    job_summary<<std::setw(20)<<it->first
      	 <<std::setw(15)<<(it->second->fWeightCalcType)
      	 <<std::setw(15)<<(it->second->fNcalls)
      	 <<std::setw(15)<<(it->second->fNmultisims)
      	 <<std::setw(15)<<(it->second->fMinWeight)
      	 <<std::setw(15)<<(it->second->fMaxWeight)
      	 <<std::setw(15)<<(it->second->fAvgWeight)
      	 <<std::endl;
  }
  for (int i=1;i<=110;i++) job_summary<<"=";
  job_summary<<std::endl;
  mf::LogInfo("")<<job_summary.str();

}

  }
}

DEFINE_ART_MODULE(sbncode::evwgh::EventWeight)

