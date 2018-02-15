////////////////////////////////////////////////////////////////////////
// Class:       EventWeight
// Module Type: producer
// File:        EventWeight_module.cc
//
// Generated at Fri Mar 20 09:36:11 2015 by Zarko Pavlovic using artmod
// from cetpkgsupport v1_08_04.
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

#include "Weight_t.h"
#include "MCEventWeight.h"
#include "WeightCalc.h"
#include "WeightCalcFactory.h"

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
  std::map<std::string, Weight_t*> fWeightCalcMap;  
  std::string fGenieModuleLabel;
};


EventWeight::EventWeight(fhicl::ParameterSet const & p)  {
   art::ServiceHandle<rndm::NuRandomService> seedservice;

  //get list of weight functions
  std::vector<std::string> rw_func=p.get<std::vector<std::string> >("weight_functions");
  fGenieModuleLabel = p.get< std::string > ("genie_module_label");

  for (auto ifunc=rw_func.begin(); ifunc!=rw_func.end(); ifunc++) {
    fhicl::ParameterSet const &ps_func=p.get<fhicl::ParameterSet> (*ifunc);
    std::string func_type=ps_func.get<std::string>("type");
    
    WeightCalc* wcalc=WeightCalcFactory::Create(func_type+"WeightCalc");
    if ( wcalc==NULL ) 
      throw cet::exception(__FUNCTION__) << "Function "<<*ifunc<<" requested in fcl file has not been registered!"<<std::endl;    
    if ( fWeightCalcMap.find(*ifunc)!=fWeightCalcMap.end() ) 
      throw cet::exception(__FUNCTION__) << "Function "<<*ifunc<<" has been requested multiple times in fcl file!"<<std::endl;
       
    mf::LogInfo("")<<"Configuring weight calculator "<<*ifunc;
    
    //create random engine for each rw function (name=*ifunc) (and seed it with random_seed set in the fcl)
    seedservice->createEngine(*this, "HepJamesRandom", *ifunc, ps_func, "random_seed");

    wcalc->SetName(*ifunc);
    wcalc->Configure(p);
    Weight_t* winfo=new Weight_t();
    winfo->fWeightCalcType=func_type;
    winfo->fWeightCalc=wcalc;
    winfo->fNmultisims=ps_func.get<int>("number_of_multisims");
      
    std::pair<std::string, Weight_t*> pwc(*ifunc,winfo);
    fWeightCalcMap.insert(pwc);
  }
 
  if ( fWeightCalcMap.size()>0 ) 
    produces<std::vector<MCEventWeight> >();
}


void EventWeight::produce(art::Event & e) {
  // Implementation of required member function here.
  std::unique_ptr<std::vector<MCEventWeight> > mcwghvec(new std::vector<MCEventWeight>);

  //get the MC generator information out of the event       
  //these are all handles to mc information.
  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;

  //actually go and get the stuff
  e.getByLabel(fGenieModuleLabel,mcTruthHandle);
  art::fill_ptr_vector(mclist, mcTruthHandle);
    

  // contains the mctruth object from genie
  for ( unsigned int inu=0; inu<mclist.size();inu++) { 
    MCEventWeight mcwgh;
    for (auto it=fWeightCalcMap.begin();it!=fWeightCalcMap.end();it++) {

      auto const & weights = it->second->GetWeight(e);
      std::cout<<"test of weights values in the module file:"<<  weights[0][0] <<std::endl; 
      if(weights.size() == 0){
	std::vector< double > empty;
	std::pair<std::string, std::vector <double> > p("empty",empty);
	mcwgh.insert(p);
      }	
      else{
	std::pair<std::string, std::vector<double> > 
	  p(it->first+"_"+it->second->fWeightCalcType,
	    weights[inu]);
	mcwgh.insert(p);
      }

    }
    (*mcwghvec).push_back(mcwgh);
  }

  e.put(std::move(mcwghvec));
}


void EventWeight::endJob() {
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
  for (auto it=fWeightCalcMap.begin();it!=fWeightCalcMap.end();it++) {
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

