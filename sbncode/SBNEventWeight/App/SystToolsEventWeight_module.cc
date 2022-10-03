////////////////////////////////////////////////////////////////////////
// Class: SystToolsEventWeight
// Module Type: producer
// File: SystToolsEventWeight_module.cc
//
// Generated at Mon Feb 21 09:36:11 2015 by Jaesung Kim (jae.sung.kim.3426@gmail.com)
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <assert.h>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightParameterSet.h"

#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interpreters/ParamHeaderHelper.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Principal/Run.h"

namespace sbn {
  namespace evwgh {

class SystToolsEventWeight : public art::EDProducer {
public:
  explicit SystToolsEventWeight(fhicl::ParameterSet const& p);

  SystToolsEventWeight(SystToolsEventWeight const &) = delete;
  SystToolsEventWeight(SystToolsEventWeight &&) = delete;
  SystToolsEventWeight& operator = (SystToolsEventWeight const&) = delete;
  SystToolsEventWeight& operator = (SystToolsEventWeight&&) = delete;

private:
  void produce(art::Event& e) override;
  void beginRun(art::Run& run) override;

private:
  std::string fGenieModuleLabel;
  bool fAllowMissingTruth;
  bool fDebugMode;

  systtools::provider_list_t fSystProviders;
  systtools::param_header_map_t fParamHeaderMap;
  systtools::ParamHeaderHelper fParamHeaderHelper;

};


SystToolsEventWeight::SystToolsEventWeight(fhicl::ParameterSet const& p)
  : EDProducer{p}
{

  fGenieModuleLabel = p.get<std::string>("generator_module_label", "generator");
  fAllowMissingTruth = p.get<bool>("AllowMissingTruth");
  fDebugMode = p.get<bool>("debugMode");

  fhicl::ParameterSet syst_provider_config = p.get<fhicl::ParameterSet>("generated_systematic_provider_configuration");

  MF_LOG_INFO("SystToolsEventWeight") << "Configuring ISystProvider";

  fSystProviders = systtools::ConfigureISystProvidersFromParameterHeaders(syst_provider_config);
  fParamHeaderMap = systtools::BuildParameterHeaders(fSystProviders);
  fParamHeaderHelper = systtools::ParamHeaderHelper(fParamHeaderMap);

  MF_LOG_INFO("SystToolsEventWeight") << "ISystProvider is configuered" << std::endl;

  produces<std::vector<sbn::evwgh::EventWeightMap> >();
  produces<art::Assns<simb::MCTruth, sbn::evwgh::EventWeightMap> >();
  produces<std::vector<sbn::evwgh::EventWeightParameterSet>, art::InRun>();

}


void SystToolsEventWeight::produce(art::Event& e) {

  auto mcwghvec = std::make_unique<std::vector<sbn::evwgh::EventWeightMap> >();
  auto wghassns = std::make_unique<art::Assns<simb::MCTruth, sbn::evwgh::EventWeightMap> >();

  art::PtrMaker<sbn::evwgh::EventWeightMap> makeWeightPtr(e);

  std::vector<art::Ptr<simb::MCTruth> > mclist;
  art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
  if(!fGenieModuleLabel.empty()) e.getByLabel(fGenieModuleLabel, mcTruthHandle);
  // Prooceed even with missing handle if we want to require the MCTruth to be
  // found, so that an exception will be thrown explaining the problem.
  if(mcTruthHandle.isValid() || !fAllowMissingTruth){

    art::fill_ptr_vector(mclist, mcTruthHandle);

    if(fDebugMode){
      std::cout << "[SystToolsEventWeight::produce] mclist.size() = " << mclist.size() << "\n"
                << "[SystToolsEventWeight::produce] fSystProviders.size() = " << fSystProviders.size() << std::endl;
    }

    for( auto &sp : fSystProviders ){

      std::unique_ptr<systtools::EventResponse> syst_resp = sp->GetEventResponse(e);
      if( !syst_resp ){
        std::cout << "[ERROR]: Got nullptr systtools::EventResponse from provider "
                  << sp->GetFullyQualifiedName();
        throw std::exception();
      }

      //==== syst_resp->size() is (Number of MCTruth)
      //==== Each index corresponds to each of MCTruth
      int nMCTruthIndex(0);
      if(fDebugMode) std::cout << "[SystToolsEventWeight::produce]   syst_resp.size() (= Number of MCTruth) of this SystProvider = " << syst_resp->size() << std::endl;

      //==== Looping over syst_resp is identical to looping over MCTruth
      for( systtools::EventResponse::iterator itResp = syst_resp->begin(); itResp != syst_resp->end(); ++itResp ){

        systtools::event_unit_response_t resp = *itResp;
        //==== resp.size() corresponds to number of knobs we altered;
        //==== e.g., MaCCQE, MaCCRES, MvCCRE -> resp.size() = 3
        if(fDebugMode){
          std::cout << "[SystToolsEventWeight::produce]     sp->GetSystMetaData().size() (expected) = " << sp->GetSystMetaData().size() << "\n"
                    << "[SystToolsEventWeight::produce]     resp.size() of this syst_resp (produced) = " << resp.size() << std::endl;
        }
        if(sp->GetSystMetaData().size()!=resp.size()){
          std::cerr << "[SystToolsEventWeight::produce]     sp->GetFullyQualifiedName() = " << sp->GetFullyQualifiedName() << std::endl;
          std::cerr << "[SystToolsEventWeight::produce]     We expect to have " << sp->GetSystMetaData().size() << " knobs for this SystProvider, but "
                    << resp.size() << " are produced. "
                    << "Probably this particular event is not relevant to this systematic variation." << std::endl;
          throw std::exception();
        }

        sbn::evwgh::EventWeightMap mcwgh;

        for( systtools::event_unit_response_t::iterator it = resp.begin(); it != resp.end(); ++it ){

          //==== responses.size() is the number of universes
          systtools::SystParamHeader const &sph = fParamHeaderHelper.GetHeader( (*it).pid );
          std::string prettyName = sph.prettyName;

          if(fDebugMode){
            std::cout << "[SystToolsEventWeight::produce]       pid of this resp = " << (*it).pid << "\n"
                      << "[SystToolsEventWeight::produce]       prettyName of this resp = " << prettyName << "\n"
                      << "[SystToolsEventWeight::produce]       paramVariations.size() of this resp (expected) = " << sph.paramVariations.size() << "\n"
                      << "[SystToolsEventWeight::produce]       responses.size() of this resp (produced) = " << (*it).responses.size() << std::endl;
          }
          if(sph.paramVariations.size()!=(*it).responses.size()){
            std::cerr << "[SystToolsEventWeight::produce]       prettyName of this resp = " << prettyName << std::endl;
            std::cerr << "[SystToolsEventWeight::produce]       We expect to have " << sph.paramVariations.size() << " universes, but "
                      << (*it).responses.size() << " are produced. "
                      << "Probably this particular event is not relevant to this systematic variation." << std::endl;
            throw std::exception();
          }

          std::vector<float> wgts;
          for( unsigned int i = 0; i < (*it).responses.size(); ++i ){

            if(fDebugMode) std::cout << "[SystToolsEventWeight::produce]         (*it).responses[i] = " << (*it).responses[i] << std::endl;
            wgts.push_back( (*it).responses[i] );

          } // END loop over universes of this knob

          mcwgh[sp->GetFullyQualifiedName()+"_"+prettyName] = wgts;

        } // END loop over knobs

        mcwghvec->push_back(std::move(mcwgh));
        art::Ptr<sbn::evwgh::EventWeightMap> wghPtr = makeWeightPtr(mcwghvec->size() - 1);
        wghassns->addSingle(mclist.at(nMCTruthIndex), wghPtr);

        nMCTruthIndex++;

      } // END loop over MCTruth particles

    } // END systprovider loop

  }

  e.put(std::move(mcwghvec));
  e.put(std::move(wghassns));

}


void SystToolsEventWeight::beginRun(art::Run& run) {

  auto p = std::make_unique<std::vector<sbn::evwgh::EventWeightParameterSet> >();
  for( auto &sp : fSystProviders ) {
    MF_LOG_INFO("SystToolsEventWeight") << "sp->GetToolType() = " << sp->GetToolType() << "\n"
                                  << "sp->GetFullyQualifiedName() = " << sp->GetFullyQualifiedName() << "\n"
                                  << "sp->GetInstanceName() = " << sp->GetInstanceName() << "\n"
                                  << "Printing each SystParamHeader of this ISystProviderTool..";
    //==== Note: typedef std::vector<SystParamHeader> SystMetaData;
    auto smd = sp->GetSystMetaData();
    for( auto &sph : smd ){
      std::cout << "  sph.prettyName = " << sph.prettyName << std::endl;

      std::string rwmode = "multisigma";
      if(sph.isRandomlyThrown) rwmode = "multisim";

      std::cout << "  rwmode = " << rwmode << std::endl;

      std::vector<float> withds;
      std::cout << "    paramVariation = ";
      for(auto pV : sph.paramVariations){
        std::cout << pV << " ";
        withds.push_back(pV);
      }
      std::cout << std::endl;

      sbn::evwgh::EventWeightParameterSet fParameterSet;
      fParameterSet.AddParameter(sph.prettyName, withds);
      fParameterSet.Configure(sp->GetFullyQualifiedName()+"_"+sph.prettyName, rwmode, sph.paramVariations.size());
      fParameterSet.FillKnobValues();

      p->push_back( std::move(fParameterSet) );

    }
  }

  run.put(std::move(p));

}

  }  // namespace evwgh
}  // namespace sbn

DEFINE_ART_MODULE(sbn::evwgh::SystToolsEventWeight)

