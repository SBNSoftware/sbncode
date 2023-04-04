////////////////////////////////////////////////////////////////////////
// Class: SystToolsEventWeight
// Module Type: producer
// File: SystToolsEventWeight_module.cc
//
// Generated at Mon Feb 21 09:36:11 2015 by Jaesung Kim (jae.sung.kim.3426@gmail.com)
//
////////////////////////////////////////////////////////////////////////

#include <iostream>

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


namespace sbn::evwgh {

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

  fGenieModuleLabel = p.get<std::string>("GeneratorModuleLabel", "generator");
  fAllowMissingTruth = p.get<bool>("AllowMissingTruth");
  fDebugMode = p.get<bool>("DebugMode", false);

  fhicl::ParameterSet syst_provider_config = p.get<fhicl::ParameterSet>("generated_systematic_provider_configuration");

  MF_LOG_DEBUG("SystToolsEventWeight") << "Configuring ISystProvider";

  fSystProviders = systtools::ConfigureISystProvidersFromParameterHeaders(syst_provider_config);
  fParamHeaderMap = systtools::BuildParameterHeaders(fSystProviders);
  fParamHeaderHelper = systtools::ParamHeaderHelper(fParamHeaderMap);

  MF_LOG_DEBUG("SystToolsEventWeight") << "ISystProvider is configuered";

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
  if(!fGenieModuleLabel.empty()) mcTruthHandle = e.getHandle<std::vector<simb::MCTruth>>(fGenieModuleLabel);

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
        throw cet::exception{ "SystToolsEventWeight" }
          << "Got nullptr systtools::EventResponse from provider "
          << sp->GetFullyQualifiedName() << "\n";
      }

      //==== syst_resp->size() is (Number of MCTruth)
      //==== Each index corresponds to each of MCTruth
      int nMCTruthIndex = 0;
      if(fDebugMode) std::cout << "[SystToolsEventWeight::produce]   syst_resp.size() (= Number of MCTruth) of this SystProvider = " << syst_resp->size() << std::endl;

      //==== Looping over syst_resp is identical to looping over MCTruth
      for(systtools::event_unit_response_t const& resp: *syst_resp) {
        //==== resp.size() corresponds to number of knobs we altered;
        //==== e.g., MaCCQE, MaCCRES, MvCCRE -> resp.size() = 3
        if(fDebugMode){
          std::cout << "[SystToolsEventWeight::produce]     sp->GetSystMetaData().size() (expected) = " << sp->GetSystMetaData().size() << "\n"
                    << "[SystToolsEventWeight::produce]     resp.size() of this syst_resp (produced) = " << resp.size() << std::endl;
        }
        if(sp->GetSystMetaData().size()!=resp.size()){
          throw cet::exception{ "SystToolsEventWeight" } 
            << "sp->GetFullyQualifiedName() = " << sp->GetFullyQualifiedName() << "\n"
            << "We expect to have " << sp->GetSystMetaData().size() << " knobs for this SystProvider, but "
                    << resp.size() << " are produced. "
                    << "Probably this particular event is not relevant to this systematic variation.\n";
        }

        sbn::evwgh::EventWeightMap mcwgh;

        for( auto const& r: resp){

          //==== responses.size() is the number of universes
          systtools::SystParamHeader const &sph = fParamHeaderHelper.GetHeader( r.pid );
          std::string prettyName = sph.prettyName;

          if(fDebugMode){
            std::cout << "[SystToolsEventWeight::produce]       pid of this resp = " << r.pid << "\n"
                      << "[SystToolsEventWeight::produce]       prettyName of this resp = " << prettyName << "\n"
                      << "[SystToolsEventWeight::produce]       paramVariations.size() of this resp (expected) = " << sph.paramVariations.size() << "\n"
                      << "[SystToolsEventWeight::produce]       responses.size() of this resp (produced) = " << r.responses.size() << std::endl;
          }
          if(sph.paramVariations.size()!=r.responses.size()){
            throw cet::exception{ "SystToolsEventWeight" }
              << "prettyName of this resp = " << prettyName << "\n"
              << "We expect to have " << sph.paramVariations.size() << " universes, but "
                      << r.responses.size() << " are produced. "
                      << "Probably this particular event is not relevant to this systematic variation.\n";
          }

          //==== r.responses : std::vector<double>
          //==== std::map<std::string, std::vector<float> > EventWeightMap
          mcwgh[sp->GetFullyQualifiedName()+"_"+prettyName].assign(r.responses.cbegin(), r.responses.cend());

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
    auto const& smd = sp->GetSystMetaData();
    for( auto &sph : smd ){
      if (fDebugMode) {
        std::cout << "  sph.prettyName = " << sph.prettyName << std::endl;
      }

      std::string rwmode = "multisigma";
      if(sph.isRandomlyThrown) rwmode = "multisim";

      if (fDebugMode) {
        std::cout << "  rwmode = " << rwmode << std::endl;
      }

      std::vector<float> widths { sph.paramVariations.begin(), sph.paramVariations.end() };

      sbn::evwgh::EventWeightParameterSet fParameterSet;
      fParameterSet.AddParameter(sph.prettyName, std::move(widths));
      fParameterSet.Configure(sp->GetFullyQualifiedName()+"_"+sph.prettyName, rwmode, sph.paramVariations.size());
      fParameterSet.FillKnobValues();

      p->push_back( std::move(fParameterSet) );

    }
  }

  run.put(std::move(p), art::fullRun());

}

}  // namespace sbn::evwgh

DEFINE_ART_MODULE(sbn::evwgh::SystToolsEventWeight)
