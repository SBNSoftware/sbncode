////////////////////////////////////////////////////////////////////////
// Class:       AssocTruthProducer
// Plugin Type: producer (art v2_11_02)
// File:        AssocTruthProducer_module.cc
//
// Generated at Wed Aug  1 12:45:20 2018 by Gray Putnam using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////
//
//
// Stub class for creating information that you want to associate with a truth object.
//
// Currently does nothing. The "in_fv" variable is not really set. Not for use in actual
// analysis, just as a way to toy around with larsoft stuff.

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "DataTypes/AssocTruthInfo.h"

namespace numuselection {
  class AssocTruthProducer;
}


class numuselection::AssocTruthProducer : public art::EDProducer {
public:
  explicit AssocTruthProducer(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  AssocTruthProducer(AssocTruthProducer const &) = delete;
  AssocTruthProducer(AssocTruthProducer &&) = delete;
  AssocTruthProducer & operator = (AssocTruthProducer const &) = delete;
  AssocTruthProducer & operator = (AssocTruthProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  // Declare member data here.
  art::InputTag _truth_tag;
  bool _verbose;

};


numuselection::AssocTruthProducer::AssocTruthProducer(fhicl::ParameterSet const & p) {
  // get truth tag
  _truth_tag = { "generator" };

  // Call appropriate produces<>() functions here.
  produces<std::vector<numuselection::AssocTruthInfo>>(); 
  produces<art::Assns<numuselection::AssocTruthInfo, simb::MCTruth>>();
}

void numuselection::AssocTruthProducer::produce(art::Event & e)
{
  // make the vector of truth info
  auto truth_info = std::make_unique<std::vector<numuselection::AssocTruthInfo>>();

  // and the associations
  auto assns = std::make_unique<art::Assns<numuselection::AssocTruthInfo, simb::MCTruth>>();

  // ptr maker
  art::PtrMaker<numuselection::AssocTruthInfo> makeTruthPtr(e, *this);

  // get truth info
  auto mctruths_handle = \
    e.getValidHandle<std::vector<simb::MCTruth>>(_truth_tag);
  std::vector<art::Ptr<simb::MCTruth>> mctruths;
  art::fill_ptr_vector(mctruths, mctruths_handle);

  for (art::Ptr<simb::MCTruth> const& mc_truth: mctruths) {
    // assoc truth object
    numuselection::AssocTruthInfo truth;

    // get FV 
    double pos[3] = {mc_truth->GetNeutrino().Nu().Vx(), mc_truth->GetNeutrino().Nu().Vy(), mc_truth->GetNeutrino().Nu().Vz()};
    truth.in_FV = true; // stub
    if (_verbose) {
      std::cout << "AssocTruth: Pos: " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
      std::cout << "AssocTruth: In FV: " <<truth.in_FV << std::endl;
    }

    // store truth object
    truth_info->push_back(std::move(truth));
    // get art pointer
    art::Ptr<numuselection::AssocTruthInfo> assoc_truth_Ptr = makeTruthPtr(truth_info->size() - 1);
    // make association
    assns->addSingle(assoc_truth_Ptr, mc_truth);
  } 

  // put stuff in the event
  e.put(std::move(truth_info));
  e.put(std::move(assns));
}

DEFINE_ART_MODULE(numuselection::AssocTruthProducer)

