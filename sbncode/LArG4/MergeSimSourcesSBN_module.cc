////////////////////////////////////////////////////////////////////////
// Class:       MergeSimSourcesSBN
// Plugin Type: producer (Unknown Unknown)
// File:        MergeSimSourcesSBN_module.cc
//
// Generated at Thu Mar 23 00:51:27 2023 by Laura Mai Hien Domine using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "larsim/MergeSimSources/MergeSimSources.h"
#include "MergeSimSourcesLiteUtility.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"
#include "lardataobj/MCBase/MCParticleLite.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
namespace sbn {
  class MergeSimSourcesSBN;
}


class sbn::MergeSimSourcesSBN : public art::EDProducer {
public:
  struct Config {

    fhicl::Sequence<art::InputTag> InputSourcesLabels{
      fhicl::Name{"InputSourcesLabels"},
      fhicl::Comment{"label of the LArG4 products to merge"}};

    fhicl::Sequence<int> TrackIDOffsets{
      fhicl::Name{"TrackIDOffsets"},
      fhicl::Comment{"offset to add to the MC particles from each source"}};

    fhicl::Atom<bool> StoreReflected{
      fhicl::Name{"StoreReflected"},
      fhicl::Comment{"whether to merge also photons from reflections"},
      false // default
    };

    fhicl::Atom<bool> FillMCParticles{
      fhicl::Name{"FillMCParticles"},
      fhicl::Comment{"whether to merge MCParticle"},
      true // default
    };

    fhicl::Atom<bool> FillMCParticlesLite{
      fhicl::Name{"FillMCParticlesLite"},
      fhicl::Comment{"whether to merge MCParticleLite"},
      true // default
    };

    fhicl::Atom<bool> FillMCParticlesAssociated{
      fhicl::Name{"FillMCParticlesAssociated"},
      fhicl::Comment{"whether to merge art::Assns<simb::MCParticle, simb::MCTruth, sim::GeneratedParticleInfo>"},
      true // default
    };

    fhicl::Atom<bool> FillSimPhotons{
      fhicl::Name{"FillSimPhotons"},
      fhicl::Comment{"whether to merge SimPhotons"},
      true // default
    };

    fhicl::Atom<bool> FillSimChannels{
      fhicl::Name{"FillSimChannels"},
      fhicl::Comment{"whether to merge SimChannels"},
      true // default
    };

    fhicl::Atom<bool> FillAuxDetSimChannels{
      fhicl::Name{"FillAuxDetSimChannels"},
      fhicl::Comment{"whether to merge AuxDetSimChannels"},
      true // default
    };

    fhicl::OptionalAtom<bool> FillSimEnergyDeposits{
      fhicl::Name{"FillSimEnergyDeposits"},
      fhicl::Comment{"merges energy deposit collection; if omitted, follows LArG4Parmeters"}};

    fhicl::Sequence<std::string> EnergyDepositInstanceLabels{
      fhicl::Name{"EnergyDepositInstanceLabels"},
      fhicl::Comment{"labels of sim::SimEnergyDeposit collections to merge"},
      std::vector<std::string>{"TPCActive", "Other"} // default
    };

    fhicl::Atom<bool> FillAuxDetHits{
      fhicl::Name{"FillAuxDetHits"},
      fhicl::Comment{"whether to merge aux det hit collections"},
      false // default
    };

    fhicl::Sequence<std::string> AuxDetHitsInstanceLabels{
      fhicl::Name{"AuxDetHitsInstanceLabels"},
      fhicl::Comment{"labels of AuxDetHits collections to merge"},
      std::vector<std::string>{"LArG4DetectorServicevolAuxDetSensitiveCRTStripY",
                               "Other"} // default
    };

  }; // struct Config

  using Parameters = art::EDProducer::Table<Config>;
 explicit MergeSimSourcesSBN(Parameters const& params);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MergeSimSourcesSBN(MergeSimSourcesSBN const&) = delete;
  MergeSimSourcesSBN(MergeSimSourcesSBN&&) = delete;
  MergeSimSourcesSBN& operator=(MergeSimSourcesSBN const&) = delete;
  MergeSimSourcesSBN& operator=(MergeSimSourcesSBN&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.


  std::vector<art::InputTag> const fInputSourcesLabels;
  std::vector<int> const fTrackIDOffsets;
  bool const fUseLitePhotons;
  bool const fStoreReflected;
  bool const fFillMCParticles;
  bool const fFillMCParticlesLite;
  bool const fFillMCParticlesAssociated;
  bool const fFillSimPhotons;
  bool const fFillSimChannels;
  bool const fFillAuxDetSimChannels;
  bool const fFillSimEnergyDeposits;
  std::vector<std::string> const fEnergyDepositionInstances;
  bool const fFillAuxDetHits;
  std::vector<std::string> const fAuxDetHitsInstanceLabels;

  static std::string const ReflectedLabel;

  void dumpConfiguration() const;
};

namespace {

  template <typename Optional>
  std::optional<typename Optional::value_type> getOptionalValue(Optional const& parameter)
  {

    using Value_t = typename Optional::value_type;

    if (!parameter.hasValue()) return std::nullopt;

    Value_t value;
    parameter(value);
    return {value};

  } // getOptionalValue(Optional& parameter)

} // local namespace

std::string const sbn::MergeSimSourcesSBN::ReflectedLabel{"Reflected"};

sbn::MergeSimSourcesSBN::MergeSimSourcesSBN(Parameters const& params)
  : EDProducer{params}
  , fInputSourcesLabels(params().InputSourcesLabels())
  , fTrackIDOffsets(params().TrackIDOffsets())
  , fUseLitePhotons(art::ServiceHandle<sim::LArG4Parameters const>()->UseLitePhotons())
  , fStoreReflected(params().StoreReflected())
  , fFillMCParticles(params().FillMCParticles())
  , fFillMCParticlesLite(params().FillMCParticlesLite())
  , fFillMCParticlesAssociated(params().FillMCParticlesAssociated())
  , fFillSimPhotons(params().FillSimPhotons())
  , fFillSimChannels(params().FillSimChannels())
  , fFillAuxDetSimChannels(params().FillAuxDetSimChannels())
  , fFillSimEnergyDeposits(
      getOptionalValue(params().FillSimEnergyDeposits)
        .value_or(art::ServiceHandle<sim::LArG4Parameters const>()->FillSimEnergyDeposits()))
  , fEnergyDepositionInstances(params().EnergyDepositInstanceLabels())
  , fFillAuxDetHits(params().FillAuxDetHits())
  , fAuxDetHitsInstanceLabels(params().AuxDetHitsInstanceLabels())
{

  if (fInputSourcesLabels.size() != fTrackIDOffsets.size()) {
    throw art::Exception(art::errors::Configuration)
      << "Unequal input vector sizes: InputSourcesLabels and TrackIDOffsets.\n";
  }

  for (art::InputTag const& tag : fInputSourcesLabels) {

    if (fFillMCParticles) {
      consumes<std::vector<simb::MCParticle>>(tag);
    }
    if (fFillMCParticlesLite) {
      consumes<std::vector<sim::MCParticleLite>>(tag);
    }

    if (fFillMCParticlesAssociated) {
      consumes<art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>(tag);
    }

    if (fFillSimChannels) { consumes<std::vector<sim::SimChannel>>(tag); }

    if (fFillAuxDetSimChannels) { consumes<std::vector<sim::AuxDetSimChannel>>(tag); }

    if (fFillSimPhotons) {
      if (!fUseLitePhotons)
        consumes<std::vector<sim::SimPhotons>>(tag);
      else
        consumes<std::vector<sim::SimPhotonsLite>>(tag);

      if (fStoreReflected) {
        art::InputTag const reflected_tag{tag.label(), ReflectedLabel};
        if (!fUseLitePhotons)
          consumes<std::vector<sim::SimPhotons>>(reflected_tag);
        else
          consumes<std::vector<sim::SimPhotonsLite>>(reflected_tag);
      }
    }

    if (fFillSimEnergyDeposits) {
      for (std::string const& edep_inst : fEnergyDepositionInstances) {
        art::InputTag const edep_tag{tag.label(), edep_inst};
        consumes<std::vector<sim::SimEnergyDeposit>>(edep_tag);
      }
    } // if fill energy deposits

    if (fFillAuxDetHits) {
      for (std::string const& auxdethit_inst : fAuxDetHitsInstanceLabels) {
        art::InputTag const auxdethit_tag{tag.label(), auxdethit_inst};
        consumes<std::vector<sim::AuxDetHit>>(auxdethit_tag);
      }
    }

  } // for input labels

  if (fFillMCParticles) {
    produces<std::vector<simb::MCParticle>>();
  }
  if (fFillMCParticlesLite) { produces<std::vector<sim::MCParticleLite>>(); }
  if (fFillMCParticlesAssociated) {
    produces<art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>();}
  if (fFillSimChannels) { produces<std::vector<sim::SimChannel>>(); }
  if (fFillAuxDetSimChannels) { produces<std::vector<sim::AuxDetSimChannel>>(); }

  if (fFillSimPhotons) {
    if (!fUseLitePhotons)
      produces<std::vector<sim::SimPhotons>>();
    else
      produces<std::vector<sim::SimPhotonsLite>>();

    if (fStoreReflected) {
      if (!fUseLitePhotons)
        produces<std::vector<sim::SimPhotons>>(ReflectedLabel);
      else
        produces<std::vector<sim::SimPhotonsLite>>(ReflectedLabel);
    }
  }

  if (fFillSimEnergyDeposits) {
    for (std::string const& edep_inst : fEnergyDepositionInstances)
      produces<std::vector<sim::SimEnergyDeposit>>(edep_inst);
  } // if

  if (fFillAuxDetHits) {
    for (std::string const& auxdethit_inst : fAuxDetHitsInstanceLabels)
      produces<std::vector<sim::AuxDetHit>>(auxdethit_inst);
  }

  dumpConfiguration();
}

void sbn::MergeSimSourcesSBN::produce(art::Event& e)
{
  auto partCol = std::make_unique<std::vector<simb::MCParticle>>();
  auto partLiteCol = std::make_unique<std::vector<sim::MCParticleLite>>();
  auto scCol = std::make_unique<std::vector<sim::SimChannel>>();
  auto PhotonCol = std::make_unique<std::vector<sim::SimPhotons>>();
  auto LitePhotonCol = std::make_unique<std::vector<sim::SimPhotonsLite>>();
  auto ReflPhotonCol = std::make_unique<std::vector<sim::SimPhotons>>();
  auto ReflLitePhotonCol = std::make_unique<std::vector<sim::SimPhotonsLite>>();
  auto tpassn =
    std::make_unique<art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>();
  auto adCol = std::make_unique<std::vector<sim::AuxDetSimChannel>>();

  using edeps_t = std::vector<sim::SimEnergyDeposit>;
  std::vector<edeps_t> edepCols;
  if (fFillSimEnergyDeposits) edepCols.resize(fEnergyDepositionInstances.size());

  using aux_det_hits_t = std::vector<sim::AuxDetHit>;
  std::vector<aux_det_hits_t> auxdethitCols;
  if (fFillAuxDetHits) auxdethitCols.resize(fAuxDetHitsInstanceLabels.size());

  sim::MergeSimSourcesUtility MergeUtility{fTrackIDOffsets};
  MergeSimSourcesLiteUtility MergeUtilityLite{fTrackIDOffsets};

  for (auto const& [i_source, input_label] : util::enumerate(fInputSourcesLabels)) {

    auto const input_partCol = e.getValidHandle<std::vector<simb::MCParticle>>(input_label);
    art::PtrMaker<simb::MCParticle> const makePartPtr{e};
    if (fFillMCParticles) {
      MergeUtility.MergeMCParticles(*partCol, *input_partCol, i_source);
    }
    if (fFillMCParticlesAssociated){
      //truth-->particle assoc stuff here
      const std::vector<size_t>& assocVectorPrimitive =
        MergeUtility.GetMCParticleListMap().at(i_source);
      art::FindOneP<simb::MCTruth, sim::GeneratedParticleInfo> mctAssn(
        input_partCol, e, input_label);
      for (auto const i_p : util::counter(mctAssn.size()))
        tpassn->addSingle(
          mctAssn.at(i_p), makePartPtr(assocVectorPrimitive[i_p]), mctAssn.data(i_p).ref());
    }
    if (fFillMCParticlesLite) {
      // MCParticleLite
      //art::PtrMaker<sim::MCParticleLite> const makePartPtrLite{e};
      auto const input_partLiteCol = e.getValidHandle<std::vector<sim::MCParticleLite>>(input_label);
      MergeUtilityLite.MergeMCParticleLites(*partLiteCol, *input_partLiteCol, i_source);
    }

    if (fFillSimChannels) {
      auto const& input_scCol = e.getProduct<std::vector<sim::SimChannel>>(input_label);
      MergeUtility.MergeSimChannels(*scCol, input_scCol, i_source);
    }

    if (fFillAuxDetSimChannels) {
      auto const& input_adCol = e.getProduct<std::vector<sim::AuxDetSimChannel>>(input_label);
      MergeUtility.MergeAuxDetSimChannels(*adCol, input_adCol, i_source);
    }

    if (fFillSimPhotons) {
      if (!fUseLitePhotons) {
        auto const& input_PhotonCol = e.getProduct<std::vector<sim::SimPhotons>>(input_label);
        MergeUtility.MergeSimPhotons(*PhotonCol, input_PhotonCol);
      }
      else {
        auto const& input_LitePhotonCol =
          e.getProduct<std::vector<sim::SimPhotonsLite>>(input_label);
        MergeUtility.MergeSimPhotonsLite(*LitePhotonCol, input_LitePhotonCol);
      }

      if (fStoreReflected) {
        art::InputTag const input_reflected_label{input_label.label(), ReflectedLabel};
        if (!fUseLitePhotons) {
          auto const& input_PhotonCol =
            e.getProduct<std::vector<sim::SimPhotons>>(input_reflected_label);
          MergeUtility.MergeSimPhotons(*ReflPhotonCol, input_PhotonCol);
        }
        else {
          auto const& input_LitePhotonCol =
            e.getProduct<std::vector<sim::SimPhotonsLite>>(input_reflected_label);
          MergeUtility.MergeSimPhotonsLite(*ReflLitePhotonCol, input_LitePhotonCol);
        }
      }
    }

    if (fFillSimEnergyDeposits) {
      for (auto const& [edep_inst, edepCol] : util::zip(fEnergyDepositionInstances, edepCols)) {
        art::InputTag const edep_tag{input_label.label(), edep_inst};
        MergeUtility.MergeSimEnergyDeposits(edepCol, e.getProduct<edeps_t>(edep_tag), i_source);
      } // for edep
    }   // if fill energy depositions

    if (fFillAuxDetHits) {
      for (auto const& [auxdethit_inst, auxdethitCol] :
           util::zip(fAuxDetHitsInstanceLabels, auxdethitCols)) {
        art::InputTag const auxdethit_tag{input_label.label(), auxdethit_inst};
        MergeUtility.MergeAuxDetHits(
          auxdethitCol, e.getProduct<aux_det_hits_t>(auxdethit_tag), i_source);
      }
    }
  }

  if (fFillMCParticles) {
    e.put(std::move(partCol));
  }
  if (fFillMCParticlesLite) { e.put(std::move(partLiteCol)); }
  if (fFillMCParticlesAssociated) { e.put(std::move(tpassn)); }
  if (fFillSimChannels) { e.put(std::move(scCol)); }
  if (fFillAuxDetSimChannels) { e.put(std::move(adCol)); }
  if (fFillSimPhotons) {
    if (!fUseLitePhotons)
      e.put(std::move(PhotonCol));
    else
      e.put(std::move(LitePhotonCol));
    if (fStoreReflected) {
      if (!fUseLitePhotons)
        e.put(std::move(ReflPhotonCol), ReflectedLabel);
      else
        e.put(std::move(ReflLitePhotonCol), ReflectedLabel);
    }
  }

  if (fFillSimEnergyDeposits) {
    for (auto&& [edep_inst, edepCol] : util::zip(fEnergyDepositionInstances, edepCols)) {
      e.put(std::make_unique<edeps_t>(move(edepCol)), edep_inst);
    } // for
  }   // if fill energy deposits

  if (fFillAuxDetHits) {
    for (auto&& [auxdethit_inst, auxdethitCol] :
         util::zip(fAuxDetHitsInstanceLabels, auxdethitCols)) {
      e.put(std::make_unique<aux_det_hits_t>(move(auxdethitCol)), auxdethit_inst);
    }
  }
}

void sbn::MergeSimSourcesSBN::dumpConfiguration() const
{

  mf::LogInfo log("MergeSimSources");
  log << "Configuration:"
      << "\n - " << fInputSourcesLabels.size() << " input sources:";
  for (auto const& [i_source, tag, offset] :
       util::enumerate(fInputSourcesLabels, fTrackIDOffsets)) {
    log << "\n   [" << i_source << "] '" << tag.encode() << "' (ID offset: " << offset << ")";
  } // for
  if (fFillMCParticles) log << "\n - filling MCParticles";

  if (fFillMCParticlesLite) log << "\n - filling MCParticlesLite";

  if (fFillMCParticlesAssociated) log << "\n - filling MCParticlesAssociated";

  if (fFillSimChannels) log << "\n - filling SimChannels";

  if (fFillAuxDetSimChannels) log << "\n - filling AuxDetSimChannels";

  if (fFillSimPhotons) {
    log << "\n - filling Simulated Photons";
    if (fUseLitePhotons)
      log << "\n   - using photon summary (`SimPhotonsLite`)";
    else
      log << "\n   - using detailed photons (`SimPhotons`)";
    if (fStoreReflected) log << "\n   - also merging reflected light";
  }

  if (fFillSimEnergyDeposits) {
    log << "\n - filling simulated energy deposits (" << fEnergyDepositionInstances.size()
        << " labels:";
    for (std::string const& label : fEnergyDepositionInstances)
      log << " '" << label << "'";
    log << ")";
  }

  if (fFillAuxDetHits) {
    log << "\n - filling auxiliary detector hits (" << fAuxDetHitsInstanceLabels.size()
        << " labels:";
    for (std::string const& label : fAuxDetHitsInstanceLabels)
      log << " '" << label << "'";
    log << ")";
  }

}

DEFINE_ART_MODULE(sbn::MergeSimSourcesSBN)
