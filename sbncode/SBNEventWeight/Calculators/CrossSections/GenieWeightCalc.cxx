// GenieWeightCalc.cxx
//
// Handles event weights for GENIE systematics studies
//
// Updated by Marco Del Tutto on Feb 18 2017
//
// Ported from uboonecode to larsim on Feb 14 2017
//   by Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
//
// Ported to/adapted for SBNCode, Dec 2020, A. Mastbaum

#include "art/Framework/Principal/Event.h"
#include "nugen/NuReweight/art/NuReweight.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"
#include "sbncode/SBNEventWeight/Base/WeightCalc.h"
#include "sbncode/SBNEventWeight/Base/WeightCalcCreator.h"

namespace sbn {
  namespace evwgh {

class GenieWeightCalc : public WeightCalc {
public:
  GenieWeightCalc() : WeightCalc() {}

  void Configure(fhicl::ParameterSet const& pset,
                 CLHEP::HepRandomEngine& engine) override;

  std::vector<float> GetWeight(art::Event& e, size_t inu) override;

private:
  std::vector<rwgt::NuReweight> rwVector;
  std::string fGenieModuleLabel;
  std::string fTuneName;

  DECLARE_WEIGHTCALC(GenieWeightCalc)
};


void GenieWeightCalc::Configure(fhicl::ParameterSet const& p,
                                CLHEP::HepRandomEngine& engine) {
  // Global config
  fGenieModuleLabel = p.get<std::string>("genie_module_label");
  const fhicl::ParameterSet& pset = p.get<fhicl::ParameterSet>(GetName());

  // Calculator config
  auto const& pars = pset.get<std::vector<std::string> >("parameter_list");
  auto const& parsigmas = pset.get<std::vector<float> >("parameter_sigma");

  if (pars.size() != parsigmas.size()) {
    throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
      << "parameter_list and parameter_sigma length mismatch."
      << std::endl;
  }

  int number_of_multisims = pset.get<int>("number_of_multisims", 1);
  std::string mode = pset.get<std::string>("mode");

  // Set up parameters
  fParameterSet.Configure(GetFullName(), mode, number_of_multisims);

  for (size_t i=0; i<pars.size(); i++) {
    fParameterSet.AddParameter(pars[i], parsigmas[i]);
  }

  fParameterSet.Sample(engine);

  // Set up GENIE with the currently active tune
  evgb::SetEventGeneratorListAndTune();
  fTuneName = evgb::ExpandEnvVar("${GENIE_XSEC_TUNE}");

  // Set up reweighters
  rwVector.resize(fParameterSet.fNuniverses);

  for (auto const& it : fParameterSet.fParameterMap) {
    std::string name = it.first.fName;
    std::cout << GetFullName() << ": Setting up " << name << std::endl;

    for (size_t i=0; i<fParameterSet.fNuniverses; i++) {
      rwgt::NuReweight& rw = rwVector[i];

      // Axial mass for NC elastic
      if (name == "NCELaxial")
        rw.ReweightNCEL(it.second[i], 0);
      // Strange axial form factor for NC elastic
      else if (name == "NCELeta")
        rw.ReweightNCEL(0, it.second[i]);
      // Axial mass for CC quasi-elastic
      else if (name == "QEMA")
        rw.ReweightQEMA(it.second[i]);
      // Choice of CCQE vector form factor (0 => BBA05, 1 => Dipole)
      else if (name == "QEVec")
        rw.ReweightQEVec(it.second[i]);
      // Axial mass for CC resonance production
      else if (name == "CCResAxial")
        rw.ReweightCCRes(it.second[i], 0);
      // Vector mass for CC resonance production
      else if (name == "CCResVector")
        rw.ReweightCCRes(0, it.second[i]);
      // Axial mass for NC resonance production
      else if (name == "NCResAxial")
        rw.ReweightNCRes(it.second[i], 0);
      // Vector mass for NC resonance production
      else if (name == "NCResVector")
        rw.ReweightNCRes(0, it.second[i]);
      // Axial mass for CC and NC coherent pion production
      else if (name == "CohMA")
        rw.ReweightCoh(it.second[i], 0);
      // Nuclear size param. controlling pi absorption in Rein-Sehgal model
      else if (name == "CohR0")
        rw.ReweightCoh(0, it.second[i]);
      // v+p and vbar + n (1 pi) type interactions (*)
      else if (name == "NonResRvp1pi")
        rw.ReweightNonResRvp1pi(it.second[i]);
      // v+n and vbar + p (1 pi) type interactions (*)
      else if (name == "NonResRvbarp1pi")
        rw.ReweightNonResRvbarp1pi(it.second[i]);
      // v+p and vbar + n (2 pi) type interactions (*)
      else if (name == "NonResRvp2pi")
        rw.ReweightNonResRvp2pi(it.second[i]);
       // v+n and vbar + p (2 pi) type interactions (*)
      else if (name == "NonResRvbarp2pi")
        rw.ReweightNonResRvbarp2pi(it.second[i]);
      // BR for radiative resonance decay
      else if (name == "ResDecayGamma")
        rw.ReweightResDecay(it.second[i], 0, 0);
      // BR for single-eta resonance decay
      else if (name == "ResDecayEta")
        rw.ReweightResDecay(0, it.second[i], 0);
      // Pion angular distibution in Delta -> pi N (0 => isotropic, 1 => RS)
      else if (name == "ResDecayTheta")
        rw.ReweightResDecay(0, 0, it.second[i]);
      // NC normalization
      else if (name == "NC")
        rw.ReweightNC(it.second[i]);
      // Ath higher twist param in BY model scaling variable xi_w
      else if (name == "DISAth")
        rw.ReweightDIS(it.second[i], 0, 0, 0);
      // Bth higher twist param in BY model scaling variable xi_w
      else if (name == "DISBth")
        rw.ReweightDIS(0, it.second[i], 0, 0);
      // Cv1u u valence GRV98 PDF correction param in BY model
      else if (name == "DISCv1u")
        rw.ReweightDIS(0, 0, it.second[i], 0);
      // Cv2u u valence GRV98 PDF correction param in BY model
      else if (name == "DISCv2u")
        rw.ReweightDIS(0, 0, 0, it.second[i]);
      // Not implemented
      else if (name == "DISnucl")
        rw.ReweightDISnucl(it.second[i]);
      // Pion Feynman x for Npi states in AGKY
      else if (name == "AGKYxF")
        rw.ReweightAGKY(it.second[i], 0);
      // Pion transverse momentum for Npi states in AGKY
      else if (name == "AGKYpT")
        rw.ReweightAGKY(0, it.second[i]);
      // Hadron Formation Zone
      else if (name == "FormZone")
        rw.ReweightFormZone(it.second[i]);
      // CCQE Pauli Suppression via changes in Fermi level kF
      else if (name == "FermiGasModelKf")
        rw.ReweightFGM(it.second[i], 0);
      // Choice of model (0 => FermiGas, 1 => SF (spectral function))
      else if (name == "FermiGasModelSf")
        rw.ReweightFGM(0, it.second[i]);
      // Nucleon mean free path (total rescattering probability)
      else if (name == "IntraNukeNmfp")
        rw.ReweightIntraNuke(rwgt::fReweightMFP_N, it.second[i]);
      // Nucleon charge exchange probability
      else if (name == "IntraNukeNcex")
        rw.ReweightIntraNuke(rwgt::fReweightFrCEx_N, it.second[i]);
      // Nucleon inelastic reaction probability
      else if (name == "IntraNukeNinel")
        rw.ReweightIntraNuke(rwgt::fReweightFrInel_N, it.second[i]);
      // Nucleon absorption probability
      else if (name == "IntraNukeNabs")
        rw.ReweightIntraNuke(rwgt::fReweightFrAbs_N, it.second[i]);
      // Nucleon pi-production probability
      else if (name == "IntraNukeNpi")
        rw.ReweightIntraNuke(rwgt::fReweightFrPiProd_N, it.second[i]);
      // Pi mean free path (total rescattering probability)
      else if (name == "IntraNukePImfp")
        rw.ReweightIntraNuke(rwgt::fReweightMFP_pi, it.second[i]);
      // Pi charge exchange probability
      else if (name == "IntraNukePIcex")
        rw.ReweightIntraNuke(rwgt::fReweightFrCEx_pi, it.second[i]);
      // Pi inelastic reaction probability
      else if (name == "IntraNukePIinel")
        rw.ReweightIntraNuke(rwgt::fReweightFrInel_pi, it.second[i]);
      // Pi absorption probability
      else if (name == "IntraNukePIabs")
        rw.ReweightIntraNuke(rwgt::fReweightFrAbs_pi, it.second[i]);
      // Pi pi-production probability
      else if (name == "IntraNukePIpi")
        rw.ReweightIntraNuke(rwgt::fReweightFrPiProd_pi, it.second[i]);
      // Unknown
      else {
        throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
          << "Unknown GENIE parameter " << name << std::endl;
      }
    }
  }

  // Configure reweight drivers
  for (auto& rw : rwVector) {
    rw.Configure();
  }
}


std::vector<float> GenieWeightCalc::GetWeight(art::Event& e, size_t inu) {
  // Get the MC generator information
  art::Handle<std::vector<simb::MCTruth> > mcTruthHandle;
  e.getByLabel(fGenieModuleLabel, mcTruthHandle);
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  art::fill_ptr_vector(mclist, mcTruthHandle);

  art::Handle<std::vector<simb::GTruth> > gTruthHandle;
  e.getByLabel(fGenieModuleLabel, gTruthHandle);
  std::vector<art::Ptr<simb::GTruth > > glist;
  art::fill_ptr_vector(glist, gTruthHandle);

  // Check tune
  const simb::MCGeneratorInfo& genInfo = mclist.at(inu)->GeneratorInfo();
  std::string tune = genInfo.generatorConfig.at("tune");
  assert(tune == fTuneName);  // FIXME cet exception

  // Calculate weights
  std::vector<float> weights(rwVector.size());
  for (size_t i=0; i<rwVector.size(); i++) {
    weights[i] = rwVector[i].CalcWeight(*(mclist.at(inu)), *(glist.at(inu)));
  }

  return weights;
}

REGISTER_WEIGHTCALC(GenieWeightCalc)

  }  // namespace evwgh
}  // namespace sbn



