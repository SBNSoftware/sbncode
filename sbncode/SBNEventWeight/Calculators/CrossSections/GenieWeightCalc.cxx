// GenieWeightCalc.cxx
//
// Handles event weights for GENIE systematics studies
//
// Updated by Marco Del Tutto on Feb 18 2017
//
// Ported from uboonecode to larsim on Feb 14 2017
//   by Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
//
// Heavily rewritten on Dec 9 2019
//   by Steven Gardiner <gardiner@fnal.gov>
// Updated for merge into larsim develop branch on Feb 22 2021 by Steven Gardiner
//
// Ported from larsim on Aug. 2021, 
//	by Keng Lin
//-----------------------------------------------------------

// Standard library includes
#include <map>
#include <memory>
#include <set>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbncode/SBNEventWeight/Base/WeightCalc.h"
#include "sbncode/SBNEventWeight/Base/WeightCalcCreator.h"

//#include "CLHEP/Random/RandGaussQ.h"

#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"

//#include "nugen/NuReweight/art/NuReweight.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

//GENIE includes
#include "GENIE/Framework/Conventions/KineVar.h"
#include "GENIE/Framework/EventGen/EventRecord.h"
#include "GENIE/Framework/Interaction/Interaction.h"
#include "GENIE/Framework/Interaction/Kinematics.h"
#include "GENIE/Framework/Messenger/Messenger.h"
#include "GENIE/Framework/Utils/AppInit.h"

#include "GENIE/RwFramework/GSystSet.h"
#include "GENIE/RwFramework/GSyst.h"
#include "GENIE/RwFramework/GReWeight.h"
#include "GENIE/RwCalculators/GReWeightNuXSecNCEL.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCCQE.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCCRES.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCOH.h"
#include "GENIE/RwCalculators/GReWeightNonResonanceBkg.h"
#include "GENIE/RwCalculators/GReWeightFGM.h"
#include "GENIE/RwCalculators/GReWeightDISNuclMod.h"
#include "GENIE/RwCalculators/GReWeightResonanceDecay.h"
#include "GENIE/RwCalculators/GReWeightFZone.h"
#include "GENIE/RwCalculators/GReWeightINuke.h"
#include "GENIE/RwCalculators/GReWeightAGKY.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCCQEaxial.h"
#include "GENIE/RwCalculators/GReWeightNuXSecCCQEvec.h"
#include "GENIE/RwCalculators/GReWeightNuXSecNCRES.h"
#include "GENIE/RwCalculators/GReWeightNuXSecDIS.h"
#include "GENIE/RwCalculators/GReWeightINukeParams.h"
#include "GENIE/RwCalculators/GReWeightNuXSecNC.h"
#include "GENIE/RwCalculators/GReWeightXSecEmpiricalMEC.h"


namespace sbn {
  namespace evwgh {

class GenieWeightCalc : public WeightCalc {
public:
  GenieWeightCalc() : WeightCalc() {}

  void Configure(fhicl::ParameterSet const& pset,
                 CLHEP::HepRandomEngine& engine) override;

  std::vector<float> GetWeight(art::Event& e, size_t inu) override;

private:
//  std::vector<rwgt::NuReweight> rwVector;
  //CHECK replace:
  std::vector< genie::rew::GReWeight > reweightVector;

  std::string fGenieModuleLabel;
  std::string fTuneName;

  //CHECK add these
  std::map< std::string, int > CheckForIncompatibleSystematics(
        const std::vector<genie::rew::GSyst_t>& knob_vec);

  void SetupWeightCalculators(genie::rew::GReWeight& rw,
        const std::map<std::string, int>& modes_to_use);

 
  bool fQuietMode;//quiet mode

  DECLARE_WEIGHTCALC(GenieWeightCalc)
};


void GenieWeightCalc::Configure(fhicl::ParameterSet const& p,
                                CLHEP::HepRandomEngine& engine) {
//  // Global config
//  fGenieModuleLabel = p.get<std::string>("genie_module_label");
//  const fhicl::ParameterSet& pset = p.get<fhicl::ParameterSet>(GetName());
//
//  // Calculator config
//  auto const& pars = pset.get<std::vector<std::string> >("parameter_list");
//  auto const& parsigmas = pset.get<std::vector<float> >("parameter_sigma");
//
//  if (pars.size() != parsigmas.size()) {
//    throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
//      << "parameter_list and parameter_sigma length mismatch."
//      << std::endl;
//  }
//
//  int number_of_multisims = pset.get<int>("number_of_multisims", 1);
//  std::string mode = pset.get<std::string>("mode");
//
//  // Set up parameters
//  fParameterSet.Configure(GetFullName(), mode, number_of_multisims);
//
//  for (size_t i=0; i<pars.size(); i++) {
//    fParameterSet.AddParameter(pars[i], parsigmas[i]);
//  }
//
//  fParameterSet.Sample(engine);
//
//  // Set up GENIE with the currently active tune
//  evgb::SetEventGeneratorListAndTune();
//  fTuneName = evgb::ExpandEnvVar("${GENIE_XSEC_TUNE}");
//
//  // Set up reweighters
//  rwVector.resize(fParameterSet.fNuniverses);
//
//  for (auto const& it : fParameterSet.fParameterMap) {
//    std::string name = it.first.fName;
//    std::cout << GetFullName() << ": Setting up " << name << std::endl;
//
//    for (size_t i=0; i<fParameterSet.fNuniverses; i++) {
//      rwgt::NuReweight& rw = rwVector[i];
//
//      // Axial mass for NC elastic
//      if (name == "NCELaxial")
//        rw.ReweightNCEL(it.second[i], 0);
//      // Strange axial form factor for NC elastic
//      else if (name == "NCELeta")
//        rw.ReweightNCEL(0, it.second[i]);
//      // Axial mass for CC quasi-elastic
//      else if (name == "QEMA")
//        rw.ReweightQEMA(it.second[i]);
//      // Choice of CCQE vector form factor (0 => BBA05, 1 => Dipole)
//      else if (name == "QEVec")
//        rw.ReweightQEVec(it.second[i]);
//      // Axial mass for CC resonance production
//      else if (name == "CCResAxial")
//        rw.ReweightCCRes(it.second[i], 0);
//      // Vector mass for CC resonance production
//      else if (name == "CCResVector")
//        rw.ReweightCCRes(0, it.second[i]);
//      // Axial mass for NC resonance production
//      else if (name == "NCResAxial")
//        rw.ReweightNCRes(it.second[i], 0);
//      // Vector mass for NC resonance production
//      else if (name == "NCResVector")
//        rw.ReweightNCRes(0, it.second[i]);
//      // Axial mass for CC and NC coherent pion production
//      else if (name == "CohMA")
//        rw.ReweightCoh(it.second[i], 0);
//      // Nuclear size param. controlling pi absorption in Rein-Sehgal model
//      else if (name == "CohR0")
//        rw.ReweightCoh(0, it.second[i]);
//      // v+p and vbar + n (1 pi) type interactions (*)
//      else if (name == "NonResRvp1pi")
//        rw.ReweightNonResRvp1pi(it.second[i]);
//      // v+n and vbar + p (1 pi) type interactions (*)
//      else if (name == "NonResRvbarp1pi")
//        rw.ReweightNonResRvbarp1pi(it.second[i]);
//      // v+p and vbar + n (2 pi) type interactions (*)
//      else if (name == "NonResRvp2pi")
//        rw.ReweightNonResRvp2pi(it.second[i]);
//       // v+n and vbar + p (2 pi) type interactions (*)
//      else if (name == "NonResRvbarp2pi")
//        rw.ReweightNonResRvbarp2pi(it.second[i]);
//      // BR for radiative resonance decay
//      else if (name == "ResDecayGamma")
//        rw.ReweightResDecay(it.second[i], 0, 0);
//      // BR for single-eta resonance decay
//      else if (name == "ResDecayEta")
//        rw.ReweightResDecay(0, it.second[i], 0);
//      // Pion angular distibution in Delta -> pi N (0 => isotropic, 1 => RS)
//      else if (name == "ResDecayTheta")
//        rw.ReweightResDecay(0, 0, it.second[i]);
//      // NC normalization
//      else if (name == "NC")
//        rw.ReweightNC(it.second[i]);
//      // Ath higher twist param in BY model scaling variable xi_w
//      else if (name == "DISAth")
//        rw.ReweightDIS(it.second[i], 0, 0, 0);
//      // Bth higher twist param in BY model scaling variable xi_w
//      else if (name == "DISBth")
//        rw.ReweightDIS(0, it.second[i], 0, 0);
//      // Cv1u u valence GRV98 PDF correction param in BY model
//      else if (name == "DISCv1u")
//        rw.ReweightDIS(0, 0, it.second[i], 0);
//      // Cv2u u valence GRV98 PDF correction param in BY model
//      else if (name == "DISCv2u")
//        rw.ReweightDIS(0, 0, 0, it.second[i]);
//      // Not implemented
//      else if (name == "DISnucl")
//        rw.ReweightDISnucl(it.second[i]);
//      // Pion Feynman x for Npi states in AGKY
//      else if (name == "AGKYxF")
//        rw.ReweightAGKY(it.second[i], 0);
//      // Pion transverse momentum for Npi states in AGKY
//      else if (name == "AGKYpT")
//        rw.ReweightAGKY(0, it.second[i]);
//      // Hadron Formation Zone
//      else if (name == "FormZone")
//        rw.ReweightFormZone(it.second[i]);
//      // CCQE Pauli Suppression via changes in Fermi level kF
//      else if (name == "FermiGasModelKf")
//        rw.ReweightFGM(it.second[i], 0);
//      // Choice of model (0 => FermiGas, 1 => SF (spectral function))
//      else if (name == "FermiGasModelSf")
//        rw.ReweightFGM(0, it.second[i]);
//      // Nucleon mean free path (total rescattering probability)
//      else if (name == "IntraNukeNmfp")
//        rw.ReweightIntraNuke(rwgt::fReweightMFP_N, it.second[i]);
//      // Nucleon charge exchange probability
//      else if (name == "IntraNukeNcex")
//        rw.ReweightIntraNuke(rwgt::fReweightFrCEx_N, it.second[i]);
//      // Nucleon inelastic reaction probability
//      else if (name == "IntraNukeNinel")
//        rw.ReweightIntraNuke(rwgt::fReweightFrInel_N, it.second[i]);
//      // Nucleon absorption probability
//      else if (name == "IntraNukeNabs")
//        rw.ReweightIntraNuke(rwgt::fReweightFrAbs_N, it.second[i]);
//      // Nucleon pi-production probability
//      else if (name == "IntraNukeNpi")
//        rw.ReweightIntraNuke(rwgt::fReweightFrPiProd_N, it.second[i]);
//      // Pi mean free path (total rescattering probability)
//      else if (name == "IntraNukePImfp")
//        rw.ReweightIntraNuke(rwgt::fReweightMFP_pi, it.second[i]);
//      // Pi charge exchange probability
//      else if (name == "IntraNukePIcex")
//        rw.ReweightIntraNuke(rwgt::fReweightFrCEx_pi, it.second[i]);
//      // Pi inelastic reaction probability
//      else if (name == "IntraNukePIinel")
//        rw.ReweightIntraNuke(rwgt::fReweightFrInel_pi, it.second[i]);
//      // Pi absorption probability
//      else if (name == "IntraNukePIabs")
//        rw.ReweightIntraNuke(rwgt::fReweightFrAbs_pi, it.second[i]);
//      // Pi pi-production probability
//      else if (name == "IntraNukePIpi")
//        rw.ReweightIntraNuke(rwgt::fReweightFrPiProd_pi, it.second[i]);
//      // Unknown
//      else {
//        throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
//          << "Unknown GENIE parameter " << name << std::endl;
//      }
//    }
//  }
//
//  // Configure reweight drivers
//  for (auto& rw : rwVector) {
//    rw.Configure();
//  }
}

//Keng:
//copied from larsim's GetWeight()
//2d weights --> 1d weights;
std::vector<float> GenieWeightCalc::GetWeight(art::Event& e, size_t inu) {
    // Get the MC generator information out of the event
    // These are both handles to MC information.
    art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
    art::Handle< std::vector<simb::GTruth> > gTruthHandle;

    // Actually go and get the stuff
    e.getByLabel( fGenieModuleLabel, mcTruthHandle );
    e.getByLabel( fGenieModuleLabel, gTruthHandle );

    std::vector< art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector( mclist, mcTruthHandle );

    std::vector< art::Ptr<simb::GTruth > > glist;
    art::fill_ptr_vector( glist, gTruthHandle );

//    size_t num_neutrinos = mclist.size();
    size_t num_knobs = reweightVector.size();

    // Calculate weight(s) here
    std::vector< float > weights( num_knobs );
//    std::vector< std::vector<double> > weights( num_neutrinos );

      // Convert the MCTruth and GTruth objects from the event
      // back into the original genie::EventRecord needed to
      // compute the weights
      std::unique_ptr< genie::EventRecord >
        genie_event( evgb::RetrieveGHEP(*mclist[inu], *glist[inu]) );

      // Set the final lepton kinetic energy and scattering cosine
      // in the owned GENIE kinematics object. This is done during
      // event generation but is not reproduced by evgb::RetrieveGHEP().
      // Several new CCMEC weight calculators developed for MicroBooNE
      // expect the variables to be set in this way (so that differential
      // cross sections can be recomputed). Failing to set them results
      // in inf and NaN weights.
      // TODO: maybe update evgb::RetrieveGHEP to handle this instead.
      genie::Interaction* interaction = genie_event->Summary();
      genie::Kinematics* kine_ptr = interaction->KinePtr();

      // Final lepton mass
      double ml = interaction->FSPrimLepton()->Mass();
      // Final lepton 4-momentum
      const TLorentzVector& p4l = kine_ptr->FSLeptonP4();
      // Final lepton kinetic energy
      double Tl = p4l.E() - ml;
      // Final lepton scattering cosine
      double ctl = p4l.CosTheta();

      kine_ptr->SetKV( kKVTl, Tl );
      kine_ptr->SetKV( kKVctl, ctl );

      // All right, the event record is fully ready. Now ask the GReWeight
      // objects to compute the weights.
      for (size_t k = 0u; k < num_knobs; ++k ) {
        weights[k] = reweightVector.at( k ).CalcWeight( *genie_event );
      }
    
    return weights;


}


std::map< std::string, int > GenieWeightCalc::CheckForIncompatibleSystematics(
		const std::vector<genie::rew::GSyst_t>& knob_vec){
		
    std::map< std::string, int > modes_to_use;
	{
	//CHECK add something

	}
    return modes_to_use;
		}

void GenieWeightCalc::SetupWeightCalculators(genie::rew::GReWeight& rw,
		const std::map<std::string, int>& modes_to_use){
		
		}



REGISTER_WEIGHTCALC(GenieWeightCalc)

  }  // namespace evwgh
}  // namespace sbn



