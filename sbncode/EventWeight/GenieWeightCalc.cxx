// GenieWeightCalc.cxx
//
// Handles event weights for GENIE systematics studies
//
// Updated by Marco Del Tutto on Feb 18 2017

#include "WeightCalcCreator.h"
#include "WeightCalc.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "nutools/NuReweight/art/NuReweight.h" //GENIEReweight.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

namespace sbncode {
namespace evwgh {
  class GenieWeightCalc : public WeightCalc
  {
  public:
    GenieWeightCalc();
    void Configure(fhicl::ParameterSet const& pset);
    std::vector<std::vector<double> > GetWeight(art::Event & e);
    
  private:
    // The reweighting utility class:
    std::vector<rwgt::NuReweight *> reweightVector;

    CLHEP::RandGaussQ *fGaussRandom;
    std::string fGenieModuleLabel;

    // What follows is the list of sereighting parameters present in LArSoft.
    // Parameters with a (*) contains more that one reweighing parameter at the same time. 
    // They can only be modified changing the relative method in $NUTOOLS_DIR/source/NuReweight/GENIEReweight.cxx
      
    enum EReweight {kNCELaxial,        // Axial mass for NC elastic
                    kNCELeta,          // Strange axial form factor for NC elastic
                    kQEMA,             // Axial mass for CC quasi-elastic
                    kQEVec,            // Choice of CCQE vector form factor (sigma = 0 => BBA05; sigma = 1 => Dipole)
                    kCCResAxial,       // Axial mass for CC resonance neutrino production
                    kCCResVector,      // Vector mass for CC resonance neutrino production
                    kResGanged,        // CC Res && NC Res (NOT ACTIVE)
                    kNCResAxial,       // Axial mass for NC resonance neutrino production
                    kNCResVector,      // Vector mass for NC resonance neutrino production
		    kCohMA,            // Axial mass for CC and NC coherent pion production
                    kCohR0,            // Nuclear size param. controlling pi absorption in Rein-Sehgal model
                    kNonResRvp1pi,     // v+p and vbar + n (1 pi) type interactions (*)
                    kNonResRvbarp1pi,  // v+n and vbar + p (1 pi) type interactions (*)
                    kNonResRvp2pi,     // v+p and vbar + n (2 pi) type interactions (*)
		    kNonResRvbarp2pi,  // v+n and vbar + p (2 pi) type interactions (*)
                    kResDecayGamma,    // BR for radiative resonance decay
                    kResDecayEta,      // BR for single-eta resonance decay
                    kResDecayTheta,    // Pion angular distibution in Delta -> pi N (sigma = 0 => isotropic; sigma = 1 => RS)
                    kNC,
                    kDISAth,           // Ath higher twist param in BY model scaling variable xi_w
                    kDISBth,           // Bth higher twist param in BY model scaling variable xi_w
                    kDISCv1u,          // Cv1u u valence GRV98 PDF correction param in BY model
                    kDISCv2u,          // Cv2u u valence GRV98 PDF correction param in BY model
                    kDISnucl,          // NOT IMPLEMENTED IN GENIE
		    kAGKYxF,           // Pion Feynman x for Npi states in AGKY
                    kAGKYpT,           // Pion transverse momentum for Npi states in AGKY
                    kFormZone,         // Hadron Formation Zone
                    kFermiGasModelKf,  // CCQE Pauli Suppression via changes in Fermi level kF
                    kFermiGasModelSf,  // Choice of model (sigma = 0 => FermiGas; sigma = 1 => SF (spectral function))
                    kIntraNukeNmfp,    // Nucleon mean free path (total rescattering probability)
                    kIntraNukeNcex,    // Nucleon charge exchange probability
                    kIntraNukeNel,     // Nucleon elastic reaction probability
                    kIntraNukeNinel,   // Nucleon inelastic reaction probability
                    kIntraNukeNabs,    // Nucleon absorption probability
                    kIntraNukeNpi,     // Nucleon pi-production probability
                    kIntraNukePImfp,   // Pi mean free path (total rescattering probability)
                    kIntraNukePIcex,   // Pi charge exchange probability
                    kIntraNukePIel,    // Pi elastic reaction probability
                    kIntraNukePIinel,  // Pi inelastic reaction probability
                    kIntraNukePIabs,   // Pi absorption probability
                    kIntraNukePIpi,    // Pi pi-production probability
                    kNReWeights};      // ?
    
    DECLARE_WEIGHTCALC(GenieWeightCalc)
  };
  GenieWeightCalc::GenieWeightCalc()
  {
  }

  void GenieWeightCalc::Configure(fhicl::ParameterSet const& p)
  {
    //global config
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");

    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    //calc config
    std::vector<std::string> pars = pset.get< std::vector<std::string> > ("parameter_list");	
    std::vector<float> parsigmas = pset.get< std::vector<float> > ("parameter_sigma");	
    std::string mode             = pset.get<std::string>("mode");

    if (pars.size() != parsigmas.size() )
      throw cet::exception(__FUNCTION__) << GetName()<<"::Bad fcl configuration. parameter_list and parameter_sigma need to have same number of parameters."<<std::endl;

    int number_of_multisims = pset.get< int > ("number_of_multisims");
      
    std::vector<EReweight> erwgh;
    for( auto & s : pars){
      if      (s == "NCELaxial") erwgh.push_back(kNCELaxial);
      else if (s == "NCELeta") erwgh.push_back(kNCELeta);
      else if (s == "QEMA") erwgh.push_back(kQEMA);
      else if (s == "QEVec") erwgh.push_back(kQEVec);
      else if (s == "CCResAxial") erwgh.push_back(kCCResAxial);
      else if (s == "CCResVector") erwgh.push_back(kCCResVector);
      else if (s == "ResGanged") erwgh.push_back(kResGanged);
      else if (s == "NCResAxial") erwgh.push_back(kNCResAxial);
      else if (s == "NCResVector") erwgh.push_back(kNCResVector);
      else if (s == "CohMA") erwgh.push_back(kCohMA);
      else if (s == "CohR0") erwgh.push_back(kCohR0);
      else if (s == "NonResRvp1pi") erwgh.push_back(kNonResRvp1pi);
      else if (s == "NonResRvbarp1pi") erwgh.push_back(kNonResRvbarp1pi);
      else if (s == "NonResRvp2pi") erwgh.push_back(kNonResRvp2pi);
      else if (s == "NonResRvbarp2pi") erwgh.push_back(kNonResRvbarp2pi);
      else if (s == "ResDecayGamma") erwgh.push_back(kResDecayGamma);
      else if (s == "ResDecayEta") erwgh.push_back(kResDecayEta);
      else if (s == "ResDecayTheta") erwgh.push_back(kResDecayTheta);
      else if (s == "NC") erwgh.push_back(kNC);
      else if (s == "DISAth") erwgh.push_back(kDISAth);
      else if (s == "DISBth") erwgh.push_back(kDISBth);
      else if (s == "DISCv1u") erwgh.push_back(kDISCv1u);
      else if (s == "DISCv2u") erwgh.push_back(kDISCv2u);
      else if (s == "DISnucl") erwgh.push_back(kDISnucl);
      else if (s == "AGKYxF") erwgh.push_back(kAGKYxF);
      else if (s == "AGKYpT") erwgh.push_back(kAGKYpT);
      else if (s == "FormZone") erwgh.push_back(kFormZone);
      else if (s == "FermiGasModelKf") erwgh.push_back(kFermiGasModelKf);
      else if (s == "FermiGasModelSf") erwgh.push_back(kFermiGasModelSf);
      else if (s == "IntraNukeNmfp") erwgh.push_back(kIntraNukeNmfp);
      else if (s == "IntraNukeNcex") erwgh.push_back(kIntraNukeNcex);
      else if (s == "IntraNukeNel") erwgh.push_back(kIntraNukeNel);
      else if (s == "IntraNukeNinel") erwgh.push_back(kIntraNukeNinel);
      else if (s == "IntraNukeNabs") erwgh.push_back(kIntraNukeNabs);
      else if (s == "IntraNukeNpi") erwgh.push_back(kIntraNukeNpi);
      else if (s == "IntraNukePImfp") erwgh.push_back(kIntraNukePImfp);
      else if (s == "IntraNukePIcex") erwgh.push_back(kIntraNukePIcex);
      else if (s == "IntraNukePIel") erwgh.push_back(kIntraNukePIel);
      else if (s == "IntraNukePIinel") erwgh.push_back(kIntraNukePIinel);
      else if (s == "IntraNukePIabs") erwgh.push_back(kIntraNukePIabs);
      else if (s == "IntraNukePIpi") erwgh.push_back(kIntraNukePIpi);
      else {
	throw cet::exception(__FUNCTION__) << GetName()<<"::Physical process "<<s<<" you requested is not available to reweight." << std::endl;
      }
    }
      
    //Prepare sigmas
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));

    std::vector<std::vector<float> > reweightingSigmas(erwgh.size());

    if (mode.find("pm1sigma") != std::string::npos ) { 
      number_of_multisims = 2; // only +-1 sigma if pm1sigma is specified
    }
    for (unsigned int i = 0; i < reweightingSigmas.size(); ++i) {
      reweightingSigmas[i].resize(number_of_multisims);
      for (int j = 0; j < number_of_multisims; j ++) {
	if (mode.find("multisim") != std::string::npos )
	  reweightingSigmas[i][j] = parsigmas[i]*fGaussRandom->shoot(&rng->getEngine(GetName()),0.,1.);
        else if (mode.find("pm1sigma") != std::string::npos )
          reweightingSigmas[i][j] = (j == 0 ? 1.: -1.); // j==0 => 1; j==1 => -1 if pm1sigma is specified
	else
	  reweightingSigmas[i][j] = parsigmas[i];
      }
    }

    reweightVector.resize(number_of_multisims);
    
    for (int weight_point = 0; 
	 weight_point < number_of_multisims;
	 weight_point++){
      
      reweightVector[weight_point] = new rwgt::NuReweight;
      
      for (unsigned int i_reweightingKnob=0;i_reweightingKnob<erwgh.size();i_reweightingKnob++) {
	std::cout<<GetName()<<"::Setting up rwgh "<<weight_point<<"\t"<<i_reweightingKnob<<"\t"<<erwgh[i_reweightingKnob]<<std::endl; 

	switch (erwgh[i_reweightingKnob]){

	case kNCELaxial:
          reweightVector[weight_point] -> ReweightNCEL(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kNCELeta:
          reweightVector[weight_point] -> ReweightNCEL(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
            
	case kQEMA:
	  reweightVector[weight_point]
	    -> ReweightQEMA(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kQEVec:
	  reweightVector[weight_point]
	    -> ReweightQEVec(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kResGanged:
	  //reweightVector[weight_point]
	  //  -> ReweightResGanged(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kCCResAxial:
          reweightVector[weight_point] -> ReweightCCRes(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kCCResVector:
          reweightVector[weight_point] -> ReweightCCRes(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;

	case kNCResAxial:
          reweightVector[weight_point] -> ReweightNCRes(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kNCResVector:
          reweightVector[weight_point] -> ReweightNCRes(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
	case kCohMA:
          reweightVector[weight_point] -> ReweightCoh(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kCohR0:
          reweightVector[weight_point] -> ReweightCoh(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
        
	case kNonResRvp1pi:
	  reweightVector[weight_point]
	    -> ReweightNonResRvp1pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvbarp1pi:
	  reweightVector[weight_point]
	    -> ReweightNonResRvbarp1pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvp2pi:
	  reweightVector[weight_point]
	    -> ReweightNonResRvp2pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvbarp2pi:
	  reweightVector[weight_point]
	    -> ReweightNonResRvbarp2pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kResDecayGamma:
          reweightVector[weight_point] -> ReweightResDecay(reweightingSigmas[i_reweightingKnob][weight_point], 0., 0.);
	  break;
        case kResDecayEta:
          reweightVector[weight_point] -> ReweightResDecay(0., reweightingSigmas[i_reweightingKnob][weight_point], 0.);
          break;
        case kResDecayTheta:
          reweightVector[weight_point] -> ReweightResDecay(0., 0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
	case kNC:
	  reweightVector[weight_point]
	    -> ReweightNC(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kDISAth:
          reweightVector[weight_point] -> ReweightDIS(reweightingSigmas[i_reweightingKnob][weight_point], 0., 0., 0.);
	  break;
        case kDISBth:
          reweightVector[weight_point] -> ReweightDIS(0., reweightingSigmas[i_reweightingKnob][weight_point], 0., 0.);
          break;
        case kDISCv1u:
          reweightVector[weight_point] -> ReweightDIS(0., 0., reweightingSigmas[i_reweightingKnob][weight_point], 0.);
          break;
        case kDISCv2u:
          reweightVector[weight_point] -> ReweightDIS(0., 0., 0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
	case kDISnucl:
	  reweightVector[weight_point]
	    -> ReweightDISnucl(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kAGKYxF:
      reweightVector[weight_point] -> ReweightAGKY(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kAGKYpT:
          reweightVector[weight_point] -> ReweightAGKY(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
      
        case kFormZone:
          reweightVector[weight_point] -> ReweightFormZone(reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
        case kFermiGasModelKf:
          reweightVector[weight_point] -> ReweightFGM(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
          break;
        case kFermiGasModelSf:
          reweightVector[weight_point] -> ReweightFGM(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
        case kIntraNukeNmfp:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightMFP_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNcex:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrCEx_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNel:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrElas_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNinel:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrInel_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNabs:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrAbs_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNpi:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrPiProd_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePImfp:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightMFP_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIcex:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrCEx_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIel:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrElas_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIinel:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrInel_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIabs:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrAbs_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIpi:
          reweightVector[weight_point] -> ReweightIntraNuke(rwgt::fReweightFrPiProd_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
	case kNReWeights:
	  break;
	}
      }
   
    } //loop over nWeights
    // Tell all of the reweight drivers to configure themselves:
    std::cout<< GetName()<<"::Setting up "<<reweightVector.size()<<" reweightcalcs"<<std::endl;
    for(auto & driver : reweightVector){
      driver -> Configure();
    }
  }

  std::vector<std::vector<double> > GenieWeightCalc::GetWeight(art::Event & e)
  { 
    //returns a vector of weights for each neutrino interaction in the event

    //get the MC generator information out of the event       
    //these are all handles to mc information.
    art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;  
    art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
    art::Handle< std::vector<simb::GTruth> > gTruthHandle;

    //actually go and get the stuff
    e.getByLabel(fGenieModuleLabel,mcTruthHandle);
    e.getByLabel(fGenieModuleLabel,mcFluxHandle);
    e.getByLabel(fGenieModuleLabel,gTruthHandle);

    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector(mclist, mcTruthHandle);

    std::vector<art::Ptr<simb::GTruth > > glist;
    art::fill_ptr_vector(glist, gTruthHandle);

    //calculate weight(s) here 
    std::vector<std::vector<double> >weight(mclist.size());
    for ( unsigned int inu=0; inu<mclist.size();inu++) {
      weight[inu].resize(reweightVector.size());    
      for (unsigned int i_weight = 0; 
	   i_weight < reweightVector.size(); 
	   i_weight ++){
	weight[inu][i_weight]= reweightVector[i_weight]-> CalcWeight(*mclist[inu],*glist[inu]);
      }
    }
    return weight;

  }
  REGISTER_WEIGHTCALC(GenieWeightCalc)
}
}
