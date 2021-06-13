// FluxUnisimWeightCalc.cxx
//
// Handles event weights for GENIE systematics studies
//
// Updated by Marco Del Tutto on Feb 18 2017
//
// Ported from uboonecode to larsim on Feb 14 2017
//   by Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
//
// Ported to/adapted for SBNCode, Dec 2020, A. Mastbaum
//
// Ported from ubcode; the structure is based on SBNEventweight/Calculator/CrossSections/GenieWeightCalc.cxx
//		by K. Lin Jun. 2021


//GetWeight() function is introduced but it has bug;

#include "art/Framework/Principal/Event.h"
#include "nugen/NuReweight/art/NuReweight.h"

#include "nusimdata/SimulationBase/MCFlux.h" //new for flux
#include "nusimdata/SimulationBase/MCTruth.h"
//#include "nusimdata/SimulationBase/GTruth.h"
//#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"

#include "sbncode/SBNEventWeight/Base/WeightCalc.h"
#include "sbncode/SBNEventWeight/Base/WeightCalcCreator.h"


#include "TH1F.h"//need these to read histograms
#include "TFile.h"

namespace sbn {
  namespace evwgh {

class FluxUnisimWeightCalc : public WeightCalc {
public:
  FluxUnisimWeightCalc() : WeightCalc() {}

  //Read FHiCL and store the settings for the reweighting environment and the calculator.
  void Configure(fhicl::ParameterSet const& pset,
                 CLHEP::HepRandomEngine& engine) override;

  //GetWeight() returns the final weights as a vector
  //	each weight calculation is carried out by WeightCalc() function;
  std::vector<float> GetWeight(art::Event& e, size_t inu) override;

	//Function for evaluating a specific weight
	//enu - neutrino energy from simb::MCTruth; 
	//ptype- charge particles label; 
	//ntype - neutrino label, 
	//noNeg - formulas for weights depending on input histograms.
  double MicroBooNEWeightCalc(double enu, int ptype, int ntype, int uni, bool noNeg);

private:
  //std::vector<rwgt::NuReweight> rwVector;//reweighter? Delete this, if this is for Genie weights
  
  std::string fGenieModuleLabel;
//  std::string fTuneName;

	//Below are from ubcode; some can actually move to the function that needs it?
  std::string fMode;//should be under fParameterSet.fRWType, tb fixed CHECK
  double fScalePos{}; 
  double fScaleNeg{}; 
  std::vector<double> fWeightArray{};//a vector of random numbers

  //contents of the input histograms;
  //[mu/pi/k-/k] x [nue,anue,numu,anumu] x [bin#]
  double fCV[4][4][200];
  double fRWpos[4][4][200];
  double fRWneg[4][4][200];
  bool PosOnly{false};

  DECLARE_WEIGHTCALC(FluxUnisimWeightCalc)
};



// Implement the above functions.

//Configure() need to:
//- Configure reweighting environment (should be done somewhere else):
//	- Set fMode;
//- Configure calculator setting: 
//	- Load histograms and save them;

void FluxUnisimWeightCalc::Configure(fhicl::ParameterSet const& p,
                                CLHEP::HepRandomEngine& engine) {
  // Global config
  // Calculator config
	//horncurrent: {
	//	<<  Reweighting  Environment >>
	//type: FluxUnisim
	//		  parameter_list: ["horncurrent"]
	//		  random_seed: 7
	//		  mode: multisim
	//		  number_of_multisims: 1000
	//
	//	<<  Calculator Setting >>
	//		  CentralValue_hist_file: "*.root"
	//		  PositiveSystematicVariation_hist_file: "*.root"
	//		  NegativeSystematicVariation_hist_file: ".root"
	//		  scale_factor_pos: 1
	//		  scale_factor_neg: 1
	//		  weight_calculator: "MicroBooNE" # "MicroBooNE" OR "MiniBooNE"
	//	<< Outdated >>
	//		  use_MiniBooNE_random_numbers: false
	//			 }
  

  fGenieModuleLabel = p.get<std::string>("genie_module_label");//use this label to get MC*Handle
  const fhicl::ParameterSet& pset = p.get<fhicl::ParameterSet>(GetName());

	//Collect FHiCL parameters:
  auto const& pars = pset.get<std::vector<std::string> >("parameter_list");
//  auto const& parsigmas = pset.get<std::vector<float> >("parameter_sigma");

  //ErrorMsg
//  if (pars.size() != parsigmas.size()) {
//    throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
//      << "parameter_list and parameter_sigma length mismatch."
//      << std::endl;
//  }

  int number_of_multisims = pset.get<int>("number_of_multisims", 1);
  std::string mode = pset.get<std::string>("mode");//3 types: multisim/pmNsigma/fixed

  // Set up parameters
  fParameterSet.Configure(GetFullName(), mode, number_of_multisims);

  for (size_t i=0; i<pars.size(); i++) {
	//Check, no sigma is used in this script.
    fParameterSet.AddParameter(pars[i], 1);//parsigmas[i]);
  }

  fParameterSet.Sample(engine);

  // Set up GENIE with the currently active tune
//  evgb::SetEventGeneratorListAndTune();
//  fTuneName = evgb::ExpandEnvVar("${GENIE_XSEC_TUNE}");

  // Set up reweighters
  fWeightArray.resize(fParameterSet.fNuniverses);

//  for (auto const& it : fParameterSet.fParameterMap) {
//    std::string name = it.first.fName;
//    std::cout << GetFullName() << ": Setting up " << name << std::endl;
//
//    for (size_t i=0; i<fParameterSet.fNuniverses; i++) {
//      rwgt::NuReweight& rw = fWeightArray[i];
//
//      // Axial mass for NC elastic
//      if (name == "NCELaxial")
//        rw.ReweightNCEL(it.second[i], 0);
//      // Strange axial form factor for NC elastic
//      else if (name == "NCELeta")
//        rw.ReweightNCEL(0, it.second[i]);
//      // Axial mass for CC quasi-elastic
//      else {
//        throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
//          << "Unknown GENIE parameter " << name << std::endl;
//      }
//    }
//  }//next weights?

  // Configure reweight drivers
//  for (auto& rw : fWeightArray) {
//    rw.Configure();
//  }
}

//K: To be updated version requires:
// GetWeight() needs:
// - art::Handle< std::vector<simb::MCFlux> > 
// - art::Handle< std::vector<simb::MCTruth> >
// - boolean (true - identical in histograms?)
//
// == Note that inu would be used in mclist.at(inu) in the following;
// [ ] mclist not working;
// [ ] fluxlist not working;

std::vector<float> FluxUnisimWeightCalc::GetWeight(art::Event& e, size_t inu) {
	std::vector<float> weights(fWeightArray.size(), 1);

//--- Copy over from ubcode
//
			//Collect the event's Flux information
			//     This specifically deals with the neutrino type and parentage
//CHECK - unable to get parent particle pdg;
//CMakeFiles/sbncode_SBNEventWeight_Calculators_BNBFlux.dir/FluxUnisimWeightCalc.cxx.o: In function `simb::MCParticle::Momentum(int) const':
///cvmfs/larsoft.opensciencegrid.org/products/nusimdata/v1_22_02/include/nusimdata/SimulationBase/MCParticle.h:220: undefined reference to `simb::MCTrajectory::Momentum(unsigned long) const'
//
//
//			art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
//			e.getByLabel(fGenieModuleLabel,mcFluxHandle);
//			std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;


			//Collect event's MC truth information
			//  This specifically deals with the neutrino energy and 
			//  counting how many interactions there are per event 
			//  (neutrino counting is CRITICALLY important for applying the 
			//   correct weights and not ending up with unphysical values)
			art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
			e.getByLabel(fGenieModuleLabel,mctruthHandle);
			std::vector<simb::MCTruth> const& mclist = *mctruthHandle;
			//Create a vector of weights for each neutrino 
//			std::vector< std::vector<double> > weight;
//			weight.resize(mclist.size());

			// No neutrinos in this event
			if(mclist.size() == 0) return weights;

			//Iterate through each neutrino in the event

				//Resize vector to the number of universes you want to generate
//				weight[inu].resize(fNuni);

				//containers for the parent and neutrino type information
				int ptype = std::numeric_limits<int>::max(); 
				int ntype = std::numeric_limits<int>::max();

//DUMMY test, make it works CHECK
//				if (      fluxlist[0].fptype==13) ptype = 0;
				// Discover the neutrino parent type
				//     This contains the neutrino's parentage information
//				if (      fluxlist[inu].fptype==13  || fluxlist[inu].fptype==-13  ) ptype = 0;
//				else if ( fluxlist[inu].fptype==211 || fluxlist[inu].fptype==-211 ) ptype = 1;
//				else if ( fluxlist[inu].fptype==130                               ) ptype = 2;
//				else if ( fluxlist[inu].fptype==321 || fluxlist[inu].fptype==-321 ) ptype = 3;                                    
//				else {
//					throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<fluxlist[0].fptype<< std::endl;
//				}
//
//				// Discover the neutrino type
//				//     This contains the neutrino's flavor information
//				if (      fluxlist[inu].fntype== 12  ) ntype=0;
//				else if ( fluxlist[inu].fntype==-12  ) ntype=1;
//				else if ( fluxlist[inu].fntype== 14  ) ntype=2;
//				else if ( fluxlist[inu].fntype==-14  ) ntype=3;
//				else {
//					throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ntype "<<fluxlist[0].fntype<< std::endl;
//				}

				// Collect neutrino energy
				//CHECK cant get mclist info; same error as fluxlist
				double enu=0;//mclist[inu].GetNeutrino().Nu().E();      
				std::cout<<"CHECK we got enu = "<<enu<<std::endl;
				std::cout<<"CHECK we got ntype = "<<ntype<<std::endl;
				std::cout<<"CHECK we got ptype = "<<ptype<<std::endl;

				//Let's make a weights based on the calculator you have requested 
				std::cout<<__LINE__<<" CHECK multisim? "<<fMode<<std::endl;
				if(1){//continue calculation with multisim
//CHECK				if((fParameterSet.GetType())==("multisim") ){//continue calculation with multisim}
					for (size_t i=0;i<fWeightArray.size();i++) {

//						if(fWeightCalc.find("MicroBooNE") != std::string::npos){
							weights[i]=1;
							//MicroBooNEWeightCalc(enu, ptype, ntype, i, PosOnly);
//						}
//						if(fWeightCalc.find("MiniBooNE") != std::string::npos){
//							weight[inu][i]=MiniBooNEWeightCalc(enu, ptype, ntype, i, PosOnly);
//						}

					}//Iterate through the number of universes      
				}
			
			std::cout<<"\n\n\n The DUMMY VERSION WORKS!********\n\n\n"<<std::endl;

//--- Copy over from ubcode

  return weights;
}

		double FluxUnisimWeightCalc::MicroBooNEWeightCalc(double enu, int ptype, int ntype, int uni, bool noNeg)
		{

			// 
			//  Largely built off the MiniBooNE code 
			//  but this is intended to expand beyond it
			//
			//  Edits from MiniBooNE code:
			//
			//  JZ (6/2017) : Remove max weight, set to max_limit of a double
			//

			double weight = 1;

			int bin = int(enu/0.05); //convert energy (in GeV) into 50 MeV bin


			//  This is based on:
			//    http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/~checkout~/...
			//          miniboone/AnalysisFramework/MultisimMatrix/src/MultisimMatrix_initialise.F?rev=1.18;content-type=text%2Fplain
			//
			//  pseudocode:
			//   Scaled Reweighting = ScaleFactor * Reweighting + ( 1 - ScaleFactor) * Central Value
			//
			double scaled_pos = fScalePos*fRWpos[ptype][ntype][bin] + 
				(1-fScalePos)*fCV[ptype][ntype][bin];

			double scaled_neg = fScaleNeg*fRWneg[ptype][ntype][bin] + 
				(1-fScaleNeg)*fCV[ptype][ntype][bin];

			// This is based on:
			//    http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/~checkout~/...
			//           miniboone/AnalysisFramework/MultisimMatrix/src/MultisimMatrix_getWeight.F?rev=1.41;content-type=text%2Fplain
			//
			//  pseudocode:
			//   Check value of Random Number array [RAND] for this universe [uni] such that:
			//   If RAND[uni] > 0
			//       Weight =  1 + ( RAND[uni] * [ { Scaled Reweighting[pos] / Central Value } - 1 ] )
			//   If RAND[uni] < 0
			//       Weight =  1 - ( RAND[uni] * [ { Scaled Reweighting[neg] / Central Value } - 1 ] )
			//
			//   if there is only one systematic histogram offered sub this:
			//   If RAND[uni] < 0
			//       Weight =  1 - ( RAND[uni] * [ < 2 - { Scaled Reweighting[pos] / Central Value } > - 1 ] )
			//

			if(fWeightArray[uni] > 0){      
				double syst = fWeightArray[uni]*((scaled_pos/fCV[ptype][ntype][bin])-1);
				weight = 1 + (syst);

				if(scaled_pos == 0) weight = 1;

			}
			else if(noNeg == true){      
				double syst = fWeightArray[uni]*( (2 - (scaled_pos/fCV[ptype][ntype][bin])) - 1);           
				weight = 1 - (syst);      

				if(scaled_pos == 0) weight = 1;

			}
			else{
				double syst = fWeightArray[uni]*((scaled_neg/fCV[ptype][ntype][bin])-1);
				weight = 1 - (syst);    

				if(scaled_neg == 0) weight = 1;

			}

			if(fabs(fCV[ptype][ntype][bin]) < 1.e-12) weight = 1;

			if(weight < 0) weight = 1; 
			if(weight > 30) weight = 30; 
			if(!(std::isfinite(weight))){
				std::cout << "UniSim " <<  fParameterSet.fName << " : Failed to get a finite weight" << std::endl;      
				weight = 30;}

			if( (ntype == 0 || ntype == 1) && ptype == 1) weight = 1;
			if( (ntype == 1 || ntype == 3) && ptype == 3) weight = 1;

			return weight;

		}

REGISTER_WEIGHTCALC(FluxUnisimWeightCalc)

  }  // namespace evwgh
}  // namespace sbn



