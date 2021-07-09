// FluxWeightCalc.cxx
// Based on  ubsim / EventWeight / Calculators / FluxUnisimWeightCalc.cxx  @ UBOONE_SUITE_v08_00_00_55 
// Ported to/adapted for SBNCode by Keng Lin July 2021

#include "art/Framework/Principal/Event.h"
#include "nugen/NuReweight/art/NuReweight.h"
#include "nusimdata/SimulationBase/MCFlux.h" 
#include "nusimdata/SimulationBase/MCTruth.h"

#include "sbncode/SBNEventWeight/Base/WeightCalc.h"
#include "sbncode/SBNEventWeight/Base/WeightCalcCreator.h"
#include <sys/stat.h> //for exit(0);

#include "TH1F.h"
#include "TFile.h"

namespace sbn {
	namespace evwgh {

		class FluxWeightCalc : public WeightCalc {
			public:
				FluxWeightCalc() : WeightCalc() {}

				//Read FHiCL and store the settings for the reweighting environment and the calculator.
				void Configure(fhicl::ParameterSet const& pset,
						CLHEP::HepRandomEngine& engine) override;

				//5 types of configurations
//				void ConfigureUnisim(fhicl::ParameterSet const& pset);//Unisim
//				void ConfigurePHFS	(fhicl::ParameterSet const& pset);//PrimaryHadronFeynmanScaling
//				void ConfigurePHSW	(fhicl::ParameterSet const& pset);//PrimaryHadronSanfordWang$
//				void ConfigurePHSWSV(fhicl::ParameterSet const& pset);//PrimaryHadronSWCentralSplineVariation
//				void ConfigurePHN	(fhicl::ParameterSet const& pset);//PrimaryHadronNormalization



				//GetWeight() returns the final weights as a vector
				//	each weight calculation is carried out by WeightCalc() function;
				//inu - the ith parameter.
				std::vector<float> GetWeight(art::Event& e, size_t inu) override;//CHECK, why fluxreader did not go through this? A: fcl setup is incorrect

				//UnisimWeightCalc() - Function for evaluating a specific weight
				//enu - neutrino energy from simb::MCTruth; 
				//ptype-  parent particles label: pi, k, k0, mu from simb:MCFlux
				//ntype - neutrino flavor label: numu, numubar, nue, nuebar from simb:MCFlux
				//randomN - input randmo number
				//noNeg - determine what formulas to use for weights depending on input histograms.
				double UnisimWeightCalc(double enu, int ptype, int ntype, double randomN, bool noNeg);//Unisim

//tool
				std::vector<double> ConvertToVector(TArrayD const* array);

			private:
//fParameterSet was prepared in `sbncode/Base/WeightManager.h`
				std::string fGenieModuleLabel;
				std::string CalcType;

				double fScalePos{}; 
				double fScaleNeg = 1; //for Unisim

//				std::vector<double> fWeightArray{};//a vector of random numbers
				//replaced by std::map<EventWeightParameter, std::vector<float> > fParameterSet.fParameterMap
				//CHECK, cannot handld the cituation that there are two sets of fWerightArray; i.e. from AddParameter()

				//for Unisim
				//Contents of the input histograms;
				//[mu/pi/k-/k] x [nue,anue,numu,anumu] x [bin#]
				double fCV[4][4][200];
				double fRWpos[4][4][200];
				double fRWneg[4][4][200];
				bool PosOnly{false};

				//for HadronsProduction
				std::vector<int> fprimaryHad = {0};
				//-- FeynmanScaling
				std::vector<double> FitVal{};//Shared to SanfordWang
				TMatrixD* FitCov{nullptr};//Shared to SanfordWang

				//-- SWCentaralSplineVariation
//				TMatrixD* HARPCov{nullptr};
				TMatrixD* HARPLowerTriangluarCov{nullptr};//reduced HARPCov matrix to be used
				TMatrixD* HARPXSec{nullptr};
				std::vector<double> HARPmomentumBounds{};
				std::vector<double> HARPthetaBounds{};
				std::vector<double> SWParam{};

				DECLARE_WEIGHTCALC(FluxWeightCalc)
		};

	}  // namespace evwgh
}  // namespace sbn
