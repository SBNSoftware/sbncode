// FluxWeightCalc.cxx
// Based on  ubsim / EventWeight / Calculators / FluxUnisimWeightCalc.cxx  @ UBOONE_SUITE_v08_00_00_55 
// Ported to/adapted for SBNCode by Keng Lin July 2021

#include "art/Framework/Principal/Event.h"
#include "nugen/NuReweight/art/NuReweight.h"
#include "nusimdata/SimulationBase/MCFlux.h" 
#include "nusimdata/SimulationBase/MCTruth.h"

#include "sbncode/SBNEventWeight/Base/WeightCalc.h"
#include "sbncode/SBNEventWeight/Base/WeightCalcCreator.h"
#include "sbncode/SBNEventWeight/Base/SmearingUtils.h"//MultiGaussianSmearing!
#include <sys/stat.h> //for exit(0);

#include "TH1F.h"
#include "TFile.h"
#include "TDecompChol.h"//for Choleskey Decomposition

namespace sbn {
	namespace evwgh {

		class FluxWeightCalc : public WeightCalc {
			public:
				FluxWeightCalc() : WeightCalc() {}

				//Read FHiCL and store the settings for the reweighting environment and the calculator.
				void Configure(fhicl::ParameterSet const& pset,
						CLHEP::HepRandomEngine& engine) override;


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
				std::pair<bool, double> PHNWeightCalc(simb::MCFlux flux, double  rand);//PrimaryHadronNormalizationWeightCalc

				std::pair<bool, double> PHFSWeightCalc(simb::MCFlux flux, std::vector<double> rand);//PrimaryHadronFeynmanScaling
				std::pair<bool, double> PHSWWeightCalc(simb::MCFlux flux, std::vector<double> rand);//PrimaryHadronSanfordWangWeightCalc
				std::pair<bool, double> PHSWCSVWeightCalc(simb::MCFlux flux, std::vector<double> rand);//PrimaryHadronSWCentralSplineVariationWeightCalc
				//tool
				std::vector<double> ConvertToVector(TArrayD const* array);

			private:
				//fParameterSet was prepared in `sbncode/Base/WeightManager.h`
				std::string fGenieModuleLabel;
				std::string CalcType;

				double fScalePos{}; 
				double fScaleNeg = 1; //for Unisim

				std::vector< std::vector< double > > fWeightArray{};//2d matrix of random numbers
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
				//				TMatrixD* HARPLowerTriangluarCov{nullptr};//reduced HARPCov matrix to be used; use FixCov for this
				TMatrixD* HARPXSec{nullptr};
				std::vector<double> HARPmomentumBounds{};
				std::vector<double> HARPthetaBounds{};
				std::vector<double> SWParam{};
				bool fIsDecomposed{false};

				int wc = 0;
				int wcn = 0;
				int wc0 = 0;
				int wc1 = 0;
				int wc30 = 0;

				DECLARE_WEIGHTCALC(FluxWeightCalc)
		};

	}  // namespace evwgh
}  // namespace sbn
