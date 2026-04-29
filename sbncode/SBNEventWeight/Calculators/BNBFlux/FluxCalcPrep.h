// Based on  ubsim / EventWeight / Calculators / FluxUnisimWeightCalc.cxx  @ UBOONE_SUITE_v08_00_00_55 
// Ported to/adapted for SBNCode by Keng Lin July 2021

#include "nugen/NuReweight/art/NuReweight.h"
#include "nusimdata/SimulationBase/MCFlux.h" 
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dk2nu/tree/NuChoice.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

#include "sbncode/SBNEventWeight/Base/WeightCalc.h"
#include "art/Framework/Principal/Event.h"
#include "sbncode/SBNEventWeight/Base/WeightCalcCreator.h"
#include "sbncode/SBNEventWeight/Base/SmearingUtils.h"//MultiGaussianSmearing!

//#include <sys/stat.h> //for exit(0); debugging purpose
#include "TH1.h"
#include "TH3D.h"
#include "TH1F.h"
#include "TFile.h"
#include "TDecompChol.h"//for Choleskey Decomposition

namespace sbn {
  namespace evwgh {

    class FluxWeightCalc : public WeightCalc {
    public:
    FluxWeightCalc() : WeightCalc() {}

      //Read FHiCL and store the settings for the reweighting environment and the calculator.
      void Configure(fhicl::ParameterSet const& p,
		     CLHEP::HepRandomEngine& engine) override;


      //GetWeight() returns the final weights as a vector
      //  each weight evaluaed by *WeightCalc() function;
      //inu - the ith parameter.
      std::vector<float> GetWeight(art::Event& e, size_t inu) override;

      //UnisimWeightCalc() - Function for evaluating a specific weight
      //enu - neutrino energy from simb::MCTruth; 
      //ptype-  parent particles label: pi, k, k0, mu from simb:MCFlux
      //ntype - neutrino flavor label: numu, numubar, nue, nuebar from simb:MCFlux
      //randomN - input randmo number
      //noNeg - determine what formulas to use for weights depending on input histograms.

      //4 *WeightCalc() functions
      double UnisimWeightCalc(double enu, int ptype, int ntype, double randomN, bool noNeg);//Unisim
      std::pair<bool, double> PHNWeightCalc  (simb::MCFlux flux, float rand);//PrimaryHadronNormalizationWeightCalc
      std::pair<bool, double> PHFSWeightCalc  (simb::MCFlux flux, std::vector<float> rand);//PrimaryHadronFeynmanScaling
      std::pair<bool, double> PHSWWeightCalc  (simb::MCFlux flux, std::vector<float> rand);//PrimaryHadronSanfordWangWeightCalc
      std::pair<bool, double> PHSWCSVWeightCalc(simb::MCFlux flux, std::vector<float> rand);//PrimaryHadronSWCentralSplineVariationWeightCalc
        
      //Handy tool
      std::vector<double> ConvertToVector(TArrayD const* array);

    private:
      //fParameterSet was prepared in `sbncode/Base/WeightManager.h`
      std::string fGeneratorModuleLabel;
      std::string CalcType;

      double fScalePos{}; 
      double fScaleNeg = 1; //for Unisim

      //for Unisim
      //load histograms values: [mu/pi/k-/k] x [nue,anue,numu,anumu] x [bin#]
      double fCV[4][4][200];
      double fRWpos[4][4][200];
      double fRWneg[4][4][200];
      bool PosOnly{false};

      //for FluxHist
      bool fUseFluxHist = false;
      //double fCVHist[7][4][200];
      //double fRWHist[7][4][200];

      //double fCVHist[14][200]{};
      //double fRWHist[14][200]{};

      //TH1* fCVHist[14][3]{};
      //TH1* fRWHist[14][3]{};

      TH3D* fCVHist3D[16]{};
      TH3D* fRWHist3D[16]{};

      int GetFluxHistIndex(int parent_pdg, int nu_pdg) const {
	if      (parent_pdg ==  211 && nu_pdg ==  14) return 0;
	else if (parent_pdg ==  211 && nu_pdg ==  12) return 1;

	else if (parent_pdg == -211 && nu_pdg == -14) return 2;
	else if (parent_pdg == -211 && nu_pdg == -12) return 3;

	else if (parent_pdg ==  321 && nu_pdg ==  14) return 4;
	else if (parent_pdg ==  321 && nu_pdg ==  12) return 5;

	else if (parent_pdg == -321 && nu_pdg == -14) return 6;
	else if (parent_pdg == -321 && nu_pdg == -12) return 7;

	else if (parent_pdg ==  130 && nu_pdg ==  14) return 8;
	else if (parent_pdg ==  130 && nu_pdg ==  12) return 9;
	else if (parent_pdg ==  130 && nu_pdg == -14) return 10;
	else if (parent_pdg ==  130 && nu_pdg == -12) return 11;

	else if (parent_pdg ==   13 && nu_pdg ==  14) return 12;
	else if (parent_pdg ==   13 && nu_pdg == -12) return 13;

	else if (parent_pdg ==  -13 && nu_pdg == -14) return 14;
	else if (parent_pdg ==  -13 && nu_pdg ==  12) return 15;

	else return -1;
      }      

      double FluxHistWeightCalc(int parent_pdg, int nu_pdg,
				double px, double py, double pz);

      //for HadronsProduction
      std::vector<int> fprimaryHad = {0};
      //-- FeynmanScaling
      std::vector<double> FitVal{};//Shared to SanfordWang
      TMatrixD* FitCov{nullptr};//Shared to SanfordWang

      //-- SWCentaralSplineVariation
      TMatrixD* HARPXSec{nullptr};
      std::vector<double> HARPmomentumBounds{};
      std::vector<double> HARPthetaBounds{};
      std::vector<double> SWParam{};
      bool fIsDecomposed{false};

      //Weight Counter
      int wc = 0;
      int wcn = 0;
      int wc0 = 0;
      int wc1 = 0;
      int wc30 = 0;

      DECLARE_WEIGHTCALC(FluxWeightCalc)
	};

  }  // namespace evwgh
}  // namespace sbn
