#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionGenerator.h"
#include "TFile.h"
#include "CAFAna/Analysis/ExpInfo.h"

#include "toysysts.h"

// Random numbers
#include "TRandom3.h"

#include <set>

using namespace ana;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

class SystWeighter
{
public:
  SystWeighter(const std::set<int>& p, unsigned int u)
    : fParams(p), fUniverse(u)
  {
  }

  double operator()(const caf::SRProxy* sr)
  {
    static std::vector<double> ws;

    if(fUniverse == 0){ // dangerous!
    //    double ret = 1;
      ws = std::vector<double>(1000, 1);
      for(auto& it: sr->truth[0].weights){
        //      if(it.second.size() >= 100) ret *= it.second[fUniverse%it.second.size()];
        //        if(it.second.size() > 50){ // what are the small ones?
        for(int i = 0; i < 1000; ++i) ws[i] *= it.second[i%it.second.size()];

      }
    }

    return ws[fUniverse];

    //    std::cout << fUniverse << " " << ret << std::endl;
    //    std::cout << ret << std::endl;
    //    return ret;
    /*
    const auto& ws = sr->truth[0].weights;

    // Learn the layout from the first record. We're assuming it's the same in
    // all files ever...
    if(fgIdxs.empty()){
      std::cout << "Searching for weight indices " << this << std::endl;
      for(unsigned int i = 0; i < ws.size(); ++i){
        const caf::SRWeight_tProxy& w = ws[i];
        if(fParams.count(w.param_idx) > 0 && w.universe == fUniverse){
          fgIdxs.push_back(i);
        }
      }
      if(fgIdxs.empty()){
        std::cout << "No matching syst entries found" << std::endl;
        abort();
      }
      std::cout << "Found " << fgIdxs.size() << " indices" << std::endl;
    }

    double weight = 1;
    for(unsigned int idx: fgIdxs) weight *= ws[idx].weight;
    return weight;
    */
  }
protected:
  std::set<int> fParams;
  unsigned int fUniverse;
  std::vector<unsigned int> fgIdxs;
};
//std::vector<unsigned int> SystWeighter::fgIdxs;

const Var kSystWeight(SystWeighter({1}, 0));

void make_state(const std::string anatype = numuStr)
{
  // use all systs
  std::set<int> systs;
  for(int i = 0; i < 1000; ++i) systs.insert(i);

  std::vector<Var*> systWs(1000);
  for(unsigned int i = 0; i < systWs.size(); ++i) systWs[i] = new Var(SystWeighter(systs/*{1}*/, i));

  Loaders loaders, loaders2;
  if (anatype == numuStr) {
    //const std::string fDir = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_reweight/processed_2.a/";
    //    const std::string fDir = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_12_05_2018/processed_1.tempwgh/";

    // const std::string fDir = "/sbnd/data/users/bckhouse/processed_1.tempwgh/";
    // const std::string fnameBeam = fDir + "output_SBNOsc_NumuSelection_Modern_SBND.root";
    // const std::string fnameBeam2 = fDir + "output_SBNOsc_NumuSelection_Modern_Icarus.root";

    const std::string dir = "/sbnd/data/users/bckhouse/sample_2.1_fitters/";
    const std::string fnameBeam = dir + "output_SBNOsc_NumuSelection_Proposal_SBND.flat.root";
    const std::string fnameBeam2 = dir + "output_SBNOsc_NumuSelection_Proposal_Icarus.flat.root";

    loaders.SetLoaderPath( fnameBeam, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders2.SetLoaderPath( fnameBeam2, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
  }
  else if (anatype == nueStr) {
    //std::cout << "Nue files not working right now!" << std::endl;
    //return;
    //const std::string fDir = "/pnfs/sbn/persistent/users/dbarker/sbnoutput/";
    const std::string fDir = "/sbnd/data/users/dbarker/sbn/selection/";

    //BNB files contain nominal non-swap beam (so numubg, nuebg, NC)
    const std::string fnameBNB = fDir + "output_SBNOsc_NueSelection_Proposal_SBND_Numu.root";
    const std::string fnameBNB2 = fDir + "output_SBNOsc_NueSelection_Proposal_Icarus_Numu.root";

    //Nue instrinsic only (to increase stats)
    const std::string fnameIntrinsic = fDir + "output_SBNOsc_NueSelection_Proposal_SBND_Int.root";
    const std::string fnameIntrinsic2 = fDir + "output_SBNOsc_NueSelection_Proposal_Icarus_Int.root";

    //Swap files are for signal
    const std::string fnameSwap = fDir + "output_SBNOsc_NueSelection_Proposal_SBND_Osc.root";
    const std::string fnameSwap2 = fDir + "output_SBNOsc_NueSelection_Proposal_Icarus_Osc.root";

    loaders.SetLoaderPath( fnameBNB,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders.SetLoaderPath( fnameIntrinsic,  Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders.SetLoaderPath( fnameSwap,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);
    loaders2.SetLoaderPath( fnameBNB2, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders2.SetLoaderPath( fnameIntrinsic2,  Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders2.SetLoaderPath( fnameSwap2,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);

  }
  else {
    std::cout << "Unrecognized analysis - use numu or nue" << std::endl;
    return;
  }

  const double sbndPOT = kPOTnominal;
  const double icarusPOT = kPOTnominal;

  // Calculator (just use no oscillations)
  osc::NoOscillations* calc = new osc::NoOscillations;

  const Var kRecoE([](const caf::SRProxy* sr)
                   {
                     // std::cout << "ENERGY" << std::endl;
                     return sr->reco.reco_energy;
                   });

  const Var kTrueE([](const caf::SRProxy* sr)
                        {
			  return sr->truth[0].neutrino.energy;
			});

  const Var kWeight([](const caf::SRProxy* sr)
		    {
		      return sr->reco.weight;
		    });

  const Var kWeighthack([](const caf::SRProxy* sr)
                        {
			  std::cout << sr->truth[0].neutrino.iscc << " " << sr->truth[0].neutrino.pdg << " " << sr->reco.weight << std::endl; 
			  if (sr->truth[0].neutrino.iscc && sr->truth[0].neutrino.pdg == 12) return 0.8;
			  if (sr->truth[0].neutrino.iscc && sr->truth[0].neutrino.pdg == 14) return 0.0058;
			  if (sr->truth[0].neutrino.isnc) return 0.058;
			  return 1.0;
			});

  const Cut kOneTrue([](const caf::SRProxy* sr)
		     {
		       return (sr->truth.size() == 1);
		     });


  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);
  const HistAxis axTrueEnergy("True energy (GeV)", binsEnergy, kTrueE);

  // List all of the systematics we'll be using
  std::cout << "\nIncluding the following systematics:" << std::endl;
  for(const ISyst* s: allSysts) std::cout << s->ShortName() << "\t\t" << s->LatexName() << std::endl;
  std::cout << "\n" << std::endl;
  std::vector<const ISyst*> noSysts{};

  std::vector<NoExtrapGenerator*> gens(systWs.size());

  for(unsigned int i = 0; i < systWs.size(); ++i){
    //Use true energy, no weights until we get new nue files
    gens[i] = new NoExtrapGenerator(anatype == numuStr ? axEnergy : axEnergy,
                                    kOneTrue,
                                    *systWs[i]*kWeight);
  }

  NoExtrapGenerator nom_gen(axEnergy, kOneTrue, kWeight);

  if (anatype == numuStr) {
    std::cout << "Using reco energy" << std::endl;
  }
  else {
    std::cout << "Using true energy" << std::endl;
  }

  std::vector<PredictionInterp*> preds_nd(systWs.size()), preds_fd(systWs.size());

  PredictionInterp pred_nom_nd(noSysts, calc, nom_gen, loaders);
  PredictionInterp pred_nom_fd(noSysts, calc, nom_gen, loaders2);

  for(unsigned int i = 0; i < systWs.size(); ++i){
    preds_nd[i] = new PredictionInterp(/*anatype == numuStr ? allSysts :*/ noSysts, calc, *gens[i], loaders);
    preds_fd[i] = new PredictionInterp(/*anatype == numuStr ? allSysts :*/ noSysts, calc, *gens[i], loaders2);
  }

  loaders.Go();
  loaders2.Go();

  std::cout << "Creating file " << ("cafe_state_smear_"+anatype+".root").c_str() << std::endl;

  TFile fout(("cafe_state_smear_"+anatype+".root").c_str(), "RECREATE");
  for(unsigned int i = 0; i < systWs.size(); ++i){
    preds_nd[i]->SaveTo(fout.mkdir(TString::Format("pred_nd_%s_%d", anatype.c_str(), i).Data()));
    preds_fd[i]->SaveTo(fout.mkdir(TString::Format("pred_fd_%s_%d", anatype.c_str(), i).Data()));
  }

  pred_nom_nd.SaveTo(fout.mkdir(TString::Format("pred_nd_%s_nom", anatype.c_str()).Data()));
  pred_nom_fd.SaveTo(fout.mkdir(TString::Format("pred_fd_%s_nom", anatype.c_str()).Data()));

}
