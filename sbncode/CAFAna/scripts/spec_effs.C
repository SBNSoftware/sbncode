// ETW May 2018
// Standard script for DUNE spectra
// Input files TensorFlow CVN training from Fall 2018 

#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Progress.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"

using namespace ana;

#include "Utilities/rootlogon.C"

#include "OscLib/func/IOscCalculator.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"

const Var kRecoE_nue = SIMPLEVAR(dune.Ev_reco_nue);
const Var kRecoE_numu = SIMPLEVAR(dune.Ev_reco_numu);
const Var kPIDFD_NUMU = SIMPLEVAR(dune.cvnnumu);
const Var kPIDFD_NUE = SIMPLEVAR(dune.cvnnue);
const Var kPID_MVA_NUMU = SIMPLEVAR(dune.mvanumu);

const Var kvtxx_truth = SIMPLEVAR(dune.vtx_x);
const Var kvtxy_truth = SIMPLEVAR(dune.vtx_y);
const Var kvtxz_truth = SIMPLEVAR(dune.vtx_z);

//All the vars for efficiency studies
const Var kNuMomX = SIMPLEVAR(dune.NuMomX);
const Var kNuMomY = SIMPLEVAR(dune.NuMomY);
const Var kNuMomZ = SIMPLEVAR(dune.NuMomZ);
const Var kEv = SIMPLEVAR(dune.Ev);
const Var kMode = SIMPLEVAR(dune.mode);
const Var kLepMomX = SIMPLEVAR(dune.LepMomX);
const Var kLepMomY = SIMPLEVAR(dune.LepMomY);
const Var kLepMomZ = SIMPLEVAR(dune.LepMomZ);
const Var kLepE = SIMPLEVAR(dune.LepE);
const Var kLepNuAngle = SIMPLEVAR(dune.LepNuAngle);
const Var kQ2 = SIMPLEVAR(dune.Q2);
const Var kW = SIMPLEVAR(dune.W);
//const Var kX = SIMPLEVAR(dune.X);
const Var kY = SIMPLEVAR(dune.Y);
const Var knP = SIMPLEVAR(dune.nP);
const Var knN = SIMPLEVAR(dune.nN);
const Var knipip = SIMPLEVAR(dune.nipip);
const Var knipim = SIMPLEVAR(dune.nipim);
const Var knipi0 = SIMPLEVAR(dune.nipi0);
//const Var knikp = SIMPLEVAR(dune.nikp);
//const Var knikm = SIMPLEVAR(dune.nikm);
//const Var knik0 = SIMPLEVAR(dune.nik0);
//const Var kniem = SIMPLEVAR(dune.niem);
//const Var kniother = SIMPLEVAR(dune.niother);
//const Var knNucleus = SIMPLEVAR(dune.nNucleus);
//const Var knUNKNOWN = SIMPLEVAR(dune.nUNKNOWN);
const Var kRecoLepEnNue = SIMPLEVAR(dune.RecoLepEnNue);
const Var kRecoHadEnNue = SIMPLEVAR(dune.RecoHadEnNue);

const Cut kPassCVN_NUE = kPIDFD_NUE>0.7 && kPIDFD_NUMU<0.5;
const Cut kPassCVN_NUMU = kPIDFD_NUMU>0.5 && kPIDFD_NUE<0.7;

// 125 MeV bins from 0.0 to 8GeV
const HistAxis axis_nue("Reconstructed energy (GeV)",
                    Binning::Simple(64, 0.0, 8.0),
                    kRecoE_nue);

const HistAxis axis_numu("Reconstructed energy (GeV)",
                    Binning::Simple(64, 0.0, 8.0),
                    kRecoE_numu);

// Set up 2D HistAxes for efficiency studies (nue only for now)


const HistAxis axis_nue_ev("Neutrino energy (GeV)",
			   Binning::Simple(64, 0.0, 8.0),
			   kEv,
			   "Reconstructed energy (GeV)",
			   Binning::Simple(64, 0.0, 8.0),
			   kRecoE_nue);
const HistAxis axis_nue_mode("Interaction Mode",
			     Binning::Simple(11, 0., 11.),
			     kMode,
			     "Reconstructed energy (GeV)",
			     Binning::Simple(64, 0.0, 8.0),
			     kRecoE_nue);

const HistAxis axis_nue_truex("True vertex x",
			       Binning::Simple(100, -450., 450.),
			       kvtxx_truth,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0., 8.0),
			       kRecoE_nue);
const HistAxis axis_nue_truey("True vertex y",
			       Binning::Simple(100, -650., 650.),
			       kvtxy_truth,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0., 8.0),
			       kRecoE_nue);
const HistAxis axis_nue_truez("True vertex z",
			       Binning::Simple(100, -50., 1500.),
			       kvtxz_truth,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0., 8.0),
			       kRecoE_nue);


const HistAxis axis_nue_numomx("Neutrino px",
			       Binning::Simple(40, -0.00001, 0.00001),
			       kNuMomX,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0., 8.0),
			       kRecoE_nue);
const HistAxis axis_nue_numomy("Neutrino py",
			       Binning::Simple(40, 0.0, 1.0),
			       kNuMomY,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0.0, 8.0),
			       kRecoE_nue);
const HistAxis axis_nue_numomz("Neutrino pz",
			       Binning::Simple(40, 0.0, 5.0),
			       kNuMomZ,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0.0, 8.0),
			       kRecoE_nue);
const HistAxis axis_nue_lepmomx("Lepton px",
			       Binning::Simple(40, -1.5, 1.5),
			       kLepMomX,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0.0, 8.0),
			       kRecoE_nue);
const HistAxis axis_nue_lepmomy("Lepton py",
			       Binning::Simple(40, -1.5, 1.5),
			       kLepMomY,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0.0, 8.0),
			       kRecoE_nue);
const HistAxis axis_nue_lepmomz("Lepton pz",
			       Binning::Simple(40, -1.0,5.0),
			       kLepMomZ,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0.0, 8.0),
			       kRecoE_nue);
const HistAxis axis_nue_lepE("Lepton energy",
			       Binning::Simple(40, 0., 5.0),
			       kLepE,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0.0, 8.0),
			       kRecoE_nue);

const HistAxis axis_nue_lepNuAngle("Lepton anlge",
			       Binning::Simple(40, 0., 2.0),
			       kLepNuAngle,
			       "Reconstructed energy (GeV)",
			       Binning::Simple(64, 0.0, 8.0),
			       kRecoE_nue);


const HistAxis axis_nue_q2("Q2",
			   Binning::Simple(40, 0., 5.),
			   kQ2,
			   "Reconstructed energy (GeV)",
			   Binning::Simple(64, 0.0, 8.0),
			   kRecoE_nue);
const HistAxis axis_nue_w("W",
			  Binning::Simple(40, 0., 5.),
			  kW,
			  "Reconstructed energy (GeV)",
			  Binning::Simple(64, 0.0, 8.0),
			  kRecoE_nue);
const HistAxis axis_nue_y("Y",
			  Binning::Simple(40, -0.1, 1.),
			  kY,
			  "Reconstructed energy (GeV)",
			  Binning::Simple(64, 0.0, 8.0),
			  kRecoE_nue);
const HistAxis axis_nue_np("N protons",
			   Binning::Simple(4, 0.0, 4.),
			   knP,
			   "Reconstructed energy (GeV)",
			   Binning::Simple(64, 0.0, 8.0),
			   kRecoE_nue);
const HistAxis axis_nue_nn("N neutrons",
			   Binning::Simple(4, 0.0, 4.),
			   knN,
			   "Reconstructed energy (GeV)",
			   Binning::Simple(64, 0.0, 8.0),
			   kRecoE_nue);
const HistAxis axis_nue_nipip("N pi+",
			      Binning::Simple(4, 0.0, 4.),
			      knipip,
			      "Reconstructed energy (GeV)",
			      Binning::Simple(64, 0.0, 8.0),
			      kRecoE_nue);
const HistAxis axis_nue_nipim("N pi-",
			      Binning::Simple(4, 0.0, 4.),
			      knipim,
			      "Reconstructed energy (GeV)",
			      Binning::Simple(64, 0.0, 8.0),
			      kRecoE_nue);
const HistAxis axis_nue_nipi0("N pi0",
			      Binning::Simple(4, 0.0, 4.),
			      knipi0,
			      "Reconstructed energy (GeV)",
			      Binning::Simple(64, 0.0, 8.0),
			      kRecoE_nue);
const HistAxis axis_nue_recolepen("Reconstruced lepton energy",
			      Binning::Simple(40, 0.0, 5.),
			      kRecoLepEnNue,
			      "Reconstructed energy (GeV)",
			      Binning::Simple(64, 0.0, 8.0),
			      kRecoE_nue);
const HistAxis axis_nue_recohaden("Reconstruced hadron energy",
			      Binning::Simple(40, 0.0, 5.),
			      kRecoHadEnNue,
			      "Reconstructed energy (GeV)",
			      Binning::Simple(64, 0.0, 8.0),
			      kRecoE_nue);



// POT/yr * 3.5yrs * mass correction
//const double potFD = 3.5 * 1.1e21 * 40/1.13;
//Run huge exposure for efficiency studies
const double potFD = 20 * 1.1e21 * 40/1.13;

const char* stateFname = "spec_state_effs.root";
const char* outputFname = "spec_hist_effs.root";

void spec_effs(bool reload = false)
{
  rootlogon(); // style
  
  if(reload || TFile(stateFname).IsZombie()){

    SpectrumLoader loaderFDFHCNonswap("/dune/data/users/marshalc/CAFs/mcc11_test/FD_FHC_nonswap.root");
    SpectrumLoader loaderFDFHCNue("/dune/data/users/marshalc/CAFs/mcc11_test/FD_FHC_nueswap.root");
    SpectrumLoaderBase& loaderFDFHCNuTau = kNullLoader;

    SpectrumLoader loaderFDRHCNonswap("/dune/data/users/marshalc/CAFs/mcc11_test/FD_RHC_nonswap.root");
    SpectrumLoader loaderFDRHCNue("/dune/data/users/marshalc/CAFs/mcc11_test/FD_RHC_nueswap.root");
    SpectrumLoaderBase& loaderFDRHCNuTau = kNullLoader;

    Loaders dummyLoaders; // PredictionGenerator insists on this

    //Selection applied
    PredictionNoExtrap predFDNumuFHC(loaderFDFHCNonswap,
                                     loaderFDFHCNue,
                                     loaderFDFHCNuTau,
                                     axis_numu,
                                     kPassCVN_NUMU && kIsTrueFV);

    PredictionNoExtrap predFDNueFHC(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue,
                                    kPassCVN_NUE && kIsTrueFV);

    PredictionNoExtrap predFDNumuRHC(loaderFDRHCNonswap,
                                     loaderFDRHCNue,
                                     loaderFDRHCNuTau,
                                     axis_numu,
                                     kPassCVN_NUMU && kIsTrueFV);

    PredictionNoExtrap predFDNueRHC(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue,
                                    kPassCVN_NUE && kIsTrueFV);


    //Fiducial Only for Efficiencies
    PredictionNoExtrap predFDNumuFHC_Fid(loaderFDFHCNonswap,
                                     loaderFDFHCNue,
                                     loaderFDFHCNuTau,
                                     axis_numu,
                                     kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_Fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNumuRHC_Fid(loaderFDRHCNonswap,
                                     loaderFDRHCNue,
                                     loaderFDRHCNuTau,
                                     axis_numu,
                                     kIsTrueFV);

    PredictionNoExtrap predFDNueRHC_Fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue,
                                    kIsTrueFV);

    //For Efficiency Studies
    PredictionNoExtrap predFDNueFHC_ev(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_ev,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_ev_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_ev,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_ev(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_ev,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_ev_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_ev,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_mode(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_mode,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_mode_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_mode,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_mode(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_mode,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_mode_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_mode,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_truex(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_truex,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_truex_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_truex,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_truex(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_truex,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_truex_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_truex,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_truey(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_truey,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_truey_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_truey,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_truey(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_truey,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_truey_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_truey,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_truez(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_truez,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_truez_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_truez,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_truez(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_truez,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_truez_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_truez,
                                    kIsTrueFV);


    PredictionNoExtrap predFDNueFHC_numomx(loaderFDFHCNonswap,
					     loaderFDFHCNue,
					     loaderFDFHCNuTau,
					     axis_nue_numomx,
					   kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_numomx_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_numomx,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_numomx(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_numomx,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_numomx_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_numomx,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_numomy(loaderFDFHCNonswap,
					     loaderFDFHCNue,
					     loaderFDFHCNuTau,
					     axis_nue_numomy,
					   kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_numomy_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_numomy,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_numomy(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_numomy,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_numomy_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_numomy,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_numomz(loaderFDFHCNonswap,
					     loaderFDFHCNue,
					     loaderFDFHCNuTau,
					     axis_nue_numomz,
					   kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_numomz_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_numomz,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_numomz(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_numomz,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_numomz_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_numomz,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_lepmomx(loaderFDFHCNonswap,
					     loaderFDFHCNue,
					     loaderFDFHCNuTau,
					     axis_nue_lepmomx,
					   kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_lepmomx_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_lepmomx,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepmomx(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_lepmomx,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepmomx_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_lepmomx,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_lepmomy(loaderFDFHCNonswap,
					     loaderFDFHCNue,
					     loaderFDFHCNuTau,
					     axis_nue_lepmomy,
					   kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_lepmomy_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_lepmomy,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepmomy(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_lepmomy,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepmomy_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_lepmomy,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_lepmomz(loaderFDFHCNonswap,
					     loaderFDFHCNue,
					     loaderFDFHCNuTau,
					     axis_nue_lepmomz,
					   kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_lepmomz_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_lepmomz,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepmomz(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_lepmomz,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepmomz_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_lepmomz,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_lepe(loaderFDFHCNonswap,
					     loaderFDFHCNue,
					     loaderFDFHCNuTau,
					     axis_nue_lepE,
					   kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_lepe_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_lepE,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepe(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_lepE,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepe_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_lepE,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_lepnuangle(loaderFDFHCNonswap,
					       loaderFDFHCNue,
					       loaderFDFHCNuTau,
					       axis_nue_lepNuAngle,
					       kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_lepnuangle_fid(loaderFDFHCNonswap,
						   loaderFDFHCNue,
						   loaderFDFHCNuTau,
						   axis_nue_lepNuAngle,
						   kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepnuangle(loaderFDRHCNonswap,
					       loaderFDRHCNue,
					       loaderFDRHCNuTau,
					       axis_nue_lepNuAngle,
					       kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_lepnuangle_fid(loaderFDRHCNonswap,
						   loaderFDRHCNue,
						   loaderFDRHCNuTau,
						   axis_nue_lepNuAngle,
						   kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_q2(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_q2,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_q2_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_q2,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_q2(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_q2,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_q2_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_q2,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_w(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_w,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_w_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_w,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_w(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_w,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_w_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_w,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_y(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_y,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_y_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_y,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_y(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_y,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_y_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_y,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_np(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_np,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_np_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_np,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_np(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_np,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_np_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_np,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_nn(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_nn,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_nn_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_nn,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_nn(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_nn,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_nn_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_nn,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_nipip(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_nipip,
                                    kPassCVN_NUE && kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_nipip_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_nipip,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_nipip(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_nipip,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_nipip_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_nipip,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_nipim(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_nipim,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_nipim_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_nipim,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_nipim(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_nipim,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_nipim_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_nipim,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_nipi0(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_nipi0,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_nipi0_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_nipi0,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_nipi0(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_nipi0,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_nipi0_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_nipi0,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_recolepen(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_recolepen,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_recolepen_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_recolepen,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_recolepen(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_recolepen,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_recolepen_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_recolepen,
                                    kIsTrueFV);

    PredictionNoExtrap predFDNueFHC_recohaden(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_recohaden,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueFHC_recohaden_fid(loaderFDFHCNonswap,
                                    loaderFDFHCNue,
                                    loaderFDFHCNuTau,
                                    axis_nue_recohaden,
                                    kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_recohaden(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_recohaden,
                                    kPassCVN_NUE && kIsTrueFV);
    PredictionNoExtrap predFDNueRHC_recohaden_fid(loaderFDRHCNonswap,
                                    loaderFDRHCNue,
                                    loaderFDRHCNuTau,
                                    axis_nue_recohaden,
                                    kIsTrueFV);


    loaderFDFHCNonswap.Go();
    loaderFDFHCNue.Go();
    loaderFDRHCNonswap.Go();
    loaderFDRHCNue.Go();

    TFile fout(stateFname, "RECREATE");
    predFDNumuFHC.SaveTo(fout.mkdir("fd_numu_fhc"));
    predFDNueFHC.SaveTo(fout.mkdir("fd_nue_fhc"));
    predFDNumuRHC.SaveTo(fout.mkdir("fd_numu_rhc"));
    predFDNueRHC.SaveTo(fout.mkdir("fd_nue_rhc"));

    predFDNumuFHC_Fid.SaveTo(fout.mkdir("fd_numu_fhc_fid"));
    predFDNueFHC_Fid.SaveTo(fout.mkdir("fd_nue_fhc_fid"));
    predFDNumuRHC_Fid.SaveTo(fout.mkdir("fd_numu_rhc_fid"));
    predFDNueRHC_Fid.SaveTo(fout.mkdir("fd_nue_rhc_fid"));

    predFDNueFHC_ev.SaveTo(fout.mkdir("fd_nue_fhc_ev"));
    predFDNueFHC_ev_fid.SaveTo(fout.mkdir("fd_nue_fhc_ev_fid"));
    predFDNueRHC_ev.SaveTo(fout.mkdir("fd_nue_rhc_ev"));
    predFDNueRHC_ev_fid.SaveTo(fout.mkdir("fd_nue_rhc_ev_fid"));
    predFDNueFHC_numomx.SaveTo(fout.mkdir("fd_nue_fhc_numomx"));
    predFDNueFHC_numomx_fid.SaveTo(fout.mkdir("fd_nue_fhc_numomx_fid"));
    predFDNueRHC_numomx.SaveTo(fout.mkdir("fd_nue_rhc_numomx"));
    predFDNueRHC_numomx_fid.SaveTo(fout.mkdir("fd_nue_rhc_numomx_fid"));
    predFDNueFHC_numomy.SaveTo(fout.mkdir("fd_nue_fhc_numomy"));
    predFDNueFHC_numomy_fid.SaveTo(fout.mkdir("fd_nue_fhc_numomy_fid"));
    predFDNueRHC_numomy.SaveTo(fout.mkdir("fd_nue_rhc_numomy"));
    predFDNueRHC_numomy_fid.SaveTo(fout.mkdir("fd_nue_rhc_numomy_fid"));
    predFDNueFHC_numomz.SaveTo(fout.mkdir("fd_nue_fhc_numomz"));
    predFDNueFHC_numomz_fid.SaveTo(fout.mkdir("fd_nue_fhc_numomz_fid"));
    predFDNueRHC_numomz.SaveTo(fout.mkdir("fd_nue_rhc_numomz"));
    predFDNueRHC_numomz_fid.SaveTo(fout.mkdir("fd_nue_rhc_numomz_fid"));
    predFDNueFHC_lepmomx.SaveTo(fout.mkdir("fd_nue_fhc_lepmomx"));
    predFDNueFHC_lepmomx_fid.SaveTo(fout.mkdir("fd_nue_fhc_lepmomx_fid"));
    predFDNueRHC_lepmomx.SaveTo(fout.mkdir("fd_nue_rhc_lepmomx"));
    predFDNueRHC_lepmomx_fid.SaveTo(fout.mkdir("fd_nue_rhc_lepmomx_fid"));
    predFDNueFHC_lepmomy.SaveTo(fout.mkdir("fd_nue_fhc_lepmomy"));
    predFDNueFHC_lepmomy_fid.SaveTo(fout.mkdir("fd_nue_fhc_lepmomy_fid"));
    predFDNueRHC_lepmomy.SaveTo(fout.mkdir("fd_nue_rhc_lepmomy"));
    predFDNueRHC_lepmomy_fid.SaveTo(fout.mkdir("fd_nue_rhc_lepmomy_fid"));
    predFDNueFHC_lepmomz.SaveTo(fout.mkdir("fd_nue_fhc_lepmomz"));
    predFDNueFHC_lepmomz_fid.SaveTo(fout.mkdir("fd_nue_fhc_lepmomz_fid"));
    predFDNueRHC_lepmomz.SaveTo(fout.mkdir("fd_nue_rhc_lepmomz"));
    predFDNueRHC_lepmomz_fid.SaveTo(fout.mkdir("fd_nue_rhc_lepmomz_fid"));

    predFDNueFHC_lepe.SaveTo(fout.mkdir("fd_nue_fhc_lepe"));
    predFDNueFHC_lepe_fid.SaveTo(fout.mkdir("fd_nue_fhc_lepe_fid"));
    predFDNueRHC_lepe.SaveTo(fout.mkdir("fd_nue_rhc_lepe"));
    predFDNueRHC_lepe_fid.SaveTo(fout.mkdir("fd_nue_rhc_lepe_fid"));
    predFDNueFHC_lepnuangle.SaveTo(fout.mkdir("fd_nue_fhc_lepnuangle"));
    predFDNueFHC_lepnuangle_fid.SaveTo(fout.mkdir("fd_nue_fhc_lepnuangle_fid"));
    predFDNueRHC_lepnuangle.SaveTo(fout.mkdir("fd_nue_rhc_lepnuangle"));
    predFDNueRHC_lepnuangle_fid.SaveTo(fout.mkdir("fd_nue_rhc_lepnuangle_fid"));
    predFDNueFHC_mode.SaveTo(fout.mkdir("fd_nue_fhc_mode"));
    predFDNueFHC_mode_fid.SaveTo(fout.mkdir("fd_nue_fhc_mode_fid"));
    predFDNueRHC_mode.SaveTo(fout.mkdir("fd_nue_rhc_mode"));
    predFDNueRHC_mode_fid.SaveTo(fout.mkdir("fd_nue_rhc_mode_fid"));
    predFDNueFHC_truex.SaveTo(fout.mkdir("fd_nue_fhc_truex"));
    predFDNueFHC_truex_fid.SaveTo(fout.mkdir("fd_nue_fhc_truex_fid"));
    predFDNueRHC_truex.SaveTo(fout.mkdir("fd_nue_rhc_truex"));
    predFDNueRHC_truex_fid.SaveTo(fout.mkdir("fd_nue_rhc_truex_fid"));
    predFDNueFHC_truey.SaveTo(fout.mkdir("fd_nue_fhc_truey"));
    predFDNueFHC_truey_fid.SaveTo(fout.mkdir("fd_nue_fhc_truey_fid"));
    predFDNueRHC_truey.SaveTo(fout.mkdir("fd_nue_rhc_truey"));
    predFDNueRHC_truey_fid.SaveTo(fout.mkdir("fd_nue_rhc_truey_fid"));
    predFDNueFHC_truez.SaveTo(fout.mkdir("fd_nue_fhc_truez"));
    predFDNueFHC_truez_fid.SaveTo(fout.mkdir("fd_nue_fhc_truez_fid"));
    predFDNueRHC_truez.SaveTo(fout.mkdir("fd_nue_rhc_truez"));
    predFDNueRHC_truez_fid.SaveTo(fout.mkdir("fd_nue_rhc_truez_fid"));


    predFDNueFHC_q2.SaveTo(fout.mkdir("fd_nue_fhc_q2"));
    predFDNueFHC_q2_fid.SaveTo(fout.mkdir("fd_nue_fhc_q2_fid"));
    predFDNueRHC_q2.SaveTo(fout.mkdir("fd_nue_rhc_q2"));
    predFDNueRHC_q2_fid.SaveTo(fout.mkdir("fd_nue_rhc_q2_fid"));
    predFDNueFHC_w.SaveTo(fout.mkdir("fd_nue_fhc_w"));
    predFDNueFHC_w_fid.SaveTo(fout.mkdir("fd_nue_fhc_w_fid"));
    predFDNueRHC_w.SaveTo(fout.mkdir("fd_nue_rhc_w"));
    predFDNueRHC_w_fid.SaveTo(fout.mkdir("fd_nue_rhc_w_fid"));
    predFDNueFHC_y.SaveTo(fout.mkdir("fd_nue_fhc_y"));
    predFDNueFHC_y_fid.SaveTo(fout.mkdir("fd_nue_fhc_y_fid"));
    predFDNueRHC_y.SaveTo(fout.mkdir("fd_nue_rhc_y"));
    predFDNueRHC_y_fid.SaveTo(fout.mkdir("fd_nue_rhc_y_fid"));
    predFDNueFHC_np.SaveTo(fout.mkdir("fd_nue_fhc_np"));
    predFDNueFHC_np_fid.SaveTo(fout.mkdir("fd_nue_fhc_np_fid"));
    predFDNueRHC_np.SaveTo(fout.mkdir("fd_nue_rhc_np"));
    predFDNueRHC_np_fid.SaveTo(fout.mkdir("fd_nue_rhc_np_fid"));
    predFDNueFHC_nn.SaveTo(fout.mkdir("fd_nue_fhc_nn"));
    predFDNueFHC_nn_fid.SaveTo(fout.mkdir("fd_nue_fhc_nn_fid"));
    predFDNueRHC_nn.SaveTo(fout.mkdir("fd_nue_rhc_nn"));
    predFDNueRHC_nn_fid.SaveTo(fout.mkdir("fd_nue_rhc_nn_fid"));
    predFDNueFHC_nipip.SaveTo(fout.mkdir("fd_nue_fhc_nipip"));
    predFDNueFHC_nipip_fid.SaveTo(fout.mkdir("fd_nue_fhc_nipip_fid"));
    predFDNueRHC_nipip.SaveTo(fout.mkdir("fd_nue_rhc_nipip"));
    predFDNueRHC_nipip_fid.SaveTo(fout.mkdir("fd_nue_rhc_nipip_fid"));
    predFDNueFHC_nipim.SaveTo(fout.mkdir("fd_nue_fhc_nipim"));
    predFDNueFHC_nipim_fid.SaveTo(fout.mkdir("fd_nue_fhc_nipim_fid"));
    predFDNueRHC_nipim.SaveTo(fout.mkdir("fd_nue_rhc_nipim"));
    predFDNueRHC_nipim_fid.SaveTo(fout.mkdir("fd_nue_rhc_nipim_fid"));
    predFDNueFHC_nipi0.SaveTo(fout.mkdir("fd_nue_fhc_nipi0"));
    predFDNueFHC_nipi0_fid.SaveTo(fout.mkdir("fd_nue_fhc_nipi0_fid"));
    predFDNueRHC_nipi0.SaveTo(fout.mkdir("fd_nue_rhc_nipi0"));
    predFDNueRHC_nipi0_fid.SaveTo(fout.mkdir("fd_nue_rhc_nipi0_fid"));
    predFDNueFHC_recolepen.SaveTo(fout.mkdir("fd_nue_fhc_recolepen"));
    predFDNueFHC_recolepen_fid.SaveTo(fout.mkdir("fd_nue_fhc_recolepen_fid"));
    predFDNueRHC_recolepen.SaveTo(fout.mkdir("fd_nue_rhc_recolepen"));
    predFDNueRHC_recolepen_fid.SaveTo(fout.mkdir("fd_nue_rhc_recolepen_fid"));
    predFDNueFHC_recohaden.SaveTo(fout.mkdir("fd_nue_fhc_recohaden"));
    predFDNueFHC_recohaden_fid.SaveTo(fout.mkdir("fd_nue_fhc_recohaden_fid"));
    predFDNueRHC_recohaden.SaveTo(fout.mkdir("fd_nue_rhc_recohaden"));
    predFDNueRHC_recohaden_fid.SaveTo(fout.mkdir("fd_nue_rhc_recohaden_fid"));

    std::cout << "All done making state..." << std::endl;

    }
  else{    
    std::cout << "Loading state from " << stateFname << std::endl; 
  }

  TFile fin(stateFname);
  PredictionNoExtrap& predFDNumuFHC = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_numu_fhc")).release();
  PredictionNoExtrap& predFDNueFHC = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc")).release();
  PredictionNoExtrap& predFDNumuRHC = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_numu_rhc")).release();
  PredictionNoExtrap& predFDNueRHC = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc")).release();
  PredictionNoExtrap& predFDNumuFHC_Fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_numu_fhc_fid")).release();
  PredictionNoExtrap& predFDNueFHC_Fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_fid")).release();
  PredictionNoExtrap& predFDNumuRHC_Fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_numu_rhc_fid")).release();
  PredictionNoExtrap& predFDNueRHC_Fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_fid")).release();

  PredictionNoExtrap& predFDNueFHC_ev = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_ev")).release();
  PredictionNoExtrap& predFDNueFHC_ev_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_ev_fid")).release();
  PredictionNoExtrap& predFDNueRHC_ev = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_ev")).release();
  PredictionNoExtrap& predFDNueRHC_ev_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_ev_fid")).release();

  PredictionNoExtrap& predFDNueFHC_NuMomX = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_numomx")).release();
  PredictionNoExtrap& predFDNueFHC_NuMomX_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_numomx_fid")).release();
  PredictionNoExtrap& predFDNueRHC_NuMomX = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_numomx")).release();
  PredictionNoExtrap& predFDNueRHC_NuMomX_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_numomx_fid")).release();
  PredictionNoExtrap& predFDNueFHC_NuMomY = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_numomy")).release();
  PredictionNoExtrap& predFDNueFHC_NuMomY_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_numomy_fid")).release();
  PredictionNoExtrap& predFDNueRHC_NuMomY = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_numomy")).release();
  PredictionNoExtrap& predFDNueRHC_NuMomY_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_numomy_fid")).release();
  PredictionNoExtrap& predFDNueFHC_NuMomZ = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_numomz")).release();
  PredictionNoExtrap& predFDNueFHC_NuMomZ_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_numomz_fid")).release();
  PredictionNoExtrap& predFDNueRHC_NuMomZ = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_numomz")).release();
  PredictionNoExtrap& predFDNueRHC_NuMomZ_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_numomz_fid")).release();

  PredictionNoExtrap& predFDNueFHC_LepMomX = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepmomx")).release();
  PredictionNoExtrap& predFDNueFHC_LepMomX_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepmomx_fid")).release();
  PredictionNoExtrap& predFDNueRHC_LepMomX = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepmomx")).release();
  PredictionNoExtrap& predFDNueRHC_LepMomX_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepmomx_fid")).release();
  PredictionNoExtrap& predFDNueFHC_LepMomY = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepmomy")).release();
  PredictionNoExtrap& predFDNueFHC_LepMomY_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepmomy_fid")).release();
  PredictionNoExtrap& predFDNueRHC_LepMomY = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepmomy")).release();
  PredictionNoExtrap& predFDNueRHC_LepMomY_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepmomy_fid")).release();
  PredictionNoExtrap& predFDNueFHC_LepMomZ = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepmomz")).release();
  PredictionNoExtrap& predFDNueFHC_LepMomZ_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepmomz_fid")).release();
  PredictionNoExtrap& predFDNueRHC_LepMomZ = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepmomz")).release();
  PredictionNoExtrap& predFDNueRHC_LepMomZ_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepmomz_fid")).release();
  PredictionNoExtrap& predFDNueFHC_LepE = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepe")).release();
  PredictionNoExtrap& predFDNueFHC_LepE_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepe_fid")).release();
  PredictionNoExtrap& predFDNueRHC_LepE = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepe")).release();
  PredictionNoExtrap& predFDNueRHC_LepE_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepe_fid")).release();
  PredictionNoExtrap& predFDNueFHC_LepNuAngle = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepnuangle")).release();
  PredictionNoExtrap& predFDNueFHC_LepNuAngle_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_lepnuangle_fid")).release();
  PredictionNoExtrap& predFDNueRHC_LepNuAngle = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepnuangle")).release();
  PredictionNoExtrap& predFDNueRHC_LepNuAngle_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_lepnuangle_fid")).release();

  PredictionNoExtrap& predFDNueFHC_mode = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_mode")).release();
  PredictionNoExtrap& predFDNueFHC_mode_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_mode_fid")).release();
  PredictionNoExtrap& predFDNueRHC_mode = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_mode")).release();
  PredictionNoExtrap& predFDNueRHC_mode_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_mode_fid")).release();

  PredictionNoExtrap& predFDNueFHC_truex = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_truex")).release();
  PredictionNoExtrap& predFDNueFHC_truex_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_truex_fid")).release();
  PredictionNoExtrap& predFDNueRHC_truex = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_truex")).release();
  PredictionNoExtrap& predFDNueRHC_truex_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_truex_fid")).release();
  PredictionNoExtrap& predFDNueFHC_truey = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_truey")).release();
  PredictionNoExtrap& predFDNueFHC_truey_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_truey_fid")).release();
  PredictionNoExtrap& predFDNueRHC_truey = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_truey")).release();
  PredictionNoExtrap& predFDNueRHC_truey_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_truey_fid")).release();
  PredictionNoExtrap& predFDNueFHC_truez = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_truez")).release();
  PredictionNoExtrap& predFDNueFHC_truez_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_truez_fid")).release();
  PredictionNoExtrap& predFDNueRHC_truez = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_truez")).release();
  PredictionNoExtrap& predFDNueRHC_truez_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_truez_fid")).release();


  PredictionNoExtrap& predFDNueFHC_q2 = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_q2")).release();
  PredictionNoExtrap& predFDNueFHC_q2_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_q2_fid")).release();
  PredictionNoExtrap& predFDNueRHC_q2 = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_q2")).release();
  PredictionNoExtrap& predFDNueRHC_q2_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_q2_fid")).release();
  PredictionNoExtrap& predFDNueFHC_w = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_w")).release();
  PredictionNoExtrap& predFDNueFHC_w_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_w_fid")).release();
  PredictionNoExtrap& predFDNueRHC_w = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_w")).release();
  PredictionNoExtrap& predFDNueRHC_w_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_w_fid")).release();
  PredictionNoExtrap& predFDNueFHC_y = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_y")).release();
  PredictionNoExtrap& predFDNueFHC_y_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_y_fid")).release();
  PredictionNoExtrap& predFDNueRHC_y = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_y")).release();
  PredictionNoExtrap& predFDNueRHC_y_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_y_fid")).release();

  PredictionNoExtrap& predFDNueFHC_nP = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_np")).release();
  PredictionNoExtrap& predFDNueFHC_nP_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_np_fid")).release();
  PredictionNoExtrap& predFDNueRHC_nP = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_np")).release();
  PredictionNoExtrap& predFDNueRHC_nP_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_np_fid")).release();
  PredictionNoExtrap& predFDNueFHC_nN = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_nn")).release();
  PredictionNoExtrap& predFDNueFHC_nN_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_nn_fid")).release();
  PredictionNoExtrap& predFDNueRHC_nN = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_nn")).release();
  PredictionNoExtrap& predFDNueRHC_nN_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_nn_fid")).release();

  PredictionNoExtrap& predFDNueFHC_nipip = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_nipip")).release();
  PredictionNoExtrap& predFDNueFHC_nipip_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_nipip_fid")).release();
  PredictionNoExtrap& predFDNueRHC_nipip = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_nipip")).release();
  PredictionNoExtrap& predFDNueRHC_nipip_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_nipip_fid")).release();
  PredictionNoExtrap& predFDNueFHC_nipim = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_nipim")).release();
  PredictionNoExtrap& predFDNueFHC_nipim_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_nipim_fid")).release();
  PredictionNoExtrap& predFDNueRHC_nipim = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_nipim")).release();
  PredictionNoExtrap& predFDNueRHC_nipim_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_nipim_fid")).release();
  PredictionNoExtrap& predFDNueFHC_nipi0 = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_nipi0")).release();
  PredictionNoExtrap& predFDNueFHC_nipi0_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_nipi0_fid")).release();
  PredictionNoExtrap& predFDNueRHC_nipi0 = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_nipi0")).release();
  PredictionNoExtrap& predFDNueRHC_nipi0_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_nipi0_fid")).release();

  PredictionNoExtrap& predFDNueFHC_RecoLepEn = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_recolepen")).release();
  PredictionNoExtrap& predFDNueFHC_RecoLepEn_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_recolepen_fid")).release();
  PredictionNoExtrap& predFDNueRHC_RecoLepEn = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_recolepen")).release();
  PredictionNoExtrap& predFDNueRHC_RecoLepEn_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_recolepen_fid")).release();

  PredictionNoExtrap& predFDNueFHC_RecoHadEn = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_recohaden")).release();
  PredictionNoExtrap& predFDNueFHC_RecoHadEn_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc_recohaden_fid")).release();
  PredictionNoExtrap& predFDNueRHC_RecoHadEn = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_recohaden")).release();
  PredictionNoExtrap& predFDNueRHC_RecoHadEn_fid = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc_recohaden_fid")).release();


  fin.Close();
  std::cout << "Done loading state" << std::endl;

  TFile* fout = new TFile(outputFname, "RECREATE");

  osc::NoOscillations noOsc;
  TH1* hnumufhc_sig_unosc = predFDNumuFHC.PredictComponent(&noOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);
  TH1* hnumufhc_sig_unosc_fid = predFDNumuFHC_Fid.PredictComponent(&noOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);

  std::string dcpnames[4] = {"0pi","piover2","pi","3piover2"};

  for(int hie = -1; hie <= +1; hie += 2){

    osc::IOscCalculatorAdjustable* inputOsc = NuFitOscCalc(hie);

    const std::string hieStr = (hie > 0) ? "nh" : "ih";

    for(int deltaIdx = 0; deltaIdx < 4; ++deltaIdx){
      inputOsc->SetdCP(deltaIdx/2.*TMath::Pi());
      const std::string dcpStr = dcpnames[deltaIdx];

      TH1* hnumufhc = predFDNumuFHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnuefhc = predFDNueFHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnumurhc = predFDNumuRHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnuerhc = predFDNueRHC.Predict(inputOsc).ToTH1(potFD);

      //FHC Dis
      TH1* hnumufhc_sig = predFDNumuFHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);
      TH1* hnumufhc_ncbg = predFDNumuFHC.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumufhc_nutaubg = predFDNumuFHC.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumufhc_wsbg = predFDNumuFHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
      TH1* hnumufhc_nuebg = predFDNumuFHC.PredictComponent(inputOsc, Flavors::kAllNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);

      
      TH1* hnumufhc_sig_fid = predFDNumuFHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);
      TH1* hnumufhc_ncbg_fid = predFDNumuFHC_Fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumufhc_nutaubg_fid = predFDNumuFHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumufhc_wsbg_fid = predFDNumuFHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
      TH1* hnumufhc_nuebg_fid = predFDNumuFHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);

      //RHC Dis
      TH1* hnumurhc_sig = predFDNumuRHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
      TH1* hnumurhc_ncbg = predFDNumuRHC.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumurhc_nutaubg = predFDNumuRHC.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumurhc_wsbg = predFDNumuRHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);

      TH1* hnumurhc_sig_fid = predFDNumuRHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
      TH1* hnumurhc_ncbg_fid = predFDNumuRHC_Fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumurhc_nutaubg_fid = predFDNumuRHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumurhc_wsbg_fid = predFDNumuRHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);

      //FHC App
      TH1* hnuefhc_sig = predFDNueFHC.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_ncbg = predFDNueFHC.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_beambg = predFDNueFHC.PredictComponent(inputOsc, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_nutaubg = predFDNueFHC.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_numubg = predFDNueFHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(potFD);

      TH1* hnuefhc_sig_fid = predFDNueFHC_Fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_ncbg_fid = predFDNueFHC_Fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_beambg_fid = predFDNueFHC_Fid.PredictComponent(inputOsc, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_nutaubg_fid = predFDNueFHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_numubg_fid = predFDNueFHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(potFD);

      //RHC App
      TH1* hnuerhc_sig = predFDNueRHC.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_ncbg = predFDNueRHC.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_beambg = predFDNueRHC.PredictComponent(inputOsc, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_nutaubg = predFDNueRHC.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_numubg = predFDNueRHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(potFD);

      TH1* hnuerhc_sig_fid = predFDNueRHC_Fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_ncbg_fid = predFDNueRHC_Fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_beambg_fid = predFDNueRHC_Fid.PredictComponent(inputOsc, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_nutaubg_fid = predFDNueRHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_numubg_fid = predFDNueRHC_Fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(potFD);

      //Efficiency studies
      TH2* hnuefhc_ev_sig = predFDNueFHC_ev.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_ev_sig_fid = predFDNueFHC_ev_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_ev_sig = predFDNueRHC_ev.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_ev_sig_fid = predFDNueRHC_ev_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_numomx_sig = predFDNueFHC_NuMomX.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomx_sig_fid = predFDNueFHC_NuMomX_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomx_sig = predFDNueRHC_NuMomX.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomx_sig_fid = predFDNueRHC_NuMomX_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomy_sig = predFDNueFHC_NuMomY.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomy_sig_fid = predFDNueFHC_NuMomY_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomy_sig = predFDNueRHC_NuMomY.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomy_sig_fid = predFDNueRHC_NuMomY_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomz_sig = predFDNueFHC_NuMomZ.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomz_sig_fid = predFDNueFHC_NuMomZ_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomz_sig = predFDNueRHC_NuMomZ.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomz_sig_fid = predFDNueRHC_NuMomZ_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_lepmomx_sig = predFDNueFHC_LepMomX.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomx_sig_fid = predFDNueFHC_LepMomX_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomx_sig = predFDNueRHC_LepMomX.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomx_sig_fid = predFDNueRHC_LepMomX_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomy_sig = predFDNueFHC_LepMomY.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomy_sig_fid = predFDNueFHC_LepMomY_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomy_sig = predFDNueRHC_LepMomY.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomy_sig_fid = predFDNueRHC_LepMomY_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomz_sig = predFDNueFHC_LepMomZ.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomz_sig_fid = predFDNueFHC_LepMomZ_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomz_sig = predFDNueRHC_LepMomZ.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomz_sig_fid = predFDNueRHC_LepMomZ_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepe_sig = predFDNueFHC_LepE.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepe_sig_fid = predFDNueFHC_LepE_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepe_sig = predFDNueRHC_LepE.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepe_sig_fid = predFDNueRHC_LepE_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepnuangle_sig = predFDNueFHC_LepNuAngle.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepnuangle_sig_fid = predFDNueFHC_LepNuAngle_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepnuangle_sig = predFDNueRHC_LepNuAngle.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepnuangle_sig_fid = predFDNueRHC_LepNuAngle_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);


      TH2* hnuefhc_mode_sig = predFDNueFHC_mode.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_mode_sig_fid = predFDNueFHC_mode_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_mode_sig = predFDNueRHC_mode.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_mode_sig_fid = predFDNueRHC_mode_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_truex_sig = predFDNueFHC_truex.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truex_sig_fid = predFDNueFHC_truex_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truex_sig = predFDNueRHC_truex.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truex_sig_fid = predFDNueRHC_truex_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truey_sig = predFDNueFHC_truey.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truey_sig_fid = predFDNueFHC_truey_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truey_sig = predFDNueRHC_truey.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truey_sig_fid = predFDNueRHC_truey_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truez_sig = predFDNueFHC_truez.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truez_sig_fid = predFDNueFHC_truez_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truez_sig = predFDNueRHC_truez.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truez_sig_fid = predFDNueRHC_truez_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);


      TH2* hnuefhc_q2_sig = predFDNueFHC_q2.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_q2_sig_fid = predFDNueFHC_q2_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_q2_sig = predFDNueRHC_q2.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_q2_sig_fid = predFDNueRHC_q2_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_w_sig = predFDNueFHC_w.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_w_sig_fid = predFDNueFHC_w_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_w_sig = predFDNueRHC_w.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_w_sig_fid = predFDNueRHC_w_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_y_sig = predFDNueFHC_y.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_y_sig_fid = predFDNueFHC_y_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_y_sig = predFDNueRHC_y.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_y_sig_fid = predFDNueRHC_y_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_np_sig = predFDNueFHC_nP.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_np_sig_fid = predFDNueFHC_nP_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_np_sig = predFDNueRHC_nP.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_np_sig_fid = predFDNueRHC_nP_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nn_sig = predFDNueFHC_nN.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nn_sig_fid = predFDNueFHC_nN_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nn_sig = predFDNueRHC_nN.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nn_sig_fid = predFDNueRHC_nN_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_nipip_sig = predFDNueFHC_nipip.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipip_sig_fid = predFDNueFHC_nipip_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipip_sig = predFDNueRHC_nipip.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipip_sig_fid = predFDNueRHC_nipip_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipim_sig = predFDNueFHC_nipim.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipim_sig_fid = predFDNueFHC_nipim_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipim_sig = predFDNueRHC_nipim.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipim_sig_fid = predFDNueRHC_nipim_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipi0_sig = predFDNueFHC_nipi0.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipi0_sig_fid = predFDNueFHC_nipi0_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipi0_sig = predFDNueRHC_nipi0.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipi0_sig_fid = predFDNueRHC_nipi0_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_recolepen_sig = predFDNueFHC_RecoLepEn.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_recolepen_sig_fid = predFDNueFHC_RecoLepEn_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_recolepen_sig = predFDNueRHC_RecoLepEn.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_recolepen_sig_fid = predFDNueRHC_RecoLepEn_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_recohaden_sig = predFDNueFHC_RecoHadEn.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_recohaden_sig_fid = predFDNueFHC_RecoHadEn_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_recohaden_sig = predFDNueRHC_RecoHadEn.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_recohaden_sig_fid = predFDNueRHC_RecoHadEn_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH2(potFD);


      TH2* hnuefhc_ev_ncbg = predFDNueFHC_ev.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_ev_ncbg_fid = predFDNueFHC_ev_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_ev_ncbg = predFDNueRHC_ev.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_ev_ncbg_fid = predFDNueRHC_ev_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomx_ncbg = predFDNueFHC_NuMomX.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomx_ncbg_fid = predFDNueFHC_NuMomX_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomx_ncbg = predFDNueRHC_NuMomX.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomx_ncbg_fid = predFDNueRHC_NuMomX_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomy_ncbg = predFDNueFHC_NuMomY.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomy_ncbg_fid = predFDNueFHC_NuMomY_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomy_ncbg = predFDNueRHC_NuMomY.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomy_ncbg_fid = predFDNueRHC_NuMomY_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomz_ncbg = predFDNueFHC_NuMomZ.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_numomz_ncbg_fid = predFDNueFHC_NuMomZ_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomz_ncbg = predFDNueRHC_NuMomZ.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_numomz_ncbg_fid = predFDNueRHC_NuMomZ_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_lepmomx_ncbg = predFDNueFHC_LepMomX.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomx_ncbg_fid = predFDNueFHC_LepMomX_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomx_ncbg = predFDNueRHC_LepMomX.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomx_ncbg_fid = predFDNueRHC_LepMomX_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomy_ncbg = predFDNueFHC_LepMomY.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomy_ncbg_fid = predFDNueFHC_LepMomY_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomy_ncbg = predFDNueRHC_LepMomY.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomy_ncbg_fid = predFDNueRHC_LepMomY_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomz_ncbg = predFDNueFHC_LepMomZ.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepmomz_ncbg_fid = predFDNueFHC_LepMomZ_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomz_ncbg = predFDNueRHC_LepMomZ.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepmomz_ncbg_fid = predFDNueRHC_LepMomZ_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepe_ncbg = predFDNueFHC_LepE.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepe_ncbg_fid = predFDNueFHC_LepE_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepe_ncbg = predFDNueRHC_LepE.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepe_ncbg_fid = predFDNueRHC_LepE_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepnuangle_ncbg = predFDNueFHC_LepNuAngle.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_lepnuangle_ncbg_fid = predFDNueFHC_LepNuAngle_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepnuangle_ncbg = predFDNueRHC_LepNuAngle.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_lepnuangle_ncbg_fid = predFDNueRHC_LepNuAngle_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_mode_ncbg = predFDNueFHC_mode.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_mode_ncbg_fid = predFDNueFHC_mode_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_mode_ncbg = predFDNueRHC_mode.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_mode_ncbg_fid = predFDNueRHC_mode_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_truex_ncbg = predFDNueFHC_truex.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truex_ncbg_fid = predFDNueFHC_truex_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truex_ncbg = predFDNueRHC_truex.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truex_ncbg_fid = predFDNueRHC_truex_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truey_ncbg = predFDNueFHC_truey.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truey_ncbg_fid = predFDNueFHC_truey_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truey_ncbg = predFDNueRHC_truey.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truey_ncbg_fid = predFDNueRHC_truey_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truez_ncbg = predFDNueFHC_truez.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_truez_ncbg_fid = predFDNueFHC_truez_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truez_ncbg = predFDNueRHC_truez.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_truez_ncbg_fid = predFDNueRHC_truez_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_q2_ncbg = predFDNueFHC_q2.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_q2_ncbg_fid = predFDNueFHC_q2_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_q2_ncbg = predFDNueRHC_q2.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_q2_ncbg_fid = predFDNueRHC_q2_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_w_ncbg = predFDNueFHC_w.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_w_ncbg_fid = predFDNueFHC_w_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_w_ncbg = predFDNueRHC_w.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_w_ncbg_fid = predFDNueRHC_w_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_y_ncbg = predFDNueFHC_y.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_y_ncbg_fid = predFDNueFHC_y_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_y_ncbg = predFDNueRHC_y.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_y_ncbg_fid = predFDNueRHC_y_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_np_ncbg = predFDNueFHC_nP.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_np_ncbg_fid = predFDNueFHC_nP_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_np_ncbg = predFDNueRHC_nP.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_np_ncbg_fid = predFDNueRHC_nP_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nn_ncbg = predFDNueFHC_nN.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nn_ncbg_fid = predFDNueFHC_nN_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nn_ncbg = predFDNueRHC_nN.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nn_ncbg_fid = predFDNueRHC_nN_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_nipip_ncbg = predFDNueFHC_nipip.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipip_ncbg_fid = predFDNueFHC_nipip_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipip_ncbg = predFDNueRHC_nipip.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipip_ncbg_fid = predFDNueRHC_nipip_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipim_ncbg = predFDNueFHC_nipim.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipim_ncbg_fid = predFDNueFHC_nipim_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipim_ncbg = predFDNueRHC_nipim.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipim_ncbg_fid = predFDNueRHC_nipim_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipi0_ncbg = predFDNueFHC_nipi0.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_nipi0_ncbg_fid = predFDNueFHC_nipi0_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipi0_ncbg = predFDNueRHC_nipi0.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_nipi0_ncbg_fid = predFDNueRHC_nipi0_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_recolepen_ncbg = predFDNueFHC_RecoLepEn.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_recolepen_ncbg_fid = predFDNueFHC_RecoLepEn_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_recolepen_ncbg = predFDNueRHC_RecoLepEn.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_recolepen_ncbg_fid = predFDNueRHC_RecoLepEn_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);

      TH2* hnuefhc_recohaden_ncbg = predFDNueFHC_RecoHadEn.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuefhc_recohaden_ncbg_fid = predFDNueFHC_RecoHadEn_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_recohaden_ncbg = predFDNueRHC_RecoHadEn.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);
      TH2* hnuerhc_recohaden_ncbg_fid = predFDNueRHC_RecoHadEn_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH2(potFD);


      hnumufhc_sig_unosc->Write("numu_fhc_unosc_sig");
      hnumufhc_sig_unosc_fid->Write("fid_numu_fhc_unosc_sig");

      hnumufhc->Write(("numu_fhc_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_sig->Write(("numu_fhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_ncbg->Write(("numu_fhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_nutaubg->Write(("numu_fhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_nuebg->Write(("numu_fhc_nuebg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_wsbg->Write(("numu_fhc_wsbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc->Write(("nue_fhc_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_sig->Write(("nue_fhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_ncbg->Write(("nue_fhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_beambg->Write(("nue_fhc_beambg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nutaubg->Write(("nue_fhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numubg->Write(("nue_fhc_numubg_"+hieStr+"_"+dcpStr).c_str());
	
      hnumurhc->Write(("numu_rhc_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_sig->Write(("numu_rhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_ncbg->Write(("numu_rhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_nutaubg->Write(("numu_rhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_wsbg->Write(("numu_rhc_wsbg_"+hieStr+"_"+dcpStr).c_str());

      hnuerhc->Write(("nue_rhc_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_sig->Write(("nue_rhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_ncbg->Write(("nue_rhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_beambg->Write(("nue_rhc_beambg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nutaubg->Write(("nue_rhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numubg->Write(("nue_rhc_numubg_"+hieStr+"_"+dcpStr).c_str());

      hnumufhc_sig_fid->Write(("fid_numu_fhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_ncbg_fid->Write(("fid_numu_fhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_nutaubg_fid->Write(("fid_numu_fhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_nuebg_fid->Write(("fid_numu_fhc_nuebg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_wsbg_fid->Write(("fid_numu_fhc_wsbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_sig_fid->Write(("fid_nue_fhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_ncbg_fid->Write(("fid_nue_fhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_beambg_fid->Write(("fid_nue_fhc_beambg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nutaubg_fid->Write(("fid_nue_fhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numubg_fid->Write(("fid_nue_fhc_numubg_"+hieStr+"_"+dcpStr).c_str());
	
      hnumurhc_sig_fid->Write(("fid_numu_rhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_ncbg_fid->Write(("fid_numu_rhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_nutaubg_fid->Write(("fid_numu_rhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_wsbg_fid->Write(("fid_numu_rhc_wsbg_"+hieStr+"_"+dcpStr).c_str());

      hnuerhc_sig_fid->Write(("fid_nue_rhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_ncbg_fid->Write(("fid_nue_rhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_beambg_fid->Write(("fid_nue_rhc_beambg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nutaubg_fid->Write(("fid_nue_rhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numubg_fid->Write(("fid_nue_rhc_numubg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_ev_sig->Write(("nue_fhc_ev_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_ev_sig_fid->Write(("fid_nue_fhc_ev_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_ev_sig->Write(("nue_rhc_ev_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_ev_sig_fid->Write(("fid_nue_rhc_ev_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomx_sig->Write(("nue_fhc_numomx_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomx_sig_fid->Write(("fid_nue_fhc_numomx_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomx_sig->Write(("nue_rhc_numomx_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomx_sig_fid->Write(("fid_nue_rhc_numomx_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomy_sig->Write(("nue_fhc_numomy_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomy_sig_fid->Write(("fid_nue_fhc_numomy_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomy_sig->Write(("nue_rhc_numomy_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomy_sig_fid->Write(("fid_nue_rhc_numomy_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomz_sig->Write(("nue_fhc_numomz_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomz_sig_fid->Write(("fid_nue_fhc_numomz_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomz_sig->Write(("nue_rhc_numomz_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomz_sig_fid->Write(("fid_nue_rhc_numomz_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomx_sig->Write(("nue_fhc_lepmomx_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomx_sig_fid->Write(("fid_nue_fhc_lepmomx_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomx_sig->Write(("nue_rhc_lepmomx_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomx_sig_fid->Write(("fid_nue_rhc_lepmomx_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomy_sig->Write(("nue_fhc_lepmomy_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomy_sig_fid->Write(("fid_nue_fhc_lepmomy_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomy_sig->Write(("nue_rhc_lepmomy_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomy_sig_fid->Write(("fid_nue_rhc_lepmomy_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomz_sig->Write(("nue_fhc_lepmomz_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomz_sig_fid->Write(("fid_nue_fhc_lepmomz_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomz_sig->Write(("nue_rhc_lepmomz_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomz_sig_fid->Write(("fid_nue_rhc_lepmomz_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepe_sig->Write(("nue_fhc_lepe_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepe_sig_fid->Write(("fid_nue_fhc_lepe_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepe_sig->Write(("nue_rhc_lepe_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepe_sig_fid->Write(("fid_nue_rhc_lepe_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepnuangle_sig->Write(("nue_fhc_lepnuangle_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepnuangle_sig_fid->Write(("fid_nue_fhc_lepnuangle_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepnuangle_sig->Write(("nue_rhc_lepnuangle_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepnuangle_sig_fid->Write(("fid_nue_rhc_lepnuangle_sig_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_mode_sig->Write(("nue_fhc_mode_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_mode_sig_fid->Write(("fid_nue_fhc_mode_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_mode_sig->Write(("nue_rhc_mode_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_mode_sig_fid->Write(("fid_nue_rhc_mode_sig_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_truex_sig->Write(("nue_fhc_truex_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truex_sig_fid->Write(("fid_nue_fhc_truex_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truex_sig->Write(("nue_rhc_truex_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truex_sig_fid->Write(("fid_nue_rhc_truex_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truey_sig->Write(("nue_fhc_truey_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truey_sig_fid->Write(("fid_nue_fhc_truey_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truey_sig->Write(("nue_rhc_truey_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truey_sig_fid->Write(("fid_nue_rhc_truey_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truez_sig->Write(("nue_fhc_truez_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truez_sig_fid->Write(("fid_nue_fhc_truez_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truez_sig->Write(("nue_rhc_truez_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truez_sig_fid->Write(("fid_nue_rhc_truez_sig_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_q2_sig->Write(("nue_fhc_q2_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_q2_sig_fid->Write(("fid_nue_fhc_q2_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_q2_sig->Write(("nue_rhc_q2_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_q2_sig_fid->Write(("fid_nue_rhc_q2_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_w_sig->Write(("nue_fhc_w_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_w_sig_fid->Write(("fid_nue_fhc_w_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_w_sig->Write(("nue_rhc_w_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_w_sig_fid->Write(("fid_nue_rhc_w_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_y_sig->Write(("nue_fhc_y_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_y_sig_fid->Write(("fid_nue_fhc_y_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_y_sig->Write(("nue_rhc_y_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_y_sig_fid->Write(("fid_nue_rhc_y_sig_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_np_sig->Write(("nue_fhc_np_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_np_sig_fid->Write(("fid_nue_fhc_np_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_np_sig->Write(("nue_rhc_np_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_np_sig_fid->Write(("fid_nue_rhc_np_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nn_sig->Write(("nue_fhc_nn_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nn_sig_fid->Write(("fid_nue_fhc_nn_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nn_sig->Write(("nue_rhc_nn_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nn_sig_fid->Write(("fid_nue_rhc_nn_sig_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_nipip_sig->Write(("nue_fhc_nipip_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipip_sig_fid->Write(("fid_nue_fhc_nipip_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipip_sig->Write(("nue_rhc_nipip_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipip_sig_fid->Write(("fid_nue_rhc_nipip_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipim_sig->Write(("nue_fhc_nipim_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipim_sig_fid->Write(("fid_nue_fhc_nipim_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipim_sig->Write(("nue_rhc_nipim_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipim_sig_fid->Write(("fid_nue_rhc_nipim_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipi0_sig->Write(("nue_fhc_nipi0_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipi0_sig_fid->Write(("fid_nue_fhc_nipi0_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipi0_sig->Write(("nue_rhc_nipi0_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipi0_sig_fid->Write(("fid_nue_rhc_nipi0_sig_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_recolepen_sig->Write(("nue_fhc_recolepen_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_recolepen_sig_fid->Write(("fid_nue_fhc_recolepen_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_recolepen_sig->Write(("nue_rhc_recolepen_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_recolepen_sig_fid->Write(("fid_nue_rhc_recolepen_sig_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_recohaden_sig->Write(("nue_fhc_recohaden_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_recohaden_sig_fid->Write(("fid_nue_fhc_recohaden_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_recohaden_sig->Write(("nue_rhc_recohaden_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_recohaden_sig_fid->Write(("fid_nue_rhc_recohaden_sig_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_ev_ncbg->Write(("nue_fhc_ev_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_ev_ncbg_fid->Write(("fid_nue_fhc_ev_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_ev_ncbg->Write(("nue_rhc_ev_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_ev_ncbg_fid->Write(("fid_nue_rhc_ev_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomx_ncbg->Write(("nue_fhc_numomx_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomx_ncbg_fid->Write(("fid_nue_fhc_numomx_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomx_ncbg->Write(("nue_rhc_numomx_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomx_ncbg_fid->Write(("fid_nue_rhc_numomx_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomy_ncbg->Write(("nue_fhc_numomy_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomy_ncbg_fid->Write(("fid_nue_fhc_numomy_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomy_ncbg->Write(("nue_rhc_numomy_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomy_ncbg_fid->Write(("fid_nue_rhc_numomy_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomz_ncbg->Write(("nue_fhc_numomz_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numomz_ncbg_fid->Write(("fid_nue_fhc_numomz_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomz_ncbg->Write(("nue_rhc_numomz_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numomz_ncbg_fid->Write(("fid_nue_rhc_numomz_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomx_ncbg->Write(("nue_fhc_lepmomx_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomx_ncbg_fid->Write(("fid_nue_fhc_lepmomx_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomx_ncbg->Write(("nue_rhc_lepmomx_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomx_ncbg_fid->Write(("fid_nue_rhc_lepmomx_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomy_ncbg->Write(("nue_fhc_lepmomy_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomy_ncbg_fid->Write(("fid_nue_fhc_lepmomy_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomy_ncbg->Write(("nue_rhc_lepmomy_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomy_ncbg_fid->Write(("fid_nue_rhc_lepmomy_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomz_ncbg->Write(("nue_fhc_lepmomz_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepmomz_ncbg_fid->Write(("fid_nue_fhc_lepmomz_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomz_ncbg->Write(("nue_rhc_lepmomz_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepmomz_ncbg_fid->Write(("fid_nue_rhc_lepmomz_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepe_ncbg->Write(("nue_fhc_lepe_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepe_ncbg_fid->Write(("fid_nue_fhc_lepe_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepe_ncbg->Write(("nue_rhc_lepe_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepe_ncbg_fid->Write(("fid_nue_rhc_lepe_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepnuangle_ncbg->Write(("nue_fhc_lepnuangle_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_lepnuangle_ncbg_fid->Write(("fid_nue_fhc_lepnuangle_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepnuangle_ncbg->Write(("nue_rhc_lepnuangle_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_lepnuangle_ncbg_fid->Write(("fid_nue_rhc_lepnuangle_ncbg_"+hieStr+"_"+dcpStr).c_str());


      hnuefhc_mode_ncbg->Write(("nue_fhc_mode_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_mode_ncbg_fid->Write(("fid_nue_fhc_mode_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_mode_ncbg->Write(("nue_rhc_mode_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_mode_ncbg_fid->Write(("fid_nue_rhc_mode_ncbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_truex_ncbg->Write(("nue_fhc_truex_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truex_ncbg_fid->Write(("fid_nue_fhc_truex_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truex_ncbg->Write(("nue_rhc_truex_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truex_ncbg_fid->Write(("fid_nue_rhc_truex_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truey_ncbg->Write(("nue_fhc_truey_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truey_ncbg_fid->Write(("fid_nue_fhc_truey_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truey_ncbg->Write(("nue_rhc_truey_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truey_ncbg_fid->Write(("fid_nue_rhc_truey_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truez_ncbg->Write(("nue_fhc_truez_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_truez_ncbg_fid->Write(("fid_nue_fhc_truez_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truez_ncbg->Write(("nue_rhc_truez_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_truez_ncbg_fid->Write(("fid_nue_rhc_truez_ncbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_q2_ncbg->Write(("nue_fhc_q2_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_q2_ncbg_fid->Write(("fid_nue_fhc_q2_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_q2_ncbg->Write(("nue_rhc_q2_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_q2_ncbg_fid->Write(("fid_nue_rhc_q2_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_w_ncbg->Write(("nue_fhc_w_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_w_ncbg_fid->Write(("fid_nue_fhc_w_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_w_ncbg->Write(("nue_rhc_w_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_w_ncbg_fid->Write(("fid_nue_rhc_w_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_y_ncbg->Write(("nue_fhc_y_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_y_ncbg_fid->Write(("fid_nue_fhc_y_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_y_ncbg->Write(("nue_rhc_y_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_y_ncbg_fid->Write(("fid_nue_rhc_y_ncbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_np_ncbg->Write(("nue_fhc_np_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_np_ncbg_fid->Write(("fid_nue_fhc_np_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_np_ncbg->Write(("nue_rhc_np_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_np_ncbg_fid->Write(("fid_nue_rhc_np_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nn_ncbg->Write(("nue_fhc_nn_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nn_ncbg_fid->Write(("fid_nue_fhc_nn_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nn_ncbg->Write(("nue_rhc_nn_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nn_ncbg_fid->Write(("fid_nue_rhc_nn_ncbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_nipip_ncbg->Write(("nue_fhc_nipip_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipip_ncbg_fid->Write(("fid_nue_fhc_nipip_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipip_ncbg->Write(("nue_rhc_nipip_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipip_ncbg_fid->Write(("fid_nue_rhc_nipip_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipim_ncbg->Write(("nue_fhc_nipim_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipim_ncbg_fid->Write(("fid_nue_fhc_nipim_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipim_ncbg->Write(("nue_rhc_nipim_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipim_ncbg_fid->Write(("fid_nue_rhc_nipim_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipi0_ncbg->Write(("nue_fhc_nipi0_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nipi0_ncbg_fid->Write(("fid_nue_fhc_nipi0_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipi0_ncbg->Write(("nue_rhc_nipi0_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nipi0_ncbg_fid->Write(("fid_nue_rhc_nipi0_ncbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_recolepen_ncbg->Write(("nue_fhc_recolepen_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_recolepen_ncbg_fid->Write(("fid_nue_fhc_recolepen_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_recolepen_ncbg->Write(("nue_rhc_recolepen_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_recolepen_ncbg_fid->Write(("fid_nue_rhc_recolepen_ncbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_recohaden_ncbg->Write(("nue_fhc_recohaden_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_recohaden_ncbg_fid->Write(("fid_nue_fhc_recohaden_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_recohaden_ncbg->Write(("nue_rhc_recohaden_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_recohaden_ncbg_fid->Write(("fid_nue_rhc_recohaden_ncbg_"+hieStr+"_"+dcpStr).c_str());

    }
  }
  fout->Close();
}

