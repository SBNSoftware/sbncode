/**
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Random/RandFlat.h"
#include "IFDH_service.h"

#include "nusimdata/SimulationBase/MCFlux.h"

// local includes
#include "IMesonGen.h"

// LArSoft includes
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/NuChoice.h"

// ROOT
#include "TVector3.h"
#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  NuMiKaonGen class definiton
 */
class NuMiKaonGen : public IMesonGen
{
public:
    /**
     *  @brief  Constructor
     */
    NuMiKaonGen(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~NuMiKaonGen();

    double GetPOT() override;
    simb::MCFlux GetNext() override;

    void configure(const fhicl::ParameterSet&) override;

    const bsim::Dk2Nu *GetNextEntry();
    std::vector<std::string> LoadFluxFiles();
    simb::MCFlux MakeMCFlux(const bsim::Dk2Nu &dk2nu);
    double LoadPOT();

    // no weights
    float MaxWeight() override { return -1.; }

private:
  // config
  std::string fSearchPath;
  std::vector<std::string> fSearchPatterns;
  unsigned long fMaxFluxFileMB;
  std::string fFluxCopyMethod;
  bool fRandomizeFiles;

  std::string fTreeName;
  std::string fMetaTreeName;

  // info for tracking files
  unsigned fFileIndex;
  bool fNewFile;
  std::vector<std::string> fFluxFiles;

  // info for tracking entry in file
  unsigned fEntry;
  unsigned fEntryStart;

  // ROOT Holders
  TTree *fFluxTree;
  TFile *fFluxFile;
  bsim::Dk2Nu *fDk2Nu;

  // count POT
  double fAccumulatedPOT;
  double fThisFilePOT;
  
};

NuMiKaonGen::NuMiKaonGen(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("NuMiKaonGen") 
{
  configure(pset);

  // copy the flux files locally
  fFluxFiles = LoadFluxFiles();

  // setup indices
  fFileIndex = 0;
  fEntry = 0;
  fEntryStart = 0;
  fNewFile = true;
  fFluxTree = NULL;
  fFluxFile = NULL;
  fDk2Nu = new bsim::Dk2Nu;

  fAccumulatedPOT = 0.;
  fThisFilePOT = 0.;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

NuMiKaonGen::~NuMiKaonGen()
{

  if (fDk2Nu) delete fDk2Nu;
}

//------------------------------------------------------------------------------------------------------------------------------------------
void NuMiKaonGen::configure(fhicl::ParameterSet const &pset)
{
  fSearchPath = pset.get<std::string>("SearchPath");
  fSearchPatterns = pset.get<std::vector<std::string>>("FluxFiles");
  fMaxFluxFileMB = pset.get<unsigned long>("MaxFluxFileMB", 2 * 1024);
  fFluxCopyMethod = pset.get<std::string>("FluxCopyMethod", "IFDH");
  fTreeName = pset.get<std::string>("TreeName");
  fMetaTreeName = pset.get<std::string>("MetaTreeName");
  fRandomizeFiles = pset.get<bool>("RandomizeFiles");

  std::cout << "Searching for flux files at path: " << fSearchPath << std::endl;
  std::cout << "With patterns:\n";
  for (const std::string &s: fSearchPatterns) std::cout << s << std::endl;
  std::cout << "With copy method: " << fFluxCopyMethod << std::endl;
}

std::vector<std::string> NuMiKaonGen::LoadFluxFiles() {
  art::ServiceHandle<IFDH> ifdhp;

  std::vector<std::pair<std::string, long>> allFiles;

  // find the flux files
  for (unsigned i = 0; i < fSearchPatterns.size(); i++) {
    std::vector<std::pair<std::string, long>> thisList = ifdhp->findMatchingFiles(fSearchPath, fSearchPatterns[i]);
    std::copy (thisList.begin(), thisList.end(), std::back_inserter(allFiles));
  }

  // first randomize the flux files
  std::vector<unsigned long> order(allFiles.size(), 0);
  if (fRandomizeFiles) {
    std::vector<double> rand(allFiles.size(), 0.);
    CLHEP::RandFlat::shootArray(fEngine, rand.size(), &rand[0]);
    TMath::Sort(allFiles.size(), &rand[0], &order[0], false);
  }
  else {
    for (unsigned i = 0; i < order.size(); i++) {
      order[i] = i;
    }
  }

  // If we are directly accessing the files, no need to copy
  if (fFluxCopyMethod == "DIRECT") {
    std::cout << "DIRECTLY ACCESSING FLUX FILES.\n";
    std::vector<std::string> files(allFiles.size());
    for (unsigned i = 0; i < order.size(); i++) {
      files[i] = allFiles[order[i]].first;
    }
    return files;
  }

  // copy over up to the provided limit
  std::vector<std::pair<std::string, long>> selected;
  unsigned long totalBytes = 0;
  unsigned ind = 0;
  while (totalBytes < (fMaxFluxFileMB * 1024 * 1024) && ind < allFiles.size()) {
    selected.push_back(allFiles[order[ind]]);
    totalBytes += allFiles[order[ind]].second;
    ind ++;
  } 

  // copy the files locally
  std::vector<std::pair<std::string, long>> localFiles = ifdhp->fetchSharedFiles(selected, fFluxCopyMethod);

  std::vector<std::string> files(localFiles.size());
  for (unsigned i = 0; i < localFiles.size(); i++) {
    files[i] = localFiles[i].first;
  }

  return files;
}

double NuMiKaonGen::LoadPOT() {
  TTreeReader metaReader(fMetaTreeName.c_str(), fFluxFile);
  TTreeReaderValue<double> pot(metaReader, "pots");
  
  double total_pot = 0.;

  while (metaReader.Next()) {
    total_pot += *pot;
  }

  return total_pot;
}

    
double NuMiKaonGen::GetPOT() {
  double ret = fAccumulatedPOT;
  fAccumulatedPOT = 0.;
  return ret;
}

const bsim::Dk2Nu *NuMiKaonGen::GetNextEntry() {
  // new file -- set the start entry 
  if (fNewFile) {
    // wrap file index around
    if (fFileIndex >= fFluxFiles.size()) {
      fFileIndex = 0;
    }
    // if (fFileIndex >= fFluxFiles.size()) {
    //  throw cet::exception("FluxReader Out of Files", 
    //                       "At file index (" + std::to_string(fFileIndex) + ") of available files (" + std::to_string(fFluxFiles.size()) + ").");
    // }

    std::cout << "New file: " << fFluxFiles[fFileIndex] << " at index: " << fFileIndex << " of: " << fFluxFiles.size() << std::endl;
    if (fFluxFile) delete fFluxFile;
    fFluxFile = new TFile(fFluxFiles[fFileIndex].c_str());
    fFluxTree = (TTree*)fFluxFile->Get(fTreeName.c_str());
    fFluxTree->SetBranchAddress("dk2nu",&fDk2Nu);

    // Start at a random index in this file
    fEntryStart = CLHEP::RandFlat::shootInt(fEngine, fFluxTree->GetEntries()-1);
    fEntry = fEntryStart;

    // load the POT in this file
    fThisFilePOT = LoadPOT();
    fNewFile = false;
  }
  else {
    fEntry = (fEntry + 1) % fFluxTree->GetEntries();
    // if this is the last entry, get ready for the next file
    if ((fEntry + 1) % fFluxTree->GetEntries() == fEntryStart) {
      fFileIndex ++;
      fNewFile = true;
    }
  }

  // count the POT
  fAccumulatedPOT += fThisFilePOT / fFluxTree->GetEntries();
    
  fFluxTree->GetEntry(fEntry);
  return fDk2Nu;
}

simb::MCFlux NuMiKaonGen::GetNext() {
  const bsim::Dk2Nu *flux = GetNextEntry();
  return MakeMCFlux(*flux);
}
  
simb::MCFlux NuMiKaonGen::MakeMCFlux(const bsim::Dk2Nu &dk2nu) {
  simb::MCFlux flux;
  
  flux.fFluxType = simb::kDk2Nu;
  flux.fntype    = dk2nu.decay.ntype;
  flux.fnimpwt   = dk2nu.decay.nimpwt;
  flux.fvx       = dk2nu.decay.vx;
  flux.fvy       = dk2nu.decay.vy;
  flux.fvz       = dk2nu.decay.vz;
  flux.fpdpx     = dk2nu.decay.pdpx;
  flux.fpdpy     = dk2nu.decay.pdpy;
  flux.fpdpz     = dk2nu.decay.pdpz;
  flux.fppdxdz   = dk2nu.decay.ppdxdz;
  flux.fppdydz   = dk2nu.decay.ppdydz;
  flux.fpppz     = dk2nu.decay.pppz;   
  flux.fppenergy  = dk2nu.decay.ppenergy;
  flux.fppmedium = dk2nu.decay.ppmedium;
  flux.fptype    = dk2nu.decay.ptype;
  flux.fndecay   = dk2nu.decay.ndecay;
  flux.fmuparpx  = dk2nu.decay.muparpx;
  flux.fmuparpy  = dk2nu.decay.muparpy;
  flux.fmuparpz  = dk2nu.decay.muparpz;
  flux.fmupare   = dk2nu.decay.mupare;
  flux.fnecm     = dk2nu.decay.necm;
  
  flux.fppvx = dk2nu.ppvx;
  flux.fppvy = dk2nu.ppvy;
  flux.fppvz = dk2nu.ppvz;
  
  flux.ftvx      = dk2nu.tgtexit.tvx;
  flux.ftvy      = dk2nu.tgtexit.tvy;
  flux.ftvz      = dk2nu.tgtexit.tvz;
  flux.ftpx      = dk2nu.tgtexit.tpx;
  flux.ftpy      = dk2nu.tgtexit.tpy;
  flux.ftpz      = dk2nu.tgtexit.tpz;
  flux.ftptype   = dk2nu.tgtexit.tptype;
  flux.ftgen     = dk2nu.tgtexit.tgen;
  
  flux.frun      = dk2nu.job;
  flux.fevtno    = dk2nu.potnum;
  flux.ftgptype  = dk2nu.ancestor[1].pdg;

  // flux.fnenergyn = flux.fnenergyf = enu;
  // flux.fnwtnear  = flux.fnwtfar = wgt; 
  // ignore variables dealing with the neutrino
  flux.fnenergyn = -1;
  flux.fnwtnear = flux.fnwtfar = -1;
  flux.fdk2gen   = -1;

  // placeholder for time
  flux.fxpoint = dk2nu.ancestor.back().startt;

  return flux;  
}

DEFINE_ART_CLASS_TOOL(NuMiKaonGen)

} // namespace ldm
} // namespace evgen
