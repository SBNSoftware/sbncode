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
#include "ifdh_art/IFDHService/IFDH_service.h"

#include "nusimdata/SimulationBase/MCFlux.h"

// local includes
#include "IMesonGen.h"
#include "boone.h"
#include "PDGCodes.h"

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
 *  @brief  BNBKaonGen class definiton
 */
class BNBKaonGen : public IMesonGen
{
public:
    /**
     *  @brief  Constructor
     */
    BNBKaonGen(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~BNBKaonGen();

    double GetPOT() override;
    simb::MCFlux GetNext() override;

    void configure(const fhicl::ParameterSet&) override;

    // const bsim::Dk2Nu *GetNextEntry();override;
    const bsim::BooNe *GetNextEntry();
    std::vector<std::string> LoadFluxFiles();
    simb::MCFlux MakeMCFlux(const bsim::Dk2Nu &dk2nu);
    simb::MCFlux MakeMCFlux(const bsim::BooNe &boone);
    double LoadPOT();

    // no weights
    double MaxWeight() override { return -1.; }

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
  bsim::BooNe *fBooNe;

  // count POT
  double fAccumulatedPOT;
  double fThisFilePOT;
  
};

BNBKaonGen::BNBKaonGen(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("BNBKaonGen") 
{
  configure(pset);

  // setup indices
  fFileIndex = 0;
  fEntry = 0;
  fEntryStart = 0;
  fNewFile = true;
  fFluxTree = NULL;
  fFluxFile = NULL;
  fDk2Nu = new bsim::Dk2Nu;
  fBooNe = new bsim::BooNe;

  fAccumulatedPOT = 0.;
  fThisFilePOT = 0.;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

BNBKaonGen::~BNBKaonGen()
{

  if (fDk2Nu) delete fDk2Nu;
  if (fBooNe) delete fBooNe;
}

//------------------------------------------------------------------------------------------------------------------------------------------
void BNBKaonGen::configure(fhicl::ParameterSet const &pset)
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
  fFluxFiles = pset.get<std::vector<std::string>>("FluxFilesFullPath");
}

std::vector<std::string> BNBKaonGen::LoadFluxFiles() {
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

double BNBKaonGen::LoadPOT() {
  TTreeReader metaReader(fMetaTreeName.c_str(), fFluxFile);
  TTreeReaderValue<double> pot(metaReader, "pots");
  
  double total_pot = 0.;

  while (metaReader.Next()) {
    total_pot += *pot;
  }

  return total_pot;
}

    
double BNBKaonGen::GetPOT() {
  double ret = fAccumulatedPOT;
  fAccumulatedPOT = 0.;
  return ret;
}

const bsim::BooNe *BNBKaonGen::GetNextEntry() {
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
    // fFluxTree->SetBranchAddress("dk2nu",&fDk2Nu);
    fBooNe = new bsim::BooNe(fFluxFiles[fFileIndex].c_str());

    // Start at a random index in this file
    fEntryStart = CLHEP::RandFlat::shootInt(fEngine, fFluxTree->GetEntries()-1);
    fEntry = fEntryStart;

    // load the POT in this file
    // fThisFilePOT = LoadPOT();
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
  fBooNe->myNtuple.run    = fEntry;
  fBooNe->myNtuple.eventn = fFileIndex;

  fAccumulatedPOT += fBooNe->GetPOT();
    
  fFluxTree->GetEntry(fEntry);
  fBooNe->GetEntry(fEntry);
  return fBooNe;
}

simb::MCFlux BNBKaonGen::GetNext() {
  // const bsim::Dk2Nu *flux = GetNextEntry();
  const bsim::BooNe *flux = GetNextEntry();
  return MakeMCFlux(*flux);
}
  
simb::MCFlux BNBKaonGen::MakeMCFlux(const bsim::BooNe &boone) {

  simb::MCFlux flux;
  auto fBooneNtp=boone.myNtuple;

  if ( fBooneNtp.ntp == 1 ) {
    flux.fntype = 12; //nue
  }
  else if ( fBooneNtp.ntp == 2 ) {
    flux.fntype = -12; //nuebar
  }
  else if ( fBooneNtp.ntp == 3 ) {
    flux.fntype = 14; //numu
  }
  else if ( fBooneNtp.ntp == 4 ) {
    flux.fntype = -14; //numubar
  }
  else{
    mf::LogWarning("BooNEInterface") << "Neutrino type not recognized! ntp = " << fBooneNtp.ntp
                                      << std::endl;
  }

  flux.fFluxType = simb::kDk2Nu;
  flux.fnimpwt   = fBooneNtp.beamwgt;
  flux.fvx       = fBooneNtp.ini_pos[0][0]; //0
  flux.fvy       = fBooneNtp.ini_pos[0][1]; //0
  flux.fvz       = fBooneNtp.ini_pos[0][2]; //0
  flux.fpdpx     = fBooneNtp.fin_mom[1][0]; //1 final
  flux.fpdpy     = fBooneNtp.fin_mom[1][1]; //1 final
  flux.fpdpz     = fBooneNtp.fin_mom[1][2]; //1 final
  flux.fpppz     = fBooneNtp.ini_mom[1][2]; //1 init
  
  // Momentum projected in dx/dz or dy/dz
  flux.fpppz     = fBooneNtp.ini_mom[1][2]; //1 init
  double pppx    = fBooneNtp.ini_mom[1][0]; //1 init
  double pppy    = fBooneNtp.ini_mom[1][1]; //1 init
  double apppz = flux.fpppz;
  if (TMath::Abs(flux.fpppz) < 1.0e-30) apppz = 1.0e-30;
  flux.fppdxdz   = pppx / apppz;
  flux.fppdydz   = pppy / apppz;

  flux.fppmedium = 0.;
  flux.fppenergy  = fBooneNtp.ini_eng[1];

  int npart      = fBooneNtp.npart;
  int ptype_input = fBooneNtp.id[1];
  int tptype_input = fBooneNtp.id[npart-2];

  if (ptype_input != 0) ptype_input =  genie::GeantToPdg(ptype_input);
  if (tptype_input != 0) tptype_input= genie::GeantToPdg(tptype_input);

  flux.fptype    = ptype_input;
  flux.ftptype   = tptype_input;

  /////
  //  Now need to calculate ndecay
  /////

  double Nenergy = fBooneNtp.ini_eng[0];
  double Ndxdz   = fBooneNtp.ini_mom[0][0] / fBooneNtp.ini_mom[0][2];
  double Ndydz   = fBooneNtp.ini_mom[0][1] / fBooneNtp.ini_mom[0][2];

  double ppenergy = fBooneNtp.ini_eng[1];
  double pdPx     = fBooneNtp.fin_mom[1][0];
  double pdPy     = fBooneNtp.fin_mom[1][1];
  double pdPz     = fBooneNtp.fin_mom[1][2];

  double ppdxdz   = fBooneNtp.ini_mom[1][0] / fBooneNtp.ini_mom[1][2];
  double ppdydz   = fBooneNtp.ini_mom[1][1] / fBooneNtp.ini_mom[1][2];
  double pppz     = fBooneNtp.ini_mom[1][2];

  // Get the neutrino energy in the parent decay cm
  double parent_mass = std::sqrt(ppenergy * ppenergy -
                                  pppz * pppz * (ppdxdz * ppdxdz +
                                  ppdydz * ppdydz +
                                  1.));

  double parent_energy = std::sqrt(pdPx * pdPx +
                                    pdPy * pdPy +
                                    pdPz * pdPz +
                                    parent_mass * parent_mass);

  double gamma = parent_energy / parent_mass;
  double beta[3];
  beta[0] = pdPx / parent_energy;
  beta[1] = pdPy / parent_energy;
  beta[2] = pdPz / parent_energy;

  double partial = fBooneNtp.ini_mom[0][2] * gamma *
                    (beta[0] * Ndxdz +
                    beta[1] * Ndydz +
                    beta[2]);

  double Necm = gamma * Nenergy - partial;

  if (fBooneNtp.id[1] == 10 && fBooneNtp.ntp == 1) flux.fndecay = 1;
  else if (fBooneNtp.id[1] == 10 && fBooneNtp.ntp == 2) flux.fndecay = 2;
  else if (fBooneNtp.id[1] == 10 && fBooneNtp.ntp == 3) flux.fndecay = 3;
  else if (fBooneNtp.id[1] == 10 && fBooneNtp.ntp == 4) flux.fndecay = 4;
  else if (fBooneNtp.id[1] == 11 && fBooneNtp.ntp == 3) {
    //check if it is a two or three body decay
    if (fabs((parent_mass*parent_mass-0.105658389*0.105658389)/(2.*parent_mass)-Necm)/Necm <= 0.001)
      //two body decay (numu + mu+)
      flux.fndecay = 5;
    else {
      //three body decay (numu + pi0 + mu+)
      flux.fndecay = 7;
    }
  }
  else if (fBooneNtp.id[1] == 11 && fBooneNtp.ntp == 1) flux.fndecay = 6;
  else if (fBooneNtp.id[1] == 12 && fBooneNtp.ntp == 4) {
    if (fabs((parent_mass*parent_mass-0.105658389*0.105658389)/(2.*parent_mass)-Necm)/Necm <= 0.001) {
      //two body decay (numu + mu+)
      flux.fndecay = 8;
    }
    else {
      //three body decay (numu + pi0 + mu+)
      flux.fndecay = 10;
    }
  }
  else if (fBooneNtp.id[1] == 12 && fBooneNtp.ntp == 2) flux.fndecay = 9;
  else if (fBooneNtp.id[1] == 5 ) flux.fndecay = 11;
  else if (fBooneNtp.id[1] == 6 ) flux.fndecay = 12;
  else if (fBooneNtp.id[1] == 8 ) flux.fndecay = 13;
  else if (fBooneNtp.id[1] == 9 ) flux.fndecay = 14;

  /////
  //  End calculation of ndecay
  /////
  double mupare;
  double muparpx;
  double muparpy;
  double muparpz;

  if ( fBooneNtp.id[1] == 5 || fBooneNtp.id[1] == 6) {
    mupare  = fBooneNtp.ini_eng[2];
    muparpx = fBooneNtp.fin_mom[2][0];
    muparpy = fBooneNtp.fin_mom[2][1];
    muparpz = fBooneNtp.fin_mom[2][2];
  } else {
    mupare  = -9999.;
    muparpx = -9999.;
    muparpy = -9999.;
    muparpz = -9999.;
  }

  flux.fmuparpx = muparpx;
  flux.fmuparpy = muparpy;
  flux.fmuparpz = muparpz;
  flux.fmupare = mupare;

  flux.fnecm     = Necm;
  
  flux.fppvx = fBooneNtp.ini_pos[1][0];
  flux.fppvy = fBooneNtp.ini_pos[1][1];
  flux.fppvz = fBooneNtp.ini_pos[1][2];
  
  flux.ftvx      = fBooneNtp.ini_pos[npart-2][0];
  flux.ftvy      = fBooneNtp.ini_pos[npart-2][1];
  flux.ftvz      = fBooneNtp.ini_pos[npart-2][2];
  flux.ftpx      = fBooneNtp.ini_mom[npart-2][0];
  flux.ftpy      = fBooneNtp.ini_mom[npart-2][1];
  flux.ftpz      = fBooneNtp.ini_mom[npart-2][2];
  // flux.ftgen     = dk2nu.tgtexit.tgen;
  
  flux.frun      = fBooneNtp.run;
  flux.fevtno    = fBooneNtp.eventn;
  flux.ftgptype  = genie::GeantToPdg(fBooneNtp.id[npart-2]);

  // flux.fnenergyn = flux.fnenergyf = enu;
  // flux.fnwtnear  = flux.fnwtfar = wgt; 
  // ignore variables dealing with the neutrino
  flux.fnenergyn = -1;
  flux.fnwtnear = flux.fnwtfar = -1;
  flux.fdk2gen   = -1;

  // // placeholder for time
  flux.fxpoint = fBooneNtp.ini_t[0];

  return flux;  
}

  

DEFINE_ART_CLASS_TOOL(BNBKaonGen)

} // namespace ldm
} // namespace evgen
