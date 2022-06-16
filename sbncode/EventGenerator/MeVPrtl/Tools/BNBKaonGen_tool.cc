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
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

// LArSoft includes
#include "Tools/Flux/GNuMIFlux.h"
#include "Tools/Flux/GSimpleNtpFlux.h"

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

    void GetNextEntry();
    std::vector<std::string> LoadFluxFiles();
    simb::MCFlux MakeMCFlux();
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
  TTree *fMetaTree;
  TFile *fFluxFile;

  genie::flux::GSimpleNtpEntry 	*fGSimpleEntry;
  genie::flux::GSimpleNtpNuMI 	*fGSimpleNuMI;
  genie::flux::GSimpleNtpAux 	*fGSimpleAux;
  genie::flux::GSimpleNtpMeta 	*fGSimpleMeta;

  // count POT
  double fAccumulatedPOT;
  double fThisFilePOT;

  // calculate Kaon total energy  
  double GetParentMass(const int pdg);
};

BNBKaonGen::BNBKaonGen(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("BNBKaonGen") 
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
  fMetaTree = NULL;
  fFluxFile = NULL;
  
  fGSimpleEntry = new genie::flux::GSimpleNtpEntry;
  fGSimpleNuMI  = new genie::flux::GSimpleNtpNuMI;
  fGSimpleAux   = new genie::flux::GSimpleNtpAux;
  fGSimpleMeta  = new genie::flux::GSimpleNtpMeta;

  fAccumulatedPOT = 0.;
  fThisFilePOT = 0.;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

BNBKaonGen::~BNBKaonGen()
{
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
  TTreeReaderValue<double> pot(metaReader, "protons");
  
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

void BNBKaonGen::GetNextEntry() {

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
    fMetaTree = (TTree*)fFluxFile->Get(fMetaTreeName.c_str());  
 
    fFluxTree->SetBranchAddress("entry",&fGSimpleEntry);
    fFluxTree->SetBranchAddress("numi" ,&fGSimpleNuMI);
    fFluxTree->SetBranchAddress("aux"  ,&fGSimpleAux);
    fMetaTree->SetBranchAddress("meta" ,&fGSimpleMeta);

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
}

simb::MCFlux BNBKaonGen::GetNext() {
  GetNextEntry(); 
  return MakeMCFlux();
}
  
simb::MCFlux BNBKaonGen::MakeMCFlux() {

  simb::MCFlux flux;
  fFluxTree->GetEntry(fEntry);

  //Assign flux type
  flux.fFluxType = simb::kSimple_Flux; 
 
  //stealing code from GenieHelper (nutools/EventGeneratorBase/GENIE/GENIEHelper.cxx)
  flux.fntype  = fGSimpleEntry->pdg;
  flux.fnimpwt = fGSimpleEntry->wgt;
  flux.fdk2gen = fGSimpleEntry->dist;
  flux.fnenergyf = fGSimpleEntry->E;

  if ( fGSimpleNuMI ) {
    flux.frun      = fGSimpleNuMI->run;
    flux.fevtno    = fGSimpleNuMI->evtno;
    flux.ftpx      = fGSimpleNuMI->tpx;
    flux.ftpy      = fGSimpleNuMI->tpy;
    flux.ftpz      = fGSimpleNuMI->tpz;
    flux.ftptype   = fGSimpleNuMI->tptype;   // converted to PDG
    flux.fvx       = fGSimpleNuMI->vx;
    flux.fvy       = fGSimpleNuMI->vy;
    flux.fvz       = fGSimpleNuMI->vz;

    flux.fndecay   = fGSimpleNuMI->ndecay;
    flux.fppmedium = fGSimpleNuMI->ppmedium;

    flux.fpdpx     = fGSimpleNuMI->pdpx;
    flux.fpdpy     = fGSimpleNuMI->pdpy;
    flux.fpdpz     = fGSimpleNuMI->pdpz;

    double apppz = fGSimpleNuMI->pppz;
    if ( TMath::Abs(fGSimpleNuMI->pppz) < 1.0e-30 ) apppz = 1.0e-30;
      flux.fppdxdz   = fGSimpleNuMI->pppx / apppz;
      flux.fppdydz   = fGSimpleNuMI->pppy / apppz;
      flux.fpppz     = fGSimpleNuMI->pppz;

      flux.fptype    = fGSimpleNuMI->ptype;
      flux.fppenergy = sqrt(	(fGSimpleNuMI->pppx)*(fGSimpleNuMI->pppx)
				+ (fGSimpleNuMI->pppy)*(fGSimpleNuMI->pppy)
				+ (fGSimpleNuMI->pppz)*(fGSimpleNuMI->pppz)
                                + GetParentMass(fGSimpleNuMI->ptype)*GetParentMass(fGSimpleNuMI->ptype));
  }

  // anything useful stuffed into vdbl or vint?
  // need to check the metadata  auxintname, auxdblname
  if ( fGSimpleAux && fGSimpleMeta ) {
    // references just for reducing complexity
    const std::vector<std::string>& auxdblname = fGSimpleMeta->auxdblname;
    const std::vector<std::string>& auxintname = fGSimpleMeta->auxintname;
    const std::vector<int>&    auxint = fGSimpleAux->auxint;
    const std::vector<double>& auxdbl = fGSimpleAux->auxdbl;

    for (size_t id=0; id<auxdblname.size(); ++id) {
      if ("muparpx"   == auxdblname[id]) flux.fmuparpx  = auxdbl[id];
      if ("muparpy"   == auxdblname[id]) flux.fmuparpy  = auxdbl[id];
      if ("muparpz"   == auxdblname[id]) flux.fmuparpz  = auxdbl[id];
      if ("mupare"    == auxdblname[id]) flux.fmupare   = auxdbl[id];
      if ("necm"      == auxdblname[id]) flux.fnecm     = auxdbl[id];
      if ("nimpwt"    == auxdblname[id]) flux.fnimpwt   = auxdbl[id];
      if ("fgXYWgt"   == auxdblname[id]) {
        flux.fnwtnear = flux.fnwtfar = auxdbl[id];
      }
    }
    for (size_t ii=0; ii<auxintname.size(); ++ii) {
      if ("tgen"      == auxintname[ii]) flux.ftgen     = auxint[ii];
      if ("tgptype"   == auxintname[ii]) flux.ftgptype  = auxint[ii];
    }
  }

  return flux;  
}

double BNBKaonGen::GetParentMass(const int pdg){
  
  double mass(0);

  switch (pdg) {
    case 130: /*K0L*/
      mass = Constants::Instance().klong_mass;
      return mass; 
    case 321:  /*K+*/
    case -321:  /*K-*/
      mass = Constants::Instance().kplus_mass;
      return mass;
    case 13:  /*mu-*/
    case -13:  /*mu+*/
      mass = Constants::Instance().muon_mass;
      return mass;
    case 211: /*pi+*/
    case -211: /*pi-*/
      mass = Constants::Instance().piplus_mass;
      return mass;
    default:
      return mass;
  } 
}
DEFINE_ART_CLASS_TOOL(BNBKaonGen)

} // namespace ldm
} // namespace evgen
