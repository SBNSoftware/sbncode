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

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  TxtFileGen class definiton
 */
class TxtFileGen : public IMesonGen
{
public:
    /**
     *  @brief  Constructor
     */
    TxtFileGen(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~TxtFileGen();

    double GetPOT() override;
    simb::MCFlux GetNext() override;

    void configure(const fhicl::ParameterSet&) override;

    std::string GetNextEntry();
    std::vector<std::string> LoadFluxFiles();
    simb::MCFlux MakeMCFlux(std::string line);

    // no weights
    double MaxWeight() override { return -1.; }

private:
  // config
  std::string fSearchPath;
  std::vector<std::string> fSearchPatterns;
  unsigned long fMaxFluxFileMB;
  std::string fFluxCopyMethod;
  bool fRandomizeFiles;
  int fNSkipLines;
  double fPOTPerMeson;
  int fPDGCode;

  // info for tracking files
  unsigned fFileIndex;
  bool fNewFile;
  std::vector<std::string> fFluxFiles;
  std::ifstream fCurrentFile;
  int fCurrentFileLines;

  // info for tracking entry in file
  unsigned fEntry;
  unsigned fEntryStart;

  // count POT
  double fAccumulatedPOT;
};

TxtFileGen::TxtFileGen(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("TxtFileGen") 
{
  configure(pset);

  // copy the flux files locally
  fFluxFiles = LoadFluxFiles();

  // setup indices
  fFileIndex = 0;
  fEntry = 0;
  fEntryStart = 0;
  fNewFile = true;
  fCurrentFileLines = 0;

  fAccumulatedPOT = 0.;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------
TxtFileGen::~TxtFileGen() {}

//------------------------------------------------------------------------------------------------------------------------------------------
void TxtFileGen::configure(fhicl::ParameterSet const &pset)
{
  fSearchPath = pset.get<std::string>("SearchPath");
  fSearchPatterns = pset.get<std::vector<std::string>>("FluxFiles");
  fMaxFluxFileMB = pset.get<unsigned long>("MaxFluxFileMB", 2 * 1024);
  fFluxCopyMethod = pset.get<std::string>("FluxCopyMethod", "IFDH");
  fRandomizeFiles = pset.get<bool>("RandomizeFiles");
  fNSkipLines = pset.get<int>("NSkipLines");
  fPOTPerMeson = pset.get<double>("POTPerMeson");
  fPDGCode = pset.get<int>("PDGCode");

  std::cout << "Searching for flux files at path: " << fSearchPath << std::endl;
  std::cout << "With patterns:\n";
  for (const std::string &s: fSearchPatterns) std::cout << s << std::endl;
  std::cout << "With copy method: " << fFluxCopyMethod << std::endl;
}

std::vector<std::string> TxtFileGen::LoadFluxFiles() {
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

double TxtFileGen::GetPOT() {
  double ret = fAccumulatedPOT;
  fAccumulatedPOT = 0.;
  return ret;
}

std::string TxtFileGen::GetNextEntry() {
  // new file -- set the start entry 
  if (fNewFile) {
    // wrap file index around
    if (fFileIndex >= fFluxFiles.size()) {
      fFileIndex = 0;
    }

    std::cout << "New file: " << fFluxFiles[fFileIndex] << " at index: " << fFileIndex << " of: " << fFluxFiles.size() << std::endl;
    fCurrentFile.open(fFluxFiles[fFileIndex]);

    // count the number of lines
    std::string line;
    int lines;
    for(lines = 0; std::getline(fCurrentFile,line); lines++);
    fCurrentFileLines = lines - fNSkipLines;

    // Seek back to start
    fCurrentFile.clear();
    fCurrentFile.seekg(0);

    // Start at a random index in this file
    fEntryStart = CLHEP::RandFlat::shootInt(fEngine, fCurrentFileLines);
    fEntry = fEntryStart;

    // skip to that line in the file
    for (lines = 0; lines < fNSkipLines + (int)fEntryStart && std::getline(fCurrentFile,line); lines++);

    fNewFile = false;
  }
  else {
    fEntry = (fEntry + 1) % fCurrentFileLines;

    // if this is the last entry, get ready for the next file
    if ((fEntry + 1) % fCurrentFileLines == fEntryStart) {
      fFileIndex ++;
      fNewFile = true;
    }
  }

  // Read the next line, if we are at the end of the file go back to the beginning
  std::string thisline;
  if (!std::getline(fCurrentFile, thisline)) {

    // seek back to start
    fCurrentFile.clear();
    fCurrentFile.seekg(0);

    // skip to the first valid line
    int lines;
    std::string line;
    for (lines = 0; lines < fNSkipLines && std::getline(fCurrentFile,line); lines++);
  }

  // count the POT
  fAccumulatedPOT += fPOTPerMeson;

  return thisline;
}

simb::MCFlux TxtFileGen::GetNext() {
  std::string meson = GetNextEntry();
  return MakeMCFlux(meson);
}
  
simb::MCFlux TxtFileGen::MakeMCFlux(std::string line) {
  std::cout << "Parsing line: " << line << std::endl;

  // read the 4-vector off the line
  double E, px, py, pz;
  std::stringstream l4vec(line);
  l4vec >> px >> py >> pz >> E;

  std::cout << "Values: " << px << " " << py << " " << pz << " " << E << std::endl;

  simb::MCFlux flux;
  
  flux.fFluxType = simb::kSimple_Flux; // number for file gen....
  flux.fnimpwt   = 1; // All the same weight!
  flux.fvx       = 0; // At the target -- TODO is this ok???
  flux.fvy       = 0;
  flux.fvz       = 0;
  flux.fpdpx     = px;
  flux.fpdpy     = py;
  flux.fpdpz     = pz;
  flux.fppdxdz   = px/pz;
  flux.fppdydz   = py/pz;
  flux.fpppz     = pz;
  flux.fppenergy  = E;
  flux.fppmedium = -1;
  flux.fptype    = fPDGCode;
  flux.fndecay   = -1;
  flux.fmuparpx  = -1;
  flux.fmuparpy  = -1;
  flux.fmuparpz  = -1;
  flux.fmupare   = -1;
  flux.fnecm     = -1;
  
  flux.fppvx = 0.;
  flux.fppvy = 0.;
  flux.fppvz = 0.;
  
  flux.ftvx      = -1;
  flux.ftvy      = -1;
  flux.ftvz      = -1;
  flux.ftpx      = -1;
  flux.ftpy      = -1;
  flux.ftpz      = -1;
  flux.ftptype   = -1;
  flux.ftgen     = -1;
  
  flux.frun      = -1;
  flux.fevtno    = -1;
  flux.ftgptype  = -1;

  // flux.fnenergyn = flux.fnenergyf = enu;
  // flux.fnwtnear  = flux.fnwtfar = wgt; 
  // ignore variables dealing with the neutrino
  flux.fnenergyn = -1;
  flux.fnwtnear = flux.fnwtfar = -1;
  flux.fdk2gen   = -1;

  // placeholder for time
  flux.fxpoint = 0.;

  return flux;  
}

DEFINE_ART_CLASS_TOOL(TxtFileGen)

} // namespace ldm
} // namespace evgen
