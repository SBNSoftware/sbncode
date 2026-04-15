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
 *  @brief  TxtFileGenNM class definiton
 */
class TxtFileGenNM : public IMesonGen
{
public:
    /**
     *  @brief  Constructor
     */
    TxtFileGenNM(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~TxtFileGenNM();

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
  bool fVerbose;

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

TxtFileGenNM::TxtFileGenNM(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("TxtFileGenNM") 
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
TxtFileGenNM::~TxtFileGenNM() {}

//------------------------------------------------------------------------------------------------------------------------------------------
void TxtFileGenNM::configure(fhicl::ParameterSet const &pset)
{
  fSearchPath = pset.get<std::string>("SearchPath");
  fSearchPatterns = pset.get<std::vector<std::string>>("FluxFiles");
  fMaxFluxFileMB = pset.get<unsigned long>("MaxFluxFileMB", 2 * 1024);
  fFluxCopyMethod = pset.get<std::string>("FluxCopyMethod", "IFDH");
  fRandomizeFiles = pset.get<bool>("RandomizeFiles");
  fNSkipLines = pset.get<int>("NSkipLines");
  fVerbose = pset.get<bool>("Verbose", false);

  std::cout << "Searching for flux files at path: " << fSearchPath << std::endl;
  std::cout << "With patterns:\n";
  for (const std::string &s: fSearchPatterns) std::cout << s << std::endl;
  std::cout << "With copy method: " << fFluxCopyMethod << std::endl;
}

std::vector<std::string> TxtFileGenNM::LoadFluxFiles() {
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

double TxtFileGenNM::GetPOT() {
  double ret = fAccumulatedPOT;
  fAccumulatedPOT = 0.;
  return ret;
}


  
std::string TxtFileGenNM::GetNextEntry() {
  // new file -- set the start entry 
  if (fNewFile) {
    // wrap file index around
    if (fFileIndex >= fFluxFiles.size()) {
      fFileIndex = 0;
    }

    if (fVerbose) std::cout << "New file: " << fFluxFiles[fFileIndex] << " at index: " << fFileIndex << " of: " << fFluxFiles.size() << std::endl;

    if (fCurrentFile.is_open()) fCurrentFile.close();
    fCurrentFile.clear();
    fCurrentFile.open(fFluxFiles[fFileIndex]);

    if (!fCurrentFile.is_open()) {
      throw cet::exception("TxtFileGenNM")
        << "Could not open file " << fFluxFiles[fFileIndex];
    }

    // count the number of lines
    std::string line;
    int lines;
    for(lines = 0; std::getline(fCurrentFile,line); lines++);

    fCurrentFileLines = lines - fNSkipLines;
    if (fCurrentFileLines <= 0) {
      throw cet::exception("TxtFileGenNM")
	<< "File " << fFluxFiles[fFileIndex]
	<< " has no usable data lines after skipping "
	<< fNSkipLines << " header lines. "
	<< "Total lines in file: " << lines;
    }
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
  /*  std::string thisline;
  if (!std::getline(fCurrentFile, thisline)) {

    // seek back to start
    fCurrentFile.clear();
    fCurrentFile.seekg(0);

    // skip to the first valid line
    int lines;
    std::string line;
    for (lines = 0; lines < fNSkipLines && std::getline(fCurrentFile,line); lines++);
  }

  return thisline;*/

  std::string thisline;
  if (!std::getline(fCurrentFile, thisline)) {

    fCurrentFile.clear();
    fCurrentFile.seekg(0);

    int lines;
    std::string line;
    for (lines = 0; lines < fNSkipLines && std::getline(fCurrentFile, line); ++lines) {}

    if (!std::getline(fCurrentFile, thisline)) {
      throw cet::exception("TxtFileGenNM")
	<< "Failed to read a valid data line from file "
	<< fFluxFiles[fFileIndex]
	<< " after rewinding.";
    }
  }
  return thisline;
}



simb::MCFlux TxtFileGenNM::GetNext() {
  std::string meson = GetNextEntry();
  return MakeMCFlux(meson);
}
  
simb::MCFlux TxtFileGenNM::MakeMCFlux(std::string line) {
  // std::cout << "Parsing line: " << line << std::endl;

  // read the 4-vector off the line
  double  px, py, pz,x, y, z, wgt, start_t;
  //  int pdgcode;
  std::stringstream l4vec(line);
  l4vec >> px >> py >> pz >> x >> y >> z >> start_t >> wgt;

  if (fVerbose) std::cout << "Values: " <<  px << py << pz << x << y << z << start_t << wgt << std::endl;

  simb::MCFlux flux;
  
  flux.fFluxType = simb::kSimple_Flux; // number for file gen....
  flux.fnimpwt   = wgt;
  flux.fvx       = x; // At the target -- TODO is this ok???
  flux.fvy       = y;
  flux.fvz       = z;
  flux.fpdpx     = px;
  flux.fpdpy     = py;
  flux.fpdpz     = pz;
  flux.fppdxdz   = px/pz;
  flux.fppdydz   = py/pz;
  flux.fpppz     = pz;

  double mpi0 = 0.1349768; 
  double p2 = px*px + py*py + pz*pz;
  flux.fppenergy  = std::sqrt(p2 + mpi0*mpi0);

  flux.fppmedium = -1;
  flux.fptype    = 111;
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
  flux.ftpx      = start_t;
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

  // Update the POT
  fAccumulatedPOT = 1000000;

  return flux;  
}

DEFINE_ART_CLASS_TOOL(TxtFileGenNM)

} // namespace ldm
} // namespace evgen
