#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"

#include <fstream>
#include <vector>
#include <map>

#include "FluxInterface.h"

class TH1D;
class TTree;
class TFile;

namespace fluxr {
  class FluxReader {
  public:
    // Required constructor
    FluxReader(fhicl::ParameterSet const &pset,
               art::ProductRegistryHelper &helper,
               art::SourceHelper const &pm);

    // Required by FileReaderSource:
    void closeCurrentFile();
    void readFile(std::string const &name,
                  art::FileBlock* &fb);
    bool readNext(art::RunPrincipal* const &inR,
                  art::SubRunPrincipal* const &inSR,
                  art::RunPrincipal* &outR,
                  art::SubRunPrincipal* &outSR,
                  art::EventPrincipal* &outE);

  private:
    art::SourceHelper const      &fSourceHelper;
    art::SubRunID                 fSubRunID;

    uint32_t                      fEventCounter;
    uint32_t                      fEntry;
    int                           fMaxEvents;    // fhicl parameter.  Maximum number of events.
    uint32_t                      fSkipEvents;   // fhicl parameter.  Number of events to skip.
    std::string                   fInputType;    // fhicl parameter.  Maximum number of events.
    float                         fPOT;
    float                         fCurrentPOT;
    bool                          fSelfIncrementRuns; // fhicl parameter. If true it increments run numbers sequentially
    art::RunNumber_t              fIncrement = 1;

    FluxInterface*                fFluxDriver;
    TFile*                        fFluxInputFile;
    int                           fNuPdgCode[4];
    TH1D*                         fHFlux[4];
    TH1D*                         fHFluxParent[4][4];
    TH1D*                         fHFluxSec[4][5];

    fhicl::ParameterSet           fConfigPS;

    art::TypeLabel                fTLmctruth;
    art::TypeLabel                fTLmcflux;
    art::TypeLabel                fTLdk2nu;
    art::TypeLabel                fTLnuchoice;
  };  // FluxReader
}


