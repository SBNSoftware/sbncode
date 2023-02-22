#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "fhiclcpp/ParameterSet.h"

#include "TRandomGen.h"
#include "TTree.h"

#include <stdlib.h>
#include <time.h> 

namespace caf {

class ICAFSelectionTool {
public:
  ICAFSelectionTool(const fhicl::ParameterSet &p):
    fExtension(p.get<std::string>("Extension")),
    fInvert(p.get<bool>("Invert", false)),
    fPreScale(p.get<float>("PreScale", -1)),
    fRandom(TRandomMT64(art::ServiceHandle<rndm::NuRandomService>()->getSeed())),
    fFirstInFile(false),
    fFirstSelectedInFile(false),
    fNEvents(0)
{
  // zero-initialize counters
  fFirstHeader.nbnbinfo = 0;
  fFirstHeader.nnumiinfo = 0;
  fFirstHeader.ngenevt = 0;
  fFirstHeader.pot = 0;
}

  virtual ~ICAFSelectionTool() noexcept = default;

  /// For external modules to call: run the actual selection
  bool DoSelect(StandardRecord &rec) {
    // Increment counters in header info
    if (fFirstInFile) {
      fFirstHeader.nbnbinfo += rec.hdr.nbnbinfo;
      fFirstHeader.nnumiinfo += rec.hdr.nnumiinfo;
      fFirstHeader.ngenevt += rec.hdr.ngenevt;
      fFirstHeader.pot += rec.hdr.pot;

      fFirstInFile = false;
    }

    // Implement prescale
    bool pass_prescale = (fPreScale < 0) || (fRandom.Uniform() < 1. / fPreScale);

    bool selected = Select(rec) == !fInvert /* implement inversion */;
    selected = selected && pass_prescale;

    if (selected) {
      // set flags
      rec.hdr.first_in_file = fFirstSelectedInFile;
      rec.hdr.first_in_subrun = fFirstSelectedInFile;

      // fix POT
      float prescale_factor = (fPreScale > 0) ? fPreScale : 1.;
      if (fFirstSelectedInFile) { 
        rec.hdr.nbnbinfo = fFirstHeader.nbnbinfo / prescale_factor;
        rec.hdr.nnumiinfo = fFirstHeader.nnumiinfo / prescale_factor;
        rec.hdr.ngenevt = fFirstHeader.ngenevt / prescale_factor;

        // scale the POT by the GetPrescale() function
        // so implementations can override it
        rec.hdr.pot = fFirstHeader.pot / GetPrescale(); 

        // reset counters in header info
        fFirstHeader.nbnbinfo = 0;
        fFirstHeader.nnumiinfo = 0;
        fFirstHeader.ngenevt = 0;
        fFirstHeader.pot = 0;
      }

      // bookkeeping
      fFirstSelectedInFile = false;
      fNEvents ++;
    }

    return selected;
  }

  // bookkeeping
  void NextFile() {
    fFirstInFile = true;
    fFirstSelectedInFile = true;

  }

  int GetTotalEvents() const { return fNEvents; }
  std::string Extension() const { return fExtension; }

  virtual float GetPrescale() const { return fPreScale; }

protected:
  /// For children to implement: Whether to select a given track
  ///
  /// The StandardRecord object is also edit-able. This can enable
  /// (e.g.) blinding certain parameters
  virtual bool Select(StandardRecord &rec) = 0;

  // config
  std::string fExtension;
  bool fInvert;
  float fPreScale;

  // random numbers
  TRandom fRandom;

  // bookkeeping
  bool fFirstInFile;
  bool fFirstSelectedInFile;
  SRHeader fFirstHeader;
  int fNEvents;
};

} // end namespace caf
