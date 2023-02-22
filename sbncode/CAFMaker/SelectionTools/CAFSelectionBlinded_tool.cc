// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "sbncode/CAFMaker/ICAFSelectionTool.h"

namespace caf {

class CAFSelectionBlind: public ICAFSelectionTool {
public:

  CAFSelectionBlind(const fhicl::ParameterSet &p);
  ~CAFSelectionBlind() {}

  bool Select(StandardRecord &rec) override;
};

CAFSelectionBlind::CAFSelectionBlind(const fhicl::ParameterSet &p):
  ICAFSelectionTool(p)
{}

bool CAFSelectionBlind::Select(StandardRecord &rec) {
  // Set blinding on
  rec.hdr.isblind = true;

  // Configuration handles prescale
  return true;
}

DEFINE_ART_CLASS_TOOL(CAFSelectionBlind)

} // end namespace sbn
