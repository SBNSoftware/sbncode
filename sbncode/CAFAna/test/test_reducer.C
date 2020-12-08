#include "CAFAna/Core/FileReducer.h"

#include "SBNAna/Cuts/Cuts.h"

using namespace ana;

void test_reducer(std::string fin,
                  std::string fout = "reduced.root")
{
  FileReducer reducer(fin, fout);

  // Only keep spills passing flash trigger
  reducer.AddSpillCut(kFlashTrigger);
  // And within those, only slices passing this cut
  reducer.AddSliceCut(kSlcNuScoreCut);

  // And when we do keep them, remove their true particle list
  reducer.AddReductionStep(ClearTrueParticles);

  reducer.Go();
}
