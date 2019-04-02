#ifndef __sbnanalysis_core_SubRun__
#define __sbnanalysis_core_SubRun__

/**
 * \file SubRun.hh
 *
 * Subrun-level information.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2019/03/01
 */

#include <map>
#include <string>
#include <vector>
#include <TTree.h>
#include <TVector3.h>
#include <sbncode/sbnanalysis/core/DataTypes/Experiment.hh>

namespace sbnanalysis {

/**
 * \class SubRun
 * \brief The standard subrun data definition.
 */
class SubRun {
public:
  SubRun()
      : runID(0), subrunID(0),
        totpot(0), totgoodpot(0),
        totspills(0), goodspills(0) {}

  SubRun(int _runID, int _subrunID,
         double _totpot, double _totgoodpot,
         int _totspills, int _goodspills)
      : runID(_runID), subrunID(_subrunID),
        totpot(_totpot), totgoodpot(_totgoodpot),
        totspills(_totspills), goodspills(_goodspills) {}

  int runID;
  int subrunID;
  double totpot;
  double totgoodpot;
  int totspills;
  int goodspills;
};

}  // namespace sbnanalysis

#endif  // __sbnanalysis_core_SubRun__

