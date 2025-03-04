#ifndef _SBNDPOTTOOLS_H
#define _SBNDPOTTOOLS_H

#include <iostream>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/TDCTimestampFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"


#include <string.h>
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <cxxabi.h>
#include <memory>
#include <bitset>
#include <tuple>
#include <algorithm>

#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/TDCTimestampFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

namespace sbn{
  typedef struct PTBInfo_t {
    double currPTBTimeStamp  = 1e20;
    double prevPTBTimeStamp  = 0;
    unsigned int GateCounter = 0;
  } PTBInfo_t;

  typedef struct TriggerInfo_t {
    double t_current_event  = 0;
    double t_previous_event = 0;
    unsigned int number_of_gates_since_previous_event = 0;
  } TriggerInfo_t;

  PTBInfo_t extractPTBInfo(art::Handle<std::vector<artdaq::Fragment> > cont_frags, int HLT);
  double extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags);
  bool BrokenClock(double time);
}

#endif
