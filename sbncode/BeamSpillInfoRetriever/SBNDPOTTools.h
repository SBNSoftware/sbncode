#ifndef _SBNDPOTTOOLS_H
#define _SBNDPOTTOOLS_H

#include <iostream>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ifdh_art/IFBeamService/IFBeam_service.h"
#include "ifbeam_c.h"

#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
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
  
  typedef struct MWRdata_t {
    std::vector< std::vector<double> > MWR_times;
    std::vector< std::vector< std::vector< int > > > unpacked_MWR;
  } MWRdata_t;

  PTBInfo_t extractPTBInfo(art::Handle<std::vector<artdaq::Fragment> > cont_frags, int HLT);
  double extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags);
  bool BrokenClock(double time, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp);
  sbn::BNBSpillInfo makeBNBSpillInfo(art::EventID const& eventID, double time, MWRdata_t const& MWRdata, std::vector<int> const& matched_MWR, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp);
}

#endif
