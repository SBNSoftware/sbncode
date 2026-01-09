#ifndef _SBNDPOTTOOLS_H
#define _SBNDPOTTOOLS_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ifdh_art/IFBeamService/IFBeam_service.h"

#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/TDCTimestampFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "sbncode/BeamSpillInfoRetriever/MWRData.h"
#include "larcorealg/CoreUtils/counter.h"
#include <vector>

namespace sbn::pot{

  typedef struct PTBInfo_t {
    std::uint64_t currPTBTimeStamp  = UINT64_MAX;  ///< Timestamp in UTC nanoseconds since Unix epoch (converted from 20ns clock ticks)
    std::uint64_t prevPTBTimeStamp  = 0; 
    unsigned int GateCounter = 0;
    bool isHLT = false;            
    uint64_t triggerWord = 0;  ///< Timestamp in s since beam extraction signal
  } PTBInfo_t;

  typedef struct TriggerInfo_t {
    int gate_type = 0; ///< Source of the spill: `1`: BNB, `2`: NuMI
    double t_current_event  = 0;  ///< Timestamp in UTC seconds since Unix epoch (converted from 20ns clock ticks)
    double t_previous_event = 0;
    unsigned int number_of_gates_since_previous_event = 0;
    std::int64_t WR_to_Spill_conversion = 0;
  } TriggerInfo_t;
  
  typedef struct MWRdata_t {
    std::vector< std::vector<double> > MWR_times;
    std::vector< std::vector< std::vector< int > > > unpacked_MWR;
  } MWRdata_t;

  /**
   * @brief Extracts information from PTB for a single HLT for use in SBND POT accounting.
   * 
   * @param cont_frags The PTB fragments to examine.
   * @param HLT The high level trigger we are searching for.
   */
  PTBInfo_t extractPTBInfo(art::Handle<std::vector<artdaq::Fragment> > cont_frags, int HLT);

  /**
   * @brief Extracts ALL PTB information for using in SBND CAF files.
   * 
   * @param cont_frags The PTB fragments to examine.
   * @return Vector of PTBInfo_t containing all triggers found.
   */
  std::vector<PTBInfo_t> extractAllPTBInfo(std::vector<artdaq::Fragment> const& cont_frags);

  /**
   * @brief Extracts information from TDC for use in SBND POT accounting.
   * 
   * @param cont_frags The TDC fragments to examine.
   */
  double extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags);

  /**
   * @brief Checks for a known IFBeam issue where TOR clocks are incorrect.
   * @return true if a time in IFBeam corresponds to a broken clock false if
   * it is a real spill. 
   * @details Check validity by seeing if all devices fire at once or
   * if time is just a single outlier. 
   * @param time The time for which we want to check for a spill.
   * @param bfp beamfolder object to check. 
   */
  bool BrokenClock(double time, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp);

  /**
   * @brief extract spill times from the IFBeam data base. Used in SBND and ICARUS. 
   * 
   * @param triggerInfo triggerInfo object to match
   * @param bfp beamfolder to match
   * @param bfp_mwr mwr devices have separate daq system, much treat separately
   * @param fTimePad Padding to use in ifbeam queries
   * @param MWRtoroidDelay Delay between MWR and toroids
   * @param mwrdata MWRData to unpack 
   */
  MWRdata_t extractSpillTimes(TriggerInfo_t const& triggerInfo, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp_mwr, double fTimePad, double MWRtoroidDelay, sbn::MWRData mwrdata );
  /**
   * @brief Compile spill information into BNBSpillInfo object 
   */
  sbn::BNBSpillInfo makeBNBSpillInfo(art::EventID const& eventID, double time, MWRdata_t const& MWRdata, std::vector<int> const& matched_MWR, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp, std::unique_ptr<ifbeam_ns::BeamFolder> const& offsets, std::unique_ptr<ifbeam_ns::BeamFolder> const& vp873);
}

#endif
