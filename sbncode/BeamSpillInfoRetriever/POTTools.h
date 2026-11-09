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
#include "sbncode/BeamSpillInfoRetriever/MWRMatching.h"
#include "larcorealg/CoreUtils/counter.h"
#include <map>
#include <string>
#include <vector>

namespace sbn::pot{

  typedef struct PTBInfo_t {
    double currPTBTimeStamp  = 1e20;
    double prevPTBTimeStamp  = 0;
    unsigned int GateCounter = 0;
  } PTBInfo_t;

  typedef struct TriggerInfo_t {
    int gate_type = 0; ///< Source of the spill: `1`: BNB, `2`: NuMI
    double t_current_event  = 0;
    double t_previous_event = 0;
    unsigned int number_of_gates_since_previous_event = 0;
    std::int64_t WR_to_Spill_conversion = 0;
  } TriggerInfo_t;
  
  typedef struct MWRdata_t {
    std::vector< std::vector<double> > MWR_times;
    std::vector< std::vector< std::vector< int > > > unpacked_MWR;
  } MWRdata_t;

  /// Last valid value of each BPM offset (device name -> offset [mm]).
  using BPMOffsetCache_t = std::map<std::string, double>;

  /**
   * @brief Extracts information from PTB for use in SBND POT accounting.
   * 
   * @param cont_frags The PTB fragments to examine.
   * @param HLT The high level trigger we are searching for.
   */
  PTBInfo_t extractPTBInfo(art::Handle<std::vector<artdaq::Fragment> > cont_frags, int HLT);

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
   * @param matched_MWR index of the multiwire reading of each device for this
   *                    spill (see `matchMWRToSpill()`); `-1` means none
   * @param offsetCache if not null, BPM offsets that cannot be read from the
   *                    database are taken from the last valid reading (they are
   *                    slowly changing settings), and valid readings update it
   *
   * Without `offsetCache` a failed offset query leaves the offset at -999, and
   * the spill gets no FOM even when all the BPM readings are present.
   */
  sbn::BNBSpillInfo makeBNBSpillInfo(art::EventID const& eventID, double time, MWRdata_t const& MWRdata, std::vector<int> const& matched_MWR, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp, std::unique_ptr<ifbeam_ns::BeamFolder> const& offsets, std::unique_ptr<ifbeam_ns::BeamFolder> const& vp873, BPMOffsetCache_t* offsetCache = nullptr);

  /**
   * @brief Fills BPM offsets still missing (-999) in `info` from `cache`.
   * @return whether any offset was filled (the FOM then needs recomputing);
   *         if so, `fomstatus::OffsetFilled` is set in `info.fom_status`
   *
   * Meant for the end of the subrun: spills recorded before the first valid
   * offset reading of the job could not use the cache when they were made.
   */
  bool fillMissingBPMOffsets(sbn::BNBSpillInfo& info, BPMOffsetCache_t const& cache);

  /**
   * @brief Fills the Booster bunch-rotation settings of a spill.
   * @param info the spill (`BRRMPS8`, `BRRMPQ1`, `BRRMPL9` are set)
   * @param time spill time [s]
   * @param folder IFBeam folder of a bundle holding the settings
   * @param devices IFBeam names of the scale factor, frequency and start time,
   *                in this order (e.g. `B_BRRMPS[8]`, `B_BRRMPQ[1]`, `B_BRRMPL[9]`)
   * @param cache if not null, last valid value of each setting: used when a query
   *              fails, and updated otherwise
   *
   * The settings are logged every 15 minutes, so the folder needs an epsilon of
   * at least that. A setting that cannot be found is -999.
   */
  void fillBunchRotation(sbn::BNBSpillInfo& info, double time, ifbeam_ns::BeamFolder& folder,
                         std::vector<std::string> const& devices, BPMOffsetCache_t* cache);
}

#endif
