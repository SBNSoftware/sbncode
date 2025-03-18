////////////////////////////////////////////////////////////////////////
// Class:       SBNDBNBZEROBIASRetriever
// Plugin Type: producer 
// File:        SBNDBNBZEROBIASRetriever_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbncode/BeamSpillInfoRetriever/POTTools.h"

namespace sbn {
  class SBNDBNBZEROBIASRetriever;
}

class sbn::SBNDBNBZEROBIASRetriever : public art::EDProducer {
public:
  explicit SBNDBNBZEROBIASRetriever(fhicl::ParameterSet const & params);
  // Required functions.
  void produce(art::Event & e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

  // Plugins should not be copied or assigned.
  SBNDBNBZEROBIASRetriever(SBNDBNBZEROBIASRetriever const &) = delete;
  SBNDBNBZEROBIASRetriever(SBNDBNBZEROBIASRetriever &&) = delete;
  SBNDBNBZEROBIASRetriever & operator = (SBNDBNBZEROBIASRetriever const &) = delete;
  SBNDBNBZEROBIASRetriever & operator = (SBNDBNBZEROBIASRetriever &&) = delete;

private:
  // Declare member data here.
  double fTimePad;
  double fBESOffset;
  std::string fDeviceUsedForTiming;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_mwr;
  sbn::MWRData mwrdata;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  static constexpr double MWRtoroidDelay = -0.035; ///< the same time point is measured _t_ by MWR and _t + MWRtoroidDelay`_ by the toroid [ms]
  std::vector< sbn::BNBSpillInfo > fOutbeamInfos;
  std::vector< sbn::BNBSpillInfo > fOutbeamInfosTotal;

  TriggerInfo_t extractTriggerInfo(art::Event const& e) const;
  void matchMultiWireData(
    art::EventID const& eventID, 
    TriggerInfo_t const& triggerInfo,
    MWRdata_t const& MWRdata,
    std::vector< sbn::BNBSpillInfo >& beamInfos
    ) const;
  unsigned int TotalBeamSpills;
};

sbn::SBNDBNBZEROBIASRetriever::SBNDBNBZEROBIASRetriever(fhicl::ParameterSet const & params)
  : EDProducer{params} {
  fTimePad = params.get<double>("TimePadding");
  fDeviceUsedForTiming = params.get<std::string>("DeviceUsedForTiming");
  fBESOffset = params.get<double>("BESOffset");
  bfp = ifbeam_handle->getBeamFolder(params.get<std::string>("Bundle"), params.get<std::string>("URL"), std::stod(params.get<std::string>("TimeWindow")));
  bfp->set_epsilon(0.02);
  bfp_mwr = ifbeam_handle->getBeamFolder(params.get<std::string>("MultiWireBundle"), params.get<std::string>("URL"), std::stod(params.get<std::string>("MWR_TimeWindow")));
  bfp_mwr->set_epsilon(0.5);
  bfp_mwr->setValidWindow(3605);
  TotalBeamSpills = 0;
  produces< std::vector< sbn::BNBSpillInfo >, art::InEvent >();
  produces< std::vector< sbn::BNBSpillInfo >, art::InSubRun >();
}

void sbn::SBNDBNBZEROBIASRetriever::produce(art::Event & e)
{
  // If this is the first event in the run, then ignore it
  // We do not currently have the ability to figure out the first
  // spill that the DAQ was sensitive to, so don't try to save any
  // spill information
  if (e.event() == 1) {
    auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
    std::swap(*p, fOutbeamInfos);
    e.put(std::move(p));
    return;
  }

  TriggerInfo_t const triggerInfo = extractTriggerInfo(e);

  if (triggerInfo.t_previous_event == 0) {
    auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
    std::swap(*p, fOutbeamInfos);
    e.put(std::move(p));
    return;
  }

  TotalBeamSpills += triggerInfo.number_of_gates_since_previous_event;
  MWRdata_t const MWRdata = extractSpillTimes(triggerInfo, bfp, bfp_mwr, fTimePad, MWRtoroidDelay, mwrdata);

  matchMultiWireData(e.id(), triggerInfo, MWRdata, fOutbeamInfos);
  fOutbeamInfosTotal.insert(fOutbeamInfosTotal.end(),fOutbeamInfos.begin(),fOutbeamInfos.end());
  auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
  std::swap(*p, fOutbeamInfos);
  e.put(std::move(p));
}

sbn::TriggerInfo_t sbn::SBNDBNBZEROBIASRetriever::extractTriggerInfo(art::Event const& e) const {
  // Using TDC for current event, but PTB for previous event
  art::InputTag PTB_itag("daq", "ContainerPTB");
  auto PTB_cont_frags = e.getHandle<artdaq::Fragments>(PTB_itag);

  art::InputTag TDC_itag("daq", "ContainerTDCTIMESTAMP");
  auto TDC_cont_frags = e.getHandle<artdaq::Fragments>(TDC_itag);

  TriggerInfo_t triggerInfo;
  PTBInfo_t PTBInfo = extractPTBInfo(PTB_cont_frags, 1);

  if (TDC_cont_frags) {
    double TDCTimeStamp = extractTDCTimeStamp(TDC_cont_frags);
    triggerInfo.t_current_event = TDCTimeStamp - fBESOffset;
  }
  else{
    // If missing TDC, use PTB instead
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Missing TDC Container Fragments!!!" << std::endl;
    triggerInfo.t_current_event = PTBInfo.currPTBTimeStamp - fBESOffset;
  }

  triggerInfo.t_previous_event = PTBInfo.prevPTBTimeStamp - fBESOffset;
  triggerInfo.number_of_gates_since_previous_event = PTBInfo.GateCounter;

  double PTBandCurrOffset = PTBInfo.currPTBTimeStamp - triggerInfo.t_current_event - fBESOffset;

  // Catch for an issue seen a few times where PTB off by a second.
  // Only need to correct prevTS because either currTS is from TDC
  // or there is no offset between currTS and PTB.
  if(abs(PTBandCurrOffset) >= 0.5){
    triggerInfo.t_previous_event-=PTBandCurrOffset;
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << "Offset between PTB and TDC, " << PTBandCurrOffset << std::endl;
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << "Corrected previous event TS is " << std::setprecision(19) << triggerInfo.t_previous_event << std::endl;
  }

  return triggerInfo;
}

void sbn::SBNDBNBZEROBIASRetriever::matchMultiWireData(
  art::EventID const& eventID,
  TriggerInfo_t const& triggerInfo,
  MWRdata_t const& MWRdata,
  std::vector< sbn::BNBSpillInfo >& beamInfos
) const {
  
  auto const& [ MWR_times, unpacked_MWR ] = MWRdata; // alias
  
  //Here we will start collecting all the other beamline devices
  // First we get the times that the beamline device fired
  //  we have to pick a specific variable to use
  std::vector<double> times_temps = bfp->GetTimeList(fDeviceUsedForTiming);

  mf::LogDebug("SBNDBNBZEROBIASRetriever") << "matchMultiWireData:: Number of time spills : " << times_temps.size() << std::endl;

  // We'll keep track of how many of these spills match to our 
  // DAQ trigger times
  int spills_removed = 0;
  std::vector<int> matched_MWR;
  matched_MWR.resize(3);

  // Iterating through each of the beamline times
 
  double best_diff = 10000000000.0; 
  double diff; 
  size_t i = 0;

  for (size_t k = 0; k < times_temps.size(); k++){
    diff = (triggerInfo.t_current_event + fTimePad) - times_temps[k];// diff is greater than zero!
    if( diff > 0 and diff < best_diff){
      if(times_temps[k] > (triggerInfo.t_current_event)+fTimePad){
        spills_removed++; 
        continue;} 
      if(times_temps[k] <= (triggerInfo.t_previous_event)+fTimePad){
        spills_removed++; 
        continue;}

      if(BrokenClock(times_temps[i], bfp)){
        continue;
      }

      best_diff = diff; 
      i = k;
    }
  }

  for(int dev = 0; dev < int(MWR_times.size()); dev++){
    //Loop through the multiwire times:
    double Tdiff = 1000000000.;
    matched_MWR[dev] = 0;

    for(int mwrt = 0;  mwrt < int(MWR_times[dev].size()); mwrt++){
      //found a candidate match! 
      if(fabs((MWR_times[dev][mwrt] - times_temps[i])) >= Tdiff){continue;}
	
      bool best_match = true;
	  
      for (size_t j = 0; j < times_temps.size(); j++) {
        //Check for a better match...
        if( j == i) continue;
        if(times_temps[j] > (triggerInfo.t_current_event+fTimePad)){continue;}
	if(times_temps[j] <= (triggerInfo.t_previous_event+fTimePad)){continue;}
	  
	//is there a better match later in the spill sequence
	if(fabs((MWR_times[dev][mwrt] - times_temps[j])) < 
	  fabs((MWR_times[dev][mwrt] - times_temps[i]))){
	  //we can have patience...
	  best_match = false;
	  break;
	}	     
      }//end better match check
	
      //Verified best match! 
      if(best_match == true){
        matched_MWR[dev] = mwrt;
	Tdiff = fabs((MWR_times[dev][mwrt] - times_temps[i]));
      }	
    }//end loop over MWR times 
  }//end loop over MWR devices
    
  sbn::BNBSpillInfo spillInfo = makeBNBSpillInfo(eventID, times_temps[i], MWRdata, matched_MWR, bfp);
  beamInfos.push_back(std::move(spillInfo));

  // We do not write these to the art::Events because 
  // we can filter events but want to keep all the POT 
  // information, so we'll write it to the SubRun
}

void sbn::SBNDBNBZEROBIASRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

void sbn::SBNDBNBZEROBIASRetriever::endSubRun(art::SubRun& sr)
{
  mf::LogDebug("SBNDBNBZEROBIASRetriever")<< "Total number of DAQ Spills : " << TotalBeamSpills << std::endl;
  mf::LogDebug("SBNDBNBZEROBIASRetriever")<< "Total number of Selected Spills : " << fOutbeamInfosTotal.size() << std::endl;
  auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
  std::swap(*p, fOutbeamInfosTotal);
  sr.put(std::move(p), art::subRunFragment());
  return;
}

DEFINE_ART_MODULE(sbn::SBNDBNBZEROBIASRetriever)
