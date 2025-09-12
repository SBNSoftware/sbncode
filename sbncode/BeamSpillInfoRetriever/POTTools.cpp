#include "POTTools.h"

using namespace std;

namespace sbn::pot{

  sbn::pot::PTBInfo_t extractPTBInfo(art::Handle<std::vector<artdaq::Fragment> > cont_frags, int HLT) {
    bool foundHLT = false;
    PTBInfo_t PTBInfo;
    for (auto const& cont : *cont_frags)
    {
      artdaq::ContainerFragment cont_frag(cont);
      for (size_t fragi = 0; fragi < cont_frag.block_count(); ++fragi)
      {
        artdaq::Fragment frag = *cont_frag[fragi];
        sbndaq::CTBFragment ctb_frag(frag);   // somehow the name CTBFragment stuck
        for(size_t word_i = 0; word_i < ctb_frag.NWords(); ++word_i)
        {
          if(ctb_frag.Trigger(word_i)){
            if (ctb_frag.Trigger(word_i)->IsHLT() && ctb_frag.Trigger(word_i)->IsTrigger(HLT))
            {
              foundHLT = true;
              uint64_t RawprevPTBTimeStamp = ctb_frag.PTBWord(word_i)->prevTS * 20;
              uint64_t RawcurrPTBTimeStamp = ctb_frag.Trigger(word_i)->timestamp * 20;
              double currTS_candidate = std::bitset<64>(RawcurrPTBTimeStamp/20).to_ullong()/50e6;
              if(currTS_candidate < PTBInfo.currPTBTimeStamp){
                PTBInfo.prevPTBTimeStamp = std::bitset<64>(RawprevPTBTimeStamp / 20).to_ullong()/50e6;
                PTBInfo.currPTBTimeStamp = currTS_candidate;
                PTBInfo.GateCounter = ctb_frag.Trigger(word_i)->gate_counter;
              }
            }
          }
        } //End of loop over the number of trigger words
      } //End of loop over the number of fragments per container
    } //End of loop over the number of containers
  
    if(foundHLT == true){
      return PTBInfo;
    }
    else{
      std::cout << "Failed to find HLT " << HLT << std::endl;
      throw std::exception();
    }
  }

  double extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags) {
  
    double TDCTimeStamp = 0;
    for (auto const& cont : *cont_frags)
    {
      artdaq::ContainerFragment cont_frag(cont);
      for (size_t fragi = 0; fragi < cont_frag.block_count(); ++fragi)
      {
        artdaq::Fragment frag = *cont_frag[fragi];
        sbndaq::TDCTimestampFragment tdc_frag(frag);
        TDCTimeStamp = static_cast<double>(tdc_frag.getTDCTimestamp()->timestamp_ns())/1e9;
      } //End of loop over the number of fragments per container
    } //End of loop over the number of containers
    return TDCTimeStamp;
  }

  sbn::BNBSpillInfo makeBNBSpillInfo
  (art::EventID const& eventID, double time, MWRdata_t const& MWRdata, std::vector<int> const& matched_MWR, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp,  const std::unique_ptr<ifbeam_ns::BeamFolder> & offsets,const  std::unique_ptr<ifbeam_ns::BeamFolder> & vp873)
  {
    
    auto const& [ MWR_times, unpacked_MWR ] = MWRdata; // alias
   
    // initializing all of our device carriers
    // device definitions can be found in BNBSpillInfo.h
    
    double TOR860 = 0; // units e12 protons
    double TOR875 = 0; // units e12 protons
    double LM875A = 0; // units R/s
    double LM875B = 0; // units R/s
    double LM875C = 0; // units R/s
    double HP873 = 0; // units mm
    double VP873 = 0; // units mm; not in the first IFBeam query bunch
    double HP875 = 0; // units mm
    double VP875 = 0; // units mm
    double HPTG1 = 0; // units mm
    double VPTG1 = 0; // units mm
    double HPTG2 = 0; // units mm
    double VPTG2 = 0; // units mm
    double BTJT2 = 0; // units Deg C
    double THCURR = 0; // units kiloAmps
    double M875HS = 0; // units mm
    double M875VS = 0; // units mm
    double M875HM = 0; // units mm
    double M875VM = 0; // units mm
    double M876HS = 0; // units mm
    double M876VS = 0; // units mm
    double M876HM = 0; // units mm
    double M876VM = 0; // units mm
    

    double HP875Offset =0;//units mm; make a another separate IFBeam query bunch for offsets
    double VP875Offset =0;//units mm
    double VP873Offset =0;//units mm
    double HPTG1Offset =0;//units mm
    double HPTG2Offset =0;//units mm
    double VPTG1Offset =0;//units mm
    double VPTG2Offset =0;//units mm

    double TOR860_time = 0; // units s

    double FOM =0; 

    // Here we request all the devices
    // since sometimes devices fail to report we'll
    // allow each to throw an exception but still move forward
    // interpreting these failures will be part of the beam quality analyses

    try{bfp->GetNamedData(time, "E:TOR860@",&TOR860,&TOR860_time);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:TOR875",&TOR875);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:LM875A",&LM875A);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:LM875B",&LM875B);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:LM875C",&LM875C);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:HP873",&HP873);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:VP873",&VP873);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:HP875",&HP875);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:VP875",&VP875);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:HPTG1",&HPTG1);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:VPTG1",&VPTG1);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:HPTG2",&HPTG2);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:VPTG2",&VPTG2);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:BTJT2",&BTJT2);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:THCURR",&THCURR);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}

    try{bfp->GetNamedData(time, "E:M875HS",&M875HS);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:M875VS",&M875VS);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:M875HM",&M875HM);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:M875VM",&M875VM);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:M876HS",&M876HS);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:M876VS",&M876VS);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:M876HM",&M876HM);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{bfp->GetNamedData(time, "E:M876VM",&M876VM);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";} 

    try{offsets->GetNamedData(time, "E_VP873S",&VP873Offset);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{offsets->GetNamedData(time, "E_HP875S",&HP875Offset);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{offsets->GetNamedData(time, "E_VP875S",&VP875Offset);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{offsets->GetNamedData(time, "E_HPTG1S",&HPTG1Offset);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{offsets->GetNamedData(time, "E_HPTG2S",&HPTG2Offset);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{offsets->GetNamedData(time, "E_VPTG1S",&VPTG1Offset);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
    try{offsets->GetNamedData(time, "E_VPTG2S",&VPTG2Offset);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}

    //crunch the times 
    unsigned long int time_closest_int = (int) TOR860_time;
    double time_closest_ns = (TOR860_time - time_closest_int)*1e9;
  
    //Store everything in our data-product
    sbn::BNBSpillInfo beamInfo;
    beamInfo.TOR860 = TOR860*1e12; //add in factor of 1e12 protons to get correct POT units
    beamInfo.TOR875 = TOR875*1e12; //add in factor of 1e12 protons to get correct POT units
    beamInfo.LM875A = LM875A;
    beamInfo.LM875B = LM875B;
    beamInfo.LM875C = LM875C;
    beamInfo.HP873 = HP873;
    beamInfo.VP873 = VP873;
    beamInfo.HP875 = HP875;
    beamInfo.VP875 = VP875;
    beamInfo.HPTG1 = HPTG1;
    beamInfo.VPTG1 = VPTG1;
    beamInfo.HPTG2 = HPTG2;
    beamInfo.VPTG2 = VPTG2;
    beamInfo.BTJT2 = BTJT2;
    beamInfo.THCURR = THCURR;
    beamInfo.M875HS = M875HS;
    beamInfo.M875VS = M875VS;
    beamInfo.M875HM = M875HM;
    beamInfo.M875VM = M875VM;
    beamInfo.M876HS = M876HS;
    beamInfo.M876VS = M876VS;
    beamInfo.M876HM = M876HM;
    beamInfo.M876VM = M876VM;
    beamInfo.spill_time_s = time_closest_int;
    beamInfo.spill_time_ns = time_closest_ns;    
    beamInfo.VP873 = VP873;
    beamInfo.VP873Offset = VP873Offset;
    beamInfo.HP875Offset = HP875Offset;
    beamInfo.VP875Offset = VP875Offset;
    beamInfo.HPTG1Offset = HPTG1Offset;
    beamInfo.HPTG2Offset = HPTG2Offset;
    beamInfo.VPTG1Offset = VPTG1Offset;
    beamInfo.VPTG2Offset = VPTG2Offset;
    beamInfo.FOM = FOM;

    
    for(auto const& MWRdata: unpacked_MWR){
      std::ignore = MWRdata;
      assert(!MWRdata.empty());
    }
  
    if(unpacked_MWR[0].empty()){
      beamInfo.M875BB.clear();
      beamInfo.M875BB_spill_time_diff = -999;//units in seconds
    }
    else{
      beamInfo.M875BB = unpacked_MWR[0][matched_MWR[0]];
      beamInfo.M875BB_spill_time_diff = (MWR_times[0][matched_MWR[0]] - time);
    }
  
   if(unpacked_MWR[1].empty()){
      beamInfo.M876BB.clear();
      beamInfo.M876BB_spill_time_diff = -999;//units in seconds
   }
   else{
     beamInfo.M876BB = unpacked_MWR[1][matched_MWR[1]];
     beamInfo.M876BB_spill_time_diff = (MWR_times[1][matched_MWR[1]] - time);
   }
  
   if(unpacked_MWR[2].empty()){
      beamInfo.MMBTBB.clear();
      beamInfo.MMBTBB_spill_time_diff = -999;//units in seconds
    }
   else{
     beamInfo.MMBTBB = unpacked_MWR[2][matched_MWR[2]];
     beamInfo.MMBTBB_spill_time_diff = (MWR_times[2][matched_MWR[2]] - time);
   }
    // We do not write these to the art::Events because 
    // we can filter events but want to keep all the POT 
    // information, so we'll write it to the SubRun
    
    beamInfo.event = eventID.event(); // the rest of ID is known by art::SubRun
    
    //for(size_t m=0; m< beamInfo.M875BB.size(); m++){
    //std::cout << "M875BB:  " <<beamInfo.M875BB[m] << " | M876BB: " << beamInfo.M876BB[m] << " | MMTBB: " << beamInfo.MMBTBB[m] <<std::endl;
    //}
    
    return beamInfo;
  }

  bool BrokenClock(double time, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp)
  {
    double TOR860 = 0; // units e12 protons
    double TOR875 = 0; // units e12 protons
    double LM875A = 0; // units R/s
    double LM875B = 0; // units R/s
    double LM875C = 0; // units R/s
    double HP875 = 0; // units mm
    double VP875 = 0; // units mm
    double HPTG1 = 0; // units mm
    double VPTG1 = 0; // units mm
    double HPTG2 = 0; // units mm
    double VPTG2 = 0; // units mm
    double BTJT2 = 0; // units Deg C
    double THCURR = 0; // units kiloAmps
    
    double TOR860_time = 0; // units s
    
    // Here we request all the devices
    // since sometimes devices fail to report we'll
    // allow each to throw an exception but still move forward
    // interpreting these failures will be part of the beam quality analyses
  
    int ExceptionCounter = 0;
  
    try{bfp->GetNamedData(time, "E:TOR860@",&TOR860,&TOR860_time);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:TOR875",&TOR875);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:LM875A",&LM875A);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:LM875B",&LM875B);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:LM875C",&LM875C);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:HP875",&HP875);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:VP875",&VP875);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:HPTG1",&HPTG1);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:VPTG1",&VPTG1);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:HPTG2",&HPTG2);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:VPTG2",&VPTG2);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:BTJT2",&BTJT2);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
    try{bfp->GetNamedData(time, "E:THCURR",&THCURR);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "got exception: " << we.what() << "\n"; ExceptionCounter++;}
  
    if(ExceptionCounter > 10){
      mf::LogDebug("SBNDBNBRetriever")<< "At time : " << std::setprecision(19) << time << " " << "found large number of device exceptions. Throwing away spill!\n";
      return true;
    }
    else{
      return false;
    }
  }

  sbn::pot::MWRdata_t extractSpillTimes(TriggerInfo_t const& triggerInfo, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp, std::unique_ptr<ifbeam_ns::BeamFolder> const& bfp_mwr, double fTimePad, double MWRtoroidDelay, sbn::MWRData mwrdata) {
    
    // These lines get everything primed within the IFBeamDB.
    try{bfp->FillCache((triggerInfo.t_current_event)+fTimePad);} catch (WebAPIException &we) {};     
    try{bfp->FillCache((triggerInfo.t_previous_event)-fTimePad);} catch (WebAPIException &we) {};      
    try{bfp_mwr->FillCache((triggerInfo.t_current_event)+fTimePad);} catch (WebAPIException &we) {};
    try{bfp_mwr->FillCache((triggerInfo.t_previous_event)-fTimePad);} catch (WebAPIException &we) {};
  
    // The multiwire chambers provide their
    // data in a vector format but we'll have 
    // to sort through it in std::string format
    // to correctly unpack it
    std::vector< std::vector< std::vector< int > > >  unpacked_MWR;
    std::vector< std::vector< double> > MWR_times;
    unpacked_MWR.resize(3);
    MWR_times.resize(3);
    std::string packed_data_str; 
    
    // Create a list of all the MWR devices with their different
    // memory buffer increments 
    // generally in the format: "E:<Device>.{Memory Block}"
    std::vector<std::string> vars = bfp_mwr->GetDeviceList();
    mf::LogDebug("SBNDBNBRetriever") << " Number of MWR Device Blocks Found : " << vars.size() << std::endl;
    // Tracking the time from the IFBeamDB
    double time_for_mwr;    
    
    // this is an iterator to track which of the 
    // three devices we will be working with
    int dev = 0;
    
    // The MWR devices are annoying and have confusing buffer
    // what we'll do is sort through all of them first and then 
    // match them to the closest spills in time
    // 
  
    int t_steps = int(((triggerInfo.t_current_event + fTimePad) - (triggerInfo.t_previous_event - fTimePad - 20.))/0.5)+25;
  
    for(int t = 0; t < t_steps; t++){//Iterate through time increments
      for (std::string const& var : vars) {// Iterate through the devices
        
        //Make sure we have a device
        if(var.empty()){ 
  	// mf::LogDebug("SBNDBNBRetriever") << " NO MWR DEVICES?!" << std::endl;
  	continue;
        }
        /// Check the device name and interate the double-vector index
        if(var.find("M875BB") != std::string::npos ) dev = 0;
        else if(var.find("M876BB") != std::string::npos ) dev = 1;
        else if(var.find("MMBTBB") != std::string::npos ) dev = 2;
        else{
  	mf::LogDebug("SBNDBNBRetriever") << " NOT matched to a MWR DEVICES?!" << var << std::endl;
  	continue;}
        
        time_for_mwr = 0;
        
        try{
  
  	//Pull the MWR data for the device
  	// these data are "packed"
  	std::vector<double> packed_MWR = bfp_mwr->GetNamedVector((triggerInfo.t_previous_event)-20.-fTimePad+double(0.5*t),var,&time_for_mwr);
  
  	//We'll convert this into a format
  	// that we can unpack doubles >> strings
  	//
  	packed_data_str.clear();
  	packed_data_str += std::to_string(int(time_for_mwr));
  	packed_data_str.append(",");
  	packed_data_str.append(var);
  	packed_data_str.append(",,");
  	
  	for(int j = 0; j < int(packed_MWR.size()); j++){
  	  packed_data_str += std::to_string(int(packed_MWR[j]));
  	  if(j < int(packed_MWR.size())-1)
  	    packed_data_str.append(",");
  	}
  	
  	// Use Zarko's unpacking function to turn this into consumeable data
  	std::vector<double> MWR_times_temp;
  	
  	// There is a 35 ms offset between the toriod and the MWR times
  	//   we'll just remove that here to match to the spill times
  	std::vector< std::vector< int > > unpacked_MWR_temp = mwrdata.unpackMWR(packed_data_str,MWR_times_temp,MWRtoroidDelay);
  	
  	//There are four events that are packed into one MWR IFBeam entry
  	for(std::size_t s: util::counter(unpacked_MWR_temp.size())){
  	  	  
  	  // If this entry has a unique time them store it for later	  
  	  if(std::find(MWR_times[dev].begin(), MWR_times[dev].end(), MWR_times_temp[s]) == MWR_times[dev].end()){
  	    unpacked_MWR[dev].push_back(unpacked_MWR_temp[s]);
  	    MWR_times[dev].push_back(MWR_times_temp[s]);
  	  }//check for unique time 
  	}//Iterate through the unpacked events
  	}//try
        catch (WebAPIException &we) {
  	//Ignore when we can't find the MWR devices
  	//   they don't always report and the timing of them can be annoying
  	}//catch
      }// Iterate over all the multiwire devices
    }// Iterate over all times
  
    mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[0] times : " << MWR_times[0].size() << std::endl;	
    mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[0]s : " << unpacked_MWR[0].size() << std::endl;	
    mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[1] times : " << MWR_times[1].size() << std::endl;	
    mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[1]s : " << unpacked_MWR[1].size() << std::endl;	
    mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[2] times : " << MWR_times[2].size() << std::endl;	
    mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[2]s : " << unpacked_MWR[2].size() << std::endl;	
    
    return { std::move(MWR_times), std::move(unpacked_MWR) };
  }

}



