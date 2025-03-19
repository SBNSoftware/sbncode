#include "sbncode/WireMod/Utility/WireModUtility.hh"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

//--- CalcROIProperties ---
sys::WireModUtility::ROIProperties_t sys::WireModUtility::CalcROIProperties(recob::Wire const& wire, size_t const& roi_idx)
{
  // get the ROI
  recob::Wire::RegionsOfInterest_t::datarange_t const& roi = wire.SignalROI().get_ranges()[roi_idx];

  // initialize the return value
  ROIProperties_t roi_vals;
  roi_vals.channel = wire.Channel();
  roi_vals.view    = wire.View();
  roi_vals.begin   = roi.begin_index();
  roi_vals.end     = roi.end_index();
  roi_vals.center  = 0;
  roi_vals.total_q = 0;
  roi_vals.sigma   = 0;

  // loop over the roi and find the charge-weighted center and total charge
  auto const& roi_data = roi.data();
  for (size_t i_t = 0; i_t < roi_data.size(); ++i_t)
  {
    roi_vals.center += roi_data[i_t]*(i_t+roi_vals.begin);
    roi_vals.total_q += roi_data[i_t];
  }
  roi_vals.center = roi_vals.center/roi_vals.total_q;

  // get the width (again charge-weighted)
  // if the ROI is only one tick set the cent to the middle of the tick and width to 0.5
  for (size_t i_t = 0; i_t<roi_data.size(); ++i_t)
  {
    roi_vals.sigma += roi_data[i_t]*(i_t+roi_vals.begin-roi_vals.center)*(i_t+roi_vals.begin-roi_vals.center);
  }
  roi_vals.sigma = std::sqrt(roi_vals.sigma/roi_vals.total_q);
  if (roi_vals.end-roi_vals.begin == 1)
  {
    roi_vals.center += 0.5;
    roi_vals.sigma   = 0.5;
  }
  
  // return the calc'd properies
  return roi_vals;
}

//--- GetTargetROIs ---
std::vector<std::pair<unsigned int, unsigned int>> sys::WireModUtility::GetTargetROIs(sim::SimEnergyDeposit const& shifted_edep, double offset)
{
  // initialize return value
  // the pairs are <channel number, tick time>
  std::vector<std::pair<unsigned int, unsigned int>> target_roi_vec;

  // if the SimEnergyDeposit isn't inside the TPC then it's not going to match a wire
  geo::TPCGeo const* curTPCGeomPtr = geometry->PositionToTPCptr(shifted_edep.MidPoint());
  if (curTPCGeomPtr == nullptr)
    return target_roi_vec;

  // iterate over planes
  for (size_t i_p = 0; i_p < curTPCGeomPtr->Nplanes(); ++i_p)
  {
    // check the wire exists in this plane
    int wireNumber = int(0.5 + curTPCGeomPtr->Plane(i_p).WireCoordinate(shifted_edep.MidPoint()));
    if ((wireNumber < 0) || (wireNumber >= (int) curTPCGeomPtr->Plane(i_p).Nwires()))
      continue;
    
    // reconstruct the wireID from the position
    geo::WireID edep_wireID = curTPCGeomPtr->Plane(i_p).NearestWireID(shifted_edep.MidPoint());
 
    // if the deposition is inside the range where it would leave a signal on the wire, look for the ROI there
    if (planeXInWindow(shifted_edep.X(), i_p, *curTPCGeomPtr, offset + tickOffset))
      target_roi_vec.emplace_back(geometry->PlaneWireToChannel(edep_wireID), std::round(planeXToTick(shifted_edep.X(), i_p, *curTPCGeomPtr, offset + tickOffset)));
  }

  return target_roi_vec;
}

//--- GetHitTargetROIs ---
std::vector<std::pair<unsigned int, unsigned int>> sys::WireModUtility::GetHitTargetROIs(recob::Hit const& hit)
{
  // initialize return value
  // the pairs are <channel number, tick time>
  std::vector<std::pair<unsigned int, unsigned int>> target_roi_vec;

  int hit_wire = hit.Channel();
  int hit_tick = int(round(hit.PeakTime()));

  if (hit_tick < tickOffset || hit_tick >= readoutWindowTicks + tickOffset)
    return target_roi_vec;


  target_roi_vec.emplace_back((unsigned int) hit_wire, (unsigned int) hit_tick);

  return target_roi_vec;
}

//--- FillROIMatchedEdepMap ---
void sys::WireModUtility::FillROIMatchedEdepMap(std::vector<sim::SimEnergyDeposit> const& edepVec, std::vector<recob::Wire> const& wireVec, double offset)
{
  // clear the map in case it was already set
  ROIMatchedEdepMap.clear();

  // get the channel from each wire and set up a map between them
  std::unordered_map<unsigned int,unsigned int> wireChannelMap;
  for (size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    wireChannelMap[wireVec[i_w].Channel()] = i_w;

  // loop over the energy deposits
  for (size_t i_e = 0; i_e < edepVec.size(); ++i_e)
  {
    // get the ROIs
    // <channel number, tick time>
    std::vector<std::pair<unsigned int, unsigned int>> target_rois = GetTargetROIs(edepVec[i_e], offset);
    
    // loop over ROI and match the energy deposits with wires
    for (auto const& target_roi : target_rois)
    {
      // if we can't find the wire, skip it
      if (wireChannelMap.find(target_roi.first) == wireChannelMap.end())
        continue;

      // get the wire
      auto const& target_wire = wireVec.at(wireChannelMap[target_roi.first]);

      // if there are no ticks in in the wire skip it
      // likewise if there's nothing in the region of interst
      if (not target_wire.SignalROI().is_valid()              ||
          target_wire.SignalROI().empty()                     ||
          target_wire.SignalROI().n_ranges() == 0             ||
          target_wire.SignalROI().size() <= target_roi.second ||
          target_wire.SignalROI().is_void(target_roi.second)   )
        continue;

      // how far into the range is the ROI?
      auto range_number = target_wire.SignalROI().find_range_iterator(target_roi.second) - target_wire.SignalROI().begin_range();

      // popluate the map
      ROIMatchedEdepMap[std::make_pair(target_wire.Channel(),range_number)].push_back(i_e);
    }
  }
}

//--- FillROIMatchedHitMap ---
void sys::WireModUtility::FillROIMatchedHitMap(std::vector<recob::Hit> const& hitVec, std::vector<recob::Wire> const& wireVec)
{
  // clear the map in case it was already set
  ROIMatchedHitMap.clear();

  // get the channel from each wire and set up a map between them
  std::unordered_map<unsigned int,unsigned int> wireChannelMap;
  for (size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    wireChannelMap[wireVec[i_w].Channel()] = i_w;

  // loop over hits
  for (size_t i_h = 0; i_h < hitVec.size(); ++i_h)
  {
    // get the ROIs
    // <channel number, tick time>
    std::vector<std::pair<unsigned int, unsigned int>> target_rois = GetHitTargetROIs(hitVec[i_h]);

    // loop over ROI and match the energy deposits with wires
    for (auto const& target_roi : target_rois)
    {
      // if we can't find the wire, skip it
      if (wireChannelMap.find(target_roi.first) == wireChannelMap.end())
        continue;

      auto const& target_wire = wireVec.at(wireChannelMap[target_roi.first]);

      // if there are no ticks in in the wire skip it
      // likewise if there's nothing in the region of interst
      if (not target_wire.SignalROI().is_valid()              ||
          target_wire.SignalROI().empty()                     ||
          target_wire.SignalROI().n_ranges() == 0             ||
          target_wire.SignalROI().size() <= target_roi.second ||
          target_wire.SignalROI().is_void(target_roi.second)   )
        continue;

      // which range is it?
      auto range_number = target_wire.SignalROI().find_range_iterator(target_roi.second) - target_wire.SignalROI().begin_range();

      // pupluate the map
      ROIMatchedHitMap[std::make_pair(target_wire.Channel(),range_number)].push_back(i_h);
    }
  }
}

//--- CalcSubROIProperties ---
std::vector<sys::WireModUtility::SubROIProperties_t> sys::WireModUtility::CalcSubROIProperties(sys::WireModUtility::ROIProperties_t const& roi_properties, std::vector<const recob::Hit*> const& hitPtrVec)
{
  std::vector<sys::WireModUtility::SubROIProperties_t> subroi_properties_vec;
  sys::WireModUtility::SubROIProperties_t subroi_properties;
  subroi_properties.channel = roi_properties.channel;
  subroi_properties.view    = roi_properties.view;

  // if this ROI doesn't contain any hits, define subROI based on ROI properities
  // otherwise, define subROIs based on hits
  if (hitPtrVec.size() == 0)
  {
    subroi_properties.key     = std::make_pair(roi_properties.key, 0);
    subroi_properties.total_q = roi_properties.total_q;
    subroi_properties.center  = roi_properties.center;
    subroi_properties.sigma   = roi_properties.sigma;
    subroi_properties_vec.push_back(subroi_properties);
  } else
  {
    for (unsigned int i_h=0; i_h < hitPtrVec.size(); ++i_h)
    {
      auto hit_ptr = hitPtrVec[i_h];
      subroi_properties.key     = std::make_pair(roi_properties.key, i_h);
      subroi_properties.total_q = hit_ptr->Integral();
      subroi_properties.center  = hit_ptr->PeakTime();
      subroi_properties.sigma   = hit_ptr->RMS();
      subroi_properties_vec.push_back(subroi_properties);
    }
  }

  return subroi_properties_vec;
}

//--- MatchEdepsToSubROIs ---
std::map<sys::WireModUtility::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> sys::WireModUtility::MatchEdepsToSubROIs(std::vector<sys::WireModUtility::SubROIProperties_t> const& subROIPropVec, 
                                                                                                                  std::vector<const sim::SimEnergyDeposit*> const& edepPtrVec, double offset)
{
  // for each TrackID, which EDeps are associated with it? keys are TrackIDs
  std::map<int, std::vector<const sim::SimEnergyDeposit*>> TrackIDMatchedEDepMap;

  // total energy of EDeps matched to the ROI (not strictly necessary, but useful for understanding/development
  double total_energy = 0.0;

  // loop over edeps, fill TrackIDMatchedEDepMap and calculate total energy
  for (auto const& edep_ptr : edepPtrVec)
  {
    TrackIDMatchedEDepMap[edep_ptr->TrackID()].push_back(edep_ptr);
    total_energy += edep_ptr->E();
  }

  // calculate EDep properties by TrackID
  std::map<int, sys::WireModUtility::TruthProperties_t> TrackIDMatchedPropertyMap;
  for (auto const& track_edeps : TrackIDMatchedEDepMap)
    TrackIDMatchedPropertyMap[track_edeps.first] = CalcPropertiesFromEdeps(track_edeps.second, offset);

  // map it all out
  std::map<unsigned int, std::vector<unsigned int>> EDepMatchedSubROIMap;        // keys are indexes of edepPtrVec, values are vectors of indexes of subROIPropVec
  std::map<int, std::unordered_set<unsigned int>>   TrackIDMatchedSubROIMap;     // keys are TrackIDs, values are sets of indexes of subROIPropVec
  std::map<unsigned int, std::vector<unsigned int>> SubROIMatchedEDepMap;        // keys are indexes of subROIPropVec, values are vectors of indexes of edepPtrVec
  std::map<unsigned int, std::map<int, double>>     SubROIMatchedTrackEnergyMap; // keys are indexes of subROIPropVec, values are maps of TrackIDs to matched energy (in MeV)

  // loop over EDeps
  for (unsigned int i_e = 0; i_e < edepPtrVec.size(); ++i_e)
  {
    // get EDep properties
    auto edep_ptr  = edepPtrVec[i_e];
    const geo::TPCGeo& curTPCGeom = geometry->PositionToTPC(edep_ptr->MidPoint());
    double ticksPercm = detPropData.GetXTicksCoefficient(); // this should be by TPCID, but isn't building right now
    double zeroTick = detPropData.ConvertXToTicks(0, curTPCGeom.Plane(0).ID());
    auto edep_tick = ticksPercm * edep_ptr->X() + (zeroTick + offset) + tickOffset;
    edep_tick = detPropData.ConvertXToTicks(edep_ptr->X(), curTPCGeom.Plane(0).ID()) + offset + tickOffset;

    // loop over subROIs
    unsigned int closest_hit = std::numeric_limits<unsigned int>::max();
    float min_dist = std::numeric_limits<float>::max();
    for (unsigned int i_h = 0; i_h < subROIPropVec.size(); ++i_h)
    {
      auto subroi_prop = subROIPropVec[i_h];
      if (edep_tick > subroi_prop.center-subroi_prop.sigma && edep_tick < subroi_prop.center+subroi_prop.sigma)
      {
        EDepMatchedSubROIMap[i_e].push_back(i_h);
        TrackIDMatchedSubROIMap[edep_ptr->TrackID()].emplace(i_h);
      }
      float hit_dist = std::abs(edep_tick - subroi_prop.center) / subroi_prop.sigma;
      if (hit_dist < min_dist)
      {
        closest_hit = i_h;
        min_dist = hit_dist;
      }
    }

    // if EDep is less than 2.5 units away from closest subROI, assign it to that subROI
    if (min_dist < 5) // try 5 for testing purposes
    {
      auto i_h = closest_hit;
      SubROIMatchedEDepMap[i_h].push_back(i_e);
      SubROIMatchedTrackEnergyMap[i_h][edep_ptr->TrackID()] += edep_ptr->E();
    }
  }

  // convert to desired format (possibly a better way to do this...?)
  std::map<SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> ReturnMap;
  for (auto it_h = SubROIMatchedEDepMap.begin(); it_h != SubROIMatchedEDepMap.end(); ++it_h)
  {
    auto key = subROIPropVec[it_h->first].key;
    for (auto const& i_e : it_h->second)
    {
      ReturnMap[key].push_back(edepPtrVec[i_e]);
    }
  }

  return ReturnMap;
}

//--- CalcPropertiesFromEdeps ---
sys::WireModUtility::TruthProperties_t sys::WireModUtility::CalcPropertiesFromEdeps(std::vector<const sim::SimEnergyDeposit*> const& edepPtrVec, double offset)
{
  //split the edeps by TrackID
  std::map< int, std::vector<const sim::SimEnergyDeposit*> > edepptrs_by_trkid;
  std::map< int, double > energy_per_trkid;
  for(auto const& edep_ptr : edepPtrVec)
  {
    edepptrs_by_trkid[edep_ptr->TrackID()].push_back(edep_ptr);
    energy_per_trkid[edep_ptr->TrackID()]+=edep_ptr->E();
  }

  int trkid_max     = std::numeric_limits<int>::min();
  double energy_max = std::numeric_limits<double>::min();
  for(auto const& e_p_id : energy_per_trkid)
  {
    if(e_p_id.second > energy_max)
    {
      trkid_max = e_p_id.first;
      energy_max = e_p_id.second;
    }
  }

  auto edepPtrVecMaxE = edepptrs_by_trkid[trkid_max];
  
  //first, let's loop over all edeps and get an average weight scale...
  sys::WireModUtility::TruthProperties_t edep_props;
  double total_energy_all = 0.0;
  sys::WireModUtility::ScaleValues_t scales_e_weighted[3];
  for(size_t i_p = 0; i_p < 3; ++i_p)
  {
    scales_e_weighted[i_p].r_Q     = 0.0;
    scales_e_weighted[i_p].r_sigma = 0.0;
  }

  for(auto const edep_ptr : edepPtrVec)
  {
    if (edep_ptr->StepLength() == 0)
      continue;

    edep_props.x = edep_ptr->X();
    edep_props.y = edep_ptr->Y();
    edep_props.z = edep_ptr->Z();

    edep_props.dxdr = (edep_ptr->EndX() - edep_ptr->StartX()) / edep_ptr->StepLength();
    edep_props.dydr = (edep_ptr->EndY() - edep_ptr->StartY()) / edep_ptr->StepLength();
    edep_props.dzdr = (edep_ptr->EndZ() - edep_ptr->StartZ()) / edep_ptr->StepLength();
    edep_props.dedr = edep_ptr->E() / edep_ptr->StepLength();

    total_energy_all += edep_ptr->E();

    const geo::TPCGeo& curTPCGeom = geometry->PositionToTPC(edep_ptr->MidPoint());

    for (size_t i_p = 0; i_p < 3; ++i_p)
    {
      auto scales = GetViewScaleValues(edep_props, curTPCGeom.Plane(i_p).View());
      scales_e_weighted[i_p].r_Q     += edep_ptr->E()*scales.r_Q;
      scales_e_weighted[i_p].r_sigma += edep_ptr->E()*scales.r_sigma; 
    }
  }

  for(size_t i_p = 0; i_p < 3; ++i_p)
  {
    if (total_energy_all > 0)
    {
      scales_e_weighted[i_p].r_Q     = scales_e_weighted[i_p].r_Q / total_energy_all;
      scales_e_weighted[i_p].r_sigma = scales_e_weighted[i_p].r_sigma / total_energy_all;
    }
    if (scales_e_weighted[i_p].r_Q == 0)
      scales_e_weighted[i_p].r_Q = 1;
    if (scales_e_weighted[i_p].r_sigma == 0)
      scales_e_weighted[i_p].r_sigma = 1;
  }

  TruthProperties_t edep_col_properties;

  //copy in the scales that were calculated before
  for(size_t i_p = 0; i_p < 3; ++i_p)
  {
    edep_col_properties.scales_avg[i_p].r_Q = scales_e_weighted[i_p].r_Q;
    edep_col_properties.scales_avg[i_p].r_sigma = scales_e_weighted[i_p].r_sigma;
  }

  // calculations happen here
  edep_col_properties.x              = 0.;
  edep_col_properties.x_rms          = 0.;
  edep_col_properties.x_rms_noWeight = 0.;
  edep_col_properties.x_min          = std::numeric_limits<float>::max();
  edep_col_properties.x_max          = std::numeric_limits<float>::min();

  edep_col_properties.y    = 0.;
  edep_col_properties.z    = 0.;
  edep_col_properties.dxdr = 0.;
  edep_col_properties.dydr = 0.;
  edep_col_properties.dzdr = 0.;
  edep_col_properties.dedr = 0.;

  double total_energy = 0.0;
  for (auto const& edep_ptr : edepPtrVecMaxE)
  {
    edep_col_properties.x += edep_ptr->X()*edep_ptr->E();
    edep_col_properties.x_min = (edep_ptr->X() < edep_col_properties.x_min) ? edep_ptr->X() : edep_col_properties.x_min;
    edep_col_properties.x_max = (edep_ptr->X() > edep_col_properties.x_max) ? edep_ptr->X() : edep_col_properties.x_max;
    total_energy += edep_ptr->E();
    edep_col_properties.y += edep_ptr->Y()*edep_ptr->E();
    edep_col_properties.z += edep_ptr->Z()*edep_ptr->E();

    if (edep_ptr->StepLength() == 0)
      continue;

    edep_col_properties.dxdr += edep_ptr->E()*(edep_ptr->EndX() - edep_ptr->StartX()) / edep_ptr->StepLength();
    edep_col_properties.dydr += edep_ptr->E()*(edep_ptr->EndY() - edep_ptr->StartY()) / edep_ptr->StepLength();
    edep_col_properties.dzdr += edep_ptr->E()*(edep_ptr->EndZ() - edep_ptr->StartZ()) / edep_ptr->StepLength();
    edep_col_properties.dedr += edep_ptr->E()*edep_ptr->E() / edep_ptr->StepLength();
  }

  if (total_energy > 0)
  {
    edep_col_properties.x    = edep_col_properties.x / total_energy;
    edep_col_properties.y    = edep_col_properties.y / total_energy;
    edep_col_properties.z    = edep_col_properties.z / total_energy;
    edep_col_properties.dxdr = edep_col_properties.dxdr / total_energy;
    edep_col_properties.dydr = edep_col_properties.dydr / total_energy;
    edep_col_properties.dzdr = edep_col_properties.dzdr / total_energy;
    edep_col_properties.dedr = edep_col_properties.dedr / total_energy;
  }

  for (auto const& edep_ptr : edepPtrVecMaxE)
  {
    edep_col_properties.x_rms          += (edep_ptr->X()-edep_col_properties.x)*(edep_ptr->X()-edep_col_properties.x)*edep_ptr->E();
    edep_col_properties.x_rms_noWeight += (edep_ptr->X()-edep_col_properties.x)*(edep_ptr->X()-edep_col_properties.x);
  }
  edep_col_properties.x_rms_noWeight = std::sqrt(edep_col_properties.x_rms_noWeight);

  if (total_energy > 0)
    edep_col_properties.x_rms = std::sqrt(edep_col_properties.x_rms/total_energy);

  const geo::TPCGeo& tpcGeom = geometry->PositionToTPC({edep_col_properties.x, edep_col_properties.y, edep_col_properties.z});
  double ticksPercm = detPropData.GetXTicksCoefficient(); // this should be by TPCID, but isn't building right now
  edep_col_properties.tick              = detPropData.ConvertXToTicks(edep_col_properties.x    , tpcGeom.Plane(0).ID()) + offset + tickOffset;
  edep_col_properties.tick_rms          = ticksPercm*edep_col_properties.x_rms;
  edep_col_properties.tick_rms_noWeight = ticksPercm*edep_col_properties.x_rms_noWeight;
  edep_col_properties.tick_min          = detPropData.ConvertXToTicks(edep_col_properties.x_min, tpcGeom.Plane(0).ID()) + offset + tickOffset;
  edep_col_properties.tick_max          = detPropData.ConvertXToTicks(edep_col_properties.x_max, tpcGeom.Plane(0).ID()) + offset + tickOffset;
  edep_col_properties.total_energy      = total_energy;
  

  return edep_col_properties; 
}

//--- GetScaleValues ---
sys::WireModUtility::ScaleValues_t sys::WireModUtility::GetScaleValues(sys::WireModUtility::TruthProperties_t const& truth_props, sys::WireModUtility::ROIProperties_t const& roi_vals)
{
  sys::WireModUtility::ScaleValues_t scales;
  sys::WireModUtility::ScaleValues_t channelScales = GetChannelScaleValues(truth_props, roi_vals.channel);
  sys::WireModUtility::ScaleValues_t viewScales    = GetViewScaleValues(truth_props, roi_vals.view);
  scales.r_Q     = channelScales.r_Q     * viewScales.r_Q;
  scales.r_sigma = channelScales.r_sigma * viewScales.r_sigma;
  return scales;
}

//--- GetChannelScaleValues ---
sys::WireModUtility::ScaleValues_t sys::WireModUtility::GetChannelScaleValues(sys::WireModUtility::TruthProperties_t const& truth_props, raw::ChannelID_t const& channel)
{
  // initialize return
  sys::WireModUtility::ScaleValues_t scales;
  scales.r_Q     = 1.0;
  scales.r_sigma = 1.0;
  
  // try to get geo
  // if not in a TPC return default values
  double const truth_coords[3] = {truth_props.x, truth_props.y, truth_props.z};
  geo::TPCGeo const* curTPCGeomPtr = geometry->PositionToTPCptr(geo::vect::makePointFromCoords(truth_coords));
  if (curTPCGeomPtr == nullptr)
    return scales;

  if (applyChannelScale)
  {
    if (spline_Charge_Channel == nullptr ||
        spline_Sigma_Channel  == nullptr  )
      throw cet::exception("WireModUtility")
        << "Tried to apply channel scale factor, but could not find splines. Check that you have set those in the utility.";
    scales.r_Q     *= spline_Charge_Channel->Eval(channel);
    scales.r_sigma *= spline_Sigma_Channel ->Eval(channel);
  }

  return scales;
}

//--- GetViewScaleValues ---
sys::WireModUtility::ScaleValues_t sys::WireModUtility::GetViewScaleValues(sys::WireModUtility::TruthProperties_t const& truth_props, geo::View_t const& view)
{
  // initialize return
  sys::WireModUtility::ScaleValues_t scales;
  scales.r_Q     = 1.0;
  scales.r_sigma = 1.0;
  
  double temp_scale=1.0;
  
  // try to get geo
  // if not in a TPC return default values
  double const truth_coords[3] = {truth_props.x, truth_props.y, truth_props.z};
  geo::TPCGeo const* curTPCGeomPtr = geometry->PositionToTPCptr(geo::vect::makePointFromCoords(truth_coords));
  if (curTPCGeomPtr == nullptr)
    return scales;

  // get the plane number by the view
  size_t plane = curTPCGeomPtr->Plane(view).ID().Plane;

  if (applyXScale)
  {
    if (splines_Charge_X[plane] == nullptr || 
        splines_Sigma_X [plane] == nullptr  )
      throw cet::exception("WireModUtility")
        << "Tried to apply X scale factor, but could not find splines. Check that you have set those in the utility.";
    scales.r_Q     *= splines_Charge_X[plane]->Eval(truth_props.x);
    scales.r_sigma *= splines_Sigma_X [plane]->Eval(truth_props.x);
  }
  
  if (applyYZScale)
  {
    if (graph2Ds_Charge_YZ[plane] == nullptr || 
        graph2Ds_Sigma_YZ [plane] == nullptr  )
      throw cet::exception("WireModUtility")
        << "Tried to apply YZ scale factor, but could not find graphs. Check that you have set those in the utility.";
    temp_scale = graph2Ds_Charge_YZ[plane]->Interpolate(truth_props.z, truth_props.y);
    if(temp_scale>0.001) scales.r_Q *= temp_scale;

    temp_scale  = graph2Ds_Sigma_YZ [plane]->Interpolate(truth_props.z, truth_props.y);
    if(temp_scale>0.001) scales.r_sigma *= temp_scale;
  }

  if (applyXZAngleScale)
  {
    if (splines_Charge_XZAngle[plane] == nullptr || 
        splines_Sigma_XZAngle [plane] == nullptr  )
      throw cet::exception("WireModUtility")
        << "Tried to apply XZ-angle scale factor, but could not find splines. Check that you have set those in the utility.";
    scales.r_Q     *= splines_Charge_XZAngle[plane]->Eval(ThetaXZ_PlaneRel(truth_props.dxdr, truth_props.dydr, truth_props.dzdr, curTPCGeomPtr->Plane(plane).ThetaZ()));
    scales.r_sigma *= splines_Sigma_XZAngle [plane]->Eval(ThetaXZ_PlaneRel(truth_props.dxdr, truth_props.dydr, truth_props.dzdr, curTPCGeomPtr->Plane(plane).ThetaZ()));
  }
  if (applyYZAngleScale)
  {
    if (splines_Charge_YZAngle[plane] == nullptr || 
        splines_Sigma_YZAngle [plane] == nullptr  )
      throw cet::exception("WireModUtility")
        << "Tried to apply YZ-angle scale factor, but could not find splines. Check that you have set those in the utility.";
    scales.r_Q     *= splines_Charge_YZAngle[plane]->Eval(ThetaYZ_PlaneRel(truth_props.dxdr, truth_props.dydr, truth_props.dzdr, curTPCGeomPtr->Plane(plane).ThetaZ()));
    scales.r_sigma *= splines_Sigma_YZAngle [plane]->Eval(ThetaYZ_PlaneRel(truth_props.dxdr, truth_props.dydr, truth_props.dzdr, curTPCGeomPtr->Plane(plane).ThetaZ()));
  }

  if(applydEdXScale)
  {
    if (splines_Charge_dEdX[plane] == nullptr || 
        splines_Sigma_dEdX [plane] == nullptr  )
      throw cet::exception("WireModUtility")
        << "Tried to apply dEdX scale factor, but could not find splines. Check that you have set those in the utility.";
    scales.r_Q     *= splines_Charge_dEdX[plane]->Eval(truth_props.dedr);
    scales.r_sigma *= splines_Sigma_dEdX [plane]->Eval(truth_props.dedr);
  }

  return scales;
}

//--- ModifyROI ---
void sys::WireModUtility::ModifyROI(std::vector<float> & roi_data,
                                    sys::WireModUtility::ROIProperties_t const& roi_prop,
                                    std::vector<sys::WireModUtility::SubROIProperties_t> const& subROIPropVec,
                                    std::map<sys::WireModUtility::SubROI_Key_t, sys::WireModUtility::ScaleValues_t> const& subROIScaleMap)
{
  // do you want a bunch of messages?
  bool verbose = false;
  
  // initialize some values
  double q_orig = 0.0;
  double q_mod  = 0.0;
  double scale_ratio = 1.0;

  // loop over the ticks
  for(size_t i_t = 0; i_t < roi_data.size(); ++i_t)
  {
    // reset your values
    q_orig = 0.0;
    q_mod  = 0.0;
    scale_ratio = 1.0;

    // loop over the subs
    for (auto const& subroi_prop : subROIPropVec)
    {
      // get your scale vals
      auto scale_vals = subROIScaleMap.find(subroi_prop.key)->second;

      q_orig += gausFunc(i_t + roi_prop.begin, subroi_prop.center,                      subroi_prop.sigma,                  subroi_prop.total_q);
      q_mod  += gausFunc(i_t + roi_prop.begin, subroi_prop.center, scale_vals.r_sigma * subroi_prop.sigma, scale_vals.r_Q * subroi_prop.total_q);

      if (verbose)
        std::cout << "    Incrementing q_orig by gausFunc(" << i_t + roi_prop.begin << ", " << subroi_prop.center << ", " <<                      subroi_prop.sigma << ", " <<                  subroi_prop.total_q << ")" << '\n'
                  << "    Incrementing q_mod  by gausFunc(" << i_t + roi_prop.begin << ", " << subroi_prop.center << ", " << scale_vals.r_sigma * subroi_prop.sigma << ", " << scale_vals.r_Q * subroi_prop.total_q << ")" << std::endl;
    }

    // do some sanity checks
    if        (isnan(q_orig))
    {
      if (verbose)
        std::cout << "WARNING: obtained q_orig = NaN... setting scale to 1" << std::endl;
      scale_ratio = 1.0;
    } else if (isnan(q_mod)) {
      if (verbose)
        std::cout << "WARNING: obtained q_mod = NaN... setting scale to 0" << std::endl;
      scale_ratio = 0.0;
    } else {
      scale_ratio = q_mod / q_orig;
    }
    if(isnan(scale_ratio) || isinf(scale_ratio))
    {
      if (verbose)
        std::cout << "WARNING: obtained scale_ratio = " << q_mod << " / " << q_orig << " = NAN/Inf... setting to 1" << std::endl;
      scale_ratio = 1.0;
    }

    roi_data[i_t] = scale_ratio * roi_data[i_t];
    if (verbose)
      std::cout << "\t tick " << i_t << ":"
                <<  " data="   << roi_data[i_t]
                << ", q_orig=" << q_orig
                << ", q_mod="  << q_mod
                << ", ratio="  << scale_ratio << std::endl;
  }

  // we're done now
  return;
}
