////////////////////////////////////////////////////////////////////////
// Class:       TrackAreaHit
// Plugin Type: producer (art v3_02_06)
// File:        TrackAreaHit_module.cc
//
// Generated at Wed Feb 19 17:38:21 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/Utilities/sparse_vector.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcorealg/Geometry/Exceptions.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include <memory>

namespace sbn {
  class TrackAreaHit;
}


class sbn::TrackAreaHit : public art::EDProducer {
public:
  explicit TrackAreaHit(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackAreaHit(TrackAreaHit const&) = delete;
  TrackAreaHit(TrackAreaHit&&) = delete;
  TrackAreaHit& operator=(TrackAreaHit const&) = delete;
  TrackAreaHit& operator=(TrackAreaHit&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  recob::Hit MakeHit(const recob::Wire &wire, const geo::WireID &wireID, const geo::GeometryCore *geo, raw::TDCtick_t startTick, raw::TDCtick_t endTick) const;
  recob::Hit MakeHit(const raw::RawDigit &digit, const geo::WireID &wireID, const geo::GeometryCore *geo, raw::TDCtick_t startTick, raw::TDCtick_t endTick) const;

private:

  art::InputTag fTrackLabel;
  art::InputTag fWireLabel;
  std::vector<art::InputTag> fDigitLabels;
  art::InputTag fT0Label;
  float fSignalShapingWindow;
  float fTransverseDiffusionScale; 
  float fLongitudinalDiffusionScale;
  bool fAppendExtraHit;
  bool fPrependExtraHit;
  bool fUseWires;
  unsigned fCryostat;

};

sbn::TrackAreaHit::TrackAreaHit(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraTrack")),
    fWireLabel(p.get<art::InputTag>("WireLabel", "calowire")),
    fDigitLabels(p.get<std::vector<art::InputTag>>("DigitLabels", {})),
    fT0Label(p.get<art::InputTag>("T0Label", "")),
    fSignalShapingWindow(p.get<float>("SignalShapingWindow")),
    fTransverseDiffusionScale(p.get<float>("TransverseDiffusionScale")),
    fLongitudinalDiffusionScale(p.get<float>("LongitudinalDiffusionScale")),
    fAppendExtraHit(p.get<bool>("AppendExtraHit", false)),
    fPrependExtraHit(p.get<bool>("PrependExtraHit", false)),
    fUseWires(p.get<bool>("UseWires", true)),
    fCryostat(p.get<unsigned>("Cryostat", 0))
{
  produces<std::vector<recob::Hit>>();
  produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
  // produces<art::Assns<recob::Hit, recob::Wire>>();
}

recob::Hit sbn::TrackAreaHit::MakeHit(const recob::Wire &wire, const geo::WireID &wireID, const geo::GeometryCore *geo, raw::TDCtick_t startTick, raw::TDCtick_t endTick) const {
  // garbage values for some fields which are not applicable
  raw::ChannelID_t hitChannel = geo->PlaneWireToChannel(wireID);
  raw::TDCtick_t hitStartTick = -1000;
  bool setstartTick = false;
  raw::TDCtick_t hitEndTick = -1000;
  bool setendTick = false;
  float hitPeakTime = (startTick + endTick) / 2.; // recob::Hit time unit is also [ticks] 
  float hitSigmaPeakTime = (endTick - startTick) / 2.; 
  float hitRMS = (endTick - startTick) / 2.;
  float hitPeakAmplitude = 0.;
  float hitSigmaPeakAmplitude = 0.;
  float hitSumADC = 0.; // we're going to set this one
  float hitIntegral = 0.;  // this one too
  float hitSigmaIntegral = 0.; // TODO: this one?
  short int hitMultiplicity = 1; // no pulse trains here
  short int hitLocalIndex = 0;
  float hitGoodnessOfFit = -999999;
  short int hitNDF = 0;
  geo::View_t hitView = geo->View(wireID);
  geo::SigType_t hitSignalType = geo->SignalType(wireID);
  geo::WireID hitWire = wireID;
  
  //std::cout << "Hit view: " << hitView << " " << " sigType: " << hitSignalType << " wireID: " << wireID << std::endl;
  
  const recob::Wire::RegionsOfInterest_t &ROI = wire.SignalROI();
  for (lar::sparse_vector<float>::datarange_t range: ROI.get_ranges()) {
    raw::TDCtick_t roiStart = range.begin_index();
    raw::TDCtick_t roiEnd = range.end_index();
    
    //std::cout << "Window from: " << startTick << " to " << endTick << std::endl;
    //std::cout << "ROI from: " << roiStart << " to " << roiEnd << std::endl;
    // check for overlap
    if (roiStart <= endTick && startTick <= roiEnd) {
      // shrink the ROI to the integration window
      if (roiStart < startTick && startTick <= roiEnd) range.move_head(startTick);
      if (roiEnd > endTick && endTick >= roiStart) range.move_tail(endTick);

      // integrate
      for (float v: range) {
        hitSumADC += v;
        hitIntegral += v;
      }

      // update start and end time
      raw::TDCtick_t thisStart = std::max(roiStart, startTick);
      raw::TDCtick_t thisEnd = std::min(roiEnd, endTick);
      
      //std::cout << "Integrating from " << thisStart << " to " << thisEnd << std::endl;
      //std::cout << "Total integral: " << hitIntegral << std::endl;
      
      if (!setstartTick || thisStart < hitStartTick) {
        hitStartTick = thisStart;
        setstartTick = true;
      }
      if (!setendTick || thisEnd > hitEndTick) {
        hitEndTick = thisEnd;
        setendTick = true;
      }
    }
  }

  // if we didn't integrate anything don't save this hit
  if (!setstartTick || !setendTick) return recob::Hit();

  // construct the hit
  recob::Hit thisHit(
    hitChannel,
    hitStartTick,
    hitEndTick,
    hitPeakTime,
    hitSigmaPeakTime,
    hitRMS,
    hitPeakAmplitude,
    hitSigmaPeakAmplitude,
    hitSumADC,
    hitIntegral,
    hitSigmaIntegral,
    hitMultiplicity,
    hitLocalIndex,
    hitGoodnessOfFit,
    hitNDF,
    hitView,
    hitSignalType,
    hitWire
  );

  return thisHit;
}

recob::Hit sbn::TrackAreaHit::MakeHit(const raw::RawDigit &digit, const geo::WireID &wireID, const geo::GeometryCore *geo, raw::TDCtick_t startTick, raw::TDCtick_t endTick) const {

  // clamp the tick range inside the digits
  startTick = std::max(0, startTick);
  endTick = std::min((int)digit.NADC()-1, endTick);

  // garbage values for some fields which are not applicable
  raw::ChannelID_t hitChannel = geo->PlaneWireToChannel(wireID);
  raw::TDCtick_t hitStartTick = startTick;
  raw::TDCtick_t hitEndTick = endTick;
  float hitPeakTime = (startTick + endTick) / 2.; // recob::Hit time unit is also [ticks] 
  float hitSigmaPeakTime = (endTick - startTick) / 2.; 
  float hitRMS = (endTick - startTick) / 2.;
  float hitPeakAmplitude = 0.;
  float hitSigmaPeakAmplitude = 0.;
  float hitSumADC = 0.; // we're going to set this one
  float hitIntegral = 0.;  // this one too
  float hitSigmaIntegral = 0.; // TODO: this one?
  short int hitMultiplicity = 1; // no pulse trains here
  short int hitLocalIndex = 0;
  float hitGoodnessOfFit = -999999;
  short int hitNDF = 0;
  geo::View_t hitView = geo->View(wireID);
  geo::SigType_t hitSignalType = geo->SignalType(wireID);
  geo::WireID hitWire = wireID;
  
  for (int i_adc = startTick; i_adc <= endTick; i_adc++) {
    hitSumADC += digit.ADC(i_adc);
    hitIntegral += digit.ADC(i_adc);
  }

  // construct the hit
  recob::Hit thisHit(
    hitChannel,
    hitStartTick,
    hitEndTick,
    hitPeakTime,
    hitSigmaPeakTime,
    hitRMS,
    hitPeakAmplitude,
    hitSigmaPeakAmplitude,
    hitSumADC,
    hitIntegral,
    hitSigmaIntegral,
    hitMultiplicity,
    hitLocalIndex,
    hitGoodnessOfFit,
    hitNDF,
    hitView,
    hitSignalType,
    hitWire
  );

  return thisHit;
}

void sbn::TrackAreaHit::produce(art::Event& e)
{
  // output data products
  std::unique_ptr<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>> assn(new art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>);
  std::unique_ptr<std::vector<recob::Hit>> outHits(new std::vector<recob::Hit>);

  // If using wires
  std::unique_ptr<art::Assns<recob::Hit, recob::Wire>> wireAssn(new art::Assns<recob::Hit, recob::Wire>);
  // If using digits
  std::unique_ptr<art::Assns<recob::Hit, raw::RawDigit>> digitAssn(new art::Assns<recob::Hit, raw::RawDigit>);


  art::PtrMaker<recob::Hit> hitPtrMaker{e};

  // input data
  art::Handle<std::vector<recob::Track>> track_handle;
  e.getByLabel(fTrackLabel, track_handle);

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, track_handle);

  std::vector<art::Ptr<recob::Wire>> wires;
  std::vector<art::Ptr<raw::RawDigit>> digits;

  if (fUseWires) {
    auto const& wire_handle = e.getValidHandle<std::vector<recob::Wire>>(fWireLabel);
    art::fill_ptr_vector(wires, wire_handle);
  }
  else {
    for (const art::InputTag &l: fDigitLabels) {
      auto const& digit_handle = e.getValidHandle<std::vector<raw::RawDigit>>(l);
      art::fill_ptr_vector(digits, digit_handle);
    }
  }

  art::FindManyP<anab::T0> fmT0(tracks, e, fT0Label);

  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();
  auto const dclock = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, dclock);

  // auto const g4param = art::ServiceHandle<sim::LArG4Parameters const>();
  art::ServiceHandle<sim::LArG4Parameters const> g4param;

  std::set<unsigned> channels;

  for (unsigned i = 0; i < tracks.size(); i++) {
    const recob::Track &track = *tracks[i];

    // Look up if the track was shifted by some T0
    double track_t0 = 0.;
    if (!fT0Label.empty()) {
      const std::vector<art::Ptr<anab::T0>> &t0 = fmT0.at(i);
      if (t0.size()) {
        track_t0 = t0.at(0)->Time() / 1e3 /* ns -> us */;
      }
      // std::cout << "T0: " << t0.at(0)->Time() << " corr: " << track_t0 << std::endl;
    }

    // get the cryostat
    geo::CryostatID cid(fCryostat);
    for (auto const &planeID: geo->IteratePlaneIDs(cid)) {
      int lastWire = -100000;
      for (size_t j = 0; j < track.NumberTrajectoryPoints()-1; j++) {
        // go to the next valid point
        j = track.NextValidPoint(j);
        // break if bad/at end of trajectory
        if (!(j < track.NumberTrajectoryPoints()-1) || j == recob::TrackTrajectory::InvalidIndex) {
          break;
        }

        // get the next valid point as well -- if the next valid point is at the end, break
        size_t next_point_ind = track.NextValidPoint(j+1);
        // break if bad/at end of trajectory
        if (!(next_point_ind < track.NumberTrajectoryPoints()) || next_point_ind == recob::TrackTrajectory::InvalidIndex) {
          break;
        }

        recob::Track::Point_t thispoint_p = track.LocationAtPoint(j);
        recob::Track::Point_t nextpoint_p = track.LocationAtPoint(next_point_ind);
        TVector3 thispoint (thispoint_p.X(), thispoint_p.Y(), thispoint_p.Z());
        TVector3 nextpoint (nextpoint_p.X(), nextpoint_p.Y(), nextpoint_p.Z());

        // if we're not in the TPC cooresponding to this plane, continue
        if (!geo->TPC(planeID).ContainsPosition(thispoint)) continue;

        // calculate overlap with wires
        double thiswirecoord = geo->WireCoordinate(thispoint_p, planeID);         
        double nextwirecoord = geo->WireCoordinate(nextpoint_p, planeID);         

        int wireStart = std::nearbyint((thiswirecoord >= nextwirecoord) ? std::floor(thiswirecoord) : std::ceil(thiswirecoord));
        int wireEnd   = std::nearbyint((thiswirecoord >= nextwirecoord) ? std::floor(nextwirecoord) : std::ceil(nextwirecoord));

        // if we're not crossing a wire continue
        if (wireStart == wireEnd) continue;

        // check the validity of the range of wires
        if (!(wireStart >= 0 && wireStart < (int)geo->Plane(planeID).Nwires())) {
          std::cout << "At end of trajectory: " <<
            "Can't find wire for track trajectory position at: " << 
            thispoint.X() << " " << thispoint.Y() << " " << thispoint.Z() <<
            ". Returned start wire is " << wireStart << " (max: " << geo->Plane(planeID).Nwires() << ")" << std::endl;
          break;
        } 
        if (!(wireEnd >= 0 && wireEnd <= (int)geo->Plane(planeID).Nwires())) {
          std::cout << "At end of trajectory: " << 
            "Can't find wire for track trajectory position at: " << 
            nextpoint.X() << " " << nextpoint.Y() << " " << nextpoint.Z() << 
            ". Returned end wire is " << wireEnd << " (max: " << geo->Plane(planeID).Nwires() << ")" << std::endl;
          break;
        } 

        // if this wire is the same as the last skip it
        if (wireStart == lastWire) {
          if (wireStart < wireEnd) wireStart ++;
          else wireStart --;
        }

        // again continue if this wouldn't cross a new wire
        if (wireStart == wireEnd) continue;

        // get the integration window 
        float x = thispoint.X();
        recob::Track::Vector_t dir_v = track.DirectionAtPoint(j);
        TVector3 dir(dir_v.X(), dir_v.Y(), dir_v.Z());
        
        // Get time from X [ticks]
        //
        // Undo any T0 correction to the track to get back to the hits
        float time = dprop.ConvertXToTicks(x, planeID) + track_t0 * geo->TPC(planeID).DriftDir()[0] / dclock.TPCClock().TickPeriod();

        // get the geometric time window [ticks]
        //
        // First get the projected length in the direction towards the wireplane
        //
        // Effective pitch as in the Calorimetry module
        double angleToVert = geo->WireAngleToVertical(geo->View(planeID), planeID) - 0.5*::util::pi<>();
        double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() + std::cos(angleToVert)*dir.Z()); 
        double effpitch = geo->WirePitch() / cosgamma;

        // Add in the contribution from diffusion
        float xdrift = abs(x - geo->TPC(planeID).PlaneLocation(0)[0]); // Copied from G4 / DriftIonizationElectrons
        float tdrift = xdrift / dprop.DriftVelocity();
        effpitch += fTransverseDiffusionScale * sqrt(2*tdrift*g4param->TransverseDiffusion());

        double proj_length = abs(effpitch * dir.X() / sqrt(1 - dir.X() * dir.X()));
        // if the projected length is too big, clamp it to the track length
        proj_length = std::min(proj_length, abs(track.Length() * dir.X()));

        // convert the length to a window in ticks
        double window = (proj_length / dprop.DriftVelocity()) / dclock.TPCClock().TickPeriod();
        
        // add in the window from signal shaping [ticks]
        window += fSignalShapingWindow;
        // And longitudinal Diffsion
        window += (fLongitudinalDiffusionScale*sqrt(2*tdrift*g4param->LongitudinalDiffusion()) / dprop.DriftVelocity()) / dclock.TPCClock().TickPeriod();
        
        // convert to a TDCTick range
        raw::TDCtick_t startTick = std::nearbyint(std::floor(time - window / 2.));
        raw::TDCtick_t endTick = std::nearbyint(std::ceil(time + window / 2.));
        
        int incl = wireStart < wireEnd ? 1 : -1; 
        for (int wire = wireStart; wire != wireEnd; wire += incl) {
          //std::cout << "Got wire: " << wire << std::endl;
          geo::WireID wireID {planeID, (geo::WireID::WireID_t) wire};
          lastWire = wire;

          // Get the channel of this wire
          unsigned channel = geo->PlaneWireToChannel(wireID);

          // If we have already seen this channel, ignore it.
          //
          // This can happen in ICARUS because TPC channels are duplicated across 
          // multiple "logical" TPC's in planes 1 and 2. We don't want to 
          // double-count the charge on those channels
          if (channels.count(channel)) continue;

          channels.insert(channel);

          // Make the hit -- using wires or digits

          recob::Hit thisHit;
          int assoc_ind;
          if (fUseWires) {
            // First find the corresponding recob::Wire
            int wire_ind = -1;
            for (unsigned i_rwire = 0; i_rwire < wires.size(); i_rwire++) {
              if (channel == wires[i_rwire]->Channel()) {
                wire_ind = (int)i_rwire;
                break;
              }
            }
            assoc_ind = wire_ind;

            // make the hit
            thisHit = MakeHit(*wires.at(wire_ind), wireID, geo, startTick, endTick);
          }
          else {
            // Find the digits
            int digit_ind = -1;
            for (unsigned i_digit = 0; i_digit < digits.size(); i_digit++) {
              if (channel == digits[i_digit]->Channel()) {
                digit_ind = (int)i_digit;
                break;
              }
            }
            assoc_ind = digit_ind;

            // Make the hit
            thisHit = MakeHit(*digits.at(digit_ind), wireID, geo, startTick, endTick);
          }

          // if invalid hit, ignore it
          if (thisHit.WireID().Wire == geo::WireID::InvalidID) continue; 

          // construct the track meta object
          recob::TrackHitMeta thisHitMeta(j);

          // save to output data
          outHits->push_back(thisHit);

          // add to the association
          art::Ptr<recob::Hit> thisHitPtr = hitPtrMaker(outHits->size()-1);
          assn->addSingle(tracks[i], thisHitPtr, thisHitMeta);
          if (fUseWires) {
            wireAssn->addSingle(thisHitPtr, wires[assoc_ind]);
          }
          else {
            digitAssn->addSingle(thisHitPtr, digits[assoc_ind]);
          }
        }
      } 
    }
  }

  e.put(std::move(outHits));
  e.put(std::move(assn));
  if (fUseWires) {
    e.put(std::move(wireAssn));
  }
  else {
    e.put(std::move(digitAssn));
  }

}

DEFINE_ART_MODULE(sbn::TrackAreaHit)
