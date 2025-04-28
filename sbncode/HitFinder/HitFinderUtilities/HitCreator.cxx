/** ****************************************************************************
 * @file   HitCreator.cxx
 * @brief  Helper functions to create a hit from ChannelROIs - implementation file
 * @date   April 20, 2025
 * @author usher@slac.stanford.edu
 * @see    Hit.h HitCreator.h
 *
 * ****************************************************************************/

// declaration header
#include "sbncode/HitFinder/HitFinderUtilities/HitCreator.h"

// C/C++ standard library
#include <algorithm> // std::accumulate(), std::max()
#include <cassert>
#include <limits>  // std::numeric_limits<>
#include <utility> // std::move()

// art libraries
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/Exception.h"

// LArSoft libraries
#include "larcore/Geometry/WireReadout.h"
#include "lardata/Utilities/MakeIndex.h"

namespace {

  /// Erases the content of an association
  template <typename Left, typename Right, typename Metadata>
  void ClearAssociations(art::Assns<Left, Right, Metadata>& assns)
  {
    art::Assns<Left, Right, Metadata> empty;
    assns.swap(empty);
  } // ClearAssociations()

} // local namespace

/// Reconstruction base classes
namespace sbn {

  //****************************************************************************
  //***  HitCreator
  //----------------------------------------------------------------------
  HitCreator::HitCreator(raw::RawDigit const& digits,
                         geo::WireID const& wireID,
                         raw::TDCtick_t start_tick,
                         raw::TDCtick_t end_tick,
                         float rms,
                         float peak_time,
                         float sigma_peak_time,
                         float peak_amplitude,
                         float sigma_peak_amplitude,
                         float hit_integral,
                         float hit_sigma_integral,
                         float ROIsummedADC,
                         float HitsummedADC,
                         short int multiplicity,
                         short int local_index,
                         float goodness_of_fit,
                         int dof)
    : hit(digits.Channel(),
          start_tick,
          end_tick,
          peak_time,
          sigma_peak_time,
          rms,
          peak_amplitude,
          sigma_peak_amplitude,
          ROIsummedADC,
          HitsummedADC,
          hit_integral,
          hit_sigma_integral,
          multiplicity,
          local_index,
          goodness_of_fit,
          dof,
          art::ServiceHandle<geo::WireReadout const>()->Get().View(digits.Channel()),
          art::ServiceHandle<geo::WireReadout const>()->Get().SignalType(digits.Channel()),
          wireID)
  {} // HitCreator::HitCreator(RawDigit)

  //----------------------------------------------------------------------
  HitCreator::HitCreator(recob::ChannelROI const& wire,
                         geo::WireID const& wireID,
                         raw::TDCtick_t start_tick,
                         raw::TDCtick_t end_tick,
                         float rms,
                         float peak_time,
                         float sigma_peak_time,
                         float peak_amplitude,
                         float sigma_peak_amplitude,
                         float hit_integral,
                         float hit_sigma_integral,
                         float ROIsummedADC,
                         float HitsummedADC,
                         short int multiplicity,
                         short int local_index,
                         float goodness_of_fit,
                         int dof)
    : hit(wire.Channel(),
          start_tick,
          end_tick,
          peak_time,
          sigma_peak_time,
          rms,
          peak_amplitude,
          sigma_peak_amplitude,
          ROIsummedADC,
          HitsummedADC,
          hit_integral,
          hit_sigma_integral,
          multiplicity,
          local_index,
          goodness_of_fit,
          dof,
          art::ServiceHandle<geo::WireReadout const>()->Get().View(wire.Channel()),
          art::ServiceHandle<geo::WireReadout const>()->Get().SignalType(wire.Channel()),
          wireID)
  {} // HitCreator::HitCreator(Wire)

  //----------------------------------------------------------------------
  HitCreator::HitCreator(recob::ChannelROI const& wire,
                         geo::WireID const& wireID,
                         raw::TDCtick_t start_tick,
                         raw::TDCtick_t end_tick,
                         float rms,
                         float peak_time,
                         float sigma_peak_time,
                         float peak_amplitude,
                         float sigma_peak_amplitude,
                         float hit_integral,
                         float hit_sigma_integral,
                         short int multiplicity,
                         short int local_index,
                         float goodness_of_fit,
                         int dof)
    : HitCreator(wire,
                 wireID,
                 start_tick,
                 end_tick,
                 rms,
                 peak_time,
                 sigma_peak_time,
                 peak_amplitude,
                 sigma_peak_amplitude,
                 hit_integral,
                 hit_sigma_integral,
                 std::accumulate(wire.SignalROI().begin() + start_tick,
                                 wire.SignalROI().begin() + end_tick,
                        0.), // sum of ADC counts between start_tick and end_tic
        std::accumulate(
          wire.SignalROI().begin() + start_tick,
          wire.SignalROI().begin() + end_tick,
          0.), // sum of ADC counts between start_tick and end_tick TO BE MODIFIED IF HIT HAS TO BE CONSIDERED
                 multiplicity,
                 local_index,
                 goodness_of_fit,
                 dof)
  {} // HitCreator::HitCreator(Wire; no summed ADC)

  //----------------------------------------------------------------------
  HitCreator::HitCreator(recob::ChannelROI const& wire,
                         geo::WireID const& wireID,
                         float rms,
                         float peak_time,
                         float sigma_peak_time,
                         float peak_amplitude,
                         float sigma_peak_amplitude,
                         float hit_integral,
                         float hit_sigma_integral,
                         float ROIsummedADC,
                         float HitsummedADC,
                         short int multiplicity,
                         short int local_index,
                         float goodness_of_fit,
                         int dof,
                         RegionOfInterest_t const& signal)
    : HitCreator(wire,
                 wireID,
                 signal.begin_index(),
                 signal.end_index(),
                 rms,
                 peak_time,
                 sigma_peak_time,
                 peak_amplitude,
                 sigma_peak_amplitude,
                 hit_integral,
                 hit_sigma_integral,
                 ROIsummedADC,
                 HitsummedADC,
                 multiplicity,
                 local_index,
                 goodness_of_fit,
                 dof)
  {} // HitCreator::HitCreator(Wire; RoI)

  //----------------------------------------------------------------------
  HitCreator::HitCreator(recob::ChannelROI const& wire,
                         geo::WireID const& wireID,
                         float rms,
                         float peak_time,
                         float sigma_peak_time,
                         float peak_amplitude,
                         float sigma_peak_amplitude,
                         float hit_integral,
                         float hit_sigma_integral,
                         float ROIsummedADC,
                         float HitsummedADC,
                         short int multiplicity,
                         short int local_index,
                         float goodness_of_fit,
                         int dof,
                         size_t iSignalRoI)
    : HitCreator(wire,
                 wireID,
                 rms,
                 peak_time,
                 sigma_peak_time,
                 peak_amplitude,
                 sigma_peak_amplitude,
                 hit_integral,
                 hit_sigma_integral,
                 ROIsummedADC,
                 HitsummedADC,
                 multiplicity,
                 local_index,
                 goodness_of_fit,
                 dof,
                 wire.SignalROI().range(iSignalRoI))
  {} // HitCreator::HitCreator(Wire; RoI index)

  HitCreator::HitCreator(recob::Hit const& from) : hit(from) {}

  // Need to do this way since we don't have friendship for the wireID
  HitCreator::HitCreator(recob::Hit const& from, geo::WireID const& wireID) 
    : hit(from.Channel(),
          from.StartTick(),
          from.EndTick(),
          from.PeakTime(),
          from.SigmaPeakTime(),
          from.RMS(),
          from.PeakAmplitude(),
          from.SigmaPeakAmplitude(),
          from.ROISummedADC(),
          from.HitSummedADC(),
          from.Integral(),
          from.SigmaIntegral(),
          from.Multiplicity(),
          from.LocalIndex(),
          from.GoodnessOfFit(),
          from.DegreesOfFreedom(),
          art::ServiceHandle<geo::WireReadout const>()->Get().View(from.Channel()),
          art::ServiceHandle<geo::WireReadout const>()->Get().SignalType(from.Channel()),
          wireID)
  { } // HitCreator::HitCreator(new wire ID)

  //****************************************************************************
  //***  HitAndAssociationsWriterBase
  //----------------------------------------------------------------------
  HitAndAssociationsWriterBase::HitAndAssociationsWriterBase(art::Event& event,
                                                             std::string instance_name,
                                                             bool doWireAssns,
                                                             bool doRawDigitAssns)
    : prod_instance(instance_name)
    , hits()
    , WireAssns(doWireAssns ? new art::Assns<recob::ChannelROI, recob::Hit> : nullptr)
    , RawDigitAssns(doRawDigitAssns ? new art::Assns<raw::RawDigit, recob::Hit> : nullptr)
    , event(&event)
    , hitPtrMaker(*(this->event), prod_instance)
  {} // HitAndAssociationsWriterBase::HitAndAssociationsWriterBase()

  //------------------------------------------------------------------------------
  void HitAndAssociationsWriterBase::declare_products(art::ProducesCollector& collector,
                                                      std::string instance_name /* = "" */,
                                                      bool doWireAssns /* = true */,
                                                      bool doRawDigitAssns /* = true */
  )
  {
    collector.produces<std::vector<recob::Hit>>(instance_name);

    // declare the other products we are creating (if any)
    if (doWireAssns) { collector.produces<art::Assns<recob::ChannelROI, recob::Hit>>(instance_name); }
    if (doRawDigitAssns) {
      collector.produces<art::Assns<raw::RawDigit, recob::Hit>>(instance_name);
    }
  } // HitAndAssociationsWriterBase::declare_products()

  //------------------------------------------------------------------------------
  void HitAndAssociationsWriterBase::put_into()
  {
    assert(event);
    if (hits) event->put(std::move(hits), prod_instance);
    if (WireAssns) event->put(std::move(WireAssns), prod_instance);
    if (RawDigitAssns) event->put(std::move(RawDigitAssns), prod_instance);
  } // HitAndAssociationsWriterBase::put_into()

  //****************************************************************************
  //***  HitCollectionCreator
  //----------------------------------------------------------------------
  HitCollectionCreator::HitCollectionCreator(art::Event& event,
                                             std::string instance_name /* = "" */,
                                             bool doWireAssns /* = true */,
                                             bool doRawDigitAssns /* = true */
                                             )
    : HitAndAssociationsWriterBase(event, instance_name, doWireAssns, doRawDigitAssns)
  {
    hits.reset(new std::vector<recob::Hit>);
  } // HitCollectionCreator::HitCollectionCreator()

  //----------------------------------------------------------------------
  void HitCollectionCreator::emplace_back(recob::Hit&& hit,
                                          art::Ptr<recob::ChannelROI> const& wire,
                                          art::Ptr<raw::RawDigit> const& digits)
  {

    // add the hit to the collection
    hits->emplace_back(std::move(hit));

    CreateAssociationsToLastHit(wire, digits);
  } // HitCollectionCreator::emplace_back(Hit&&)

  //----------------------------------------------------------------------
  void HitCollectionCreator::emplace_back(recob::Hit const& hit,
                                          art::Ptr<recob::ChannelROI> const& wire,
                                          art::Ptr<raw::RawDigit> const& digits)
  {

    // add the hit to the collection
    hits->push_back(hit);

    CreateAssociationsToLastHit(wire, digits);
  } // HitCollectionCreator::emplace_back(Hit)

  //----------------------------------------------------------------------
  void HitCollectionCreator::put_into()
  {
    if (!hits) {
      throw art::Exception(art::errors::LogicError)
        << "HitCollectionCreator is trying to put into the event"
           " a hit collection that was never created!\n";
    }
    HitAndAssociationsWriterBase::put_into();
  } // HitCollectionCreator::put_into()

  //----------------------------------------------------------------------
  void HitCollectionCreator::CreateAssociationsToLastHit(art::Ptr<recob::ChannelROI> const& wire,
                                                         art::Ptr<raw::RawDigit> const& digits)
  {
    // if no association is required, we are done
    if (!WireAssns && !RawDigitAssns) return;

    // art pointer to the hit we just created
    HitPtr_t hit_ptr(CreatePtrToLastHit());

    // association with wires
    if (WireAssns && wire.isNonnull())
      WireAssns->addSingle(wire, hit_ptr); // if it fails, it throws

    // association with wires
    if (RawDigitAssns && digits.isNonnull())
      RawDigitAssns->addSingle(digits, hit_ptr); // if it fails, it throws

  } // HitCollectionCreator::CreateAssociationsToLastHit()

  //****************************************************************************
  //***  HitCollectionAssociator
  //----------------------------------------------------------------------
  HitCollectionAssociator::HitCollectionAssociator(art::Event& event,
                                                   std::string instance_name,
                                                   art::InputTag const& WireModuleLabel,
                                                   art::InputTag const& RawDigitModuleLabel)
    : HitAndAssociationsWriterBase(event,
                                   instance_name,
                                   WireModuleLabel != "",
                                   RawDigitModuleLabel != "")
    , wires_label(WireModuleLabel)
    , digits_label(RawDigitModuleLabel)
  {
    hits.reset(new std::vector<recob::Hit>);
  } // HitCollectionAssociator::HitCollectionAssociator()

  //----------------------------------------------------------------------
  sbn::HitCollectionAssociator::HitCollectionAssociator(art::Event& event,
                                                        std::string instance_name,
                                                        art::InputTag const& WireModuleLabel,
                                                        bool doRawDigitAssns /* = false */
                                                        )
    : HitAndAssociationsWriterBase(event, instance_name, WireModuleLabel != "", doRawDigitAssns)
    , wires_label(WireModuleLabel)
    , digits_label()
  {
    if (RawDigitAssns && !WireAssns) {
      throw art::Exception(art::errors::LogicError)
        << "HitCollectionAssociator can't create hit <--> raw digit"
           " associations through wires, without wires!\n";
    }
    hits.reset(new std::vector<recob::Hit>);
  } // HitCollectionAssociator::HitCollectionAssociator()

  //----------------------------------------------------------------------
  void HitCollectionAssociator::use_hits(std::unique_ptr<std::vector<recob::Hit>>&& srchits)
  {
    hits = std::move(srchits);
  } // HitCollectionAssociator::use_hits()

  //----------------------------------------------------------------------
  void HitCollectionAssociator::put_into()
  {
    prepare_associations();
    HitAndAssociationsWriterBase::put_into();
  } // HitCollectionAssociator::put_into()

  //----------------------------------------------------------------------
  void HitCollectionAssociator::prepare_associations(std::vector<recob::Hit> const& srchits)
  {
    if (!RawDigitAssns && !WireAssns) return; // no associations needed
    assert(event);

    // we make the associations anew
    if (RawDigitAssns) ClearAssociations(*RawDigitAssns);
    if (WireAssns) ClearAssociations(*WireAssns);

    // the following is true is we want associations with digits
    // but we don't know where digits are; in that case, we try to use wires
    const bool bUseWiresForDigits = RawDigitAssns && (digits_label == "");

    if (WireAssns || bUseWiresForDigits) {
      // do we use wires for digit associations too?

      // get the wire collection
      art::ValidHandle<std::vector<recob::ChannelROI>> hWires =
        event->getValidHandle<std::vector<recob::ChannelROI>>(wires_label);

      // fill a map of wire index vs. channel number
      std::vector<size_t> WireMap = util::MakeIndex(*hWires, std::mem_fn(&recob::ChannelROI::Channel));

      // use raw rigit - wire association, assuming they have been produced
      // by the same producer as the wire and with the same instance name;
      // we don't check whether the data product is found, but the following
      // code will have FindOneP throw if that was not the case
      // (that's what we would do here anyway, maybe with a better message...)
      std::unique_ptr<art::FindOneP<raw::RawDigit>> WireToDigit;
      if (bUseWiresForDigits) {
        WireToDigit.reset(new art::FindOneP<raw::RawDigit>(hWires, *event, wires_label));
      }

      // add associations, hit by hit:
      for (size_t iHit = 0; iHit < srchits.size(); ++iHit) {

        // find the channel
        size_t iChannel = size_t(srchits[iHit].Channel()); // forcibly converted

        // find the wire associated to that channel
        size_t iWire = std::numeric_limits<size_t>::max();
        if (iChannel < WireMap.size()) iWire = WireMap[iChannel];
        if (iWire == std::numeric_limits<size_t>::max()) {
          throw art::Exception(art::errors::LogicError)
            << "No wire associated to channel #" << iChannel << " whence hit #" << iHit
            << " comes!\n";
        } // if no channel

        // make the association with wires
        if (WireAssns) {
          art::Ptr<recob::ChannelROI> wire(hWires, iWire);
          WireAssns->addSingle(wire, CreatePtr(iHit));
        }

        if (bUseWiresForDigits) {
          // find the digit associated to that channel
          art::Ptr<raw::RawDigit> const& digit = WireToDigit->at(iWire);
          if (digit.isNull()) {
            throw art::Exception(art::errors::LogicError)
              << "No raw digit associated to channel #" << iChannel << " whence hit #" << iHit
              << " comes!\n";
          } // if no channel

          // make the association
          RawDigitAssns->addSingle(digit, CreatePtr(iHit));
        } // if create digit associations through wires
      }   // for hit

    } // if wire associations

    if (RawDigitAssns && !bUseWiresForDigits) {
      // get the digit collection
      art::ValidHandle<std::vector<raw::RawDigit>> hDigits =
        event->getValidHandle<std::vector<raw::RawDigit>>(digits_label);

      // fill a map of wire index vs. channel number
      std::vector<size_t> DigitMap =
        util::MakeIndex(*hDigits, std::mem_fn(&raw::RawDigit::Channel));

      // add associations, hit by hit:
      for (size_t iHit = 0; iHit < srchits.size(); ++iHit) {

        // find the channel
        size_t iChannel = size_t(srchits[iHit].Channel()); // forcibly converted

        // find the digit associated to that channel
        size_t iDigit = std::numeric_limits<size_t>::max();
        if (iChannel < DigitMap.size()) iDigit = DigitMap[iChannel];
        if (iDigit == std::numeric_limits<size_t>::max()) {
          throw art::Exception(art::errors::LogicError)
            << "No raw digit associated to channel #" << iChannel << " whence hit #" << iHit
            << " comes!\n";
        } // if no channel

        // make the association
        art::Ptr<raw::RawDigit> digit(hDigits, iDigit);
        RawDigitAssns->addSingle(digit, CreatePtr(iHit));

      } // for hit
    }   // if we have rawdigit label

  } // HitCollectionAssociator::put_into()

  //****************************************************************************
  //***  HitRefinerAssociator
  //----------------------------------------------------------------------
  HitRefinerAssociator::HitRefinerAssociator(art::Event& event,
                                             art::InputTag const& HitModuleLabel,
                                             std::string instance_name /* = "" */,
                                             bool doWireAssns /* = true */,
                                             bool doRawDigitAssns /* = true */
                                             )
    : HitAndAssociationsWriterBase(event, instance_name, doWireAssns, doRawDigitAssns)
    , hits_label(HitModuleLabel)
  {
    hits.reset(new std::vector<recob::Hit>);
  } // HitRefinerAssociator::HitRefinerAssociator()

  //----------------------------------------------------------------------
  void HitRefinerAssociator::use_hits(std::unique_ptr<std::vector<recob::Hit>>&& srchits)
  {
    hits = std::move(srchits);
  } // HitRefinerAssociator::use_hits()

  //----------------------------------------------------------------------
  void HitRefinerAssociator::put_into()
  {
    prepare_associations();
    HitAndAssociationsWriterBase::put_into();
  } // HitRefinerAssociator::put_into()

  //----------------------------------------------------------------------
  void HitRefinerAssociator::prepare_associations(std::vector<recob::Hit> const& srchits)
  {
    if (!RawDigitAssns && !WireAssns) return; // no associations needed
    assert(event);

    // we make the associations anew
    if (RawDigitAssns) ClearAssociations(*RawDigitAssns);

    // read the hits; this is going to hurt performances...
    // no solution to that until there is a way to have a lazy read
    art::ValidHandle<std::vector<recob::Hit>> hHits =
      event->getValidHandle<std::vector<recob::Hit>>(hits_label);

    // now get the associations
    if (WireAssns) {
      // we make the associations anew
      ClearAssociations(*WireAssns);

      // find the associations between the hits and the wires
      art::FindOneP<recob::ChannelROI> HitToWire(hHits, *event, hits_label);
      if (!HitToWire.isValid()) {
        throw art::Exception(art::errors::ProductNotFound)
          << "Can't find the associations between hits and wires produced by '" << hits_label
          << "'!\n";
      } // if no association

      // fill a map of wire vs. channel number
      std::vector<art::Ptr<recob::ChannelROI>> WireMap;
      for (size_t iAssn = 0; iAssn < HitToWire.size(); ++iAssn) {
        art::Ptr<recob::ChannelROI> wire = HitToWire.at(iAssn);
        if (wire.isNull()) continue;
        size_t channelID = (size_t)wire->Channel();
        if (WireMap.size() <= channelID) // expand the map of necessary
          WireMap.resize(std::max(channelID + 1, 2 * WireMap.size()), {});
        WireMap[channelID] = std::move(wire);
      } // for

      // now go through all the hits...
      for (size_t iHit = 0; iHit < srchits.size(); ++iHit) {
        recob::Hit const& hit = srchits[iHit];
        size_t channelID = (size_t)hit.Channel();

        // no association if there is no wire to associate with
        if ((channelID >= WireMap.size()) || !WireMap[channelID]) continue;

        // create an association using the same wire pointer
        WireAssns->addSingle(WireMap[channelID], CreatePtr(iHit));
      } // for hits
    }   // if wire associations

    // now get the associations
    if (RawDigitAssns) {
      // we make the associations anew
      ClearAssociations(*RawDigitAssns);

      // find the associations between the hits and the raw digits
      art::FindOneP<raw::RawDigit> HitToDigits(hHits, *event, hits_label);
      if (!HitToDigits.isValid()) {
        throw art::Exception(art::errors::ProductNotFound)
          << "Can't find the associations between hits and raw digits"
          << " produced by '" << hits_label << "'!\n";
      } // if no association

      // fill a map of digits vs. channel number
      std::vector<art::Ptr<raw::RawDigit>> DigitMap;
      for (size_t iAssn = 0; iAssn < HitToDigits.size(); ++iAssn) {
        art::Ptr<raw::RawDigit> digits = HitToDigits.at(iAssn);
        if (digits.isNull()) continue;
        size_t channelID = (size_t)digits->Channel();
        if (DigitMap.size() <= channelID) // expand the map of necessary
          DigitMap.resize(std::max(channelID + 1, 2 * DigitMap.size()), {});
        DigitMap[channelID] = std::move(digits);
      } // for

      // now go through all the hits...
      for (size_t iHit = 0; iHit < srchits.size(); ++iHit) {
        recob::Hit const& hit = srchits[iHit];
        size_t channelID = (size_t)hit.Channel();

        // no association if there is no digits to associate with
        if ((channelID >= DigitMap.size()) || !DigitMap[channelID]) continue;

        // create an association using the same digits pointer
        RawDigitAssns->addSingle(DigitMap[channelID], CreatePtr(iHit));
      } // for hits
    }   // if digit associations

  } // HitRefinerAssociator::put_into()

} // namespace recob
