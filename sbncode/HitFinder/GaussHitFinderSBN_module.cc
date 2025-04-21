////////////////////////////////////////////////////////////////////////
//
// GaussHitFinder class
//
// jaasaadi@syr.edu
//
//  This algorithm is designed to find hits on wires after deconvolution.
// -----------------------------------
// This algorithm is based on the FFTHitFinder written by Brian Page,
// Michigan State University, for the ArgoNeuT experiment.
//
//
// The algorithm walks along the wire and looks for pulses above threshold
// The algorithm then attempts to fit n-gaussians to these pulses where n
// is set by the number of peaks found in the pulse
// If the Chi2/NDF returned is "bad" it attempts to fit n+1 gaussians to
// the pulse. If this is a better fit it then uses the parameters of the
// Gaussian fit to characterize the "hit" object
//
// To use this simply include the following in your producers:
// gaushit:     @local::microboone_GaussHitFinderSBN
// gaushit:	@local::argoneut_GaussHitFinderSBN
////////////////////////////////////////////////////////////////////////

// C/C++ standard library
#include <algorithm> // std::accumulate()
#include <atomic>
#include <memory> // std::unique_ptr()
#include <string>
#include <utility> // std::move()

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/Globals.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft Includes
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larreco/HitFinder/HitFilterAlg.h"

#include "sbnobj/ICARUS/TPC/ChannelROI.h"
#include "sbnobj/Common/Utilities//HitCreator.h"
#include "sbnobj/ICARUS/TPC/Hit.h"

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"

// ROOT Includes
#include "TH1F.h"
#include "TMath.h"

#include "tbb/concurrent_vector.h"
#include "tbb/parallel_for.h"

namespace hit {
  class GaussHitFinderSBN : public art::SharedProducer {
  public:
    explicit GaussHitFinderSBN(fhicl::ParameterSet const& pset, art::ProcessingFrame const&);

  private:
    void produce(art::Event& evt, art::ProcessingFrame const&) override;
    void beginJob(art::ProcessingFrame const&) override;

    std::vector<double> FillOutHitParameterVector(const std::vector<double>& input);

    const bool fFilterHits;
    const bool fFillHists;

    const std::string fCalDataModuleLabel;
    const std::string fAllHitsInstanceName;

    const std::vector<int> fLongMaxHitsVec;    ///<Maximum number hits on a really long pulse train
    const std::vector<int> fLongPulseWidthVec; ///<Sets width of hits used to describe long pulses

    const size_t fMaxMultiHit; ///<maximum hits for multi fit
    const int fAreaMethod;     ///<Type of area calculation
    const std::vector<double>
      fAreaNormsVec;       ///<factors for converting area to same units as peak height
    const double fChi2NDF; ///maximum Chisquared / NDF allowed for a hit to be saved

    const std::vector<float> fPulseHeightCuts;
    const std::vector<float> fPulseWidthCuts;
    const std::vector<float> fPulseRatioCuts;

    std::atomic<size_t> fEventCount{0};

    //only Standard and Morphological implementation is threadsafe.
    std::vector<std::unique_ptr<reco_tool::ICandidateHitFinder>>
      fHitFinderToolVec; ///< For finding candidate hits
    // only Marqdt implementation is threadsafe.
    std::unique_ptr<reco_tool::IPeakFitter> fPeakFitterTool; ///< Perform fit to candidate peaks
    //HitFilterAlg implementation is threadsafe.
    std::unique_ptr<HitFilterAlg> fHitFilterAlg; ///< algorithm used to filter out noise hits

    //only used when fFillHists is true and in single threaded mode.
    TH1F* fFirstChi2;
    TH1F* fChi2;

  }; // class GaussHitFinderSBN

  //-------------------------------------------------
  //-------------------------------------------------
  GaussHitFinderSBN::GaussHitFinderSBN(fhicl::ParameterSet const& pset, art::ProcessingFrame const&)
    : SharedProducer{pset}
    , fFilterHits(pset.get<bool>("FilterHits", false))
    , fFillHists(pset.get<bool>("FillHists", false))
    , fCalDataModuleLabel(pset.get<std::string>("CalDataModuleLabel"))
    , fAllHitsInstanceName(pset.get<std::string>("AllHitsInstanceName", ""))
    , fLongMaxHitsVec(pset.get<std::vector<int>>("LongMaxHits", std::vector<int>() = {25, 25, 25}))
    , fLongPulseWidthVec(
        pset.get<std::vector<int>>("LongPulseWidth", std::vector<int>() = {16, 16, 16}))
    , fMaxMultiHit(pset.get<int>("MaxMultiHit"))
    , fAreaMethod(pset.get<int>("AreaMethod"))
    , fAreaNormsVec(FillOutHitParameterVector(pset.get<std::vector<double>>("AreaNorms")))
    , fChi2NDF(pset.get<double>("Chi2NDF"))
    , fPulseHeightCuts(
        pset.get<std::vector<float>>("PulseHeightCuts", std::vector<float>() = {3.0, 3.0, 3.0}))
    , fPulseWidthCuts(
        pset.get<std::vector<float>>("PulseWidthCuts", std::vector<float>() = {2.0, 1.5, 1.0}))
    , fPulseRatioCuts(
        pset.get<std::vector<float>>("PulseRatioCuts", std::vector<float>() = {0.35, 0.40, 0.20}))
  {
    if (fFillHists && art::Globals::instance()->nthreads() > 1u) {
      throw art::Exception(art::errors::Configuration)
        << "Cannot fill histograms when multiple threads configured, please set fFillHists to "
           "false or change number of threads to 1\n";
    }
    async<art::InEvent>();
    if (fFilterHits) {
      fHitFilterAlg = std::make_unique<HitFilterAlg>(pset.get<fhicl::ParameterSet>("HitFilterAlg"));
    }

    // recover the tool to do the candidate hit finding
    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& hitFinderTools = pset.get<fhicl::ParameterSet>("HitFinderToolVec");

    fHitFinderToolVec.resize(hitFinderTools.get_pset_names().size());

    for (const std::string& hitFinderTool : hitFinderTools.get_pset_names()) {
      const fhicl::ParameterSet& hitFinderToolParamSet =
        hitFinderTools.get<fhicl::ParameterSet>(hitFinderTool);
      size_t planeIdx = hitFinderToolParamSet.get<size_t>("Plane");

      fHitFinderToolVec.at(planeIdx) =
        art::make_tool<reco_tool::ICandidateHitFinder>(hitFinderToolParamSet);
    }

    // Recover the peak fitting tool
    fPeakFitterTool =
      art::make_tool<reco_tool::IPeakFitter>(pset.get<fhicl::ParameterSet>("PeakFitter"));

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // We want the option to output two hit collections, one filtered
    // and one with all hits. The key to doing this will be a non-null
    // instance name for the second collection
    // (with no particular product label)
    sbn::HitCollectionCreator::declare_products(
      producesCollector(), fAllHitsInstanceName, true, false); //fMakeRawDigitAssns);

    // and now the filtered hits...
    if (fAllHitsInstanceName != "")
      sbn::HitCollectionCreator::declare_products(
        producesCollector(), "", true, false); //fMakeRawDigitAssns);

    return;
  } // GaussHitFinderSBN::GaussHitFinderSBN()

  //-------------------------------------------------
  //-------------------------------------------------
  std::vector<double> GaussHitFinderSBN::FillOutHitParameterVector(const std::vector<double>& input)
  {
    if (input.size() == 0)
      throw std::runtime_error(
        "GaussHitFinderSBN::FillOutHitParameterVector ERROR! Input config vector has zero size.");

    std::vector<double> output;
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    const unsigned int N_PLANES = wireReadoutGeom.Nplanes();

    if (input.size() == 1)
      output.resize(N_PLANES, input[0]);
    else if (input.size() == N_PLANES)
      output = input;
    else
      throw std::runtime_error("GaussHitFinderSBN::FillOutHitParameterVector ERROR! Input config "
                               "vector size !=1 and !=N_PLANES.");
    return output;
  }

  //-------------------------------------------------
  //-------------------------------------------------
  void GaussHitFinderSBN::beginJob(art::ProcessingFrame const&)
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;

    // ======================================
    // === Hit Information for Histograms ===
    if (fFillHists) {
      fFirstChi2 = tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 5000);
      fChi2 = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 5000);
    }
  }

  //  This algorithm uses the fact that deconvolved signals are very smooth
  //  and looks for hits as areas between local minima that have signal above
  //  threshold.
  //-------------------------------------------------
  void GaussHitFinderSBN::produce(art::Event& evt, art::ProcessingFrame const&)
  {
    unsigned int count = fEventCount.fetch_add(1);
    //==================================================================================================

    TH1::AddDirectory(kFALSE);

    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

    // ###############################################
    // ### Making a ptr vector to put on the event ###
    // ###############################################
    // this contains the hit collection
    // and its associations to wires and raw digits
    sbn::HitCollectionCreator allHitCol(evt, fAllHitsInstanceName, true, false);

    // Handle the filtered hits collection...
    sbn::HitCollectionCreator hcol(evt, "", true, false);
    sbn::HitCollectionCreator* filteredHitCol = 0;

    if (fFilterHits) filteredHitCol = &hcol;

    //store in a thread safe way
    struct hitstruct {
      sbn::Hit hit_tbb;
      art::Ptr<recob::ChannelROI> channelROI_tbb;
    };

    tbb::concurrent_vector<hitstruct> hitstruct_vec;
    tbb::concurrent_vector<hitstruct> filthitstruct_vec;

    // ##########################################
    // ### Reading in the Wire List object(s) ###
    // ##########################################
    art::Handle<std::vector<recob::ChannelROI>> channelROIVecHandle;
    evt.getByLabel(fCalDataModuleLabel, channelROIVecHandle);

    //#################################################
    //###    Set the charge determination method    ###
    //### Default is to compute the normalized area ###
    //#################################################
    std::function<double(double, double, double, double, int, int)> chargeFunc =
      [](double /* peakMean */,
         double peakAmp,
         double peakWidth,
         double areaNorm,
         int /* low */,
         int /* hi */) { return std::sqrt(2 * TMath::Pi()) * peakAmp * peakWidth / areaNorm; };

    //##############################################
    //### Alternative is to integrate over pulse ###
    //##############################################
    if (fAreaMethod == 0)
      chargeFunc = [](double peakMean,
                      double peakAmp,
                      double peakWidth,
                      double /* areaNorm */,
                      int low,
                      int hi) {
        double charge(0);
        for (int sigPos = low; sigPos < hi; sigPos++)
          charge += peakAmp * TMath::Gaus(sigPos, peakMean, peakWidth);
        return charge;
      };

    //##############################
    //### Looping over the wires ###
    //##############################
    tbb::parallel_for(
      static_cast<std::size_t>(0),
      channelROIVecHandle->size(),
      [&](size_t& channelROIIter) {
        // ####################################
        // ### Getting this particular channelROI ###
        // ####################################
        art::Ptr<recob::ChannelROI> channelROI(channelROIVecHandle, channelROIIter);

        // --- Setting Channel Number and Signal type ---

        raw::ChannelID_t channel = channelROI->Channel();

        // get the WireID for this hit
        std::vector<geo::WireID> wids = wireReadoutGeom.ChannelToWire(channel);
        // for now, just take the first option returned from ChannelToWire
        geo::WireID wid = wids[0];
        // We need to know the plane to look up parameters
        geo::PlaneID::PlaneID_t plane = wid.Plane;

        // ----------------------------------------------------------
        // -- Setting the appropriate signal widths and thresholds --
        // --    for the right plane.      --
        // ----------------------------------------------------------

        // #################################################
        // ### Set up to loop over ROI's for this wire   ###
        // #################################################
        const recob::ChannelROI::RegionsOfInterest_f signalROI = channelROI->SignalROIF();

        tbb::parallel_for(
          static_cast<std::size_t>(0),
          signalROI.n_ranges(),
          [&](size_t& rangeIter) {
            const auto& range = signalROI.range(rangeIter);
            // ROI start time
            raw::TDCtick_t roiFirstBinTick = range.begin_index();

            // ###########################################################
            // ### Scan the waveform and find candidate peaks + merge  ###
            // ###########################################################

            reco_tool::ICandidateHitFinder::HitCandidateVec hitCandidateVec;
            reco_tool::ICandidateHitFinder::MergeHitCandidateVec mergedCandidateHitVec;

            fHitFinderToolVec.at(plane)->findHitCandidates(
              range, 0, channel, count, hitCandidateVec);
            fHitFinderToolVec.at(plane)->MergeHitCandidates(
              range, hitCandidateVec, mergedCandidateHitVec);

            // #######################################################
            // ### Lets loop over the pulses we found on this wire ###
            // #######################################################

            for (auto& mergedCands : mergedCandidateHitVec) {
              int startT = mergedCands.front().startTick;
              int endT = mergedCands.back().stopTick;

              // ### Putting in a protection in case things went wrong ###
              // ### In the end, this primarily catches the case where ###
              // ### a fake pulse is at the start of the ROI           ###
              if (endT - startT < 5) continue;

              // #######################################################
              // ### Clearing the parameter vector for the new Pulse ###
              // #######################################################

              // === Setting The Number Of Gaussians to try ===
              int nGausForFit = mergedCands.size();

              // ##################################################
              // ### Calling the function for fitting Gaussians ###
              // ##################################################
              double chi2PerNDF(0.);
              int NDF(1);
              /*stand alone
                reco_tool::IPeakFitter::PeakParamsVec peakParamsVec(nGausForFit);
                */
              reco_tool::IPeakFitter::PeakParamsVec peakParamsVec;

              // #######################################################
              // ### If # requested Gaussians is too large then punt ###
              // #######################################################
              if (mergedCands.size() <= fMaxMultiHit) {
                fPeakFitterTool->findPeakParameters(
                  range.data(), mergedCands, peakParamsVec, chi2PerNDF, NDF);

                // If the chi2 is infinite then there is a real problem so we bail
                if (!(chi2PerNDF < std::numeric_limits<double>::infinity())) {
                  chi2PerNDF = 2. * fChi2NDF;
                  NDF = 2;
                }

                if (fFillHists) fFirstChi2->Fill(chi2PerNDF);
              }

              // #######################################################
              // ### If too large then force alternate solution      ###
              // ### - Make n hits from pulse train where n will     ###
              // ###   depend on the fhicl parameter fLongPulseWidth ###
              // ### Also do this if chi^2 is too large              ###
              // #######################################################
              if (mergedCands.size() > fMaxMultiHit || nGausForFit * chi2PerNDF > fChi2NDF) {
                int longPulseWidth = fLongPulseWidthVec.at(plane);
                int nHitsThisPulse = (endT - startT) / longPulseWidth;

                if (nHitsThisPulse > fLongMaxHitsVec.at(plane)) {
                  nHitsThisPulse = fLongMaxHitsVec.at(plane);
                  longPulseWidth = (endT - startT) / nHitsThisPulse;
                }

                if (nHitsThisPulse * longPulseWidth < endT - startT) nHitsThisPulse++;

                int firstTick = startT;
                int lastTick = std::min(firstTick + longPulseWidth, endT);

                peakParamsVec.clear();
                nGausForFit = nHitsThisPulse;
                NDF = 1.;
                chi2PerNDF = chi2PerNDF > fChi2NDF ? chi2PerNDF : -1.;

                for (int hitIdx = 0; hitIdx < nHitsThisPulse; hitIdx++) {
                  // This hit parameters
                  double ROIsumADC =
                    std::accumulate(range.begin() + firstTick, range.begin() + lastTick, 0.);
                  double peakSigma = (lastTick - firstTick) / 3.;  // Set the width...
                  double peakAmp = 0.3989 * ROIsumADC / peakSigma; // Use gaussian formulation
                  double peakMean = (firstTick + lastTick) / 2.;

                  // Store hit params
                  reco_tool::IPeakFitter::PeakFitParams_t peakParams;

                  peakParams.peakCenter = peakMean;
                  peakParams.peakCenterError = 0.1 * peakMean;
                  peakParams.peakSigma = peakSigma;
                  peakParams.peakSigmaError = 0.1 * peakSigma;
                  peakParams.peakAmplitude = peakAmp;
                  peakParams.peakAmplitudeError = 0.1 * peakAmp;

                  peakParamsVec.push_back(peakParams);

                  // set for next loop
                  firstTick = lastTick;
                  lastTick = std::min(lastTick + longPulseWidth, endT);
                }
              }

              // #######################################################
              // ### Loop through returned peaks and make recob hits ###
              // #######################################################

              int numHits(0);

              // Make a container for what will be the filtered collection
              std::vector<sbn::Hit> filteredHitVec;

              float nextpeak(0);
              float prevpeak(0);
              float nextpeakSig(0);
              float prevpeakSig(0);
              float nsigmaADC(2.0);
              float newright(0);
              float newleft(0);
              for (const auto& peakParams : peakParamsVec) {
                // Extract values for this hit
                float peakAmp = peakParams.peakAmplitude;
                float peakMean = peakParams.peakCenter;
                float peakWidth = peakParams.peakSigma;

                //std::cout<<" ans hits "<<numHits<<" gaus "<<nGausForFit<<std::endl;

                //ANS get prev and next
                if (numHits == 0) {
                  newleft = -9999;
                  newright = 9999;
                  nextpeak = 0;
                  prevpeak = 0;
                  nextpeakSig = 0;
                  prevpeakSig = 0;
                }
                if (numHits < nGausForFit - 1) {
                  nextpeak = (peakParamsVec.at(numHits + 1)).peakCenter;
                  nextpeakSig = (peakParamsVec.at(numHits + 1)).peakSigma;
                  //std::cout<<" ans size "<<peakParamsVec.size()<<" hit "<<numHits<<" next peak "<<nextpeak<<" sig "<<nextpeakSig<<std::endl;
                }
                if (numHits > 0) {
                  prevpeak = (peakParamsVec.at(numHits - 1)).peakCenter;
                  prevpeakSig = (peakParamsVec.at(numHits - 1)).peakSigma;
                  //std::cout<<" ans size "<<peakParamsVec.size()<<"hit "<<numHits<<" prev peak "<<prevpeak<<" sig "<<prevpeakSig<<std::endl;
                }

                // Place one bit of protection here
                if (std::isnan(peakAmp)) {
                  std::cout << "**** hit peak amplitude is a nan! Channel: " << channel
                            << ", start tick: " << startT << std::endl;
                  continue;
                }

                // Extract errors
                float peakAmpErr = peakParams.peakAmplitudeError;
                float peakMeanErr = peakParams.peakCenterError;
                float peakWidthErr = peakParams.peakSigmaError;

                // ### Charge ###
                float charge =
                  chargeFunc(peakMean, peakAmp, peakWidth, fAreaNormsVec[plane], startT, endT);
                ;
                float chargeErr =
                  std::sqrt(TMath::Pi()) * (peakAmpErr * peakWidthErr + peakWidthErr * peakAmpErr);

                // ### limits for getting sums
                std::vector<float>::const_iterator sumStartItr = range.begin() + startT;
                std::vector<float>::const_iterator sumEndItr = range.begin() + endT;

                //### limits for the sum of the Hit based on the gaussian peak and sigma
                std::vector<float>::const_iterator HitsumStartItr =
                  range.begin() + peakMean - nsigmaADC * peakWidth;
                std::vector<float>::const_iterator HitsumEndItr =
                  range.begin() + peakMean + nsigmaADC * peakWidth;

                if (nGausForFit > 1) {
                  if (numHits > 0) {
                    if ((peakMean - nsigmaADC * peakWidth) < (prevpeak + nsigmaADC * prevpeakSig)) {
                      float difPeak = peakMean - prevpeak;
                      float weightpeak = prevpeakSig / (prevpeakSig + peakWidth);
                      HitsumStartItr = range.begin() + prevpeak + difPeak * weightpeak;
                      newleft = prevpeak + difPeak * weightpeak;
                    }
                  }

                  if (numHits < nGausForFit - 1) {
                    if ((peakMean + nsigmaADC * peakWidth) > (nextpeak - nsigmaADC * nextpeakSig)) {
                      float difPeak = nextpeak - peakMean;
                      float weightpeak = peakWidth / (nextpeakSig + peakWidth);
                      HitsumEndItr = range.begin() + peakMean + difPeak * weightpeak;
                      newright = peakMean + difPeak * weightpeak;
                    }
                  }
                }

                //protection to avoid negative ranges
                if (newright - newleft < 0) continue;
                if (HitsumStartItr > HitsumEndItr) continue;

                //avoid ranges out of ROI if it happens
                if (HitsumStartItr < sumStartItr) HitsumStartItr = sumStartItr;

                if (HitsumEndItr > sumEndItr) HitsumEndItr = sumEndItr;

                // ### Sum of ADC counts
                float ROIsumADC = std::accumulate(sumStartItr, sumEndItr, 0.);
                float HitsumADC = std::accumulate(HitsumStartItr, HitsumEndItr, 0.);
                float peakTime  = peakMean + roiFirstBinTick;
                float chi2ByNDF = chi2PerNDF;  // convert to float

                // ok, now create the hit
                sbn::HitCreator hitcreator(
                  *channelROI,                // channelROI reference
                  wid,                        // wire ID
                  startT + roiFirstBinTick,   // start_tick TODO check
                  endT + roiFirstBinTick,     // end_tick TODO check
                  peakWidth,                  // rms
                  peakTime,                   // peak_time
                  peakMeanErr,                // sigma_peak_time
                  peakAmp,                    // peak_amplitude
                  peakAmpErr,                 // sigma_peak_amplitude
                  charge,                     // hit_integral
                  chargeErr,                  // hit_sigma_integral
                  ROIsumADC,                  // summedADC of the ROI
                  HitsumADC,                  // summedADC of the Hit
                  nGausForFit,                // multiplicity
                  numHits,                    // local_index TODO check that the order is correct
                  chi2ByNDF,                  // goodness_of_fit
                  NDF                         // dof
                );

                if (filteredHitCol) filteredHitVec.push_back(hitcreator.copy());

                const sbn::Hit hit(hitcreator.move());

                // This loop will store ALL hits
                hitstruct tmp{std::move(hit), channelROI};
                hitstruct_vec.push_back(std::move(tmp));

                numHits++;
              } // <---End loop over gaussians

              // Should we filter hits?
              if (filteredHitCol && !filteredHitVec.empty()) {
                // #######################################################################
                // Is all this sorting really necessary?  Would it be faster to just loop
                // through the hits and perform simple cuts on amplitude and width on a
                // hit-by-hit basis, either here in the module (using fPulseHeightCuts and
                // fPulseWidthCuts) or in HitFilterAlg?
                // #######################################################################

                // Sort in ascending peak height
                std::sort(filteredHitVec.begin(),
                          filteredHitVec.end(),
                          [](const auto& left, const auto& right) {
                            return left.PeakAmplitude() > right.PeakAmplitude();
                          });

                // Reject if the first hit fails the PH/wid cuts
                if (filteredHitVec.front().PeakAmplitude() < fPulseHeightCuts.at(plane) ||
                    filteredHitVec.front().RMS() < fPulseWidthCuts.at(plane))
                  filteredHitVec.clear();

                // Now check other hits in the snippet
                if (filteredHitVec.size() > 1) {
                  // The largest pulse height will now be at the front...
                  float largestPH = filteredHitVec.front().PeakAmplitude();

                  // Find where the pulse heights drop below threshold
                  float threshold(fPulseRatioCuts.at(plane));

                  std::vector<sbn::Hit>::iterator smallHitItr =
                    std::find_if(filteredHitVec.begin(),
                                 filteredHitVec.end(),
                                 [largestPH, threshold](const auto& hit) {
                                   return hit.PeakAmplitude() < 8. &&
                                          hit.PeakAmplitude() / largestPH < threshold;
                                 });

                  // Shrink to fit
                  if (smallHitItr != filteredHitVec.end())
                    filteredHitVec.resize(std::distance(filteredHitVec.begin(), smallHitItr));

                  // Resort in time order
                  std::sort(filteredHitVec.begin(),
                            filteredHitVec.end(),
                            [](const auto& left, const auto& right) {
                              return left.PeakTime() < right.PeakTime();
                            });
                }

                // Copy the hits we want to keep to the filtered hit collection
//                for (const auto& filteredHit : filteredHitVec)
//                  if (!fHitFilterAlg || fHitFilterAlg->IsGoodHit(filteredHit)) {
//                    hitstruct tmp{std::move(filteredHit), channelROI};
//                    filthitstruct_vec.push_back(std::move(tmp));
//                  }
//
                if (fFillHists) fChi2->Fill(chi2PerNDF);
              }
            } //<---End loop over merged candidate hits
          }   //<---End looping over ROI's
        );    //end tbb parallel for
      }       //<---End looping over all the wires
    );        //end tbb parallel for

    for (size_t i = 0; i < hitstruct_vec.size(); i++) {
      allHitCol.emplace_back(hitstruct_vec[i].hit_tbb, hitstruct_vec[i].channelROI_tbb);
    }

    for (size_t j = 0; j < filthitstruct_vec.size(); j++) {
      filteredHitCol->emplace_back(filthitstruct_vec[j].hit_tbb, filthitstruct_vec[j].channelROI_tbb);
    }

    //==================================================================================================
    // End of the event -- move the hit collection and the associations into the event

    if (filteredHitCol) {

      // If we filtered hits but no instance name was
      // specified for the "all hits" collection, then
      // only save the filtered hits to the event
      if (fAllHitsInstanceName == "") {
        filteredHitCol->put_into(evt);

        // otherwise, save both
      }
      else {
        filteredHitCol->put_into(evt);
        allHitCol.put_into(evt);
      }
    }
    else {
      allHitCol.put_into(evt);
    }
  } // End of produce()

  DEFINE_ART_MODULE(GaussHitFinderSBN)

} // end of hit namespace
