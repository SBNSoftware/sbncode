//////////////////////////////////////////////////////////////////
// \file     CAFMaker_module.cc
/// \brief   This module creates Common Analysis Files.
//           Inspired by the NOvA CAFMaker package
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////

// ---------------- TO DO ----------------
//
// - Add in cycle and batch to params
// - Move this list some place useful
// - Add reco.CRT branch
// ---------------------------------------


#include "sbncode/CAFMaker/CAFMakerParams.h"
#include "sbncode/CAFMaker/FillFlashMatch.h"
#include "sbncode/CAFMaker/FillTrue.h"
#include "sbncode/CAFMaker/FillReco.h"
#include "sbncode/CAFMaker/FillExposure.h"
#include "sbncode/CAFMaker/FillTrigger.h"
#include "sbncode/CAFMaker/Utils.h"

// C/C++ includes
#include <fenv.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <array>

#ifdef DARWINBUILD
#include <libgen.h>
#endif

#include "ifdh_art/IFDHService/IFDH_service.h"

// ROOT includes
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TRandomGen.h"
#include "TObjString.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "cetlib_except/exception.h"
#include "cetlib_except/demangle.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "fhiclcpp/ParameterSetRegistry.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlTruth.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightParameterSet.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"
#include "sbnobj/Common/POTAccounting/NuMISpillInfo.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"


#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "larcoreobj/SummaryData/POTSummary.h"

// StandardRecord
#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"

#include "sbnanaobj/StandardRecord/Flat/FlatRecord.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/TriggerData.h"

// // CAFMaker
#include "sbncode/CAFMaker/AssociationUtil.h"
// #include "sbncode/CAFMaker/Blinding.h"

// Metadata
#include "sbncode/Metadata/MetadataSBN.h"

namespace sbn{
  namespace evwgh{
    std::ostream& operator<<(std::ostream& os, const sbn::evwgh::EventWeightParameterSet& p)
    {
      // TODO proper implementation of this should be added in sbnobj
      os << p.fName << " " << p.fRWType << std::endl;
      for(const auto& it: p.fParameterMap){
        os << it.first.fName << " " << it.first.fMean << " " << it.first.fWidth << std::endl << " ";
        for(float v: it.second) os << " " << v;
        os << std::endl;
      }
      return os;
    }
  }
}

namespace caf {

/// Function to calculate a timestamp from the spill info product
template <typename SpillInfo>
double spillInfoToTimestamp(SpillInfo const& info) {
  return static_cast<double>(info.spill_time_s) +
         static_cast<double>(info.spill_time_ns)*1.0e-9;
}

/// Module to create Common Analysis Files from ART files
class CAFMaker : public art::EDProducer {
 public:
  // Allows 'nova --print-description' to work
  using Parameters = art::EDProducer::Table<CAFMakerParams>;

  explicit CAFMaker(const Parameters& params);
  virtual ~CAFMaker();

  void produce(art::Event& evt) noexcept;

  void respondToOpenInputFile(const art::FileBlock& fb);

  void beginJob();
  void endJob();
  virtual void beginRun(art::Run& r);
  virtual void beginSubRun(art::SubRun& sr);
  virtual void endSubRun(art::SubRun& sr);

 protected:
  CAFMakerParams fParams;

  std::string fCafFilename;
  std::string fCafBlindFilename;
  std::string fCafPrescaleFilename;

  std::string fFlatCafFilename;
  std::string fFlatCafBlindFilename;
  std::string fFlatCafPrescaleFilename;
  
  std::string fSourceFile;

  bool fFirstInSubRun;
  unsigned int fIndexInFile = SRHeader::NoSourceIndex();
  bool fFirstBlindInFile;
  bool fFirstPrescaleInFile;
  int fFileNumber;
  double fTotalPOT;
  double fSubRunPOT;
  unsigned int fOffbeamBNBGates;
  unsigned int fOffbeamNuMIGates;
  double fTotalSinglePOT;
  double fTotalEvents;
  double fBlindEvents;
  double fPrescaleEvents;
  std::vector<caf::SRBNBInfo> fBNBInfo; ///< Store detailed BNB info to save into the first StandardRecord of the output file
  std::vector<caf::SRNuMIInfo> fNuMIInfo; ///< Store detailed NuMI info to save into the first StandardRecord of the output file
  std::map<unsigned int,sbn::BNBSpillInfo> fBNBInfoEventMap; ///< Store detailed BNB info to save for the particular spills of events
  std::map<unsigned int,sbn::NuMISpillInfo> fNuMIInfoEventMap; ///< Store detailed NuMI info to save for the particular spills of events
  bool fHasBNBInfo;
  bool fHasNuMIInfo;

  // int fCycle;
  // int fBatch;

  TFile* fFile = 0;
  TFile* fFileb = 0;
  TFile* fFilep = 0;

  TTree* fRecTree = 0;
  TTree* fRecTreeb = 0;
  TTree* fRecTreep = 0;

  TFile* fFlatFile = 0;
  TFile* fFlatFileb = 0;
  TFile* fFlatFilep = 0;

  TTree* fFlatTree = 0;
  TTree* fFlatTreeb = 0;
  TTree* fFlatTreep = 0;

  flat::Flat<caf::StandardRecord>* fFlatRecord = 0;
  flat::Flat<caf::StandardRecord>* fFlatRecordb = 0;
  flat::Flat<caf::StandardRecord>* fFlatRecordp = 0;

  Det_t fDet;  ///< Detector ID in caf namespace typedef

  // volumes
  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;

  // random number generator for fake reco
  TRandom *fFakeRecoTRandom;

  // random number generator for prescaling
  TRandom *fBlindTRandom;

  /// What position in the vector each parameter set take
  std::map<std::string, unsigned int> fWeightPSetIndex;
  /// Map from parameter labels to previously seen parameter set configuration
  std::map<std::string, std::vector<sbn::evwgh::EventWeightParameterSet>> fPrevWeightPSet;

  std::string DeriveFilename(const std::string& inname,
                             const std::string& ext) const;

  static std::string Basename(const std::string& path);

  void AddEnvToFile(TFile* f);
  void AddMetadataToFile(TFile* f,
                         const std::map<std::string, std::string>& metadata);
  void AddGlobalTreeToFile(TFile* outfile, caf::SRGlobal& global) const;
  void AddHistogramsToFile(TFile* outfile,bool isBlindPOT, bool isPrescalePOT) const;

  void InitializeOutfiles();

  void BlindEnergyParameters(StandardRecord* brec);
  double GetBlindPOTScale() const;

  void InitVolumes(); ///< Initialize volumes from Gemotry service

  void FixPMTReferenceTimes(StandardRecord &rec, double PMT_reference_time);
  void FixCRTReferenceTimes(StandardRecord &rec, double CRTT0_reference_time, double CRTT1_reference_time);

  /// Equivalent of FindManyP except a return that is !isValid() prints a
  /// messsage and aborts if StrictMode is true.
  template <class T, class U>
  art::FindManyP<T> FindManyPStrict(const U& from, const art::Event& evt,
                                    const art::InputTag& label) const;

  template <class T, class D, class U>
  art::FindManyP<T, D> FindManyPDStrict(const U& from,
                                        const art::Event& evt,
                                        const art::InputTag& tag) const;

  /// Equivalent of FindOneP except a return that is !isValid() prints a
  /// messsage and aborts if StrictMode is true.
  template <class T, class U>
  art::FindOneP<T> FindOnePStrict(const U& from, const art::Event& evt,
				  const art::InputTag& label) const;


  /// \brief Retrieve an object from an association, with error handling
  ///
  /// This can go wrong in two ways: either the FindManyP itself is
  /// invalid, or the result for the requested index is empty. In most
  /// cases these have the same response, so conflating them here
  /// saves redundancy elsewhere.
  ///
  /// \param      fm  The FindManyP object describing the association
  /// \param      idx Which element of the FindManyP to look it
  /// \param[out] ret The product retrieved
  /// \return          Whether \a ret was filled
  template <class T>
  bool GetAssociatedProduct(const art::FindManyP<T>& fm, int idx, T& ret) const;

  /// Equivalent of evt.getByLabel(label, handle) except failedToGet
  /// prints a message and aborts if StrictMode is true.
  template <class EvtT, class T>
  void GetByLabelStrict(const EvtT& evt, const std::string& label,
                        art::Handle<T>& handle) const;

  /// Equivalent of evt.getByLabel(label, handle) except failedToGet
  /// prints a message.
  template <class T>
  void GetByLabelIfExists(const art::Event& evt, const std::string& label,
                          art::Handle<T>& handle) const;

  /// \param      pset The parameter set
  /// \param      name Pass "foo.bar.baz" as {"foo", "bar", "baz"}
  /// \param[out] ret  Value of the key, not set if we return false
  /// \return          Whether the key was found
  template <class T>
  bool GetPsetParameter(const fhicl::ParameterSet& pset,
                        const std::vector<std::string>& name, T& ret) const;

  static bool EssentiallyEqual(double a, double b, double precision = 0.0001) {
    return a <= (b + precision) && a >= (b - precision);
  }

  static bool sortRBTrkLength(const art::Ptr<recob::Track>& a,
                                const art::Ptr<recob::Track>& b) {
    return a->Length() > b->Length();
  }
  // static bool sortTrackLength(const SRTrack& a, const SRTrack& b) {
  //   return a.len > b.len;
  // }
//.......................................................................
}; //Producer

//.......................................................................

  CAFMaker::CAFMaker(const Parameters& params)
  : art::EDProducer{params},
    fParams(params()), fFile(0)
  {
  // Note: we will define isRealData on a per event basis in produce function [using event.isRealData()], at least for now.

  fCafFilename = fParams.CAFFilename();
  fFlatCafFilename = fParams.FlatCAFFilename();

  // Normally CAFMaker is run wit no output ART stream, so these go
  // nowhere, but can be occasionally useful for filtering in ART

  produces<std::vector<caf::StandardRecord>>();
  //produces<art::Assns<caf::StandardRecord, recob::Slice>>();

  // setup volume definitions
  InitVolumes();

  // setup random number generators
  fFakeRecoTRandom = new TRandomMT64(art::ServiceHandle<rndm::NuRandomService>()->getSeed());
  if (fParams.CreateBlindedCAF()) {
    fBlindTRandom = new TRandomMT64(art::ServiceHandle<rndm::NuRandomService>()->getSeed());
  }
}

//......................................................................
  double CAFMaker::GetBlindPOTScale() const {
   std::string bstring = std::to_string(fParams.POTBlindSeed());
   int slen = bstring.length();
   std::string s1 = bstring.substr(0,int(slen/2));
   std::string s2 = bstring.substr(int(slen/2));
   double rat = stod(s1)/stod(s2);
   while (abs(rat)>1){
     rat = -1 * (abs(rat) - 1);
   }
   return 1 + rat*0.3;
   
  }
//......................................................................
void CAFMaker::BlindEnergyParameters(StandardRecord* brec) {

  //simple cuts for trk and shower variables
  //blind events with a potential lepton with momentum > 0.6 that starts in fiducial volume
  for (caf::SRPFP& pfp: brec->reco.pfp) {
    const caf::SRVector3D start = pfp.trk.start;
    if ( ((start.x < -71.1 - 25 && start.x > -369.33 + 25 ) ||
	  (start.x > 71.1 + 25 && start.x < 369.33 - 25 )) &&
	 (start.y > -181.7 + 25 && start.y < 134.8 - 25 ) &&
	 (start.z  > -895.95 + 30 && start.z < 895.95 - 50)) {

      if (pfp.trk.mcsP.fwdP_muon > 0.6) {
	pfp.trk.mcsP.fwdP_muon = TMath::QuietNaN();    
      }
      if (pfp.trk.rangeP.p_muon > 0.6) {
	pfp.trk.rangeP.p_muon = TMath::QuietNaN();
      }
    }
  }

  //Note shower energy may not be currently very functional
  for (caf::SRPFP& pfp: brec->reco.pfp) {
    const caf::SRVector3D start = pfp.shw.start;
    if ( ((start.x < -71.1 - 25 && start.x > -369.33 + 25 ) ||
	  (start.x > 71.1 + 25 && start.x < 369.33 - 25 )) &&
	 (start.y > -181.7 + 25 && start.y < 134.8 - 25 ) &&
	 (start.z  > -895.95 + 30 && start.z < 895.95 - 50)) {
      if (pfp.shw.bestplane_energy > 0.6) {
	pfp.shw.bestplane_energy = TMath::QuietNaN();
	pfp.shw.plane[0].energy = TMath::QuietNaN();
	pfp.shw.plane[1].energy = TMath::QuietNaN();
	pfp.shw.plane[2].energy = TMath::QuietNaN();
      }
    }
  }

  // And for slices, check vertex in FV and then check tracks and showers
  for (caf::SRSlice& slc: brec->slc) {
    const caf::SRVector3D vtx = slc.vertex;
    if ( ((vtx.x < -71.1 - 25 && vtx.x > -369.33 + 25 ) ||
	  (vtx.x > 71.1 + 25 && vtx.x < 369.33 - 25 )) &&
	 (vtx.y > -181.7 + 25 && vtx.y < 134.8 - 25 ) &&
	 (vtx.z  > -895.95 + 30 && vtx.z < 895.95 - 50)) {

      for (caf::SRPFP& pfp: slc.reco.pfp) {
	if (pfp.trk.mcsP.fwdP_muon > 0.6) {
	  pfp.trk.mcsP.fwdP_muon = TMath::QuietNaN();    
	}
	if (pfp.trk.rangeP.p_muon > 0.6) {
	  pfp.trk.rangeP.p_muon = TMath::QuietNaN();
	}
      }
      for (caf::SRPFP& pfp: slc.reco.pfp) {
	if (pfp.shw.bestplane_energy > 0.6) {
	  pfp.shw.bestplane_energy = TMath::QuietNaN();
	  pfp.shw.plane[0].energy = TMath::QuietNaN();
	  pfp.shw.plane[1].energy = TMath::QuietNaN();
	  pfp.shw.plane[2].energy = TMath::QuietNaN();
	}
      }
    }
  }
}

void CAFMaker::FixPMTReferenceTimes(StandardRecord &rec, double PMT_reference_time) {
  // Fix the flashes
  for (SROpFlash &f: rec.opflashes) {
    f.time += PMT_reference_time;
    f.timemean += PMT_reference_time;
    f.firsttime += PMT_reference_time;
  }

  // Fix the flash matches
  for (SRSlice &s: rec.slc) {
    s.fmatch.time += PMT_reference_time;
    s.fmatch_a.time += PMT_reference_time;
    s.fmatch_b.time += PMT_reference_time;
  }

  // TODO: fix more?

}

void CAFMaker::FixCRTReferenceTimes(StandardRecord &rec, double CRTT0_reference_time, double CRTT1_reference_time) {
  // Fix the hits

  double crttime_to_shift = fParams.CRTUseTS0() ? CRTT0_reference_time : CRTT1_reference_time;

  // As discussed/described in https://github.com/SBNSoftware/sbncode/pull/251,
  // we added CRTHit::t0 and CRTHit::t1 in addition to CRTHit::time,
  // not to break any existing studies that still use "CRTHit::time"
  for (SRCRTHit &h: rec.crt_hits) {
    h.t0 += CRTT0_reference_time;
    h.t1 += CRTT1_reference_time;
    h.time += crttime_to_shift;
  }

  // Fix the hit matches
  for (SRSlice &s: rec.slc) {
    for (SRPFP &pfp: s.reco.pfp) {
      pfp.trk.crthit.hit.t0 += CRTT0_reference_time; 
      pfp.trk.crthit.hit.t1 += CRTT1_reference_time;
      pfp.trk.crthit.hit.time += crttime_to_shift;
    }
  }
  for (SRPFP &pfp: rec.reco.pfp) {
    pfp.trk.crthit.hit.t0 += CRTT0_reference_time;
    pfp.trk.crthit.hit.t1 += CRTT1_reference_time;
    pfp.trk.crthit.hit.time += crttime_to_shift;
  }

  // TODO: fix more?
  // Tracks?

}

void CAFMaker::InitVolumes() {
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // first the TPC volumes
  for (auto const &cryo: geometry->Iterate<geo::CryostatGeo>()) {
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryo.ID())) {
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
    }
     fTPCVolumes.push_back(std::move(this_tpc_volumes));
  }

  // then combine them into active volumes
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: fTPCVolumes) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
  }
}

//......................................................................
CAFMaker::~CAFMaker()
{
  delete fRecTree;
  delete fFile;

  delete fFlatRecord;
  delete fFlatTree;
  delete fFlatFile;

  if (fParams.CreateBlindedCAF()) {
    delete fRecTreeb;
    delete fRecTreep;
    delete fFileb;
    delete fFilep;
    delete fFlatRecordb;
    delete fFlatRecordp;
    delete fFlatTreeb;
    delete fFlatTreep;
    delete fFlatFileb;
    delete fFlatFilep;
    delete fBlindTRandom;
  }

  delete fFakeRecoTRandom;

}

//......................................................................
std::string CAFMaker::Basename(const std::string& path)
{
  // C++17: use filesystem library (Clang 7 still not compliant)
  constexpr char sep = '/';
  std::size_t const iSep = path.rfind(sep);
  return (iSep == std::string::npos)? path: path.substr(iSep + 1);
}

//......................................................................
std::string CAFMaker::DeriveFilename(const std::string& inname,
                                     const std::string& ext) const
{
  char* temp = new char[inname.size()+1];
  std::strcpy(temp, inname.c_str());
  std::string ret = basename(temp);
  delete[] temp;
  const size_t dotpos = ret.rfind('.'); // Find last dot
  assert(dotpos != std::string::npos);  // Must have a dot, surely?
  ret.resize(dotpos); // Truncate everything after dot
  ret += ext;
  return ret;
}

//......................................................................
void CAFMaker::respondToOpenInputFile(const art::FileBlock& fb) {
  
  std::string const inputBasename = Basename(fb.fileName()); // includes suffix
  
  if ((fParams.CreateCAF() && !fFile) ||
      (fParams.CreateFlatCAF() && !fFlatFile) ||
      (fParams.CreateBlindedCAF() && (!fFileb || !fFilep))) {
    // If Filename wasn't set in the FCL, and this is the
    // first file we've seen
    if(fParams.CreateCAF() && fCafFilename.empty()){
      if (fParams.CreateBlindedCAF()){
	fCafFilename = DeriveFilename(fb.fileName(), fParams.UnblindFileExtension());
	fCafFilename = DeriveFilename(fCafFilename, fParams.FileExtension());
	fCafBlindFilename = DeriveFilename(fb.fileName(), fParams.BlindFileExtension());
	fCafBlindFilename = DeriveFilename(fCafBlindFilename, fParams.FileExtension());
	fCafPrescaleFilename = DeriveFilename(fb.fileName(), fParams.PrescaleFileExtension());
	fCafPrescaleFilename = DeriveFilename(fCafPrescaleFilename, fParams.FileExtension());
      }
      else {
	fCafFilename = DeriveFilename(fb.fileName(), fParams.FileExtension());
      }
    }
    if(fParams.CreateFlatCAF() && fFlatCafFilename.empty()){
      if (fParams.CreateBlindedCAF()){
	fFlatCafFilename = DeriveFilename(fb.fileName(), fParams.UnblindFileExtension());
	fFlatCafFilename = DeriveFilename(fFlatCafFilename, fParams.FlatCAFFileExtension());
	fFlatCafBlindFilename = DeriveFilename(fb.fileName(), fParams.BlindFileExtension());
	fFlatCafBlindFilename = DeriveFilename(fFlatCafBlindFilename, fParams.FlatCAFFileExtension());
	fFlatCafPrescaleFilename = DeriveFilename(fb.fileName(), fParams.PrescaleFileExtension());
	fFlatCafPrescaleFilename = DeriveFilename(fFlatCafPrescaleFilename, fParams.FlatCAFFileExtension());
      }
      else {
	fFlatCafFilename = DeriveFilename(fb.fileName(), fParams.FlatCAFFileExtension());
      }
    }
    if (fParams.CreateBlindedCAF() && fCafBlindFilename.empty()) {
      const std::string basename = fCafFilename;
      fCafBlindFilename = DeriveFilename(basename, fParams.BlindFileExtension());
      fCafBlindFilename = DeriveFilename(fCafBlindFilename, fParams.FileExtension());
      fCafPrescaleFilename = DeriveFilename(basename, fParams.PrescaleFileExtension());
      fCafPrescaleFilename = DeriveFilename(fCafPrescaleFilename, fParams.FileExtension());
    }
    if (fParams.CreateBlindedCAF() && fParams.CreateFlatCAF() && fFlatCafBlindFilename.empty()) {
      const std::string basename = fFlatCafFilename;
      fFlatCafBlindFilename = DeriveFilename(basename, fParams.BlindFileExtension());
      fFlatCafBlindFilename = DeriveFilename(fFlatCafBlindFilename, fParams.FlatCAFFileExtension());
      fFlatCafPrescaleFilename = DeriveFilename(basename, fParams.PrescaleFileExtension());
      fFlatCafPrescaleFilename = DeriveFilename(fFlatCafPrescaleFilename, fParams.FlatCAFFileExtension());
    }
    InitializeOutfiles();
  }

  fFileNumber ++;
  fIndexInFile = 0;
  fFirstBlindInFile = true;
  fFirstPrescaleInFile = true;
  fSourceFile = inputBasename;

}

//......................................................................
void CAFMaker::beginJob()
{
}

//......................................................................
void CAFMaker::AddGlobalTreeToFile(TFile* outfile, caf::SRGlobal& global) const
{
  outfile->cd();

  TTree* globalTree = new TTree("globalTree", "globalTree");
  SRGlobal* pglobal = &global;
  TBranch* br = globalTree->Branch("global", "caf::SRGlobal", &pglobal);
  if(!br) abort();
  globalTree->Fill();
  globalTree->Write();
}

//......................................................................
void CAFMaker::beginRun(art::Run& run) {
  fDet = kUNKNOWN;

  caf::Det_t override = kUNKNOWN;
  if(fParams.DetectorOverride() == "sbnd") override = kSBND;
  if(fParams.DetectorOverride() == "icarus") override = kICARUS;
  if(!fParams.DetectorOverride().empty() && override == kUNKNOWN){
    std::cout << "CAFMaker: unrecognized value for DetectorOverride parameter: '" << fParams.DetectorOverride() << "'" << std::endl;
    abort();
  }

  // Heuristic method to determine the detector ID
  const geo::GeometryCore* geom = lar::providerFrom<geo::Geometry>();

  std::string gdml = geom->GDMLFile();
  gdml = basename(gdml.c_str()); // strip directory part
  std::cout << "CAFMaker: Attempting to deduce detector from GDML file name: '" << gdml
            << "' and configured detector name: '" << geom->DetectorName() << "'. ";
  // Lowercase filename, in case it contains "SBND" or "Icarus" etc
  for(unsigned int i = 0; i < gdml.size(); ++i) gdml[i] = std::tolower(gdml[i]);
  // Do we find the string in either of the names?
  const bool hasSBND = ((gdml.find("sbnd") != std::string::npos) ||
                        (geom->DetectorName().find("sbnd") != std::string::npos));
  const bool hasIcarus = ((gdml.find("icarus") != std::string::npos) ||
                          (geom->DetectorName().find("icarus") != std::string::npos));

  // Either no evidence, or ambiguous evidence
  if(hasSBND == hasIcarus){
    std::cout << "Unable to automatically determine detector!" << std::endl;
    if(override == kUNKNOWN) abort();
  }
  // Now must be one or the other
  if(hasSBND){
    fDet = kSBND;
    std::cout << "Detected SBND" << std::endl;
  }
  if(hasIcarus){
    fDet = kICARUS;
    std::cout << "Detected Icarus" << std::endl;
  }

  if(override != kUNKNOWN){
    std::cout << "Detector set to ";
    std::cout << ((override == kSBND) ? "SBND" : "Icarus");
    std::cout << " based on user configuration." << std::endl;
    if(fDet == override){
      std::cout << "  This was redundant with the auto-detection. Suggest to not specify DetectorOverride" << std::endl;
    }
    else if(fDet != kUNKNOWN){
      std::cout << "  This OVERRODE the auto-detection. Are you sure this is what you wanted?" << std::endl;
    }
    fDet = override;
  }


  if(fParams.SystWeightLabels().empty()) return; // no need for globalTree

  SRGlobal global;

  for(const std::string& label: fParams.SystWeightLabels()){
    art::Handle<std::vector<sbn::evwgh::EventWeightParameterSet>> wgt_params;
    GetByLabelStrict(run, label, wgt_params);

    if(fPrevWeightPSet.count(label)){
      if(fPrevWeightPSet[label] != *wgt_params){
        std::cout << "CAFMaker: Run-level EventWeightParameterSet mismatch."
                  << std::endl;
        std::cout << "Previous parameter sets:";
        for(const sbn::evwgh::EventWeightParameterSet& p: fPrevWeightPSet[label]){
          std::cout << p << std::endl;
        }
        std::cout << "\nNew parameter sets:";
        for(const sbn::evwgh::EventWeightParameterSet& p: *wgt_params){
          std::cout << p << std::endl;
        }
        abort();
      }
      return; // Match, no need to refill into tree
    }

    // If there were no weights available, return
    if (!wgt_params.isValid()){
      std::cout << "CAFMaker: no EventWeightParameterSet found under label '" << label << "'" << std::endl;
      return;
    }

    fPrevWeightPSet[label] = *wgt_params;

    for(const sbn::evwgh::EventWeightParameterSet& pset: *wgt_params){
      FillSRGlobal(pset, global, fWeightPSetIndex);
    } // end for pset
  } // end for label

  if(fFile) AddGlobalTreeToFile(fFile, global);
  if(fParams.CreateBlindedCAF() && fFileb) AddGlobalTreeToFile(fFileb, global);
  if(fParams.CreateBlindedCAF() && fFilep) AddGlobalTreeToFile(fFilep, global);
  if(fFlatFile) AddGlobalTreeToFile(fFlatFile, global);
  if(fParams.CreateBlindedCAF() && fFlatFileb) AddGlobalTreeToFile(fFlatFileb, global);
  if(fParams.CreateBlindedCAF() && fFlatFilep) AddGlobalTreeToFile(fFlatFilep, global);
}

//......................................................................
void CAFMaker::beginSubRun(art::SubRun& sr) {

  // get POT information
  fBNBInfo.clear();
  fNuMIInfo.clear();

  fBNBInfoEventMap.clear();
  fNuMIInfoEventMap.clear();
  fHasBNBInfo = false;
  fHasNuMIInfo = false;

  fSubRunPOT = 0;
  fOffbeamBNBGates = 0;
  fOffbeamNuMIGates = 0;

  auto bnb_spill          = sr.getHandle<std::vector<sbn::BNBSpillInfo>>(fParams.BNBPOTDataLabel());
  auto numi_spill         = sr.getHandle<std::vector<sbn::NuMISpillInfo>>(fParams.NuMIPOTDataLabel());
  auto bnb_offbeam_spill  = sr.getHandle<std::vector<sbn::EXTCountInfo>>(fParams.OffbeamBNBCountDataLabel());
  auto numi_offbeam_spill = sr.getHandle<std::vector<sbn::EXTCountInfo>>(fParams.OffbeamNuMICountDataLabel());

  if(bool(bnb_spill) + bool(numi_spill) + bool(bnb_offbeam_spill) + bool(numi_offbeam_spill) > 1) {
    std::cout << "Expected at most one of " << fParams.BNBPOTDataLabel() << ", "
              << fParams.NuMIPOTDataLabel() << ", " << fParams.OffbeamBNBCountDataLabel() << ", and "
              << fParams.OffbeamNuMICountDataLabel() << ". Found ";
    if(bnb_spill) std::cout << fParams.BNBPOTDataLabel() << " ";
    if(numi_spill) std::cout << fParams.NuMIPOTDataLabel() << " ";
    if(bnb_offbeam_spill) std::cout << fParams.OffbeamBNBCountDataLabel() << " ";
    if(numi_offbeam_spill) std::cout << fParams.OffbeamNuMICountDataLabel();
    std::cout << std::endl;
    abort();
  }

  if(bnb_spill){
    FillExposure(*bnb_spill, fBNBInfo, fSubRunPOT);
    fTotalPOT += fSubRunPOT;

    // Find the spill for each event and fill the event map:
    // We take the latest spill for a given event number to be the one to keep
    fHasBNBInfo = true;
    for(const sbn::BNBSpillInfo& info: *bnb_spill)
    {
      auto& storedInfo = fBNBInfoEventMap[info.event]; // creates if needed
      if ( (storedInfo.event == UINT_MAX) || spillInfoToTimestamp(info) > spillInfoToTimestamp(storedInfo) ) {
        storedInfo = std::move(info);
      }
    }
  }
  else if (numi_spill) {
    FillExposureNuMI(*numi_spill, fNuMIInfo, fSubRunPOT);
    fTotalPOT += fSubRunPOT;

    // Find the spill for each event and fill the event map:
    // We take the latest spill for a given event number to be the one to keep
    fHasNuMIInfo = true;
    for(const sbn::NuMISpillInfo& info: *numi_spill)
    {
      auto& storedInfo = fNuMIInfoEventMap[info.event]; // creates if needed
      if ( (storedInfo.event == UINT_MAX) || spillInfoToTimestamp(info) > spillInfoToTimestamp(storedInfo) ) {
        storedInfo = std::move(info);
      }
    }
  }
  else if (bnb_offbeam_spill){
    for(const auto& spill: *bnb_offbeam_spill) {
      fOffbeamBNBGates += spill.gates_since_last_trigger;
    }
  }
  else if (numi_offbeam_spill){
    for(const auto& spill: *numi_offbeam_spill) {
      fOffbeamNuMIGates += spill.gates_since_last_trigger;
    }
  }
  else if(auto pot_handle = sr.getHandle<sumdata::POTSummary>(fParams.GenLabel())){
    fSubRunPOT = pot_handle->totgoodpot;
    fTotalPOT += fSubRunPOT;
  }
  else{
    if(!fParams.BNBPOTDataLabel().empty() || !fParams.GenLabel().empty() || !fParams.NuMIPOTDataLabel().empty() ||
       !fParams.OffbeamBNBCountDataLabel().empty() || !fParams.OffbeamNuMICountDataLabel().empty()){
      std::cout << "Found neither BNB data POT info under '"
                << fParams.BNBPOTDataLabel()
                << "' nor NuMIdata POT info under '"
                << fParams.NuMIPOTDataLabel()
                << "' nor BNB EXT Count info under '"
                << fParams.OffbeamBNBCountDataLabel()
                << "' nor NuMI EXT Count info under '"
                << fParams.OffbeamNuMICountDataLabel()
                << "' nor MC POT info under '"
                << fParams.GenLabel() << "'"
                << std::endl;
      if(fParams.StrictMode()) abort();
    }

    // Otherwise, if one label is blank, maybe no POT was the expected result
  }

  std::cout << "POT: " << fSubRunPOT << std::endl;

  fFirstInSubRun = true;
}

//......................................................................
void CAFMaker::AddEnvToFile(TFile* outfile)
{
  // Global information about the processing details:
  std::map<std::string, std::string> envmap;

  // Environ comes from unistd.h
  // environ is not present on OSX for some reason, so just use getenv to
  // grab the variables we care about.
#ifdef DARWINBUILD
  std::set<TString> variables;
  variables.insert("USER");
  variables.insert("HOSTNAME");
  variables.insert("PWD");
  for (auto var : variables) if(getenv(var)) envmap[var] = getenv(var);
#else
  for (char** penv = environ; *penv; ++penv) {
    const std::string pair = *penv;
    const size_t split = pair.find("=");
    if(split == std::string::npos) continue;  // Huh?
    const std::string key = pair.substr(0, split);
    const std::string value = pair.substr(split + 1);
    envmap[key] = value;
  }
#endif

  // Default constructor is "now"
  envmap["date"] = TTimeStamp().AsString();
  envmap["output"] = fCafFilename;

  // Get the command-line we were invoked with. What I'd really like is
  // just the fcl script and list of input filenames in a less hacky
  // fashion. I'm not sure that's actually possible in ART.
  // TODO: ask the artists.
  FILE* cmdline = fopen("/proc/self/cmdline", "rb");
  char* arg = 0;
  size_t size = 0;
  std::string cmd;
  while (getdelim(&arg, &size, 0, cmdline) != -1) {
    cmd += arg;
    cmd += " ";
  }
  free(arg);
  fclose(cmdline);

  envmap["cmd"] = cmd;

  outfile->mkdir("env")->cd();

  TTree* trenv = new TTree("envtree", "envtree");
  std::string key, value;
  trenv->Branch("key", &key);
  trenv->Branch("value", &value);
  for(const auto& keyval: envmap){
    key = keyval.first;
    value = keyval.second;
    trenv->Fill();
  }
  trenv->Write();
}

//......................................................................
void CAFMaker::AddMetadataToFile(TFile* outfile, const std::map<std::string, std::string>& metadata)
{
  outfile->mkdir("metadata")->cd();

  TTree* trmeta = new TTree("metatree", "metatree");
  std::string key, value;
  trmeta->Branch("key", &key);
  trmeta->Branch("value", &value);
  for(const auto& keyval: metadata){
    key = keyval.first;
    value = keyval.second;
    trmeta->Fill();
  }
  trmeta->Write();
}

//......................................................................
void CAFMaker::InitializeOutfiles()
{
  if(fParams.CreateCAF()){

    mf::LogInfo("CAFMaker") << "Output filename is " << fCafFilename;

    fFile = new TFile(fCafFilename.c_str(), "RECREATE");
    
    fRecTree = new TTree("recTree", "records");

    // Tell the tree it's expecting StandardRecord objects
    StandardRecord* rec = 0;
    fRecTree->Branch("rec", "caf::StandardRecord", &rec);

    AddEnvToFile(fFile);
 
    if (fParams.CreateBlindedCAF()) {
      mf::LogInfo("CAFMaker") << "Blinded output filenames are " << fCafBlindFilename << ", and " << fCafPrescaleFilename;
      fFileb = new TFile(fCafBlindFilename.c_str(), "RECREATE");
      fRecTreeb = new TTree("recTree", "records");
      fRecTreeb->Branch("rec", "caf::StandardRecord", &rec);

      fFilep = new TFile(fCafPrescaleFilename.c_str(), "RECREATE");
      fRecTreep = new TTree("recTree", "records");
      fRecTreep->Branch("rec", "caf::StandardRecord", &rec);

      AddEnvToFile(fFileb);
      AddEnvToFile(fFilep);
    }

  }     

  if(fParams.CreateFlatCAF()){
    mf::LogInfo("CAFMaker") << "Output flat filename is " << fFlatCafFilename;

    // LZ4 is the fastest format to decompress. I get 3x faster loading with
    // this compared to the default, and the files are only slightly larger.
    fFlatFile = new TFile(fFlatCafFilename.c_str(), "RECREATE", "",
                          ROOT::CompressionSettings(ROOT::kLZ4, 1));

    fFlatTree = new TTree("recTree", "recTree");

    fFlatRecord = new flat::Flat<caf::StandardRecord>(fFlatTree, "rec", "", 0);

    AddEnvToFile(fFlatFile);

    if (fParams.CreateBlindedCAF()) {

      mf::LogInfo("CAFMaker") << "Blinded output flat filename are " << fFlatCafBlindFilename << ", and " << fFlatCafPrescaleFilename;

      // LZ4 is the fastest format to decompress. I get 3x faster loading with
      // this compared to the default, and the files are only slightly larger.
      fFlatFileb = new TFile(fFlatCafBlindFilename.c_str(), "RECREATE", "",
			     ROOT::CompressionSettings(ROOT::kLZ4, 1));
      fFlatTreeb = new TTree("recTree", "recTree");
      fFlatRecordb = new flat::Flat<caf::StandardRecord>(fFlatTreeb, "rec", "", 0);
      AddEnvToFile(fFlatFileb);

      fFlatFilep = new TFile(fFlatCafPrescaleFilename.c_str(), "RECREATE", "",
			     ROOT::CompressionSettings(ROOT::kLZ4, 1));
      fFlatTreep = new TTree("recTree", "recTree");
      fFlatRecordp = new flat::Flat<caf::StandardRecord>(fFlatTreep, "rec", "", 0);
      AddEnvToFile(fFlatFilep);
    }

  }

  fFileNumber = -1;
  fTotalPOT = 0;
  fSubRunPOT = 0;
  fTotalSinglePOT = 0;
  fTotalEvents = 0;
  fBlindEvents = 0;
  fPrescaleEvents = 0;
  fIndexInFile = SRHeader::NoSourceIndex();
  fFirstInSubRun = false;
  // fCycle = -5;
  // fBatch = -5;
}


//......................................................................
template <class T, class U>
art::FindManyP<T> CAFMaker::FindManyPStrict(const U& from,
                                            const art::Event& evt,
                                            const art::InputTag& tag) const {
  art::FindManyP<T> ret(from, evt, tag);

  if (!tag.label().empty() && !ret.isValid() && fParams.StrictMode()) {
    std::cout << "CAFMaker: No Assn from '"
              << cet::demangle_symbol(typeid(from).name()) << "' to '"
              << cet::demangle_symbol(typeid(T).name())
              << "' found under label '" << tag << "'. "
              << "Set 'StrictMode: false' to continue anyway." << std::endl;
    abort();
  }

  return ret;
}

//......................................................................
template <class T, class D, class U>
art::FindManyP<T, D> CAFMaker::FindManyPDStrict(const U& from,
                                            const art::Event& evt,
                                            const art::InputTag& tag) const {
  art::FindManyP<T, D> ret(from, evt, tag);

  if (!tag.label().empty() && !ret.isValid() && fParams.StrictMode()) {
    std::cout << "CAFMaker: No Assn from '"
              << cet::demangle_symbol(typeid(from).name()) << "' to '"
              << cet::demangle_symbol(typeid(T).name())
              << "' found under label '" << tag << "'. "
              << "Set 'StrictMode: false' to continue anyway." << std::endl;
    abort();
  }

  return ret;
}

//......................................................................
template <class T, class U>
art::FindOneP<T> CAFMaker::FindOnePStrict(const U& from,
					  const art::Event& evt,
					  const art::InputTag& tag) const {
  art::FindOneP<T> ret(from, evt, tag);

  if (!tag.label().empty() && !ret.isValid() && fParams.StrictMode()) {
    std::cout << "CAFMaker: No Assn from '"
              << cet::demangle_symbol(typeid(from).name()) << "' to '"
              << cet::demangle_symbol(typeid(T).name())
              << "' found under label '" << tag << "'. "
              << "Set 'StrictMode: false' to continue anyway." << std::endl;
    abort();
  }

  return ret;
}

//......................................................................
template <class T>
bool CAFMaker::GetAssociatedProduct(const art::FindManyP<T>& fm, int idx,
                                    T& ret) const {
  if (!fm.isValid()) return false;

  const std::vector<art::Ptr<T>> prods = fm.at(idx);

  if (prods.empty()) return false;

  ret = *prods[0];

  return true;
}

//......................................................................
template <class EvtT, class T>
void CAFMaker::GetByLabelStrict(const EvtT& evt, const std::string& label,
                                art::Handle<T>& handle) const {
  evt.getByLabel(label, handle);
  if (!label.empty() && handle.failedToGet() && fParams.StrictMode()) {
    std::cout << "CAFMaker: No product of type '"
              << cet::demangle_symbol(typeid(*handle).name())
              << "' found under label '" << label << "'. "
              << "Set 'StrictMode: false' to continue anyway." << std::endl;
    abort();
  }
}

//......................................................................
template <class T>
void CAFMaker::GetByLabelIfExists(const art::Event& evt,
                                  const std::string& label,
                                  art::Handle<T>& handle) const {
  evt.getByLabel(label, handle);
  if (!label.empty() && handle.failedToGet() && fParams.StrictMode()) {
    std::cout << "CAFMaker: No product of type '"
              << cet::demangle_symbol(typeid(*handle).name())
              << "' found under label '" << label << "'. "
              << "Continuing without it." << std::endl;
  }
}

//......................................................................
template <class T>
bool CAFMaker::GetPsetParameter(const fhicl::ParameterSet& pset,
                                const std::vector<std::string>& name,
                                T& ret) const {
  fhicl::ParameterSet p = pset;
  for (unsigned int i = 0; i < name.size() - 1; ++i) {
    if (!p.has_key(name[i])) return false;
    p = p.get<fhicl::ParameterSet>(name[i]);
  }
  if (!p.has_key(name.back())) return false;
  ret = p.get<T>(name.back());
  return true;
}

//......................................................................
void CAFMaker::produce(art::Event& evt) noexcept {

  bool const firstInFile = (fIndexInFile++ == 0);
  
  // is this event real data?
  bool isRealData = evt.isRealData();

  std::unique_ptr<std::vector<caf::StandardRecord>> srcol(
      new std::vector<caf::StandardRecord>);

  std::unique_ptr<art::Assns<caf::StandardRecord, recob::Slice>> srAssn(
      new art::Assns<caf::StandardRecord, recob::Slice>);

  fTotalEvents += 1;

  // get all the truth's
  art::Handle<std::vector<simb::MCTruth>> mctruth_handle;
  GetByLabelStrict(evt, fParams.GenLabel(), mctruth_handle);

  std::vector<art::Ptr<simb::MCTruth>> mctruths;
  if (mctruth_handle.isValid()) {
    art::fill_ptr_vector(mctruths, mctruth_handle);
  }

  // And associated GTruth objects
  art::FindManyP<simb::GTruth> fmp_gtruth = FindManyPStrict<simb::GTruth>(mctruths, evt, fParams.GenLabel());

  art::Handle<std::vector<simb::MCTruth>> cosmic_mctruth_handle;
  evt.getByLabel(fParams.CosmicGenLabel(), cosmic_mctruth_handle);

  art::Handle<std::vector<simb::MCTruth>> pgun_mctruth_handle;
  evt.getByLabel(fParams.ParticleGunGenLabel(), pgun_mctruth_handle);

  // use the MCTruth to determine the simulation type
  caf::MCType_t mctype = caf::kMCUnknown;
  if (mctruth_handle.isValid() && cosmic_mctruth_handle.isValid()) {
    mctype = caf::kMCOverlay;
  }
  else if (mctruth_handle.isValid()) {
    mctype = caf::kMCNeutrino;
  }
  else if (cosmic_mctruth_handle.isValid()) {
    mctype = caf::kMCCosmic;
  }
  else if (pgun_mctruth_handle.isValid()) {
    mctype = caf::kMCParticleGun;
  }

  // Lookup the MeV-Portal info if it is there
  //
  // Don't be "strict" because this will only be true for a subset of MC
  art::Handle<std::vector<evgen::ldm::MeVPrtlTruth>> mevprtltruth_handle;
  evt.getByLabel(fParams.GenLabel(), mevprtltruth_handle);

  std::vector<art::Ptr<evgen::ldm::MeVPrtlTruth>> mevprtl_truths;
  if (mevprtltruth_handle.isValid()) art::fill_ptr_vector(mevprtl_truths, mevprtltruth_handle);

  // prepare map of track ID's to energy depositions
  art::Handle<std::vector<sim::SimChannel>> simchannel_handle;
  GetByLabelStrict(evt, fParams.SimChannelLabel(), simchannel_handle);

  std::vector<art::Ptr<sim::SimChannel>> simchannels;
  if (simchannel_handle.isValid()) {
    art::fill_ptr_vector(simchannels, simchannel_handle);
  }

  art::Handle<std::vector<simb::MCFlux>> mcflux_handle;
  GetByLabelStrict(evt, std::string("generator"), mcflux_handle);

  std::vector<art::Ptr<simb::MCFlux>> mcfluxes;
  if (mcflux_handle.isValid()) {
    art::fill_ptr_vector(mcfluxes, mcflux_handle);
  }

  // get the MCReco for the fake-reco
  art::Handle<std::vector<sim::MCTrack>> mctrack_handle;
  GetByLabelStrict(evt, std::string("mcreco"), mctrack_handle);
  std::vector<art::Ptr<sim::MCTrack>> mctracks;
  if (mctrack_handle.isValid()) {
    art::fill_ptr_vector(mctracks, mctrack_handle);
  }

  // get all of the true particles from G4
  std::vector<caf::SRTrueParticle> true_particles;
  art::Handle<std::vector<simb::MCParticle>> mc_particles;
  GetByLabelStrict(evt, fParams.G4Label(), mc_particles);

  // collect services
  // Moved ParticleInventory and BackTracker services definition as needed elsewhere (BH)
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  auto const *sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  // Collect the input TPC reco tags
  std::vector<std::string> pandora_tag_suffixes;
  fParams.PandoraTagSuffixes(pandora_tag_suffixes);
  if (pandora_tag_suffixes.size() == 0) pandora_tag_suffixes.push_back("");

  // collect the TPC hits
  std::vector<art::Ptr<recob::Hit>> hits;
  for (unsigned i_tag = 0; i_tag < pandora_tag_suffixes.size(); i_tag++) {
    const std::string &pandora_tag_suffix = pandora_tag_suffixes[i_tag];
    art::Handle<std::vector<recob::Hit>> thisHits;
    GetByLabelStrict(evt, fParams.HitLabel() + pandora_tag_suffix, thisHits);
    if (thisHits.isValid()) {
      art::fill_ptr_vector(hits, thisHits);
    }
  }

  // Prep truth-to-reco-matching info
  std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map;
  std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map;
  std::map<int, caf::HitsEnergy> id_to_hit_energy_map;

  if ( !isRealData ) {
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    id_to_ide_map = PrepSimChannels(simchannels, *geometry);
    id_to_truehit_map = PrepTrueHits(hits, clock_data, *bt_serv);
    id_to_hit_energy_map = SetupIDHitEnergyMap(hits, clock_data, *bt_serv);
  }

  //#######################################################
  // Fill truths & fake reco
  //#######################################################

  caf::SRTruthBranch                  srtruthbranch;

  if (mc_particles.isValid()) {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    for (const simb::MCParticle part: *mc_particles) {
      true_particles.emplace_back();

      FillTrueG4Particle(part,
                         fActiveVolumes,
                         fTPCVolumes,
                         id_to_ide_map,
                         id_to_truehit_map,
                         *bt_serv,
                         *pi_serv,
                         mctruths,
                         true_particles.back());
    }
  }

  std::vector<art::FindManyP<sbn::evwgh::EventWeightMap>> fmpewm;

  // holder for invalid MCFlux
  simb::MCFlux badflux; // default constructor gives nonsense values

  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);
    const simb::MCFlux &mcflux = (mcfluxes.size()) ? *mcfluxes.at(i) : badflux;

    simb::GTruth gtruth;
    bool ok = GetAssociatedProduct(fmp_gtruth, i, gtruth);
    if(!ok){
      std::cout << "Failed to get GTruth object!" << std::endl;
    }

    srtruthbranch.nu.push_back(SRTrueInteraction());
    srtruthbranch.nnu ++;

    if ( !isRealData ) FillTrueNeutrino(mctruth, mcflux, gtruth, true_particles, id_to_truehit_map, srtruthbranch.nu.back(), i, fActiveVolumes);

    // Don't check for syst weight assocations until we have something (MCTruth
    // corresponding to a neutrino) that could plausibly be reweighted. This
    // avoids the need for special configuration for cosmics or single particle
    // simulation, and real data.
    if(fmpewm.empty() && mctruth->NeutrinoSet()){
      for(const std::string& label: fParams.SystWeightLabels()){
        fmpewm.push_back(FindManyPStrict<sbn::evwgh::EventWeightMap>(mctruths, evt, label));
      }
    }

    // For each of the sources of systematic weights
    for(auto& fm: fmpewm){
      if (!fm.isValid()) continue; // Don't crash if StrictMode==false

      // Find the weights associated with this particular interaction
      const std::vector<art::Ptr<sbn::evwgh::EventWeightMap>> wgts = fm.at(i);

      // For all the weights associated with this MCTruth
      for(const art::Ptr<sbn::evwgh::EventWeightMap>& wgtmap: wgts){
        FillEventWeight(*wgtmap, srtruthbranch.nu.back(), fWeightPSetIndex);
      } // end for wgtmap
    } // end for fm
  } // end for i (mctruths)

  // get the number of events generated in the gen stage
  unsigned n_gen_evt = 0;
  for (const art::ProcessConfiguration &process: evt.processHistory()) {
    fhicl::ParameterSet gen_config;
    bool success = evt.getProcessParameterSet(process.processName(), gen_config);
    if (success && gen_config.has_key("source") && gen_config.has_key("source.maxEvents") && gen_config.has_key("source.module_type") ) {
      int max_events = gen_config.get<int>("source.maxEvents");
      std::string module_type = gen_config.get<std::string>("source.module_type");
      if (module_type == "EmptyEvent") {
        n_gen_evt += max_events;
      }
    }
  }

  std::vector<caf::SRFakeReco> srfakereco;
  FillFakeReco(mctruths, true_particles, mctracks, fActiveVolumes, *fFakeRecoTRandom, srfakereco);

  // Fill the MeVPrtl stuff
  for (unsigned i_prtl = 0; i_prtl < mevprtl_truths.size(); i_prtl++) {
    srtruthbranch.prtl.emplace_back();
    FillMeVPrtlTruth(*mevprtl_truths[i_prtl], fActiveVolumes, srtruthbranch.prtl.back());
    srtruthbranch.nprtl = srtruthbranch.prtl.size();
  } 

  //#######################################################
  // Fill detector & reco
  //#######################################################

  //Beam gate and Trigger info
  art::Handle<sbn::ExtraTriggerInfo> extratrig_handle;
  GetByLabelStrict(evt, fParams.TriggerLabel().encode(), extratrig_handle);

  art::Handle<std::vector<raw::Trigger>> trig_handle;
  GetByLabelStrict(evt, fParams.TriggerLabel().encode(), trig_handle);

  caf::SRTrigger srtrigger; 
  if (extratrig_handle.isValid() && trig_handle.isValid() && trig_handle->size() == 1) {
    FillTrigger(*extratrig_handle, trig_handle->at(0), srtrigger);
  }
  // If not real data, fill in enough of the SRTrigger to make (e.g.) the CRT 
  // time referencing work. TODO: add more stuff to a "MC"-Trigger?
  else if(!isRealData) {
    FillTriggerMC(fParams.CRTSimT0Offset(), srtrigger);
  }

  // try to find the result of the Flash trigger if it was run
  bool pass_flash_trig = false;
  art::Handle<bool> flashtrig_handle;
  GetByLabelStrict(evt, fParams.FlashTrigLabel(), flashtrig_handle);

  if (flashtrig_handle.isValid()) {
    pass_flash_trig = *flashtrig_handle;
  }

  // Fill various detector information associated with the event
  //
  // Get all of the CRT hits
  std::vector<caf::SRCRTHit> srcrthits;

  art::Handle<std::vector<sbn::crt::CRTHit>> crthits_handle;
  GetByLabelStrict(evt, fParams.CRTHitLabel(), crthits_handle);
  // fill into event
  //int64_t CRT_T0_reference_time = fParams.ReferenceCRTT0ToBeam() ? -srtrigger.beam_gate_time_abs : 0; // ns, signed
  //double CRT_T1_reference_time = fParams.ReferenceCRTT1FromTriggerToBeam() ? srtrigger.trigger_within_gate : 0.;
  int64_t CRT_T0_reference_time = isRealData ?  -srtrigger.beam_gate_time_abs : -fParams.CRTSimT0Offset();
  double CRT_T1_reference_time = isRealData ? srtrigger.trigger_within_gate : -fParams.CRTSimT0Offset();

  if (crthits_handle.isValid()) {
    const std::vector<sbn::crt::CRTHit> &crthits = *crthits_handle;
    for (unsigned i = 0; i < crthits.size(); i++) {
      srcrthits.emplace_back();
      FillCRTHit(crthits[i], fParams.CRTUseTS0(), CRT_T0_reference_time, CRT_T1_reference_time, srcrthits.back());
    }
  }

  // Get all of the CRT Tracks
  std::vector<caf::SRCRTTrack> srcrttracks;

  art::Handle<std::vector<sbn::crt::CRTTrack>> crttracks_handle;
  GetByLabelStrict(evt, fParams.CRTTrackLabel(), crttracks_handle);
  // fill into event
  if (crttracks_handle.isValid()) {
    const std::vector<sbn::crt::CRTTrack> &crttracks = *crttracks_handle;
    for (unsigned i = 0; i < crttracks.size(); i++) {
      srcrttracks.emplace_back();
      FillCRTTrack(crttracks[i], fParams.CRTUseTS0(), srcrttracks.back());
    }
  }

  // Get all of the CRTPMT Matches .. 
  std::vector<caf::SRCRTPMTMatch> srcrtpmtmatches;
  std::cout << "srcrtpmtmatches.size = " << srcrtpmtmatches.size() << "\n";
  art::Handle<std::vector<sbn::crt::CRTPMTMatching>> crtpmtmatch_handle;
  GetByLabelStrict(evt, fParams.CRTPMTLabel(), crtpmtmatch_handle);
  if(crtpmtmatch_handle.isValid()){
    std::cout << "valid handle! label: " << fParams.CRTPMTLabel() << "\n";
    const std::vector<sbn::crt::CRTPMTMatching> &crtpmtmatches = *crtpmtmatch_handle;
    for (unsigned i = 0; i < crtpmtmatches.size(); i++) {
      srcrtpmtmatches.emplace_back();
      FillCRTPMTMatch(crtpmtmatches[i],srcrtpmtmatches.back());
    }
  }
  else{
    std::cout << "crtpmtmatch_handle.isNOTValid!\n";
  }

  // Get all of the OpFlashes
  std::vector<caf::SROpFlash> srflashes;

  for (const std::string& pandora_tag_suffix : pandora_tag_suffixes) {
    art::Handle<std::vector<recob::OpFlash>> flashes_handle;
    GetByLabelStrict(evt, fParams.OpFlashLabel() + pandora_tag_suffix, flashes_handle);
    // fill into event
    if (flashes_handle.isValid()) {
      const std::vector<recob::OpFlash> &opflashes = *flashes_handle;
      int cryostat = ( pandora_tag_suffix.find("W") != std::string::npos ) ? 1 : 0;

      // get associated OpHits for each OpFlash
      art::FindMany<recob::OpHit> findManyHits(flashes_handle, evt, fParams.OpFlashLabel() + pandora_tag_suffix);

      int iflash=0;
      for (const recob::OpFlash& flash : opflashes) {

        std::vector<recob::OpHit const*> const& ophits = findManyHits.at(iflash);

        srflashes.emplace_back();
        FillOpFlash(flash, ophits, cryostat, srflashes.back());
        iflash++;
      }
    }
  }

  // collect the TPC slices
  std::vector<art::Ptr<recob::Slice>> slices;
  std::vector<std::string> slice_tag_suffixes;
  std::vector<unsigned> slice_tag_indices;
  for (unsigned i_tag = 0; i_tag < pandora_tag_suffixes.size(); i_tag++) {
    const std::string &pandora_tag_suffix = pandora_tag_suffixes[i_tag];
    // Get a handle on the slices
    art::Handle<std::vector<recob::Slice>> thisSlices;
    GetByLabelStrict(evt, fParams.PFParticleLabel() + pandora_tag_suffix, thisSlices);
    if (thisSlices.isValid()) {
      art::fill_ptr_vector(slices, thisSlices);
      for (unsigned i = 0; i < thisSlices->size(); i++) {
        slice_tag_suffixes.push_back(pandora_tag_suffix);
        slice_tag_indices.push_back(i_tag);
      }
    }
  }

  // The Standard Record
  // Branch entry definition -- contains list of slices, CRT information, and truth information
  StandardRecord rec;

  //#######################################################
  // Loop over slices
  //#######################################################
  for (unsigned sliceID = 0; sliceID < slices.size(); sliceID++) {
    // Holder for information on this slice
    caf::SRSlice recslc;
    recslc.truth.det = fDet;

    art::Ptr<recob::Slice> slice = slices[sliceID];
    const std::string &slice_tag_suff = slice_tag_suffixes[sliceID];
    unsigned producer = slice_tag_indices[sliceID];

    // Get tracks & showers here
    std::vector<art::Ptr<recob::Slice>> sliceList {slice};
    art::FindManyP<recob::PFParticle> findManyPFParts =
       FindManyPStrict<recob::PFParticle>(sliceList, evt,  fParams.PFParticleLabel() + slice_tag_suff);

    std::vector<art::Ptr<recob::PFParticle>> fmPFPart;
    if (findManyPFParts.isValid()) {
      fmPFPart = findManyPFParts.at(0);
    }

    art::FindManyP<recob::Hit> fmSlcHits =
      FindManyPStrict<recob::Hit>(sliceList, evt, fParams.PFParticleLabel() + slice_tag_suff);

    std::vector<art::Ptr<recob::Hit>> slcHits;
    if (fmSlcHits.isValid()) {
      slcHits = fmSlcHits.at(0);
    }

    art::FindOneP<sbn::CRUMBSResult> foSlcCRUMBS =
      FindOnePStrict<sbn::CRUMBSResult>(sliceList, evt,
          fParams.CRUMBSLabel() + slice_tag_suff);
    const sbn::CRUMBSResult *slcCRUMBS = nullptr;
    if (foSlcCRUMBS.isValid()) {
      slcCRUMBS = foSlcCRUMBS.at(0).get();
    }

    art::FindManyP<sbn::SimpleFlashMatch> fm_sFM =
      FindManyPStrict<sbn::SimpleFlashMatch>(fmPFPart, evt,
                                             fParams.FlashMatchLabel() + slice_tag_suff);

    art::FindManyP<larpandoraobj::PFParticleMetadata> fmPFPMeta =
      FindManyPStrict<larpandoraobj::PFParticleMetadata>(fmPFPart, evt,
               fParams.PFParticleLabel() + slice_tag_suff);

    art::FindManyP<recob::SpacePoint> fmSpacePoint =
      FindManyPStrict<recob::SpacePoint>(slcHits, evt, fParams.PFParticleLabel() + slice_tag_suff);

    std::vector<art::Ptr<recob::SpacePoint>> slcSpacePoints;
    if (fmSpacePoint.isValid()) {
      for (unsigned i = 0; i < fmSpacePoint.size(); i++) {
        const std::vector<art::Ptr<recob::SpacePoint>> &thisSpacePoints = fmSpacePoint.at(i);
        if (thisSpacePoints.size() == 0) {
          slcSpacePoints.emplace_back(); // nullptr
        }
        else if (thisSpacePoints.size() == 1) {
          slcSpacePoints.push_back(fmSpacePoint.at(i).at(0));
        }
        else abort();
      }
    }

    art::FindManyP<recob::PFParticle> fmSpacePointPFPs =
      FindManyPStrict<recob::PFParticle>(slcSpacePoints, evt, fParams.PFParticleLabel() + slice_tag_suff);

    art::FindManyP<recob::Shower> fmShower =
      FindManyPStrict<recob::Shower>(fmPFPart, evt, fParams.RecoShowerLabel() + slice_tag_suff);

    // make Ptr's to showers for shower -> other object associations
    std::vector<art::Ptr<recob::Shower>> slcShowers;
    if (fmShower.isValid()) {
      for (unsigned i = 0; i < fmShower.size(); i++) {
        const std::vector<art::Ptr<recob::Shower>> &thisShowers = fmShower.at(i);
        if (thisShowers.size() == 0) {
          slcShowers.emplace_back(); // nullptr
        }
        else if (thisShowers.size() == 1) {
          slcShowers.push_back(fmShower.at(i).at(0));
        }
        else assert(false); // bad
      }
    }

    art::FindManyP<float> fmShowerCosmicDist =
      FindManyPStrict<float>(slcShowers, evt, fParams.ShowerCosmicDistLabel() + slice_tag_suff);

    art::FindManyP<float> fmShowerResiduals =
      FindManyPStrict<float>(slcShowers, evt, fParams.RecoShowerSelectionLabel() + slice_tag_suff);

    art::FindManyP<sbn::ShowerTrackFit> fmShowerTrackFit =
      FindManyPStrict<sbn::ShowerTrackFit>(slcShowers, evt, fParams.RecoShowerSelectionLabel() + slice_tag_suff);

    art::FindManyP<sbn::ShowerDensityFit> fmShowerDensityFit =
      FindManyPStrict<sbn::ShowerDensityFit>(slcShowers, evt, fParams.RecoShowerSelectionLabel() + slice_tag_suff);

    art::FindManyP<recob::Track> fmTrack =
      FindManyPStrict<recob::Track>(fmPFPart, evt,
            fParams.RecoTrackLabel() + slice_tag_suff);

    art::FindOneP<anab::T0> f1PFPT0 =
      FindOnePStrict<anab::T0>(fmPFPart, evt,
            fParams.PFParticleLabel() + slice_tag_suff);

    // make Ptr's to tracks for track -> other object associations
    std::vector<art::Ptr<recob::Track>> slcTracks;
    if (fmTrack.isValid()) {
      for (unsigned i = 0; i < fmTrack.size(); i++) {
        const std::vector<art::Ptr<recob::Track>> &thisTracks = fmTrack.at(i);
        if (thisTracks.size() == 0) {
          slcTracks.emplace_back(); // nullptr
        }
        else if (thisTracks.size() == 1) {
          slcTracks.push_back(fmTrack.at(i).at(0));
        }
        else assert(false); // bad
      }
    }

    // Get the stubs!
    art::FindManyP<sbn::Stub> fmSlcStubs =
      FindManyPStrict<sbn::Stub>(sliceList, evt,
          fParams.StubLabel() + slice_tag_suff);

    std::vector<art::Ptr<sbn::Stub>> fmStubs;
    if (fmSlcStubs.isValid()) {
      fmStubs = fmSlcStubs.at(0);
    } 

    // Lookup stubs to overlaid PFP
    art::FindManyP<recob::PFParticle> fmStubPFPs =
      FindManyPStrict<recob::PFParticle>(fmStubs, evt,
          fParams.StubLabel() + slice_tag_suff);
    // and get the stub hits for truth matching
    art::FindManyP<recob::Hit> fmStubHits =
      FindManyPStrict<recob::Hit>(fmStubs, evt,
          fParams.StubLabel() + slice_tag_suff);

    art::FindManyP<anab::Calorimetry> fmCalo =
      FindManyPStrict<anab::Calorimetry>(slcTracks, evt,
           fParams.TrackCaloLabel() + slice_tag_suff);

    art::FindManyP<anab::ParticleID> fmChi2PID =
      FindManyPStrict<anab::ParticleID>(slcTracks, evt,
          fParams.TrackChi2PidLabel() + slice_tag_suff);

    art::FindManyP<sbn::ScatterClosestApproach> fmScatterClosestApproach =
      FindManyPStrict<sbn::ScatterClosestApproach>(slcTracks, evt,
          fParams.TrackScatterClosestApproachLabel() + slice_tag_suff);

    art::FindManyP<sbn::StoppingChi2Fit> fmStoppingChi2Fit =
      FindManyPStrict<sbn::StoppingChi2Fit>(slcTracks, evt,
          fParams.TrackStoppingChi2FitLabel() + slice_tag_suff);

    art::FindManyP<sbn::MVAPID> fmTrackDazzle =
      FindManyPStrict<sbn::MVAPID>(slcTracks, evt,
          fParams.TrackDazzleLabel() + slice_tag_suff);

    art::FindManyP<sbn::MVAPID> fmShowerRazzle =
      FindManyPStrict<sbn::MVAPID>(slcShowers, evt,
          fParams.ShowerRazzleLabel() + slice_tag_suff);

    art::FindManyP<recob::Vertex> fmVertex =
      FindManyPStrict<recob::Vertex>(fmPFPart, evt,
             fParams.PFParticleLabel() + slice_tag_suff);

    art::FindManyP<recob::Hit> fmTrackHit =
      FindManyPStrict<recob::Hit>(slcTracks, evt,
          fParams.RecoTrackLabel() + slice_tag_suff);

    art::FindManyP<recob::Hit> fmShowerHit =
      FindManyPStrict<recob::Hit>(slcShowers, evt,
          fParams.RecoShowerLabel() + slice_tag_suff);

    // NOTE: The sbn::crt::CRTHit is associated to the T0. It's a bit awkward to 
    // access that here, so we do it per-track (see code where fmCRTHitMatch is accessed below)
    art::FindManyP<anab::T0> fmCRTHitMatch =
      FindManyPStrict<anab::T0>(slcTracks, evt,
               fParams.CRTHitMatchLabel() + slice_tag_suff);

    // TODO: also save the sbn::crt::CRTTrack in the matching so that CAFMaker has access to it
    art::FindManyP<anab::T0> fmCRTTrackMatch =
      FindManyPStrict<anab::T0>(slcTracks, evt,
               fParams.CRTTrackMatchLabel() + slice_tag_suff);

    std::vector<art::FindManyP<recob::MCSFitResult>> fmMCSs;
    static const std::vector<std::string> PIDnames {"muon", "pion", "kaon", "proton"};
    for (std::string pid: PIDnames) {
      art::InputTag tag(fParams.TrackMCSLabel() + slice_tag_suff, pid);
      fmMCSs.push_back(FindManyPStrict<recob::MCSFitResult>(slcTracks, evt, tag));
    }

    std::vector<art::FindManyP<sbn::RangeP>> fmRanges;
    static const std::vector<std::string> rangePIDnames {"muon", "pion", "proton"};
    for (std::string pid: rangePIDnames) {
      art::InputTag tag(fParams.TrackRangeLabel() + slice_tag_suff, pid);
      fmRanges.push_back(FindManyPStrict<sbn::RangeP>(slcTracks, evt, tag));
    }

    //    if (slice.IsNoise() || slice.NCell() == 0) continue;
    // Because we don't care about the noise slice and slices with no hits.

    // get the primary particle
    size_t iPart;
    for (iPart = 0; iPart < fmPFPart.size(); ++iPart ) {
      const recob::PFParticle &thisParticle = *fmPFPart[iPart];
      if (thisParticle.IsPrimary()) break;
    }
    // primary particle and meta-data
    const recob::PFParticle *primary = (iPart == fmPFPart.size()) ? NULL : fmPFPart[iPart].get();
    const larpandoraobj::PFParticleMetadata *primary_meta = (iPart == fmPFPart.size()) ? NULL : fmPFPMeta.at(iPart).at(0).get();
    // get the flash match
    const sbn::SimpleFlashMatch* fmatch = nullptr;
    if (fm_sFM.isValid() && primary != NULL) {
      std::vector<art::Ptr<sbn::SimpleFlashMatch>> fmatches = fm_sFM.at(iPart);
      if (fmatches.size() != 0) {
        assert(fmatches.size() == 1);
        fmatch = fmatches[0].get();
      }
    }
    // get the primary vertex
    const recob::Vertex *vertex = (iPart == fmPFPart.size() || !fmVertex.at(iPart).size()) ? NULL : fmVertex.at(iPart).at(0).get();

    //#######################################################
    // Add slice info.
    //#######################################################
    FillSliceVars(*slice, primary, producer, recslc);
    FillSliceMetadata(primary_meta, recslc);
    FillSliceFlashMatch(fmatch, recslc);
    FillSliceFlashMatchA(fmatch, recslc);
    FillSliceVertex(vertex, recslc);
    FillSliceCRUMBS(slcCRUMBS, recslc);
    FillSliceBarycenter(slcHits, slcSpacePoints, recslc);

    // select slice
    if (!SelectSlice(recslc, fParams.CutClearCosmic())) continue;

    // Whether Pandora thinks this slice is a neutrino
    //
    // This requirement is used to determine whether to save additional
    // per-hit information about the slice.
    bool NeutrinoSlice = !recslc.is_clear_cosmic;

    // Fill truth info after decision on selection is made
    if ( !isRealData ) {
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

      FillSliceTruth(slcHits, mctruths, srtruthbranch,
		     *pi_serv, clock_data, recslc);

      FillSliceFakeReco(slcHits, mctruths, srtruthbranch,
			*pi_serv, clock_data, recslc, true_particles, mctracks, 
                        fActiveVolumes, *fFakeRecoTRandom);
    }

    //#######################################################
    // Add detector dependent slice info.
    //#######################################################
    // if (fDet == kSBND) {
    //   rec.sel.contain.nplanestofront = rec.slc.firstplane - (plnfirst - 1);
    //   rec.sel.contain.nplanestoback = (plnlast) - 1 - rec.slc.lastplane;
    // }

    //#######################################################
    // Add stub reconstructed objects.
    //#######################################################
    for (size_t iStub = 0; iStub < fmStubs.size(); iStub++) {
      const sbn::Stub &thisStub = *fmStubs[iStub];

      art::Ptr<recob::PFParticle> thisStubPFP;
      if (!fmStubPFPs.at(iStub).empty()) thisStubPFP = fmStubPFPs.at(iStub).at(0);

      rec.reco.stub.emplace_back();
      FillStubVars(thisStub, thisStubPFP, rec.reco.stub.back());
      if ( !isRealData ) FillStubTruth(fmStubHits.at(iStub), id_to_hit_energy_map, true_particles, clock_data, rec.reco.stub.back());
      rec.reco.nstub = rec.reco.stub.size();

      // Duplicate stub reco info in the srslice
      recslc.reco.stub.push_back(rec.reco.stub.back());
      recslc.reco.nstub = recslc.reco.stub.size();
    }
 
    if (fParams.FillHits()) {
      for ( size_t iHit = 0; iHit < slcHits.size(); ++iHit ) {
        const recob::Hit &thisHit = *slcHits[iHit];
      
        std::vector<art::Ptr<recob::PFParticle>> thisParticle;
        if (fmSpacePointPFPs.isValid()) {
          thisParticle = fmSpacePointPFPs.at(iHit);
        }
        std::vector<art::Ptr<recob::SpacePoint>> thisPoint;
        if (fmSpacePoint.isValid()) {
          thisPoint = fmSpacePoint.at(iHit); 
        }
        if (!thisParticle.empty() && !thisPoint.empty()) {
          assert(thisParticle.size() == 1);
          assert(thisPoint.size() == 1);
          rec.reco.nhit++;
          rec.reco.hit.push_back(SRHit());

          FillHitVars(thisHit, producer, *thisPoint[0], *thisParticle[0], rec.reco.hit.back());
          recslc.reco.hit.push_back(rec.reco.hit.back());
          recslc.reco.nhit = recslc.reco.hit.size();
        }
      }
    }

    //#######################################################
    // Add track/shower reconstructed objects.
    //#######################################################
    // Reco objects have assns to the slice PFParticles
    // This depends on the findMany object created above.
    for ( size_t iPart = 0; iPart < fmPFPart.size(); ++iPart ) {
      const recob::PFParticle &thisParticle = *fmPFPart[iPart];
      std::vector<art::Ptr<recob::Track>> thisTrack;
      if (fmTrack.isValid()) {
        thisTrack = fmTrack.at(iPart);
      }
      std::vector<art::Ptr<recob::Shower>> thisShower;
      if (fmShower.isValid()) {
        thisShower = fmShower.at(iPart);
      }
     
      SRPFP pfp;

      art::Ptr<anab::T0> thisPFPT0;
      if (f1PFPT0.isValid()) {
        thisPFPT0 = f1PFPT0.at(iPart);
      }

      const larpandoraobj::PFParticleMetadata *pfpMeta = (fmPFPMeta.at(iPart).empty()) ? NULL : fmPFPMeta.at(iPart).at(0).get();
      FillPFPVars(thisParticle, primary, pfpMeta, thisPFPT0, pfp);

      if (!thisTrack.empty())  { // it has a track!
        assert(thisTrack.size() == 1);

        // collect all the stuff
        std::array<std::vector<art::Ptr<recob::MCSFitResult>>, 4> trajectoryMCS;
        for (unsigned index = 0; index < 4; index++) {
          if (fmMCSs[index].isValid()) {
            trajectoryMCS[index] = fmMCSs[index].at(iPart);
          }
          else {
            trajectoryMCS[index] = std::vector<art::Ptr<recob::MCSFitResult>>();
          }
        }

        std::array<std::vector<art::Ptr<sbn::RangeP>>, 3> rangePs;
        for (unsigned index = 0; index < 3; index++) {
          if (fmRanges[index].isValid()) {
            rangePs[index] = fmRanges[index].at(iPart);
          }
          else {
            rangePs[index] = std::vector<art::Ptr<sbn::RangeP>>();
          }
        }

        // fill all the stuff
        SRTrack& trk = pfp.trk;
        FillTrackVars(*thisTrack[0], producer, trk);
        FillTrackMCS(*thisTrack[0], trajectoryMCS, trk);
        FillTrackRangeP(*thisTrack[0], rangePs, trk);

        if (fmChi2PID.isValid()) {
           FillTrackChi2PID(fmChi2PID.at(iPart), lar::providerFrom<geo::Geometry>(), trk);
        }
        if (fmScatterClosestApproach.isValid() && fmScatterClosestApproach.at(iPart).size()==1) {
           FillTrackScatterClosestApproach(fmScatterClosestApproach.at(iPart).front(), trk);
        }
        if (fmStoppingChi2Fit.isValid() && fmStoppingChi2Fit.at(iPart).size()==1) {
           FillTrackStoppingChi2Fit(fmStoppingChi2Fit.at(iPart).front(), trk);
        }
        if (fmTrackDazzle.isValid() && fmTrackDazzle.at(iPart).size()==1) {
           FillTrackDazzle(fmTrackDazzle.at(iPart).front(), trk);
        }
        if (fmCalo.isValid()) {
          FillTrackCalo(fmCalo.at(iPart), fmTrackHit.at(iPart),
              (fParams.FillHitsNeutrinoSlices() && NeutrinoSlice) || fParams.FillHitsAllSlices(), 
              fParams.TrackHitFillRRStartCut(), fParams.TrackHitFillRREndCut(),
              lar::providerFrom<geo::Geometry>(), dprop, trk);
        }
        if (fmCRTHitMatch.isValid()) {
          art::FindManyP<sbn::crt::CRTHit> CRTT02Hit = FindManyPStrict<sbn::crt::CRTHit>
              (fmCRTHitMatch.at(iPart), evt, fParams.CRTHitMatchLabel() + slice_tag_suff);
         
          std::vector<art::Ptr<sbn::crt::CRTHit>> crthitmatch;
          if (CRTT02Hit.isValid() && CRTT02Hit.size() == 1) crthitmatch = CRTT02Hit.at(0);

          FillTrackCRTHit(fmCRTHitMatch.at(iPart), crthitmatch, fParams.CRTUseTS0(), CRT_T0_reference_time, CRT_T1_reference_time, trk);
        }
        // NOTE: SEE TODO AT fmCRTTrackMatch
        if (fmCRTTrackMatch.isValid()) {
          FillTrackCRTTrack(fmCRTTrackMatch.at(iPart), trk);
        }
        // Truth matching
        if (fmTrackHit.isValid()) {
          if ( !isRealData ) {
            // Track -> particle matching
            FillTrackTruth(fmTrackHit.at(iPart), id_to_hit_energy_map, true_particles, clock_data, trk);
            // Hit truth information corresponding to Calo-Points
            // Assumes truth matching and calo-points are filled
            if (mc_particles.isValid() && fParams.FillTrackCaloTruth()) FillTrackCaloTruth(id_to_ide_map, *mc_particles, geometry, clock_data, sce, trk);
          }
        }
      } // thisTrack exists

      if (!thisShower.empty()) { // it has shower!
        assert(thisShower.size() == 1);
	
        SRShower& shw = pfp.shw;
        FillShowerVars(*thisShower[0], vertex, fmShowerHit.at(iPart), lar::providerFrom<geo::Geometry>(), producer, shw);

        // We may have many residuals per shower depending on how many showers ar in the slice
        if (fmShowerRazzle.isValid() && fmShowerRazzle.at(iPart).size()==1) {
           FillShowerRazzle(fmShowerRazzle.at(iPart).front(), shw);
        }
        if (fmShowerCosmicDist.isValid() && fmShowerCosmicDist.at(iPart).size() != 0) {
          FillShowerCosmicDist(fmShowerCosmicDist.at(iPart), shw);
        }
        if (fmShowerResiduals.isValid() && fmShowerResiduals.at(iPart).size() != 0) {
          FillShowerResiduals(fmShowerResiduals.at(iPart), shw);
        }
        if (fmShowerTrackFit.isValid() && fmShowerTrackFit.at(iPart).size()  == 1) {
          FillShowerTrackFit(*fmShowerTrackFit.at(iPart).front(), shw);
        }
        if (fmShowerDensityFit.isValid() && fmShowerDensityFit.at(iPart).size() == 1) {
          FillShowerDensityFit(*fmShowerDensityFit.at(iPart).front(), shw);
        }
        if (fmShowerHit.isValid()) {
          if ( !isRealData ) FillShowerTruth(fmShowerHit.at(iPart), id_to_hit_energy_map, true_particles, clock_data, shw);
        }

      } // thisShower exists
      
      recslc.reco.pfp.push_back(std::move(pfp));
      recslc.reco.npfp = recslc.reco.pfp.size();
    }// end for pfparts



    //#######################################################
    // Fill slice in rec tree
    //#######################################################

    // // Set mc branch values to default
    // rec.mc.setDefault();
    // if (fParams.EnableBlindness()) BlindThisRecord(&rec);
    //util::CreateAssn(*this, evt, *srcol, art::Ptr<recob::Slice>(slices, sliceID),
    //                 *srAssn);

    rec.slc.push_back(recslc);

  }  // end loop over slices

  //#######################################################
  //  Fill rec Tree
  //#######################################################
  rec.nslc            = rec.slc.size();
  rec.mc              = srtruthbranch;
  rec.fake_reco       = srfakereco;
  rec.nfake_reco      = srfakereco.size();
  rec.pass_flashtrig  = pass_flash_trig;  // trigger result
  rec.crt_hits        = srcrthits;
  rec.ncrt_hits       = srcrthits.size();
  rec.crt_tracks        = srcrttracks;
  rec.ncrt_tracks       = srcrttracks.size();
  rec.opflashes       = srflashes;
  rec.nopflashes      = srflashes.size();
  if (fParams.FillTrueParticles()) {
    rec.true_particles  = true_particles;
  }
  rec.ntrue_particles = true_particles.size();
  rec.crtpmt_matches = srcrtpmtmatches;
  rec.ncrtpmt_matches = srcrtpmtmatches.size();

  // Fix the Reference time
  //
  // We want MC and Data to have the same reference time.
  // In MC/LArSoft the "reference time" is canonically defined
  // as the time when the start of the beam spill reaches the detector.
  //
  // In data it may be defined differently for different subsystems. In 
  // particular, some sub-systems define the reference time as the time
  // of the trigger. We want to correct those to the universal reference
  // time from MC.
  //
  // PMT's:
  double PMT_reference_time = fParams.ReferencePMTFromTriggerToBeam() ? srtrigger.trigger_within_gate : 0.;
  FixPMTReferenceTimes(rec, PMT_reference_time);

  // TODO: TPC?

  // Get metadata information for header
  unsigned int run = evt.run();
  unsigned int subrun = evt.subRun();
  unsigned int evtID = evt.event();
  //   unsigned int spillNum = evt.id().event();

  rec.hdr = SRHeader();

  // Get the Process and Cluser number
  const char *process_str = std::getenv("PROCESS");
  if (process_str) {
    try {
      rec.hdr.proc = std::stoi(process_str);
    }
    catch (...) {}
  }

  const char *cluster_str = std::getenv("CLUSTER");
  if (cluster_str) {
    try {
      rec.hdr.cluster = std::stoi(cluster_str);
    }
    catch (...) {}
  }

  rec.hdr.run     = run;
  rec.hdr.subrun  = subrun;
  rec.hdr.evt     = evtID;
  // rec.hdr.subevt = sliceID;
  rec.hdr.ismc    = !isRealData;
  rec.hdr.det     = fDet;
  rec.hdr.fno     = fFileNumber;
  if(firstInFile)
  {
    rec.hdr.nbnbinfo = fBNBInfo.size();
    rec.hdr.bnbinfo = fBNBInfo;
    rec.hdr.nnumiinfo = fNuMIInfo.size();
    rec.hdr.numiinfo = fNuMIInfo;
    rec.hdr.noffbeambnb = fOffbeamBNBGates;
    rec.hdr.noffbeamnumi = fOffbeamNuMIGates;
    rec.hdr.pot   = fSubRunPOT;
  }
  
  rec.hdr.ngenevt = n_gen_evt;
  rec.hdr.mctype  = mctype;
  rec.hdr.sourceName = fSourceFile;
  rec.hdr.sourceIndex = fIndexInFile;
  rec.hdr.first_in_file = firstInFile;
  rec.hdr.first_in_subrun = fFirstInSubRun;
  rec.hdr.triggerinfo = srtrigger;
  // rec.hdr.cycle = fCycle;
  // rec.hdr.batch = fBatch;
  // rec.hdr.blind = 0;
  // rec.hdr.filt = rb::IsFiltered(evt, slices, sliceID);

  // Fill the header info for the given event's spill quality info
  if ( fHasBNBInfo && fHasNuMIInfo ) {
    std::cout << "Found > 0 BNBInfo size and NuMIInfo size, which seems strange. Throwing..." << std::endl;
    abort();
  }
  else if ( fHasBNBInfo ) {
    if ( fBNBInfoEventMap.find(evt.id().event()) == fBNBInfoEventMap.end() )
      std::cout << "We think (since fHasBNBInfo) that this event should be BNB, but did not find this event in the spill info map." << std::endl;
    else rec.hdr.spillbnbinfo = makeSRBNBInfo(fBNBInfoEventMap.at(evt.id().event()));
  }
  else if ( fHasNuMIInfo ) {
    if ( fNuMIInfoEventMap.find(evt.id().event()) == fNuMIInfoEventMap.end() )
      std::cout << "We think (since fHasNuMIInfo) that this event should be NuMI, but did not find this event in the spill info map." << std::endl;
    else rec.hdr.spillnumiinfo = makeSRNuMIInfo(fNuMIInfoEventMap.at(evt.id().event()));
  }

  if(fRecTree){
    // Save the standard-record
    StandardRecord* prec = &rec;
    fRecTree->SetBranchAddress("rec", &prec);
    fRecTree->Fill();

    if(fFlatTree){
      fFlatRecord->Clear();
      fFlatRecord->Fill(rec);
      fFlatTree->Fill();
    }

    //Generate random number to decide if event is saved in prescale or blinded file
    if (fParams.CreateBlindedCAF()) {
      const bool keepprescale = fBlindTRandom->Uniform() < 1/fParams.PrescaleFactor();
      rec.hdr.evt = 0;
      rec.hdr.isblind = true;
      if (keepprescale) {
      	StandardRecord* precp = new StandardRecord (*prec);
	if (fFirstPrescaleInFile) {
	  precp->hdr.pot = fSubRunPOT*(1/fParams.PrescaleFactor());
	  precp->hdr.first_in_file = true;
	  precp->hdr.first_in_subrun = true;
	  precp->hdr.nbnbinfo = fBNBInfo.size()*(1/fParams.PrescaleFactor());
	  precp->hdr.nnumiinfo = fNuMIInfo.size()*(1/fParams.PrescaleFactor());
	}
	precp->hdr.ngenevt = n_gen_evt*(1/fParams.PrescaleFactor());
	precp->hdr.evt = evtID;
	fRecTreep->SetBranchAddress("rec", &precp);
      	fRecTreep->Fill();
	fPrescaleEvents += 1;
	if (fFlatTreep) {
	  fFlatRecordp->Clear();
	  fFlatRecordp->Fill(*precp);
	  fFlatTreep->Fill();
	}
	fFirstPrescaleInFile = false;
      }
      else {
	StandardRecord* precb = new StandardRecord (*prec);
	BlindEnergyParameters(precb);
	if (fFirstBlindInFile) {
	  precb->hdr.pot = fSubRunPOT*(1-(1/fParams.PrescaleFactor()))*GetBlindPOTScale();
	  precb->hdr.first_in_file = true;
	  precb->hdr.first_in_subrun = true;
	  precb->hdr.nbnbinfo = fBNBInfo.size()*(1 - (1/fParams.PrescaleFactor()));
	  precb->hdr.nnumiinfo = fNuMIInfo.size()*(1-(1/fParams.PrescaleFactor()));
	}
	precb->hdr.ngenevt = n_gen_evt*(1 - (1/fParams.PrescaleFactor()));
	precb->hdr.evt = evtID;
	fRecTreeb->SetBranchAddress("rec", &precb);
	fRecTreeb->Fill();
	fBlindEvents += 1;
	if (fFlatTreeb) {
	  fFlatRecordb->Clear();
	  fFlatRecordb->Fill(*precb);
	  fFlatTreeb->Fill();
	}
	fFirstBlindInFile = false;
      }
    }
  }

// reset
  fFirstInSubRun = false;
  srcol->push_back(rec);
  evt.put(std::move(srcol));

  fBNBInfo.clear();
  fNuMIInfo.clear();
  rec.hdr.pot = 0;
}

void CAFMaker::endSubRun(art::SubRun& sr) {

}

//......................................................................
  void CAFMaker::AddHistogramsToFile(TFile* outfile,bool isBlindPOT = false, bool isPrescalePOT = false) const
{

  outfile->cd();

  TH1* hPOT = new TH1D("TotalPOT", "TotalPOT;; POT", 1, 0, 1);
  TH1* hEvents = new TH1D("TotalEvents", "TotalEvents;; Events", 1, 0, 1);

  if (isBlindPOT) {
    hPOT->Fill(0.5,fTotalPOT*(1-(1/fParams.PrescaleFactor()))*GetBlindPOTScale());
  }
  else if (isPrescalePOT) {
    hPOT->Fill(0.5,fTotalPOT*(1/fParams.PrescaleFactor()));
  }
  else {
    hPOT->Fill(0.5,fTotalPOT);
  }
  hEvents->Fill(0.5,fTotalEvents);

  hPOT->Write();
  hEvents->Write();

  if (fParams.CreateBlindedCAF()) {
    TH1*hBlindEvents = new TH1D("BlindEvents", "BlindEvents;; Events", 1, 0, 1);
    TH1* hPrescaleEvents = new TH1D("PrescaleEvents", "PrescaleEvents;; Events", 1, 0, 1);
    hBlindEvents->Fill(0.5, fBlindEvents);
    hPrescaleEvents->Fill(0.5, fPrescaleEvents);
    hBlindEvents->Write();
    hPrescaleEvents->Write();
  }
}

//......................................................................
void CAFMaker::endJob() {
  if (fTotalEvents == 0) {

    std::cerr << "No events processed in this file. Aborting rather than "
                 "produce an empty CAF."
              << std::endl;
    // n.b. changed abort() to return so that eny exceptions thrown during startup
    // still get printed to the user by art
    return;
  }



  if(fFile){

    AddHistogramsToFile(fFile);
    fRecTree->SetDirectory(fFile);
    if (fParams.CreateBlindedCAF()) {
      fRecTreeb->SetDirectory(fFileb);
      fRecTreep->SetDirectory(fFilep);
    }
    fFile->cd();
    fFile->Write();
    if (fParams.CreateBlindedCAF()) {
      AddHistogramsToFile(fFileb,true,false);
      fFileb->cd();
      fFileb->Write();
      AddHistogramsToFile(fFilep,false,true);
      fFilep->cd();
      fFilep->Write();
    }

  }

  if(fFlatFile){

    AddHistogramsToFile(fFlatFile);
    fFlatTree->SetDirectory(fFlatFile);
    if (fParams.CreateBlindedCAF() && fFlatFileb) {
      fFlatTreeb->SetDirectory(fFlatFileb);
      fFlatTreep->SetDirectory(fFlatFilep);
    }
    fFlatFile->cd();
    fFlatFile->Write();
    if (fParams.CreateBlindedCAF()) {
      AddHistogramsToFile(fFlatFileb,true,false);
      fFlatFileb->cd();
      fFlatFileb->Write();
      AddHistogramsToFile(fFlatFilep,false,true);
      fFlatFilep->cd();
      fFlatFilep->Write();
    }
  }

  std::map<std::string, std::string> metamap;

  try{
    art::ServiceHandle<util::MetadataSBN> meta;

    std::map<std::string, std::string> strs;
    std::map<std::string, int> ints;
    std::map<std::string, double> doubles;
    std::map<std::string, std::string> objs;
    meta->GetMetadataMaps(strs, ints, doubles, objs);

    for(auto it: strs) metamap[it.first] = "\""+it.second+"\"";
    for(auto it: ints) metamap[it.first] = std::to_string(it.second);
    for(auto it: doubles) metamap[it.first] = std::to_string(it.second);
    for(auto it: objs) metamap[it.first] = it.second;
  }
  catch(art::Exception& e){//(art::errors::ServiceNotFound)
    // I don't know any way to detect this apart from an exception, unfortunately
    std::cout << "\n\nCAFMaker: TFileMetadataSBN service not configured -- this CAF will not have any metadata saved.\n" << std::endl;
  }

  if(fFile) AddMetadataToFile(fFile, metamap);
  if(fParams.CreateBlindedCAF() && fFileb) AddMetadataToFile(fFileb, metamap);
  if(fParams.CreateBlindedCAF() && fFilep) AddMetadataToFile(fFilep, metamap);
  if(fFlatFile) AddMetadataToFile(fFlatFile, metamap);
  if(fParams.CreateBlindedCAF() && fFlatFileb) AddMetadataToFile(fFlatFileb, metamap);
  if(fParams.CreateBlindedCAF() && fFlatFilep) AddMetadataToFile(fFlatFilep, metamap);
}


}  // end namespace caf
DEFINE_ART_MODULE(caf::CAFMaker)
////////////////////////////////////////////////////////////////////////
