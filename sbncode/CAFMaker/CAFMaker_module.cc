//////////////////////////////////////////////////////////////////
// \file     CAFMaker_module.cc
/// \brief   This module creates Common Analysis Files.
//           Inspired by the NOvA CAFMaker package
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////

// ---------------- TO DO ----------------
//
// - Give real fDet values
// - Add in cycle and batch to params
// - Fill reco tree a bit more
// - Move this list some place useful
// - Add Truth branch
// - Add reco.CRT branch
//
// ---------------------------------------



#include "CAFMakerParams.h"
#include "FillReco.h"
#include "FillTrue.h"
#include "Utils.h"

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

#ifdef DARWINBUILD
#include <libgen.h>
#endif

#include <IFDH_service.h>

// ROOT includes
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TTimeStamp.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "cetlib_except/exception.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

// StandardRecord
#include "sbncode/StandardRecord/StandardRecord.h"

// // CAFMaker 
#include "AssociationUtil.h"
// #include "Blinding.h"

namespace caf {

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

  bool   fIsRealData;  // use this instead of evt.isRealData(), see init in
                     // produce(evt)
  double fTotalPOT;
  double fTotalSinglePOT;
  double fTotalEvents;
  // int fCycle;
  // int fBatch;

  TFile* fFile;
  TTree* fRecTree;

  TH1D* hPOT;
  TH1D* hSinglePOT;
  TH1D* hEvents;

  Det_t fDet;  ///< Detector ID in caf namespace typedef

  // algorithms
  trkf::TrajectoryMCSFitter fMCSCalculator; 
  trkf::TrackMomentumCalculator fRangeCalculator;

  // volumes
  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;

  void InitializeOutfile();

  /// Equivalent of FindManyP except a return that is !isValid() prints a
  /// messsage and aborts if StrictMode is true.
  template <class T, class U>
  art::FindManyP<T> FindManyPStrict(const U& from, const art::Event& evt,
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
  template <class T>
  void GetByLabelStrict(const art::Event& evt, const std::string& label,
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
    fParams(params()), fIsRealData(false), fFile(0),
    fMCSCalculator(fParams.MCSConfig),
    fRangeCalculator(fParams.RangePMinTrackLength())
  {
  fCafFilename = fParams.CAFFilename();

  // Normally CAFMaker is run wit no output ART stream, so these go
  // nowhere, but can be occasionally useful for filtering in ART

  produces<std::vector<caf::StandardRecord>>();
  //produces<art::Assns<caf::StandardRecord, recob::Slice>>();

  // setup volume definitions
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // first the TPC volumes 
  for (auto const &cryo: geometry->IterateCryostats()) {
    geo::GeometryCore::TPC_iterator iTPC = geometry->begin_TPC(cryo.ID()),
                                    tend = geometry->end_TPC(cryo.ID());
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    while (iTPC != tend) {
      geo::TPCGeo const& TPC = *iTPC;
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
      iTPC++;
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
CAFMaker::~CAFMaker() {}

//......................................................................
void CAFMaker::respondToOpenInputFile(const art::FileBlock& fb) {
  if (!fFile) {
    // If Filename wasn't set in the FCL, and this is the
    // first file we've seen
    const int sizef = fb.fileName().size() + 1;
    char* temp  = new char[sizef];
    std::strcpy(temp, fb.fileName().c_str());
    fCafFilename = basename(temp);
    const size_t dotpos = fCafFilename.find('.');
    assert(dotpos != std::string::npos);  // Must have a dot, surely?
    fCafFilename.resize(dotpos);
    fCafFilename += fParams.FileExtension();

    InitializeOutfile();
  }
}

//......................................................................
void CAFMaker::beginJob() {
  if (!fCafFilename.empty()) InitializeOutfile();
}

//......................................................................
void CAFMaker::beginRun(art::Run& r) {
  // fDetID = geom->DetId();
  fDet = (Det_t)1;//(Det_t)fDetID;

}

//......................................................................
void CAFMaker::beginSubRun(art::SubRun& sr) {

}

//......................................................................
void CAFMaker::InitializeOutfile() {
  assert(!fFile);
  assert(!fCafFilename.empty());

  mf::LogInfo("CAFMaker") << "Output filename is " << fCafFilename;

  fFile = new TFile(fCafFilename.c_str(), "RECREATE");

  hPOT = new TH1D("TotalPOT", "TotalPOT;; POT", 1, 0, 1);
  hSinglePOT =
      new TH1D("TotalSinglePOT", "TotalSinglePOT;; Single POT", 1, 0, 1);
  hEvents = new TH1D("TotalEvents", "TotalEvents;; Events", 1, 0, 1);

  fRecTree = new TTree("recTree", "records");

  // Tell the tree it's expecting StandardRecord objects
  StandardRecord* rec = 0;
  fRecTree->Branch("rec", "caf::StandardRecord", &rec);

  fTotalPOT = 0;
  fTotalSinglePOT = 0;
  fTotalEvents = 0;
  // fCycle = -5;
  // fBatch = -5;

  // Global information about the processing details:
  std::map<TString, TString> envmap;
  std::string envstr;
  // Environ comes from unistd.h
  // environ is not present on OSX for some reason, so just use getenv to
  // grab the variables we care about.
#ifdef DARWINBUILD
  std::set<TString> variables;
  variables.insert("USER");
  variables.insert("HOSTNAME");
  variables.insert("PWD");
  variables.insert("SRT_PUBLIC_CONTEXT");
  variables.insert("SRT_PRIVATE_CONTEXT");
  for (auto var : variables) envmap[var] = getenv(var);
#else
  for (char** penv = environ; *penv; ++penv) {
    const std::string pair = *penv;
    envstr += pair;
    envstr += "\n";
    const size_t split = pair.find("=");
    if (split == std::string::npos) continue;  // Huh?
    const std::string key = pair.substr(0, split);
    const std::string value = pair.substr(split + 1);
    envmap[key] = value;
  }
#endif

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

  fFile->mkdir("env")->cd();

  TObjString(envmap["USER"]).Write("user");
  TObjString(envmap["HOSTNAME"]).Write("hostname");
  TObjString(envmap["PWD"]).Write("pwd");
  TObjString(envmap["SRT_PUBLIC_CONTEXT"]).Write("publiccontext");
  TObjString(envmap["SRT_PRIVATE_CONTEXT"]).Write("privatecontext");
  // Default constructor is "now"
  TObjString(TTimeStamp().AsString()).Write("date");
  TObjString(cmd.c_str()).Write("cmd");
  TObjString(fCafFilename.c_str()).Write("output");
  TObjString(envstr.c_str()).Write("env");
}

//......................................................................
template <class T, class U>
art::FindManyP<T> CAFMaker::FindManyPStrict(const U& from,
                                            const art::Event& evt,
                                            const art::InputTag& tag) const {
  art::FindManyP<T> ret(from, evt, tag);

  if (!tag.label().empty() && !ret.isValid() && fParams.StrictMode()) {
    std::cout << "CAFMaker: No Assn from '"
              << abi::__cxa_demangle(typeid(from).name(), 0, 0, 0) << "' to '"
              << abi::__cxa_demangle(typeid(T).name(), 0, 0, 0)
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
template <class T>
void CAFMaker::GetByLabelStrict(const art::Event& evt, const std::string& label,
                                art::Handle<T>& handle) const {
  evt.getByLabel(label, handle);
  if (!label.empty() && handle.failedToGet() && fParams.StrictMode()) {
    std::cout << "CAFMaker: No product of type '"
              << abi::__cxa_demangle(typeid(*handle).name(), 0, 0, 0)
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
              << abi::__cxa_demangle(typeid(*handle).name(), 0, 0, 0)
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

  std::unique_ptr<std::vector<caf::StandardRecord>> srcol(
      new std::vector<caf::StandardRecord>);

  std::unique_ptr<art::Assns<caf::StandardRecord, recob::Slice>> srAssn(
      new art::Assns<caf::StandardRecord, recob::Slice>);

  fTotalEvents += 1;

  // get all of the true particles from G4
  std::vector<caf::SRTrueParticle> true_particles;
  art::Handle<std::vector<simb::MCParticle>> mc_particles;
  GetByLabelStrict(evt, "largeant", mc_particles);

  art::Handle<std::vector<simb::MCTruth>> neutrinos;
  GetByLabelStrict(evt, "generator", neutrinos);

  if (mc_particles.isValid()) {
    for (const simb::MCParticle part: *mc_particles) {
      true_particles.emplace_back();

      ParticleData pdata;
      pdata.AV = &fActiveVolumes;
      pdata.TPCVolumes = &fTPCVolumes;
      art::fill_ptr_vector(pdata.neutrinos, neutrinos);

      art::ServiceHandle<cheat::BackTrackerService> bt_serv;
      pdata.backtracker = bt_serv.get();
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      pdata.inventory_service = pi_serv.get();

      // TODO: fix weird crashes
      FillTrueG4Particle(part, pdata, true_particles.back());
    }
  }

  // keep track of ID offsets
  unsigned pfp_index_offset = 0;
  unsigned slice_index_offset = 0;

  std::vector<std::string> slice_tag_suffixes;
  fParams.SliceTagSuffixes(slice_tag_suffixes);
  if (slice_tag_suffixes.size() == 0) slice_tag_suffixes.push_back("");

  // Loop over TPC reco tags
  for (const std::string &slice_tag_suffix: slice_tag_suffixes) {
    std::cout << "Getting slices with label: " << (fParams.ClusterLabel() + slice_tag_suffix) << std::endl;

    // Get a handle on the slices
    art::Handle<std::vector<recob::Slice>> slices;
    GetByLabelStrict(evt, fParams.ClusterLabel() + slice_tag_suffix, slices);
    // TO DO - Create 3 labels (SliceLabel a,b,&c = pandora, 
    //         pandorawtvCryo01, etc.) Then add them to the same 
    //         vector with the associated CRYO var for the tree

    std::cout << "Got N slices: " << slices->size() << std::endl;

    // Get tracks & showers here
    art::FindManyP<recob::PFParticle> fmPFPart = 
      FindManyPStrict<recob::PFParticle>(slices, evt, "pandora" + slice_tag_suffix);//fParams.PFPartLabel());

    art::Handle<std::vector<recob::PFParticle>> pfparticles;
    GetByLabelStrict(evt, "pandora" + slice_tag_suffix, pfparticles);

    art::FindManyP<anab::T0> fmT0 =
      FindManyPStrict<anab::T0>(pfparticles, evt, "fmatch" + slice_tag_suffix);

    art::Handle<std::vector<recob::Track>> tracks;
    GetByLabelStrict(evt, "pandoraTrack" + slice_tag_suffix, tracks);

    art::FindManyP<larpandoraobj::PFParticleMetadata> fmPFPMeta =
      FindManyPStrict<larpandoraobj::PFParticleMetadata>(pfparticles, evt, "pandora" + slice_tag_suffix);

    art::FindManyP<recob::Track> fmTrack = 
      FindManyPStrict<recob::Track>(pfparticles, evt, "pandoraTrack" + slice_tag_suffix);

    art::FindManyP<recob::Shower> fmShower =
      FindManyPStrict<recob::Shower>(pfparticles, evt, "pandoraShower" + slice_tag_suffix);

    art::FindManyP<anab::Calorimetry> fmCalo = 
      FindManyPStrict<anab::Calorimetry>(tracks, evt, "pandoraCalo" + slice_tag_suffix);

    art::FindManyP<anab::ParticleID> fmPID = 
      FindManyPStrict<anab::ParticleID>(tracks, evt, "pandoraPid" + slice_tag_suffix);

    art::FindManyP<recob::Hit> fmHit =
      FindManyPStrict<recob::Hit>(tracks, evt, "pandoraTrack" + slice_tag_suffix);

    //#######################################################
    // Loop over slices 
    //#######################################################

    for (unsigned int sliceId = 0; sliceId < slices->size(); ++sliceId) {
      const recob::Slice& slice = (*slices)[sliceId];

      //    if (slice.IsNoise() || slice.NCell() == 0) continue;
      // Because we don't care about the noise slice and slices with no hits.
      StandardRecord rec;
      StandardRecord* prec = &rec;  // TTree wants a pointer-to-pointer
      fRecTree->SetBranchAddress("rec", &prec);

      // fill up the true particles
      rec.true_particles = true_particles;

      //#######################################################
      // Fill slice header.
      //#######################################################
      // Get metadata information for header
      unsigned int run = evt.run();
      unsigned int subrun = evt.subRun();
      //   unsigned int spillNum = evt.id().event();

      rec.hdr = SRHeader();

      rec.hdr.run = run;
      rec.hdr.subrun = subrun;
      rec.hdr.subevt = sliceId + slice_index_offset;
      rec.hdr.ismc = !fIsRealData;
      rec.hdr.det = fDet;
      // rec.hdr.cycle = fCycle;
      // rec.hdr.batch = fBatch;
      // rec.hdr.blind = 0;
      // rec.hdr.filt = rb::IsFiltered(evt, slices, sliceId);

      // get the primary particle associated with this slice
      std::vector<art::Ptr<recob::PFParticle>> slcPFParts = fmPFPart.at(sliceId);
      size_t iPart;
      for (iPart = 0; iPart < slcPFParts.size(); ++iPart ) {
        const recob::PFParticle &thisParticle = *slcPFParts[iPart];
        if (thisParticle.IsPrimary()) break;
      }
      // primary particle and meta-data
      SliceData sdata;
      sdata.particle_id_offset = pfp_index_offset;
      sdata.slice_id_offset = slice_index_offset;
      sdata.primary = (iPart == slcPFParts.size()) ? NULL : slcPFParts[iPart].get();
      sdata.primary_meta = (iPart == slcPFParts.size()) ? NULL : fmPFPMeta.at(slcPFParts[iPart]->Self()).at(0).get();
      // get the flash match
      sdata.fmatch = NULL; 
      if (fmT0.isValid() && sdata.primary != NULL) {
        std::vector<art::Ptr<anab::T0>> fmatch = fmT0.at(sdata.primary->Self());
        if (fmatch.size() != 0) {
          assert(fmatch.size() == 1);
          sdata.fmatch = fmatch[0].get(); 
        }
      }

      //#######################################################
      // Add slice info.
      //#######################################################
      FillSliceVars(slice, sdata, rec.slc);

      // select slice
      if (!SelectSlice(rec.slc, fParams.CutClearCosmic())) continue;

      //#######################################################
      // Add detector dependent slice info.
      //#######################################################
      // if (fDet == kSBND) {
      //   rec.sel.contain.nplanestofront = rec.slc.firstplane - (plnfirst - 1);
      //   rec.sel.contain.nplanestoback = (plnlast) - 1 - rec.slc.lastplane;
      // }

      //#######################################################
      // Add reconstructed objects.
      //#######################################################
      // Reco objects have assns to the slice PFParticles 
      // This depends on the findMany object created above.

      if ( fmPFPart.isValid() ) {
        for ( size_t iPart = 0; iPart < slcPFParts.size(); ++iPart ) {
          const recob::PFParticle &thisParticle = *slcPFParts[iPart];

          const std::vector<art::Ptr<recob::Track>> &thisTrack = fmTrack.at(thisParticle.Self());
          const std::vector<art::Ptr<recob::Shower>> &thisShower = fmShower.at(thisParticle.Self());

          if (thisTrack.size())  { // it's a track!
            assert(thisTrack.size() == 1);
            assert(thisShower.size() == 0);
            rec.reco.ntrk ++;
	    rec.reco.trk.push_back(SRTrack()); 

            // setup the data the track will need
            TrackData tdata;

            // set the ID
            // use the PFParticle ID for the id so that it is
            // global across tracks and showers
            // also include the index offset in the case of multiple
            // reconstruction being run (as in ICARUS)
            tdata.particle_index_offset = pfp_index_offset;

            // set the tdata
            tdata.particleIDs = fmPID.at(thisTrack[0]->ID());
            tdata.calos = fmCalo.at(thisTrack[0]->ID());
            tdata.hits = fmHit.at(thisTrack[0]->ID());

            // set the algorithms
            tdata.mcs_calculator = &fMCSCalculator;
            tdata.range_calculator = &fRangeCalculator;
            tdata.geom = lar::providerFrom<geo::Geometry>();
              
	    FillTrackVars(*thisTrack[0], thisParticle, tdata, rec.reco.trk.back());
	    
	  } // thisTrack exists
          else if (thisShower.size()) {
            assert(thisTrack.size() == 0);
            assert(thisShower.size() == 1);
            // TODO: fill shower vars
          } // thisShower exists
          else {}
	
        }// end for pfparts
      } // fmPFPart ok
    


      //#######################################################
      // Fill truth information
      //#######################################################

      // // Set mc branch values to default
      // rec.mc.setDefault();
      // if (fParams.EnableBlindness()) BlindThisRecord(&rec);

      fRecTree->Fill();
      srcol->push_back(rec);
      //util::CreateAssn(*this, evt, *srcol, art::Ptr<recob::Slice>(slices, sliceId),
      //                 *srAssn);
    }  // end loop over slices
    // update id offsets
    pfp_index_offset += pfparticles->size();
    slice_index_offset += slices->size();
  } // end loop over TPC-reco tags

  evt.put(std::move(srcol));
}

void CAFMaker::endSubRun(art::SubRun& sr) {

}

//......................................................................
void CAFMaker::endJob() {
  if (fTotalEvents == 0) {

    std::cerr << "No events processed in this file. Aborting rather than "
                 "produce an empty CAF."
              << std::endl;
    abort();
  }

  // Make sure the recTree is in the file before filling other items 
  // for debugging. 
  fFile->Write();

  fFile->cd();
  hEvents->Fill(.5, fTotalEvents);

  hEvents->Write();
  fFile->Write();

}


}  // end namespace caf
DEFINE_ART_MODULE(caf::CAFMaker)
////////////////////////////////////////////////////////////////////////
