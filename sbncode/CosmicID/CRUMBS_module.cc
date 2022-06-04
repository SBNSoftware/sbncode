////////////////////////////////////////////////////////////////////////
// Class:       CRUMBS
// Plugin Type: producer
// File:        CRUMBS_module.cc
//
// Generated at Wed Jan  5 08:25:29 2022 by Henry Lay using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art_root_io/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include "sbncode/GeometryTools/TPCGeoAlg.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"
#include "sbncode/LArRecoProducer/TrackStoppingChi2Alg.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"

#include "TTree.h"
#include "TMVA/Reader.h"

class CRUMBS;

namespace sbn {
  class CRUMBS : public art::EDProducer {
  public:
    explicit CRUMBS(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRUMBS(CRUMBS const&) = delete;
    CRUMBS(CRUMBS&&) = delete;
    CRUMBS& operator=(CRUMBS const&) = delete;
    CRUMBS& operator=(CRUMBS&&) = delete;

    // Required functions.
    void produce(art::Event& e) override;

    void ResetVars();
    void GetMaps(art::Event const& e, std::map<int, int> &trackIDToGenMap, std::map<int, std::string> &genTypeMap);

    art::Ptr<recob::PFParticle> GetSlicePrimary(art::Event const& e, 
                                                const art::Ptr<recob::Slice> &slice, 
                                                const art::ValidHandle<std::vector<recob::Slice> > &handleSlices);

    std::vector<art::Ptr<anab::T0> > GetCRTTrackT0s(art::Event const& e, const art::Ptr<recob::Slice> &slice, 
                                                    const art::ValidHandle<std::vector<recob::PFParticle> > &handlePFPs,
                                                    const art::ValidHandle<std::vector<recob::Slice> > &handleSlices);

    std::vector<art::Ptr<anab::T0> > GetCRTHitT0s(art::Event const& e, const art::Ptr<recob::Slice> &slice, 
                                                  const art::ValidHandle<std::vector<recob::PFParticle> > &handlePFPs,
                                                  const art::ValidHandle<std::vector<recob::Slice> > &handleSlices);

    float GetLongestTrackStoppingChi2Ratio(art::Event const& e, const art::Ptr<recob::Slice> &slice, 
                                           const art::ValidHandle<std::vector<recob::PFParticle> > &handlePFPs,
                                           const art::ValidHandle<std::vector<recob::Slice> > &handleSlices);

    void FillCRTVars(const std::vector<art::Ptr<anab::T0> > &trackT0s, const std::vector<art::Ptr<anab::T0> > &hitT0s);

    void FillPandoraNuScoreVars(std::map<std::string, float> &propertiesMap);

    std::vector<art::Ptr<recob::Hit> > GetAllSliceHits(art::Event const& e, 
                                                       const art::Ptr<recob::Slice> &slice, 
                                                       const art::ValidHandle<std::vector<recob::Slice> > &handleSlices);

    void GetTruthMatching(art::Event const& e, const std::vector<art::Ptr<recob::Hit> > &sliceHits, const std::vector<art::Ptr<recob::Hit> > &allHits, 
                          std::map<int, int> &trackIDToGenMap, int &matchedID, double &purity, double &completeness);

    int SliceTruthId(std::map<int, float> &purities);

  private:

    // Bools to control training
    bool fTrainingMode, fProcessNeutrinos, fProcessCosmics;

    // Module labels
    std::string fMCParticleModuleLabel, fGeneratorModuleLabel, fCosmicModuleLabel, fPFParticleModuleLabel, fHitModuleLabel, fTrackModuleLabel, fSliceModuleLabel, 
      fFlashMatchModuleLabel, fCRTTrackMatchModuleLabel, fCRTHitMatchModuleLabel, fCalorimetryModuleLabel;

    // MVA location and type for loading
    std::string fMVAName, fMVAFileName;

    // Parameter set to pass to the stopping chi2 alg
    fhicl::ParameterSet fChi2FitParams;

    // Tree for storing training information
    TTree *fSliceTree;
    
    // TMVA reader for calculating CRUMBS score
    TMVA::Reader *fMVAReader;

    // Other useful information for training tree
    float tpc_NuScore;
    unsigned eventID, subRunID, runID, slicePDG;
    std::string matchedType;
    double matchedPurity, matchedCompleteness;

    // Algorithms used for calculating variables
    sbn::TrackStoppingChi2Alg fTrackStoppingChi2Alg;
    sbn::TPCGeoAlg fTpcGeo;

    // ======================== //
    //  CRUMBS INPUT VARIABLES  //
    // ======================== //

    // Pandora Cosmic Hypothesis Variables
    float tpc_CRFracHitsInLongestTrack;   // fraction of slice’s space points in longest track
    float tpc_CRLongestTrackDeflection;   // 1 - the cosine of the angle between the starting and finishing directions of the longest track
    float tpc_CRLongestTrackDirY;         // relative direction of the longest track in Y
    float tpc_CRNHitsMax;                 // the number of space points in the largest pfp

    // Pandora Neutrino Hypothesis Variables
    float tpc_NuEigenRatioInSphere;       // the ratio between the first and second eigenvalues from a PCA of spacepoints within 10cm of the vertex
    float tpc_NuNFinalStatePfos;          // the number of final state pfos
    float tpc_NuNHitsTotal;               // the total number of space points
    float tpc_NuNSpacePointsInSphere;     // the total number of space points within 10cm of the vertex
    float tpc_NuVertexY;                  // the vertex position in Y [cm]
    float tpc_NuWeightedDirZ;             // the Z component of the space-point weighted direction of the final state pfos

    // Other TPC Variables
    float tpc_StoppingChi2CosmicRatio;    // a ratio of chi2 values intended to find Bragg peaks in stopping muon tracks

    // SBN Simple Flash Match Variables
    float pds_FMTotalScore;               // the total score
    float pds_FMPE;                       // the total number of photoelectrons in the associated flash
    float pds_FMTime;                     // the time associated with the flash [us]

    // CRT Track and Hit Matching Variables
    float crt_TrackScore;                // a combination of the DCA and angle between the best matched TPC & CRT tracks
    float crt_HitScore;                  // the best distance from an extrapolated TPC track to a CRT hit [cm]
    float crt_TrackTime;                 // the time associated with the matched CRT track [us]
    float crt_HitTime;                   // the time associated with the matched CRT hit [us]
  };


  CRUMBS::CRUMBS(fhicl::ParameterSet const& p)
    : EDProducer{p},
    fTrainingMode                 (p.get<bool>("TrainingMode",false)),
    fProcessNeutrinos             (p.get<bool>("ProcessNeutrinos",true)),
    fProcessCosmics               (p.get<bool>("ProcessCosmics",true)),
    fMCParticleModuleLabel        (p.get<std::string>("MCParticleModuleLabel","")),
    fGeneratorModuleLabel         (p.get<std::string>("GeneratorModuleLabel","")),
    fCosmicModuleLabel            (p.get<std::string>("CosmicModuleLabel","")),
    fPFParticleModuleLabel        (p.get<std::string>("PFParticleModuleLabel")),
    fHitModuleLabel               (p.get<std::string>("HitModuleLabel")),
    fTrackModuleLabel             (p.get<std::string>("TrackModuleLabel")),
    fSliceModuleLabel             (p.get<std::string>("SliceModuleLabel")),
    fFlashMatchModuleLabel        (p.get<std::string>("FlashMatchModuleLabel")),
    fCRTTrackMatchModuleLabel     (p.get<std::string>("CRTTrackMatchModuleLabel")),
    fCRTHitMatchModuleLabel       (p.get<std::string>("CRTHitMatchModuleLabel")),
    fCalorimetryModuleLabel       (p.get<std::string>("CalorimetryModuleLabel")),
    fMVAName                      (p.get<std::string>("MVAName")),
    fMVAFileName                  (p.get<std::string>("MVAFileName")),
    fChi2FitParams                (p.get<fhicl::ParameterSet>("Chi2FitParams")),
    fTrackStoppingChi2Alg(fChi2FitParams)
    {
      produces<std::vector<CRUMBSResult>>();
      produces<art::Assns<recob::Slice, CRUMBSResult>>();

      fMVAReader = new TMVA::Reader("V");

      fMVAReader->AddVariable("tpc_CRFracHitsInLongestTrack",&tpc_CRFracHitsInLongestTrack);
      fMVAReader->AddVariable("tpc_CRLongestTrackDeflection",&tpc_CRLongestTrackDeflection);
      fMVAReader->AddVariable("tpc_CRLongestTrackDirY",&tpc_CRLongestTrackDirY);
      fMVAReader->AddVariable("tpc_CRNHitsMax",&tpc_CRNHitsMax);
      fMVAReader->AddVariable("tpc_NuEigenRatioInSphere",&tpc_NuEigenRatioInSphere);
      fMVAReader->AddVariable("tpc_NuNFinalStatePfos",&tpc_NuNFinalStatePfos);
      fMVAReader->AddVariable("tpc_NuNHitsTotal",&tpc_NuNHitsTotal);
      fMVAReader->AddVariable("tpc_NuNSpacePointsInSphere",&tpc_NuNSpacePointsInSphere);
      fMVAReader->AddVariable("tpc_NuVertexY",&tpc_NuVertexY);
      fMVAReader->AddVariable("tpc_NuWeightedDirZ",&tpc_NuWeightedDirZ);
      fMVAReader->AddVariable("tpc_StoppingChi2CosmicRatio",&tpc_StoppingChi2CosmicRatio);

      fMVAReader->AddVariable("pds_FMTotalScore",&pds_FMTotalScore);
      fMVAReader->AddVariable("pds_FMPE",&pds_FMPE);
      fMVAReader->AddVariable("pds_FMTime",&pds_FMTime);

      fMVAReader->AddVariable("crt_TrackScore",&crt_TrackScore);
      fMVAReader->AddVariable("crt_HitScore",&crt_HitScore);
      fMVAReader->AddVariable("crt_TrackTime",&crt_TrackTime);
      fMVAReader->AddVariable("crt_HitTime",&crt_HitTime);

      cet::search_path searchPath("FW_SEARCH_PATH");
      std::string weightFileFullPath;
      if (!searchPath.find_file(fMVAFileName, weightFileFullPath))
        throw cet::exception("CRUMBS") << "Unable to find weight file: " << fMVAFileName << " in FW_SEARCH_PATH: " << searchPath.to_string();

      fMVAReader->BookMVA(fMVAName, weightFileFullPath);

      art::ServiceHandle<art::TFileService> tfs;
      if(fTrainingMode)
        {
          fSliceTree = tfs->make<TTree>("SliceTree","Slice data TTree");

          fSliceTree->Branch("tpc_NuScore",&tpc_NuScore);
          fSliceTree->Branch("tpc_CRFracHitsInLongestTrack",&tpc_CRFracHitsInLongestTrack);
          fSliceTree->Branch("tpc_CRLongestTrackDeflection",&tpc_CRLongestTrackDeflection);
          fSliceTree->Branch("tpc_CRLongestTrackDirY",&tpc_CRLongestTrackDirY);
          fSliceTree->Branch("tpc_CRNHitsMax",&tpc_CRNHitsMax);
          fSliceTree->Branch("tpc_NuEigenRatioInSphere",&tpc_NuEigenRatioInSphere);
          fSliceTree->Branch("tpc_NuNFinalStatePfos",&tpc_NuNFinalStatePfos);
          fSliceTree->Branch("tpc_NuNHitsTotal",&tpc_NuNHitsTotal);
          fSliceTree->Branch("tpc_NuNSpacePointsInSphere",&tpc_NuNSpacePointsInSphere);
          fSliceTree->Branch("tpc_NuVertexY",&tpc_NuVertexY);
          fSliceTree->Branch("tpc_NuWeightedDirZ",&tpc_NuWeightedDirZ);
          fSliceTree->Branch("tpc_StoppingChi2CosmicRatio",&tpc_StoppingChi2CosmicRatio);

          fSliceTree->Branch("pds_FMTotalScore",&pds_FMTotalScore);
          fSliceTree->Branch("pds_FMPE",&pds_FMPE);
          fSliceTree->Branch("pds_FMTime",&pds_FMTime);

          fSliceTree->Branch("crt_TrackScore",&crt_TrackScore);
          fSliceTree->Branch("crt_HitScore",&crt_HitScore);
          fSliceTree->Branch("crt_TrackTime",&crt_TrackTime);
          fSliceTree->Branch("crt_HitTime",&crt_HitTime);

          fSliceTree->Branch("eventID",&eventID);
          fSliceTree->Branch("subRunID",&subRunID);
          fSliceTree->Branch("runID",&runID);
          fSliceTree->Branch("slicePDG",&slicePDG);
          fSliceTree->Branch("matchedType",&matchedType);
          fSliceTree->Branch("matchedPurity",&matchedPurity);
          fSliceTree->Branch("matchedCompleteness",&matchedCompleteness);
        }
    }

  void CRUMBS::ResetVars()
  {
    tpc_NuScore = -999999.; tpc_CRFracHitsInLongestTrack = -999999.; tpc_CRLongestTrackDeflection = -999999.; tpc_CRLongestTrackDirY = -999999.; tpc_CRNHitsMax = -999999.;
    tpc_NuEigenRatioInSphere = -999999.; tpc_NuNFinalStatePfos = -999999.; tpc_NuNHitsTotal = -999999.; tpc_NuNSpacePointsInSphere = -999999.; tpc_NuVertexY = -999999.;
    tpc_NuWeightedDirZ = -999999.; tpc_StoppingChi2CosmicRatio = -4.;

    pds_FMTotalScore = -999999.; pds_FMPE = -999999.; pds_FMTime = -500.;

    crt_TrackScore = -4.; crt_HitScore = -4.; crt_TrackTime = -3000; crt_HitTime = -3000;

    slicePDG = 999999;
    matchedType = "";
    matchedPurity = -999999.; matchedCompleteness = -999999.;
  }

  void CRUMBS::GetMaps(art::Event const& e, std::map<int, int> &trackIDToGenMap, std::map<int, std::string> &genTypeMap)
  {

    unsigned nNu(0), nCos(0);

    if(fProcessNeutrinos)
      {
        art::Handle<std::vector<simb::MCTruth> > handleMCTruthNu;
        e.getByLabel(fGeneratorModuleLabel, handleMCTruthNu);
        art::FindManyP<simb::MCParticle> truthNuMCPAssn(handleMCTruthNu,e,fMCParticleModuleLabel);

        for (unsigned int i = 0; i < handleMCTruthNu->size(); ++i){
          const art::Ptr<simb::MCTruth> mcTruth(handleMCTruthNu, i);
          const simb::MCParticle nu = mcTruth->GetNeutrino().Nu();

          if(!fTpcGeo.InVolume(nu))
            genTypeMap[i] = "DirtNu";
          else
            genTypeMap[i] = "Nu";
    
          const std::vector<art::Ptr<simb::MCParticle> > particles = truthNuMCPAssn.at(mcTruth.key());
    
          for (auto const& particle : particles)
            {
              trackIDToGenMap[particle->TrackId()] = i;
            }
          ++nNu;
        }
      }

    if(fProcessCosmics)
      {
        art::Handle<std::vector<simb::MCTruth> > handleMCTruthCosmic;
        e.getByLabel(fCosmicModuleLabel, handleMCTruthCosmic);

        art::FindManyP<simb::MCParticle> truthCosmicMCPAssn(handleMCTruthCosmic,e,fMCParticleModuleLabel);

        for (unsigned int i = 0; i < handleMCTruthCosmic->size(); ++i){
          const art::Ptr<simb::MCTruth> mcTruth(handleMCTruthCosmic, i);

          genTypeMap[i + nNu] = "Cosmic";
    
          const std::vector<art::Ptr<simb::MCParticle> > particles = truthCosmicMCPAssn.at(mcTruth.key());
    
          for (auto const& particle : particles)
            {
              trackIDToGenMap[particle->TrackId()] = i + nNu;
            }
          ++nCos;
        }
      }
  
    eventID = e.event();
    subRunID = e.subRun();
    runID = e.run();
  }

  void CRUMBS::produce(art::Event& e)
  {
    std::map<int, int> trackIDToGenMap;
    std::map<int, std::string> genTypeMap;

    if(fTrainingMode)
      this->GetMaps(e, trackIDToGenMap, genTypeMap);

    auto resultsVec = std::make_unique<std::vector<CRUMBSResult>>();
    auto sliceAssns = std::make_unique<art::Assns<recob::Slice, CRUMBSResult>>();

    auto const handleSlices(e.getValidHandle<std::vector<recob::Slice>>(fSliceModuleLabel));
    std::vector<art::Ptr<recob::Slice>> slices;
    art::fill_ptr_vector(slices, handleSlices);

    auto const handlePFPs(e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleModuleLabel));
    std::vector<art::Ptr<recob::PFParticle>> pfps;
    art::fill_ptr_vector(pfps, handlePFPs);

    auto const handleHits(e.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel));
    std::vector<art::Ptr<recob::Hit>> allHits;
    art::fill_ptr_vector(allHits, handleHits);

    art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssoc(handlePFPs, e, fPFParticleModuleLabel);
    art::FindManyP<sbn::SimpleFlashMatch> pfpFMAssoc(handlePFPs, e, fFlashMatchModuleLabel);

    for(auto const &slice : slices)
      {
        this->ResetVars();

        auto const primary = this->GetSlicePrimary(e, slice, handleSlices);

        if(primary.isNull())
          continue;

        if(primary->PdgCode() == 13 || primary->PdgCode() == 11)
          continue;

        const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = pfpMetadataAssoc.at(primary.key());
        const std::vector<art::Ptr<sbn::SimpleFlashMatch> > pfpFMVec = pfpFMAssoc.at(primary.key());
        const std::vector<art::Ptr<anab::T0> > sliceCRTTrackT0s = this->GetCRTTrackT0s(e, slice, handlePFPs, handleSlices);
        const std::vector<art::Ptr<anab::T0> > sliceCRTHitT0s = this->GetCRTHitT0s(e, slice, handlePFPs, handleSlices);

        this->FillCRTVars(sliceCRTTrackT0s, sliceCRTHitT0s);

        const art::Ptr<larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetaVec.front();
        std::map<std::string, float> propertiesMap = pfpMeta->GetPropertiesMap();
      
        this->FillPandoraNuScoreVars(propertiesMap);

        tpc_StoppingChi2CosmicRatio = this->GetLongestTrackStoppingChi2Ratio(e, slice, handlePFPs, handleSlices);
    
        const art::Ptr<sbn::SimpleFlashMatch> flashmatch = pfpFMVec.front();
        pds_FMTotalScore = flashmatch->score.total;
        pds_FMPE = flashmatch->light.pe;
        pds_FMTime = std::max(flashmatch->time, -100.);
      
        const float score = fMVAReader->EvaluateMVA(fMVAName);

        resultsVec->emplace_back(score, tpc_CRFracHitsInLongestTrack, tpc_CRLongestTrackDeflection, tpc_CRLongestTrackDirY, std::round(tpc_CRNHitsMax),
                                 tpc_NuEigenRatioInSphere, std::round(tpc_NuNFinalStatePfos), std::round(tpc_NuNHitsTotal), std::round(tpc_NuNSpacePointsInSphere), 
                                 tpc_NuVertexY, tpc_NuWeightedDirZ, tpc_StoppingChi2CosmicRatio, pds_FMTotalScore, pds_FMPE, pds_FMTime, crt_TrackScore, crt_HitScore, 
                                 crt_TrackTime, crt_HitTime);
        
        util::CreateAssn(*this, e, *resultsVec, slice, *sliceAssns);

        if(fTrainingMode)
          {
            std::vector<art::Ptr<recob::Hit> > sliceHits = this->GetAllSliceHits(e, slice, handleSlices);

            int matchedID(-1);
            this->GetTruthMatching(e, sliceHits, allHits, trackIDToGenMap, matchedID, matchedPurity, matchedCompleteness);

            slicePDG = primary->PdgCode();
            matchedType = genTypeMap[matchedID];
      
            fSliceTree->Fill();
          }
      }

    e.put(std::move(resultsVec));
    e.put(std::move(sliceAssns));
  }

  void CRUMBS::FillCRTVars(const std::vector<art::Ptr<anab::T0> > &trackT0s, const std::vector<art::Ptr<anab::T0> > &hitT0s)
  {
    if (!trackT0s.empty()){
      crt_TrackScore = std::numeric_limits<float>::max();
      for(auto const crttrackmatcht0 : trackT0s)
        {
          if(crttrackmatcht0->TriggerConfidence() < crt_TrackScore)
            {
              crt_TrackScore = crttrackmatcht0->TriggerConfidence();
              crt_TrackTime = crttrackmatcht0->Time() * 1e-3;
            }
        }
    }
  
    if (!hitT0s.empty()){
      crt_HitScore = std::numeric_limits<float>::max();
      for(auto const crthitmatcht0 : hitT0s)
        {
          if(crthitmatcht0->TriggerConfidence() < crt_HitScore)
            {
              crt_HitScore = crthitmatcht0->TriggerConfidence();
              crt_HitTime = crthitmatcht0->Time() * 1e-3;
            }
        }
    }
  }

  void CRUMBS::FillPandoraNuScoreVars(std::map<std::string, float> &propertiesMap)
  {
    auto propertiesMapIter = propertiesMap.find("NuScore");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- NuScore" << std::endl;
      abort();
    }
    tpc_NuScore = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("CRFracHitsInLongestTrack");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- CRFracHitsInLongestTrack" << std::endl;
      abort();
    }
    tpc_CRFracHitsInLongestTrack = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("CRLongestTrackDeflection");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- CRLongestTrackDeflection" << std::endl;
      abort();
    }
    tpc_CRLongestTrackDeflection = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("CRLongestTrackDirY");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- CRLongestTrackDirY" << std::endl;
      abort();
    }
    tpc_CRLongestTrackDirY = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("CRNHitsMax");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- CRNHitsMax" << std::endl;
      abort();
    }
    tpc_CRNHitsMax = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("NuEigenRatioInSphere");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- NuEigenRatioInSphere" << std::endl;
      abort();
    }
    tpc_NuEigenRatioInSphere = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("NuNFinalStatePfos");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- NuNFinalStatePfos" << std::endl;
      abort();
    }
    tpc_NuNFinalStatePfos = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("NuNHitsTotal");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- NuNHitsTotal" << std::endl;
      abort();
    }
    tpc_NuNHitsTotal = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("NuNSpacePointsInSphere");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- NuNSpacePointsInSphere" << std::endl;
      abort();
    }
    tpc_NuNSpacePointsInSphere = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("NuVertexY");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- NuVertexY" << std::endl;
      abort();
    }
    tpc_NuVertexY = propertiesMapIter->second;

    propertiesMapIter = propertiesMap.find("NuWeightedDirZ");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout << "CRUMBS_module: Error finding variable -- NuWeightedDirZ" << std::endl;
      abort();
    }
    tpc_NuWeightedDirZ = propertiesMapIter->second;
  }

  std::vector<art::Ptr<recob::Hit> > CRUMBS::GetAllSliceHits(art::Event const& e, const art::Ptr<recob::Slice> &slice, const art::ValidHandle<std::vector<recob::Slice> > &handleSlices)
  {
    art::FindManyP<recob::Hit> sliceHitAssn(handleSlices,e,fSliceModuleLabel);
    return sliceHitAssn.at(slice.key());
  }

  art::Ptr<recob::PFParticle> CRUMBS::GetSlicePrimary(art::Event const& e, const art::Ptr<recob::Slice> &slice, const art::ValidHandle<std::vector<recob::Slice> > &handleSlices)
  {
    art::FindManyP<recob::PFParticle> slicePfpAssn(handleSlices,e,fSliceModuleLabel);
    std::vector<art::Ptr<recob::PFParticle> > pfps = slicePfpAssn.at(slice.key());

    for(auto const &pfp : pfps)
      {
        if(pfp->IsPrimary())
          return pfp;
      }
  
    art::Ptr<recob::PFParticle> nullReturn;
    return nullReturn;
  }

  void CRUMBS::GetTruthMatching(art::Event const& e, const std::vector<art::Ptr<recob::Hit> > &sliceHits, const std::vector<art::Ptr<recob::Hit> > &allHits, 
                                std::map<int, int> &trackIDToGenMap, int &matchedID, double &purity, double &completeness)
  {
    std::map<int, int> sliceHitMap;
    std::map<int, float> slicePurityMap;

    auto clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    for (auto const& hit : sliceHits)
      {
        ++sliceHitMap[trackIDToGenMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]];
      }

    for (auto const& [id, nHits] : sliceHitMap)
      {
        slicePurityMap[id] = (float) nHits / (float) sliceHits.size();
      }

    for (auto const& [id, pur] : slicePurityMap)
      {
        if(pur > purity) 
          {
            matchedID = id;
            purity = pur;
          }
      }

    int totalTrueHits(0);

    for (auto const& hit : allHits)
      {
        if(trackIDToGenMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)] == matchedID)
          ++totalTrueHits;
      }
  
    if(totalTrueHits == 0) 
      completeness = 0;
    else
      completeness = sliceHitMap[matchedID] / (float) totalTrueHits;
  }

  std::vector<art::Ptr<anab::T0> > CRUMBS::GetCRTTrackT0s(art::Event const& e, const art::Ptr<recob::Slice> &slice, const art::ValidHandle<std::vector<recob::PFParticle> > &handlePFPs,
                                                          const art::ValidHandle<std::vector<recob::Slice> > &handleSlices)
  {
    std::vector<art::Ptr<anab::T0> > t0Vec;

    art::Handle<std::vector<recob::Track> > handleTracks;
    e.getByLabel(fTrackModuleLabel, handleTracks);

    art::FindManyP<recob::PFParticle> slicePFPAssn(handleSlices,e,fSliceModuleLabel);
    art::FindManyP<recob::Track> pfpTrackAssn(handlePFPs,e,fTrackModuleLabel);
    art::FindManyP<anab::T0> trackT0Assn(handleTracks,e,fCRTTrackMatchModuleLabel);

    const std::vector<art::Ptr<recob::PFParticle> > pfps = slicePFPAssn.at(slice.key());
  
    for(auto const& pfp : pfps)
      {
        if(pfp->PdgCode() != 13)
          continue;

        const std::vector<art::Ptr<recob::Track> > tracks = pfpTrackAssn.at(pfp.key());

        if(tracks.size() != 1)
          continue;

        const art::Ptr<recob::Track> track = tracks.front();

        const std::vector<art::Ptr<anab::T0> > t0s = trackT0Assn.at(track.key());
        t0Vec.insert(t0Vec.end(), t0s.begin(), t0s.end());
      }
  
    return t0Vec;
  }

  std::vector<art::Ptr<anab::T0> > CRUMBS::GetCRTHitT0s(art::Event const& e, const art::Ptr<recob::Slice> &slice, const art::ValidHandle<std::vector<recob::PFParticle> > &handlePFPs,
                                                        const art::ValidHandle<std::vector<recob::Slice> > &handleSlices)
  {
    std::vector<art::Ptr<anab::T0> > t0Vec;

    art::Handle<std::vector<recob::Track> > handleTracks;
    e.getByLabel(fTrackModuleLabel, handleTracks);

    art::FindManyP<recob::PFParticle> slicePFPAssn(handleSlices,e,fSliceModuleLabel);
    art::FindManyP<recob::Track> pfpTrackAssn(handlePFPs,e,fTrackModuleLabel);
    art::FindManyP<anab::T0> trackT0Assn(handleTracks,e,fCRTHitMatchModuleLabel);

    const std::vector<art::Ptr<recob::PFParticle> > pfps = slicePFPAssn.at(slice.key());
  
    for(auto const& pfp : pfps)
      {
        if(pfp->PdgCode() != 13)
          continue;

        const std::vector<art::Ptr<recob::Track> > tracks = pfpTrackAssn.at(pfp.key());

        if(tracks.size() != 1)
          continue;

        const art::Ptr<recob::Track> track = tracks.front();

        const std::vector<art::Ptr<anab::T0> > t0s = trackT0Assn.at(track.key());
        t0Vec.insert(t0Vec.end(), t0s.begin(), t0s.end());
      }
  
    return t0Vec;
  }

  float CRUMBS::GetLongestTrackStoppingChi2Ratio(art::Event const& e, const art::Ptr<recob::Slice> &slice, const art::ValidHandle<std::vector<recob::PFParticle> > &handlePFPs,
                                                 const art::ValidHandle<std::vector<recob::Slice> > &handleSlices)
  {
    art::Ptr<anab::Calorimetry> longestTrackCalo;
    float maxLength = -std::numeric_limits<float>::max();

    art::Handle<std::vector<recob::Track> > handleTracks;
    e.getByLabel(fTrackModuleLabel, handleTracks);

    art::FindManyP<recob::PFParticle> slicePFPAssn(handleSlices,e,fSliceModuleLabel);
    art::FindOneP<recob::Track> pfpTrackAssn(handlePFPs,e,fTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> trackCaloAssn(handleTracks,e,fCalorimetryModuleLabel);

    const std::vector<art::Ptr<recob::PFParticle> > pfps = slicePFPAssn.at(slice.key());
  
    for(auto const& pfp : pfps)
      {
        if(pfp->PdgCode() != 13)
          continue;

        const art::Ptr<recob::Track> track = pfpTrackAssn.at(pfp.key());

        if(track.isNull())
          continue;

        const std::vector<art::Ptr<anab::Calorimetry> > calos = trackCaloAssn.at(track.key());

        const unsigned int maxHits(std::max({ calos[0]->dEdx().size(), calos[1]->dEdx().size(), calos[2]->dEdx().size() }));
        const int bestPlane((calos[2]->dEdx().size() == maxHits) ? 2 : (calos[0]->dEdx().size() == maxHits) ? 0 : (calos[1]->dEdx().size() == maxHits) ? 1 : -1);

        if (bestPlane == -1)
          continue;

        const art::Ptr<anab::Calorimetry> calo = calos.at(bestPlane);

        if(track->Length() > maxLength)
          {
            maxLength = track->Length();
            longestTrackCalo = calo;
          }
      }

    sbn::StoppingChi2Fit longestTrackFit = longestTrackCalo.isNull() ? StoppingChi2Fit() : fTrackStoppingChi2Alg.RunFitForCosmicID(*longestTrackCalo);

    if(longestTrackFit.pol0Chi2 < 0 || longestTrackFit.expChi2 <= 0)
      return -4.f;
  
    return longestTrackFit.pol0Chi2 / longestTrackFit.expChi2;
  }
}

DEFINE_ART_MODULE(sbn::CRUMBS)
