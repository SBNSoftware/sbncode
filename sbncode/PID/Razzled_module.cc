////////////////////////////////////////////////////////////////////////
// Class:       Razzled
// Plugin Type: producer (art v3_05_01)
// File:        Razzled_module.cc
//
// Created by Henry Lay, May 2023 to combine the work done by Edward
// Tyley in the Razzle and Dazzle modules into a single tool for
// PFP PID regardless of track or shower decision.
//
// Module that attempts to classify each shower as either:
// - Electron (11)
// - Muon (13)
// - Photon (22)
// - Pion (211)
// - Proton (2212)
// - Other (0)
//
// If a PFP is not classified it will be assigned a pdg of -5
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/DataViewImpl.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"

#include "TTree.h"
#include "TMVA/Reader.h"

#include <numeric>

namespace sbn {
  class Razzled : public art::EDProducer {
  public:
    explicit Razzled(fhicl::ParameterSet const& p);
  
    Razzled(Razzled const&) = delete;
    Razzled(Razzled&&) = delete;
    Razzled& operator=(Razzled const&) = delete;
    Razzled& operator=(Razzled&&) = delete;

    void produce(art::Event& e) override;
    void beginJob() override;

  private:
    art::ServiceHandle<art::TFileService> tfs;

    art::InputTag fPFPLabel, fClusterLabel, fTrackLabel, fShowerLabel, fCaloLabel, 
      fMCSLabel, fChi2Label, fRangeLabel, fClosestApproachLabel, fStoppingChi2Label;

    const float fMinTrackLength, fMinShowerEnergy;
    const bool fMakeTree, fRunMVA;
    const std::string fMethodName, fWeightFile;
    const float fXMin, fXMax, fYMin, fYMax, fZMin, fZMax;

    float pfp_numDaughters;    // PFP's number of daughters in particle hierarchy
    float pfp_maxDaughterHits; // Max number of hits from any PFP daughters
    float pfp_trackScore;
    float pfp_chargeEndFrac;
    float pfp_chargeFracSpread;
    float pfp_linearFitDiff;
    float pfp_linearFitLength;
    float pfp_linearFitGapLength;
    float pfp_linearFitRMS;
    float pfp_openAngleDiff;
    float pfp_secondaryPCARatio;
    float pfp_tertiaryPCARatio;
    float pfp_vertexDist;

    float trk_length;                // The length of the track [cm]
    float trk_chi2PIDMuon;           // Chi2 PID score for muon hypothesis
    float trk_chi2PIDProton;         // Chi2 PID score for proton hypothesis
    float trk_chi2PIDMuonPionDiff;   // The muon and pion scores are very correlated, instead take the relative difference
    float trk_mcsScatterMean;        // Mean scattering angle from MCS
    float trk_mcsScatterMaxRatio;    // Scattering angles used in MCS calculation [mrad]
    float trk_momDiff;               // Relative momentum agreement between range and MCS
    float trk_meanDCA;               // Distance of closest approach to interpolated track [cm]
    float trk_stoppingdEdxChi2Ratio; // Ratio of exp/pol0 chi2 fits to end of track
    float trk_chi2Pol0dEdxFit;       // Fitted pol0 to find the dE/dx of the track [MeV/cm]

    float shw_bestdEdx;          // The dE/dx at the start of the shower (in the best plane) [MeV/cm]
    float shw_convGap;           // The gap between the shower start and parent vertex [cm]
    float shw_openAngle;         // Opening Angle of the shower, defined as atan(width/length) [deg]
    float shw_modHitDensity;     // The hit density corrected for the pitch (hits/wire)
    float shw_sqrtEnergyDensity; // Sqrt of the energy divided by the length [cm^-1]


    float electronScore, muonScore, photonScore, pionScore, protonScore, otherScore, bestScore;
    int bestPDG;

    TMVA::Reader* reader;
    TTree* pfpTree;

    int truePdg;
    std::string trueType, trueEndProcess;
    float energyComp, energyPurity;

    bool recoPrimary, trackContained, showerContained;
    float trackStartX, trackStartY, trackStartZ, trackEndX, trackEndY, trackEndZ,
      showerStartX, showerStartY, showerStartZ, showerEndX, showerEndY, showerEndZ;


    void ClearTreeValues();

    void FillTrueParticleMetrics(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit>> &hits);

    void FillPFPMetrics(const art::Ptr<recob::PFParticle> &pfp, const std::map<size_t, art::Ptr<recob::PFParticle>> &pfpMap,
			const art::FindManyP<recob::Cluster> &pfpsToClusters, const art::FindManyP<recob::Hit> &clustersToHits,
			const art::FindOneP<larpandoraobj::PFParticleMetadata> &pfpsToMetadata);

    void FillPandoraTrackIDScoreVars(const std::map<std::string, float> &propertiesMap);

    void FillPandoraTrackIDScoreVar(const std::map<std::string, float> &propertiesMap, float &var, std::string name);

    void FillTrackMetrics(const art::Ptr<recob::Track> &track);

    void FillMCSMetrics(const art::Ptr<recob::MCSFitResult> &mcs);

    void FillMomDiff(const art::Ptr<RangeP> &rangeP, const art::Ptr<recob::MCSFitResult> &mcs);

    void FillChi2PIDMetrics(const art::Ptr<anab::ParticleID> &pid);

    void FillStoppingChi2Metrics(const art::Ptr<StoppingChi2Fit> &stoppingChi2);

    void FillShowerMetrics(const art::Ptr<recob::Shower> &shower, const std::vector<art::Ptr<recob::Hit>> &hits);

    void FillShowerConversionGap(const art::Ptr<recob::PFParticle> &pfp, const std::map<size_t, art::Ptr<recob::PFParticle>> &pfpMap,
				 const art::Ptr<recob::Shower> &shower, const art::FindOneP<recob::Vertex> &pfpsToVertices);

    bool InFV(const TVector3 &pos);

    std::map<size_t, art::Ptr<recob::PFParticle>> GetPFPMap(const std::vector<art::Ptr<recob::PFParticle>> &pfps);

    std::string PDGString(const int pdg);

    MVAPID RunMVA();
  };

  Razzled::Razzled(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fPFPLabel(p.get<std::string>("PFPLabel"))
    , fClusterLabel(p.get<std::string>("ClusterLabel"))
    , fTrackLabel(p.get<std::string>("TrackLabel"))
    , fShowerLabel(p.get<std::string>("ShowerLabel"))
    , fCaloLabel(p.get<std::string>("CaloLabel"))
    , fMCSLabel(p.get<std::string>("MCSLabel"), std::string("muon"))
    , fChi2Label(p.get<std::string>("Chi2Label"))
    , fRangeLabel(p.get<std::string>("RangeLabel"), std::string("muon"))
    , fClosestApproachLabel(p.get<std::string>("ClosestApproachLabel"))
    , fStoppingChi2Label(p.get<std::string>("StoppingChi2Label"))
    , fMinTrackLength(p.get<float>("MinTrackLength"))    
    , fMinShowerEnergy(p.get<float>("MinShowerEnergy"))
    , fMakeTree(p.get<bool>("MakeTree"))
    , fRunMVA(p.get<bool>("RunMVA"))
    , fMethodName(p.get<std::string>("MethodName", ""))
    , fWeightFile(p.get<std::string>("WeightFile", ""))
    , fXMin(p.get<float>("XMin"))
    , fXMax(p.get<float>("XMax"))
    , fYMin(p.get<float>("YMin"))
    , fYMax(p.get<float>("YMax"))
    , fZMin(p.get<float>("ZMin"))
    , fZMax(p.get<float>("ZMax"))
    {
      if(!fMakeTree && !fRunMVA)
        throw cet::exception("Razzled") << "Configured to do nothing";

      if(fRunMVA) 
        {
          if(fMethodName == "" || fWeightFile == "")
            throw cet::exception("Razzled") << "Trying to run MVA with inputs not set: MethodName: " << fMethodName << " and WeightFile: " << fWeightFile;

          cet::search_path searchPath("FW_SEARCH_PATH");
          std::string fWeightFileFullPath;
          if(!searchPath.find_file(fWeightFile, fWeightFileFullPath))
            throw cet::exception("Razzled") << "Unable to find weight file: " << fWeightFile << " in FW_SEARCH_PATH: " << searchPath.to_string();

          reader = new TMVA::Reader("V");

          reader->AddVariable("pfp_numDaughters", &pfp_numDaughters);
          reader->AddVariable("pfp_maxDaughterHits", &pfp_maxDaughterHits);

          reader->AddVariable("trk_length", &trk_length);
          reader->AddVariable("trk_chi2PIDMuon", &trk_chi2PIDMuon);
          reader->AddVariable("trk_chi2PIDProton", &trk_chi2PIDProton);
          reader->AddVariable("trk_chi2PIDMuonPionDiff", &trk_chi2PIDMuonPionDiff);
          reader->AddVariable("trk_mcsScatterMean", &trk_mcsScatterMean);
          reader->AddVariable("trk_mcsScatterMaxRatio", &trk_mcsScatterMaxRatio);
          reader->AddVariable("trk_momDiff", &trk_momDiff);
          reader->AddVariable("trk_meanDCA", &trk_meanDCA);
          reader->AddVariable("trk_stoppingdEdxChi2Ratio", &trk_stoppingdEdxChi2Ratio);
          reader->AddVariable("trk_chi2Pol0dEdxFit", &trk_chi2Pol0dEdxFit);

          reader->AddVariable("shw_bestdEdx", &shw_bestdEdx);
          reader->AddVariable("shw_convGap", &shw_convGap);
          reader->AddVariable("shw_openAngle", &shw_openAngle);
          reader->AddVariable("shw_modHitDensity", &shw_modHitDensity);
          reader->AddVariable("shw_sqrtEnergyDensity", &shw_sqrtEnergyDensity);

          reader->BookMVA(fMethodName, fWeightFileFullPath);
        }

      produces<std::vector<MVAPID>>();
      produces<art::Assns<recob::PFParticle, MVAPID>>();
    }

  void Razzled::beginJob()
  {
    if(fMakeTree)
      {
        pfpTree = tfs->make<TTree>("pfpTree", "Tree filled per PFP with  PID variables");

        if(fRunMVA)
          {
            pfpTree->Branch("electronScore", &electronScore);
            pfpTree->Branch("muonScore", &muonScore);
            pfpTree->Branch("photonScore", &photonScore);
            pfpTree->Branch("pionScore", &pionScore);
            pfpTree->Branch("protonScore", &protonScore);
            pfpTree->Branch("otherScore", &otherScore);
            pfpTree->Branch("bestScore", &bestScore);
            pfpTree->Branch("bestPDG", &bestPDG);
          }

        pfpTree->Branch("truePdg", &truePdg);
        pfpTree->Branch("trueType", &trueType);
        pfpTree->Branch("trueEndProcess", &trueEndProcess);
        pfpTree->Branch("energyComp", &energyComp);
        pfpTree->Branch("energyPurity", &energyPurity);

        pfpTree->Branch("recoPrimary", &recoPrimary);

        pfpTree->Branch("trackStartX", &trackStartX);
        pfpTree->Branch("trackStartY", &trackStartY);
        pfpTree->Branch("trackStartZ", &trackStartZ);
        pfpTree->Branch("trackEndX", &trackEndX);
        pfpTree->Branch("trackEndY", &trackEndY);
        pfpTree->Branch("trackEndZ", &trackEndZ);
        pfpTree->Branch("trackContained", &trackContained);

        pfpTree->Branch("showerStartX", &showerStartX);
        pfpTree->Branch("showerStartY", &showerStartY);
        pfpTree->Branch("showerStartZ", &showerStartZ);
        pfpTree->Branch("showerEndX", &showerEndX);
        pfpTree->Branch("showerEndY", &showerEndY);
        pfpTree->Branch("showerEndZ", &showerEndZ);
        pfpTree->Branch("showerContained", &showerContained);

        pfpTree->Branch("pfp_numDaughters", &pfp_numDaughters);
        pfpTree->Branch("pfp_maxDaughterHits", &pfp_maxDaughterHits);
        pfpTree->Branch("pfp_trackScore", &pfp_trackScore);
        pfpTree->Branch("pfp_chargeEndFrac", &pfp_chargeEndFrac);
        pfpTree->Branch("pfp_chargeFracSpread", &pfp_chargeFracSpread);
        pfpTree->Branch("pfp_linearFitDiff", &pfp_linearFitDiff);
        pfpTree->Branch("pfp_linearFitLength", &pfp_linearFitLength);
        pfpTree->Branch("pfp_linearFitGapLength", &pfp_linearFitGapLength);
        pfpTree->Branch("pfp_linearFitRMS", &pfp_linearFitRMS);
        pfpTree->Branch("pfp_openAngleDiff", &pfp_openAngleDiff);
        pfpTree->Branch("pfp_secondaryPCARatio", &pfp_secondaryPCARatio);
        pfpTree->Branch("pfp_tertiaryPCARatio", &pfp_tertiaryPCARatio);
        pfpTree->Branch("pfp_vertexDist", &pfp_vertexDist);

        pfpTree->Branch("trk_length", &trk_length);
        pfpTree->Branch("trk_chi2PIDMuon", &trk_chi2PIDMuon);
        pfpTree->Branch("trk_chi2PIDProton", &trk_chi2PIDProton);
        pfpTree->Branch("trk_chi2PIDMuonPionDiff", &trk_chi2PIDMuonPionDiff);
        pfpTree->Branch("trk_mcsScatterMean", &trk_mcsScatterMean);
        pfpTree->Branch("trk_mcsScatterMaxRatio", &trk_mcsScatterMaxRatio);
        pfpTree->Branch("trk_meanDCA", &trk_meanDCA);
        pfpTree->Branch("trk_stoppingdEdxChi2Ratio", &trk_stoppingdEdxChi2Ratio);
        pfpTree->Branch("trk_chi2Pol0dEdxFit", &trk_chi2Pol0dEdxFit);
        pfpTree->Branch("trk_momDiff", &trk_momDiff);

        pfpTree->Branch("shw_bestdEdx", &shw_bestdEdx);
        pfpTree->Branch("shw_convGap", &shw_convGap);
        pfpTree->Branch("shw_openAngle", &shw_openAngle);
        pfpTree->Branch("shw_modHitDensity", &shw_modHitDensity);
        pfpTree->Branch("shw_sqrtEnergyDensity", &shw_sqrtEnergyDensity);
      }
  }

  void Razzled::produce(art::Event& e)
  {
    auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));

    auto const pfpHandle(e.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel));
    auto const clusterHandle(e.getValidHandle<std::vector<recob::Cluster>>(fClusterLabel));
    auto const trackHandle(e.getValidHandle<std::vector<recob::Track>>(fTrackLabel));
    auto const showerHandle(e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel));

    std::vector<art::Ptr<recob::PFParticle>> pfps;
    art::fill_ptr_vector(pfps, pfpHandle);

    art::FindOneP<recob::Vertex> pfpsToVertices(pfpHandle, e, fPFPLabel);
    art::FindOneP<larpandoraobj::PFParticleMetadata> pfpsToMetadata(pfpHandle, e, fPFPLabel);
    art::FindManyP<recob::Cluster> pfpsToClusters(pfpHandle, e, fPFPLabel);
    art::FindManyP<recob::Hit> clustersToHits(clusterHandle, e, fPFPLabel);

    art::FindOneP<recob::Track> pfpsToTracks(pfpHandle, e, fTrackLabel);
    art::FindManyP<anab::Calorimetry> tracksToCalos(trackHandle, e, fCaloLabel);
    art::FindOneP<recob::MCSFitResult> tracksToMCSs(trackHandle, e, fMCSLabel);
    art::FindManyP<anab::ParticleID> tracksToChi2s(trackHandle, e, fChi2Label);
    art::FindOneP<RangeP> tracksToRangePs(trackHandle, e, fRangeLabel);
    art::FindOneP<ScatterClosestApproach> tracksToClosestApproaches(trackHandle, e, fClosestApproachLabel);
    art::FindOneP<StoppingChi2Fit> tracksToStoppingChi2s(trackHandle, e, fStoppingChi2Label);

    art::FindOneP<recob::Shower> pfpsToShowers(pfpHandle, e, fShowerLabel);
    art::FindManyP<recob::Hit> showersToHits(showerHandle, e, fShowerLabel);

    auto mvaPIDVec = std::make_unique<std::vector<MVAPID>>();
    auto pfpAssns  = std::make_unique<art::Assns<recob::PFParticle, MVAPID>>();

    const std::map<size_t, art::Ptr<recob::PFParticle>> pfpMap = this->GetPFPMap(pfps);

    for(auto const& pfp : pfps) 
      {
        this->ClearTreeValues();

        this->FillPFPMetrics(pfp, pfpMap, pfpsToClusters, clustersToHits, pfpsToMetadata);

        const art::Ptr<recob::Track> track   = pfpsToTracks.at(pfp.key());
        const art::Ptr<recob::Shower> shower = pfpsToShowers.at(pfp.key());

        if(track.isNull() && shower.isNull())
          continue;

        if(track.isNonnull())
          {
            if(track->Length() < fMinTrackLength)
              continue;

            this->FillTrackMetrics(track);

            const std::vector<art::Ptr<anab::Calorimetry>> calos = tracksToCalos.at(track.key());
            if(calos.size() != 3)
              continue;
          
            const unsigned int maxHits = std::max({ calos[0]->dEdx().size(), calos[1]->dEdx().size(), calos[2]->dEdx().size()});
            const int bestPlane        = (calos[2]->dEdx().size() == maxHits) ? 2 : (calos[0]->dEdx().size() == maxHits) ? 0 : (calos[1]->dEdx().size() == maxHits) ? 1 : -1;
          
            if(bestPlane < 0 || bestPlane > 3)
              throw cet::exception("Dazzle") << "Best plane: " << bestPlane;
          
            const std::vector<art::Ptr<anab::ParticleID>> chi2s = tracksToChi2s.at(track.key());
            if(chi2s.size() == 3)
              this->FillChi2PIDMetrics(chi2s[bestPlane]);

            const art::Ptr<recob::MCSFitResult> mcs = tracksToMCSs.at(track.key());
            if(mcs.isNonnull())
              this->FillMCSMetrics(mcs);
          
            const art::Ptr<RangeP> rangeP = tracksToRangePs.at(track.key());
            if(rangeP.isNonnull() && mcs.isNonnull())
              this->FillMomDiff(rangeP, mcs);

            const art::Ptr<ScatterClosestApproach> closestApproach = tracksToClosestApproaches.at(track.key());
            if(closestApproach.isNonnull())
              trk_meanDCA = closestApproach->mean;

            const art::Ptr<StoppingChi2Fit> stoppingChi2 = tracksToStoppingChi2s.at(track.key());
            if(stoppingChi2.isNonnull())
              this->FillStoppingChi2Metrics(stoppingChi2);
          }
         
        if(shower.isNonnull())
          {
            if(shower->best_plane() < 0 || shower->Energy().at(shower->best_plane()) < fMinShowerEnergy)
              continue;

            const std::vector<art::Ptr<recob::Hit>> showerHits = showersToHits.at(shower.key());

            this->FillShowerMetrics(shower, showerHits);
            this->FillShowerConversionGap(pfp, pfpMap, shower, pfpsToVertices);
          }
    
        if(fRunMVA)
          {
            MVAPID mvaPID = this->RunMVA();
            mvaPIDVec->push_back(mvaPID);
            util::CreateAssn(*this, e, *mvaPIDVec, pfp, *pfpAssns);
          }

        if(fMakeTree)
          {
            std::vector<art::Ptr<recob::Hit>> hits;
          
            const std::vector<art::Ptr<recob::Cluster>> clusters = pfpsToClusters.at(pfp.key());
            for(auto const& cluster : clusters)
              {
                const std::vector<art::Ptr<recob::Hit>> clusterHits = clustersToHits.at(cluster.key());
                hits.insert(hits.end(), clusterHits.begin(), clusterHits.end());
              }

            this->FillTrueParticleMetrics(clockData, hits);
            pfpTree->Fill();
          }
      }

    e.put(std::move(mvaPIDVec));
    e.put(std::move(pfpAssns));
  }

  void Razzled::ClearTreeValues()
  {
    muonScore   = -5.f; electronScore = -5.f; photonScore = -5.f; pionScore = -5.f;
    protonScore = -5.f; otherScore    = -5.f; bestScore   = -5.f;
    bestPDG     = -5;

    truePdg = -5; trueType = ""; trueEndProcess = ""; energyComp = -5.f; energyPurity = -5.f;

    recoPrimary = false;

    trackStartX    = -999.f; trackStartY = -999.f; trackStartZ = -999.f;
    trackEndX      = -999.f; trackEndY   = -999.f; trackEndZ   = -999.f;
    trackContained = false;

    showerStartX    = -999.f; showerStartY = -999.f; showerStartZ = -999.f;
    showerEndX      = -999.f; showerEndY   = -999.f; showerEndZ   = -999.f;
    showerContained = false;

    pfp_numDaughters    = -5.f; pfp_maxDaughterHits    = -5.f; pfp_trackScore       = -5.f;
    pfp_chargeEndFrac   = -5.f; pfp_chargeFracSpread   = -5.f; pfp_linearFitDiff    = -5.f;
    pfp_linearFitLength = -5.f; pfp_linearFitGapLength = -5.f; pfp_linearFitRMS     = -5.f;
    pfp_openAngleDiff   = -5.f; pfp_secondaryPCARatio  = -5.f; pfp_tertiaryPCARatio = -5.f;
    pfp_vertexDist      = -5.f; 

    trk_length              = -5.f; trk_chi2PIDMuon           = -5.f; trk_chi2PIDProton      = -5.f;
    trk_chi2PIDMuonPionDiff = -5.f; trk_mcsScatterMean        = -5.f; trk_mcsScatterMaxRatio = -5.f;
    trk_meanDCA             = -5.f; trk_stoppingdEdxChi2Ratio = -5.f; trk_chi2Pol0dEdxFit    = -5.f;
    trk_momDiff             = -5.f;

    shw_bestdEdx      = -5.f; shw_convGap           = -5.f; shw_openAngle = -5.f;
    shw_modHitDensity = -5.f; shw_sqrtEnergyDensity = -5.f;
  }

  void Razzled::FillTrueParticleMetrics(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit>> &hits)
  {
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

    const int bestMatch = TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, hits, true) ;

    if(!TruthMatchUtils::Valid(bestMatch))
      return;

    float totalHitEnergy = 0.f, totalTrueHitEnergy = 0.f;

    for(auto const& hit : hits)
      {
        const std::vector<sim::TrackIDE> ides = bt_serv->HitToTrackIDEs(clockData, hit);

        totalHitEnergy = std::accumulate(ides.cbegin(), ides.cend(), totalHitEnergy,
                                         [](float sum, auto const& ide) { return sum + ide.energy; });
        totalTrueHitEnergy = std::accumulate(ides.cbegin(), ides.cend(), totalTrueHitEnergy,
                                             [bestMatch](float sum, auto const& ide) { return (std::abs(ide.trackID) == bestMatch) ? sum + ide.energy : sum; });
      }
    
    const std::vector<const sim::IDE*> allIDEs = bt_serv->TrackIdToSimIDEs_Ps(bestMatch);
    
    float totalTrueEnergy = std::accumulate(allIDEs.cbegin(), allIDEs.cend(), 0.f,
                                            [](float sum, auto const& ide) { return sum + ide->energy; });

    energyComp   = totalTrueHitEnergy / totalTrueEnergy;
    energyPurity = totalTrueHitEnergy / totalHitEnergy;

    const simb::MCParticle* trueParticle = particleInv->TrackIdToParticle_P(bestMatch);

    if(trueParticle == NULL)
      return;
    
    truePdg        = trueParticle->PdgCode();
    trueType       = this->PDGString(truePdg);
    trueEndProcess = trueParticle->EndProcess();
  }

  void Razzled::FillPFPMetrics(const art::Ptr<recob::PFParticle> &pfp, const std::map<size_t, art::Ptr<recob::PFParticle>> &pfpMap,
                               const art::FindManyP<recob::Cluster> &pfpsToClusters, const art::FindManyP<recob::Hit> &clustersToHits,
                               const art::FindOneP<larpandoraobj::PFParticleMetadata> &pfpsToMetadata)
  {
    pfp_numDaughters = pfp->Daughters().size();

    for(auto const& daughterID : pfp->Daughters())
      {
        auto const &daughterIter = pfpMap.find(daughterID);
        if(daughterIter == pfpMap.end())
          continue;

        int hits = 0;

        const std::vector<art::Ptr<recob::Cluster>> clusters = pfpsToClusters.at(daughterIter->second.key());
        for(auto const& cluster : clusters)
          {
            const std::vector<art::Ptr<recob::Hit>> clusterHits = clustersToHits.at(cluster.key());
            hits += clusterHits.size();
          }

        if(hits > pfp_maxDaughterHits)
          pfp_maxDaughterHits = hits;
      }

    const art::Ptr<larpandoraobj::PFParticleMetadata> metadata(pfpsToMetadata.at(pfp.key()));

    if(metadata.isNull())
      return;

    const std::map<std::string, float> propertiesMap = metadata->GetPropertiesMap();
    this->FillPandoraTrackIDScoreVars(propertiesMap);

    const int parentId     = pfp->Parent();
    auto const& parentIter = pfpMap.find(parentId);

    if(parentIter == pfpMap.end())
      return;
    
    const art::Ptr<recob::PFParticle> parent = parentIter->second;

    recoPrimary = parent->IsPrimary();
  }

  void Razzled::FillPandoraTrackIDScoreVars(const std::map<std::string, float> &propertiesMap)
  {
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_trackScore, "TrackScore");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_chargeEndFrac, "LArThreeDChargeFeatureTool_EndFraction");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_chargeFracSpread, "LArThreeDChargeFeatureTool_FractionalSpread");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_linearFitDiff, "LArThreeDLinearFitFeatureTool_DiffStraightLineMean");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_linearFitLength, "LArThreeDLinearFitFeatureTool_Length");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_linearFitGapLength, "LArThreeDLinearFitFeatureTool_MaxFitGapLength");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_linearFitRMS, "LArThreeDLinearFitFeatureTool_SlidingLinearFitRMS");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_openAngleDiff, "LArThreeDOpeningAngleFeatureTool_AngleDiff");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_secondaryPCARatio, "LArThreeDPCAFeatureTool_SecondaryPCARatio");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_tertiaryPCARatio, "LArThreeDPCAFeatureTool_TertiaryPCARatio");
    this->FillPandoraTrackIDScoreVar(propertiesMap, pfp_vertexDist, "LArThreeDVertexDistanceFeatureTool_VertexDistance");
  }

  void Razzled::FillPandoraTrackIDScoreVar(const std::map<std::string, float> &propertiesMap, float &var, std::string name)
  {
    auto propertiesMapIter = propertiesMap.find(name);
    if(propertiesMapIter == propertiesMap.end())
      throw cet::exception("Razzled") << "Error finding variable -- " << name;
    var = propertiesMapIter->second;
  }


  void Razzled::FillTrackMetrics(const art::Ptr<recob::Track> &track)
  {
    trk_length = track->Length();

    const TVector3 start = track->Start<TVector3>();
    const TVector3 end   = track->End<TVector3>();
    trackContained       = this->InFV(start) && this->InFV(end);

    if(!fMakeTree)
      return;

    trackStartX = start.X();
    trackStartY = start.Y();
    trackStartZ = start.Z();

    trackEndX = end.X();
    trackEndY = end.Y();
    trackEndZ = end.Z();
  }

  void Razzled::FillMCSMetrics(const art::Ptr<recob::MCSFitResult> &mcs)
  {
    if(mcs->scatterAngles().empty())
      return;
    
    unsigned int counter = 0;
    float maxScatter = 0.f, meanScatter = 0.f;

    for(auto const& angle : mcs->scatterAngles())
      {
        if(angle < 0)
          continue;

        maxScatter = std::max(maxScatter, angle);
        meanScatter += angle;
        counter++;
      }
    
    if(!counter)
      return;
    
    trk_mcsScatterMean     = meanScatter / counter;
    trk_mcsScatterMaxRatio = maxScatter / meanScatter;
  }

  void Razzled::FillMomDiff(const art::Ptr<RangeP> &rangeP, const art::Ptr<recob::MCSFitResult> &mcs)
  {
    const float rangeMom = rangeP->range_p;
    const float mcsMom   = mcs->fwdMomentum();

    trk_momDiff = (rangeMom > 0 && mcsMom > 0) ? (mcsMom - rangeMom) / rangeMom : -5.f;
  }

  void Razzled::FillChi2PIDMetrics(const art::Ptr<anab::ParticleID> &pid)
  {
    trk_chi2PIDMuon   = 0.;
    trk_chi2PIDProton = 0.;

    float chi2PIDPion = 0.;

    const std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pid->ParticleIDAlgScores();

    for(size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
      {
        const anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);

        if(AlgScore.fAlgName == "Chi2")
          {
            switch(TMath::Abs(AlgScore.fAssumedPdg))
              {
              case 13:
                trk_chi2PIDMuon = AlgScore.fValue;
              case 211:
                chi2PIDPion = AlgScore.fValue;
              case 2212:
                trk_chi2PIDProton = AlgScore.fValue;
              }
          }
      }

    trk_chi2PIDMuonPionDiff = trk_chi2PIDMuon - chi2PIDPion;
  }

  void Razzled::FillStoppingChi2Metrics(const art::Ptr<StoppingChi2Fit> &stoppingChi2)
  {
    const float pol0Chi2 = stoppingChi2->pol0Chi2;
    const float expChi2  = stoppingChi2->expChi2;

    trk_stoppingdEdxChi2Ratio = (pol0Chi2 > 0.f && expChi2 > 0.f) ? pol0Chi2 / expChi2 : -5.f;
    trk_chi2Pol0dEdxFit       = stoppingChi2->pol0Fit;
  }

  void Razzled::FillShowerMetrics(const art::Ptr<recob::Shower> &shower, const std::vector<art::Ptr<recob::Hit>> &hits)
  {
    const geo::GeometryCore* geom = lar::providerFrom<geo::Geometry>();

    const float length = shower->Length();
    shw_openAngle      = TMath::RadToDeg() * shower->OpenAngle();

    const TVector3 start = shower->ShowerStart();
    const TVector3 end   = start + length * shower->Direction();

    showerContained = this->InFV(start) && this->InFV(end);

    std::array<int, 3> showerPlaneHits = { 0, 0, 0 };

    for(auto const& hit : hits)
      showerPlaneHits[hit->WireID().Plane]++;

    std::array<float, 3> showerPlanePitches = { -1.f, -1.f, -1.f };

    for(geo::PlaneGeo const& plane : geom->Iterate<geo::PlaneGeo>())
      {
        const float angleToVert = geom->WireAngleToVertical(plane.View(), plane.ID()) - 0.5 * M_PI;
        const float cosgamma    = std::abs(std::sin(angleToVert) * shower->Direction().Y() + std::cos(angleToVert) * shower->Direction().Z());

        showerPlanePitches[plane.ID().Plane] = plane.WirePitch() / cosgamma;
      }

    const int bestPlane = shower->best_plane();

    if(bestPlane < 0 || bestPlane > 3)
      throw cet::exception("Razzled") << "Best plane: " << bestPlane;

    shw_bestdEdx = shower->dEdx()[bestPlane];
    shw_bestdEdx = std::min(shw_bestdEdx, 20.f);
    shw_bestdEdx = std::max(shw_bestdEdx, -5.f);

    const float bestEnergy  = shower->Energy()[bestPlane];
    const int bestPlaneHits = showerPlaneHits[bestPlane];
    const float bestPitch   = showerPlanePitches[bestPlane];
    const float wiresHit    = bestPitch > std::numeric_limits<float>::epsilon() ? length / bestPitch : -5.f;

    shw_sqrtEnergyDensity = (length > 0 && bestEnergy > 0) ? std::sqrt(bestEnergy) / length : -5.f;
    shw_modHitDensity     = wiresHit > 1.f ? bestPlaneHits / wiresHit : -5.f;
    shw_modHitDensity     = std::min(shw_modHitDensity, 40.f);

    if(!fMakeTree)
      return;

    showerStartX = start.X();
    showerStartY = start.Y();
    showerStartZ = start.Z();

    showerEndX = end.X();
    showerEndY = end.Y();
    showerEndZ = end.Z();
  }

  void Razzled::FillShowerConversionGap(const art::Ptr<recob::PFParticle> &pfp, const std::map<size_t, art::Ptr<recob::PFParticle>> &pfpMap,
					const art::Ptr<recob::Shower> &shower, const art::FindOneP<recob::Vertex> &pfpsToVertices)
  {
    const int parentId     = pfp->Parent();
    auto const& parentIter = pfpMap.find(parentId);

    if(parentIter == pfpMap.end())
      return;
    
    const art::Ptr<recob::PFParticle> parent = parentIter->second;

    const art::Ptr<recob::Vertex> vertex = pfpsToVertices.at(parent.key());

    if(vertex.isNull())
      return;

    const geo::Point_t vertexPoint = vertex->position();
    const TVector3 vertexTV3       = { vertexPoint.X(), vertexPoint.Y(), vertexPoint.Z() };

    shw_convGap = (shower->ShowerStart() - vertexTV3).Mag();
    shw_convGap = std::min(shw_convGap, 50.f);
  }

  bool Razzled::InFV(const TVector3 &pos)
  {
    return (pos.X() > fXMin && pos.X() < fXMax && pos.Y() > fYMin && pos.Y() < fYMax && pos.Z() > fZMin && pos.Z() < fZMax);
  }

  std::map<size_t, art::Ptr<recob::PFParticle>> Razzled::GetPFPMap(const std::vector<art::Ptr<recob::PFParticle>> &pfps)
  {
    std::map<size_t, art::Ptr<recob::PFParticle>> pfpMap;
    for(auto const& pfp : pfps)
      pfpMap[pfp->Self()] = pfp;

    return pfpMap;
  }

  std::string Razzled::PDGString(const int pdg)
  {
    switch (std::abs(pdg))
      {
      case 11:
        return "Electron";
      case 13:
        return "Muon";
      case 22:
        return "Photon";
      case 211:
        return "Pion";
      case 2212:
        return "Proton";
      default:
        return "Other";
      }
  }

  MVAPID Razzled::RunMVA()
  {
    const std::vector<float> mvaScores = reader->EvaluateMulticlass(fMethodName);

    MVAPID pidResults;

    pidResults.AddScore(11, mvaScores.at(0));
    pidResults.AddScore(13, mvaScores.at(1));
    pidResults.AddScore(22, mvaScores.at(2));
    pidResults.AddScore(211, mvaScores.at(3));
    pidResults.AddScore(2212, mvaScores.at(4));
    pidResults.AddScore(0, mvaScores.at(5));

    if(!fMakeTree)
      return pidResults;

    electronScore = mvaScores.at(0);
    muonScore     = mvaScores.at(1);
    photonScore   = mvaScores.at(2);
    pionScore     = mvaScores.at(3);
    protonScore   = mvaScores.at(4);
    otherScore    = mvaScores.at(5);

    bestScore = pidResults.BestScore();
    bestPDG   = pidResults.BestPDG();

    return pidResults;
  }
}

DEFINE_ART_MODULE(sbn::Razzled)
