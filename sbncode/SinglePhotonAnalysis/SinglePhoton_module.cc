#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//#include "ubobj/Optical/UbooneOpticalFilter.h"

// Helper function for PID stuff
//#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "sbnobj/Common/CRT/CRTHit.hh"

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TRandom3.h"
#include "TGeoPolygon.h"

//#include "Pandora/PdgTable.h"
#include <chrono>

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <sys/stat.h>

#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_gadget.h"
#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_PandoraPFParticles.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/init_branches.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/Processors.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/analyze_PandoraReco.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/analyze_MC.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/reco_truth_matching.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/fiducial_volume.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/second_shower_search.h"
//#include "sbncode/SinglePhotonAnalysis/Libraries/DBSCAN.h"
#include "sbncode/SinglePhotonAnalysis/SEAview/SEAviewer.h"

namespace single_photon
{

  /**
   *  @brief  SinglePhoton class
   */
  class SinglePhoton : public art::EDFilter
  {
    public:
      var_all vars;
      para_all paras;
      /**
       *  @brief  Constructor
       *
       *  @param  pset the set of input fhicl parameters
       */
      SinglePhoton(fhicl::ParameterSet const &pset);

      /**
       *  @brief  Configure memeber variables using FHiCL parameters
       *
       *  @param  pset the set of input fhicl parameters
       */
      void reconfigure(fhicl::ParameterSet const &pset);

      /**
       *  @brief  Analyze an event!
       *
       *  @param  evt the art event to analyze
       */
      bool filter(art::Event &evt) override;

      /**
       *  @brief  Begin the job, setting up !
       *
       */
      void beginJob() override;

      /**
       *  @brief  End the job, setting down !
       *
       */
      void endJob() override;
      /**
       * @brief: grab run, subrun number, and subrun POT, fill the TTree */
      bool beginSubRun(art::SubRun& sr) override;
      bool endSubRun(art::SubRun& sr) override;

  };

  //Constructor from .fcl parameters
  SinglePhoton::SinglePhoton(fhicl::ParameterSet const &pset) : art::EDFilter(pset)
  {
    std::cout<<"SinglePhoton::"<<__FUNCTION__<<" kicks off ------------------------------"<<std::endl;
    this->reconfigure(pset);
    //Set up some detector, timing, spacecharge and geometry services
    //Keng, grab theDetector and detClocks in each event.
    //        theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    //        detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    paras.s_SCE = lar::providerFrom<spacecharge::SpaceChargeService>();//Get space charge service
    paras.s_geom =  lar::providerFrom<geo::Geometry>();
    std::cout<<"SinglePhoton::"<<__FUNCTION__<<" now can start jobs --------------------"<<std::endl;

  }

  //Reconfigure the internal class parameters from .fcl parameters
  void SinglePhoton::reconfigure(fhicl::ParameterSet const &pset)
  {
    //Get Geometry
    std::cout<<"SinglePhoton::"<<__FUNCTION__<<" starts ------------------------------"<<std::endl;


    //input parameters for what file/mode were running in
    paras.s_print_out_event = pset.get<bool>("PrintOut", false);
    g_is_verbose =    pset.get<bool>("Verbose",false);
    paras.s_is_data =      pset.get<bool>("isData",false);
    paras.s_is_overlayed =  pset.get<bool>("isOverlayed",false);
    paras.s_is_textgen =    pset.get<bool>("isTextGen",false);

    //some specific additonal info, default not include
    paras.s_use_PID_algorithms =  pset.get<bool>("usePID",false);
    paras.s_use_delaunay =    pset.get<bool>("useDelaunay",false);
    paras.s_delaunay_max_hits =  pset.get<int>("maxDelaunayHits",1000);

    //When we moved to filter instead of analyzer, kept ability to run 2g filter. A bit depreciated (not in code, but in our use case)
    paras.s_fill_trees          = pset.get<bool>("FillTrees",true);
    paras.s_run_pi0_filter      = pset.get<bool>("RunPi0Filter",false);
    paras.s_run_pi0_filter_2g1p = pset.get<bool>("FilterMode2g1p",false);
    paras.s_run_pi0_filter_2g0p = pset.get<bool>("FilterMode2g0p",false);

    if(paras.s_run_pi0_filter) paras.s_is_data = true;// If running in filter mode, treat all as data

    //Some output for logging
    std::cout<<"SinglePhoton::reconfigure || whats configured? "<<std::endl;
    std::cout<<"SinglePhoton::reconfigure || paras.s_is_data: "<<paras.s_is_data<<std::endl;
    std::cout<<"SinglePhoton::reconfigure || paras.s_is_overlayed: "<<paras.s_is_overlayed<<std::endl;
    std::cout<<"SinglePhoton::reconfigure || paras.s_is_textgen: "<<paras.s_is_textgen<<std::endl;


    //Instead of running over ALL evnts in a file, can pass in a file to run over just the run:subrun:event listed
    paras.s_runSelectedEvent    = pset.get<bool>("SelectEvent", false);
    paras.s_selected_event_list = pset.get<std::string>("SelectEventList", "");

    //Studies for photo nuclear EventWeights
    paras.s_runPhotoNuTruth = pset.get<bool>("RunPhotoNu",false); 

    //Ability to save some FULL eventweight components, rather than run later. Useful for systematic studies. Harcoded to two currently (TODO)
    paras.s_runTrueEventweight = pset.get<bool>("RunTrueEventWeight",false); 
    paras.s_true_eventweight_label = pset.get<std::string>("true_eventweight_label","eventweight");
    paras.s_Spline_CV_label = pset.get<std::string>("SplineCVLabel", "eventweight");//4to4aFix");


    //Input ArtRoot data products 
    paras.s_pandoraLabel = pset.get<std::string>("PandoraLabel");
    paras.s_trackLabel = pset.get<std::string>("TrackLabel");
    //KENG        paras.s_sliceLabel = pset.get<std::string>("SliceLabel","pandora");
    paras.s_showerLabel = pset.get<std::string>("ShowerLabel");
    paras.s_pidLabel = pset.get<std::string>("ParticleIDLabel","pandoracalipidSCE");
    paras.s_caloLabel = pset.get<std::string>("CaloLabel");
    paras.s_flashLabel = pset.get<std::string>("FlashLabel");
    paras.s_potLabel = pset.get<std::string>          ("POTLabel");
    paras.s_hitfinderLabel = pset.get<std::string>      ("HitFinderModule", "gaushit");
    //KENG  no such labels in sbnd?
    //KENG  paras.s_badChannelLabel = pset.get<std::string>      ("BadChannelLabel","badmasks");
    //KENG  paras.s_showerKalmanLabel = pset.get<std::string>      ("ShowerTrackFitter","pandoraKalmanShower");
    //KENG  paras.s_showerKalmanCaloLabel =  pset.get<std::string>  ("ShowerTrackFitterCalo","pandoraKalmanShowercali");
    paras.s_generatorLabel = pset.get<std::string>      ("GeneratorLabel","generator");
    //KENG        paras.s_mcTrackLabel = pset.get<std::string>        ("MCTrackLabel","mcreco");
    //KENG        paras.s_mcShowerLabel = pset.get<std::string>        ("MCShowerLabel","mcreco");
    paras.s_geantModuleLabel = pset.get<std::string>      ("GeantModule","largeant");
    //KENG        paras.s_backtrackerLabel = pset.get<std::string>      ("BackTrackerModule","gaushitTruthMatch");
    paras.s_hitMCParticleAssnsLabel = pset.get<std::string>  ("HitMCParticleAssnLabel","gaushitTruthMatch");
    paras.s_shower3dLabel = pset.get<std::string>        ("Shower3DLabel","shrreco3d");



    //Flash related variables. 
    paras.s_beamgate_flash_start = pset.get<double>("beamgateStartTime",3.2); //Defaults to MC for now. Probably should change
    paras.s_beamgate_flash_end = pset.get<double>("beamgateEndTime",4.8);

    // paras.s_badChannelProducer = pset.get<std::string>("BadChannelProducer","nfspl1");
    //        paras.s_badChannelProducer = pset.get<std::string>("BadChannelProducer","simnfspl1");

    //CRT related variables, should run only for RUN3+ enabled
    paras.s_runCRT =      pset.get<bool>("runCRT",false);
    paras.s_CRTTzeroLabel =  pset.get<std::string>("CRTTzeroLabel","crttzero");
    paras.s_CRTVetoLabel =  pset.get<std::string>("CRTVetoLabel","crtveto");
    paras.s_CRTHitProducer =  pset.get<std::string>("CRTHitProducer", "crthitcorr");
    paras.s_DTOffset =    pset.get<double>("DTOffset" , 68600.); //us, taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter.fcl
    paras.s_Resolution =    pset.get<double>("Resolution" ,  1.0); //us, taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter.fcl
    //        paras.s_DAQHeaderProducer = pset.get<std::string>("DAQHeaderProducer" ,  "daq");

    //Some track calorimetry parameters
    paras.s_track_calo_min_dEdx = pset.get<double>("Min_dEdx",0.005);
    paras.s_track_calo_max_dEdx = pset.get<double>("Max_dEdx", 30);
    paras.s_track_calo_min_dEdx_hits = pset.get<double>("Min_dEdx_hits",5); //might be good?
    paras.s_track_calo_trunc_fraction = pset.get<double>("TruncMeanFraction",20.0);

    //Some shower calorimetry parameters
    paras.s_work_function =      pset.get<double>("work_function");
    paras.s_recombination_factor =  pset.get<double>("recombination_factor");
    paras.s_gain_mc =          pset.get<std::vector<double>>("gain_mc");
    paras.s_gain_data =        pset.get<std::vector<double>>("gain_data");
    paras.s_wire_spacing =      pset.get<double>("wire_spacing");
    paras.s_width_dqdx_box =      pset.get<double>("width_box");
    paras.s_length_dqdx_box =      pset.get<double>("length_box");
    paras.s_truthmatching_signaldef = pset.get<std::string>("truthmatching_signaldef");

    //A seperate mode to run over AllPFPs and not just slice particles 
    paras.s_run_all_pfps = pset.get<bool>("runAllPFPs",false);


    //Some paramaters for counting protons & photons
    paras.s_exiting_photon_energy_threshold = pset.get<double>("exiting_photon_energy");
    paras.s_exiting_proton_energy_threshold = pset.get<double>("exiting_proton_energy");
    paras.s_max_conv_dist =          pset.get<double>("convention_distance_cutoff", 999);
    paras.s_mass_pi0_mev =  139.57;

    //SEAviwer Settings for shower clustering and proton stub finding
    //Have two sets:
    //Base SEAview is for Second Shower Veto
    paras.s_runSEAview =      pset.get<bool>("runSEAviewShower", false);
    paras.s_SEAviewHitThreshold = pset.get<double>("SEAviewShowerHitThreshold",25);
    paras.s_SEAviewPlotDistance = pset.get<double>("SEAviewShowerPlotDistance",80);
    paras.s_SEAviewDbscanMinPts = pset.get<double>("SEAviewShowerDBSCANMinPts",8);
    paras.s_SEAviewDbscanEps =  pset.get<double>("SEAviewShowerDBSCANEps",4);
    paras.s_SEAviewMaxPtsLinFit = pset.get<double>("SEAviewShowerMaxHitsLinFit",20.0);
    paras.s_SEAviewMakePDF =    pset.get<bool>("SEAviewShowerMakePDF",false);
    paras.s_SEAviewNumRecoShower= pset.get<int>("SEAviewShowerNumRecoShower", -1);
    paras.s_SEAviewNumRecoTrack = pset.get<int>("SEAviewShowerNumRecoTrack", -1);

    // Second set is for Proton Stub finding
    paras.s_runSEAviewStub = pset.get<bool>("runSEAviewStub", false);
    paras.s_SEAviewStubHitThreshold = pset.get<double>("SEAviewStubHitThreshold",25);
    paras.s_SEAviewStubPlotDistance = pset.get<double>("SEAviewStubPlotDistance",80);
    paras.s_SEAviewStubDbscanMinPts = pset.get<double>("SEAviewStubDBSCANMinPts",1);
    paras.s_SEAviewStubDbscanEps =    pset.get<double>("SEAviewStubDBSCANEps",1);
    paras.s_SEAviewStubMakePDF =      pset.get<bool>("SEAviewStubMakePDF",false);
    paras.s_SEAviewStubNumRecoShower =pset.get<int>("SEAviewStubNumRecoShower", -1);
    paras.s_SEAviewStubNumRecoTrack = pset.get<int>("SEAviewStubNumRecoTrack", -1);

    paras.s_make_sss_plots = true;

    //Misc setup  CHECK
    setTPCGeom(paras);//set activeTPC
    paras.rangen = new TRandom3(22);

    std::cout<<"SinglePhoton::"<<__FUNCTION__<<" finishes ------------------------------"<<std::endl;

  }

  //------------------------------------------------------------------------------------------------------------------------------------------



  //--------------------------------------- Primary Filter------------------------------------------------------------------------------------
  // Runs over every artroot event
  bool SinglePhoton::filter(art::Event &evt)
  {
    //    //Grab services 
    auto detClocks = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);//it is detinfo::DetectorClocksData 
    auto theDetector = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, detClocks);//it is detinfo::DetectorPropertiesData 
    std::cout<<"SinglePhoton::"<<__FUNCTION__<<" a new event ------------------------------"<<std::endl;

    //Clear all output branches 
    ClearMeta(vars);
    ClearIsolation(vars);
    ClearSecondShowers(vars);
    ClearSecondShowers3D(vars);
    ClearStubs(vars);
    ClearFlashes(vars);
    ClearTracks(vars);
    ClearShowers(vars);
    ClearMCTruths(vars);
    ClearEventWeightBranches(vars);
    ClearGeant4Branches(vars);
    ClearSlices(vars);

    Save_EventMeta( evt, vars);

    //if module is run in selected-event mode, and current event is not in the list, skip it
    if(paras.s_runSelectedEvent && !IsEventInList(vars.m_run_number, vars.m_subrun_number, vars.m_event_number, vars)){
      std::cout << "SinglePhoton::analyze()\t||\t event " << vars.m_run_number << "/" << vars.m_subrun_number << "/" << vars.m_event_number << " is not in the list, skip it" << std::endl;
      return true;
    }

    //Timing and TPC info
    //        auto const TPC = (*geom).begin_TPC();
    //        auto ID = TPC.ID();
    //        paras.s_Cryostat = ID.Cryostat;
    //        paras.s_TPC = ID.TPC;

    //use this in calculating dQdx
    paras._time2cm = sampling_rate(detClocks) / 1000.0 * theDetector.DriftVelocity( theDetector.Efield(), theDetector.Temperature() );//found in ProtoShowerPandora_tool.cc


    //******************************Setup*****************Setup**************************************/
    //Use a new PandoraPFParticle class! 
    //Each class has a core PFParticle member and other associated recob objects;


    //        std::vector<art::Ptr<recob::Track>> trackVector;
    //        art::fill_ptr_vector(trackVector,trackHandle);

    //Collect the PFParticles from the event. This is the core!
    art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(paras.s_pandoraLabel);
    std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
    art::fill_ptr_vector(pfParticleVector,pfParticleHandle);
    //So a cross check
    if (!pfParticleHandle.isValid())
    { mf::LogDebug("SinglePhoton") << "  Failed to find the PFParticles.\n";
      return (paras.s_run_pi0_filter ? false : true) ;
    }

    //get the cluster handle for the dQ/dx calc
    art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(paras.s_pandoraLabel);
    std::vector< art::Ptr<recob::Cluster> > clusterVector;
    art::fill_ptr_vector(clusterVector,clusterHandle);

    // This is another pandora helper. I don't like PFParticle ID lookups but I guess lets keep for now;
    // typedef std::map< size_t, art::Ptr<recob::PFParticle>>
    // Produce a map of the PFParticle IDs for fast navigation through the hierarchy
    //        PFParticleIdMap pfParticleMap;
    //        this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);


    //And some verticies.        
    art::ValidHandle<std::vector<recob::Vertex>> const & vertexHandle = evt.getValidHandle<std::vector<recob::Vertex>>(paras.s_pandoraLabel);
    std::vector<art::Ptr<recob::Vertex>> vertexVector;
    art::fill_ptr_vector(vertexVector,vertexHandle);
    if(vertexVector.size()>0)  vars.m_number_of_vertices++;

    //PFParticle to Vertices
    art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, paras.s_pandoraLabel);
    std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticlesToVerticesMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
      auto pfp = pfParticleVector[i];
      pfParticlesToVerticesMap[pfp] =vertices_per_pfparticle.at(pfp.key());
    }

    //------- 3D showers
    art::FindOneP<recob::Shower> showerreco3D_per_pfparticle(pfParticleHandle, evt, paras.s_shower3dLabel);

    if(g_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get PandoraMetadata"<<std::endl;
    //add the associaton between PFP and metadata, this is important to look at the slices and scores
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt,  paras.s_pandoraLabel);

    //PFPartciles -->Showers/Tracks
    art::FindManyP< recob::Track     > pfPartToTrackAssoc(pfParticleHandle, evt, paras.s_trackLabel);
    art::FindManyP< recob::Shower    > pfPartToShowerAssoc(pfParticleHandle, evt, paras.s_showerLabel);
    //Assign Cluster here;
    art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfParticleHandle, evt, paras.s_pandoraLabel);
    art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, paras.s_pandoraLabel);

    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > pfParticleToMetadataMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
      const art::Ptr<recob::PFParticle> pfp = pfParticleVector[i];
      pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
    }

    //try the new class, Keng

    std::vector<PandoraPFParticle> allPPFParticles;
    for(size_t index=0; index< pfParticleVector.size(); ++index){
      auto pfp = pfParticleVector[index];
      int temp_key = pfp.key();

      PandoraPFParticle temp_pf(
          pfp, 
          pfPartToMetadataAssoc.at(temp_key),
          vertices_per_pfparticle.at(temp_key),
          clusters_per_pfparticle.at(temp_key),
          pfPartToShowerAssoc.at(temp_key),
          pfPartToTrackAssoc.at(temp_key),
          hits_per_cluster
          );

      allPPFParticles.push_back(temp_pf);

    }
    PPFP_FindAncestor(allPPFParticles);


    //Add slices & hits info.
    //Slices
    art::ValidHandle<std::vector<recob::Slice>> const & sliceHandle  = evt.getValidHandle<std::vector<recob::Slice>>(paras.s_pandoraLabel);
    std::vector<art::Ptr<recob::Slice>> sliceVector;
    art::fill_ptr_vector(sliceVector,sliceHandle);

    //And some associations
    art::FindManyP<recob::PFParticle> pfparticles_per_slice(sliceHandle, evt, paras.s_pandoraLabel);
    art::FindManyP<recob::Hit> hits_per_slice(sliceHandle, evt, paras.s_pandoraLabel);

    //Slice to PFParticle
    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::PFParticle>> > sliceToPFParticlesMap;
    std::map<int, std::vector<art::Ptr<recob::PFParticle>> > sliceIDToPFParticlesMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
      auto slice = sliceVector[i];
      PPFP_FindSliceIDandHits(allPPFParticles, slice, pfparticles_per_slice.at(slice.key()), hits_per_slice.at(slice.key()) );
      sliceToPFParticlesMap[slice] =pfparticles_per_slice.at(slice.key());
      sliceIDToPFParticlesMap[slice->ID()] = pfparticles_per_slice.at(slice.key());
    }

    //Mark all Pandora neutrino slices as neutrino slices and assign the highes neutrino score.
    vars.pfp_w_bestnuID = DefineNuSlice(allPPFParticles);

    Save_PFParticleInfo( allPPFParticles, vars, paras);
    ResizeSlices(vars.m_reco_slice_num, vars);

    //Slice to Hits
    std::map<int, std::vector<art::Ptr<recob::Hit>> > sliceIDToHitsMap;
    for(size_t i=0; i< sliceVector.size(); ++i){

      auto slice = sliceVector[i];
      sliceIDToHitsMap[slice->ID()] = hits_per_slice.at(slice.key());
    }


    // Once we have actual verticies, lets concentrate on JUST the neutrino PFParticles for now:
    //--------------------------------
    // Produce two PFParticle vectors containing final-state particles:
    // 1. Particles identified as cosmic-rays - recontructed under cosmic-hypothesis
    // 2. Daughters of the neutrino PFParticle - reconstructed under the neutrino hypothesis
    std::vector< art::Ptr<recob::PFParticle> > nuParticles;
    std::vector< art::Ptr<recob::PFParticle> > crParticles;

    for(size_t jndex=0; jndex< allPPFParticles.size(); ++jndex){
      PandoraPFParticle* temp_pf = &allPPFParticles[jndex];
      int temp_id = temp_pf->get_AncestorID();//indentify the jndex for the ancestor PPFParticle
      PandoraPFParticle* temp_ancestor_ppfp = PPFP_GetPPFPFromPFID( allPPFParticles, temp_id);

      if ( temp_ancestor_ppfp->get_IsNeutrino() || paras.s_run_all_pfps){
        nuParticles.push_back( temp_pf->pPFParticle);
        //        std::cout<<"Get nuParticles from pfpID "<<temp_pf->get_PFParticleID()<<std::endl;
      }else{
        crParticles.push_back( temp_pf->pPFParticle);
      }
    }


    if(g_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Spacepoints"<<std::endl;
    //Spacepoint associaitions
    art::FindManyP<recob::SpacePoint> spacePoints_per_pfparticle(pfParticleHandle, evt, paras.s_pandoraLabel);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
    for(size_t i=0; i< nuParticles.size(); ++i){
      const art::Ptr<recob::PFParticle> pfp = nuParticles[i];
      pfParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());
    }

    if(g_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Clusters"<<std::endl;
    //Get a map between the PFP's and the clusters  they're imporant for the shower dQ/dx
    //Also need a map between clusters and hits
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
    std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;
    //fill map PFP to Clusters
    for(size_t i=0; i< nuParticles.size(); ++i){
      auto pfp = nuParticles[i];
      pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
    }
    //fill map Cluster to Hits
    for(size_t i=0; i< clusterVector.size(); ++i){
      auto cluster = clusterVector[i];
      clusterToHitsMap[cluster] = hits_per_cluster.at(cluster.key());
    }
    if(g_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Build hits to PFP Maps"<<std::endl;


    //OK Here we build two IMPORTANT maps for the analysis, (a) given a PFParticle get a vector of hits..
    //and (b) given a single hit, get the PFParticle it is in (MARK: is it only one? always? RE-MARK: Yes)
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;


    //use pfp->cluster and cluster->hit to build pfp->hit map
    //for each PFP
    for(size_t i=0; i<nuParticles.size(); ++i){
      auto pfp = nuParticles[i];

      //get the associated clusters
      std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp] ;

      //make empty vector to store hits
      std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};


      //for each cluster, get the associated hits
      for (art::Ptr<recob::Cluster> cluster: clusters_vec){
        std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];

        //insert hits into vector
        hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
      }

      //fill the map
      pfParticleToHitsMap[pfp] = hits_for_pfp;

    }//for each pfp



    /**************************************************************************
     * For SEAview: grab cosmic-related PFPaticles and recob::Hits
     *
     **************************************************************************/
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > cr_pfParticleToClustersMap;
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > cr_pfParticleToHitsMap;

    //first, collect all daughters of primary cosmic
    //        int nuvars.m_primary_cosmic_particle = crParticles.size();
    //        for(int i =0; i!=nuvars.m_primary_cosmic_particle; ++i){
    //            auto& pParticle = crParticles[i];
    //            for(const size_t daughterId : pParticle->Daughters())
    //            {
    //                if (pfParticleMap.find(daughterId) == pfParticleMap.end())
    //                    throw cet::exception("SinglePhoton") << "  Invalid PFParticle collection!";
    //
    //                crParticles.push_back(pfParticleMap.at(daughterId));
    //            }
    //        }

    //second, build PFP to hits map for cosmic-related PFParticles
    for(size_t i=0; i< crParticles.size(); ++i){
      auto pfp = crParticles[i];
      cr_pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
    }

    for(size_t i=0; i< crParticles.size(); ++i){
      auto pfp = crParticles[i];

      // std::cout<<"starting to match to hits for pfp "<<pfp->Self()<<std::endl;
      //get the associated clusters
      std::vector<art::Ptr<recob::Cluster>> clusters_vec  = cr_pfParticleToClustersMap[pfp] ;

      //make empty vector to store hits
      std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

      // std::cout<<"-- there are "<<clusters_vec.size()<<" associated clusters"<<std::endl;

      //for each cluster, get the associated hits
      for (art::Ptr<recob::Cluster> cluster: clusters_vec){
        std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];

        //   std::cout<<"looking at cluster in pfp "<<pfp->Self()<<" with "<<hits_vec.size() <<" hits"<<std::endl;
        //insert hits into vector
        hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
      }

      //fill the map
      cr_pfParticleToHitsMap[pfp] = hits_for_pfp;
      //std::cout<<"saving a total of "<<hits_for_pfp.size()<<" hits for pfp "<<pfp->Self()<<std::endl;

    }//for each pfp
    //        std::cout << "Guanqun SEAview test: initial crParticle size: " << nuvars.m_primary_cosmic_particle << ", current size: " << crParticles.size() <<" (including daughters)"<< std::endl;
    /********************************************************************
     * End of SEAview: grab cosmic-related PFPaticles and recob::Hits
     *
     **************************************************************************/


    //******************************  Collecting Info and Analysis  **************************************/

    // These are the vectors to hold the tracks and showers for the final-states of the reconstructed neutrino
    //At this point, nuParticles is a std::vector< art::Ptr<recon::PFParticle>> of the PFParticles that we are interested in.
    //tracks is a vector of recob::Tracks and same for showers.
    //Implicitly, tracks.size() + showers.size() =  nuParticles.size(); At this point I would like two things.
    std::vector< art::Ptr<recob::Track> > tracks;
    std::vector< art::Ptr<recob::Shower> > showers;
    std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle >> trackToNuPFParticleMap; 
    std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap;

    if(g_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Tracks and Showers"<<std::endl;

    for(size_t jndex=0; jndex< allPPFParticles.size(); ++jndex){
      PandoraPFParticle* temp_pf = &allPPFParticles[jndex];
      if(!temp_pf->get_IsNuSlice() ) continue;//pass slice not defined as nu-slice in above lines
      if(temp_pf->get_HasTrack()){
        tracks.push_back(temp_pf->pTrack);
        trackToNuPFParticleMap[temp_pf->pTrack] = temp_pf->pPFParticle;
      }
      if(temp_pf->get_HasShower()){ 
        showers.push_back(temp_pf->pShower);
        showerToNuPFParticleMap[temp_pf->pShower] = temp_pf->pPFParticle;
      }
    }

    //Helper function (can be found below) to collect tracks and showers in neutrino slice
    //this->CollectTracksAndShowers(nuParticles, pfParticleMap,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);

    //Track Calorimetry. Bit odd here but bear with me, good to match and fill here

    //tracks
    art::ValidHandle<std::vector<recob::Track>> const & trackHandle  = evt.getValidHandle<std::vector<recob::Track>>(paras.s_trackLabel);
    //Keng, use pandoraCalo to get the correct Calorimetry;
    art::FindManyP<anab::Calorimetry> calo_per_track(trackHandle, evt, paras.s_caloLabel);
    //        std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<anab::Calorimetry>> > trackToCalorimetryMap;
    //So a cross check
    if (!calo_per_track.isValid())
    {
      mf::LogDebug("SinglePhoton") << "  Failed to get Assns between recob::Track and anab::Calorimetry.\n";
      return (paras.s_run_pi0_filter ? false : true);
    }

//    art::FindOneP<anab::ParticleID> pid_per_track(trackHandle, evt, paras.s_pidLabel);
    art::FindOneP<anab::ParticleID>* pid_per_track = nullptr;

	if(paras.s_use_PID_algorithms){
		std::cout<<"Using PID stuff"<<__LINE__<<std::endl;
		pid_per_track = new art::FindOneP<anab::ParticleID>(trackHandle, evt, paras.s_pidLabel);
	}

    std::map<art::Ptr<recob::Track>, art::Ptr<anab::ParticleID> > trackToPIDMap;
    //Loop over all tracks we have to fill calorimetry map
    for(size_t i=0; i< tracks.size(); ++i){
      if(calo_per_track.at(tracks[i].key()).size() ==0){
        std::cerr<<"Track Calorimetry Breaking! the vector of calo_per_track is of length 0 at this track."<<std::endl;
      }
      //size_t calo_size = calo_per_track.at(tracks[i].key()).size();
      //std::cout<<"Track Calo from producer: "<<paras.s_caloLabel<<" has "<<calo_size<<" anab::Calorimetry objects associaed."<<std::endl;
      //trackToCalorimetryMap[tracks[i]] = calo_per_track.at(tracks[i].key());
      PandoraPFParticle * temp_ppfp = PPFP_GetPPFPFromTrack(allPPFParticles, tracks[i]);
      temp_ppfp->set_Calorimetries(calo_per_track.at(tracks[i].key()));

      if(paras.s_use_PID_algorithms){
        temp_ppfp->set_HasPID(true);
        temp_ppfp->set_ParticleID(pid_per_track->at(tracks[i].key()));
      }
    }


    // If we want PID algorithms to run. do so here
    // Build a map to get PID from PFParticles, then call PID collection function
    if(paras.s_use_PID_algorithms){//yes, as set in the FHiCL;
      for(size_t i=0; i< tracks.size(); ++i){
        trackToPIDMap[tracks[i]] = pid_per_track->at(tracks[i].key());
      }
    }



    //**********************************************************************************************/
    //**********************************************************************************************/
    //---------------------------------- MC TRUTH, MC Only---------------------------
    //**********************************************************************************************/
    //**********************************************************************************************/

    //Get the MCtruth handles and vectors
    std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

    //Then build a map from MCparticles to Hits and vice versa
    std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  mcParticleToHitsMap;
    std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  hitToMCParticleMap;

    //Apparrently a MCParticle doesn't know its origin (thanks Andy!)
    //I would also like a map from MCparticle to MCtruth and then I will be done.  and Vice Versa
    //Note which map is which!       //First  is one-to-many.         //Second is one-to-one
    std::map< art::Ptr<simb::MCTruth>,    std::vector<art::Ptr<simb::MCParticle>>>  MCTruthToMCParticlesMap;
    std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>                  MCParticleToMCTruthMap;
    std::map<int, art::Ptr<simb::MCParticle> >                                     MCParticleToTrackIdMap;

    std::vector<art::Ptr<sim::MCTrack>> mcTrackVector;
    std::vector<art::Ptr<sim::MCShower>> mcShowerVector;

    std::vector<art::Ptr<simb::MCParticle>> matchedMCParticleVector;
    std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > trackToMCParticleMap;
    std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > showerToMCParticleMap;

    //Given a simb::MCParticle we would like a map to either a sim::MCTrack or sim::MCShower
    std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCTrack> > MCParticleToMCTrackMap;
    std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCShower> > MCParticleToMCShowerMap;

    if(g_is_verbose){
      std::cout << "SinglePhoton::analyze()\t||\t Consolidated event summary:" << "\n";
      std::cout << "SinglePhoton::analyze()\t||\t - Number of primary cosmic-ray PFParticles   : " << crParticles.size() << "\n";
      std::cout << "SinglePhoton::analyze()\t||\t - Number of neutrino final-state PFParticles : " << nuParticles.size() << "\n";
      std::cout << "SinglePhoton::analyze()\t||\t    ... of which are track-like   : " << tracks.size() << "\n";
      std::cout << "SinglePhoton::analyze()\t||\t    ... of which are showers-like : " << showers.size() << "\n";
    }

    //**********************************************************************************************/
    //**********************************************************************************************/

    //and now get the simb::MCparticle to both MCtrack and MCshower maps (just for the MCparticles matched ok).
    //        if(!paras.s_run_pi0_filter) badChannelMatching<art::Ptr<recob::Track>>(badChannelVector, tracks, trackToNuPFParticleMap, pfParticleToHitsMap,geom,bad_channel_list_fixed_mcc9);
    //        badChannelMatching<art::Ptr<recob::Track>>(badChannelVector, tracks, trackToNuPFParticleMap, pfParticleToHitsMap,geom,bad_channel_list_fixed_mcc9);


    //*******************************Slices***************************************************************/

    //these are all filled in analyze slice
    //        std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind
    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> primaryPFPSliceIdVec; //stores a pair of only the primary PFP's in the event and the slice ind
    std::map<int, double> sliceIdToNuScoreMap; //map between a slice Id and neutrino score
    std::map<art::Ptr<recob::PFParticle>, bool> PFPToClearCosmicMap; //returns true for clear cosmic, false otherwise
    std::map<art::Ptr<recob::PFParticle>,bool> PFPToNuSliceMap;
    std::map<art::Ptr<recob::PFParticle>,double> PFPToTrackScoreMap;
    std::map<int, int> sliceIdToNumPFPsMap;

    if(g_is_verbose) std::cout<<"\nSinglePhoton::AnalyzeSlice()\t||\t Starting"<<std::endl;



    //******************************* CRT CRT***************************************************************/
    //Optical Flashes
    art::ValidHandle<std::vector<recob::OpFlash>> const & flashHandle  = evt.getValidHandle<std::vector<recob::OpFlash>>(paras.s_flashLabel);
    std::vector<art::Ptr<recob::OpFlash>> flashVector;
    art::fill_ptr_vector(flashVector,flashHandle);

    std::map<art::Ptr<recob::OpFlash>, std::vector< art::Ptr<sbn::crt::CRTHit>>> crtvetoToFlashMap;

    if(paras.s_runCRT){
      art::FindManyP<sbn::crt::CRTHit> crtveto_per_flash(flashHandle, evt,   paras.s_CRTVetoLabel);
      for(size_t i=0; i< flashVector.size(); ++i){
        crtvetoToFlashMap[flashVector[i]] = crtveto_per_flash.at(flashVector[i].key());
      }
    }

    art::Handle<std::vector<sbn::crt::CRTHit>> crthit_h; //only filled when there are hits, otherwise empty
    //art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;//Keng, this is to track CRT hit time;
    double evt_timeGPS_nsec = -999 ;
    if(paras.s_runCRT){
      //evt.getByLabel(paras.s_DAQHeaderProducer, rawHandle_DAQHeader);


      //       raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
      //Keng, time() is set to the trigger time.
      art::Timestamp evtTimeGPS = evt.time();
      evt_timeGPS_nsec = evtTimeGPS.timeLow();

      evt.getByLabel(paras.s_CRTHitProducer, crthit_h);
      std::cout<<"SinglePhoton::analyze \t||\t Got CRT hits"<<std::endl;
    }

    //Analize the CRT flashes and store them in the vertex_tree.
    //Found in analyze_OpFlashes.h
    ResizeFlashes(flashVector.size(), vars);
    AnalyzeFlashes(flashVector, crthit_h, evt_timeGPS_nsec, crtvetoToFlashMap, vars, paras);


    //******************************* Common Optical Filter **************************************************************/
    //Raw Optical fltr

    //        art::Handle<uboone::UbooneOpticalFilter> uBooNE_common_optFltr;
    //        if(evt.getByLabel("opfiltercommon", uBooNE_common_optFltr)){
    //            vars.m_flash_optfltr_pe_beam     = uBooNE_common_optFltr->PE_Beam();
    //            vars.m_flash_optfltr_pe_beavars.m_tot = uBooNE_common_optFltr->PE_Beavars.m_Total();
    //            vars.m_flash_optfltr_pe_veto     = uBooNE_common_optFltr->PE_Veto();
    //            vars.m_flash_optfltr_pe_veto_tot = uBooNE_common_optFltr->PE_Veto_Total();
    //        }else{
    //            vars.m_flash_optfltr_pe_beam     = -999;
    //            vars.m_flash_optfltr_pe_beavars.m_tot = -999;
    //            vars.m_flash_optfltr_pe_veto     = -999;
    //            vars.m_flash_optfltr_pe_veto_tot = -999;
    //            std::cout<<"No opfiltercommon product:"<<std::endl;
    //        }


    //*******************************  Tracks **************************************************************/
    std::cout<<"\nSinglePhoton::analyze \t||\t Start on Track Analysis "<<std::endl;

    //found in analyze_Tracks.h
    vars.m_reco_asso_tracks=tracks.size();
    ResizeTracks(vars.m_reco_asso_tracks, vars);
    AnalyzeTracks(
        allPPFParticles, 
        tracks, 
        pfParticleToSpacePointsMap,  MCParticleToTrackIdMap, sliceIdToNuScoreMap,
        vars, paras);
    AnalyzeTrackCalo(tracks, allPPFParticles, vars, paras);


    //Run over PID?
    if(paras.s_use_PID_algorithms)  CollectPID(tracks, allPPFParticles, vars);


    //*******************************  Showers **************************************************************/
    std::cout<<"\nSinglePhoton::analyze \t||\t Start on Shower Analysis "<<std::endl;


    //found in analyze_Showers.h
    vars.m_reco_asso_showers=showers.size();
    ResizeShowers(vars.m_reco_asso_showers, vars);
    AnalyzeShowers(allPPFParticles, showers,
        clusterToHitsMap,
        trigger_offset(detClocks), 
        theDetector, 
        vars,
        paras);

    std::cout<<"\nSinglePhoton::analyze \t||\t Finish Shower Analysis "<<std::endl;

    std::cout<<std::endl;
    std::cout<<"SinglePhoton::Reconstruction Summary ---------------------------------------------- "<<std::endl;
    std::vector<int> spacers = Printer_header({" pfpID", "           Parent        ", "Vertex(x,  ","   y,      ",",      z  )", "slice ","    nu_score","  trk_score "," prim? "," nu? "," nuSlice?"," cos? "," S# "," T# "});
    for(size_t jndex=0; jndex< allPPFParticles.size(); ++jndex){
      PandoraPFParticle temp_pf = allPPFParticles[jndex];
      art::Ptr<recob::PFParticle> pfp = temp_pf.pPFParticle;

      Printer_content( 
          {std::to_string(pfp->Self()), 
          std::to_string(pfp->Parent()) , 
          std::to_string(temp_pf.get_Vertex_pos()[0]),
          std::to_string(temp_pf.get_Vertex_pos()[1]),
          std::to_string(temp_pf.get_Vertex_pos()[2]),
          std::to_string(temp_pf.get_SliceID()),
          std::to_string(temp_pf.get_NuScore()),
          std::to_string(temp_pf.get_TrackScore()),
          (pfp->IsPrimary())?      "T":" ", 
          (temp_pf.get_IsNeutrino() )?    "T":" ",
          (temp_pf.get_IsNuSlice() )?    "T":" ",
          (temp_pf.get_IsClearCosmic() )?  "T":" ", 
          (temp_pf.get_HasShower() >0)?    "1":" ",
          (temp_pf.get_HasTrack() >0)?    "1":" "
          }, spacers);

    }

    //KENG no KalmanShowers in SBND?
    //    this->AnalyzeKalmanShowers(showers,showerToNuPFParticleMap,pfParticlesToShowerKalmanMap, kalmanTrackToCaloMap, pfParticleToHitsMap, theDetector);


    //Some misc things thrown in here rather than in a proper helper function. TODO. fix
    //Calc a fake shower "end" distance. How to define an end distance? good question
    for(size_t i_shr = 0; i_shr<showers.size();i_shr++){
      const art::Ptr<recob::Shower> s = showers[i_shr];
      const art::Ptr<recob::PFParticle> pfp = showerToNuPFParticleMap[s];
      const std::vector< art::Ptr<recob::SpacePoint> > shr_spacepoints = pfParticleToSpacePointsMap[pfp];

      vars.m_reco_shower_end_dist_to_active_TPC[i_shr] = 99999;
      vars.m_reco_shower_end_dist_to_SCB[i_shr] = 99999;

      for(auto &sp: shr_spacepoints){
        std::vector<double> tmp_spt = {sp->XYZ()[0],sp->XYZ()[1] , sp->XYZ()[2]};
        vars.m_reco_shower_end_dist_to_active_TPC[i_shr] = std::min(vars.m_reco_shower_end_dist_to_active_TPC[i_shr], distToTPCActive(tmp_spt, paras));
        double tmo;
        distToSCB(tmo,tmp_spt, paras);
        vars.m_reco_shower_end_dist_to_SCB[i_shr] = std::min(vars.m_reco_shower_end_dist_to_SCB[i_shr],tmo);

        //This section runs for only 1 shower events for purpose of testing delta specifics 
        if(showers.size()==1){
          vars.m_reco_shower_spacepoint_x.push_back(sp->XYZ()[0]);
          vars.m_reco_shower_spacepoint_y.push_back(sp->XYZ()[1]);
          vars.m_reco_shower_spacepoint_z.push_back(sp->XYZ()[2]);
        }

      }
    }

    //*******************************   MCTruth  **************************************************************/
    // Collect all the hits. We will need these. Lets grab both the handle as well as a vector of art::Ptr as I like both. 
    art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(paras.s_hitfinderLabel); 
    std::vector<art::Ptr<recob::Hit>> hitVector;
    art::fill_ptr_vector(hitVector,hitHandle);


    //Grab the backtracker info for MCTruth Matching
    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcparticles_per_hit(hitHandle, evt, paras.s_hitMCParticleAssnsLabel);

    // MCTruth, MCParticle, MCNeutrino information all comes directly from GENIE.
    // MCShower and MCTrack come from energy depositions in GEANT4

    //Only run if its not data, i.e. MC events :)
    if(!paras.s_is_data){
//GTruth information, don't need this at the moment; 
//      if(!paras.s_is_textgen){
//        std::vector<art::Ptr<simb::GTruth>> gTruthVector;
//        std::cout<<"\n Get some GTruth info."<<std::endl;
//        // if Text is in the generator label, skip it. TODO this is a bit simple but works, maybe add a boolean
//        art::ValidHandle<std::vector<simb::GTruth>> const & gTruthHandle= evt.getValidHandle<std::vector<simb::GTruth>>(paras.s_generatorLabel);
//        art::fill_ptr_vector(gTruthVector,gTruthHandle);
//        if(g_is_verbose){
//          for(size_t p=0; p< gTruthVector.size();p++) std::cout<<gTruthVector[p]<<" "<<*gTruthVector[p]<<std::endl;
//        }
//        std::cout<<"done"<<std::endl;
//      }

      //get MCTruth (GENIE)
      art::ValidHandle<std::vector<simb::MCTruth>> const & mcTruthHandle= evt.getValidHandle<std::vector<simb::MCTruth>>(paras.s_generatorLabel);
      art::fill_ptr_vector(mcTruthVector,mcTruthHandle);

      //get MCPartilces (GEANT4)
      art::ValidHandle<std::vector<simb::MCParticle>> const & mcParticleHandle= evt.getValidHandle<std::vector<simb::MCParticle>>(paras.s_geantModuleLabel);
      art::fill_ptr_vector(mcParticleVector,mcParticleHandle);

      //Found inanalyze_Geant4.h
      //Currently just saves first 2 particles. TODO have a input list of things to save. Dont want to save everything!!

      if(g_is_verbose) std::cout<<"SinglePhoton::AnalyzeGeant4()\t||\t Begininning recob::Geant4 analysis suite"<<std::endl;;
      AnalyzeGeant4(mcParticleVector, vars);


      //Get the MCParticles (move to do this ourselves later)
      CollectMCParticles(evt, paras.s_geantModuleLabel, MCTruthToMCParticlesMap, MCParticleToMCTruthMap, MCParticleToTrackIdMap,vars);


      //mcc9 march miniretreat fix
      std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing
      vars.m_test_matched_hits = 0;

      for(size_t j=0; j<hitVector.size();j++){
        const art::Ptr<recob::Hit> hit = hitVector[j];

        particle_vec.clear(); match_vec.clear(); //only store per hit
        mcparticles_per_hit.get(hit.key(), particle_vec, match_vec);
        if(particle_vec.size() > 0){
          vars.m_test_matched_hits++;
        }
      }


      //Important map, given a MCparticle, whats the "hits" associated
      BuildMCParticleHitMaps(evt, paras.s_geantModuleLabel, hitVector,  mcParticleToHitsMap, hitToMCParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters,  MCParticleToTrackIdMap, vars);


      //The recoMC was originally templated for any track shower, but sufficient differences in showers emerged to have stand alone sadly
      std::cout<<"\nSinglePhoton\t||\t Starting backtracker on recob::track"<<std::endl;
      std::vector<double> trk_overlay_vec = trackRecoMCmatching( tracks, 
          trackToMCParticleMap, 
          trackToNuPFParticleMap, 
          pfParticleToHitsMap, 
          mcparticles_per_hit, 
          matchedMCParticleVector,
          vars);

      std::cout<<"\nSinglePhoton\t||\t Starting backtracker on recob::shower"<<std::endl;
      if(g_is_verbose) std::cout<<"Matching "<<showers.size()<<" showers"<<std::endl;
      showerRecoMCmatching(
          allPPFParticles,
          showers, 
          showerToMCParticleMap, 
          mcparticles_per_hit, 
          matchedMCParticleVector, 
          MCParticleToTrackIdMap,
          vars);

      //photoNuclearTesting(matchedMCParticleVector);
      std::cout<<"\nSinglePhoton\t||\tStarting RecoMCTracks "<<std::endl;
      RecoMCTracks(
          allPPFParticles,
          tracks, 
          trackToMCParticleMap, 
          MCParticleToMCTruthMap,
          mcParticleVector, 
          MCParticleToTrackIdMap, 
          trk_overlay_vec,
          vars);

      std::cout<<"\nSinglePhoton\t||\tStarting AnalyzeMCTruths"<<std::endl;

      ResizeMCTruths(mcTruthVector[0]->NParticles(), vars);//assume only one MCTruth 
      AnalyzeMCTruths(mcTruthVector, mcParticleVector, vars, paras);


      if(!paras.s_is_textgen){// if Text is in the generator label, skip it. wont contain any

        //This is important for rebuilding MCFlux and GTruth to do eventweight 
        std::cout<<"\nSinglePhoton\t||\tStarting AnalyzeEventWeight"<<std::endl;
        AnalyzeEventWeight(evt, vars);


        std::cout<<"\nSinglePhoton\t||\tStarting AnalyzeRecoMCSlices"<<std::endl;
        AnalyzeRecoMCSlices(paras.s_truthmatching_signaldef, allPPFParticles,  MCParticleToTrackIdMap,  showerToMCParticleMap,  trackToMCParticleMap, vars, paras);

        std::cout<<"\nSinglePhoton\t||\tFinish AnalyzeRecoMCSlices"<<std::endl;

        //*******************************   PhotoNuclear Absorption  **************************************************************/
        vars.m_photonu_weight_low = -999;
        vars.m_photonu_weight_high = -999;
        if(paras.s_runPhotoNuTruth){
          art::Handle<std::vector<evwgh::MCEventWeight>>  ev_evw_ph ;
          if(    evt.getByLabel("eventweight",ev_evw_ph)){
            std::map<std::string, std::vector<double>> const & weight_map = ev_evw_ph->front().fWeight;
            for (auto const& x : weight_map){
              std::cout << x.first  // string (key)
                << ':' 
                << x.second.size() << std::endl ;
              if(x.first == "photonuclear_photon_PhotoNuclear"){
                auto vec  = x.second;
                double ph_low = vec[1];
                double ph_high = vec[0];
                std::cout<<"PhotoNuBit: "<<ph_low<<" "<<ph_high<<std::endl;
                vars.m_photonu_weight_low = ph_low;
                vars.m_photonu_weight_high = ph_high;
              }
            }

          }
        }


        //*******************************   True EventWeight  **************************************************************/
        //Remove this
        //        if(paras.s_runTrueEventweight){
        //
        //          art::ValidHandle<std::vector<evwgh::MCEventWeight>> const & ev_evw_true =  evt.getValidHandle<std::vector<evwgh::MCEventWeight>>(paras.s_true_eventweight_label);
        //          std::map<std::string, std::vector<double>> const & weight_map = ev_evw_true->front().fWeight;
        //          if(ev_evw_true->size() > 1) {
        //            std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
        //              << "WARNING: eventweight has more than one entry\n";
        //          }
        //          fmcweight=weight_map;
        //        }
        //
      }//end NOT textgen 



      std::cout<<"SinglePhoton::analyze\t||\t finnished loop for this event"<<std::endl;
    }//end MC loop

    //*******************************   Isolation (SSV precursor)  **************************************************************/
    if(!paras.s_run_all_pfps && ! paras.s_run_pi0_filter){ 
      std::cout<<"\nSinglePhoton::analyze \t||\t Start IsolationStudy "<<std::endl;
      IsolationStudy(allPPFParticles, tracks,  showers, theDetector, vars, paras);
      std::cout<<"\nSinglePhoton::analyze \t||\t Finish IsolationStudy "<<std::endl;
    }



    // ################################################### SEAview SEAview #########################################################
    // ################################################### Proton Stub ###########################################
    // ------------- stub clustering ---------------------------
    std::cout << "----------------- Stub clustering --------------------------- " << std::endl;
    //    std::cout << "SEAview Stub formation: " << (paras.s_runSEAviewStub ? "true" : "false" ) << " nshower requirement: " << paras.s_SEAviewStubNumRecoShower << ", actual num shower: " << showers.size() << " | ntrack requirement: " << paras.s_SEAviewStubNumRecoTrack << ", actual num track: " << tracks.size() << std::endl;

    if(  !paras.s_run_pi0_filter && 
        paras.s_runSEAviewStub && 
        (paras.s_SEAviewStubNumRecoShower== -1 || (int)showers.size()== paras.s_SEAviewStubNumRecoShower) && 
        (paras.s_SEAviewStubNumRecoTrack == -1 || (int)tracks.size() == paras.s_SEAviewStubNumRecoTrack)){

      // grab all hits in the slice of the reco shower
      art::Ptr<recob::Shower> p_shr = showers.front();
      PandoraPFParticle* ppfp = PPFP_GetPPFPFromShower(allPPFParticles, p_shr);
      art::Ptr<recob::PFParticle> p_pfp = ppfp->pPFParticle;//showerToNuPFParticleMap[p_shr];
      std::vector<art::Ptr<recob::Hit>> p_hits = pfParticleToHitsMap[p_pfp];

      int p_sliceid = ppfp->get_SliceID();
      auto p_slice_hits =    sliceIDToHitsMap[p_sliceid];

      std::string uniq_tag = "HitThres_"+ std::to_string(static_cast<int>(paras.s_SEAviewStubHitThreshold)) + "_" + std::to_string(vars.m_run_number)+"_"+std::to_string(vars.m_subrun_number)+"_"+std::to_string(vars.m_event_number);

      //Setup seaviewr object
      ////CHECK undefined reference? Found in header, but not defined;
      //Solution: add MODULE_LIBRARIES in cmake;
      seaview::SEAviewer sevd("Stub_"+uniq_tag,paras.s_geom, theDetector);
      //Pass in any bad channels you like
      //            sevd.setBadChannelList(bad_channel_list_fixed_mcc9);
      //Give it a vertex to center around
      sevd.loadVertex(vars.m_vertex_pos_x,vars.m_vertex_pos_y, vars.m_vertex_pos_z);

      //Add the hits from just this slice, as well as hits within 150cm of the vertex 
      sevd.addHitsToConsider(hitVector);   // std::vector<art::Ptr<recob::Hit>> 
      sevd.filterConsideredHits(150); //remve hits that're not within 150cm of the vertex on 2D view
      sevd.addHitsToConsider(p_slice_hits);

      sevd.addAllHits(hitVector); // std::vector<art::Ptr<recob::Hit>> 
      sevd.setHitThreshold(paras.s_SEAviewStubHitThreshold); 


      //Add all the "nice " PFParticle Hits, as well as what to label
      //sevd.addPFParticleHits(p_hits, "Shower");  //std::vector<art::Ptr<recob::Hit>> and std::string
      sevd.addPFParticleHits(p_hits, "Shower", vars.m_reco_shower_energy_max[0], vars.m_reco_shower_conversion_distance[0], vars.m_reco_shower_impact_parameter[0]);  //std::vector<art::Ptr<recob::Hit>> and std::string

      //and add the SingleShower we like
      sevd.addShower(p_shr); // art::Ptr<recob::Shower>

      //Add all track PFP
      int i_trk = 0;
      for(auto &trk: tracks){
        art::Ptr<recob::PFParticle> p_pfp_trk = trackToNuPFParticleMap[trk];
        std::vector<art::Ptr<recob::Hit>> p_hits_trk = pfParticleToHitsMap[p_pfp_trk];
        //sevd.addPFParticleHits(p_hits_trk,"track");
        sevd.addPFParticleHits(p_hits_trk,"track", vars.m_reco_track_length[i_trk], vars.m_reco_track_spacepoint_principal0[i_trk]);
        sevd.addTrack(trk);
        ++i_trk;
      }

      //Add all cosmic-relatd PFP
      for(auto &cr: crParticles){
        std::vector<art::Ptr<recob::Hit>> p_hits_cr = cr_pfParticleToHitsMap[cr];
        sevd.addPFParticleHits(p_hits_cr,"cosmic");
      }


      //We then calculate Unassociated hits, i.e the hits not associated to the "Shower" or tracksyou passed in. 
      auto vnh= sevd.calcUnassociatedHits();
      vars.m_trackstub_num_unassociated_hits =vnh[1]+vnh[2];
      vars.m_trackstub_unassociated_hits_below_threshold = vnh[2];
      vars.m_trackstub_associated_hits = vnh[0]-vnh[1]-vnh[2];

      //Recluster, group unassociated hits into different clusters
      sevd.runseaDBSCAN(paras.s_SEAviewStubDbscanMinPts, paras.s_SEAviewStubDbscanEps);

      //And some plotting
      // If we want to plot pdfs again later, then we can't plot here
      //if(vars.m_SEAviewStubMakePDF) sevd.Print(vars.m_SEAviewStubPlotDistance);

      //Analyze formed clusters and save info
      std::vector<seaview::cluster> vec_SEAclusters ;
      sevd.analyzeTrackLikeClusters(paras.s_SEAviewStubDbscanEps, showerToNuPFParticleMap, pfParticleToHitsMap, vec_SEAclusters);


      //And save to file.
      std::cout<<"After SEAview we have "<<vec_SEAclusters.size()<<" Stub clusters to chat about"<<std::endl;

      vars.m_trackstub_num_candidates = 0;
      for(size_t c=0; c< vec_SEAclusters.size(); c++){
        auto& clu = vec_SEAclusters.at(c); //type: seaview::cluster
        int pl = clu.getPlane();
        auto hitz = clu.getHits();
        double Ep = CalcEShowerPlane(hitz,pl, paras); 
        int remerge = clu.getShowerRemerge();
        seaview::cluster_score * ssscorz = clu.getScore();

        std::cout<<c<<" "<<pl<<" "<<Ep<<" "<<clu.getMinHitImpactParam()<<" "<<clu.getMinHitConvDist()<<" "<<clu.getMinHitIOC()<<" "<<clu.getMeanADC()<<" "<<clu.getADCrms()<<" "<< clu.getLinearChi() << " " << remerge<<std::endl;

        //if the cluster is too close to the recob::shower, then do not include it
        if(remerge>=0 && remerge< (int)vars.m_reco_shower_reclustered_energy_plane2.size()){
          //decide not to add energy of the cluster to reco shower if it's matched
          //
          //if(pl==0)vars.m_reco_shower_reclustered_energy_plane0[remerge]+=Ep;
          //if(pl==1)vars.m_reco_shower_reclustered_energy_plane1[remerge]+=Ep;
          //if(pl==2)vars.m_reco_shower_reclustered_energy_plane2[remerge]+=Ep;

          continue;// Dont include this as a viable cluster!
        }

        ++vars.m_trackstub_num_candidates;
        //determine if this cluster is in neutrino slice
        vars.m_trackstub_candidate_in_nu_slice.push_back(clu.InNuSlice(sliceIDToHitsMap, p_sliceid));

        //Fill All the bits
        vars.m_trackstub_candidate_num_hits.push_back((int)hitz.size());
        vars.m_trackstub_candidate_num_wires.push_back((int)ssscorz->n_wires);
        vars.m_trackstub_candidate_num_ticks.push_back((int)ssscorz->n_ticks);
        vars.m_trackstub_candidate_plane.push_back(pl);
        vars.m_trackstub_candidate_PCA.push_back(ssscorz->pca_0);
        vars.m_trackstub_candidate_mean_tick.push_back(ssscorz->mean_tick);
        vars.m_trackstub_candidate_max_tick.push_back(ssscorz->max_tick);
        vars.m_trackstub_candidate_min_tick.push_back(ssscorz->min_tick);
        vars.m_trackstub_candidate_min_wire.push_back(ssscorz->min_wire);
        vars.m_trackstub_candidate_max_wire.push_back(ssscorz->max_wire);
        vars.m_trackstub_candidate_mean_wire.push_back(ssscorz->mean_wire);
        vars.m_trackstub_candidate_min_dist.push_back(ssscorz->min_dist);
        vars.m_trackstub_candidate_min_impact_parameter_to_shower.push_back(clu.getMinHitImpactParam());
        vars.m_trackstub_candidate_min_conversion_dist_to_shower_start.push_back(clu.getMinHitConvDist());
        vars.m_trackstub_candidate_min_ioc_to_shower_start.push_back(clu.getMinHitIOC());
        vars.m_trackstub_candidate_ioc_based_length.push_back(clu.getIOCbasedLength());
        vars.m_trackstub_candidate_wire_tick_based_length.push_back(clu.getWireTickBasedLength());
        vars.m_trackstub_candidate_mean_ADC_first_half.push_back(clu.getMeanADCFirstHalf());
        vars.m_trackstub_candidate_mean_ADC_second_half.push_back(clu.getMeanADCSecondHalf());
        vars.m_trackstub_candidate_mean_ADC_first_to_second_ratio.push_back(clu.getMeanADCRatio());
        vars.m_trackstub_candidate_track_angle_wrt_shower_direction.push_back(clu.getTrackAngleToShowerDirection());
        vars.m_trackstub_candidate_linear_fit_chi2.push_back(clu.getLinearChi());
        vars.m_trackstub_candidate_mean_ADC.push_back(clu.getMeanADC());
        vars.m_trackstub_candidate_ADC_RMS.push_back(clu.getADCrms());
        vars.m_trackstub_candidate_energy.push_back(Ep);
        vars.m_trackstub_candidate_remerge.push_back(remerge);


        //MCTruth matching for pi0's
        if(paras.s_is_data){
          vars.m_trackstub_candidate_matched.push_back(-1);
          vars.m_trackstub_candidate_pdg.push_back(-1);
          vars.m_trackstub_candidate_parent_pdg.push_back(-1);
          vars.m_trackstub_candidate_trackid.push_back(-1);
          vars.m_trackstub_candidate_true_energy.push_back(-1);
          vars.m_trackstub_candidate_overlay_fraction.push_back(-1);
          vars.m_trackstub_candidate_matched_energy_fraction_best_plane.push_back(-1);
        }else{

          auto ssmatched = SecondShowerMatching(hitz, mcparticles_per_hit, mcParticleVector,  MCParticleToTrackIdMap, vars);
          vars.m_trackstub_candidate_matched.push_back(ssmatched[0]);
          vars.m_trackstub_candidate_pdg.push_back(ssmatched[1]);
          vars.m_trackstub_candidate_parent_pdg.push_back(ssmatched[2]);
          vars.m_trackstub_candidate_trackid.push_back(ssmatched[3]);
          vars.m_trackstub_candidate_true_energy.push_back(ssmatched[4]);
          vars.m_trackstub_candidate_overlay_fraction.push_back(ssmatched[5]);
          vars.m_trackstub_candidate_matched_energy_fraction_best_plane.push_back(ssmatched[6]);

          //Guanqun: print out (best-matched) truth information of the cluster
          std::cout << "Cluster: " << vars.m_trackstub_num_candidates-1  << " plane: " << vars.m_trackstub_candidate_plane.back() << ", energy: " << vars.m_trackstub_candidate_energy.back() << ", min IOC of hit(wrt shower): " << vars.m_trackstub_candidate_min_ioc_to_shower_start.back() << "\n";
          std::cout << "Cluster is matched: " << vars.m_trackstub_candidate_matched.back() << ", matched PDG: " << vars.m_trackstub_candidate_pdg.back() << " track ID: " << vars.m_trackstub_candidate_trackid.back() << " overlay fraction: " << vars.m_trackstub_candidate_overlay_fraction.back() << std::endl; 
          std::cout << "===============================================================" << std::endl;
        }
        sevd.SetClusterLegend(c, vars.m_trackstub_candidate_energy.back(),  vars.m_trackstub_candidate_matched.back(), vars.m_trackstub_candidate_pdg.back() , vars.m_trackstub_candidate_overlay_fraction.back() );


      } //end of cluster loop

      // Plot the event
      if(paras.s_SEAviewStubMakePDF){
        sevd.Print(paras.s_SEAviewStubPlotDistance);
      }

      //group clusters HERE
      std::pair<int, std::pair<std::vector<std::vector<double>>, std::vector<double>> > group_result = GroupClusterCandidate(vars.m_trackstub_num_candidates,  vars.m_trackstub_candidate_plane, vars.m_trackstub_candidate_max_tick, vars.m_trackstub_candidate_min_tick);
      vars.m_trackstub_num_candidate_groups = group_result.first;
      vars.m_grouped_trackstub_candidate_indices = group_result.second.first;
      vars.m_trackstub_candidate_group_timeoverlap_fraction = group_result.second.second;
    }

    // --------------- shower clustering --------------------------
    std::cout << "------------- Shower clustering --------------------" << std::endl;
    std::cout << "SEAview Shower cluster formation: " << (paras.s_runSEAview ? "true" : "false" ) << " nshower requirement: " << paras.s_SEAviewNumRecoShower << ", actual num shower: " << showers.size() << " | ntrack requirement: " << paras.s_SEAviewNumRecoTrack << ", actual num track: " << tracks.size() << std::endl;

    if(!paras.s_run_pi0_filter &&  
        paras.s_runSEAview && 
        (paras.s_SEAviewNumRecoShower == -1 || (int)showers.size()== paras.s_SEAviewNumRecoShower) && 
        (paras.s_SEAviewNumRecoTrack == -1 || (int)tracks.size() == paras.s_SEAviewNumRecoTrack) ){

      art::Ptr<recob::Shower> p_shr = showers.front();

      PandoraPFParticle* ppfp = PPFP_GetPPFPFromShower(allPPFParticles, p_shr);
      art::Ptr<recob::PFParticle> p_pfp = ppfp->pPFParticle;//showerToNuPFParticleMap[p_shr];
      std::vector<art::Ptr<recob::Hit>> p_hits = pfParticleToHitsMap[p_pfp];


      int p_sliceid = ppfp->get_SliceID();
      auto p_slice_hits =    sliceIDToHitsMap[p_sliceid];

      std::string uniq_tag = "HitThres_"+ std::to_string(static_cast<int>(paras.s_SEAviewHitThreshold)) + "_" + std::to_string(vars.m_run_number)+"_"+std::to_string(vars.m_subrun_number)+"_"+std::to_string(vars.m_event_number);

      //Setup seaviewr object
      seaview::SEAviewer sevd("Shower_"+uniq_tag, paras.s_geom, theDetector );
      //Pass in any bad channels you like
      //            sevd.setBadChannelList(bad_channel_list_fixed_mcc9);
      //Give it a vertex to center around
      sevd.loadVertex(vars.m_vertex_pos_x,vars.m_vertex_pos_y, vars.m_vertex_pos_z);

      //Add hits to consider for clustering 
      //sevd.addHitsToConsider(hitVector);// DONT do this yet, need something smarter for SSV
      //sevd.filterConsideredHits(150); 
      sevd.addHitsToConsider(p_slice_hits);

      //Add all hits in the events
      sevd.addAllHits(hitVector); // std::vector<art::Ptr<recob::Hit>> 
      sevd.setHitThreshold(paras.s_SEAviewHitThreshold); 

      //Add all the "nice " PFParticle Hits, as well as what to label
      //sevd.addPFParticleHits(p_hits, "Shower");  //std::vector<art::Ptr<recob::Hit>> and std::string
      sevd.addPFParticleHits(p_hits, "Shower", vars.m_reco_shower_energy_max[0], vars.m_reco_shower_conversion_distance[0], vars.m_reco_shower_impact_parameter[0]);  //std::vector<art::Ptr<recob::Hit>> and std::string

      //and add the SingleShower we like
      sevd.addShower(p_shr); // art::Ptr<recob::Shower>

      //Add all track PFP
      int i_trk = 0;
      for(auto &trk: tracks){
        art::Ptr<recob::PFParticle> p_pfp_trk = trackToNuPFParticleMap[trk];
        std::vector<art::Ptr<recob::Hit>> p_hits_trk = pfParticleToHitsMap[p_pfp_trk];
        //sevd.addPFParticleHits(p_hits_trk,"track");
        sevd.addPFParticleHits(p_hits_trk,"track", vars.m_reco_track_length[i_trk], vars.m_reco_track_spacepoint_principal0[i_trk]);
        sevd.addTrack(trk);
        ++i_trk;
      }

      //Add all cosmic-relatd PFP // DONT do this yet, see line 1206
      /*for(auto &cr: crParticles){
        std::vector<art::Ptr<recob::Hit>> p_hits_cr = cr_pfParticleToHitsMap[cr];
        sevd.addPFParticleHits(p_hits_cr,"cosmic");
        }
        */

      //We then calculate Unassociated hits, i.e the hits not associated to the "Shower" or tracksyou passed in. 
      auto vnh= sevd.calcUnassociatedHits();
      vars.m_sss_num_unassociated_hits =vnh[1]+vnh[2];
      vars.m_sss_num_unassociated_hits_below_threshold = vnh[2];
      vars.m_sss_num_associated_hits = vnh[0]-vnh[1]-vnh[2];

      //Recluster, group unassociated hits into different clusters
      sevd.runseaDBSCAN(paras.s_SEAviewDbscanMinPts, paras.s_SEAviewDbscanEps);


      //This is the place I will put the new Second Shower Search
      std::vector<seaview::cluster> vec_SEAclusters ;
      sevd.analyzeShowerLikeClusters(paras.s_SEAviewDbscanEps, showerToNuPFParticleMap, pfParticleToHitsMap, vec_SEAclusters);

      //And save to file.
      std::cout<<"After SEAview we have "<<vec_SEAclusters.size()<<" Shower clusters to chat about"<<std::endl;

      vars.m_sss_num_candidates = 0;
      for(size_t c=0; c< vec_SEAclusters.size(); c++){
        auto clu = vec_SEAclusters.at(c); //type: seaview::cluster
        int pl = clu.getPlane();
        auto hitz = clu.getHits();
        double Ep = CalcEShowerPlane(hitz,pl, paras); 
        int remerge = clu.getShowerRemerge();
        seaview::cluster_score * ssscorz = clu.getScore();

        std::cout<<c<<" "<<pl<<" "<<Ep<<" "<<clu.getImpactParam()<<" "<<clu.getFitSlope()<<" "<<clu.getFitCons()<<" "<<clu.getMeanADC() << " " << clu.getADCrms() << " "<<clu.getAngleWRTShower()<<" "<<remerge<<std::endl;

        //if the cluster is too close to the recob::shower, then do not include it
        if(remerge>=0 && remerge< (int)vars.m_reco_shower_reclustered_energy_plane2.size()){
          if(pl==0)vars.m_reco_shower_reclustered_energy_plane0[remerge]+=Ep;
          if(pl==1)vars.m_reco_shower_reclustered_energy_plane1[remerge]+=Ep;
          if(pl==2)vars.m_reco_shower_reclustered_energy_plane2[remerge]+=Ep;

          continue;// Dont include this as a viable cluster!
        }

        ++vars.m_sss_num_candidates;

        //determine if this cluster is in neutrino slice
        vars.m_sss_candidate_in_nu_slice.push_back(clu.InNuSlice(sliceIDToHitsMap, p_sliceid));

        //Fill All the bits
        vars.m_sss_candidate_num_hits.push_back((int)hitz.size());
        vars.m_sss_candidate_num_wires.push_back((int)ssscorz->n_wires);
        vars.m_sss_candidate_num_ticks.push_back((int)ssscorz->n_ticks);
        vars.m_sss_candidate_plane.push_back(pl);
        vars.m_sss_candidate_PCA.push_back(ssscorz->pca_0);
        vars.m_sss_candidate_impact_parameter.push_back(clu.getImpactParam());
        vars.m_sss_candidate_fit_slope.push_back(clu.getFitSlope());
        vars.m_sss_candidate_fit_constant.push_back(clu.getFitCons());
        vars.m_sss_candidate_mean_tick.push_back(ssscorz->mean_tick);
        vars.m_sss_candidate_max_tick.push_back(ssscorz->max_tick);
        vars.m_sss_candidate_min_tick.push_back(ssscorz->min_tick);
        vars.m_sss_candidate_min_wire.push_back(ssscorz->min_wire);
        vars.m_sss_candidate_max_wire.push_back(ssscorz->max_wire);
        vars.m_sss_candidate_mean_wire.push_back(ssscorz->mean_wire);
        vars.m_sss_candidate_min_dist.push_back(ssscorz->min_dist);
        vars.m_sss_candidate_wire_tick_based_length.push_back(clu.getWireTickBasedLength());
        vars.m_sss_candidate_mean_ADC.push_back(clu.getMeanADC());
        vars.m_sss_candidate_ADC_RMS.push_back(clu.getADCrms());
        vars.m_sss_candidate_energy.push_back(Ep);
        vars.m_sss_candidate_angle_to_shower.push_back(clu.getAngleWRTShower());
        vars.m_sss_candidate_remerge.push_back(remerge);


        //MCTruth matching for pi0's
        if(paras.s_is_data){
          vars.m_sss_candidate_matched.push_back(-1);
          vars.m_sss_candidate_pdg.push_back(-1);
          vars.m_sss_candidate_parent_pdg.push_back(-1);
          vars.m_sss_candidate_trackid.push_back(-1);
          vars.m_sss_candidate_true_energy.push_back(-1);
          vars.m_sss_candidate_overlay_fraction.push_back(-1);
          vars.m_sss_candidate_matched_energy_fraction_best_plane.push_back(-1);
        }else{

          auto ssmatched = SecondShowerMatching(hitz, mcparticles_per_hit, mcParticleVector,  MCParticleToTrackIdMap, vars);
          vars.m_sss_candidate_matched.push_back(ssmatched[0]);
          vars.m_sss_candidate_pdg.push_back(ssmatched[1]);
          vars.m_sss_candidate_parent_pdg.push_back(ssmatched[2]);
          vars.m_sss_candidate_trackid.push_back(ssmatched[3]);
          vars.m_sss_candidate_true_energy.push_back(ssmatched[4]);
          vars.m_sss_candidate_overlay_fraction.push_back(ssmatched[5]);
          vars.m_sss_candidate_matched_energy_fraction_best_plane.push_back(ssmatched[6]);

          //Guanqun: print out (best-matched) truth information of the cluster
          std::cout << "Cluster: " << vars.m_sss_num_candidates-1  << " plane: " << vars.m_sss_candidate_plane.back() << ", energy: " << vars.m_sss_candidate_energy.back() << "\n";
          std::cout << "Cluster is matched: " << vars.m_sss_candidate_matched.back() << ", matched PDG: " << vars.m_sss_candidate_pdg.back() << " track ID: " << vars.m_sss_candidate_trackid.back() << " overlay fraction: " << vars.m_sss_candidate_overlay_fraction.back() << std::endl; 
          std::cout << "===============================================================" << std::endl;
        }


        sevd.SetClusterLegend(c, vars.m_sss_candidate_energy.back(),  vars.m_sss_candidate_matched.back(), vars.m_sss_candidate_pdg.back() , vars.m_sss_candidate_overlay_fraction.back() );

      } //end of cluster loop

      // Plot the event
      if(paras.s_SEAviewMakePDF){
        sevd.Print(paras.s_SEAviewPlotDistance);
      }

    }

    for(int i =0; i<(int)showers.size(); i++){
      vars.m_reco_shower_reclustered_energy_max[i] = std::max(vars.m_reco_shower_reclustered_energy_plane1[i],std::max(vars.m_reco_shower_reclustered_energy_plane0[i],vars.m_reco_shower_reclustered_energy_plane2[i]));
    }

    // ################################################### END SEAview END SEAview #########################################################
    // #####################################################################################################################################



    // PandoraAllOutComes
    // I.e This runs over all 3D reco showers in the whole event and find second shower candidates
    if(!paras.s_run_pi0_filter){
      std::cout<<"------------ Shower3D --------------"<<std::endl;
      /*for(auto &s : showers){
        std::cout<<"shower pfp key : "<<showerToNuPFParticleMap[s].key()<<" self: "<<showerToNuPFParticleMap[s]->Self()<<std::endl;
        }
        for(auto &s : tracks){
        std::cout<<"track pfp key : "<<trackToNuPFParticleMap[s].key()<<" self: "<<trackToNuPFParticleMap[s]->Self()<<std::endl;
        }*/

      SecondShowerSearch3D(showers, showerToNuPFParticleMap, tracks,trackToNuPFParticleMap,evt, vars, paras);


      //And cluster the 2d and 3d second showers. Very simple TODO
      SimpleSecondShowerCluster( vars,paras);

    }



    //---------------------- END OF LOOP, fill vertex ---------------------
    bool filter_pass_2g1p = Pi0PreselectionFilter(vars,paras);
    bool filter_pass_2g0p = Pi0PreselectionFilter2g0p(vars,paras);

    if (paras.s_fill_trees &&  (  (filter_pass_2g1p && paras.s_run_pi0_filter_2g1p) || (filter_pass_2g0p && paras.s_run_pi0_filter_2g0p) || !paras.s_run_pi0_filter ) ) {
      vars.vertex_tree->Fill();
      vars.ncdelta_slice_tree->Fill();
      vars.eventweight_tree->Fill();
      vars.true_eventweight_tree->Fill();
      vars.geant4_tree->Fill();
    }

    //Rest the vertex after filling
    //        this->ClearMeta();
    if(paras.s_run_pi0_filter_2g1p)  return filter_pass_2g1p;
    else if(paras.s_run_pi0_filter_2g0p)  return filter_pass_2g0p;

    //if not in filter mode pass all
    return true;

  }// end filter_module main class



  //-------------------------------------------------------------------------------------------

  //This runs ONCE at the start of the job and sets up all the necessary services and TTrees
  void SinglePhoton::beginJob()
  {
    mf::LogDebug("SinglePhoton") << " *** beginJob() *** " << "\n";

    art::ServiceHandle<art::TFileService> tfs;//output ROOT 

    vars.run_subrun_tree = tfs->make<TTree>("run_subrun_tree","run_subrun_tree");
    vars.true_eventweight_tree = tfs->make<TTree>("true_eventweight_tree", "true_eventweight_tree");
    vars.pot_tree = tfs->make<TTree>("pot_tree", "pot_tree");

    vars.eventweight_tree = tfs->make<TTree>("eventweight_tree", "eventweight_tree");
    vars.ncdelta_slice_tree = tfs->make<TTree>("ncdelta_slice_tree", "ncdelta_slice_tree");
    vars.geant4_tree = tfs->make<TTree>("geant4_tree","geant4_tree");
    vars.vertex_tree = tfs->make<TTree>("vertex_tree", "vertex_tree");

    //run_subrun_tree, reset some POT
    vars.m_run = 0;
    vars.m_subrun = 0;
    vars.m_subrun_pot = 0;

    // --------------------- POT Releated variables -----------------
    vars.m_number_of_events = 0;
    vars.m_number_of_vertices = 0;
    vars.m_pot_count=0;
    vars.m_pot_per_event = 0;
    vars.m_pot_per_subrun = 0;
    vars.m_number_of_events_in_subrun=0;


    CreateMetaBranches(vars);
    CreateIsolationBranches(vars);
    CreateSecondShowerBranches(vars);
    CreateSecondShowerBranches3D(vars);
    CreateStubBranches(vars);
    CreateFlashBranches(vars);
    CreateShowerBranches(vars);
    CreateMCTruthBranches(vars);
    CreateTrackBranches(vars);

    CreateEventWeightBranches(vars);
    CreateSliceBranches(vars);
    CreateGeant4Branches(vars);


    //------------------- List of Selected Events to run --------
    if(paras.s_runSelectedEvent){
      std::cout << "SinglePhoton \t||\t Running in selected-event only mode " << std::endl;

      std::ifstream infile(paras.s_selected_event_list);
      if(!infile){
        std::cerr << "Fail to open file: " <<paras.s_selected_event_list << std::endl;
        return;
      }

      //read from file, run number, subrun number ,event number that should be run
      vars.m_selected_set.clear();
      std::string line;
      while(std::getline(infile, line)){
        std::istringstream ss(line);

        std::vector<int> event_info; 
        for(int i; ss >> i; ) event_info.push_back(i);

        vars.m_selected_set.insert(event_info);
      }

      infile.close();

      if(g_is_verbose){
        std::cout << "Selected Events: " << std::endl;
        std::cout << "Run \t SubRun \t Event" << std::endl;
        for(auto & v: vars.m_selected_set){
          std::for_each(v.begin(), v.end(), [](int n){std::cout << n<<" \t "; });
          std::cout << std::endl;
        }
      }
    }
  }



  bool SinglePhoton::beginSubRun(art::SubRun& sr) {

    vars.m_run = sr.run();
    vars.m_subrun = sr.subRun();

    double this_pot = 0;

    //reset subrun count 
    vars.m_subrun_counts = 0;


    if(paras.s_potLabel != ""){
      if(paras.s_potLabel == "generator"){

        art::Handle<sumdata::POTSummary> gen_pot_hand;
        if(sr.getByLabel(paras.s_potLabel,gen_pot_hand)){
          this_pot =  gen_pot_hand->totgoodpot;
          vars.m_pot_count += this_pot;
          std::cout<<"SinglePhoton::beginSubRun()\t||\t SubRun POT: "<<this_pot<<" . Current total POT this file: "<<vars.m_pot_count<<" (label) "<<paras.s_potLabel<<std::endl;
        }
      }else{

        art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR875;
        if (sr.getByLabel("beamdata","bnbETOR875",potSummaryHandlebnbETOR875)){
          this_pot =potSummaryHandlebnbETOR875->totpot; 
          vars.m_pot_count += this_pot;
          std::cout<<"SinglePhoton::beginSubRun()\t||\t SubRun POT: "<<potSummaryHandlebnbETOR875->totpot<<" . Current total POT this file: "<<vars.m_pot_count<<" (label) "<<paras.s_potLabel<<std::endl;
        }
      }
    }
    vars.m_subrun_pot = this_pot; 
    return true;
  }


  bool SinglePhoton::endSubRun(art::SubRun&sr){

    vars.run_subrun_tree->Fill();
    return true;
  }

  void SinglePhoton::endJob()
  {
    //        if (vars.m_print_out_event){
    //            out_stream.close();
    //        }
    vars.pot_tree->Fill();
  }

  DEFINE_ART_MODULE(SinglePhoton)
} //namespace single_photon
