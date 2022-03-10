#include "SinglePhoton_module.h"
#include "analyze_OpFlashes.h"
#include "analyze_Tracks.h"
#include "analyze_Showers.h"
#include "analyze_Template.h"
#include "analyze_MCTruth.h"
#include "analyze_EventWeight.h"
#include "analyze_Slice.h"
#include "second_shower_search.h"
#include "isolation.h"
#include "BobbyVertexBuilder.h"
#include "Atlas.h"

using namespace std;

namespace single_photon
{
	//constructor
    SinglePhoton::SinglePhoton(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
    {
        this->reconfigure(pset);
        theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
        detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
        SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
        geom = lar::providerFrom<geo::Geometry>();

    }

    void SinglePhoton::reconfigure(fhicl::ParameterSet const &pset)
    {
		//	cout<<"\n\n\n\n\n\n Reconfigure!"<<endl;
        m_print_out_event = pset.get<bool>("PrintOut", true);
        m_is_verbose = pset.get<bool>("Verbose",false);
        m_use_PID_algorithms = pset.get<bool>("usePID",false);
        m_use_delaunay = pset.get<bool>("useDelaunay",false);
        m_is_data = pset.get<bool>("isData",false);
        m_is_overlayed = pset.get<bool>("isOverlayed",false);

        m_pandoraLabel = pset.get<std::string>("PandoraLabel");
        m_trackLabel = pset.get<std::string>("TrackLabel");
        m_sliceLabel = pset.get<std::string>("SliceLabel","pandora");
        m_showerLabel = pset.get<std::string>("ShowerLabel");
        m_caloLabel = pset.get<std::string>("CaloLabel");
        m_flashLabel = pset.get<std::string>("FlashLabel");
        m_potLabel = pset.get<std::string>("POTLabel");
        m_beamgate_flash_start = pset.get<double>("beamgateStartTime",3.2); //Defaults to MC for now. Probably should change
        m_beamgate_flash_end = pset.get<double>("beamgateEndTime",4.8);
        m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
        m_badChannelLabel = pset.get<std::string>("BadChannelLabel","badmasks");
        // m_badChannelProducer = pset.get<std::string>("BadChannelProducer","nfspl1");
        m_badChannelProducer = pset.get<std::string>("BadChannelProducer","simnfspl1");
        m_showerKalmanLabel = pset.get<std::string>("ShowerTrackFitter","pandoraKalmanShower");
        m_showerKalmanCaloLabel =  pset.get<std::string>("ShowerTrackFitterCalo","pandoraKalmanShowercali");

        m_generatorLabel = pset.get<std::string>("GeneratorLabel","generator");
        m_mcTrackLabel = pset.get<std::string>("MCTrackLabel","mcreco");
        m_mcShowerLabel = pset.get<std::string>("MCShowerLabel","mcreco");
        m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
        m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","gaushitTruthMatch");
        m_hitMCParticleAssnsLabel = pset.get<std::string>("HitMCParticleAssnLabel","gaushitTruthMatch");


        m_CRTTzeroLabel = pset.get<std::string>("CRTTzeroLabel","crttzero");
		m_runCRT = pset.get<bool>("runCRT",false);
		m_CRTHitProducer = pset.get<std::string>("CRTHitProducer", "crthitcorr");

		m_DTOffset = pset.get<double>("DTOffset" , 68600.); //us, taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter.fcl
		m_Resolution = pset.get<double>("Resolution" ,  1.0); //us, taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter.fcl
		m_DAQHeaderProducer = pset.get<std::string>("DAQHeaderProducer" ,  "daq");

		//Some track calorimetry parameters
        m_track_calo_min_dEdx = pset.get<double>("Min_dEdx",0.005);
        m_track_calo_max_dEdx = pset.get<double>("Max_dEdx", 30);
        m_track_calo_min_dEdx_hits = pset.get<double>("Min_dEdx_hits",5); //might be good?
        m_track_calo_trunc_fraction = pset.get<double>("TruncMeanFraction",20.0);

        //Some shower calorimetry parameters
        m_work_function = pset.get<double>("work_function");
        m_recombination_factor =pset.get<double>("recombination_factor");
        //m_gain =pset.get<double>("gain");
        m_gain_mc =pset.get<std::vector<double>>("gain_mc");
        m_gain_data =pset.get<std::vector<double>>("gain_data");
        m_wire_spacing = pset.get<double>("wire_spacing");
        m_width_dqdx_box = pset.get<double>("width_box");
        m_length_dqdx_box = pset.get<double>("length_box");
        m_truthmatching_signaldef = pset.get<std::string>("truthmatching_signaldef");
        m_pidLabel = pset.get<std::string>("ParticleIDLabel","particleid");
        m_shower3dLabel = pset.get<std::string>("Shower3DLabel","shrreco3d");

        m_run_all_pfps = pset.get<bool>("runAllPFPs",false); //See .fcl file!
        m_exiting_photon_energy_threshold = pset.get<double>("exiting_photon_energy");
        m_exiting_proton_energy_threshold = pset.get<double>("exiting_proton_energy");

		//------------ BobbyVertexBuilder ---------------
		m_bobbyvertexing_more = pset.get<bool>("BobbyVertex_more_object", false);
		//------------ BobbyVertexBuilder ---------------

        rangen = new TRandom3(22);
        bool_make_sss_plots = true;

        std::vector<std::string> delta_names = {"Delta++","Delta+","Delta-","Delta0"};
        std::vector<int> delta_pdg_list = {2224,2214,1114,2114};
        for(size_t i=0; i< delta_pdg_list.size(); ++i){
            is_delta_map[delta_pdg_list[i]] = delta_names[i];
            is_delta_map[-delta_pdg_list[i]] ="Anti-"+delta_names[i];
        }

        //output text file of events can be set to true
        //initialize outfile
        //io stream to write to .txt file for EVD      
        //  std::ofstream out_stream;

        if (m_print_out_event && false){

            out_stream.open("v12_ncdelta_missing_trackshower_events.txt");
            if (!out_stream.is_open()){
                std::cout<<"ERROR output file not open"<<std::endl;
                exit(0);
            }
        }
 
                std::vector<std::string> inputVars = { "sss_candidate_num_hits", "sss_candidate_num_wires", "sss_candidate_num_ticks", "sss_candidate_PCA", "log10(sss_candidate_impact_parameter)", "log10(sss_candidate_min_dist)", "sss_candidate_impact_parameter/sss_candidate_min_dist", "sss_candidate_energy*0.001", "cos(sss_candidate_angle_to_shower)", "sss_candidate_fit_slope", "sss_candidate_fit_constant", "sss_candidate_plane", "sss_reco_shower_energy*0.001", "2*0.001*0.001*sss_reco_shower_energy*sss_candidate_energy*(1-cos(sss_candidate_angle_to_shower))", "log10(2*0.001*0.001*sss_reco_shower_energy*sss_candidate_energy*(1-cos(sss_candidate_angle_to_shower)))", "sss_candidate_energy*0.001/(sss_reco_shower_energy*0.001)", "sss_candidate_closest_neighbour" };
                sssVetov1 = new ReadBDT(inputVars);


    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void SinglePhoton::analyze(const art::Event &evt)
	{//analyzing one event per run! 

		m_is_verbose = true;
		std::cout<<"---------------------------------------------------------------------------------"<<std::endl;
		std::cout<<"SinglePhoton::analyze()\t||\t On entry: "<<m_number_of_events<<std::endl;

		auto const TPC = (*geom).begin_TPC();
		auto ID = TPC.ID();
		m_Cryostat = ID.Cryostat;
		m_TPC = ID.TPC;

		_time2cm = theDetector->SamplingRate() / 1000.0 * theDetector->DriftVelocity( theDetector->Efield(), theDetector->Temperature() );//found in ProtoShowerPandora_tool.cc

		this->ClearVertex();

		//******************************Setup*******************************************************/
		//******************************************************************************************/
		// OK in this section we will get a bunch of stuff we need later, general associations and maps. These are either from pandora helper stuff or from direct associations. 
		// Make sure under the hood you understand this!
		// ------------------------
		// The basic idea is that here we get every possible vector of data products we might need, along with all maps. e.g 
		// tracks->pfparticles->hits
		// tracks->pfparticles->spacepoints ..etc..
		//
		// And then later we just write pseudo-independant code assuming you have every objecets you want (see analyze_Tracks.h) and assume you have access to everything. 


		//BadChannels
		art::Handle<std::vector<int> > badChannelHandle;
		evt.getByLabel(m_badChannelProducer, m_badChannelLabel, badChannelHandle);
		std::vector<int> badChannelVector = *(badChannelHandle);



		//And some verticies.        
		art::ValidHandle<std::vector<recob::Vertex>> const & vertexHandle = evt.getValidHandle<std::vector<recob::Vertex>>(m_pandoraLabel);
		std::vector<art::Ptr<recob::Vertex>> vertexVector;
		art::fill_ptr_vector(vertexVector,vertexHandle);

	
		//HERE COMES THE Atlas!!!
		Atlas object_container(evt, {m_trackLabel, m_showerLabel, m_hitfinderLabel, m_flashLabel, m_pandoraLabel,m_shower3dLabel,m_showerKalmanLabel,m_showerKalmanCaloLabel, m_generatorLabel, m_geantModuleLabel}, m_is_data);
		
		//still need these handls for now;
		art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);//This is useful for FindManyP< reco::Track/Shower>
		art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitfinderLabel); //yea, this should be gone;
		art::ValidHandle<std::vector<recob::Track>> const & trackHandle  = evt.getValidHandle<std::vector<recob::Track>>(m_trackLabel);
		std::vector<art::Ptr<recob::OpFlash>> flashVector = object_container.all_opflashes;

		//add new FindManyP
		art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
		art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt,  m_pandoraLabel);
	
		//Get back to original structure;
		//typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;
		std::vector<art::Ptr<recob::PFParticle>> pfParticleVector = object_container.all_pfparticles;
		std::vector<art::Ptr<recob::Hit>> hitVector = object_container.all_hits;
		std::vector< art::Ptr<recob::Cluster> > clusterVector = object_container.all_clusters;
		std::vector<art::Ptr<recob::Slice>> sliceVector = object_container.sliceVector;


		PFParticleIdMap pfParticleMap = object_container.IDToPFParticleMap;
		std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticlesToVerticesMap = object_container.PFParticlesToVerticesMap;
		std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > pfParticleToMetadataMap = object_container.PFParticleToMetadataMap;
//		std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>> pfParticlesToShowerReco3DMap = object_container.PFParticlesToShowerReco3DMap;
        std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Track>> pfParticlesToShowerKalmanMap = object_container.PFParticlesToShowerKalmanMap;
        std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>>> kalmanTrackToCaloMap = object_container.kalmanTrackToCaloMap;
		std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap = object_container.PFParticleToSpacePointsMap;
		std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::PFParticle>> > sliceToPFParticlesMap = object_container.sliceToPFParticlesMap;
		std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::Hit>> > sliceToHitsMap = object_container.sliceToHitsMap;
		std::map<int, std::vector<art::Ptr<recob::PFParticle>> > sliceIDToPFParticlesMap = object_container.sliceIDToPFParticlesMap;
		std::map<int, std::vector<art::Ptr<recob::Hit>> > sliceIDToHitsMap = object_container.sliceIDToHitsMap;
		std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap = object_container.PFParticleToClustersMap;
//		std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap = object_container.ClusterToHitsMap;


		//Once we have actual verticies, lets concentrate on JUST the neutrino PFParticles for now:
		//--------------------------------
		// Produce two PFParticle vectors containing final-state particles:
		// 1. Particles identified as cosmic-rays - recontructed under cosmic-hypothesis
		// 2. Daughters of the neutrino PFParticle - reconstructed under the neutrino hypothesis
		std::vector< art::Ptr<recob::PFParticle> > crParticles;
		std::vector< art::Ptr<recob::PFParticle> > nuParticles;
		this->GetFinalStatePFParticleVectors(pfParticleMap, pfParticlesToVerticesMap, crParticles, nuParticles);

		//if not running over neutrino slice only, use all pfp's in event
		if (m_run_all_pfps){
			nuParticles = pfParticleVector;
		}


		//taking out the Larpandora helper functions here because they don't match to non-neutrino slice hits for some reason

		//OK Here we build two IMPORTANT maps for the analysis, (a) given a PFParticle get a vector of hits..
		//and (b) given a single hit, get the PFParticle it is in (MARK: is it only one? always? RE-MARK: Yes)
		std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
		//        std::map<art::Ptr<recob::Hit>, art::Ptr<recob::PFParticle>>                hitToPFParticleMap;
		//Using a pandora helper here, but to be honest we should probably just build using normal associations so keep independant if pssoble
		// lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(evt, m_pandoraLabel, pfParticleToHitsMap, hitToPFParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters);

		//use pfp->cluster and cluster->hit to build pfp->hit map
		//for each PFP
		for(size_t i=0; i<nuParticles.size(); ++i){
			auto pfp = nuParticles[i];

			std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp] ;

			//make empty vector to store hits
			std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

			//for each cluster, get the associated hits
			for (art::Ptr<recob::Cluster> cluster: clusters_vec){
				std::vector<art::Ptr<recob::Hit>> hits_vec =  object_container.ClusterToHitsMap[cluster];

				//insert hits into vector
				hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
			}

			//fill the map
			pfParticleToHitsMap[pfp] = hits_for_pfp;

		}//for each pfp

        //      Test ground for some slice stuff
        std::cout<<"SliceTest: there are "<<sliceVector.size()<<" slices in this event"<<std::endl;
        for(size_t s =0; s<sliceVector.size(); s++){
            auto slice = sliceVector[s];
            auto pfps = sliceToPFParticlesMap[slice]; 

            std::cout<<"SliceTest: On Slice "<<s<<" it has "<<pfps.size()<<" pfparticles"<<std::endl;
            std::vector<float> nu_scores;
            bool isSelectedSlice = false;
            int primaries = 0;
            int primary_pdg = 0;

            for(auto &pfp: pfps){
                std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadatas = pfParticleToMetadataMap[pfp];
                for(auto &meta: metadatas){
                    std::map<std::string, float> propertiesmap  = meta->GetPropertiesMap();
                    //for each of the things in the list
                    if(propertiesmap.count("NuScore")==1){
                        nu_scores.push_back(propertiesmap["NuScore"]);
                    }
                    if(propertiesmap.count("IsNeutrino")==1){
                       isSelectedSlice = true; 
                    }
                }

                if (pfp->IsPrimary()) {
                    primaries++;
                    primary_pdg = (pfp->PdgCode());    
                }
                /*if (!pfp->IsPrimary()) continue;
                // Check if this particle is identified as the neutrino
                const int pdg(pfp->PdgCode());
                const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
                if(isNeutrino){
                    isSelectedSlice = true; 
                }*/
            }

            if(nu_scores.size()>0){
                double mean  = std::accumulate(nu_scores.begin(), nu_scores.end(), 0.0)/(double)nu_scores.size();
                if(mean!=nu_scores.front()){
                    std::cout<<"ERROR! Somehow the pfp's in this slice have different nu-scores? IMpossible."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                std::cout<<"SliceTest: -- and has a nu_score of "<<nu_scores.front()<<std::endl;
                std::cout<<"SliceTest: -- with "<<primaries<<" primaries: pdg last: "<<primary_pdg<<std::endl;
            }else{
                std::cout<<"SliceTest: -- and does not have a nu_score of. "<<std::endl;
            }
            if(isSelectedSlice) std::cout<<"SliceTest: -- -- And is the Selected Neutrino Slice"<<std::endl;

        }




		// These are the vectors to hold the tracks and showers for the final-states of the reconstructed neutrino
		//At this point, nuParticles is a std::vector< art::Ptr<recon::PFParticle>> of the PFParticles that we are interested in.
		//tracks is a vector of recob::Tracks and same for showers.
		//Implicitly, tracks.size() + showers.size() =  nuParticles.size(); At this point I would like two things.

		if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Tracks and Showers"<<std::endl;

		//these are all filled in analyze slice
//		std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind
		std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> primaryPFPSliceIdVec; //stores a pair of only the primary PFP's in the event and the slice ind
//		std::map<int, double> sliceIdToNuScoreMap; //map between a slice Id and neutrino score
//		std::map<art::Ptr<recob::PFParticle>, bool> object_container.PFPToClearCosmicMap; //returns true for clear cosmic, false otherwise
//		std::map<art::Ptr<recob::PFParticle>, int> PFPToSliceIdMap; //returns the slice id for all PFP's
//		std::map<art::Ptr<recob::PFParticle>,bool> PFPToNuSliceMap;
//		std::map<art::Ptr<recob::PFParticle>,double> PFPToTrackScoreMap;
//		std::map<int, int> sliceIdToNumPFPsMap;
		//above are all filled in analyze slice

		//Need the following to identify a PFParticle is a track or a shower;
		//Keng is going to pack these FindManyP up to a function using the pfParticleHandle;


//	cout<<"CHECK! create object\n\n\n\n"<<endl;
		//Make an object and pack up stuffs we need.
		object_container.particles				= nuParticles;//no cosmic ray PFParticles

		//The following 4 lines are going to be transferred to atlas;
		art::FindManyP< recob::Track     > pfPartToTrackAssoc(pfParticleHandle, evt, m_trackLabel);
		art::FindManyP< recob::Shower    > pfPartToShowerAssoc(pfParticleHandle, evt, m_showerLabel);
		object_container.PFParticleAsATrack		= &pfPartToTrackAssoc;
		object_container.PFParticleAsAShower		= &pfPartToShowerAssoc;
		//Oooohu, here we go!

//		int run_count = 0;
//redo_event:
		this->AnalyzeSlices(
				object_container.PFParticleToMetadataMap, 
				object_container.IDToPFParticleMap,  
				object_container.primaryPFPSliceIdVec, 
				object_container.sliceIdToNuScoreMap, 
				object_container.PFPToClearCosmicMap, 
				object_container.PFPToSliceIdMap, 
				object_container.PFPToNuSliceMap, 
				object_container.PFPToTrackScoreMap);
		
		this->CollectTracksAndShowers_v2(evt, object_container); //This tells what showers and tracks to use.
		std::vector< art::Ptr<recob::Track> > tracks;
		tracks.reserve( object_container.selected_tracks.size() + object_container.more_tracks.size() );
		tracks.insert( tracks.end(), object_container.selected_tracks.begin(), object_container.selected_tracks.end() );
		tracks.insert( tracks.end(), object_container.more_tracks.begin(), object_container.more_tracks.end() );

		std::vector< art::Ptr<recob::Shower> > showers;
		showers.reserve( object_container.selected_showers.size() + object_container.more_showers.size() );
		showers.insert( showers.end(), object_container.selected_showers.begin(), object_container.selected_showers.end() );
		showers.insert( showers.end(), object_container.more_showers.begin(), object_container.more_showers.end() );



		//		std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle >> trackToNuPFParticleMap = object_container.trackToNuPFParticleMap; //give access to the PFParticle via track/shower
		//		std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap = object_container.showerToNuPFParticleMap;

		//		std::map<int, art::Ptr<simb::MCParticle>> MCParticleToTrackIdMap = object_container.MCParticleToTrackIdMap;

		//m_vertex_pos_x / y/ z are ready to be used now;


		 geoalgo::Point_t pvertex(m_vertex_pos_x, m_vertex_pos_y, m_vertex_pos_z);
		//use object_container.trackToDistMap/showerToDistMap;
		m_is_verbose = false;
		this->AnalyzeTracks(
				object_container,
				pvertex,
				tracks,
				object_container.trackToNuPFParticleMap,
				pfParticleToSpacePointsMap,
				object_container.MCParticleToTrackIdMap,//disabled
				object_container.sliceIdToNuScoreMap,
				object_container.PFPToClearCosmicMap,
				object_container.PFPToSliceIdMap,
				object_container.PFPToTrackScoreMap,
				object_container.PFPToNuSliceMap,
				pfParticleMap);

		this->AnalyzeShowers(
				object_container,
				pvertex,
				showers,
				object_container.showerToNuPFParticleMap,
				pfParticleToHitsMap, 
				pfParticleToClustersMap, 
				object_container.ClusterToHitsMap,
				object_container.sliceIdToNuScoreMap, 
				object_container.PFPToClearCosmicMap,
				object_container.PFPToSliceIdMap, 
				object_container.PFPToNuSliceMap, 
				object_container.PFPToTrackScoreMap,
				pfParticleMap,
				object_container.PFParticlesToShowerReco3DMap); 
		m_is_verbose = true;

		if(!m_is_data){
		//Borrow the following from MCTruth
		art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcparticles_per_hit(hitHandle, evt, m_hitMCParticleAssnsLabel);


		//mcc9 march miniretreat fix
		std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
		std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing

		m_test_matched_hits = 0;

		for(size_t j=0; j<hitVector.size();j++){
			const art::Ptr<recob::Hit> hit = hitVector[j];

			particle_vec.clear(); match_vec.clear(); //only store per hit

			mcparticles_per_hit.get(hit.key(), particle_vec, match_vec);

			if(particle_vec.size() > 0){
				m_test_matched_hits++;
			}

		}
		this->CollectMCParticles_v2(evt, object_container);
//		this->CollectMCParticles(
//				evt, 
//				m_geantModuleLabel, 
//				object_container.MCTruthToMCParticlesMap, 
//				object_container.MCParticleToMCTruthMap, 
//				object_container.MCParticleToTrackIdMap);//created here;

		this->ResizeShowers(showers.size());
		this->ResizeTracks(tracks.size());

		this->showerRecoMCmatching(//Shower MCTruth here?
				showers,
				object_container.showerToMCParticleMap, //created here
				object_container.showerToNuPFParticleMap, 
				pfParticleToHitsMap, 
				mcparticles_per_hit, //see above
				object_container.matchedMCParticleVector, //created here
				pfParticleMap,  
				object_container.MCParticleToTrackIdMap,  //input
				object_container.sliceIdToNuScoreMap, 
				object_container.PFPToClearCosmicMap, 
				object_container.PFPToSliceIdMap, 
				object_container.PFPToNuSliceMap);

		std::vector<double> trk_overlay_vec = recoMCmatching<art::Ptr<recob::Track>>( 
				tracks, 
				object_container.trackToMCParticleMap, //created here
				object_container.trackToNuPFParticleMap, 
				pfParticleToHitsMap, 
				mcparticles_per_hit, 
				object_container.matchedMCParticleVector);

		this->RecoMCTracks(//Find Track MCTruth here?
				tracks,
				object_container,
				object_container.trackToNuPFParticleMap,
				object_container.trackToMCParticleMap,
				object_container.MCParticleToMCTruthMap,
				object_container.matchedMCParticleVector, 
				//mcParticleVector,
				object_container.MCParticleToTrackIdMap,
				object_container.sliceIdToNuScoreMap,
				object_container.PFPToClearCosmicMap,
				object_container.PFPToSliceIdMap,
				trk_overlay_vec);
		}
		//-------------------------------


		//---------- VertexBuilder--------------
		//use the new the new class for variables and vertexing.
		//		ParticleAssociations_all const & bobby_particle_associations = BobbyVertexBuilder_ext(object_container, m_bobbyvertexing_more );
		BobbyVertexBuilder(object_container, m_bobbyvertexing_more );
		//introduce a for loop for all particle associations identified by Bobby's VertexBuilder

/*CHECK

std::cout<<"Filling in Bobby's Vertex info. with "<<bobby_particle_associations.GetSelectedAssociations().size()<<" Vertex candidates."<<std::endl;
		bool reset_bobbyvertex = true;
		if(bobby_particle_associations.GetSelectedAssociations().size()==0){
		cout<<"No vertex is reconstructed."<<endl;
		}

		for(size_t const nth_associations : bobby_particle_associations.GetSelectedAssociations()) {//Loop over all associations, which is a vector
			ParticleAssociation const & particle_associated = bobby_particle_associations.GetAssociations().at(nth_associations);//grab the "pn"th association;
			geoalgo::Point_t const & reco_vertex = particle_associated.GetRecoVertex();//Grab the vertec of the "pn"th association.
			if(reset_bobbyvertex){
				m_bobbyvertex_pos_xv.clear();
				m_bobbyvertex_pos_yv.clear();
				m_bobbyvertex_pos_zv.clear();
				m_bobbytracksv.clear();
				m_bobbyshowersv.clear();

				m_bobbyprotontrack.clear();
				m_bobbyphotonshower.clear();
				m_bobbypi0daughter.clear();
				reset_bobbyvertex = false;
			}

			m_bobbyvertex_pos_xv.push_back( reco_vertex.at(0));
			m_bobbyvertex_pos_yv.push_back( reco_vertex.at(1));
			m_bobbyvertex_pos_zv.push_back( reco_vertex.at(2));

			cout<<"Vertex Coordinates found by Bobby Vertex Builder: ";
			cout<<reco_vertex.at(0)<<", ";
			cout<<reco_vertex.at(1)<<", ";
			cout<<reco_vertex.at(2)<<endl;

			//calculate the # of tracks/showers;
			DetectorObjects_all const & detos = bobby_particle_associations.GetDetectorObjects();
			int temp_num_tracks = 0;
			int temp_num_showers = 0;
			
			int get_a_proton = 0;
			int get_a_photon = 0;
			int get_a_pi0daughter = 0;
			for(size_t const n : particle_associated.GetObjectIndices()) {
				if(detos.GetRecoType(n) == detos.ftrack_reco_type) {
					++temp_num_tracks;
					
					int trackindex = detos.GetTrackIndexFromObjectIndex(n);
					art::Ptr<simb::MCParticle> temp_MCtrack = object_container.trackToMCParticleMap.find(tracks[trackindex])->second;
					if(temp_MCtrack->PdgCode()==2212){
						get_a_proton++;
					}
				}
				if(detos.GetRecoType(n) == detos.fshower_reco_type) {

					++temp_num_showers;
					int showerindex = detos.GetShowerIndexFromObjectIndex(n);//CHECK
					art::Ptr<simb::MCParticle> temp_MCshower = object_container.showerToMCParticleMap.find(showers[showerindex])->second;
					if(temp_MCshower->PdgCode()==22){
						get_a_photon++;
					}
					art::Ptr<simb::MCParticle> amother = object_container.MCParticleToTrackIdMap[temp_MCshower->Mother()];
					if(amother){//sometime Mother is unknown..
						if(amother->PdgCode() == 111){
							get_a_pi0daughter++;
						}
					}
				}
			}
			cout<<"# of showers: "<<temp_num_showers<<endl;
			cout<<"# of tracks : "<<temp_num_tracks<<endl;

			m_bobbytracksv.push_back(temp_num_tracks);
			m_bobbyshowersv.push_back(temp_num_showers);

			m_bobbyprotontrack.push_back(get_a_proton);
			m_bobbyphotonshower.push_back(get_a_photon);
			m_bobbypi0daughter.push_back(get_a_pi0daughter);
		}

//			double best_vertex_dist = SIZE_MAX;
			for( size_t index = 0; index< m_bobbyvertex_pos_xv.size() ; index++){
//				double temp_dist =  pow(m_vertex_pos_x - m_bobbyvertex_pos_xv[index],2) + pow(m_vertex_pos_y - m_bobbyvertex_pos_yv[index],2)+ pow(m_vertex_pos_z - m_bobbyvertex_pos_zv[index],2);

//				if(temp_dist < best_vertex_dist){}//update when find a closer vertex to the pandora vertex.
				if((m_bobbyshowersv[index] > 0 && m_bobbytracksv[index] > 0)|| index == 0){
				//	best_vertex_dist = temp_dist;
					
					m_bobbyvertex_pos_x = m_bobbyvertex_pos_xv[index];
					m_bobbyvertex_pos_y = m_bobbyvertex_pos_yv[index];
					m_bobbyvertex_pos_z = m_bobbyvertex_pos_zv[index];
					m_bobbytracks =  m_bobbytracksv[index];
					m_bobbyshowers = m_bobbyshowersv[index];
					break;
				}
			}

//			if( m_run_all_pfps && m_bobbyvertexing_more &&m_bobbytracks + m_bobbyshowers < 2&& run_count < 2 ){//repeat with 0,1
//				run_count++;
//				cout<<"Need more objects. Now run the "<<run_count<<" time the vertexing."<<endl;
//				goto redo_event;
//			}

			std::cout<<"Got Bobby's info.!\n"<<std::endl;
*/

			//----------------------------------------------------
			//Track Calorimetry
		art::FindManyP<anab::Calorimetry> calo_per_track(trackHandle, evt, m_caloLabel);
		std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<anab::Calorimetry>> > trackToCalorimetryMap;
		//So a cross check
		if (!calo_per_track.isValid())
		{
			mf::LogDebug("SinglePhoton") << "  Failed to get Assns between recob::Track and anab::Calorimetry.\n";
			return;
		}
		for(size_t i=0; i< tracks.size(); ++i){
			if(calo_per_track.at(tracks[i].key()).size() ==0){
				std::cerr<<"Track Calorimetry Breaking!  the vector of calo_per_track is of length 0 at this track."<<std::endl;
			}

			size_t calo_size = calo_per_track.at(tracks[i].key()).size();
			//std::cout<<"Track Calo from producer: "<<m_caloLabel<<" has "<<calo_size<<" anab::Calorimetry objects associaed."<<std::endl;
			trackToCalorimetryMap[tracks[i]] = calo_per_track.at(tracks[i].key());
			for(size_t k=0; k<calo_size; k++){
				//std::cout<<"Calo "<<k<<" PlaneID: "<<calo_per_track.at(tracks[i].key())[k]->PlaneID()<<std::endl;
			}
		}

		art::FindOneP<anab::ParticleID> pid_per_track(trackHandle, evt, m_pidLabel);
		std::map<art::Ptr<recob::Track>, art::Ptr<anab::ParticleID> > trackToPIDMap;

		if(m_use_PID_algorithms){
			// Build a map to get PID from PFParticles, then call PID collection function
			for(size_t i=0; i< tracks.size(); ++i){
				art::Ptr<recob::Track> track = tracks[i];
				trackToPIDMap[track] = pid_per_track.at(track.key());
			}
		}
	
			//-------------------------------------

			//CRT 
			/*
			   if(m_has_CRT){
			   art::ValidHandle<std::vector<crt::CRTTzero>> const & crtHandle  = evt.getValidHandle<std::vector<crt::CRTTzero>>(m_CRTTzeroLabel);
			   std::vector<art::Ptr<crt::CRTTzero>> crtVector;
			   art::fill_ptr_vector(crtVector,crtHandle);
			   }
			   */


			//**********************************************************************************************/
			//**********************************************************************************************/
			//---------------------------------- MC TRUTH Data Only---------------------------
			//**********************************************************************************************/
			//**********************************************************************************************/

			//Get the MCtruth handles and vectors
			std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;
//			std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

			//Then build a map from MCparticles to Hits and vice versa
			std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  mcParticleToHitsMap;
			std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> > hitToMCParticleMap;

			//Apparrently a MCParticle doesn't know its origin (thanks Andy!)
			//I would also like a map from MCparticle to MCtruth and then I will be done.  and Vice Versa
			//Note which map is which!       //First  is one-to-many.         //Second is one-to-one
//			std::map< art::Ptr<simb::MCTruth>,    std::vector<art::Ptr<simb::MCParticle>>>  MCTruthToMCParticlesMap;
//			std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>                  MCParticleToMCTruthMap;


//			std::vector<art::Ptr<simb::MCParticle>> matchedMCParticleVector;
//			std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > trackToMCParticleMap;
//			std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > showerToMCParticleMap;

			//Given a simb::MCParticle we would like a map to either a sim::MCTrack or sim::MCShower
			std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCTrack> > MCParticleToMCTrackMap;
			std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCShower> > MCParticleToMCShowerMap;


			//**********************************************************************************************/
			//**********************************************************************************************/
			//Some event based properties

			m_number_of_events++;

			m_run_number = evt.run();
			m_subrun_number = evt.subRun();
			m_event_number = evt.id().event();

			if(vertexVector.size()>0){
				m_number_of_vertices++;
			}


			//and now get the simb::MCparticle to both MCtrack and MCshower maps (just for the MCparticles matched ok).

			badChannelMatching<art::Ptr<recob::Track>>(badChannelVector, tracks, object_container.trackToNuPFParticleMap, pfParticleToHitsMap,geom,bad_channel_list_fixed_mcc9);

			if(m_is_verbose){
				std::cout << "SinglePhoton::analyze()\t||\t Consolidated event summary:" << "\n";
				std::cout << "SinglePhoton::analyze()\t||\t - Number of primary cosmic-ray PFParticles   : " << crParticles.size() << "\n";
				std::cout << "SinglePhoton::analyze()\t||\t - Number of neutrino final-state PFParticles : " << nuParticles.size() << "\n";
				std::cout << "SinglePhoton::analyze()\t||\t    ... of which are track-like   : " << tracks.size() << "\n";
				std::cout << "SinglePhoton::analyze()\t||\t    ... of which are showers-like : " << showers.size() << "\n";
			}

			std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t Starting"<<std::endl;
//			this->AnalyzeSlices(pfParticleToMetadataMap, pfParticleMap,  primaryPFPSliceIdVec, object_container.sliceIdToNuScoreMap, object_container.PFPToClearCosmicMap, object_container.PFPToSliceIdMap, object_container.PFPToNuSliceMap, object_container.PFPToTrackScoreMap);
			//std::cout<<"There are "<< allPFPSliceIdVec.size()<<" pfp-slice id matches stored in the vector"<<std::endl;
			std::cout<<"the number of PPF's with stored clear cosmic info is "<<object_container.PFPToClearCosmicMap.size()<<std::endl;
			std::cout<<"the number of PFP's stored in the object_container.PFPToSliceIdMap is "<<object_container.PFPToSliceIdMap.size()<<std::endl;
			if (object_container.PFPToSliceIdMap.size() < 1){
				std::cout<<"ERROR, not storing PFP's in object_container.PFPToSliceIdMap"<<std::endl;
			}

			for (auto pair:object_container.PFPToNuSliceMap){
				auto pfp = pair.first;
				auto is_nuslice = pair.second;
				if (is_nuslice){
					std::cout<<"pfp in nuslice "<<pfp->Self()<<std::endl;
				}

			}

			for (auto pair:sliceIDToPFParticlesMap){ 
				std::vector<art::Ptr<recob::PFParticle>> pfp_vec = pair.second;
				int slice_id = pair.first;
				//if (slice_vec[0]->Slice() != object_container.PFPToSliceIdMap[pfp] )
				for(auto pfp: pfp_vec){
					if (slice_id != object_container.PFPToSliceIdMap[pfp] && object_container.PFPToSliceIdMap[pfp]>=0){
						std::cout<<"sliceIDToPFParticlesMap[slice->ID()] for pfp "<<pfp->Self()<<" is slice "<< slice_id<< "but object_container.PFPToSliceIdMap[pfp] = "<< object_container.PFPToSliceIdMap[pfp]<<std::endl;
					}
				}

			}


			
			//if CRT info, get CRT hits
			art::Handle<std::vector<crt::CRTHit>> crthit_h; //only filled when there are hits, otherwise empty
			art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
			double evt_timeGPS_nsec = -999 ;
			if(m_runCRT){
				evt.getByLabel(m_DAQHeaderProducer, rawHandle_DAQHeader);

				evt.getByLabel(m_CRTHitProducer, crthit_h);
				raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
				art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
				evt_timeGPS_nsec = evtTimeGPS.timeLow();

					std::cout<<"SinglePhoton::analyze \t||\t Got CRT hits"<<std::endl;
			}

        this->AnalyzeFlashes(flashVector, crthit_h, evt_timeGPS_nsec);
        //   this->AnalyzeFlashes(flashVector, crthit_h);

//			std::cout<<"start track"<<std::endl;
//			this->AnalyzeTracks(tracks, object_container.trackToNuPFParticleMap, pfParticleToSpacePointsMap,  object_container.MCParticleToTrackIdMap, object_container.sliceIdToNuScoreMap, object_container.PFPToClearCosmicMap,  object_container.PFPToSliceIdMap,  object_container.PFPToTrackScoreMap, object_container.PFPToNuSliceMap,pfParticleMap);
//			this->AnalyzeShowers(showers,object_container.showerToNuPFParticleMap, pfParticleToHitsMap, pfParticleToClustersMap, clusterToHitsMap,sliceIdToNuScoreMap, object_container.PFPToClearCosmicMap,  object_container.PFPToSliceIdMap, object_container.PFPToNuSliceMap, object_container.PFPToTrackScoreMap,pfParticleMap,pfParticlesToShowerReco3DMap); 
			this->AnalyzeKalmanShowers(showers,object_container.showerToNuPFParticleMap,pfParticlesToShowerKalmanMap, kalmanTrackToCaloMap, pfParticleToHitsMap);


			this->AnalyzeTrackCalo(tracks,   trackToCalorimetryMap);

			if(m_use_PID_algorithms)  this->CollectPID(tracks, trackToPIDMap);
			// MCTruth, MCParticle, MCNeutrino information all comes directly from GENIE.
			// MCShower and MCTrack come from energy depositions in GEANT4
			if(!m_is_data){

				art::ValidHandle<std::vector<simb::MCTruth>> const & mcTruthHandle= evt.getValidHandle<std::vector<simb::MCTruth>>(m_generatorLabel);
				art::fill_ptr_vector(mcTruthVector,mcTruthHandle);

				art::ValidHandle<std::vector<simb::MCParticle>> const & mcParticleHandle= evt.getValidHandle<std::vector<simb::MCParticle>>(m_geantModuleLabel);
//				art::fill_ptr_vector(mcParticleVector,mcParticleHandle);
				art::fill_ptr_vector(object_container.matchedMCParticleVector,mcParticleHandle);



				//Get the MCParticles (move to do this ourselves later)
//				this->CollectMCParticles(evt, m_geantModuleLabel, MCTruthToMCParticlesMap, MCParticleToMCTruthMap, object_container.MCParticleToTrackIdMap);

				//OK lets get all set up with sim::MCTrack and sim::MCShower .

				//   art::ValidHandle<std::vector<sim::MCTrack>> const & mcTrackHandle  = evt.getValidHandle<std::vector<sim::MCTrack>>(m_mcTrackLabel);
				// art::ValidHandle<std::vector<sim::MCShower>> const & mcShowerHandle  = evt.getValidHandle<std::vector<sim::MCShower>>(m_mcShowerLabel);
				//  art::fill_ptr_vector(mcTrackVector,mcTrackHandle);
				//  art::fill_ptr_vector(mcShowerVector,mcShowerHandle);

/*CHECK: move this to abvoe
				art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcparticles_per_hit(hitHandle, evt, m_hitMCParticleAssnsLabel);


				//mcc9 march miniretreat fix
				std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
				std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing

				m_test_matched_hits = 0;

				for(size_t j=0; j<hitVector.size();j++){
					const art::Ptr<recob::Hit> hit = hitVector[j];

					particle_vec.clear(); match_vec.clear(); //only store per hit

					mcparticles_per_hit.get(hit.key(), particle_vec, match_vec);

					if(particle_vec.size() > 0){
						m_test_matched_hits++;
					}

				}
*/

				this->BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector,  mcParticleToHitsMap, hitToMCParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters,  object_container.MCParticleToTrackIdMap);

				std::cout<<"SinglePhoton\t||\t Starting backtracker on recob::track"<<std::endl;
//				std::vector<double> trk_overlay_vec = recoMCmatching<art::Ptr<recob::Track>>( tracks, object_container.trackToMCParticleMap, object_container.trackToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, object_container.matchedMCParticleVector);


				std::cout<<"SinglePhoton\t||\t Starting backtracker on recob::shower"<<std::endl;
//				this->showerRecoMCmatching(showers, showerToMCParticleMap, object_container.showerToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, matchedMCParticleVector, pfParticleMap,  object_container.MCParticleToTrackIdMap, object_container.sliceIdToNuScoreMap, object_container.PFPToClearCosmicMap,  object_container.PFPToSliceIdMap, object_container.PFPToNuSliceMap);


				//showerRecoMCmatching( showers, showerToMCParticleMap, object_container.showerToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, matchedMCParticleVector, pfParticleMap,  object_container.MCParticleToTrackIdMap);

				//looking at metadata
				//std::map<art::Ptr<recob::PFParticle>, double >  pfParticleToNuScoreMap;//is filled during analyze slices
				/*std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind
				  std::cout<<"SinglePhoton\t||\t Analyze Metadata"<<std::endl;
				  this->AnalyzeSlices( pfParticleToMetadataMap, pfParticleMap, allPFPSliceIdVec);
				  std::cout<<"There are "<< allPFPSliceIdVec.size()<<" pfp-slice id matches stored in the vector"<<std::endl;
				  if (showers.size()>0){
				  std::cout<<"the shower at 0 is in slice "<<this->GetShowerSlice(showers[0], object_container.showerToNuPFParticleMap, allPFPSliceIdVec)<<std::endl;
				  }

				  this->FindSignalSlice( m_truthmatching_signaldef, object_container.MCParticleToTrackIdMap);

				  for(auto & track: tracks){
				  std::cout<<"CHECKTRACK 0: "<<trackToMCParticleMap.count(track)<<std::endl;
				  }*/

				//  perfectRecoMatching<art::Ptr<sim::MCTrack>>(matchedMCParticleVector, mcTrackVector, MCParticleToMCTrackMap);
				// perfectRecoMatching<art::Ptr<sim::MCShower>>(matchedMCParticleVector, mcShowerVector, MCParticleToMCShowerMap);
				//OK a really wierd bug in which by accessing the map here in line 355, everything breaks.. but commenting it out is OK


				//for(auto & shower: showers){
				//    auto mp = showerToMCParticleMap[shower];
				//    std::cout<<"CHECKSHOWER: count trackmap: "<<MCParticleToMCTrackMap.count(mp)<<" "<< MCParticleToMCShowerMap.count(mp)<<std::endl;
				//}




//				this->RecoMCTracks(tracks, object_container.trackToNuPFParticleMap, trackToMCParticleMap, object_container.MCParticleToMCTruthMap,mcParticleVector, object_container.MCParticleToTrackIdMap, object_container.sliceIdToNuScoreMap, object_container.PFPToClearCosmicMap,  object_container.PFPToSliceIdMap,trk_overlay_vec);



				//Obsolete function
				//this->RecoMCShowers(showers, object_container.showerToNuPFParticleMap, showerToMCParticleMap, MCParticleToMCTruthMap,mcParticleVector);
//				this->AnalyzeMCTruths(mcTruthVector, mcParticleVector);
				this->AnalyzeMCTruths(mcTruthVector, object_container.matchedMCParticleVector);
				this->AnalyzeEventWeight(evt);

				//added since last time?
				std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind
            //this one was for testing, leaving out for now
            // this->FindSignalSlice( m_truthmatching_signaldef, object_container.MCParticleToTrackIdMap, object_container.showerToNuPFParticleMap , allPFPSliceIdVec, showerToMCParticleMap, object_container.trackToNuPFParticleMap, trackToMCParticleMap);
				if(m_is_verbose)std::cout<<"Starting SecondShowerSearch"<<std::endl;
//				this->SecondShowerSearch(tracks,  object_container.trackToNuPFParticleMap, showers, object_container.showerToNuPFParticleMap, pfParticleToHitsMap, object_container.PFPToSliceIdMap, sliceIDToHitsMap,mcparticles_per_hit, object_container.matchedMCParticleVector, pfParticleMap,  object_container.MCParticleToTrackIdMap);

            std::cout<<"filling info in ncdelta slice tree"<<std::endl;
            this->AnalyzeRecoMCSlices( m_truthmatching_signaldef, object_container.MCParticleToTrackIdMap, object_container.showerToNuPFParticleMap , allPFPSliceIdVec, object_container.showerToMCParticleMap, object_container.trackToNuPFParticleMap, object_container.trackToMCParticleMap,  object_container.PFPToSliceIdMap);

            if (m_print_out_event){
                if (m_matched_signal_shower_num != 1 || m_matched_signal_track_num != 1){
                    out_stream <<"run subrunevent "<<m_run_number<<" "<<m_subrun_number<<" "<<m_event_number<<"\n";
                }

            }
            std::cout<<"Going to grab eventweightSplines for CCQE genie fix, won't be necessary long term"<<std::endl;
            art::Handle<std::vector<evwgh::MCEventWeight>>  ev_evw ;
            if(    evt.getByLabel("eventweightSplines",ev_evw)){

                std::map<std::string, std::vector<double>> const & weight_map = ev_evw->front().fWeight;
                if(ev_evw->size() > 1) std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"<< "WARNING: eventweight slice genie fix has more than one entry\n";
                //m_genie_spline_weight=weight_map;
                for (auto const& x : weight_map){
                    std::cout << x.first  // string (key)
                        << ':' 
                        << x.second.size() << std::endl ;
                    if(x.second.size()==1){
                        m_genie_spline_weight = x.second.front();
                    }
                }

            }else{
                std::cout<<"No data producet called eventweightSplines"<<std::endl;
                m_genie_spline_weight =1.0;
            }


            std::cout<<"SinglePhoton::analyze\t||\t finnished loop for this event"<<std::endl;
        }else{

//            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> * tmp_mcparticles_per_hit = NULL;
            std::vector<art::Ptr<simb::MCParticle>> tmp_matchedMCParticleVector;
            
//            this->SecondShowerSearch(tracks,  object_container.trackToNuPFParticleMap, showers, object_container.showerToNuPFParticleMap, pfParticleToHitsMap, object_container.PFPToSliceIdMap, sliceIDToHitsMap,*tmp_mcparticles_per_hit, tmp_matchedMCParticleVector, pfParticleMap,  object_container.MCParticleToTrackIdMap);


        }

        //Second Shower Search-Pandora style
        if(!m_run_all_pfps){
			//Isolation
            this-> IsolationStudy(tracks,  object_container.trackToNuPFParticleMap, showers, object_container.showerToNuPFParticleMap, pfParticleToHitsMap, object_container.PFPToSliceIdMap, sliceIDToHitsMap);
        }


        //This is a quick check 

		size_t n_neutrino_slice=0;
        size_t n_neutrino_candidate_pfp_id=0;

        for(size_t s=0; s< sliceVector.size(); s++){
            auto slice = sliceVector[s];
            std::vector<art::Ptr<recob::PFParticle>> pfps = sliceToPFParticlesMap[slice]; 

            int primaries=0;
            int n_dau=0;
            int found = 0;
            //std::cout<<"Starting a loop over "<<pfps.size()<<" pfparticles"<<std::endl;
            for(auto &pfp: pfps){
                //std::cout<<pfp->Self()<<" Primary: "<<pfp->IsPrimary()<<" PDG "<<pfp->PdgCode()<<" NDau: "<<pfp->NumDaughters()<<" Parent: "<<pfp->Parent()<<std::endl;

                if (!pfp->IsPrimary()) continue;
                // Check if this particle is identified as the neutrino
                const int pdg(pfp->PdgCode());
                const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
                primaries++;
                // If it is, lets get the vertex position
                if(isNeutrino){
                    found++;
                    //Ok this is neutrino candidate. 
                    
                    std::cout<<"Found Neutrinoi Slice "<<s<<std::endl;
                    for(auto &pfp: pfps){
                        std::cout<<pfp->Self()<<" Primary: "<<pfp->IsPrimary()<<" PDG "<<pfp->PdgCode()<<" NDau: "<<pfp->NumDaughters()<<" Parent: "<<pfp->Parent()<<std::endl;
                    }
                    std::cout<<"************   Printing hierarcy "<<m_run_number<<" "<<m_subrun_number<<" "<<m_event_number<<" **************"<<std::endl;
                    n_neutrino_candidate_pfp_id = pfp->Self();
                    for (const size_t daughterId : pfp->Daughters()){
                        n_dau++;
                        auto dau = pfParticleMap[daughterId];
                        std::cout<<"---> gen1 --->"<<daughterId<<" trkScore: "<<object_container.PFPToTrackScoreMap[dau]<<" PDG: "<<dau->PdgCode()<<" NumDau: "<<dau->NumDaughters()<<std::endl;
                        auto tmp = dau;
                        int n_gen = 2;
                            for (const size_t granDaughterId : tmp->Daughters()){
                                while(tmp->NumDaughters()>0 && n_gen < 4){
                                    for(int k=0; k< n_gen; k++){
                                        std::cout<<"---> ";
                                    }
                                    auto grandau = pfParticleMap[granDaughterId];
 	                                 std::cout<<"gen"<<n_gen<<"  --->"<<granDaughterId<<" trkScore: "<< object_container.PFPToTrackScoreMap[grandau]<<" PDG: "<<grandau->PdgCode()<<" NumDau: "<<grandau->NumDaughters()<<std::endl;
                                    tmp = grandau;    
                                    n_gen++;
                                }
                            if(n_gen >=4) break;
                            }

                    }
                    std::cout<<"************   Finished hierarcy **************"<<std::endl;

                }
            }

            if(found==1){
                n_neutrino_slice = s;
                std::cout<<"Found a neutrino slice @ slice "<<n_neutrino_slice<<" ID "<<slice->ID()<<" key "<<slice.key()<<" pdfID "<<n_neutrino_candidate_pfp_id<<std::endl;
                std::cout<<"And there is "<<pfps.size()<<" PFParticles of which "<<primaries<<" are primary and "<<n_dau<<" are daughters of the Neutrino."<<std::endl;
                if((int)pfps.size() > n_dau+1){
                    std::cout<<"We're Missing Something!."<<std::endl;
                }
                m_reco_slice_objects = (int)pfps.size();
            }else if(found >1){
                throw cet::exception("DetachedVertexFinder") << "  This event contains multiple reconstructed neutrinos! Size: "<<found<<std::endl;
            }else if(found ==0){

            }
        }


		//---------------------- END OF LOOP, fill vertex ---------------------

		vertex_tree->Fill();
		ncdelta_slice_tree->Fill();

		std::cout<<"---------------------------------------------------------------------------------"<<std::endl;

	}




    //-------------------------------------------------------------------------------------------

    void SinglePhoton::endJob()
    {
        if (m_print_out_event){
            out_stream.close();
        }
        pot_tree->Fill();
    }

    //-------------------------------------------------------------------------------------------


    void SinglePhoton::beginJob()
    {
        mf::LogDebug("SinglePhoton") << " *** beginJob() *** " << "\n";

        art::ServiceHandle<art::TFileService> tfs;

        vertex_tree = tfs->make<TTree>("vertex_tree", "vertex_tree");
        pot_tree = tfs->make<TTree>("pot_tree", "pot_tree");
        eventweight_tree = tfs->make<TTree>("eventweight_tree", "eventweight_tree");
        ncdelta_slice_tree = tfs->make<TTree>("ncdelta_slice_tree", "ncdelta_slice_tree");


        // --------------------- POT Releated variables -----------------
        m_number_of_events = 0;
        m_number_of_vertices = 0;
        m_pot_count=0;
        pot_tree->Branch("number_of_events",&m_number_of_events,"number_of_events/I");
        pot_tree->Branch("number_of_vertices",&m_number_of_vertices,"number_of_vertices/I");
        pot_tree->Branch("POT",&m_pot_count,"POT/D");

        // --------------------- Event Related variables ------------
        vertex_tree->Branch("run_number", &m_run_number, "run_number/I");
        vertex_tree->Branch("subrun_number", &m_subrun_number, "subrun_number/I");
        vertex_tree->Branch("event_number", &m_event_number, "event_number/I");

        vertex_tree->Branch("genie_spline_weight", &m_genie_spline_weight, "genie_spline_weight/D");


        vertex_tree->Branch("test_matched_hits", &m_test_matched_hits, "test_matched_hits/I");
        // --------------------- Vertex Related variables ------------
        vertex_tree->Branch("reco_vertex_size", &m_reco_vertex_size);
        vertex_tree->Branch("reco_vertex_x", &m_vertex_pos_x);
        vertex_tree->Branch("reco_vertex_y", &m_vertex_pos_y);
        vertex_tree->Branch("reco_vertex_z", &m_vertex_pos_z);
        vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane0",&m_reco_vertex_to_nearest_dead_wire_plane0);
        vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane1",&m_reco_vertex_to_nearest_dead_wire_plane1);
        vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane2",&m_reco_vertex_to_nearest_dead_wire_plane2);
		
		//---------------------- BobbyVertexBuilder -----------------
        vertex_tree->Branch("reco_bobbyvertex_x", &m_bobbyvertex_pos_x);
        vertex_tree->Branch("reco_bobbyvertex_y", &m_bobbyvertex_pos_y);
        vertex_tree->Branch("reco_bobbyvertex_z", &m_bobbyvertex_pos_z);
        vertex_tree->Branch("reco_bobbyshowers", &m_bobbyshowers);
        vertex_tree->Branch("reco_bobbytracks", &m_bobbytracks);
		//vector
		
        vertex_tree->Branch("reco_bobbyvertex_xv", &m_bobbyvertex_pos_xv);
        vertex_tree->Branch("reco_bobbyvertex_yv", &m_bobbyvertex_pos_yv);
        vertex_tree->Branch("reco_bobbyvertex_zv", &m_bobbyvertex_pos_zv);
        vertex_tree->Branch("reco_bobbytracksv", &m_bobbytracksv);
		vertex_tree->Branch("reco_bobbyshowersv", &m_bobbyshowersv);
		vertex_tree->Branch("reco_bobbyvertexradiusv", &m_bobbyvertexradiusv);
		vertex_tree->Branch("reco_bobbyvertexradius", &m_bobbyvertexradius);
		//MCTruth Matrching;
		vertex_tree->Branch("mctruth_bobbysameslicev", &m_bobbysameslicev);
		vertex_tree->Branch("mctruth_bobbyprotontrackv", &m_bobbyprotontrackv);
		vertex_tree->Branch("mctruth_bobbyphotonshowerv", &m_bobbyphotonshowerv);
		vertex_tree->Branch("mctruth_bobbypi0daughterv", &m_bobbypi0daughterv);
		vertex_tree->Branch("mctruth_bobbydeltaradppdaughterv", &m_bobbydeltaradppdaughterv);
		vertex_tree->Branch("mctruth_bobbydeltaradmdaughterv", &m_bobbydeltaradmdaughterv);
		vertex_tree->Branch("mctruth_bobbydeltaradpdaughterv", &m_bobbydeltaradpdaughterv);
		vertex_tree->Branch("mctruth_bobbydeltarad0daughterv", &m_bobbydeltarad0daughterv);
		vertex_tree->Branch("mctruth_bobbyotherdaughterv", &m_bobbyotherdaughterv);
		vertex_tree->Branch("mctruth_bobbyoverlayv", &m_bobbyoverlayv);
//		vertex_tree->Branch("mctruth_bobbytrackdaughter_pdg", & m_bobbytrackdaughter_pdg);
//		vertex_tree->Branch("mctruth_bobbyshowerdaughter_pdg", & m_bobbyshowerdaughter_pdg);

		vertex_tree->Branch("mctruth_bobbyprotontrack", &m_bobbyprotontrack);
		vertex_tree->Branch("mctruth_bobbyphotonshower", &m_bobbyphotonshower);
		vertex_tree->Branch("mctruth_bobbypi0daughter", &m_bobbypi0daughter);
		vertex_tree->Branch("mctruth_bobbyotherdaughter", &m_bobbyotherdaughter);
		vertex_tree->Branch("mctruth_bobbyoverlay", &m_bobbyoverlay);
		vertex_tree->Branch("mctruth_bobbydeltaradppdaughter", &m_bobbydeltaradppdaughter);
		vertex_tree->Branch("mctruth_bobbydeltaradpdaughter", &m_bobbydeltaradpdaughter);
		vertex_tree->Branch("mctruth_bobbydeltaradmdaughter", &m_bobbydeltaradmdaughter);
		vertex_tree->Branch("mctruth_bobbydeltarad0daughter", &m_bobbydeltarad0daughter);

		vertex_tree->Branch("parameter_dist_tt",&m_dist_tt);
		vertex_tree->Branch("parameter_dist_sx",&m_dist_sx);
		vertex_tree->Branch("parameter_dist_st",&m_dist_st);
		vertex_tree->Branch("parameter_dist_sst",&m_dist_sst);
		//-----------------------------------------------------------
		
        vertex_tree->Branch("reco_slice_objects", &m_reco_slice_objects, "reco_slice_objects/I");

        this->CreateIsolationBranches();

        this->CreateSecondShowerBranches();
        // --------------------- Flash Related Variables ----------------------
        this->CreateFlashBranches();

        // --------------------- Track Related variables ------------
        this->CreateTrackBranches();

        std::string gpvm_location ="/pnfs/uboone/resilient/users/markross/tars/";


        //Get the info for length->energy conversion from PSTAR database.
        TFile *fileconv;
        struct stat buffer;   

        if(stat("proton_conversion.root", &buffer) == 0){
            fileconv = new TFile("proton_conversion.root", "read");
        }else{
            fileconv = new TFile((gpvm_location+"proton_conversion.root").c_str(), "read");
        }

        proton_length2energy_tgraph = *(TGraph*)fileconv->Get("Graph");
        proton_length2energy_tgraph.GetMean();
        fileconv->Close();

        // --------------------- Shower Related variables ------------
        this->CreateShowerBranches();


        //Metadata Branches
        this->CreateSliceBranches();
        //this->CreateMatchedSliceBranches();


        // ---------------------- MCTruth Related Variables ----------
        this->CreateMCTruthBranches();

        // ---------------------- Eventweight CTruth Related Variables ---------
        this->CreateEventWeightBranches();


        //std::string bad_channel_file = "/pnfs/uboone/resilient/users/markross/tars/MCC9_channel_list.txt";

        std::string bad_channel_file = "MCC9_channel_list.txt";

        if(stat(bad_channel_file.c_str(), &buffer) != 0){
            bad_channel_file = gpvm_location+bad_channel_file;
        }

        std::ifstream bc_file(bad_channel_file);

        if (bc_file.is_open())
        {
            std::string line;
            while ( getline (bc_file,line) )
            {
                std::vector<int> res;
                std::istringstream iss(line);
                for(std::string s; iss >> s; )
                    res.push_back( std::stof(s));

                std::pair<int,int> t(res[0],res[1]);
                bad_channel_list_fixed_mcc9.push_back(t);
            }
            bc_file.close();
        }


        std::cout<<"SinglePhoton \t||\t beginJob() is complete"<<std::endl;

    }




    //-------------------------------------------------------------------------------------------
    void SinglePhoton::ClearVertex(){

        //------------ Event related Variables -------------
        m_event_number = -99;
        m_subrun_number = -99;
        m_run_number = -99;
        m_test_matched_hits = 0;

        m_genie_spline_weight = 1.0;

        //------------ Vertex related Variables -------------
        m_reco_vertex_size = 0;
        m_vertex_pos_x=-99999;
        m_vertex_pos_y=-99999;
        m_vertex_pos_z=-99999;
        m_vertex_pos_tick=-9999;
        m_vertex_pos_wire_p0=-9999;
        m_vertex_pos_wire_p1=-9999;
        m_vertex_pos_wire_p2=-9999;

		//---------------------- BobbyVertexBuilder -----------------
		m_bobbyvertex_pos_x=-9999;
		m_bobbyvertex_pos_y=-9999;
		m_bobbyvertex_pos_z=-9999;
        m_bobbyshowers = 0;
        m_bobbytracks = 0;
		m_bobbyvertexradiusv = {999};
		m_bobbyvertexradius = 999;
		m_bobbyvertex_pos_xv={-9999};
		m_bobbyvertex_pos_yv={-9999};
		m_bobbyvertex_pos_zv={-9999};

		m_bobbysameslicev = {true};
        m_bobbyshowersv = {0};
        m_bobbytracksv = {0};
		m_bobbyprotontrackv = {0};
		m_bobbyphotonshowerv = {0};
		m_bobbypi0daughterv = {0};
		m_bobbydeltaradppdaughterv = {0};
		m_bobbydeltaradmdaughterv = {0};
		m_bobbydeltaradpdaughterv = {0};
		m_bobbydeltarad0daughterv = {0};
		m_bobbyotherdaughterv = {0};
		m_bobbyoverlayv = {0};
		m_bobbytrackdaughter_pdg = {-999};
		m_bobbyshowerdaughter_pdg = {-999};

		m_bobbyprotontrack = 0;
		m_bobbyphotonshower = 0;
		m_bobbypi0daughter = 0;
		m_bobbyotherdaughter = 0;
		m_bobbyoverlay = 0;
		m_bobbydeltaradppdaughter = 0;
		m_bobbydeltaradmdaughter = 0;
		m_bobbydeltaradpdaughter = 0;
		m_bobbydeltarad0daughter = 0;
		m_dist_tt = {999};
		m_dist_sx = {999};
		m_dist_st = {999};
		m_dist_sst ={999};
		//-----------------------------------------------------------

        m_reco_vertex_to_nearest_dead_wire_plane0=-99999;
        m_reco_vertex_to_nearest_dead_wire_plane1=-99999;
        m_reco_vertex_to_nearest_dead_wire_plane2=-99999;

        m_reco_slice_objects = 0;


        this->ClearIsolation();

        this->ClearSecondShowers();
        //------------- Flash related Variables ------------------
        this->ClearFlashes();

        //------------- Track Related Variables -----------------
        this->ClearTracks();

        //------------- Track Related Variables -----------------
        this->ClearShowers();
        this->ClearMCTruths();

        //------------- EventWeight Related Variables -----------------
        this->ClearEventWeightBranches();



        //MetaData Related Varibles
        this->ClearSlices();


    }


    void SinglePhoton::beginSubRun(art::SubRun const & sr) {

        if(m_potLabel != ""){
            if(m_potLabel == "generator"){
                double this_pot =  sr.getValidHandle<sumdata::POTSummary>(m_potLabel)->totgoodpot;
                m_pot_count += this_pot;
                std::cout<<"SinglePhoton::beginSubRun()\t||\t SubRun POT: "<<this_pot<<" . Current total POT this file: "<<m_pot_count<<std::endl;
            }else{
                art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR875;
                if (sr.getByLabel("beamdata","bnbETOR875",potSummaryHandlebnbETOR875)){
                    m_pot_count += potSummaryHandlebnbETOR875->totpot;
                }
            }
        }

    }









    //-----------------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------

    void SinglePhoton::GetVertex(const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, const art::Ptr<recob::PFParticle> & particle ){

        if(m_is_verbose) std::cout<<"SinglePhoton::Getvertex()\t||\t Starting to analyze recob::Vertex\n";
        int n_vert =0;

        //std::cout<<"There are "<<pfParticlesToVerticesMap.count(particle)<<" verticies associated with this particle"<<std::endl;

        lar_pandora::PFParticlesToVertices::const_iterator vIter = pfParticlesToVerticesMap.find(particle);
        if (pfParticlesToVerticesMap.end() != vIter)
        {
            const lar_pandora::VertexVector &vertexVector = vIter->second;
            if (!vertexVector.empty())
            {
                if (vertexVector.size() !=1)
                    std::cout << " Warning: Found particle with more than one associated vertex " << "\n";

                const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
                double xyz[3] = {0.0, 0.0, 0.0} ;
                vertex->XYZ(xyz);

                n_vert++;
                //std::cout<<"Vertex!"<<"\t "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<"\n";

                m_vertex_pos_x = xyz[0];
                m_vertex_pos_y = xyz[1];
                m_vertex_pos_z = xyz[2];

                m_reco_vertex_to_nearest_dead_wire_plane0 = distanceToNearestDeadWire(0, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);
                m_reco_vertex_to_nearest_dead_wire_plane1 = distanceToNearestDeadWire(1, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);
                m_reco_vertex_to_nearest_dead_wire_plane2 = distanceToNearestDeadWire(2, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);


            }else{
                std::cout << " Error: vertexVector associated with this particle is empty " << "\n";
                std::cerr << " Error: vertexVector associated with this particle is empty " << "\n";
                //exit(0);

            }
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::Getvertex()\t||\t Finished. Found "<<n_vert<<" vertices.\n";
    }

	//DUMP
    void SinglePhoton::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap)
    {
        //std::cout<<"Filling pfParticleMap with from the handle with total number "<<pfParticleHandle->size()<<std::endl;
        for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
        {
            const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
            // std::cout<<"Adding PFP to pfParticleMap with pfp id  "<<pParticle->Self()<<std::endl;
            if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second)
            {
                throw cet::exception("SinglePhoton") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
            }
        }
    }


	/***********************
	 *
	 * GetFinalStatePFParticleVectors() - fill in crParticles and nuParticles.
	 *	m_reco_vertex_size (filled) - the size of neutrino vertex is determined here.
	 * *********************/
    void SinglePhoton::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, PFParticleVector &crParticles, PFParticleVector &nuParticles )
    {

        int found = 0;
        int primaries = 0;
        int full = 0;
		//everything from the pfParticleMap, and look inside the element one by one.
        for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
        {
            const art::Ptr<recob::PFParticle> pParticle(it->second);

            full++;
            // Only look for primary particles
            if (!pParticle->IsPrimary()) continue;

            // Check if this particle is identified as the neutrino
            const int pdg(pParticle->PdgCode());
            const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);


            primaries++;
            // If it is, lets get the vertex position
            if(isNeutrino){
                found++;
                this->GetVertex(pfParticlesToVerticesMap, pParticle );

            } else{ // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
                crParticles.push_back(pParticle);
                continue;
            }

            // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
            //       If this is not the case please handle accordingly
            if (!nuParticles.empty())
            {
                throw cet::exception("SinglePhoton") << "  This event contains multiple reconstructed neutrinos!";
            }

            // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
            for (const size_t daughterId : pParticle->Daughters())
            {
                if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                    throw cet::exception("SinglePhoton") << "  Invalid PFParticle collection!";

                nuParticles.push_back(pfParticleMap.at(daughterId));
            }
        }
        std::cout<<"SinglePhoton::GetFinalStatePFParticleVectors()\t||\t Found "<<primaries<<" primary PFParticles (out of "<<full<<") of which: "<<found<<" were neutrinos."<<std::endl;
        m_reco_vertex_size = found;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

	    double SinglePhoton::triangle_area(double a1, double a2, double b1, double b2, double c1, double c2){
        return fabs((a1*(b2-c2)+b1*(c2-a2)+c1*(a2-b2))/2.0);
    }

    int SinglePhoton::quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area){

        std::vector<double> z(n,0.0);

        TGraph2D *g = new TGraph2D(n,X,Y,&z[0]);
        TGraphDelaunay delan(g);
        delan.SetMarginBinsContent(0);
        delan.ComputeZ(0,0);
        delan.FindAllTriangles();
        (*num_triangles)=delan.GetNdt();

        //Grab the locations of all the trianges. These will be intergers referencing to position in X,Y arrays
        Int_t *MT = delan.GetMTried();
        Int_t *NT = delan.GetNTried();
        Int_t *PT = delan.GetPTried();

        (*area)=0.0;
        for(int i = 0; i<delan.GetNdt(); i++){
            (*area)+=triangle_area(X[MT[i]-1],Y[MT[i]-1],X[NT[i]-1],Y[NT[i]-1],X[PT[i]-1],Y[PT[i]-1]);
        }

        delete g;
        return 0;
    }

    int SinglePhoton::delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area){

        int n = hits.size();
        std::vector<double> C0,T0;
        std::vector<double> C1,T1;
        std::vector<double> C2,T2;
        size_t n_0=0;
        size_t n_1=0;
        size_t n_2=0;

        for(int i=0;i<n; i++){
            const art::Ptr<recob::Hit> hit = hits[i];
            switch(hit->View()){
                case 0:
                    C0.push_back((double)hit->Channel());         
                    T0.push_back(hit->PeakTime());         
                    n_0++;
                    break;
                case 1:
                    C1.push_back((double)hit->Channel());         
                    T1.push_back(hit->PeakTime());         
                    n_1++;
                    break;
                case 2:
                    C2.push_back((double)hit->Channel());         
                    T2.push_back(hit->PeakTime());         
                    n_2++;
                    break;
                default:
                    break;
            }
        }
        if(m_use_delaunay){
            if(n_0>0) this->quick_delaunay_fit(n_0, &C0[0]  , &T0[0]  , &num_triangles[0],&area[0]);
            if(n_1>0) this->quick_delaunay_fit(n_1, &C1[0]  , &T1[0]  , &num_triangles[1],&area[1]);
            if(n_2>0) this->quick_delaunay_fit(n_2, &C2[0]  , &T2[0]  , &num_triangles[2],&area[2]);
        }
        num_hits[0] = n_0;
        num_hits[1] = n_1;
        num_hits[2] = n_2;

        //std::cout<<"Plane 0: "<<n_0<<" hits with "<<num_triangles[0]<<" triangles of area: "<< area[0]<<std::endl;
        //std::cout<<"Plane 1: "<<n_1<<" hits with "<<num_triangles[1]<<" triangles of area: "<< area[1]<<std::endl;
        //std::cout<<"Plane 2: "<<n_2<<" hits with "<<num_triangles[2]<<" triangles of area: "<< area[2]<<std::endl;

        return 0;
    }

    int SinglePhoton::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input){
        corrected.resize(3);

        double kx = input[0];
        double ky = input[1];
        double kz = input[2];

        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();

        double xtimeoffset = theDetector->ConvertTicksToX(g4Ticks,0,0,0);

//        double xOffset = -scecorr.X() +xtimeoffset+0.6;
        double yOffset = scecorr.Y();
        double zOffset = scecorr.Z();

        corrected[0]=kx - scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
        corrected[1]=ky+yOffset;
        corrected[2]=kz+zOffset;

        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<kx<<" "<<xOffset<<" "<<theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)<<" "<<scecorr.X()<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<xOffset<<" "<<yOffset<<" "<<zOffset<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: mcp->T(): "<<mcparticle->T()<<" TPCG4Time2Tick(): "<<detClocks->TPCG4Time2Tick(mcparticle->T())<<". "<<theDetector->GetXTicksOffset(0,0,0)<<" "<<theDetector->TriggerOffset()<<std::endl;
        return 0;
    }





    int SinglePhoton::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected){
        corrected.resize(3);

        double kx = mcparticle->Vx();
        double ky = mcparticle->Vy();
        double kz = mcparticle->Vz();

        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();

        double xtimeoffset = theDetector->ConvertTicksToX(g4Ticks,0,0,0);

        //double xOffset = -scecorr.X() +xtimeoffset+0.6;
        double yOffset = scecorr.Y();
        double zOffset = scecorr.Z();

        corrected[0]=kx - scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
        corrected[1]=ky+yOffset;
        corrected[2]=kz+zOffset;

        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<kx<<" "<<xOffset<<" "<<theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)<<" "<<scecorr.X()<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<xOffset<<" "<<yOffset<<" "<<zOffset<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: mcp->T(): "<<mcparticle->T()<<" TPCG4Time2Tick(): "<<detClocks->TPCG4Time2Tick(mcparticle->T())<<". "<<theDetector->GetXTicksOffset(0,0,0)<<" "<<theDetector->TriggerOffset()<<std::endl;
        return 0;
    }










    int SinglePhoton::spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected){
        corrected.resize(3);
        //Space Charge Effect! functionize this soon.
        double kx = mcparticle.Vx();
        double ky = mcparticle.Vy();
        double kz = mcparticle.Vz();
        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks->TPCG4Time2Tick(mcparticle.T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();

        double xtimeoffset = theDetector->ConvertTicksToX(g4Ticks,0,0,0);
        
        corrected[0]=kx - scecorr.X() +xtimeoffset+0.6;
        corrected[1]=ky + scecorr.Y();
        corrected[2]=kz + scecorr.Z();
        return 0;
    }

    void SinglePhoton::CollectMCParticles(
		const art::Event &evt, 
		const std::string &label, //m_geantModuleLabel
		std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,
		std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> &particlesToTruth, 
		std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap){

        //    if (evt.isRealData())
        //      throw cet::exception("LArPandora") << " PandoraCollector::CollectMCParticles --- Trying to access MC truth from real data ";

        art::Handle< std::vector< simb::MCParticle>  > theParticles;
        evt.getByLabel(label, theParticles);

        if (!theParticles.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles " << std::endl;
        }

        art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, label);

        for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i)
        {//theParticles are pointers to all MCParticles.
		//theTruthAssns are MCTruth with geantModuleLabel
            const art::Ptr<simb::MCParticle> particle(theParticles, i);
            const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));

            truthToParticles[truth].push_back(particle);
            particlesToTruth[particle] = truth;
            MCParticleToTrackIdMap[particle->TrackId()] = particle;
        }

        std::cout<<"SinglePhoton::CollectMCParticles() \t||\t the number of MCParticles in the event is "<<theParticles->size()<<std::endl;
    }

    void SinglePhoton::CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector)
    {
        //    if (evt.isRealData())
        //      throw cet::exception("LArPandora") << " PandoraCollector::CollectSimChannels --- Trying to access MC truth from real data ";

        art::Handle< std::vector<sim::SimChannel> > theSimChannels;
        evt.getByLabel(label, theSimChannels);

        if (!theSimChannels.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find sim channels... " << std::endl;
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theSimChannels->size() << " SimChannels " << std::endl;
        }

        for (unsigned int i = 0; i < theSimChannels->size(); ++i)
        {
            const art::Ptr<sim::SimChannel> channel(theSimChannels, i);
            simChannelVector.push_back(channel);
        }
    }


    void SinglePhoton::BuildMCParticleHitMaps(
			const art::Event &evt, 
			const std::string &label, 
			const std::vector<art::Ptr<recob::Hit>> &hitVector,   
			std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         
			std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> > &hitsToParticles, 
			const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, 
			std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap)
    {
        std::vector< art::Ptr<sim::SimChannel> >   simChannelVector;
        std::map< art::Ptr<simb::MCTruth>,     std::vector<art::Ptr<simb::MCParticle>>  >    truthToParticles;
        std::map< art::Ptr<simb::MCParticle>,  art::Ptr<simb::MCTruth> > particlesToTruth;
        std::map< art::Ptr<recob::Hit>,    std::vector< sim::TrackIDE >    >               hitsToTrackIDEs;

        this->CollectSimChannels(evt, label, simChannelVector);
        this->CollectMCParticles(evt, label, truthToParticles, particlesToTruth, MCParticleToTrackIdMap);
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(hitVector, simChannelVector, hitsToTrackIDEs);
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(hitsToTrackIDEs, truthToParticles, particlesToHits, hitsToParticles, daughterMode);


    }


} //namespace
