#ifndef __ATLAS_H__
#define __ATLAS_H__

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "SinglePhoton_module.h"

/*
 * Purpose of this header file, get the skeleton of all the maps and
 *	leave the contents to be filled inside the actual functions in
 *	other  header files.
 *
 */

namespace single_photon
{		/***************
	 * A class that designed for storing addresses for all associated (to an event) tracks, showers, 
	 *		and their cooresponding PFParticles.
	 *	evt (input) - the event that we currently look at.
	 * *************/

//	typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
//	typedef std::vector< art::Ptr<recob::Track> > TrackVector;
//	typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
//	typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;

	class Atlas{//Initialize this with event address;
		friend class SinglePhoton;//private as default

		public:
		//Constructor
		Atlas();//this might be useless;
		//Overflow constructor 1;
		Atlas( const art::Event &evt,
			std::vector< std::string > labels,
			bool is_data);//this initialize vectors and maps below

//		~Atlas(){
//			MCParticleToTrackIdMap.clear();
//		};
//		//main stuffs that we feed into the vertex builder.

/*
 * The constructor takes care of the following basic recob objects.
 *
 */
		//the following recob objects are created through the overflow constructor 1 (see above);
		std::vector< art::Ptr<recob::PFParticle> >	particles;//this is loaded depends on the option m_run_all_pfps, configurable in the .fcl file.
		std::vector< art::Ptr<recob::PFParticle> >	all_pfparticles;
		std::vector< art::Ptr<recob::Track> >		all_tracks;
		std::vector< art::Ptr<recob::Shower> >		all_showers;
		std::vector< art::Ptr<recob::Hit> >			all_hits;
		std::vector< art::Ptr<recob::OpFlash> >		all_opflashes;
		//get the cluster handle for the dQ/dx calc
		std::vector< art::Ptr<recob::Cluster> >		all_clusters;
        std::vector<art::Ptr<recob::Track>>			kalmanTrackVector;
		std::vector<art::Ptr<recob::Slice>>			sliceVector;

		//MCTruth (only initialized when the sample is not data)
		std::vector<art::Ptr<simb::MCTruth>>   mcTruthVector;
		std::vector<art::Ptr<simb::MCParticle>> matchedMCParticleVector;

/*
 * The overload constructor 1 takes care of the following maps.
 *
 */
		std::map< size_t, art::Ptr<recob::PFParticle>>	IDToPFParticleMap;
//		std::map< art::Ptr<recob::PFParticle, size_t>>	PFParticleToIDMap;//This makse more consistant, but it needs works!
		std::map< art::Ptr<recob::PFParticle> , std::vector<art::Ptr<recob::Vertex>> > PFParticlesToVerticesMap;//old name pfParticlesToVerticesMap;
		std::map< art::Ptr<recob::PFParticle> , std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > PFParticleToMetadataMap;//old name pfParticleToMetadataMap;
		std::map< art::Ptr<recob::PFParticle> , std::vector<art::Ptr<recob::SpacePoint>> > PFParticleToSpacePointsMap;
		std::map< art::Ptr<recob::PFParticle> , std::vector<art::Ptr<recob::Cluster>> > PFParticleToClustersMap;
		std::map< art::Ptr<recob::Cluster>	  , std::vector<art::Ptr<recob::Hit>> > ClusterToHitsMap;
        std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>> PFParticlesToShowerReco3DMap;
        std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Track>> PFParticlesToShowerKalmanMap;
        std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>>> kalmanTrackToCaloMap;
		std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::PFParticle>> > sliceToPFParticlesMap;
		std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::Hit>> > sliceToHitsMap;
		std::map<int, std::vector<art::Ptr<recob::PFParticle>> > sliceIDToPFParticlesMap;
		std::map<int, std::vector<art::Ptr<recob::Hit>> > sliceIDToHitsMap;
/*
 * Initially empty variables to be filled from other parts of the code.
 *
 */
//The followings are taken care by the CollectTracksAndShowers_v2() in BobbyVertexBuilder.h
		std::vector< art::Ptr<recob::Track> >							selected_tracks;
		std::vector< art::Ptr<recob::Shower> >							selected_showers;
		std::vector< art::Ptr<recob::Track> >							more_tracks;//non-cosmic objects, but not selected nu objects.
		std::vector< art::Ptr<recob::Shower> >							more_showers;
		//Maps for more pandora objects.
		std::map< art::Ptr<recob::Track>  , art::Ptr<recob::PFParticle>> trackToNuPFParticleMap;
		std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap;

		std::map< art::Ptr<recob::Track>  , double > trackToDistMap;
		std::map< art::Ptr<recob::Shower> , double > showerToDistMap;

//The followings are taken care by the AnalyzeSlices() in analyze_Slice.h
		std::map<int, double>											sliceIdToNuScoreMap;
		//Pairs that connect PFParticle to sliceID.
		std::vector<std::pair<art::Ptr<recob::PFParticle>,int>>			primaryPFPSliceIdVec;
		std::map< art::Ptr<recob::PFParticle>, int> 					PFPToSliceIdMap;
		std::map< art::Ptr<recob::PFParticle>, bool>					PFPToClearCosmicMap;
		std::map< art::Ptr<recob::PFParticle>, bool> 					PFPToNuSliceMap;
		std::map< art::Ptr<recob::PFParticle>, double>					PFPToTrackScoreMap;

//Maps for simb objects; for MC sample, not applied to data

	
		//Filled in the CollectMCParticles_v2() in Singlephoton_module.cc
		std::map<int, art::Ptr<simb::MCParticle>>									MCParticleToTrackIdMap;
		std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> MCTruthToMCParticlesMap;
		std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>              MCParticleToMCTruthMap;
		std::map< art::Ptr<simb::MCParticle>, int>            MCParticleToAncestorPdgMap;

		//Filled in the showerRecoMCmatching() in reco_truth_matching.h
//		std::vector<art::Ptr<simb::MCParticle>>										matchedMCParticleVector;
		std::map< art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> >				showerToMCParticleMap;

		//Filled in the RecoMCTracks() in analyze_Tracks.h
		std::map< art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> >				trackToMCParticleMap;


		//FindManyP's!
		//specially for the number of coorresponding recob (pandora_objects) to a PFParticle;
		// example: PFParticleIsATrack[particles.key()] gives the vector that containts all 
		//		cooresponding tracks;
		art::FindManyP< recob::Track	>*								PFParticleAsATrack;
		art::FindManyP< recob::Shower	>*								PFParticleAsAShower;

//		art::FindManyP< recob::Track	> PFParticleAsATrack(PFParticleHandle, evt, m_trackLabel);
//		art::FindManyP< recob::Shower	> PFParticleAsAShower(PFParticleHandle, evt, m_showerLabel);
		
		private:
		/***********************
		 *
		 * Wigets for constructing the different types of varaibels.in this class
		 * 1. HandleToVector;
		 *
		 * **********************/
		//1. Scratch from Handle, and return a equivalently useful vector.
		//sample usage:
		//		art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitfinderLabel); //yea, this should be gone;
		//		recob::Hit dummy_hit;//This is to specify the template;
		//		std::vector<art::Ptr<recob::Hit>> hitVector = HandleToVector(dummy_hit, evt, m_hitfinderLabel);
		template <typename recob_object>//A helper template that allows you to make compliated types.
			struct temporary_types{ 
				using type1 = std::vector<art::Ptr<recob_object>>;
				using type2 = art::ValidHandle<std::vector<recob_object>>;
				using type3 = std::vector<recob_object>;
			};
		template <class recob_object>//ref_type is only used to identify the temporary_types late.
			typename temporary_types<recob_object>::type1 HandleToVector(recob_object ref_type, const art::Event &evt, std::string &label){

				typename temporary_types<recob_object>::type2 const & Handle = evt.getValidHandle<typename temporary_types<recob_object>::type3>(label);
				typename temporary_types<recob_object>::type1 Vector;
				art::fill_ptr_vector(Vector,Handle);
				return Vector;
			}

		//2. Wiget 2 comes here..

	};
	

//sth similar to class DetectorObjects

//create struct for Track and Shower to add more detail to them, i.e. the geoalgo class;

//struct Track
//struct Shower








//THE HEADER FILE IS THE ABOVE



	//Constructor
	Atlas::Atlas (){}
	//Overloaded Constructor 1, initialize the essential variables
	Atlas::Atlas ( const art::Event &evt,
				std::vector<std::string > labels,
				bool is_data){

		//PREPARE some recob objects;
		//vector<string> labels = {m_trackLabel, m_showerLabel, m_hitfinderLabel, m_flashLabel, m_pandoraLabel,m_shower3dLabel,m_showerKalmanLabel,m_showerKalmanCaloLabel,m_generatorLabel, m_geantModuleLabel}
		recob::PFParticle dummy_PFParticle;
		all_pfparticles = HandleToVector(dummy_PFParticle, evt, labels[4]);
		recob::Track dummy_track;//This is to specify the template;
		all_tracks = HandleToVector(dummy_track, evt, labels[0]);//m_trackLabel
		kalmanTrackVector = HandleToVector(dummy_track,evt,labels[6]);//m_showerKalmanLabel

		recob::Shower dummy_shower;//This is to specify the template;
		all_showers = HandleToVector(dummy_shower, evt, labels[1]);//m_showerLabel
		recob::Hit	dummy_hit;
		all_hits = HandleToVector(dummy_hit, evt, labels[2]);//m_hitfinderLabel

		recob::OpFlash	dummy_opflash;
		all_opflashes = HandleToVector(dummy_opflash, evt, labels[3]);//m_flashLabel

		recob::Cluster	dummy_cluster;
		all_clusters = HandleToVector(dummy_cluster, evt, labels[4]);//m_pandoraLabel

		recob::Slice dummy_slice;
		sliceVector = HandleToVector(dummy_slice, evt, labels[4]);
		
		if(!is_data){
			//MCTruth Handle
			simb::MCTruth dummy_geant; //for Geant info.
			mcTruthVector = HandleToVector(dummy_geant, evt, labels[8]);

			simb::MCParticle dummy_genie; //for Genie info.
			matchedMCParticleVector = HandleToVector(dummy_genie, evt, labels[9]);
		}
		//CREATE maps!
		//Ingredient 1: Handles; I temporary define it here for mapping purpose;
		art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(labels[4]);//This is useful for FindManyP< reco::Track/Shower>
		art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(labels[4]);
        art::ValidHandle<std::vector<recob::Track>> const & kalmanTrackHandle  = evt.getValidHandle<std::vector<recob::Track>>(labels[6]);
		art::ValidHandle<std::vector<recob::Slice>> const & sliceHandle  = evt.getValidHandle<std::vector<recob::Slice>>(labels[4]);
		//a cross check
		if (!pfParticleHandle.isValid())
		{
			mf::LogDebug("SinglePhoton") << "  Failed to find the PFParticles.\n";
			return;
		}

	
		//Ingredient 2: FindManyPs; these will be gone when construction finished
		art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt,  labels[4]);
		art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, labels[4]);
		art::FindManyP<recob::SpacePoint> spacePoints_per_pfparticle(pfParticleHandle, evt, labels[4]);
		art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfParticleHandle, evt, labels[4]);
		art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, labels[4]);
        art::FindOneP<recob::Shower> showerreco3D_per_pfparticle(pfParticleHandle, evt, labels[5]);
        art::FindOneP<recob::Track> showerKalman_per_pfparticle(pfParticleHandle, evt, labels[6]);
        art::FindManyP<anab::Calorimetry> cali_per_kalmantrack(kalmanTrackHandle, evt, labels[7]);
		art::FindManyP<recob::PFParticle> pfparticles_per_slice(sliceHandle, evt, labels[4]);
		art::FindManyP<recob::Hit> hits_per_slice(sliceHandle, evt, labels[4]);
		
		//make maps here;
		for(size_t i=0; i< all_pfparticles.size(); ++i){
			const art::Ptr<recob::PFParticle> pfp = all_pfparticles[i];
			PFParticlesToVerticesMap[pfp] = vertices_per_pfparticle.at(pfp.key());//old name: pfParticlesToVerticesMap;

			PFParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
			IDToPFParticleMap[pfp->Self()] = pfp;
			PFParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());
			PFParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());

			//3D Showers
			if(!showerreco3D_per_pfparticle.at(pfp.key()).isNull()){
                PFParticlesToShowerReco3DMap[pfp] = showerreco3D_per_pfparticle.at(pfp.key());
             }
			//---------Kalman Track Showers
             if(!showerKalman_per_pfparticle.at(pfp.key()).isNull()){ 
                 PFParticlesToShowerKalmanMap[pfp] =showerKalman_per_pfparticle.at(pfp.key());
             }
		}
		//other maps not on pfParticles;
		for(size_t i=0; i< all_clusters.size(); ++i){
			auto cluster = all_clusters[i];
			ClusterToHitsMap[cluster] = hits_per_cluster.at(cluster.key());
		}

        //----- kalmon Cali
		for(size_t i=0; i< kalmanTrackVector.size(); ++i){
			auto trk = kalmanTrackVector[i];
			if(cali_per_kalmantrack.at(trk.key()).size()!=0){
				kalmanTrackToCaloMap[trk] =cali_per_kalmantrack.at(trk.key());
			}
		}

		for(size_t i=0; i< sliceVector.size(); ++i){
			auto slice = sliceVector[i];
			sliceToPFParticlesMap[slice] = pfparticles_per_slice.at(slice.key());
			sliceIDToPFParticlesMap[slice->ID()] = pfparticles_per_slice.at(slice.key());
			sliceToHitsMap[slice] =hits_per_slice.at(slice.key());
			sliceIDToHitsMap[slice->ID()] = hits_per_slice.at(slice.key());
		}

//------- trackToMCParticleMap --------
//		if(!is_data){
//			//Fill in trackToMCParticleMap;
//			//https://indico.fnal.gov/event/20453/session/6/contribution/4/material/slides/0.pdf
//			//BackTrackerService and ParticleInventoryService
//			art::ServiceHandle<cheat::BackTrackerService> bt_serv;
//			art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
//			//Get tracks
//			art::Handle< std::vector<recob::Track> > trackListHandle;
//			std::vector<art::Ptr<recob::Track> > tracklist;
//			if (evt.getByLabel(labels[0],trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
//			//Get hit-track association
//			art::FindManyP<recob::Hit> fmth(trackListHandle, evt, labels[0]);
//			//Loop over all tracks
//			for(size_t i=0; i<tracklist.size();++i){ 
//				if (fmth.isValid()){
//					// Find true track for each reconstructed track
//					int TrackID = 0;
//					//Get all hits associated with the track 
//					std::vector< art::Ptr<recob::Hit> > allHits = fmth.at(i);
//					std::map<int,double> trkide;
//					cout<<"Size of hits "<<allHits.size()<<endl;
//					for(size_t h = 0; h < allHits.size(); ++h){ 
//						art::Ptr<recob::Hit> hit = allHits[h]; 
//						//TrackIDE saves the energy deposition for each Geant particle ID
//						std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(hit);
//						cout<<"Size of TrackIDs "<<TrackIDs.size()<<endl;
//						for(size_t e = 0; e < TrackIDs.size(); ++e){ 
//						trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
//						}
//					}
//					// Work out which IDE despoited the most charge in the hit if there was more than one.
//					double maxe = -1;
//					double tote = 0;
//					for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
//							cout<<"Possible IDs "<<ii->first<<endl;
//						tote += ii->second;
//						if ((ii->second)>maxe){
//							TrackID = ii->first; //TrackID maxe = ii->second; //Energy
//						} }
//					// Now have trackID, so get PdG code.
//					
//				//	const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
//
//					cout<<"Ptr::MCP CHECK "<<TrackID;//NO working, TrackID is always 0;
//					
//				//	if(particle){
//					for(const art::Ptr<simb::MCParticle> mcp : matchedMCParticleVector){
//						//cout<<"Match ID "<<mcp->TrackId()<<endl;
//						if(mcp->TrackId() == TrackID){ cout<<"Yes"<<endl;}
//						else{ cout<<"NO"<<endl;}
//					}
//				//	}
//					cout<<"Ptr::MCP CHECK Finished"<<endl;
//
////				 	trackToMCParticleMap.emplace(tracklist[i], particle);
////					if (particle){
////						std::cout<<"Pdgcode = "<<particle-
////							>PdgCode()<<std::endl; }
//				} 
//			}
//		}
//	//CHECK
//		std::vector<simb::MCParticle const *> allhit_particle_vec;
//		std::vector<anab::BackTrackerHitMatchingData const *> allhit_match_vec;
//		std::vector<std::unordered_map<int, double>> allhit_trkide(3);
//
//		for (size_t i = 0; i < hitlist.size(); ++i) {
//			art::Ptr<recob::Hit> hit = hitlist[i];
//			allhit_particle_vec.clear();
//			allhit_match_vec.clear();
//			particles_per_hit.get(hit.key(), allhit_particle_vec, allhit_match_vec);
//			for (size_t i_p = 0; i_p < allhit_particle_vec.size(); ++i_p) {
//				allhit_trkide[hit->WireID().Plane][pi_serv->TrackIdToEveTrackId(allhit_particle_vec.at(i_p)->TrackId())] += allhit_match_vec[i_p]->energy; // match_vec[i_p]->numElectrons;
//			}
//		} // end of loop hitlist i

//-------------------

	}
//		//3D Showers
//        std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>> pfParticlesToShowerReco3DMap;
//        for(size_t i=0; i< pfParticleVector.size(); ++i){
//            auto pfp = pfParticleVector[i];
//             
//
//        }
//        std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Track>> pfParticlesToShowerKalmanMap;
//        for(size_t i=0; i< pfParticleVector.size(); ++i){
//            auto pfp = pfParticleVector[i];
//
//        }

//        art::ValidHandle<std::vector<recob::Track>> const & kalmanTrackHandle  = evt.getValidHandle<std::vector<recob::Track>>(m_showerKalmanLabel);
//        std::vector<art::Ptr<recob::Track>> kalmanTrackVector;
//        art::fill_ptr_vector(kalmanTrackVector,kalmanTrackHandle);
//        
//        art::FindManyP<anab::Calorimetry> cali_per_kalmantrack(kalmanTrackHandle, evt, m_showerKalmanCaloLabel);
//        std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>>> kalmanTrackToCaloMap;


}




#endif
