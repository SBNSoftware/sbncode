

#ifndef __BOBBYVERTEXBUILDER_H__
#define __BOBBYVERTEXBUILDER_H__

#include "SinglePhoton_module.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "Atlas.h"

//---VertexBuilder---
#include "VertexBuilder/VertexBuilder.h"
#include "VertexBuilder/ParticleAssociations.h" 
//#include "larsim/MCCheater/ParticleInventoryService.h"

//#include "VertexBuilder/DetectorObjects.h"
//
//#include "VertexQuality.h"
//#include "FillTreeVariables.h"
//#include "RecoMCMatching.h"
//-------------------

//#include "TCanvas.h"
//#include "TGraph.h"
//#include "TFile.h"
//#include "TAxis.h"
//#include "TLine.h"
//#include "TLatex.h"
//#include "TLegend.h"
//#include "TPrincipal.h"
//#include "TVectorD.h"
//#include "TMatrixD.h"
//#include "TF1.h"
//#include "TEllipse.h"

//the following is not working yet
//#include <curses.h> //Use LINES and COLS to get the dimensions of the screen

using namespace std;

namespace single_photon
{
	//Printout Gadgets
	void verbose_check( bool fverbose, bool debug_message ){

		std::cout<<" Verbose output is ";
		if(fverbose){
			std::cout<<"on;";
		}else{
			std::cout<<"off;";
		}
		std::cout<<" Debug message is ";
		if(debug_message){
			std::cout<<"on.";
		}else{
			std::cout<<"off.";
		}
		std::cout<<std::endl;
	}
	/*****************************
	 * CollectTracksAndShowers_v2() - Fill in package.selected_showers and
	 *	package.selected_tracks as well as other nu-like objects (more_tracks/showers)
	 *
	 *  packahe (to be modified) - contains all vectors and maps that are
	 *		needed to access pandora_objects; it will be modified as
	 *		the code runs.
	 *  //not using run_count here now;
	 *	run_count (input) - has values 0,1,2; each number indicates differnt collections of slices: 
	 *		0 - selected nu_slice; 
	 *		1 - not-selected nu_slices; 
	 *		2 - cosmic slices (non-nu_slices).  NOT DOING THIS
	 *
	 * **************************/


	void SinglePhoton::CollectTracksAndShowers_v2(
			const art::Event &evt,
			class Atlas &package){
		bool fverbose = true;//m_is_verbose;
		bool debug_message = true;

		std::cout<<"SinglePhoton::CollectTracksAndShowers_v2() \t||\t";
		verbose_check( fverbose, debug_message);

		std::vector< art::Ptr<recob::PFParticle> > candidate_particles(package.particles);//CHECK, always vertex the most nu-like slice.

		//we dont need the run_all_particle boolean, because it is considered under package.particles;
		if(candidate_particles.size() == 0){
			mf::LogDebug("SinglePhoton") << "  No PFParticles can be considered for vertexing.\n";
		}

		for(const art::Ptr<recob::PFParticle> &pParticle : candidate_particles){ //candidate_particles are determined in above if-else statement.
			const std::vector< art::Ptr<recob::Track> > ToBeAddedTracks((package.PFParticleAsATrack)->at(pParticle.key()));
			const std::vector< art::Ptr<recob::Shower> > ToBeAddedShowers((package.PFParticleAsAShower)->at(pParticle.key()));

			const unsigned int nTracks(ToBeAddedTracks.size());
			const unsigned int nShowers(ToBeAddedShowers.size());

			if( nTracks + nShowers == 0 ){//Um... the PFParticle is not identified as track or shower;
				mf::LogDebug("SinglePhoton") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << "\n";
			}else if( nTracks + nShowers > 1 ){
				//Check, do I need throw??
				throw cet::exception("SinglePhoton") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
				//Ok, if going through the below if-elses, it means we have a shower/ a track!

			}else if( nTracks == 1 ){ //Add a Track
				if(package.PFPToNuSliceMap[pParticle]){//add the package.selected nu_slice particle;
					(package.selected_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(package.selected_tracks).back()]= pParticle;

				}else{//add other nu-like PFParticle;
					(package.more_tracks).push_back(ToBeAddedTracks.front());
					(package.trackToNuPFParticleMap)[(package.more_tracks).back()]= pParticle;
				}

				if(fverbose){
					std::cout<<"SinglePhoton::CollectTracksAndShowers_v2() \t||\t";
					std::cout<<"adding to trackToNuPFParticleMap this track with id "<<  ToBeAddedTracks.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;
				}

			}else if( nShowers == 1 ){ //Add a Shower
				if(package.PFPToNuSliceMap[pParticle]) {//add the package.selected nu_slice particle;
					(package.selected_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(package.selected_showers).back()] = pParticle;
				}else{//add other nu-like PFParticle
					(package.more_showers).push_back(ToBeAddedShowers.front());
					(package.showerToNuPFParticleMap)[(package.more_showers).back()]= pParticle;
				}

				if(fverbose){
					std::cout<<"SinglePhoton::CollectTracksAndShowers_v2() \t||\t";
					std::cout<<"adding to showerToNuPFParticleMap this shower with id "<<  ToBeAddedShowers.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;
				}
			} //repeat this loop for another PFParticle in the candidate_particles.
		}
	}

	/*****************************
	 * TraceMothers - based on input MCParticle, return its mother() or mother's mother,depends on loops
	 * particle (input) - input MCParticle;
	 * map (input) - TrackId, MCPartile Map
	 * target_pdg (input) - look for target pdg;
	 * loops (modified) - if loops goes to 0, then stop searching for parent.
	 * **************************/

	simb::MCParticle TraceMothers( simb::MCParticle particle, std::map<int, simb::MCParticle >& map, vector <int> target_pdg, int loops){
		int current_pdg = (particle).PdgCode();
		int mother_id = (particle).Mother();
		simb::MCParticle mother_particle = map[mother_id];
		bool match_targetpdg = false;

		for(int& i :target_pdg){
			if(current_pdg == i){ 
				match_targetpdg = true;
				break;
			}
		}
		loops--;
		if(loops < 0 || (current_pdg != mother_particle.PdgCode() && match_targetpdg) || mother_id < 0){//stop searching mother;
			return particle;
		}

		return TraceMothers(mother_particle, map, target_pdg,loops);
	}

	art::Ptr< simb::MCParticle > TraceMothers( art::Ptr< simb::MCParticle > particle, std::map<int, art::Ptr< simb::MCParticle > >& map, vector <int> target_pdg, int loops){
		int current_pdg = particle->PdgCode();
		int mother_id = particle->Mother();
		art::Ptr< simb::MCParticle > mother_particle = map[mother_id];
		bool match_targetpdg = false;

		for(int& i :target_pdg){
			if(current_pdg == i){ 
				match_targetpdg = true;
				break;
			}
		}
		loops--;
		if(loops < 0 || (current_pdg != mother_particle->PdgCode() && match_targetpdg) || mother_id < 0){//stop searching mother;
			return particle;
		}

		return TraceMothers(mother_particle, map, target_pdg,loops);
	}


    void SinglePhoton::CollectMCParticles_v2(
		const art::Event &evt, 
		class Atlas & package){

		int screen_width = 86;
		bool debug_message = true;
		bool fverbose = m_is_verbose;

		std::cout<<"SinglePhoton::CollectMCParticles_v2() \t||\t";
		verbose_check( fverbose, debug_message);

        art::Handle< std::vector< simb::MCParticle>  > theParticles;
        evt.getByLabel(m_geantModuleLabel, theParticles);

        if (!theParticles.isValid()){
            mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
            return;
        } else {
            mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles " << std::endl;
        }

        art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, m_geantModuleLabel);
		package.MCParticleToTrackIdMap.clear();

		// Strategy: -----------------------
		// 1. Prepare store_MCTruths, which contain MCParticles from MCTruths 
		//		that could be possibly mapped to the MCParticles from Geant;
		// 2. Examine all MCParticles at Geant level, and find their ancestor Pdg, 
		//		i.e. Pdg obtained from MCParticles from Genie;
		//		fill in MCParticleToAncestorPdgMap;
		
		//doing 1. Prepare store_MCTruths.
		std::map < int, simb::MCParticle > MCT_trackIdToDaughterMap;//MCTruth that can be mapped by Geant MCP;
		std::map < int, simb::MCParticle > MCT_trackIdToMCParticleMap; //MCTruth that can provide ancestor Pdg to MCPmother;
		vector<int> looking_for_mothers = {111, 1114, 2114, 2214, 2224, -111, -1114, -2114, -2214, -2224};

		for(size_t i = 0; i<package.mcTruthVector.size(); i++){
			bool found_ancestor = false; //find pi0 or delta.
			for(int j = 0; j<package.mcTruthVector[i]->NParticles(); j++){
				//collect MCTruth candidate;
				simb::MCParticle check_mct = package.mcTruthVector[i]->GetParticle(j);
				MCT_trackIdToMCParticleMap[ check_mct.TrackId()] = check_mct;

				if(debug_message){
					std::cout<<" \n MCTruth Track Id:"<<setw(4)<<check_mct.TrackId(); //gives the id of the MCParticle
					std::cout<<" Mother ID:"<<setw(3)<<check_mct.Mother(); 
					std::cout<<" Pdg Code:"<<setw(10)<<check_mct.PdgCode(); 
					std::cout<<" Status Code:"<<setw(4)<<check_mct.StatusCode(); 
					cout<<"; An ancestor MCParticle with energy: "<<setw(8)<<check_mct.E();
				}
				int current_pdg = abs(check_mct.PdgCode());
				switch(current_pdg){
					case 12:
					case 14:
					case 16:
						continue;//do nothing on neutrinos
						break;
					case 111://pi0
						{
							if(!found_ancestor){
								if(debug_message) cout<<" Identify this Pi0 Ancestor. Store this Pi0 Daughter.";
								found_ancestor = true;
							}
							MCT_trackIdToDaughterMap.emplace(check_mct.TrackId() ,check_mct );
							simb::MCParticle mother_mct = TraceMothers(check_mct, MCT_trackIdToMCParticleMap, {111,-111} , 4);
							if( check_mct.TrackId() != mother_mct.TrackId() && check_mct.E() - mother_mct.E()<0.005){//oops, duplicated pion;
								MCT_trackIdToDaughterMap.erase(mother_mct.TrackId());
							}else{
								if(debug_message) cout<<" Store this Pi0 Daughter";
							}

							break;
						}
					case 1114://delta- 
					case 2114://delta0
					case 2214://delta+
					case 2224://delta++
						//							MCT_trackIdToAncestorMap.emplace(check_mct.TrackId() ,check_mct );
						if(debug_message) cout<<"Identify this Ancestor ";
						found_ancestor = true;
						break;
					case 22://photon
						//						case 211://pi+
					case 2112://neutron
					case 2212://proton
						if(found_ancestor && check_mct.StatusCode() == 1){
							MCT_trackIdToDaughterMap.emplace(check_mct.TrackId() ,check_mct );
							if(debug_message) cout<<" Store this Daughter ";
						}
						break;
					default:
						{}
				}
			}
		}//End of looping all MCParticles from MCTruth; the map store_MCTruths are prepared.
		cout<<endl;

		//doing 2. Map MCParticle to its mother(), Geant level;
		if(debug_message) std::cout<<"Identifying MCParticles ancestors:"<<endl;

//		std::vector< simb::MCParticle> store_InvisibleMCTruths;//MCTruth that can be mapped to MCPmother;
//		std::vector< simb::MCParticle> store_VisibleMCTruths;//MCTruth that can be mapped to MCPmother;

		std::map < int, int > MCP_trackIdToPi0PdgMap; //Record trackId of ancestor Pi0 from MCPrticle;
		
		if(false){
			cout<<"This is a preview"<<endl;
			for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i){//theParticles are pointers to all MCParticles.
				const art::Ptr<simb::MCParticle> particle(theParticles, i);
				if(debug_message){
					std::cout<<"\n MCParticle Track Id:"<<setw(4)<<particle->TrackId(); //gives the id of the MCParticle
					std::cout<<" Mother ID:"<<setw(3)<<particle->Mother(); 
					std::cout<<" Pdg Code:"<<setw(10)<<particle->PdgCode(); 
					std::cout<<" Status Code:"<<setw(4)<<particle->StatusCode(); 
					cout<<"; energy: "<<setw(14)<<particle->E();
				}
			}
			cout<<"End of Previse"<<endl;
		}

		for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i){//theParticles are pointers to all MCParticles.
			//theTruthAssns are MCTruth with geantModuleLabel
			const art::Ptr<simb::MCParticle> particle(theParticles, i);
			const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));

			package.MCTruthToMCParticlesMap[truth].push_back(particle);
			package.MCParticleToMCTruthMap[particle] = truth;//truth ->GetParticles() will give a bunch of particles that has not been linked.
			package.MCParticleToTrackIdMap[particle->TrackId()] = particle;

			//---------- Find out MCParticles that would be mapped to MCTruth Pdg, i.e. particles produced by Genie
			int mother_track = particle->Mother();	
			int ancestor_pdg = particle->PdgCode();//itself pdg by default
			//pdg==22 is for pi0 search; Mother() is for both pi0 and Delta;
//			cout<<"\n   # of daughter is left: "<<MCT_trackIdToDaughterMap.size()<<" Id: ";
//			for(auto it = MCT_trackIdToDaughterMap.begin(); it!=MCT_trackIdToDaughterMap.end();it++){
//				cout<<it->first<<", ";
//			}
			if(false){
				std::cout<<"\n MCParticle Track Id:"<<setw(4)<<particle->TrackId(); //gives the id of the MCParticle
				std::cout<<" Mother ID:"<<setw(3)<<particle->Mother(); 
				std::cout<<" Pdg Code:"<<setw(10)<<particle->PdgCode(); 
				std::cout<<" Status Code:"<<setw(4)<<particle->StatusCode(); 
				cout<<"; energy: "<<setw(14)<<particle->E();
			}

			//true, fill the MCParticleToAncestorPdgMap with ancestor pdg; 
			//false, fill with self pdg;

			if(particle->Mother() == 0){//the goal is to change ancestor_pdg, and erase an element in MCT_trackIdToDaughterMap.
				for(auto const & [a,b] : MCT_trackIdToDaughterMap){//map the ancestor pdg to each MCParticle
					bool match_pdg = ( b.PdgCode() == ancestor_pdg );
					bool match_energy = ( b.E() - particle->E() < 0.0005);
					if(match_pdg && match_energy){//woo, match to the MCTruth! now look for ancestor!

//						cout<<"\n Trace (candidate MCTruth daughter) TrackId "<<particle->TrackId()<<endl;
						simb::MCParticle mother_mct = TraceMothers( b, MCT_trackIdToMCParticleMap, looking_for_mothers, 4);
						cout<<"Trace found mother, TrackId "<<mother_mct.TrackId()<<" with Pdg "<<mother_mct.PdgCode()<<endl;
						int temp_mcpdg = abs(particle->PdgCode());
						cout<<"this pgd is "<<temp_mcpdg<<endl;
						switch(temp_mcpdg){//the following particles need new ancestor pdg;
							case 111://pi0
								MCP_trackIdToPi0PdgMap.emplace(particle->TrackId(), particle->PdgCode());
								if(debug_message) cout<<" Mark down pi0 Track Id "<<particle->TrackId()<<endl;
								break;
							case 22://gamma
								mother_mct = TraceMothers( b, MCT_trackIdToMCParticleMap, {111, -111}, 1);
							case 211://pi+
							case 2112://neutron
							case 2212://proton
								ancestor_pdg = mother_mct.PdgCode();
								break;
							default://ancestor_pdg is unchanged, the particle is not listed above;
								{}
						}

						MCT_trackIdToDaughterMap.erase(a);
						break;
					}
				}
			}else if( MCP_trackIdToPi0PdgMap.find(mother_track) != MCP_trackIdToPi0PdgMap.end() ){//a photon that could come from the pi0!
				if(debug_message) cout<<"Identify photons from pi0"<<endl;
				ancestor_pdg = MCP_trackIdToPi0PdgMap[mother_track];
				//				cout<<"\n # Pi0 ancestor: "<<MCT_trackIdToPi0AncestorMap.size()<<" Trace (Photon) TrackId "<<particle->TrackId()<< " Pdg: "<<particle->PdgCode();
				//				art::Ptr< simb::MCParticle > mother_mcp = TraceMothers(particle, package.MCParticleToTrackIdMap, {111, -111} , 1);
				//				cout<<" Found mother, TrackId "<<mother_mcp->TrackId()<<" with Pdg "<<mother_mcp->PdgCode()<<" energy: "<<mother_mcp->E();
				//				int mother_pdg = mother_mcp->PdgCode();
				//				for(auto iterator = MCT_trackIdToPi0AncestorMap.begin(); iterator != MCT_trackIdToPi0AncestorMap.end(); ++iterator){
				//					bool match_pdg = ( (iterator->second).PdgCode() == mother_pdg );
				//					bool match_energy = ( (iterator->second).E() - mother_mcp->E() < 0.0005);
				//					if(match_pdg && match_energy){//woo, match to the MCTruth! this photon is from MCTruth pi0!
				//						ancestor_pdg = mother_pdg;
				//						break;
				//					}
				//				}
			}

			package.MCParticleToAncestorPdgMap.emplace(particle, ancestor_pdg);
		}//End of looping all MCparticles; store_MCPcandidates are prepared.

		//---
		//Below is to check the result;
		//---

		if(fverbose){
			std::cout<<"\nSinglePhoton::CollectMCParticles_v2() \t||\t ";
			std::cout<<"the number of all MCParticles (not necessarily in use) in the event is "<<theParticles->size()<<std::endl;
		}

		if(debug_message){
			cout<<"Take a look at the AncestorMap size:"<<package.MCParticleToAncestorPdgMap.size()<<endl;
			for(int i = 0; i<screen_width; i++)
				cout<<"-";
			cout<<endl;

			cout<<setw(10)<<" Track Id ";
			cout<<setw(14)<<left<<" Pdg ";
			cout<<setw(9)<<left<<" Mother ";
			cout<<setw(15)<<"| Ancestor Pdg "<<endl;

			for(int i = 0; i<screen_width; i++) 
				cout<<"-";
			cout<<endl;

			for(auto const & [a,b] : package.MCParticleToAncestorPdgMap){//map the ancestor pdg to each MCParticle

				cout<<setw(10)<<a->TrackId();
				cout<<" "<<setw(12)<<left<<a->PdgCode();
				cout<<" "<<setw(7)<<left<<a->Mother()<<" | ";
				cout<<setw(5)<<left<<b<<endl;
			}

			for(int i = 0; i<screen_width; i++) 
				cout<<"-";
			cout<<endl;
			cout<<"Finish taking a look at the AncestorMap"<<endl;
		}

		if(fverbose) std::cout<<"SinglePhoton::CollectMCParticles_v2() \t||\t Finish!"<<endl;
	}


	/****************************
	 *
	 * BobbyVertexBuilder_ext() - Bobby's vertexbuilder! Find vertex from given tracks and showers.
	 *
	 * **************************/

	void SinglePhoton::BobbyVertexBuilder(class Atlas &package, bool more_objects){
		bool debug_message = true;
		bool fverbose = m_is_verbose;

		std::cout<<"SinglePhoton::BobbyVertexBuilder() \t||\t";
		verbose_check( fverbose, debug_message);

		//		initscr();//initialize the COLS and LINES for the screen size
		//		int screen_width = 86;//ncurses::COLS;//get the width of the screen!

		//Initialize criteria for vertexing,  unit: cm.
		//{4,30,100,25} from bobby's thesis
		double start_prox = 0.5;//4;//8;//Max. track proximity threshold (t_max)
		double shower_prox = 5;//30;//12;//Max. shower bkw point proximity threshold (s_max)
		double cpoa_vert_prox = 55;//7;//Max. distance btw shower start & cloest approach (dp_max), when no vertex, look at shower bkw projection and track distance
		double cpoa_trackend_prox = 5;//11;//Max. distance between shower bkw projection to a track end point; 
		double nearby_vertex_dist = 80;
		//		for (auto it:package.trackToDistMap){
		//			if(it.second > temp_max_value) temp_max_value = it.second;
		//		}
		//		for (auto it:package.showerToDistMap){
		//			if(it.second > temp_max_value) temp_max_value = it.second;
		//		}//ok, now find the furthest distance;


		//Declear two classes to kick off vertexing
		VertexBuilder vbuilder;//it was named vb
		ParticleAssociations_all candidates;// it was named pas

		//Initialize creteria for reconstruction.
		vbuilder.SetParameters({start_prox, shower_prox, cpoa_vert_prox, cpoa_trackend_prox});

		//		if(fvbuildert.ftree) vbuilder.SetVBT(&fvbuildert);
		//prepare parameters for pre-check # of tracks and showers;
		std::vector< art::Ptr<recob::Track> > use_tracks(package.selected_tracks);//CHECK all tracks
		std::vector< art::Ptr<recob::Shower> > use_showers(package.selected_showers);


		auto trackmap = package.trackToDistMap;
		auto showermap = package.showerToDistMap;
		//-- Get a preview of # of tracks and showers first	
		std::vector<int>	tem_bobbyshowersv;
		std::vector<int> 	tem_bobbytracksv;
//		double temp_max_value = 200;//objects that is too far away from the Pandora Reco...
		bool successful_vertex = false;//true, do them all together; false, do them by step;
		bool one_for_all = false;//true, means load more_tracks/showers all at once; sucessful_vertex has to be true to make it work;
		//		std::vector< int > target_topo = {1,1};//{1,2}<->1s2t signal;
		bool reset_bobbyvertex = true;
		//		int loop_tracker = 1;
		//		int max_loops = 4;//consider all objects within temp_max_value distance after the 4th loop.

		//update objects?
		if(more_objects && one_for_all) {
			for (auto object : package.more_tracks){
				use_tracks.push_back(object);
			}
			for (auto object : package.more_showers){
				use_showers.push_back(object);
			}
		}
		//sort out values of both maps first
		vector<double> all_dist;
		for (auto it:trackmap){
			if(it.second > nearby_vertex_dist){
				trackmap.erase(it.first);
			} else{
				all_dist.push_back(it.second);
			}
		}
		for (auto it:showermap){
			if(it.second > nearby_vertex_dist){
				showermap.erase(it.first);
			} else{
				all_dist.push_back(it.second);
			}
		}//ok, now find the furthest distance;
		std::sort(all_dist.begin(), all_dist.end());// distance from low to high;

		if(debug_message){
			cout<<"\nPreview Bobby Reco. Preformance"<<endl;
			cout<<setw(9)<<right<<"Attempts ";
			cout<<setw(14)<<right<<"|Input shower ";
			cout<<setw(14)<<left<< "|Input track";
			cout<<setw(14)<<right<<"|Reco. showers ";
			cout<<setw(14)<<left<< "|Reco. of tracks "<<endl;
		}
		while(!successful_vertex&&all_dist.size()>0){
			ParticleAssociations_all candidates_copy(candidates);
			candidates_copy.SetVerbose(false);
			vbuilder.SetVerbose(false);
			//CHECK control what tracks and showers to look at
			if(more_objects){
				//nearby_vertex_dist = all_dist[(all_dist.size()-1)*loop_tracker/max_loops] + 0.01;
//				nearby_vertex_dist = 80;//all_dist[(all_dist.size()-1)*loop_tracker/max_loops] + 0.01;
				cout<<"within the distance: "<<nearby_vertex_dist<<endl;
				//				cout<<"Search for additional tracks"<<endl;
				for (auto it:trackmap){
					if(it.second < nearby_vertex_dist){
						use_tracks.push_back(it.first);
						trackmap.erase(it.first);//eliminate for next round of consideration
					}
				}
				//				cout<<"Search for additional showers"<<endl;
				for (auto it:showermap){
					if(it.second < nearby_vertex_dist){
						use_showers.push_back(it.first);
						showermap.erase(it.first);
					}
				}
			} 
			break;
			//Disable the following, because we dont need expand range of slices step by step;
			/*
			   candidates_copy.GetDetectorObjects().AddShowers(use_showers);//load tracks
			   candidates_copy.GetDetectorObjects().AddTracks( use_tracks);//load showers
			   vbuilder.Run(candidates_copy);//here deals with the candidates_copy and find the vertex.

			//Coming update:
			//1. Label slice with Pandora reco. position
			//2. Add 4 slices at m.
			//3. Improve algo. later;


			//Fill in Bobby varaibles;
			//			std::cout<<"Preview Topological Reco. result status: "<<candidates_copy.GetSelectedAssociations().size()<<" vertex candidates_copy."<<std::endl;


			if(candidates_copy.GetSelectedAssociations().size() == 0){
			cout<<"No vertex is reconstructed."<<endl;
			}


			for(size_t const nth_associations : candidates_copy.GetSelectedAssociations()) {//Loop over all associations, which is a vector
			ParticleAssociation const & particle_associated = candidates_copy.GetAssociations().at(nth_associations);//grab the "pn"th association;

			//calculate the # of tracks/showers;
			DetectorObjects_all const & detos = candidates_copy.GetDetectorObjects();
			int temp_num_tracks = 0;
			int temp_num_showers = 0;

			for(size_t const n : particle_associated.GetObjectIndices()) {
			if(detos.GetRecoType(n) == detos.ftrack_reco_type) {
			++temp_num_tracks;
			}
			if(detos.GetRecoType(n) == detos.fshower_reco_type) {
			++temp_num_showers;
			}
			}
			if(( temp_num_showers >= target_topo[0]  && temp_num_tracks >= target_topo[1])||all_dist.size()<1) successful_vertex = true;

			if(fverbose){//Over view of the inputs
			cout<<setw(9)<<left<<loop_tracker<<"|";
			cout<<setw(13)<<use_showers.size()<<"|";
			cout<<setw(13)<<use_tracks.size()<<"|";
			cout<<setw(13)<<left<<temp_num_showers<<"|";
			cout<<setw(13)<<left<<temp_num_tracks<<endl;
			}
			}
			if(fverbose){
			for(int i = 0; i<screen_width; i++) 
			cout<<"-";
			cout<<endl;
			}

			loop_tracker++;
			if(loop_tracker > max_loops){
			break;
			}
			*/		
		}//end loop of testing diff. inputs for vertexing.

		candidates.SetVerbose(debug_message);
		vbuilder.SetVerbose(debug_message);
		//Looks good, then proceed to really fill in trees;

		//			cout<<endl;//printout 
		//			for(int i = 0 ; i < screen_width/2 ;i++)
		//				cout<<"/*";
		//			cout<<"\nBobby Revertexing is finished. Now start to fill in the TTree."<<endl;
		////			for(int i = 0 ; i < screen_width/2 ;i++)
		//				cout<<"/*";
		//			cout<<endl;
		//			cout<<" Reports from the ParticleAssociations_all"<<endl;	

		//			if(fverbose)candidates.PrintAssociations_all();
		//
		//			cout<<endl;//printout 
		//			for(int i = 0 ; i < screen_width/2 ;i++)
		//				cout<<"/*";
		//			cout<<"\n"<<endl;
		//			if(fverbose)
		//				candidates.PrintNodes();
		//			cout<<endl;//printout 
		//			for(int i = 0 ; i < screen_width/2 ;i++)
		//				cout<<"/*";
		//			cout<<"\n"<<endl;

		if(fverbose){
			std::cout<<"SinglePhoton::BobbyVertexBuilder() \t||\t";
			std::cout<<"Loading VB with "<<use_showers.size()<< " showers, ";
			std::cout<<use_tracks.size()<< " tracks."<<std::endl;
		}
		candidates.GetDetectorObjects().AddShowers(use_showers);//load showers
		candidates.GetDetectorObjects().AddTracks( use_tracks);//load tracks
		vbuilder.Run(candidates);//here deals with the candidates_copy and find the vertex.

		m_dist_tt = vbuilder.f_dist_tt;
		m_dist_sx = vbuilder.f_dist_sx;
		m_dist_st = vbuilder.f_dist_st;
		m_dist_sst = vbuilder.f_dist_sst;


		//temporary varaibles
		std::vector<double> tem_bobbyvertex_pos_xv;
		std::vector<double> tem_bobbyvertex_pos_yv;
		std::vector<double> tem_bobbyvertex_pos_zv;
		std::vector<int> 	tem_bobbyphotonshowerv;
		std::vector<int> 	tem_bobbypi0daughterv;
		std::vector<int> 	tem_bobbydeltaradppdaughterv;
		std::vector<int> 	tem_bobbydeltaradmdaughterv;
		std::vector<int> 	tem_bobbydeltaradpdaughterv;
		std::vector<int> 	tem_bobbydeltarad0daughterv;
		std::vector<int> 	tem_bobbyotherdaughterv;
		std::vector<int> 	tem_bobbyoverlayv;
		std::vector<int> 	tem_bobbyprotontrackv;
		std::vector<bool>	tem_same_slicev;

		std::vector<std::vector <int>> tem_track_daughter_pdgv;
		std::vector<std::vector <int>> tem_shower_daughter_pdgv;
		double min_bobbyvertexradius = 999;
		size_t min_index = 0;
		size_t temp_counter = 0;

		if(debug_message){
			cout<<"-- Check the reco-MC track Map (x-position - TrackId):"<<endl;
			for(auto const & [a,b]:package.trackToMCParticleMap){
				cout<<a->Start().X()<<" "<<b->TrackId()<<endl;
			}

			cout<<"-- Check the reco-MC shower Map (Shower Id - TrackId):"<<endl;
			for(auto const & [a,b]:package.showerToMCParticleMap){
				cout<<a->ShowerStart().X()<<" "<<b->TrackId()<<endl;
			}
			cout<<"-----------"<<endl;
		}

		for(size_t const nth_associations : candidates.GetSelectedAssociations()) {//Loop over all associations, which is a vector
			ParticleAssociation const & particle_associated = candidates.GetAssociations().at(nth_associations);//grab the "pn"th association;
			geoalgo::Point_t const & reco_vertex = particle_associated.GetRecoVertex();//Grab the vertec of the "pn"th association.
			if(reset_bobbyvertex){
				m_bobbysameslicev.clear();
				m_bobbyvertex_pos_xv.clear();
				m_bobbyvertex_pos_yv.clear();
				m_bobbyvertex_pos_zv.clear();
				m_bobbytracksv.clear();
				m_bobbyshowersv.clear();

				m_bobbyprotontrackv.clear();
				m_bobbyphotonshowerv.clear();
				m_bobbypi0daughterv.clear();
				m_bobbydeltaradppdaughterv.clear();
				m_bobbydeltaradmdaughterv.clear();
				m_bobbydeltaradpdaughterv.clear();
				m_bobbydeltarad0daughterv.clear();
				m_bobbyotherdaughterv.clear();
				m_bobbyoverlayv.clear();
				m_bobbyvertexradiusv.clear();
				m_bobbytrackdaughter_pdg.clear();
				m_bobbyshowerdaughter_pdg.clear();
				reset_bobbyvertex = false;
			}

			tem_bobbyvertex_pos_xv.push_back( reco_vertex.at(0));
			tem_bobbyvertex_pos_yv.push_back( reco_vertex.at(1));
			tem_bobbyvertex_pos_zv.push_back( reco_vertex.at(2));

			if(fverbose){
				std::cout<<"SinglePhoton::BobbyVertexBuilder() \t||\t";
				cout<<"BB Vertex: (";
				cout<<reco_vertex.at(0)<<", ";
				cout<<reco_vertex.at(1)<<", ";
				cout<<reco_vertex.at(2)<<") "<<endl;
			}

			//calculate the # of tracks/showers;
			DetectorObjects_all const & detos = candidates.GetDetectorObjects();
			int temp_num_tracks = 0;
			int temp_num_showers = 0;

			int get_a_proton = 0;
			int get_a_photon = 0;
			int get_a_pi0daughter = 0;
			int get_a_deltaradppdaughter = 0;
			int get_a_deltaradmdaughter = 0;
			int get_a_deltaradpdaughter = 0;
			int get_a_deltarad0daughter = 0;
			int get_a_otherdaughter = 0;
			int get_a_overlay = 0;
			bool same_slice = true;
			bool first_time = true;
			int temp_sliceid = 0;


			std::vector <int> get_a_track_daughter_pdg;
			std::vector <int> get_a_shower_daughter_pdg;
			bool find_MC = true;

			for(size_t const n : particle_associated.GetObjectIndices()) {
				int index;
				art::Ptr<simb::MCParticle> temp_mcp;

				m_bobbyvertexradius = particle_associated.GetGoodness();
				//mark down the smallest radius;
				if(m_bobbyvertexradius < min_bobbyvertexradius && m_bobbyvertexradius > 0){
					min_bobbyvertexradius = m_bobbyvertexradius;
					min_index = temp_counter;
				}

				if(detos.GetRecoType(n) == detos.ftrack_reco_type) {//it is a track
					++temp_num_tracks;
					index = detos.GetTrackIndexFromObjectIndex(n);
					if(debug_message)cout<<"A track "<< n <<" feed-in order "<<index;
					//---same slice check---
					art::Ptr<recob::PFParticle> tempp = package.trackToNuPFParticleMap.find(use_tracks[index])->second;	
					for(auto const &[a,b]:package.sliceIDToPFParticlesMap){
						if(std::find(b.begin(),b.end(),tempp)!=b.end() ) {
							if(debug_message) cout<<"SLICE ID "<<a<<endl;

							if(first_time){
								temp_sliceid = a;
								first_time = false;
							} else if(temp_sliceid != a){
								same_slice = false;
							}
							break;
						}
					}
					//---------------------

					if(package.trackToMCParticleMap.find(use_tracks[index])==package.trackToMCParticleMap.end()){
						//Overlays!
						get_a_track_daughter_pdg.push_back(-999);
						if(debug_message) cout<<" No match to any MCParticles."<<endl;
						get_a_overlay++;
						find_MC = false;
						continue;
					}else{ 
						temp_mcp = package.trackToMCParticleMap.find(use_tracks[index])->second;//CHECK THIS
						get_a_track_daughter_pdg.push_back(temp_mcp->PdgCode());
					}
				}

				if(detos.GetRecoType(n) == detos.fshower_reco_type) {//it is a shower

					++temp_num_showers;
					index = detos.GetShowerIndexFromObjectIndex(n);
					if(debug_message) cout<<"A shower index "<< n <<" feed-in order "<<index;
					//--------- same slice check --------
					art::Ptr<recob::PFParticle> tempp = package.showerToNuPFParticleMap.find(use_showers[index])->second;	
					for(auto const &[a,b]:package.sliceIDToPFParticlesMap){
						if(std::find(b.begin(),b.end(),tempp)!=b.end() ) {
							//						cout<<"SLICE ID "<<a<<endl;
							if(first_time){
								temp_sliceid = a;
								first_time = false;
							} else if(temp_sliceid != a){
								same_slice = false;
							}
							break;
						}
					}
					// ---------------------------

					if(package.showerToMCParticleMap.find(use_showers[index])==package.showerToMCParticleMap.end()){
						//Overlays!
						get_a_shower_daughter_pdg.push_back(-999);
						if(debug_message) cout<<" No match to any MCParticles."<<endl;
						get_a_overlay++;
						find_MC = false;
						continue;
					}else { 
						temp_mcp = package.showerToMCParticleMap.find(use_showers[index])->second;
						get_a_shower_daughter_pdg.push_back(temp_mcp->PdgCode());
					}
				}
				//identify shower/track MCTruth info.
				if(debug_message) cout<<" TrackId "<<temp_mcp->TrackId()<<" PdgCode "<<temp_mcp->PdgCode();
				switch (temp_mcp->PdgCode()){
					case 2212:
						get_a_proton++;
						break;
					case 22:
						get_a_photon++;
						break;
					default:
						break;
				}
				//identify the ancestor!
//				art::Ptr<simb::MCParticle> temp_mcpmother= package.MCParticleToAncestorMap.find(temp_mcp)->second;

//				if(debug_message) cout<<" Ancestor PdgCode "<<package.AncestorToPdgMap.find(temp_mcpmother)->second;
				switch( package.MCParticleToAncestorPdgMap.find(temp_mcp)->second){
//				switch(111){
				//second+ generation MCParticles lose connection to the ancestor, because I didnt trace all the way up in the AncestorToPdgMap;
					case 2224://delta++
						if(debug_message) cout<<" from delta++"<<endl;
						get_a_deltaradppdaughter++;
						break;
					case 1114://delta- (not likely)
						if(debug_message) cout<<" from delta-"<<endl;
						get_a_deltaradmdaughter++;
						break;
					case 2214://delta+
						if(debug_message) cout<<" from delta+"<<endl;
						get_a_deltaradpdaughter++;
						break;
					case 2114://delta0
						if(debug_message) cout<<" from delta0"<<endl;
						get_a_deltarad0daughter++;
						break;
					case 111:
						if(debug_message) cout<<" from pi0"<<endl;
						get_a_pi0daughter++;
						break;
					default:
						get_a_otherdaughter++;
						if(debug_message) cout<<" from other"<<endl;
						break;
				}

				temp_counter++;
			}
			if(find_MC&debug_message) cout<<" All MC, no overlays!"<<endl;
			if(fverbose){
				std::cout<<"SinglePhoton::BobbyVertexBuilder() \t||\t";
				cout<<"Summary: reconstruct a vertex with ";
				cout<<temp_num_showers<<" showers, ";
				cout<<temp_num_tracks<<" tracks. Vertex radius ";
				cout<<min_bobbyvertexradius;
			}

			if(debug_message){
				if(same_slice){
					cout<<" In the same Slice"<<endl;
				}else{
					cout<<" Not in the same Slice"<<endl;
				}
			}
			tem_same_slicev.push_back(same_slice);
			tem_bobbytracksv.push_back(temp_num_tracks);
			tem_bobbyshowersv.push_back(temp_num_showers);

			tem_bobbyprotontrackv.push_back(get_a_proton);
			tem_bobbyphotonshowerv.push_back(get_a_photon);
			tem_bobbypi0daughterv.push_back(get_a_pi0daughter);
			tem_bobbydeltaradppdaughterv.push_back(get_a_deltaradppdaughter);
			tem_bobbydeltaradmdaughterv.push_back(get_a_deltaradmdaughter);
			tem_bobbydeltaradpdaughterv.push_back(get_a_deltaradpdaughter);
			tem_bobbydeltarad0daughterv.push_back(get_a_deltarad0daughter);
			tem_bobbyotherdaughterv.push_back(get_a_otherdaughter);
			tem_bobbyoverlayv.push_back(get_a_overlay);
			//	CHECK, need to figure out how to push vector into a vector<vector>
			//			tem_track_daughter_pdgv.push_back(get_a_track_daughter_pdg);
			tem_shower_daughter_pdgv.push_back(get_a_shower_daughter_pdg);

			m_bobbyvertexradiusv.push_back(m_bobbyvertexradius);
		}


		//Finalize Bobby info. and pick a vertex set.
		m_bobbysameslicev = tem_same_slicev;
		m_bobbyvertex_pos_xv = tem_bobbyvertex_pos_xv;
		m_bobbyvertex_pos_yv = tem_bobbyvertex_pos_yv;
		m_bobbyvertex_pos_zv = tem_bobbyvertex_pos_zv;
		m_bobbyprotontrackv  = tem_bobbyprotontrackv;
		m_bobbyphotonshowerv = tem_bobbyphotonshowerv;
		m_bobbypi0daughterv  = tem_bobbypi0daughterv;
		m_bobbydeltaradppdaughterv  = tem_bobbydeltaradppdaughterv;
		m_bobbydeltaradmdaughterv  = tem_bobbydeltaradmdaughterv;
		m_bobbydeltaradpdaughterv  = tem_bobbydeltaradpdaughterv;
		m_bobbydeltarad0daughterv  = tem_bobbydeltarad0daughterv;
		m_bobbyotherdaughterv  = tem_bobbyotherdaughterv;
		m_bobbyoverlayv  = tem_bobbyoverlayv;
		m_bobbytracksv = tem_bobbytracksv;
		m_bobbyshowersv= tem_bobbyshowersv;

		//		for( size_t index = 0; index< m_bobbyvertex_pos_xv.size() ; index++){//loop over bobbyvertex
		//metric for vertex:
		//Topology:
		//VertexSize/doogness:

		//				if(temp_dist < best_vertex_dist){}//update when find a closer vertex to the pandora vertex.
		if(true){//maybe pick the smallest vertex radius (for now),
			//	best_vertex_dist = temp_dist;
			m_bobbyvertex_pos_x = m_bobbyvertex_pos_xv[min_index];
			m_bobbyvertex_pos_y = m_bobbyvertex_pos_yv[min_index];
			m_bobbyvertex_pos_z = m_bobbyvertex_pos_zv[min_index];
			m_bobbytracks =  m_bobbytracksv[min_index];
			m_bobbyshowers = m_bobbyshowersv[min_index];
			m_bobbyphotonshower = m_bobbyphotonshowerv[min_index];
			m_bobbypi0daughter = m_bobbypi0daughterv[min_index];
			m_bobbydeltaradppdaughter = m_bobbydeltaradppdaughterv[min_index];
			m_bobbydeltaradmdaughter = m_bobbydeltaradmdaughterv[min_index];
			m_bobbydeltaradpdaughter = m_bobbydeltaradpdaughterv[min_index];
			m_bobbydeltarad0daughter = m_bobbydeltarad0daughterv[min_index];
			m_bobbyotherdaughter = m_bobbyotherdaughterv[min_index];
			m_bobbyoverlay = m_bobbyoverlayv[min_index];
			m_bobbyprotontrack = m_bobbyprotontrackv[min_index];
			//				m_bobbytrackdaughter_pdg = tem_track_daughter_pdgv[min_index];
			//				m_bobbyshowerdaughter_pdg = tem_shower_daughter_pdgv[min_index];
			//CHECK add shower/track daughter here;
		}
		//		}

		//if( m_run_all_pfps && m_bobbyvertexing_more &&m_bobbytracks + m_bobbyshowers < 2&& run_count < 2 ){//repeat with 0,1
		//	run_count++;
		//	cout<<"Need more objects. Now run the "<<run_count<<" time the vertexing."<<endl;
		//	goto redo_event;
		//}
		if(fverbose){
			std::cout<<"SinglePhoton::BobbyVertexBuilder() \t||\t";
			std::cout<<"Finish!\n"<<std::endl;
		}

	}
}
#endif
