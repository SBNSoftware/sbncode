#ifndef VERTEXBUILDER_H
#define VERTEXBUILDER_H

#include "ParticleAssociations.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TTree.h"
#include <fstream>//read and write txt file.

/**********************
 *stuct VertexBuilderTree: it store values of variables extracted 
	from the artObjects? (values from the raw root file)
 *
 *class VertexBuilder: it contains functions for determining the vertex of objects
	assigned by the class ParticleAssociations_all.
 *
 *********************/

using namespace std;

namespace single_photon{
	struct VertexBuilderTree {

		TTree * ftree;
		int frun_number;
		int fsubrun_number;
		int fevent_number;
		int ftrack_number;
		int fshower_number;
		//	int fassociation_track_number; //CHECK,number of track associations;
		//	int fassociation_shower_number; //CHECK,number of shower associations;
		int fassociation_final_number; //number of total associations;

		VertexBuilderTree() :
			ftree(nullptr){}
		/*********************************************************************
		 *	Setup() - Obtain the addresses of event variables 
		 *		(indices for particular events, showers, and tracks).
		 *********************************************************************/
		void Setup() {
			art::ServiceHandle< art::TFileService > tfs;
			ftree = tfs->make<TTree>("VertexBuilder", "");
			ftree->Branch("run_number", &frun_number, "run_number/I");
			ftree->Branch("subrun_number", &fsubrun_number, "subrun_number/I");
			ftree->Branch("event_number", &fevent_number, "event_number/I");  
			ftree->Branch("track_number", &ftrack_number, "track_number/I");
			ftree->Branch("shower_number", &fshower_number, "shower_number/I");
			//		ftree->Branch("association_track_number", &fassociation_track_number, "association_track_number/I");
			//		ftree->Branch("association_shower_number", &fassociation_shower_number, "association_shower_number/I");
			ftree->Branch("association_final_number", &fassociation_final_number, "association_final_number/I");
		}

	};


	class VertexBuilder {

		geoalgo::GeoAlgo const falgo;

		size_t fobject_id;

		double fstart_prox;//the maximum track proximity threshold (in cm)
		double fshower_prox;//the maximum shower proximity threshold (in cm)
		double fcpoa_vert_prox;//the maximum distance btw shower start to a track end point(in cm); 
		double fcpoa_trackend_prox;//the maximum distance btw mid point of impact parameter of shower to the vertex (in cm);

		DetectorObjects_all const * fdetos;

		VertexBuilderTree * fvbt;

		bool fverbose;

		void CheckSetVariables();// exit if the criteria are not set.

		void Erase(std::multimap<size_t, geoalgo::Point_t const *> & pn,
				std::multimap<size_t, geoalgo::Point_t const *>::iterator const best_it,
				geoalgo::Point_t const & sv);

		/*****************************
		 *
		 * AssociateTracks() - this associates tracks to either end points of other tracks
		 *	whose are nearby.
		 *
		 * ***************************/
		void AssociateTracks(ParticleAssociations_all & pas);

		/*
		 * FindClosestApproach() - find the center points of two points along 
		 *		the backward projection that give minimum impact parameter; this
		 *		also returns the distance of two of such points.
		 *shr1, shr2 - two showers with start points and direction
		 *vtx - a vertex to be updated.
		 */
		double FindClosestApproach(const geoalgo::HalfLine_t & shr1,
				const geoalgo::HalfLine_t & shr2,
				geoalgo::Point_t & vtx) const;

		/*****************************
		 *
		 * AssociateShowers() - Run after AssociateTracks(); 
		 *	this associates showers to either end points of a track
		 *	that is considered to be from the same vertex.
		 *
		 * ***************************/
		void AssociateShowers(ParticleAssociations_all & pas);


		void AddLoneTracks(ParticleAssociations_all & pas);
		void AddLoneShowers(ParticleAssociations_all & pas);


		/*****************************
		 *
		 * AssociateShowers() - Run after AssociateTracks(); 
		 *	this associates showers to either end points of a track
		 *	that is considered to be from the same vertex.
		 *
		 * ***************************/
		void FillVBT(ParticleAssociations_all & pas);

		void dist_marker ( vector <double>& dist_field, double input) {
			if(dist_field[0] == 999) {
				dist_field[0] = input;
			}else {
				dist_field.push_back(input);
			}
			cout<<"Mark a distance."<<endl;
		}

		public:
		vector <double> f_dist_tt;//track&track, start_prox
		vector <double> f_dist_sx;//shower&anything, shower_prox
		vector <double> f_dist_st;//shower&track, cpoa_vert_prox
		vector <double> f_dist_sst;//shower&shower&track, cpoa_trackend_prox


		VertexBuilder();

		void SetVerbose(bool const verbose = true) {//allow output info. to the terminal
			fverbose = verbose;
		}

		void SetParameters (std::vector< double > p){

			fstart_prox			= p[0];
			fshower_prox		= p[1];
			fcpoa_vert_prox		= p[2];
			fcpoa_trackend_prox = p[3];

			f_dist_tt.push_back(999);//track&track, start_prox
			f_dist_sx.push_back(999);//shower&anything, shower_prox
			f_dist_st.push_back(999);//shower&track, cpoa_vert_prox
			f_dist_sst.push_back(999);//shower&shower&track, cpoa_trackend_prox

			if(fverbose){
				cout<<right<<setw(82)<<" Max. track proximity threshold (t_max)= "<<fstart_prox<<endl;
				cout<<right<<setw(82)<<" Max. shower proximity threshold (s_max)= "<<fshower_prox<<endl;
				cout<<right<<setw(82)<<" Max. distance btw shower start & cloest approach (dp_max)= "<<fcpoa_vert_prox<<endl;
				cout<<right<<setw(82)<<" Max. distance btw midway point of impact parameter to a precandidate vertex (a_max)= "<<fcpoa_trackend_prox<<endl;
			}
		}


		//Output Gadget
		void print_gadget ( const std::string& input, const std::string& symbol){
			double screen_width = 86;
			for(int i = 0; i<screen_width; i++) 
				cout<<symbol;
			cout<<"\n"<<input<<endl;
			for(int i = 0; i<screen_width; i++) 
				cout<<symbol;
			cout<<endl;
		}

		//	void SetVBT(VertexBuilderTree * vbt) {fvbt = vbt;}

		//  void AddTracks(art::ValidHandle<std::vector<recob::Track>> const & ev_t);
		//  void AddShowers(art::ValidHandle<std::vector<recob::Shower>> const & ev_s);
		//	void AddTracks(std::vector<art::Ptr<recob::Track>> const & ev_t);
		//	void AddShowers(std::vector<art::Ptr<recob::Shower>> const & ev_t);

		// The main code is here; it calls up all the necessary process.
		void Run(ParticleAssociations_all & pas);

	};















	//HEADER FILE ARE ABOVE














	VertexBuilder::VertexBuilder() :
		fobject_id(0),
		fstart_prox(-1),
		fshower_prox(-1),
		fcpoa_vert_prox(-1),
		fcpoa_trackend_prox(-1),
		fdetos(nullptr),
		fvbt(nullptr),
		fverbose(true) {}


	void VertexBuilder::CheckSetVariables() { //this was initialized 

		if(fstart_prox == -1) {//It start_prox was not previously defined, exit.
			std::cout << "fstart_prox, the maximum track proximity threshold, is not set\n";
			exit(1);
		}

		if(fshower_prox == -1) {
			std::cout << "fshower_prox, the maximum shower proximity threshold, is  not set\n";
			exit(1);
		}

		if(fcpoa_vert_prox == -1) {
			std::cout << "fcpoa_vert_prox, the maximum distance btw shower start to the cloest approach, is not set\n";
			exit(1);
		}

		if(fcpoa_trackend_prox == -1) {
			std::cout << "fcpoa_trackend_prox, the maximum distance btw mid point of impact parameter of shower to the vertex, is not set\n";
			exit(1);
		}

	}


	void VertexBuilder::Erase(std::multimap<size_t, geoalgo::Point_t const *> & pn,
			std::multimap<size_t, geoalgo::Point_t const *>::iterator const best_it,
			geoalgo::Point_t const & sv) {//Remove candidates that have been analyzed.

		size_t const index = best_it->first;
		pn.erase(best_it);

		if(fdetos->GetRecoType(index) == fdetos->ftrack_reco_type) {
			auto const pn_it = pn.find(index);
			if(fdetos->GetTrack(index).ftrajectory.Length() < fstart_prox ||
					(pn_it != pn.end() && pn_it->second->Dist(sv) < fstart_prox)) {
				pn.erase(pn_it);
			}
		}

	}


	void VertexBuilder::AssociateTracks(ParticleAssociations_all & pas) {//Group tracks

		//Ingredients: Key for tracks, Point coordinates along the track (recob:: Trajectory)
		int screen_width = 86;

		std::multimap<size_t, geoalgo::Point_t const *> pn;//This maps an # to a pointer to the track start/end points
		/* Structure of pn:
		 *	Iterator Key	Value	
		 *	0		Key[0]	point1
		 *	1		Key[0]	point2
		 *	2		Key[1]	pointa
		 *	3		Key[1]	pointb
		 *	.
		 *	example of calling first row:
		 *	pn.begin(), pn[0].first, pn[0].second
		 */

		//#1 load up tracks to pn: 1 id for 2 end points;
		for(size_t const i : fdetos->GetTrackIndices()) {

			Track const & t = fdetos->GetTrack(i);
			geoalgo::Vector ref_point(1,1,1);//use this to pick valid end points; i.e. kill (-999,-999,-999)

			if(t.ftrajectory.front().Dot(ref_point) > -2996){
				pn.emplace(t.fid, &t.ftrajectory.front());//track ID and point;
			}
			if(t.ftrajectory.back().Dot(ref_point) > -2996){
				pn.emplace(t.fid, &t.ftrajectory.back());
			}
			/* appended map of one of this loop:
			 * < track_id, {track beginning point, track ending point}>
			 */
		}

		if(fverbose) {
			std::cout << "STEP I: Load " << pn.size() << " track end points to be evaluated.\n\n";

			if(pn.size()) {//Print the title
				std::stringstream output_buf;
				output_buf<<setw(9)<<right<<"Track ID ";
				output_buf<<setw(15)<<left<<"| Track Length ";
				output_buf<<setw(37)<<left<<"| Track End Point Coordinates (x,y,z) [cm] from Pandora"<<endl;
				print_gadget ( output_buf.str(), "-");

				for(auto p : pn){
					cout<<setw(9)<<left<< p.first<<"| ";//Index
					cout<<setw(13)<<left<<fdetos->GetTrack(p.first).ftrajectory.Length();//Track Length
					//ftrajectory is the obejct in art library
					cout<<"| "<<*p.second<<endl;//Coordinates
				}
				for(int i = 0; i<screen_width; i++) 
					cout<<"-";
				cout<<"\n\nSTEP II: Look for a candidate vertex among tracks:"<<endl;
			}
		}

		//#2 evaluate all loaded end points of tracks; 
		while(pn.size() > 1) {//while there are tracks to be evaluated, loop them and try to make candidate vertex.

			auto best_match = pn.end(); //best match (iterator) is assumed to be the iterator of the last track;
			auto best_candidate = best_match; //best candidate;

			double best_dist = fstart_prox;//Max. track proximity threshold (t_max)
			if(fverbose)
				std::cout <<  "\tRemain "<< pn.size() <<" track end points for evaluation," << std::endl;

			for(auto match_this = pn.begin(); match_this != pn.end(); ++match_this) {//iterator to indicate the point we are evaluating

				if(std::next(match_this) == pn.end()) break;//stop when there is only one track end point left;
				if(fverbose)
					std::cout << "\t\tEnd point from the track (ID: " << match_this->first <<"):"<< std::endl;

				double temp_save_this_dist = 999;
				for(auto compare_this = std::next(match_this); compare_this != pn.end(); ++compare_this) {
					if(fverbose)
						std::cout << "\t\t\tCompare to a end point from a track (ID: " << compare_this->first<<"); ";
					if(match_this->first == compare_this->first){//if two end points from the same track, skip this round.
						if(fverbose) std::cout<<"\t\t\t\tSame track, skip!"<<endl;
						continue;
					}

					double const dist = compare_this->second->Dist(*match_this->second);//distance between match_this point to compare_this point;


					if(fverbose) std::cout << "distance btw them: " << dist << ", the shortest: "<< best_dist << std::endl;
					//				if(dist==0) std::cout<<" Reject this point."<<endl;
					if(dist < best_dist) {
						if(fverbose)
							std::cout << "\t\t\t\t>>Take this track end point as a better candidate vertex!\n\n";	
						best_match = match_this;
						best_candidate = compare_this;
						best_dist = dist;
					}
					if(temp_save_this_dist>dist) temp_save_this_dist = dist;
				}
				dist_marker(f_dist_tt, temp_save_this_dist);
			}

			if(best_match == pn.end() || pn.size()<2) {
				if(fverbose)
					std::cout << "\tNah, tracks are all single. No vertex candidate for tracks.\n";
				return;
			}
			//at this moment, we have:
			//	best_match - the iterator of an end-point of one track
			//	best_candidate - the iterator of an end-point of another track
			//	best_dist - the cloest distance between these two points

			std::vector<size_t> track_IDs;
			track_IDs.push_back(best_match->first);
			track_IDs.push_back(best_candidate->first);

			std::vector<geoalgo::Point_t> points;
			points.push_back(*best_match->second);
			points.push_back(*best_candidate->second);

			geoalgo::Sphere sphere(falgo.boundingSphere(points));//falgo is an GeoAlgo object, this sets the mid point of two track endpoints as new candidate verex.

			Erase(pn, best_match, sphere.Center());//remove the points that has been added as candidate vertex.
			Erase(pn, best_candidate, sphere.Center());

			//#3 add more tracks around the vertex candidate
			if(fverbose){
				std::cout << "\n\nSTEP III: Add tracks around the vertex candidate. ";
				std::cout << "ID: "<< best_match->first<<" at (";
				cout<<falgo.boundingSphere(points).Center().at(0)<<",";
				cout<<falgo.boundingSphere(points).Center().at(1)<<",";
				cout<<falgo.boundingSphere(points).Center().at(2)<<")"<<endl;
			}

			do {//have to break to leave.
				auto best_o = pn.end();//best object to add;
				double sbest_dist = fstart_prox;//Max. track proximity threshold (t_max)

				geoalgo::Point_t const & sphere_center = sphere.Center();//this sets the candidate vertex;

				if(fverbose)
					std::cout << "\t"<< pn.size()<<" track end points to be considered: "<<endl;

				for(auto this_iterator = pn.begin(); this_iterator != pn.end(); ++this_iterator) {

					double const dist = this_iterator->second->Dist(sphere_center);
					if(fverbose){
						std::cout << "\t\tTrack (ID: " << this_iterator->first<< ") end point ";
						std::cout << "is at distant: " << dist << " maximum allowed distance: " << sbest_dist << std::endl;
					}

					if(std::find(track_IDs.begin(), track_IDs.end(), this_iterator->first) != track_IDs.end()) {
						if(fverbose)
							std::cout << " This is not the last track?(CHECK) Skip\n";
						continue;
					}

					if(dist < sbest_dist) {//when find a close track end point, attch it to the candidate vertex.
						if(fverbose)
							std::cout << "\t\t\tTake this track. Update the shortest distance!\n";
						best_o = this_iterator;
						sbest_dist = dist;
					}
				}

				if(best_o == pn.end()) {//exit the loop only when 
					if(fverbose)
						std::cout << "\tNo more tracks can be added.\n";
					break;
				}

				track_IDs.push_back(best_o->first);
				points.push_back(*best_o->second);

				Erase(pn, best_o, *best_o->second);

				//s = algo.boundingSphere(points);

			} while(true);//if false, then it only loops once;
			//points - point objects that contains positions;
			//boundingSphere(points) is a minimum sphere that contains points.
			pas.AddAssociation(track_IDs, points, sphere.Center(), falgo.boundingSphere(points).Radius());

			if(fverbose){//no worry, the code would not run this part is there were no candidate vertex.
				cout<<"\nSummary: add "<<track_IDs.size();
				cout<<" tracks (IDs: ";
				for(size_t i = 0; i < track_IDs.size(); i++){
					cout<<track_IDs[i]<<" ";
				}
				cout<<") to the candidate vertex: (";
				cout<<falgo.boundingSphere(points).Center().at(0)<<",";
				cout<<falgo.boundingSphere(points).Center().at(1)<<",";
				cout<<falgo.boundingSphere(points).Center().at(2)<<")"<<endl;
			}
			//move to next track end point, but actually, only one vertex candidate is considered here, so no more end point here;
			//Keng, current version takes only the best candidate vertex, maybe take all under consideration is better?
		}
		//  return;
	}

	double VertexBuilder::FindClosestApproach(
			const geoalgo::HalfLine_t & shr1,
			const geoalgo::HalfLine_t & shr2,
			geoalgo::Point_t & vtx) const {
		// Find mininum impact parameter between a two showers
		// flip their directions and look backwards

		// Create a half-line pointing backwards from the shower
		geoalgo::HalfLine_t shr1Bkwd(shr1.Start(), shr1.Dir() * (-1));
		geoalgo::HalfLine_t shr2Bkwd(shr2.Start(), shr2.Dir() * (-1));

		// coordinates for closest points on the two objects
		geoalgo::Point_t PtShr1(3);
		geoalgo::Point_t PtShr2(3);
		double IP = falgo.SqDist(shr1Bkwd, shr2Bkwd, PtShr1, PtShr2);//give the square of the cloesest distance in 3D
		//see https://nusoft.fnal.gov/larsoft/doxsvn/html/classgeoalgo_1_1GeoAlgo.html#a605d7e7a2736237727cad00c80c27eb8 for more;

		// build candidate vertex
		vtx = (PtShr1 + PtShr2) / 2.;

		return sqrt(IP);
	}


	void VertexBuilder::AssociateShowers(ParticleAssociations_all & pas) {

		//for creterias, cr1=10, cr2, cr3..
		//#1 load showers
		//#2 find precandidate
		//- distance shower_i to shower_j < cr1
		//- distance shower_i to track_k
		// ---> next shower_i
		//
		int screen_width = 86;

		std::map<size_t, Shower const *> shower_map;

		//#1 load showers
		for(size_t const i : fdetos->GetShowerIndices()) { 
			shower_map.emplace(i, &pas.GetDetectorObjects().GetShower(i));
		}

		std::vector<ParticleAssociation> const & associations = pas.GetAssociations();//Load vertex candidates, that contains vertices, radius of the vertex (bounding phere)

		while(shower_map.size()>0) {
			if(fverbose) {
				std::cout << "STEP I: Load " << shower_map.size() << " showers for evaluation.\n\n";

				if(shower_map.size()) {
					//Print the title
					std::stringstream output_buf;
					output_buf<<setw(10)<<right<<"Shower ID ";
					output_buf<<"| Shower Start Point Coordinates (x,y,z) [cm] from Pandora ";
					output_buf<<"| Shower Direction from Pandora"<<endl;
					print_gadget ( output_buf.str(), "-");

					for(auto const &c : shower_map){
						cout<<setw(10)<<left<< c.first<<"| ";//Index
						cout<<left<<c.second->fcone.Start();//Shower Coordinates
						cout<<" | "<<c.second->fcone.Dir()<<endl;//Direction
					}
					for(int i = 0; i<screen_width; i++) 
						cout<<"-";
					cout<<"\n\nSTEP II: Look for a precandidate vertex between a shower to other showers/tracks/candidate vertices."<<endl;
					cout<<"Note: 'bkw point' means the point along the shower backward projection that gives the shortest distance to other objects."<<endl;
				}
			}

			if(fverbose)	std::cout << "\tRemain "<< shower_map.size()<<" showers to be evaluated: "<<endl;
			size_t best_shower_id = SIZE_MAX;
			size_t best_other_id = SIZE_MAX;//the object that is cloest to the shower bkw prjection line;
			size_t index = SIZE_MAX;//updated when the vertex candidate is associated as the last action.
			bool AddShowerToVertex = false;//true means a shower is attached to a vertex candidate.

			Double_t best_dist = fshower_prox;//save the fshower_prox, distance form shower to shower/track/candidate vertex
			geoalgo::Point_t best_vert(2000, 2000, 2000);

			//#2 find precandidate vertex for any showers to 1) other showers, 2) tracks, and 3) candidate vertex.
			for(auto this_iterator = shower_map.begin(); this_iterator != shower_map.end(); ++ this_iterator) {//look at each shower candidate, 

				//#2.1 shower to shower
				if(fverbose){
					std::cout << "\t\tShower (ID: " << this_iterator->first <<")"<< std::endl;
					std::cout << "\t\t\t Compare to other showers:"<<std::endl;
				}

				geoalgo::Point_t const & c_start = this_iterator->second->fcone.Start();
				geoalgo::Vector_t const & c_dir = this_iterator->second->fcone.Dir();

				for(auto that_iterator = std::next(this_iterator); that_iterator != shower_map.end(); ++ that_iterator) {//look at each shower candidate, 

					if(fverbose) std::cout << "\t\t\t Shower (ID: " << that_iterator->first <<") ";

					if( that_iterator->first == this_iterator->first) {
						if(fverbose) std::cout << "is the same shower, skip (CHECK)\n";
						continue;
					}
					geoalgo::Point_t temp_vert;

					double dist = FindClosestApproach( that_iterator->second->fcone, this_iterator->second->fcone, temp_vert);//minimum impact parameter
					//this is the shortest distance between two backward projection in 3D.
					//update temp_vert as the mid point of these two points that gives the shortest distance.

					if(fverbose)
						std::cout << "has bkw-projection distant to the current shower's bkw-projection, " << dist << ", and the shortest-dist so far is "
							<< best_dist << ".\n";

					double temp_dist = best_dist;//fshower_prox now is temp_dist

					if(dist < temp_dist) {//min. impact parameter < max shower proximity threshold

						if(fverbose) std::cout << "\t\t\t\t>>Take their mid point as the candidate vertex and update the shortest-dist!\n";

						best_shower_id = this_iterator->first;
						best_other_id = that_iterator->first;
						best_vert = temp_vert;
						best_dist = dist;
						index = SIZE_MAX;
						AddShowerToVertex = false;
					}
				}

				//			if(fverbose) std::cout << "\t\t >>Finish comparing showers to the shower with ID: "<<this_iterator->first<<endl;
				if(fverbose){ 
					std::cout << "\t\t\t Finish comparison between shower and shower;\n"<<endl;
					std::cout << "\t\t\t Compare to tracks:"<<std::endl;
				}


				//#2.2 shower to track
				for(size_t const i : fdetos->GetTrackIndices()) {//compare tracks

					if(fverbose)
						std::cout << "\t\t\t Track (ID: " << i <<") ";
					Track const & t = fdetos->GetTrack(i);

					geoalgo::Point_t dont_care;
					geoalgo::Point_t temp_vert;

					Double_t dist =
						sqrt(falgo.SqDist(t.ftrajectory, 
									geoalgo::HalfLine(c_start,
										c_dir*-1),
									temp_vert,//this point will be changed along the track; id does not matter what value it is;
									dont_care));//this is along the backward projection
							//CHECK, temp_vert will becomes the first point of the ftrajectory, if parallel? this happens !

					if(fverbose)
						std::cout << "is distant to the bkw point, "
						<< dist << ", and the shortest_dist so far is " << best_dist << ".\n";

					if(dist < best_dist) {

						if(fverbose) std::cout << "\t\t\t\t>>Take the bkw point ("<<temp_vert.at(0)<<","<<temp_vert.at(1)<<","<<temp_vert.at(2)
									<<") as the precandidate vertex and update the shortest-dist!\n";
						best_shower_id = this_iterator->first;
						best_other_id = i;
						best_vert = temp_vert;
						best_dist = dist;
						index  = SIZE_MAX;
						AddShowerToVertex = false;
					}
				}


				//if(fverbose) std::cout << "\t\t >>Finish comparing tracks to the shower with ID: "<<this_iterator->first<<endl;
				if(fverbose) {
					std::cout << "\t\t\t Finish comparison between shower and track;\n"<<endl;
					std::cout << "\t\t\t Compare to candidate vertex:"<<std::endl;
				}

				//#2.3 shower to vertex candidate
				for(size_t i = 0; i < associations.size(); ++i) {//compare vertex candidates

					if(fverbose)
						std::cout << "\t\t\t Vertex candidate vertex (ID: "<< i << " ) ";

					ParticleAssociation const & pa = associations.at(i);

					double dist = sqrt(falgo.SqDist(
								pa.GetRecoVertex(),
								geoalgo::HalfLine(c_start, c_dir*-1)));

					if(fverbose) std::cout << "is distant to the bkw point, " << dist << ", and the shortest_dist so far is "<< best_dist<<".\n";

					if(dist < best_dist) {

						if(fverbose) std::cout << "\t\t\t\t>>Take this candidate vertex as the precandidate vertex and update the shortest-dist!\n";

						best_shower_id = this_iterator->first;
						best_other_id = SIZE_MAX;//shower matches no shower nor track.
						best_vert = falgo.ClosestPt(pa.GetRecoVertex(), this_iterator->second->fcone);//c.second->fcone is a Trajectory_t
						best_dist = dist;
						index = i;//the ith association
						AddShowerToVertex = true;
					}
				}


				if(fverbose) std::cout << "\t\t\t Finish comparison between shower and candidate vertex; \n"<<endl;
			}// >>Finish looking at all showers with other showers, tracks, and vertex candidates.
			//it normally ends up with a shower with a small distance to another shower/track/vertex candidate


			dist_marker(f_dist_sx, best_dist);
			//#3 precandidate vertex promotion
			if(fverbose){
				cout << "STEP III: Try to promote the above precandidate vertex."<<endl;
				std::cout << "The shortest_dist btw the bkw point to other showers/tracks/candidate vertex is " << best_dist<<","<<endl; 
			}

			if(best_dist >= fshower_prox) {//ds>=s_{max}
				if(fverbose) std::cout << "\t larger than "<<fshower_prox<< ", so no shower can be added to any vertex.\n";
				return;
			} else{
				if(fverbose) std::cout << "\t smaller than "<<fshower_prox<< ", so attempt promoting this precandidate vertex.\n";
			}

			//#3 Promote precandidate to candidate vertex, if it is 
			//		1) candidate vertex already, promote it directly
			//		2) shower or track, then create a candidate vertex if needed.

			if(AddShowerToVertex) {
			//#3.1 Precandidate vertex is a candidate vertex, add it and look at next shower;
				if(fverbose) std::cout << "\t\tThe precandidate vertex is a candidate vertex. Take it."<<endl;
				pas.AddShower(index, best_shower_id, best_vert);
				shower_map.erase(best_shower_id);
				continue;
			}
			if(fverbose)std::cout << "\t\tThe precandidate vertex is: ";

			//#3.2 The precandidate vertex is not a candidate vertex, but it can be
			//		1) from a shower and possibly promoted to candiate vertex.
			//		2) from a track and possibly promoted to candiate vertex.
			switch( fdetos->GetRecoType(best_other_id) ){
				case 1://#3.2.1 shower that gives the precandidate vertex; if no candidate vertex is found nearby, create one
					{
						size_t association_index = SIZE_MAX;//a large number;
						double best_association_dist = fcpoa_vert_prox;//dist between precandidate vertex and candidate vertex
						bool found_candidate_vertex = false;

						if(fverbose) {
							std::cout << "A shower bkw point;"<<endl;
							std::cout << "\t\tSearch for a nearby candidate vertex among "<< associations.size()<< " existent candidate vertices:" <<endl;
						}

						for(size_t i = 0; i < associations.size(); ++i) {
							double const dist = best_vert.Dist(associations.at(i).GetRecoVertex());


							if(fverbose){
								std::cout << "\t\t\tVertex candidate (ID: " << i <<") with distant "<< dist;
								std::cout<< ", and the shortest_dist so far is " << best_association_dist << std::endl;
							}

							if(dist < best_association_dist) {//<fcpoa_vert_prox, Max. distance btw shower start & cloest approach (dp_max);
								if(fverbose)
									std::cout << "\t\t\t\t>>Find a nearby candidate vertex, update the shortest_dist.\n";

								association_index = i;//get the index of association that gives the best association distance.
								best_association_dist = dist;
								found_candidate_vertex = true;
							}
						}


						if(fverbose) std::cout << "\t\t\t Finish candidate vertex search.\n"<<endl;


						if(fverbose){
							//							<< "Compare the shortest_association_dist, "<< best_association_dist
							//							<< ", which is required to be less than "
							//							<< fcpoa_vert_prox << ".\n";//max. shower projection distance;
						}

						//if(best_association_dist >= fcpoa_vert_prox){} //this is a yes when the precandidate vertex is close to a candidate vertex
						if(found_candidate_vertex){//candidate vertex is nearby

							if(fverbose) {
								std::cout << "\t\t\t\t>>Find a candidate vertex (ID: "<< association_index; 
								std::cout<<") near the bkw point, so add both showers (ID: ";
								std::cout<<best_shower_id<<" "<<best_other_id<<") to that candidate vertex."<< std::endl;
							}
							pas.AddShower(association_index, best_shower_id, best_vert);
							pas.AddShower(association_index, best_other_id, best_vert);

						} else{//Didnt find any candidate vertex nearby, so try to see if any track end point nearby.
							if(fverbose) std::cout << "\t\tSearch for nearby track end points:"<<endl;

							size_t best_track = SIZE_MAX;
							geoalgo::Point_t const * best_tp = nullptr;
							geoalgo::Point_t const * best_other_tp = nullptr;
							double best_showerp_dist = fcpoa_vert_prox;
							bool promote_the_precandidate = false;

							for(size_t const i : fdetos->GetTrackIndices()) {//look at all tracks;

								geoalgo::Trajectory const & t = fdetos->GetTrack(i).ftrajectory;
								//look at two end points of a track, and see if one of them is close to the precandidate vertex
								double const trackend_dist = t.back().Dist(best_vert);
								if(fverbose){
									std::cout << "\t\t\tStart point from the track (ID: " << i <<"): ";
									std::cout << "distance btw them: " << trackend_dist << ", the shortest: "<< best_showerp_dist << std::endl;
								}
								if(trackend_dist < best_showerp_dist) {
									best_track = i;
									best_tp = &t.back();
									best_other_tp = &t.front();
									best_showerp_dist = trackend_dist;
									promote_the_precandidate = true;
									if(fverbose) std::cout << "\t\t\t\t>>Take this track end point as the closest!\n\n";	
								}

								double const trackstart_dist = t.front().Dist(best_vert);

								if(fverbose){
									std::cout << "\t\t\tEnd point from the track (ID: " << i <<"): ";
									std::cout << "distance btw them: " << trackstart_dist << ", the shortest: "<< best_showerp_dist << std::endl;
								}
								if(trackstart_dist < best_showerp_dist) {
									best_track = i;
									best_tp = &t.front();
									best_other_tp = &t.back();
									best_showerp_dist = trackstart_dist;
									promote_the_precandidate = true;
									if(fverbose) std::cout << "\t\t\t\t>>Take this track end point as the closest!\n\n";	
								}
							}

							dist_marker(f_dist_st, best_dist);//distance between shower&shower to a track end point, if updated. 
							if(fverbose){ 
								if(promote_the_precandidate){
									std::cout << "\t\t\t A nearby track end point is used as a candidate vertex."<<endl;
								} else{
									std::cout<<"\t\t\t No track end point is found, so going to create a new candidate vertex for these two showers."<<endl;
								}
							}

							if(promote_the_precandidate) {//it means at least one track end points is close to the precandidate vertex

								std::vector<size_t> const & index_positions =
									pas.GetAssociationIndicesFromObject(best_track);//iterator that gives the best_track inside pas;
								switch(index_positions.size()){
									case 0://no track end point is found from the candidate vertices;
										{
											if(fverbose) {
												std::cout << "\t\t\t\t>>This is a new track! Create a candidate vertex from ";
												std::cout<<" the mid point btw the bkw point and track end point.\n";
												if(best_tp == nullptr) std::cout << "The track-end point does not exist!\n";//tp = track point
											}

											std::vector<size_t> objects;
											objects.push_back(best_shower_id);//a shower
											objects.push_back(best_other_id);//other shower 
											objects.push_back(best_track);//the new track
											std::vector<geoalgo::Point_t> verts(2, best_vert);//best_vert is the bkw point (of 2 showers);

											verts.push_back(*best_tp);
											pas.AddAssociation(objects,
													verts,
													geoalgo::Sphere(*best_tp, best_dist).Center(),//bounding sphere
													best_dist);//"goodness"
										}
										break;
									case 1://1 end point is found from the candidate vertices;
										{
											size_t const index =
												pas.GetAssociationIndices().at(index_positions.front());

											geoalgo::Point_t const & added_point =
												associations.at(index).GetVertexFromNode(best_track);//add the vertex candidate (obtained from two track end points) from pas;

											double const point_dist = added_point.Dist(*best_tp);//bounding sphere radius, distance btw track end point and vertex candidate;
											double const otherpoint_dist = added_point.Dist(*best_other_tp);//same as above, but with the other track end point;

											if(otherpoint_dist < point_dist) {//another end point is cloeser than the chosen one;

												if(associations.at(index).GetRecoVertex().
														Dist(*best_tp) < fstart_prox) {//mid point of both showers, track end point, and candidate are nearby.
													//if the vertex candidate is close enough to the track end point,
													//		then we add both showers (the track was already associated) to this candidate vertex;
													pas.AddShower(index, best_shower_id, best_vert);
													pas.AddShower(index, best_other_id, best_vert);

												} else {//the candidate vertex is far away from the track end point, 
													//		then make a new candidate vertex for these 2 showers and the track;
													if(fverbose) {
														std::cout << "\t\t\t\t>>This is a new track! Create a candidate vertex from ";
														std::cout<<" the mid point btw the bkw point and track end point.\n";
														if(best_tp == nullptr) std::cout << "The track-end point does not exist!\n";//tp = track point
													}

													std::vector<size_t> showers;
													showers.push_back(best_shower_id);
													showers.push_back(best_other_id);
													showers.push_back(best_track);
													std::vector<geoalgo::Point_t> verts(2, best_vert);
													verts.push_back(*best_tp);
													pas.AddAssociation(showers,
															verts,
															geoalgo::Sphere(*best_tp, best_dist).Center(),
															best_dist);
												}
											} else {//this is cloesd enought; associate both showers to the candidate vertex
												pas.AddShower(index, best_shower_id, best_vert);
												pas.AddShower(index, best_other_id, best_vert);
											}
										}
										break;
									case 2://2 end points are found from the candidate vertices;
										{//for the track end point that is closer to the candidate vertex, 
											//		add associate both showers respecting to that end point to the candidate vertex.
											size_t const indexa =
												pas.GetAssociationIndices().at(index_positions.front());
											double dista = 
												associations.at(indexa).GetRecoVertex().Dist(best_vert);

											size_t const indexb =
												pas.GetAssociationIndices().at(index_positions.back());   
											double distb = 
												associations.at(indexb).GetRecoVertex().Dist(best_vert);

											if(dista < distb) {
												pas.AddShower(indexa, best_shower_id, best_vert);
												pas.AddShower(indexa, best_other_id, best_vert);
											} else {
												pas.AddShower(indexb, best_shower_id, best_vert);      
												pas.AddShower(indexb, best_other_id, best_vert);
											}
										}
										break;
									default: {
												 cout<<"WARNING, if you see this line, it means you find more than two end points on a track."<<endl;	
												 cout<<__FILE__<<" at line "<<__LINE__<<endl;
											 }
								}
							} else {//no track end points are near the precandidate vertex
								//	if(fverbose) std::cout << "\t\t\t\t>>promote the precandidate connecting two showers as a vertex candidate: ";

								std::vector<size_t> showers;
								showers.push_back(best_shower_id);
								showers.push_back(best_other_id);
								std::vector<geoalgo::Point_t> verts(2, best_vert);
								geoalgo::Point_t temp_vert2 = geoalgo::Sphere(best_vert, best_dist).Center();
								pas.AddAssociation(showers,
										verts,
										//geoalgo::Sphere(best_vert, best_dist).Center(),
										temp_vert2,
										best_dist);

								if(fverbose) {
									std::cout << "\t\t\t\t>>No candidate vertex or track end points near the bkw point,"<<endl;
									std::cout<<"\t\t\t\t   so add both showers (ID: ";
									std::cout<<best_shower_id<<" "<<best_other_id<<") to a new candidate vertex: (";
									cout<<temp_vert2.at(0)<<",";
									cout<<temp_vert2.at(1)<<",";
									cout<<temp_vert2.at(2)<<")"<<endl;
								}
							}
						}

						//	if(fverbose)	std::cout << "\t\t\t >>Finish comparing to the shower with ID: " << best_other_id << std::endl;
						shower_map.erase(best_other_id);//When a shower is analyzed, eliminate it;
					}
					break;

				case 2:
					{//#3.2.2 track that gives the precandidate vertex; if the track end point is not a candidate vertex, create one; 
						//best_vert is the shower start.
						if(fverbose) std::cout << "A track end point with ID "<<best_other_id<<endl;

						geoalgo::Trajectory const & t = fdetos->GetTrack(best_other_id).ftrajectory;

						double best_trackend_dist = t.front().Dist(best_vert);//dist btw track end point and precandidate vertex;
						geoalgo::Point_t const * point = &t.front();//precandidate vertex, i.e. track end point;
						geoalgo::Point_t const * otherpoint = &t.back();//another end of the track;

						double const trackend_dist = t.back().Dist(best_vert);

						if(trackend_dist != 0 && trackend_dist < best_trackend_dist) {//compare edges of q track to best_vert point non-zero distance
							if(fverbose) std::cout<<"\t\t\tanother end of the track is closer to the bkw point, so take this."<<endl;
							best_trackend_dist = trackend_dist;
							point = &t.back();
							otherpoint = &t.front();
						}

						if(fverbose)
							std::cout << "\t\t\tShower bkw point is distant to the track end point: "
								<< best_trackend_dist
								<< " compared to the allowed distance: "
								<< fcpoa_trackend_prox << " \n";

						dist_marker(f_dist_sst, best_trackend_dist);
						if(best_trackend_dist < fcpoa_trackend_prox) {//track is close enough to be used as a candidate vertex.

							std::vector<size_t> const index_positions =
								pas.GetAssociationIndicesFromObject(best_other_id);
							//this helps to find the object with best_other_id as value
							// in the pas.

							if(fverbose)
								std::cout << "\t\t\tAdding shower to the existing candidate vertex, "
									<< "# of track end points as the candidate vertex: "
									<< index_positions.size() << std::endl;

							std::vector<size_t> objects;
							std::vector<geoalgo::Point_t> verts;

							switch(index_positions.size()){
								case 0://the track has no end points as cnadidate vertex;
									{
										if(fverbose)
											std::cout << "\t\t\t\tPromote the bkw point as the candidate vertex\n";

										objects.push_back(best_shower_id);
										objects.push_back(best_other_id);
										verts.push_back(best_vert);
										verts.push_back(*point);
										pas.AddAssociation(objects,
												verts,
												geoalgo::Sphere(*point, best_dist).Center(),
												best_dist);
									}
									break;
								case 1://the track has 1 end point as candidate vertex
									{
										if(fverbose)
											std::cout << "\t\t\t\tAdd shower to the candidate vertex or create a new candidate vertex with the track end point\n";

										size_t const index =
											pas.GetAssociationIndices().at(index_positions.front());//get the end point as candidate vertex

										geoalgo::Point_t const & added_point =
											associations.at(index).GetVertexFromNode(best_other_id);

										double const point_dist = added_point.Dist(*point);//point is the close track end point
										double const otherpoint_dist = added_point.Dist(*otherpoint);//other point is the far track end point, introduced because of the candidate vertex;

										if(fverbose)
											std::cout << "\t\t\t\tdist of the candidate to related track end point: "<< otherpoint_dist
											<<" or the other end of that track: "<< point_dist << "\n";

										if(otherpoint_dist < point_dist) {

											if(fverbose)
												std::cout << "\t\t\t\t\tanother end is closer\n"
													<< "\t\t\t\t\tdistance between candidate vertex and the track end point (as used to find the candidate vetex): "
													<< associations.at(index).GetRecoVertex().Dist(*point)
													<< " < track-track distance: " << fstart_prox << " ?\n";

											if(associations.at(index).GetRecoVertex().
													Dist(*point) < fstart_prox) {
												if(fverbose) std::cout << "\t\t\t\t\t\tYes, add the shower to the candidate vertex\n";
												pas.AddShower(index, best_shower_id, best_vert);//Add Shower to the candidate vertex
											} else {

												if(fverbose) std::cout << "\t\t\t\t\t\tNo, add the shower and the other track end point as candidate vertex."<<std::endl;
												objects.push_back(best_shower_id);
												objects.push_back(best_other_id);
												//		std::vector<geoalgo::Point_t> verts;
												verts.push_back(best_vert);
												verts.push_back(*point);
												pas.AddAssociation(objects,
														verts,
														geoalgo::Sphere(*point, best_dist).Center(),
														best_dist);//make a bounding sphere by merging different vertices
											}
										} else{

											if(fverbose)
												std::cout << "\t\t\t\t\tthe candidate vertex is indeed closer\n"
													<< "\t\t\t\t\tAdd the shower with ID: " << best_shower_id<<std::endl;
//													<< " to association: " << index << std::endl;

											pas.AddShower(index, best_shower_id, best_vert);//Add Shower to the candidate vertex
										}
									}
									break;
								case 2://both track end points are candidate vertex
									{
										if(fverbose)
											std::cout << "\t\t\t\tAdd shower the cloest (distant to the bkw projection point) track end point.\n";

										size_t const indexa =
											pas.GetAssociationIndices().at(index_positions.front());
										double dista = 
											associations.at(indexa).GetRecoVertex().Dist(best_vert);

										size_t const indexb = pas.GetAssociationIndices().at(index_positions.back());   
										double distb = 
											associations.at(indexb).GetRecoVertex().Dist(best_vert);

										if(dista < distb) {
											pas.AddShower(indexa,
													best_shower_id,
													best_vert);
										} else {
											pas.AddShower(indexb,
													best_shower_id,
													best_vert);
										}
									}
									break;
								default:
									std::cout << "Warning: more than two indices found, node: "
										<< best_other_id << std::endl;
							}
						} else {

							if(fverbose) std::cout << "\t\t\tToo far. Not adding the shower the existing candidate vertex; create a new candidate vertex for this shower.\n";

							pas.GetDetectorObjects().SetAssociated(best_shower_id);
						}
					}
					break;
				default:{
							cout<<"If you see this line, it means you find an object that is neither track nor shower."<<endl;	
							cout<<__FILE__<<" at line "<<__LINE__<<endl;
						}
			}
			shower_map.erase(best_shower_id);
		}

	}


	void VertexBuilder::AddLoneTracks(ParticleAssociations_all & pas) {
		//fdetos, detectorobjects_all
		for(size_t const gn : fdetos->GetTrackIndices()) {
			//goes over recob::track index, gn

			//Track is a struct defined in DetectorObjects.h
			Track const & t = fdetos->GetTrack(gn);//t is track

			if(t.fis_associated) continue;//this will jump to the end of the loop and goes with next gn; triger this when the track is associated.

			geoalgo::Point_t const * track_end = nullptr;
			double zmin = 2000;//the unit .. 2000 cm??? it doesnt matter though; it changes after first if-else

			geoalgo::Point_t const & front = t.ftrajectory.front();//first point along the trajectory.
			if(front.at(2) < zmin) {//first point's z coordinate compared to zmin
				track_end = &front;
				zmin = front.at(2);
			}

			geoalgo::Point_t const & back = t.ftrajectory.back();
			if(back.at(2) < zmin) {
				track_end = &back;
				zmin = back.at(2);
			}

			if(track_end) {

				pas.AddAssociation(std::vector<size_t>(1, gn),
						std::vector<geoalgo::Point_t>(1, *track_end),
						geoalgo::Sphere(*track_end, 0).Center(),
						0);

			} else{
				std::cout << "Warning: No track end pointer\n";
			}
		}

	}


	void VertexBuilder::AddLoneShowers(ParticleAssociations_all & pas) {

		for(size_t const gn : fdetos->GetShowerIndices()) {

			Shower const & s = fdetos->GetShower(gn);

			if(s.fis_associated) continue;

			geoalgo::Point_t const p = s.fcone.Start();

			pas.AddAssociation(std::vector<size_t>(1, gn),
					std::vector<geoalgo::Point_t>(1, p),
					geoalgo::Sphere(p, 0).Center(),
					0);

		}

	}

	// Fill the TTree with vertex info.
	void VertexBuilder::FillVBT(ParticleAssociations_all & pas) {

		fvbt->ftrack_number = fdetos->GetTrackIndices().size();
		fvbt->fshower_number = fdetos->GetShowerIndices().size();
		fvbt->ftree->Fill();

	}


	void VertexBuilder::Run(ParticleAssociations_all & pas) {//Analysis the tracks & showers objects.
//		int screen_width = 86;

		CheckSetVariables();//Variables are set in SinglePhoton_module.cc.

		fdetos = &pas.GetDetectorObjects();//initialize an empty object;

		//Make two associations, for tracks and showers.
		if(fverbose) print_gadget(">>>>> Associate tracks","=");

		AssociateTracks(pas);//Gives candidate vertex from tracks.

		//fvbt is the object, it means.. if the object is not empty, do something.
		//  if(fvbt) fvbt->fassociation_track_number = pas.GetAssociations().size();


		if(fverbose) print_gadget(">>>> Associate showers","=");

		AssociateShowers(pas);//select associated candidates for showers

		//if(fvbt) fvbt->fassociation_shower_number = pas.GetAssociations().size();

		//CHECK, repeat for long objects? Add them into consideration even if they are not associated to a vertex?
		if(fverbose) print_gadget(">>>> Add lone tracks","=");
		AddLoneTracks(pas);

		if(fverbose) print_gadget(">>>> Add lone showers","=");
		AddLoneShowers(pas);

		if(fverbose) print_gadget(">>>> Get shower associations","=");
		pas.GetShowerAssociations();//CHECK

		if(fvbt) fvbt->fassociation_final_number = pas.GetSelectedAssociations().size();
		if(fverbose)cout<<"\n\n # of Vertex found: "<<pas.GetSelectedAssociations().size()<<endl;

		if(fvbt) {//after the association is finished (found the vertex), fill in in the tree.
			if(fverbose) std::cout << "Fill VBT\n";
			FillVBT(pas);
		}

		pas.NodeCheck();

	}
}


#endif
