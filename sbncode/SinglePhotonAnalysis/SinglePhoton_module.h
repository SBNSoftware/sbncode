#ifndef SINGLE_PHOTON_ANALYSIS
#define SINGLE_PHOTON_ANALYSIS

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//KENG, should use these in sbndcode
//#include "art/Framework/Services/Optional/TFileService.h"
//#include "art/Framework/Services/Optional/TFileDirectory.h"
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
//KENG
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include  "nusimdata/SimulationBase/GTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"

//#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

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

#include "HelperFunctions/helper_math.h"
#include "HelperFunctions/helper_gadget.h"
//#include "Libraries/Atlas.h"

//#include "Libraries/bad_channel_matching.h"
#include "Libraries/DBSCAN.h"

#include "SEAview/SEAviewer.h"


//------------------------------------------------------------------------------------------------------------------------------------------

namespace single_photon
{

	//Create a class based on PFParticle, connected to other Pandora objects
	//CHECK, blending this to the code; not active yet;
	class PandoraPFParticle{

	public:
		//constructor:
		PandoraPFParticle( 
				art::Ptr<recob::PFParticle> input_PFParticle,
				std::vector< art::Ptr< larpandoraobj::PFParticleMetadata > > input_MetaData,
				std::vector< art::Ptr<recob::Vertex > > input_Vertex,
				std::vector< art::Ptr<recob::Cluster> > input_Clusters,
				std::vector< art::Ptr<recob::Shower > > input_Showers,
				std::vector< art::Ptr<recob::Track  > > input_Tracks
		)
		:
		pPFParticle(input_PFParticle),
		pMetaData(input_MetaData),
		pVertex(input_Vertex),
		pClusters(input_Clusters)
		{
			pPFParticleID = pPFParticle->Self();
			pPdgCode = pPFParticle->PdgCode();

			//Get recob::Shower/Track
			const unsigned int nShowers(input_Showers.size());
			const unsigned int nTracks(input_Tracks.size());
			if(nShowers+nTracks != 1) mf::LogDebug("SinglePhoton") << "  No just one shower/track is associated to PFParticle " << pPFParticleID << "\n";

			pHasShower = nShowers;
			pHasTrack = nTracks;
			if(pHasShower == 1) pShower=input_Showers.front();
			if(pHasTrack == 1)  pTrack=input_Tracks.front();

			//set vertex
			if (!pVertex.empty()){
				const art::Ptr<recob::Vertex> vertex = * (pVertex.begin());
				vertex->XYZ(pVertex_pos);
			} 
			//else{
			//	throw art::Exception(art::errors::StdException)
			//		<< " Pandor did not have an associated vertex for a particle. ";
			//}


			//(pVertex.begin()).XYZ(pVertex_pos);
//			std::cout<<"CHECK Vertex position: ("<<pVertex_pos[0]<<","<<pVertex_pos[1]<<","<<pVertex_pos[2]<<")"<<std::endl;
	
			//get ancestor for a 1st generation PFParticle
			if(pPFParticle->IsPrimary()){ 
				pAncestor = pPFParticle;
				pAncestorID = -1;
			}

			//fill in some booleans
			for(auto &meta: pMetaData){
				std::map<std::string, float> propertiesmap  = meta->GetPropertiesMap();
				if(propertiesmap.count("NuScore")==1) pNuScore = propertiesmap["NuScore"];
				if(propertiesmap.count("TrackScore")==1) pTrackScore = propertiesmap["TrackScore"];
				if(propertiesmap.count("IsNeutrino")==1) pIsNeutrino = true; 
				if(propertiesmap.count("IsClearCosmic")==1) pIsClearCosmic = true; 
			}

		};

		art::Ptr< recob::PFParticle > pPFParticle;//d
		art::Ptr< recob::Shower>	pShower;//d with 0 or 1 element
		art::Ptr< recob::Track >	pTrack;//d with 0 or 1 element
		art::Ptr< recob::Slice >	pSlice;//d in helper_connector.h get the id from pSlice->key()
		art::Ptr< recob::PFParticle > pAncestor;//d found by tracing Parent()
		art::Ptr< simb::MCTruth >	pMCTruth;

		std::vector< art::Ptr< larpandoraobj::PFParticleMetadata > > pMetaData;//d
		std::vector< art::Ptr< recob::Vertex > > pVertex;//d
		std::vector< art::Ptr< recob::Hit > >	pHits;//d
		std::vector< art::Ptr< recob::Cluster > > pClusters;//d
		std::vector< art::Ptr< recob::SpacePoint > > pSpacePoints;
		std::vector< art::Ptr< simb::MCParticle > > pMCParticles;
		
		double pVertex_pos[3] = {-9999,-9999,-9999};//d

		double pNuScore = -999 ;//d
		double pTrackScore = -999;//d

		bool pIsNeutrino = false;//d
		bool pIsClearCosmic = false;//d
		bool pIsNuSlice = false;//d
		
		int pHasShower = 0;//d
		int pHasTrack = 0;//d
		int pPdgCode = -999;//d
		int pPFParticleID = -9;//d
		int pAncestorID = -9;//d
		int pSliceID = -9;//d
	};


	//this is for second_shower_search.h
    struct sss_score{
        int plane;
        int cluster_label;
        double point_score;
        int n_hits;

        double mean_wire;
        double max_wire;
        double min_wire;
        double mean_tick;
        double max_tick;
        double min_tick;

        double close_tick;
        double close_wire; /* wire of hit that's closest to vertex */
        double angle;//w.r.t shower primary

        double impact_parameter;

        double max_dist_tick;
        double mean_dist_tick;
        double min_dist_tick;
        double max_dist_wire;
        double mean_dist_wire;
        double min_dist_wire;

        double mean_dist;
        double max_dist;
        double min_dist;

        double pca_0;
        double pca_1;
        double pca_theta;

        int n_wires; /* number of wires hits correspond to */
        int n_ticks;

        bool pass;

        sss_score(int ip, int cl): plane(ip), cluster_label(cl){};
    }; //end of class sss_score

	//this works with SEAview/SEAviewer.h
    class cluster {

        public:

            cluster(int ID, int plane, std::vector<std::vector<double>> &pts, std::vector<art::Ptr<recob::Hit>> &hits) 
				:f_ID(ID), f_plane(plane), f_pts(pts), f_hits(hits) {
				
				std::cout<<"\n\n\nCHECK !! Do we need this?"<<std::endl;
                f_npts = f_pts.size();
                if(pts.size() != hits.size()){
                    std::cerr<<"seaviewer::cluster, input hits and pts not alligned"<<std::endl;
                }
                std::vector<double> wires(f_npts);
                std::vector<double> ticks(f_npts);
                for(int p =0; p< f_npts; ++p){
                    wires[p]=f_pts[p][0];
                    ticks[p]=f_pts[p][1];
                }
                TGraph af_graph(f_npts,&wires[0],&ticks[0]);
                f_graph = af_graph;

            };

//            int getID() {return f_ID;}
//            int getN() {return f_npts;}
//            int getPlane(){ return f_plane;}
//            TGraph * getGraph(){ return &f_graph;}
//            std::vector<art::Ptr<recob::Hit>>  getHits(){return f_hits;}
//            int setSSScore(sss_score & scorein){ f_SSScore = &scorein; return 0;}
//            sss_score * getSSScore(){return f_SSScore;}
//        private:
//
            int f_ID;
            int f_npts;
            int f_plane;
            std::vector<std::vector<double>> f_pts;
            std::vector<art::Ptr<recob::Hit>> f_hits;
            TGraph f_graph;
//            sss_score *f_SSScore;
    }; // end of class cluster



    /**
     *  @brief  SinglePhoton class
     */
    class SinglePhoton : public art::EDFilter
    {
        public:
            // name alias from pandora
            typedef art::ValidHandle< std::vector<recob::PFParticle> > PFParticleHandle;
            typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
            typedef std::vector< art::Ptr<recob::Track> > TrackVector;
            typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
            typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;

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

        private:
            /**
             * @brief: reset/clear data members
             */
            void ClearMeta();
			void CreateMetaBranches();
            /**
             *  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
             *
             *  @param  pfParticleHandle the handle for the PFParticle collection
             *  @param  pfParticleMap the mapping from ID to PFParticle
             */
            void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

            /**
             * @brief Print out scores in PFParticleMetadata
             *
             * @param evt the art event to analyze
             * @param pfParticleHandle the handle for the PFParticle collection
             */
            void PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const;

            /**
             *  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
             *
             *  @param  pfParticleMap the mapping from ID to PFParticle
             *  @param  crParticles a vector to hold the top-level PFParticles reconstructed under the cosmic hypothesis
             *  @param  nuParticles a vector to hold the final-states of the reconstruced neutrino
             */
            void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap,
			const lar_pandora::PFParticlesToVertices &particlesToVertices, 
			PFParticleVector &crParticles, 
			PFParticleVector &nuParticles,
			size_t fpfp_w_bestnuID);

            /**
             *  @brief  Collect associated tracks and showers to particles in an input particle vector
             *
             *  @param  particles a vector holding PFParticles from which to find the associated tracks and showers
             *  @param  pfParticleHandle the handle for the PFParticle collection
             *  @param  evt the art event to analyze
             *  @param  tracks a vector to hold the associated tracks
             *  @param  showers a vector to hold the associated showers
             */
            void CollectTracksAndShowers(const PFParticleVector &particles,const PFParticleIdMap pfParticleMap,  const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap);

            /**
             * @brief: analyze associated tracks/showers for an PFParticle
             * @param: pParticle: PFParticle to be analyzed
             * @param: associatedTracks: a vector of asso track for pParticle
             * @param: associatedShowers: a vector of asso shower for pParticle
             * @param: tracks: associated track will be added into tracks
             * @param: showers: associated shower will be added into showers
             * @param: trackToNuPFParticleMap/showerToNuPFParticleMap: map of associated track/shower to the PFParticle
             */
            void FillTracksAndShowers( const std::vector< art::Ptr<recob::Track> > & associatedTracks, const std::vector< art::Ptr<recob::Shower> > & associatedShowers, const art::Ptr<recob::PFParticle> &pParticle , const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap);

            /**
             * @brief: get vertex for particle
             * @param: particle: a primary neutrino
             */
            void GetVertex(const lar_pandora::PFParticlesToVertices & particlestoVertices, const art::Ptr<recob::PFParticle> & particle );



            double CalcEShowerPlane(const std::vector<art::Ptr<recob::Hit>>& hits, int plane); /* returns energy sum of hits on a certain plane */
            /*
             *@brief Calculated the shower energy by looping over all the hits and summing the charge
             *@param hits -  an art pointer of all the hits in a shower
             * 
			 */
            double CalcEShower(const std::vector<art::Ptr<recob::Hit>> &hits); /* returns max energy on any of the planes */

            /**
             *@brief Takes a hit and multiplies the charge by the gain
             *@param thishitptr art pointer to a hit
             *@param plane the plane the hit is on
             **/
            double GetQHit(art::Ptr<recob::Hit> thishitptr, int plane);


            int getNHitsPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane);
            double getMeanHitWidthPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane);


            /**
             * @brief Calculate the E value in MeV for a given hit
             * @param thishit - an individual hit 
             * 
             *
             * */

            /**
             * @brief Calculate the E value in MeV from a given Q value
             * @param q - the charge value
             * 
             * */
            double QtoEConversion(double q);


            /**
             *@brief Takes a vector of dQ/dx values and converts to dE/dx
             *@param dqdx - vector of dqdx points
             *
             * */
            std::vector<double> CalcdEdxFromdQdx(std::vector<double> dqdx);

            /**
             *
             *@brief For a single shower, calculates the dQdx for each hit in the clusters in the shower for a single plane
             *@param shower - a Pandora shower
             *@param clusters - all of the clusters in the shower
             *@param clusterToHitMap - a map between each cluster and all of the hits in the cluster
             *@param plane - a single plane
             * * */

            std::vector<double> CalcdQdxShower(
                    const art::Ptr<recob::Shower>& shower,
                    const std::vector<art::Ptr<recob::Cluster>> & clusters, 
                    std::map<art::Ptr<recob::Cluster>,    std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,int plane,
					double triggeroffset,
					detinfo::DetectorPropertiesData const & theDetector);
            /**
             *@brief Gets the pitch between the 3D reconstructed shower direction and the wires for a given plane (the dx in dQdx)
             *@param shower_dir - the 3D shower direction
             *@param plane - a single plane
             * */
            double getPitch(TVector3 shower_dir, int plane);  /* distance between hit in shower direction projected on plane */


            /**
             *
             *@brief For a 2d point on a plane in cm and a rectangle, returns true if ponint is inside or on the boundary
             *uses triangle area check
             *
             * */
            bool isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle);



            //----------------  Templatees ----------------------------
            void AnalyzeTemplates();
            void ClearTemplates();
            void ResizeTemplates(size_t);
            void CreateTemplateBranches();

            //---------------- Potential Track Stub --------------------
            void ClearStubs();
            void CreateStubBranches();

            /* @brief: given indices of clusters, determine if they overlap in time
             * @ arguments: cluster_planes, cluster_max_ticks, cluster_min_ticks are full vector containing plane, max/min tick information of all clusters
             * 		    candidate_indices provided the indices of clusters of which we'd like to check the overlap
             */ 
            std::pair<bool, std::vector<double>> clusterCandidateOverlap(const std::vector<int> & candidate_indices, const std::vector<int>& cluster_planes, const std::vector<double>& cluster_max_ticks, const std::vector<double>& cluster_min_ticks);


            /* @brief: given all clusters, and their plane, tick information, find all possible matching clusters using time information
             * @brief: candidate clusters on different plane that overlap in time tick will be grouped together
             * @return: return.first -> number of possible matches
             * 		return.second.first -> 2D vector, indices of clusters in every possible match
             * 		return.second.second -> 1D vector, time overlap fraction of clusters in every possible match
             */	   
            std::pair<int, std::pair<std::vector<std::vector<double>>, std::vector<double>>> GroupClusterCandidate(int num_clusters,  const std::vector<int>& cluster_planes, const std::vector<double>& cluster_max_ticks, const std::vector<double>& cluster_min_ticks);



            //---------------- SecondShower----
            void ClearSecondShowers(); /* reset and clear variable/vectors related to second shower */
            void ResizeSecondShowers(size_t size); /* currently does nothing */ 
            void CreateSecondShowerBranches(); /*create branches in vertex tree for second shower related variables */

            void ClearSecondShowers3D(); /* reset and clear variables/vectors realted to second shower 3D */
            void CreateSecondShowerBranches3D(); /*create branches in vertex tree for second shower 3D */

            void SimpleSecondShowerCluster();


            void SecondShowerSearch3D(std::vector<art::Ptr<recob::Shower>> & showers,std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & NormalShowerToPFParticleMap,  std::vector<art::Ptr<recob::Track>> & tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & normaltrkmap,art::Event const & evt);



            /* brief: analyze hits (second shower), find out which primary MCParticle is the best-match for these hits
             * and return a vector of 7 elements: 
             * {has_found_match, PDG for the ancestor match, PDG for the mother particle of ancestor match, track ID of ancestor match, true energy of ancestor match, fraction of overlay in hits, fraction of energy deposited by the matched ancestor particle on the best-plane}
             */
            std::vector<double>SecondShowerMatching(std::vector<art::Ptr<recob::Hit>>& hitz,
                    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
                    std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
                    std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
                    std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap);


            /* analyze a cluster of hits, and return corresponding sss_score */
            sss_score ScoreCluster(int,int,std::vector<art::Ptr<recob::Hit>>&,double,double, const art::Ptr<recob::Shower>&);

            /* get the nearest N hits surrounding given position, and return a TGraph of wire-tick for these hits 
             * This function is currently used in function 'SecondShowerSearch'
             * @parameter: plane, cluster are not in use
             * @note: need to make sure all the hits in hitz are on the same plane as vertex_wire
             */
            TGraph* GetNearestNpts(int plane,int cluter,std::vector<art::Ptr<recob::Hit>>& hitz,double vertex_wire,double vertex_tick,int Npt);


            /* brief: returns index of the first shower in showers which is close enough to one of hit in hitz */
            int CompareToShowers(int,int,std::vector<art::Ptr<recob::Hit>>& hitz,double,double,
                    const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & showertopfparticlemap,      const   std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfparticletohitsmap,                    double eps);
            //---------------- Isolation ----------------- 

            void ClearIsolation();  /* clear vector members related to isolation */
            void CreateIsolationBranches();  /* create branches for vectors related to isolation in vertex_tree */

            /* for a slice, study the min distance between track hits, shower hits, and unassociated hits */
            void IsolationStudy(
                    const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap,
                    const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,
                    const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,  
                    const std::map<art::Ptr<recob::PFParticle>, int> & pfParticleToSliceIDMap, const std::map<int, std::vector<art::Ptr<recob::Hit>>>& sliceIDToHitsMap,
					detinfo::DetectorPropertiesData const & theDetector,
					int slice_w_bestnuID);



            //----------------  Flashes ----------------------------
            void AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes,  art::Handle<std::vector<sbn::crt::CRTHit>> crthit_h, double evt_timeGPS_nsec,  std::map<art::Ptr<recob::OpFlash>, std::vector< art::Ptr<sbn::crt::CRTHit>>> crtvetoToFlashMap);


            void ClearFlashes();  /* clear and reset all the flash-related vectors/variables */
            void ResizeFlashes(size_t); /* resize flash-related vectors */
            void CreateFlashBranches(); /* create branches for flashes in vertex_tree */

            //----------------  Tracks ----------------------------
            /* @brief: analyze each reco track in vector 'tracks', and save info to track-related data members */
            void AnalyzeTracks(const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & tracktopfparticlemap,
                    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitsMap,
                    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>>> & pfparticletospacepointmap , std::map<int, art::Ptr<simb::MCParticle> > &  MCParticleToTrackIdMap, std::map<int, double> &sliceIdToNuScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool> &PFPToClearCosmicMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
                    std::map<art::Ptr<recob::PFParticle>,double> &PFPToTrackScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool> &PFPToNuSliceMap,
                    PFParticleIdMap &pfParticleMap
                    );

            void ClearTracks();  /* clear track related variable and vectors */
            void ResizeTracks(size_t);  /* resize track related vectors */
            void CreateTrackBranches(); /* create track related branch in vertex tree */
            void AnalyzeTrackCalo(const std::vector<art::Ptr<recob::Track>> &tracks, std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<anab::Calorimetry>>> &trackToCaloMap);


            /* @brief: analyze MCParticle related to recob::Track if it has one 
             * variables starting with 'm_sim_track_' will be updated
             * */
            void RecoMCTracks(const std::vector<art::Ptr<recob::Track>>& tracks,  std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle>> & trackToPFParticleMap, std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > & trackToMCParticleMap,  std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector, std::map< int, art::Ptr<simb::MCParticle> > &      MCParticleToTrackIdMap, 
                    std::map<int, double>& sliceIdToNuScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
                    std::vector<double>& vec);


            /* collect information from anab::sParticleIDAlgScores of reco track */
            void CollectPID(std::vector<art::Ptr<recob::Track>> & tracks,std::map< art::Ptr<recob::Track>, art::Ptr<anab::ParticleID>> & trackToPIDMap);
            TGraph proton_length2energy_tgraph;

            //----------------  Showers ----------------------------
			void ClearShowers();  /* clear and reset shower related vectors/variables, for reco, and sim shower */
			void ResizeShowers(size_t);   /* resize vectors related to reco, and sim showers */
			void CreateShowerBranches();  /* create branches for shower-related vec/vars in vertex_tree */

	void AnalyzeShowers(
			std::vector<PandoraPFParticle> all_PPFPs,
			const std::vector<art::Ptr<recob::Shower>>& showers,  
//			std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap, 
//			std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap, 
//			std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > & pfParticleToClusterMap,
			std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> >  & clusterToHitMap , 
//			std::map<int, double>& sliceIdToNuScoreMap,
//			std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
//			std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap, 
//			std::map<art::Ptr<recob::PFParticle>,bool> &PFPToNuSliceMap, 
//			std::map<art::Ptr<recob::PFParticle>,double> &PFPToTrackScoreMap,
//			PFParticleIdMap &pfParticleMap,
			std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>>& PFPtoShowerReco3DMap,
			double triggeroffset,
			detinfo::DetectorPropertiesData const & theDetector
			);

			void AnalyzeKalmanShowers(
					const std::vector<art::Ptr<recob::Shower>>& showers,
					std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> &showerToPFParticleMap,
					std::map<art::Ptr<recob::PFParticle>,art::Ptr<recob::Track>> & pfParticlesToShowerKalmanMap,
					std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>>>&  kalmanTrackToCaloMap, 
					std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap,
					detinfo::DetectorPropertiesData const & theDetector);

            /**
             * @brief: match showers to MCParticles 
             * @arguments filled during function execution:
             * 		mcParticleVector: vector of mother particles of showers
             * 		objectToMCParticleMap: map of shower to its mother particle
             */
            void showerRecoMCmatching(std::vector<art::Ptr<recob::Shower>>& objectVector,
                    std::map<art::Ptr<recob::Shower>,art::Ptr<simb::MCParticle>>& objectToMCParticleMap,
                    std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
                    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
                    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
                    std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
                    std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
                    std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap,
                    std::map<int, double> & sliceIdToNuScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToNuSliceMap);


            /* tranverse through mcParticleVector, and print out infos for photons */
            int   photoNuclearTesting(std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector);

			//KENG Geometry dimensions; Fiducial volume and SCB (no SCB yet?)
			std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
			std::vector<geo::BoxBoundedGeo> fActiveVolumes;
			double m_tpc_active_XMin;
			double m_tpc_active_YMin;
			double m_tpc_active_ZMin;
			double m_tpc_active_XMax;
			double m_tpc_active_YMax;
			double m_tpc_active_ZMax;

            int setTPCGeom();
            int isInTPCActive(std::vector<double>& vec); /* if point is in active TPC volume */
            double distToTPCActive(std::vector<double>&vec); 
            int distToSCB(double & dist, std::vector<double> &vec); 
					/* if point in active TPC, returns distance from point to closest TPC wall
                                                              * otherwise, returns -999 */


            // ------------ Fid Volume and SCB------------------------- //
//CHECK, no SCB stuff yet;
//            double m_tpc_active_x_low;
//            double m_tpc_active_x_high;
//            double m_tpc_active_y_low;
//            double m_tpc_active_y_high;
//            double m_tpc_active_z_low ;
//            double m_tpc_active_z_high;
//
//            double m_SCB_YX_TOP_y1_array;
//            std::vector<double> m_SCB_YX_TOP_x1_array;
//            std::vector<double> m_SCB_YX_TOP_y2_array;
//            double  m_SCB_YX_TOP_x2_array;
//            double m_SCB_YX_BOT_y1_array;
//            std::vector<double> m_SCB_YX_BOT_x1_array;
//            std::vector<double> m_SCB_YX_BOT_y2_array;
//            double m_SCB_YX_BOT_x2_array;
//
//            double m_SCB_ZX_Up_z1_array ;
//            double m_SCB_ZX_Up_x1_array ;
//            double m_SCB_ZX_Up_z2_array ;
//            double m_SCB_ZX_Up_x2_array ;
//
//            double m_SCB_ZX_Dw_z1_array;
//            std::vector<double> m_SCB_ZX_Dw_x1_array;
//            std::vector<double> m_SCB_ZX_Dw_z2_array;
//            double m_SCB_ZX_Dw_x2_array;


			//Fiducial stuff?
//            int isInSCB(std::vector<double>&);  /* if point is inside SCB */
//            int isInSCB(double cut,std::vector<double>&);
		/* calc the minimum distance from point to the SC boundary,save to dist. 
                                                                     * return value (0, or 1) indicates whether the point is in SCB */

            //---------------- MCTruths ----------------------------

            /**
             * @brief: analyze simb::MCTruth neutrino interaction, update truth related variable ('m_mctruth_***' )
             */
            void AnalyzeMCTruths(std::vector<art::Ptr<simb::MCTruth>> & mcTruthVector,  std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector );
            void ClearMCTruths();
            void ResizeMCTruths(size_t);  /* resize mctruth daughters vectors */
            void CreateMCTruthBranches();

            std::map<int,std::string> is_delta_map;

            //---------------- EventWeight ----------------------------

            /**
             * @brief: fill event weight related variables */
            void AnalyzeEventWeight(art::Event const & e );
            void ClearEventWeightBranches();  /* reset eventweight related variable */
            void CreateEventWeightBranches();  /* create branches for eventweight related variable in eventweight_tree */

            //These three are shameless steals from LArPandorHelper But overlays dont work so this is a direct clone. We will filter out later.


            //---------------- Geant4 ----------------------------

            /**
             * @brief: fill event weight related variables */
            void ClearGeant4Branches();  /* reset eventweight related variable */
            void CreateGeant4Branches();  /* create branches for eventweight related variable in eventweight_tree */
            void AnalyzeGeant4( const    std::vector<art::Ptr<simb::MCParticle>> &mcParticleVector);    





            /**
             * @brief: given an event and a label, collect all the SimChannel with that label
             * @ param: simChannelVector: a vector of SimChannel [to be filled]
             */
            void CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector);

            /**
             * @brief: given a event, and a label, grab all MCParticle with that label, and fill the corresponding map for future use
             * @param: evt: event, label: given label
             * @param: truthToParticles: a map of MCTruth to a vector of MCParticle [to be filled]
             * @param: particlesToTruth: a map of MCParticle to MCTruth [to be filled]
             * @param: MCParticleToTrackIdMap: a map pf MCParticle track ID to itself [to be filled]
             */
            void CollectMCParticles(const art::Event &evt, const std::string &label, std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>              &particlesToTruth, std::map< int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIdMap);
            void BuildMCParticleHitMaps(const art::Event &evt, const std::string &label, const std::vector<art::Ptr<recob::Hit>> &hitVector,   std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  &hitsToParticles, const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap);


            //-------------- Slices/Pandora Metadata ---------------//
            void  ClearSlices(); /* reset and clear variables/vectors related to slice */
            void  ResizeSlices(size_t size);   /* resize vectors related to slice */ 
            void CreateSliceBranches();  /* create slice branches in ncdelta_slice_tree and vertex_tree */


            /** 
             * brief: analyze metadata of PFParticles, and fill in all these maps
             * argument: primaryPFPSliceIdVec, sliceIdToNuScoreMap, PFPToClearCosmicMap,PFPToSliceIdMap
             * PFPToNuSliceMap, PFPToTrackScoreMap will be filled in the function boby
             */
            void AnalyzeSlices(std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > & pfParticleToMetadataMap,
                    PFParticleIdMap &pfParticleMap,
                    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & primaryPFPSliceIdVec,   
                    std::map<int, double> & sliceIdToNuScoreMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,   
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
                    std::map<art::Ptr<recob::PFParticle>,bool>& PFPToNuSliceMap,
                    std::map<art::Ptr<recob::PFParticle>,double>& PFPToTrackScoreMap);

            // void GetPFPsPerSlice( std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
            //        std::map<int, int>& sliceIdToNumPFPsMap );

            std::vector<int>  GetPFPsPerSlice( std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap ); /* get number of PFParticles per slice */



            /* returns numbr of PFParticles that correspond to showers (and not cosmic) per slice */
            std::vector<int> GetNumShowersPerSlice(std::map< art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& showerToPFParticleMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap );
            //        void GetNumShowersPerSlice(std::map< art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& showerToPFParticleMap,
            //                 std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
            //                std::map<int, int>& sliceIdToNumShowersMap );


            /* returns numbr of PFParticles that correspond to tracks  (and not cosmic) per slice */
            std::vector<int> GetNumTracksPerSlice(std::map< art::Ptr<recob::Track>,art::Ptr<recob::PFParticle>>& trackToPFParticleMap,   
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap);

            /** @brief: look at reco showers and reco tracks in the event, together with MCParticle info
             * to determine how many eligible tracks and showers there are in the event 
             */
            void AnalyzeRecoMCSlices(std::string signal_def, std::map<int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIDMap,
                    std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle> > & showerToPFParticleMap, 
                    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec, 
                    std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,
                    std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle> > & trackToNuPFParticleMap,
                    std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > &trackToMCParticleMap,
                    std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap);


            int  m_reco_slice_num; //total number of slices in the event
            std::vector<double> m_reco_slice_nuscore; //vector of the neutrino score for each slice in an event
            int m_reco_slice_shower_num_matched_signal; //the number of sim showers matched an MCP in the signal def
            int m_reco_slice_track_num_matched_signal; //the number of sim showers matched an MCP in the signal def
            std::vector<int> m_reco_slice_shower_matched_sliceId; //the slice id for each matched shower
            std::vector<int> m_reco_slice_track_matched_sliceId; //the slice id for each matched track

            std::vector<int> m_reco_slice_num_pfps; //the total number of PFP's per slice
            std::vector<int> m_reco_slice_num_showers; //the subset of PFP's that are showers, ie number of showers per slice 
            std::vector<int> m_reco_slice_num_tracks; //the subset of PFP's that are tracks

            std::vector<double> m_reco_slice_shower_matched_energy; //the energy for each matched shower
            std::vector<double> m_reco_slice_track_matched_energy; //the energy for each matched track
            std::vector<double> m_reco_slice_shower_matched_conversion; //the conversion distance for each matched shower
            std::vector<double> m_reco_slice_shower_matched_overlay_frac; //fraction of overlay hits for each matched shower
            //std::map<art::Ptr<recob::PFParticle>, double > & pfParticleToNuScoreMap;//is filled during analyze slices


            //-------  matched shower: reco shower that matches to a primary photon + max energy of 3 plane > 20 + definition being ncdelta----- 
            std::vector<double> m_matched_signal_shower_overlay_fraction;
            //std::vector<double> m_matched_signal_shower_conversion_length;
            std::vector<double> m_matched_signal_shower_true_E;  /* energy of the best-matched MCparticle for the shower */
            std::vector<double> m_matched_signal_shower_nuscore; /* the neutrino score of the slice containing the reco shower */
            std::vector<int> m_matched_signal_shower_sliceId;    /* reco shower slice ID */
            std::vector<bool> m_matched_signal_shower_is_clearcosmic;
            int m_matched_signal_shower_num = 0;  /* number of match showers (that has unique best-matched primary photon ?)  */
            std::vector<bool> m_matched_signal_shower_is_nuslice;
            std::vector<int> m_matched_signal_shower_tracks_in_slice; /* number of showers in the same slice as of this reco shower */
            std::vector<int> m_matched_signal_shower_showers_in_slice; /* number of tracks in the same slice as of this reco shower */

            //-------  matched shower -------------------------------------



            //-------- for reco tracks that match to a primary proton ---------
            std::vector<double> m_matched_signal_track_true_E; /*  the true energy of matched MCparticle (proton) */
            std::vector<double> m_matched_signal_track_nuscore;  /* nu score of the slice containing the reco track */
            std::vector<int> m_matched_signal_track_sliceId;
            std::vector<bool> m_matched_signal_track_is_clearcosmic; /* if reco track is in clear cosmic slice */
            //  std::vector<bool> m_matched_signal_track_is_nuslice;
            std::vector<bool> m_matched_signal_track_is_nuslice;
            std::vector<int> m_matched_signal_track_tracks_in_slice; /* num of PFP that are tracks in the slice this reco track is in */
            std::vector<int> m_matched_signal_track_showers_in_slice;

            int m_matched_signal_track_num = 0;  /* num of reco tracks matched to primary proton */ 

            //int m_matched_signal_total_num_slices;

            //---------for reco tracks that match to a primary proton ---------

            bool m_reco_1g1p_is_same_slice;
            bool m_reco_1g1p_is_multiple_slices;
            bool m_reco_1g1p_is_nuslice;
            bool m_reco_1g0p_is_nuslice;
            double m_reco_1g1p_nuscore;
            double  m_reco_1g0p_nuscore;
            bool m_is_matched_1g1p;
            bool m_is_matched_1g0p;
            bool m_no_matched_showers;
            bool m_multiple_matched_showers; //if there is more than 1 eligible shower (match to primary photon, pass energy threshold)
            bool m_multiple_matched_tracks; /* if there is more than 1 eligible track (match to primary proton) */



            //------------------ Delaunay triangle tools -----------//

            double triangle_area(double a1, double a2, double b1, double b2, double c1, double c2); /* returns area of triangles */
            int quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area); /* get number of Delaunay triangle found
                                                                                                     * and total area of these triangles, 
                                                                                                     * save to num_triangles & area */
            int delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area); /* given hits, calc number of hits, Delaunay triangles and total areas on each plane */

            // given a MCParticle, get its corrected vertex
			//CHECK can merge three in one?
            int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected);
            int spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected);
            // given a particle, and input location calculate its corrected true position, so we can compare it to reco
            int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input);

            //databased http://dbdata0vm.fnal.gov:8186/uboonecon_prod/app/data?f=channelstatus_data&t=357812824
            std::vector<std::pair<int,int>> bad_channel_list_fixed_mcc9;
            std::map<int,bool> bad_channel_map_fixed_mcc9;


            /* @brief: given run/subrun/event number, determine if this event is in the selected event list */
            bool IsEventInList(int run, int subrun, int event);

            TRandom3 *rangen;
            std::string m_shower3dLabel;
            std::string m_showerKalmanLabel;
            std::string m_showerKalmanCaloLabel;
            std::string m_pandoraLabel;         ///< The label for the pandora producer
            std::string m_trackLabel;           ///< The label for the track producer from PFParticles
//            std::string m_sliceLabel;          
            std::string m_showerLabel;          ///< The label for the shower producer from PFParticles
            std::string m_caloLabel;            ///< The label for calorimetry associations producer
            std::string m_flashLabel;
            std::string m_geantModuleLabel;
//            std::string m_backtrackerLabel;
            std::string m_hitfinderLabel;
            std::string m_hitMCParticleAssnsLabel;
            std::string m_potLabel;
            std::string m_generatorLabel;
//            std::string m_badChannelLabel;
//            std::string m_badChannelProducer;
//            std::string m_mcTrackLabel;
//            std::string m_mcShowerLabel;
			


            std::string m_pidLabel;            ///< For PID stuff
            std::string m_CRTVetoLabel;
            std::string m_CRTTzeroLabel;
            std::string m_CRTHitProducer;
            std::string m_true_eventweight_label;

            bool m_use_PID_algorithms;
            bool m_use_delaunay;
            int  m_delaunay_max_hits;
            bool m_is_verbose;
            bool m_print_out_event;
            bool m_is_data; // value provided by pset
            bool m_is_overlayed;
            bool m_is_textgen;
            bool m_run_all_pfps;
            bool m_has_CRT;
            bool m_fill_trees;
            bool m_run_pi0_filter; //value provided by pset
            bool m_run_pi0_filter_2g1p;
            bool m_run_pi0_filter_2g0p;

            bool m_runPhotoNuTruth;
            bool m_runTrueEventweight;

            bool m_runSelectedEvent;  //if it should run only selected events
            std::string m_selected_event_list; //full path for the file containing run/subrun/event number of selected events
            std::set<std::vector<int>> m_selected_set;  //set of selected events  	 

            //SEAviwer bits
            bool m_runSEAview;
            double m_SEAviewPlotDistance;   //parameters related to shower-like object finding
            double m_SEAviewHitThreshold;
            double  m_SEAviewDbscanMinPts;
            double m_SEAviewDbscanEps;
            double m_SEAviewMaxPtsLinFit;
            bool   m_SEAviewMakePDF;
            int m_SEAviewNumRecoShower;
            int m_SEAviewNumRecoTrack;

            bool m_runSEAviewStub;
            double m_SEAviewStubHitThreshold; //parameters related to track-like object finding
            double m_SEAviewStubPlotDistance;
            double m_SEAviewStubDbscanMinPts;
            double m_SEAviewStubDbscanEps;
            bool m_SEAviewStubMakePDF;
            int m_SEAviewStubNumRecoShower;
            int m_SEAviewStubNumRecoTrack;

            std::string m_Spline_CV_label; //"eventweight4to4aFix"

            bool m_runCRT;
            double m_DTOffset;
            double  m_Resolution;
            std::string  m_DAQHeaderProducer;//"daq"
            std::ofstream out_stream;

			//SSS parameters
			double m_max_conv_dist;
            double m_mass_pi0_mev;

            double m_exiting_photon_energy_threshold ;
            double m_exiting_proton_energy_threshold ;

            double m_track_calo_min_dEdx;
            double m_track_calo_max_dEdx;
            double m_track_calo_min_dEdx_hits;
            double m_track_calo_trunc_fraction;

			//Keng, DetectorClocks now is declared in each event
//            detinfo::DetectorProperties const * theDetector;// = lar::providerFrom<detinfo::DetectorPropertiesService>();
//            detinfo::DetectorClocks    const *  detClocks ;//= lar::providerFrom<detinfo::DetectorClocksService>();
            spacecharge::SpaceCharge const * SCE;
            geo::GeometryCore const * geom;
            double m_work_function;  //value provided by pset
            double m_recombination_factor; // value provided by pset
            //double m_gain;
            std::vector<double> m_gain_mc; // value provided by pset 
            std::vector<double> m_gain_data; 
            double m_wire_spacing;

            int m_Cryostat;
            int m_TPC;

            double m_width_dqdx_box; // value provided by pset
            double m_length_dqdx_box;

            TTree* run_subrun_tree;
            TTree* pot_tree;
            TTree* vertex_tree;
            TTree* eventweight_tree;
            TTree* ncdelta_slice_tree;

            TTree* geant4_tree;

            TTree* true_eventweight_tree;
            std::map<std::string, std::vector<double>> fmcweight;

            //------------ POT related variables --------------
            int m_number_of_events;
            int m_number_of_events_in_subrun;
            double m_pot_count;
            int m_number_of_vertices;

            int m_run;
            int m_subrun;
            double m_subrun_pot;
            int m_subrun_counts;

            //------------ Event Related Variables -------------
            int m_run_number;
            int m_subrun_number;
            int m_event_number;
            double m_pot_per_event;
            double m_pot_per_subrun;

            int m_test_matched_hits;
            int m_reco_slice_objects;

            //------- Potential Unreconstructed Track Stub related variables ----
            int m_trackstub_num_unassociated_hits; /* number of hits in the slice that're associated with neither shower nor tracks */
            int m_trackstub_unassociated_hits_below_threshold; /*number of unassociated hits that also didn't pass hit threshold,in the slice*/
            int m_trackstub_associated_hits; /* total number of hits from showers and tracks in the slice */


            int m_trackstub_num_candidates; /* number of unasso hit clusters which are not close enough to reco showers */
	    std::vector<int> m_trackstub_candidate_in_nu_slice; /* check if candidate is in neutrino slice: 1->YES, 0->Parts in neutrino slice, -1->Not at all */
            std::vector<int> m_trackstub_candidate_num_hits;
            std::vector<int> m_trackstub_candidate_num_wires; //number of wires spanned by the candidate cluster
            std::vector<int>  m_trackstub_candidate_num_ticks;
            std::vector<int>  m_trackstub_candidate_plane; /* on which plan the unasso cluster is */
            std::vector<double> m_trackstub_candidate_PCA;
            std::vector<double> m_trackstub_candidate_mean_ADC;
            std::vector<double> m_trackstub_candidate_ADC_RMS;
            std::vector<double> m_trackstub_candidate_veto_score;
            std::vector<double> m_trackstub_candidate_mean_tick;
            std::vector<double> m_trackstub_candidate_max_tick;
            std::vector<double> m_trackstub_candidate_min_tick;
            std::vector<double> m_trackstub_candidate_min_wire;
            std::vector<double> m_trackstub_candidate_max_wire;
            std::vector<double> m_trackstub_candidate_mean_wire;
            std::vector<double> m_trackstub_candidate_min_dist;  // min distance from unasso cluter to the vertex */
            std::vector<double> m_trackstub_candidate_min_impact_parameter_to_shower; //min impact parameter of all hits in cluster to the recob::shower direction line (on 2D plane)
            std::vector<double> m_trackstub_candidate_min_conversion_dist_to_shower_start;  //min distance between hits and recob::shower start (on 2D plane)
            std::vector<double> m_trackstub_candidate_min_ioc_to_shower_start;        //min ratio of impact_parameter_to_shower/conversion_dist_to_shower_start of all hits in the cluster
            std::vector<double> m_trackstub_candidate_ioc_based_length;		//length of the cluster, calculated based on the IOC of hit
            std::vector<double> m_trackstub_candidate_wire_tick_based_length;		//length of the cluster, calculated based on the wire & tick span of the cluster
            std::vector<double> m_trackstub_candidate_mean_ADC_first_half;		// mean ADC per hit for the first half of cluster (cluster divided into halves based on hit IOC)
            std::vector<double> m_trackstub_candidate_mean_ADC_second_half;
            std::vector<double> m_trackstub_candidate_mean_ADC_first_to_second_ratio; // ratio of the mean ADC per hit, first half of cluster over second half.
            std::vector<double> m_trackstub_candidate_track_angle_wrt_shower_direction;   //treat cluster as a track, angle between track direction and the shower direction
            std::vector<double> m_trackstub_candidate_linear_fit_chi2;		// chi2 from linear fit of the  {wire, tick} distribution of the cluster
            std::vector<double> m_trackstub_candidate_energy;
            std::vector<int>    m_trackstub_candidate_remerge; // index of the recob::shower candidate cluster is close to (expect it to be -1)
            std::vector<int>    m_trackstub_candidate_matched; /* has matched this unasso cluter to a primary MCParticle: 0-No, 1-Yes */
            std::vector<double> m_trackstub_candidate_matched_energy_fraction_best_plane; /* matched energy fraction of the best-matched MCParticle on best-plane */ 
            std::vector<int>    m_trackstub_candidate_pdg;   /* pdg of the matched MCParticle */
            std::vector<int>    m_trackstub_candidate_parent_pdg;
            std::vector<int>    m_trackstub_candidate_trackid; /* track ID of the matched MCParticle */
            std::vector<double> m_trackstub_candidate_true_energy;  /* true energy of the matched MCParticle */
            std::vector<double> m_trackstub_candidate_overlay_fraction; /* fraction of overlay in the unasso cluster hits */

            //------- grouped stub clusters --------------
            int m_trackstub_num_candidate_groups;	         /* number of groups */ 
            std::vector<std::vector<double>> m_grouped_trackstub_candidate_indices; /* indices of stub clusters that are matched as a group */
            std::vector<double> m_trackstub_candidate_group_timeoverlap_fraction;   /* minimum fraction of the time overlap of grouped stub clusters */



            //------- Second shower related variables ----
            int m_sss_num_unassociated_hits; /* number of hits in the slice that're associated with neither shower nor tracks */
            int m_sss_num_unassociated_hits_below_threshold; /*number of unassociated hits that also didn't pass hit threshold,in the slice*/
            int m_sss_num_associated_hits; /* total number of hits from showers and tracks in the slice */


            //currently commenting this out for speed as its not used
            //ReadBDT * sssVetov1;

            int m_sss_num_candidates; /* number of unasso hit clusters which are not close enough to reco showers */
			std::vector<int> m_sss_candidate_in_nu_slice;
            std::vector<int> m_sss_candidate_num_hits;
            std::vector<int> m_sss_candidate_num_wires; //number of wires spanned by the candidate cluster
            std::vector<int>  m_sss_candidate_num_ticks;
            std::vector<int>  m_sss_candidate_plane; /* on which plan the unasso cluster is */
            std::vector<double> m_sss_candidate_PCA;
            std::vector<double> m_sss_candidate_mean_ADC;
            std::vector<double> m_sss_candidate_ADC_RMS;
            std::vector<double> m_sss_candidate_impact_parameter;
            std::vector<double> m_sss_candidate_fit_slope; //slope of the cluster direction
            std::vector<double> m_sss_candidate_veto_score;
            std::vector<double> m_sss_candidate_fit_constant; //intercept of the cluster direction
            std::vector<double> m_sss_candidate_mean_tick;
            std::vector<double> m_sss_candidate_max_tick;
            std::vector<double> m_sss_candidate_min_tick;
            std::vector<double> m_sss_candidate_min_wire;
            std::vector<double> m_sss_candidate_max_wire;
            std::vector<double> m_sss_candidate_mean_wire;
            std::vector<double> m_sss_candidate_min_dist;  // min distance from unasso cluter to the vertex */
            std::vector<double> m_sss_candidate_wire_tick_based_length;		//length of the cluster, calculated based on the wire & tick span of the cluster
            std::vector<double> m_sss_candidate_energy;
            std::vector<double> m_sss_candidate_angle_to_shower;
            std::vector<double> m_sss_candidate_closest_neighbour;
            std::vector<int>    m_sss_candidate_remerge; // index of the recob::shower candidate cluster is close to (expect it to be -1)
            std::vector<int>    m_sss_candidate_matched; /* has matched this unasso cluter to a primary MCParticle: 0-No, 1-Yes */
            std::vector<double> m_sss_candidate_matched_energy_fraction_best_plane; /* matched energy fraction of the best-matched MCParticle on best-plane */ 
            std::vector<int>    m_sss_candidate_pdg;   /* pdg of the matched MCParticle */
            std::vector<int>    m_sss_candidate_parent_pdg;
            std::vector<int>    m_sss_candidate_trackid; /* track ID of the matched MCParticle */
            std::vector<double> m_sss_candidate_true_energy;
            std::vector<double> m_sss_candidate_overlay_fraction; /* fraction of overlay in the unasso cluster hits */



            //------------ sss3d_showers variables are for reco::showers which are in the events, but not in the slice ----

            int m_sss3d_num_showers;  /* number of showers in the event but not in the slice */
            std::vector<double> m_sss3d_shower_start_x; /* shower start in X axis, for all showers in the event but not in the slice*/
            std::vector<double> m_sss3d_shower_start_y;
            std::vector<double> m_sss3d_shower_start_z;
            std::vector<double> m_sss3d_shower_dir_x; /* shower direction projection on X axis */
            std::vector<double> m_sss3d_shower_dir_y;
            std::vector<double> m_sss3d_shower_dir_z;
            std::vector<double> m_sss3d_shower_length;
            std::vector<double> m_sss3d_shower_conversion_dist; /* dist between shower start and vertex*/

            std::vector<double> m_sss3d_shower_invariant_mass; /* invariant mass of primary recob::shower, and each shower in the event, 
                                                                * calculated assuming vertex is where their mother particle decays */

            std::vector<double> m_sss3d_shower_implied_invariant_mass; /* similar to invariance mass, except this invariant mass  
                                                                        * is calced direclty using shower direction of two showers */

            std::vector<double> m_sss3d_shower_impact_parameter; /* dist between vertex and shower direction line */
            std::vector<double> m_sss3d_shower_ioc_ratio; /* ratio of impact parameter over conversion dist 
                                                           * 0 if the conversion distance is 0*/
            std::vector<double> m_sss3d_shower_energy_max; /* max energy of all planes (E summed from hits) */
            std::vector<double> m_sss3d_shower_score;
            std::vector<int> m_sss3d_slice_nu;
            std::vector<int> m_sss3d_slice_clear_cosmic;

            bool bool_make_sss_plots;


            //------ max_energy, conversion dist, ioc of the sss3d shower that has the smallest ioc parameter ----
            double m_sss3d_ioc_ranked_en;
            double m_sss3d_ioc_ranked_conv;
            double m_sss3d_ioc_ranked_invar;
            double m_sss3d_ioc_ranked_implied_invar;
            double m_sss3d_ioc_ranked_ioc;
            double m_sss3d_ioc_ranked_opang;
            double m_sss3d_ioc_ranked_implied_opang;
            int m_sss3d_ioc_ranked_id; //index of the sss3d_shower that has the smallest ioc.

            // --- same parameters, of the sss3d shower whose implied invariant mass together with primary recob::shower is closest to pi0 mass --
            double m_sss3d_invar_ranked_en;
            double m_sss3d_invar_ranked_conv;
            double m_sss3d_invar_ranked_invar;
            double m_sss3d_invar_ranked_implied_invar;
            double m_sss3d_invar_ranked_ioc;
            double m_sss3d_invar_ranked_opang;
            double m_sss3d_invar_ranked_implied_opang;
            int m_sss3d_invar_ranked_id;


            //--------------- sss2d showers are essentially group of cluters on 3 planes, that have the potential to be a shower -------
            //--------------- they are not recob::showers --------------------------

            // sss2d_ioc_ranked variables are the varaibles (mean shower energy, conv. dist, ioc, etc) of the sss2d shower that has the smallest ioc 
            // sss2d_conv_ranked variables are the varaibles (energy, conv. dist, ioc, etc) of the sss2d shower that has the smallest conv. distance 
            // sss2d_invar_ranked variables are the varaibles of the sss2d shower whose invariant mass together with primary shower is closest to pi0. 
            double m_sss2d_ioc_ranked_en;
            double m_sss2d_ioc_ranked_conv;
            double m_sss2d_ioc_ranked_ioc;
            double m_sss2d_ioc_ranked_pca;
            double m_sss2d_ioc_ranked_invar;
            double m_sss2d_ioc_ranked_angle_to_shower;
            int m_sss2d_ioc_ranked_num_planes;

            double m_sss2d_conv_ranked_en;
            double m_sss2d_conv_ranked_conv;
            double m_sss2d_conv_ranked_ioc;
            double m_sss2d_conv_ranked_pca;
            double m_sss2d_conv_ranked_invar;
            double m_sss2d_conv_ranked_angle_to_shower;
            int m_sss2d_conv_ranked_num_planes;

            double m_sss2d_invar_ranked_en;
            double m_sss2d_invar_ranked_conv;
            double m_sss2d_invar_ranked_ioc;
            double m_sss2d_invar_ranked_pca;
            double m_sss2d_invar_ranked_invar;
            double m_sss2d_invar_ranked_angle_to_shower;
            int m_sss2d_invar_ranked_num_planes;







            //------------ Vertex Related variables -------------
            int m_reco_vertex_size;
            double m_vertex_pos_x;
            double m_vertex_pos_y;
            double m_vertex_pos_z;
            double m_vertex_pos_tick; /* time tick of vertex pos */
            double m_vertex_pos_wire_p0;
            double m_vertex_pos_wire_p2;
            double m_vertex_pos_wire_p1;
            int m_reco_vertex_in_SCB; /* is vertex in SCB: 0- No, 1- Yes */
            double m_reco_vertex_dist_to_SCB; /* dist between vertex to SCB */
            double m_reco_vertex_dist_to_active_TPC; /* dist from vertex to closest active TPC wall, -999 if not in active TPC */


            int m_reco_asso_showers;

            double m_reco_vertex_to_nearest_dead_wire_plane0;
            double m_reco_vertex_to_nearest_dead_wire_plane1;
            double m_reco_vertex_to_nearest_dead_wire_plane2;

            //added eventweight
            //-------------- EventWeight related variables -------------
            static const int k_max_mc_particles=100;

            int m_run_number_eventweight;
            int m_subrun_number_eventweight;
            int m_event_number_eventweight;

            double m_mcflux_nu_pos_x;
            double m_mcflux_nu_pos_y;
            double m_mcflux_nu_pos_z;
            double m_mcflux_nu_mom_x;
            double m_mcflux_nu_mom_y;
            double m_mcflux_nu_mom_z;
            double m_mcflux_nu_mom_E;
            int m_mcflux_ntype;
            int m_mcflux_ptype;
            double m_mcflux_nimpwt;
            double m_mcflux_dk2gen;
            double m_mcflux_nenergyn;
            double m_mcflux_tpx;
            double m_mcflux_tpy;
            double m_mcflux_tpz;
            double m_mcflux_vx;
            double m_mcflux_vy;
            double m_mcflux_vz;
            int m_mcflux_tptype;
            int m_mctruth_nparticles;
            int m_mctruth_particles_track_Id[k_max_mc_particles];
            int m_mctruth_particles_pdg_code[k_max_mc_particles];
            int m_mctruth_particles_mother[k_max_mc_particles];
            int m_mctruth_particles_status_code[k_max_mc_particles];
            int m_mctruth_particles_num_daughters[k_max_mc_particles]; //other similar variables
            int m_mctruth_particles_daughters[100][100];
            double m_mctruth_particles_Gvx[k_max_mc_particles];
            double m_mctruth_particles_Gvy[k_max_mc_particles];
            double m_mctruth_particles_Gvz[k_max_mc_particles];
            double m_mctruth_particles_Gvt[k_max_mc_particles];
            double m_mctruth_particles_px0[k_max_mc_particles];
            double m_mctruth_particles_py0[k_max_mc_particles];
            double m_mctruth_particles_pz0[k_max_mc_particles];
            double m_mctruth_particles_e0[k_max_mc_particles];
            int m_mctruth_particles_rescatter[k_max_mc_particles];
            double m_mctruth_particles_polx[k_max_mc_particles];
            double m_mctruth_particles_poly[k_max_mc_particles];
            double m_mctruth_particles_polz[k_max_mc_particles];
            int m_mctruth_neutrino_ccnc;
            int m_mctruth_neutrino_mode;
            int m_mctruth_neutrino_interaction_type;
            int m_mctruth_neutrino_target;
            int m_mctruth_neutrino_nucleon;
            int m_mctruth_neutrino_quark;
            double m_mctruth_neutrino_w;
            double m_mctruth_neutrino_x;
            double m_mctruth_neutrino_y;
            double m_mctruth_neutrino_qsqr;
            bool m_gtruth_is_sea_quark;
            int m_gtruth_tgt_pdg;
            int m_gtruth_tgt_Z;
            int m_gtruth_tgt_A;
            double m_gtruth_tgt_p4_x;
            double m_gtruth_tgt_p4_y;
            double m_gtruth_tgt_p4_z;
            double m_gtruth_tgt_p4_E;
            double m_gtruth_weight;
            double m_gtruth_probability;
            double m_gtruth_xsec;
            double m_gtruth_diff_xsec;
            int m_gtruth_gphase_space;
            double m_gtruth_vertex_x;
            double m_gtruth_vertex_y;
            double m_gtruth_vertex_z;
            double m_gtruth_vertex_T;
            int m_gtruth_gscatter;
            int m_gtruth_gint;
            int m_gtruth_res_num;
            int m_gtruth_num_piplus;
            int m_gtruth_num_pi0;
            int m_gtruth_num_piminus;
            int m_gtruth_num_proton;
            int m_gtruth_num_neutron;
            bool m_gtruth_is_charm;
            bool m_gtruth_is_strange;
            int m_gtruth_charm_hadron_pdg;
            int m_gtruth_strange_hadron_pdg;
            int m_gtruth_decay_mode;
            double m_gtruth_gx;
            double m_gtruth_gy;
            double m_gtruth_gt;
            double m_gtruth_gw;
            double m_gtruth_gQ2;
            double m_gtruth_gq2;
            int m_gtruth_probe_pdg;
            double m_gtruth_probe_p4_x;
            double m_gtruth_probe_p4_y;
            double m_gtruth_probe_p4_z;
            double m_gtruth_probe_p4_E;
            double m_gtruth_hit_nuc_p4_x;
            double m_gtruth_hit_nuc_p4_y;
            double m_gtruth_hit_nuc_p4_z;
            double m_gtruth_hit_nuc_p4_E;
            double m_gtruth_hit_nuc_pos;
            double m_gtruth_fs_had_syst_p4_x;
            double m_gtruth_fs_had_syst_p4_y;
            double m_gtruth_fs_had_syst_p4_z;
            double m_gtruth_fs_had_syst_p4_E;

            //-------------- Flash related variables -------------
            int m_reco_num_templates;
            std::vector<double> m_reco_template;  /* temp comment: does not seem to be used */


            //-------------- Flash related variables -------------
            std::vector<double> m_reco_flash_total_pe;
            std::vector<double> m_reco_flash_time;
            std::vector<double> m_reco_flash_time_width;
            std::vector<double> m_reco_flash_abs_time;
            std::vector<int>    m_reco_flash_frame;
            std::vector<double> m_reco_flash_ycenter;
            std::vector<double> m_reco_flash_ywidth;
            std::vector<double> m_reco_flash_zcenter;
            std::vector<double> m_reco_flash_zwidth;
            std::vector<double> m_reco_flash_total_pe_in_beamgate;
            std::vector<double> m_reco_flash_time_in_beamgate;
            std::vector<double> m_reco_flash_ycenter_in_beamgate;
            std::vector<double> m_reco_flash_zcenter_in_beamgate;

            int m_reco_num_flashes;
            int m_reco_num_flashes_in_beamgate;

            double m_beamgate_flash_start;
            double m_beamgate_flash_end;


            //----------- CRT related variables -----------------

            //for crt hits from the CRT veto product
            int m_CRT_veto_nhits;  /* number of CRT veto hits */
            std::vector<double> m_CRT_veto_hit_PE;  

            //fields storing information about the CRT hit closest to the flash
            double m_CRT_min_hit_time;
            double m_CRT_min_hit_PE;
            double m_CRT_min_hit_x;
            double m_CRT_min_hit_y;
            double m_CRT_min_hit_z;

            //Fields storing information about all CRT hits in event
            std::vector<double> m_CRT_hits_time;
            std::vector<double> m_CRT_hits_PE;
            std::vector<double> m_CRT_hits_x;
            std::vector<double> m_CRT_hits_y;
            std::vector<double> m_CRT_hits_z;

            double m_CRT_dt; //time between flash and nearest CRT hit

            //------------ Track related Variables -------------
            int m_reco_asso_tracks;  /* number of track. (temp: will figure out what associate means later) */
            std::vector<int> m_reco_track_num_daughters;
            std::vector<double> m_reco_track_daughter_trackscore;  /* track score of this reco track's first daughter */
            std::vector<double> m_reco_track_length;  /* whole length of the reco track */
            std::vector<double> m_reco_track_dirx;    /* need to understand what the pair track->Direction() returns means*/
            std::vector<double> m_reco_track_diry;
            std::vector<double> m_reco_track_dirz;
            std::vector<double> m_reco_track_startx;  /* start pos of the track in cartesian X */
            std::vector<double> m_reco_track_starty;
            std::vector<double> m_reco_track_startz;
            std::vector<double> m_reco_track_endx;    /* end of the track in cartesian X */
            std::vector<double> m_reco_track_endy;
            std::vector<double> m_reco_track_endz;
            std::vector<double> m_reco_track_end_dist_to_active_TPC;  /* min distance from track end to TPC active volume boundaries */
            std::vector<double> m_reco_track_start_dist_to_active_TPC; /* min dist from trk start to TPC active boundaries */
            std::vector<double> m_reco_track_end_dist_to_SCB;          /* min dist from track end to SCB */
            std::vector<double> m_reco_track_start_dist_to_SCB;
            std::vector<int> m_reco_track_end_in_SCB;   /* if track end is in SCB boundary, 1- yes, 0- no */
            std::vector<int> m_reco_track_start_in_SCB;
            std::vector<double> m_reco_track_calo_energy_plane0;  /* energy sum of hits on plane 0 that correspond to the reco track */
            std::vector<double> m_reco_track_calo_energy_plane1;
            std::vector<double> m_reco_track_calo_energy_plane2;
            std::vector<double> m_reco_track_calo_energy_max;	  /* max energy of 3 plane for the reco track */


            std::vector<double>   m_reco_track_theta_yz; /* theta, phi of the track */
            std::vector<double>   m_reco_track_phi_yx;


            std::vector<int> m_reco_track_num_trajpoints;	/* number of valid points in the track */
            std::vector<int> m_reco_track_num_spacepoints;	/* number of recob::spacepoints coresponding to the reco track */
            std::vector<double> m_reco_track_proton_kinetic_energy; /* energy of the track, under the asssumption it's a proton track 
                                                                     * set to -9999 if m_run_pi0_filter is set to true */

            std::vector<size_t>  m_reco_track_ordered_energy_index; /* index of m_reco_track_proton_kinetic_energy such that element values are in descending order */
            std::vector<size_t>  m_reco_track_ordered_displacement_index; /* index of m_reco_track_length so that track length are in descending order */
            std::vector<double> m_reco_track_spacepoint_principal0; /* PCA of reco track (in 3D spacepoint) */
            std::vector<double> m_reco_track_spacepoint_principal1;
            std::vector<double> m_reco_track_spacepoint_principal2;

            std::vector<double> m_reco_track_spacepoint_chi;  /* essentially sum of square of distances between spacepoint and the track line*/
            std::vector<double> m_reco_track_spacepoint_max_dist; /* max distance between a track and its coresponding spacepoints */


            // ---- corresponding variables on the best plane of reco track, which is defined as such------
            // if plane 2 have good hits, then plane 2 is the best-plane
            // otherwise, which plane of plane 0 and 1 has more good hits will be best plane
            // if none of 3 planes has good hits, then best-plane is set to -1
            std::vector<int> m_reco_track_best_calo_plane;
            std::vector<double> m_reco_track_mean_dEdx_best_plane;
            std::vector<double> m_reco_track_mean_dEdx_start_half_best_plane;
            std::vector<double> m_reco_track_mean_dEdx_end_half_best_plane;
            std::vector<int> m_reco_track_good_calo_best_plane;
            std::vector<std::vector<double>> m_reco_track_trunc_dEdx_best_plane;
            std::vector<double> m_reco_track_mean_trunc_dEdx_best_plane;
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_best_plane;
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_best_plane;
            std::vector<double> m_reco_track_trunc_PIDA_best_plane;
            std::vector<std::vector<double>> m_reco_track_resrange_best_plane;
            std::vector<std::vector<double>> m_reco_track_dEdx_best_plane;



            std::vector<double> m_reco_track_mean_dEdx_p0;  /* mean dEdx of hits on plane 0 of the reco track */
            std::vector<double> m_reco_track_mean_dEdx_start_half_p0; /* mean dEdx of first half of the track */
            std::vector<double> m_reco_track_mean_dEdx_end_half_p0;
            std::vector<int> m_reco_track_good_calo_p0; /* number of good dEdx hits on plane 0 of track calorimetry */
            std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p0;
            std::vector<double> m_reco_track_mean_trunc_dEdx_p0;  /* mean of truncated dEdx's of good hits */
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p0; /*mean of first half of trucated dEdx's of good hits */
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p0;
            std::vector<double> m_reco_track_trunc_PIDA_p0; /* mean of constant A in residual range formula, calc'd from good hits */
            std::vector<std::vector<double>> m_reco_track_resrange_p0; /* vec of residual range of good hits per reco track */
            std::vector<std::vector<double>> m_reco_track_dEdx_p0;     /* vec of dEdx of good hits per reco track */

            std::vector<double> m_reco_track_mean_dEdx_p1;
            std::vector<double> m_reco_track_mean_dEdx_start_half_p1;
            std::vector<double> m_reco_track_mean_dEdx_end_half_p1;
            std::vector<int> m_reco_track_good_calo_p1;
            std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p1;
            std::vector<double> m_reco_track_mean_trunc_dEdx_p1;
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p1;
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p1;
            std::vector<double> m_reco_track_trunc_PIDA_p1;
            std::vector<std::vector<double>> m_reco_track_resrange_p1;
            std::vector<std::vector<double>> m_reco_track_dEdx_p1;

            std::vector<double> m_reco_track_mean_dEdx_p2;
            std::vector<double> m_reco_track_mean_dEdx_start_half_p2;
            std::vector<double> m_reco_track_mean_dEdx_end_half_p2;
            std::vector<int> m_reco_track_good_calo_p2;
            std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p2;
            std::vector<double> m_reco_track_mean_trunc_dEdx_p2;
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p2;
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p2;
            std::vector<double> m_reco_track_trunc_PIDA_p2;
            std::vector<std::vector<double>> m_reco_track_resrange_p2;
            std::vector<std::vector<double>> m_reco_track_dEdx_p2;


            std::vector<int> m_reco_track_num_calo_hits_p0; /* number of hits in calorimetry on plane 0 of each reco track */
            std::vector<int> m_reco_track_num_calo_hits_p1;
            std::vector<int> m_reco_track_num_calo_hits_p2;


            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane0; /* distance between track end and the nearest dead wire on plane*/
            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane1;
            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane2;

            std::vector<int> m_reco_track_sliceId; //the slice id for the slice continaing the reco track
            std::vector<double> m_reco_track_nuscore; //the neutrino score of the slice containing the reco track
            std::vector<bool> m_reco_track_isclearcosmic;//true if reco track is in a clear cosmic slice
            std::vector<double> m_reco_track_trackscore; /* track score of reco track, -999 if track is not found in PFPToTrackScoreMap */
            std::vector<int> m_reco_track_pfparticle_pdg; /* PDG of track's corresponding PFParticle, -999 if track is not found in PFPToTrackScoreMap*/
            std::vector<bool> m_reco_track_is_nuslice;  /* if reco track is in a neutrino slice */




            std::vector<int> m_sim_track_matched;  /* if reco track has been matched to a MCParticle, 1-YES, 0-NO */

            //-------- energy, mass, pdg ..etc.. of the matched MCParticle of reco track -----
            std::vector<double> m_sim_track_overlay_fraction;
            std::vector<double> m_sim_track_energy;
            std::vector<double> m_sim_track_mass;
            std::vector<double> m_sim_track_kinetic_energy;
            std::vector<int> m_sim_track_pdg;
            std::vector<int> m_sim_track_parent_pdg;

            /* event origin types:
             * kUnknown: ???	
             * kBeamNeutrino: Beam neutrinos.
             * kCosmicRay: Cosmic rays.
             * kSuperNovaNeutrino: Supernova neutrinos.
             * kSingleParticle: single particles thrown at the detector
             */
            std::vector<int> m_sim_track_origin;   /* truth origin of the matched MCParticle */
            std::vector<std::string> m_sim_track_process;
            std::vector<double> m_sim_track_startx;  /* space-charge corrected start point of the match MCParticle */
            std::vector<double> m_sim_track_starty;
            std::vector<double> m_sim_track_startz;
            std::vector<double> m_sim_track_px;
            std::vector<double> m_sim_track_py;
            std::vector<double> m_sim_track_pz;
            std::vector<double> m_sim_track_endx;  /* space-charge corrected end-point of the matched MCParticle */
            std::vector<double> m_sim_track_endy;
            std::vector<double> m_sim_track_endz;
            std::vector<double> m_sim_track_length; /* track length calculated based on the SC-corrected start and end point of the matched MCParticle */

            std::vector<int> m_sim_track_trackID;
            //-------- energy, mass, pdg ..etc.. of the matched MCParticle of reco track -----



            std::vector<int> m_sim_track_sliceId; //the slice id for the slice continaing the sim track, based on corresponding recob:PFP
            std::vector<double> m_sim_track_nuscore; //the neutrino score of the slice containing the sim track
            std::vector<bool> m_sim_track_isclearcosmic;//true if sim track is in a clear cosmic slice

            /*-------------------------------------------------------------------------------------*/
            std::vector<double> m_isolation_min_dist_trk_shr; /* minimum distance betwee shower hits and track hits on each plane 
                                                               * if there is no shower hits, set to 999
                                                               * if there is shower hits but no track hits, set to -999
                                                               */	
            std::vector<double> m_isolation_nearest_shr_hit_to_trk_wire; /* the wire number of shower hit closest to track hits */	
            std::vector<double> m_isolation_nearest_shr_hit_to_trk_time; /* the time tick of shower hit closest to track hits in the slice */


            std::vector<double> m_isolation_num_shr_hits_win_1cm_trk; /* number of shower hits whose min distance to track hits <= 1cm 
                                                                       * of each plane (this is a 3 element vector) */	    
            std::vector<double> m_isolation_num_shr_hits_win_2cm_trk;	    
            std::vector<double> m_isolation_num_shr_hits_win_5cm_trk;	    
            std::vector<double> m_isolation_num_shr_hits_win_10cm_trk;	    


            std::vector<double> m_isolation_min_dist_trk_unassoc; /* of all unassociated hits, min distance to closest track hits 
                                                                   * set to -999 if there is no unassociated hits or track hits on plane
                                                                   */
            std::vector<double> m_isolation_nearest_unassoc_hit_to_trk_wire;/* wire number of the unassociated hit that of all is nearest to track hits in the slice */	
            std::vector<double> m_isolation_nearest_unassoc_hit_to_trk_time; /* time tick of the unasso hit that is nearest to track hits in the slice */
            std::vector<double> m_isolation_num_unassoc_hits_win_1cm_trk; /* number of unasso hits whose min distance to track hits <= 1cm
                                                                           * on each plane (this vector has 3 elements) */ 
            std::vector<double> m_isolation_num_unassoc_hits_win_2cm_trk;	    
            std::vector<double> m_isolation_num_unassoc_hits_win_5cm_trk;	    
            std::vector<double> m_isolation_num_unassoc_hits_win_10cm_trk;	    


            /*-------------------------------------------------------------------------------------*/

            //------------ Shower related Variables  -------------

            std::vector<int> m_reco_shower_num_daughters;
            std::vector<double> m_reco_shower_daughter_trackscore;

            std::vector<int>   m_reco_shower3d_exists;
            std::vector<double>   m_reco_shower3d_startx;
            std::vector<double>   m_reco_shower3d_starty;
            std::vector<double>   m_reco_shower3d_startz;
            std::vector<double>   m_reco_shower3d_dirx;
            std::vector<double>   m_reco_shower3d_diry;
            std::vector<double>   m_reco_shower3d_dirz;
            std::vector<double>   m_reco_shower3d_theta_yz; /* theta, phi of the 3D shower (direction) */
            std::vector<double>   m_reco_shower3d_phi_yx;

            std::vector<double> m_reco_shower3d_openingangle;
            std::vector<double> m_reco_shower3d_length;
            std::vector<double> m_reco_shower3d_conversion_distance;

            std::vector<double>   m_reco_shower3d_impact_parameter;  /* distance between vertex and 3D shower direction */
            std::vector<double>    m_reco_shower3d_implied_dirx; /* X component of the unit vector point from vertex to 3D shower start */
            std::vector<double>     m_reco_shower3d_implied_diry;
            std::vector<double>     m_reco_shower3d_implied_dirz;

            std::vector<double> m_reco_shower3d_energy_plane0;
            std::vector<double> m_reco_shower3d_energy_plane1;
            std::vector<double> m_reco_shower3d_energy_plane2;

            std::vector<double> m_reco_shower3d_dEdx_plane0;
            std::vector<double> m_reco_shower3d_dEdx_plane1;
            std::vector<double> m_reco_shower3d_dEdx_plane2;



            std::vector<double>   m_reco_shower_startx;
            std::vector<double>   m_reco_shower_starty;
            std::vector<double>   m_reco_shower_startz;
            std::vector<double> m_reco_shower_start_dist_to_active_TPC; /* distance from shower start to closest TPC wall */
            std::vector<double> m_reco_shower_start_dist_to_SCB;
            std::vector<int> m_reco_shower_start_in_SCB;
            std::vector<double> m_reco_shower_end_dist_to_active_TPC;
            std::vector<double> m_reco_shower_end_dist_to_SCB;

            std::vector<double>   m_reco_shower_dirx; /* X component of shower direction */
            std::vector<double>   m_reco_shower_diry;
            std::vector<double>   m_reco_shower_dirz;
            std::vector<double>   m_reco_shower_theta_yz; /* theta, phi of the shower direction */
            std::vector<double>   m_reco_shower_phi_yx;

            std::vector<double> m_reco_shower_openingangle;
            std::vector<double> m_reco_shower_length;
            std::vector<double> m_reco_shower_conversion_distance;  /* distance between shower start and vertex */

            std::vector<double>   m_reco_shower_impact_parameter; /* distance from vertex to the shower direction */
            std::vector<double>    m_reco_shower_implied_dirx; /* the X component of the unit vector pointing from vertex to shower start */
            std::vector<double>     m_reco_shower_implied_diry;
            std::vector<double>     m_reco_shower_implied_dirz;

            std::vector<int> m_reco_shower_delaunay_num_triangles_plane0; /* num of delaunay triangles found on plane 0 for each shower */
            std::vector<int> m_reco_shower_delaunay_num_triangles_plane1;
            std::vector<int> m_reco_shower_delaunay_num_triangles_plane2;

            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane0;/* dist from shower start to nearest dead wire on plane 0 */
            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane1;
            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane2;


            //shower flash matching

            std::vector<double> m_reco_shower_flash_shortest_distz;
            std::vector<double> m_reco_shower_flash_shortest_disty;
            std::vector<double> m_reco_shower_flash_shortest_distyz;

            std::vector<int> m_reco_shower_flash_shortest_index_z;
            std::vector<int> m_reco_shower_flash_shortest_index_y;
            std::vector<int> m_reco_shower_flash_shortest_index_yz;

            double  m_flash_optfltr_pe_beam;
            double  m_flash_optfltr_pe_beam_tot;
            double  m_flash_optfltr_pe_veto;
            double  m_flash_optfltr_pe_veto_tot;

            //end flash matching
            std::vector<int> m_reco_shower_num_hits_plane0; /* number of hits on plane 0 for each shower */
            std::vector<int> m_reco_shower_num_hits_plane1;
            std::vector<int> m_reco_shower_num_hits_plane2;


            std::vector<double> m_reco_shower_delaunay_area_plane0; /* total area of delaunay triangles found on plane 0 for each shower */
            std::vector<double> m_reco_shower_delaunay_area_plane1;
            std::vector<double> m_reco_shower_delaunay_area_plane2;

            std::vector<int> m_reco_shower_sliceId; //the slice id for the slice continaing the reco shower
            std::vector<double> m_reco_shower_nuscore; //the neutrino score of the slice containing the reco shower
            std::vector<bool> m_reco_shower_isclearcosmic;//true if reco shower is in a clear cosmic slice
            std::vector<bool> m_reco_shower_is_nuslice;//true if reco shower is in a clear cosmic slice
            std::vector<double> m_reco_shower_trackscore;
            std::vector<double> m_reco_shower_pfparticle_pdg;

            std::vector<double> m_reco_shower_kalman_exists; /* if there is a kalman track and reco::Calo related to this shower - 0, 1 */
            std::vector<double>   m_reco_shower_kalman_median_dEdx_plane0;
            std::vector<double>     m_reco_shower_kalman_median_dEdx_plane1;
            std::vector<double>   m_reco_shower_kalman_median_dEdx_plane2;
            std::vector<double>   m_reco_shower_kalman_median_dEdx_allplane;
            std::vector<double>      m_reco_shower_kalman_mean_dEdx_plane0;
            std::vector<double>    m_reco_shower_kalman_mean_dEdx_plane1;
            std::vector<double>    m_reco_shower_kalman_mean_dEdx_plane2;




            std::vector<int> m_sim_shower_matched;  /* whether shower has been matched to a MCParticle, 0 - False, 1 - True */

            // ----- energy, mass, pdg ... of the best-matched MCParticle for the shower ------
            std::vector<double> m_sim_shower_energy;
            std::vector<double> m_sim_shower_kinetic_energy;
            std::vector<double> m_sim_shower_mass;
            std::vector<int> m_sim_shower_pdg;
            std::vector<int> m_sim_shower_trackID;
            std::vector<int> m_sim_shower_parent_pdg;
            std::vector<int> m_sim_shower_parent_trackID;
            std::vector<int> m_sim_shower_origin;
            std::vector<std::string> m_sim_shower_process;
            std::vector<std::string> m_sim_shower_end_process;
            // ----- energy, mass, pdg ... of the best-matched MCParticle for the shower ------


            std::vector<double> m_sim_shower_start_x;  /* space charge corrected shower starting point */
            std::vector<double> m_sim_shower_start_y;
            std::vector<double> m_sim_shower_start_z;
            std::vector<double> m_sim_shower_vertex_x;  /* spacecharge corrected shower vertex */
            std::vector<double> m_sim_shower_vertex_y;
            std::vector<double> m_sim_shower_vertex_z;

            std::vector<double> m_sim_shower_px;
            std::vector<double> m_sim_shower_py;
            std::vector<double> m_sim_shower_pz;


            std::vector<int> m_sim_shower_is_true_shower;
            std::vector<int> m_sim_shower_best_matched_plane;
            std::vector<double> m_sim_shower_matched_energy_fraction_plane0; /* fraction of energy of the best-matched mother for shower on 
                                                                              * plane 0 over all energy deposited on plane 0 by the shower */
            std::vector<double> m_sim_shower_matched_energy_fraction_plane1;
            std::vector<double> m_sim_shower_matched_energy_fraction_plane2;
            std::vector<double> m_sim_shower_overlay_fraction; /* fraction of hits from overlay over all hits in the shower */

            std::vector<int> m_sim_shower_sliceId; //the slice id for the slice continaing the sim shower matched to reco
            std::vector<double> m_sim_shower_nuscore; //the neutrino score of the slice containing the sim shower matched to reco
            std::vector<bool> m_sim_shower_isclearcosmic;//true if sim shower matched to reco is in a clear cosmic slice
            std::vector<bool> m_sim_shower_is_nuslice;//true if sim shower matched to reco is in a clear cosmic slice



            //------------ MCTruth related Variables  -------------
            int m_mctruth_num;
            int m_mctruth_origin;
            double m_mctruth_nu_E;
            double m_mctruth_nu_vertex_x;
            double m_mctruth_nu_vertex_y;
            double m_mctruth_nu_vertex_z;
            double m_mctruth_reco_vertex_dist;
            double m_mctruth_lepton_E;
            int m_mctruth_nu_pdg;
            int m_mctruth_lepton_pdg;
            int m_mctruth_mode ;
            int m_mctruth_interaction_type ;
            int m_mctruth_ccnc;
            double m_mctruth_qsqr;

            int m_mctruth_num_daughter_particles;
            std::vector<int> m_mctruth_daughters_pdg;
            std::vector<double> m_mctruth_daughters_E;

            std::vector<int> m_mctruth_daughters_status_code;
            std::vector<int> m_mctruth_daughters_trackID;
            std::vector<int> m_mctruth_daughters_mother_trackID;
            std::vector<double> m_mctruth_daughters_px;
            std::vector<double> m_mctruth_daughters_py;
            std::vector<double> m_mctruth_daughters_pz;
            std::vector<double> m_mctruth_daughters_startx;
            std::vector<double> m_mctruth_daughters_starty;
            std::vector<double> m_mctruth_daughters_startz;
            std::vector<double> m_mctruth_daughters_time;
            std::vector<double> m_mctruth_daughters_endx;
            std::vector<double> m_mctruth_daughters_endy;
            std::vector<double> m_mctruth_daughters_endz;
            std::vector<double> m_mctruth_daughters_endtime;
            std::vector<std::string> m_mctruth_daughters_process;
            std::vector<std::string> m_mctruth_daughters_end_process;


            int     m_mctruth_num_exiting_photons ;
            int      m_mctruth_num_exiting_protons ;
            int    m_mctruth_num_exiting_pi0 ;
            int   m_mctruth_num_exiting_pipm ;
            int   m_mctruth_num_exiting_neutrons; 
            int   m_mctruth_num_exiting_delta0; 
            int   m_mctruth_num_exiting_deltapm; 
            int   m_mctruth_num_exiting_deltapp; 

            double m_mctruth_leading_exiting_proton_energy;

            int m_mctruth_is_delta_radiative;
            int m_mctruth_delta_radiative_1g1p_or_1g1n;
            double m_mctruth_delta_photon_energy;
            double m_mctruth_delta_proton_energy;
            double m_mctruth_delta_neutron_energy;
            std::vector<int> m_mctruth_exiting_delta0_num_daughters;

            std::vector<int> m_mctruth_exiting_photon_trackID;
            std::vector<int> m_mctruth_exiting_photon_mother_trackID;
            std::vector<int> m_mctruth_exiting_photon_from_delta_decay;
            std::vector<double> m_mctruth_exiting_photon_energy;
            std::vector<double> m_mctruth_exiting_photon_px;
            std::vector<double> m_mctruth_exiting_photon_py;
            std::vector<double> m_mctruth_exiting_photon_pz;

            std::vector<int> m_mctruth_exiting_proton_trackID;
            std::vector<int> m_mctruth_exiting_proton_mother_trackID;
            std::vector<int> m_mctruth_exiting_proton_from_delta_decay;
            std::vector<double> m_mctruth_exiting_proton_energy;
            std::vector<double> m_mctruth_exiting_proton_px;
            std::vector<double> m_mctruth_exiting_proton_py;
            std::vector<double> m_mctruth_exiting_proton_pz;

            std::vector<int> m_mctruth_exiting_neutron_trackID;
            std::vector<int> m_mctruth_exiting_neutron_mother_trackID;
            std::vector<int> m_mctruth_exiting_neutron_from_delta_decay;
            std::vector<double> m_mctruth_exiting_neutron_energy;
            std::vector<double> m_mctruth_exiting_neutron_px;
            std::vector<double> m_mctruth_exiting_neutron_py;
            std::vector<double> m_mctruth_exiting_neutron_pz;


            int  m_mctruth_num_reconstructable_protons;

            bool  m_mctruth_is_reconstructable_1g1p;
            bool  m_mctruth_is_reconstructable_1g0p;


            std::vector<double>        m_mctruth_exiting_pi0_E;
            std::vector<double>        m_mctruth_exiting_pi0_mom;
            std::vector<double>        m_mctruth_exiting_pi0_px;
            std::vector<double>        m_mctruth_exiting_pi0_py;
            std::vector<double>        m_mctruth_exiting_pi0_pz;

            double m_mctruth_pi0_leading_photon_energy;
            std::string m_mctruth_pi0_leading_photon_end_process;
            double m_mctruth_pi0_subleading_photon_energy;
            std::string m_mctruth_pi0_subleading_photon_end_process;
            std::vector<double> m_mctruth_pi0_subleading_photon_end;
            std::vector<double> m_mctruth_pi0_subleading_photon_start;
            std::vector<double> m_mctruth_pi0_leading_photon_end;
            std::vector<double> m_mctruth_pi0_leading_photon_start;
            int    m_mctruth_pi0_leading_photon_exiting_TPC;
            int    m_mctruth_pi0_subleading_photon_exiting_TPC;
            std::vector<double> m_mctruth_pi0_leading_photon_mom;
            std::vector<double> m_mctruth_pi0_subleading_photon_mom;

            std::string  m_truthmatching_signaldef;

            //the calo calculated quantities 
            std::vector<double> m_reco_shower_energy_max; //for each hit in a shower, converts Q->E, and sums. The max energy of all planes
            std::vector<double> m_reco_shower_energy_plane0; /* shower energy (summed hit energy) on plan 0 */
            std::vector<double> m_reco_shower_energy_plane1;
            std::vector<double> m_reco_shower_energy_plane2;

            std::vector<double> m_reco_shower_reclustered_energy_max;
            std::vector<double> m_reco_shower_reclustered_energy_plane0; /* total energy of the reco shower, and unassociated hit clusters 
                                                                          * close enough to it */
            std::vector<double> m_reco_shower_reclustered_energy_plane1;
            std::vector<double> m_reco_shower_reclustered_energy_plane2;


            std::vector<double> m_reco_shower_plane0;
            std::vector<double> m_reco_shower_plane1;
            std::vector<double> m_reco_shower_plane2;

            std::vector<double> m_reco_shower_plane0_nhits; /* num of shower hits on plane 0 */
            std::vector<double> m_reco_shower_plane1_nhits;
            std::vector<double> m_reco_shower_plane2_nhits;

            std::vector<double> m_reco_shower_plane0_meanRMS; /* the mean of RMS of the shower hit shape (in tick unit) on plane 0 */
            std::vector<double> m_reco_shower_plane1_meanRMS;
            std::vector<double> m_reco_shower_plane2_meanRMS;

            std::vector<int> m_reco_shower_hit_wire;
            std::vector<int> m_reco_shower_hit_plane;
            std::vector<double> m_reco_shower_hit_tick;
            std::vector<double> m_reco_shower_spacepoint_x;
            std::vector<double> m_reco_shower_spacepoint_z;
            std::vector<double> m_reco_shower_spacepoint_y;


            std::vector<size_t>  m_reco_shower_ordered_energy_index; /* indices of 'm_reco_shower_energy_max' such that energy max is in descending order */
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane0; //for each shower, looks at the hits for all clusters in the plane, stores the dQ/dx for each hit 
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane1;
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane2;
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane0; //dE/dx from the calculated dQ/dx for each hit of all clusters in shower on plane 	
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane1;
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane2;

            std::vector<double> m_reco_shower_dEdx_plane0_mean; /* mean of dE/dx of each hit in shower */
            std::vector<double> m_reco_shower_dEdx_plane1_mean;
            std::vector<double> m_reco_shower_dEdx_plane2_mean;
            std::vector<double> m_reco_shower_dEdx_plane0_max;
            std::vector<double> m_reco_shower_dEdx_plane1_max;
            std::vector<double> m_reco_shower_dEdx_plane2_max;
            std::vector<double> m_reco_shower_dEdx_plane0_min;
            std::vector<double> m_reco_shower_dEdx_plane1_min;
            std::vector<double> m_reco_shower_dEdx_plane2_min;
            std::vector<double> m_reco_shower_dEdx_plane0_median;/* median of dE/dx of each hit in shower (median of vector element of m_reco_shower_dEdx_plane0) */
            std::vector<double> m_reco_shower_dEdx_plane1_median;
            std::vector<double> m_reco_shower_dEdx_plane2_median;

            std::vector<double>  m_reco_shower_angle_wrt_wires_plane0; /* angle between shower direction and wire dir on plane, in radian*/
            std::vector<double>  m_reco_shower_angle_wrt_wires_plane1;
            std::vector<double>  m_reco_shower_angle_wrt_wires_plane2;

            std::vector<double>  m_reco_shower_dEdx_amalgamated;
            std::vector<int>  m_reco_shower_dEdx_amalgamated_nhits;

            std::vector<double> m_reco_shower_dQdx_plane0_median;/* median of dQ/dx of each hit in shower (median of m_reco_shower_dQdx_plane0) */
            std::vector<double> m_reco_shower_dQdx_plane1_median;
            std::vector<double> m_reco_shower_dQdx_plane2_median;

            std::vector<double> m_reco_shower_dEdx_plane0_nhits; /* number of hits of all clusters of the shower on plane 0 */
            std::vector<double> m_reco_shower_dEdx_plane1_nhits;
            std::vector<double> m_reco_shower_dEdx_plane2_nhits;

            double _time2cm;//value modeled from David's shower code

            // PID-related variables
            std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane0;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane1;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane2;
            std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane0;
            std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane1;
            std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane2;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane0;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane1;
            std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane2;
            std::vector<double> m_reco_track_pid_pida_plane0;
            std::vector<double> m_reco_track_pid_pida_plane1;
            std::vector<double> m_reco_track_pid_pida_plane2;
            std::vector<double> m_reco_track_pid_chi2_mu_plane0;
            std::vector<double> m_reco_track_pid_chi2_mu_plane1;
            std::vector<double> m_reco_track_pid_chi2_mu_plane2;
            std::vector<double> m_reco_track_pid_chi2_p_plane0;
            std::vector<double> m_reco_track_pid_chi2_p_plane1;
            std::vector<double> m_reco_track_pid_chi2_p_plane2;
            std::vector<double> m_reco_track_pid_three_plane_proton_pid;


            //Geant4
            std::vector<int> m_geant4_pdg;
            std::vector<int>          m_geant4_trackid;
            std::vector<int>          m_geant4_mother;
            std::vector<int>         m_geant4_statuscode;
            std::vector<double>          m_geant4_E;
            std::vector<double>          m_geant4_mass;
            std::vector<double>          m_geant4_px;
            std::vector<double>          m_geant4_py;
            std::vector<double>          m_geant4_pz;
            std::vector<double>          m_geant4_vx;
            std::vector<double>          m_geant4_vy;
            std::vector<double>          m_geant4_vz;
            std::vector<double>          m_geant4_dx;
            std::vector<double>          m_geant4_dy;
            std::vector<double>          m_geant4_dz;
            std::vector<std::string>          m_geant4_process;
            std::vector<std::string>          m_geant4_end_process;

            std::vector<double>          m_geant4_costheta;




            double m_genie_spline_weight;
            double m_genie_CV_tune_weight;

            double m_photonu_weight_low;
            double m_photonu_weight_high;

            bool Pi0PreselectionFilter(); /* returns whether the event pass pi0 pre-selection for 2g1p */
            bool Pi0PreselectionFilter2g0p(); /* returns whether the event pass pre-selction for 2g0p */
    };

    DEFINE_ART_MODULE(SinglePhoton)

} // namespace lar_pandora
#endif
//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
