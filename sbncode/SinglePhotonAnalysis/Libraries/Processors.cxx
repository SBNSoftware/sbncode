#include "sbncode/SinglePhotonAnalysis/Libraries/variables.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/Processors.h"
#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_math.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"


#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"

//May move this to Libaries
namespace single_photon
{


  int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input){
    corrected.resize(3);

    //CHECK        double kx = input[0];
    double ky = input[1];
    double kz = input[2];
    //CHECK
    //        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
    // CHECK
    // double g4Ticks =  detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->trigger_offset();

    //CHECK       double xtimeoffset = 0;//CHECK  theDetector->ConvertTicksToX(g4Ticks,0,0,0);

    //        double xOffset = -scecorr.X() +xtimeoffset+0.6;
    double yOffset = 0;//CHECK scecorr.Y();
    double zOffset = 0;//CHECK scecorr.Z();

    corrected[0]=0;//CHECK kx - scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
    corrected[1]=ky+yOffset;
    corrected[2]=kz+zOffset;
    return 0;
  }





  int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected){
    corrected.resize(3);

    double kx = mcparticle->Vx();
    double ky = mcparticle->Vy();
    double kz = mcparticle->Vz();

    //CHECK        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));  // to get position offsets to be used in ionization electron drift
    //CHECK         double g4Ticks =  detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->trigger_offset();

    //CHECK        double xtimeoffset =  theDetector->ConvertTicksToX(g4Ticks,0,0,0);

    //double xOffset = -scecorr.X() +xtimeoffset+0.6;
    double yOffset = 0;//CHECK scecorr.Y();
    double zOffset = 0;//CHECK scecorr.Z();

    corrected[0]= kx;//CHECK- scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
    corrected[1]=ky+yOffset;
    corrected[2]=kz+zOffset;

    return 0;
  }





  int spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected){
    corrected.resize(3);
    //Space Charge Effect! functionize this soon.
    double kx = mcparticle.Vx();
    double ky = mcparticle.Vy();
    double kz = mcparticle.Vz();
    //CHECK         auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
    //CHECK        double g4Ticks =  detClocks->TPCG4Time2Tick(mcparticle.T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->trigger_offset();

    //CHECK        double xtimeoffset = 0;//CHECK theDetector->ConvertTicksToX(g4Ticks,0,0,0);

    corrected[0]=kx;//CHECK  - scecorr.X() +xtimeoffset+0.6;
    corrected[1]=ky;//CHECK  + scecorr.Y();
    corrected[2]=kz;//CHECK  + scecorr.Z();
    return 0;
  }

  void CollectMCParticles(
      const art::Event &evt,
      const std::string &label,
      std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,
      std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> &particlesToTruth,
      std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap,
	  var_all& vars){

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
    {
      const art::Ptr<simb::MCParticle> particle(theParticles, i);
      const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
      truthToParticles[truth].push_back(particle);
      particlesToTruth[particle] = truth;
      MCParticleToTrackIdMap[particle->TrackId()] = particle;
    }

    //        std::cout<<"CollectMCParticles() \t||\t the number of MCParticles in the event is "<<theParticles->size()<<std::endl;
  }

  void CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector)
  {
        if (evt.isRealData())
          throw cet::exception("LArPandora") << " PandoraCollector::CollectSimChannels --- Trying to access MC truth from real data ";

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


  void BuildMCParticleHitMaps(
      const art::Event &evt, 
      const std::string &label, 
      const std::vector<art::Ptr<recob::Hit>> &hitVector,   
      std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         
      std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  &hitsToParticles, 
      const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, 
      std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap,   
      var_all& vars){
    std::vector< art::Ptr<sim::SimChannel> >   simChannelVector;
    std::map< art::Ptr<simb::MCTruth>,     std::vector<art::Ptr<simb::MCParticle>>  >    truthToParticles;
    std::map< art::Ptr<simb::MCParticle>,  art::Ptr<simb::MCTruth> > particlesToTruth;
    std::map< art::Ptr<recob::Hit>,    std::vector< sim::TrackIDE >    >               hitsToTrackIDEs;

    CollectSimChannels(evt, label, simChannelVector);
    CollectMCParticles(evt, label, truthToParticles, particlesToTruth, MCParticleToTrackIdMap, vars);

    //Collect the links from reconstructed hits to their true energy deposits.
    lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(evt, hitVector, simChannelVector, hitsToTrackIDEs);
    //Build mapping between Hits and MCParticles, starting from Hit/TrackIDE/MCParticle information
    lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(hitsToTrackIDEs, truthToParticles, particlesToHits, hitsToParticles, daughterMode);


  }

  bool Pi0PreselectionFilter(var_all& vars, para_all& paras){

    if(vars.m_vertex_pos_x < paras.s_tpc_active_XMin||  vars.m_vertex_pos_x > paras.s_tpc_active_XMax) return false;
    if(vars.m_vertex_pos_y < paras.s_tpc_active_YMin || vars.m_vertex_pos_y > paras.s_tpc_active_YMax) return false;
    if(vars.m_vertex_pos_z < paras.s_tpc_active_ZMin||  vars.m_vertex_pos_z > paras.s_tpc_active_ZMax) return false;

    if(vars.m_reco_asso_showers!=2) return false;
    if(vars.m_reco_asso_tracks!=1) return false;
    if(vars.m_reco_vertex_size<1) return false;

    if(vars.m_reco_shower_conversion_distance.size()!=2) return false;
    if(vars.m_reco_shower_conversion_distance[0]<1. || vars.m_reco_shower_conversion_distance[1]<1.) return false;

    return true;
  }



  bool Pi0PreselectionFilter2g0p(var_all& vars, para_all& paras){
    if(vars.m_vertex_pos_x < paras.s_tpc_active_XMin||  vars.m_vertex_pos_x > paras.s_tpc_active_XMax) return false;
    if(vars.m_vertex_pos_y < paras.s_tpc_active_YMin || vars.m_vertex_pos_y > paras.s_tpc_active_YMax) return false;
    if(vars.m_vertex_pos_z < paras.s_tpc_active_ZMin||  vars.m_vertex_pos_z > paras.s_tpc_active_ZMax) return false;

    if(vars.m_reco_asso_showers!=2) return false;
    if(vars.m_reco_asso_tracks!=0) return false;
    if(vars.m_reco_vertex_size<1) return false;

    if(vars.m_reco_shower_energy_max.size()!=2) return false;
    //if the maximum energy of all showers on all planes is smaller than 30
    if(vars.m_reco_shower_energy_max[vars.m_reco_shower_ordered_energy_index[0]]<30.) return false;

    return true;
  }

  bool IsEventInList(int run, int subrun, int event, var_all& vars){
    if(vars.m_selected_set.find( {run, subrun, event} ) == vars.m_selected_set.end()){
      if(vars.m_selected_set.find({run, subrun})  == vars.m_selected_set.end() ){
        if(vars.m_selected_set.find({run}) == vars.m_selected_set.end())
          return false;
      }
    }
    return true;
  }

  //----------- Above are migrated from Singlephoton_module.cc

  //determines if a point is inside the rectangle by summing the areas of the four triangles made by 
  //if the point is inside, the sum of the triangles should exactly equal the area of the rectangle
  //also returns true if the point is on the boundary
  bool isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle, para_all& paras){
    int n_vertices = (int)rectangle.size();
    //bool inside = false;
    int i, j = 0;
    double areas = 0;
    //for each pair of vertices
    for (i = 0, j = n_vertices-1; i < n_vertices; j = i++) {
      //calculate the area of a triangle with the point and two vertices
      double this_area = 0.5*abs(rectangle[i][0]*(rectangle[j][1] - thishit_pos[1]) 
          + rectangle[j][0]*(thishit_pos[1] - rectangle[i][1]) 
          + thishit_pos[0]*(rectangle[i][1] - rectangle[j][1]));
      areas += this_area;
    }        
    //calc area of the rectangle
    double area_rectangle = paras.s_width_dqdx_box* paras.s_length_dqdx_box;

    //check the sum of areas match
    if (abs(areas - area_rectangle) <= 0.001 ){
      return true;
    }
    return false;
  }


  //helpers for calculating calometry
//  double CalcEShower(const std::vector<art::Ptr<recob::Hit>> &hits, var_all& vars){
//    double energy[3] = {0., 0., 0.};
//
//    //std::cout<<"AnalyzeShowers() \t||\t Looking at shower with "<<hits.size() <<" hits on all planes"<<std::endl;
//
//    //for each hit in the shower
//    for (auto &thishitptr : hits){
//      //check the plane
//      int plane= thishitptr->View();
//
//      //skip invalid planes       
//      if (plane > 2 || plane < 0)  continue;
//
//      //calc the energy of the hit
//      double E = QtoEConversion(GetQHit(thishitptr, plane, vars), vars);
//
//      //add the energy to the plane
//      energy[plane] += E;
//    }//for each hiti
//
//    //find the max energy on a single plane
//    double max = energy[0];
//    for (double en: energy){
//      if( en > max){
//        max = en;
//      }
//    }
//    // std::cout<<"AnalyzeShowers() \t||\t The energy on each plane for this shower is "<<energy[0]<<", "<<energy[1]<<", "<<energy[2]<<std::endl;
//
//    //return the highest energy on any of the planes
//    return max;
//
//  }

  double CalcEShowerPlane(const std::vector<art::Ptr<recob::Hit>>& hits, int this_plane, para_all& paras){
    double energy = 0.;

    //for each hit in the shower
    for (auto &thishitptr : hits){
      //check the plane
      int plane= thishitptr->View();

      //skip invalid planes       
      if (plane != this_plane )  continue;

      //calc the energy of the hit
      double E = QtoEConversion(GetQHit(thishitptr, plane, paras), paras);
      //add the energy to the plane
      energy += E;
    }//for each hit

    return energy;

  }




  double GetQHit(art::Ptr<recob::Hit> thishitptr, int plane, para_all& paras){
    double gain;
    //choose gain based on whether data/mc and by plane
    if (paras.s_is_data == false &&  paras.s_is_overlayed == false){
      gain = paras.s_gain_mc[plane] ;
      //if (g_is_verbose) std::cout<<"the gain for mc on plane "<<plane<<" is "<<gain<<std::endl;
    } if (paras.s_is_data == true ||  paras.s_is_overlayed == true){
      gain = paras.s_gain_data[plane] ;
      //if (g_is_verbose) std::cout<<"the gain for data on plane "<<plane<<" is "<<gain<<std::endl;

    }

    double Q = thishitptr->Integral()*gain;
    return Q;
  }


  double QtoEConversion(double Q, para_all& paras){
    //return the energy value converted to MeV (the factor of 1e-6)
    double E = Q* paras.s_work_function *1e-6 /paras.s_recombination_factor;
    return E;

  }


  std::vector<double> CalcdEdxFromdQdx(std::vector<double> dqdx, para_all& paras){
    int n = dqdx.size();
    std::vector<double> dedx(n,0.0);
    for (int i = 0; i < n; i++){
      //std::cout<<"The dQ/dx is "<<dqdx[i]<<std::endl;
      dedx[i] = QtoEConversion(dqdx[i], paras);
      //std::cout<<"The dE/dx is "<<dedx[i]<<std::endl;
    }
    return dedx;
  }


  std::vector<double> CalcdQdxShower(
      const art::Ptr<recob::Shower>& shower,
      const std::vector<art::Ptr<recob::Cluster>> & clusters, 
      std::map<art::Ptr<recob::Cluster>,    std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,
      int plane,
      double triggeroffset,
      detinfo::DetectorPropertiesData const & theDetector,
      para_all & paras){

    //if(g_is_verbose) std::cout<<"AnalyzeShowers() \t||\t The number of clusters in this shower is "<<clusters.size()<<std::endl;
    std::vector<double> dqdx;

    //get the 3D shower direction
    //note: in previous versions of the code there was a check to fix cases where the shower direction was inverted - this hasn't been implemented
    TVector3 shower_dir(shower->Direction().X(), shower->Direction().Y(),shower->Direction().Z());

    //calculate the pitch for this plane
    //CHECK upgradable, see ./Calibration/TrackCaloSkimmer_module.cc line 746
    double pitch = getPitch(shower_dir, paras)[plane];  

    if(g_is_verbose) std::cout<<"AnalyzeShowers() \t||\t The pitch between the shower and plane "<<plane<<" is "<<pitch<<std::endl;

    //for all the clusters in the shower
    for (const art::Ptr<recob::Cluster> &thiscluster: clusters){
      //keep only clusters on the plane
      if(thiscluster->View() != plane) continue;

      //calculate the cluster direction
      std::vector<double> cluster_axis = {cos(thiscluster->StartAngle()), sin(thiscluster->StartAngle())};    

      //get the cluster start and and in CM
      //std::cout<<"for plane/tpc/cryo:"<<plane<<"/"<<paras.s_TPC<<"/"<<paras.s_Cryostat<<", fXTicksOffset: "<<theDetector->GetXTicksOffset(plane, paras.s_TPC, paras.s_Cryostat)<<" fXTicksCoefficient: "<<theDetector->GetXTicksCoefficient(paras.s_TPC, paras.s_Cryostat)<<std::endl;

      //convert the cluster start and end positions to time and wire coordinates
      std::vector<double> cluster_start = {thiscluster->StartWire() * paras.s_wire_spacing,(thiscluster->StartTick() - triggeroffset)* paras._time2cm};
      std::vector<double> cluster_end = {thiscluster->EndWire() * paras.s_wire_spacing,(thiscluster->EndTick() - triggeroffset)* paras._time2cm };

      //check that the cluster has non-zero length
      double length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) + pow(cluster_end[1] - cluster_start[1], 2));
      //if(g_is_verbose) std::cout<<"AnalyzeShowers() \t||\t The cluster length is "<<length<<std::endl;
      if (length <= 0){ 
        std::cout<<"skipping cluster on plane "<<plane<<", length = "<<length<<std::endl;
        continue;
      }


      //draw a rectangle around the cluster axis 
      std::vector<std::vector<double>> rectangle = buildRectangle(cluster_start, cluster_axis, paras.s_width_dqdx_box, paras.s_length_dqdx_box);  

      //get all the hits for this cluster
      std::vector<art::Ptr<recob::Hit>> hits =  clusterToHitMap[thiscluster];

      //for each hit in the cluster
      for (art::Ptr<recob::Hit> &thishit: hits){  
        //get the hit position in cm from the wire and time
        std::vector<double> thishit_pos ={thishit->WireID().Wire * paras.s_wire_spacing, (thishit->PeakTime() - triggeroffset)* paras._time2cm};

        //check if inside the box
        bool v2 = isInsidev2(thishit_pos, rectangle, paras);
        if (v2){
          double q = GetQHit(thishit, plane, paras); 
          double this_dqdx = q/pitch; 
          dqdx.push_back(this_dqdx);
        }//if hit falls inside the box

      }//for each hit inthe cluster
    }//for each cluster
    return dqdx;
  }

  std::vector<double> getPitch(TVector3 shower_dir, para_all& paras){
		
	std::vector<double> pitches;
        for (geo::PlaneGeo const& plane: paras.s_wireReadout->Iterate<geo::PlaneGeo>()) {
        //6 planes in SBND
        //WireAngleToVertical  : 30 ,150,90,150,30 ,90
        //ub wire angles    : 30 ,150,90  (respected to beam,z)
        //Pitch        : 0.3,0.3,0.3,0.3,0.3,0.3

        const double angToVert(paras.s_wireReadout->WireAngleToVertical(plane.View(), plane.ID())+0.5*M_PI);//wire angle respected to z + pi/2

        TVector3 wire_vector;  
        if(abs(angToVert) < 1e-9 ){
          wire_vector = {0,0,1};
        } else{
          wire_vector = { 0 , sin(angToVert) ,cos(angToVert) };
        }

        //    std::cout<<" Angle "<<angToVert<<" Get Vec y="<<wire_vector[1]<< " z= "<<wire_vector[2]<<std::endl;
		double cos = abs(wire_vector.Dot(shower_dir))/(wire_vector.Mag()*shower_dir.Mag());

		pitches.push_back((cos==0)? std::numeric_limits<double>::max() : paras.s_wire_spacing/cos );

        if(pitches.size()==3) break;
      }
	  return pitches;
  }



  double getMeanHitWidthPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane){
    int nhits = 0;
    double widths = 0;
    for (art::Ptr<recob::Hit> thishitptr : hits){
      //check the plane
      int plane= thishitptr->View();

      //skip invalid planes       
      if (plane != this_plane) continue;

      widths += thishitptr->RMS(); // recob::Hit->RMS() returns RMS of the hit shape in tick units
      nhits++;
    }//for each hiti
    return   widths/(double)nhits;
  }



  int getNHitsPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane){
    int nhits = 0;
    for (art::Ptr<recob::Hit> thishitptr : hits){
      //check the plane
      int plane= thishitptr->View();

      //skip invalid planes       
      if (plane != this_plane) continue;

      nhits++;
    }//for each hiti
    return nhits;
  }



  double triangle_area(double a1, double a2, double b1, double b2, double c1, double c2){
    double m1 = 0.3;
    double m2 = 1.0/25.0;

    return fabs((a1*m1*(b2*m2-c2*m2)+b1*m1*(c2*m2-a2*m2)+c1*m1*(a2*m2-b2*m2))/2.0);
  }

  int quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area){

    std::vector<double> z(n,0.0);

    TGraph2D *g = new TGraph2D(n,X,Y,&z[0]);
    TGraphDelaunay delan(g);
    delan.SetMarginBinsContent(0);
    delan.ComputeZ(0,0);
    delan.FindAllTriangles();
    (*num_triangles)=delan.GetNdt(); // number of Delaunay triangles found

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

  int delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area, para_all& paras){

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
          C0.push_back((double)hit->WireID().Wire);         
          T0.push_back(hit->PeakTime());         
          n_0++;
          break;
        case 1:
          C1.push_back((double)hit->WireID().Wire);         
          T1.push_back(hit->PeakTime());         
          n_1++;
          break;
        case 2:
          C2.push_back((double)hit->WireID().Wire);         
          T2.push_back(hit->PeakTime());         
          n_2++;
          break;
        default:
          break;
      }
    }
    if(paras.s_use_delaunay){
      if(n_0>0 && (int)n_0 < paras.s_delaunay_max_hits) quick_delaunay_fit(n_0, &C0[0]  , &T0[0]  , &num_triangles[0], &area[0]);
      if(n_1>0 && (int)n_1 < paras.s_delaunay_max_hits) quick_delaunay_fit(n_1, &C1[0]  , &T1[0]  , &num_triangles[1], &area[1]);
      if(n_2>0 && (int)n_2 < paras.s_delaunay_max_hits) quick_delaunay_fit(n_2, &C2[0]  , &T2[0]  , &num_triangles[2], &area[2]);
    }
    num_hits[0] = n_0;
    num_hits[1] = n_1;
    num_hits[2] = n_2;

    //std::cout<<"Plane 0: "<<n_0<<" hits with "<<num_triangles[0]<<" triangles of area: "<< area[0]<<std::endl;
    //std::cout<<"Plane 1: "<<n_1<<" hits with "<<num_triangles[1]<<" triangles of area: "<< area[1]<<std::endl;
    //std::cout<<"Plane 2: "<<n_2<<" hits with "<<num_triangles[2]<<" triangles of area: "<< area[2]<<std::endl;

    return 0;
  }
}
