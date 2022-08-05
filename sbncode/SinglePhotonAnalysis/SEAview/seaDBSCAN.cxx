#include "seaDBSCAN.h"

namespace seaview{

  std::vector<int> seaDBSCAN::Scan2D(std::vector<std::vector<double>> &pts){

    int cluster_count = 0;//we have no clusters
    size_t N = pts.size();
    int l_undef = -99;
    int l_noise = 0;
    std::vector<int> label(N,l_undef);

    for(size_t i=0; i<N; i++){
      if(label[i]!=l_undef) continue;
      std::vector<std::vector<double>> neighbours = this->GetNeighbours(i,pts,false);
      //std::cout<<i<<" has neightbours "<<neighbours.size()<<std::endl;

      if((int)neighbours.size()+1 < m_minpts){ // if there is less than minpts, its a noise point
        label[i]= l_noise;
        //  std::cout<<i<<" thats less than min, noise"<<std::endl;
        continue;
      }

      cluster_count+=1;
      label[i] = cluster_count;


      std::vector<std::vector<double>> seed_set = neighbours;
      for(size_t q=0; q<seed_set.size(); q++){
        size_t iq = (size_t)seed_set[q][2]; //This is original 

        if(label[iq]==l_noise){
          label[iq] = cluster_count;//Change noise to a border point;
        }

        if(label[iq]!=l_undef){
          continue; // previously processed, already identified to a cluster
        }

        // if label[iq] is l_undef
        label[iq]=cluster_count;// wasn't noise, new point, add to cluster

        std::vector<std::vector<double>> new_neighbours = this->GetNeighbours(iq,pts,true);//Get neighbours of this point, including itslef eh.

        if((int)new_neighbours.size() >= m_minpts ){ //expand the seed set
          //Guanqun: change vector while looping over its elements
          //new elements are pushed back to seed_set
          this->UnionSets(seed_set, new_neighbours); 
        }
      }
    }
    return label;

  }



  std::vector<std::vector<double>> seaDBSCAN::GetNeighbours(size_t ipoint, std::vector<std::vector<double>>&pts,bool include_self){
    std::vector<std::vector<double>> neighbours;
    std::vector<double> point = pts[ipoint];

    //VERY simple, will update soon to a DB 


    for(size_t ip=0; ip<pts.size(); ip++){
      std::vector<double> p = pts[ip];

      double dist = this->SimpleDist( p[0], p[1], point[0], point[1]);

      if(include_self){
        if(dist <= m_eps){
          std::vector<double> tp = {p[0],p[1],(double)ip};//Push back original index too
          neighbours.push_back(tp);
        }
      }else{
        if(dist <= m_eps && p != point ){
          std::vector<double> tp = {p[0],p[1],(double)ip};//Push back original index too
          neighbours.push_back(tp);
        }

      }
    }
    return neighbours;
  }

  int seaDBSCAN::UnionSets(std::vector<std::vector<double>> &seed, std::vector<std::vector<double>> &pts){

    //VERY simple, will update soon if it works
    for(auto &p:pts){

      bool is_in=false;
      for(auto &s:seed){
        if(s==p){
          is_in = true;
          break;
        }
      }

      if(is_in == false){
        seed.push_back(p);
      }

    }



    return 0;
  }

  double seaDBSCAN::SimpleDist(double w1, double t1, double w2, double t2){
    // Guanqun: wire and tick conversion to distance??, if so should only be done on plane 2
    double wire_con = 0.3;
    double tick_con = 1.0/25.0;

    return  sqrt(pow(w1*wire_con-w2*wire_con,2)+pow(t1*tick_con-t2*tick_con,2));

  }


}
