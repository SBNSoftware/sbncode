/**
 * \file DBSCAN.h
 *
 * 
 * \brief Class def header for a class DBSCAN
 *
 * @author mark ross-lonergan markrl@nevis.columbia.edu
 * Written 20th May 2019.
 */

#ifndef DBSCAN_H
#define DBSCAN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <climits>
#include <limits>


class DBSCAN{

    public:
        double m_eps;
        int m_minpts;

        /// Default constructor
        DBSCAN(double in_eps, int in_minpts): m_eps(in_eps), m_minpts(in_minpts) {}

        /// Default destructor
        ~DBSCAN(){}

        std::vector<int> Scan2D(std::vector<std::vector<double>> &pts);
        std::vector<std::vector<double>> GetNeighbours(size_t i, std::vector<std::vector<double>> &pts,bool);
        int UnionSets(std::vector<std::vector<double>> &seed, std::vector<std::vector<double>> &pts);
};

std::vector<int> DBSCAN::Scan2D(std::vector<std::vector<double>> &pts){

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
                continue; // previously processed
            }

            label[iq]=cluster_count;// wasn't noise, new point, add to cluster

            std::vector<std::vector<double>> new_neighbours = this->GetNeighbours(iq,pts,true);//Get neighbours of this point, including itslef eh.

            if((int)new_neighbours.size() >= m_minpts ){ //expand the seed set
                this->UnionSets(seed_set, new_neighbours); 
            }
        }
    }
    return label;

}



std::vector<std::vector<double>> DBSCAN::GetNeighbours(size_t ipoint, std::vector<std::vector<double>>&pts,bool include_self){
    std::vector<std::vector<double>> neighbours;
    std::vector<double> point = pts[ipoint];

    //VERY simple, will update soon to a DB 


    for(size_t ip=0; ip<pts.size(); ip++){
        std::vector<double> p = pts[ip];
        
        double dist = sqrt(pow(p[0]*0.3-point[0]*0.3,2)+pow(p[1]/25.0-point[1]/25.0,2));

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

int DBSCAN::UnionSets(std::vector<std::vector<double>> &seed, std::vector<std::vector<double>> &pts){

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


#endif
