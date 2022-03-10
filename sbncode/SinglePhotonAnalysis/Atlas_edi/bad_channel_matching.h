#include "SinglePhoton_module.h"


namespace single_photon
{

    //line between x1 and x2, point x0;
    double dist_line_point( std::vector<double>&X1, std::vector<double>& X2, std::vector<double>& point){
        double x1 =X1.at(0);
        double y1 =X1.at(1);
        double z1 =X1.at(2);

        double x2 =X2.at(0);
        double y2 =X2.at(1);
        double z2 =X2.at(2);

        double x0 =point.at(0);
        double y0 =point.at(1);
        double z0 =point.at(2);

        double x10 = x1-x0;
        double y10 = y1-y0;
        double z10 = z1-z0;

        double x21 = x2-x1;
        double y21 = y2-y1;
        double z21 = z2-z1;

        double t = -(x10*x21+y10*y21+z10*z21)/fabs(x21*x21+y21*y21+z21*z21 );

        double d2 = pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2)+2*t*((x2-x1)*(x1-x0)+(y2-y1)*(y1-y0)+(z2-z1)*(z1-z0))+t*t*( pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));


        return sqrt(d2);

    }

    double distanceToNearestDeadWire(int plane, double Ypoint, double Zpoint,
                const geo::GeometryCore * geom,
                std::vector<std::pair<int,int>> & bad_channel_list_fixed_mcc9 ){

             double min_dist = 999999;
            
             for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
                   int channel = bad_channel_list_fixed_mcc9[i].first;                                       
                   int is_ok = bad_channel_list_fixed_mcc9[i].second;       
                        if(is_ok>1)continue;

                   auto wireids = geom->ChannelToWire(channel);
                   auto result = geom->WireEndPoints(wireids[0]);
                        
                    //std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl; 
//                    std::cout<<wireids[0].Plane<<" "<<result.start().X()<<std::endl;
                    if(plane != (int)wireids[0].Plane) continue;

                    std::vector<double> start = {0.0,result.start().Y(),result.start().Z()};
                    std::vector<double> end = {0.0,result.end().Y(),result.end().Z()};
                    std::vector<double> point = {0.0,Ypoint,Zpoint};
                    double dist = dist_line_point(start,end,point);
                    min_dist = std::min(dist,min_dist);
             }
                 
    return min_dist;

    }


    //Typenamed for recob::Track and recob::Shower
    template<typename T>
        int badChannelMatching(std::vector<int>& badchannels, 
                std::vector<T>& objects, 
                std::map< T, art::Ptr<recob::PFParticle> > & objectToPFParticleMap,
                std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,
                const geo::GeometryCore * geom,
                std::vector<std::pair<int,int>> & bad_channel_list_fixed_mcc9){




            for(size_t i_trk=0; i_trk<objects.size(); i_trk++){


                const T object = objects[i_trk];
                const art::Ptr<recob::PFParticle> pfp = objectToPFParticleMap[object];
                const std::vector<art::Ptr<recob::Hit>> hits = pfParticleToHitsMap[pfp];

                int min_dist_from_bad_channel = 99999;

                for(size_t h=0; h<hits.size(); h++){
                const art::Ptr<recob::Hit> hit = hits[h];
                    
            
                    //int nch = (int)badchannels.size()/3;
                    //for(int i=0; i<nch; i++){
                    for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
                      //  const int offset = 3*i;
                      //  int bc = badchannels[offset+0];                                       
                        int bc = bad_channel_list_fixed_mcc9[i].first;                                       
                        int ok = bad_channel_list_fixed_mcc9[i].second;       
                        if(ok>1)continue;
                        int dist =hit->Channel()-bc;
                        auto hs = geom->ChannelToWire(bc);
                        //std::cout<<"AG: "<<hs.size()<<"  BC("<<bc<<"): "<<hs[0]<<" ours: ("<<hit->Channel()<<"): "<<hit->WireID()<<std::endl;
                        //this is the right format for my plotting routine
                        //std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<< badchannels[offset+1]<<" "<<badchannels[offset+2]<<std::endl;
                        std::vector<double> start(3);
                        std::vector<double> end(3);
                        auto result = geom->WireEndPoints(hs[0]);
                        
                //        std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl; 

                        

                        if(fabs(dist) < min_dist_from_bad_channel) min_dist_from_bad_channel = fabs(dist);
                   
            //            std::cout<<"BD: "<<badchannels[offset+0]<<" "<<badchannels[offset+1]<<" "<<badchannels[offset+2]<<" "<< fabs(badchannels[offset+1]-badchannels[offset+2])<<" "<<dist<<std::endl;
                    }
                }//hitloop
            }

            return 0;
        }

}//namespace end
