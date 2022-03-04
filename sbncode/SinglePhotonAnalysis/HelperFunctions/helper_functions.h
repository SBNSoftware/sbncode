
namespace single_photon
{
//-----------------HELPER FUNCTIONS -----------
	//CHECK Copy from bad_channel_matching.h    
	////line between x1 and x2, point x0;
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

	// right, but can be simplified
        double d2 = pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2)+2*t*((x2-x1)*(x1-x0)+(y2-y1)*(y1-y0)+(z2-z1)*(z1-z0))+t*t*( pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));


        return sqrt(d2);

    }
    // minimum distance between point and dead wire
    double distanceToNearestDeadWire(int plane, double Ypoint, double Zpoint,
                const geo::GeometryCore * geom,
                std::vector<std::pair<int,int>> & bad_channel_list_fixed_mcc9 ){

             double min_dist = 999999;
            
             for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
                   int channel = bad_channel_list_fixed_mcc9[i].first;                                       
                   int is_ok = bad_channel_list_fixed_mcc9[i].second;       
                        if(is_ok>1)continue;

                   auto wireids = geom->ChannelToWire(channel); //type of wireids: IDs of all the connected wires to the channel
                   auto result = geom->WireEndPoints(wireids[0]);
                        
                    //std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl; 
//                    std::cout<<wireids[0].Plane<<" "<<result.start().X()<<std::endl;
                    if(plane != (int)wireids[0].Plane) continue;

		    // for the dead wires on the same plane 
                    std::vector<double> start = {0.0,result.start().Y(),result.start().Z()};
                    std::vector<double> end = {0.0,result.end().Y(),result.end().Z()};
                    std::vector<double> point = {0.0,Ypoint,Zpoint};
                    double dist = dist_line_point(start,end,point);
                    min_dist = std::min(dist,min_dist);
             }
                 
    return min_dist;

    }

//--------------- end of copying from bad_channel_matching.h

    double impact_paramater_shr(double x, double y, double z, art::Ptr<recob::Shower> & shr){

        std::vector<double> vert = {x,y,z}; 
        std::vector<double> start = {shr->ShowerStart().X(), shr->ShowerStart().Y(),shr->ShowerStart().Z()};
        std::vector<double> abit = {shr->ShowerStart().X() + shr->Direction().X(),  shr->ShowerStart().Y()+shr->Direction().Y(),  shr->ShowerStart().Z()+shr->Direction().Z()};

        return dist_line_point(start, abit, vert);

    }

    // invariant mass of a particle that decays to two showers
    double  implied_invar_mass(double vx, double vy, double vz, art::Ptr<recob::Shower> & s1, double E1,  art::Ptr<recob::Shower> &s2, double E2){

        double s1x = s1->ShowerStart().X()-vx;
        double s1y = s1->ShowerStart().Y()-vy;
        double s1z = s1->ShowerStart().Z()-vz;
        double norm1  = sqrt(pow(s1x,2)+pow(s1y,2)+pow(s1z,2));
        s1x = s1x/norm1; //unit vector pointing to shower start from point (vx, vy,vz)
        s1y = s1y/norm1;
        s1z = s1z/norm1;

        double s2x = s2->ShowerStart().X()-vx;
        double s2y = s2->ShowerStart().Y()-vy;
        double s2z = s2->ShowerStart().Z()-vz;
        double norm2  = sqrt(pow(s2x,2)+pow(s2y,2)+pow(s2z,2));
        s2x = s2x/norm2; // unit vector pointing to shower start from point (vx, vy, vz)
        s2y = s2y/norm2;
        s2z = s2z/norm2;

        return sqrt(2.0*E1*E2*(1.0-(s1x*s2x+s1y*s2y+s1z*s2z)));


    }

    // invariant mass of two showers, calculated directly from shower directions
    double  invar_mass(art::Ptr<recob::Shower> & s1, double E1,  art::Ptr<recob::Shower> &s2, double E2){

        double s1x = s1->Direction().X();
        double s1y = s1->Direction().Y();
        double s1z = s1->Direction().Z();

        double s2x = s2->Direction().X();
        double s2y = s2->Direction().Y();
        double s2z = s2->Direction().Z();

        return sqrt(2.0*E1*E2*(1.0-(s1x*s2x+s1y*s2y+s1z*s2z)));

    }



    // sort indices in descending order
    template <typename T>
        std::vector<size_t> sort_indexes(const std::vector<T> &v) {

            std::vector<size_t> idx(v.size());
            std::iota(idx.begin(), idx.end(), 0);

            // sort indexes based on comparing values in v
            std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

            return idx;
        }

    // sort indices such that elements in v are in ascending order
    template <typename T>
        std::vector<size_t> sort_indexes_rev(const std::vector<T> &v) {

            std::vector<size_t> idx(v.size());
            std::iota(idx.begin(), idx.end(), 0);

            // sort indexes based on comparing values in v
            std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

            return idx;
        }


    // check if two vectors have same elements (regardless of the order), and arrange their elements in order
    template<typename T>
        bool marks_compare_vec_nonsense(std::vector<T>& v1, std::vector<T>& v2)
        {
            std::sort(v1.begin(), v1.end());
            std::sort(v2.begin(), v2.end());
            return v1 == v2;
        }



    double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
        double wire = geo.WireCoordinate(Y, Z, plane, fTPC, fCryostat);
        return wire;
    }


//    double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorPropertiesData const& detpropdata){
//		
//        double time = detpropdata.ConvertXToTicks(X, plane, fTPC,fCryostat);
//        return time;
//    }


//----------- END OF HELPER FUNCTIONS -----------


}
