#include <vector>

namespace single_photon
{
    //ask YJ
    int setTPCGeom();


  /* inside TPC or not? */
    int isInTPCActive(std::vector<double> & vec);


    /* returns minimum distance to the TPC active boundary; returns -999 if the point is not in TPC active volume */
    double distToTPCActive(std::vector<double>&vec);


    /* returns minimum distance to the TPCActive boundary around the Cathode Plane Assemble; returns -999 if the point is not in TPC active volume */
    double distToCPA(std::vector<double>&vec);



    int distToSCB(double & dist, std::vector<double> &vec);

}
