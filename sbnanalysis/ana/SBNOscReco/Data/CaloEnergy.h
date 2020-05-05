#ifndef _sbnumurecodata_CaloEnergy_hh
#define _sbnumurecodata_CaloEnergy_hh
namespace numu {

struct CaloEnergy {
  float best_plane;
  float coll_plane;
  int coll_nhit;
  int best_nhit;

  CaloEnergy():
    best_plane(-1),
    coll_plane(-1),
    coll_nhit(-1),
    best_nhit(-1)
  {}
};

} // end namespace
#endif
