/**
 * \file QWeightPoint.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class QWeightPoint
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_QWEIGHTPOINT_H
#define OPT0FINDER_QWEIGHTPOINT_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseFlashFilter.h"
#include "flashmatch/Base/FlashMatchFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseFlashFilter.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashMatchFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif

#include <cmath>
#include <sstream>
#include <numeric>


namespace flashmatch {
  
  /**
     \class QWeightPoint
     Implementation of flashmatch::BaseFlashHypothesis algorithm class. \n
     Given a TPC object (flashmatch::QCluster_t), it calcultes a list of flash hypothesis \n
     points based on charge deposition and its geometrical position. Each energy deposition \n
     point is weighted by its charge and inverse-squared-x position. As the absolute \n
     x-position is not known by a TPC object, it uses a relative position for each point \n
     w.r.t. the closest point to the wire plane (x=0). The algorithm then assigns an overall \n
     absolute x-position offset in a successive step of _x_step_size value, assigned by a user, \n
     to compute possible flash hypothesis points.\n
  */
  class QWeightPoint : public BaseFlashMatch {
    
  public:
    
    /// Default constructor
    QWeightPoint(const std::string name="QWeightPoint");
    
    /// Default destructor
    ~QWeightPoint(){}

    FlashMatch_t Match(const QCluster_t&, const Flash_t&);

    void SetTPCCryo(int tpc, int cryo) {};

  protected:

    void _Configure_(const Config_t &pset);

  private:
    double _x_step_size; ///< step size in x-direction
    double _zdiff_max;   ///< allowed diff in z-direction to be considered as a match
    flashmatch::QCluster_t _tpc_qcluster;
    flashmatch::Flash_t    _vis_array;
  };

  /**
     \class flashmatch::QWeightPointFactory
  */
  class QWeightPointFactory : public FlashMatchFactoryBase {
  public:
    /// ctor
    QWeightPointFactory() { FlashMatchFactory::get().add_factory("QWeightPoint",this); }
    /// dtor
    ~QWeightPointFactory() {}
    /// creation method
    BaseFlashMatch* create(const std::string instance_name) { return new QWeightPoint(instance_name); }
  };
  
}
#endif
/** @} */ // end of doxygen group 

