/**
 * \file NPtFilter.h
 *
 * \ingroup Algorithms
 *
 * \brief Class def header for a class NPtFilter
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_NPTFILTER_H
#define OPT0FINDER_NPTFILTER_H

#include "sbncode/OpT0Finder/flashmatch/Base/BaseTPCFilter.h"
#include "sbncode/OpT0Finder/flashmatch/Base/TPCFilterFactory.h"

namespace flashmatch {
  /**
     \class NPtFilter
     Implementation of flashana::BaseTPCFilter abstract algorithm class. \n
     It applies a _very_ simple filtering: excludes TPC objects (flashana::QCluster_t) \n
     that contains less than specified number of points. \n
  */
  class NPtFilter : public BaseTPCFilter{

  public:

    /// Default constructor
    NPtFilter(const std::string name="NPtFilter");

    /// Default destructor
    ~NPtFilter(){}

    /// Implementation of virtualfunction
    IDArray_t Filter(const QClusterArray_t&);

    /// set minimum number of point in TPC track
    void SetMinNumPoints(size_t n) { _min_num_pt = n; }

  protected:

    void _Configure_(const Config_t &pset);

  private:

    size_t _min_num_pt; ///< mininum number of QPoint_t to pass the filter

  };

  /**
     \class flashana::NPtFilterFactory
  */
  class NPtFilterFactory : public TPCFilterFactoryBase {
  public:
    /// ctor
    NPtFilterFactory() { TPCFilterFactory::get().add_factory("NPtFilter",this); }
    /// dtor
    ~NPtFilterFactory() {}
    /// creation method
    BaseTPCFilter* create(const std::string instance_name) { return new NPtFilter(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group

