/**
 * \file CommonAmps.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class CommonAmps
 *
 * @author ariana Hackenburg 
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_COMMONAMPS_H
#define OPT0FINDER_COMMONAMPS_H

#include "flashmatch/Base/BaseFlashMatch.h"
#include "flashmatch/Base/FlashMatchFactory.h"
#include "TH1D.h"

namespace flashmatch {
  
  /**
     \class CommonAmps
     Implementation of flashmatch::BaseFlashHypothesis algorithm class. \n
     The goal of this algorithm is to compare the most prominent pieces of opflash spectra, 
     and spectra produced by photon library given QCluster points. 
     Tpc_point calculated as in QWeightPoint and stored based on the best match amplitudes.

     Only works with Photon Library currently
  */
  class CommonAmps : public BaseFlashMatch {
    
  public:
    
    /// Default constructor
    CommonAmps(const std::string name="CommonAmps");
    
    /// Default destructor
    ~CommonAmps(){}

    FlashMatch_t Match(const QCluster_t&, const Flash_t&);

  protected:

    void _Configure_(const Config_t &pset);

  private:

    float _percent;
    float _score ;
    float _x_step_size;
    flashmatch::QCluster_t _tpc_qcluster;
    flashmatch::Flash_t    _vis_array;

  };

  /**
     \class flashmatch::CommonAmpsFactory
  */
  class CommonAmpsFactory : public FlashMatchFactoryBase {
  public:
    /// ctor
    CommonAmpsFactory() { FlashMatchFactory::get().add_factory("CommonAmps",this); }
    /// dtor
    ~CommonAmpsFactory() {}
    /// creation method
    BaseFlashMatch* create(const std::string instance_name) { return new CommonAmps(instance_name); }
  };
}
#endif
/** @} */ // end of doxygen group 

