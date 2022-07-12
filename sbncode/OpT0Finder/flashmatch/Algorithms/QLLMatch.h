/**
 * \file QLLMatch.h
 *
 * \ingroup Algorithms
 *
 * \brief Class def header for a class QLLMatch
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef QLLMATCH_H
#define QLLMATCH_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseFlashMatch.h"
#include "flashmatch/Base/FlashMatchFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseFlashMatch.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashMatchFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif


#include <iostream>
#include <TMinuit.h>
#include <TF2.h>
namespace flashmatch {
  /**
     \class QLLMatch
     User defined class QLLMatch ... these comments are used to generate
     doxygen documentation!
  */
  class QLLMatch : public BaseFlashMatch {

  public:

    enum QLLMode_t { kChi2, kLLHD, kSimpleLLHD, kWeightedLLHD, kIntegralLLHD, kZIP, kPEWeightedLLHD };

  private:
    /// Valid ctor hidden (singleton)
    QLLMatch(const std::string);

  public:

    /// Default ctor throws exception (singleton)
    QLLMatch();

    /// Default destructor
    ~QLLMatch(){}

    /// Singleton shared instance getter
    static QLLMatch* GetME(std::string name="")
    {
      if(!_me) _me = new QLLMatch(name);
      else if(!name.empty() && name != _me->AlgorithmName()) {
	std::cerr << "QLLMatch instance must be uniquely named. Requested: "
		  << name << " vs. Existing: " << _me->AlgorithmName() << std::endl;
	throw std::exception();
      }
      return _me;
    }

    /// Core function: execute matching
    void Match(const QCluster_t&, const Flash_t&, FlashMatch_t& match);

    const Flash_t& ChargeHypothesis(const double);
    const Flash_t& Measurement() const;

    double QLL(const flashmatch::Flash_t&,
	       const flashmatch::Flash_t&);

    void Record(const double x)
    {
      if(_record) {
      	_minimizer_record_chi2_v.push_back(_current_chi2);
      	_minimizer_record_llhd_v.push_back(_current_llhd);
      	_minimizer_record_x_v.push_back(x);
        _minimizer_record_pe_v.push_back(_current_pe);
      }
    }

    void OneStep(const double x) {
        _num_steps = _num_steps + 1;
        if (x < _minimizer_min_x) _minimizer_min_x = x;
        if (x > _minimizer_max_x) _minimizer_max_x = x;
    }

    double CallMinuit(const Flash_t& pmt, const double x0);
      
    const std::vector<double>& HistoryLLHD() const { return _minimizer_record_llhd_v; }
    const std::vector<double>& HistoryChi2() const { return _minimizer_record_chi2_v; }
    const std::vector<double>& HistoryX()    const { return _minimizer_record_x_v;    }
    const std::vector<double>& HistoryPE()   const { return _minimizer_record_pe_v;   }

    void DumpHistory() const;

  protected:

    void _Configure_(const Config_t &pset);

  private:

    //FlashMatch_t TouchingTrack(const QCluster_t &pt_v, const Flash_t & flash, double score, bool tpc0);
    void PESpectrumMatch(const Flash_t &flash, const double x0, FlashMatch_t& match);
    std::vector<double> CalculateX0(const Flash_t &pmt);
    void OnePMTMatch(const Flash_t &flash,FlashMatch_t& match);

    static QLLMatch* _me;

    QLLMode_t _mode;   ///< Minimizer mode
    bool _record;      ///< Boolean switch to record minimizer history
    double _normalize; ///< Noramalize hypothesis PE spectrum
    bool _check_touching_track; ///< Whether to match immediately touching track with flash if timing coincides.
    double _touching_track_window; ///< Time(us) such that we use this tolerance T to find touching tracks

    std::vector<double>  _penalty_threshold_v;
    std::vector<double>  _penalty_value_v;
    double _pe_hypothesis_threshold;
    double _pe_observation_threshold;

    flashmatch::QCluster_t _raw_trk;
    QPoint_t _raw_xmin_pt;
    QPoint_t _raw_xmax_pt;
    flashmatch::QCluster_t _var_trk;
    flashmatch::Flash_t    _hypothesis;  ///< Hypothesis PE distribution over PMTs
    flashmatch::Flash_t    _measurement; ///< Flash PE distribution over PMTs

    double _current_pe;
    double _current_chi2;
    double _current_llhd;
    std::vector<double> _minimizer_record_chi2_v; ///< Minimizer record chi2 value
    std::vector<double> _minimizer_record_llhd_v; ///< Minimizer record llhd value
    std::vector<double> _minimizer_record_x_v;    ///< Minimizer record X values
    std::vector<double> _minimizer_record_pe_v;   ///< Minimizer record PE values
    double _minimizer_min_x;
    double _minimizer_max_x;

    double _reco_x_offset;     ///< reconstructed X offset (from wire-plane to min-x point)
    double _reco_x_offset_err; ///< reconstructed X offset w/ error
    double _qll;               ///< Minimizer return value

    double _minuit_x_buffer; ///< a buffer along x (drift) direction for the range in which minuit runs

    TF2 _poisson;

    bool _converged;
    TMinuit* _minuit_ptr;
    double _migrad_tolerance;
    int _num_steps;
    double _offset;
		double _time_shift;

    double _recox_penalty_threshold;
    double _recoz_penalty_threshold;

    double _onepmt_score_threshold;
    double _onepmt_xdiff_threshold;
    double _onepmt_pesum_threshold;
    double _onepmt_pefrac_threshold;

    double _vol_xmax, _vol_xmin;
    std::vector<double> _xpos_v, _ypos_v, _zpos_v;
    std::vector<double> _exp_frac_v, _exp_tau_v;
  };

  /**
     \class flashmatch::QLLMatchFactory
  */
  class QLLMatchFactory : public FlashMatchFactoryBase {
  public:
    /// ctor
    QLLMatchFactory() { FlashMatchFactory::get().add_factory("QLLMatch",this); }
    /// dtor
    ~QLLMatchFactory() {}
    /// creation method
    BaseFlashMatch* create(const std::string instance_name) { return QLLMatch::GetME(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
