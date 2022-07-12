#ifndef QLLMATCH_CXX
#define QLLMATCH_CXX

#include "QLLMatch.h"
#include <TMinuit.h>
#include <cmath>
#include <numeric>
#include <TMath.h>
#include <cassert>
#include <chrono>
#include <climits>
using namespace std::chrono;
namespace flashmatch {

  static QLLMatchFactory __global_QLLMatchFactory__;

  QLLMatch *QLLMatch::_me = nullptr;

  void MIN_vtx_qll(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);

  QLLMatch::QLLMatch(const std::string name)
    : BaseFlashMatch(name), _mode(kChi2), _record(false), _normalize(false), _poisson("QLLMatch_poisson","TMath::Poisson(x,y)"), _minuit_ptr(nullptr)
  { _current_llhd = _current_chi2 = _current_pe = -1.0; }

  QLLMatch::QLLMatch()
  { throw OpT0FinderException("Use QLLMatch::GetME() to obtain singleton pointer!"); }

  void QLLMatch::_Configure_(const Config_t &pset) {
    _record = pset.get<bool>("RecordHistory");
    _check_touching_track = pset.get<bool>("CheckTouchingTrack");
    _normalize = pset.get<bool>("NormalizeHypothesis");
    _mode   = (QLLMode_t)(pset.get<unsigned short>("QLLMode"));
    _pe_observation_threshold = pset.get<double>("PEObservationThreshold", 1.e-6);
    _pe_hypothesis_threshold  = pset.get<double>("PEHypothesisThreshold", 1.e-6);
    _migrad_tolerance         = pset.get<double>("MIGRADTolerance", 0.1);
    _offset                   = pset.get<double>("Offset", 0.0);
		_time_shift               = pset.get<double>("BeamTimeShift", 0.0);
    _touching_track_window    = pset.get<double>("TouchingTrackWindow", 5.0);
    _minuit_x_buffer          = pset.get<double>("MinuitXBuffer", 10.0);
    // _custom_algo              = pset.get<std::string>("CustomAlgo", "");
    // if (!_custom_algo.empty()) _alg_custom_algo = CustomAlgoFactory::get().create(_custom_algo, _custom_algo);
    // if (_alg_custom_algo) {
    //
    // }
    _penalty_threshold_v = pset.get<std::vector<double> >("PEPenaltyThreshold");
    _penalty_value_v = pset.get<std::vector<double> >("PEPenaltyValue");

    _recox_penalty_threshold = pset.get<double>("XPenaltyThreshold");
    _recoz_penalty_threshold = pset.get<double>("ZPenaltyThreshold");

    _exp_frac_v.clear(); _exp_frac_v.push_back(0.23); _exp_frac_v.push_back(0.77);
    _exp_tau_v.clear();  _exp_tau_v.push_back(0.002); _exp_tau_v.push_back(1.5);

    _exp_frac_v = pset.get<std::vector<double> >("PhotonDecayFractions",_exp_frac_v);
    _exp_tau_v  = pset.get<std::vector<double> >("PhotonDecayTimes",_exp_tau_v);

    _xpos_v.resize(DetectorSpecs::GetME().NOpDets(),0.);
    _ypos_v.resize(DetectorSpecs::GetME().NOpDets(),0.);
    _zpos_v.resize(DetectorSpecs::GetME().NOpDets(),0.);
    for(size_t ch=0; ch<DetectorSpecs::GetME().NOpDets(); ++ch) {
      auto const& pmt_pos = DetectorSpecs::GetME().PMTPosition(ch);
      _xpos_v[ch] = pmt_pos[0];
      _ypos_v[ch] = pmt_pos[1];
      _zpos_v[ch] = pmt_pos[2];
    }
    auto const& bbox = DetectorSpecs::GetME().ActiveVolume();
    _vol_xmax = bbox.Max()[0];
    _vol_xmin = bbox.Min()[0];
  }


  std::vector<double> QLLMatch::CalculateX0(const Flash_t &pmt) {
    // scenarios include: vol_min, using flash time + tpc0 assumption, using flash time + tpc1 assumption
    FLASH_INFO() << "Setting initial X positions. Volume min " << _vol_xmin << " Flash time " << pmt.time << std::endl;
    std::vector<double> x0_v;

    // default: include the vol minimum
    x0_v.push_back(_vol_xmin);

    // if the trajectory is touching both TPC planes, really nothing to minimize (i.e. initial position suffices)
    bool touching_both_planes = (_vol_xmax - _vol_xmin) <= (_raw_xmax_pt.x - _raw_xmin_pt.x);
    if(touching_both_planes) return x0_v;


    // Assume this is the right flash, derive the true x position under this assumption
    double reco_x_tpc0 = _raw_xmin_pt.x - (pmt.time-_time_shift) * DetectorSpecs::GetME().DriftVelocity(); // assuming that the original track is in tpc0
    double reco_x_tpc1 = _raw_xmin_pt.x + (pmt.time-_time_shift) * DetectorSpecs::GetME().DriftVelocity(); // assuming that the original track is in tpc1

    double tolerance = _touching_track_window/2. * DetectorSpecs::GetME().DriftVelocity();
    // Inspect, in either assumption (original track is in tpc0 or tpc1), the track is contained in the whole active volume or not
    bool contained_tpc0 = reco_x_tpc0 >= _vol_xmin - tolerance && (reco_x_tpc0 + _raw_xmax_pt.x - _raw_xmin_pt.x) <= _vol_xmax + tolerance; 
    bool contained_tpc1 = reco_x_tpc1 >= _vol_xmin - tolerance && (reco_x_tpc1 + _raw_xmax_pt.x - _raw_xmin_pt.x) <= _vol_xmax + tolerance;
    FLASH_INFO() << " tpc0 " << reco_x_tpc0 << " " << (reco_x_tpc0 + _raw_xmax_pt.x - _raw_xmin_pt.x) << " contained? " << contained_tpc0 << std::endl;
    FLASH_INFO() << " tpc1 " << reco_x_tpc1 << " " << (reco_x_tpc1 + _raw_xmax_pt.x - _raw_xmin_pt.x) << " contained? " << contained_tpc1 << std::endl;
    if (contained_tpc0 && reco_x_tpc0 > _vol_xmin) {
      x0_v.push_back(std::max(reco_x_tpc0, _vol_xmin + _offset));
    }
    if (contained_tpc1) {
      x0_v.push_back(std::min(reco_x_tpc1, _vol_xmax - _offset));
    }

    return x0_v;
  }

  void QLLMatch::DumpHistory() const {
    FLASH_INFO() << "Dumping minimizer history..." << std::endl;
    for(size_t i=0; i<_minimizer_record_x_v.size(); ++i) {
      FLASH_INFO() << "Step " << i 
      << " X " << _minimizer_record_x_v[i] 
      << " PE " << _minimizer_record_pe_v[i]
      << " LLHD " << _minimizer_record_llhd_v[i] << std::endl;
    }
  }

  void QLLMatch::Match(const QCluster_t &pt_v, const Flash_t &flash, FlashMatch_t& match) {
    FLASH_INFO() << "Starting a match... "
    << "MC info: tpc true time " << pt_v.time_true 
    << " Flash true time " << flash.time_true << std::endl;    
    //
    // Prepare TPC
    //
    _raw_trk.resize(pt_v.size());
    double min_x =  1e20;
    double max_x = -1e20;
    for (size_t i = 0; i < pt_v.size(); ++i) {
      auto const &pt = pt_v[i];
      _raw_trk[i] = pt;
      if (pt.x < min_x) { min_x = pt.x; _raw_xmin_pt = pt; }
      if (pt.x > max_x) { max_x = pt.x; _raw_xmax_pt = pt; }
    }
    for (auto &pt : _raw_trk) pt.x -= min_x;

    // Calculate initial x positions
    auto x0_v = this->CalculateX0(flash);

    std::vector<FlashMatch_t> res_v;
    for(auto const& x0 : x0_v) {
      FLASH_INFO() << "Flash match with x0 = " << x0 << std::endl;
      FlashMatch_t res = match;
      PESpectrumMatch(flash,x0,res);
      
      if(_record)
        this->DumpHistory();
      
      res_v.emplace_back(std::move(res));
    }

    size_t best_res_idx = 0;
    double best_res_score = -1;
    for(size_t i=0; i<x0_v.size(); ++i) {
      FLASH_INFO() << "Using x0 = " << x0_v[i] 
      << " maximized 1/param Score=" << res_v[i].score 
      << " @ X=" << res_v[i].tpc_point.x << " [cm]" << std::endl;
      if(res_v[i].score > best_res_score) {
        best_res_idx = i;
        best_res_score = res_v[i].score;
      }
    }

    match = res_v[best_res_idx];
  }


  void QLLMatch::PESpectrumMatch(const Flash_t &flash, const double x0, FlashMatch_t& match) {

    match.num_steps = _num_steps;
    match.minimizer_min_x = _minimizer_min_x;
    match.minimizer_max_x = _minimizer_max_x;
    match.tpc_point.x = match.tpc_point.y = match.tpc_point.z = 0;
    match.score = 0;

    _hypothesis.time_width = flash.time_width;

    this->CallMinuit(flash, x0);

    // Estimate position
    if (std::isnan(_qll)) return;

    double weight = 0;

    for (size_t pmt_index = 0; pmt_index < DetectorSpecs::GetME().NOpDets(); ++pmt_index) {

      match.tpc_point.y += _ypos_v.at(pmt_index) * _hypothesis.pe_v[pmt_index];
      match.tpc_point.z += _zpos_v.at(pmt_index) * _hypothesis.pe_v[pmt_index];

      weight += _hypothesis.pe_v[pmt_index];
    }

    match.tpc_point.y /= weight;
    match.tpc_point.z /= weight;

    match.tpc_point.x = _reco_x_offset;
    match.tpc_point_err.x = _reco_x_offset_err;

    match.hypothesis  = _hypothesis.pe_v;

    //
    // Compute score
    //
    if(_mode == kSimpleLLHD)
      match.score = _qll * -1.;
    else
      match.score = 1. / _qll;

    // Compute X-weighting
    /*
    double x0 = _raw_xmin_pt.x - flash.time * DetectorSpecs::GetME().DriftVelocity();
    if( fabs(_reco_x_offset - x0) > _recox_penalty_threshold )
      match.score *= 1. / (1. + fabs(_reco_x_offset - x0) - _recox_penalty_threshold);
    // Compute Z-weighting
    double z0 = 0;
    weight = 0;
    for (size_t pmt_index = 0; pmt_index < DetectorSpecs::GetME().NOpDets(); ++pmt_index) {
      z0 += _zpos_v.at(pmt_index) * flash.pe_v[pmt_index];
      weight += flash.pe_v[pmt_index];
    }
    z0 /= weight;
    if( fabs(match.tpc_point.z - z0) > _recoz_penalty_threshold )
      match.score *= 1. / (1. + fabs(match.tpc_point.z - z0) - _recoz_penalty_threshold);
    */

    return;
  }

  const Flash_t &QLLMatch::ChargeHypothesis(const double xoffset) {
    //auto start = high_resolution_clock::now();
    if (_hypothesis.pe_v.empty()) _hypothesis.pe_v.resize(DetectorSpecs::GetME().NOpDets(), 0.);
    if (_hypothesis.pe_v.size() != DetectorSpecs::GetME().NOpDets()) {
      throw OpT0FinderException("Hypothesis vector length != PMT count");
    }

    for (auto &v : _hypothesis.pe_v) v = 0;

    // Apply xoffset
    _var_trk.resize(_raw_trk.size());
    for (size_t pt_index = 0; pt_index < _raw_trk.size(); ++pt_index) {
      //std::cout << "x point : " << _raw_trk[pt_index].x << "\t offset : " << xoffset << std::endl;
      _var_trk[pt_index].x = _raw_trk[pt_index].x + xoffset;
      _var_trk[pt_index].y = _raw_trk[pt_index].y;
      _var_trk[pt_index].z = _raw_trk[pt_index].z;
      _var_trk[pt_index].q = _raw_trk[pt_index].q;
    }

    //auto end = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration ChargeHypothesis 1 = " << duration.count() << "us" << std::endl;

    //start = high_resolution_clock::now();
    FillEstimate(_var_trk, _hypothesis);
    //end = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration ChargeHypothesis 2 = " << duration.count() << "us" << std::endl;

    //start = high_resolution_clock::now();
    if (_normalize) {
      double qsum = std::accumulate(std::begin(_hypothesis.pe_v),
				    std::end(_hypothesis.pe_v),
				    0.0);
      for (auto &v : _hypothesis.pe_v) v /= qsum;
    }
    //end = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration ChargeHypothesis 3 = " << duration.count() << "us" << std::endl;

    return _hypothesis;
  }

  const Flash_t &QLLMatch::Measurement() const { return _measurement; }

  double QLLMatch::QLL(const Flash_t &hypothesis,
		       const Flash_t &measurement) 
  {

    double nvalid_pmt = 0;
    _current_chi2 = _current_llhd = _current_pe = 0.;

    if (measurement.pe_v.size() != hypothesis.pe_v.size())
      throw OpT0FinderException("Cannot compute QLL for unmatched length!");

    double integral_factor = 0;
    for(size_t i=0; i<_exp_frac_v.size(); ++i) {
      integral_factor += _exp_frac_v[i] * (1 - exp(-1 * measurement.time_width / _exp_tau_v[i]));
    }
    assert(integral_factor > 0);

    double O, H, Error;
    const double epsilon = 1.e-300;

    int count_observation = 0;
    int count_hypothesis = 0;
    for (size_t pmt_index = 0; pmt_index < hypothesis.pe_v.size(); ++pmt_index) {

      O = measurement.pe_v[pmt_index] / integral_factor; // observation
      H = hypothesis.pe_v[pmt_index];  // hypothesis
      _current_pe += H;
      if( H < 0 ) throw OpT0FinderException("Cannot have hypothesis value < 0!");

      if(O < _pe_observation_threshold) ++count_observation;
      if(H < _pe_hypothesis_threshold ) ++count_hypothesis;
    }
    FLASH_DEBUG() << count_observation << " (Reco) v.s. " << count_hypothesis << " (Hypothesis) PMTs below PE threshold " << std::endl;

    for (size_t pmt_index = 0; pmt_index < hypothesis.pe_v.size(); ++pmt_index) {

      O = measurement.pe_v[pmt_index] / integral_factor; // observation
      H = hypothesis.pe_v[pmt_index];  // hypothesis

      if( H < 0 ) throw OpT0FinderException("Cannot have hypothesis value < 0!");

      // O(bservation) must be above threshold if set
      if(O < _pe_observation_threshold) {
        if (!_penalty_value_v.empty()) {
          O = _penalty_value_v[pmt_index];
        }
        else {
          O = _pe_observation_threshold;
        }
      }

      // H(ypothesis) must be above threshold if set
      if (H < _pe_hypothesis_threshold) {
        if(!_penalty_threshold_v.empty()) {
          H = _penalty_threshold_v[pmt_index];
        }
        else {
          H = _pe_hypothesis_threshold;
        }
      }

      //_current_pe += H;

      if(_mode == kLLHD) {
      	assert(H>0);
      	double arg = TMath::Poisson(O,H) + epsilon;
      	if(!std::isnan(arg) && !std::isinf(arg)) {
      	  _current_llhd -= std::log10(arg);
      	  nvalid_pmt += 1;
      	  if(_converged) FLASH_DEBUG() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " LHD "<<arg << " -LLHD " << -1 * std::log10(arg) << std::endl;
      	}
      }
      else if(_mode == kWeightedLLHD) {
      	assert(H>0);
      	double arg = TMath::Poisson(O,H) + epsilon;
      	if(!std::isnan(arg) && !std::isinf(arg)) {
      	  _current_llhd -= std::log10(arg * sqrt(std::max(H,epsilon)));
      	  nvalid_pmt += 1;
      	  if(_converged) FLASH_DEBUG() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " LHD "<<arg << " -LLHD " << -1 * std::log10(arg * sqrt(std::max(H,epsilon))) << std::endl;
      	}
      }
      else if(_mode == kZIP) {
        double pzero = H > _pe_hypothesis_threshold ? 0. : 1.;
        double arg = O > _pe_observation_threshold ? (TMath::Poisson(O,H) * (1-pzero) + epsilon) : (pzero + (1 - pzero)*TMath::Exp(-H) + epsilon);
        if(!std::isnan(arg) && !std::isinf(arg)) {
            _current_llhd -= std::log10(arg);
            nvalid_pmt += 1;
            if(_converged) FLASH_DEBUG() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " LHD "<<arg << " -LLHD " << -1 * std::log10(arg) << std::endl;
        }
      }
      else if(_mode == kPEWeightedLLHD) {
        assert(H>0);
      	double arg = TMath::Poisson(O,H) + epsilon;
      	if(!std::isnan(arg) && !std::isinf(arg)) {
      	  _current_llhd -= std::log10(arg * sqrt(std::max(H,epsilon))) * std::max(H,epsilon)/hypothesis.TotalPE() ;
      	  nvalid_pmt += 1;
      	  if(_converged) FLASH_DEBUG() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " LHD "<<arg << " -LLHD " << -1 * std::log10(arg * sqrt(std::max(H,epsilon))) * std::max(H,epsilon)/hypothesis.TotalPE() << std::endl;
        }
      }
      else if(_mode == kIntegralLLHD) {
      	double hmin = H-0.5;
      	double hmax = H+0.5;
      	double omin = O-sqrt(O);
      	double omax = O+sqrt(O);
      	if(hmin<0.) hmin = 0.;
      	if(omin<0.) omin = 0.;
      	if(hmax<hmin+1) hmax = hmin+1;
      	if(omax<omin+1) omax = omin+1;
      	double arg = _poisson.Integral(omin,omax,hmin,hmax) + epsilon;
      	if(!std::isnan(arg) && !std::isinf(arg)) {
      	  _current_llhd -= std::log10(arg);
      	  nvalid_pmt += 1;
      	  if(_converged) FLASH_DEBUG() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " LHD "<<arg << " -LLHD " << -1 * std::log10(arg) << std::endl;
      	}
      }else if (_mode == kSimpleLLHD) {

      	double arg = (H - O * std::log(H));
      	_current_llhd += arg;
      	if(_converged) FLASH_DEBUG() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " ... -LLHD " << arg << std::endl;
      	//nvalid_pmt += 1;

      } else if (_mode == kChi2) {

      	Error = O;
      	if( Error < 1.0 ) Error = 1.0;
      	_current_chi2 += std::pow((O - H), 2) / (Error);
      	nvalid_pmt += 1;

      } else {
      	FLASH_ERROR() << "Unexpected mode" << std::endl;
      	throw OpT0FinderException();
      }
    }
    //FLASH_DEBUG() <<"Mode " << (int)(_mode) << " Chi2 " << _current_chi2 << " LLHD " << _current_llhd << " nvalid " << nvalid_pmt << std::endl;
    _current_chi2 /= nvalid_pmt;
    _current_llhd /= (nvalid_pmt +1);

    if(_converged)
      FLASH_INFO() << "Combined LLHD: " << _current_llhd << " (divided by nvalid_pmt+1 = " << nvalid_pmt+1<<")"<<std::endl;

    return (_mode == kChi2 ? _current_chi2 : _current_llhd);
  }

  void MIN_vtx_qll(Int_t & /*Npar*/, // Number of parameters
		   Double_t * /*Grad*/, // Partial derivatives (return values)
		   Double_t &Fval, // Function value (return value)
		   Double_t *Xval, // Parameter values
		   Int_t) /*Flag*/{ // flag word
    //std::cout << "minuit offset : " << Fval << std::endl;
    //std::cout << "minuit Xval?? : " << *Xval << std::endl;

    //auto start = high_resolution_clock::now();
    auto const &hypothesis = QLLMatch::GetME()->ChargeHypothesis(*Xval);
    //std::cout << hypothesis.TotalPE() << " @x=" << Xval[0] << std::endl;
    //auto end = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration ChargeHypothesis = " << duration.count() << "us" << std::endl;

    //start = high_resolution_clock::now();
    auto const &measurement = QLLMatch::GetME()->Measurement();
    //end = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration Measurement = " << duration.count() << "us" << std::endl;

    //start = high_resolution_clock::now();
    Fval = QLLMatch::GetME()->QLL(hypothesis, measurement);
    //std::cout << "Fval " << Fval << " x = " << Xval[0] << std::endl;
    //end = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration QLL = " << duration.count() << "us" << std::endl;

    QLLMatch::GetME()->Record(Xval[0]);
    QLLMatch::GetME()->OneStep(Xval[0]);

    return;
  }

  double QLLMatch::CallMinuit(const Flash_t &pmt, const double x0) {

    if (_measurement.pe_v.empty()) {
      _measurement.pe_v.resize(DetectorSpecs::GetME().NOpDets(), 0.);
    }
    if (_measurement.pe_v.size() != pmt.pe_v.size()) {
      FLASH_CRITICAL() << _measurement.pe_v.size() << " v.s. " << pmt.pe_v.size() << std::endl;
      throw OpT0FinderException("PMT dimension has changed!");
    }

    if (!_penalty_threshold_v.empty() && _penalty_threshold_v.size() != pmt.pe_v.size()) {
      throw OpT0FinderException("Penalty threshold array has a different size than PMT array size!");
    }

    if (!_penalty_value_v.empty() && _penalty_value_v.size() != pmt.pe_v.size()) {
      throw OpT0FinderException("Penalty value array has a different size than PMT array size!");
    }

    _converged = false;

    //
    // Prepare PMT
    //
    double max_pe = 1.;

    // Debug: Print out expected PE spectrum
    //for(size_t i=0; i<pmt.pe_v.size(); ++i) {
    //std::cout << "PE meas: " << i << " " << pmt.pe_v[i] << std::endl;
    //}

    if (_normalize) {
      max_pe = 0;
      for (auto const &v : pmt.pe_v) if (v > max_pe) max_pe = v;
    }

    for (size_t i = 0; i < pmt.pe_v.size(); ++i)  _measurement.pe_v[i] = pmt.pe_v[i] / max_pe;

    _minimizer_record_chi2_v.clear();
    _minimizer_record_llhd_v.clear();
    _minimizer_record_x_v.clear();
    _minimizer_record_pe_v.clear();
    _num_steps = 0;
    _minimizer_min_x = std::numeric_limits<double>::max();
    _minimizer_max_x = -std::numeric_limits<double>::max();

    //double reco_x = _vol_xmin + _offset;
    double reco_x = x0 + _offset;

    double reco_x_err = ((_vol_xmax - _vol_xmin) - (_raw_xmax_pt.x - _raw_xmin_pt.x)) / 2.;
    double xmin = std::max(_vol_xmin, _vol_xmin - _minuit_x_buffer);
    double xmax = std::min(_vol_xmax, (_vol_xmax - _vol_xmin) - (_raw_xmax_pt.x - _raw_xmin_pt.x) + _vol_xmin + _minuit_x_buffer);
    FLASH_INFO() << _raw_xmax_pt.x << " " << _raw_xmin_pt.x << " " << (xmin < xmax) << " " << (xmin == xmax) << std::endl;
    FLASH_INFO() << "Running Minuit x: " << xmin << " => " << xmax
     << " (x buffer set at " << _minuit_x_buffer << ")"
		 << " ... initial state x=" <<reco_x <<" x_err=" << reco_x_err << std::endl;


    if (!_minuit_ptr) _minuit_ptr = new TMinuit(4);
    double MinFval;
    int ierrflag, npari, nparx, istat;
    double arglist[4], Fmin, Fedm, Errdef;
    ierrflag = npari = nparx = istat = 0;

    assert(this == QLLMatch::GetME());

    _minuit_ptr->SetPrintLevel(-1);
    arglist[0] = 2.0;  // set strategy level
    _minuit_ptr->mnexcm("SET STR", arglist, 1, ierrflag);

    _minuit_ptr->SetFCN(MIN_vtx_qll);

    _minuit_ptr->DefineParameter(0, "X", reco_x, reco_x_err, xmin, xmax);

    _minuit_ptr->Command("SET NOW");

    // use Migrad minimizer

    arglist[0] = 5000;  // maxcalls
    arglist[1] = _migrad_tolerance; // tolerance*1e-3 = convergence condition
    _minuit_ptr->mnexcm("MIGRAD", arglist, 2, ierrflag);

    _converged = true;

    //arglist[0]   = 5.0e+2;
    //arglist[1]   = 1.0e-6;
    //_minuit_ptr->mnexcm ("simplex",arglist,2,ierrflag);
    FLASH_INFO() << " reco x before " << reco_x << std::endl;
    _minuit_ptr->GetParameter(0, reco_x, reco_x_err);
    FLASH_INFO() << " reco x after " << reco_x << std::endl;

    _minuit_ptr->mnstat(Fmin, Fedm, Errdef, npari, nparx, istat);

    // use this for debugging, maybe uncomment the actual minimzing function (MIGRAD / simplex calls)
    // scanning the parameter set
    //arglist[0] = 0;       // Parameter No (in our case x is the only and therefore the 0th parameter)
    //arglist[1] = 500;    // Number of points
    //arglist[2] = 0;       // Start point of scan
    //arglist[3] = 256;     // End point of scan
    //_minuit_ptr->mnexcm("scan", arglist,4, ierrflag);

    double *grad = 0;
    int nPar = 1;
    double fValue[1];
    fValue[0] = reco_x;

    // Transfer the minimization variables:
    MIN_vtx_qll(nPar, grad, Fmin, fValue, ierrflag);
    MinFval = Fmin;
    _reco_x_offset = reco_x;
    _reco_x_offset_err = reco_x_err;
    _qll = MinFval;

    // Clear:
    _minuit_ptr->mnexcm("clear", arglist, 0, ierrflag);

    if (_minuit_ptr) delete _minuit_ptr;
    _minuit_ptr = 0;

    return _qll;
  }

}
#endif
