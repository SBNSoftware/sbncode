#ifndef PHOTONLIBHYPOTHESIS_CXX
#define PHOTONLIBHYPOTHESIS_CXX

#include "PhotonLibHypothesis.h"

#include <omp.h>
#define NUM_THREADS 4

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

using namespace std::chrono;
namespace flashmatch {

  static PhotonLibHypothesisFactory __global_PhotonLibHypothesisFactory__;

  PhotonLibHypothesis::PhotonLibHypothesis(const std::string name)
    : BaseFlashHypothesis(name)
    , _vis(art::ServiceHandle<phot::PhotonVisibilityService const>().get())
  {}

  void PhotonLibHypothesis::_Configure_(const Config_t &pset)
  {
    #if USING_LARSOFT == 0
    omp_set_num_threads(NUM_THREADS);
    #endif

    _global_qe = pset.get<double>("GlobalQE");
    _global_qe_refl = pset.get<double>("GlobalQERefl", 0);
    _use_semi_analytical = pset.get<bool>("UseSemiAnalytical", 0);
    _calc_recombination  = pset.get<bool>("CalcRecombination", 1);

    _qe_v.clear();
    _qe_v = pset.get<std::vector<double> >("CCVCorrection",_qe_v);
    if(_qe_v.empty()) _qe_v.resize(DetectorSpecs::GetME().NOpDets(),1.0);
    if(_qe_v.size() != DetectorSpecs::GetME().NOpDets()) {
      FLASH_CRITICAL() << "CCVCorrection factor array has size " << _qe_v.size()
                       << " != number of opdet (" << DetectorSpecs::GetME().NOpDets() << ")!" << std::endl;
      throw OpT0FinderException();
    }

    // By default, add all opdets to the channel mask
    // Note that this may be overridden by the manager
    // via the SetChannelMask() method.
    _channel_mask.clear();
    _channel_mask.reserve(DetectorSpecs::GetME().NOpDets());
    for (size_t i = 0; i < DetectorSpecs::GetME().NOpDets(); i++) {
      _channel_mask[i] = i;
    }
  }

  //
  // LArSoft
  //
  #if USING_LARSOFT == 1

  void PhotonLibHypothesis::FillEstimate(const QCluster_t& trk, Flash_t &flash) const
  {

    for ( auto& v : flash.pe_v ) v = 0;

    if (_use_semi_analytical) {
      FillEstimateSemiAnalytical(trk, flash);
    } else {
      FillEstimateLibrary(trk, flash);
    }
  }


  void PhotonLibHypothesis::FillEstimateSemiAnalytical(const QCluster_t& trk, Flash_t &flash) const
  {
    // double test_pitch  = trk.pitch();
    // std::cout << "pitch: " << test_pitch << std::endl;

    for ( size_t ipt = 0; ipt < trk.size(); ++ipt) {

      /// Get the 3D point in space from where photons should be propagated
      auto const& pt = trk[ipt];

      geo::Point_t const xyz = {pt.x, pt.y, pt.z};
      auto const track = pt.trk; // track==1 if point belongs to a track, != 1 if belongs to a shower/non-track/unassociated

      // recomb calculation
      double n_original_photons = 0;
      if (track!=1){
        // attempt to correct for attenuation
        double drift_time = abs(200.0 - abs(pt.x))/(0.16); // in us, drift velocity = 0.16 cm/us 
        double atten_correction = std::exp(drift_time/10e3); // electron lifetime = 10e3 us, or 10 ms 

        // Get the number of photons produced in such point
        double charge = atten_correction*pt.q;
      // end attempt 
        // set up recomb variables
        // double E_field = 0.5;
        double W_ion = 2.36016*1e-5; 
        double W_ph  = 19.5*1e-6;
        double ds = 0.5;
        // double ds = (track==1)? test_pitch : 0.5; 
        // if (track==1 && (std::isnan(test_pitch) || std::isinf(test_pitch)))
          // ds = 0.5; 
        // Modified Box constants 
        double ModBoxA = 0.93;
        double ModBoxB = 0.305344; 

        // Birks Model constants 
        // double A_b = 0.8;
        // double k_b = 0.0349993;
        
        double e_dep = (ds/ModBoxB)*(exp(ModBoxB*W_ion*charge/ds)-ModBoxA); // energy deposition 
        // if (charge > 26000) // use modbox
        // else // use birks 
        //   e_dep = (charge)/((A_b/W_ion) - (k_b/E_field)*(charge/ds));
        n_original_photons = (e_dep/W_ph) - charge;
      }
      if (track==1){
        n_original_photons = pt.q;
      }

      std::vector<double> direct_visibilities;
      _semi_model->detectedDirectVisibilities(direct_visibilities, xyz);

      std::vector<double> reflected_visibilities;
      _semi_model->detectedReflectedVisibilities(reflected_visibilities, xyz);

      //
      // Fill Estimate with Direct light
      //
      for (size_t op_det=0; op_det<direct_visibilities.size(); ++op_det) {

        const double visibility = direct_visibilities[op_det];

        double q = n_original_photons * visibility * _global_qe / _qe_v[op_det];

        // Coated PMTs don't see direct photons
        if (std::find(_uncoated_pmt_list.begin(), _uncoated_pmt_list.end(), op_det) != _uncoated_pmt_list.end()) {
          q = 0;
        }

        // std::cout << "OpDet: " << op_det << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q << "]" << std::endl;

        if (std::find(_channel_mask.begin(), _channel_mask.end(), op_det) != _channel_mask.end()) {
          flash.pe_v[op_det] += q;
        } else {
          flash.pe_v[op_det] = 0;
        }
      }

      //
      // Fill Estimate with Reflected light
      //
      for (size_t op_det=0; op_det<reflected_visibilities.size(); ++op_det) {

        const double visibility = reflected_visibilities[op_det];

        double q = n_original_photons * visibility * _global_qe_refl / _qe_v[op_det];

        if (std::find(_channel_mask.begin(), _channel_mask.end(), op_det) != _channel_mask.end()) {
          flash.pe_v[op_det] += q;
        } else {
          flash.pe_v[op_det] = 0;
        }

      }

    }
  }



  void PhotonLibHypothesis::FillEstimateLibrary(const QCluster_t& trk, Flash_t &flash) const
  {

    static double xyz[3] = {0.};

    size_t n_pmt = DetectorSpecs::GetME().NOpDets();

    for ( size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {

      if (std::find(_channel_mask.begin(), _channel_mask.end(), ipmt) == _channel_mask.end()) {
        continue;
      }

      bool is_uncoated = false;
      if (std::find(_uncoated_pmt_list.begin(), _uncoated_pmt_list.end(), ipmt) != _uncoated_pmt_list.end()) {
        is_uncoated = true;
      }

      for (size_t ipt = 0; ipt < trk.size(); ++ipt) {
        auto const& pt = trk[ipt];

        double q = pt.q;

        xyz[0] = pt.x;
        xyz[1] = pt.y;
        xyz[2] = pt.z;

        double q_direct = 0, q_refl = 0;

        // Direct light
        if (is_uncoated) {
          q_direct = 0;
        } else {
          q_direct = q * _vis->GetVisibility(xyz, ipmt) * _global_qe / _qe_v[ipmt];
        }
        // std::cout << "PMT : " << ipmt << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q << "]" << std::endl;

        // Reflected light
        q_refl = q * _vis->GetVisibility(xyz, ipmt, true) * _global_qe_refl / _qe_v[ipmt];

        flash.pe_v[ipmt] += q_direct;
        flash.pe_v[ipmt] += q_refl;
      }
    }
    return;
  }


  //
  // Standalone
  //
  #else

  void PhotonLibHypothesis::FillEstimate(const QCluster_t& trk, Flash_t &flash) const
  {

    static double xyz[3] = {0.};

    size_t n_pmt = DetectorSpecs::GetME().NOpDets();

    for ( auto& v : flash.pe_v ) v = 0;


    size_t n_pmt = DetectorSpecs::GetME().NOpDets();//n_pmt returns 0 now, needs to be fixed
    if(flash.pe_v.empty()) flash.pe_v.resize(n_pmt);
    if(flash.pe_err_v.empty()) flash.pe_err_v.resize(n_pmt);

    assert(flash.pe_v.size()     == n_pmt);
    assert(flash.pe_err_v.size() == n_pmt);

    for (auto& v : flash.pe_v     ) v = 0;
    for (auto& v : flash.pe_err_v ) v = 0;

    auto det = DetectorSpecs::GetME();

    auto const& lib_data = DetectorSpecs::GetME().GetPhotonLibraryData();

    //start = high_resolution_clock::now();
    #pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      size_t num_threads = omp_get_num_threads();
      size_t num_pts = trk.size() / num_threads;
      size_t start_pt = num_pts * thread_id;
      if(thread_id+1 == num_threads) num_pts += (trk.size() % num_threads);

      auto const& vox_def = DetectorSpecs::GetME().GetVoxelDef();
      // auto s = vox_def.GetVoxelSize();
      // auto s1 = vox_def.GetRegionLowerCorner();
      // auto s2 = vox_def.GetRegionUpperCorner();
      // std::cout << s[0] << " " << s[1] << " " << s[2] << std::endl;
      // std::cout << s1[0] << " " << s2[1] << " " << s1[2] << std::endl;
      // std::cout << s2[0] << " " << s2[1] << " " << s2[2] << std::endl;
      std::vector<double> local_pe_v(n_pmt,0);
      int vox_id;
      for( size_t ipt = start_pt; ipt < start_pt + num_pts; ++ipt) {
        auto const& pt = trk[ipt];
        vox_id = vox_def.GetVoxelID(pt.x,pt.y,pt.z);
        if (vox_id < 0) continue;
        auto const& vis_pmt = lib_data[vox_id];
        for ( size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {
          local_pe_v[ipmt] += pt.q * vis_pmt[ipmt];
        }
      }
      #pragma omp critical
      for(size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {
        flash.pe_v[ipmt] += local_pe_v[ipmt] * _global_qe / _qe_v[ipmt];
      }

    }
    return;


  }

  void PhotonLibHypothesis::FillEstimateSemiAnalytical(const QCluster_t& trk, Flash_t &flash) const
  {}

  void PhotonLibHypothesis::FillEstimateLibrary(const QCluster_t& trk, Flash_t &flash) const
  {}

  #endif

}
#endif
