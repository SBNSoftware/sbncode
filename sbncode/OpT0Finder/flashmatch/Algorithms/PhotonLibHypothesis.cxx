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
    : BaseFlashHypothesis(name),
    _opfast_scintillation(new larg4::OpFastScintillation())
  {}

  void PhotonLibHypothesis::_Configure_(const Config_t &pset)
  {
    #if USING_LARSOFT == 0
    omp_set_num_threads(NUM_THREADS);
    #endif
    _global_qe = pset.get<double>("GlobalQE");
    _global_qe_refl = pset.get<double>("GlobalQERefl", 0);
    _use_semi_analytical = pset.get<bool>("UseSemiAnalytical", 0);
    _uncoated_pmt_list = pset.get<std::vector<int>>("UncoatedPMTList"); // FIXME This should be provided by the geometry service, eventually

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

    // static double xyz[3] = {0.};

    // size_t n_pmt = DetectorSpecs::GetME().NOpDets();

    for ( size_t ipt = 0; ipt < trk.size(); ++ipt) {

      /// Get the 3D point in space from where photons should be propagated
      auto const& pt = trk[ipt];

      // Get the number of photons produced in such point
      double n_original_photons = pt.q;

      // std::map<size_t, int> temp;
      // geo::Point_t const temp_xyz = {-100, 100, 200};
      // _opfast_scintillation->detectedReflecHits(temp, 50000, temp_xyz);
      // std::cout << "*** PhotonLibHypothesis Number of reflected photons to ch 118: " << temp[118] << std::endl;
      // std::map<size_t, int> temp_2;
      // _opfast_scintillation->detectedDirectHits(temp_2, 50000, temp_xyz);
      // std::cout << "*** PhotonLibHypothesis Number of direct photons to ch 118: " << temp_2[118] << std::endl;

      geo::Point_t const xyz = {pt.x, pt.y, pt.z};

      std::map<size_t, int> direct_photons;
      _opfast_scintillation->detectedDirectHits(direct_photons, n_original_photons, xyz);

      std::map<size_t, int> reflected_photons;
      _opfast_scintillation->detectedReflecHits(reflected_photons, n_original_photons, xyz);

      //
      // Fill Estimate with Direct light
      //
      for (auto it = direct_photons.begin(); it != direct_photons.end(); ++it) {

        const size_t op_det = it->first;
        const int n_photons = it->second;

        double q = n_photons * _global_qe / _qe_v[op_det];

        // Coated PMTs don't see direct photons
        if (std::find(_uncoated_pmt_list.begin(), _uncoated_pmt_list.end(), op_det) != _uncoated_pmt_list.end()) {
          q = 0;
        }

        // std::cout << "OpDet: " << op_det << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q << std::endl;

        if (std::find(_channel_mask.begin(), _channel_mask.end(), op_det) != _channel_mask.end()) {
          flash.pe_v[op_det] += q;
        } else {
          flash.pe_v[op_det] = 0;
        }
      }

      //
      // Fill Estimate with Reflected light
      //
      for (auto it = reflected_photons.begin(); it != reflected_photons.end(); ++it) {

        const size_t op_det = it->first;
        const int n_photons = it->second;

        double q = n_photons * _global_qe_refl / _qe_v[op_det];

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

    art::ServiceHandle<phot::PhotonVisibilityService const> vis;
    for ( size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {

      if (std::find(_channel_mask.begin(), _channel_mask.end(), ipmt) == _channel_mask.end()) {
        continue;
      }

      for ( size_t ipt = 0; ipt < trk.size(); ++ipt) {
        auto const& pt = trk[ipt];

        double q = pt.q;

        // q *= ::phot::PhotonVisibilityService::GetME().GetVisibility( pt.x, pt.y, pt.z, ipmt) * _global_qe / _qe_v[ipmt];
        xyz[0] = pt.x;
        xyz[1] = pt.y;
        xyz[2] = pt.z;

        // Direct light
        q *= vis->GetVisibility(xyz, ipmt) * _global_qe / _qe_v[ipmt];
        flash.pe_v[ipmt] += q;
        // std::cout << "direct -- ipmt: " << ipmt << " - xyx: " << " - visibility: " << vis->GetVisibility(xyz, ipmt) << " - q: " << q << std::endl;
        // std::cout << "PMT : " << ipmt << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q << std::endl;

        // Reflected light
        q *= vis->GetVisibility(xyz, ipmt, true) * _global_qe_refl / _qe_v[ipmt];
        flash.pe_v[ipmt] += q;
        // std::cout << "reflected -- ipmt: " << ipmt << " - xyx: " << " - visibility: " << vis->GetVisibility(xyz, ipmt, true) << " - q: " << q << std::endl;
      }
    }
    return;
  }


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
