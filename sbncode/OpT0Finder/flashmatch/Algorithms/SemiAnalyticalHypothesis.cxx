#ifndef SEMIANALYTICALHYPOTHESIS_CXX
#define SEMIANALYTICALHYPOTHESIS_CXX

#include "SemiAnalyticalHypothesis.h"
#include <chrono>

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

using namespace std::chrono;
namespace flashmatch {

  static SemiAnalyticalHypothesisFactory __global_SemiAnalyticalHypothesisFactory__;

  SemiAnalyticalHypothesis::SemiAnalyticalHypothesis(const std::string name)
   : PhotonLibHypothesis(name)
  {
#if USING_LARSOFT == 1
    _opfast_scintillation = new larg4::OpFastScintillation();
#endif
  }

  //
  // LArSoft
  //
  #if USING_LARSOFT == 1


  void SemiAnalyticalHypothesis::BuildHypothesis(const QCluster_t& trk, Flash_t &flash) const
  {

    for ( size_t ipt = 0; ipt < trk.size(); ++ipt) {

      /// Get the 3D point in space from where photons should be propagated
      auto const& pt = trk[ipt];

      // Get the number of photons produced in such point
      double n_original_photons = pt.q;

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

  #else

  // Not implemented outside larsoft

  void SemiAnalyticalHypothesis::BuildHypothesis(const QCluster_t& trk, Flash_t &flash) const
  {}

  #endif

}
#endif
