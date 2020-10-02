////////////////////////////////////////////////////////////////////////
// Class:       CRTHitConverter
// Plugin Type: producer (art v3_02_06)
// File:        CRTHitConverter_module.cc
//
// Generated at Wed Feb 19 17:38:21 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

//#include "Products/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTHit_Legacy.hh"

#include <memory>

namespace sbn {
  class CRTHitConverter;
}


class sbn::CRTHitConverter : public art::EDProducer {
public:
  explicit CRTHitConverter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTHitConverter(CRTHitConverter const&) = delete;
  CRTHitConverter(CRTHitConverter&&) = delete;
  CRTHitConverter& operator=(CRTHitConverter const&) = delete;
  CRTHitConverter& operator=(CRTHitConverter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fCRTHitLabel;
  std::string fExperiment;

};

sbn::crt::CRTHit SBNDCRTHit(const sbnd::crt::CRTHit &inp) {
  sbn::crt::CRTHit ret;

  ret.peshit = inp.peshit;
  ret.ts0_s = inp.ts0_s;
  ret.ts0_s_corr = inp.ts0_s_corr;
  ret.ts0_ns = inp.ts0_ns;
  ret.ts0_ns_corr = inp.ts0_ns_corr; 
  ret.ts1_ns = inp.ts1_ns;
  ret.plane = inp.plane;
  ret.x_pos = inp.x_pos;
  ret.x_err = inp.x_err;
  ret.y_pos = inp.y_pos;
  ret.y_err = inp.y_err;
  ret.z_pos = inp.z_pos;
  ret.z_err = inp.z_err;
  ret.tagger = inp.tagger;

  return ret;
} 

sbn::crt::CRTHit ICARUSCRTHit(const icarus::crt::CRTHit& inp) {
  sbn::crt::CRTHit ret;
  ret.feb_id = inp.feb_id;
  ret.pesmap = inp.pesmap;
  // convert ADC -> PE's
  // TODO: fix -- hardcoded for now as temporary hack
  unsigned n_strip = 2;
  double baseline = 63.6; // ADC
  double gain = 70; // ADC / PE
  ret.peshit = (inp.peshit - n_strip*baseline) / (gain * n_strip);
  ret.ts0_s = inp.ts0_s;
  ret.ts0_s_corr = inp.ts0_s_corr;
  ret.ts0_ns = inp.ts0_ns;
  ret.ts0_ns_corr = inp.ts0_ns_corr;

  ret.ts1_ns = inp.ts1_ns;

  ret.plane = inp.plane;
  ret.x_pos = inp.x_pos;
  ret.x_err = inp.x_err;
  ret.y_pos = inp.y_pos;
  ret.y_err = inp.y_err;
  ret.z_pos = inp.z_pos;
  ret.z_err = inp.z_err;
  ret.tagger = inp.tagger;
  return ret;
}


sbn::CRTHitConverter::CRTHitConverter(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fCRTHitLabel(p.get<std::string>("CRTHitLabel", "crthit")),
    fExperiment(p.get<std::string>("Experiment"))
{
  if (fExperiment != "SBND" && fExperiment != "ICARUS") {
    std::cerr << "ERROR: Bad experiment configuration (" << fExperiment << "). Experiment must be one of 'SBND' or 'ICARUS'";
  }

  produces< std::vector<sbn::crt::CRTHit> >();
}

void sbn::CRTHitConverter::produce(art::Event& e)
{

  std::unique_ptr<std::vector<sbn::crt::CRTHit>> crthits(new std::vector<sbn::crt::CRTHit>);

  if (fExperiment == "SBND") {
    art::Handle<std::vector<sbnd::crt::CRTHit>> crthit_handle;
    e.getByLabel(fCRTHitLabel, crthit_handle);
    if (crthit_handle.isValid()) {
      for (const sbnd::crt::CRTHit &sbndhit: *crthit_handle) {
        crthits->push_back(SBNDCRTHit(sbndhit));
      }
    }
  }
  else if (fExperiment == "ICARUS") {
    art::Handle<std::vector<icarus::crt::CRTHit>> crthit_handle;
    e.getByLabel(fCRTHitLabel, crthit_handle);
    if (crthit_handle.isValid()) {
      for (const icarus::crt::CRTHit &icarushit: *crthit_handle) {
        crthits->push_back(ICARUSCRTHit(icarushit));
      }
    }
  }

  e.put(std::move(crthits));

}

DEFINE_ART_MODULE(sbn::CRTHitConverter)
