#include "CRTHisto.h"

#include "TH1D.h"
#include "TH2D.h"

void ana::SBNOsc::CRTHistos::Initialize(const std::string &postfix, const std::vector<double> &tagger_volume) {
#define CRT_HISTO2D(name, binx, lo_x, hi_x, biny, lo_y, hi_y)  name = new TH2D((#name"_" + postfix).c_str(), #name, binx, lo_x, hi_x, biny, lo_y, hi_y); fAllHistos.push_back(name) 

  unsigned n_bins = 200;
  CRT_HISTO2D(crt_hits_xy, n_bins, tagger_volume[0], tagger_volume[3], n_bins, tagger_volume[1], tagger_volume[4]);
  CRT_HISTO2D(crt_hits_xz, n_bins, tagger_volume[0], tagger_volume[3], n_bins, tagger_volume[2], tagger_volume[5]);
  CRT_HISTO2D(crt_hits_yz, n_bins, tagger_volume[1], tagger_volume[4], n_bins, tagger_volume[2], tagger_volume[5]);
#undef CRT_HISTO2D
}

void ana::SBNOsc::CRTHistos::Get(TFile &f, const std::string &postfix) {
#define GET_CRT_HISTO2D(name, postfix) name = (TH2D *) f.Get((#name "_" + postfix).c_str()); fAllHistos.push_back(name);
  GET_CRT_HISTO2D(crt_hits_xy, postfix);
  GET_CRT_HISTO2D(crt_hits_xz, postfix);
  GET_CRT_HISTO2D(crt_hits_yz, postfix);
#undef GET_CRT_HISTO2D
}

void ana::SBNOsc::CRTHistos::Fill(const numu::CRTHit &hit) {
  crt_hits_xy->Fill(hit.location.X(), hit.location.Y()); 
  crt_hits_xz->Fill(hit.location.X(), hit.location.Z()); 
  crt_hits_yz->Fill(hit.location.Y(), hit.location.Z()); 
}


void ana::SBNOsc::CRTHistos::Fill(const sbnd::crt::CRTHit &hit) {
  crt_hits_xy->Fill(hit.x_pos, hit.y_pos);  
  crt_hits_xz->Fill(hit.x_pos, hit.z_pos);  
  crt_hits_yz->Fill(hit.y_pos, hit.z_pos);  
}

