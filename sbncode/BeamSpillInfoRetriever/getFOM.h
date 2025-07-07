#ifndef _GETFOM_H
#define _GETFOM_H

//#include "datatypes/ub_BeamHeader.h"
//#include "datatypes/ub_BeamData.h"
//#include "lardataobj/RawData/BeamInfo.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
//#include "bnbAutoTune.h"
#include <string>
#include "TMath.h"
#include "TH1.h"
#include "TGraph.h"
//#include "TH1D.h"
//#include "TF1.h"




namespace sbn
{
  //typedef std::vector<bnb::bnbAutoTune> autoTunes;
  //float getFOM2(std::string beam, const gov::fnal::uboone::datatypes::ub_BeamHeader& bh, const std::vector<gov::fnal::uboone::datatypes::ub_BeamData>& bd, const bnb::bnbAutoTune settings = bnb::bnbAutoTune(), bool useAutoTune=true);
  float getFOM(BNBSpillInfo& spill);

  /* autoTunes cacheAutoTuneHistory(); */
  /* bnb::bnbAutoTune getSettings(const autoTunes& history, const gov::fnal::uboone::datatypes::ub_BeamHeader& bh); */
  /* bnb::bnbAutoTune getSettings(const autoTunes& history, const raw::BeamInfo& bi); */
  /* bnb::bnbAutoTune getSettings(const autoTunes& history, const uint64_t utctstamp); */
  double calcFOM(double horpos,double horang,double verpos,double verang,double tor,double tgtsx,double tgtsy);
  void swimBNB(const double centroid1[6], const double sigma1[6][6], 
               const double xferc[6][6], const double xfers[6][6],
               double &cx, double& cy, double &sx, double &sy, double &rho);
  double func_intbivar(const double cx, const double cy, const double sx, const double sy, const double rho );
  void processBNBprofile(const double* mwdata, double &x, double& sx, double& chi2); 
  
}
#endif
