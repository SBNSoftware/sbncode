#ifndef _WEIGHT_T_H_
#define _WEIGHT_T_H_

#include "WeightCalc.h"

namespace sbncode {
namespace evwgh {
  struct Weight_t {
    
  Weight_t() :
    fMinWeight(std::numeric_limits<double>::max()),
    fMaxWeight(std::numeric_limits<double>::min()),
    fAvgWeight(0),
    fNcalls(0)
    {}

    std::vector<std::vector<double> > GetWeight(art::Event& e) {
      std::vector<std::vector<double> >wgh=fWeightCalc->GetWeight(e);
      for (unsigned int inu=0;inu<wgh.size();inu++) {
	double avgwgh=std::accumulate(wgh[inu].begin(),wgh[inu].end(),0.0)/wgh[inu].size();
	fAvgWeight=(fAvgWeight*fNcalls+avgwgh)/float(fNcalls+1);
	fMinWeight=std::min(fMinWeight,
			    *std::min_element(wgh[inu].begin(),wgh[inu].end()));
	fMaxWeight=std::max(fMaxWeight,
			    *std::max_element(wgh[inu].begin(),wgh[inu].end()));
	fNcalls++;
      }
      
      return wgh;
    }
    std::string fName;
    WeightCalc* fWeightCalc;
    std::string fWeightCalcType;
    double      fMinWeight;
    double      fMaxWeight;
    double      fAvgWeight;
    long        fNcalls;
    int         fNmultisims;
  };

}
}

#endif // _WEIGHT_T_H_
