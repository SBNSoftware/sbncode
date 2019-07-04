#include "CAFAna/Core/GenieWeightList.h"

#include <algorithm>
#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  std::vector<std::string> GetGenieWeightNames()
  {
    // I wonder if we can share this somehow with the code that generates these?
    return {"MaCCQE","VecFFCCQEshape",
	"MaNCEL","EtaNCEL","MaCCRES","MvCCRES",
	"MaNCRES","MvNCRES","RDecBR1gamma","RDecBR1eta",
	"Theta_Delta2Npi","AhtBY","BhtBY","CV1uBY",
	"CV2uBY","FormZone","MFP_pi","FrCEx_pi",
	"FrElas_pi","FrInel_pi","FrAbs_pi","FrPiProd_pi",
	"MFP_N","FrCEx_N","FrElas_N","FrInel_N",
	"FrAbs_N","FrPiProd_N","CCQEPauliSupViaKF","Mnv2p2hGaussEnhancement",
	"MKSPP_ReWeight","E2p2h_A_nu","E2p2h_B_nu","E2p2h_A_nubar",
	"E2p2h_B_nubar","NR_nu_n_CC_2Pi","NR_nu_n_CC_3Pi","NR_nu_p_CC_2Pi",
	"NR_nu_p_CC_3Pi","NR_nu_np_CC_1Pi","NR_nu_n_NC_1Pi","NR_nu_n_NC_2Pi",
	"NR_nu_n_NC_3Pi","NR_nu_p_NC_1Pi","NR_nu_p_NC_2Pi","NR_nu_p_NC_3Pi",
	"NR_nubar_n_CC_1Pi","NR_nubar_n_CC_2Pi","NR_nubar_n_CC_3Pi","NR_nubar_p_CC_1Pi",
	"NR_nubar_p_CC_2Pi","NR_nubar_p_CC_3Pi","NR_nubar_n_NC_1Pi","NR_nubar_n_NC_2Pi",
	"NR_nubar_n_NC_3Pi","NR_nubar_p_NC_1Pi","NR_nubar_p_NC_2Pi","NR_nubar_p_NC_3Pi",
	"BeRPA_A","BeRPA_B","BeRPA_D","BeRPA_E",
	"C12ToAr40_2p2hScaling_nu","C12ToAr40_2p2hScaling_nubar",
	"nuenuebar_xsec_ratio","nuenumu_xsec_ratio","SPPLowQ2Suppression"
	};
  }

  std::pair<double,double> GetGenieDialLimits(const std::string& name){

    static std::map<std::string,std::pair<double,double>> genieMap = {
      // Regular dials
      {"MaCCQE",{-3,3}},
      {"MaNCEL",{-3,3}},
      {"EtaNCEL",{-3,3}},
      {"MaCCRES",{-3,3}},
      {"MvCCRES",{-3,3}},
      {"MaNCRES",{-3,3}},
      {"MvNCRES",{-3,3}},
      {"RDecBR1gamma",{-3,3}},
      {"RDecBR1eta",{-3,3}},
      {"AhtBY",{-3,3}},
      {"BhtBY",{-3,3}},
      {"CV1uBY",{-3,3}},
      {"CV2uBY",{-3,3}},
      {"BeRPA_A",{-3,3}},
      {"BeRPA_B",{-3,3}},
      {"BeRPA_D",{-3,3}},
      {"BeRPA_E",{-3,3}},
      
      // NRpi dials
      {"NR_nu_n_CC_2Pi",{-2,3}},
      {"NR_nu_n_CC_3Pi",{-2,3}},
      {"NR_nu_p_CC_2Pi",{-2,3}},
      {"NR_nu_p_CC_3Pi",{-2,3}},
      {"NR_nu_np_CC_1Pi",{-2,3}},
      {"NR_nu_n_NC_1Pi",{-2,3}},
      {"NR_nu_n_NC_2Pi",{-2,3}},
      {"NR_nu_n_NC_3Pi",{-2,3}},
      {"NR_nu_p_NC_1Pi",{-2,3}},
      {"NR_nu_p_NC_2Pi",{-2,3}},
      {"NR_nu_p_NC_3Pi",{-2,3}},
      {"NR_nubar_n_CC_1Pi",{-2,3}},
      {"NR_nubar_n_CC_2Pi",{-2,3}},
      {"NR_nubar_n_CC_3Pi",{-2,3}},
      {"NR_nubar_p_CC_1Pi",{-2,3}},
      {"NR_nubar_p_CC_2Pi",{-2,3}},
      {"NR_nubar_p_CC_3Pi",{-2,3}},
      {"NR_nubar_n_NC_1Pi",{-2,3}},
      {"NR_nubar_n_NC_2Pi",{-2,3}},
      {"NR_nubar_n_NC_3Pi",{-2,3}},
      {"NR_nubar_p_NC_1Pi",{-2,3}},
      {"NR_nubar_p_NC_2Pi",{-2,3}},
      {"NR_nubar_p_NC_3Pi",{-2,3}},
      
      // FSI
      {"FormZone",{-2,2}},
      {"MFP_pi",{-2,2}},
      {"FrCEx_pi",{-2,2}},
      {"FrElas_pi",{-2,2}},
      {"FrInel_pi",{-2,2}},
      {"FrAbs_pi",{-2,2}},
      {"FrPiProd_pi",{-2,2}},
      {"MFP_N",{-2,2}},
      {"FrCEx_N",{-2,2}},
      {"FrElas_N",{-2,2}},
      {"FrInel_N",{-2,2}},
      {"FrAbs_N",{-2,2}},
      {"FrPiProd_N",{-2,2}},
      
      // On/Off dials
      {"VecFFCCQEshape",{-1,1}},
      {"Theta_Delta2Npi",{-1,1}},
      {"CCQEPauliSupViaKF",{-1,1}},
      {"Mnv2p2hGaussEnhancement",{-1,1}},
      {"MKSPP_ReWeight",{-1,1}},
      {"E2p2h_A_nu",{-1,1}},
      {"E2p2h_B_nu",{-1,1}},
      {"E2p2h_A_nubar",{-1,1}},
      {"E2p2h_B_nubar",{-1,1}},
      {"C12ToAr40_2p2hScaling_nu",{-1,1}},
      {"C12ToAr40_2p2hScaling_nubar",{-1,1}},
      {"nuenuebar_xsec_ratio",{-1,1}},
      {"nuenumu_xsec_ratio",{-1,1}},
      {"SPPLowQ2Suppression",{-1,1}}
    };
    
    if (genieMap.find(name) == genieMap.end()){
      std::cout << "Warning, no known systematic called " << name << " returning nominal" << std::endl;
      return {-3,3};
    }
    return genieMap[name];
  }
  
  
  //----------------------------------------------------------------------
  int GetGenieIndex(const std::string& name, bool quiet)
  {
    const std::vector<std::string> names = GetGenieWeightNames();

    auto it = std::find(names.begin(), names.end(), name);

    if(it == names.end()){
      if(!quiet){
        std::cerr << "Warning: couldn't find " << name
                  << " in list of genie systs" << std::endl;
      }
      return -1;
    }

    return it-names.begin();
  }

  std::string GetGenieWeightName(int index){
    const std::vector<std::string> names = GetGenieWeightNames();
    return names[index];
  }

  double GetGenieMin(int index){
    static const std::vector<std::string> names = GetGenieWeightNames();
    return GetGenieDialLimits(names[index]).first;
  }

  double GetGenieMax(int index){
    static const std::vector<std::string> names = GetGenieWeightNames();
    return GetGenieDialLimits(names[index]).second;
  }  
  
}
