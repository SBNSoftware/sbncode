#pragma once

#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Core/GenieWeightList.h"
#include "StandardRecord/StandardRecord.h"
#include <cmath>
#include <cassert>

namespace ana
{
  // class GenieSyst: public ISyst
  // {
  // public:
  //   virtual ~GenieSyst(){};

  //   void Shift(double sigma,
	 //       Restorer& restore,
	 //       caf::StandardRecord* sr,
	 //       double& weight) const override{

  //     assert(std::abs(sigma) <= 3 && "GENIE XSECs only valid up to +/-3 sigma!");

  //     // How far apart are the points
  //     const double spacing = 1;

  //     // Get the top and bottom values in the array
  //     int low_index  = std::floor(sigma/spacing) + 3;
  //     int high_index = std::ceil(sigma/spacing) + 3;

  //     double diff = (sigma-double(low_index))/spacing;

  //     double low_weight  = sr->dune.genie_wgt[fID][low_index];
  //     double high_weight = sr->dune.genie_wgt[fID][high_index];

  //     weight *= low_weight + (high_weight - low_weight)*diff;

  //   }

  // protected:
  // GenieSyst(int genie_id, bool applyPenalty=true) :
  //   ISyst(GetGenieWeightName(genie_id),
	 //  GetGenieWeightName(genie_id),
	 //  applyPenalty,
	 //  GetGenieMin(genie_id),
	 //  GetGenieMax(genie_id)),
	 //  fID(genie_id) {}
  //   friend std::vector<const ISyst*> GetGenieSysts(std::vector<std::string> names, bool applyPenalty);
	  
  //   int fID;
  // };

  // std::vector<const ISyst*> GetGenieSysts(std::vector<std::string> names = {}, bool applyPenalty=true){
  //   static std::vector<const ISyst*> ret;

  //   if (names.empty()) names = GetGenieWeightNames();

  //   if(ret.empty()){
  //     for (auto & it : names){
  //       ret.push_back(new GenieSyst(GetGenieIndex(it), applyPenalty));
  //     }
  //   }
  //   return ret;
  // }

}
