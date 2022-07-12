#ifndef OPT0FINDER_SelectionGreedy_CXX
#define OPT0FINDER_SelectionGreedy_CXX

#include "SelectionGreedy.h"

namespace flashmatch {
  
  static SelectionGreedyFactory __global_SelectionGreedyFactory__;  

  SelectionGreedy::SelectionGreedy(const std::string name)
  : BaseMatchSelection(name)
  {}

  void SelectionGreedy::_Configure_(const Config_t &pset)
  {
    _score_max_threshold = pset.get<double>("TouchMatchMaxThreshold",2.0);
    _score_min_threshold = pset.get<double>("FlashScoreMinThreshold",0.);
    _score_max_ceiling = pset.get<double>("FlashScoreMaxCeiling",1.0);
    _allow_reuse_flash = pset.get<bool>("AllowReuseFlash",false);
  }

  std::vector<FlashMatch_t> SelectionGreedy::Select(const std::vector<std::vector<FlashMatch_t> >& match_data)
  {

    std::vector<FlashMatch_t> result;
    result.reserve(match_data.size());

    // 
    // Step 1: select a touch match from the lowest to the highest touch match score 
    // Step 2: select a flash match from the highest to the lowest flash match score
    //

    // Ordered score container
    std::multimap<double,FlashMatch_t> score_map;

    // Step 1.1 register matches with a touch-match score
    for(auto const& match_v : match_data) {
      for(auto const& match : match_v) {
        if(match.touch_match == flashmatch::kNoTouchMatch) continue;
        if(std::fabs(match.touch_score) > _score_max_threshold) continue;
        score_map.insert(std::make_pair(std::fabs(match.touch_score),match));
      }
    }

    // Step 1.2 select matches w/o re-using tpc/flash object (unless allow to reuse flash)
    std::set<ID_t> tpc_used, flash_used;
    for (auto& score_info : score_map) {

      auto&       match_info  = score_info.second;   // match information
      auto const& tpc_index   = match_info.tpc_id;   // matched tpc original id
      auto const& flash_index = match_info.flash_id; // matched flash original id


      // If this tpc object is already assigned (=better match found), ignore
      if (tpc_used.find(tpc_index) != tpc_used.end()) continue;

      // If this flash object is already assigned + re-use is not allowed, ignore
      if (!_allow_reuse_flash && flash_used.find(flash_index) != flash_used.end()) continue;

      // Reaching this point means a new match. Yay!
      FLASH_INFO () << "Concrete Match: " << " TPC=" << tpc_index << " Flash=" << flash_index
        << " Score=" << match_info.score
        << std::endl;

      // Register to a list of a "used" flash and tpc info
      tpc_used.insert(tpc_index);
      flash_used.insert(flash_index);

      // std::move matched info from the map to result vector
      result.push_back( match_info );
    }

    score_map.clear();

    // Step 2.1 register matches with a flash-match score (inverse in order to use multimap)
    for(auto const& match_v : match_data) {
      for(auto const& match : match_v) {
        if(match.score < _score_min_threshold) continue;
        double score = std::min(match.score, _score_max_ceiling);
        score_map.insert(std::make_pair(1./(score - _score_min_threshold) ,match));
      }
    }

    // Loop over score map created with matching algorithm
    for (auto& score_info : score_map) {

      auto&       match_info  = score_info.second;   // match information
      auto const& tpc_index   = match_info.tpc_id;   // matched tpc original id
      auto const& flash_index = match_info.flash_id; // matched flash original id


      // If this tpc object is already assigned (=better match found), ignore
      if (tpc_used.find(tpc_index) != tpc_used.end()) continue;

      // If this flash object is already assigned + re-use is not allowed, ignore
      if (!_allow_reuse_flash && flash_used.find(flash_index) != flash_used.end()) continue;

      // Reaching this point means a new match. Yay!
      FLASH_INFO () << "Concrete Match: " << " TPC=" << tpc_index << " Flash=" << flash_index
        << " Score=" << match_info.score
        << std::endl;

      // Register to a list of a "used" flash and tpc info
      tpc_used.insert(tpc_index);
      flash_used.insert(flash_index);

      // std::move matched info from the map to result vector
      result.push_back( match_info );
    }

    return result;
  }


}

#endif
