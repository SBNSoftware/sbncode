#ifndef OPT0FINDER_FLASHMATCHMANAGER_CXX
#define OPT0FINDER_FLASHMATCHMANAGER_CXX

#include <sstream>
#include <map>
#include <set>
#include "FlashMatchManager.h"
#include "OpT0FinderException.h"
#include "FlashFilterFactory.h"
#include "TPCFilterFactory.h"
#include "FlashMatchFactory.h"
#include "FlashHypothesisFactory.h"
#include "FlashProhibitFactory.h"
#include "CustomAlgoFactory.h"
#include <chrono>

using namespace std::chrono;
namespace flashmatch {

  FlashMatchManager::FlashMatchManager(const std::string name)
    : LoggerFeature(name)
    , _alg_flash_filter(nullptr)
    , _alg_tpc_filter(nullptr)
    , _alg_match_prohibit(nullptr)
    , _alg_flash_match(nullptr)
    , _alg_flash_hypothesis(nullptr)
    , _configured(false)
    , _name(name)
  {
    _allow_reuse_flash = true;
  }

  const std::string& FlashMatchManager::Name() const
  { return _name; }

  /*
  void FlashMatchManager::SetAlgo(BaseAlgorithm* alg)
  {
    _configured = false;
    // Figure out the type of a provided algorithm
    switch (alg->AlgorithmType()) {

    // TPC filter
    case kTPCFilter:
      _alg_tpc_filter   = (BaseTPCFilter*)alg; break;

    // Flash filter
    case kFlashFilter:
      _alg_flash_filter = (BaseFlashFilter*)alg; break;

    // Match prohibit algo
    case kMatchProhibit:
      _alg_match_prohibit = (BaseProhibitAlgo*)alg; break;

    // Flash matching
    case kFlashMatch:
      _alg_flash_match  = (BaseFlashMatch*)alg; break;

    // Flash hypothesis
    case kFlashHypothesis:
      _alg_flash_hypothesis = (BaseFlashHypothesis*)alg; break;

    // Other algorithms
    case kCustomAlgo:
      if(_custom_alg_m.find(alg->AlgorithmName()) != _custom_alg_m.end()) {
	std::stringstream ss;
	ss << "Duplicate name: " << alg->AlgorithmName() << std::endl;
	throw OpT0FinderException(ss.str());
      }
      _custom_alg_m[alg->AlgorithmName()] = alg;
      break;
    // Fuck it
    default:
      std::stringstream ss;
      ss << "Unsupported algorithm type: " << alg->AlgorithmType();
      throw OpT0FinderException(ss.str());
    }
  }

  */

  void FlashMatchManager::AddCustomAlgo(BaseAlgorithm* alg)
  {
    if(_custom_alg_m.find(alg->AlgorithmName()) != _custom_alg_m.end()) {
      std::stringstream ss;
      ss << "Duplicate name: " << alg->AlgorithmName() << std::endl;
      throw OpT0FinderException(ss.str());
      }
    _custom_alg_m[alg->AlgorithmName()] = alg;
  }

  void FlashMatchManager::Configure(const Config_t& main_cfg)
  {
    /*
    ::fcllite::ConfigManager cfg_mgr("FlashMatchManager");

    cfg_mgr.AddCfgFile(_config_file);

    auto const& main_cfg = cfg_mgr.Config();
    */
    auto const& mgr_cfg = main_cfg.get<flashmatch::Config_t>(Name());

    _allow_reuse_flash = mgr_cfg.get<bool>("AllowReuseFlash");
    this->set_verbosity((msg::Level_t)(mgr_cfg.get<unsigned int>("Verbosity")));
    _store_full = mgr_cfg.get<bool>("StoreFullResult");

    auto const flash_filter_name = mgr_cfg.get<std::string>("FlashFilterAlgo","");
    auto const tpc_filter_name   = mgr_cfg.get<std::string>("TPCFilterAlgo","");
    auto const prohibit_name     = mgr_cfg.get<std::string>("ProhibitAlgo","");
    auto const hypothesis_name   = mgr_cfg.get<std::string>("HypothesisAlgo","");
    auto const match_name        = mgr_cfg.get<std::string>("MatchAlgo","");
    std::vector<std::string> custom_algo_v;
    custom_algo_v = mgr_cfg.get<std::vector<std::string> >("CustomAlgo",custom_algo_v);

    if(!flash_filter_name.empty()) _alg_flash_filter     = FlashFilterFactory::get().create(flash_filter_name,flash_filter_name);
    if(!tpc_filter_name.empty()  ) _alg_tpc_filter       = TPCFilterFactory::get().create(tpc_filter_name,tpc_filter_name);
    if(!prohibit_name.empty()    ) _alg_match_prohibit   = FlashProhibitFactory::get().create(prohibit_name,prohibit_name);
    if(!hypothesis_name.empty()  ) _alg_flash_hypothesis = FlashHypothesisFactory::get().create(hypothesis_name,hypothesis_name);
    if(!match_name.empty()       ) _alg_flash_match      = FlashMatchFactory::get().create(match_name,match_name);
    for(auto const& name : custom_algo_v)
      if(!name.empty()) AddCustomAlgo(CustomAlgoFactory::get().create(name,name));

    // checks
    if (_alg_match_prohibit)
      _alg_match_prohibit->Configure(main_cfg.get<flashmatch::Config_t>(_alg_match_prohibit->AlgorithmName()));

    if (_alg_flash_hypothesis)
      _alg_flash_hypothesis->Configure(main_cfg.get<flashmatch::Config_t>(_alg_flash_hypothesis->AlgorithmName()));

    if (_alg_flash_match) {
      _alg_flash_match->SetFlashHypothesis(_alg_flash_hypothesis);
      _alg_flash_match->Configure(main_cfg.get<flashmatch::Config_t>(_alg_flash_match->AlgorithmName()));
    }

    for (auto& name_ptr : _custom_alg_m)
      name_ptr.second->Configure(main_cfg.get<flashmatch::Config_t>(name_ptr.first));

    _configured = true;
  }

  BaseAlgorithm* FlashMatchManager::GetAlgo(flashmatch::Algorithm_t type)
  {
    if (!_configured)
      FLASH_WARNING() << "Algorithm may be not configured yet!" << std::endl;

    // Figure out the type of a provided algorithm
    switch (type) {

    // TPC filter
    case kTPCFilter:
      return _alg_tpc_filter;

    // Flash filter
    case kFlashFilter:
      return _alg_flash_filter;

    // Match prohibit algo
    case kMatchProhibit:
      return _alg_match_prohibit;

    // Flash matching
    case kFlashMatch:
      return _alg_flash_match;

    // Flash hypothesis
    case kFlashHypothesis:
      return _alg_flash_hypothesis;

    // Fuck it
    default:
      std::stringstream ss;
      ss << "Unsupported algorithm type: " << type;
      throw OpT0FinderException(ss.str());
    }
    return nullptr;
  }

  flashmatch::BaseAlgorithm* FlashMatchManager::GetCustomAlgo(std::string name)
  {
    if(_custom_alg_m.find(name) == _custom_alg_m.end()) {
      FLASH_ERROR() << "Algorithm name " << name << " not found!" << std::endl;
      throw OpT0FinderException();
    }
    return _custom_alg_m[name];
  }

  void FlashMatchManager::Add(flashmatch::QCluster_t& obj)
  { _tpc_object_v.push_back(obj); }

  void FlashMatchManager::Emplace(flashmatch::QCluster_t&& obj)
  { _tpc_object_v.emplace_back(std::move(obj)); }

  void FlashMatchManager::Add(flashmatch::Flash_t& obj)
  {
    if(!obj.Valid()) throw OpT0FinderException("Invalid Flash_t object cannot be registered!");
    _flash_v.push_back(obj);
  }

  void FlashMatchManager::Emplace(flashmatch::Flash_t&& obj)
  {
    if(!obj.Valid()) throw OpT0FinderException("Invalid Flash_t object cannot be registered!");
    _flash_v.emplace_back(std::move(obj));
  }

  // CORE FUNCTION
  std::vector<FlashMatch_t> FlashMatchManager::Match()
  {
    // Clear some history variables
    _res_tpc_flash_v.clear();
    _res_flash_tpc_v.clear();
    if(_store_full) {
      _res_tpc_flash_v.resize(_tpc_object_v.size(),std::vector<flashmatch::FlashMatch_t>(_flash_v.size()));
      _res_flash_tpc_v.resize(_flash_v.size(),std::vector<flashmatch::FlashMatch_t>(_tpc_object_v.size()));
    }

    // Create also a result container
    std::vector<FlashMatch_t> result;

    if (!_alg_flash_match)
      throw OpT0FinderException("Flash matching algorithm is reuqired! (not attached)");
    if (!_alg_flash_hypothesis)
      throw OpT0FinderException("Flash hypothesis algorithm is required! (not attached)");

    if (!_configured)
      throw OpT0FinderException("Have not configured yet!");

    if(_tpc_object_v.empty() || _flash_v.empty()) return result;

    //
    // Filter stage: for both TPC and Flash
    //

    // IDArray_t to store candidate list of tpc/flash to be used for matching
    IDArray_t tpc_index_v;
    IDArray_t flash_index_v;

    // Figure out which tpc object to use: if algorithm provided, ask it. Else use all.
    if (_alg_tpc_filter)
      tpc_index_v = _alg_tpc_filter->Filter(_tpc_object_v);
    else {
      tpc_index_v.reserve(_tpc_object_v.size());
      for (size_t i = 0; i < _tpc_object_v.size(); ++i) tpc_index_v.push_back(i);
    }

    FLASH_INFO() << "TPC Filter: " << _tpc_object_v.size() << " => " << tpc_index_v.size() << std::endl;

    // Figure out which flash to use: if algorithm provided, ask it. Else use all
    if (_alg_flash_filter)
      flash_index_v = _alg_flash_filter->Filter(_flash_v);
    else {
      flash_index_v.reserve(_flash_v.size());
      for (size_t i = 0; i < _flash_v.size(); ++i) flash_index_v.push_back(i);
    }
    FLASH_INFO() << "Flash Filter: " << _flash_v.size() << " => " << flash_index_v.size() << std::endl;

    //
    // Flash matching stage
    //

    // use multi-map for possible equally-scored matches
    std::multimap<double, FlashMatch_t> score_map;

    // Double loop over a list of tpc object & flash
    // Call matching function to inspect the compatibility.
    for (size_t tpc_index = 0; tpc_index < tpc_index_v.size(); ++tpc_index) {
      // Loop over flash list
      for (auto const& flash_index : flash_index_v) {
        FLASH_INFO() << "TPC index " << tpc_index << ", Flash index " << flash_index << std::endl;
        auto const& tpc   = _tpc_object_v[tpc_index_v[tpc_index]]; // Retrieve TPC object
        auto const& flash = _flash_v[flash_index];    // Retrieve flash

        if (tpc.size() == 0 )
          continue;

        // run the match-prohibit algo first
        if (_alg_match_prohibit) {
          bool compat = _alg_match_prohibit->MatchCompatible( tpc, flash);
          if (compat == false) {
            FLASH_INFO() << "Match not compatible. " << std::endl;
            continue;
          }
        }
        auto start = high_resolution_clock::now();
        auto res = _alg_flash_match->Match( tpc, flash ); // Run matching
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<nanoseconds>(end - start);
        FLASH_INFO() << "Match duration = " << duration.count() << "ns" << std::endl;

        // ignore this match if the score is <= 0
        if (res.score <= 0) continue;

        // Else we store this match. Assign TPC & flash index info
        res.tpc_id = tpc_index_v[tpc_index];//_index;
        res.flash_id = flash_index;//_index;
        res.duration = duration.count();


        if(_store_full) {
          _res_tpc_flash_v[res.tpc_id][res.flash_id] = res;
          _res_flash_tpc_v[res.flash_id][res.tpc_id] = res;
        }
        // For ordering purpose, take an inverse of the score for sorting
        score_map.emplace( 1. / res.score, res);

        FLASH_DEBUG() << "Candidate Match: "
		      << " TPC=" << tpc_index << " (" << tpc.min_x() << " min x)" << " @ " << tpc.time
		      << " with Flash=" << flash_index << " @ " << flash.time
		      << " ... Score=" << res.score
		      << " ... PE=" << flash.TotalPE()
		      << std::endl;
      }
    }

    // We have a score-ordered list of match information at this point.
    // Prepare return match information by respecting a score of each possible match.
    // Note _allow_reuse_flash becomes relevant here as well.

    // Create a std::set of tpc/flash IDs to keep track of already-matched tpc/flash input.
    std::set<ID_t> tpc_used, flash_used;
    result.reserve(tpc_index_v.size());
    // Loop over score map created with matching algorithm
    for (auto& score_info : score_map) {

      auto&       match_info  = score_info.second;   // match information
      auto const& tpc_index   = match_info.tpc_id;   // matched tpc original id
      auto const& flash_index = match_info.flash_id; // matched flash original id

//      std::cout<<"tpc_index and flash_index : "<<tpc_index<<", "<<flash_index<<std::endl ;

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
      result.emplace_back( match_info );

    }
    // Return result
    return result;

  }

  void FlashMatchManager::PrintConfig() {

    std::cout << "---- FLASH MATCH MANAGER PRINTING CONFIG     ----" << std::endl
	      << "_allow_reuse_flash = " << _allow_reuse_flash << std::endl
	      << "_name = " << _name << std::endl
	      << "_alg_flash_filter?" << std::endl;
    if (_alg_flash_filter)
      std::cout << "\t" << _alg_flash_filter->AlgorithmName() << std::endl;
    std::cout << "_alg_tpc_filter?" << std::endl;
    if (_alg_tpc_filter)
      std::cout << "\t" << _alg_tpc_filter->AlgorithmName() << std::endl;
    std::cout << "_alg_match_prohibit?" << std::endl;
    if (_alg_match_prohibit)
      std::cout << "\t" << _alg_match_prohibit->AlgorithmName() << std::endl;
    std::cout << "_alg_flash_hypothesis?" << std::endl;
    if (_alg_flash_hypothesis)
      std::cout << "\t" << _alg_flash_hypothesis->AlgorithmName() << std::endl;
    std::cout << "_alg_flash_match?" << std::endl;
    if (_alg_flash_match)
      std::cout << "\t" << _alg_flash_match->AlgorithmName() << std::endl;
    std::cout << "_custom_alg_m?" << std::endl;
    for (auto& name_ptr : _custom_alg_m)
      std::cout << "\t" << name_ptr.first << std::endl;
    std::cout << "---- END FLASH MATCH MANAGER PRINTING CONFIG ----" << std::endl;
  }

  void FlashMatchManager::SetChannelMask(std::vector<int> ch_mask) {

    if (!_alg_flash_hypothesis) {
      throw OpT0FinderException("Flash hypothesis algorithm is required to set channel mask!");
    }

    _alg_flash_hypothesis->SetChannelMask(ch_mask);
    
    if (_alg_flash_match) {
      _alg_flash_match->SetChannelMask(ch_mask);
    }
  }

  void FlashMatchManager::SetTPCCryo(int tpc, int cryo) {
    _tpc = tpc;
    _cryo = cryo;

    if (_alg_flash_match) {
      _alg_flash_match->SetTPCCryo(_tpc, _cryo);
    }
  }

  void FlashMatchManager::SetUncoatedPMTs(std::vector<int> ch_uncoated) {
    if (_alg_flash_hypothesis) {
      _alg_flash_hypothesis->SetUncoatedPMTs(ch_uncoated);
    }
  }

  // void FlashMatchManager::SetVisibleChannels(std::vector<int> ch_visible) {
  //   if (_alg_flash_hypothesis) {
  //     _alg_flash_hypothesis->SetVisibleChannels(ch_visible);
  //   }
  // }

  #if USING_LARSOFT == 1
  void FlashMatchManager::SetSemiAnalyticalModel(std::unique_ptr<phot::SemiAnalyticalModel> model) {
    if (_alg_flash_hypothesis) {
      _alg_flash_hypothesis->SetSemiAnalyticalModel(std::move(model));
    }
  }
  #endif

}

#endif
