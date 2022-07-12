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
#include "TouchMatchFactory.h"
#include "MatchSelectionFactory.h"
#include <chrono>

//using namespace std::chrono;
namespace flashmatch {

  FlashMatchManager::FlashMatchManager(const std::string name)
    : LoggerFeature(name)
    , _alg_flash_filter(nullptr)
    , _alg_tpc_filter(nullptr)
    , _alg_match_prohibit(nullptr)
    , _alg_flash_match(nullptr)
    , _alg_flash_hypothesis(nullptr)
    , _alg_touch_match(nullptr)
    , _alg_match_select(nullptr)
    , _configured(false)
    , _name(name)
  {}

  const std::string& FlashMatchManager::Name() const
  { return _name; }

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
    
    auto const& mgr_cfg = main_cfg.get<flashmatch::Config_t>(Name());

    this->set_verbosity((msg::Level_t)(mgr_cfg.get<unsigned int>("Verbosity")));
    _store_full = mgr_cfg.get<bool>("StoreFullResult");

    auto const flash_filter_name = mgr_cfg.get<std::string>("FlashFilterAlgo","");
    auto const tpc_filter_name   = mgr_cfg.get<std::string>("TPCFilterAlgo",  "");
    auto const prohibit_name     = mgr_cfg.get<std::string>("ProhibitAlgo",   "");
    auto const hypothesis_name   = mgr_cfg.get<std::string>("HypothesisAlgo"    ); // required
    auto const match_name        = mgr_cfg.get<std::string>("MatchAlgo"         ); // required
    auto const touch_match_name  = mgr_cfg.get<std::string>("TouchMatchAlgo", "");
    auto const match_select_name = mgr_cfg.get<std::string>("MatchSelectionAlgo"); // required
    std::vector<std::string> custom_algo_v;
    custom_algo_v = mgr_cfg.get<std::vector<std::string> >("CustomAlgo",custom_algo_v);

    if(!flash_filter_name.empty()) _alg_flash_filter     = FlashFilterFactory::get().create(flash_filter_name,flash_filter_name);
    if(!tpc_filter_name.empty()  ) _alg_tpc_filter       = TPCFilterFactory::get().create(tpc_filter_name,tpc_filter_name);
    if(!prohibit_name.empty()    ) _alg_match_prohibit   = FlashProhibitFactory::get().create(prohibit_name,prohibit_name);
    if(!hypothesis_name.empty()  ) _alg_flash_hypothesis = FlashHypothesisFactory::get().create(hypothesis_name,hypothesis_name);
    if(!match_name.empty()       ) _alg_flash_match      = FlashMatchFactory::get().create(match_name,match_name);
    if(!touch_match_name.empty() ) _alg_touch_match      = TouchMatchFactory::get().create(touch_match_name,touch_match_name);
    if(!match_select_name.empty()) _alg_match_select     = MatchSelectionFactory::get().create(match_select_name,match_select_name);
    for(auto const& name : custom_algo_v)
      if(!name.empty()) AddCustomAlgo(CustomAlgoFactory::get().create(name,name));

    // checks

    if (_alg_flash_filter) {
      FLASH_INFO() << "Configuring FlashFilter " << _alg_flash_filter->AlgorithmName() << std::endl;
      _alg_flash_filter->Configure(main_cfg.get<flashmatch::Config_t>(_alg_flash_filter->AlgorithmName()));
    }

    if (_alg_tpc_filter) {
      FLASH_INFO() << "Configuring TPCFilter " << _alg_tpc_filter->AlgorithmName() << std::endl;
      _alg_tpc_filter->Configure(main_cfg.get<flashmatch::Config_t>(_alg_tpc_filter->AlgorithmName()));
    }

    if (_alg_match_prohibit) {
      FLASH_INFO() << "Configuring MatchProhibit " << _alg_match_prohibit->AlgorithmName() << std::endl;
      _alg_match_prohibit->Configure(main_cfg.get<flashmatch::Config_t>(_alg_match_prohibit->AlgorithmName()));
    }

    if (_alg_flash_hypothesis) {
      FLASH_INFO() << "Configuring Hypothesis " << _alg_flash_hypothesis->AlgorithmName() << std::endl;
      _alg_flash_hypothesis->Configure(main_cfg.get<flashmatch::Config_t>(_alg_flash_hypothesis->AlgorithmName()));
    }

    for (auto& name_ptr : _custom_alg_m) {
      FLASH_INFO() << "Configuring CustomAlgo " << name_ptr.first << std::endl;
      name_ptr.second->Configure(main_cfg.get<flashmatch::Config_t>(name_ptr.first));
    }

    if (_alg_flash_match) {
      FLASH_INFO() << "Configuring FlashMatch " << _alg_flash_match->AlgorithmName() << std::endl;
      _alg_flash_match->SetFlashHypothesis(_alg_flash_hypothesis);
      _alg_flash_match->Configure(main_cfg.get<flashmatch::Config_t>(_alg_flash_match->AlgorithmName()));
    }

    if (_alg_touch_match) {
      FLASH_INFO() << "Configuring TouchMatch " << _alg_touch_match->AlgorithmName() << std::endl;
      _alg_touch_match->Configure(main_cfg.get<flashmatch::Config_t>(_alg_touch_match->AlgorithmName()));
    }

    if(_alg_match_select) {
      FLASH_INFO() << "Configuring MatchSelect " << _alg_match_select->AlgorithmName() << std::endl;
      _alg_match_select->Configure(main_cfg.get<flashmatch::Config_t>(_alg_match_select->AlgorithmName()));
    }

    _configured = true;

    this->PrintConfig();

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

    // Touch match
    case kTouchMatch:
      return _alg_touch_match;

    // Match select
    case kMatchSelect:
      return _alg_match_select;

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

  void FlashMatchManager::Add(flashmatch::QCluster_t obj)
  { _tpc_object_v.push_back(obj); }

  void FlashMatchManager::Emplace(flashmatch::QCluster_t&& obj)
  { _tpc_object_v.emplace_back(std::move(obj)); }

  void FlashMatchManager::Add(flashmatch::Flash_t obj)
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

    // Report what's put into matching
    if(this->logger().level() == flashmatch::msg::kINFO) {
      for(size_t idx=0; idx<tpc_index_v.size(); ++idx) {
        auto const& tpc_index = tpc_index_v[idx];
        auto const& qcluster = _tpc_object_v[tpc_index];
        FLASH_INFO() << "Input QCluster " << idx << " (ID=" << tpc_index << ") ... "
        << qcluster.size() << " pts ... point qsum " << qcluster.sum()
        << " ... X span " << qcluster.min_x() << " => " << qcluster.max_x() << std::endl << std::endl; 
      }

      for(size_t idx=0; idx<flash_index_v.size(); ++idx) {
        auto const& flash_index = flash_index_v[idx];
        auto const& flash = _flash_v[flash_index];
        FLASH_INFO() << "Input Flash " << idx << " (ID=" << flash_index << ") ... "
        << "Reco " << flash.time << " [us] " << flash.TotalPE() << " p.e. "
        << "... True " << flash.time_true << " [us] " << flash.TotalTruePE() << " p.e. " << std::endl << std::endl;
      }
    }

    // use multi-map for possible equally-scored matches
    std::vector<std::vector<FlashMatch_t> > match_result;
    match_result.resize(tpc_index_v.size());
    for(auto& match_v : match_result) match_v.resize(flash_index_v.size());

    // Double loop over a list of tpc object & flash
    // Call matching function to inspect the compatibility.

    for (size_t tpc_index = 0; tpc_index < tpc_index_v.size(); ++tpc_index) {

      auto const& tpc = _tpc_object_v[tpc_index_v[tpc_index]]; // Retrieve TPC object

      // Loop over flash list
      for (size_t flash_index=0; flash_index < flash_index_v.size(); ++flash_index) {

        FLASH_INFO() << "TPC/Flash index " << tpc_index << "/" << flash_index
          << " ID " << tpc_index_v[tpc_index] << "/" << flash_index_v[flash_index] << std::endl;
        auto const& flash = _flash_v[flash_index_v[flash_index]];    // Retrieve flash

        auto& match = match_result[tpc_index][flash_index];
        match.tpc_id = tpc_index_v[tpc_index];
        match.flash_id = flash_index_v[flash_index];

        if (tpc.size() == 0 )
          continue;

        // run the match-prohibit algo first
        if (_alg_match_prohibit) {
          if(! _alg_match_prohibit->MatchCompatible( tpc, flash) )
            continue;
        }

        // run touch match algorithm
        if (_alg_touch_match)
          _alg_touch_match->Match(tpc,flash,match);

        if(match.tpc_id != tpc_index_v[tpc_index]){
          FLASH_CRITICAL() << "TPC ID changed by FlashMatch algorithm. Not supposed to happen..." << std::endl;
          throw OpT0FinderException();
        }
        if(match.flash_id != flash_index_v[flash_index]){
          FLASH_CRITICAL() << "Flash ID changed by FlashMatch algorithm. Not supposed to happen..." << std::endl;
          throw OpT0FinderException();
        }
        auto start = std::chrono::high_resolution_clock::now();
        _alg_flash_match->Match( tpc, flash, match ); // Run matching
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        FLASH_INFO() << "Match duration = " << duration.count() << "ns" << std::endl;

        match.duration = duration.count();

        if(_store_full) {
          _res_tpc_flash_v[match.tpc_id][match.flash_id] = match;
          _res_flash_tpc_v[match.flash_id][match.tpc_id] = match;
      	}

        FLASH_DEBUG() << "Candidate Match: "
          << " TPC=" << tpc_index_v[tpc_index] << " @ " << tpc.time
          << " with Flash=" << flash_index << " @ " << flash.time
          << " ... Score=" << match.score
          << " ... PE=" << flash.TotalPE()
          << " /hyp. PE=" << std::accumulate(match.hypothesis.begin(), match.hypothesis.end(), 0)
          << std::endl;
      }
    }

    // We have a score-ordered list of match information at this point.

    result = _alg_match_select->Select(match_result);

    for(size_t idx=0; idx < result.size(); ++idx) {
      auto const& match = result[idx];
      auto const& tpc   = _tpc_object_v[match.tpc_id];
      auto const& flash = _flash_v[match.flash_id];

      FLASH_INFO() << "Match " << idx << " ... TPC " << match.tpc_id << " with PMT " << match.flash_id << std::endl
      << "  Match score       : " << match.score << " touch? " << match.touch_match << " (" << match.touch_score << ")" << std::endl 
      << "  X position        : true " << tpc.min_x_true << " hypo " << match.tpc_point.x << " touch " << match.touch_point.x << std::endl
      << "  Time              : true " << tpc.time_true << " flash " << flash.time << std::endl
      << "  PEs               : true " << flash.TotalTruePE() << " hypo " << std::accumulate(match.hypothesis.begin(),match.hypothesis.end(),0.) << " flash " << flash.TotalPE() << std::endl
      << std::endl;
    }

    // Return result
    return result;

  }

  void FlashMatchManager::PrintConfig() {

    std::cout << "---- FLASH MATCH MANAGER PRINTING CONFIG     ----" << std::endl
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
}

#endif
