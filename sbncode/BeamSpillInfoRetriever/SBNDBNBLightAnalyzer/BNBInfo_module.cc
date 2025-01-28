////////////////////////////////////////////////////////////////////////
// Class:       BNBInfo
// Plugin Type: analyzer (Unknown Unknown)
// File:        BNBInfo_module.cc
//
// Generated at Tue Jan 28 14:29:25 2025 by Max Dubnowski using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h" 


using namespace art;
using namespace std;
namespace test {
  class BNBInfo;
}


class test::BNBInfo : public art::EDAnalyzer {
public:
  explicit BNBInfo(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BNBInfo(BNBInfo const&) = delete;
  BNBInfo(BNBInfo&&) = delete;
  BNBInfo& operator=(BNBInfo const&) = delete;
  BNBInfo& operator=(BNBInfo&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void beginSubRun(art::SubRun const& sr) override;

private:
  std::vector<double> fsrGates;
  // Declare member data here.

};


test::BNBInfo::BNBInfo(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::BNBInfo::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

void test::BNBInfo::beginSubRun(art::SubRun const& sr)

{

  double totGates = 0;

  // InputTag bnb_tag { "sbndbnbextinfo" };
  InputTag bnb_tag { "sbndbnbinfo" };

  // auto const& bnb_handle = sr.getValidHandle<vector<sbn::EXTCountInfo>>(bnb_tag);
  auto const& bnb_handle = sr.getValidHandle<vector<sbn::BNBSpillInfo>>(bnb_tag);
  auto const& bnb_vec(*bnb_handle);

  for (size_t j = 0; j < bnb_vec.size(); j++){
    // totGates += bnb_vec.at(j).gates_since_last_trigger;
    double real_time = bnb_vec.at(j).spill_time_s + bnb_vec.at(j).spill_time_ns*1e-9;
    std::cout << "spill_time: " << std::setprecision(19) << real_time << std::endl;
    std::cout << "TOR860: " << std::setprecision(12) << bnb_vec.at(j).TOR860 << std::endl;
    std::cout << "TOR875: " << std::setprecision(12) << bnb_vec.at(j).TOR875 << std::endl;
    std::cout << "THCURR: " << std::setprecision(12) << bnb_vec.at(j).THCURR << std::endl;
  }

  fsrGates.push_back(totGates);
}



void test::BNBInfo::beginJob()
{
  // Implementation of optional member function here.
}

void test::BNBInfo::endJob()
{
  for(size_t k = 0; k < fsrGates.size(); k++){
    std::cout << "Gates: " << fsrGates.at(k) << std::endl;
  }
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::BNBInfo)
