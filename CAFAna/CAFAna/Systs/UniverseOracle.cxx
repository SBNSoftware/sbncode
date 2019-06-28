#include "CAFAna/Systs/UniverseOracle.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  UniverseOracle::UniverseOracle()
  {
    const std::string dir = "/sbnd/data/users/bckhouse/sample_2.1_fitters/";

    for(const std::string prefix: {"genie", "fluxunisim"}){
      // Must be in this order to match the indices in the CAFs
      for(const std::string ab: {"a", "b"}){
        const std::string fname = dir+prefix+"_params_"+ab+".txt";
        std::ifstream ifin(fname);
        if(ifin.fail()){
          std::cout << "UniverseOracle: Couldn't open file '" << fname << "'" << std::endl;
          abort();
        }

        while(true){
          ifin.clear(); // clear error flag
          std::string name;
          ifin >> name;
          if(!ifin.good()) break; // probably EOF
          if(prefix == "genie"){
            std::string junk;
            ifin >> junk; // second version of the name...
            assert(ifin.good());
          }

          while(true){
            double shift;
            ifin >> shift;

            if(!ifin.good()) break; // either the next label or EOF
            fData[name].push_back(shift);
          } // end loop over values
        } // end loop over systs
      } // end for ab
    } // end for prefix

    // Check we didn't screw up the file reading
    for(auto it: fData){
      const unsigned int N = it.second.size();
      if(N != 1000 && N != 50){
        std::cout << "UniverseOracle: Unexpected number of shifts (" << N << ") for '" << it.first << "'" << std::endl;
      }
    }
  }

  // --------------------------------------------------------------------------
  bool UniverseOracle::SystExists(const std::string& name) const
  {
    return fData.find(name) != fData.end();
  }

  // --------------------------------------------------------------------------
  std::vector<std::string> UniverseOracle::Systs() const
  {
    std::vector<std::string> ret;
    ret.reserve(fData.size());
    for(auto it: fData) ret.push_back(it.first);
    return ret;
  }

  // --------------------------------------------------------------------------
  const std::vector<double>& UniverseOracle::ShiftsForSyst(const std::string& name) const
  {
    assert(SystExists(name));
    return fData.find(name)->second;
  }

  // --------------------------------------------------------------------------
  unsigned int UniverseOracle::ClosestIndex(const std::string& name,
                                            double shift,
                                            double* trueShift) const
  {
    const std::vector<double>& v = ShiftsForSyst(name);
    int bestIdx = -1;
    double bestDist;
    for(unsigned int i = 0; i < v.size(); ++i){
      const double dv = fabs(v[i]-shift);
      if(bestIdx == -1 || dv < bestDist){
        bestIdx = i;
        bestDist = dv;
        if(trueShift) *trueShift = v[i];
      }
    }
    return unsigned(bestIdx);
  }
}
