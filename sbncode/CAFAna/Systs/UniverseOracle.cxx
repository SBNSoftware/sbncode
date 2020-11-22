#include "CAFAna/Systs/UniverseOracle.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  UniverseOracle& UniverseOracle::Instance()
  {
    static UniverseOracle uo;
    return uo;
  }

  // --------------------------------------------------------------------------
  UniverseOracle::UniverseOracle()
  {
    const std::string dir = "/sbnd/data/users/bzamoran/workshop-game-0320/weights/";

	const std::string fname = dir+"combined_params_2020a.txt";
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
	  if(name.find("genie")!=std::string::npos){
            std::string junk;
            ifin >> junk; // second version of the name...
            assert(ifin.good());
          }
          else{
            // Fix up the flux names
            if(name == "expskin" ||
               name == "horncurrent" ||
               name == "nucleoninexsec" ||
               name == "nucleonqexsec" ||
               name == "nucleontotxsec" ||
               name == "pioninexsec" ||
               name == "pionqexsec" ||
               name == "piontotxsec") name += "_FluxUnisim";
          }

          while(true){
            double shift;
            ifin >> shift;

            if(!ifin.good()) break; // either the next label or EOF
            fData[name].push_back(shift);
          } // end loop over values
        } // end loop over systs

    // Check we didn't screw up the file reading
    for(auto it: fData){
      const unsigned int N = it.second.size();
      if(N != 500){
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
  std::vector<SystShifts> UniverseOracle::ShiftsForSysts(const std::vector<const ISyst*>& systs, int nUniv) const
  {
    std::vector<SystShifts> ret(nUniv);

    for(const ISyst* s: systs){
      const std::vector<double>& xs = ShiftsForSyst(s->ShortName());
      for(int i = 0; i < nUniv; ++i) ret[i].SetShift(s, xs[i%xs.size()]);
    }

    return ret;
  }

  // --------------------------------------------------------------------------
  unsigned int UniverseOracle::ClosestIndex(const std::string& name,
                                            double shift,
                                            ESide side,
                                            double* trueShift) const
  {
    const std::vector<double>& v = ShiftsForSyst(name);
    int bestIdx = -1;
    double bestDist;
    for(unsigned int i = 0; i < v.size(); ++i){
      if(side == ESide::kBelow && v[i] > shift) continue;
      if(side == ESide::kAbove && v[i] < shift) continue;
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
