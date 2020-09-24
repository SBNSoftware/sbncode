#include "CAFAna/Systs/SBNWeightSysts.h"

#include "CAFAna/Systs/UniverseOracle.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  UniverseWeight::UniverseWeight(const std::vector<std::string>& systs, int univIdx)
    : fNames(systs), fUnivIdx(univIdx)
  {
  }

  // --------------------------------------------------------------------------
  UniverseWeight::UniverseWeight(const std::vector<const ISyst*>& systs, int univIdx)
    : fUnivIdx(univIdx)
  {
    for(const ISyst* s: systs) fNames.push_back(s->ShortName());
  }

  // --------------------------------------------------------------------------
  double UniverseWeight::operator()(const caf::SRProxy* sr) const
  {
    if(sr->truth.empty()) return 1;
    const auto& ws = sr->truth[0].weights;

    // This hack improves throughput vastly
    if(fUnivIdx == 0){
      for(unsigned int i = 0; i < fNames.size(); ++i){
        for(const auto& b: ws[GetIdx(ws, i)].second) (void)((float)b);
      }
    }

    double w = 1;

    for(unsigned int i = 0; i < fNames.size(); ++i){
      // TODO: might want to "wrap around" differently in different systs to
      // avoid unwanted correlations between systs with the same number of
      // universes.
      const unsigned int unividx = fUnivIdx % ws[GetIdx(ws, i)].second.size();

      w *= ws[GetIdx(ws, i)].second[unividx];
    }

    return w;
  }

  // --------------------------------------------------------------------------
  int UniverseWeight::GetIdx(const caf::VectorProxy<caf::PairProxy>& ws,
                             int isyst) const
  {
    if(fSystIdxs.empty()){
      for(const std::string& name: fNames){
        bool any = false;
        for(unsigned int i = 0; i < ws.size(); ++i){
          if(ws[i].first == name){
            fSystIdxs.push_back(i);
            any = true;
            break;
          }
        }
        if(!any){
          std::cout << "UniverseWeight: Failed to find syst "
                    << name << " in file" << std::endl;
          abort();
        }
      } // end for name
    }

    return fSystIdxs[isyst];
  }


  // --------------------------------------------------------------------------
  SBNWeightSyst::SBNWeightSyst(const std::string& name)
    : ISyst(name, name), fIdx(-1)
  {
    assert(UniverseOracle::Instance().SystExists(name));
  }

  // --------------------------------------------------------------------------
  void SBNWeightSyst::Shift(double x, caf::SRProxy* sr, double& weight) const
  {
    if(sr->truth.empty()) return;

    const auto& ws = sr->truth[0].weights;

    const int i = GetIdx(ws);
    const Univs u = GetUnivs(x);

    const double y0 = ws[i].second[u.i0];
    const double y1 = ws[i].second[u.i1];

    weight *= u.w0*y0 + u.w1*y1;
  }

  // --------------------------------------------------------------------------
  int SBNWeightSyst::GetIdx(const caf::VectorProxy<caf::PairProxy>& ws) const
  {
    if(fIdx == -1){
      for(unsigned int i = 0; i < ws.size(); ++i){
        if(ws[i].first == ShortName()) fIdx = i;
      }

      if(fIdx == -1){
        std::cout << "SBNWeightSyst: '" << ShortName() << "' not found in record" << std::endl;
        std::cout << "Available:";
        for(auto& it: ws) std::cout << " " << std::string(it.first);
        std::cout << std::endl;
        abort();
      }
    }

    return fIdx;
  }

  // --------------------------------------------------------------------------
  SBNWeightSyst::Univs SBNWeightSyst::GetUnivs(double x) const
  {
    auto it = fUnivs.find(x);
    if(it != fUnivs.end()) return it->second;

    Univs u;
    const UniverseOracle& uo = UniverseOracle::Instance();
    // Neighbours
    double x0, x1;
    u.i0 = uo.ClosestIndex(ShortName(), x, ESide::kBelow, &x0);
    u.i1 = uo.ClosestIndex(ShortName(), x, ESide::kAbove, &x1);
    // Interpolation weights
    u.w0 = (x1-x)/(x1-x0);
    u.w1 = (x-x0)/(x1-x0);

    //      std::cout << ShortName() << " " << x << " sigma, found indices " << u.i0 << " and " << u.i1 << " at " << x0 << " and " << x1 << ", will use weights " << u.w0 << " and " << u.w1 << std::endl;

    // If one of the neighbours wasn't found, we fall back to just using the
    // neighbour we did find. It would probably be better to find two
    // neighbours on the same side and extrapolate.
    if(u.i0 == -1){u.i0 = u.i1; u.w0 = u.w1 = 0.5;}
    if(u.i1 == -1){u.i1 = u.i0; u.w0 = u.w1 = 0.5;}

    fUnivs[x] = u;
    return u;
  }

  // --------------------------------------------------------------------------
  const std::vector<const ISyst*>& GetSBNGenieWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(ret.empty()){
      const UniverseOracle& uo = UniverseOracle::Instance();
      ret.reserve(uo.Systs().size());
      for(const std::string& name: uo.Systs()){
        if(name.find("Alt") != std::string::npos) continue; // these seem to be duplicates
        if(name.find("genie") == std::string::npos) continue;
        ret.push_back(new SBNWeightSyst(name));
      }
    }
    return ret;
  }

  // --------------------------------------------------------------------------
  const std::vector<const ISyst*>& GetSBNFluxWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(ret.empty()){
      const UniverseOracle& uo = UniverseOracle::Instance();
      ret.reserve(uo.Systs().size());
      for(const std::string& name: uo.Systs()){
        if(name.find("genie") != std::string::npos) continue;
        ret.push_back(new SBNWeightSyst(name));
      }
    }
    return ret;
  }

  // --------------------------------------------------------------------------
  const std::vector<const ISyst*>& GetSBNWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(ret.empty()){
      const std::vector<const ISyst*>& g = GetSBNGenieWeightSysts();
      const std::vector<const ISyst*>& f = GetSBNFluxWeightSysts();
      ret.reserve(g.size()+f.size());
      ret.insert(ret.end(), g.begin(), g.end());
      ret.insert(ret.end(), f.begin(), f.end());
    }
    return ret;
  }
}
