#include "OscCalculatorSterile.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

namespace osc
{
  OscCalculatorSterile::OscCalculatorSterile()
    : fPMNS_Sterile(0), fNFlavors(3), fDirty(true), fPrevE(0), fPrevAnti(0), fPrevFlavBefore(0)
  {
    fPMNS_Sterile = new PMNS_Sterile(fNFlavors);
  }

  //---------------------------------------------------------------------------
  std::vector<double> OscCalculatorSterile::GetState() const
  {
    std::vector<double> state;
    state.push_back((double)fNFlavors);
    state.push_back(fL);
    state.push_back(fRho);
    for(int i = 2; i <= fNFlavors; ++i) state.push_back(GetDm(i));
    for(int j = 2; j <= fNFlavors; ++j) {
      for(int i = 1; i < j; ++i) {
        state.push_back(GetAngle(i, j));
        if(i+1 != j) state.push_back(GetDelta(i, j));
      }
    }
    return state;
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetState(std::vector<double> state)
  {
    int iState(0);
    fDirty = true ;
    SetNFlavors((int)state[iState++]);
    SetL(state[iState++]);
    SetRho(state[iState++]);
    for (int i = 2; i <= fNFlavors; ++i) SetDm(i, state[iState++]);
    for(int j = 2; j <= fNFlavors; ++j) {
      for(int i = 1; i < j; ++i) {
        SetAngle(i, j, state[iState++]);
        if(i+1 != j) SetDelta(i, j, state[iState++]);
      }
    }
  }

  //---------------------------------------------------------------------------
  OscCalculatorSterile::OscCalculatorSterile(const OscCalculatorSterile& calc)
    : OscCalculatorSterile()
  {
    std::vector<double> state = calc.GetState();
    SetState(state);
  }

  //---------------------------------------------------------------------------
  TMD5* OscCalculatorSterile::GetParamsHash() const
  {
    TMD5* ret = new TMD5;
    std::string txt = "PMNSSterile";
    ret->Update((unsigned char*)txt.c_str(), txt.size());
    std::vector<double> buf = GetState();
    ret->Update((unsigned char*)&buf[0], sizeof(double)*buf.size());
    ret->Final();
    return ret;
  }

  //---------------------------------------------------------------------------
  OscCalculatorSterile::~OscCalculatorSterile()
  {
    delete fPMNS_Sterile;
  }

  //---------------------------------------------------------------------------
  IOscCalculatorAdjustable* OscCalculatorSterile::Copy() const
  {
    return new OscCalculatorSterile(*this);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetNFlavors(int nflavors)
  {
    fDirty = true;
    delete fPMNS_Sterile;
    fNFlavors = nflavors;
    fPMNS_Sterile = new PMNS_Sterile(fNFlavors);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetAngle(int i, int j, double th) 
  {
    fDirty = true; 
    fPMNS_Sterile->SetAngle(i, j, th);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetDelta(int i, int j, double delta)
  {
    fDirty = true; 
    fPMNS_Sterile->SetDelta(i, j, delta);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetDm(int i, double dm) 
  {
    fDirty = true; 
    fPMNS_Sterile->SetDm(i, dm);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetDmsq21(double dmsq21)
  {
    std::cerr << "Must use SetDm!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetDmsq32(double dms32)
  {
    std::cerr << "Must use SetDm!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetTh12(double th12)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetTh13(double th13)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetTh23(double th23)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterile::SetdCP(double dCP)
  {
    std::cerr << "Must use SetDelta!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  double OscCalculatorSterile::P(int flavBefore, int flavAfter, double E)
  {
    const int anti = (flavBefore > 0) ? +1 : -1;
    //anti must be +/- 1 but flavAfter can be zero
    assert(flavAfter/anti >= 0);
    if (anti != fPrevAnti) fDirty = true;
    fPrevAnti = anti;

    if (flavBefore != fPrevFlavBefore) fDirty = true;
    fPrevFlavBefore = flavBefore;

    if (abs(fPrevE - E) > 1e-8) fDirty = true;
    fPrevE = E;
    
    int i = -1, j = -1;
    if(abs(flavBefore) == 12) i = 0;
    if(abs(flavBefore) == 14) i = 1;
    if(abs(flavBefore) == 16) i = 2;
    if(abs(flavAfter) == 12) j = 0;
    if(abs(flavAfter) == 14) j = 1;
    if(abs(flavAfter) == 16) j = 2;
    if(abs(flavAfter) == 0)  j = 3;
    assert(i >= 0 && j >= 0);

    if (fDirty) {
      fPMNS_Sterile->ResetToFlavour(i);      
      // Assume Z/A=0.5
      const double Ne = fRho/2;
      fPMNS_Sterile->PropMatter(fL, E, Ne, anti);
    }
    
    // Return the active fraction
    if (j == 3) return fPMNS_Sterile->P(0) + fPMNS_Sterile->P(1) + fPMNS_Sterile->P(2);
    else        return fPMNS_Sterile->P(j);
  }

  //---------------------------------------------------------------------------
  OscCalculatorSterileTrivial::OscCalculatorSterileTrivial()
    : OscCalculatorSterile()
    {}

  //---------------------------------------------------------------------------
  OscCalculatorSterileTrivial::OscCalculatorSterileTrivial(
    const OscCalculatorSterile& calc)
    : OscCalculatorSterileTrivial()
  {
    std::vector<double> state = calc.GetState();
    SetState(state);
  }

  //---------------------------------------------------------------------------
  OscCalculatorSterileTrivial::OscCalculatorSterileTrivial(
    const OscCalculatorSterileTrivial& calc)
    : OscCalculatorSterileTrivial()
  {
    std::vector<double> state = calc.GetState();
    SetState(state);
  }

  //---------------------------------------------------------------------------
  IOscCalculatorAdjustable* OscCalculatorSterileTrivial::Copy() const
  {
    return new OscCalculatorSterileTrivial(*this);
  }

  //---------------------------------------------------------------------------
  double OscCalculatorSterileTrivial::P(
    int flavBefore, int flavAfter, double E)
  {
    return 1;
  }

  //---------------------------------------------------------------------------
  const OscCalculatorSterile* DowncastToSterile(const IOscCalculator* calc)
  {
    const OscCalculatorSterileTrivial* calc_trivial
            = dynamic_cast<const OscCalculatorSterileTrivial*>(calc);
    if (calc_trivial) return calc_trivial; 
    const OscCalculatorSterile* calc_sterile = dynamic_cast<const OscCalculatorSterile*>(calc);
    if(calc_sterile) return calc_sterile;
    else             std::cout << "Input calculator was not of type OscCalculatorSterile." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }

  //---------------------------------------------------------------------------
  OscCalculatorSterile* DowncastToSterile(IOscCalculator* calc)
  {
    OscCalculatorSterileTrivial* calc_trivial
            = dynamic_cast<OscCalculatorSterileTrivial*>(calc);
    if (calc_trivial) return calc_trivial;
    OscCalculatorSterile* calc_sterile = dynamic_cast<OscCalculatorSterile*>(calc);
    if(calc_sterile) return calc_sterile;
    else             std::cout << "Input calculator was not of type OscCalculatorSterile." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }
} // namespace
