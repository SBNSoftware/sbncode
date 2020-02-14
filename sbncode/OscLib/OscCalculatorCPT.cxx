#include "OscLib/OscCalculatorCPT.h"

#include <cassert>
#include <cstdlib>

namespace osc
{

  OscCalculatorCPT::OscCalculatorCPT()
  {
    fCalc    = new OscCalculatorPMNSOpt;
    fBarCalc = new OscCalculatorPMNSOpt;
    fSigDel  = {};
  }

  OscCalculatorCPT::OscCalculatorCPT(IOscCalculatorAdjustable* calc,
                                     IOscCalculatorAdjustable* barcalc,
                                     SDMap sigdel)
  {
    fCalc    = calc;
    fBarCalc = barcalc;
    fSigDel  = sigdel;
  }

  OscCalculatorCPT::~OscCalculatorCPT()
  {
    delete fCalc;
    delete fBarCalc;
  }

  IOscCalculatorAdjustable* OscCalculatorCPT::Copy() const
  {
    return new OscCalculatorCPT(fCalc->Copy(),
                                fBarCalc->Copy(),
                                fSigDel);
  }

  double OscCalculatorCPT::P(int flavBefore, int flavAfter, double E)
  {
    assert(flavBefore*flavAfter > 0); // check for matter<->anti-matter
    return ( (flavBefore > 0) ?    fCalc->P(flavBefore, flavAfter, E)
                              : fBarCalc->P(flavBefore, flavAfter, E) );
  }

  // asymmetric setters
  void OscCalculatorCPT::SetL(double L, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetL(L) :
                           fBarCalc->SetL(L) ;
  }

  void OscCalculatorCPT::SetRho(double rho, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetRho(rho) :
                           fBarCalc->SetRho(rho) ;
  }

  void OscCalculatorCPT::SetDmsq21(double dmsq21, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetDmsq21(dmsq21) :
                           fBarCalc->SetDmsq21(dmsq21) ;
  }

  void OscCalculatorCPT::SetDmsq32(double dmsq32, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetDmsq32(dmsq32) :
                           fBarCalc->SetDmsq32(dmsq32) ;
  }

  void OscCalculatorCPT::SetTh12(double th12, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetTh12(th12) :
                           fBarCalc->SetTh12(th12) ;
  }

  void OscCalculatorCPT::SetTh13(double th13, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetTh13(th13) :
                           fBarCalc->SetTh13(th13) ;
  }

  void OscCalculatorCPT::SetTh23(double th23, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetTh23(th23) :
                           fBarCalc->SetTh23(th23) ;
  }

  void OscCalculatorCPT::SetdCP(double dCP, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetdCP(dCP) :
                           fBarCalc->SetdCP(dCP) ;
  }

  // symmetric getters
  double OscCalculatorCPT::GetL() const
  {
    double L = GetL(ENuSign::kNu);
    assert( L == GetL(ENuSign::kNuBar) );
    return L;
  }

  double OscCalculatorCPT::GetRho() const
  {
    double rho = GetRho(ENuSign::kNu);
    assert( rho == GetRho(ENuSign::kNuBar) );
    return rho;
  }

  double OscCalculatorCPT::GetDmsq21() const
  {
    double dmsq21 = GetDmsq21(ENuSign::kNu);
    assert( dmsq21 == GetDmsq21(ENuSign::kNuBar) );
    return dmsq21;
  }

  double OscCalculatorCPT::GetDmsq32() const
  {
    double dmsq32 = GetDmsq32(ENuSign::kNu);
    assert( dmsq32 == GetDmsq32(ENuSign::kNuBar) );
    return dmsq32;
  }

  double OscCalculatorCPT::GetTh12() const
  {
    double th12 = GetTh12(ENuSign::kNu);
    assert( th12 == GetTh12(ENuSign::kNuBar) );
    return th12;
  }

  double OscCalculatorCPT::GetTh13() const
  {
    double th13 = GetTh13(ENuSign::kNu);
    assert( th13 == GetTh13(ENuSign::kNuBar) );
    return th13;
  }

  double OscCalculatorCPT::GetTh23() const
  {
    double th23 = GetTh23(ENuSign::kNu);
    assert( th23 == GetTh23(ENuSign::kNuBar) );
    return th23;
  }

  double OscCalculatorCPT::GetdCP() const
  {
    double dCP = GetdCP(ENuSign::kNu);
    assert( dCP == GetdCP(ENuSign::kNuBar) );
    return dCP;
  }

  // asymmetric getters
  double OscCalculatorCPT::GetL(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetL()
                                  : fBarCalc->GetL() );
  }

  double OscCalculatorCPT::GetRho(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetRho()
                                  : fBarCalc->GetRho() );
  }

  double OscCalculatorCPT::GetDmsq21(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetDmsq21()
                                  : fBarCalc->GetDmsq21() );
  }

  double OscCalculatorCPT::GetDmsq32(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetDmsq32()
                                  : fBarCalc->GetDmsq32() );
  }

  double OscCalculatorCPT::GetTh12(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetTh12()
                                  : fBarCalc->GetTh12() );
  }

  double OscCalculatorCPT::GetTh13(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetTh13()
                                  : fBarCalc->GetTh13() );
  }

  double OscCalculatorCPT::GetTh23(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetTh23()
                                  : fBarCalc->GetTh23() );
  }

  double OscCalculatorCPT::GetdCP(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetdCP()
                                  : fBarCalc->GetdCP() );
  }


  TMD5* OscCalculatorCPT::GetParamsHash() const
  {
    TMD5* hash    = fCalc->GetParamsHash();
    TMD5* barhash = fBarCalc->GetParamsHash();

    // don't provide hash unless both sub-cals do
    if ( !(hash && barhash) )
    {
      delete hash;
      delete barhash;
      return 0;
    }

    TMD5* ret = new TMD5;

    //hash together sub-hashes
    ret->Update( (unsigned char*)hash->AsString(), 16);
    ret->Update( (unsigned char*)barhash->AsString(), 16);

    delete hash;
    delete barhash;

    //also hash in class name in case another CPT class has same idea
    ret->Update( (unsigned char*)"OscCalculatorCPT", 16);

    ret->Final();
    return ret;
  } 

  //---------------------------------------------------------------------------

  const OscCalculatorCPT* DowncastToCPT(const IOscCalculatorAdjustable* osc)
  {
    const OscCalculatorCPT* cpt = dynamic_cast<const OscCalculatorCPT*>(osc);
    assert( cpt && "Must use OscCalculatorCPT with CPT FitVars." );
    return cpt;
  }

  OscCalculatorCPT* DowncastToCPT(IOscCalculatorAdjustable* osc)
  {
    OscCalculatorCPT* cpt = dynamic_cast<OscCalculatorCPT*>(osc);
    assert( cpt && "Must use OscCalculatorCPT with CPT FitVars." );
    return cpt;
  }

} // namespace
