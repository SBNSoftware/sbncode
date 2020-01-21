#include "CAFAna/Core/Ratio.h"

#include "CAFAna/Core/HistCache.h"
#include "CAFAna/Core/Utilities.h"

#include "TH1.h"

#include <cassert>

namespace ana
{
  //----------------------------------------------------------------------
  Ratio::Ratio(const Spectrum& num, const Spectrum& denom,
	       bool purOrEffErrs)
  {
    // Scale to same arbitrary POT
    fHist = num.ToTH1(1e20);
    TH1D* temp = denom.ToTH1(1e20);
    if(purOrEffErrs){
      fHist->Divide(fHist, temp, 1, 1, "B");
    }
    else{
      fHist->Divide(temp);
    }
    HistCache::Delete(temp);

    fHist->GetYaxis()->SetTitle("Ratio");

    // TODO: set error bars smartly
  }

  //----------------------------------------------------------------------
  Ratio::Ratio(TH1* h, std::string varName)
  {
    if(!h){
      fHist = 0;
      return;
    }

    DontAddDirectory guard;

    const TString className = h->ClassName();

    if(className == "TH1D"){
      // Shortcut if types match
      fHist = HistCache::Copy((TH1D*)h);
    }
    else{
      fHist = HistCache::New(UniqueName(), h->GetXaxis());
      fHist->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
      fHist->Add(h);
    }

    if(!varName.empty()) fHist->GetXaxis()->SetTitle(varName.c_str());
  }

  //----------------------------------------------------------------------
  Ratio::~Ratio()
  {
    HistCache::Delete(fHist);
  }

  //----------------------------------------------------------------------
  Ratio::Ratio(const Ratio& rhs)
  {
    DontAddDirectory guard;

    assert(rhs.fHist);
    fHist = HistCache::Copy(rhs.fHist);
  }

  //----------------------------------------------------------------------
  Ratio& Ratio::operator=(const Ratio& rhs)
  {
    if(this == &rhs) return *this;

    DontAddDirectory guard;

    HistCache::Delete(fHist);
    assert(rhs.fHist);
    fHist = HistCache::Copy(rhs.fHist);
    return *this;
  }

  //----------------------------------------------------------------------
  Ratio& Ratio::operator*=(const Ratio& rhs)
  {
    fHist->Multiply(rhs.fHist);
    return *this;
  }

  //----------------------------------------------------------------------
  Ratio Ratio::operator*(const Ratio& rhs) const
  {
    Ratio ret = *this;
    ret *= rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  Ratio& Ratio::operator/=(const Ratio& rhs)
  {
    fHist->Divide(rhs.fHist);
    return *this;
  }

  //----------------------------------------------------------------------
  Ratio Ratio::operator/(const Ratio& rhs) const
  {
    Ratio ret = *this;
    ret /= rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  TH1D* Ratio::ToTH1(Color_t col, Style_t style) const
  {
    // Could have a file temporarily open
    DontAddDirectory guard;

    TH1D* ret = HistCache::Copy(fHist);
    ret->SetLineColor(col);
    ret->SetLineStyle(style);
    return ret;
  }
} // namespace
