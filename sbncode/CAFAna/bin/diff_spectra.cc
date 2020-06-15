#include <cassert>
#include <set>
#include <iostream>

#include "TClass.h"
#include "TFile.h"
#include "TH1.h"
#include "TObject.h"
#include "TObjString.h"
#include "TPad.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TVector3.h"

// Set by command-line switch
bool gEnablePdf = true;

// Has the output pdf been started?
bool gPdfStarted = false;

// ----------------------------------------------------------------------------
std::string ConcatPath(const std::string& a, const std::string& b)
{
  if(a.empty()) return b;
  return a+"/"+b;
}

// ----------------------------------------------------------------------------
bool operator==(const TObjString& a, const TObjString& b)
{
  return a.GetString() == b.GetString();
}

bool operator==(const TVectorD& a, const TVectorD& b)
{
  return (a-b).Norm1() == 0;
}

bool operator==(const TH1& a, const TH1& b)
{
  if(a.GetNcells() != b.GetNcells()) return false;
  // TODO - test bin edges
  for(int i = 0; i <= a.GetNcells(); ++i)
    if(a.GetBinContent(i) != b.GetBinContent(i))
      return false;

  return true;
}

// ----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const TObjString& a)
{
  return os << "\"" << a.GetString() << "\"";
}

std::ostream& operator<<(std::ostream& os, const TH1& a)
{
  return os << "integral " << a.Integral(0, -1);
}

std::ostream& operator<<(std::ostream& os, const TObject& a)
{
  a.Print();
  return os;
}

// ----------------------------------------------------------------------------
void Plots(const TH1& a, const TH1& b, const std::string& key)
{
  TH1* ac = (TH1*)a.Clone();
  TH1* bc = (TH1*)b.Clone();

  ac->SetLineColor(kRed);
  bc->SetLineColor(kBlue);

  if(ac->Integral() < bc->Integral()) std::swap(ac, bc);

  ac->SetTitle(key.c_str());

  ac->Draw("hist");
  bc->Draw("hist same");

  if(!gPdfStarted){
    gPdfStarted = true;
    gPad->Print("diff.pdf[");
  }

  gPad->Print("diff.pdf");
}

// ----------------------------------------------------------------------------
int gNChecked = 0;
template<class T> bool Compare(const T& a, const T& b, const std::string& key)
{
  if(a == b){
    std::cout << "\rChecked " << ++gNChecked << " keys" << std::flush;
    return true;
  }

  std::cout << "\n'" << key << "' differs: " << a << " vs " << b << std::endl;

  if constexpr(std::is_same_v<T, TH1>)
    if(gEnablePdf)
      Plots(a, b, key);

  return false;
}

// ----------------------------------------------------------------------------
bool CompareObjects(const TObject& a, const TObject& b, const std::string& key)
{
  if(a.ClassName() != std::string(b.ClassName())){
    std::cout << key << " is of unlike type: "
              << a.ClassName() << " vs "
              << b.ClassName() << std::endl;
    return false;
  }

  if(a.ClassName() == std::string("TObjString")){
    return Compare((TObjString&)a, (TObjString&)b, key);
  }

  if(a.ClassName() == std::string("TVectorT<double>")){
    return Compare((TVectorD&)a, (TVectorD&)b, key);
  }

  if(a.ClassName() == std::string("TVector3")){
    return Compare((TVector3&)a, (TVector3&)b, key);
  }

  if(a.InheritsFrom(TH1::Class()) && b.InheritsFrom(TH1::Class())){
    return Compare((TH1&)a, (TH1&)b, key);
  }

  std::cout << "Don't know how to compare '" << key << "' of type "
            << a.ClassName() << std::endl;
  return true; // give them the benefit of the doubt
}

// ----------------------------------------------------------------------------
void ReportUnique(const std::vector<std::string>& only, TDirectory* dir)
{
  if(only.empty()) return;

  std::cout << "The following keys exist only in " << dir->GetFile()->GetName() << std::endl;
  for(const std::string& k: only) std::cout << "  " << k << std::endl;
}

// ----------------------------------------------------------------------------
bool CheckDirectory(TDirectory* a, TDirectory* b, std::string path = "")
{
  std::set<std::string> keysA, keysB;
  TIter nextA(a->GetListOfKeys()), nextB(b->GetListOfKeys());
  while(TObject* x = nextA()) keysA.insert(x->GetName());
  while(TObject* x = nextB()) keysB.insert(x->GetName());

  std::vector<std::string> onlyA, onlyB, bothAB;
  std::set_difference  (keysA.begin(), keysA.end(), keysB.begin(), keysB.end(), std::back_inserter(onlyA));
  std::set_difference  (keysB.begin(), keysB.end(), keysA.begin(), keysA.end(), std::back_inserter(onlyB));
  std::set_intersection(keysA.begin(), keysA.end(), keysB.begin(), keysB.end(), std::back_inserter(bothAB));

  ReportUnique(onlyA, a);
  ReportUnique(onlyB, b);

  bool ok = onlyA.empty() && onlyB.empty();

  for(const std::string& key: bothAB){
    TObject* objA = a->Get(key.c_str());
    TObject* objB = b->Get(key.c_str());

    TDirectory* dirA = dynamic_cast<TDirectoryFile*>(objA);
    TDirectory* dirB = dynamic_cast<TDirectoryFile*>(objB);
    if(dirA && dirB){
      if(!CheckDirectory(dirA, dirB, ConcatPath(path, key))) ok = false;
    }
    else{
      if(!CompareObjects(*objA, *objB, ConcatPath(path, key))) ok = false;
    }

    delete objA;
    delete objB;
  } // end for key

  return ok;
}


// ----------------------------------------------------------------------------
void usage()
{
  std::cout << "Usage: diff_spectra [-np] input1.root input2.root" << std::endl;
  std::cout << "  -np Disable plotting functionality" << std::endl;

  exit(2);
}

// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  if(argc != 3 && argc != 4) usage();
  if(argc == 4 && argv[1] != std::string("-np")) usage();

  if(argc == 4){
    // Must have been -np
    gEnablePdf = false;
    // Shift to maintain filenames as [1] and [2]
    ++argv;
  }

  const std::string fnameA = argv[1];
  const std::string fnameB = argv[2];

  std::cout << "Comparing " << fnameA << " and " << fnameB << std::endl;

  if(CheckDirectory(TFile::Open(fnameA.c_str()),
                    TFile::Open(fnameB.c_str()))){
    std::cout << "\nFiles match!" << std::endl;
    return 0;
  }
  else{
    if(gPdfStarted) gPad->Print("diff.pdf]");

    std::cout << "\nFiles do not match" << std::endl;
    return 1;
  }
}
