#pragma once

#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "TDirectory.h"
#include "TString.h"
#include "TVector3.h"

class TFormLeafInfo;
class TBranch;
class TLeaf;
class TTreeFormula;
class TTree;

namespace caf
{
  class SRProxySystController
  {
  public:
    static void ResetSysts()
    {
      if(fgAnyShifted){
        ++fgSeqNo;
        fgAnyShifted = false;
      }
    }
    static bool AnyShifted() {return fgAnyShifted;}
  protected:
    template<class T> friend class Proxy;
    friend class VectorProxyBase;

    static long CurrentSeqNo() {return fgSeqNo;};
    static void SetShifted() {fgAnyShifted = true;}

    static long fgSeqNo;
    static bool fgAnyShifted;
  };

  class SRBranchRegistry
  {
  public:
    static void AddBranch(const std::string& b){fgBranches.insert(b);}
    static const std::set<std::string>& GetBranches(){return fgBranches;}
    static void clear() {fgBranches.clear();}

    static void Print(bool abbrev = true);
    static void ToFile(const std::string& fname);
  protected:
    static std::set<std::string> fgBranches;
  };

  template<class T> class Proxy
  {
  public:
    Proxy(TDirectory* d, TTree* tr, const std::string& name, const long& base, int offset);

    // Need to be copyable because Vars return us directly
    Proxy(const Proxy&);
    Proxy(const Proxy&&);
    // No need to be assignable though
    Proxy& operator=(const Proxy&) = delete;

    // Somehow including this helps us not get automatically converted to a
    // type we might not want to be in ternary expressions (we now get a type
    // error instead).
    Proxy(T v) = delete;

    ~Proxy();

    operator T() const {return GetValue();}

    bool operator==(const T& x) const {return T(*this) == x;}

    T GetValue() const;

    // In practice these are the only operations that systematic shifts use
    Proxy<T>& operator=(T x);
    Proxy<T>& operator+=(T x);
    Proxy<T>& operator*=(T x);

    std::string Name() const {return fName;}

  protected:
    T GetValueFlat() const;
    T GetValueNested() const;

    void SetShifted();

    // The type to fetch from the TLeaf - get template errors inside of ROOT
    // for enums.
    typedef typename std::conditional<std::is_enum<T>::value, int, T>::type U;

    // Shared
    std::string fName;
    mutable TLeaf* fLeaf;
    mutable T fVal;
    TTree* fTree;

    // Flat
    TDirectory* fDir;
    const long& fBase;
    int fOffset;

    // Nested
    mutable TFormLeafInfo* fLeafInfo;
    mutable TBranch* fBranch;
    mutable TTreeFormula* fTTF;
    mutable long fEntry;
    mutable int fSubIdx;

    // Syst
    T fSystOverrideValue;
    mutable long fSystOverrideEntry;
    mutable long fSystOverrideSeqNo;
  };

  // Helper functions that don't need to be templated
  class VectorProxyBase
  {
  public:
    VectorProxyBase(TDirectory* d, TTree* tr, const std::string& name,
                    const long& base, int offset);

    std::string Name() const {return fName;}

    size_t size() const;
    bool empty() const;
    void resize(size_t i);
  protected:
    void CheckIndex(size_t i) const;

    // Used by nested variant
    std::string AtSize() const;
    std::string Subscript(int i) const;

    TDirectory* fDir;
    TTree* fTree;
    std::string fName;
    const long& fBase;
    int fOffset;
    Proxy<int> fSize;
    Proxy<long long> fIdxP;
    mutable long fIdx;

    size_t fSystOverrideSize;
    mutable long fSystOverrideEntry;
    mutable long fSystOverrideSeqNo;
  };

  template<class T> class VectorProxy: public VectorProxyBase
  {
  public:
    VectorProxy(TDirectory* d, TTree* tr, const std::string& name, const long& base, int offset)
      : VectorProxyBase(d, tr, name, base, offset)
    {
    }

    ~VectorProxy(){for(T* e: fElems) delete e;}

    T& at(size_t i) const {EnsureSize(i); return *fElems[i];}
    T& at(size_t i)       {EnsureSize(i); return *fElems[i];}

    T& operator[](size_t i) const {return at(i);}
    T& operator[](size_t i)       {return at(i);}

    // U should be either T or const T
    template<class U> class iterator
    {
    public:
      U& operator*() {return (*fParent)[fIdx];}
      iterator<U>& operator++(){++fIdx; return *this;}
      bool operator!=(const iterator<U>& it) const {return fIdx != it.fIdx;}
      bool operator==(const iterator<U>& it) const {return fIdx == it.fIdx;}
    protected:
      friend class VectorProxy;
      iterator(const VectorProxy<T>* p, int i) : fParent(p), fIdx(i) {}

      const VectorProxy<T>* fParent;
      size_t fIdx;
    };

    iterator<const T> begin() const {return iterator<const T>(this, 0);}
    iterator<T> begin() {return iterator<T>(this, 0);}
    iterator<const T> end() const {return iterator<const T>(this, size());}
    iterator<T> end() {return iterator<T>(this, size());}

  protected:
    /// Implies CheckIndex()
    void EnsureSize(size_t i) const
    {
      CheckIndex(i);
      if(i >= fElems.size()) fElems.resize(i+1);

      if(fDir){
        // Flat
        fIdx = fIdxP; // store into an actual value we can point to
        if(!fElems[i]){
          TTree* tr = (TTree*)fDir->Get(fName.c_str());
          if(!tr){
            std::cout << "Couldn't find TTree " << fName << " in " << fDir->GetName() << std::endl;
            abort();
          }
          fElems[i] = new T(fDir, tr, fName, fIdx, i);
        }
      }
      else{
        // Nested
        if(!fElems[i]) fElems[i] = new T(0, fTree, Subscript(i), 0, 0);
      }
    }

    mutable std::vector<T*> fElems;
  };

  // class TVector3Proxy
  // {
  // public:
  //   TVector3Proxy(TDirectory* d, TTree* tr, const std::string& name, const long& base, int offset);

  //   float X() const {return x;}
  //   float Y() const {return y;}
  //   float Z() const {return z;}

  //   operator TVector3() const {return TVector3(x, y, z);}

  //   float Mag2() const {return x*x+y*y+z*z;}
  //   float Mag() const {return sqrt(Mag2());}
  //   float Dot(const TVector3Proxy& v) const {return x*v.x + y*v.y + z*v.z;}
  //   float Dot(const TVector3& v) const {return x*v.X() + y*v.Y() + z*v.Z();}
  //   TVector3 Unit() const
  //   {
  //     const float m = Mag();
  //     return TVector3(x/m, y/m, z/m);
  //   }

  //   Proxy<double> x, y, z;
  // };

  /// Used in comparison of GENIE version numbers
  template<class T> bool operator<(const VectorProxy<Proxy<T>>& a,
                                   const std::vector<T>& b)
  {
    const size_t N = a.size();
    if(N != b.size()) return N < b.size();
    for(size_t i = 0; i < N; ++i){
      if(a[i] != b[i]) return a[i] < b[i];
    }
    return false;
  }

  template <class T, unsigned int N> class ArrayProxy
  {
  public:
    ArrayProxy(TDirectory* d, TTree* tr, const std::string& name, const long& base, int offset)
      : fIdxP(Proxy<long long>(d, tr, name+"_idx", base, offset))
    {
      fElems.reserve(N);
      for(unsigned int i = 0; i < N; ++i){
        if(d){
          // Flat
          fElems.emplace_back(d, tr, name, fIdx, i);
        }
        else{
          // Nested
          fElems.emplace_back(nullptr, tr, TString::Format("%s[%d]", name.c_str(), i).Data(), 0, 0);
        }
      }
    }

    const Proxy<T>& operator[](size_t i) const {fIdx = fIdxP; return fElems[i];}
          Proxy<T>& operator[](size_t i)       {fIdx = fIdxP; return fElems[i];}
  protected:
    std::vector<Proxy<T>> fElems;

    // Flat
    Proxy<long long> fIdxP;
    mutable long fIdx;
  };

} // namespace

namespace std
{
  template<class T> T min(const caf::Proxy<T>& a, T b)
  {
    return std::min(a.GetValue(), b);
  }

  template<class T> T min(T a, const caf::Proxy<T>& b)
  {
    return std::min(a, b.GetValue());
  }

  template<class T> T max(const caf::Proxy<T>& a, T b)
  {
    return std::max(a.GetValue(), b);
  }

  template<class T> T max(T a, const caf::Proxy<T>& b)
  {
    return std::max(a, b.GetValue());
  }
}
