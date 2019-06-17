#include "StandardRecord/Proxy/BasicTypesProxy.h"
#include "StandardRecord/SREnums.h"

#include "TFormLeafInfo.h"
#include "TTreeFormula.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

namespace caf
{
  long SRProxySystController::fgSeqNo = 0;
  bool SRProxySystController::fgAnyShifted = false;

  std::set<std::string> SRBranchRegistry::fgBranches;

  //----------------------------------------------------------------------
  void SRBranchRegistry::Print(bool abbrev)
  {
    std::string prev;
    for(std::string b: fgBranches){
      if(abbrev){
        unsigned int cutto = 0;
        for(unsigned int i = 0; i < std::min(b.size(), prev.size()); ++i){
          if(b[i] != prev[i]) break;
          if(b[i] == '.') cutto = i;
        }
        prev = b;
        for(unsigned int i = 0; i < cutto; ++i) b[i] = ' ';
      }
      std::cout << b << std::endl;
    }
  }

  //----------------------------------------------------------------------
  void SRBranchRegistry::ToFile(const std::string& fname)
  {
    std::ofstream fout(fname);
    for(const std::string& b: fgBranches) fout << b << std::endl;
  }

  //----------------------------------------------------------------------
  template<class T> Proxy<T>::Proxy(TDirectory* d, TTree* tr, const std::string& name, const long& base, int offset)
    : fName(name), fLeaf(0), fTree(tr),
      fDir(d), fBase(base), fOffset(offset),
      fLeafInfo(0), fBranch(0), fTTF(0), fEntry(-1), fSubIdx(0),
      fSystOverrideEntry(-1), fSystOverrideSeqNo(-1)
  {
  }

  //----------------------------------------------------------------------
  template<class T> Proxy<T>::Proxy(const Proxy<T>& p)
    : fName("copy of "+p.fName), fLeaf(0), fTree(p.fTree),
      fDir(0), fBase(p.fBase), fOffset(p.fOffset),
      fLeafInfo(0), fBranch(0), fTTF(0), fSubIdx(-1),
      fSystOverrideValue(p.fSystOverrideValue),
      fSystOverrideEntry(p.fSystOverrideEntry),
      fSystOverrideSeqNo(p.fSystOverrideSeqNo)
  {
    // Ensure that the value is evaluated and baked in in the parent object, so
    // that fTTF et al aren't re-evaluated in every single copy.
    fVal = p.GetValue();
    fEntry = p.fEntry;
  }

  //----------------------------------------------------------------------
  template<class T> Proxy<T>::Proxy(const Proxy&& p)
    : fName("move of "+p.fName), fLeaf(0), fTree(p.fTree),
      fDir(0), fBase(p.fBase), fOffset(p.fOffset),
      fLeafInfo(0), fBranch(0), fTTF(0), fSubIdx(-1),
      fSystOverrideValue(p.fSystOverrideValue),
      fSystOverrideEntry(p.fSystOverrideEntry),
      fSystOverrideSeqNo(p.fSystOverrideSeqNo)
  {
    // Ensure that the value is evaluated and baked in in the parent object, so
    // that fTTF et al aren't re-evaluated in every single copy.
    fVal = p.GetValue();
    fEntry = p.fEntry;
  }

  //----------------------------------------------------------------------
  template<class T> Proxy<T>::~Proxy()
  {
    // The other pointers aren't ours
    delete fTTF;
  }

  //----------------------------------------------------------------------
  template<class T> T Proxy<T>::GetValue() const
  {
    if(fDir) return GetValueFlat(); else return GetValueNested();
  }

  template<class T> void GetTypedValueWrapper(TLeaf* leaf, T& x, int subidx)
  {
    x = leaf->GetTypedValue<T>(subidx);
  }

  void GetTypedValueWrapper(TLeaf* leaf, std::string& x, int subidx)
  {
    std::cout << "Reading string from flat tree unsupported" << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  template<class T> T Proxy<T>::GetValueFlat() const
  {
    // If there's a valid systematic override value in place, give that
    if((!fTree || fSystOverrideEntry == fBase+fOffset) &&
       fSystOverrideSeqNo == SRProxySystController::CurrentSeqNo()){
      return fSystOverrideValue;
    }

    if(!fTree) return fVal; // copied or moved?

    if(!fLeaf){
      fLeaf = fTree->GetLeaf(fName.c_str());
      if(!fLeaf){
        std::cout << std::endl << "BasicTypeProxy: Branch '" << fName
                  << "' not found in tree '" << fTree->GetName() << "'."
                  << std::endl;
        abort();
      }

      fEntry = -1;
      fSystOverrideSeqNo = -1;

      if(fName.find("_idx") == std::string::npos &&
         fName.find("_length") == std::string::npos &&
         fName.find(".size()") == std::string::npos){ // specific to "nested"
        SRBranchRegistry::AddBranch(fName);
      }
    }

    if(fEntry == fBase+fOffset) return fVal; // cache up-to-date

    fLeaf->GetBranch()->GetEntry(fBase+fOffset);

    //    fVal = (T)fLeaf->GetTypedValue<U>(0);
    U tmp;
    GetTypedValueWrapper(fLeaf, tmp, 0);
    fVal = (T)tmp;

    fEntry = fBase+fOffset;
    fSystOverrideSeqNo = -1;

    return fVal;
  }

  template<class T> void EvalInstanceWrapper(TTreeFormula* ttf, T& x)
  {
    // TODO is this the safest way to cast?
    x = (T)ttf->EvalInstance(0);
  }

  void EvalInstanceWrapper(TTreeFormula* ttf, std::string& x)
  {
    x = ttf->EvalStringInstance(0);
  }

  //----------------------------------------------------------------------
  template<class T> T Proxy<T>::GetValueNested() const
  {
    // If there's a valid systematic override value in place, give that
    if((!fTree || fSystOverrideEntry == fTree->GetReadEntry()) &&
       fSystOverrideSeqNo == SRProxySystController::CurrentSeqNo()){
      return fSystOverrideValue;
    }

    // If there's no tree set but we've reached this point, either we were
    // constructed wrong, or CopyRecord() hasn't initialized our SystOverride
    // properly.
    assert(fTree);

    // If we have the value cached already, give that
    if(fEntry == fTree->GetReadEntry()) return fVal;

    // We're about to update to the current entry, and we won't be appropriate
    // for any systematic shift.
    fEntry = fTree->GetReadEntry();
    fSystOverrideSeqNo = -1;

    // First time calling, set up the branches etc
    if(!fTTF){
      SRBranchRegistry::AddBranch(fName);

      // Leaves are attached to the TTF, must keep it
      fTTF = new TTreeFormula(("TTFProxy-"+fName).c_str(), fName.c_str(), fTree);
      fLeafInfo = fTTF->GetLeafInfo(0); // Can fail (for a regular branch?)
      fLeaf = fTTF->GetLeaf(0);
      fBranch = fLeaf->GetBranch();

      if(!fLeaf || !fBranch){
        // If it ends in .x or .y or .z it's probably a SRVector3D. And maybe
        // we're running on one of the weird MR files where some of these are
        // still TVector3. Make an attempt to patch things up.
        if(fName.size() > 2 && fName[fName.size()-2] == '.'){
          const char lastChar = fName[fName.size()-1];
          if(lastChar >= 'x' && lastChar <= 'z'){
            // ".x" -> ".fX" etc
            std::string name = fName;
            name.resize(name.size()-1);
            name += 'f';
            name += std::toupper(lastChar);

            std::cout << "Branch " << fName << " not found. "
                      << "Trying " << name << " instead." << std::endl;

            delete fTTF;
            delete fLeaf;
            delete fBranch;

            fTTF = new TTreeFormula(("TTFProxy-"+name).c_str(), name.c_str(), fTree);
            fLeafInfo = fTTF->GetLeafInfo(0);
            fLeaf = fTTF->GetLeaf(0);
            fBranch = fLeaf->GetBranch();
          }
        }

        // If we still didn't find the leaf or branch, abort
        if(!fLeaf || !fBranch){
          std::cout << "Couldn't find " << fName << " in tree. Abort."
                    << std::endl;

          abort();
        }
      }

      // TODO - parsing the array indices out sucks - pass in as an int somehow
      const size_t open_idx = fName.find_first_of('[');
      // Do we have exactly one set of [] in the name?
      if(open_idx != std::string::npos && open_idx == fName.find_last_of('[')){
	const size_t close_idx = fName.find_first_of(']');

	std::string numPart = fName.substr(open_idx+1, close_idx-open_idx-1);
	fSubIdx = atoi(numPart.c_str());
      }
    }

    if(fLeafInfo){
      // Using TTreeFormula always works, and is sometimes necessary

      fTTF->GetNdata(); // for some reason this is necessary for fTTF to work
                        // in all cases.

      EvalInstanceWrapper(fTTF, fVal);
    }
    else{
      // But when this is possible the hope is it might be faster

      if(fBranch->GetReadEntry() != fEntry){
        fBranch->GetEntry(fEntry);
      }

      // This check is much quicker than what CheckIndex() does, which winds up
      // calling a TTF, but I can't figure out a safe way to automatically
      // elide that check.
      if(fSubIdx > fLeaf->GetLen()){
        std::cout << std::endl << fName << " out of range (size() == " << fLeaf->GetLen() << "). Aborting." << std::endl;
        abort();
      }

      //      fVal = (T)fLeaf->GetTypedValue<U>(fSubIdx);
      U tmp;
      GetTypedValueWrapper(fLeaf, tmp, fSubIdx);
      fVal = (T)tmp;
    }

    return fVal;
  }

  //----------------------------------------------------------------------
  template<class T> void Proxy<T>::SetShifted()
  {
    if(fDir){
      // Flat
      fSystOverrideEntry = fBase+fOffset;
    }
    else{
      // Nested
      fSystOverrideEntry = fTree ? fTree->GetReadEntry() : 0;
    }
    fSystOverrideSeqNo = SRProxySystController::CurrentSeqNo();
    SRProxySystController::SetShifted();
  }

  //----------------------------------------------------------------------
  template<class T> Proxy<T>& Proxy<T>::operator=(T x)
  {
    fSystOverrideValue = x;
    SetShifted();
    return *this;
  }

  //----------------------------------------------------------------------
  template<class T> Proxy<T>& Proxy<T>::operator+=(T x)
  {
    fSystOverrideValue = T(GetValue() + x);
    SetShifted();
    return *this;
  }

  //----------------------------------------------------------------------
  template<class T> Proxy<T>& Proxy<T>::operator*=(T x)
  {
    fSystOverrideValue = T(GetValue() * x);
    SetShifted();
    return *this;
  }

  //----------------------------------------------------------------------
  template<> Proxy<std::string>& Proxy<std::string>::operator*=(std::string x)
  {
    std::cout << "BasicTypesProxy.cxx: Multiplying strings makes no sense..." << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  VectorProxyBase::VectorProxyBase(TDirectory* d, TTree* tr,
                                   const std::string& name,
                                   const long& base, int offset)
    : fDir(d), fTree(tr), fName(name), fBase(base), fOffset(offset), fSize(d, tr, fDir ? fName+"_length" : AtSize(), base, offset), fIdxP(d, tr, name+"_idx", base, offset),
      fSystOverrideSize(-1),
      fSystOverrideEntry(-1),
      fSystOverrideSeqNo(-1)
  {
  }

  //----------------------------------------------------------------------
  void VectorProxyBase::CheckIndex(size_t i) const
  {
    // This is the only way to get good error messages. But it also seems like
    // the call to size() here is necessary in Nested mode to trigger some
    // side-effect within ROOT, otherwise we get some bogus index out-of-range
    // crashes.
    if(i >= size()){
      std::cout << std::endl << fName << "[" << i << "] out of range (size() == " << size() << "). Aborting." << std::endl;
      abort();
    }
  }

  //----------------------------------------------------------------------
  // Used by nested variant
  std::string VectorProxyBase::AtSize() const
  {
    // These fields don't have a corresponding "n" so we have no choice
    static const std::unordered_set<std::string> missing = {
      // None for SBN for now!
    };

    const int idx = fName.find_last_of('.');

    if(missing.count(fName.substr(idx+1))){
      // foo.bar.baz -> foo.bar.@baz.size()
      return fName.substr(0, idx+1)+"@"+fName.substr(idx+1)+".size()";
    }
    else{
      // foo.bar.baz -> foo.bar.nbaz
      return fName.substr(0, idx)+".n"+fName.substr(idx+1);
    }
  }

  //----------------------------------------------------------------------
  // Used by nested variant
  std::string VectorProxyBase::Subscript(int i) const
  {
    // Only have to do the at() business for subscripts from the 3rd one on
    const int nSubs = std::count(fName.begin(), fName.end(), '[');
    if(nSubs < 2) return TString::Format("%s[%d]", fName.c_str(), i).Data();

    const int idx = fName.find_last_of('.');
    return TString::Format("%s.@%s.at(%d)",
                           fName.substr(0, idx).c_str(),
                           fName.substr(idx+1).c_str(),
                           i).Data();
  }

  //----------------------------------------------------------------------
  size_t VectorProxyBase::size() const
  {
    return 1234; // HACK HACK HACK

    // If there's a valid systematic override value in place, give that
    if(fDir){
      // Flat
      if((!fTree || fSystOverrideEntry == fBase+fOffset) &&
         fSystOverrideSeqNo == SRProxySystController::CurrentSeqNo()){
        return fSystOverrideSize;
      }
    }
    else{
      // Nested
      if((!fTree || fSystOverrideEntry == fTree->GetReadEntry()) &&
         fSystOverrideSeqNo == SRProxySystController::CurrentSeqNo()){
        return fSystOverrideSize;
      }
    }

    return fSize;
  }

  //----------------------------------------------------------------------
  bool VectorProxyBase::empty() const
  {
    return size() == 0;
  }

  //----------------------------------------------------------------------
  void VectorProxyBase::resize(size_t i)
  {
    fSystOverrideSize = i;
    if(fDir){
      // Flat
      fSystOverrideEntry = fBase+fOffset;
    }
    else{
      // Nested
      fSystOverrideEntry = fTree ? fTree->GetReadEntry() : 0;
    }
    fSystOverrideSeqNo = SRProxySystController::CurrentSeqNo();
    SRProxySystController::SetShifted();
  }

  //----------------------------------------------------------------------
  TVector3Proxy::TVector3Proxy(TDirectory* d, TTree* tr, const std::string& name, const long& base, int offset)
    : x(d, tr, name+".fX", base, offset),
      y(d, tr, name+".fY", base, offset),
      z(d, tr, name+".fZ", base, offset)
  {
  }


  // Enumerate all the variants we expect
  template class Proxy<int>;
  template class Proxy<long>;
  template class Proxy<unsigned int>;
  template class Proxy<unsigned short>;
  template class Proxy<unsigned long>;
  template class Proxy<long long>;
  template class Proxy<unsigned long long>;
  template class Proxy<short>;
  template class Proxy<float>;
  template class Proxy<double>;
  template class Proxy<bool>;
  template class Proxy<unsigned char>;

  template class Proxy<std::string>;

  template class Proxy<SRExperiment>;
  //  template class Proxy<View_t>;
  //  template class Proxy<generator_>;
  //  template class Proxy<mode_type_>;
  //  template class Proxy<int_type_>;

} // namespace
