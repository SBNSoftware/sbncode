#include "CAFAna/Core/Cut.h"

#include <iostream>
#include <map>

namespace ana
{
  // Make sure all three versions get generated
  template class GenericCut<caf::StandardRecord>;
  template class GenericCut<caf::SRSpill>;
  template class GenericCut<caf::SRSpillTruthBranch>;

  template<class T> int GenericCut<T>::fgNextID = 0;

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>::
  GenericCut(const std::set<std::string>& /*reqs*/,
             const std::function<CutFunc_t>& func,
	     const std::function<ExposureFunc_t>& liveFunc,
	     const std::function<ExposureFunc_t>& potFunc)
    : fFunc(func), fLiveFunc(liveFunc), fPOTFunc(potFunc),
      fID(fgNextID++)
  {
  }

  //----------------------------------------------------------------------
  std::function<ExposureFunc_t>
  CombineExposures(const std::function<ExposureFunc_t>& a,
                   const std::function<ExposureFunc_t>& b)
  {
    if(!a && !b) return 0;
    if(!a) return b;
    if(!b) return a;

    return [a, b](const caf::SRSpill* spill){
      const double va = a(spill);
      const double vb = b(spill);

      if(va >= 0 && vb >= 0){
        std::cout << "Inconsistent pot/livetime values of "
                  << va << " and " << vb
                  << " from two cuts being combined." << std::endl;
        abort();
      }

      return std::max(va, vb);
    };
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T> operator&&(const GenericCut<T>& a,
                                             const GenericCut<T>& b)
  {
    // The same pairs of cuts are frequently and-ed together. Make sure those
    // duplicates get the same IDs by remembering what we've done in the past.
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return GenericCut<T>({},
			   [a, b](const T* sr){return a(sr) && b(sr);},
                           CombineExposures(a.fLiveFunc, b.fLiveFunc),
                           CombineExposures(a.fPOTFunc, b.fPOTFunc),
			   ids[key]);
    }
    else{
      const GenericCut<T> ret({},
                              [a, b](const T* sr){return a(sr) && b(sr);},
                              CombineExposures(a.fLiveFunc, b.fLiveFunc),
                              CombineExposures(a.fPOTFunc, b.fPOTFunc));
      ids[key] = ret.ID();
      return ret;
    }
  }

  // Make sure all versions get generated
  template Cut operator&&<caf::StandardRecord>(const Cut& a, const Cut& b);
  template SpillCut operator&&<caf::SRSpill>(const SpillCut& a, const SpillCut& b);
  template SpillTruthCut operator&&<caf::SRSpillTruthBranch>(const SpillTruthCut& a, const SpillTruthCut& b);

  //----------------------------------------------------------------------
  template<class T> GenericCut<T> operator||(const GenericCut<T>& a,
                                             const GenericCut<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());
    if(ids.count(key)){
      return GenericCut<T>({},
			   [a, b](const T* sr){return a(sr) || b(sr);},
                           CombineExposures(a.fLiveFunc, b.fLiveFunc),
                           CombineExposures(a.fPOTFunc, b.fPOTFunc),
			   ids[key]);
    }
    else{
      const GenericCut<T> ret({},
                              [a, b](const T* sr){return a(sr) || b(sr);},
                              CombineExposures(a.fLiveFunc, b.fLiveFunc),
                              CombineExposures(a.fPOTFunc, b.fPOTFunc));
      ids[key] = ret.ID();
      return ret;
    }
  }

  // Make sure all versions get generated
  template Cut operator||<caf::StandardRecord>(const Cut& a, const Cut& b);
  template SpillCut operator||<caf::SRSpill>(const SpillCut& a, const SpillCut& b);
  template SpillTruthCut operator||<caf::SRSpillTruthBranch>(const SpillTruthCut& a, const SpillTruthCut& b);

  //----------------------------------------------------------------------
  template<class T> GenericCut<T> operator!(const GenericCut<T>& a)
  {
    static std::map<int, int> ids;
    if(ids.count(a.ID())){
      return GenericCut<T>({},
			   [a](const T* sr){return !a(sr);},
			   0, 0, ids[a.ID()]);
    }
    else{
      const GenericCut<T> ret({},
                              [a](const T* sr){return !a(sr);});
      ids[a.ID()] = ret.ID();
      return ret;
    }
  }

  // Make sure all versions get generated
  template Cut operator!<caf::StandardRecord>(const Cut& a);
  template SpillCut operator!<caf::SRSpill>(const SpillCut& a);
  template SpillTruthCut operator!<caf::SRSpillTruthBranch>(const SpillTruthCut& a);


  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator>(const GenericVar<T>& v, double c)
  {
    return GenericCut<T>({},
                         [v, c](const T* sr){return v(sr) > c;});
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator>=(const GenericVar<T>& v, double c)
  {
    return GenericCut<T>({},
                         [v, c](const T* sr){return v(sr) >= c;});
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator<(const GenericVar<T>& v, double c)
  {
    return GenericCut<T>({},
                         [v, c](const T* sr){return v(sr) < c;});
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator<=(const GenericVar<T>& v, double c)
  {
    return GenericCut<T>({},
                         [v, c](const T* sr){return v(sr) <= c;});
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator==(const GenericVar<T>& v, double c)
  {
    return GenericCut<T>({},
                         [v, c](const T* sr){return v(sr) == c;});
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator!=(const GenericVar<T>& v, double c)
  {
    return !(v == c);
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator>(const GenericVar<T>& a, const GenericVar<T>& b)
  {
    return GenericCut<T>({},
                         [a, b](const T* sr){return a(sr) > b(sr);});
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator>=(const GenericVar<T>& a, const GenericVar<T>& b)
  {
    return GenericCut<T>({},
                         [a, b](const T* sr){return a(sr) >= b(sr);});
  }

  //----------------------------------------------------------------------
  template<class T> GenericCut<T>
  operator==(const GenericVar<T>& a, const GenericVar<T>& b)
  {
    return GenericCut<T>({},
                         [a, b](const T* sr){return a(sr) == b(sr);});
  }

  // Build the rest up through simple logic

  template<class T> GenericCut<T> operator>(double c, const GenericVar<T>& v){return v < c;}
  template<class T> GenericCut<T> operator<(double c, const GenericVar<T>& v){return v > c;}
  template<class T> GenericCut<T> operator>=(double c, const GenericVar<T>& v){return v <= c;}
  template<class T> GenericCut<T> operator<=(double c, const GenericVar<T>& v){return v >= c;}
  template<class T> GenericCut<T> operator!=(double c, const GenericVar<T>& v){return v != c;}

  template<class T> GenericCut<T> operator<(const GenericVar<T>& a, const GenericVar<T>& b){return !(b >= a);}
  template<class T> GenericCut<T> operator<=(const GenericVar<T>& a, const GenericVar<T>& b){return !(b > a);}
  template<class T> GenericCut<T> operator!=(const GenericVar<T>& a, const GenericVar<T>& b){return !(b == a);}


  // explicitly instantiate the templates for the types we know we have
  template Cut operator>(const Var&, double);
  template Cut operator<(const Var&, double);
  template Cut operator>=(const Var&, double);
  template Cut operator<=(const Var&, double);
  template Cut operator==(const Var&, double);
  template Cut operator!=(const Var&, double);

  template Cut operator>(const Var&, const Var&);
  template Cut operator<(const Var&, const Var&);
  template Cut operator>=(const Var&, const Var&);
  template Cut operator<=(const Var&, const Var&);
  template Cut operator==(const Var&, const Var&);
  template Cut operator!=(const Var&, const Var&);

  template Cut operator>(double, const Var&);
  template Cut operator<(double, const Var&);
  template Cut operator>=(double, const Var&);
  template Cut operator<=(double, const Var&);

  template SpillCut operator>(const SpillVar&, double);
  template SpillCut operator<(const SpillVar&, double);
  template SpillCut operator>=(const SpillVar&, double);
  template SpillCut operator<=(const SpillVar&, double);
  template SpillCut operator==(const SpillVar&, double);

  template SpillCut operator>(const SpillVar&, const SpillVar&);
  template SpillCut operator<(const SpillVar&, const SpillVar&);
  template SpillCut operator>=(const SpillVar&, const SpillVar&);
  template SpillCut operator<=(const SpillVar&, const SpillVar&);
  template SpillCut operator==(const SpillVar&, const SpillVar&);

  template SpillCut operator>(double, const SpillVar&);
  template SpillCut operator<(double, const SpillVar&);
  template SpillCut operator>=(double, const SpillVar&);
  template SpillCut operator<=(double, const SpillVar&);

  template SpillTruthCut operator>(const SpillTruthVar&, double);
  template SpillTruthCut operator<(const SpillTruthVar&, double);
  template SpillTruthCut operator>=(const SpillTruthVar&, double);
  template SpillTruthCut operator<=(const SpillTruthVar&, double);
  template SpillTruthCut operator==(const SpillTruthVar&, double);
  template SpillTruthCut operator!=(const SpillTruthVar&, double);

  template SpillTruthCut operator>(const SpillTruthVar&, const SpillTruthVar&);
  template SpillTruthCut operator<(const SpillTruthVar&, const SpillTruthVar&);
  template SpillTruthCut operator>=(const SpillTruthVar&, const SpillTruthVar&);
  template SpillTruthCut operator<=(const SpillTruthVar&, const SpillTruthVar&);
  template SpillTruthCut operator==(const SpillTruthVar&, const SpillTruthVar&);
  template SpillTruthCut operator!=(const SpillTruthVar&, const SpillTruthVar&);

  template SpillTruthCut operator>(double, const SpillTruthVar&);
  template SpillTruthCut operator<(double, const SpillTruthVar&);
  template SpillTruthCut operator>=(double, const SpillTruthVar&);
  template SpillTruthCut operator<=(double, const SpillTruthVar&);
}
