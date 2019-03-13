#include "CAFAna/Core/Var.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace ana
{
  // explicitly instantiate the template for the types we know we have
  template class GenericVar<caf::StandardRecord>;
  template class GenericVar<caf::SRSpill>;
  template class GenericVar<caf::SRSpillTruthBranch>;

  template<class T> int GenericVar<T>::fgNextID = 0;

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>::
  GenericVar(const std::set<std::string>& /*reqs*/,
             const std::function<VarFunc_t>& fun)
    : fFunc(fun), fID(fgNextID++)
  {
  }

  //----------------------------------------------------------------------
  /// Helper for \ref Var2D
  template<class T> class Var2DFunc
  {
  public:
    Var2DFunc(const GenericVar<T>& a, const Binning binsa,
              const GenericVar<T>& b, const Binning binsb)
      : fA(a), fBinsA(binsa),
	fB(b), fBinsB(binsb)
    {
    }

    double operator()(const T* sr) const
    {
      // Calculate current values of the variables in StandardRecord once
      const double va = fA(sr);
      const double vb = fB(sr);

      // Since there are no overflow/underflow bins, check the range
      if(va < fBinsA.Min() || vb < fBinsB.Min()) return -1;
      if(va > fBinsA.Max() || vb > fBinsB.Max()) return fBinsA.NBins() * fBinsB.NBins();

      // FindBin uses root convention, first bin is bin 1, bin 0 is underflow
      const int ia = fBinsA.FindBin(va) - 1;
      const int ib = fBinsB.FindBin(vb) - 1;

      const int i = ia*fBinsB.NBins()+ib;

      return i+.5;
    }

  protected:
    const GenericVar<T> fA;
    const Binning fBinsA;
    const GenericVar<T> fB;
    const Binning fBinsB;
  };

  /// Helper for \ref Var3D
  template<class T> class Var3DFunc
  {
  public:
    Var3DFunc(const GenericVar<T>& a, const Binning binsa,
              const GenericVar<T>& b, const Binning binsb,
              const GenericVar<T>& c, const Binning binsc)
      : fA(a), fBinsA(binsa),
	fB(b), fBinsB(binsb),
	fC(c), fBinsC(binsc)
    {
    }

    double operator()(const T* sr) const
    {
      /// Calculate current values of the variables in StandardRecord once
      const double va = fA(sr);
      const double vb = fB(sr);
      const double vc = fC(sr);

      /// Since there are no overflow/underflow bins, check the range
      if(va < fBinsA.Min() || vb < fBinsB.Min() || vc < fBinsC.Min()){
        return -1.0;
      }

      if(va > fBinsA.Max() || vb > fBinsB.Max() || vc > fBinsC.Max()){
        return fBinsA.NBins() * fBinsB.NBins() * fBinsC.NBins();
      }

      const int ia = fBinsA.FindBin(va) - 1;
      const int ib = fBinsB.FindBin(vb) - 1;
      const int ic = fBinsC.FindBin(vc) - 1;
      const int i = ia*fBinsB.NBins()*fBinsC.NBins() + ib*fBinsC.NBins() + ic;

      return i+.5;
    }

  protected:
    const GenericVar<T> fA;
    const Binning fBinsA;
    const GenericVar<T> fB;
    const Binning fBinsB;
    const GenericVar<T> fC;
    const Binning fBinsC;
  };

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>
  Var2D(const GenericVar<T>& a, const Binning& binsa,
        const GenericVar<T>& b, const Binning& binsb)
  {
    return GenericVar<T>({},
                         Var2DFunc<T>(a, binsa, b, binsb));
  }

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>
  Var2D(const GenericVar<T>& a, int na, double a0, double a1,
        const GenericVar<T>& b, int nb, double b0, double b1)
  {
    return Var2D(a, Binning::Simple(na, a0, a1),
                 b, Binning::Simple(nb, b0, b1));
  }

  // explicitly instantiate the template for the types we know we have
  template Var Var2D(const Var&, const Binning&, const Var&, const Binning&);
  template SpillVar Var2D(const SpillVar&, const Binning&, const SpillVar&, const Binning&);
  template SpillTruthVar Var2D(const SpillTruthVar&, const Binning&, const SpillTruthVar&, const Binning&);

  template Var Var2D(const Var&, int, double, double, const Var&, int, double, double);
  template SpillVar Var2D(const SpillVar&, int, double, double, const SpillVar&, int, double, double);
  template SpillTruthVar Var2D(const SpillTruthVar&, int, double, double, const SpillTruthVar&, int, double, double);

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>
  Var3D(const GenericVar<T>& a, const Binning& binsa,
        const GenericVar<T>& b, const Binning& binsb,
        const GenericVar<T>& c, const Binning& binsc)
  {
    return GenericVar<T>({},
                         Var3DFunc<T>(a, binsa, b, binsb, c, binsc));
  }

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>
  Var3D(const GenericVar<T>& a, int na, double a0, double a1,
        const GenericVar<T>& b, int nb, double b0, double b1,
        const GenericVar<T>& c, int nc, double c0, double c1)
  {
    return Var3D(a, Binning::Simple(na, a0, a1),
                 b, Binning::Simple(nb, b0, b1),
                 c, Binning::Simple(nc, c0, c1));
  }

  // explicitly instantiate the template for the types we know we have
  template Var Var3D(const Var&, const Binning&, const Var&, const Binning&, const Var&, const Binning&);
  template SpillVar Var3D(const SpillVar&, const Binning&, const SpillVar&, const Binning&, const SpillVar&, const Binning&);
  template SpillTruthVar Var3D(const SpillTruthVar&, const Binning&, const SpillTruthVar&, const Binning&, const SpillTruthVar&, const Binning&);

  template Var Var3D(const Var&, int, double, double, const Var&, int, double, double, const Var&, int, double, double);
  template SpillVar Var3D(const SpillVar&, int, double, double, const SpillVar&, int, double, double, const SpillVar&, int, double, double);
  template SpillTruthVar Var3D(const SpillTruthVar&, int, double, double, const SpillTruthVar&, int, double, double, const SpillTruthVar&, int, double, double);

  //----------------------------------------------------------------------
  Var Scaled(const Var& v, double s)
  {
    return Var({},
               [v, s](const caf::StandardRecord* sr){return s*v(sr);});
  }

  //----------------------------------------------------------------------
  Var Constant(double c)
  {
    return Var({}, [c](const caf::StandardRecord*){return c;});
  }

  //--------------------------------------------------------------------

  Var Sqrt(const Var& v)
  {
    return Var({},
               [v](const caf::StandardRecord* sr){return sqrt(v(sr));});
  }

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>
  operator*(const GenericVar<T>& a, const GenericVar<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return GenericVar<T>({},
                           [a, b](const T* sr){return a(sr) * b(sr);},
                           ids[key]);
    }
    else{
      const GenericVar<T> ret({},
                              [a, b](const T* sr){return a(sr) * b(sr);});
      ids[key] = ret.ID();
      return ret;
    }
  }

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>
  operator/(const GenericVar<T>& a, const GenericVar<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return GenericVar<T>({},
                           [a, b](const T* sr)
                           {
                             const double denom = b(sr);
                             if(denom != 0)
                               return a(sr) / denom;
                             else
                               return 0.0;
                           },
                           ids[key]);
    }
    else{
      const GenericVar<T> ret({},
                              [a, b](const T* sr)
                              {
                                const double denom = b(sr);
                                if(denom != 0)
                                  return a(sr) / denom;
                                else
                                  return 0.0;
                              });
      ids[key] = ret.ID();
      return ret;
    }
  }

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>
  operator+(const GenericVar<T>& a, const GenericVar<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return GenericVar<T>({},
                           [a, b](const T* sr){return a(sr) + b(sr);},
                           ids[key]);
    }
    else{
      const GenericVar<T> ret({},
                              [a, b](const T* sr){return a(sr) + b(sr);});
      ids[key] = ret.ID();
      return ret;
    }
  }

  //----------------------------------------------------------------------
  template<class T> GenericVar<T>
  operator-(const GenericVar<T>& a, const GenericVar<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return GenericVar<T>({},
                           [a, b](const T* sr){return a(sr) - b(sr);},
                           ids[key]);
    }
    else{
      const GenericVar<T> ret({},
                              [a, b](const T* sr){return a(sr) - b(sr);});
      ids[key] = ret.ID();
      return ret;
    }
  }


  // explicitly instantiate the templates for the types we know we have
  template Var operator*(const Var&, const Var&);
  template Var operator/(const Var&, const Var&);
  template Var operator+(const Var&, const Var&);
  template Var operator-(const Var&, const Var&);
  template SpillVar operator*(const SpillVar&, const SpillVar&);
  template SpillVar operator/(const SpillVar&, const SpillVar&);
  template SpillVar operator+(const SpillVar&, const SpillVar&);
  template SpillVar operator-(const SpillVar&, const SpillVar&);
  template SpillTruthVar operator*(const SpillTruthVar&, const SpillTruthVar&);
  template SpillTruthVar operator/(const SpillTruthVar&, const SpillTruthVar&);
  template SpillTruthVar operator+(const SpillTruthVar&, const SpillTruthVar&);
  template SpillTruthVar operator-(const SpillTruthVar&, const SpillTruthVar&);
}
