#include "CAFAna/Core/Var.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace ana
{
  //----------------------------------------------------------------------
  template<class T> _Var<T>::
  _Var(const std::function<VarFunc_t>& fun)
    : fFunc(fun), fID(fgNextID++)
  {
  }

  //----------------------------------------------------------------------
  template<class T> int _Var<T>::MaxID()
  {
    return fgNextID-1;
  }

  //----------------------------------------------------------------------
  /// Helper for \ref Var2D
  template<class T> class Var2DFunc
  {
  public:
    Var2DFunc(const _Var<T>& a, const Binning binsa,
              const _Var<T>& b, const Binning binsb)
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
    const _Var<T> fA;
    const Binning fBinsA;
    const _Var<T> fB;
    const Binning fBinsB;
  };

  /// Helper for \ref Var3D
  template<class T> class Var3DFunc
  {
  public:
    Var3DFunc(const _Var<T>& a, const Binning binsa,
              const _Var<T>& b, const Binning binsb,
              const _Var<T>& c, const Binning binsc)
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
    const _Var<T> fA;
    const Binning fBinsA;
    const _Var<T> fB;
    const Binning fBinsB;
    const _Var<T> fC;
    const Binning fBinsC;
  };

  //----------------------------------------------------------------------
  template<class T> _Var<T>
  Var2D(const _Var<T>& a, const Binning& binsa,
        const _Var<T>& b, const Binning& binsb)
  {
    return _Var<T>(Var2DFunc<T>(a, binsa, b, binsb));
  }

  //----------------------------------------------------------------------
  template<class T> _Var<T>
  Var2D(const _Var<T>& a, int na, double a0, double a1,
        const _Var<T>& b, int nb, double b0, double b1)
  {
    return Var2D(a, Binning::Simple(na, a0, a1),
                 b, Binning::Simple(nb, b0, b1));
  }

  // explicitly instantiate the template for the types we know we have
  template Var Var2D(const Var&, const Binning&, const Var&, const Binning&);
  template SpillVar Var2D(const SpillVar&, const Binning&, const SpillVar&, const Binning&);

  template Var Var2D(const Var&, int, double, double, const Var&, int, double, double);
  template SpillVar Var2D(const SpillVar&, int, double, double, const SpillVar&, int, double, double);

  //----------------------------------------------------------------------
  template<class T> _Var<T>
  Var3D(const _Var<T>& a, const Binning& binsa,
        const _Var<T>& b, const Binning& binsb,
        const _Var<T>& c, const Binning& binsc)
  {
    return _Var<T>(Var3DFunc<T>(a, binsa, b, binsb, c, binsc));
  }

  //----------------------------------------------------------------------
  template<class T> _Var<T>
  Var3D(const _Var<T>& a, int na, double a0, double a1,
        const _Var<T>& b, int nb, double b0, double b1,
        const _Var<T>& c, int nc, double c0, double c1)
  {
    return Var3D(a, Binning::Simple(na, a0, a1),
                 b, Binning::Simple(nb, b0, b1),
                 c, Binning::Simple(nc, c0, c1));
  }

  // explicitly instantiate the template for the types we know we have
  template Var Var3D(const Var&, const Binning&, const Var&, const Binning&, const Var&, const Binning&);
  template SpillVar Var3D(const SpillVar&, const Binning&, const SpillVar&, const Binning&, const SpillVar&, const Binning&);

  template Var Var3D(const Var&, int, double, double, const Var&, int, double, double, const Var&, int, double, double);
  template SpillVar Var3D(const SpillVar&, int, double, double, const SpillVar&, int, double, double, const SpillVar&, int, double, double);

  //----------------------------------------------------------------------
  Var Scaled(const Var& v, double s)
  {
    return Var([v, s](const caf::SRSliceProxy* sr){return s*v(sr);});
  }

  //----------------------------------------------------------------------
  Var Constant(double c)
  {
    return Var([c](const caf::SRSliceProxy*){return c;});
  }

  //--------------------------------------------------------------------

  Var Sqrt(const Var& v)
  {
    return Var([v](const caf::SRSliceProxy* sr){return sqrt(v(sr));});
  }

  //----------------------------------------------------------------------
  template<class T> _Var<T>
  operator*(const _Var<T>& a, const _Var<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return _Var<T>([a, b](const T* sr){return a(sr) * b(sr);},
                           ids[key]);
    }
    else{
      const _Var<T> ret([a, b](const T* sr){return a(sr) * b(sr);});
      ids[key] = ret.ID();
      return ret;
    }
  }

  //----------------------------------------------------------------------
  template<class T> _Var<T>
  operator/(const _Var<T>& a, const _Var<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return _Var<T>([a, b](const T* sr)
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
      const _Var<T> ret([a, b](const T* sr)
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
  template<class T> _Var<T>
  operator+(const _Var<T>& a, const _Var<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return _Var<T>([a, b](const T* sr){return a(sr) + b(sr);},
                           ids[key]);
    }
    else{
      const _Var<T> ret([a, b](const T* sr){return a(sr) + b(sr);});
      ids[key] = ret.ID();
      return ret;
    }
  }

  //----------------------------------------------------------------------
  template<class T> _Var<T>
  operator-(const _Var<T>& a, const _Var<T>& b)
  {
    static std::map<std::pair<int, int>, int> ids;
    const std::pair<int, int> key(a.ID(), b.ID());

    if(ids.count(key)){
      return _Var<T>([a, b](const T* sr){return a(sr) - b(sr);},
                           ids[key]);
    }
    else{
      const _Var<T> ret([a, b](const T* sr){return a(sr) - b(sr);});
      ids[key] = ret.ID();
      return ret;
    }
  }

  // explicitly instantiate the templates for the types we know we have
  template class _Var<caf::SRSpillProxy>;
  template class _Var<caf::SRSliceProxy>;

  template<class T> int _Var<T>::fgNextID = 0;

  template Var operator*(const Var&, const Var&);
  template Var operator/(const Var&, const Var&);
  template Var operator+(const Var&, const Var&);
  template Var operator-(const Var&, const Var&);
  template SpillVar operator*(const SpillVar&, const SpillVar&);
  template SpillVar operator/(const SpillVar&, const SpillVar&);
  template SpillVar operator+(const SpillVar&, const SpillVar&);
  template SpillVar operator-(const SpillVar&, const SpillVar&);
}
