#pragma once

#include <functional>
#include <set>
#include <string>
#include <vector>

namespace caf{class StandardRecord;}

namespace ana
{
  /// A Var that returns multiple results for each slice. eg the properties of
  /// multiple prongs. All results will be filled into the Spectrum.
  class MultiVar
  {
  public:
    /// The type of the function part of a var
    typedef std::vector<double> (VarFunc_t)(const caf::StandardRecord* sr);

    /// std::function can wrap a real function, function object, or lambda
    MultiVar(const std::set<std::string>& reqs,
             const std::function<VarFunc_t>& fun)
      : fFunc(fun), fID(fgNextID--)
    {
    }

    /// Allows a variable to be called with double value = myVar(sr) syntax
    std::vector<double> operator()(const caf::StandardRecord* sr) const
    {
      return fFunc(sr);
    }

    /// Vars with the same definition will have the same ID
    int ID() const {return fID;}

    static int MaxID() {return fgNextID-1;}
  protected:
    std::function<VarFunc_t> fFunc;

    int fID;
    /// The next ID that hasn't yet been assigned
    static int fgNextID;
  };

} // namespace
