#pragma once

#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"

#include <set>

namespace ana
{
  // Example reduction steps
  void ClearTrueParticles(caf::StandardRecord* sr);

  /// \brief Create smaller CAFs
  ///
  /// This class produces new CAFs, removing entries that fail a cut. It also
  /// allows the event record to be edited in custom ways.
  class FileReducer: protected SpectrumLoaderBase
  {
  public:
    FileReducer(const std::string& wildcard,
		const std::string& outfile);
    FileReducer(const std::vector<std::string>& fnames,
                const std::string& outfile);
    virtual ~FileReducer();

    /// Only copy records to the output file if they pass this cut
    void AddSpillCut(const SpillCut& cut);
    void AddSliceCut(const SliceCut& cut);

    /// \brief If called, only events whose run/subrun/event occur in \a fname
    /// will be retained.
    void SetEventList(const std::string& fname);

    typedef void (ReductionFunc)(caf::StandardRecord*);
    //    typedef void (ReductionFuncWithProxy)(caf::StandardRecord*,
    //                                          const caf::Proxy<caf::StandardRecord>*);

    /// Run the specified reduction function over each event
    void AddReductionStep(const std::function<ReductionFunc> &f) {fReductionFuncs.push_back(f);}
    //    void AddReductionStep(const std::function<ReductionFuncWithProxy> &f) {fReductionFuncsWithProxy.push_back(f);}

    /// Override any metadata key in the output file
    /*
    void SetMetadata(const std::string& key, const std::string& val)
    {
      fMetaMap[key] = val;
    }
    */

    //    void CopyMetadata(bool copy=true){fCopyMetadata = copy;};

    virtual void Go() override;

    // required by the interface, but not needed for anything done by FileReducer
    virtual void AccumulateExposures(const caf::SRSpill* spill) override {};

  protected:
    //    void UpdateMetadata(std::map<std::string, std::string>& meta,
    //                        const std::set<std::string>& mask,
    //                        const std::vector<std::string>& fnames) const;

    std::string fOutfile;
    SpillCut*   fSpillCut;
    SliceCut*   fSliceCut;

    std::set<std::tuple<int, int, int>> fEventList;

    std::vector<std::function<ReductionFunc>> fReductionFuncs;
    //    std::vector<std::function<ReductionFuncWithProxy>> fReductionFuncsWithProxy;

    std::map<std::string, std::string> fMetaMap;
    bool        fCopyMetadata;
  };
}
