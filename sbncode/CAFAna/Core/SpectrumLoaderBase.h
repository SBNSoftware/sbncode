#pragma once

#include <functional>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/IFileSource.h"
#include "CAFAna/Core/MultiVar.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Core/Var.h"

namespace caf{class StandardRecord;}

class TFile;
class TH1;
class TTree;

namespace ana
{
  class Spectrum;
  class ReweightableSpectrum;

  /// Is this data-file representing beam spills or cosmic spills?
  enum DataSource{
    // TODO there are no longer any actions taken as a result of this
    // distinction. Remove after SA.
    kBeam,
    kCosmic
  };

  /// Base class for the various types of spectrum loader
  class SpectrumLoaderBase
  {
  public:

    friend class ReweightableSpectrum;
    friend class NDOscillatableSpectrum;
    friend class OscillatableSpectrum;
    friend class Spectrum;

    virtual ~SpectrumLoaderBase();

    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const Var& var,
                             const Cut& cut,
                             const SystShifts& shift,
                             const Var& wei = kUnweighted);

    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const MultiVar& var,
                             const Cut& cut,
                             const SystShifts& shift,
                             const Var& wei = kUnweighted);

    /// For use by the constructors of \ref ReweightableSpectrum subclasses
    virtual void AddReweightableSpectrum(ReweightableSpectrum& spect,
                                         const Var& var,
                                         const Cut& cut,
                                         const SystShifts& shift,
                                         const Var& wei);

    /// Load all the registered spectra
    virtual void Go() = 0;

    /// Indicate whether or not \ref Go has been called
    virtual bool Gone() const {return fGone;}

  protected:
    /// Component of other constructors
    SpectrumLoaderBase(DataSource src = kBeam);
    /// Construct from a filename, wildcard, SAM definition, or SAM query
    SpectrumLoaderBase(const std::string& wildcard, DataSource src = kBeam);
    /// Construct from an explicit list of files
    SpectrumLoaderBase(const std::vector<std::string>& fnames, DataSource src = kBeam);

    // Move operations
    SpectrumLoaderBase(SpectrumLoaderBase&&) = default;
    SpectrumLoaderBase& operator=(SpectrumLoaderBase&&) = default;

    // No copy operations because I don't want to deal with pointers
    SpectrumLoaderBase(const SpectrumLoaderBase&) = delete;
    SpectrumLoaderBase& operator=(const SpectrumLoaderBase&) = delete;

    /// Figure out if \a str is a wildcard or SAM query and return a source
    IFileSource* WildcardOrSAMQuery(const std::string& str) const;

    friend class SpectrumLoaderMockData;
    virtual void RemoveSpectrum(Spectrum*);
    virtual void RemoveReweightableSpectrum(ReweightableSpectrum*);

    virtual void AccumulateExposures(const caf::SRSpill* spill) = 0;

    /// Forwards to \ref fFileSource
    int NFiles() const;

    /// Forwards to \ref fFileSource but also accumulates POT and livetime
    TFile* GetNextFile();

    std::string fWildcard;
    std::unique_ptr<IFileSource> fFileSource;

    DataSource fSource;

    bool fGone; ///< Has Go() been called? Can't add more histograms after that

    double fPOT; ///< Accumulated by calls to \ref GetNextFile

    /// \brief Helper class for \ref SpectrumLoaderBase
    ///
    /// List of Spectrum and OscillatableSpectrum, some utility functions
    struct SpectList
    {
      void Erase(Spectrum* s);
      void Erase(ReweightableSpectrum* os);
      void RemoveLoader(SpectrumLoaderBase* l);
      size_t TotalSize() const;
      void GetSpectra(std::vector<Spectrum*>& ss);
      void GetReweightableSpectra(std::vector<ReweightableSpectrum*>& ss);

      std::vector<Spectrum*> spects;
      std::vector<ReweightableSpectrum*> rwSpects;
    };

    /// \brief Helper class for \ref SpectrumLoaderBase
    ///
    /// Functions like std::map<T, U> except it should be faster to iterate
    /// through the elements (while slower to fill) and it knows to compare Ts
    /// via their ID() function. Various methods that forward through to the
    /// \ref SpectList at the end of the chain.
    template<class T, class U> struct IDMap
    {
      U& operator[](const T& key);

      // Make class iterable. Keep inline for speed
      typedef typename std::vector<std::pair<T, U>>::iterator it_t;
      inline it_t begin(){return fElems.begin();}
      inline it_t end(){return fElems.end();}

      template<class V> void Erase(const V& v);
      void RemoveLoader(SpectrumLoaderBase* l);
      void Clear();
      size_t TotalSize();
      void GetSpectra(std::vector<Spectrum*>& ss);
      void GetReweightableSpectra(std::vector<ReweightableSpectrum*>& ss);
    protected:
      std::vector<std::pair<T, U>> fElems;
    };

    class VarOrMultiVar
    {
    public:
      // v could easily be a temporary, have to make a copy
      VarOrMultiVar(const Var& v) : fVar(new Var(v)), fMultiVar(0) {}
      VarOrMultiVar(const MultiVar& v) : fVar(0), fMultiVar(new MultiVar(v)) {}
      ~VarOrMultiVar() {delete fVar; delete fMultiVar;}

      VarOrMultiVar(const VarOrMultiVar& v)
        : fVar(v.fVar ? new Var(*v.fVar) : 0),
          fMultiVar(v.fMultiVar ? new MultiVar(*v.fMultiVar) : 0)
      {
      }

      VarOrMultiVar(VarOrMultiVar&& v)
      {
        fVar = v.fVar;
        fMultiVar = v.fMultiVar;
        v.fVar = 0;
        v.fMultiVar = 0;
      }

      bool IsMulti() const {return fMultiVar;}
      const Var& GetVar() const {assert(fVar); return *fVar;}
      const MultiVar& GetMultiVar() const {assert(fMultiVar); return *fMultiVar;}

      int ID() const {return fVar ? fVar->ID() : fMultiVar->ID();}

    protected:
      const Var* fVar;
      const MultiVar* fMultiVar;
    };

    /// \brief All the spectra that need to be filled
    ///
    /// [shift][cut][wei][var]
    IDMap<SystShifts, IDMap<Cut, IDMap<Var, IDMap<VarOrMultiVar, SpectList>>>> fHistDefs;
  };

  /// \brief Dummy loader that doesn't load any files
  ///
  /// Useful when a loader is required for a component you want to ignore
  class NullLoader: public SpectrumLoaderBase
  {
  public:
    NullLoader() {}
    ~NullLoader();

    virtual void Go() override;
    void AddSpectrum(Spectrum& spect,
                     const Var& var,
                     const Cut& cut,
                     const SystShifts& shift,
                     const Var& wei = kUnweighted) override {}
    void AddSpectrum(Spectrum& spect,
                     const MultiVar& var,
                     const Cut& cut,
                     const SystShifts& shift,
                     const Var& wei = kUnweighted) override {}

    void AddReweightableSpectrum(ReweightableSpectrum& spect,
                                 const Var& var,
                                 const Cut& cut,
                                 const SystShifts& shift,
                                 const Var& wei) override {}

    void AccumulateExposures(const caf::SRSpill* spill) override {};
  };
  /// \brief Dummy loader that doesn't load any files
  ///
  /// Useful when a loader is required for a component you want to ignore
  static NullLoader kNullLoader;
}

