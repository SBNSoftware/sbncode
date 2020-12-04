#include "CAFAna/Core/EventList.h"

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"

namespace ana
{
  /// Helper class for \ref MakeTextListFile
  class ASCIIMaker: public SpectrumLoader
  {
  public:
    ASCIIMaker(const std::string& wildcard,
               const std::vector<SpillCut>& cut,
               const std::vector<FILE*>& f,
               const std::vector<const SpillVar*>& floatVars,
               const std::vector<const SpillVar*>& intVars)
      : SpectrumLoader(wildcard),
        fCut(cut),
        fFile(f),
        fFloatVars(floatVars),
        fIntVars(intVars),
        fNPassed(0)
    {
    }

    ASCIIMaker(const std::string& wildcard,
               const std::vector<SpillCut>& cut,
               const std::vector<FILE*>& f,
               const std::vector<const SpillMultiVar*>& multivars,
               const std::vector<const SpillMultiVar*>& /*multiIntVars*/)
      : SpectrumLoader(wildcard),
        fCut(cut),
        fFile(f),
        fMultiVars(multivars),
        fNPassed(0)
    {
    }

    ASCIIMaker(const std::vector<std::string>& fnames,
               const std::vector<SpillCut>& cut,
               const std::vector<FILE*>& f,
               const std::vector<const SpillVar*>& floatVars,
               const std::vector<const SpillVar*>& intVars)
      : SpectrumLoader(fnames),
        fCut(cut),
        fFile(f),
        fFloatVars(floatVars),
        fIntVars(intVars),
        fNPassed(0)
    {
    }

    ~ASCIIMaker()
    {
      std::cout << "Selected " << fNPassed << " slices." << std::endl;
    }

    void HandleRecord(caf::Proxy<caf::StandardRecord>* sr) override
    {
      for(unsigned int ic = 0; ic < fCut.size(); ic++){
        if(!fCut[ic](sr))
          continue;

        ++fNPassed;

        for(const SpillVar* v: fIntVars){
          fprintf(fFile[ic], "%d ", int(std::round((*v)(sr))));
        }

        for(const SpillVar* v: fFloatVars){
          fprintf(fFile[ic], "%g ", (*v)(sr));
        }

	for(const SpillMultiVar* mv: fMultiVars){
          for(const double md: (*mv)(sr)){
	    fprintf(fFile[ic], "%g ", md);
          }
	}

        fprintf(fFile[ic], "\n");
      }// end loop over cuts
    }
  protected:
    std::vector<SpillCut> fCut;
    std::vector<FILE*> fFile;
    std::vector<const SpillVar*> fFloatVars;
    std::vector<const SpillVar*> fIntVars;
    std::vector<const SpillMultiVar*> fMultiVars;
    int fNPassed;
  };

  /// Helper class for \ref MakeEventTTreeFile
  class TreeMaker: public SpectrumLoader
  {
  public:
    TreeMaker(const std::string& wildcard,
              const std::vector<std::pair<std::string, SpillCut>>& cuts,
              const std::vector<std::pair<std::string, SpillVar>>& floatVars,
              const std::vector<std::pair<std::string, SpillVar>>& intVars)
      : SpectrumLoader(wildcard),
        fFloats(floatVars.size()), fInts(intVars.size())
    {
      for(unsigned int ic = 0; ic < cuts.size(); ++ic){
        fCuts.push_back(cuts[ic].second);
        fTrees.push_back(new TTree(cuts[ic].first.c_str(),
                                   cuts[ic].first.c_str()));
        for(unsigned int iv = 0; iv < floatVars.size(); ++iv){
          fFloatVars.push_back(floatVars[iv].second);
          fTrees.back()->Branch(floatVars[iv].first.c_str(), &fFloats[iv]);
        }
        for(unsigned int iv = 0; iv < intVars.size(); ++iv){
          fIntVars.push_back(intVars[iv].second);
          fTrees.back()->Branch(intVars[iv].first.c_str(), &fInts[iv]);
        }
      }
    }

    void HandleRecord(caf::Proxy<caf::StandardRecord>* sr) override
    {
      bool any = false;

      for(unsigned int ic = 0; ic < fCuts.size(); ++ic){
        if(!fCuts[ic](sr)) continue;

        if(!any){
          any = true;
          for(unsigned int iv = 0; iv < fFloatVars.size(); ++iv){
            fFloats[iv] = fFloatVars[iv](sr);
          }
          for(unsigned int iv = 0; iv < fIntVars.size(); ++iv){
            fInts[iv] = fIntVars[iv](sr);
          }
        }

        fTrees[ic]->Fill();
      }
    }

    void Write()
    {
      for(TTree* tr: fTrees) tr->Write();
    }

  protected:
    std::vector<SpillCut> fCuts;
    std::vector<SpillVar> fFloatVars;
    std::vector<SpillVar> fIntVars;
    std::vector<TTree*> fTrees;
    std::vector<float> fFloats;
    std::vector<int> fInts;
  };

  //----------------------------------------------------------------------
  //
  // All the variants below have the same bodies (here), but slightly different
  // arguments. We have T=wildcard/file list and U=Var/MultiVar
  template<class T, class U>
  void MakeTextListFileHelper(const T& source,
                              const std::vector<SpillCut>& cut,
                              const std::vector<std::string>& output,
                              const std::vector<U>& floatVars,
                              const std::vector<U>& intVars)
  {
    assert(output.size() == cut.size());

    std::vector<FILE*> files;
    for(const std::string& out: output)
      files.push_back(fopen(out.c_str(), "w"));

    ASCIIMaker maker(source, cut, files, floatVars, intVars);
    maker.Go();

    for(const std::string& out: output)
      std::cout << "Wrote text list file " << out << std::endl;

    for(FILE* f: files) fclose(f);
  }

  //----------------------------------------------------------------------
  void MakeTextListFile(const std::string& wildcard,
			const std::vector<SpillCut>& cut,
			const std::vector<std::string>& output,
			const std::vector<const SpillVar*>& floatVars,
			const std::vector<const SpillVar*>& intVars)
  {
    MakeTextListFileHelper(wildcard, cut, output, floatVars, intVars);
  }

  //----------------------------------------------------------------------

  void MakeTextListFile(const std::vector<std::string>& fnames,
			const std::vector<SpillCut>& cut,
			const std::vector<std::string>& output,
			const std::vector<const SpillVar*>& floatVars,
			const std::vector<const SpillVar*>& intVars)
  {
    MakeTextListFileHelper(fnames, cut, output, floatVars, intVars);
  }

  //----------------------------------------------------------------------
  void MakeTextListFile(const std::string& wildcard,
                        const std::vector<SpillCut>& cut,
                        const std::vector<std::string>& output,
                        const std::vector<const SpillMultiVar*>& multivars)
  {
    MakeTextListFileHelper(wildcard, cut, output, multivars, {});
  }

  //----------------------------------------------------------------------
  //
  // T can be a wildcard or list of file names
  template<class T>
  void MakeEventListFileHelper(const T& source,
                               const std::vector<SpillCut>& cuts,
                               const std::vector<std::string>& outputs,
                               bool includeSliceIndex,
                               bool includeSliceTime,
                               bool includeCycleNumber,
                               bool includeBatchNumber)
  {
    assert(cuts.size() == outputs.size());

    std::vector<const SpillVar*> intVars = {new SIMPLESPILLVAR(hdr.run),
                                            new SIMPLESPILLVAR(hdr.subrun)};

    // TODO TODO TODO
    //    if(includeCycleNumber) intVars.push_back(new SIMPLESPILLVAR(hdr.cycle));
    //    if(includeBatchNumber) intVars.push_back(new SIMPLESPILLVAR(hdr.batch));
    intVars.push_back(new SIMPLESPILLVAR(hdr.evt));
    if(includeSliceIndex) intVars.push_back(new SIMPLESPILLVAR(hdr.subevt));

    std::vector<const SpillVar*> floatVars;
    // TODO TODO TODO
    //    if(includeSliceTime) floatVars.push_back(new SIMPLESPILLVAR(hdr.subevtmeantime));

    MakeTextListFile(source, cuts, outputs, floatVars, intVars);
  }

  //----------------------------------------------------------------------
  void MakeEventListFile(const std::string& wildcard,
                         const std::vector<SpillCut>& cuts,
                         const std::vector<std::string>& outputs,
                         bool includeSliceIndex,
                         bool includeSliceTime,
                         bool includeCycleNumber,
			 bool includeBatchNumber)
  {
    MakeEventListFileHelper(wildcard, cuts, outputs,
                            includeSliceIndex, includeSliceTime, includeCycleNumber, includeBatchNumber);
  }

  //----------------------------------------------------------------------
  void MakeEventListFile(const std::vector<std::string>& fnames,
                         const std::vector<SpillCut>& cuts,
                         const std::vector<std::string>& outputs,
                         bool includeSliceIndex,
                         bool includeSliceTime,
                         bool includeCycleNumber,
			 bool includeBatchNumber)
  {
    MakeEventListFileHelper(fnames, cuts, outputs,
                            includeSliceIndex,
                            includeSliceTime,
                            includeCycleNumber,
                            includeBatchNumber);
  }

  //----------------------------------------------------------------------
  void MakeEventListFile(const std::string& wildcard,
                         const SpillCut& cut,
                         const std::string& output,
                         bool includeSliceIndex,
                         bool includeSliceTime,
                         bool includeCycleNumber,
			 bool includeBatchNumber)
  {
    MakeEventListFileHelper(wildcard, {cut}, {output},
                            includeSliceIndex,
                            includeSliceTime,
                            includeCycleNumber,
                            includeBatchNumber);
  }

  //----------------------------------------------------------------------
  void MakeEventListFile(const std::vector<std::string>& fnames,
                         const SpillCut& cut,
                         const std::string& output,
                         bool includeSliceIndex,
                         bool includeSliceTime,
                         bool includeCycleNumber,
			 bool includeBatchNumber)
  {
    MakeEventListFileHelper(fnames, {cut}, {output},
                            includeSliceIndex,
                            includeSliceTime,
                            includeCycleNumber,
                            includeBatchNumber);

  }

  //----------------------------------------------------------------------
  void MakeEventTTreeFile(const std::string& wildcard,
                          const std::string& output,
                          const std::vector<std::pair<std::string, SpillCut>>& cuts,
                          const std::vector<std::pair<std::string, SpillVar>>& floatVars,
                          const std::vector<std::pair<std::string, SpillVar>>& intVars)
  {
    // Scope all the TTrees within this
    TFile fout(output.c_str(), "RECREATE");

    TreeMaker maker(wildcard, cuts, floatVars, intVars);

    maker.Go();

    fout.cd();
    maker.Write();
    fout.Close();
    std::cout << "Wrote TTree file " << output << std::endl;
  }

}
