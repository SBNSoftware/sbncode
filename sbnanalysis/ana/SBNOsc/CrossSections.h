#ifndef __sbnanalysis_ana_SBNOsc_CrossSections__
#define __sbnanalysis_ana_SBNOsc_CrossSections__

/**
 * \file CrossSections.h
 */

#include "fhiclcpp/ParameterSet.h"
#include "core/PostProcessorBase.hh"

#include <string>
#include <vector>
#include <map>
#include <string>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom2.h>

class TTree;

namespace ana {
namespace SBNOsc {

class CrossSections: public core::PostProcessorBase {
    public:
        // Constructor
        CrossSections() {}

        // implementing PostProcessor
        void FileCleanup(TTree *eventTree);
        void Initialize(fhicl::ParameterSet *config);
        void ProcessEvent(const Event *event);
        void Finalize(){Write();}
        
        // API Functions
        void Write();
        double SmearMcsMomentum(double momentum);
        double SmearRangeMomentum(double momentum);

        // config
        std::string fOutputFile;

        double fMinX;
        double fMaxX;
        double fMinY;
        double fMaxY;
        double fMinZ;
        double fMaxZ;

        double fMinContainedLength;
        double fMinExitingLength;

        double fProtonThreshold;
        double fPionThreshold;
        double fMuonThreshold;
        double fPi0Threshold;

        double fProtonPidEff;
        double fPionPidEff;

        double fProtonRecoEff;
        double fMuonRecoEff;
        double fPionRecoEff;
        double fPi0RecoEff;

        TRandom2* fRandom;

        // Hists
        std::map<std::string, TH1D*> hTrueNuE;
        std::map<std::string, TH1D*> hTrueMuP;
        std::map<std::string, TH1D*> hTrueMuTheta;
        std::map<std::string, TH1D*> hRecoNuE;
        std::map<std::string, TH1D*> hRecoMuP;
        std::map<std::string, TH1D*> hRecoMuTheta;
        TH1D* hNuETrue;
        TH1D* hNuEReco;
        
};

}   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_CrossSections__
