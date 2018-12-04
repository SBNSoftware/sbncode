#ifndef __sbnanalysis_ana_SBNOsc_Covariance__
#define __sbnanalysis_ana_SBNOsc_Covariance__

/**
 * \file Covariance.h
 */

#include "json/json.h"
#include "core/PostProcessorBase.hh"

#include <string>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cassert>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMatrixDSym.h>

class TTree;

namespace ana {
namespace SBNOsc {

class Covariance: public core::PostProcessorBase {
    public:
        // Constructor
        Covariance() {}

        // implementing PostProcessor
        void FileCleanup(TTree *eventTree);
        void Initialize(Json::Value *config);
        void ProcessEvent(const Event *event);
        void Finalize() { GetCovs(); Write(); }

        // API Functions
        void GetCovs();
        void Write();

        // build the covariance matrix
        TMatrixDSym CovarianceMatrix();
        
        // Output
        TH2D *cov; //!< Covariance Matrix
        TH2D *fcov; //!< Fractional Covariance Matrix
        TH2D *corr; //!< Correlation Matrix

    
    private:
        class EventSample {
          public:
    
	    /** Constructors. */
	    EventSample(const Json::Value &config, unsigned nUniverses);
	    
	    double fScaleFactor;         //!< Factor for POT (etc.) scaling
	    std::vector <double> fBins; //!< Energy bin limits
	    TH1D *fCentralValue; //!< central value histogram
	    std::vector<TH1D *> fUniverses; //!< histogram per systematic universe
	    std::string fName; //!< Name for the sample
        };

        // config
        std::vector<std::string> fWeightKeys;
        std::vector<std::string> fUniformWeights;
        int fNumAltUnis;
        std::string fEnergyType;
        
        double fSelectionEfficiency;
        double fBackgroundRejection;
        std::string fOutputFile;

        bool fSaveCentralValue;
        bool fSaveUniverses;

        // file counter
        unsigned fSampleIndex;

        // Stored Event Samples
        std::vector<EventSample> fEventSamples;
        
};

}   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Covariance__
