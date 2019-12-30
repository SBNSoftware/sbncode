#ifndef __sbnanalysis_ana_SBNOsc_Covariance__
#define __sbnanalysis_ana_SBNOsc_Covariance__

/**
 * \file Covariance.h
 */

#include "fhiclcpp/ParameterSet.h"
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
        void Initialize(fhicl::ParameterSet *config);
        void ProcessEvent(const event::Event *event, unsigned thread_index);
        void ProcessSubRun(const SubRun *subrun);
        void Finalize() { Scale(); GetCovs(); Write(); }

        // API Functions
        void GetCovs();
        void Write();
        void Scale();

        // build the covariance matrix
        TMatrixDSym CovarianceMatrix();
        
        // Output
        std::vector<TH2D> cov; //!< Covariance Matrix per variation
        std::vector<TH2D> fcov; //!< Fractional Covariance Matrix per variation
        std::vector<TH2D> corr; //!< Correlation Matrix per variation
    
    private:
        class EventSample {
          public:
    
	    /** Constructors. */
	    EventSample(const fhicl::ParameterSet &config, unsigned nUniverses, unsigned nVariations);
	    
	    double fScaleCC;
	    double fScaleNC;
	    double fScaleEnergy;
            double fPOT;
	    double fScalePOT;         //!< Factor for POT (etc.) scaling
	    std::vector <double> fBins; //!< Energy bin limits
	    TH1D *fCentralValue; //!< central value histogram
	    std::vector<std::vector<TH1D *>> fUniverses; //!< List of histogram per systematic universe. 
              // List has index per covariance matrix to be generated.
	    std::string fName; //!< Name for the sample
             std::vector<double> fEnergyBinScale;
        };
        void GetCovPerVariation(unsigned variation);

        // config
        std::vector<std::vector<std::string>> fWeightKeys;
        std::vector<std::vector<std::string>> fWeightKeysCC;
        std::vector<std::vector<std::string>> fWeightKeysNC;
        int fNumAltUnis;
        int fAltUniOffset;
        std::string fEnergyType;


        unsigned fNVariations;
        
        double fSelectionEfficiency;
        double fBackgroundRejection;
        std::string fOutputFile;

        bool fSaveCentralValue;
        bool fSaveUniverses;


        double fWeightMax;

        // file counter
        unsigned fSampleIndex;

        // Stored Event Samples
        std::vector<EventSample> fEventSamples;
};

}   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Covariance__
