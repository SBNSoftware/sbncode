#ifndef __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__
#define __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__

/**
 * \file Chi2Sensitivity.h
 */

#include "core/PostProcessorBase.hh"

#include "Covariance.h"

#include <vector>
#include <string>

#include <TGraph2D.h>
#include <TGraph.h>

namespace ana {
namespace SBNOsc {

class Chi2Sensitivity: public core::PostProcessorBase {
    
    public:
        
        // Constructor
        Chi2Sensitivity() {}
        // Implement post-processor
        void Initialize(fhicl::ParameterSet* config);
        void ProcessEvent(const Event* event);
        void FileCleanup(TTree* eventTree);

        void Finalize() { fCovariance.Finalize(); Scale(); GetChi2(); GetContours(); Write(); }

        // API Functions
        void Scale();
        void GetChi2();
        void GetContours();
        void Write();

        TH1D *Oscillate(double sinth, double dm2);
        
        // Output
        
        TGraph2D *chisqplot;
        TGraph *contour_90pct, *contour_3sigma, *contour_5sigma;
    
    private:
        // struct for additional config on event sample detectors
        class EventSample {
          public:
            // Constructor
            EventSample(const fhicl::ParameterSet& config);
            // Oscillate Signal counts into 1D Histogram
            std::vector<double> Oscillate(double sinth, double dm2) const;
            // Get Vector of Signal counts
            std::vector<double> Signal() const;
            // Get vector of Background counts
            std::vector<double> Background() const;
 
            // Config
            std::string fName;
            double fDistance; //!< Distance in km
            std::array<double, 2> fXlim; //!< Detector size in cm
            std::array<double, 2> fYlim; //!< Detector size in cm
            std::array<double, 2> fZlim; //!< Detector size in cm
            double fScaleFactor; //!< Factor for POT (etc.) scaling   
            int fOscType; //!< Oscilaltion type: 0 == None, 1 == numu -> nue, 2 == numu -> numu
             std::vector<double> fEnergyBinScale;

            // Storage
            TH1D *fBkgCounts; //!< Background count Histogram
            TH3D *fSignalCounts; //!< Signal Count Histogram 
  
            // bins
            std::vector<double> fBins; //!< Reco Energy bin limits
            std::vector<double> fTrueEBins; //!< True energy bin limits
            std::vector<double> fDistBins; //!< Distance bin limits 
        };
        
        // From config file
        std::string fEnergyType;
        
        double fSelectionEfficiency, fBackgroundRejection;
        
        int fNumDm2;
        int fNumSin;
        std::vector <double> fLogDm2Lims, fLogSinLims;
        std::vector<std::string> fUniformWeights;
        
        std::string fOutputFile;
        // whether to save stuff
        int fSavePDFs;
        bool fSaveSignal;
        bool fSaveBackground;
        std::vector<std::array<double, 2>> fSaveOscillations;

        // index into sample
        unsigned fSampleIndex;
        
        // keep internal covariance
        Covariance fCovariance;
        // keep own Event Samples
        std::vector<EventSample> fEventSamples;
        
        // containers for sin2theta and dm values
        std::vector<double> sin2theta;
        std::vector<double> dm2;

        std::vector <std::vector <double> > chisq_diffs; //!< Container for chi2 values
    
};

}  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__
