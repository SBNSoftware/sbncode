////////////////////////////////////////////////////////////////////////
// Class:       opHitFinderSBND
// Module Type: producer
// File:        opHitFinderSBND_module.cc
//
// This module produces an OpHit object for light analysis
// Created by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
//#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
//#include "lardataobj/Simulation/BeamGateInfo.h"

#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
//#include "larsim/MCCheater/PhotonBackTracker.h"

#include <memory>
#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TF1.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.h" 

namespace opdet{

  class opHitFinderSBND {
  public:
    opHitFinderSBND(fhicl::ParameterSet const & p, const detinfo::DetectorClocks *timeService);

    std::vector<recob::OpHit> MakeHits(const std::vector<raw::OpDetWaveform> &waveforms);
     
    opdet::sbndPDMapAlg map; //map for photon detector types

  private:

  // Declare member data here.
    std::string fInputModuleName;
  //  art::ServiceHandle<cheat::PhotonBackTracker> pbt;
    double fSampling; //in GHz
    double fBaselineSample; //in ticks
    double fUseDenoising; 
    double fPulsePolarityPMT; 
    double fPulsePolarityArapuca; 
    double fSaturation; //in number of p.e.
    double fArea1pePMT; //area of 1 pe in ADC*ns for PMTs
    double fArea1peSiPM; //area of 1 pe in ADC*ns for Arapucas
    int fThresholdPMT; //in ADC
    int fThresholdArapuca; //in ADC
    int fEvNumber;
    int fChNumber;
    //int fSize;
    //int fTimePMT;         //Start time of PMT signal
    //int fTimeMax;         //Time of maximum (minimum) PMT signal
    TH1D* wvfHist;          //processed waveform histogram 
    //TH1D* wvfHistPrint;          //processed waveform histogram 
    void subtractBaseline(TH1D* hist, std::string pdtype, double& rms);
    bool findPeak(TH1D* h, size_t& time, double& Area, double rms, double& amplitude, std::string type);
    void denoise(TH1D* h);
    void TV1D_denoise(float* input, float*& output, const int width, const float lambda);
    std::stringstream histname;
  };

}
