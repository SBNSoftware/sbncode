//C++ Includes 
#include <iostream>
#include <vector>

//Root Includes 
#include <TH2D.h>

//Framework Includes 
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"

//Larsoft Includes 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"

//SBN Includes 
#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

    NueSelection::NueSelection() : SelectionBase(), EventCounter(0), NuCount(0) {}

    
    void NueSelection::Initialize(fhicl::ParameterSet* config) {

      rand.SetSeed(0);

      hello();

      fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("NueSelection");

      // setup active volume bounding boxes
      std::vector<fhicl::ParameterSet> AVs =				\
        pconfig.get<std::vector<fhicl::ParameterSet> >("active_volumes");
      for (auto const& AV : AVs) {
        double xmin = AV.get<double>("xmin");
        double ymin = AV.get<double>("ymin");
        double zmin = AV.get<double>("zmin");
        double xmax = AV.get<double>("xmax");
        double ymax = AV.get<double>("ymax");
        double zmax = AV.get<double>("zmax");
        fConfig.active_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
      }

      std::vector<fhicl::ParameterSet> FVs =				\
        pconfig.get<std::vector<fhicl::ParameterSet> >("fiducial_volumes");
      for (auto const& FV : FVs) {
        double xmin = FV.get<double>("xmin");
        double ymin = FV.get<double>("ymin");
        double zmin = FV.get<double>("zmin");
        double xmax = FV.get<double>("xmax");
        double ymax = FV.get<double>("ymax");
        double zmax = FV.get<double>("zmax");
        fConfig.fiducial_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
      }


      //Set up the fcl parameters.
      fConfig.minLengthExitingTrack           = pconfig.get<double>("minLengthExitingTrack",100); 
      fConfig.trackVisibleEnergyThreshold     = pconfig.get<double>("trackVisibleEnergyThreshold",0.1);
      fConfig.showerVisibleEnergyThreshold    = pconfig.get<double>("showerVisibleEnergyThreshold",200); 
      fConfig.showerEnergyDistortion          = pconfig.get<double>("showerEnergyDistortion",1); 
      fConfig.trackEnergyDistortion           = pconfig.get<double>("trackEnergyDistortion",1); 
      fConfig.photonVisibleEnergyThreshold    = pconfig.get<double>("photonVisibleEnergyThreshold",100);
      fConfig.leptonEnergyDistortionContained = pconfig.get<double>("leptonEnergyDistortionContained",1); 
      fConfig.nueEfficency                    = pconfig.get<double>("nueEfficency",0.8);
      fConfig.vtxEnergyCut                    = pconfig.get<double>("vtxEnergyCut",20); 
      fConfig.photonConvLenghCut              = pconfig.get<double>("photonConvLenghCut",2); 
      fConfig.dEdxPhotonCut                   = pconfig.get<double>("dEdxPhotonCut",0.06);
      fConfig.GlobalWeight                    = pconfig.get<double>("GlobalWeight",1);
      fConfig.POTWeight                       = pconfig.get<double>("POTWeight",1);
      fConfig.CosmicWeight                    = pconfig.get<double>("CosmicWeight",0.05);
      fConfig.CosmicGlobalWeight              = pconfig.get<double>("CosmicWeight",1);
      fConfig.CosmicVolumeRadius              = pconfig.get<double>("CosmicVolumeRadius",15);
      fConfig.UniformWeights = pconfig.get<std::vector<std::string>>("UniformWeights", {});
      fConfig.ElectronWeight                  = pconfig.get<double>("ElectronWeight",1);
      fConfig.PhotonWeight                    = pconfig.get<double>("PhotonWeight",1);

      fConfig.Detector               = pconfig.get<std::string>("Detector"," "); 
      
      fConfig.ApplyNueEfficiency          = pconfig.get<bool>("ApplyNueEfficiency",false);
      fConfig.ApplyElectronEnergyCut      = pconfig.get<bool>("ApplyElectronEnergyCut",false);
      fConfig.ApplyFVCut                  = pconfig.get<bool>("ApplyFVCut",false);
      fConfig.ApplyAVCut                  = pconfig.get<bool>("ApplyAVCut",false);
      fConfig.ApplyPhotonEnergyCut        = pconfig.get<bool>("ApplyPhotonEnergyCut",false); 
      fConfig.ApplyConversionGapCut       = pconfig.get<bool>("ApplyConversionGapCut",false);
      fConfig.ApplydEdxCut                = pconfig.get<bool>("ApplydEdxCut",false);
      fConfig.ApplyMuonLenghtCut          = pconfig.get<bool>("ApplyMuonLenghtCut",false);
      fConfig.ApplyKMECCut                = pconfig.get<bool>("ApplyKMECCut",false);
      fConfig.ApplyCosmicCylinderCut      = pconfig.get<bool>("ApplyCosmicCylinderCut",false);
      fConfig.ApplyCosmicFVCut            = pconfig.get<bool>("ApplyCosmicFVCut",false);
      fConfig.ApplyCosmicInSpillWindowCut = pconfig.get<bool>("ApplyCosmicInSpillWindowCut",false);
      fConfig.UseAllCosmics               =  pconfig.get<bool>("UseAllCosmics",false);

      fConfig.ApplyGENIEVersionWeight     = pconfig.get<bool>("ApplyGENIEVersionWeight",false);

      fConfig.IgnoreNeutrinoDepsInCosmic  = pconfig.get<bool>("IgnoreNeutrinoDepsInCosmic",false);
      fConfig.UseGenieHists               = pconfig.get<bool>("UseGenieHists",false);
	
      fConfig.IncludeCosmics         = pconfig.get<bool>("IncludeCosmics",true);
      fConfig.IncludeDirt            = pconfig.get<bool>("IncludeDirt",true);  
      fConfig.CosmicsOnly            = pconfig.get<bool>("Cosmics Only",false);   
      fConfig.DirtOnly               = pconfig.get<bool>("DirtOnly",false); 
      fConfig.IntrinsicOnly          = pconfig.get<bool>("IntrinsicOnly",false);   
      fConfig.NuMuOnly               = pconfig.get<bool>("NuMuOnly",false); 
      fConfig.OscOnly                = pconfig.get<bool>("OscOnly",false);  

      fConfig.DontUseSimChannels     = pconfig.get<bool>("DontUseSimChannels",false);
      fConfig.Verbose                = pconfig.get<bool>("Verbose",false);  
      
      fConfig.GlobalTimeOffset       = pconfig.get<float>("GlobalTimeOffset",0);
      fConfig.SpillTime              = pconfig.get<float>("SpillTime",1600);
      fConfig.ReadoutWindowSize      = pconfig.get<float>("ReadoutWindowSize",1);

      fConfig.FillHistograms         = pconfig.get<bool>("FillHistograms",false);
      fConfig.FillModeHistograms     = pconfig.get<bool>("FillModeHistograms",false);
      fConfig.FillIntTypeHistograms  = pconfig.get<bool>("FillIntTypeHistograms",false);
      
      //Set up the selection histograms 
      fOutputFile->cd();

      //Initilise the histograms 
      InitialiseHistograms();

      //Time the CPU.
      //c_start = std::clock();

      MCTruthCounter =0;

      //Initialise The Genie Map
      GENIEWeight_map[1002] = {1.03493,0.0337303};
      GENIEWeight_map[1006] = {1.0319,-0.118565};
      GENIEWeight_map[1007] = {1.01131,-0.0855016};
      GENIEWeight_map[1008] = {0.964424,-0.0382905};
      GENIEWeight_map[1009] = {1.11478,-0.181611};
      GENIEWeight_map[1092] = {1.12048,-0.0810433};
      GENIEWeight_map[1096] = {1.00436,-0.190635};
      GENIEWeight_map[1098] = {-0.184961,1.37043};
      GENIEWeight_map[1001] = {0.99774,-0.0742519};
      GENIEWeight_map[1003] = {0.830306,0.0623186};
      GENIEWeight_map[1004] = {0.860698,0.00769489};
      GENIEWeight_map[1005] = {0.841456,0.0485801};
      GENIEWeight_map[1091] = {1.17868,-0.0717679};
      GENIEWeight_map[1097] = {0.854546,-0.0204999};

      //Initialising the MEC Weight map.
      // MECWeight_map["ICARUS"]["InNuE"][1001] = {1.00754,1.18254,3.51601};
      // MECWeight_map["ICARUS"]["InNuE"][1003] = {1.03578,1.64466,2.12899};
      // MECWeight_map["ICARUS"]["InNuE"][1004] = {1.06253,1.83729,2.63582};
      // MECWeight_map["ICARUS"]["InNuE"][1005] = {1.18223,-3.0788,7.77411};
      // MECWeight_map["ICARUS"]["InNuE"][1091] = {0.98294,2.63614,8.04431};
      // MECWeight_map["ICARUS"]["InNuE"][1097] = {0.979699,1.3552,0.982815};
      // MECWeight_map["ICARUS"]["InNuE"][1098] = {0.633218,2.09945,1.0884};
      // MECWeight_map["ICARUS"]["NCInNuE"][1002] = {1.02806,1.20218,2.32323};
      // MECWeight_map["ICARUS"]["NCInNuE"][1006] = {1.03922,-2.47018,7.96965};
      // MECWeight_map["ICARUS"]["NCInNuE"][1007] = {1.023,2.617,3.12625};
      // MECWeight_map["ICARUS"]["NCInNuE"][1008] = {0.988719,1.2917,1.74901};
      // MECWeight_map["ICARUS"]["NCInNuE"][1009] = {0.971067,1.52246,1.28223};
      // MECWeight_map["ICARUS"]["NCInNuE"][1013] = {0.387976,3.56248,8.32928};
      // MECWeight_map["ICARUS"]["NCInNuE"][1014] = {3.86577,-9.84088,5.15899};
      // MECWeight_map["ICARUS"]["NCInNuE"][1015] = {17.51,-14.3757,6.21929};
      // MECWeight_map["ICARUS"]["NCInNuE"][1016] = {0.425494,-0.694368,5.46637};
      // MECWeight_map["ICARUS"]["NCInNuE"][1092] = {1.01505,1.91317,1.76786};
      // MECWeight_map["ICARUS"]["NCInNuE"][1096] = {0.482165,1.84919,1.10869};
      // MECWeight_map["ICARUS"]["NCNuMu"][1002] = {1.05194,1.10104,1.74686};
      // MECWeight_map["ICARUS"]["NCNuMu"][1006] = {1.07158,1.16942,1.18191};
      // MECWeight_map["ICARUS"]["NCNuMu"][1007] = {1.0165,1.05931,1.00825};
      // MECWeight_map["ICARUS"]["NCNuMu"][1008] = {0.986577,1.2544,1.43143};
      // MECWeight_map["ICARUS"]["NCNuMu"][1009] = {1.073,0.671748,1.48902};
      // MECWeight_map["ICARUS"]["NCNuMu"][1092] = {1.04116,4.57795,7.69638};
      // MECWeight_map["ICARUS"]["NCNuMu"][1096] = {0.852231,1.22317,0.748372};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1002] = {1.00343,1.02725,1.87281};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1006] = {0.958863,1.1561,1.462};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1007] = {1.03635,1.17081,0.876149};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1008] = {1.02585,1.18397,1.17684};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1009] = {1.19696,2.89249,0.00583123};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1013] = {0.697982,2.03713,3.81921};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1014] = {0.858117,1.00461,0.630546};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1015] = {1.06073,0.796518,0.367197};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1016] = {0.843982,1.27052,0.22399};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1092] = {1.02029,1.30137,2.1216};
      // MECWeight_map["ICARUS"]["NCOscNuE"][1096] = {1.04173,1.13053,0.759037};
      // MECWeight_map["ICARUS"]["NuMu"][1001] = {1.09267,-0.974541,5.30914};
      // MECWeight_map["ICARUS"]["NuMu"][1003] = {1.02895,1.33454,2.56503};
      // MECWeight_map["ICARUS"]["NuMu"][1004] = {1.08513,1.52893,1.4254};
      // MECWeight_map["ICARUS"]["NuMu"][1005] = {1.03796,1.36873,1.22084};
      // MECWeight_map["ICARUS"]["NuMu"][1091] = {0.996049,1.38321,2.1598};
      // MECWeight_map["ICARUS"]["NuMu"][1097] = {1.02223,1.33239,0.986173};
      // MECWeight_map["ICARUS"]["OscNuE"][1001] = {1.01627,1.20791,2.14031};
      // MECWeight_map["ICARUS"]["OscNuE"][1003] = {1.01746,1.11931,3.23307};
      // MECWeight_map["ICARUS"]["OscNuE"][1004] = {1.00403,1.06331,1.31472};
      // MECWeight_map["ICARUS"]["OscNuE"][1005] = {1.02676,1.16789,1.53706};
      // MECWeight_map["ICARUS"]["OscNuE"][1091] = {1.09584,-2.74039,8.08869};
      // MECWeight_map["ICARUS"]["OscNuE"][1097] = {1.0642,1.03104,0.878908};
      // MECWeight_map["ICARUS"]["OscNuE"][1098] = {0.787406,-0.372674,3.44994};
      // MECWeight_map["SBND"]["InNuE"][1001] = {0.909349,1.33457,2.87285};
      // MECWeight_map["SBND"]["InNuE"][1003] = {0.914933,1.32031,3.94143};
      // MECWeight_map["SBND"]["InNuE"][1004] = {0.940545,1.79466,2.06287};
      // MECWeight_map["SBND"]["InNuE"][1005] = {0.866657,1.71607,3.53104};
      // MECWeight_map["SBND"]["InNuE"][1091] = {0.927632,1.08278,3.92851};
      // MECWeight_map["SBND"]["InNuE"][1097] = {0.969957,-0.450089,2.92167};
      // MECWeight_map["SBND"]["InNuE"][1098] = {0.414376,0.601776,3.9717};
      // MECWeight_map["SBND"]["NCInNuE"][1002] = {0.916797,1.34155,2.24621};
      // MECWeight_map["SBND"]["NCInNuE"][1006] = {0.89958,1.25565,1.96131};
      // MECWeight_map["SBND"]["NCInNuE"][1007] = {0.900102,1.51508,1.41904};
      // MECWeight_map["SBND"]["NCInNuE"][1008] = {0.915634,1.67886,1.68061};
      // MECWeight_map["SBND"]["NCInNuE"][1009] = {0.939751,1.4115,1.72854};
      // MECWeight_map["SBND"]["NCInNuE"][1013] = {0.720146,2.12013,1.12051};
      // MECWeight_map["SBND"]["NCInNuE"][1014] = {2.90794,12.5589,5.28556};
      // MECWeight_map["SBND"]["NCInNuE"][1015] = {0.912586,1.29938,1.06656};
      // MECWeight_map["SBND"]["NCInNuE"][1016] = {0.625069,-1.09846,4.02656};
      // MECWeight_map["SBND"]["NCInNuE"][1092] = {0.881963,2.15713,2.52658};
      // MECWeight_map["SBND"]["NCInNuE"][1096] = {0.488307,1.32607,0.939895};
      // MECWeight_map["SBND"]["NCNuMu"][1002] = {1.21721,1.02646,2.05462};
      // MECWeight_map["SBND"]["NCNuMu"][1006] = {3.43867,-11.0586,8.20568};
      // MECWeight_map["SBND"]["NCNuMu"][1007] = {1.62552,-2.46786,3.91621};
      // MECWeight_map["SBND"]["NCNuMu"][1008] = {1.12444,0.993556,1.50678};
      // MECWeight_map["SBND"]["NCNuMu"][1009] = {1.19479,0.822479,1.67179};
      // MECWeight_map["SBND"]["NCNuMu"][1092] = {1.18019,0.899514,2.45275};
      // MECWeight_map["SBND"]["NCNuMu"][1096] = {1.53533,4.05814,2.07738};
      // MECWeight_map["SBND"]["NCOscNuE"][1002] = {1.04273,1.22903,1.78029};
      // MECWeight_map["SBND"]["NCOscNuE"][1006] = {1.50084,-6.46427,8.38485};
      // MECWeight_map["SBND"]["NCOscNuE"][1007] = {0.981944,1.26284,1.42235};
      // MECWeight_map["SBND"]["NCOscNuE"][1008] = {0.988806,1.14034,1.87636};
      // MECWeight_map["SBND"]["NCOscNuE"][1009] = {1.00258,0.909995,1.25386};
      // MECWeight_map["SBND"]["NCOscNuE"][1013] = {2.47185,-6.23185,4.96047};
      // MECWeight_map["SBND"]["NCOscNuE"][1014] = {1.58532,1.498,0.433534};
      // MECWeight_map["SBND"]["NCOscNuE"][1015] = {19.3425,14.9531,5.15792};
      // MECWeight_map["SBND"]["NCOscNuE"][1016] = {1.24229,1.40831,0.461593};
      // MECWeight_map["SBND"]["NCOscNuE"][1092] = {1.01002,2.02769,2.95598};
      // MECWeight_map["SBND"]["NCOscNuE"][1096] = {1.15813,1.0215,0.408322};
      // MECWeight_map["SBND"]["NuMu"][1001] = {1.18542,1.53117,5.72593};
      // MECWeight_map["SBND"]["NuMu"][1003] = {1.22574,1.35337,3.21991};
      // MECWeight_map["SBND"]["NuMu"][1004] = {1.2316,1.00696,1.89215};
      // MECWeight_map["SBND"]["NuMu"][1005] = {1.14269,0.763421,2.73513};
      // MECWeight_map["SBND"]["NuMu"][1091] = {1.20171,1.45091,2.48061};
      // MECWeight_map["SBND"]["NuMu"][1097] = {1048.87,-16.0473,4.49884};
      // MECWeight_map["SBND"]["OscNuE"][1001] = {1.01624,1.51105,3.15648};
      // MECWeight_map["SBND"]["OscNuE"][1003] = {1.02042,1.50292,2.43135};
      // MECWeight_map["SBND"]["OscNuE"][1004] = {1.03144,1.41465,2.22496};
      // MECWeight_map["SBND"]["OscNuE"][1005] = {1.01276,0.963411,2.07326};
      // MECWeight_map["SBND"]["OscNuE"][1091] = {0.978842,1.9671,8.20263};
      // MECWeight_map["SBND"]["OscNuE"][1097] = {0.755367,0.787758,2.32375};
      // MECWeight_map["SBND"]["OscNuE"][1098] = {1.3418,1.32831,0.561326};
      // MECWeight_map["UBOONE"]["InNuE"][1001] = {0.967302,1.3232,4.26436};
      // MECWeight_map["UBOONE"]["InNuE"][1003] = {0.974279,1.05649,4.78063};
      // MECWeight_map["UBOONE"]["InNuE"][1004] = {1.00283,-0.0699605,7.8231};
      // MECWeight_map["UBOONE"]["InNuE"][1005] = {0.92884,1.07184,3.11373};
      // MECWeight_map["UBOONE"]["InNuE"][1091] = {0.962806,1.81841,2.59309};
      // MECWeight_map["UBOONE"]["InNuE"][1097] = {0.843007,0.650223,2.61095};
      // MECWeight_map["UBOONE"]["InNuE"][1098] = {0.377019,1.83915,2.09811};
      // MECWeight_map["UBOONE"]["NCInNuE"][1002] = {0.989346,1.72277,3.04771};
      // MECWeight_map["UBOONE"]["NCInNuE"][1006] = {1.03636,1.60584,1.78212};
      // MECWeight_map["UBOONE"]["NCInNuE"][1007] = {0.925134,1.85083,2.40083};
      // MECWeight_map["UBOONE"]["NCInNuE"][1008] = {0.966648,1.41976,1.90753};
      // MECWeight_map["UBOONE"]["NCInNuE"][1009] = {0.936313,1.05297,1.70373};
      // MECWeight_map["UBOONE"]["NCInNuE"][1013] = {2.06117,0.621816,0.0822705};
      // MECWeight_map["UBOONE"]["NCInNuE"][1014] = {0.590339,-1.29939,6.13645};
      // MECWeight_map["UBOONE"]["NCInNuE"][1015] = {0.738685,2.36891,2.29632};
      // MECWeight_map["UBOONE"]["NCInNuE"][1016] = {0.789607,1.96217,0.277538};
      // MECWeight_map["UBOONE"]["NCInNuE"][1092] = {0.936107,3.58765,7.65049};
      // MECWeight_map["UBOONE"]["NCInNuE"][1096] = {8.63686,18.4178,7.80626};
      // MECWeight_map["UBOONE"]["NCNuMu"][1002] = {1.02679,1.05534,2.186};
      // MECWeight_map["UBOONE"]["NCNuMu"][1006] = {1.02898,1.10175,1.30697};
      // MECWeight_map["UBOONE"]["NCNuMu"][1007] = {1.07184,1.14692,0.897739};
      // MECWeight_map["UBOONE"]["NCNuMu"][1008] = {1.06476,1.21771,1.13764};
      // MECWeight_map["UBOONE"]["NCNuMu"][1009] = {1.84531,-4.23085,4.78537};
      // MECWeight_map["UBOONE"]["NCNuMu"][1092] = {0.971414,1.34581,2.44091};
      // MECWeight_map["UBOONE"]["NCNuMu"][1096] = {0.953025,5.01469,5.24142};
      // MECWeight_map["UBOONE"]["NCNuMu"][1098] = {0.829187,0.913197,2.42932};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1002] = {1.01367,1.1796,1.85584};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1006] = {1.00531,1.19958,1.17791};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1007] = {1.10668,1.2309,1.39937};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1008] = {1.02413,1.14009,1.27173};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1009] = {1.03496,1.02342,1.18895};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1013] = {2.50422,7.36642,3.60874};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1014] = {1.01556,1.03634,0.870108};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1015] = {12.9448,-2.00325,1.23789};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1016] = {1.29207,1.54445,0.530082};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1092] = {1.0643,1.53987,1.73028};
      // MECWeight_map["UBOONE"]["NCOscNuE"][1096] = {1.17772,0.943034,0.777501};
      // MECWeight_map["UBOONE"]["NuMu"][1001] = {1.01588,0.54904,8.39826};
      // MECWeight_map["UBOONE"]["NuMu"][1003] = {1.02961,0.983058,2.29481};
      // MECWeight_map["UBOONE"]["NuMu"][1004] = {1.04928,1.43963,1.72202};
      // MECWeight_map["UBOONE"]["NuMu"][1005] = {1.02698,0.341929,2.63204};
      // MECWeight_map["UBOONE"]["NuMu"][1091] = {0.983417,1.81046,7.58807};
      // MECWeight_map["UBOONE"]["NuMu"][1097] = {0.734648,1.07835,0.711813};
      // MECWeight_map["UBOONE"]["OscNuE"][1001] = {1.02071,1.05916,3.31857};
      // MECWeight_map["UBOONE"]["OscNuE"][1003] = {1.04012,1.4165,2.54485};
      // MECWeight_map["UBOONE"]["OscNuE"][1004] = {1.00368,1.56953,1.80087};
      // MECWeight_map["UBOONE"]["OscNuE"][1005] = {0.989574,1.43803,1.39383};
      // MECWeight_map["UBOONE"]["OscNuE"][1091] = {1.03442,1.36656,1.75237};
      // MECWeight_map["UBOONE"]["OscNuE"][1097] = {1.04903,1.05837,0.829909};
      // MECWeight_map["UBOONE"]["OscNuE"][1098] = {75.6033,-13.9954,4.93722};


      if(fConfig.UseGenieHists){
	std::string GenieDiffString = "root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr/sbnd/persistent/users/dbarker/sensitivity/RatioPlots.root";
	const char* GenieDiffName = GenieDiffString.c_str();
	GenieDiffFile = TFile::Open(GenieDiffName);

	for(int int_type=1000; int_type<1101; ++int_type){
	  std::string  Hist_String  ="genieratio_" + std::to_string(int_type);
	  const char*  Hist_Name    = Hist_String.c_str();
	  if(gDirectory->Get(Hist_Name) == NULL){continue;}
	  TH1D* GenieHist = (TH1D*)(gDirectory->Get(Hist_Name))->Clone();
	  GenieHist->SetDirectory(0);
	  GenieHists[int_type] = GenieHist;
	} 
	GenieDiffFile->Close();
	fOutputFile->cd();
      }
    }
    
    
    void NueSelection::Finalize(){

      //Read the End time
      //      c_end = std::clock(); 
      //long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
      //std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

      std::cout << "Neutrinos selected: " <<  NuCount << " in " << EventCounter << " readout events" << std::endl;

      std::cout << "MCTruth Events count was: " <<  MCTruthCounter << std::endl;

      //Set up the selection histograms                                                  
      fOutputFile->cd();    
      gDirectory->mkdir("Histograms");
      gDirectory->mkdir("ModeHistograms");
      gDirectory->mkdir("InteractionTypeHistograms");
      fOutputFile->cd("Histograms");

      fEventHist->Write();

      if(fConfig.FillHistograms){
	
	fRootHists.XZCosmic->Write(); 
	fRootHists.YZCosmic->Write(); 
	fRootHists.XYCosmic->Write(); 
	fRootHists.DCACosmic->Write(); 
	fRootHists.XZPhot->Write(); 
	fRootHists.YZPhot->Write(); 
	fRootHists.XYPhot->Write(); 
	fRootHists.DCAPhot->Write(); 
	fRootHists.XZNu->Write(); 
	fRootHists.YZNu->Write(); 
	fRootHists.XYNu->Write(); 
	fRootHists.DCANu->Write();
	fRootHists.VisibleEnergy_SpillWindow_Hist->Write(); 

	
	fRootHists.VisibleEnergy_FNKPDecay_Hist->Write();
	fRootHists.VisibleEnergy_FNKMDecay_Hist->Write();
	fRootHists.VisibleEnergy_MuDecays_Hist->Write();
	
	for(int i=0; i<fRootHists.HistTypes.size(); ++i){

	  fOutputFile->cd("Histograms");
	  fRootHists.TrueNumber_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.TrueEnergy_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.TrueEnergyAll_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.CCQEEnergy_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_AVCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_FVBefore_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_FVCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_EnergyCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_PhotonEnergyCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_ConversionGapCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_MuLenghtCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_NCCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_LeptonPlusPhotonCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_Selection_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_PiZero_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_Photon_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_PhotonSmall_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.LowNCEnergy_Hist[fRootHists.HistTypes[i]]->Write();

	  fRootHists.VisibleEnergy_NoShowerCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_TwoPhotonCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_TwoShowerCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_PhotonFVCut_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_TwoLeptons_Hist[fRootHists.HistTypes[i]]->Write();

	  fRootHists.VisibleEnergyAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_AVCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_FVBeforeAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_FVCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_EnergyCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_PhotonEnergyCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_ConversionGapCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_MuLenghtCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_NCCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_LeptonPlusPhotonCutAfter_Hist[fRootHists.HistTypes[i]]->Write();

	  fRootHists.VisibleEnergy_NoShowerCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_TwoPhotonCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_TwoShowerCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_PhotonFVCutAfter_Hist[fRootHists.HistTypes[i]]->Write();


	  fRootHists.Weights_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.ProtonE_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.PionE_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.KaonE_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.ProtonN_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.PionN_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.KaonN_Hist[fRootHists.HistTypes[i]]->Write();

	  fRootHists.HadronE_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.PhotonCon_Hist[fRootHists.HistTypes[i]]->Write();

	  if(fConfig.FillModeHistograms){
	    for(int mode=-1; mode<6; ++mode){
	      fOutputFile->cd("ModeHistograms");
	      fRootHists.TrueNumber_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.TrueEnergy_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.TrueEnergyAll_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.CCQEEnergy_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_AVCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_FVCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_FVBefore_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_EnergyCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_PhotonEnergyCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ConversionGapCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MuLenghtCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_NCCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_LeptonPlusPhotonCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
 	      fRootHists.VisibleEnergy_Selection_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_PiZero_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_Photon_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_PhotonSmall_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.LowNCEnergy_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      
	      fRootHists.Weights_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      
	      fRootHists.ProtonE_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.PionE_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.KaonE_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.ProtonN_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.PionN_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.KaonN_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.HadronE_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.PhotonCon_HistMode[fRootHists.HistTypes[i]][mode]->Write();

	      fRootHists.VisibleEnergyAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_AVCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_FVCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_FVBeforeAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_EnergyCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_PhotonEnergyCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ConversionGapCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MuLenghtCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_NCCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_LeptonPlusPhotonCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();


	    }
	  }

	  if(fConfig.FillIntTypeHistograms){
	    for(int inttype=1000; inttype<1101; ++inttype){
	      fOutputFile->cd("InteractionTypeHistograms");

	      if(fRootHists.TrueNumber_HistIntType[fRootHists.HistTypes[i]][inttype]->GetEntries() == 0){continue;} 

	      fRootHists.TrueNumber_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.TrueEnergy_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.TrueEnergyAll_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.CCQEEnergy_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_AVCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_FVCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_FVBefore_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_EnergyCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_PhotonEnergyCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ConversionGapCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MuLenghtCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_NCCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_LeptonPlusPhotonCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_Selection_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_PiZero_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_Photon_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_PhotonSmall_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.LowNCEnergy_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();

	      fRootHists.Weights_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      
	      fRootHists.ProtonE_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.PionE_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.KaonE_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.ProtonN_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.PionN_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.KaonN_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.HadronE_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.PhotonCon_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();

	      fRootHists.VisibleEnergyAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_AVCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_FVCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_FVBeforeAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_EnergyCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_PhotonEnergyCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ConversionGapCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MuLenghtCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_NCCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_LeptonPlusPhotonCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      

	    }
	  }
	}

	if(fConfig.FillHistograms){
	
	  for(int i=0; i<fRootHists.HistTypes.size(); ++i){
	  
	    fOutputFile->cd("Histograms");

	    fRootHists.VisibleEnergyLep_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_AVCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_FVBefore_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_FVCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_EnergyCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_PhotonEnergyCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_ConversionGapCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_MuLenghtCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_NCCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_LeptonPlusPhotonCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_Selection_Hist[fRootHists.HistTypes[i]]->Write();


	    fRootHists.VisibleEnergyLep_TwoShowerCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_PhotonFVCut_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_TwoLeptons_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_TwoPhotonCut_Hist[fRootHists.HistTypes[i]]->Write();

	    fRootHists.VisibleEnergyLepAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_AVCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_FVBeforeAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_FVCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_EnergyCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_PhotonEnergyCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_ConversionGapCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_MuLenghtCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_NCCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_LeptonPlusPhotonCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_SelectionAfter_Hist[fRootHists.HistTypes[i]]->Write();

	    fRootHists.VisibleEnergyLep_NoShowerCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_TwoPhotonCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_TwoShowerCutAfter_Hist[fRootHists.HistTypes[i]]->Write();
	    fRootHists.VisibleEnergyLep_PhotonFVCutAfter_Hist[fRootHists.HistTypes[i]]->Write();

	    if(fConfig.FillModeHistograms){
	      for(int mode=-1; mode<6; ++mode){
		fOutputFile->cd("ModeHistograms");

		fRootHists.VisibleEnergyLep_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_AVCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_FVCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_FVBefore_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_EnergyCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_PhotonEnergyCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_ConversionGapCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_MuLenghtCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_NCCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_LeptonPlusPhotonCut_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_Selection_HistMode[fRootHists.HistTypes[i]][mode]->Write();

		fRootHists.VisibleEnergyLepAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_AVCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_FVCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_FVBeforeAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_EnergyCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_PhotonEnergyCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_ConversionGapCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_MuLenghtCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_NCCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_LeptonPlusPhotonCutAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();
		fRootHists.VisibleEnergyLep_SelectionAfter_HistMode[fRootHists.HistTypes[i]][mode]->Write();




	      }
	    }

	    if(fConfig.FillIntTypeHistograms){
	      for(int inttype=1000; inttype<1101; ++inttype){
		fOutputFile->cd("InteractionTypeHistograms");

		if(fRootHists.TrueNumber_HistIntType[fRootHists.HistTypes[i]][inttype]->GetEntries() == 0){continue;} 

		fRootHists.VisibleEnergyLep_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_AVCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_FVCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_FVBefore_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_EnergyCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_PhotonEnergyCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_ConversionGapCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_MuLenghtCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_NCCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_LeptonPlusPhotonCut_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_Selection_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();

		fRootHists.VisibleEnergyLepAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_AVCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_FVCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_FVBeforeAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_EnergyCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_PhotonEnergyCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_ConversionGapCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_MuLenghtCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_NCCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_LeptonPlusPhotonCutAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
		fRootHists.VisibleEnergyLep_SelectionAfter_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();


	      }
	    }
	  }
	}

	
	fOutputFile->cd("Histograms");
	fRootHists.VisibleEnergy_CosmicFVCut_Hist->Write();
	fRootHists.VisibleEnergy_CosmicAVCut_Hist->Write();
	fRootHists.VisibleEnergy_CosmicClyinderCut_Hist->Write();
	fRootHists.VisibleEnergy_CosmicdEdxCut_Hist->Write();
	fRootHists.VisibleEnergy_CosmicWeightCut_Hist->Write();
	fRootHists.VisibleEnergy_CosmicEnergyCut_Hist->Write();
	fRootHists.VisibleEnergy_CosmicSelection_Hist->Write();
	fRootHists.StartXZCosmic->Write();
	fRootHists.StartYZCosmic->Write(); 
	fRootHists.StartXYCosmic->Write(); 
	fRootHists.StartMCPXZCosmic->Write(); 
	fRootHists.StartMCPYZCosmic->Write(); 
	fRootHists.StartMCPXYCosmic->Write(); 
	fRootHists.EndMCPXZCosmic->Write(); 
	fRootHists.EndMCPYZCosmic->Write(); 
	fRootHists.EndMCPXYCosmic->Write(); 
	
	fRootHists.NeutrinoT0->Write();
	fRootHists.CosmicShowerT0->Write();
      }
    }

    bool NueSelection::ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco){
     
      if(fConfig.Verbose){
	if (EventCounter % 10 == 0) {
	  std::cout << "NueSelection: Processing event " << EventCounter << " "
		    << "(" << NuCount << " neutrinos selected)"
		    << std::endl;
	}

	std::cout << "New Event" << std::endl;
      }

      bool selected = false;

      EventCounter++;

      //Lets keep track 
      fEventHist->Fill(1);

      int truth_int = 0; 
      int truth_cosint = 0;

      float cosweight_tot = 0;

      //Get tracks and showers 
      auto const& mctracks  = *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
      auto const& mcshowers = *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

      //Get the flux info 
      auto const& mcflux = ev.getValidHandle<std::vector<simb::MCFlux> >(fFluxTag);

      //Get this list of track ides that actually deposited energy.
      std::map<int,double> visible_mcparticles;
      
      if(!fConfig.DontUseSimChannels){

      	auto const& mcsimchannels =  *ev.getValidHandle<std::vector<sim::SimChannel> >(fMCParticleTag);
      	for(auto const& mcsimchannel: mcsimchannels){
	  
      	  //Get the TDCIDEMap (Charge vs Time)
      	  auto tdc_ide_map = mcsimchannel.TDCIDEMap();
	  
      	  //Loop through the map     
      	  for(auto const& tdc_ide_pair : tdc_ide_map){ 
	    
      	    //Get the IDEs associated to the TDC? 
      	    auto const& ide_v = tdc_ide_pair.second;

      	    //Loop over the IDEs and add the energy. Only count from the collection plane. 
      	    for(auto const& ide : ide_v) {
      	      //account for three planes.
      	      visible_mcparticles[TMath::Abs(ide.trackID)] += ide.energy/(3*1000);
      	    }
      	  }
      	}
      }

      // //Get a map of the all the particles (not just the final state ones).
      std::map<int,const simb::MCParticle*> mcparticles;
      auto const& mcparticle_list = *ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);
      for(auto const &mcparticle:  mcparticle_list){
	
      	mcparticles[mcparticle.TrackId()] = &mcparticle;

	// if(TMath::Abs(mcparticle.PdgCode()) == 11 || TMath::Abs(mcparticle.PdgCode()) == 13){continue;}
	
	// if(mcparticle.PdgCode() == 22){
	//   if(mcparticles.find(mcparticle.Mother()) != mcparticles.end()){
	//     if(TMath::Abs(mcparticles[mcparticle.Mother()]->PdgCode()) == 13){
	//       continue;
	//     }
	//   }
	// }

	if(mcparticle.E() < 0.03){continue;}

      	if(fConfig.Verbose){
      	  std::cout << "MC Particle with track ID: " << mcparticle.TrackId() << " has pdgcode: " << mcparticle.PdgCode() << " with energy: " << mcparticle.E() << " and mother: " << mcparticle.Mother() << " start position: " <<  mcparticle.Vx() << ", " <<  mcparticle.Vy() << ", " <<  mcparticle.Vz() << std::endl;
      	}
      }

      if(fConfig.DontUseSimChannels){
      	for(auto const &mcparticle:  mcparticle_list){
      	  //int numtrajpoints = mcparticle.NumberTrajectoryPoints() - 2;
      	  //if(numtrajpoints<0){numtrajpoints = 0;}
      	  visible_mcparticles[mcparticle.TrackId()] = mcparticle.E();
      	}
      }

      if(!fConfig.DontUseSimChannels){
	  
      	//      Roll up the visible energy for shower like particles.
      	std::vector<int> rm_particles;
      	for(auto const& visible_mcparticle: visible_mcparticles){

      	  if(mcparticles.find(visible_mcparticle.first) == mcparticles.end()){
      	    rm_particles.push_back(visible_mcparticle.first);
      	    continue;
      	  }

      	  const simb::MCParticle*  mcparticle =  mcparticles[visible_mcparticle.first];

      	  if(TMath::Abs(mcparticle->PdgCode()) != 11 && TMath::Abs(mcparticle->PdgCode()) != 22){continue;}
	  
      	  //Find the mother
      	  int track_id =  mcparticle->TrackId();

      	  //Find the initial  mother
      	  int mother_id = track_id;
      	  while(mother_id != 0){
      	    if(mcparticles.find(mother_id) != mcparticles.end()){
      	      if(TMath::Abs(mcparticles[mother_id]->PdgCode()) != 11 && TMath::Abs(mcparticles[mother_id]->PdgCode()) != 22){break;}
      	      track_id = mother_id;
      	      mother_id = mcparticles[mother_id]->Mother();
      	    }
      	    else{
      	      break;
      	    }
      	  }
      	  if(track_id != mcparticle->TrackId()){
      	    visible_mcparticles[track_id] += visible_mcparticle.second;
      	    rm_particles.push_back(visible_mcparticle.first);
      	  }
      	}
      	for(auto const& rm_particle: rm_particles)
      	  visible_mcparticles.erase(rm_particle);
      }
      
      if(fConfig.Verbose){
	//for(auto const& visible_mcparticle: visible_mcparticles)
	//  std::cout << "visible track id: " << visible_mcparticle.first << " Energy: " << visible_mcparticle.second  << std::endl;
      }

      float totalbnbweight = 1;

      std::vector<simb::MCTruth> mctruth_neutrino;
      
      //Grab a neutrino datat product from the event
      for(int fTruthTag=0; fTruthTag<fTruthTags.size(); ++fTruthTag){
	
	auto const& mctruths =  *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTags[fTruthTag]);

	//Iterate through the neutrinos
	for (size_t i=0; i<mctruths.size(); i++) {

	  //lazy pea ranch
	  if(fConfig.CosmicsOnly){continue;}  
	  
	  ++MCTruthCounter;

	  auto const& mctruth = mctruths.at(i);

	  //This part is for neutrinos only.
	  if (!mctruth.NeutrinoSet()) {++truth_int; continue;}
		  
	  mctruth_neutrino.push_back(mctruth);

	  int showers=0;

	  if(fConfig.Verbose){
	    //Remove for debugging perposes
	    for (auto const &mct: mctracks) {
	      double mass = PDGMass(mct.PdgCode());
	      TLorentzVector nuVtx     = mctruths[i].GetNeutrino().Nu().Trajectory().Position(0);
	      TLorentzVector partstart = mct.Start().Position();
	      std::cout << "Track with id: " << mct.TrackID() << " has pdgcode: " << mct.PdgCode() << " and energy: " << (mct.Start().E() - mass) / 1000  << " distance from vertex: " << TMath::Abs((partstart - nuVtx).Mag()) << std::endl;
	    }
	    for (auto const &mcs: mcshowers) {
	      TLorentzVector nuVtx     = mctruths[i].GetNeutrino().Nu().Trajectory().Position(0);
	      TLorentzVector partstart = mcs.Start().Position();
	      double mass = PDGMass(mcs.PdgCode());

	      if(mcs.PdgCode() == 22 && mcs.Start().E() > 100){++showers;}

	      std::cout << "Shower with id: " << mcs.TrackID() << " has pdgcode: " << mcs.PdgCode() << " and energy: " << (mcs.Start().E() - mass) / 1000 <<  " distance from vertex: " << TMath::Abs((partstart - nuVtx).Mag()) << " and time: " << mcparticles[mcs.TrackID()]->T() << std::endl;
	    }
	  }
	  
	  
	  if(fConfig.Verbose){
	    std::cout << "##########################################################################################################" << std::endl;
	    std::cout << mctruth << std::endl;
	  }
	
	  
	  const simb::MCNeutrino& nu = mctruth.GetNeutrino();

	  //Print out the decay type. Find enums here: http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/
	  int initnu =  mcflux->at(i).fntype;
	  int fnd = mcflux->at(i).fndecay;
	
	  if(fConfig.IntrinsicOnly && (initnu != nu.Nu().PdgCode() || nu.Nu().PdgCode() != 12)){
	    // //build a temp holder so incase there is a one to one matching.
	    // event::Interaction interaction = truth[truth_int]; 
	    // event::RecoInteraction reco_interaction(truth_int);
	    // reco_interaction.reco_energy = -9999;
	    // reco_interaction.weight = -9999;
	    // reco.push_back(reco_interaction);
	    ++truth_int; 
	    continue;
	  }

	  if(fConfig.NuMuOnly && (nu.Nu().PdgCode() != 14 || initnu != nu.Nu().PdgCode())){
	    //build a temp holder so incase there is a one to one matching.
	    // event::Interaction interaction = truth[truth_int]; 
	    // event::RecoInteraction reco_interaction(truth_int);
	    // reco_interaction.reco_energy = -9999;
	    // reco_interaction.weight = -9999;
	    // reco.push_back(reco_interaction);
	    ++truth_int; 
	    continue;
	  }
	  
	  if(fConfig.OscOnly && initnu == nu.Nu().PdgCode()){

	    //build a temp holder so incase there is a one to one matching.
	    // event::Interaction interaction = truth[truth_int]; 
	    // event::RecoInteraction reco_interaction(truth_int);
	    // reco_interaction.reco_energy = -9999;
	    // reco_interaction.weight = -9999;
	    // reco.push_back(reco_interaction);
	    ++truth_int; 
	    continue;
	  }
	
	
	  //KMec didn't exist in the past might not be correct but remove for proposal
	  bool pass_kMEC = !(fConfig.ApplyKMECCut && nu.Mode() == simb::kMEC); 
	  if(!pass_kMEC){
	    if(fConfig.Verbose){std::cout << "Event Was kMEC. Event was not selected" << std::endl;}
	    //build a temp holder so incase there is a one to one matching.
	    // event::Interaction interaction = truth[truth_int]; 
	    // event::RecoInteraction reco_interaction(truth_int);
	    // reco_interaction.reco_energy = -9999;
	    // reco_interaction.weight = -9999;
	    // reco.push_back(reco_interaction);
	    ++truth_int; 
	    continue; 
	  }

	  fRootHists.NeutrinoT0->Fill(nu.Nu().T());
	  
	  
	  //holder for the electron track ID
	  std::vector<int> leptontrackIDs;
	  float lepton_energy = -99999;
	  int leptontrackID = -1; 
	
	  //Calculate the lepton track ID
	  for(auto const &mcparticle_data: visible_mcparticles){
	    const simb::MCParticle*  mcparticle = mcparticles[mcparticle_data.first];
	    float Energy = mcparticle_data.second;
	    if(TMath::Abs(mcparticle->PdgCode()) == 11 && Energy*1000 > fConfig.photonVisibleEnergyThreshold){

	      //Don't consider daughters of showering particles
	      if(mcparticles.find(mcparticle->Mother()) != mcparticles.end()){
		if(TMath::Abs(mcparticle->Mother()) == 11 || TMath::Abs(mcparticle->Mother()) == 22){continue;}
	      }

	      //Check its contained in the AV.
	      if(!containedInAV(nu.Nu().Position().Vect())){continue;}
	      
	      //Check to see if its from the vertex.
	      if(isFromNuVertex(mctruth,mcparticle)){
		
		//Say the lepton is the one with the highest energy.
		if(Energy > lepton_energy){
		  lepton_energy = Energy;
		  leptontrackID = mcparticle->TrackId();
		}
		leptontrackIDs.push_back(mcparticle->TrackId());
	      }
	    }
	  }
	  
	  if(leptontrackIDs.size() == 0){
	    leptontrackIDs.push_back(-1);
	  }
	  
	  if(fConfig.Verbose){
	    std::cout << "leptontrackID: " << leptontrackID << " lepton_energy: " << lepton_energy << std::endl;
	  }

	  //Check if the neutrino is contained.
	  bool containedinfv = containedInFV(nu.Nu().Position().Vect());
	  
	  //Set the values required to calculate the track and shower energys + distortion 
	  VisibleEnergyCalculator calculator;
	  calculator.lepton_pdgid = 11;
	  calculator.track_threshold =  fConfig.trackVisibleEnergyThreshold;
	  calculator.shower_energy_distortion = fConfig.showerEnergyDistortion;
	  calculator.track_energy_distortion = fConfig.trackEnergyDistortion;
	  calculator.lepton_energy_distortion_contained = fConfig.leptonEnergyDistortionContained;
	  calculator.lepton_energy_distortion_leaving_A = 1; //Not required here.
	  calculator.lepton_energy_distortion_leaving_B = 1; //Not required here. 
	  calculator.lepton_contained = containedinfv; // We don't worrry at the moment if the electron is contained at this moment.
	  calculator.lepton_index = leptontrackID;
	  
	  if(fConfig.Verbose){
	    if(mctruth.GetNeutrino().CCNC() == simb::kCC){
	      std::cout << "is charged current" << std::endl;
	    }
	    else if(mctruth.GetNeutrino().CCNC() == simb::kNC) {
	      std::cout << "is neutral current" << std::endl; 
	    }
	    else{
	      std::cout << "Somehow other. Maybe cosmic if you have some cosmics" << std::endl;
	    }
	  }

	  
	  //Calculate the Energy 
	  const std::vector<double> visible_energy = FlavourEnergyDeposition(rand, mctruth, mcparticles,visible_mcparticles,fConfig.active_volumes,calculator);
	  double Hadronic_energy = visible_energy[0];
	  double Other_energy    = visible_energy[1];
	  double Leptonic_energy = visible_energy[2];
	  double Shower_energy   = 0;
	  //	double Nue_energy = Hadronic_energy + Shower_energy + Leptonic_energy;

	  //Weighting for efficiency of selection.
	  double weight = 1;
	  
	  //Make a nue interaction 
	  NueSelection::NueInteraction intInfo({Hadronic_energy, Shower_energy, Leptonic_energy, Other_energy, weight,initnu,leptontrackID,11,fnd}); 

	  ///Weight with the POT
	  intInfo.weight *= fConfig.POTWeight;
	  
	  //Add the Global weightings if any.
	  intInfo.weight *= fConfig.GlobalWeight;

	  //Apply Uniform weight for bnb correction 
	  if(fConfig.Detector != "SBND"){
	    // for (auto const &key: fConfig.UniformWeights) {
	    //   //	      intInfo.weight *= interaction.weightmap.at(key)[0];
	    //   totalbnbweight *=  interaction.weightmap.at(key)[0];     
	    //   if(fConfig.Verbose){
	    // 	std::cout << "Flux correction is: " <<  interaction.weightmap.at(key)[0] << std::endl;
	    //   }

	    // }
	  }

	  //Add the weight for removing MEC events if add.
	  if(fConfig.ApplyKMECCut){
	    //intInfo.weight *= MECWeight(nu,fConfig.Detector,intInfo);
	  }

	  //Add the weight to get back to genie v02_08
	  if(fConfig.ApplyGENIEVersionWeight){
	    intInfo.weight *= GENIEWeight(nu);
	  }

	  cosweight_tot +=  intInfo.weight;

	  dirtevent = false;

	  	  int nprotons = 0;
	  int npions   = 0; 
	  int nkaons   = 0;
	  // //Get a map of the all the particles (not just the final state ones).
	  double HadronE = 0;
	  for(auto const &mcparticleit:  mcparticles){
	    const simb::MCParticle*  mcparticle = mcparticleit.second;

	    //Count the protons kaons and pions and there energy. 
	    if(isFromNuVertex(mctruth,mcparticle)){

	      double mass = PDGMass(mcparticle->PdgCode());
	      double KE = (mcparticle->E()*1000-mass)/1000;
	      if(TMath::Abs(mcparticle->PdgCode()) == 2212){
		++nprotons;
		FillHistograms(fRootHists.ProtonE_Hist,fRootHists.ProtonE_HistMode,fRootHists.ProtonE_HistIntType,nu,intInfo,KE,dirtevent);
		if(KE > 0.021) HadronE += KE;
	      }
	      else if(TMath::Abs(mcparticle->PdgCode()) == 211){
		++npions;
		FillHistograms(fRootHists.PionE_Hist,fRootHists.PionE_HistMode,fRootHists.PionE_HistIntType,nu,intInfo,KE,dirtevent);
		HadronE += KE;
	      }
	      else if(TMath::Abs(mcparticle->PdgCode()) == 321){
		++nkaons;
		FillHistograms(fRootHists.KaonE_Hist,fRootHists.KaonE_HistMode,fRootHists.KaonE_HistIntType,nu,intInfo,mcparticle->E(),dirtevent);
		HadronE += mcparticle->E();
	      }
	    }
	  }
	  FillHistograms(fRootHists.ProtonN_Hist,fRootHists.ProtonN_HistMode,fRootHists.ProtonN_HistIntType,nu,intInfo,nprotons,dirtevent);
	  FillHistograms(fRootHists.PionN_Hist,fRootHists.PionN_HistMode,fRootHists.PionN_HistIntType,nu,intInfo,npions,dirtevent);
	  FillHistograms(fRootHists.KaonN_Hist,fRootHists.KaonN_HistMode,fRootHists.KaonN_HistIntType,nu,intInfo,nkaons,dirtevent);
	  FillHistograms(fRootHists.HadronE_Hist,fRootHists.HadronE_HistMode,fRootHists.HadronE_HistIntType,nu,intInfo,HadronE,dirtevent);
	  

	  //Fill the number graphs
	  FillHistograms(fRootHists.TrueNumber_Hist,fRootHists.TrueNumber_HistMode,fRootHists.TrueNumber_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	  FillHistograms(fRootHists.TrueEnergyAll_Hist,fRootHists.TrueEnergyAll_HistMode,fRootHists.TrueEnergyAll_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	  FillHistograms(fRootHists.Weights_Hist,fRootHists.Weights_HistMode,fRootHists.Weights_HistIntType,nu,intInfo,intInfo.weight,dirtevent);  

	  //Selection Criteria 
	  //Check for Tracks - Remove any track that come from the neutrino if 1m 
	  //Find electron vertex candiates - Remove 80% due to fake dEdx cut 
	  //If pi0 is made 2 photons come out as candiates in both are above 100 MeV remove.
	  //If hadrons produce have more than 50 MeV this is visible if all photons pair produce more than 3cm away from the vertex remove.
	  //Remove 94% pion of events as dEdx cut
	  //Event could be CC0Pi Muons Interaction with photon background which mimics the CC1Pi Electron interatactions. these need to go throught the cuts above as well.  
	  //Run the Selection
	  
	  //Selection. 
	  bool selection = Select(ev, mctruth, truth_int, mcparticles, visible_mcparticles, intInfo, leptontrackIDs);
	  
	  if(selection){

	    if((int) intInfo.leptonpdgID == 11){
	      intInfo.weight *= fConfig.ElectronWeight;
	    }
	    else if ((int)intInfo.leptonpdgID == 22){
	      intInfo.weight *= fConfig.PhotonWeight;
	    }

	    
	    //Make the reco interaction event. 
	    event::RecoInteraction reco_interaction(truth_int);
	    reco_interaction.reco_energy = intInfo.GetNueEnergy();

	    // reco_interaction.reco_energy = nu.Nu().E();
	    
	    reco_interaction.weight = intInfo.weight;
	    if(dirtevent){ 
	      reco_interaction.wasDirt = true;
	    }

	    reco.push_back(reco_interaction);
	    selected =  selection;

	    //Fill the histograms
	    FillHistograms(fRootHists.VisibleEnergy_Selection_Hist,fRootHists.VisibleEnergy_Selection_HistMode,fRootHists.VisibleEnergy_Selection_HistIntType,nu,intInfo,dirtevent);
	    FillHistograms(fRootHists.VisibleEnergy_Selection_Hist,fRootHists.VisibleEnergy_Selection_HistMode,fRootHists.VisibleEnergy_Selection_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);

	    FillHistograms(fRootHists.TrueEnergy_Hist,fRootHists.TrueEnergy_HistMode,fRootHists.TrueEnergy_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	    //	    FillHistograms(fRootHists.Weights_Hist,nu,intInfo,intInfo.weight,dirtevent);
	    
	    NuCount++;
	    
	    if(nu.Nu().PdgCode() == intInfo.initnu){
	      if(intInfo.fnd == 5 || intInfo.fnd == 6 || intInfo.fnd == 7){fRootHists.VisibleEnergy_FNKPDecay_Hist->Fill(intInfo.GetNueEnergy());}
	      if(intInfo.fnd == 8 || intInfo.fnd == 9 || intInfo.fnd == 10){fRootHists.VisibleEnergy_FNKMDecay_Hist->Fill(intInfo.GetNueEnergy());}  
	      if(intInfo.fnd == 11 || intInfo.fnd == 12){fRootHists.VisibleEnergy_MuDecays_Hist->Fill(intInfo.GetNueEnergy());} 
	    }
	  }
	  else{
	    // event::RecoInteraction reco_interaction(truth_int);
	    // reco_interaction.reco_energy = -9999;
	    // reco_interaction.weight = -9999;
	    // reco.push_back(reco_interaction);
	  }
	  
	  if(fConfig.Verbose){
	    
	    if(selection){
	      std::cout << "Event Was Selected with energy: " << intInfo.GetNueEnergy()<< "and weight" << intInfo.weight << std::endl;
	      if(mctruth.GetNeutrino().CCNC() == simb::kNC){
		if(intInfo.weight > 0.1){std::cout << "dodgy weight" << std::endl;}
		if(intInfo.GetNueEnergy() < 0.35){
		  FillHistograms(fRootHists.LowNCEnergy_Hist,fRootHists.LowNCEnergy_HistMode,fRootHists.LowNCEnergy_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
		}
	      }
	      if(mctruth.GetNeutrino().CCNC() == simb::kNC && showers>2){std::cout << "More than one shower" << std::endl;}
	    }
	    
	    PrintInformation(mctruth, intInfo);
	  }
	  FillHistograms(fRootHists.TrueEnergyAll_Hist,fRootHists.TrueEnergyAll_HistMode,fRootHists.TrueEnergyAll_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	  ++truth_int;
	}
	
	//########################
	//### Cosmic Selection ###
	//########################

	//I really don't want to split but lets split  
	if(fConfig.IncludeCosmics){

	
	  //Iterate through the neutrinos
	  for (int j=0; j<mctruths.size(); ++j) {
	    
	    //This is for cosmics only 
	    auto const& mctruth = mctruths.at(j);
	    if (mctruth.NeutrinoSet()){ ++truth_cosint; continue;} 
	    
	    if(cosweight_tot == 0){ ++truth_cosint; continue;}

	    for(auto const &mcs: mcshowers){ 

	      //Weighting for efficiency of selection. Start with weights on the neutrino as we are scaling with the POT.
	      double weightcos = cosweight_tot;

	      //Add the Global weightings if any.
	      weightcos *= fConfig.CosmicGlobalWeight;
	      
	      //Account for bnb weights?
	      //	      weightcos *= totalbnbweight;
	      
	      //	      fRootHists.CosmicShowerT0->Fill(mcparticles[mcs.TrackID()]->T());

	      //Make a nue interaction 
	      NueSelection::NueInteraction intInfo({0, 0, 0, 0, weightcos,-9999,-9999,-9999,-9999}); 
	      bool selection = SelectCosmic(mcs, mcparticles, intInfo, mctruth_neutrino);

	      if(selection){
	    	//Make the reco interaction event. 
	    	event::RecoInteraction reco_interaction(truth_cosint);
		reco_interaction.wasCosmic = true;
	    	reco_interaction.reco_energy = intInfo.GetNueEnergy();
	    	reco_interaction.weight = intInfo.weight; 
	    	std::cout << " Cosmic Selected with track id: " << mcs.TrackID() << " energy: " << reco_interaction.reco_energy << " weight: " << reco_interaction.weight  << std::endl;
	    	reco.push_back(reco_interaction);
	    	selected =  selection;
	    	NuCount++;
	      }
	    }
	    ++truth_cosint;
	  }
	}
      }
      return selected;
    }
    
    bool NueSelection::Select(const gallery::Event& ev, const simb::MCTruth& mctruth,  unsigned i, std::map<int, const simb::MCParticle*>& mcparticles,std::map<int,double>& visible_mcparticles, NueSelection::NueInteraction& intInfo, std::vector<int>& leptontrackIDs){


      //Neutrino Info 
      const simb::MCNeutrino& nu = mctruth.GetNeutrino();

      // bool isCC;
      // if(mctruth.GetNeutrino().CCNC() == simb::kCC){
      // 	isCC = true;
      // }
      // else{
      // 	isCC = false;
      // }
      
      //######################
      //### Reco Efficency ###
      //######################

      //if(isCC){
      //Check to see if the event passes the 80% ID efficiency cut from the proposal
      //bool pass_NueIDEfficiencyCut = NueIDEfficiencyCut(nu);
      // if(!pass_NueIDEfficiencyCut){
      intInfo.weight *= fConfig.nueEfficency;

      //#########################
      //### Active Volume Cut ###
      //#########################

      //Resetting direvent bool
      dirtevent = false;
      // pass activevolume cut                               
      FillHistograms(fRootHists.VisibleEnergy_AVCut_Hist,fRootHists.VisibleEnergy_AVCut_HistMode,fRootHists.VisibleEnergy_AVCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
      FillHistograms(fRootHists.VisibleEnergyLep_AVCut_Hist,fRootHists.VisibleEnergyLep_AVCut_HistMode,fRootHists.VisibleEnergyLep_AVCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 

      fRootHists.XZNu->Fill(nu.Nu().Position().Vect().X(),nu.Nu().Position().Vect().Z());
      fRootHists.YZNu->Fill(nu.Nu().Position().Vect().Y(),nu.Nu().Position().Vect().Z());
      fRootHists.XYNu->Fill(nu.Nu().Position().Vect().X(),nu.Nu().Position().Vect().Y());
      fRootHists.DCANu->Fill(DistanceToClosestSurface(nu.Nu().Position().Vect()));

      bool pass_AV = passAV(nu.Nu().Position().Vect());
      if(!pass_AV){
	if(fConfig.Verbose){std::cout << "Failed the AV Cut" << std::endl;}

	if(fConfig.IncludeDirt){dirtevent = true;}
	else{return false;}
      }
      else{
	if(fConfig.DirtOnly){return false;}
      }
      if(fConfig.Verbose){std::cout << "Passed the AV Cut" << std::endl;}
      FillHistograms(fRootHists.VisibleEnergy_AVCutAfter_Hist,fRootHists.VisibleEnergy_AVCutAfter_HistMode,fRootHists.VisibleEnergy_AVCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
      FillHistograms(fRootHists.VisibleEnergyLep_AVCutAfter_Hist,fRootHists.VisibleEnergyLep_AVCutAfter_HistMode,fRootHists.VisibleEnergyLep_AVCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 


      FillHistograms(fRootHists.VisibleEnergy_TwoLeptons_Hist,fRootHists.VisibleEnergy_TwoLeptons_HistMode,fRootHists.VisibleEnergy_TwoLeptons_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
      FillHistograms(fRootHists.VisibleEnergyLep_TwoLeptons_Hist,fRootHists.VisibleEnergyLep_TwoLeptons_HistMode,fRootHists.VisibleEnergyLep_TwoLeptons_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 

      //If there two or more lepton candiates we remove the events.
      if(leptontrackIDs.size() > 1){
	if(fConfig.Verbose){
	  std::cout << "Too many lepton candiadates. removing" << std::endl;
	}
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_TwoLeptonsAfter_Hist,fRootHists.VisibleEnergy_TwoLeptonsAfter_HistMode,fRootHists.VisibleEnergy_TwoLeptonsAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
      FillHistograms(fRootHists.VisibleEnergyLep_TwoLeptonsAfter_Hist,fRootHists.VisibleEnergyLep_TwoLeptonsAfter_HistMode,fRootHists.VisibleEnergyLep_TwoLeptonsAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 

      
      //#######################
      //### Beam Backgrounds ##
      //#######################
      
      //Check the Neutral current cuts. Get the photons from pi0s
      std::vector<int> pi_zeros =  findNeutralPions(mcparticles, mctruth);

      if(pi_zeros.size() > 0){
	for(int pi=0; pi<pi_zeros.size(); ++pi){

	  const simb::MCParticle* pizero = mcparticles[pi_zeros[pi]];
	  FillHistograms(fRootHists.VisibleEnergy_PiZero_Hist,fRootHists.VisibleEnergy_PiZero_HistMode,fRootHists.VisibleEnergy_PiZero_HistIntType,nu,intInfo,pizero->E(),dirtevent);
	}
      }


      //Check the number of photons 
      std::vector<int> photons  =  findPhotons(pi_zeros, mcparticles, mctruth, visible_mcparticles);

      float EBig   = -999;
      float ESmall = -999;
      float BigID = -1;
      float convdist = -999;
      //      float SmallID = -1;
      for(int ph=0; ph<photons.size(); ++ph){
	if(mcparticles[photons[ph]]->E() > EBig){
	  EBig   = mcparticles[photons[ph]]->E(); 
	  BigID  = mcparticles[photons[ph]]->TrackId();
	  convdist =  (mcparticles[BigID]->EndPosition().Vect() -  nu.Nu().Trajectory().Position(0).Vect()).Mag(); 
	}

	TVector3 EndPos   = mcparticles[photons[ph]]->EndPosition().Vect(); 
	for(int i=0; i<mcparticles[photons[ph]]->NumberTrajectoryPoints();++i){
	  if(mcparticles[photons[ph]]->E(i) < 0.1*mcparticles[photons[ph]]->E()){
	    EndPos   = mcparticles[photons[ph]]->Position(i).Vect();
	    break;
	  }
	}

	fRootHists.XZPhot->Fill(EndPos.X(),EndPos.Z());
	fRootHists.YZPhot->Fill(EndPos.Y(),EndPos.Z());
	fRootHists.XYPhot->Fill(EndPos.X(),EndPos.Y());
	fRootHists.DCAPhot->Fill(DistanceToClosestSurface(EndPos));

      }
      for(int ph=0; ph<photons.size(); ++ph){
	if(mcparticles[photons[ph]]->TrackId() == BigID){continue;}
	if(mcparticles[photons[ph]]->E() > ESmall){
	  ESmall   = mcparticles[photons[ph]]->E(); 
	  //	  SmallID = mcparticles[photons[ph]]->TrackId();
	}
      }


      FillHistograms(fRootHists.VisibleEnergy_Photon_Hist,fRootHists.VisibleEnergy_Photon_HistMode,fRootHists.VisibleEnergy_Photon_HistIntType,nu,intInfo,EBig,dirtevent);
      FillHistograms(fRootHists.VisibleEnergy_PhotonSmall_Hist,fRootHists.VisibleEnergy_PhotonSmall_HistMode,fRootHists.VisibleEnergy_PhotonSmall_HistIntType,nu,intInfo,ESmall,dirtevent);
      FillHistograms(fRootHists.PhotonCon_Hist,fRootHists.PhotonCon_HistMode,fRootHists.PhotonCon_HistIntType,nu,intInfo,convdist,dirtevent);


      //Calculate the visible photon energy
      //intInfo.shower_energy = PhotonVisibleEnergy(mcparticles,photons); 
      intInfo.shower_energy = 0;

      if(fConfig.Verbose){std::cout << "Number of pi_zeros: " << pi_zeros.size() << " Number of photons: " << photons.size() << std::endl;}

      
      FillHistograms(fRootHists.VisibleEnergy_NoShowerCut_Hist,fRootHists.VisibleEnergy_NoShowerCut_HistMode,fRootHists.VisibleEnergy_NoShowerCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
      FillHistograms(fRootHists.VisibleEnergyLep_NoShowerCut_Hist,fRootHists.VisibleEnergyLep_NoShowerCut_HistMode,fRootHists.VisibleEnergyLep_NoShowerCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 


      //Check to see if there is actually a shower
      if(photons.size() == 0 && intInfo.leptontrackID < 0){
	if(fConfig.Verbose){std::cout << "No lepton or photon = no shower. Event Not Selected" << std::endl;}
	return false;
      }

      FillHistograms(fRootHists.VisibleEnergy_NoShowerCutAfter_Hist,fRootHists.VisibleEnergy_NoShowerCutAfter_HistMode,fRootHists.VisibleEnergy_NoShowerCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
      FillHistograms(fRootHists.VisibleEnergyLep_NoShowerCutAfter_Hist,fRootHists.VisibleEnergyLep_NoShowerCutAfter_HistMode,fRootHists.VisibleEnergyLep_NoShowerCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 


      //Single photons can mimic the Nue CC events they come from NC events and misidentified Numu CC events.
      if(photons.size() > 0){

	
	FillHistograms(fRootHists.VisibleEnergy_TwoPhotonCut_Hist,fRootHists.VisibleEnergy_TwoPhotonCut_HistMode,fRootHists.VisibleEnergy_TwoPhotonCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
	FillHistograms(fRootHists.VisibleEnergyLep_TwoPhotonCut_Hist,fRootHists.VisibleEnergyLep_TwoPhotonCut_HistMode,fRootHists.VisibleEnergyLep_TwoPhotonCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 
	
	int photon_trackID = -1;
	//Check to see if there are two photons are above 200 MeV and so are detectable.
	bool pass_photonEnergyCut = passPhotonEnergyCut(photons,mcparticles,photon_trackID,visible_mcparticles);


	if(!pass_photonEnergyCut){
	  if(fConfig.Verbose){std::cout << "Failed the Photon Energy Cut. There are more than one photons in the event. Event not Selected" << std::endl;}
	  return false;
	}
      	if(fConfig.Verbose){std::cout << "Passed the Photon Energy Cut" << std::endl;}

	FillHistograms(fRootHists.VisibleEnergy_TwoPhotonCutAfter_Hist,fRootHists.VisibleEnergy_TwoPhotonCutAfter_HistMode,fRootHists.VisibleEnergy_TwoPhotonCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
	FillHistograms(fRootHists.VisibleEnergyLep_TwoPhotonCutAfter_Hist,fRootHists.VisibleEnergyLep_TwoPhotonCutAfter_HistMode,fRootHists.VisibleEnergyLep_TwoPhotonCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 


	FillHistograms(fRootHists.VisibleEnergy_TwoShowerCut_Hist,fRootHists.VisibleEnergy_TwoShowerCut_HistMode,fRootHists.VisibleEnergy_TwoShowerCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
	FillHistograms(fRootHists.VisibleEnergyLep_TwoShowerCut_Hist,fRootHists.VisibleEnergyLep_TwoShowerCut_HistMode,fRootHists.VisibleEnergyLep_TwoShowerCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 


	//Unfortunatly we must remove CC events where a photon exists (at this stage one photon exists)
	if(intInfo.leptontrackID > 0 && intInfo.leptonic_energy*1000 > fConfig.showerVisibleEnergyThreshold){
	  FillHistograms(fRootHists.VisibleEnergy_LeptonPlusPhotonCut_Hist,fRootHists.VisibleEnergy_LeptonPlusPhotonCut_HistMode,fRootHists.VisibleEnergy_LeptonPlusPhotonCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	  if(fConfig.Verbose){std::cout << "Failed becuase lepton + shower exists. Event Not Selected" << std::endl;}
	  return false;
	}
	if(fConfig.Verbose){std::cout << "Passed the lepton + photon cut" << std::endl;}

	FillHistograms(fRootHists.VisibleEnergy_TwoShowerCutAfter_Hist,fRootHists.VisibleEnergy_TwoShowerCutAfter_HistMode,fRootHists.VisibleEnergy_TwoShowerCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
	FillHistograms(fRootHists.VisibleEnergyLep_TwoShowerCutAfter_Hist,fRootHists.VisibleEnergyLep_TwoShowerCutAfter_HistMode,fRootHists.VisibleEnergyLep_TwoShowerCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 

	
	FillHistograms(fRootHists.VisibleEnergy_PhotonFVCut_Hist,fRootHists.VisibleEnergy_PhotonFVCut_HistMode,fRootHists.VisibleEnergy_PhotonFVCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
	FillHistograms(fRootHists.VisibleEnergyLep_PhotonFVCut_Hist,fRootHists.VisibleEnergyLep_PhotonFVCut_HistMode,fRootHists.VisibleEnergyLep_PhotonFVCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 


	//No photon was matched in the FV the track id is not set remove event.
	if(photon_trackID < 0 && intInfo.leptontrackID < 0){
	  if(fConfig.Verbose){std::cout << "No Photon in FV. Event not selected" << std::endl;}
	  return false;
	}

	FillHistograms(fRootHists.VisibleEnergy_PhotonFVCutAfter_Hist,fRootHists.VisibleEnergy_PhotonFVCutAfter_HistMode,fRootHists.VisibleEnergy_PhotonFVCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
	FillHistograms(fRootHists.VisibleEnergyLep_PhotonFVCutAfter_Hist,fRootHists.VisibleEnergyLep_PhotonFVCutAfter_HistMode,fRootHists.VisibleEnergyLep_PhotonFVCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 

	//	std::cout << "photon_trackID: " << photon_trackID << " intInfo.leptontrackID: " << intInfo.leptontrackID << std::endl;
   
	if(photons.size() == 0){std::cout << "no photons when there should be" << std::endl;}

	//Continue if one photon has been matched if we havn't matched then we are a CC event with no visible photons in the FV.
	if(photon_trackID !=-99999){

	  //If we have one photon candandiate choose that to be the energy of the lepton 
	  VisibleEnergyCalculator calculator;
	  calculator.lepton_pdgid = 22;
	  calculator.track_threshold =  fConfig.trackVisibleEnergyThreshold;
	  calculator.shower_energy_distortion = fConfig.showerEnergyDistortion;
	  calculator.track_energy_distortion = fConfig.trackEnergyDistortion;
	  calculator.lepton_energy_distortion_contained = fConfig.leptonEnergyDistortionContained;
	  calculator.lepton_energy_distortion_leaving_A = 1; //Not required here.
	  calculator.lepton_energy_distortion_leaving_B = 1; //Not required here. 
	  calculator.lepton_contained = true; // We don't worry at the moment if the electron is contained at this moment.
	  calculator.lepton_index = photon_trackID;

	  //	  int numtrajpoints = mcparticles[photon_trackID]->NumberTrajectoryPoints() - 2;
	  //if(numtrajpoints<0){numtrajpoints = 0;}
	  //double photonenergy = 0;
	  //	  photonenergy = mcparticles[photon_trackID]->E(numtrajpoints);
	   
	  double photonenergy = visible_mcparticles[photon_trackID];
	  
	  
	  if(fConfig.Verbose){std::cout << "Photon Visible energy: " << photonenergy << std::endl;}

	  //Add the energy 
	  intInfo.leptonic_energy = smearLeptonEnergy(rand,photonenergy,calculator.lepton_pdgid,calculator);
	  intInfo.leptonpdgID = 22; 

	  FillHistograms(fRootHists.VisibleEnergy_PhotonEnergyCut_Hist,fRootHists.VisibleEnergy_PhotonEnergyCut_HistMode,fRootHists.VisibleEnergy_PhotonEnergyCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	  FillHistograms(fRootHists.VisibleEnergyLep_PhotonEnergyCut_Hist,fRootHists.VisibleEnergyLep_PhotonEnergyCut_HistMode,fRootHists.VisibleEnergyLep_PhotonEnergyCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);

	  if(!dirtevent){

	    FillHistograms(fRootHists.VisibleEnergy_ConversionGapCut_Hist,fRootHists.VisibleEnergy_ConversionGapCut_HistMode,fRootHists.VisibleEnergy_ConversionGapCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	    FillHistograms(fRootHists.VisibleEnergyLep_ConversionGapCut_Hist,fRootHists.VisibleEnergyLep_ConversionGapCut_HistMode,fRootHists.VisibleEnergyLep_ConversionGapCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);


	    //If hadrons produce have more than 50 MeV vertex is visible. If all photons pair produce more than 3cm away from the vertex remove.
	    bool pass_conversionGapCut  = passConversionGapCut(mcparticles,photon_trackID,intInfo.hadronic_energy,nu);
	    if(!pass_conversionGapCut){
	      if(fConfig.Verbose){std::cout << "Failed the photon conversion gap cut. Event not selected." << std::endl;}
	      return false;
	    }
	    if(fConfig.Verbose){std::cout << "Passed the Photon Conversion gap Cut" << std::endl;}
	    
	    FillHistograms(fRootHists.VisibleEnergy_ConversionGapCutAfter_Hist,fRootHists.VisibleEnergy_ConversionGapCutAfter_HistMode,fRootHists.VisibleEnergy_ConversionGapCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	    FillHistograms(fRootHists.VisibleEnergyLep_ConversionGapCutAfter_Hist,fRootHists.VisibleEnergyLep_ConversionGapCutAfter_HistMode,fRootHists.VisibleEnergyLep_ConversionGapCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);

	 
	    FillHistograms(fRootHists.VisibleEnergy_MuLenghtCut_Hist,fRootHists.VisibleEnergy_MuLenghtCut_HistMode,fRootHists.VisibleEnergy_MuLenghtCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
	    FillHistograms(fRootHists.VisibleEnergyLep_MuLenghtCut_Hist,fRootHists.VisibleEnergyLep_MuLenghtCut_HistMode,fRootHists.VisibleEnergyLep_MuLenghtCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 

   
	    //Check if it is a numu CC where the muon is less than 1m
	    bool pass_muLenghtCut = passMuLengthCut(mcparticles, mctruth);
	    if(!pass_muLenghtCut){
	      if(fConfig.Verbose){std::cout << "Failed the muon length cut. Event not selected" << std::endl;}
	      return false;
	    }
	    if(fConfig.Verbose){std::cout << "Passed the muon min length Cut" << std::endl;}
	    FillHistograms(fRootHists.VisibleEnergy_MuLenghtCutAfter_Hist,fRootHists.VisibleEnergy_MuLenghtCutAfter_HistMode,fRootHists.VisibleEnergy_MuLenghtCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent); 
	    FillHistograms(fRootHists.VisibleEnergyLep_MuLenghtCutAfter_Hist,fRootHists.VisibleEnergyLep_MuLenghtCutAfter_HistMode,fRootHists.VisibleEnergyLep_MuLenghtCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent); 

	  }
	
	
	  //Remove 94% on a dEdx cut                                       
	  FillHistograms(fRootHists.VisibleEnergy_NCCut_Hist,fRootHists.VisibleEnergy_NCCut_HistMode,fRootHists.VisibleEnergy_NCCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	  FillHistograms(fRootHists.VisibleEnergyLep_NCCut_Hist,fRootHists.VisibleEnergyLep_NCCut_HistMode,fRootHists.VisibleEnergyLep_NCCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);

	  bool pass_dEdxCut = passdEdxCut(22);
	  if(!pass_dEdxCut){intInfo.weight *= fConfig.dEdxPhotonCut;}
	  FillHistograms(fRootHists.VisibleEnergy_NCCutAfter_Hist,fRootHists.VisibleEnergy_NCCutAfter_HistMode,fRootHists.VisibleEnergy_NCCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
	  FillHistograms(fRootHists.VisibleEnergyLep_NCCutAfter_Hist,fRootHists.VisibleEnergyLep_NCCutAfter_HistMode,fRootHists.VisibleEnergyLep_NCCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);

	}
      }

      //##########################
      //### Fiducal Volume cut ###
      //##########################
      
      FillHistograms(fRootHists.VisibleEnergy_FVBefore_Hist,fRootHists.VisibleEnergy_FVBefore_HistMode,fRootHists.VisibleEnergy_FVBefore_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
      FillHistograms(fRootHists.VisibleEnergy_FVCut_Hist,fRootHists.VisibleEnergy_FVCut_HistMode,fRootHists.VisibleEnergy_FVCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
      FillHistograms(fRootHists.VisibleEnergyLep_FVCut_Hist,fRootHists.VisibleEnergyLep_FVCut_HistMode,fRootHists.VisibleEnergyLep_FVCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);

    // pass fiducial volume cut (already passed on photon events. )                                             
      if(intInfo.leptonpdgID == 11){
	bool pass_FV = passFV(nu.Nu().Position().Vect());
	if(!pass_FV){
	  if(fConfig.Verbose){std::cout << "Failed the FV cut. Event not selected wit position X:" << nu.Nu().Position().Vect().X() << " Y: " << nu.Nu().Position().Vect().Y() << " Z: " << nu.Nu().Position().Vect().Z() << std::endl;}
	  return false;
	}
	if(fConfig.Verbose){std::cout << "Passed the FV Cut" << std::endl;}
	//Should have no dirts here.

      }
      FillHistograms(fRootHists.VisibleEnergy_FVCutAfter_Hist,fRootHists.VisibleEnergy_FVCutAfter_HistMode,fRootHists.VisibleEnergy_FVCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
      FillHistograms(fRootHists.VisibleEnergyLep_FVCutAfter_Hist,fRootHists.VisibleEnergyLep_FVCutAfter_HistMode,fRootHists.VisibleEnergyLep_FVCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);

      //##################
      //### Energy Cut ###
      //##################

      FillHistograms(fRootHists.VisibleEnergy_EnergyCut_Hist,fRootHists.VisibleEnergy_EnergyCut_HistMode,fRootHists.VisibleEnergy_EnergyCut_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
      FillHistograms(fRootHists.VisibleEnergyLep_EnergyCut_Hist,fRootHists.VisibleEnergyLep_EnergyCut_HistMode,fRootHists.VisibleEnergyLep_EnergyCut_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);


      //Check the electron is not removed due to the 200 MeV cut. Just check we havn't matched 
      bool pass_eEnergyCut = passeEnergyCut(intInfo.leptonic_energy);
      if(!pass_eEnergyCut){
	if(fConfig.Verbose){std::cout << "Shower was too low in energy. Event not Selected" << std::endl;}
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_EnergyCutAfter_Hist,fRootHists.VisibleEnergy_EnergyCutAfter_HistMode,fRootHists.VisibleEnergy_EnergyCutAfter_HistIntType,nu,intInfo,nu.Nu().E(),dirtevent);
      FillHistograms(fRootHists.VisibleEnergyLep_EnergyCutAfter_Hist,fRootHists.VisibleEnergyLep_EnergyCutAfter_HistMode,fRootHists.VisibleEnergyLep_EnergyCutAfter_HistIntType,nu,intInfo,nu.Lepton().E(),dirtevent);

      if(fConfig.Verbose){std::cout << "Passed the Energy Cut." << std::endl;}
	
	return true;
    
    }
   
    //Check if the point is the fiducal volume.
    bool NueSelection::containedInFV(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& FV: fConfig.fiducial_volumes) {
	if (FV.Contain(p)) return true;
      }
      return false;
    }

    //Check if the point is in the Active volume.
    bool NueSelection::containedInAV(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& AV: fConfig.active_volumes) {
	if (AV.Contain(p)) return true;
      }
      return false;
    }


    //Is the neutrino a electron neutrino if so apply ID efficiency. 
    bool NueSelection::NueIDEfficiencyCut(const simb::MCNeutrino& nu){
      
      if(nu.Nu().PdgCode() == 12){
	  return true;
      }
      
      return false;
    }


    //Remove Events where the lepton is below 200 MeV in energy.
    bool NueSelection::eEnergyCut(double leptonenergy){
      
      if(fConfig.ApplyElectronEnergyCut){
	if(leptonenergy*1000 < fConfig.showerVisibleEnergyThreshold){
	  return false;
	}
      }
      return true;
    }

    //Weight events by the dEdx cut.
    bool NueSelection::dEdxCut(int pdgcode){
      
      //Apply weight for photon evenets 
      if(pdgcode==22){return false;}
      
      return true;
    }

    //Function to find the neutral pions from a neutral current interaction.
    std::vector<int> NueSelection::findNeutralPions(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth){
      
      //Pions
      std::vector<int> pions;
      
      //Loop over the map and identify the pions associated to the vertex. 
      for(std::map<int, const simb::MCParticle*>::iterator mcparticle=mcparticles.begin(); mcparticle!=mcparticles.end(); ++mcparticle){

	//Get the data
	const simb::MCParticle* mcparticle_data = mcparticle->second; 

	//Check the particle is a pion  
	if(mcparticle_data->PdgCode() != 111){continue;}

	//Check to see if the pion came from the vertex
	if(isFromNuVertex(mctruth,mcparticle_data)){
	  pions.push_back(mcparticle_data->TrackId());
	}
      }
      return pions;
    }

    //Function to get photons from neutral pions in a neutral current interaction 
    std::vector<int> NueSelection::findPhotons(std::vector<int>& pi_zeros, std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth, std::map<int,double>& visible_mcparticles){

      //Photons
      std::vector<int> photons;

      //Loop over the pions and see if they have any daughters. 
      for(int pion=0; pion<pi_zeros.size(); ++pion){
	
      	//Get the pion
      	const simb::MCParticle* mcparticle = mcparticles[pi_zeros.at(pion)]; 

      	if(mcparticle->NumberDaughters() == 0){continue;}

      	//Loop of the daughters and id photons. 
      	for(int d=0; d<mcparticle->NumberDaughters(); ++d){

      	  const simb::MCParticle* daughter = mcparticles[mcparticle->Daughter(d)];

      	  //	  int numtrajpoints = daughter->NumberTrajectoryPoints() - 2;
      	  //if(numtrajpoints<0){numtrajpoints = 0;}
      	  if(visible_mcparticles.find(daughter->TrackId()) == visible_mcparticles.end()){continue;}

	  TVector3 EndPos   = daughter->EndPosition().Vect(); 
	  for(int i=0; i<daughter->NumberTrajectoryPoints();++i){
	    if(daughter->E(i) < 0.1*daughter->E()){
	      EndPos   = daughter->Position(i).Vect();
	      break;
	    }
	  }
	  if(!containedInAV(EndPos)){continue;}


      	  //Are we a photon? (very rare we are not I think) or energy has to be above threshold
      	  if(daughter->PdgCode() == 22 && visible_mcparticles[daughter->TrackId()]*1000 > fConfig.photonVisibleEnergyThreshold){
      	    photons.push_back(daughter->TrackId());
      	  }

	}
      }

      //Look for any misc photons.
      for(std::map<int, const simb::MCParticle*>::iterator mcparticle=mcparticles.begin(); mcparticle!=mcparticles.end(); ++mcparticle){

	const simb::MCParticle* mcparticle_data = mcparticle->second;

	// For Comparision Purposes
	//	if(mcparticle_data->PdgCode() == 22 && isFromNuVertex(mctruth,mcparticle_data,1)){
	//  photons.push_back(mcparticle_data->TrackId());
	//}

	//Check the particle is a photons and is above threshold.  
	//	int numtrajpoints = mcparticle_data->NumberTrajectoryPoints() - 2;
	//	if(numtrajpoints<0){numtrajpoints = 0;}

	if(visible_mcparticles.find(mcparticle_data->TrackId()) == visible_mcparticles.end()){continue;}


	if(mcparticle_data->PdgCode() != 22 ||  visible_mcparticles[mcparticle_data->TrackId()]*1000 < fConfig.photonVisibleEnergyThreshold){continue;}

	std::cout << "Potential Photon Track ID: " << mcparticle_data->TrackId() << std::endl;


	//Check we are not a daughter product of a shower.
	if(mcparticles.find(mcparticle_data->Mother()) != mcparticles.end()){
	  if(mcparticles[mcparticle_data->Mother()]->PdgCode() == 22 || TMath::Abs(mcparticles[mcparticle_data->Mother()]->PdgCode()) == 11){continue;}
	}

	//Check we have not added the photon already.
	if(std::find(photons.begin(),photons.end(),mcparticle_data->TrackId()) != photons.end()){std::cout << "photon used" << std::endl; continue;}

	std::cout << "test" << std::endl;
	
	//See if the photon came from the vertex.
	int track_id  = mcparticle_data->TrackId();
	int mother_id = track_id;
	bool motherphoton = false;
	while(mother_id != 0){
	  if(mcparticles.find(mother_id) != mcparticles.end()){
	    if(track_id != mother_id && (mcparticles[mother_id]->PdgCode() == 22 || TMath::Abs(mcparticles[mother_id]->PdgCode()) == 11)){motherphoton=true; break;}
	    track_id = mother_id;
	    mother_id = mcparticles[mother_id]->Mother();
	  }
	  else{
	    mother_id = track_id;
	    break;
	  }
	}

	//Add the photon 
      	if(motherphoton == true){continue;}


	if(visible_mcparticles.find(mcparticle_data->TrackId()) == visible_mcparticles.end()){continue;}
	
	TVector3 EndPos   = mcparticle_data->EndPosition().Vect(); 
	for(int i=0; i<mcparticle_data->NumberTrajectoryPoints();++i){
	  if(mcparticle_data->E(i) < 0.1*mcparticle_data->E()){
	    EndPos   = mcparticle_data->Position(i).Vect();
	    break;
	  }
	}
	if(!containedInAV(EndPos)){continue;}


	if(isFromNuVertex(mctruth,mcparticles[track_id])){   
	  std::cout << "Photon Selected" << std::endl;
	  photons.push_back(mcparticle_data->TrackId());
	}
      }

      return photons;
    }
    
    //Look to see if there any photons in the fiducial volume in above the energy cut. Check there is not a complementary one in the active volume.
    bool NueSelection::PhotonEnergyCut(std::vector<int>& photons, std::map<int, const simb::MCParticle*>& mcparticles, int &photonTrackID, std::map<int,double>& visible_mcparticles){
      
      //Initialise the photons above the cut 
      std::vector<int> photons_abovecut;
      int              photons_inAV = 0;

      for(int p=0; p<photons.size(); ++p){
	const simb::MCParticle* photon = mcparticles[photons[p]];
	
	
	//	int numtrajpoints = photon->NumberTrajectoryPoints() - 2;
	//if(numtrajpoints<0){numtrajpoints = 0;}
	if(visible_mcparticles[photon->TrackId()]*1000 < fConfig.photonVisibleEnergyThreshold){continue;}

	std::cout << "Photon End Position: " << photon->EndPosition().Vect().X() << " " << photon->EndPosition().Vect().Y() << " " << photon->EndPosition().Vect().Z() << std::endl;

	TVector3 EndPos   = photon->EndPosition().Vect(); 
	for(int i=0; i<photon->NumberTrajectoryPoints();++i){
	  if(photon->E(i) < 0.1*photon->E()){
	    EndPos   = photon->Position(i).Vect();
	    break;
	  }
	}



	//Check the photons in the active volume/
	if(containedInAV(EndPos)){
	  std::cout << "not contained" << std::endl;
	  ++photons_inAV;
	}
	else{
	  std::cout << "contained" << std::endl;
	  continue;
	}
		
	//Check the photons energy.  
	//if(photon->E(photon->NumberTrajectoryPoints())*1000 < fConfig.showerVisibleEnergyThreshold){continue;}

	//Check if the photon is within the fiducal volume.
	if(containedInFV(EndPos)){
	  std::cout << "in FV" << std::endl;
	  photons_abovecut.push_back(photons[p]);
	}
      }

      //Count the number of photons above threshold from the neutrino in the FV. 
      if(photons_abovecut.size() > 1){
	std::cout << "too many photons" << std::endl;
	photonTrackID = -99999;
	return false;
      }

      //If there is no photon in the FV it still might be a nue CC event so pass -- Personally think it should have been AV 
      if(photons_abovecut.size() == 0){
	std::cout << "no photons in fv" << std::endl;
	photonTrackID = -99999;
	return true;
      }

      //There might be one in the active volume as well we can remove the event if this is the case 
      if(photons_inAV > 1){
	std::cout << "too many in AV" << std::endl;
	photonTrackID = -99999;
        return false;
      }

      //Return the boolean and the track ID 
      photonTrackID = photons_abovecut[0];
      return true;
    }

    //Function to check if the photon passes the conversion gap cut. We can only get to this stage with one photon        
    bool NueSelection::ConversionGapCut(std::map<int, const simb::MCParticle*>& mcparticles, int photonTrackID, const float& hadronic_energy, const simb::MCNeutrino& nu){
      
      //Calculate the conversion lengh of the photon.First Get the end position of the photon
      TVector3 EndPos   = mcparticles[photonTrackID]->EndPosition().Vect(); 
      for(int i=0; i<mcparticles[photonTrackID]->NumberTrajectoryPoints();++i){
	if(mcparticles[photonTrackID]->E(i) < 0.1*mcparticles[photonTrackID]->E()){
	  EndPos   = mcparticles[photonTrackID]->Position(i).Vect();
	  break;
	}
      }
	
      //Get the Vertex position
      TVector3 VtxPos = nu.Nu().Trajectory().Position(0).Vect();

      //Get the conversion length.
      double ConvLength = (EndPos - VtxPos).Mag();

      if(fConfig.Verbose){
	std::cout << "Photon Conversation length is "<< ConvLength << " cm cut is:"  <<  fConfig.photonConvLenghCut << std::endl;
	std::cout << "hadronic_energy is: " << hadronic_energy*1000 << " MeV cut is: " << fConfig.vtxEnergyCut  << std::endl;
      }

      //Check the cuts.
      if(hadronic_energy*1000 > fConfig.vtxEnergyCut && ConvLength > fConfig.photonConvLenghCut){
	return false;
      }

      return true;
    }

    //Check to see if a muon exists 
    bool NueSelection::MuLengthCut(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth){

      //Find the track index which is associated to the muon 
      int track_ind = -1;
      for(auto const& mcparticle_iter: mcparticles){
	const simb::MCParticle* mcparticle =  mcparticle_iter.second;
	if (isFromNuVertex(mctruth,mcparticle) && abs(mcparticle->PdgCode()) == 13 && mcparticle->Process() == "primary") {
	  if (track_ind == -1 || mcparticle->E() > mcparticles[track_ind]->E()) {
	    track_ind = mcparticle_iter.first;
	  }
	}
      }

      //No Muon. No problem pass the event. 
      if(track_ind == -1){return true;}

      //Get the track in question.
      const simb::MCParticle*mcparticle = mcparticles[track_ind];

      TVector3 Start = mcparticle->Position().Vect();  
      TVector3 End   = mcparticle->Position().Vect();

      //Get the End of the track in the tpc
      for(int traj=0; traj< mcparticle->NumberTrajectoryPoints(); ++traj){
	TVector3 pos = mcparticle->Position(traj).Vect();
	if(!containedInAV(pos)){End = pos; break;}
	End = pos;
      }

      //Find the length of the muon.
      double track_length = (End - Start).Mag();

      if(fConfig.Verbose){
	std::cout << "Muon Track Lengh is: " << track_length  << " cm" << std::endl;
      }

      //check if it is contained
      if(track_length > fConfig.minLengthExitingTrack){return false;}
      
      return true;
    }

    //Currently not used 
    double NueSelection::PhotonVisibleEnergy(std::map<int, const simb::MCParticle*>& mcparticles, std::vector<int>& photons){

      double energy=0;

      //Loop over the photons, check that they start in the FV then add the energy. 
      for(int p=0; p<photons.size(); ++p){
	
	const simb::MCParticle* photon = mcparticles[photons[p]];
	
	if(containedInFV(photon->EndPosition().Vect())){
	  energy += photon->E();
	}
      }

      return energy;
    }
  
    //Fill the histograms.... can't think of a better comment.
    void NueSelection::FillHistograms(std::map<std::string,TH1D*>& HistMap, std::map<std::string,std::map<int,TH1D*> >& HistMapMode, std::map<std::string,std::map<int,TH1D*> >& HistMapIntType, const simb::MCNeutrino& nu, NueSelection::NueInteraction& intInfo, bool& booldirtevent){
      
      double Energy = intInfo.GetNueEnergy();
      
      FillHistograms(HistMap,HistMapMode,HistMapIntType,nu,intInfo,Energy,booldirtevent);

    }

    void NueSelection::FillHistograms(std::map<std::string,TH1D*>& HistMap, std::map<std::string,std::map<int,TH1D*> >& HistMapMode, std::map<std::string,std::map<int,TH1D*> >& HistMapIntType, const simb::MCNeutrino& nu, NueSelection::NueInteraction& intInfo,double Energy, bool& booldirtevent){
      
      if(!fConfig.FillHistograms){return;}

      int mode = -1;
      if(nu.Mode() > -1 && nu.Mode() < 6){mode = nu.Mode();}

      if(nu.Nu().PdgCode() != intInfo.initnu && nu.Nu().PdgCode() == 11 && nu.CCNC() == 0  && dirtevent == false){
	HistMap["AllSignal"]->Fill(Energy,intInfo.weight);
      }
      else{
	HistMap["AllBackground"]->Fill(Energy,intInfo.weight);
      }

      if(dirtevent == true){

	HistMap["AllDirt"]->Fill(Energy,intInfo.weight);

	if(nu.CCNC() == simb::kCC){
	  
	  //Check to see if we are Numu background
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["DirtNuMu"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtNuMu"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)    HistMapMode["DirtNuMu"][mode]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["DirtInNuE"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtInNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)    HistMapMode["DirtInNuE"][mode]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["DirtOscNuE"]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtOscNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillModeHistograms)    HistMapMode["DirtOscNuE"][mode]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
	else{ 
	  //Neutral current histogtam fill
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["DirtNCNuMu"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtNCNuMu"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)    HistMapMode["DirtNCNuMu"][mode]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["DirtNCInNuE"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtNCInNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)    HistMapMode["DirtNCInNuE"][mode]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["DirtNCOscNuE"]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtNCOscNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillModeHistograms)     HistMapMode["DirtNCOscNuE"][mode]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
      }
      else{
	//Check if we are charged current
	if(nu.CCNC() == simb::kCC){
	  //Check to see if we are Numu background
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["NuMu"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["NuMu"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)     HistMapMode["NuMu"][mode]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["InNuE"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["InNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)     HistMapMode["InNuE"][mode]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["OscNuE"]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillIntTypeHistograms) HistMapIntType["OscNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillModeHistograms)     HistMapMode["OscNuE"][mode]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
	else{ 
	  //Neutral current histogtam fill
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["NCNuMu"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["NCNuMu"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)     HistMapMode["NCNuMu"][mode]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["NCInNuE"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["NCInNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)     HistMapMode["NCInNuE"][mode]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["NCOscNuE"]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillIntTypeHistograms) HistMapIntType["NCOscNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillModeHistograms)     HistMapMode["NCOscNuE"][mode]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
      }
      return;
    }

    //Returns the scale factor to get back to v2_08.
    double NueSelection::GENIEWeight(const simb::MCNeutrino& nu){
      
 
      double mode = nu.InteractionType();
      double weight = 1;

      if(fConfig.UseGenieHists){
	if(GenieHists.find(mode) == GenieHists.end()){return 1;}
	double E = nu.Nu().E();
	TAxis *xaxis = GenieHists[mode]->GetXaxis();
	Int_t binx = xaxis->FindBin(E);
	weight = GenieHists[mode]->GetBinContent(binx);
      }
      else{
	if(GENIEWeight_map.find(mode) == GENIEWeight_map.end()){return 1;}	
	weight = GENIEWeight_map[mode].at(0) + nu.Nu().E()*GENIEWeight_map[mode].at(1);
      }

      if(fConfig.Verbose){
	std::cout << "GENIE Weight is: " << weight << std::endl;
      }
      
      return weight;

    }


    double NueSelection::MECWeight(const simb::MCNeutrino& nu,std::string& Detector, NueSelection::NueInteraction& intInfo){
      
      double mode = nu.InteractionType();
      //      if(nu.Mode() < -1 && nu.Mode() > 6){mode =-1;}

      std::string flavour = "";

      //Get the parameters.
      if(nu.CCNC() == simb::kCC){
	
	//Check to see if we are Numu background
	if(nu.Nu().PdgCode() == 14){
	  flavour = "NuMu";
	}
	//Check if we are an intrinsic electron or oscillated 
	else if(nu.Nu().PdgCode() == intInfo.initnu){
	  flavour = "InNuE";
	}
	else if(nu.Nu().PdgCode() != intInfo.initnu){
	  flavour = "OscNuE";
	}
      }
      else{ 
	  //Neutral current histogtam fill
	  if(nu.Nu().PdgCode() == 14){
	    flavour = "NCNuMu";
	  }
	  //Check if we are an intrinsic electron or oscillated 
	  else if(nu.Nu().PdgCode() == intInfo.initnu){
	    flavour = "NCInNuE";
	  }
	  else if(nu.Nu().PdgCode() != intInfo.initnu){
	    flavour = "NCOscNuE";
	  }
      }

      if(fConfig.Verbose){
	std::cout << "Detector: " << Detector << " flavour: " << flavour << " mode: " << mode << std::endl;
      }

      if(flavour == ""){return 1;}

      if(MECWeight_map.find(Detector) == MECWeight_map.end()){return 1;}
      if(MECWeight_map[Detector].find(flavour) == MECWeight_map[Detector].end()){return 1;}
      if(MECWeight_map[Detector][flavour].find(mode) == MECWeight_map[Detector][flavour].end()){return 1;}
      

      double scale = MECWeight_map[Detector][flavour][mode].at(0);
      double mean  = MECWeight_map[Detector][flavour][mode].at(1);
      double rms   = MECWeight_map[Detector][flavour][mode].at(2);
      double x     = nu.Nu().E();
      double MECWeight = scale*TMath::Exp(-((x-mean)*(x-mean))/(2*rms*rms));

      if(fConfig.Verbose){
	std::cout << "MEC Weight is: " << MECWeight << std::endl;

	//	std::cout << "MEC Weight values are: " << MECWeight_map[Detector][flavour][mode].at(0) << + " 1: " << MECWeight_map[Detector][flavour][mode].at(1) << " = " <<  MECWeight_map[Detector][flavour][mode].at(0) + nu.Nu().E()*MECWeight_map[Detector][flavour][mode].at(1) << std::endl;
      }


      return MECWeight;

    }

    
    void NueSelection::PrintInformation(const simb::MCTruth& mctruth, NueSelection::NueInteraction& intInfo){
      
      std::cout << "#############################################" << std::endl;
      std::cout << "######### Truth Neutrino Information ########" << std::endl; 

      std::cout << mctruth.GetNeutrino() << std::endl;
      std::cout << "Positon X: " << mctruth.GetNeutrino().Nu().Vx() << " Y: " << mctruth.GetNeutrino().Nu().Vy() << " Z: " << mctruth.GetNeutrino().Nu().Vz() << " T: " << mctruth.GetNeutrino().Nu().T() << std::endl;
      
      std::cout << "######### Reco Neutrino Information ########"        << std::endl;
      std::cout << "Hadronic Energy:    "     << intInfo.hadronic_energy << std::endl;
      std::cout << "Shower Like Energy: "     << intInfo.shower_energy   << std::endl;
      std::cout << "Lepton Energy:          " << intInfo.leptonic_energy << std::endl;
      std::cout << "Other Energy: "           << intInfo.other_energy    << std::endl;
      std::cout << "Neutrino Visible Energy " << intInfo.GetNueEnergy()  << std::endl;
      std::cout << "#############################################"       << std::endl;
      
    } 

    void NueSelection::InitialiseHistograms(){



      fEventHist = new TH1I("EventHist", "EventHist",2,0,1);
      fRootHists.VisibleEnergy_FNKPDecay_Hist = new TH1D("VisibleEnergy_FNKPDecay_Hist","VisibleEnergy_FNKPDecay_Hist",fRootHists.ebins,fRootHists.emin,fRootHists.emax); 
      fRootHists.VisibleEnergy_FNKMDecay_Hist = new TH1D("VisibleEnergy_FNKMDecay_Hist","VisibleEnergy_FNKMDecay_Hist",fRootHists.ebins,fRootHists.emin,fRootHists.emax); 
      fRootHists.VisibleEnergy_MuDecays_Hist = new TH1D("VisibleEnergy_MuDecays_Hist","VisibleEnergy_MuDecays_Hist",fRootHists.ebins,fRootHists.emin,fRootHists.emax);


      //Loop over the types of interactions and set the bins
      if(fConfig.FillHistograms){
	for(auto& Type: fRootHists.HistTypes){

	  std::string  TrueNumber_String                     = Type + " TrueNumber";
	  std::string  TrueEnergy_String                     = Type + " TrueEnergy";
	  std::string  TrueEnergyAll_String                  = Type + " TrueEnergyAll";
	  std::string  CCQEEnergy_String                     = Type + " CCQEEnergy";
	  std::string  VisibleEnergy_String                  = Type + " VisibleEnergy";
	  std::string  VisibleEnergy_AVCut_String            = Type + " VisibleEnergy_AVCut";
	  std::string  VisibleEnergy_FVCut_String            = Type + " VisibleEnergy_FVCut";
	  std::string  VisibleEnergy_FVBefore_String         = Type + " VisibleEnergy_FVBefore";
	  std::string  VisibleEnergy_EnergyCut_String        = Type + " VisibleEnergy_EnergyCut";
	  std::string  VisibleEnergy_PhotonEnergyCut_String  = Type + " VisibleEnergy_PhotonEnergyCut";
	  std::string  VisibleEnergy_ConversionGapCut_String = Type + " VisibleEnergy_ConversionGapCut";
	  std::string  VisibleEnergy_MuLenghtCut_String      = Type + " VisibleEnergy_MuLenghtCut";
	  std::string  VisibleEnergy_NCCut_String            = Type + " VisibleEnergy_NCCut";
	  std::string  VisibleEnergy_Selection_String        = Type + " VisibleEnergy_Selection";
	  std::string  VisibleEnergy_PiZero_String           = Type + " VisibleEnergy_PiZero";
	  std::string  VisibleEnergy_Photon_String           = Type + " VisibleEnergy_Photon";
	  std::string  VisibleEnergy_PhotonSmall_String      = Type + " VisibleEnergy_PhotonSmall";
	  std::string  LowNCEnergy_String                    = Type + " LowNCEnergy";
	  std::string  Weights_String                        = Type + " Weights";

	  std::string  VisibleEnergy_NoShowerCut_String      = Type + " VisibleEnergy_NoShowerCut";
	  std::string  VisibleEnergy_TwoPhotonCut_String     = Type + " VisibleEnergy_TwoPhotonCut";
	  std::string   VisibleEnergy_TwoShowerCut_String     = Type + " VisibleEnergy_TwoShowerCut";
	  std::string  VisibleEnergy_PhotonFVCut_String      = Type + " VisibleEnergy_PhotonFVCut";
	  std::string  VisibleEnergy_TwoLeptons_String      = Type + " VisibleEnergy_TwoLeptons";

	  std::string  VisibleEnergy_LeptonPlusPhotonCut_String        = Type + " VisibleEnergy_LeptonPlusPhotonCut";

	  
	  std::string  TrueNumber_StringAfter                     = Type + " TrueNumberAfter";
	  std::string  TrueEnergy_StringAfter                     = Type + " TrueEnergyAfter";
	  std::string  TrueEnergyAll_StringAfter                  = Type + " TrueEnergyAllAfter";
	  std::string  CCQEEnergy_StringAfter                     = Type + " CCQEEnergyAfter";
	  std::string  VisibleEnergy_StringAfter                  = Type + " VisibleEnergyAfter";
	  std::string  VisibleEnergy_AVCut_StringAfter            = Type + " VisibleEnergy_AVCutAfter";
	  std::string  VisibleEnergy_FVCut_StringAfter            = Type + " VisibleEnergy_FVCutAfter";
	  std::string  VisibleEnergy_FVBefore_StringAfter            = Type + " VisibleEnergy_FVBeforeAfter";
	  std::string  VisibleEnergy_EnergyCut_StringAfter        = Type + " VisibleEnergy_EnergyCutAfter";
	  std::string  VisibleEnergy_PhotonEnergyCut_StringAfter  = Type + " VisibleEnergy_PhotonEnergyCutAfter";
	  std::string  VisibleEnergy_ConversionGapCut_StringAfter = Type + " VisibleEnergy_ConversionGapCutAfter";
	  std::string  VisibleEnergy_MuLenghtCut_StringAfter      = Type + " VisibleEnergy_MuLenghtCutAfter";
	  std::string  VisibleEnergy_NCCut_StringAfter            = Type + " VisibleEnergy_NCCutAfter";
	  std::string  VisibleEnergy_Selection_StringAfter        = Type + " VisibleEnergy_SelectionAfter";
	  
	  std::string  VisibleEnergy_NoShowerCut_StringAfter      = Type + " VisibleEnergy_NoShowerCutAfter";
	  std::string  VisibleEnergy_TwoPhotonCut_StringAfter     = Type + " VisibleEnergy_TwoPhotonCutAfter";
	  std::string   VisibleEnergy_TwoShowerCut_StringAfter     = Type + " VisibleEnergy_TwoShowerCutAfter";
	  std::string  VisibleEnergy_PhotonFVCut_StringAfter      = Type + " VisibleEnergy_PhotonFVCutAfter";
	  std::string  VisibleEnergy_TwoLeptons_StringAfter      = Type + " VisibleEnergy_TwoLeptonsAfter";
	  std::string VisibleEnergy_LeptonPlusPhotonCut_StringAfter = Type + " VisibleEnergy_LeptonPlusPhotonCut_After";



	  std::string ProtonE_String = Type + " ProtonE";
	  std::string PionE_String = Type + " PionE";
	  std::string KaonE_String = Type + " KaonE";
	  std::string ProtonN_String = Type + " ProtonN";
	  std::string PionN_String = Type + " PionN";
	  std::string KaonN_String = Type + " KaonN";
	  std::string HadronE_String = Type + " HadronE";
	  std::string PhotonCon_String = Type + " PhotonCon";


	  const char* TrueNumber_Name                     = TrueNumber_String.c_str();
	  const char* TrueEnergy_Name                     = TrueEnergy_String.c_str();
	  const char* TrueEnergyAll_Name                  = TrueEnergyAll_String.c_str();
	  const char* CCQEEnergy_Name                     = CCQEEnergy_String.c_str();
	  const char* VisibleEnergy_Name                  = VisibleEnergy_String.c_str();
	  const char* VisibleEnergy_AVCut_Name            = VisibleEnergy_AVCut_String.c_str();
	  const char* VisibleEnergy_FVCut_Name            = VisibleEnergy_FVCut_String.c_str();
	  const char* VisibleEnergy_FVBefore_Name            = VisibleEnergy_FVBefore_String.c_str();
	  const char* VisibleEnergy_EnergyCut_Name        = VisibleEnergy_EnergyCut_String.c_str();
	  const char* VisibleEnergy_PhotonEnergyCut_Name  = VisibleEnergy_PhotonEnergyCut_String.c_str();
	  const char* VisibleEnergy_ConversionGapCut_Name = VisibleEnergy_ConversionGapCut_String.c_str();
	  const char* VisibleEnergy_MuLenghtCut_Name      = VisibleEnergy_MuLenghtCut_String.c_str();
	  const char* VisibleEnergy_NCCut_Name            = VisibleEnergy_NCCut_String.c_str();
	  const char* VisibleEnergy_Selection_Name        = VisibleEnergy_Selection_String.c_str();
	  const char* VisibleEnergy_NoShowerCut_Name      = VisibleEnergy_NoShowerCut_String.c_str();
	  const char* VisibleEnergy_TwoPhotonCut_Name     = VisibleEnergy_TwoPhotonCut_String.c_str();
	  const char* VisibleEnergy_TwoShowerCut_Name     = VisibleEnergy_TwoShowerCut_String.c_str();
	  const char* VisibleEnergy_PhotonFVCut_Name      = VisibleEnergy_PhotonFVCut_String.c_str();
	  const char* VisibleEnergy_TwoLeptons_Name      = VisibleEnergy_TwoLeptons_String.c_str();

	  const char* TrueNumber_NameAfter                     = TrueNumber_StringAfter.c_str();
	  const char* TrueEnergy_NameAfter                     = TrueEnergy_StringAfter.c_str();
	  const char* TrueEnergyAll_NameAfter                  = TrueEnergyAll_StringAfter.c_str();
	  const char* CCQEEnergy_NameAfter                     = CCQEEnergy_StringAfter.c_str();
	  const char* VisibleEnergy_NameAfter                  = VisibleEnergy_StringAfter.c_str();
	  const char* VisibleEnergy_AVCut_NameAfter            = VisibleEnergy_AVCut_StringAfter.c_str();
	  const char* VisibleEnergy_FVCut_NameAfter            = VisibleEnergy_FVCut_StringAfter.c_str();
	  const char* VisibleEnergy_FVBefore_NameAfter            = VisibleEnergy_FVBefore_StringAfter.c_str();
	  const char* VisibleEnergy_EnergyCut_NameAfter        = VisibleEnergy_EnergyCut_StringAfter.c_str();
	  const char* VisibleEnergy_PhotonEnergyCut_NameAfter  = VisibleEnergy_PhotonEnergyCut_StringAfter.c_str();
	  const char* VisibleEnergy_ConversionGapCut_NameAfter = VisibleEnergy_ConversionGapCut_StringAfter.c_str();
	  const char* VisibleEnergy_MuLenghtCut_NameAfter      = VisibleEnergy_MuLenghtCut_StringAfter.c_str();
	  const char* VisibleEnergy_NCCut_NameAfter            = VisibleEnergy_NCCut_StringAfter.c_str();
	  const char* VisibleEnergy_Selection_NameAfter        = VisibleEnergy_Selection_StringAfter.c_str();
	  const char* VisibleEnergy_NoShowerCut_NameAfter      = VisibleEnergy_NoShowerCut_StringAfter.c_str();
	  const char* VisibleEnergy_TwoPhotonCut_NameAfter     = VisibleEnergy_TwoPhotonCut_StringAfter.c_str();
	  const char* VisibleEnergy_TwoShowerCut_NameAfter     = VisibleEnergy_TwoShowerCut_StringAfter.c_str();
	  const char* VisibleEnergy_PhotonFVCut_NameAfter      = VisibleEnergy_PhotonFVCut_StringAfter.c_str();
	  const char* VisibleEnergy_TwoLeptons_NameAfter      = VisibleEnergy_TwoLeptons_StringAfter.c_str();


	  const char* VisibleEnergy_LeptonPlusPhotonCut_NameAfter = VisibleEnergy_LeptonPlusPhotonCut_StringAfter.c_str();

	  const char* VisibleEnergy_PiZero_Name           = VisibleEnergy_PiZero_String.c_str();
	  const char* VisibleEnergy_Photon_Name           = VisibleEnergy_Photon_String.c_str();
	  const char* VisibleEnergy_PhotonSmall_Name      = VisibleEnergy_PhotonSmall_String.c_str();
	  const char* LowNCEnergy_Name                    = LowNCEnergy_String.c_str();

	  const char* Weights_Name        = Weights_String.c_str();

	  const char* VisibleEnergy_LeptonPlusPhotonCut_Name        = VisibleEnergy_LeptonPlusPhotonCut_String.c_str();
	  
	  const char* ProtonE_Name = ProtonE_String.c_str();
	  const char* PionE_Name = PionE_String.c_str();
	  const char* KaonE_Name = KaonE_String.c_str();
	  const char* ProtonN_Name = ProtonN_String.c_str();
	  const char* PionN_Name = PionN_String.c_str();
	  const char* KaonN_Name = KaonN_String.c_str();
	  const char* HadronE_Name = HadronE_String.c_str();
	  const char* PhotonCon_Name = PhotonCon_String.c_str();


	  fRootHists.TrueNumber_Hist[Type] = new TH1D(TrueNumber_Name,TrueNumber_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.TrueEnergy_Hist[Type] = new TH1D(TrueEnergy_Name,TrueEnergy_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.TrueEnergyAll_Hist[Type] = new TH1D(TrueEnergyAll_Name,TrueEnergyAll_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.CCQEEnergy_Hist[Type] = new TH1D(CCQEEnergy_Name,CCQEEnergy_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_Hist[Type] = new TH1D(VisibleEnergy_Name,VisibleEnergy_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_AVCut_Hist[Type] = new TH1D(VisibleEnergy_AVCut_Name,VisibleEnergy_AVCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_FVCut_Hist[Type] = new TH1D(VisibleEnergy_FVCut_Name,VisibleEnergy_FVCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_FVBefore_Hist[Type] = new TH1D(VisibleEnergy_FVBefore_Name,VisibleEnergy_FVBefore_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergy_EnergyCut_Hist[Type] = new TH1D(VisibleEnergy_EnergyCut_Name,VisibleEnergy_EnergyCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_PhotonEnergyCut_Hist[Type] = new TH1D(VisibleEnergy_PhotonEnergyCut_Name,VisibleEnergy_PhotonEnergyCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_ConversionGapCut_Hist[Type] = new TH1D(VisibleEnergy_ConversionGapCut_Name,VisibleEnergy_ConversionGapCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_MuLenghtCut_Hist[Type] = new TH1D(VisibleEnergy_MuLenghtCut_Name,VisibleEnergy_MuLenghtCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_NCCut_Hist[Type] = new TH1D(VisibleEnergy_NCCut_Name,VisibleEnergy_NCCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_Selection_Hist[Type] = new TH1D(VisibleEnergy_Selection_Name,VisibleEnergy_Selection_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_LeptonPlusPhotonCut_Hist[Type] = new TH1D(VisibleEnergy_LeptonPlusPhotonCut_Name,VisibleEnergy_LeptonPlusPhotonCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergy_NoShowerCut_Hist[Type] = new TH1D(VisibleEnergy_NoShowerCut_Name,VisibleEnergy_NoShowerCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_TwoPhotonCut_Hist[Type] = new TH1D(VisibleEnergy_TwoPhotonCut_Name,VisibleEnergy_TwoPhotonCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_TwoShowerCut_Hist[Type] = new TH1D(VisibleEnergy_TwoShowerCut_Name,VisibleEnergy_TwoShowerCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_PhotonFVCut_Hist[Type] = new TH1D(VisibleEnergy_PhotonFVCut_Name,VisibleEnergy_PhotonFVCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_TwoLeptons_Hist[Type] = new TH1D(VisibleEnergy_TwoLeptons_Name,VisibleEnergy_TwoLeptons_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);




	  fRootHists.VisibleEnergyAfter_Hist[Type] = new TH1D(VisibleEnergy_NameAfter,VisibleEnergy_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_AVCutAfter_Hist[Type] = new TH1D(VisibleEnergy_AVCut_NameAfter,VisibleEnergy_AVCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_FVCutAfter_Hist[Type] = new TH1D(VisibleEnergy_FVCut_NameAfter,VisibleEnergy_FVCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_FVBeforeAfter_Hist[Type] = new TH1D(VisibleEnergy_FVBefore_NameAfter,VisibleEnergy_FVBefore_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergy_EnergyCutAfter_Hist[Type] = new TH1D(VisibleEnergy_EnergyCut_NameAfter,VisibleEnergy_EnergyCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_PhotonEnergyCutAfter_Hist[Type] = new TH1D(VisibleEnergy_PhotonEnergyCut_NameAfter,VisibleEnergy_PhotonEnergyCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_ConversionGapCutAfter_Hist[Type] = new TH1D(VisibleEnergy_ConversionGapCut_NameAfter,VisibleEnergy_ConversionGapCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_MuLenghtCutAfter_Hist[Type] = new TH1D(VisibleEnergy_MuLenghtCut_NameAfter,VisibleEnergy_MuLenghtCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_NCCutAfter_Hist[Type] = new TH1D(VisibleEnergy_NCCut_NameAfter,VisibleEnergy_NCCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_LeptonPlusPhotonCutAfter_Hist[Type] = new TH1D(VisibleEnergy_LeptonPlusPhotonCut_NameAfter,VisibleEnergy_LeptonPlusPhotonCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergy_NoShowerCutAfter_Hist[Type] = new TH1D(VisibleEnergy_NoShowerCut_NameAfter,VisibleEnergy_NoShowerCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_TwoPhotonCutAfter_Hist[Type] = new TH1D(VisibleEnergy_TwoPhotonCut_NameAfter,VisibleEnergy_TwoPhotonCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_TwoShowerCutAfter_Hist[Type] = new TH1D(VisibleEnergy_TwoShowerCut_NameAfter,VisibleEnergy_TwoShowerCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_PhotonFVCutAfter_Hist[Type] = new TH1D(VisibleEnergy_PhotonFVCut_NameAfter,VisibleEnergy_PhotonFVCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_TwoLeptonsAfter_Hist[Type] = new TH1D(VisibleEnergy_TwoLeptons_NameAfter,VisibleEnergy_TwoLeptons_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);



	  fRootHists.VisibleEnergy_PiZero_Hist[Type] = new TH1D(VisibleEnergy_PiZero_Name,VisibleEnergy_PiZero_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_Photon_Hist[Type] = new TH1D(VisibleEnergy_Photon_Name,VisibleEnergy_Photon_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergy_PhotonSmall_Hist[Type] = new TH1D(VisibleEnergy_PhotonSmall_Name,VisibleEnergy_PhotonSmall_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.LowNCEnergy_Hist[Type] = new TH1D(LowNCEnergy_Name,LowNCEnergy_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	  fRootHists.Weights_Hist[Type] = new TH1D(Weights_Name,Weights_Name,100,0,2);
	  fRootHists.ProtonE_Hist[Type] = new TH1D(ProtonE_Name,ProtonE_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);	
	  fRootHists.PionE_Hist[Type] = new TH1D(PionE_Name,PionE_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.KaonE_Hist[Type] = new TH1D(KaonE_Name,KaonE_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.ProtonN_Hist[Type] = new TH1D(ProtonN_Name,ProtonN_Name,20,0,20);	
	  fRootHists.PionN_Hist[Type] = new TH1D(PionN_Name,PionN_Name,20,0,20);
	  fRootHists.KaonN_Hist[Type] = new TH1D(KaonN_Name,KaonN_Name,20,0,20);
	  fRootHists.HadronE_Hist[Type] = new TH1D(HadronE_Name,HadronE_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.PhotonCon_Hist[Type] = new TH1D(PhotonCon_Name,PhotonCon_Name,100,0,40);



	  if(fConfig.FillModeHistograms){

	    for(int mode=-1; mode<6; ++mode){
	    
	      std::string  TrueNumber_StringMode                     = TrueNumber_String + " " +  std::to_string(mode);
	      std::string  TrueEnergy_StringMode                     = TrueEnergy_String + " " +  std::to_string(mode);
	      std::string  TrueEnergyAll_StringMode                  = TrueEnergyAll_String + " " +  std::to_string(mode);  
	      std::string  CCQEEnergy_StringMode                     = CCQEEnergy_String + " " +  std::to_string(mode);                    
	      std::string  VisibleEnergy_StringMode                  = VisibleEnergy_String + " " +  std::to_string(mode);                  
	      std::string  VisibleEnergy_AVCut_StringMode            = VisibleEnergy_AVCut_String + " " +  std::to_string(mode);           
	      std::string  VisibleEnergy_FVCut_StringMode            = VisibleEnergy_FVCut_String + " " +  std::to_string(mode);            
	      std::string  VisibleEnergy_FVBefore_StringMode            = VisibleEnergy_FVBefore_String + " " +  std::to_string(mode);            

	      std::string  VisibleEnergy_EnergyCut_StringMode        = VisibleEnergy_EnergyCut_String + " " +  std::to_string(mode);        
	      std::string  VisibleEnergy_PhotonEnergyCut_StringMode  = VisibleEnergy_PhotonEnergyCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_ConversionGapCut_StringMode = VisibleEnergy_ConversionGapCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_MuLenghtCut_StringMode      = VisibleEnergy_MuLenghtCut_String + " " +  std::to_string(mode);      
	      std::string  VisibleEnergy_NCCut_StringMode            = VisibleEnergy_NCCut_String + " " +  std::to_string(mode);            
	      std::string  VisibleEnergy_Selection_StringMode        = VisibleEnergy_Selection_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_PiZero_StringMode           = VisibleEnergy_PiZero_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_Photon_StringMode           = VisibleEnergy_Photon_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_PhotonSmall_StringMode           = VisibleEnergy_PhotonSmall_String + " " +  std::to_string(mode);
	      std::string  LowNCEnergy_StringMode                     = LowNCEnergy_String + " " +  std::to_string(mode);

	      std::string  VisibleEnergy_NoShowerCut_StringMode      = VisibleEnergy_NoShowerCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_TwoPhotonCut_StringMode     = VisibleEnergy_TwoPhotonCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_TwoShowerCut_StringMode     = VisibleEnergy_TwoShowerCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_PhotonFVCut_StringMode      = VisibleEnergy_PhotonFVCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_TwoLeptons_StringMode      = VisibleEnergy_PhotonFVCut_String + " " +  std::to_string(mode);
	  

	      std::string  VisibleEnergy_LeptonPlusPhotonCut_StringMode = VisibleEnergy_LeptonPlusPhotonCut_String + " " +  std::to_string(mode);
	      std::string  Weights_StringMode                           = Weights_String  + " " +  std::to_string(mode);
	      std::string  ProtonE_StringMode                     = ProtonE_String + " " +  std::to_string(mode);
	      std::string  PionE_StringMode                       = PionE_String + " " +  std::to_string(mode);
	      std::string  KaonE_StringMode                       = KaonE_String + " " +  std::to_string(mode);
	      std::string  ProtonN_StringMode                     = ProtonN_String + " " +  std::to_string(mode);
	      std::string  PionN_StringMode                       = PionN_String + " " +  std::to_string(mode);
	      std::string  KaonN_StringMode                       = KaonN_String + " " +  std::to_string(mode);
	      std::string  HadronE_StringMode                       = HadronE_String + " " +  std::to_string(mode);
	      std::string  PhotonCon_StringMode                       = PhotonCon_String + " " +  std::to_string(mode);


	      std::string  VisibleEnergy_StringAfterMode                  = VisibleEnergy_StringAfter + " " +  std::to_string(mode);                  
	      std::string  VisibleEnergy_AVCut_StringAfterMode            = VisibleEnergy_AVCut_StringAfter + " " +  std::to_string(mode);           
	      std::string  VisibleEnergy_FVCut_StringAfterMode            = VisibleEnergy_FVCut_StringAfter + " " +  std::to_string(mode);            
	      std::string  VisibleEnergy_FVBefore_StringAfterMode            = VisibleEnergy_FVBefore_StringAfter + " " +  std::to_string(mode);            

	      std::string  VisibleEnergy_EnergyCut_StringAfterMode        = VisibleEnergy_EnergyCut_StringAfter + " " +  std::to_string(mode);        
	      std::string  VisibleEnergy_PhotonEnergyCut_StringAfterMode  = VisibleEnergy_PhotonEnergyCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_ConversionGapCut_StringAfterMode = VisibleEnergy_ConversionGapCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_MuLenghtCut_StringAfterMode      = VisibleEnergy_MuLenghtCut_StringAfter + " " +  std::to_string(mode);      
	      std::string  VisibleEnergy_NCCut_StringAfterMode            = VisibleEnergy_NCCut_StringAfter + " " +  std::to_string(mode);            
	      std::string  VisibleEnergy_Selection_StringAfterMode        = VisibleEnergy_Selection_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_NoShowerCut_StringAfterMode      = VisibleEnergy_NoShowerCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_TwoPhotonCut_StringAfterMode     = VisibleEnergy_TwoPhotonCut_StringAfter + " " +  std::to_string(mode);
	      std::string   VisibleEnergy_TwoShowerCut_StringAfterMode     = VisibleEnergy_TwoShowerCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_PhotonFVCut_StringAfterMode      = VisibleEnergy_PhotonFVCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergy_TwoLeptons_StringAfterMode      = VisibleEnergy_PhotonFVCut_StringAfter + " " +  std::to_string(mode);

	      std::string  VisibleEnergy_LeptonPlusPhotonCut_StringAfterMode = VisibleEnergy_LeptonPlusPhotonCut_StringAfter + " " +  std::to_string(mode);


	      const char* TrueNumber_NameMode                     = TrueNumber_StringMode.c_str();
	      const char* TrueEnergy_NameMode                     = TrueEnergy_StringMode.c_str();
	      const char* TrueEnergyAll_NameMode                  = TrueEnergyAll_StringMode.c_str();
	      const char* CCQEEnergy_NameMode                     = CCQEEnergy_StringMode.c_str();
	      const char* VisibleEnergy_NameMode                  = VisibleEnergy_StringMode.c_str();
	      const char* VisibleEnergy_AVCut_NameMode            = VisibleEnergy_AVCut_StringMode.c_str();
	      const char* VisibleEnergy_FVCut_NameMode            = VisibleEnergy_FVCut_StringMode.c_str();
	      const char* VisibleEnergy_FVBefore_NameMode            = VisibleEnergy_FVBefore_StringMode.c_str();

	      const char* VisibleEnergy_EnergyCut_NameMode        = VisibleEnergy_EnergyCut_StringMode.c_str();
	      const char* VisibleEnergy_PhotonEnergyCut_NameMode  = VisibleEnergy_PhotonEnergyCut_StringMode.c_str();
	      const char* VisibleEnergy_ConversionGapCut_NameMode = VisibleEnergy_ConversionGapCut_StringMode.c_str();
	      const char* VisibleEnergy_MuLenghtCut_NameMode      = VisibleEnergy_MuLenghtCut_StringMode.c_str();
	      const char* VisibleEnergy_NCCut_NameMode            = VisibleEnergy_NCCut_StringMode.c_str();
	      const char* VisibleEnergy_Selection_NameMode        = VisibleEnergy_Selection_StringMode.c_str();
	      const char* VisibleEnergy_PiZero_NameMode           = VisibleEnergy_PiZero_StringMode.c_str();
	      const char* VisibleEnergy_Photon_NameMode           = VisibleEnergy_Photon_StringMode.c_str();
	      const char* VisibleEnergy_PhotonSmall_NameMode           = VisibleEnergy_PhotonSmall_StringMode.c_str();
	      const char* LowNCEnergy_NameMode                     = LowNCEnergy_StringMode.c_str();

	      const char*  VisibleEnergy_NoShowerCut_NameMode      = VisibleEnergy_NoShowerCut_StringMode.c_str();
	      const char*  VisibleEnergy_TwoPhotonCut_NameMode     = VisibleEnergy_TwoPhotonCut_StringMode.c_str();
	      const char*  VisibleEnergy_TwoShowerCut_NameMode     = VisibleEnergy_TwoShowerCut_StringMode.c_str();
	      const char*  VisibleEnergy_PhotonFVCut_NameMode      = VisibleEnergy_TwoShowerCut_StringMode.c_str();
	      const char*  VisibleEnergy_TwoLeptons_NameMode      = VisibleEnergy_TwoLeptons_StringMode.c_str();


	      const char* VisibleEnergy_LeptonPlusPhotonCut_NameMode = VisibleEnergy_LeptonPlusPhotonCut_StringMode.c_str();
	      const char* Weights_NameMode                           = Weights_StringMode.c_str();
	      const char* ProtonE_NameMode                     = ProtonE_StringMode.c_str();
	      const char* PionE_NameMode                       = PionE_StringMode.c_str();
	      const char* KaonE_NameMode                       = KaonE_StringMode.c_str();
	      const char* ProtonN_NameMode                     = ProtonN_StringMode.c_str();
	      const char* PionN_NameMode                       = PionN_StringMode.c_str();
	      const char* KaonN_NameMode                       = KaonN_StringMode.c_str();
	      const char* HadronE_NameMode                       = HadronE_StringMode.c_str();
	      const char* PhotonCon_NameMode                       = PhotonCon_StringMode.c_str();


	      const char* VisibleEnergy_NameAfterMode                  = VisibleEnergy_StringAfterMode.c_str();
	      const char* VisibleEnergy_AVCut_NameAfterMode            = VisibleEnergy_AVCut_StringAfterMode.c_str();
	      const char* VisibleEnergy_FVCut_NameAfterMode            = VisibleEnergy_FVCut_StringAfterMode.c_str();
	      const char* VisibleEnergy_FVBefore_NameAfterMode            = VisibleEnergy_FVBefore_StringAfterMode.c_str();

	      const char* VisibleEnergy_EnergyCut_NameAfterMode        = VisibleEnergy_EnergyCut_StringAfterMode.c_str();
	      const char* VisibleEnergy_PhotonEnergyCut_NameAfterMode  = VisibleEnergy_PhotonEnergyCut_StringAfterMode.c_str();
	      const char* VisibleEnergy_ConversionGapCut_NameAfterMode = VisibleEnergy_ConversionGapCut_StringAfterMode.c_str();
	      const char* VisibleEnergy_MuLenghtCut_NameAfterMode      = VisibleEnergy_MuLenghtCut_StringAfterMode.c_str();
	      const char* VisibleEnergy_NCCut_NameAfterMode            = VisibleEnergy_NCCut_StringAfterMode.c_str();
	      const char* VisibleEnergy_Selection_NameAfterMode        = VisibleEnergy_Selection_StringAfterMode.c_str();
	      const char* VisibleEnergy_LeptonPlusPhotonCut_NameAfterMode = VisibleEnergy_LeptonPlusPhotonCut_StringAfterMode.c_str();

	      const char*  VisibleEnergy_NoShowerCut_NameAfterMode      = VisibleEnergy_NoShowerCut_StringAfterMode.c_str();
	      const char*  VisibleEnergy_TwoPhotonCut_NameAfterMode     = VisibleEnergy_TwoPhotonCut_StringAfterMode.c_str();
	      const char*  VisibleEnergy_TwoShowerCut_NameAfterMode     = VisibleEnergy_TwoShowerCut_StringAfterMode.c_str();
	      const char*  VisibleEnergy_PhotonFVCut_NameAfterMode      = VisibleEnergy_TwoShowerCut_StringAfterMode.c_str();
	      const char*  VisibleEnergy_TwoLeptons_NameAfterMode      = VisibleEnergy_TwoLeptons_StringAfterMode.c_str();


	      
	      fRootHists.TrueNumber_HistMode[Type][mode] = new TH1D(TrueNumber_NameMode,TrueNumber_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.TrueEnergy_HistMode[Type][mode] = new TH1D(TrueEnergy_NameMode,TrueEnergy_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.TrueEnergyAll_HistMode[Type][mode] = new TH1D(TrueEnergyAll_NameMode,TrueEnergyAll_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.CCQEEnergy_HistMode[Type][mode] = new TH1D(CCQEEnergy_NameMode,CCQEEnergy_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_HistMode[Type][mode] = new TH1D(VisibleEnergy_NameMode,VisibleEnergy_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_AVCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_AVCut_NameMode,VisibleEnergy_AVCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_FVCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_FVCut_NameMode,VisibleEnergy_FVCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_FVBefore_HistMode[Type][mode] = new TH1D(VisibleEnergy_FVBefore_NameMode,VisibleEnergy_FVBefore_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      fRootHists.VisibleEnergy_EnergyCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_EnergyCut_NameMode,VisibleEnergy_EnergyCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonEnergyCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_PhotonEnergyCut_NameMode,VisibleEnergy_PhotonEnergyCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_ConversionGapCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_ConversionGapCut_NameMode,VisibleEnergy_ConversionGapCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_MuLenghtCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_MuLenghtCut_NameMode,VisibleEnergy_MuLenghtCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_NCCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_NCCut_NameMode,VisibleEnergy_NCCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_Selection_HistMode[Type][mode] = new TH1D(VisibleEnergy_Selection_NameMode,VisibleEnergy_Selection_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PiZero_HistMode[Type][mode] = new TH1D(VisibleEnergy_PiZero_NameMode,VisibleEnergy_PiZero_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_Photon_HistMode[Type][mode] = new TH1D(VisibleEnergy_Photon_NameMode,VisibleEnergy_Photon_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonSmall_HistMode[Type][mode] = new TH1D(VisibleEnergy_PhotonSmall_NameMode,VisibleEnergy_PhotonSmall_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	    
	      fRootHists.Weights_HistMode[Type][mode]  = new TH1D(Weights_NameMode,Weights_NameMode,100,0,2); 
	 
	      fRootHists.VisibleEnergy_NoShowerCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_NoShowerCut_NameMode,VisibleEnergy_NoShowerCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoPhotonCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_TwoPhotonCut_NameMode,VisibleEnergy_TwoPhotonCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoShowerCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_TwoShowerCut_NameMode,VisibleEnergy_TwoShowerCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonFVCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_PhotonFVCut_NameMode,VisibleEnergy_PhotonFVCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoLeptons_HistMode[Type][mode] = new TH1D(VisibleEnergy_TwoLeptons_NameMode,VisibleEnergy_TwoLeptons_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      
 
	      fRootHists.VisibleEnergy_LeptonPlusPhotonCut_HistMode[Type][mode] = new TH1D(VisibleEnergy_LeptonPlusPhotonCut_NameMode,VisibleEnergy_LeptonPlusPhotonCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.LowNCEnergy_HistMode[Type][mode] = new TH1D(LowNCEnergy_NameMode,LowNCEnergy_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.ProtonE_HistMode[Type][mode] = new TH1D(ProtonE_NameMode,ProtonE_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.PionE_HistMode[Type][mode] = new TH1D(PionE_NameMode,PionE_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.KaonE_HistMode[Type][mode] = new TH1D(KaonE_NameMode,KaonE_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.HadronE_HistMode[Type][mode] = new TH1D(HadronE_NameMode,HadronE_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.PhotonCon_HistMode[Type][mode] = new TH1D(PhotonCon_NameMode,PhotonCon_NameMode,100,0,40);

	      fRootHists.ProtonN_HistMode[Type][mode] = new TH1D(ProtonN_NameMode,ProtonN_NameMode,20,0,20);
	      fRootHists.PionN_HistMode[Type][mode] = new TH1D(PionN_NameMode,PionN_NameMode,20,0,20);
	      fRootHists.KaonN_HistMode[Type][mode] = new TH1D(KaonN_NameMode,KaonN_NameMode,20,0,20);

	       fRootHists.VisibleEnergyAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_NameAfterMode,VisibleEnergy_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_AVCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_AVCut_NameAfterMode,VisibleEnergy_AVCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_FVCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_FVCut_NameAfterMode,VisibleEnergy_FVCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_FVBeforeAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_FVBefore_NameAfterMode,VisibleEnergy_FVBefore_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      fRootHists.VisibleEnergy_EnergyCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_EnergyCut_NameAfterMode,VisibleEnergy_EnergyCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonEnergyCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_PhotonEnergyCut_NameAfterMode,VisibleEnergy_PhotonEnergyCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_ConversionGapCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_ConversionGapCut_NameAfterMode,VisibleEnergy_ConversionGapCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_MuLenghtCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_MuLenghtCut_NameAfterMode,VisibleEnergy_MuLenghtCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_NCCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_NCCut_NameAfterMode,VisibleEnergy_NCCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_SelectionAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_Selection_NameAfterMode,VisibleEnergy_Selection_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_LeptonPlusPhotonCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_LeptonPlusPhotonCut_NameAfterMode,VisibleEnergy_LeptonPlusPhotonCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      fRootHists.VisibleEnergy_NoShowerCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_NoShowerCut_NameAfterMode,VisibleEnergy_NoShowerCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoPhotonCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_TwoPhotonCut_NameAfterMode,VisibleEnergy_TwoPhotonCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoShowerCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_TwoShowerCut_NameAfterMode,VisibleEnergy_TwoShowerCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonFVCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_PhotonFVCut_NameAfterMode,VisibleEnergy_PhotonFVCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoLeptonsAfter_HistMode[Type][mode] = new TH1D(VisibleEnergy_TwoLeptons_NameAfterMode,VisibleEnergy_TwoLeptons_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);



	    }
	  }

	  if(fConfig.FillIntTypeHistograms){
	    for(int inttype=1000; inttype<1101; ++inttype){

	      std::string  TrueNumber_StringIntType                     = TrueNumber_String + " " +  std::to_string(inttype);
	      std::string  TrueEnergy_StringIntType                     = TrueEnergy_String + " " +  std::to_string(inttype);
	      std::string  TrueEnergyAll_StringIntType                  = TrueEnergyAll_String + " " +  std::to_string(inttype);  
	      std::string  CCQEEnergy_StringIntType                     = CCQEEnergy_String + " " +  std::to_string(inttype);                    
	      std::string  VisibleEnergy_StringIntType                  = VisibleEnergy_String + " " +  std::to_string(inttype);                  
	      std::string  VisibleEnergy_AVCut_StringIntType            = VisibleEnergy_AVCut_String + " " +  std::to_string(inttype);           
	      std::string  VisibleEnergy_FVCut_StringIntType            = VisibleEnergy_FVCut_String + " " +  std::to_string(inttype);            
	      std::string  VisibleEnergy_FVBefore_StringIntType            = VisibleEnergy_FVBefore_String + " " +  std::to_string(inttype);            

	      std::string  VisibleEnergy_EnergyCut_StringIntType        = VisibleEnergy_EnergyCut_String + " " +  std::to_string(inttype);        
	      std::string  VisibleEnergy_PhotonEnergyCut_StringIntType  = VisibleEnergy_PhotonEnergyCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_ConversionGapCut_StringIntType = VisibleEnergy_ConversionGapCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_MuLenghtCut_StringIntType      = VisibleEnergy_MuLenghtCut_String + " " +  std::to_string(inttype);      
	      std::string  VisibleEnergy_NCCut_StringIntType            = VisibleEnergy_NCCut_String + " " +  std::to_string(inttype);            
	      std::string  VisibleEnergy_Selection_StringIntType        = VisibleEnergy_Selection_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_PiZero_StringIntType           = VisibleEnergy_PiZero_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_Photon_StringIntType           = VisibleEnergy_Photon_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_PhotonSmall_StringIntType           = VisibleEnergy_PhotonSmall_String + " " +  std::to_string(inttype);
	      std::string  LowNCEnergy_StringIntType                     = LowNCEnergy_String + " " +  std::to_string(inttype);

	      std::string  VisibleEnergy_LeptonPlusPhotonCut_StringIntType = VisibleEnergy_LeptonPlusPhotonCut_String + " " +  std::to_string(inttype);
	      std::string  Weights_StringIntType                           = Weights_String  + " " +  std::to_string(inttype);
	      std::string  ProtonE_StringIntType                     = ProtonE_String + " " +  std::to_string(inttype);
	      std::string  PionE_StringIntType                       = PionE_String   + " " +  std::to_string(inttype);
	      std::string  KaonE_StringIntType                       = KaonE_String   + " " +  std::to_string(inttype);
	      std::string  ProtonN_StringIntType                     = ProtonN_String + " " +  std::to_string(inttype);
	      std::string  PionN_StringIntType                       = PionN_String   + " " +  std::to_string(inttype);
	      std::string  KaonN_StringIntType                       = KaonN_String   + " " +  std::to_string(inttype);
	      std::string  HadronE_StringIntType                       = HadronE_String   + " " +  std::to_string(inttype);
	      std::string  PhotonCon_StringIntType                       = PhotonCon_String   + " " +  std::to_string(inttype);
	      
	      std::string  VisibleEnergy_NoShowerCut_StringIntType      = VisibleEnergy_NoShowerCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_TwoPhotonCut_StringIntType     = VisibleEnergy_TwoPhotonCut_String + " " +  std::to_string(inttype);
	      std::string   VisibleEnergy_TwoShowerCut_StringIntType     = VisibleEnergy_TwoShowerCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_PhotonFVCut_StringIntType      = VisibleEnergy_PhotonFVCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_TwoLeptons_StringIntType      = VisibleEnergy_PhotonFVCut_String + " " +  std::to_string(inttype);



	      const char* TrueNumber_NameIntType                     = TrueNumber_StringIntType.c_str();
	      const char* TrueEnergy_NameIntType                     = TrueEnergy_StringIntType.c_str();
	      const char* TrueEnergyAll_NameIntType                  = TrueEnergyAll_StringIntType.c_str();
	      const char* CCQEEnergy_NameIntType                     = CCQEEnergy_StringIntType.c_str();
	      const char* VisibleEnergy_NameIntType                  = VisibleEnergy_StringIntType.c_str();
	      const char* VisibleEnergy_AVCut_NameIntType            = VisibleEnergy_AVCut_StringIntType.c_str();
	      const char* VisibleEnergy_FVCut_NameIntType            = VisibleEnergy_FVCut_StringIntType.c_str();
	      const char* VisibleEnergy_FVBefore_NameIntType            = VisibleEnergy_FVBefore_StringIntType.c_str();

	      const char* VisibleEnergy_EnergyCut_NameIntType        = VisibleEnergy_EnergyCut_StringIntType.c_str();
	      const char* VisibleEnergy_PhotonEnergyCut_NameIntType  = VisibleEnergy_PhotonEnergyCut_StringIntType.c_str();
	      const char* VisibleEnergy_ConversionGapCut_NameIntType = VisibleEnergy_ConversionGapCut_StringIntType.c_str();
	      const char* VisibleEnergy_MuLenghtCut_NameIntType      = VisibleEnergy_MuLenghtCut_StringIntType.c_str();
	      const char* VisibleEnergy_NCCut_NameIntType            = VisibleEnergy_NCCut_StringIntType.c_str();
	      const char* VisibleEnergy_Selection_NameIntType        = VisibleEnergy_Selection_StringIntType.c_str();
	      const char* VisibleEnergy_PiZero_NameIntType           = VisibleEnergy_PiZero_StringIntType.c_str();
	      const char* VisibleEnergy_Photon_NameIntType           = VisibleEnergy_Photon_StringIntType.c_str();
	      const char* VisibleEnergy_PhotonSmall_NameIntType           = VisibleEnergy_PhotonSmall_StringIntType.c_str();
	      const char* LowNCEnergy_NameIntType                     = LowNCEnergy_StringIntType.c_str();

	      const char* VisibleEnergy_LeptonPlusPhotonCut_NameIntType = VisibleEnergy_LeptonPlusPhotonCut_StringIntType.c_str();
	      const char* Weights_NameIntType                           = Weights_StringIntType.c_str();
	      const char* ProtonE_NameIntType                     = ProtonE_StringIntType.c_str();
	      const char* PionE_NameIntType                       = PionE_StringIntType.c_str();
	      const char* KaonE_NameIntType                       = KaonE_StringIntType.c_str();
	      const char* ProtonN_NameIntType                     = ProtonN_StringIntType.c_str();
	      const char* PionN_NameIntType                       = PionN_StringIntType.c_str();
	      const char* KaonN_NameIntType                       = KaonN_StringIntType.c_str();
	      const char* HadronE_NameIntType                       = HadronE_StringIntType.c_str();
	      const char* PhotonCon_NameIntType                       = PhotonCon_StringIntType.c_str();

	      const char*  VisibleEnergy_NoShowerCut_NameIntType      = VisibleEnergy_NoShowerCut_StringIntType.c_str();
	      const char*  VisibleEnergy_TwoPhotonCut_NameIntType     = VisibleEnergy_TwoPhotonCut_StringIntType.c_str();
	      const char*  VisibleEnergy_TwoShowerCut_NameIntType     = VisibleEnergy_TwoShowerCut_StringIntType.c_str();
	      const char*  VisibleEnergy_PhotonFVCut_NameIntType      = VisibleEnergy_TwoShowerCut_StringIntType.c_str();
	      const char*  VisibleEnergy_TwoLeptons_NameIntType      = VisibleEnergy_TwoLeptons_StringIntType.c_str();


	      fRootHists.TrueNumber_HistIntType[Type][inttype] = new TH1D(TrueNumber_NameIntType,TrueNumber_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.TrueEnergy_HistIntType[Type][inttype] = new TH1D(TrueEnergy_NameIntType,TrueEnergy_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.TrueEnergyAll_HistIntType[Type][inttype] = new TH1D(TrueEnergyAll_NameIntType,TrueEnergyAll_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.CCQEEnergy_HistIntType[Type][inttype] = new TH1D(CCQEEnergy_NameIntType,CCQEEnergy_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NameIntType,VisibleEnergy_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_AVCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_AVCut_NameIntType,VisibleEnergy_AVCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_FVCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_FVCut_NameIntType,VisibleEnergy_FVCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_FVBefore_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_FVBefore_NameIntType,VisibleEnergy_FVBefore_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      fRootHists.VisibleEnergy_EnergyCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_EnergyCut_NameIntType,VisibleEnergy_EnergyCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonEnergyCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_PhotonEnergyCut_NameIntType,VisibleEnergy_PhotonEnergyCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_ConversionGapCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ConversionGapCut_NameIntType,VisibleEnergy_ConversionGapCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_MuLenghtCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MuLenghtCut_NameIntType,VisibleEnergy_MuLenghtCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_NCCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NCCut_NameIntType,VisibleEnergy_NCCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_Selection_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_Selection_NameIntType,VisibleEnergy_Selection_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PiZero_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_PiZero_NameIntType,VisibleEnergy_PiZero_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_Photon_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_Photon_NameIntType,VisibleEnergy_Photon_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonSmall_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_PhotonSmall_NameIntType,VisibleEnergy_PhotonSmall_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	    
	      fRootHists.Weights_HistIntType[Type][inttype]  = new TH1D(Weights_NameIntType,Weights_NameIntType,100,0,2); 
	  
	  
	      fRootHists.VisibleEnergy_LeptonPlusPhotonCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_LeptonPlusPhotonCut_NameIntType,VisibleEnergy_LeptonPlusPhotonCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.LowNCEnergy_HistIntType[Type][inttype] = new TH1D(LowNCEnergy_NameIntType,LowNCEnergy_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.ProtonE_HistIntType[Type][inttype] = new TH1D(ProtonE_NameIntType,ProtonE_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.PionE_HistIntType[Type][inttype] = new TH1D(PionE_NameIntType,PionE_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.KaonE_HistIntType[Type][inttype] = new TH1D(KaonE_NameIntType,KaonE_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.HadronE_HistIntType[Type][inttype] = new TH1D(HadronE_NameIntType,HadronE_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.PhotonCon_HistIntType[Type][inttype] = new TH1D(PhotonCon_NameIntType,PhotonCon_NameIntType,100,0,40);

	      fRootHists.ProtonN_HistIntType[Type][inttype] = new TH1D(ProtonN_NameIntType,ProtonN_NameIntType,20,0,20);
	      fRootHists.PionN_HistIntType[Type][inttype] = new TH1D(PionN_NameIntType,PionN_NameIntType,20,0,20);
	      fRootHists.KaonN_HistIntType[Type][inttype] = new TH1D(KaonN_NameIntType,KaonN_NameIntType,20,0,20);
	      
	      fRootHists.VisibleEnergy_NoShowerCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NoShowerCut_NameIntType,VisibleEnergy_NoShowerCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoPhotonCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_TwoPhotonCut_NameIntType,VisibleEnergy_TwoPhotonCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoShowerCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_TwoShowerCut_NameIntType,VisibleEnergy_TwoShowerCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonFVCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_PhotonFVCut_NameIntType,VisibleEnergy_PhotonFVCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoLeptons_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_TwoLeptons_NameIntType,VisibleEnergy_TwoLeptons_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      std::string  VisibleEnergy_StringAfterIntType                  = VisibleEnergy_StringAfter + " " +  std::to_string(inttype);                  
	      std::string  VisibleEnergy_AVCut_StringAfterIntType            = VisibleEnergy_AVCut_StringAfter + " " +  std::to_string(inttype);           
	      std::string  VisibleEnergy_FVCut_StringAfterIntType            = VisibleEnergy_FVCut_StringAfter + " " +  std::to_string(inttype);            
	      std::string  VisibleEnergy_FVBefore_StringAfterIntType            = VisibleEnergy_FVBefore_StringAfter + " " +  std::to_string(inttype);            

	      std::string  VisibleEnergy_EnergyCut_StringAfterIntType        = VisibleEnergy_EnergyCut_StringAfter + " " +  std::to_string(inttype);        
	      std::string  VisibleEnergy_PhotonEnergyCut_StringAfterIntType  = VisibleEnergy_PhotonEnergyCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_ConversionGapCut_StringAfterIntType = VisibleEnergy_ConversionGapCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_MuLenghtCut_StringAfterIntType      = VisibleEnergy_MuLenghtCut_StringAfter + " " +  std::to_string(inttype);      
	      std::string  VisibleEnergy_NCCut_StringAfterIntType            = VisibleEnergy_NCCut_StringAfter + " " +  std::to_string(inttype);            
	      std::string  VisibleEnergy_Selection_StringAfterIntType        = VisibleEnergy_Selection_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_LeptonPlusPhotonCut_StringAfterIntType = VisibleEnergy_LeptonPlusPhotonCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_NoShowerCut_StringAfterIntType      = VisibleEnergy_NoShowerCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_TwoPhotonCut_StringAfterIntType     = VisibleEnergy_TwoPhotonCut_StringAfter + " " +  std::to_string(inttype);
	      std::string   VisibleEnergy_TwoShowerCut_StringAfterIntType     = VisibleEnergy_TwoShowerCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_PhotonFVCut_StringAfterIntType      = VisibleEnergy_PhotonFVCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergy_TwoLeptons_StringAfterIntType      = VisibleEnergy_PhotonFVCut_StringAfter + " " +  std::to_string(inttype);




	      const char* VisibleEnergy_NameAfterIntType                  = VisibleEnergy_StringAfterIntType.c_str();
	      const char* VisibleEnergy_AVCut_NameAfterIntType            = VisibleEnergy_AVCut_StringAfterIntType.c_str();
	      const char* VisibleEnergy_FVCut_NameAfterIntType            = VisibleEnergy_FVCut_StringAfterIntType.c_str();
	      const char* VisibleEnergy_FVBefore_NameAfterIntType            = VisibleEnergy_FVBefore_StringAfterIntType.c_str();

	      const char* VisibleEnergy_EnergyCut_NameAfterIntType        = VisibleEnergy_EnergyCut_StringAfterIntType.c_str();
	      const char* VisibleEnergy_PhotonEnergyCut_NameAfterIntType  = VisibleEnergy_PhotonEnergyCut_StringAfterIntType.c_str();
	      const char* VisibleEnergy_ConversionGapCut_NameAfterIntType = VisibleEnergy_ConversionGapCut_StringAfterIntType.c_str();
	      const char* VisibleEnergy_MuLenghtCut_NameAfterIntType      = VisibleEnergy_MuLenghtCut_StringAfterIntType.c_str();
	      const char* VisibleEnergy_NCCut_NameAfterIntType            = VisibleEnergy_NCCut_StringAfterIntType.c_str();
	      const char* VisibleEnergy_Selection_NameAfterIntType        = VisibleEnergy_Selection_StringAfterIntType.c_str();
	      const char* VisibleEnergy_LeptonPlusPhotonCut_NameAfterIntType = VisibleEnergy_LeptonPlusPhotonCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergy_NoShowerCut_NameAfterIntType      = VisibleEnergy_NoShowerCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergy_TwoPhotonCut_NameAfterIntType     = VisibleEnergy_TwoPhotonCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergy_TwoShowerCut_NameAfterIntType     = VisibleEnergy_TwoShowerCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergy_PhotonFVCut_NameAfterIntType      = VisibleEnergy_TwoShowerCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergy_TwoLeptons_NameAfterIntType      = VisibleEnergy_TwoLeptons_StringAfterIntType.c_str();


	      fRootHists.VisibleEnergyAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NameAfterIntType,VisibleEnergy_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_AVCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_AVCut_NameAfterIntType,VisibleEnergy_AVCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_FVCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_FVCut_NameAfterIntType,VisibleEnergy_FVCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_FVBeforeAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_FVBefore_NameAfterIntType,VisibleEnergy_FVBefore_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      fRootHists.VisibleEnergy_EnergyCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_EnergyCut_NameAfterIntType,VisibleEnergy_EnergyCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonEnergyCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_PhotonEnergyCut_NameAfterIntType,VisibleEnergy_PhotonEnergyCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_ConversionGapCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ConversionGapCut_NameAfterIntType,VisibleEnergy_ConversionGapCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_MuLenghtCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MuLenghtCut_NameAfterIntType,VisibleEnergy_MuLenghtCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_NCCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NCCut_NameAfterIntType,VisibleEnergy_NCCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_SelectionAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_Selection_NameAfterIntType,VisibleEnergy_Selection_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  
	      fRootHists.VisibleEnergy_LeptonPlusPhotonCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_LeptonPlusPhotonCut_NameAfterIntType,VisibleEnergy_LeptonPlusPhotonCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      fRootHists.VisibleEnergy_NoShowerCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NoShowerCut_NameAfterIntType,VisibleEnergy_NoShowerCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoPhotonCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_TwoPhotonCut_NameAfterIntType,VisibleEnergy_TwoPhotonCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoShowerCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_TwoShowerCut_NameAfterIntType,VisibleEnergy_TwoShowerCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_PhotonFVCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_PhotonFVCut_NameAfterIntType,VisibleEnergy_PhotonFVCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergy_TwoLeptonsAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_TwoLeptons_NameAfterIntType,VisibleEnergy_TwoLeptons_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	    }
	  }
	}
      }


       //Loop over the types of interactions and set the bins
      if(fConfig.FillHistograms){
	for(auto& Type: fRootHists.HistTypes){

	  std::string  VisibleEnergyLep_String                  = Type + " VisibleEnergyLep";
	  std::string  VisibleEnergyLep_AVCut_String            = Type + " VisibleEnergyLep_AVCut";
	  std::string  VisibleEnergyLep_FVCut_String            = Type + " VisibleEnergyLep_FVCut";
	  std::string  VisibleEnergyLep_FVBefore_String         = Type + " VisibleEnergyLep_FVBefore";
	  std::string  VisibleEnergyLep_EnergyCut_String        = Type + " VisibleEnergyLep_EnergyCut";
	  std::string  VisibleEnergyLep_PhotonEnergyCut_String  = Type + " VisibleEnergyLep_PhotonEnergyCut";
	  std::string  VisibleEnergyLep_ConversionGapCut_String = Type + " VisibleEnergyLep_ConversionGapCut";
	  std::string  VisibleEnergyLep_MuLenghtCut_String      = Type + " VisibleEnergyLep_MuLenghtCut";
	  std::string  VisibleEnergyLep_NCCut_String            = Type + " VisibleEnergyLep_NCCut";
	  std::string  VisibleEnergyLep_Selection_String        = Type + " VisibleEnergyLep_Selection";
	  std::string  VisibleEnergyLep_PiZero_String           = Type + " VisibleEnergyLep_PiZero";
	  std::string  VisibleEnergyLep_Photon_String           = Type + " VisibleEnergyLep_Photon";
	  std::string  VisibleEnergyLep_PhotonSmall_String      = Type + " VisibleEnergyLep_PhotonSmall";
	  std::string  VisibleEnergyLep_NoShowerCut_String      = Type + " VisibleEnergyLep_NoShowerCut";
	  std::string  VisibleEnergyLep_TwoPhotonCut_String     = Type + " VisibleEnergyLep_TwoPhotonCut";
	  std::string   VisibleEnergyLep_TwoShowerCut_String     = Type + " VisibleEnergyLep_TwoShowerCut";
	  std::string  VisibleEnergyLep_PhotonFVCut_String      = Type + " VisibleEnergyLep_PhotonFVCut";
	  std::string  VisibleEnergyLep_TwoLeptons_String      = Type + " VisibleEnergyLep_TwoLeptons";
	  std::string  VisibleEnergyLep_LeptonPlusPhotonCut_String        = Type + " VisibleEnergyLep_LeptonPlusPhotonCut";

	  std::string  VisibleEnergyLep_StringAfter                  = Type + " VisibleEnergyLepAfter";
	  std::string  VisibleEnergyLep_AVCut_StringAfter            = Type + " VisibleEnergyLep_AVCutAfter";
	  std::string  VisibleEnergyLep_FVCut_StringAfter            = Type + " VisibleEnergyLep_FVCutAfter";
	  std::string  VisibleEnergyLep_FVBefore_StringAfter            = Type + " VisibleEnergyLep_FVBeforeAfter";
	  std::string  VisibleEnergyLep_EnergyCut_StringAfter        = Type + " VisibleEnergyLep_EnergyCutAfter";
	  std::string  VisibleEnergyLep_PhotonEnergyCut_StringAfter  = Type + " VisibleEnergyLep_PhotonEnergyCutAfter";
	  std::string  VisibleEnergyLep_ConversionGapCut_StringAfter = Type + " VisibleEnergyLep_ConversionGapCutAfter";
	  std::string  VisibleEnergyLep_MuLenghtCut_StringAfter      = Type + " VisibleEnergyLep_MuLenghtCutAfter";
	  std::string  VisibleEnergyLep_NCCut_StringAfter            = Type + " VisibleEnergyLep_NCCutAfter";
	  std::string  VisibleEnergyLep_Selection_StringAfter        = Type + " VisibleEnergyLep_SelectionAfter";
	  std::string  VisibleEnergyLep_NoShowerCut_StringAfter      = Type + " VisibleEnergyLep_NoShowerCutAfter";
	  std::string  VisibleEnergyLep_TwoPhotonCut_StringAfter     = Type + " VisibleEnergyLep_TwoPhotonCutAfter";
	  std::string   VisibleEnergyLep_TwoShowerCut_StringAfter     = Type + " VisibleEnergyLep_TwoShowerCutAfter";
	  std::string  VisibleEnergyLep_PhotonFVCut_StringAfter      = Type + " VisibleEnergyLep_PhotonFVCutAfter";
	  std::string  VisibleEnergyLep_TwoLeptons_StringAfter      = Type + " VisibleEnergyLep_TwoLeptonsAfter";
	  std::string VisibleEnergyLep_LeptonPlusPhotonCut_StringAfter = Type + "VisibleEnergyLep_LeptonPlusPhotonCut_StringAfter"; 

	  const char* VisibleEnergyLep_Name                  = VisibleEnergyLep_String.c_str();
	  const char* VisibleEnergyLep_AVCut_Name            = VisibleEnergyLep_AVCut_String.c_str();
	  const char* VisibleEnergyLep_FVCut_Name            = VisibleEnergyLep_FVCut_String.c_str();
	  const char* VisibleEnergyLep_FVBefore_Name            = VisibleEnergyLep_FVBefore_String.c_str();
	  const char* VisibleEnergyLep_EnergyCut_Name        = VisibleEnergyLep_EnergyCut_String.c_str();
	  const char* VisibleEnergyLep_PhotonEnergyCut_Name  = VisibleEnergyLep_PhotonEnergyCut_String.c_str();
	  const char* VisibleEnergyLep_ConversionGapCut_Name = VisibleEnergyLep_ConversionGapCut_String.c_str();
	  const char* VisibleEnergyLep_MuLenghtCut_Name      = VisibleEnergyLep_MuLenghtCut_String.c_str();
	  const char* VisibleEnergyLep_NCCut_Name            = VisibleEnergyLep_NCCut_String.c_str();
	  const char* VisibleEnergyLep_Selection_Name        = VisibleEnergyLep_Selection_String.c_str();
	  const char* VisibleEnergyLep_NoShowerCut_Name      = VisibleEnergyLep_NoShowerCut_String.c_str();
	  const char* VisibleEnergyLep_TwoPhotonCut_Name     = VisibleEnergyLep_TwoPhotonCut_String.c_str();
	  const char* VisibleEnergyLep_TwoShowerCut_Name     = VisibleEnergyLep_TwoShowerCut_String.c_str();
	  const char* VisibleEnergyLep_PhotonFVCut_Name      = VisibleEnergyLep_PhotonFVCut_String.c_str();
	  const char* VisibleEnergyLep_TwoLeptons_Name      = VisibleEnergyLep_TwoLeptons_String.c_str();

	  const char* VisibleEnergyLep_NameAfter                  = VisibleEnergyLep_StringAfter.c_str();
	  const char* VisibleEnergyLep_AVCut_NameAfter            = VisibleEnergyLep_AVCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_FVCut_NameAfter            = VisibleEnergyLep_FVCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_FVBefore_NameAfter            = VisibleEnergyLep_FVBefore_StringAfter.c_str();
	  const char* VisibleEnergyLep_EnergyCut_NameAfter        = VisibleEnergyLep_EnergyCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_PhotonEnergyCut_NameAfter  = VisibleEnergyLep_PhotonEnergyCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_ConversionGapCut_NameAfter = VisibleEnergyLep_ConversionGapCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_MuLenghtCut_NameAfter      = VisibleEnergyLep_MuLenghtCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_NCCut_NameAfter            = VisibleEnergyLep_NCCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_Selection_NameAfter        = VisibleEnergyLep_Selection_StringAfter.c_str();
	  const char* VisibleEnergyLep_NoShowerCut_NameAfter      = VisibleEnergyLep_NoShowerCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_TwoPhotonCut_NameAfter     = VisibleEnergyLep_TwoPhotonCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_TwoShowerCut_NameAfter     = VisibleEnergyLep_TwoShowerCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_PhotonFVCut_NameAfter      = VisibleEnergyLep_PhotonFVCut_StringAfter.c_str();
	  const char* VisibleEnergyLep_TwoLeptons_NameAfter      = VisibleEnergyLep_TwoLeptons_StringAfter.c_str();


	  const char* VisibleEnergyLep_PiZero_Name           = VisibleEnergyLep_PiZero_String.c_str();
	  const char* VisibleEnergyLep_Photon_Name           = VisibleEnergyLep_Photon_String.c_str();
	  const char* VisibleEnergyLep_LeptonPlusPhotonCut_Name        = VisibleEnergyLep_LeptonPlusPhotonCut_String.c_str();

	  const char* VisibleEnergyLep_LeptonPlusPhotonCut_NameAfter        = VisibleEnergyLep_LeptonPlusPhotonCut_StringAfter.c_str();

	  fRootHists.VisibleEnergyLep_Hist[Type] = new TH1D(VisibleEnergyLep_Name,VisibleEnergyLep_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_AVCut_Hist[Type] = new TH1D(VisibleEnergyLep_AVCut_Name,VisibleEnergyLep_AVCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_FVCut_Hist[Type] = new TH1D(VisibleEnergyLep_FVCut_Name,VisibleEnergyLep_FVCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_FVBefore_Hist[Type] = new TH1D(VisibleEnergyLep_FVBefore_Name,VisibleEnergyLep_FVBefore_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergyLep_EnergyCut_Hist[Type] = new TH1D(VisibleEnergyLep_EnergyCut_Name,VisibleEnergyLep_EnergyCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_PhotonEnergyCut_Hist[Type] = new TH1D(VisibleEnergyLep_PhotonEnergyCut_Name,VisibleEnergyLep_PhotonEnergyCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_ConversionGapCut_Hist[Type] = new TH1D(VisibleEnergyLep_ConversionGapCut_Name,VisibleEnergyLep_ConversionGapCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_MuLenghtCut_Hist[Type] = new TH1D(VisibleEnergyLep_MuLenghtCut_Name,VisibleEnergyLep_MuLenghtCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_NCCut_Hist[Type] = new TH1D(VisibleEnergyLep_NCCut_Name,VisibleEnergyLep_NCCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_Selection_Hist[Type] = new TH1D(VisibleEnergyLep_Selection_Name,VisibleEnergyLep_Selection_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_LeptonPlusPhotonCut_Hist[Type] = new TH1D(VisibleEnergyLep_LeptonPlusPhotonCut_Name,VisibleEnergyLep_LeptonPlusPhotonCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergyLep_NoShowerCut_Hist[Type] = new TH1D(VisibleEnergyLep_NoShowerCut_Name,VisibleEnergyLep_NoShowerCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_TwoPhotonCut_Hist[Type] = new TH1D(VisibleEnergyLep_TwoPhotonCut_Name,VisibleEnergyLep_TwoPhotonCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_TwoShowerCut_Hist[Type] = new TH1D(VisibleEnergyLep_TwoShowerCut_Name,VisibleEnergyLep_TwoShowerCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_PhotonFVCut_Hist[Type] = new TH1D(VisibleEnergyLep_PhotonFVCut_Name,VisibleEnergyLep_PhotonFVCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_TwoLeptons_Hist[Type] = new TH1D(VisibleEnergyLep_TwoLeptons_Name,VisibleEnergyLep_TwoLeptons_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergyLepAfter_Hist[Type] = new TH1D(VisibleEnergyLep_NameAfter,VisibleEnergyLep_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_AVCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_AVCut_NameAfter,VisibleEnergyLep_AVCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_FVCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_FVCut_NameAfter,VisibleEnergyLep_FVCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_FVBeforeAfter_Hist[Type] = new TH1D(VisibleEnergyLep_FVBefore_NameAfter,VisibleEnergyLep_FVBefore_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergyLep_EnergyCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_EnergyCut_NameAfter,VisibleEnergyLep_EnergyCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_PhotonEnergyCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_PhotonEnergyCut_NameAfter,VisibleEnergyLep_PhotonEnergyCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_ConversionGapCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_ConversionGapCut_NameAfter,VisibleEnergyLep_ConversionGapCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_MuLenghtCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_MuLenghtCut_NameAfter,VisibleEnergyLep_MuLenghtCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_NCCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_NCCut_NameAfter,VisibleEnergyLep_NCCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_SelectionAfter_Hist[Type] = new TH1D(VisibleEnergyLep_Selection_NameAfter,VisibleEnergyLep_Selection_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_LeptonPlusPhotonCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_LeptonPlusPhotonCut_NameAfter,VisibleEnergyLep_LeptonPlusPhotonCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  fRootHists.VisibleEnergyLep_NoShowerCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_NoShowerCut_NameAfter,VisibleEnergyLep_NoShowerCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_TwoPhotonCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_TwoPhotonCut_NameAfter,VisibleEnergyLep_TwoPhotonCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_TwoShowerCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_TwoShowerCut_NameAfter,VisibleEnergyLep_TwoShowerCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_PhotonFVCutAfter_Hist[Type] = new TH1D(VisibleEnergyLep_PhotonFVCut_NameAfter,VisibleEnergyLep_PhotonFVCut_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  fRootHists.VisibleEnergyLep_TwoLeptonsAfter_Hist[Type] = new TH1D(VisibleEnergyLep_TwoLeptons_NameAfter,VisibleEnergyLep_TwoLeptons_NameAfter,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	  if(fConfig.FillModeHistograms){

	    for(int mode=-1; mode<6; ++mode){
	      
	      std::string  VisibleEnergyLep_StringMode                  = VisibleEnergyLep_String + " " +  std::to_string(mode);                  
	      std::string  VisibleEnergyLep_AVCut_StringMode            = VisibleEnergyLep_AVCut_String + " " +  std::to_string(mode);           
	      std::string  VisibleEnergyLep_FVCut_StringMode            = VisibleEnergyLep_FVCut_String + " " +  std::to_string(mode);            
	      std::string  VisibleEnergyLep_FVBefore_StringMode            = VisibleEnergyLep_FVBefore_String + " " +  std::to_string(mode);            

	      std::string  VisibleEnergyLep_EnergyCut_StringMode        = VisibleEnergyLep_EnergyCut_String + " " +  std::to_string(mode);        
	      std::string  VisibleEnergyLep_PhotonEnergyCut_StringMode  = VisibleEnergyLep_PhotonEnergyCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_ConversionGapCut_StringMode = VisibleEnergyLep_ConversionGapCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_MuLenghtCut_StringMode      = VisibleEnergyLep_MuLenghtCut_String + " " +  std::to_string(mode);      
	      std::string  VisibleEnergyLep_NCCut_StringMode            = VisibleEnergyLep_NCCut_String + " " +  std::to_string(mode);            
	      std::string  VisibleEnergyLep_Selection_StringMode        = VisibleEnergyLep_Selection_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_NoShowerCut_StringMode      = VisibleEnergyLep_NoShowerCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_TwoPhotonCut_StringMode     = VisibleEnergyLep_TwoPhotonCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_TwoShowerCut_StringMode     = VisibleEnergyLep_TwoShowerCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_PhotonFVCut_StringMode      = VisibleEnergyLep_PhotonFVCut_String + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_TwoLeptons_StringMode      = VisibleEnergyLep_PhotonFVCut_String + " " +  std::to_string(mode);
	  

	      std::string  VisibleEnergyLep_LeptonPlusPhotonCut_StringMode = VisibleEnergyLep_LeptonPlusPhotonCut_String + " " +  std::to_string(mode);
	      
	      std::string  VisibleEnergyLep_StringAfterMode                  = VisibleEnergyLep_StringAfter + " " +  std::to_string(mode);                  
	      std::string  VisibleEnergyLep_AVCut_StringAfterMode            = VisibleEnergyLep_AVCut_StringAfter + " " +  std::to_string(mode);           
	      std::string  VisibleEnergyLep_FVCut_StringAfterMode            = VisibleEnergyLep_FVCut_StringAfter + " " +  std::to_string(mode);            
	      std::string  VisibleEnergyLep_FVBefore_StringAfterMode            = VisibleEnergyLep_FVBefore_StringAfter + " " +  std::to_string(mode);            

	      std::string  VisibleEnergyLep_EnergyCut_StringAfterMode        = VisibleEnergyLep_EnergyCut_StringAfter + " " +  std::to_string(mode);        
	      std::string  VisibleEnergyLep_PhotonEnergyCut_StringAfterMode  = VisibleEnergyLep_PhotonEnergyCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_ConversionGapCut_StringAfterMode = VisibleEnergyLep_ConversionGapCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_MuLenghtCut_StringAfterMode      = VisibleEnergyLep_MuLenghtCut_StringAfter + " " +  std::to_string(mode);      
	      std::string  VisibleEnergyLep_NCCut_StringAfterMode            = VisibleEnergyLep_NCCut_StringAfter + " " +  std::to_string(mode);            
	      std::string  VisibleEnergyLep_Selection_StringAfterMode        = VisibleEnergyLep_Selection_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_NoShowerCut_StringAfterMode      = VisibleEnergyLep_NoShowerCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_TwoPhotonCut_StringAfterMode     = VisibleEnergyLep_TwoPhotonCut_StringAfter + " " +  std::to_string(mode);
	      std::string   VisibleEnergyLep_TwoShowerCut_StringAfterMode     = VisibleEnergyLep_TwoShowerCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_PhotonFVCut_StringAfterMode      = VisibleEnergyLep_PhotonFVCut_StringAfter + " " +  std::to_string(mode);
	      std::string  VisibleEnergyLep_TwoLeptons_StringAfterMode      = VisibleEnergyLep_PhotonFVCut_StringAfter + " " +  std::to_string(mode);

	      std::string  VisibleEnergyLep_LeptonPlusPhotonCut_StringAfterMode = VisibleEnergyLep_LeptonPlusPhotonCut_StringAfter + " " +  std::to_string(mode);


	      const char* VisibleEnergyLep_NameMode                  = VisibleEnergyLep_StringMode.c_str();
	      const char* VisibleEnergyLep_AVCut_NameMode            = VisibleEnergyLep_AVCut_StringMode.c_str();
	      const char* VisibleEnergyLep_FVCut_NameMode            = VisibleEnergyLep_FVCut_StringMode.c_str();
	      const char* VisibleEnergyLep_FVBefore_NameMode            = VisibleEnergyLep_FVBefore_StringMode.c_str();

	      const char* VisibleEnergyLep_EnergyCut_NameMode        = VisibleEnergyLep_EnergyCut_StringMode.c_str();
	      const char* VisibleEnergyLep_PhotonEnergyCut_NameMode  = VisibleEnergyLep_PhotonEnergyCut_StringMode.c_str();
	      const char* VisibleEnergyLep_ConversionGapCut_NameMode = VisibleEnergyLep_ConversionGapCut_StringMode.c_str();
	      const char* VisibleEnergyLep_MuLenghtCut_NameMode      = VisibleEnergyLep_MuLenghtCut_StringMode.c_str();
	      const char* VisibleEnergyLep_NCCut_NameMode            = VisibleEnergyLep_NCCut_StringMode.c_str();
	      const char* VisibleEnergyLep_Selection_NameMode        = VisibleEnergyLep_Selection_StringMode.c_str();
	      const char*  VisibleEnergyLep_NoShowerCut_NameMode      = VisibleEnergyLep_NoShowerCut_StringMode.c_str();
	      const char*  VisibleEnergyLep_TwoPhotonCut_NameMode     = VisibleEnergyLep_TwoPhotonCut_StringMode.c_str();
	      const char*  VisibleEnergyLep_TwoShowerCut_NameMode     = VisibleEnergyLep_TwoShowerCut_StringMode.c_str();
	      const char*  VisibleEnergyLep_PhotonFVCut_NameMode      = VisibleEnergyLep_TwoShowerCut_StringMode.c_str();
	      const char*  VisibleEnergyLep_TwoLeptons_NameMode      = VisibleEnergyLep_TwoLeptons_StringMode.c_str();


	      const char* VisibleEnergyLep_LeptonPlusPhotonCut_NameMode = VisibleEnergyLep_LeptonPlusPhotonCut_StringMode.c_str();
	      const char* VisibleEnergyLep_NameAfterMode                  = VisibleEnergyLep_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_AVCut_NameAfterMode            = VisibleEnergyLep_AVCut_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_FVCut_NameAfterMode            = VisibleEnergyLep_FVCut_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_FVBefore_NameAfterMode            = VisibleEnergyLep_FVBefore_StringAfterMode.c_str();

	      const char* VisibleEnergyLep_EnergyCut_NameAfterMode        = VisibleEnergyLep_EnergyCut_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_PhotonEnergyCut_NameAfterMode  = VisibleEnergyLep_PhotonEnergyCut_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_ConversionGapCut_NameAfterMode = VisibleEnergyLep_ConversionGapCut_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_MuLenghtCut_NameAfterMode      = VisibleEnergyLep_MuLenghtCut_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_NCCut_NameAfterMode            = VisibleEnergyLep_NCCut_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_Selection_NameAfterMode        = VisibleEnergyLep_Selection_StringAfterMode.c_str();
	      const char* VisibleEnergyLep_LeptonPlusPhotonCut_NameAfterMode = VisibleEnergyLep_LeptonPlusPhotonCut_StringAfterMode.c_str();

	      const char*  VisibleEnergyLep_NoShowerCut_NameAfterMode      = VisibleEnergyLep_NoShowerCut_StringAfterMode.c_str();
	      const char*  VisibleEnergyLep_TwoPhotonCut_NameAfterMode     = VisibleEnergyLep_TwoPhotonCut_StringAfterMode.c_str();
	      const char*  VisibleEnergyLep_TwoShowerCut_NameAfterMode     = VisibleEnergyLep_TwoShowerCut_StringAfterMode.c_str();
	      const char*  VisibleEnergyLep_PhotonFVCut_NameAfterMode      = VisibleEnergyLep_TwoShowerCut_StringAfterMode.c_str();
	      const char*  VisibleEnergyLep_TwoLeptons_NameAfterMode      = VisibleEnergyLep_TwoLeptons_StringAfterMode.c_str();


	      

	      fRootHists.VisibleEnergyLep_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_NameMode,VisibleEnergyLep_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_AVCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_AVCut_NameMode,VisibleEnergyLep_AVCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_FVCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_FVCut_NameMode,VisibleEnergyLep_FVCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_FVBefore_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_FVBefore_NameMode,VisibleEnergyLep_FVBefore_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      fRootHists.VisibleEnergyLep_EnergyCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_EnergyCut_NameMode,VisibleEnergyLep_EnergyCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_PhotonEnergyCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_PhotonEnergyCut_NameMode,VisibleEnergyLep_PhotonEnergyCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_ConversionGapCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_ConversionGapCut_NameMode,VisibleEnergyLep_ConversionGapCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_MuLenghtCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_MuLenghtCut_NameMode,VisibleEnergyLep_MuLenghtCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_NCCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_NCCut_NameMode,VisibleEnergyLep_NCCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_Selection_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_Selection_NameMode,VisibleEnergyLep_Selection_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	 
	      fRootHists.VisibleEnergyLep_NoShowerCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_NoShowerCut_NameMode,VisibleEnergyLep_NoShowerCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoPhotonCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_TwoPhotonCut_NameMode,VisibleEnergyLep_TwoPhotonCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoShowerCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_TwoShowerCut_NameMode,VisibleEnergyLep_TwoShowerCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_PhotonFVCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_PhotonFVCut_NameMode,VisibleEnergyLep_PhotonFVCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoLeptons_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_TwoLeptons_NameMode,VisibleEnergyLep_TwoLeptons_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      
 
	      fRootHists.VisibleEnergyLep_LeptonPlusPhotonCut_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_LeptonPlusPhotonCut_NameMode,VisibleEnergyLep_LeptonPlusPhotonCut_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLepAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_NameAfterMode,VisibleEnergyLep_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_AVCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_AVCut_NameAfterMode,VisibleEnergyLep_AVCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_FVCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_FVCut_NameAfterMode,VisibleEnergyLep_FVCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_FVBeforeAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_FVBefore_NameAfterMode,VisibleEnergyLep_FVBefore_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      fRootHists.VisibleEnergyLep_EnergyCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_EnergyCut_NameAfterMode,VisibleEnergyLep_EnergyCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_PhotonEnergyCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_PhotonEnergyCut_NameAfterMode,VisibleEnergyLep_PhotonEnergyCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_ConversionGapCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_ConversionGapCut_NameAfterMode,VisibleEnergyLep_ConversionGapCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_MuLenghtCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_MuLenghtCut_NameAfterMode,VisibleEnergyLep_MuLenghtCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_NCCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_NCCut_NameAfterMode,VisibleEnergyLep_NCCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_SelectionAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_Selection_NameAfterMode,VisibleEnergyLep_Selection_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_LeptonPlusPhotonCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_LeptonPlusPhotonCut_NameAfterMode,VisibleEnergyLep_LeptonPlusPhotonCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      fRootHists.VisibleEnergyLep_NoShowerCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_NoShowerCut_NameAfterMode,VisibleEnergyLep_NoShowerCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoPhotonCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_TwoPhotonCut_NameAfterMode,VisibleEnergyLep_TwoPhotonCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoShowerCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_TwoShowerCut_NameAfterMode,VisibleEnergyLep_TwoShowerCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_PhotonFVCutAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_PhotonFVCut_NameAfterMode,VisibleEnergyLep_PhotonFVCut_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoLeptonsAfter_HistMode[Type][mode] = new TH1D(VisibleEnergyLep_TwoLeptons_NameAfterMode,VisibleEnergyLep_TwoLeptons_NameAfterMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);



	    }
	  }

	  if(fConfig.FillIntTypeHistograms){
	    for(int inttype=1000; inttype<1101; ++inttype){

	      std::string  VisibleEnergyLep_StringIntType                  = VisibleEnergyLep_String + " " +  std::to_string(inttype);                  
	      std::string  VisibleEnergyLep_AVCut_StringIntType            = VisibleEnergyLep_AVCut_String + " " +  std::to_string(inttype);           
	      std::string  VisibleEnergyLep_FVCut_StringIntType            = VisibleEnergyLep_FVCut_String + " " +  std::to_string(inttype);            
	      std::string  VisibleEnergyLep_FVBefore_StringIntType            = VisibleEnergyLep_FVBefore_String + " " +  std::to_string(inttype);            

	      std::string  VisibleEnergyLep_EnergyCut_StringIntType        = VisibleEnergyLep_EnergyCut_String + " " +  std::to_string(inttype);        
	      std::string  VisibleEnergyLep_PhotonEnergyCut_StringIntType  = VisibleEnergyLep_PhotonEnergyCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_ConversionGapCut_StringIntType = VisibleEnergyLep_ConversionGapCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_MuLenghtCut_StringIntType      = VisibleEnergyLep_MuLenghtCut_String + " " +  std::to_string(inttype);      
	      std::string  VisibleEnergyLep_NCCut_StringIntType            = VisibleEnergyLep_NCCut_String + " " +  std::to_string(inttype);            
	      std::string  VisibleEnergyLep_Selection_StringIntType        = VisibleEnergyLep_Selection_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_LeptonPlusPhotonCut_StringIntType = VisibleEnergyLep_LeptonPlusPhotonCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_NoShowerCut_StringIntType      = VisibleEnergyLep_NoShowerCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_TwoPhotonCut_StringIntType     = VisibleEnergyLep_TwoPhotonCut_String + " " +  std::to_string(inttype);
	      std::string   VisibleEnergyLep_TwoShowerCut_StringIntType     = VisibleEnergyLep_TwoShowerCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_PhotonFVCut_StringIntType      = VisibleEnergyLep_PhotonFVCut_String + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_TwoLeptons_StringIntType      = VisibleEnergyLep_PhotonFVCut_String + " " +  std::to_string(inttype);


	      const char* VisibleEnergyLep_NameIntType                  = VisibleEnergyLep_StringIntType.c_str();
	      const char* VisibleEnergyLep_AVCut_NameIntType            = VisibleEnergyLep_AVCut_StringIntType.c_str();
	      const char* VisibleEnergyLep_FVCut_NameIntType            = VisibleEnergyLep_FVCut_StringIntType.c_str();
	      const char* VisibleEnergyLep_FVBefore_NameIntType            = VisibleEnergyLep_FVBefore_StringIntType.c_str();

	      const char* VisibleEnergyLep_EnergyCut_NameIntType        = VisibleEnergyLep_EnergyCut_StringIntType.c_str();
	      const char* VisibleEnergyLep_PhotonEnergyCut_NameIntType  = VisibleEnergyLep_PhotonEnergyCut_StringIntType.c_str();
	      const char* VisibleEnergyLep_ConversionGapCut_NameIntType = VisibleEnergyLep_ConversionGapCut_StringIntType.c_str();
	      const char* VisibleEnergyLep_MuLenghtCut_NameIntType      = VisibleEnergyLep_MuLenghtCut_StringIntType.c_str();
	      const char* VisibleEnergyLep_NCCut_NameIntType            = VisibleEnergyLep_NCCut_StringIntType.c_str();
	      const char* VisibleEnergyLep_Selection_NameIntType        = VisibleEnergyLep_Selection_StringIntType.c_str();
	      const char*  VisibleEnergyLep_NoShowerCut_NameIntType      = VisibleEnergyLep_NoShowerCut_StringIntType.c_str();
	      const char*  VisibleEnergyLep_TwoPhotonCut_NameIntType     = VisibleEnergyLep_TwoPhotonCut_StringIntType.c_str();
	      const char*  VisibleEnergyLep_TwoShowerCut_NameIntType     = VisibleEnergyLep_TwoShowerCut_StringIntType.c_str();
	      const char*  VisibleEnergyLep_PhotonFVCut_NameIntType      = VisibleEnergyLep_TwoShowerCut_StringIntType.c_str();
	      const char*  VisibleEnergyLep_TwoLeptons_NameIntType      = VisibleEnergyLep_TwoLeptons_StringIntType.c_str();

	      const char* VisibleEnergyLep_LeptonPlusPhotonCut_NameIntType = VisibleEnergyLep_LeptonPlusPhotonCut_StringIntType.c_str();

	      fRootHists.VisibleEnergyLep_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_NameIntType,VisibleEnergyLep_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_AVCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_AVCut_NameIntType,VisibleEnergyLep_AVCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_FVCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_FVCut_NameIntType,VisibleEnergyLep_FVCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_FVBefore_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_FVBefore_NameIntType,VisibleEnergyLep_FVBefore_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      fRootHists.VisibleEnergyLep_EnergyCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_EnergyCut_NameIntType,VisibleEnergyLep_EnergyCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_PhotonEnergyCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_PhotonEnergyCut_NameIntType,VisibleEnergyLep_PhotonEnergyCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_ConversionGapCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_ConversionGapCut_NameIntType,VisibleEnergyLep_ConversionGapCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_MuLenghtCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_MuLenghtCut_NameIntType,VisibleEnergyLep_MuLenghtCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_NCCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_NCCut_NameIntType,VisibleEnergyLep_NCCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_Selection_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_Selection_NameIntType,VisibleEnergyLep_Selection_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_LeptonPlusPhotonCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_LeptonPlusPhotonCut_NameIntType,VisibleEnergyLep_LeptonPlusPhotonCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      
	      fRootHists.VisibleEnergyLep_NoShowerCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_NoShowerCut_NameIntType,VisibleEnergyLep_NoShowerCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoPhotonCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_TwoPhotonCut_NameIntType,VisibleEnergyLep_TwoPhotonCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoShowerCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_TwoShowerCut_NameIntType,VisibleEnergyLep_TwoShowerCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_PhotonFVCut_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_PhotonFVCut_NameIntType,VisibleEnergyLep_PhotonFVCut_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoLeptons_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_TwoLeptons_NameIntType,VisibleEnergyLep_TwoLeptons_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      std::string  VisibleEnergyLep_StringAfterIntType                  = VisibleEnergyLep_StringAfter + " " +  std::to_string(inttype);                  
	      std::string  VisibleEnergyLep_AVCut_StringAfterIntType            = VisibleEnergyLep_AVCut_StringAfter + " " +  std::to_string(inttype);           
	      std::string  VisibleEnergyLep_FVCut_StringAfterIntType            = VisibleEnergyLep_FVCut_StringAfter + " " +  std::to_string(inttype);            
	      std::string  VisibleEnergyLep_FVBefore_StringAfterIntType            = VisibleEnergyLep_FVBefore_StringAfter + " " +  std::to_string(inttype);            

	      std::string  VisibleEnergyLep_EnergyCut_StringAfterIntType        = VisibleEnergyLep_EnergyCut_StringAfter + " " +  std::to_string(inttype);        
	      std::string  VisibleEnergyLep_PhotonEnergyCut_StringAfterIntType  = VisibleEnergyLep_PhotonEnergyCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_ConversionGapCut_StringAfterIntType = VisibleEnergyLep_ConversionGapCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_MuLenghtCut_StringAfterIntType      = VisibleEnergyLep_MuLenghtCut_StringAfter + " " +  std::to_string(inttype);      
	      std::string  VisibleEnergyLep_NCCut_StringAfterIntType            = VisibleEnergyLep_NCCut_StringAfter + " " +  std::to_string(inttype);            
	      std::string  VisibleEnergyLep_Selection_StringAfterIntType        = VisibleEnergyLep_Selection_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_LeptonPlusPhotonCut_StringAfterIntType = VisibleEnergyLep_LeptonPlusPhotonCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_NoShowerCut_StringAfterIntType      = VisibleEnergyLep_NoShowerCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_TwoPhotonCut_StringAfterIntType     = VisibleEnergyLep_TwoPhotonCut_StringAfter + " " +  std::to_string(inttype);
	      std::string   VisibleEnergyLep_TwoShowerCut_StringAfterIntType     = VisibleEnergyLep_TwoShowerCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_PhotonFVCut_StringAfterIntType      = VisibleEnergyLep_PhotonFVCut_StringAfter + " " +  std::to_string(inttype);
	      std::string  VisibleEnergyLep_TwoLeptons_StringAfterIntType      = VisibleEnergyLep_PhotonFVCut_StringAfter + " " +  std::to_string(inttype);




	      const char* VisibleEnergyLep_NameAfterIntType                  = VisibleEnergyLep_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_AVCut_NameAfterIntType            = VisibleEnergyLep_AVCut_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_FVCut_NameAfterIntType            = VisibleEnergyLep_FVCut_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_FVBefore_NameAfterIntType            = VisibleEnergyLep_FVBefore_StringAfterIntType.c_str();

	      const char* VisibleEnergyLep_EnergyCut_NameAfterIntType        = VisibleEnergyLep_EnergyCut_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_PhotonEnergyCut_NameAfterIntType  = VisibleEnergyLep_PhotonEnergyCut_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_ConversionGapCut_NameAfterIntType = VisibleEnergyLep_ConversionGapCut_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_MuLenghtCut_NameAfterIntType      = VisibleEnergyLep_MuLenghtCut_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_NCCut_NameAfterIntType            = VisibleEnergyLep_NCCut_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_Selection_NameAfterIntType        = VisibleEnergyLep_Selection_StringAfterIntType.c_str();
	      const char* VisibleEnergyLep_LeptonPlusPhotonCut_NameAfterIntType = VisibleEnergyLep_LeptonPlusPhotonCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergyLep_NoShowerCut_NameAfterIntType      = VisibleEnergyLep_NoShowerCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergyLep_TwoPhotonCut_NameAfterIntType     = VisibleEnergyLep_TwoPhotonCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergyLep_TwoShowerCut_NameAfterIntType     = VisibleEnergyLep_TwoShowerCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergyLep_PhotonFVCut_NameAfterIntType      = VisibleEnergyLep_TwoShowerCut_StringAfterIntType.c_str();
	      const char*  VisibleEnergyLep_TwoLeptons_NameAfterIntType      = VisibleEnergyLep_TwoLeptons_StringAfterIntType.c_str();


	      fRootHists.VisibleEnergyLepAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_NameAfterIntType,VisibleEnergyLep_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_AVCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_AVCut_NameAfterIntType,VisibleEnergyLep_AVCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_FVCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_FVCut_NameAfterIntType,VisibleEnergyLep_FVCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_FVBeforeAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_FVBefore_NameAfterIntType,VisibleEnergyLep_FVBefore_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      fRootHists.VisibleEnergyLep_EnergyCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_EnergyCut_NameAfterIntType,VisibleEnergyLep_EnergyCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_PhotonEnergyCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_PhotonEnergyCut_NameAfterIntType,VisibleEnergyLep_PhotonEnergyCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_ConversionGapCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_ConversionGapCut_NameAfterIntType,VisibleEnergyLep_ConversionGapCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_MuLenghtCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_MuLenghtCut_NameAfterIntType,VisibleEnergyLep_MuLenghtCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_NCCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_NCCut_NameAfterIntType,VisibleEnergyLep_NCCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_SelectionAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_Selection_NameAfterIntType,VisibleEnergyLep_Selection_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      fRootHists.VisibleEnergyLep_LeptonPlusPhotonCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_LeptonPlusPhotonCut_NameAfterIntType,VisibleEnergyLep_LeptonPlusPhotonCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      fRootHists.VisibleEnergyLep_NoShowerCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_NoShowerCut_NameAfterIntType,VisibleEnergyLep_NoShowerCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoPhotonCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_TwoPhotonCut_NameAfterIntType,VisibleEnergyLep_TwoPhotonCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoShowerCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_TwoShowerCut_NameAfterIntType,VisibleEnergyLep_TwoShowerCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_PhotonFVCutAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_PhotonFVCut_NameAfterIntType,VisibleEnergyLep_PhotonFVCut_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      fRootHists.VisibleEnergyLep_TwoLeptonsAfter_HistIntType[Type][inttype] = new TH1D(VisibleEnergyLep_TwoLeptons_NameAfterIntType,VisibleEnergyLep_TwoLeptons_NameAfterIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	    }
	  }
	}
      }
      
      
      
      fRootHists.VisibleEnergy_CosmicFVCut_Hist = new TH1D("VisibleEnergy_CosmicFVCut","VisibleEnergy_CosmicFVCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicAVCut_Hist = new TH1D("VisibleEnergy_CosmicAVCut","VisibleEnergy_CosmicAVCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);

      fRootHists.VisibleEnergy_CosmicClyinderCut_Hist = new TH1D("VisibleEnergy_CosmicClyinderCut","VisibleEnergy_CosmicClyinderCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicdEdxCut_Hist = new TH1D("VisibleEnergy_CosmicdEdxCut","VisibleEnergy_CosmicdEdxCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicWeightCut_Hist = new TH1D("VisibleEnergy_CosmicWeightCut","VisibleEnergy_CosmicWeightCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicEnergyCut_Hist = new TH1D("VisibleEnergy_CosmicEnergyCut","VisibleEnergy_CosmicEnergyCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicSelection_Hist = new TH1D("Cosmic VisibleEnergy_Selection","Cosmic VisibleEnergy_Selection",fRootHists.ebins,fRootHists.emin,fRootHists.emax);

     fRootHists.VisibleEnergy_SpillWindow_Hist = new TH1D("Cosmic VisibleEnergy_SpillWindow","Cosmic VisibleEnergy_SpillWindow",fRootHists.ebins,fRootHists.emin,fRootHists.emax);

     fRootHists.NeutrinoT0 = new TH1D("Neutrino T0","Neutrino T0",100,0,5000);
     fRootHists.CosmicShowerT0 = new TH1D("Cosic T0","Cosmic T0",200,-500000,500000);

    }

    bool NueSelection::SelectCosmic(const sim::MCShower& mcs, std::map<int, const simb::MCParticle*>& mcparticles, NueSelection::NueInteraction& intInfo, std::vector<simb::MCTruth>& mcneutrinotruth ){

      dirtevent = false;

      //Check the shower is not a neutrino shower. 
      if(mcs.Origin() ==  simb::kBeamNeutrino){return false;}
      
      if(fConfig.Verbose){std::cout << "Is a Cosmic" << std::endl;}

      if(CosmicInAV(mcs,mcparticles[mcs.TrackID()])){
	fRootHists.VisibleEnergy_CosmicAVCut_Hist->Fill(mcparticles[mcs.TrackID()]->E(),intInfo.weight);
      }


      //Vertex for the cosmic shower start.
      TVector3 vertex;
      double   time;

      //Check if it is in the fiducal volume.
      bool Pass_CosmicInFV = PassCosmicInFV(mcs,vertex,mcparticles[mcs.TrackID()],time);

      if(mcparticles[mcs.TrackID()]->E() > 0.1){

	if(CosmicInAV(mcs,mcparticles[mcs.TrackID()])){

	  fRootHists.XZCosmic->Fill(vertex.X(),vertex.Z());
	  fRootHists.YZCosmic->Fill(vertex.Y(),vertex.Z());
	  fRootHists.XYCosmic->Fill(vertex.X(),vertex.Y());
	  fRootHists.DCACosmic->Fill(DistanceToClosestSurface(vertex));
		
	  fRootHists.StartXZCosmic->Fill(mcs.Start().Position().Vect().X(),mcs.Start().Position().Vect().Z());
	  fRootHists.StartYZCosmic->Fill(mcs.Start().Position().Vect().Y(),mcs.Start().Position().Vect().Z());
	  fRootHists.StartXYCosmic->Fill(mcs.Start().Position().Vect().X(),mcs.Start().Position().Vect().Y());
		
	  fRootHists.StartMCPXZCosmic->Fill(mcparticles[mcs.TrackID()]->Position().Vect().X(),mcparticles[mcs.TrackID()]->Position().Vect().Z());
	  fRootHists.StartMCPYZCosmic->Fill(mcparticles[mcs.TrackID()]->Position().Vect().Y(),mcparticles[mcs.TrackID()]->Position().Vect().Z());
	  fRootHists.StartMCPXYCosmic->Fill(mcparticles[mcs.TrackID()]->Position().Vect().X(),mcparticles[mcs.TrackID()]->Position().Vect().Y());
		
	  fRootHists.EndMCPXZCosmic->Fill(mcparticles[mcs.TrackID()]->EndPosition().Vect().X(),mcparticles[mcs.TrackID()]->Position().Vect().Z());
	  fRootHists.EndMCPYZCosmic->Fill(mcparticles[mcs.TrackID()]->EndPosition().Vect().Y(),mcparticles[mcs.TrackID()]->Position().Vect().Z());
	  fRootHists.EndMCPXYCosmic->Fill(mcparticles[mcs.TrackID()]->EndPosition().Vect().X(),mcparticles[mcs.TrackID()]->Position().Vect().Y());
	}



      }	

      if(!Pass_CosmicInFV){
	if(fConfig.Verbose){std::cout << "Cosmic shower not in the FV removed" << std::endl;}
	return false;
      }
      fRootHists.VisibleEnergy_CosmicFVCut_Hist->Fill(mcparticles[mcs.TrackID()]->E(),intInfo.weight);


      fRootHists.CosmicShowerT0->Fill(time);



      //Check if the cosmic is within the spill time window.
      bool pass_cosmic_spillwindow_cut  = PassCosmicInSpillWindow(mcs,mcparticles,intInfo,mcneutrinotruth,time);
      if(!pass_cosmic_spillwindow_cut){
	if(fConfig.Verbose){
	  std::cout << "Cosmic Not in spill window" << std::endl;
	}
	return false;
      }
      fRootHists.VisibleEnergy_SpillWindow_Hist->Fill(mcparticles[mcs.TrackID()]->E(),intInfo.weight);


      //If we have one photon candndiate choose that to be the energy of the lepton 
      VisibleEnergyCalculator calculator;
      calculator.lepton_pdgid = 22;
      calculator.track_threshold =  fConfig.trackVisibleEnergyThreshold;
      calculator.shower_energy_distortion = fConfig.showerEnergyDistortion;
      calculator.track_energy_distortion = fConfig.trackEnergyDistortion;
      calculator.lepton_energy_distortion_contained = fConfig.leptonEnergyDistortionContained;
      calculator.lepton_energy_distortion_leaving_A = 1; //Not required here.
      calculator.lepton_energy_distortion_leaving_B = 1; //Not required here. 
      calculator.lepton_contained = true; // We don't worry at the moment if the electron is contained at this moment.
      calculator.lepton_index = mcs.TrackID();
      
      //Add the energy 
      intInfo.leptonic_energy = smearLeptonEnergy(rand,mcparticles[mcs.TrackID()],calculator);
      intInfo.leptonpdgID = mcs.PdgCode();

      
      bool pass_eEnergyCut = passeEnergyCut(intInfo.leptonic_energy);
      if(!pass_eEnergyCut){
       	if(fConfig.Verbose){std::cout << "Cosmic Shower was too low in energy. Event not Selected" << std::endl;}
       	return false;
      }
      
      fRootHists.VisibleEnergy_CosmicEnergyCut_Hist->Fill(mcparticles[mcs.TrackID()]->E(),intInfo.weight);
      if(fConfig.Verbose){std::cout << "Passed the Cosmic Energy Cut." << std::endl;}
           
      //Check if the cosmic passes the cylinder cut.
      bool Pass_CosmicCylinderCut = PassCosmicCylinderCut(mcs,vertex,mcparticles);
      if(!Pass_CosmicCylinderCut){
	if(fConfig.Verbose){std::cout << "Cosmic shower in the cyclinder" << std::endl;}
	return false;
      }
      fRootHists.VisibleEnergy_CosmicClyinderCut_Hist->Fill(mcparticles[mcs.TrackID()]->E(),intInfo.weight);

      //dEdx weight 
      bool pass_dEdxCut = passdEdxCut(mcs.PdgCode());
      if(!pass_dEdxCut){intInfo.weight *= fConfig.dEdxPhotonCut;}
      fRootHists.VisibleEnergy_CosmicdEdxCut_Hist->Fill(mcparticles[mcs.TrackID()]->E(),intInfo.weight);

      //Cosmic global weight for time and CRTs 
      intInfo.weight *= fConfig.CosmicWeight;
      fRootHists.VisibleEnergy_CosmicWeightCut_Hist->Fill(mcparticles[mcs.TrackID()]->E(),intInfo.weight);
                
      return true;
    }

    //Check the cosmic is in the spill window. 
    bool NueSelection::CosmicInSpillWindow(const sim::MCShower& mcs, std::map<int, const simb::MCParticle*>& mcparticles, NueSelection::NueInteraction& intInfo,std::vector<simb::MCTruth>& mcneutrinotruths, double& time){
      
      //      std::cout << "new cosmic " << mcparticles.size() << std::endl;

      //If we want to weight the the spillwindow/windowtime for better stats lets do that.
      if(fConfig.UseAllCosmics){
	//	std::cout << "should not be here" << std::endl;
	intInfo.weight *= fConfig.SpillTime/(1000*fConfig.ReadoutWindowSize);
	return true;
      }
	
      
      //Get a rough shower time (photon should decay within a few ns at most, which is within errro of the light system,
      float interaction_time = time;

      //See if its in the window
      if(interaction_time > fConfig.GlobalTimeOffset && interaction_time < fConfig.GlobalTimeOffset + fConfig.SpillTime){
	std::cout << "in time: " << std::endl;
	return true;
      }
      
      

      //If its out of time check that not another cosmic was in time nn
      for(auto const& mcparticle: mcparticles){
	
	//Particle must be charged 
	if(!abs(PDGCharge(mcparticle.second->PdgCode())) > 1e-4){continue;}

	//Particle must deposit some sort of energy.
	double mass = PDGMass(mcparticle.second->PdgCode());
	if(mcparticle.second->E()*1000 - mass < 8.3e-6){continue;}

	//Option to ignore neutrino 
	if(fConfig.IgnoreNeutrinoDepsInCosmic){
	  if(isFromNuVertex(mcneutrinotruths,mcparticles,mcparticle.second->TrackId())){
	    std::cout << "From Vertex" << std::endl;
	    continue;
	  }
	}

	const TLorentzVector StartPosition =  mcparticle.second->Position(); 
	const TLorentzVector EndPosition = mcparticle.second->EndPosition();
	
	//Check the times
	if(EndPosition.T() < fConfig.GlobalTimeOffset){
	  //	  std::cout << " endpos time: " << EndPosition.T() << "  fConfig.GlobalTimeOffset: " <<  fConfig.GlobalTimeOffset << std::endl;
	  continue;
	}
	if(StartPosition.T() > fConfig.GlobalTimeOffset + fConfig.SpillTime){
	  //	  std::cout << " startpos time: " << EndPosition.T() << "  fConfig.GlobalTimeOffset: " <<  fConfig.GlobalTimeOffset << std::endl;

	  continue;
	}
	//	std::cout << "cosmic in time" << std::endl;
	
	//Check if any of the traj points are in the tpc 
	int numpoints = mcparticle.second->NumberTrajectoryPoints();
	int traj =0;
	for(;traj<numpoints; ++traj){
	  //	  std::cout << "position: mcparticle.second->Position(traj).Vect().X(): " << mcparticle.second->Position(traj).Vect().X() << " " << mcparticle.second->Position(traj).Vect().Y() << " "  << mcparticle.second->Position(traj).Vect().Z() << std::endl;
	  if(containedInAV(mcparticle.second->Position(traj).Vect())){break;}
	}
	//Never was in the tpc
	if(traj = numpoints-1){
	  //	  std::cout << "concient cosmic not in AV" << std::endl;
	  continue;
	}

	std::cout << "particle in the AV at least time: " <<  mcparticle.second->Position(traj).T() << " shower time: " << interaction_time  << std::endl;

	//Now check the time of the traj point
	float part_interaction_time = mcparticle.second->Position(traj).T();
	if(part_interaction_time > fConfig.GlobalTimeOffset && part_interaction_time < fConfig.GlobalTimeOffset + fConfig.SpillTime){
	  std::cout << "cosmic was in splill time in the AV" << std::endl;
	  return true;
	}
      }
      return false;
    }
    
    //See if a shower is in the FV.
    bool NueSelection::CosmicInFV(const sim::MCShower& mcs, TVector3& vertex, const simb::MCParticle*& mcparticle, double& time){

      if(TMath::Abs(mcs.PdgCode()) == 11){
	vertex = mcs.Start().Position().Vect();
	time   = mcparticle->T();
      }
      else if(mcs.PdgCode() == 22){
	if (mcparticle->EndProcess()=="conv"){
	  vertex = mcs.End().Position().Vect();
	  time   = mcparticle->EndT();
	}
	else{
	  vertex = mcs.End().Position().Vect();
	  time   = mcparticle->EndT();
	  for(int i=0; i<mcparticle->NumberTrajectoryPoints();++i){
	    if(mcparticle->E(i) < 0.1*mcparticle->E()){
	      vertex = mcparticle->Position(i).Vect();
	      time   = mcparticle->T(i);
	      break;
	    }
	  }
	}
      }
      else{
	std::cerr << "what kind of shower is this? :S" << std::endl;
	return false;
      }
      return containedInFV(vertex);
    }


    bool NueSelection::CosmicInAV(const sim::MCShower& mcs, const simb::MCParticle*& mcparticle){
      
      TVector3 vertex;
      if(TMath::Abs(mcs.PdgCode()) == 11){
	vertex = mcs.Start().Position().Vect();
      }
      else if(mcs.PdgCode() == 22){
	if (mcparticle->EndProcess()=="conv"){
	  vertex = mcs.End().Position().Vect();
	}
	else{
	  vertex = mcs.End().Position().Vect();
	  for(int i=0; i<mcparticle->NumberTrajectoryPoints();++i){
	    if(mcparticle->E(i) < 0.1*mcparticle->E()){
	      vertex = mcparticle->Position(i).Vect();
	      break;
	    }
	  }
	}
      }
      else{
	std::cerr << "what kind of shower is this? :S" << std::endl;
	return false;
      }
      return containedInAV(vertex);
    }


    bool NueSelection::CosmicCylinderCut(const sim::MCShower& mcs, TVector3& vertex, std::map<int, const simb::MCParticle*>& mcparticles){
      
      //Get the mother track. If it doesn't have a mother we cannot do this. 
      if(mcs.MotherTrackID() == -1){return true;}

      if(mcparticles.find(mcs.MotherTrackID()) == mcparticles.end()){
	if(fConfig.Verbose){std::cout << "Warning! Shower is not associated to mother we cannot track this" << std::endl;}
	return true;
      }

      const simb::MCParticle* mother = mcparticles[mcs.MotherTrackID()];
      if(fConfig.Verbose){std::cout << "cosmic mother is pdgcode: " << mother->PdgCode() << " track id: " << mother->TrackId() << " and energy: " << mother->E() << std::endl;} 
      
      //Mother might have been neutral but cam from a muon track 
      while(!abs(PDGCharge(mother->PdgCode())) > 1e-4){
	if(mcparticles.find(mother->Mother()) == mcparticles.end()){break;}
	mother = mcparticles[mother->Mother()];
      }

      //Cant do this if its not a muon 
      if(TMath::Abs(mother->PdgCode()) != 13){return true;}
      
      //Muon will never be not in the FV with the shower being in this case
      //Check if the vertex is within a cylinder of the muon. Assuming the muon is just a straight line then work out mag from start to vector and star and end
      int numpoints = mother->NumberTrajectoryPoints();

      TVector3 StartPos   = mother->Position().Vect();
      TVector3 EndPos     = mother->EndPosition().Vect();

      TVector3 UnitDir = mother->Momentum().Vect().Unit();
      double vertex_proj = (vertex-StartPos).Dot(UnitDir);
      

      //Loop over the trajectory points and find the two which sandwich the vertex
      for(int traj=0; traj<numpoints-1; ++traj){
      	TVector3 TrajPointBefore   = mother->Position(traj).Vect();
      	TVector3 TrajPointAfter    = mother->Position(traj+1).Vect();

      	double trajbefore_proj =  (TrajPointBefore-StartPos).Dot(UnitDir);
      	double trajafter_proj =  (TrajPointAfter-StartPos).Dot(UnitDir);

      	TrajPointBefore.Dot(UnitDir);

      	if(trajafter_proj > vertex_proj){
	  StartPos = TrajPointBefore;
      	  EndPos = TrajPointAfter;
      	  break;
      	}
      }

      TVector3 Mag = (StartPos - EndPos);
      double DistFromMother = ((vertex-StartPos) - (vertex-StartPos).Dot(Mag)*Mag).Mag();
      
      // double   StartVtxMag = (StartPos - vertex).Mag();
      // double   StartEndMag = (StartPos - EndPos).Mag();
	   
      // //What is the angle
      // double AngleStart = (vertex - StartPos).Angle(EndPos - StartPos);
      // if(AngleStart > TMath::Pi()){AngleStart = 2*TMath::Pi() - AngleStart;}
      // if(AngleStart > TMath::Pi()/2){
      // 	//angle is above the start of the track
      // 	std::cout << "Dom Check the angle" << std::endl;
      // 	if(fConfig.Verbose){std::cout << "Shower above start of the track" << std::endl;}
      // 	return true;
      // }
	  
      // double AngleEnd = (vertex - EndPos).Angle(StartPos - EndPos);
      // if(AngleEnd > TMath::Pi()){AngleEnd = 2*TMath::Pi() - AngleEnd;}
      // if(AngleEnd > TMath::Pi()/2){
      // 	//angle is below the end of the track    
      // 	if(fConfig.Verbose){std::cout << "Shower below end of the track" << std::endl;}
      // 	std::cout << "Dom Check the angle for end" << std::endl;
      // 	return true;
      // }
      
      // //Simple trig to work out the distance from the vertex to the track now.
      // double DistFromMother = StartVtxMag*TMath::Sin(AngleStart);
      // if(fConfig.Verbose){std::cout << "Angle from the start is: " << AngleStart << " Hypot is: " << StartVtxMag << std::endl;} 
      
      if(DistFromMother > fConfig.CosmicVolumeRadius){return true;}
      
      return false;
    }

    //Check if the point is in the Active volume.
    double  NueSelection::DistanceToClosestSurface(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& AV: fConfig.active_volumes) {
	if (AV.Contain(p)){
	  const geoalgo::Point_t min_p =  AV.Min();
	  const geoalgo::Point_t max_p =  AV.Max();
	
	  const TVector3 min = {min_p.at(0), min_p.at(1), min_p.at(2)}; 
	  const TVector3 max = {max_p.at(0), max_p.at(1), max_p.at(2)}; 
 	
	  const double x_length = TMath::Abs((v.X()-min.X()));
	  const double y_length =  TMath::Abs((v.Y()-min.Y()));
	  const double z_length =  TMath::Abs((v.Z()-min.Z()));

	  const double x_length_max =  TMath::Abs((v.X()-max.X()));
	  const double y_length_max =  TMath::Abs((v.Y()-max.Y()));
	  const double z_length_max =  TMath::Abs((v.Z()-max.Z()));

	  double x_pos = x_length_max;
	  double y_pos = y_length_max;  
	  double z_pos = z_length_max;  
	
	  if(x_length < x_length_max){
	    x_pos = x_length;
	  }
	  if(y_length < y_length_max){
	    y_pos = y_length;
	  }
	  if(z_length < z_length_max){
	    z_pos = z_length;
	  }


	  if(x_pos < y_pos){
	    if(x_pos < z_pos){
	      return x_pos;
	    }
	    else{
	      return z_pos;
	    }
	  }
	  else{
	    if(y_pos < z_pos){
	      return y_pos;
	    }
	    else{
	      return z_pos;
	    }
	  }
	}
      }
      return  -999;
    }


    
  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelection)

