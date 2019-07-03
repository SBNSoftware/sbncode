#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

void NueRatioPlots(){

  Double_t xerr[11] =  {0.075,0.075,0.075,0.075,0.075,0.075,0.075,0.1,0.125,0.125,0.5};
  Double_t nue_bins[11] = {0.275,0.425,0.575,0.725,0.875,1.025,1.2,1.400,1.625,1.875,2.5};

  //Int
  Double_t sbnd_int_new[11] = {
    4244.44,
    7670.06,
    9193.13,
    9373.18,
    9065.31,
    8185.25,
    7159.19,
    6165.64,
    4960.89,
    3927.56,
    2087.04,
  };
  Double_t sbnd_int_prop[11] = {
    5157.0102,
    8748.3414,
    11012.8262,
    11065.9, 
    10092.8793,
    9606.3689,
    8084.9182,
    6775.7629,
    5749.6683,
    4015.9222,
    1813.3569
  };
  Double_t sbnd_int_error[11] = {
    59.3959,
    79.8538,
    87.4263,
    88.2801,
    86.819,
    82.4983,
    66.8181,
    62.0087,
    49.7497,
    44.2664,
    16.1341,
  };

  Double_t sbnd_int_ratio[11];
  for(int i=0; i<11; ++i){
    sbnd_int_prop[i] = sbnd_int_prop[i];
    
    if(sbnd_int_prop[i] > 1){
      sbnd_int_error[i] = sbnd_int_error[i]/sbnd_int_prop[i];   
      sbnd_int_ratio[i] = sbnd_int_new[i]/sbnd_int_prop[i];
    }
    else if(sbnd_int_prop[i] < 1 && sbnd_int_new[i] < 1){ 
      sbnd_int_ratio[i] = 1;
    }
    else{
      sbnd_int_ratio[i] = 0;
    }
  }

    
  auto sbnd_int_canvas = new TCanvas("sbnd_int_canvas","sbnd_int_canvas",800,500);
  auto sbnd_int_graph= new TGraphErrors(11,nue_bins,sbnd_int_ratio,xerr,sbnd_int_error);
  sbnd_int_graph->SetTitle("Option ACP example");
  sbnd_int_graph->GetXaxis()->SetTitle("X title");
  sbnd_int_graph->GetYaxis()->SetTitle("Y title");
  sbnd_int_graph->Draw("AP");
  sbnd_int_canvas->Draw();
    
    
  //NC
  Double_t sbnd_nc_new[11] = { 
    6359.2,
    3471.71, 
    1340.07,
    712.489,
    148.821,
    84.8202,
    50.8921,
    42.4101,
    20.3568,
    16.964,
    8.4820,
  };

  Double_t sbnd_nc_prop[11] = {
    8553.7372,
    10880.1415,
    11791.2428, 
    11587.793, 
    10375.9398,
    9783.2817, 
    8164.529, 
    6961.5215, 
    5926.5812, 
    4130.9155, 
    1813.3569
  };
  Double_t sbnd_nc_ratio[11];

  Double_t sbnd_nc_error[11] = {
    247.131,
    307.918,
    202.454,
    170.518,
    29.2394,
    21.9005,
    14.6913,
    13.4112,
    8.31064,
    7.58654,
    2.68225,
  };

  for(int i=0; i<11; ++i){
    sbnd_nc_prop[i] = sbnd_nc_prop[i] - sbnd_int_prop[i];
    std::cout << "sbnd_nc_prop[i]: " << sbnd_nc_prop[i] << std::endl;
    if(sbnd_nc_prop[i] > 1){
      sbnd_nc_ratio[i] = sbnd_nc_new[i]/sbnd_nc_prop[i]; 
      sbnd_nc_error[i] = sbnd_nc_error[i]/sbnd_nc_prop[i];
    }
    else if(sbnd_nc_prop[i] < 1 && sbnd_nc_new[i] <1){ 
       std::cout << "test2" << std::endl;
      sbnd_nc_ratio[i] = 1;
    }
    else{
       std::cout << "test3" << std::endl;
      sbnd_nc_ratio[i] = -999;
    }
    sbnd_nc_prop[i] = sbnd_nc_prop[i] + sbnd_int_prop[i];
  };
    
  auto sbnd_nc_canvas = new TCanvas("sbnd_nc_canvas","sbnd_nc_canvas",800,500);
  auto sbnd_nc_graph= new TGraphErrors(11,nue_bins,sbnd_nc_ratio,xerr,sbnd_nc_error);
  sbnd_nc_graph->SetTitle("Option ACP example");
  sbnd_nc_graph->GetXaxis()->SetTitle("X title");
  sbnd_nc_graph->GetYaxis()->SetTitle("Y title");
  sbnd_nc_graph->Draw("AP");
  sbnd_nc_canvas->Draw();


  //Numu 
  Double_t sbnd_numu_new[11] = {
    141.367,
    906.633,
    1121.43,
    1366.55,
    422.216,
    365.669,
    352.004,
    84.8202,
    68.9356,
    90.4748,
    112.528,
  };

  Double_t sbnd_numu_prop[11] = {
    9340.9996, 
    11587.793,
    12242.3706,
    11817.7797,
    10535.1614,
    9783.2817, 
    8164.529, 
    6961.5215,
    6050.4202,
    4290.1371,
    1813.3569,
  };
  Double_t sbnd_numu_ratio[11];
  Double_t sbnd_numu_error[11] = {
    28.2734,
    147.733,
    177.206,
    257.238,
    103.618,
    102.063,
    124.828,
    18.9664,
    15.4503,
    57.5556,
    37.5557,
  };

  for(int i=0; i<11; ++i){
    sbnd_numu_prop[i] = sbnd_numu_prop[i] - sbnd_nc_prop[i];
    std::cout << "sbnd_numu_prop[i]: " << sbnd_numu_prop[i] << " sbnd_numu_new[i]: " << sbnd_numu_new[i]  << std::endl; 
    if(sbnd_numu_prop[i] > 1){
      sbnd_numu_ratio[i] = sbnd_numu_new[i]/sbnd_numu_prop[i];
      sbnd_numu_error[i] = sbnd_numu_error[i]/sbnd_numu_prop[i];
    }
    else if(sbnd_numu_prop[i] < 1 && sbnd_numu_new[i] <1){ 
      sbnd_numu_ratio[i] = 1;
    }
    else{
      sbnd_numu_ratio[i] = -999;
    }
  }
    
  auto sbnd_numu_canvas = new TCanvas("sbnd_numu_canvas","sbnd_numu_canvas",800,500);
  auto sbnd_numu_graph = new TGraphErrors(11,nue_bins,sbnd_numu_ratio,xerr,sbnd_numu_error);
  sbnd_numu_graph->SetTitle("Option ACP example");
  sbnd_numu_graph->GetXaxis()->SetTitle("X title");
  sbnd_numu_graph->GetYaxis()->SetTitle("Y title");
  sbnd_numu_graph->Draw("AP");
  sbnd_numu_canvas->Draw();
  

}


// Bin Entry is: 202.136 +- 3.13764
// Bin Entry is: 384.656 +- 4.32843
// Bin Entry is: 488.247 +- 4.87682
// Bin Entry is: 512.31 +- 4.99557
// Bin Entry is: 507.897 +- 4.97414
// Bin Entry is: 486.991 +- 4.87072
// Bin Entry is: 443.213 +- 4.02412
// Bin Entry is: 367.064 +- 3.66216
// Bin Entry is: 300.761 +- 2.96495
// Bin Entry is: 243.613 +- 2.66844
// Bin Entry is: 120.731 +- 0.939267
// NC F
// Bin Entry is: 276.025 +- 21.2407
// Bin Entry is: 209.232 +- 32.3669
// Bin Entry is: 105.121 +- 29.6212
// Bin Entry is: 53.951 +- 28.1721
// Bin Entry is: 14.9533 +- 4.931
// Bin Entry is: 6.67262 +- 3.28734
// Bin Entry is: 1.28843 +- 1.23279
// Bin Entry is: 3.72717 +- 2.13518
// Bin Entry is: 1.0058 +- 0.986207
// Bin Entry is: 0.0178233 +- 0.00563623
// Bin Entry is: 0.249666 +- 0.24655
// NuMu F
// Bin Entry is: 31.2293 +- 7.1645
// Bin Entry is: 41.0912 +- 8.21825
// Bin Entry is: 27.942 +- 6.77694
// Bin Entry is: 14.7928 +- 4.93095
// Bin Entry is: 14.7928 +- 4.93095
// Bin Entry is: 1.64365 +- 1.64365
// Bin Entry is: 7.39642 +- 3.01958
// Bin Entry is: 2.46547 +- 1.74335
// Bin Entry is: 0 +- 0
// Bin Entry is: 0.986189 +- 0.986189
// Bin Entry is: 0 +- 0


// Bin Entry is: 337.875 +- 11.708
// Bin Entry is: 726.722 +- 17.1729
// Bin Entry is: 928.808 +- 19.4152
// Bin Entry is: 991.243 +- 20.058
// Bin Entry is: 930.838 +- 19.4364
// Bin Entry is: 910.064 +- 19.219
// Bin Entry is: 764.115 +- 15.2515
// Bin Entry is: 691.376 +- 14.5072
// Bin Entry is: 580.606 +- 11.891
// Bin Entry is: 428.883 +- 10.2199
// Bin Entry is: 219.069 +- 3.65204
// NC F
// Bin Entry is: 542.562 +- 63.5844
// Bin Entry is: 279.364 +- 45.5813
// Bin Entry is: 173.503 +- 35.9376
// Bin Entry is: 60.1732 +- 21.1947
// Bin Entry is: 22.8566 +- 12.9794
// Bin Entry is: 30.0239 +- 14.9869
// Bin Entry is: 11.353 +- 7.94811
// Bin Entry is: 0.0564468 +- 0.0325896
// Bin Entry is: 0.0451575 +- 0.0260717
// Bin Entry is: 0 +- 0
// Bin Entry is: 0.0112894 +- 0.00651792
// NuMu F
// Bin Entry is: 44.9606 +- 18.3551
// Bin Entry is: 59.9474 +- 21.1946
// Bin Entry is: 37.4671 +- 16.7558
// Bin Entry is: 37.4671 +- 16.7558
// Bin Entry is: 22.4803 +- 12.979
// Bin Entry is: 14.9869 +- 10.5973
// Bin Entry is: 0 +- 0
// Bin Entry is: 0 +- 0
// Bin Entry is: 0 +- 0
// Bin Entry is: 4.49606 +- 4.49606
// Bin Entry is: 0 +- 0
