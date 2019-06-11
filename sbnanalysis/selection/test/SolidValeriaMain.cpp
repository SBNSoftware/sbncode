 #include "../include/CC0piAnalysisHelper.h"                                                    
 #include "../include/GeneralAnalysisHelper.h"                                                  
 #include "../include/EventSelectionTool.h"                                                     
 #include "../include/Event.h"                                                                  
 #include <iostream>                                                                            
 #include <sstream>                                                                             
 #include <numeric>                                                                             
 #include <time.h>                                                                              
 #include <algorithm>                                                                           
 #include "TVector3.h"                                                                          
 #include "TH1.h"                                                                               
 #include "TH2.h"                                                                               
#include "TGraph.h"
#include "TCanvas.h"                                                                           
 #include "TLegend.h"                                                                           
 #include "TLatex.h"                                                                            
 #include "TStyle.h"                                                                            
 #include "TColor.h"                                                                            
 #include "TObjArray.h"                                                                         
 #include "THStack.h"
 #include "TProfile.h"
 #include "TGraph.h"
 #include <math.h>
#include <vector>
using namespace selection;                                                                     
                                                                                                
 int MainTest(){                                                                                
                                                                                                
   time_t rawtime;                                                                              
   struct tm * timeinfo;                                                                        
   time (&rawtime);                                                                             
   timeinfo = localtime (&rawtime);                                                             
   std::cout << "-----------------------------------------------------------" << std::endl;     
   std::cout << " Start local time and date:  " << asctime(timeinfo)          << std::endl;     
   std::cout << "-----------------------------------------------------------" << std::endl;     
                                                                                                
   // Output file location                                                                      
   std::string angle = "../Output_Selection_Tool/plots/solid/angle_distributions/";   
   std::string momenta = "../Output_Selection_Tool/plots/solid/momenta_distributions/";
   std::string stats_location= "../Output_Selection_Tool/statistics/";
                                                                                                
   //Load events                                         
   
   // Initialise event list and the topology maps                                               
   EventSelectionTool::EventList events;
   int start = static_cast<int>(time(NULL));                                                    
   unsigned int total = 500 ;                                                                    
                                                                                               
   // Load the events into the event list                                                       
   for( unsigned int i = 0; i < total; ++i ){                                                   
                                                                                                
     std::string name;                                                                          
     name.clear();       
     char file_name[1024];                                                                      
     name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
     //name = "/pnfs/sbnd/persistent/users/rsjones/analysis_sample/120918_ana_files/11509725_"+std::to_string(i)+"/output_file.root";
     strcpy( file_name, name.c_str() );                                                         
    
     EventSelectionTool::LoadEventList(file_name, events, i);                                                                                                                              
     EventSelectionTool::GetTimeLeft(start,total,i);                                            
   }                                                                                            
   std::cout << std::endl;  

   TopologyMap cc_signal_map         = GeneralAnalysisHelper::GetCCIncTopologyMap(); 
   TopologyMap cc0pi2p_signal_map    = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap(); 
int tottrue      = 0; // counter of total CC true events
int totreco      = 0; // counter of total CC reconstructed events
int cc0pi2pmc    = 0; // counter of cc0pi2p true events
int cc0pi2pre    = 0; // counter of cc0pi2p reconstructed events  
int cc0pi2ps2    = 0; // counter of cc0pi2p signal events by protons
int mutrue       = 0; // counter of true muons
int mureco       = 0; // counter of reconstructed particles with escaping tracks 
int musignal     = 0; // counter of signal muons
int nproton      = 0; // counter for comparing true proton momenta
int nprotonreco  = 0; // counter for comparing reco proton momenta
int nprotonsig   = 0; // counter for comparing signal proton momenta
int mucount      = 0; // counter for true muons
int mucountr     = 0; // counter for reco muons
int mucounts     = 0; // counter for signal muons
int ptotreco     = 0; // counter of reconstructed particles with not escaping & w/ at least 5 hits tracks 
int ptottrue     = 0; // counter of true protons 
int ptotsig      = 0; // counter of signal protons
int truev        = 0; // counter of true events considering both protons and muons  
int ntotevents   = 0; // counter of total number of events 
int recopevent   = 0; // counter of total number of reconstructed events
Double_t momento; 
Double_t momento2;
Double_t momento3; 
std:: vector <float> mo;
std:: vector <Double_t> res, res2;
std:: vector <Double_t> zezi;
std:: vector <Double_t> mores;
std:: vector <Double_t> true_momentum_mu, true_momentum_p1, true_momentum_p2;
std:: vector <Double_t> reco_momentum_mu, reco_momentum_p1, reco_momentum_p2;
std:: vector <Double_t> signal_momentum_mu, signal_momentum_p1, signal_momentum_p2;

TVector3 momentum, momentum2;
TVector3 p1, p2, p1r, p2r, p1s, p2s; 

std:: vector <Double_t> angle_true;
std:: vector <Double_t> angle_reco;
std:: vector <Double_t> angle_sign;

std:: vector <int> interactions; 
interactions.push_back(1); // QE 
interactions.push_back(3); // DIS 
interactions.push_back(4); // RESONANCE
interactions.push_back(5); // COHERENT 
interactions.push_back(10);// MEC 
std:: vector < vector <Double_t> >  angle_int_sign;
std:: vector <Double_t> angle_QE_sign, angle_DIS_sign, angle_RES_sign, angle_COH_sign, angle_MEC_sign;
angle_int_sign.push_back(angle_QE_sign);
angle_int_sign.push_back(angle_DIS_sign);
angle_int_sign.push_back(angle_RES_sign);
angle_int_sign.push_back(angle_COH_sign);
angle_int_sign.push_back(angle_MEC_sign);
std:: vector < vector <Double_t> >  angle_int_reco;
std:: vector <Double_t> angle_QE_reco, angle_DIS_reco, angle_RES_reco, angle_COH_reco, angle_MEC_reco;
angle_int_reco.push_back(angle_QE_reco);
angle_int_reco.push_back(angle_DIS_reco);
angle_int_reco.push_back(angle_RES_reco);
angle_int_reco.push_back(angle_COH_reco);
angle_int_reco.push_back(angle_MEC_reco); 
std:: vector < vector <Double_t> >  angle_int_true;
std:: vector <Double_t> angle_QE_true, angle_DIS_true, angle_RES_true, angle_COH_true, angle_MEC_true;
angle_int_true.push_back(angle_QE_true);
angle_int_true.push_back(angle_DIS_true);
angle_int_true.push_back(angle_RES_true);
angle_int_true.push_back(angle_COH_true);
angle_int_true.push_back(angle_MEC_true);

Double_t angletrue, anglereco, anglesign;
ofstream file;
file.open(stats_location+"solid.txt");

// loop over the events 
for(const Event &e : events){
  if(e.IsSBNDTrueFiducial()){
    if(bool max = GeneralAnalysisHelper::MaxOneEscapingTrack(e)){
    ++ntotevents;
    ParticleList true_particles = e.GetMCParticleList();  
    ParticleList reco_particles = e.GetRecoParticleList();                                    
    if(e.CheckRecoTopology(cc_signal_map)){
      ++totreco;
    }
    if(e.CheckMCTopology(cc_signal_map)){
    ++tottrue;
    }

// cc0pi2P true events

  
    if(e.CheckMCTopology(cc0pi2p_signal_map)){
    ++cc0pi2pmc;  
    //loop over the true particles
    if(unsigned int j=e.CountMCParticlesWithPdg(2212)==2){ 
      for(Particle &p_true : true_particles){
        if(p_true.GetPdgCode() == 13){
        ++mutrue; 
        ++mucount;
        momento3=p_true.GetModulusMomentum();
        }
        else if(p_true.GetPdgCode()==2212 && p_true.GetKineticEnergy()>=0.021){
        if(p_true.GetModulusMomentum()>0.250 && p_true.GetModulusMomentum() < 100) {
          ++nproton; 
          ++ptottrue; 
          momentum=p_true.GetMomentum();
          mo.push_back(p_true.GetModulusMomentum());       
          for(int i=0; i<3; i++){
           res.push_back(momentum(i));
        }
        }
      }
  }
 }

      true_momentum_mu.push_back(momento3);
      if(mo[0]<mo[1] && mo[1]!=0 && mo[0]!=0){
      true_momentum_p1.push_back(mo[0]);
      true_momentum_p2.push_back(mo[1]);
      for(unsigned int i=0; i<3; i++){
        p1(i)=res[i];
      }
      for(unsigned int i=0; i<3; i++){
        p2(i)=res[i+3];
      }      
      }
      else if(mo[1]<mo[0] && mo[1]!=0 && mo[0]!= 0){
      true_momentum_p1.push_back(mo[1]);
      true_momentum_p2.push_back(mo[0]);
      for(unsigned int i=0; i<3; i++){
        p2(i)=res[i];
      }
      for(unsigned int i=0; i<3; i++){
        p1(i)=res[i+3];
      }
      }

      angle_true.push_back(p1.Angle(p2));
      angletrue=p1.Angle(p2);
      for(unsigned int i=0; i<interactions.size(); i++){
        if (unsigned int j=e.GetPhysicalProcess()==interactions[i]){
          angle_int_true[i].push_back(angletrue); 
        } 
      }      
      mo.clear(); 
      res.clear();
      nproton=0;
      mucount=0;
    }
     
   // CC0Pi2P reconstructed events
  
    if(e.CheckRecoTopology(cc0pi2p_signal_map)){
    ++cc0pi2pre;
    if(unsigned int k=e.CountRecoParticlesWithPdg(2212)==2){
      // loop over the reconstructed particles  
      for(Particle &p_reco : reco_particles){        
        if(p_reco.GetFromRecoTrack()){
          if(p_reco.GetPdgCode()==13){         
            ++mureco;
            ++mucountr;
            momento2=p_reco.GetModulusMomentum();
            //get the MC particle by hits 
            int phits = p_reco.GetMCParticleIdHits();           
            for(Particle &p_true : true_particles){
              int mcid = p_true.GetMCId();  
              if(p_true.GetPdgCode()==13){ 
                if(mcid==phits){             
                  ++musignal;
                  ++mucounts;
                  momento=p_reco.GetModulusMomentum();
                }
              }
            }
          }
          else if(p_reco.GetPdgCode()==2212 && p_reco.GetKineticEnergy()>=0.021){
            ++ptotreco;
            ++nprotonreco; 
           if(p_reco.GetModulusMomentum()>=0.250){
     
            mo.push_back(p_reco.GetModulusMomentum());
            momentum=p_reco.GetMomentum();
            for(int i=0; i<3; i++){
              res.push_back(momentum(i));
            }
             
               //get the MC particle by hits 
            
            int  phits = p_reco.GetMCParticleIdHits();           
            for(Particle &p_true : true_particles){
              int mcid = p_true.GetMCId();  
              if(p_true.GetPdgCode()==2212 && p_true.GetKineticEnergy()>=0.021 && p_true.GetModulusMomentum()>=0.250){ 
                if(mcid==phits){
                  ++nprotonsig;
                  ++ptotsig;
                  zezi.push_back(p_reco.GetModulusMomentum());         
                  momentum2=p_reco.GetMomentum();
                  for(int i=0; i<3; i++){
                    res2.push_back(momentum(i));
                  }
                } 
              }              
            }
            }
          }
          }
    } //closes the loop on reco particles
  
      reco_momentum_mu.push_back(momento2); 
      if(mo[0]<mo[1]){
        reco_momentum_p1.push_back(mo[0]);
        reco_momentum_p2.push_back(mo[1]);
        for(unsigned int i=0; i<3; i++){
          p1r(i)=res[i];
        }
        for(unsigned int i=0; i<3; i++){
          p2r(i)=res[i+3];
        }
      }
      else if(mo[1]<=mo[0]){
        reco_momentum_p1.push_back(mo[1]);
        reco_momentum_p2.push_back(mo[0]);
        for(unsigned int i=0; i<3; i++){
          p2r(i)=res[i];
        }
        for(unsigned int i=0; i<3; i++){
          p1r(i)=res[i+3];
        }
      }   
    }
    
    angle_reco.push_back(p1r.Angle(p2r));
    anglereco=p1r.Angle(p2r);
    for(unsigned int i=0; i<interactions.size(); i++){
      if (unsigned int j=e.GetPhysicalProcess()==interactions[i]){
        angle_int_reco[i].push_back(anglereco); 
      } 
    }
    
    if(nprotonsig==2 && mucounts==1){
      ++cc0pi2ps2;
      signal_momentum_mu.push_back(momento);
      if(zezi[0]<=zezi[1]){
        signal_momentum_p1.push_back(zezi[0]);                                                
        signal_momentum_p2.push_back(zezi[1]); 
        for(unsigned int i=0; i<3; i++){
          p1s(i)=res2[i];
        }
        for(unsigned int i=0; i<3; i++){
          p2s(i)=res2[i+3];
        }      
      }
      else if(zezi[1]<zezi[0]){
        signal_momentum_p1.push_back(zezi[1]);
        signal_momentum_p2.push_back(zezi[0]);
        for(unsigned int i=0; i<3; i++){
          p2s(i)=res2[i];
        }
        for(unsigned int i=0; i<3; i++){
          p1s(i)=res2[i+3];
        }
      }
      
      angle_sign.push_back(p1s.Angle(p2s));
      anglesign=p1s.Angle(p2s);
      for(unsigned int i=0; i<interactions.size(); i++){
        if (unsigned int j=e.GetPhysicalProcess()==interactions[i]){
          angle_int_sign[i].push_back(anglesign); 
        } 
      }
    }
    res.clear(); 
    nprotonreco=0;
    mucountr=0;
    zezi.clear();
    nprotonsig=0;
    mucounts=0;
    res2.clear();
    //}
    }  
  } // end of the condition on reconstructed topology map
}
} // end of the loop on the event         

// ============================================
//               Fill histograms
// ============================================

// __________ momenta distributions ___________


TH1D *h_reco_mu = new TH1D("h_reco_mu", "Muon momentum distribution", 40, 0 ,5);                                        
TH1D *h_true_mu = new TH1D("h_true_mu", "Muon momentum distribution" , 40, 0, 5);                                    
TH1D *h_true_p1 = new TH1D("h_true_p1", "p1 momentum distribution", 40, 0, 2);                                       
TH1D *h_reco_p1 = new TH1D("h_reco_p1", "p1 momentum distribution", 40, 0, 2);                                         
TH1D *h_true_p2 = new TH1D("h_true_p2", "p2 momentum distribution", 40, 0, 2);                                       
TH1D *h_reco_p2 = new TH1D("h_reco_p2", "p2 momentum distribution", 40, 0, 2);                                       
TH1D *h_sign_mu = new TH1D("h_sign_mu","Muon momentum distribution", 40,0,5);                                            
TH1D *h_sign_p1 = new TH1D("h_sign_p1", "p1 momentum distribution", 40, 0, 2);                                       
TH1D *h_sign_p2 = new TH1D("h_sign_p2", "p2 momentum distribution", 40, 0, 2);                                       
for(unsigned int i=0; i<reco_momentum_mu.size(); i++){                                         
  h_reco_mu->Fill(reco_momentum_mu[i]);                                                        
}                                                                                              
for(unsigned int i=0; i<true_momentum_mu.size(); i++){                                        
  h_true_mu->Fill(true_momentum_mu[i]);                                                       
}                                                                                             
for(unsigned int i=0; i<true_momentum_p1.size(); i++){                                         
  h_true_p1->Fill(true_momentum_p1[i]);                                                         
}                                                                                             
for(unsigned int i=0; i<reco_momentum_p1.size(); i++){                                         
  h_reco_p1->Fill(reco_momentum_p1[i]);                                                         
}                                                                                              
for(unsigned int i=0; i<true_momentum_p2.size(); i++){                                         
  h_true_p2->Fill(true_momentum_p2[i]);                                                         
}                                                                                              
for(unsigned int i=0; i<reco_momentum_p2.size(); i++){                                         
  h_reco_p2->Fill(reco_momentum_p2[i]);                                                         
}                                                                                              
for(unsigned int i=0; i<signal_momentum_mu.size(); i++){                                      
  h_sign_mu->Fill(signal_momentum_mu[i]);                                                      
}                                                                                              
for(unsigned int i=0; i<signal_momentum_p1.size(); i++){                                      
  h_sign_p1->Fill(signal_momentum_p1[i]);                                                      
}                                                                                 
for(unsigned int i=0; i<signal_momentum_p2.size(); i++){                                       
  h_sign_p2->Fill(signal_momentum_p2[i]);                                                      
}
file << "=====================================================================" << endl;
file << "enties for the histogram of true mu " << h_true_mu->GetEntries() << endl;
file << "entries for the histogram of reco mu " << h_reco_mu->GetEntries() << endl;
file << "entries for the histogram of signal mu " << h_sign_mu->GetEntries() << endl;
file << "=====================================================================" << endl;
file << "=====================================================================" << endl;
file << "entries for the histogram of true p1 " << h_true_p1->GetEntries() << endl;
file << " entries for the histogram of reco p1 " << h_reco_p1->GetEntries() << endl;
file << "entries for the histogram of signal p1 " << h_sign_p1->GetEntries() << endl; 
file << "=====================================================================" << endl;
file << "entries for the histogram of true p2 " << h_true_p1->GetEntries() << endl;
file << " entries for the histogram of reco p2 " << h_reco_p1->GetEntries() << endl;
file << "entries for the histogram of signal p2 " << h_sign_p1->GetEntries() << endl; 
file << "=====================================================================" << endl;


//________ cosine of solid angle vs proton momenta plots ____________


Int_t n1 = angle_true.size();
std:: vector <Double_t> cos1; 
for(unsigned int i=0; i<angle_true.size(); i++){
  cos1.push_back(cos(angle_true[i]));
}
TGraph* g11=new TGraph(n1, &true_momentum_p1[0], &cos1[0]);
TGraph* g12=new TGraph(n1, &true_momentum_p2[0], &cos1[0]);

Int_t n2 = angle_reco.size();
std:: vector <Double_t> cos2; 
for(unsigned int i=0; i<angle_reco.size(); i++){
  cos2.push_back(cos(angle_reco[i]));
}
TGraph* g21=new TGraph(n2, &reco_momentum_p1[0], &cos2[0]);
TGraph* g22=new TGraph(n2, &reco_momentum_p2[0], &cos2[0]);

Int_t n3 = angle_sign.size();
std:: vector <Double_t> cos3; 
for(unsigned int i=0; i<angle_sign.size(); i++){
  cos3.push_back(cos(angle_sign[i]));
}
TGraph* g31=new TGraph(n3, &signal_momentum_p1[0], &cos3[0]);
TGraph* g32=new TGraph(n3, &signal_momentum_p2[0], &cos3[0]);



// ____________ cosine of solid angle distributions _______________


TH1D *a_reco = new TH1D("a_reco", " angle distribution", 40, -1,1);                                        
TH1D *a_true = new TH1D("a_true", " angle distribution" , 40, -1, 1);                                    
TH1D *a_sign = new TH1D("a_sign", " angle distribution", 40, -1, 1);                                       

for(unsigned int i=0; i<angle_true.size(); i++){                                       
  a_true->Fill(cos(angle_true[i]));                                                      
}

for(unsigned int i=0; i<angle_reco.size(); i++){                                       
  a_reco->Fill(cos(angle_reco[i]));                                                      
}

for(unsigned int i=0; i<angle_sign.size(); i++){                                       
  a_sign->Fill(cos(angle_sign[i]));                                                      
}

file << "=====================================================================" << endl;
file << "enties for the histogram of true mu " << a_true->GetEntries() << endl;
file << "entries for the histogram of reco mu " << a_reco->GetEntries() << endl;
file << "entries for the histogram of signal mu " << a_sign->GetEntries() << endl;
file << "=====================================================================" << endl;


//____________p2 momentum vs p1 momentum________________

Int_t dim_true = true_momentum_p1.size();
TGraph *p1p2_true = new TGraph(dim_true, &true_momentum_p1[0], &true_momentum_p2[0]);

Int_t dim_reco = reco_momentum_p1.size();
TGraph *p1p2_reco = new TGraph(dim_reco, &reco_momentum_p1[0], &reco_momentum_p2[0]);

Int_t dim_sign = signal_momentum_p1.size();
TGraph *p1p2_sign = new TGraph(dim_sign, &signal_momentum_p1[0], &signal_momentum_p2[0]);

//__________ costheta12 interactions __________________

//QE
TH1D *qe_reco = new TH1D("qe_reco", " angle distribution", 40, -1,1);                                        
TH1D *qe_true = new TH1D("qe_true", " angle distribution" , 40, -1, 1);                                    
TH1D *qe_sign = new TH1D("qe_sign", " angle distribution", 40, -1, 1);                                       

for(unsigned int i=0; i<angle_int_true[0].size(); i++){                                       
  qe_true->Fill(cos(angle_int_true[0][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[0].size(); i++){                                       
  qe_reco->Fill(cos(angle_int_reco[0][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[0].size(); i++){                                       
  qe_sign->Fill(cos(angle_int_sign[0][i]));                                                      
}




// DIS
TH1D *dis_reco = new TH1D("dis_reco", " angle distribution", 40, -1,1);                                        
TH1D *dis_true = new TH1D("dis_true", " angle distribution" , 40, -1, 1);                                    
TH1D *dis_sign = new TH1D("dis_sign", " angle distribution", 40, -1, 1);                                       

for(unsigned int i=0; i<angle_int_true[1].size(); i++){                                       
  dis_true->Fill(cos(angle_int_true[1][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[1].size(); i++){                                       
  dis_reco->Fill(cos(angle_int_reco[1][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[1].size(); i++){                                       
  dis_sign->Fill(cos(angle_int_sign[1][i]));                                                      
}

// RES
TH1D *res_reco = new TH1D("res_reco", " angle distribution", 40, -1,1);                                        
TH1D *res_true = new TH1D("res_true", " angle distribution" , 40, -1, 1);                                    
TH1D *res_sign = new TH1D("res_sign", " angle distribution", 40, -1, 1);                                       


for(unsigned int i=0; i<angle_int_true[2].size(); i++){                                       
  res_true->Fill(cos(angle_int_true[2][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[2].size(); i++){                                       
  res_reco->Fill(cos(angle_int_reco[2][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[2].size(); i++){                                       
 res_sign->Fill(cos(angle_int_sign[2][i]));                                                      
}


// COH
TH1D *coh_reco = new TH1D("coh_reco", " angle distribution", 40, -1,1);                                        
TH1D *coh_true = new TH1D("coh_true", " angle distribution" , 40, -1, 1);                                    
TH1D *coh_sign = new TH1D("coh_sign", " angle distribution", 40, -1, 1);                                       


for(unsigned int i=0; i<angle_int_true[3].size(); i++){                                       
  coh_true->Fill(cos(angle_int_true[3][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[3].size(); i++){                                       
  coh_reco->Fill(cos(angle_int_reco[3][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[3].size(); i++){                                       
  coh_sign->Fill(cos(angle_int_sign[3][i]));                                                      
}


// MEC

TH1D *mec_reco = new TH1D("mec_reco", " angle distribution", 40, -1,1);                                        
TH1D *mec_true = new TH1D("mec_true", " angle distribution" , 40, -1, 1);                                    
TH1D *mec_sign = new TH1D("mec_sign", " angle distribution", 40, -1, 1);                                       


for(unsigned int i=0; i<angle_int_true[4].size(); i++){                                       
  mec_true->Fill(cos(angle_int_true[4][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[4].size(); i++){                                       
  mec_reco->Fill(cos(angle_int_reco[4][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[4].size(); i++){                                       
  mec_sign->Fill(cos(angle_int_sign[4][i]));                                                      
}



//=========================================================
//                        Drawing
//=========================================================

// cosine of solid angle in LAB frame VS protons momenta 


TCanvas *f1= new TCanvas();
g11->SetTitle("cos(\\theta_{p1p2}) \\ vs \\ p1 \\ true");
g11->SetMarkerSize(1);
g11->SetMarkerStyle(7);
g11->SetMarkerColor(kRed);
g11->GetYaxis()->SetRangeUser(-1, 1);                                                   
g11->GetXaxis()->SetRangeUser(0, 2);                                                                                                      
g11->GetYaxis()->SetTitle("Cos(\\theta_{p1p2})");                            
g11->GetXaxis()->SetTitle("p1(GeV/c)");
g11->Draw("AP");
f1->SaveAs((angle+"cos_theta12_p1_true.root").c_str());                                       
f1->SaveAs((angle+"cos_theta12_p1_true.png").c_str());                                        


TCanvas *f2= new TCanvas();
g12->SetTitle("cos(\\theta_{p1p2} \\ vs \\ p2 \\ true");
g12->SetMarkerSize(1);
g12->SetMarkerStyle(7);
g12->SetMarkerColor(kRed);
g12->GetXaxis()->SetRangeUser(0, 2);                                                   
g12->GetYaxis()->SetRangeUser(-1, 1);                                                                                
g12->GetYaxis()->SetTitle("Cos(\\theta_{p1p2})");                            
g12->GetXaxis()->SetTitle("p2(GeV/c)");
g12->Draw("AP");
f2->SaveAs((angle+"cos_theta12_p2_true.root").c_str());                                       
f2->SaveAs((angle+"cos_theta12_p2_true.png").c_str());                                        


TCanvas *f3= new TCanvas();
g21->SetTitle("cos(\\theta_{p1p2}) \\ vs \\ p1 \\ reco");
g21->SetMarkerSize(1);
g21->SetMarkerStyle(7);
g21->SetMarkerColor(kBlue);
g21->GetYaxis()->SetRangeUser(-1, 1);                                                   
g21->GetXaxis()->SetRangeUser(0, 2);                                                                                
g21->GetYaxis()->SetTitle("Cos(\\theta_{p1p2})");                            
g21->GetXaxis()->SetTitle("p1 (GeV/c)");
g21->Draw("AP");
f3->SaveAs((angle+"cos_theta12_p1_reco.root").c_str());                                       
f3->SaveAs((angle+"cos_theta12_p1_reco.png").c_str());                                        


TCanvas *f4= new TCanvas();
g22->SetTitle("cos(\\theta_{p1p2}) \\ vs \\ p2 \\ RECO");
g22->SetMarkerSize(1);
g22->SetMarkerStyle(7);
g22->SetMarkerColor(kBlue);
g22->GetYaxis()->SetRangeUser(-1, +1);                                                   
g22->GetXaxis()->SetRangeUser(0, 2);                                                                                
g22->GetYaxis()->SetTitle("Cos(\\theta_{p1p2})");                            
g22->GetXaxis()->SetTitle("p2(GeV/c)") ; 
g22->Draw("AP");
f4->SaveAs((angle+"cos_theta12_p2_reco.root").c_str());                                       
f4->SaveAs((angle+"cos_theta12_p2_reco.png").c_str());                                        


TCanvas *f5 = new TCanvas();
g31->SetTitle("cos(\\theta_{p1p2}) \\ vs \\ p1 \\ sign");
g31->SetMarkerSize(1);
g31->SetMarkerStyle(7);
g31->SetMarkerColor(kRed);
g31->GetXaxis()->SetRangeUser(0, 2);                                                   
g31->GetYaxis()->SetRangeUser(-1, 1);                                                                                
g31->GetYaxis()->SetTitle("Cos(\\theta_{p1p2})");                            
g31->GetXaxis()->SetTitle("p1 (GeV/c)");
g31->Draw("AP");
g31->SetMarkerSize(.6);
f5->SaveAs((angle+"cos_theta12_p1_sign.root").c_str());                                       
f5->SaveAs((angle+"cos_theta12_p1_sign.png").c_str());                                        


TCanvas *f6 = new TCanvas();
g32->SetTitle("cos(\\theta_{p1p2}) \\ vs \\ p2 \\ sign");
g32->SetMarkerSize(1);
g32->SetMarkerStyle(7);
g32->SetMarkerColor(kRed);
g32->GetYaxis()->SetRangeUser(-1,1);                                                   
g32->GetXaxis()->SetRangeUser(0, 2);                                                                                
g32->GetYaxis()->SetTitle("Cos(\\theta_{p1p2})");                            
g32->GetXaxis()->SetTitle("p2 (GeV/c)");
g32->Draw("AP");
f6->SaveAs((angle+"cos_theta12_p2_sign.root").c_str());                                       
f6->SaveAs((angle+"cos_theta12_p2_sign.png").c_str());                                        



// p2 vs p1 graphs


TCanvas *f7 = new TCanvas();
p1p2_true->SetTitle("p2 \\ vs \\ p1 \\ true");
p1p2_true->SetMarkerSize(1);
p1p2_true->SetMarkerStyle(7);
p1p2_true->SetMarkerColor(kBlue);
p1p2_true->GetYaxis()->SetRangeUser(0,2);                                                   
p1p2_true->GetXaxis()->SetRangeUser(0, 2);                                                                                
p1p2_true->GetYaxis()->SetTitle("p2 (GeV/c)");                            
p1p2_true->GetXaxis()->SetTitle("p1 (GeV/c)");
p1p2_true->Draw("AP");
f7->SaveAs((momenta+"p2_vs_p1_true.root").c_str());                                       
f7->SaveAs((momenta+"p2_vs_p1_true.png").c_str());                                        

TCanvas *f8 = new TCanvas();
p1p2_reco->SetTitle("p2 \\ vs \\ p1 \\ reco");
p1p2_reco->SetMarkerSize(1);
p1p2_reco->SetMarkerStyle(7);
p1p2_reco->SetMarkerColor(kBlue);
p1p2_reco->GetYaxis()->SetRangeUser(0,2);                                                   
p1p2_reco->GetXaxis()->SetRangeUser(0, 2);                                                                                
p1p2_reco->GetXaxis()->SetTitle("p1(GeV/c)");                            
p1p2_reco->GetYaxis()->SetTitle("p2 (GeV/c)");
p1p2_reco->Draw("AP");
f8->SaveAs((momenta+"p2_vs_p1_reco.root").c_str());                                       
f8->SaveAs((momenta+"p2_vs_p1_reco.png").c_str());                                        

TCanvas *f9 = new TCanvas();
p1p2_sign->SetTitle("p2 \\ vs \\ p1 \\ sign");
p1p2_sign->SetMarkerSize(1);
p1p2_sign->SetMarkerStyle(7);
p1p2_sign->SetMarkerColor(kBlue);
p1p2_sign->GetYaxis()->SetRangeUser(0,2);                                                   
p1p2_sign->GetXaxis()->SetRangeUser(0, 2);                                                                                
p1p2_sign->GetXaxis()->SetTitle("p1(GeV/c)");                            
p1p2_sign->GetYaxis()->SetTitle("p2 (GeV/c)");
p1p2_sign->Draw("AP");
f9->SaveAs((momenta+"p2_vs_p1_sign.root").c_str());                                       
f9->SaveAs((momenta+"p2_vs_p1_sign.png").c_str());                                        


// momenta distributions 

TCanvas *c1 = new TCanvas();                                                                    
gStyle->SetOptStat(0);
h_reco_mu->SetLineColorAlpha(kBlue, 0.30);
h_reco_mu->SetLineWidth(2);
h_reco_mu->SetFillColorAlpha(kBlue, 0.20);
h_reco_mu->SetFillStyle(3001);
h_true_mu->SetLineColorAlpha(kRed, 0.30);
h_true_mu->SetLineWidth(2);
h_true_mu->SetFillColorAlpha(kRed, 0.20);
h_true_mu->SetFillStyle(3001);
h_sign_mu->SetLineColorAlpha(kGreen, 0.30);
h_sign_mu->SetLineWidth(2);
h_sign_mu->SetFillColorAlpha(kGreen, 0.20);
h_sign_mu->SetFillStyle(3001);
h_true_mu->Scale(1/(h_true_mu->GetEntries()));                                                 
h_reco_mu->Scale(1/(h_reco_mu->GetEntries()));                                                 
h_sign_mu->Scale(1/(h_sign_mu->GetEntries()));                                                 
h_sign_mu->GetYaxis()->SetRangeUser(0, 0.35);                                                   
h_reco_mu->GetYaxis()->SetRangeUser(0, 0.35);                                                   
h_true_mu->GetYaxis()->SetRangeUser(0, 0.35);                                                   
h_sign_mu->GetYaxis()->SetTitle("Arbitrary units");                                                     
h_sign_mu->GetXaxis()->SetTitle("p_{\\mu} (GeV/c)");                                             
h_reco_mu->GetYaxis()->SetTitle("Arbitrary units");                                                     
h_reco_mu->GetXaxis()->SetTitle("p_{\\mu} (GeV/c)");                                             
h_true_mu->GetYaxis()->SetTitle("Arbitrary units");                                                     
h_true_mu->GetXaxis()->SetTitle("p_{\\mu} (GeV/c)");                                             
h_reco_mu->Draw("HIST");                                                                       
h_true_mu->Draw("HIST SAME");                                                                  
h_sign_mu->Draw("HIST SAME"); 
TLegend *legend = new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend->AddEntry(h_true_mu,"true","f");                                                            
legend->AddEntry(h_reco_mu, "reco", "f");
legend->AddEntry(h_sign_mu, "signal", "f");
legend->Draw();  
c1->Update();
c1->SaveAs((momenta+"Muon_Histogram_momentum.root").c_str());                                       
c1->SaveAs((momenta+"Muon_Histogram_momentum.png").c_str());                                        


TCanvas *c2 = new TCanvas();                                                                    
gStyle->SetOptStat(0);                                                                     
h_reco_p1->SetLineColorAlpha(kBlue, 0.30);
h_reco_p1->SetLineWidth(2);
h_reco_p1->SetFillColorAlpha(kBlue, 0.20);
h_reco_p1->SetFillStyle(3001);
h_true_p1->SetLineColorAlpha(kRed, 0.30);
h_true_p1->SetLineWidth(2);
h_true_p1->SetFillColorAlpha(kRed, 0.20);
h_true_p1->SetFillStyle(3001);
h_sign_p1->SetLineColorAlpha(kGreen, 0.30);
h_sign_p1->SetLineWidth(2);
h_sign_p1->SetFillColorAlpha(kGreen, 0.20);
h_sign_p1->SetFillStyle(3001);
h_true_p1->Scale(1/(h_true_p1->GetEntries()));                                                 
h_reco_p1->Scale(1/(h_reco_p1->GetEntries()));                                                
h_sign_p1->Scale(1/(h_sign_p1->GetEntries()));                                             
h_sign_p1->GetYaxis()->SetRangeUser(0, 0.25);                                                   
h_reco_p1->GetYaxis()->SetRangeUser(0, 0.25);                                                   
h_true_p1->GetYaxis()->SetRangeUser(0, 0.25);                   
h_sign_p1->GetYaxis()->SetTitle("Arbitrary units");                                       
h_sign_p1->GetXaxis()->SetTitle("p_{p_{1}} (GeV/c)");                                             
h_reco_p1->GetYaxis()->SetTitle("Arbitrary units");                                                     
h_reco_p1->GetXaxis()->SetTitle("p_{p_{1}} (GeV/c)");                                             
h_true_p1->GetYaxis()->SetTitle("Arbitrary units");                                                     
h_true_p1->GetXaxis()->SetTitle("p_{p_{1}} (GeV/c)");                                             
h_reco_p1->Draw("HIST");                                                                       
h_true_p1->Draw("HIST SAME");                                                                  
h_sign_p1->Draw("HIST SAME"); 
TLegend *legend2 = new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend2->AddEntry(h_true_p1,"true","f");                                                            
legend2->AddEntry(h_reco_p1, "reco", "f");
legend2->AddEntry(h_sign_p1, "signal", "f");
legend2->Draw();  
c2->Update();
c2->SaveAs((momenta+"p1_Histogram_momentum.root").c_str());                                       
c2->SaveAs((momenta+"p1_Histogram_momentum.png").c_str());                                        

TCanvas *c3 = new TCanvas();                                                                    
gStyle->SetOptStat(0);                                                                      
h_reco_p2->SetLineColorAlpha(kBlue, 0.30);
h_reco_p2->SetLineWidth(2);
h_reco_p2->SetFillColorAlpha(kBlue, 0.20);
h_reco_p2->SetFillStyle(3001);
h_true_p2->SetLineColorAlpha(kRed, 0.30);
h_true_p2->SetLineWidth(2);
h_true_p2->SetFillColorAlpha(kRed, 0.20);
h_true_p2->SetFillStyle(3001);
h_sign_p2->SetLineColorAlpha(kGreen, 0.30);
h_sign_p2->SetLineWidth(2);
h_sign_p2->SetFillColorAlpha(kGreen, 0.20);
h_sign_p2->SetFillStyle(3001);
h_true_p2->Scale(1/(h_true_p2->GetEntries()));                                                 
h_reco_p2->Scale(1/(h_reco_p2->GetEntries()));                                                 
h_sign_p2->Scale(1/(h_sign_p2->GetEntries()));                                                 
h_sign_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                   
h_reco_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                   
h_true_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                  
h_sign_p2->GetYaxis()->SetTitle("Arbitrary units");                                                     
h_sign_p2->GetXaxis()->SetTitle("p_{p_{2}} (GeV/c)");                                             
h_reco_p2->GetYaxis()->SetTitle("Arbitrary units");                                                     
h_reco_p2->GetXaxis()->SetTitle("p_{p_{2}} (GeV/c)");                                             
h_true_p2->GetYaxis()->SetTitle("Arbitrary units");                                                     
h_true_p2->GetXaxis()->SetTitle("p_{p_{2}} (GeV/c)");                                             
h_reco_p2->Draw("HIST");                                                                       
h_true_p2->Draw("HIST SAME");                                                                  
h_sign_p2->Draw("HIST SAME"); 
TLegend *legend3 = new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend3->AddEntry(h_true_p2,"true","f");                                                            
legend3->AddEntry(h_reco_p2, "reco", "f");
legend3->AddEntry(h_sign_p2, "signal", "f");
legend3->Draw();  
c3->Update();
c3->SaveAs((momenta+"p2_Histogram_momentum.root").c_str());                                       
c3->SaveAs((momenta+"p2_Histogram_momentum.png").c_str());                                        

// cosine of solid angle in LAB frame distribution 

TCanvas *a = new TCanvas();                                                                    
gStyle->SetOptStat(0);
a_true->SetLineColorAlpha(kBlue, 0.30);
a_true->SetLineWidth(2);
a_true->SetFillColorAlpha(kBlue, 0.20);
a_true->SetFillStyle(3001);
a_reco->SetLineColorAlpha(kRed, 0.30);
a_reco->SetLineWidth(2);
a_reco->SetFillColorAlpha(kRed, 0.20);
a_reco->SetFillStyle(3001);
a_sign->SetLineColorAlpha(kGreen, 0.30);
a_sign->SetLineWidth(2);
a_sign->SetFillColorAlpha(kGreen, 0.20);
a_sign->SetFillStyle(3001);
a_true->Scale(1/(a_true->GetEntries()));                                                 
a_reco->Scale(1/(a_reco->GetEntries()));                                                 
a_sign->Scale(1/(a_sign->GetEntries()));                                                 
a_sign->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_reco->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_true->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_sign->GetYaxis()->SetTitle("Arbitrary units");                                                     
a_sign->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             
a_reco->GetYaxis()->SetTitle("Arbitrary units");                                                     
a_reco->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             
a_true->GetYaxis()->SetTitle("Arbitrary units");                                                     
a_true->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             
a_reco->Draw("HIST");                                                                       
a_true->Draw("HIST SAME");                                                                  
a_sign->Draw("HIST SAME"); 
TLegend *legendn= new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legendn->AddEntry(a_true,"true","f");                                                            
legendn->AddEntry(a_reco, "reco", "f");
legendn->AddEntry(a_sign, "signal", "f");
legendn->Draw();  
a->Update();
a->SaveAs((angle+"Angle_Histogram.root").c_str());                                       
a->SaveAs((angle+"Angle_Histogram.png").c_str());                                        

// true events 
//
TCanvas *b = new TCanvas();                                                                    
gStyle->SetOptStat(0); 
qe_true->SetFillColorAlpha(4, 0.5);
qe_true->Scale(1/(a_true->GetEntries()));                                                 
qe_true->GetYaxis()->SetRangeUser(0, 0.15);                                                   
qe_true->GetYaxis()->SetTitle("Arbitrary units");                                                     
qe_true->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             

dis_true->SetFillColorAlpha(2, 0.5);
dis_true->Scale(1/(a_true->GetEntries()));                                                 
dis_true->GetYaxis()->SetRangeUser(0, 0.15);                                                   
dis_true->GetYaxis()->SetTitle("Arbitrary units");                                                     
dis_true->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             

res_true->SetFillColorAlpha(5, 0.5);
res_true->Scale(1/(a_true->GetEntries()));                                                 
res_true->GetYaxis()->SetRangeUser(0, 0.15);                                                   
res_true->GetYaxis()->SetTitle("Arbitrary units");                                                     
res_true->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             

coh_true->SetFillColorAlpha(6, 0.5);
//coh_true->Scale(1/(coh_true->GetEntries()));                                                 
coh_true->GetYaxis()->SetRangeUser(0, 0.15);                                                   
coh_true->GetYaxis()->SetTitle("Arbitrary units");                                                     
coh_true->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             

mec_true->SetFillColorAlpha(3, 0.5);
mec_true->Scale(1/(a_true->GetEntries()));                                                 
mec_true->GetYaxis()->SetRangeUser(0, 0.15);                                                   
mec_true->GetYaxis()->SetTitle("Arbitrary units");                                                     
mec_true->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             


THStack *hs_true=new THStack("hs_true", "cos(\\theta_{p1p2}) true" ); 
hs_true->Add(qe_true);
hs_true->Add(res_true);
hs_true->Add(coh_true);
hs_true->Add(dis_true);
hs_true->Add(mec_true);
hs_true->Write();
hs_true->Draw("HIST");
hs_true->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");
hs_true->GetYaxis()->SetTitle("Arbitrary units");
hs_true->GetYaxis()->SetRangeUser(0, 0.15);

TLegend *legend_true= new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend_true->AddEntry(qe_true,"QE","f");                                                            
legend_true->AddEntry(res_true, "RES", "f");
legend_true->AddEntry(coh_true, "COH", "f");
legend_true->AddEntry(dis_true, "DIS" , "f");
legend_true->AddEntry(mec_true, "MEC", "f");
legend_true->Draw("SAME");  
b->Update();
b->SaveAs((angle+"Angle_p1p2_interactions_true.root").c_str());                                       
b->SaveAs((angle+"Angle_p1p2_interactions_true.png").c_str());                                        

// reco events 

TCanvas *c = new TCanvas();
gStyle->SetOptStat(0);
qe_reco->SetLineWidth(2);
qe_reco->SetFillColorAlpha(4, 0.5);
qe_reco->Scale(1/(a_reco->GetEntries()));                                                 
qe_reco->GetYaxis()->SetRangeUser(0, 0.15);                                                   
qe_reco->GetYaxis()->SetTitle("Arbitrary units");                                                     
qe_reco->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             

dis_reco->SetFillColorAlpha(2, 0.5);
dis_reco->Scale(1/(a_reco->GetEntries()));                                                 
dis_reco->GetYaxis()->SetRangeUser(0, 0.15);                                                   
dis_reco->GetYaxis()->SetTitle("Arbitrary units");                                                     
dis_reco->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             

res_reco->SetFillColorAlpha(5, 0.5);
res_reco->Scale(1/(a_reco->GetEntries()));                                                 
res_reco->GetYaxis()->SetRangeUser(0, 0.15);                                                   
res_reco->GetYaxis()->SetTitle("Arbitrary units");                                                     
res_reco->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             

coh_reco->SetFillColorAlpha(6, 0.5);
//coh_reco->Scale(1/(coh_reco->GetEntries()));                                                 
coh_reco->GetYaxis()->SetRangeUser(0, 0.15);                                                   
coh_reco->GetYaxis()->SetTitle("Arbitrary units");                                                     
coh_reco->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             

mec_reco->SetFillColorAlpha(3, 0.5);
mec_reco->Scale(1/(a_reco->GetEntries()));                                                 
mec_reco->GetYaxis()->SetRangeUser(0, 0.15);                                                   
mec_reco->GetYaxis()->SetTitle("Arbitrary units");                                                     
mec_reco->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             




THStack *hs_reco=new THStack("hs_reco","cos(\\theta_{p1p2}) reco" );
hs_reco->Add(qe_reco);
hs_reco->Add(res_reco);
hs_reco->Add(coh_reco);
hs_reco->Add(dis_reco);
hs_reco->Add(mec_reco);
hs_reco->Draw("HIST"); 
hs_reco->GetYaxis()->SetRangeUser(0, 0.15);
hs_reco->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");
hs_reco->GetYaxis()->SetTitle("Arbitrary units");


TLegend *legend_reco= new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend_reco->AddEntry(qe_reco,"QE","f");                                                            
legend_reco->AddEntry(res_reco, "RES", "f");
legend_reco->AddEntry(coh_reco, "COH", "f");
legend_reco->AddEntry(dis_reco, "DIS" , "f");
legend_reco->AddEntry(mec_reco, "MEC", "f");
legend_reco->Draw("SAME");  
c->Update();
c->SaveAs((angle+"Angle_p1p2_interactions_reco.root").c_str());                                       
c->SaveAs((angle+"Angle_p1p2_interactions_reco.png").c_str());                                        


// reco events 

TCanvas *d = new TCanvas();
gStyle->SetOptStat(0);
qe_sign->SetLineWidth(2);
qe_sign->SetFillColorAlpha(4, 0.5);
qe_sign->Scale(1/(a_sign->GetEntries()));                                                 
qe_sign->GetYaxis()->SetRangeUser(0, 0.15);                                                   
qe_sign->GetYaxis()->SetTitle("Arbitrary units");                                                     
qe_sign->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             
qe_sign->Draw("HIST");

dis_sign->SetFillColorAlpha(2, 0.5);
dis_sign->Scale(1/(a_sign->GetEntries()));                                                 
dis_sign->GetYaxis()->SetRangeUser(0, 0.15);                                                   
dis_sign->GetYaxis()->SetTitle("Arbitrary units");                                                     
dis_sign->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             
dis_sign->Draw("HIST SAME");

res_sign->SetFillColorAlpha(5, 0.5);
res_sign->Scale(1/(a_sign->GetEntries()));                                                 
res_sign->GetYaxis()->SetRangeUser(0, 0.15);                                                   
res_sign->GetYaxis()->SetTitle("Arbitrary units");                                                     
res_sign->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             
res_sign->Draw("HIST SAME");

coh_sign->SetFillColorAlpha(6, 0.5);
//coh_sign->Scale(1/(coh_sign->GetEntries()));                                                 
coh_sign->GetYaxis()->SetRangeUser(0, 0.15);                                                   
coh_sign->GetYaxis()->SetTitle("Arbitrary units");                                                     
coh_sign->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             



mec_sign->SetFillColorAlpha(3, 0.5);
mec_sign->Scale(1/(a_sign->GetEntries()));                                                 
mec_sign->GetYaxis()->SetRangeUser(0, 0.15);                                                   
mec_sign->GetYaxis()->SetTitle("Arbitrary units");                                                     
mec_sign->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");                                             
mec_sign->Draw("HIST SAME");



THStack *hs_sign=new THStack("hs_sign","cos(\\theta_{p1p2}) signal");
hs_sign->Add(qe_sign);
hs_sign->Add(res_sign);
hs_sign->Add(coh_sign);
hs_sign->Add(dis_sign);
hs_sign->Add(mec_sign);
hs_sign->Draw("HIST");
hs_sign->GetXaxis()->SetTitle("cos(\\theta_{p1p2})");
hs_sign->GetYaxis()->SetTitle("Arbitrary units");
hs_sign->GetYaxis()->SetRangeUser(0, 0.15);

TLegend *legend_sign= new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend_sign->AddEntry(qe_sign,"QE","f");                                                            
legend_sign->AddEntry(res_sign, "RES", "f");
legend_sign->AddEntry(coh_sign, "COH", "f");
legend_sign->AddEntry(dis_sign, "DIS" , "f");
legend_sign->AddEntry(mec_sign, "MEC", "f");
legend_sign->Draw("SAME");  
d->Update();
d->SaveAs((angle+"Angle_p1p2_interactions_sign.root").c_str());                                       
d->SaveAs((angle+"Angle_p1p2_interactions_sign.png").c_str());                                        





file.close();

return 0;
}
