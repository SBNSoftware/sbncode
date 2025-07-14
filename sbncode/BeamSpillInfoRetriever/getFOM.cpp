#include "getFOM.h"
#include <iostream>
#include <math.h>
#include "TH1D.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include "TCanvas.h"


using namespace std;


namespace sbn {
  
  bool onePlot = true;

  float getFOM(BNBSpillInfo & spill )
  //float getFOM2( const ub_BeamHeader& bh, const std::vector<ub_BeamData>& bd, const bnb::bnbAutoTune settings, bool useAutoTune)
  {
    double fom=0;
    
    
    
    
    
    // ========================================= I believe all of this can be replaced with the BNBSpillInfo which I fill already ============================================
    
    
    
    
    


    double hp875_offset= spill.HP875Offset;
    double vp875_offset= spill.VP875Offset;
    double hptg1_offset= spill.HPTG1Offset;
    double vptg1_offset= spill.VPTG1Offset;
    double hptg2_offset= spill.HPTG2Offset;
    double vptg2_offset= spill.VPTG2Offset;
    
    
    //int useHTG = bnb::kTG1;
    //int useVTG = bnb::kTG1;
    
    int useHTG = 0;
    int useVTG = 0;
    
    
    //if(useAutoTune){
    // bnb::bnbOffsets data = settings.getData();
    // hp875_offset = data.hp875;
    // vp875_offset = data.vp875;
    // hptg1_offset = data.hptg1;
    // vptg1_offset = data.vptg1;
    // hptg2_offset = data.hptg2;
    // vptg2_offset = data.vptg2;
    //useHTG = settings.getHTG();
    //useVTG = settings.getVTG();
    //}
    
    
    
    //double hptg2_offset= +0.79;
    //double vptg2_offset= +1.02;
    double hp875_zpos= 202.116104;
    double vp875_zpos= 202.3193205;
    double hptg1_zpos= 204.833267;
    double vptg1_zpos= 204.629608;
    double hptg2_zpos= 205.240662;
    double vptg2_zpos= 205.036835;
    double target_center_zpos= 206.870895;
    double p875x[]={0.431857, 0.158077, 0.00303551};
    double p875y[]={0.279128, 0.337048, 0};
    double p876x[]={0.166172, 0.30999, -0.00630299};
    double p876y[]={0.13425, 0.580862, 0};
    
    
    std::vector<double> tor860;
    std::vector<double> tor875;
    std::vector<double> hp875;
    std::vector<double> vp875;
    std::vector<double> hptg1;
    std::vector<double> vptg1;
    std::vector<double> hptg2;
    std::vector<double> vptg2;
    
    std::vector<double> mw875(spill.M875BB.begin(), spill.M875BB.end());
    std::vector<double> mw876(spill.M876BB.begin(), spill.M876BB.end());
    std::vector<double> mwtgt(spill.MMBTBB.begin(), spill.MMBTBB.end());
    
    //std::cout << "MWR 875 Size: " << mw875.size() << "  | MWR876 Size:  " << mw876.size() << "  |  MWTGT Size: " << mwtgt.size() << std::endl;  
    //for (size_t m =0; m< mw875.size() ; m++){
    //std::cout << "MWR 875: " << mw875[m] << "  | MWR876: " << mw876[m] << "  |  MWTGT: " << mwtgt[m] << std::endl;  
    //}
    
    /*
    if (onePlot){
      TH1D* MWRx =new TH1D("MWRx","MWRx",48,-12,12);
      TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

      TH1D* MWRy =new TH1D("MWRy","MWRy",48,-12,12);
      TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);

      
      for(int w =0 ; w< 48; w++){
	MWRx->SetBinContent(w+1, -mw875[w]);
	MWRx->SetBinContent(w+1, -mw875[w+48]);
      }

      c1->cd();
      MWRx->Draw("HIST");
      c1->Draw();

      c2->cd();
      MWRy->Draw("HIST");
      c2->Draw();

      onePlot =false;
    }
    */
    
    //std::vector<double> mw875 = spill.M875BB;
    //std::vector<double> mw876 = spill.M876BB;
    //std::vector<double> mwtgt = spill.MMBTBB;
    
    tor860.push_back(spill.TOR860);
    tor875.push_back(spill.TOR875);
    hp875.push_back(spill.HP875);
    vp875.push_back(spill.VP875);
    hptg1.push_back(spill.HPTG1);
    vptg1.push_back(spill.VPTG1);
    hptg2.push_back(spill.HPTG2);
    vptg2.push_back(spill.VPTG2);
    


    //std::cout << spill.HP875 << "," <<spill.HP875Offset << "," << spill.HPTG2 << "," << spill.HPTG2Offset << "," << spill.VP873 << "," << spill.VP873Offset << ","  << spill.VP875 << "," << spill.VP875Offset  << "," <<  spill.VPTG1 << "," << spill.VPTG1Offset << "," << spill.VPTG2 << "," << spill.VPTG2Offset << " | Figure of Merit:  " << spill.FOM  << std::endl;


    /*
      for(auto& bdata : bd) {        // get toroid, BPM, and multiwire data
      if(bdata.getDeviceName().find("E:TOR860") != std::string::npos) {        // get E:TOR860 reading
      tor860 = bdata.getData();
      continue;
      } else if(bdata.getDeviceName().find("E:TOR875") != std::string::npos) {        // get E:TOR875 reading
      tor875 = bdata.getData();
      continue;
      }        else if(bdata.getDeviceName().find("E:HP875") != std::string::npos) {        // get E:HP875 reading
      hp875 = bdata.getData();
      continue;
      }        else if(bdata.getDeviceName().find("E:VP875") != std::string::npos) {        // get E:VP875 reading
      vp875 = bdata.getData();
      continue;
      }        else if(bdata.getDeviceName().find("E:HPTG1") != std::string::npos) {        // get E:HPTG1 reading
      hptg1 = bdata.getData();
      continue;
      } else if(bdata.getDeviceName().find("E:VPTG1") != std::string::npos) {        // get E:VPTG1 reading
      vptg1 = bdata.getData();
      continue;
      }        else if(bdata.getDeviceName().find("E:HPTG2") != std::string::npos) {        // get E:HPTG1 reading
      hptg2 = bdata.getData();
      continue;
      } else if(bdata.getDeviceName().find("E:VPTG2") != std::string::npos) {        // get E:VPTG1 reading
      vptg2 = bdata.getData();
      continue;
      } else if(bdata.getDeviceName().find("E:M875BB") != std::string::npos) {        // get horizontal and vertical target multiwire reading
      mw875.resize(96);
      for (int ii=0;ii<96;ii++)
      mw875[ii]=bdata.getData()[ii];
      continue;
      } else if(bdata.getDeviceName().find("E:M876BB") != std::string::npos) {        // get horizontal and vertical target multiwire reading
      mw876.resize(96);
      for (int ii=0;ii<96;ii++)
      mw876[ii]=bdata.getData()[ii];
      continue;
      } else if(bdata.getDeviceName().find("E:MMBTBB") != std::string::npos) {        // get horizontal and vertical target multiwire reading
      mwtgt.resize(96);
      for (int ii=0;ii<96;ii++)
      mwtgt[ii]=bdata.getData()[ii];
      continue;
      }
      }
    */
    
    double tor;
    if (tor860.size()>0)
      tor=tor860[0];
    else if (tor875.size()>0)
      tor=tor875[0];
    else
      return -1;
    
    
    
    
    // ==================================================================== Replace with BNBSpillInfo ====================================================================================
    
    
    
    
    
    //when creating ntuples for pot counting script the variables are filled with -999
    //this could create a difference when passing events with FOM>1 since
    //events with missing BPM data would get FOM=2, while events with BPM set to -999 will get FOM=0
    //bad or missing MWR data gets FOM 4 in either case
    if (hptg2.size()==0) hptg2.push_back(-999);
    if (hptg1.size()==0) hptg1.push_back(-999);
    if (hp875.size()==0) hp875.push_back(-999);
    if (vptg2.size()==0) vptg2.push_back(-999);
    if (vptg1.size()==0) vptg1.push_back(-999);
    if (vp875.size()==0) vp875.push_back(-999);
    double horang,horpos;
    //if(useHTG == bnb::kTG1){
    if(useHTG == 1){
      if (hptg1.size()>0 && hp875.size()>0) {
	horang=((hptg1[0]-hptg1_offset)-(hp875[0]-hp875_offset))/(hptg1_zpos-hp875_zpos);
	horpos=(hp875[0]-hp875_offset)+horang*(target_center_zpos-hp875_zpos);
      } else if (hptg2.size()>0 && hp875.size()>0) {
	horang=((hptg2[0]-hptg2_offset)-(hp875[0]-hp875_offset))/(hptg2_zpos-hp875_zpos);
	horpos=(hp875[0]-hp875_offset)+horang*(target_center_zpos-hp875_zpos);
      } else {
	//missing horizontal BPM data
	return 2;
      }
    }
    //else if(useHTG == bnb::kTG2){
    else if(useHTG == 0){
      if (hptg2.size()>0 && hp875.size()>0) {
	horang=((hptg2[0]-hptg2_offset)-(hp875[0]-hp875_offset))/(hptg2_zpos-hp875_zpos);
	horpos=(hp875[0]-hp875_offset)+horang*(target_center_zpos-hp875_zpos);
      } else if (hptg1.size()>0 && hp875.size()>0) {
	horang=((hptg1[0]-hptg1_offset)-(hp875[0]-hp875_offset))/(hptg1_zpos-hp875_zpos);
	horpos=(hp875[0]-hp875_offset)+horang*(target_center_zpos-hp875_zpos);
      } else {
	//missing horizontal BPM data
	return 2;
      }
    }
    else { return 2; }
    double verang,verpos;
    if(useVTG == 1){
      //if(useVTG == bnb::kTG1){
      if (vptg1.size()>0 && vp875.size()>0) {
	verang=((vptg1[0]-vptg1_offset)-(vp875[0]-vp875_offset))/(vptg1_zpos-vp875_zpos);
	verpos=(vp875[0]-vp875_offset)+verang*(target_center_zpos-vp875_zpos);
      } else if (vptg2.size()>0 && vp875.size()>0) {
	verang=((vptg2[0]-vptg2_offset)-(vp875[0]-vp875_offset))/(vptg2_zpos-vp875_zpos);
	verpos=(vp875[0]-vp875_offset)+verang*(target_center_zpos-vp875_zpos);
      } else {
	//missing vertical BPM data
	return 3;
      }
    }
    else if(useVTG == 0){
      //else if(useVTG == bnb::kTG2){
      if (vptg2.size()>0 && vp875.size()>0) {
	verang=((vptg2[0]-vptg2_offset)-(vp875[0]-vp875_offset))/(vptg2_zpos-vp875_zpos);
	verpos=(vp875[0]-vp875_offset)+verang*(target_center_zpos-vp875_zpos);
      } else if (vptg1.size()>0 && vp875.size()>0) {
	verang=((vptg1[0]-vptg1_offset)-(vp875[0]-vp875_offset))/(vptg1_zpos-vp875_zpos);
	verpos=(vp875[0]-vp875_offset)+verang*(target_center_zpos-vp875_zpos);
      } else {
	//missing vertical BPM data
	return 3;
      }
    }
    else { return 3; }
    horang=atan(horang);
    verang=atan(verang);
    double xx,yy,sx,sy,chi2x,chi2y;
    double tgtsx, tgtsy;
    bool good_tgt=false;
    bool good_876=false;
    bool good_875=false;
    if (mwtgt.size()>0) {
      //if (bh.getSeconds()<1463677200 && mwtgt.size()>0) {
      //target multiwire started failing after May 20 2016 (1463677200)
      //skip this if data from when tgt multiwire was already dead
      processBNBprofile(&mwtgt[0], xx, sx,chi2x);
      processBNBprofile(&mwtgt[48], yy, sy, chi2y);
      if (sx>0.5 && sx<10 && sy>0.3 && sy<10 && chi2x<20 && chi2y<20 && mwtgt.size()>0) {
	tgtsx=sx;
	tgtsy=sy;
	good_tgt=true;
      }
    }
    if (!good_tgt && mw876.size()>0) {
      processBNBprofile(&mw876[0], xx,sx,chi2x);
      processBNBprofile(&mw876[48], yy,sy,chi2y);
      double tgtsx876=p876x[0]+p876x[1]*sx+p876x[2]*sx*sx;
      double tgtsy876=p876y[0]+p876y[1]*sy+p876y[2]*sy*sy;
      if (tgtsx876>0.5 && tgtsx876<10 && tgtsy876>0.3 && tgtsy876<10 && chi2x<20 && chi2y<20) {
	tgtsx=tgtsx876;
	tgtsy=tgtsy876;
	good_876=true;
      }
    }
    if (!good_tgt && !good_876 && mw875.size()>0){
      processBNBprofile(&mw875[0], xx,sx,chi2x);
      processBNBprofile(&mw875[48], yy,sy,chi2y);
      double tgtsx875=p875x[0]+p875x[1]*sx+p875x[2]*sx*sx;
      double tgtsy875=p875y[0]+p875y[1]*sy+p875y[2]*sy*sy;
      if (tgtsx875>0.5 && tgtsx875<10 && tgtsy875>0.3 && tgtsy875<10 && chi2x<20 && chi2y<20) {
	tgtsx=tgtsx875;
	tgtsy=tgtsy875;
	good_875=true;
      }
    }
    if (!good_tgt && !good_876 && !good_875) {
      //failed getting  multiwire data
      return 4;
    }
    //std::cout << "FOM before power: " << sbn::calcFOM(horpos,horang,verpos,verang,tor,tgtsx,tgtsy)  << ", " << horpos << ", " << horang << ", " << verpos << ", " << verang << ", " << tor << ", " << tgtsx << ", " << tgtsy<< std::endl;
    
    //std::cout << "Tgtsx: " << tgtsx << " | Tgtsy: " << tgtsy << std::endl;
    fom=1-pow(10,sbn::calcFOM(horpos,horang,verpos,verang,tor,tgtsx,tgtsy));
    //std::cout << "FOM: " <<fom << std::endl;
    return fom;
  }
  
  /*
    bmd::autoTunes bmd::cacheAutoTuneHistory()
    {
    bmd::autoTunes history;
    bnb::bnbAutoTune item;
    // run1-3 parameters
    item.setEntry(202435254, 132, 313, bnb::bnbOffsets{-3.4, 1.48, 0.46, 0.39, -1., -1.}, bnb::kTG1, bnb::kTG1);
    history.push_back(item);
    // run4a (misaligned beam)
    item.setEntry(203988466, 174, 315, bnb::bnbOffsets{-3.4, 1.48, 0.46, -1., -1., 1.84}, bnb::kTG1, bnb::kTG2);
    history.push_back(item);
    item.setEntry(204211306, 174, 316, bnb::bnbOffsets{-3.4, 1.48, -1., -1., 1.07, 1.84}, bnb::kTG2, bnb::kTG2);
    history.push_back(item);
    // run4b-4d parameters (same as run1-3)
    item.setEntry(215107286, 410, 313, bnb::bnbOffsets{-3.4, 1.48, 0.46, 0.39, -1., -1.}, bnb::kTG1, bnb::kTG1);
    history.push_back(item);
    // item.setEntry(232368446, 268, 316, bnb::bnbOffsets{-3.4, 1.48, -1., -1., 1.07, 1.84}, bnb::kTG2, bnb::kTG2);
  // history.push_back(item);
  // item.setEntry(232368986, 268, 313, bnb::bnbOffsets{-3.4, 1.48, 0.46, 0.39, -1., -1.}, bnb::kTG1, bnb::kTG1);
  // history.push_back(item);
  // before run 5
  item.setEntry(247063273, 276, 316, bnb::bnbOffsets{-3.4, 1.48, -1., -1., 1.07, 1.84}, bnb::kTG2, bnb::kTG2);
  history.push_back(item);
  // run5a parameters used on Nov 4th
  item.setEntry(247501553, 130, 317, bnb::bnbOffsets{-3.4, 1.48, -1., -1., 1.07, -3.15}, bnb::kTG2, bnb::kTG2);
  history.push_back(item);
  // parameters being tested in the autotune from Nov5th - Nov 10th (not used because not parameters not reflected in BPM readings)
  // item.setEntry(247573532, 305, 318, bnb::bnbOffsets{-3.4, 1.48, -1., -1., 1.07, -0.03}, bnb::kTG2, bnb::kTG2);
  // history.push_back(item);
  // item.setEntry(247576952, 304, 319, bnb::bnbOffsets{-3.4, 1.48, -1., -0.03, 1.07, -1.}, bnb::kTG2, bnb::kTG1);
  // history.push_back(item);
  // item.setEntry(247577492, 304, 320, bnb::bnbOffsets{-3.4, 1.48, -1., -3.15, 1.07, -1.}, bnb::kTG2, bnb::kTG1);
  // history.push_back(item);
  // item.setEntry(247587752, 304, 322, bnb::bnbOffsets{-3.4, 1.48, -1., -1., 1.07, 0.}, bnb::kTG2, bnb::kTG2);
  // history.push_back(item);
  // item.setEntry(248020832, 304, 320, bnb::bnbOffsets{-3.4, 1.48, -1., -3.15, 1.07, -1.}, bnb::kTG2, bnb::kTG1);
  // history.push_back(item);
  // item.setEntry(248458772, 304, 322, bnb::bnbOffsets{-3.4, 1.48, -1., -1., 1.07, 0.}, bnb::kTG2, bnb::kTG2);
  // history.push_back(item);
  // later run5 data (after target scan)
  item.setEntry(248792844, 573, 323, bnb::bnbOffsets{-1.20, 1.30, -1., -1., -0.4, 0.1}, bnb::kTG2, bnb::kTG2);
  history.push_back(item);
  return history;
}


bnb::bnbAutoTune bmd::getSettings(const bmd::autoTunes& history, const ub_BeamHeader& bh)
{
  if(*(history.end()-1) <= bh) return *(history.end()-1);
  if(*(history.begin()) >= bh) return *(history.begin());
  for(auto it = history.begin(); it != history.end()-1; ++it){
    auto curr = *it;
    auto following = *(it+1);
    if((curr <= bh) && (following > bh)) return curr;
  }
  bnb::bnbAutoTune ret = bnb::bnbAutoTune();
  return ret;
}

bnb::bnbAutoTune bmd::getSettings(const bmd::autoTunes& history, const raw::BeamInfo& bi)
{
  if(*(history.end()-1) <= bi) return *(history.end()-1);
  if(*(history.begin()) >= bi) return *(history.begin());
  for(auto it = history.begin(); it != history.end()-1; ++it){
    auto curr = *it;
    auto following = *(it+1);
    if((curr <= bi) && (following > bi)) return curr;
  }
  bnb::bnbAutoTune ret = bnb::bnbAutoTune();
  return ret;
}

bnb::bnbAutoTune bmd::getSettings(const bmd::autoTunes& history, const uint64_t utctstamp)
{
  if(*(history.end()-1) <= utctstamp) return *(history.end()-1);
  if(*(history.begin()) >= utctstamp) return *(history.begin());
  for(auto it = history.begin(); it != history.end()-1; ++it){
    auto curr = *it;
    auto following = *(it+1);
    if((curr <= utctstamp) && (following > utctstamp)) return curr;
  }
  bnb::bnbAutoTune ret = bnb::bnbAutoTune();
  return ret;
}
*/

    void processBNBprofile(const double* mwdata, double &x, double& sx, double& chi2)
  {
    /*
      function takes multiwire data as input (48 channels horizontal or vertical,
      finds the min and max
      finds the first and last bin where amplitude is greater than 20%
      fits the peak between first and last bin with gaussian (assuming 2% relative errors)
      returns:
      x  - mean position from the fit
      sx - sigma x from the fit
      chi2 - chi2/NDF for the gaussian fit
    */
    double minx=0;
    double maxx=0;
    for (unsigned int i=0;i<48;i++) {
      minx=std::min(-mwdata[i],minx);
      maxx=std::max(-mwdata[i],maxx);
    }
    int first_x=-1;
    int last_x=-1;
    static int entry = 1;
    TH1D* hProf=new TH1D("hProf"+TString(std::to_string(entry)),"",48,-12,12);
    for (unsigned int i=0;i<48;i++) {
      hProf->SetBinContent(i+1,-mwdata[i]-minx);
      if (-mwdata[i]-minx    > (maxx-minx)*0.2 && first_x==-1) first_x=i;
      if (-mwdata[i]-minx    > (maxx-minx)*0.2)                last_x=i+1;
      hProf->SetBinError(i+1,(maxx-minx)*0.02);
    }
    if (hProf->GetSumOfWeights()>0) {
      hProf->Fit("gaus","Q0","",-12+first_x*0.5,-12+last_x*0.5);
      x   =hProf->GetFunction("gaus")->GetParameter(1);
      sx  =hProf->GetFunction("gaus")->GetParameter(2);
      //std::cout << "sx value: " << sx << std::endl;
      // if ( sx > 50.){
      // 	for (size_t m =0; m< 48 ; m++){
      // 	  std::cout << "Bad MWData: " << mwdata[m]  << std::endl;  
      // 	}
      //}

      chi2=hProf->GetFunction("gaus")->GetChisquare()/hProf->GetFunction("gaus")->GetNDF();
    } else {
      x=99999;
      sx=99999;
      chi2=99999;
    }
    entry += 1;
  }
  
  
  double calcFOM(double horpos, double horang, double verpos, double verang, double ppp, double tgtsx, double tgtsy)
  {
    ppp = ppp / pow(10,12);


    //code from MiniBooNE AnalysisFramework with the addition of scaling the beam profile to match tgtsx, tgtsy
    //form DQ_BeamLine_twiss_init.F
    double bx  =  4.68;
    double ax  =  0.0389;
    double gx  = (1+ax*ax)/bx;
    double nx  =  0.0958;
    double npx = -0.0286;
    double by  = 59.12;
    double ay  =  2.4159;
    double gy  = (1+ay*ay)/by;
    double ny  =  0.4577;
    double npy = -0.0271;
    //from DQ_BeamLine_make_tgt_fom2.F
    double ex = 0.1775E-06 + 0.1827E-07*ppp;
    double ey = 0.1382E-06 + 0.2608E-08*ppp;
    double dp = 0.4485E-03 + 0.6100E-04*ppp;
    double tex = ex;
    double tey = ey;
    double tdp = dp;
    //from DQ_BeamLine_beam_init.F
    double sigma1[6][6]={{0}};
    double centroid1[6]={0};
    centroid1[0] = horpos;
    centroid1[1] = horang;
    centroid1[2] = verpos;
    centroid1[3] = verang;
    centroid1[4] = 0.0;
    centroid1[5] = 0.0;
    sigma1[5][5] =  tdp*tdp;
    sigma1[0][0] =  tex*bx+ nx*nx *tdp*tdp;
    sigma1[0][1] = -tex*ax+ nx*npx*tdp*tdp;
    sigma1[1][1] =  tex*gx+npx*npx*tdp*tdp;
    sigma1[0][5] =  nx*tdp*tdp;
    sigma1[1][5] =  npx*tdp*tdp;
    sigma1[1][0] =  sigma1[0][1];
    sigma1[5][0] =  sigma1[0][5];
    sigma1[5][1] =  sigma1[1][5];
    sigma1[2][2] =  tey*by+ny*ny*tdp*tdp;
    sigma1[2][3] = -tey*ay+ny*npy*tdp*tdp;
    sigma1[3][3] =  tey*gy+npy*npy*tdp*tdp;
    sigma1[2][5] =  ny*tdp*tdp;
    sigma1[3][5] =  npy*tdp*tdp;
    sigma1[3][2] =  sigma1[2][3];
    sigma1[5][2] =  sigma1[2][5];
    sigma1[5][3] =  sigma1[3][5];

    

    for (int fi =0; fi <6; fi++){
      for(int la =0; la<6; la++){
	//std::cout << "sigma1[" << fi << "][" << la << "] = " << sigma1[fi][la] << "  | ppp: " << ppp << std::endl;
      }
    }

    double begtocnt[6][6]={
      { 0.65954,  0.43311,  0.00321,  0.10786, 0.00000,  1.97230},
      { 0.13047,  1.60192,  0.00034,  0.00512, 0.00000,  1.96723},
      {-0.00287, -0.03677, -0.35277, -4.68056, 0.00000,  0.68525},
      {-0.00089, -0.00430, -0.17722, -5.18616, 0.00000,  0.32300},
      {-0.00104,  0.00232, -0.00001, -0.00224, 1.00000, -0.00450},
      { 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  1.00000}
    };
    double cnttoups[6][6]={{0}};
    double cnttodns[6][6]={{0}};
    double identity[6][6]={{0}};
    double begtoups[6][6]={{0}};
    double begtodns[6][6]={{0}};
    for (int i=0;i<6;i++) {
      for (int j=0;j<6;j++) {
	if (i==j) {
	  cnttoups[i][j] = 1.0;
	  cnttodns[i][j] = 1.0;
	  identity[i][j] = 1.0;
	} else {
	  cnttoups[i][j] = 0.0;
	  cnttodns[i][j] = 0.0;
	  identity[i][j] = 0.0;
	}
      }
    }
    cnttoups[0][1] = -0.35710;
    cnttoups[2][3] = -0.35710;
    cnttodns[0][1] = +0.35710;
    cnttodns[2][3] = +0.35710;
    for (int i=0;i<6;i++) {
      for (int j=0;j<6;j++) {
	for (int k=0;k<6;k++) {
	  begtoups[i][k] = begtoups[i][k] + cnttoups[i][j]*begtocnt[j][k];
	  begtodns[i][k] = begtodns[i][k] + cnttodns[i][j]*begtocnt[j][k];
	}
      }
    }
    //swim to upstream of target
    double cx, cy, sx, sy, rho;
    sbn::swimBNB(centroid1,sigma1,
		 cnttoups, begtoups,
		 cx, cy, sx, sy, rho);
    double scalex=tgtsx/sx;
    double scaley=tgtsy/sy;
    //std::cout << "Scalex: " << scalex << "  | Scaley: " << scaley << " | sx: " << sx << " | sy: " << sy << " | tgtsx: " << tgtsx << " | tgtsy: " << tgtsy << std::endl;
    double fom_a=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    //swim to center of target
    sbn::swimBNB(centroid1,sigma1,
		 identity, begtocnt,
		 cx, cy, sx, sy, rho);
    double fom_b=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    //std::cout << "After Scalex: " << scalex << "  | Scaley: " << scaley << " | sx: " << sx << " | sy: " << sy << " | tgtsx: " << tgtsx << " | tgtsy: " << tgtsy  << " | rho: " << rho << std::endl;
    //swim to downstream of target
    sbn::swimBNB(centroid1,sigma1,
		 cnttodns, begtodns,
		 cx, cy, sx, sy, rho);
    double fom_c=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    // add a guard for double precision
    //std::cout << "Foma,b,c: " << fom_a << ", " << fom_b << ", " << fom_c << std::endl;
    
    if(fom_a <= -10000. || fom_b <= -10000. || fom_c <= -10000) return -10000.;
    double fom2=fom_a*0.6347 +
      fom_b*0.2812 +
      fom_c*0.0841;
    return fom2;
  }
  
  void swimBNB(const double centroid1[6], const double sigma1[6][6],
	       const double xferc[6][6], const double xfers[6][6],
	       double &cx, double& cy, double &sx, double &sy, double &rho)
  {
    //centroid
    double centroid2[6]={0};
    for (int i=0;i<6;i++) {
      for (int j=0;j<6;j++) {
	centroid2[i] = centroid2[i] + xferc[i][j]*centroid1[j];
      }
    }
    cx = centroid2[0];
    cy = centroid2[2];
    //sigma
    double sigma2[6][6]={{0}};
    for (int i = 0; i<6;i++) {
      for (int j = 0;j<6;j++) {
	for (int k = 0;k<6;k++) {
	  for (int m = 0;m<6;m++) {
	    sigma2[i][m] = sigma2[i][m] + xfers[i][j]*sigma1[j][k]*xfers[m][k];
	    //std::cout << "sigma2: " << sigma2[i][m] << " = " << sigma2[i][m]<< " + " << xfers[i][j] << " * " << sigma1[j][k] << " * " << xfers[m][k] << std::endl; 
	  }
	}
      }
    }
    //get beam sigma
    sx  = sqrt(sigma2[0][0])*1000.0;
    sy  = sqrt(sigma2[2][2])*1000.0;
    rho = sigma2[0][2]/sqrt(sigma2[0][0]*sigma2[2][2]);
    //cout<<"Swim "<<sx<<"\t"<<sy<<"\t"<<rho<<endl;
  }
  
  
  
  double func_intbivar(const double cx, const double cy, const double sx, const double sy, const double rho )
  {
    //std::cout << cx << ", " << cy << ", " << sx << ", " << sy << ", " << rho << std::endl;
    //integrate beam overlap with target cylinder
    double x0  =  cx;
    double y0  =  cy;
    double dbin = 0.1;
    double dx = dbin;
    double dy = dbin;
    double r    = 4.75;
    double rr   = r*r;
    double rho2 = rho*rho;
    double xmin = -r;
    double ymin = -r;
    int imax = round((2.0*r)/dx);
    int jmax = round((2.0*r)/dy);
    double sum =  0.0;
    double x = xmin;
    //cout <<sx<<"\t"<<sy<<"\t"<<imax<<"\t"<<jmax<<endl;
    for (int i=0;i<=imax;i++) {
      double y = ymin;
      for (int j=0;j<=jmax;j++) {
	if ( (x*x+y*y)<rr ) {
	  double tx = (x+x0)/sx;
	  double ty = (y+y0)/sy;
	  double z = tx*tx - 2.0*rho*tx*ty + ty*ty;
	  double t = exp(-z/(2.0*(1.0-rho2)));
	  sum = sum + t;
	  //std::cout << " t: " << t << " z: " << z << std::endl;
	}
	//      cout << i<<"\t"<<j<<"\t"<<x<<"\t"<<y<<"\t"<<sum<<endl;
	y = y + dy;
      }
      x = x + dx;
    }
    sum = sum*dx*dy/(2.0*3.14159*sx*sy*sqrt(1.0-rho2));

    //std::cout << sum << ", " << dx << ", " << dy  << ", " << rho2<< std::endl;

    // add a guard for double precision
    if(sum >= 1.) return -10000.;
    return log10(1-sum);
  }
  
}
