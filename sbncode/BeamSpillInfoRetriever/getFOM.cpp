//#############################################################################################
// Script: Figure of Merit for BNB Spills using SBND information adapted from MicroBooNE FOM
//
// Author: Max Dubnowski
// Contact: maxdub@sas.upenn.edu or sbn slack
//
//###############################################################################################

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
  {
    double fom=0;

    double hp875_offset= spill.HP875Offset;
    double vp875_offset= spill.VP875Offset;
    double hptg1_offset= spill.HPTG1Offset;
    double vptg1_offset= spill.VPTG1Offset;
    double hptg2_offset= spill.HPTG2Offset;
    double vptg2_offset= spill.VPTG2Offset;
        
    //Decides which position monitor to check first
    int useHTG = 0;
    int useVTG = 0;
    
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
    
    
    tor860.push_back(spill.TOR860);
    tor875.push_back(spill.TOR875);
    hp875.push_back(spill.HP875);
    vp875.push_back(spill.VP875);
    hptg1.push_back(spill.HPTG1);
    vptg1.push_back(spill.VPTG1);
    hptg2.push_back(spill.HPTG2);
    vptg2.push_back(spill.VPTG2);
    
    double tor;
    if (tor860.size()>0)
      tor=tor860[0];
    else if (tor875.size()>0)
      tor=tor875[0];
    else
      return -1;
        
    
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
    
    fom=1-pow(10,sbn::calcFOM(horpos,horang,verpos,verang,tor,tgtsx,tgtsy));
    return fom;
  }
  


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
    double fom_a=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    //swim to center of target
    sbn::swimBNB(centroid1,sigma1,
		 identity, begtocnt,
		 cx, cy, sx, sy, rho);
    double fom_b=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    //swim to downstream of target
    sbn::swimBNB(centroid1,sigma1,
		 cnttodns, begtodns,
		 cx, cy, sx, sy, rho);
    double fom_c=sbn::func_intbivar(cx, cy, sx*scalex, sy*scaley, rho);
    // add a guard for double precision

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
	  }
	}
      }
    }
    //get beam sigma
    sx  = sqrt(sigma2[0][0])*1000.0;
    sy  = sqrt(sigma2[2][2])*1000.0;
    rho = sigma2[0][2]/sqrt(sigma2[0][0]*sigma2[2][2]);
  }
  
  
  
  double func_intbivar(const double cx, const double cy, const double sx, const double sy, const double rho )
  {
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
    for (int i=0;i<=imax;i++) {
      double y = ymin;
      for (int j=0;j<=jmax;j++) {
	if ( (x*x+y*y)<rr ) {
	  //Looking up Bivariate Normal Distribution, tx and ty should be x-x0 not x+x0; doesn't change much for a symmetric spill
	  double tx = (x-x0)/sx;
	  double ty = (y-y0)/sy;
	  //Code used in the MicroBooNE code, possible error in the formula
	  //double tx = (x+x0)/sx;
	  //double ty = (y+y0)/sy;
	  double z = tx*tx - 2.0*rho*tx*ty + ty*ty;
	  double t = exp(-z/(2.0*(1.0-rho2)));
	  sum = sum + t;
	}
	y = y + dy;
      }
      x = x + dx;
    }
    sum = sum*dx*dy/(2.0*3.14159*sx*sy*sqrt(1.0-rho2));


    // add a guard for double precision
    if(sum >= 1.) return -10000.;
    return log10(1-sum);
  }
  
}
