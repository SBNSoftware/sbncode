////////////////////////////////////////////////////////////////////////
// Class:       opHitFinderSBND
// Module Type: producer
// File:        opHitFinderSBND_module.cc
//
// This module produces an OpHit object for light analysis
// Created by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

#include "opHitFinderSBND.hh"

namespace opdet{
  opHitFinderSBND::opHitFinderSBND(fhicl::ParameterSet const & p, const detinfo::DetectorClocks *timeService)
// Initialize member data here.
  {
    fInputModuleName = p.get< std::string >("InputModule" );
    fBaselineSample  = p.get< int    >("BaselineSample"); //in ticks
    fSaturation      = p.get< double >("Saturation"   ); //in number of p.e.
    fArea1pePMT      = p.get< double >("Area1pePMT"   ); //in ADC*ns for PMTs
    fArea1peSiPM     = p.get< double >("Area1peSiPM"  ); //in ADC*ns for SiPMs
    fThresholdPMT    = p.get< double >("ThresholdPMT" ); //in ADC
    fThresholdArapuca = p.get< double >("ThresholdArapuca"); //in ADC
    fUseDenoising     = p.get< int   >("UseDenoising"); 
    fPulsePolarityPMT = p.get< int   >("PulsePolarityPMT");
    fPulsePolarityArapuca = p.get< int   >("PulsePolarityArapuca");

    fSampling = (timeService->OpticalClock().Frequency())*500./64; //in GHz. This number is wrong! Therefore the hard coded value
  }

  std::vector<recob::OpHit> opHitFinderSBND::MakeHits(const std::vector<raw::OpDetWaveform> &waveforms) {
    std::vector< recob::OpHit > ret;

    size_t timebin=0;
    double FWHM=1, Area=1, phelec, fasttotal=3./4., rms=0, amplitude=0, time=0;
    unsigned short frame=1;
    int histogram_number = 0;

    for(auto const& wvf : waveforms) {
	fChNumber = wvf.ChannelNumber();
	histname.str(std::string());
        histname << "opchannel_" << fChNumber << "_histo_" << histogram_number;
	wvfHist = new TH1D(histname.str().c_str(), "Histogram", wvf.size(),0, double(wvf.size()));
	for(unsigned int i=0;i<wvf.size();i++){
	  wvfHist->SetBinContent(i,wvf[i]);
	}
	subtractBaseline(wvfHist, map.pdName(fChNumber), rms);
	if((map.pdName(fChNumber)=="pmt") || (map.pdName(fChNumber)== "barepmt")){
	}else{
	  if(fUseDenoising==1) denoise(wvfHist);    
	}
        int i=1;
        while(findPeak(wvfHist,timebin,Area,rms,amplitude,map.pdName(fChNumber))){
          time = wvf.TimeStamp() + (double)timebin/fSampling;
	  if(map.pdName(fChNumber)=="pmt" || map.pdName(fChNumber) == "barepmt"){
	    phelec=Area/fArea1pePMT;
        //    std::cout << 0 << " " << time << " " << Area << " " << phelec << std::endl;
	  }else{
	    phelec=Area/fArea1peSiPM;
          //  std::cout << 1 << " " << time << " " << Area << " " << phelec << std::endl;
	  }
          i++;
	  recob::OpHit opHit(fChNumber, time, time, frame, FWHM, Area, amplitude, phelec, fasttotal);//including hit info: OpChannel, PeakTime, PeakTimeAbs, Frame, Width, Area, PeakHeight, PE, FastToTotal
          ret.push_back(opHit);
        }
        histogram_number += 1;
	delete wvfHist;
    }
    return ret;
  }

  void opHitFinderSBND::subtractBaseline(TH1D* h, std::string pdtype, double& rms){
    double baseline = 0.0; 
    rms=0.0;
    int cnt = 0;
    for(int i=0; i<fBaselineSample; i++){
	baseline+=h->GetBinContent(i);
	rms+=pow((h->GetBinContent(i)),2.0);
        cnt++;
    }
    baseline=baseline/cnt;
    rms=sqrt(rms/cnt-baseline*baseline);
    rms=rms/sqrt(cnt-1);

    if(pdtype=="pmt" || pdtype == "barepmt"){
        for(int i=0; i<h->GetNbinsX(); i++)h->SetBinContent(i, (fPulsePolarityPMT*(h->GetBinContent(i)-baseline)));
    }else{
        for(int i=0; i<h->GetNbinsX(); i++)h->SetBinContent(i, (fPulsePolarityArapuca*(h->GetBinContent(i)-baseline)));
    }
  }

  bool opHitFinderSBND::findPeak(TH1D* h, size_t& time, double& Area, double rms, double& amplitude, std::string type){

    //Gets info from highest peak and suppress it in histogram
    double aux = h->GetMaximum();
    double max;
    size_t time_end, bin, binmax = h->GetMaximumBin();
    int threshold;

    if(type=="pmt" || type == "barepmt"){
      threshold=fThresholdPMT;
    }else{
      threshold=fThresholdArapuca;
    }

    bin=binmax;
    amplitude=aux;
    max=aux;

    if(aux<threshold)return false;

    while(aux>=rms){
        bin++;	
        aux = h->GetBinContent(bin);
    }
    time_end=bin-1; //looking for the length of the peak

    aux=max;
    while(aux>=rms){
        bin--;
        aux = h->GetBinContent(bin);
    }
    time = bin+1; //for rise time

    Area=(h->Integral(time,time_end))/fSampling;

    bin = time;
    aux = h->GetBinContent(time);  

   while(aux>=rms){
        h->SetBinContent(bin,0.0);
        bin++;
        aux = h->GetBinContent(bin);
    }
    //std::cout << time << " " << time_end << " " << (time_end - time);
    time=binmax; //returning the peak time
    return true;		
  }

  void opHitFinderSBND::denoise(TH1D* h){

    int wavelength = (int)h->GetNbinsX();
    float lambda = 10.0;
    float* input;
    float* output;

    //std::cout <<"denoise wavl:" << wavelength << std::endl;

    input = (float*)malloc(sizeof(*input)*wavelength);
    output = (float*)malloc(sizeof(*output)*wavelength);

    for(int i=0; i<wavelength; i++){
      input[i] = h->GetBinContent(i);
      output[i] = input[i];
    }

    //std::cout <<"denoise full array" << std::endl;

    TV1D_denoise(input, output, (const int)wavelength, (const float)lambda);

    //std::cout <<"denoise denoise" << std::endl;

    for(int i=0; i<wavelength; i++){
      if(output[i])h->SetBinContent(i,output[i]);
    }

 }

  void opHitFinderSBND::TV1D_denoise(float* input, float*& output, const int width, const float lambda) {
	if (width>0) {				/*to avoid invalid memory access to input[0]*/
		int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
		float umin=lambda, umax=-lambda;	/*u is the dual variable*/
		float vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
		int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const float twolambda=2.0*lambda;	/*auxiliary variable*/
		const float minlambda=-lambda;		/*auxiliary variable*/
		for (;;) {				/*simple loop, the exit test is inside*/
			while (k==width-1) {	/*we use the right boundary condition*/
				if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
					do output[k0++]=vmin; while (k0<=kminus);
					umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
				} else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
					do output[k0++]=vmax; while (k0<=kplus);
					umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
				} else {
					vmin+=umin/(k-k0+1); 
					do output[k0++]=vmin; while(k0<=k); 
					return;
				}
			}
			if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
				do output[k0++]=vmin; while (k0<=kminus);
				vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
				umin=lambda; umax=minlambda;
			} else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
				do output[k0++]=vmax; while (k0<=kplus);
				vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
				umin=lambda; umax=minlambda;
			} else { 	/*no jump necessary, we continue*/
				k++;
				if (umin>=lambda) {		/*update of vmin*/
					vmin+=(umin-lambda)/((kminus=k)-k0+1);
					umin=lambda;
				} 
				if (umax<=minlambda) {	/*update of vmax*/
					vmax+=(umax+lambda)/((kplus=k)-k0+1);
					umax=minlambda;
				} 	
			}
		}
	}
  }

}
