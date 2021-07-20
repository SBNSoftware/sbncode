//Built from PrimaryHadronSWCentralSplineVariation.cxx in ubcode/EventWeight/Calculator/Flux/BNBPrimaryHadron directory

#include "FluxCalcPrep.h"
#include "TSpline.h"


namespace sbn {
	namespace evwgh {

		std::pair<bool, double> FluxWeightCalc::PHSWCSVWeightCalc(simb::MCFlux flux, std::vector<float> rand){

			// 
			//  Largely built off the MiniBooNE code 
			//  but this is intended to expand beyond it
			//
			//  Edits from MiniBooNE code:
			//
			//  JZ (6/2017) : Remove max weight, set to max_limit of a double
			//  JZ (6/2017) : Changed guards on c9 to be double point precision 
			//   

			bool parameters_pass = true;

			double c1 = SWParam[0];
			double c2 = SWParam[1];
			double c3 = SWParam[2];
			double c4 = SWParam[3];
			double c5 = SWParam[4];
			double c6 = SWParam[5];
			double c7 = SWParam[6];
			double c8 = SWParam[7];
			double c9 = 1.0; // This isn't in the table but it is described in the text

			//  Lay out the event kinimatics 
			double HadronMass = 0.13957010;

			//      if(fabs(fprimaryHad) == 211) HadronMass = 0.13957010; //Charged Pion
			//      else{ 
			//	throw art::Exception(art::errors::StdException)
			//	  << "sanford-wang is only configured for charged pions ";
			//      }

			TLorentzVector HadronVec; 
			double HadronPx = flux.ftpx;
			double HadronPy = flux.ftpy;
			double HadronPz = flux.ftpz;
			double HadronE  = sqrt(HadronPx*HadronPx + 
					HadronPy*HadronPy + 
					HadronPz*HadronPz + 
					HadronMass*HadronMass);
			HadronVec.SetPxPyPzE(HadronPx,HadronPy,HadronPz,HadronE);

			////////
			//  
			//   Based on MiniBooNE code to evaluate the theta value within below the maximum 
			//    HARP theta coverage. This helps keep the splines well formed and constrains  
			//    the uncertainties at very low neutrino energy
			//
			////////

			double ThetaOfInterest;
			if(HadronVec.Theta() > 0.195){
				ThetaOfInterest = 0.195;
			}
			else{
				ThetaOfInterest = HadronVec.Theta();
			}


			// Get Initial Proton Kinitmatics 
			//   CURRENTLY GSimple flux files drop information about 
			//   the initial state proton, but that this 
			TLorentzVector ProtonVec;
			double ProtonMass = 0.9382720;
			double ProtonPx = 0;
			double ProtonPy = 0;
			double ProtonPz = 8.89; //GeV
			double ProtonE  = sqrt(ProtonPx*ProtonPx +
					ProtonPy*ProtonPy +
					ProtonPz*ProtonPz +
					ProtonMass*ProtonMass);
			ProtonVec.SetPxPyPzE(ProtonPx,ProtonPy,ProtonPz,ProtonE);

			//Sanford-Wang Parameterization 
			//  Eq 11 from PhysRevD.79.072002

			double CV = c1 * pow(HadronVec.P(), c2) * 
				(1. - HadronVec.P()/(ProtonVec.P() - c9)) *
				exp(-1. * c3 * pow(HadronVec.P(), c4) / pow(ProtonVec.P(), c5)) *
				exp(-1. * c6 * ThetaOfInterest *(HadronVec.P() - c7 * ProtonVec.P() * pow(cos(ThetaOfInterest), c8)));

			// Check taken from MiniBooNE code
			if((HadronVec.P()) > ((ProtonVec.P()) - (c9))){
				CV = 0;} 

			//////
			//
			// Now that we have our central value based on the Sanford-Wang parameterization we 
			// need to create variations around this value. MiniBooNE did this by performing 
			// spline fits to the HARP data. This is a reproduction of that code:
			// 
			////

			double* HARPmomentumBins = HARPmomentumBounds.data(); // Convert std::vector to array
			double* HARPthetaBins = HARPthetaBounds.data();       // Convert std::vector to array

			int Ntbins = int(HARPthetaBounds.size()) - 1;
			int Npbins = int(HARPmomentumBounds.size()) - 1;

			//
			//  Using the HARP cross section and covariance matrices 
			//  we will now create variations around the measured HARP
			//  meson production cross sections. 
			//

			// Important notes of the HARP data: 
			//  The cross section matrix has 13 momentum bins and 6 theta bins 
			//  The covariance matrix contains the correlated uncertainties across
			//  all 78 cross section measurements. 
			//
			//  The covariance matrix encodes the uncertainty on a given HARP 
			//  ANALYSIS BIN instead of being a 3D matrix. 
			//
			//  A HARP analysis bin is defined as 
			//  bin = momentum[theta[]] meaning that:
			//      analysis bin 0 is the zeroth theta and zeroth momentum bin
			//  but analysis bin 27 is the 3rd theta and 5th momentum bin
			// 
			// The first thing to do is convert our cross section matrix into 
			// an std::vector for the analysis bins
			//
			std::vector< double > HARPCrossSectionAnalysisBins;
			HARPCrossSectionAnalysisBins.resize(int(Ntbins*Npbins));

			int anaBin = 0;
			for(int pbin = 0; pbin < Npbins; pbin++){
				for(int tbin = 0; tbin < Ntbins; tbin++){	
					HARPCrossSectionAnalysisBins[anaBin] = HARPXSec[0][pbin][tbin];
					anaBin++;
				}
			}

			//
			// Now using this we can vary the cross section based on the HARP covariance matrix 
			// this will allow us to reweigh each cross section measurement based on
			// the multigaussian smearing of this matrix. 
			//

			std::vector< double > smearedHARPCrossSectionAnalysisBins = 
				MultiGaussianSmearing(HARPCrossSectionAnalysisBins, FitCov, fIsDecomposed,rand); 
			HARPCrossSectionAnalysisBins.clear();

			//
			//   Check all the smeared cross sections, if any come out to be negative then 
			//   we will not pass this given parameter set.
			//

			for(int check = 0; check < int(smearedHARPCrossSectionAnalysisBins.size()); check++){
				if(smearedHARPCrossSectionAnalysisBins[check] < 0){ parameters_pass = false;}
			}

			//  
			// With a checked set of smeared cross sections we can convert our analysis bin
			// std::vector back into a TMatrixD, though this is a bit annoying since TMatrices 
			// are an annoying format. 
			//
			TMatrixD*  smearedHARPXSec = new TMatrixD(Npbins,Ntbins);

			anaBin = 0;
			for(int pbin = 0; pbin < Npbins; pbin++){
				for(int tbin = 0; tbin < Ntbins; tbin++){
					smearedHARPXSec[0][pbin][tbin] = smearedHARPCrossSectionAnalysisBins[anaBin];	  
					anaBin++;
				}
			}

			//
			// We want to now make a vector of histograms that are 
			//  projections across the cross section matrix, 
			//  this means we want to have vectors where one histogram 
			//  is made for each bin on variable X where the histogram
			//  then contains bins for variable Y. 
			//
			//  This will allow us to do a 2D interpolation across the
			//  full parameter space using cubic spline fits 
			//

			std::vector< TH1F > MomentumBins;
			MomentumBins.resize(Ntbins);
			for(int bin = 0; bin < Ntbins; bin++){
				MomentumBins[bin] = TH1F("HARPp",";;;", Npbins, HARPmomentumBins);
				MomentumBins[bin].SetBit(kCanDelete);
			} 

			//Setup vectors of splines for fitting
			std::vector< TSpline3 > SplinesVsMomentum; 
			SplinesVsMomentum.resize(Ntbins);

			for(int tbin = 0; tbin < Ntbins; tbin++){
				for(int pbin = 0; pbin < Npbins; pbin++){

					//
					// We want to create spline at fixed theta values, that way we can 
					// study the cross section at the given meson momentum
					MomentumBins[tbin].SetBinContent(pbin+1, smearedHARPXSec[0][pbin][tbin]);		    
					// This will give us the cross section in discreet slices of theta (tbin)
					// but each histogram in the vector will be the cross section in bins of
					// momentum that can then be cublic splined  
					//
				}

				//
				// in a given theta slice we can now spline over the cross section entries
				// in bins of momentum
				//
				// Important note about Splines, MiniBooNE used constraints on the
				// second derivative of the first and final knot point and required 
				// that they both equal to zero, this is not naturally the case in 
				// in TSpline3 but it is default in DCSPLC (the Fortran CERN library spline function)
				// this helps to control the smoothness of the spline and minimizes the variation bin to bin
				// For TSpine3 this is controlled by:
				//
				//  b1 = constrain first knot 1st derivative 
				//  b2 = constrain first knot 2nd derivative 
				//  e1 = constrain final knot 1st derivative 
				//  e2 = constrain final knot 2nd derivative 
				//
				//  the numbers that follow are the values you constrain those conditions to
				//
				SplinesVsMomentum[tbin] = TSpline3(&(MomentumBins[tbin]),"b2e2",0,0); 

			}
			MomentumBins.clear();
			delete smearedHARPXSec;
			//
			// Now that we have the 1D Splines over momentum we want to know
			// what the cross section is for the meson's momentum in bins of
			// theta that way we can extract the cross section at the meson's
			// theta. 
			//
			TH1F* ThetaBins = new TH1F("HARPt",";;;", Ntbins, HARPthetaBins);
			ThetaBins->SetBit(kCanDelete);
			for(int tbin = 0; tbin < Ntbins; tbin++){
				ThetaBins->SetBinContent(tbin+1, SplinesVsMomentum[tbin].Eval(HadronVec.P()));
			} 

			////      
			// Now we can spline across the theta bins and we can extract the exact value of the 
			// cross section for the meson's theta value
			//
			//  Options described above.
			/////

			TSpline3* FinalSpline = new TSpline3(ThetaBins,"b2e2",0,0);

			double RW = FinalSpline->Eval(ThetaOfInterest);

			SplinesVsMomentum.clear();
			delete ThetaBins;
			delete FinalSpline;

			//
			// These guards are inherited from MiniBooNE code
			//
			// These were defined here: 
			// 
			//   cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/~checkout~/ ... 
			//              miniboone/AnalysisFramework/MultisimMatrix/src/MultisimMatrix.inc    
			//
			//  This forces any negative spline fit to be 1     

			////////
			//  Possible Bug.
			////////     
			// This seems to be a feature in the MiniBooNE code
			//  It looks like the intension is to set this to zero
			//  but it is set to 1 before it is set to zero  

			double weight = 1; 

			if(RW < 0 || CV < 0){
				weight = 1;
			}
			else if(fabs(CV) < 1.e-12){
				weight = 1;
			}
			else{
				weight *= RW/CV;
			}

			if(weight < 0) weight = 0; 
			if(weight > 30) weight = 30; 
			if(!(std::isfinite(weight))){
				std::cout << "SW+Splines : Failed to get a finite weight" << std::endl; 	
				weight = 30;
			}

			std::pair<bool, double> output(parameters_pass, weight);
			return output; 

		}// Done with the WeigthCalc function

	}  // namespace evwgh
}  // namespace sbn



