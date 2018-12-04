#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>
#include "Covariance.h"

#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>

#include <TFile.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TROOT.h>
// #include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCanvas.h>
#include <core/Event.hh>

namespace ana {
namespace SBNOsc {

Covariance::EventSample::EventSample(const Json::Value &config, unsigned nUniverses) {
    // get configuration stuff
    Json::Value bins = config.get("binlims", "");
    for (auto const& bin: bins) {
      fBins.push_back(bin.asDouble());
    }
    fName = config.get("name", "").asString(); 
    fScaleFactor = config.get("scalefactor", 0.).asDouble();

    // setup histograms
    std::string cv_title = fName + " Central Value";
    fCentralValue = new TH1D(cv_title.c_str(), cv_title.c_str(), fBins.size() - 1, &fBins[0]);

    for (unsigned i = 0; i < nUniverses; i++) {
      std::string uni_title = fName + " Universe " + std::to_string(i);
      fUniverses.push_back(
        new TH1D(uni_title.c_str(), uni_title.c_str(), fBins.size() - 1, &fBins[0])
      );
    }

}


// Gets scale factors (weights) for different universes
std::vector <double> GetUniWeights(const std::map <std::string, std::vector <double> > &weights, const std::vector<std::string> &keys, int n_unis) {
    
    // Tentative format: universe u scale factor is the product of the u-th entries on each vector 
    // inside the map. For vectors with less than u entries, use the (u - vec_size)-th entry
    
    std::vector <double> uweights;
    
    for (int u = 0; u < n_unis; u++) {
        double weight = 1.;
        for (auto const &key: keys) {
            const std::vector<double>& this_weights = weights.at(key);
            int wind = u % this_weights.size();
            weight *= this_weights.at(wind);
        }
        uweights.push_back(weight);
    }
    return uweights;
}


void Covariance::Initialize(Json::Value *config) {
    // must have config
    assert(config != NULL);
    // Weight and universe stuff
    // WeightKey can either be the name of a single weight, or a name signifying a list of weights,
    // or a list of weights to be oscillated
    
    Json::Value configWeightKey = (*config)["Covariance"].get("WeightKey", ""); 
    for (auto const& keyName: configWeightKey) {
        fWeightKeys.push_back(keyName.asString());
    }

    // uniformly applied weights
    if ((*config)["Sensitivity"].isMember("UniformWeights") &&
        (*config)["Sensitivity"]["UniformWeights"].isArray()) {
        for (auto const &key: (*config)["Sensitivity"]["UniformWeights"]) {
            fUniformWeights.push_back(key.asString());
        }
    }

    // number of universes to be used
    fNumAltUnis = (*config)["Covariance"].get("NumAltUnis", 0).asInt();
    
    // Type of energy
    fEnergyType = (*config)["Covariance"].get("EnergyType", "").asString();
    
    // Further selection and rejection 'efficiencies'
    fSelectionEfficiency = (*config)["Covariance"].get("SelectionEfficiency", 1.0).asDouble();
    fBackgroundRejection = (*config)["Covariance"].get("BackgroundRejection", 0.0).asDouble();
    
    // get event samples
    Json::Value configEventSamples = (*config)["EventSamples"];
    for (auto const& sample: configEventSamples) {
        fEventSamples.emplace_back(sample, fNumAltUnis);
    }
    
    // set output directory
    fOutputFile = (*config)["Covariance"].get("OutputFile", "").asString();
    // whether to save things
    fSaveCentralValue = (*config)["Covariance"].get("SaveCentralValue", false).asBool();
    fSaveUniverses = (*config)["Covariance"].get("SaveUniverses", false).asBool();
    // start out at the zeroth sample
    fSampleIndex = 0;
}

void Covariance::ProcessEvent(const Event *event) {
    // iterate over each interaction in the event
    for (int n = 0; n < event->reco.size(); n++) {
        unsigned truth_ind = event->reco[n].truth_index;

        // Get energy
        double true_nuE = event->reco[n].truth.neutrino.energy;
        double nuE; 
        if (fEnergyType == "CCQE") {
            nuE = event->truth[truth_ind].neutrino.eccqe;
        } else if (fEnergyType == "True") {
            nuE = true_nuE;
        } else if (fEnergyType == "Reco") {
            nuE = event->reco[n].reco_energy;
        }
    
        // Apply selection (or rejection) efficiencies
        int isCC = event->truth[truth_ind].neutrino.iscc;
        double wgt = isCC*(fSelectionEfficiency) + (1-isCC)*(1 - fBackgroundRejection);
        // apply scaling from fScaleFactor
        wgt *= fEventSamples[fSampleIndex].fScaleFactor;
        // apply uniform weights
        for (auto const &key: fUniformWeights) {
            wgt *= event->truth[truth_ind].weights.at(key)[0];
        }
    
        // Get weights for each alternative universe
        std::vector <double> uweights = GetUniWeights(event->truth[truth_ind].weights, fWeightKeys, fNumAltUnis);
    
        // Add to base and bkg count histogram
        fEventSamples[fSampleIndex].fCentralValue->Fill(nuE, wgt);
    
        // Fill alternative universe histograms
        for (int u = 0; u < uweights.size(); u++) {
            fEventSamples[fSampleIndex].fUniverses[u]->Fill(nuE, wgt*uweights[u]);
        }
    }
}


void Covariance::FileCleanup(TTree *eventTree) {
    // onto the next sample
    fSampleIndex ++;
}

void Covariance::GetCovs() {
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Getting covs..." << std::endl;

    // get the total number of bins across all samples
    unsigned num_bins = 0;
    for (auto const &sample: fEventSamples) {
       num_bins += sample.fCentralValue->GetNbinsX();
    }
    
    // Covariance and fractional covariance
    cov = new TH2D("Covariance", "Covariance", num_bins, 0, num_bins, num_bins, 0, num_bins);
    fcov = new TH2D("Fractional Covariance", "Fractional Covariance", num_bins, 0, num_bins, num_bins, 0, num_bins);

    // indexes into covariance
    unsigned cov_i = 0;
    unsigned cov_j = 0;
    // iterate over samples and bins
    for (auto const &sample_i: fEventSamples) {
        for (unsigned sample_i_bin_index = 0; sample_i_bin_index < sample_i.fCentralValue->GetNbinsX(); sample_i_bin_index++) {
            // Get Central Value
            double i_central = sample_i.fCentralValue->GetBinContent(sample_i_bin_index+1);
            for (auto const &sample_j: fEventSamples) {
                for (unsigned sample_j_bin_index = 0; sample_j_bin_index < sample_j.fCentralValue->GetNbinsX(); sample_j_bin_index++) {
                    // Get Central Value
                    double j_central = sample_j.fCentralValue->GetBinContent(sample_j_bin_index+1);
                    // calculate covariance
                    double cov_value = 0.; 
                    for (int u = 0; u < fNumAltUnis; u++) {
                        // get variations
                        double i_variation = sample_i.fUniverses[u]->GetBinContent(sample_i_bin_index+1);
                        double j_variation = sample_j.fUniverses[u]->GetBinContent(sample_j_bin_index+1);
                        cov_value += (i_central - i_variation) * (j_central - j_variation);
                    }
                    // average
                    if (fNumAltUnis != 0) cov_value /= fNumAltUnis;

                    // calculate fractional covariance
                    double fcov_value;
                    if (i_central * j_central < 1e-6) {
                        fcov_value = 0.;
                    }
                    else {
		        fcov_value = cov_value / (i_central * j_central);
                    }

                    // save values
                    cov->SetBinContent(cov_i + 1, cov_j + 1, cov_value);
                    fcov->SetBinContent(cov_i + 1, cov_j + 1, fcov_value);
                    // update covariance position
                    cov_j ++;
                }
            }   
	    // update covariance position
            cov_i ++;
            // reset j
            cov_j = 0;
        }
    } 

    // Pearson Correlation Coefficients
    corr = new TH2D("Correlation", "Correlation", num_bins, 0, num_bins, num_bins, 0, num_bins);
    
    for (int i = 0; i < num_bins; i++) {
        double covii = cov->GetBinContent(i+1, i+1);
        for (int j = 0; j < num_bins; j++) {
            double covjj = cov->GetBinContent(j+1, j+1);
            double covij = cov->GetBinContent(i+1, j+1);

            double corrij;
            // handle case where covariance is 0
            if (covii*covjj < 1e-6) {
                corrij = 0.;
            }
            else {
                corrij = covij / TMath::Sqrt(covii*covjj);
            }
            corr->SetBinContent(i+1, j+1, corrij);
        }
    }
    
    // Add bin labels
    std::vector <TH2D*> hists = {cov, fcov, corr};
    unsigned bin = 0;
    for (auto const &sample: fEventSamples) {
        // Set label
        for (TH2D* hist : hists) {
            hist->GetXaxis()->SetBinLabel(1+bin, sample.fName.c_str());
            hist->GetYaxis()->SetBinLabel(1+bin, sample.fName.c_str());
        }
        // update bin positions
        bin += sample.fBins.size() - 1;

    }

    for (TH2D* hist : hists) { hist->GetXaxis()->LabelsOption("h"); hist->GetYaxis()->LabelsOption("v"); }
    
    std::cout << std::endl << "  Got covs." << std::endl;
    
}
    
void Covariance::Write() {
    // if output file not set, don't write
    if (fOutputFile == "") return;

    // Write Covariances
    TFile *output = TFile::Open(fOutputFile.c_str(), "recreate");
    assert(output && output->IsOpen());
    
    cov->Write();
    fcov->Write();
    corr->Write();
    // write histos
    if (fSaveCentralValue) {
      for (auto const &sample: fEventSamples) {
        sample.fCentralValue->Write();
      }
    }
    if (fSaveUniverses) {
        for (auto const &sample: fEventSamples) {
            for (TH1D *h: sample.fUniverses) {
                h->Write();
            }
        }
    }
    // close file
    output->Close(); 
}

// Turn the covariance TH2D into a matrix -- include stat uncertainty
TMatrixDSym Covariance::CovarianceMatrix() {
    TMatrixDSym E_mat(cov->GetNbinsX());
    for (int i = 0; i < cov->GetNbinsX(); i++) {
        for (int j = 0; j < cov->GetNbinsY(); j++) {
            E_mat[i][j] = cov->GetBinContent(i+1, j+1);
            // add statistical uncertainty
            if (i == j) {
                unsigned stat_ind = i;
                double stat_uncertainty;
                // get index into CV bin to get stat uncertainty
                for (auto const &sample: fEventSamples) {
                    if (stat_ind < sample.fCentralValue->GetNbinsX()) {
                        stat_uncertainty = sample.fCentralValue->GetBinContent(1+stat_ind);
                        break;
                    }
                    else stat_ind -= sample.fCentralValue->GetNbinsX();
                }
                E_mat[i][i] += stat_uncertainty;
            }
        }
    }
    return E_mat;
}
    
}   // namespace SBNOsc
}   // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Covariance);

