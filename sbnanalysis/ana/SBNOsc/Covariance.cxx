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
#include <TCanvas.h>
#include <TMath.h>
#include <TCanvas.h>
#include <core/Event.hh>

namespace ana {
namespace SBNOsc {

Covariance::EventSample::EventSample(const fhicl::ParameterSet &config, unsigned nUniverses, unsigned nVariations) {
    // get configuration stuff
    fBins = config.get<std::vector<double> >("binlims");
    fName = config.get<std::string>("name", "");
    fScalePOT = config.get<double>("ScalePOT", -1);
    fScaleCC = config.get<double>("ScaleCC", 1);
    fScaleNC = config.get<double>("ScaleNC", 1);
    fScaleEnergy = config.get<double>("ScaleEnergy", 1);
    fPOT = 0.;
    fEnergyBinScale = config.get<std::vector<double>>("energy_bin_scale", {});
    if (fEnergyBinScale.size() == 0) {
      fEnergyBinScale = std::vector<double>(fBins.size() - 1, 1.0);
    }
    else assert(fEnergyBinScale.size() == fBins.size() - 1);

    // setup histograms
    std::string cv_title = fName + " Central Value";
    fCentralValue = new TH1D(cv_title.c_str(), cv_title.c_str(), fBins.size() - 1, &fBins[0]);

    for (unsigned var_i = 0; var_i < nVariations; var_i++) {
      fUniverses.emplace_back();
      for (unsigned i = 0; i < nUniverses; i++) {
        std::string uni_title = fName + " Variation " + std::to_string(var_i) + " Universe " + std::to_string(i);
        fUniverses[var_i].push_back(
          new TH1D(uni_title.c_str(), uni_title.c_str(), fBins.size() - 1, &fBins[0])
        );
      }
    }

}


// Gets scale factors (weights) for different universes
std::vector <float> GetUniWeights(const std::map <std::string, std::vector <float> > &weights, const std::vector<std::string> &keys, int n_unis, int uni_offset) {
    std::vector <float> uweights;
    for (int u = 0; u < n_unis; u++) {
        float weight = 1.;
        for (auto const &key: keys) {
            const std::vector<float>& this_weights = weights.at(key);
            int wind = u + uni_offset;
            weight *= this_weights.at(wind);
        }
        uweights.push_back(weight);
    }
    return uweights;
}


void Covariance::Initialize(fhicl::ParameterSet* config) {
    // must have config
    assert(config != NULL);
    // Weight and universe stuff
    // WeightKey can either be the name of a single weight, or a name signifying a list of weights,
    // or a list of weights to be oscillated
    
    fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("Covariance");

    fWeightKeys = pconfig.get<std::vector<std::vector<std::string>>>("WeightKey");
    fWeightKeysCC = pconfig.get<std::vector<std::vector<std::string>>>("WeightKeyCC", {});
    fWeightKeysNC = pconfig.get<std::vector<std::vector<std::string>>>("WeightKeyNC", {});
    fNVariations = fWeightKeys.size();

    // all should have the same size
    assert(fWeightKeysCC.size() == 0 || fWeightKeysCC.size() == fWeightKeys.size());
    assert(fWeightKeysNC.size() == 0 || fWeightKeysNC.size() == fWeightKeys.size());

    // merge weight keys into CC/NC keys
    if (fWeightKeysCC.size() == 0) {
      fWeightKeysCC = fWeightKeys;
    }
    else {
      for (unsigned i = 0; i < fNVariations; i++) {
        fWeightKeysCC[i].insert(fWeightKeysCC[i].end(), fWeightKeys[i].begin(), fWeightKeys[i].end());
      }
    }
    if (fWeightKeysNC.size() == 0) {
      fWeightKeysNC = fWeightKeys;
    }
    else {
      for (unsigned i = 0; i < fNVariations; i++) {
        fWeightKeysNC[i].insert(fWeightKeysNC[i].end(), fWeightKeys[i].begin(), fWeightKeys[i].end());
      }
    }

    // number of universes to be used
    fNumAltUnis = pconfig.get<int>("NumAltUnis", 0);
    // and offset into the eventweight vector
    fAltUniOffset = pconfig.get<int>("AltUniOffset", 0);
    
    // Type of energy
    fEnergyType = pconfig.get<std::string>("EnergyType", "Reco");
    
    // Further selection and rejection 'efficiencies'
    fSelectionEfficiency = pconfig.get<double>("SelectionEfficiency", 1.0);
    fBackgroundRejection = pconfig.get<double>("BackgroundRejection", 0.0);

    // maximum weight for covariance calculation
    // If an event has any weight larger than this value, it will be 
    // thrown away for the purposes of covariance construction
    fWeightMax = pconfig.get<double>("WeightMax", -1.);
     
    // get event samples
    std::vector<fhicl::ParameterSet> configEventSamples = \
      config->get<std::vector<fhicl::ParameterSet> >("EventSamples");
    for (const fhicl::ParameterSet& sample : configEventSamples) {
      fEventSamples.emplace_back(sample, fNumAltUnis, fNVariations);
    }
    
    // set output directory
    fOutputFile = pconfig.get<std::string>("OutputFile", "");

    // whether to save things
    fSaveCentralValue = pconfig.get<bool>("SaveCentralValue", false);
    fSaveUniverses = pconfig.get<bool>("SaveUniverses", false);

    // start out at the zeroth sample
    fSampleIndex = 0;
}

void Covariance::ProcessSubRun(const SubRun *subrun) {
  fEventSamples[fSampleIndex].fPOT += subrun->totgoodpot;
}

void Covariance::ProcessEvent(const event::Event *event) {
    // iterate over each interaction in the event
    for (int n = 0; n < event->reco.size(); n++) {
        unsigned truth_ind = event->reco[n].truth_index;

        // Get energy
        int idx = event->reco[n].truth_index;
        if (idx < 0 || idx > event->truth.size()) {
          continue;
        }
        double true_nuE = event->truth[idx].neutrino.energy;
        double nuE; 
        if (fEnergyType == "CCQE") {
            nuE = event->truth[truth_ind].neutrino.eccqe;
        } else if (fEnergyType == "True") {
            nuE = true_nuE;
        } else if (fEnergyType == "Reco") {
            nuE = event->reco[n].reco_energy;
        }
        nuE *= fEventSamples[fSampleIndex].fScaleEnergy;
    
        // Apply selection (or rejection) efficiencies
        int isCC = event->truth[truth_ind].neutrino.iscc;
        // double wgt = isCC*(fSelectionEfficiency) + (1-isCC)*(1 - fBackgroundRejection);
        // apply scaling from fScaleFactor
        // wgt *= fEventSamples[fSampleIndex].fScaleFactor;
        // apply uniform weights
        // for (auto const &key: fUniformWeights) {
        //     wgt *= event->truth[truth_ind].weights.at(key)[0];
        // }
        double wgt = event->reco[n].weight;
        // apply POT scaling if configured
        if (fEventSamples[fSampleIndex].fScalePOT > 0) {
          wgt *= fEventSamples[fSampleIndex].fScalePOT / fEventSamples[fSampleIndex].fPOT;
        }

        if (event->truth[truth_ind].neutrino.iscc) {
          wgt *= fEventSamples[fSampleIndex].fScaleCC;
        }
        else {
          wgt *= fEventSamples[fSampleIndex].fScaleNC;
        }
    
        std::vector<std::vector <float>> uweights; 
        // Get weights for each alternative universe
        if (isCC) {
          for (std::vector<std::string> &weight_keys: fWeightKeysCC) {
            uweights.push_back(
              GetUniWeights(event->truth[truth_ind].weights, weight_keys, fNumAltUnis, fAltUniOffset)
            );
          }
        }
        else {
          for (std::vector<std::string> &weight_keys: fWeightKeysNC) {
            uweights.push_back(
              GetUniWeights(event->truth[truth_ind].weights, weight_keys, fNumAltUnis, fAltUniOffset)
            );
          }
        }

        // see if weight is too big
        if (fWeightMax > 0) {
          double max_weight = 0;
          for (std::vector<float> &weights: uweights) {
            float this_max_weight = *std::max_element(weights.begin(), weights.end());
            if (this_max_weight > max_weight) max_weight = this_max_weight;
          }
          if (max_weight > fWeightMax) {
            std::cout << "Weight (" << max_weight << ") over cap (" << fWeightMax << 
              ") for event with energy (" << nuE << ") in sample " << fEventSamples[fSampleIndex].fName << std::endl;
            continue;
          }
        }
    
        // Add to base and bkg count histogram
        fEventSamples[fSampleIndex].fCentralValue->Fill(nuE, wgt);
    
        // Fill alternative universe histograms
        for (unsigned variation = 0; variation < uweights.size(); variation++) {
          for (int u = 0; u < uweights[variation].size(); u++) {
              fEventSamples[fSampleIndex].fUniverses[variation][u]->Fill(nuE, wgt*uweights[variation][u]);
          }
        }
    }
}

void Covariance::Scale() {
  for (int sample_i = 0; sample_i < fEventSamples.size(); sample_i++) {
    for (int i = 1; i < fEventSamples[sample_i].fCentralValue->GetNbinsX()+1; i++) {
      //std::cout << "pre sample: " << sample_i << " bin content: " << fEventSamples[sample_i].fCentralValue->GetBinContent(i);
      fEventSamples[sample_i].fCentralValue->SetBinContent(i,
        fEventSamples[sample_i].fEnergyBinScale[i-1] *  fEventSamples[sample_i].fCentralValue->GetBinContent(i));
      //std::cout << "post sample: " << sample_i << " bin content: " << fEventSamples[sample_i].fCentralValue->GetBinContent(i);
    }
    for (int variation = 0; variation < fNVariations; variation++) {
      for (int uni = 0; uni < fNumAltUnis; uni++) {
        for (int i = 1; i < fEventSamples[sample_i].fUniverses[variation][uni]->GetNbinsX(); i++) {
          double content = fEventSamples[sample_i].fUniverses[variation][uni]->GetBinContent(i);
          fEventSamples[sample_i].fUniverses[variation][uni]->SetBinContent(i, content * fEventSamples[sample_i].fEnergyBinScale[i-1]);
        }
      }
    }
  }
}


void Covariance::FileCleanup(TTree *eventTree) {
    // onto the next sample
    fSampleIndex ++;
}

void Covariance::GetCovs() {
    for (unsigned i = 0; i < fNVariations; i++) GetCovPerVariation(i);
}

void Covariance::GetCovPerVariation(unsigned variation) {
    
    //// Get covariances, fractional covariances and correlation coefficients
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::cout << std::endl << "Getting covs..." << std::endl;

    // get the total number of bins across all samples
    unsigned num_bins = 0;
    for (auto const &sample: fEventSamples) {
       num_bins += sample.fCentralValue->GetNbinsX();
    }
    

    std::string covariance_name = "Covariance "  + std::to_string(variation);
    std::string fractional_covariance_name = "Fractional Covariance "  + std::to_string(variation);
    // Covariance and fractional covariance
    cov.emplace_back(covariance_name.c_str(), covariance_name.c_str(), num_bins, 0, num_bins, num_bins, 0, num_bins);
    fcov.emplace_back(fractional_covariance_name.c_str(), fractional_covariance_name.c_str(), num_bins, 0, num_bins, num_bins, 0, num_bins);

    // indexes into covariance
    unsigned cov_i = 0;
    unsigned cov_j = 0;
    // iterate over samples and bins
    for (auto const &sample_i: fEventSamples) {
        for (unsigned sample_i_bin_index = 0; sample_i_bin_index < sample_i.fCentralValue->GetNbinsX(); sample_i_bin_index++) {
            // Get Central Value
            float i_central = sample_i.fCentralValue->GetBinContent(sample_i_bin_index+1);
            for (auto const &sample_j: fEventSamples) {
                for (unsigned sample_j_bin_index = 0; sample_j_bin_index < sample_j.fCentralValue->GetNbinsX(); sample_j_bin_index++) {
                    // Get Central Value
                    float j_central = sample_j.fCentralValue->GetBinContent(sample_j_bin_index+1);
                    // calculate covariance
                    float cov_value = 0.; 
                    for (int u = 0; u < fNumAltUnis; u++) {
                        // get variations
                        float i_variation = sample_i.fUniverses[variation][u]->GetBinContent(sample_i_bin_index+1);
                        float j_variation = sample_j.fUniverses[variation][u]->GetBinContent(sample_j_bin_index+1);
                        cov_value += (i_central - i_variation) * (j_central - j_variation);
                    }
                    // average
                    if (fNumAltUnis != 0) cov_value /= fNumAltUnis;

                    // calculate fractional covariance
                    float fcov_value;
                    if (i_central * j_central < 1e-6) {
                        fcov_value = 0.;
                    }
                    else {
		        fcov_value = cov_value / (i_central * j_central);
                    }

                    // save values
                    cov[variation].SetBinContent(cov_i + 1, cov_j + 1, cov_value);
                    fcov[variation].SetBinContent(cov_i + 1, cov_j + 1, fcov_value);
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

    std::string correlation_name = "Correlation " + std::to_string(variation);

    // Pearson Correlation Coefficients
    corr.emplace_back(correlation_name.c_str(), correlation_name.c_str(), num_bins, 0, num_bins, num_bins, 0, num_bins);
    
    for (int i = 0; i < num_bins; i++) {
        float covii = cov[variation].GetBinContent(i+1, i+1);
        for (int j = 0; j < num_bins; j++) {
            float covjj = cov[variation].GetBinContent(j+1, j+1);
            float covij = cov[variation].GetBinContent(i+1, j+1);

            float corrij;
            // handle case where covariance is 0
            if (covii*covjj < 1e-6) {
                corrij = 0.;
            }
            else {
                corrij = covij / TMath::Sqrt(covii*covjj);
            }
            corr[variation].SetBinContent(i+1, j+1, corrij);
        }
    }
    
    // Add bin labels
    std::vector <TH2D*> hists = {&cov[variation], &fcov[variation], &corr[variation]};
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
    
    for (TH2D &h: cov) h.Write();
    for (TH2D &h: fcov) h.Write();
    for (TH2D &h: corr) h.Write();

    // write histos
    if (fSaveCentralValue) {
      for (auto const &sample: fEventSamples) {
        sample.fCentralValue->Write();
      }
    }
    if (fSaveUniverses) {
        for (auto const &sample: fEventSamples) {
            for (unsigned i =0; i < fNVariations; i++) {
              for (TH1D *h: sample.fUniverses[i]) {
                  h->Write();
              }
            }
        }
    }
    // close file
    output->Close(); 
}

// Turn the covariance TH2D into a matrix -- include stat uncertainty
TMatrixDSym Covariance::CovarianceMatrix() {
    TMatrixDSym E_mat(cov[0].GetNbinsX());
    for (int i = 0; i < cov[0].GetNbinsX(); i++) {
        for (int j = 0; j < cov[0].GetNbinsY(); j++) {
            E_mat[i][j] = 0.; // zero-initialize
            for (TH2D &cov_histo: cov)  E_mat[i][j] += cov_histo.GetBinContent(i+1, j+1); // add all systematics
            // add statistical uncertainty
            if (i == j) {
                unsigned stat_ind = i;
                float stat_uncertainty;
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

