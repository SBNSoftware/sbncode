#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <string>
#include "core/Event.hh"
#include "core/SubRun.hh"
void read_numu_file(const char *input_file, const char *output_file) {
  // open the input file
  TFile f(input_file);
  // and output file
  TFile fout(output_file, "CREATE");

  TTree *ana_tree = new TTree("sbnana", "SBN Analysis Tree");
  event::Event *event = new event::Event();
  ana_tree->Branch("events", &event);
  
  TTree *reco_tree = new TTree("sbnreco", "SBN Reco Event Tree");
  event::RecoEvent *reco_event = new event::RecoEvent();
  reco_tree->Branch("reco_events", &reco_event);
  
  TTree *subrun_tree = new TTree("sbnsubrun", "SBN Analysis Subrun Tree");
  // subrun_tree->Autosave("overwrite");
  SubRun *subrun = new SubRun();
  subrun_tree->Branch("subruns", &subrun);

  // read from input file
  TTreeReader myReader("sbnana", &f);
  TTreeReaderArray<float> reco_energies(myReader, "reco.reco_energy");
  TTreeReaderArray<float> reco_weights(myReader, "reco.weight");
  TTreeReaderArray<float> reco_had_energy(myReader, "reco.hadronic_energy");
  TTreeReaderArray<float> reco_lep_energy(myReader, "reco.lepton_energy");
  TTreeReaderArray<float> reco_lep_costh(myReader, "reco.lepton_costh");
  TTreeReaderArray<int> reco_npion(myReader, "reco.npion");
  TTreeReaderArray<int> reco_npi0(myReader, "reco.npi0");
  TTreeReaderArray<int> reco_nproton(myReader, "reco.nproton");
  TTreeReaderValue<int> experiment(myReader, "experiment");
  TTreeReaderValue<int> m_run(myReader, "metadata.run");
  TTreeReaderValue<int> m_eventID(myReader, "metadata.eventID");
  TTreeReaderValue<int> m_subrun(myReader, "metadata.subrun");

  int max_event = -1;
  int event_ind = 0;
  
  // Loop over the sbnana Tree
  while (myReader.Next()) {
    // break if we are done
    if (max_event >= 0 && event_ind == max_event) break;
 
    // fill each neutrino enrgy in each interaction
    for (unsigned i = 0; i < reco_energies.GetSize(); i++) {

      reco_event->reco.weight = reco_weights[i];
      reco_event->reco.reco_energy = reco_energies[i];
      reco_event->reco.hadronic_energy = reco_had_energy[i];
      reco_event->reco.lepton_energy = reco_lep_energy[i];
      reco_event->reco.lepton_costh = reco_lep_costh[i];
      reco_event->reco.npion = reco_npion[i];
      reco_event->reco.npi0 = reco_npi0[i];
      reco_event->reco.nproton = reco_nproton[i];
      reco_event->reco.truth_index = -1;
      reco_event->reco.index = i;

      reco_event->metadata.run = *m_run;
      reco_event->metadata.subrun = *m_subrun;
      reco_event->metadata.eventID = *m_eventID;
      reco_event->experiment = (Experiment)*experiment;

      event->reco.push_back(reco_event->reco);
      reco_tree->Fill();
    }

    event->metadata.run = *m_run;
    event->metadata.subrun = *m_subrun;
    event->metadata.eventID = *m_eventID;
    event->experiment = (Experiment)*experiment;
    event->nreco = event->reco.size();

    ana_tree->Fill();
    event->reco.clear();

    event_ind += 1;
  }

  // now loop over the subrun TTree
  TTreeReader readSubRun("sbnsubrun", &f);
  TTreeReaderValue<int> runID(readSubRun, "runID");
  TTreeReaderValue<int> subrunID(readSubRun, "subrunID");
  TTreeReaderValue<float> totpot(readSubRun, "totpot");
  TTreeReaderValue<float> totgoodpot(readSubRun, "totgoodpot");
  TTreeReaderValue<int> totspills(readSubRun, "totspills");
  TTreeReaderValue<int> goodspills(readSubRun, "goodspills");
  while (readSubRun.Next()) {
    subrun->runID = *runID;
    subrun->subrunID = *subrunID;
    subrun->totpot = *totpot;
    subrun->totgoodpot = *totgoodpot;
    subrun->totspills = *totspills;
    subrun->goodspills = *goodspills;
    subrun_tree->Fill();
  }

  fout.cd();
  ana_tree->Write();
  reco_tree->Write();
  subrun_tree->Write();
}
