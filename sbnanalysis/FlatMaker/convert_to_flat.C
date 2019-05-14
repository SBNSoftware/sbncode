#include "sbnanalysis/core/Event.hh"

#include "FlatRecord.cxx" // don't need to figure out how to build libraries yet...

#include "TFile.h"

#include <iostream>

void convert_to_flat()
{
  TFile* fin = TFile::Open("/pnfs/sbnd/persistent/users/gputnam/numu_simulation_12_05_2018/processed/output_SBNOsc_NumuSelection_Modern_SBND.root");

  TTree* tr = (TTree*)fin->Get("sbnana");

  Event* event = 0;
  tr->SetBranchAddress("events", &event);

  TFile fout("flat.root", "RECREATE");

  TTree* trout = new TTree("rec", "rec");
  // Have trouble with memory usage (because several trees are open at
  // once?). Set the maximum buffer usage (per tree) to 3MB (10x less than
  // default)
  trout->SetAutoFlush(-3*1000*1000);

  flat::FlatEvent* rec = new flat::FlatEvent("rec.", trout, 0);//policy);

  for(int i = 0; i < tr->GetEntries(); ++i){
    // TODO use CAFAna progress bar
    if(i%1000 == 0) std::cout << i/1000 << " / " << tr->GetEntries()/1000 << std::endl;

    tr->GetEntry(i);

    rec->Fill(*event);
    trout->Fill();
  }

  trout->Write();
  delete rec;
}
