#include "sbnanalysis/core/Event.hh"

#include "FlatRecord.cxx"

void convert_to_flat()
{
  //  gSystem->Load("srcs/sbncode/sbnanalysis/build/core/libsbnanalysis_Event.so");

  //  TFile* f = new TFile("/sbnd/data/users/bckhouse/processed_1.tempwgh/output_SBNOsc_NumuSelection_Modern_SBND.root");

  TFile* f = new TFile("/pnfs/sbnd/persistent/users/gputnam/numu_simulation_12_05_2018/processed/output_SBNOsc_NumuSelection_Modern_SBND.root");

  TTree* tr = (TTree*)f->Get("sbnana");

  Event* event = 0;
  tr->SetBranchAddress("events", &event);

  //  TH1* hwei = new TH1F("", "", 100, 0, 2);

  for(int i = 0; i < tr->GetEntries(); ++i){
    std::cout << i << " / " << tr->GetEntries() << std::endl;

    tr->GetEntry(i);

    // //    std::cout << event->truth.size() << std::endl;
    // if(!event->truth.empty()){
    //   //      std::cout << event->truth[0].neutrino.pdg << std::endl;
    //   double weiUniv0 = 1;
    //   for(auto it: event->truth[0].weights){
    //     weiUniv0 *= it.second[0];
    //     //        std::cout << it.first << " " << it.second.size() << " " << it.second[0] << std::endl;
    //   }
    //   //      std::cout << weiUniv0 << std::endl;
    //   hwei->Fill(weiUniv0);
    // }
  }

  //  hwei->Draw();
}
