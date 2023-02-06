#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void setpotweight(std::string &fName, const int &nScale) {
  
  // Open the input file and get the subruntree
  TFile f(fName.c_str(),"UPDATE");
  TTree *t = static_cast<TTree*>(f.Get("sbnsubrun"));

  // Make sure there is only 1 entry in the file
  if(t->GetEntries() != 1){
    std::cerr << "Error: Expected a single sunrun entry but got " << t->GetEntries() << std::endl;
    std::exit(1);
  }
  std::cout << " Old POT: " << t->GetBranch("totpot")->GetLeaf("totpot")->GetValue(0) << std::endl;

  // Now access the subrun branches for reading, scaling and writing
  double totpot, totgoodpot;
  t->SetBranchAddress("totpot", &totpot);
  t->SetBranchAddress("totgoodpot", &totgoodpot);

  // Since we only have a single entry, don't need a loop
  t->GetEntry(0);

  totpot = totpot*nScale;
  totgoodpot = totgoodpot*nScale;

  std::cout << " New POT: " << t->GetBranch("totpot")->GetLeaf("totpot")->GetValue(0) << std::endl;
  t->Write();
  f.Close();

}
