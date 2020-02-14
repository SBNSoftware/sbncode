// This macro loads a CAF with regular Reflex-based StandardRecord and with
// SRProxy. It recurses through all branches in every event checking for
// differences and printing a message if it finds any. Run me with cafe.

#include "TFile.h"
#include "TTree.h"

#include "StandardRecord/StandardRecord.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "StandardRecord/Proxy/CheckEquals.h"

void check_proxy(std::string fname, int N = -1)
{
  TFile* f = new TFile(fname.c_str());
  assert(!f->IsZombie());
  TTree* recTree = (TTree*)f->Get("recTree");
  assert(recTree);

  // Use an entirely seperate tree to prevent any interference between the two
  // methods.
  TFile* f2 = new TFile(fname.c_str());
  assert(!f2->IsZombie());
  TTree* recTree2 = (TTree*)f2->Get("recTree");
  assert(recTree2);


  caf::StandardRecord* sr = 0;
  recTree->SetBranchAddress("rec", &sr);

  caf::SRProxy srProxy(recTree2, "rec");

  if(N < 0 || N > recTree->GetEntries()) N = recTree->GetEntries();
  for(int i = 0; i < N; ++i){
    std::cout << i << " / " << N << std::endl;

    recTree->GetEntry(i);

    recTree2->LoadTree(i);

    CheckEquals(srProxy, *sr);
  }
}

