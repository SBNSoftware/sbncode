// This macro loads a CAF with regular Reflex-based StandardRecord and with
// SRProxy. It recurses through all branches in every event checking for
// differences and printing a message if it finds any. Run me with cafe.

#include "TFile.h"
#include "TTree.h"

#include "StandardRecord/StandardRecord.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include <iostream>

void check_proxy(std::string fname, std::string fname2 = "", int N = -1)
{
  // Default to comparing a file to itself accessed two ways
  if(fname2.empty()) fname2 = fname;

  TFile* f = new TFile(fname.c_str());
  assert(!f->IsZombie());
  TTree* recTree = (TTree*)f->Get("recTree");
  assert(recTree);

  TFile* f2 = TFile::Open(fname2.c_str());
  assert(!f2->IsZombie());

  // This is the logic from SpectrumLoader than also handles multitree flatcafs
  TDirectory* dir = 0;
  TTree* tr = 0;

  TObject* obj = f2->Get("recTree");
  assert(obj);
  if(obj->ClassName() == std::string("TTree")){
    tr = (TTree*)obj;
  }
  else{
    dir = (TDirectory*)obj;
    tr = (TTree*)dir->Get("rec");
    assert(tr);
  }

  const caf::CAFType type = caf::GetCAFType(dir, tr);

  // Read one directly via dictionary
  caf::StandardRecord* sr = 0;
  recTree->SetBranchAddress("rec", &sr);

  // And the other via proxy
  long i;
  caf::Proxy<caf::StandardRecord> srProxy(dir, tr, "rec", i, 0);

  if(N < 0 || N > recTree->GetEntries()) N = recTree->GetEntries();
  for(i = 0; i < N; ++i){
    std::cout << i << " / " << N << std::endl;

    recTree->GetEntry(i);

    if(type != caf::kFlatMultiTree) tr->LoadTree(i);

    srProxy.CheckEquals(*sr);
  }
}

