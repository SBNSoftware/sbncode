#include "FlatMaker/BranchPolicy.h"

#include "TTree.h"
#include <string>

// Helper functions for FlatRecord.cxx

namespace
{
  //----------------------------------------------------------------------
  TTree* make_tree(const std::string& name, const std::string& title, TTree* model)
  {
    TTree* ret = new TTree(name.c_str(), title.c_str());
    ret->SetAutoSave(model->GetAutoSave());
    ret->SetAutoFlush(model->GetAutoFlush());
    return ret;
  }

  //----------------------------------------------------------------------
  inline char root_code(char){return 'B';}
  inline char root_code(unsigned char){return 'b';}
  inline char root_code(short){return 'S';}
  inline char root_code(short unsigned int){return 's';}
  inline char root_code(int){return 'I';}
  inline char root_code(unsigned int){return 'i';}
  inline char root_code(float){return 'F';}
  inline char root_code(double){return 'D';}
  inline char root_code(long){return 'L';}
  inline char root_code(long unsigned int){return 'l';}
  inline char root_code(bool){return 'O';}

  //----------------------------------------------------------------------
  template<class T> void
  branch(TTree* tr, const std::string& name, T* target,
         const flat::IBranchPolicy* policy)
  {
    if(policy && !policy->Include(name)) return;
    tr->Branch(name.c_str(), target, (name+"/"+root_code(*target)).c_str());
  }

  //----------------------------------------------------------------------
  void branch(TTree* tr, const std::string& name, std::string* target,
              const flat::IBranchPolicy* policy)
  {
    if(policy && !policy->Include(name)) return;
    // we really don't want the string to ever reallocate its storage
    target->reserve(1024*1024);
    tr->Branch(name.c_str(), (void*)target->c_str(), (name+"/C").c_str());
  }
}
