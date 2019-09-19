#include "CAFAna/Core/SystRegistry.h"

#include <cxxabi.h>

#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  std::map<std::string, const ISyst*>& SystRegistry::Map()
  {
    // https://isocpp.org/wiki/faq/ctors#static-init-order
    static auto m = new std::map<std::string, const ISyst*>;
    return *m;
  }

  //----------------------------------------------------------------------
  void SystRegistry::Register(const ISyst* s)
  {
    const std::string name = s->ShortName();

    if(Map().count(name)){
      std::cout << "SystRegistry: Warning: systematic '" << name
                << "' registered multiple times." << std::endl
                << "Check you declared it with 'extern' in the .h and "
                << "instantiated it in the .cxx and check for other ISysts "
                << "that happen to have the same name." << std::endl;
      // Store an entry there, so we can detect further duplicates, but make it
      // invalid.
      Map()[name] = 0;
      return;
    }

    Map()[name] = s;
  }

  //----------------------------------------------------------------------
  void SystRegistry::UnRegister(const ISyst* s)
  {
    auto it = Map().find(s->ShortName());

    if(it == Map().end()){
      std::cout << "SystRegistry: Error: unregistering ISyst '"
                << s->ShortName() << "' that was never registered!"
                << std::endl;
      return;
    }

    // If was multiply registered, leave that signal
    if(it->second) Map().erase(it);
  }

  //----------------------------------------------------------------------
  const ISyst* SystRegistry::ShortNameToSyst(const std::string& s,
                                             bool allowFail)
  {
    auto it = Map().find(s);
    if(it == Map().end()){
      if(allowFail) return 0;
      std::cout << "SystRegistry: Error: Syst '" << s << "' not found. "
                << "Pass allowFail=true to return NULL in this case."
                << std::endl;
      abort();
    }

    const ISyst* ret = it->second;

    if(!ret){
      std::cout << "SystRegistry: Error: Syst '" << s
                << "' was registered multiple times. Refusing to return a"
                << " random instance. Go fix the syst registration."
                << std::endl;
      abort();
    }

    return ret;
  }

  //----------------------------------------------------------------------
  void SystRegistry::Print()
  {
    std::cout << Map().size() << " ISysts:" << std::endl;
    for(auto it: Map()){
      std::cout << it.first << " :\t";

      if(it.second){
        std::cout << abi::__cxa_demangle(typeid(*it.second).name(), 0, 0, 0)
                  << " at " << it.second << std::endl;
      }
      else{
        std::cout << "MULTIPLY REGISTERED!" << std::endl;
      }
    }
  }
}
