#ifndef uscript_api_h
#define uscript_api_h

#include <cassert>
#include <functional>
#include <string>

#include "common.h"
#include "compile.h"
#include "vm.h"

namespace uscript {

// convert any random template type into a char
template<typename T>
using JustAString = std::string;

template<typename T>
using JustAStringRef = const std::string&;

void inline SetAll(VM &vm, std::initializer_list<std::string> names) {
  for (const std::string &name: names) {
    vm.AddGlobal(name.c_str());
  }
}

template <typename... Ts>
void _swallow(Ts&&...) {}

template<typename... TObjs>
void RegisterAll() {
  _swallow(uscript::Compiler::Register<TObjs>...);
}

template<typename... TObjs>
void AddAll(VM &vm) {}

template<typename T, typename... TObjs>
void AddAll(VM &vm, JustAStringRef<T> name, JustAStringRef<TObjs> ...names, const T *&obj, const TObjs *&...objs) { 
  vm.AddGlobal(name.c_str(), obj); 
  AddAll<TObjs...>(vm, names..., objs...);
}


template<typename... TObjs>
Chunk compileChunk(const char *source) {
  Chunk chunk;
  RegisterAll<TObjs...>();

  uscript::Compiler::Compile(source, &chunk);

  return std::move(chunk);
}

template<typename... TObjs>
std::function<uscript::Value (const TObjs*...)> compile(JustAString<TObjs> ...names, const char *source) {
  Chunk chunk = compileChunk<TObjs...>(source);

  VM vm;
  SetAll(vm, {names...}); 
  return [chunk, vm, names...](const TObjs *... data) mutable {
    vm.SetChunk(&chunk);
    AddAll<TObjs...>(vm, names..., data...);
    Value value;
    vm.Run(&value);
    return value;
  }; 
}

} // end namespace

#endif
