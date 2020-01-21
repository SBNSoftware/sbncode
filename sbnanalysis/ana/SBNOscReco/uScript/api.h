#ifndef uscript_api_h
#define uscript_api_h

#include <cassert>
#include <functional>

#include "common.h"
#include "compile.h"
#include "vm.h"

namespace uscript {

// convert any random template type into a char
template<typename T>
using JustAChar = const char *;

template<typename T>
using JustACharRef = const char *&;

void inline SetAll(VM &vm, std::initializer_list<const char *> names) {
  for (const char *name: names) {
    vm.AddGlobal(name);
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
void AddAll(VM &vm, JustACharRef<T> name, JustACharRef<TObjs> ...names, const T *&obj, const TObjs *&...objs) { 
  vm.AddGlobal(name, obj); 
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
std::function<uscript::Value (const TObjs*...)> compile(JustAChar<TObjs> ...names, const char *source) {
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
