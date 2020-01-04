#ifndef uscript_api_h
#define uscript_api_h

#include <cassert>
#include <functional>

#include "common.h"
#include "chunk.h"
#include "vm.h"

namespace uscript {

// convert any random template type into a char
template<typename T>
using JustAChar = const char *;

template<typename... TObjs>
void AddAll(VM &vm, JustAChar<TObjs> ...names, const TObjs *...objs); 

template<typename none = void>
void RegisterAll() {}

template<typename T, typename... TObjs>
void RegisterAll() {
  Compiler::Register<T>();
  RegisterAll<TObjs...>();
}


template<typename... TObjs>
void AddAll(VM &vm) {}

template<typename T, typename... TObjs>
void AddAll(VM &vm, const char *name, JustAChar<TObjs> ...names, const T *obj, const TObjs *...objs) { 
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
std::function<bool (std::reference_wrapper<const TObjs>...)> compile(JustAChar<TObjs> ...names, const char *source) {
  Chunk chunk = compile<TObjs...>(source);

  return [chunk](std::reference_wrapper<const TObjs>... data) {
    VM vm;
    vm.SetChunk(&chunk);
    AddAll(vm, &data...);
    Value value;
    vm.Run(&value);
    return !!value;
  }; 
}

} // end namespace

#endif
