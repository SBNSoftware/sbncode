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
void RegisterAll(Chunk &chunk, const char *name, JustAChar<TObjs> ...names);

template<typename... TObjs>
void RegisterAll(VM &vm, JustAChar<TObjs> ...names, const TObjs *...objs); 

template<typename... TObjs>
void RegisterAll(Chunk &chunk) {}

template<typename T, typename... TObjs>
void RegisterAll(Chunk &chunk, const char *name, JustAChar<TObjs> ...names) {
  chunk.Register<T>(name);
  RegisterAll<TObjs...>(chunk, names...);
}


template<typename... TObjs>
void RegisterAll(VM &vm) {}

template<typename T, typename... TObjs>
void RegisterAll(VM &vm, const char *name, JustAChar<TObjs> ...names, const T *obj, const TObjs *...objs) { 
  vm.Register(name, obj); 
  RegisterAll<TObjs...>(vm, names..., objs...);
}


template<typename... TObjs>
Chunk compileChunk(JustAChar<TObjs> ...names, const char *source) {
  Scanner scanner(source);
  Chunk chunk;
  RegisterAll<TObjs...>(chunk, names...); 

  Compiler compiler;
  assert(compiler.Compile(source, &chunk));

  return std::move(chunk);
}

template<typename... TObjs>
std::function<bool (std::reference_wrapper<const TObjs>...)> compile(JustAChar<TObjs> ...names, const char *source) {
  Chunk chunk = compile<TObjs...>(source, names...);

  return [chunk](std::reference_wrapper<const TObjs>... data) {
    VM vm;
    vm.SetChunk(chunk);
    RegisterAll(vm, &data...);
    Value value;
    vm.Run(&value);
    return !!value;
  }; 
}

} // end namespace

#endif
