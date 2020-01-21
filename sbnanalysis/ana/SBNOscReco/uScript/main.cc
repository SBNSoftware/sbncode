#include <stdio.h> 
#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "common.h"
#include "chunk.h"
#include "vm.h"
#include "compile.h"

#include "api.h"

#include "../Data/RecoEvent.h"
#include "../Data/FlatInteraction.h"

static void repl() {
  numu::RecoEvent event;
  event.type = (numu::MCType)1;
  numu::CRTHit hit;
  hit.pes = -5.6;
  event.in_time_crt_hits.push_back(hit);

  numu::RecoSlice slice;
  slice.flash_match.time = 2.;

  numu::TrueParticle particle;

  numu::flat::FlatInteraction flat;
  flat.ptrack.start[1] = 5.;

  uscript::Compiler::Register<numu::TrueParticle>();
  uscript::Compiler::Register<numu::RecoEvent>();
  uscript::Compiler::Register<numu::RecoSlice>();
  uscript::Compiler::Register<numu::flat::FlatInteraction>();

  char line[1024];
  uscript::VM vm;
  while (1) {
    std::cout << "> ";
    if (!fgets(line, sizeof(line), stdin)) {
      std::cout << "\n";
      break;
    }
    uscript::Chunk chunk = uscript::compileChunk(line);
    if (uscript::Compiler::HadError()) continue;
    vm.SetChunk(&chunk);
    vm.AddGlobal("flat", &flat);
    vm.AddGlobal("event", &event);
    vm.AddGlobal("slice", &slice);
    vm.AddGlobal("particle", &particle);
    int global = 5;
    vm.AddGlobal("int", &global);
    uscript::Value v;
    vm.Run(&v);
    if (!IS_NIL(v)) {
      v.Print();
      std::cout << std::endl;
    }
  }
}

static char* readFile(const char* path) {                        
  FILE* file = fopen(path, "rb");
  if (file == NULL) {
    std::cerr << "Could not open file: " << path << std::endl; 
    return NULL;
  }

  fseek(file, 0L, SEEK_END); 
  size_t fileSize = ftell(file); 
  rewind(file); 

  char* buffer = (char*)malloc(fileSize + 1); 
  if (buffer == NULL) {
    std::cerr << "Not enough memory to read file\n";
    return NULL;
  }
  size_t bytesRead = fread(buffer, sizeof(char), fileSize, file);
  if (bytesRead < fileSize) {
    std::cerr << "Could not read file.\n";
    free(buffer);
    return NULL;
  }
  buffer[bytesRead] = '\0';

  fclose(file);
  return buffer; 
}   

static int runFile(const char* path) { 
  char* source = readFile(path);            
  if (source == NULL) return 60;
  uscript::VM vm;
  uscript::InterpretResult result = vm.Interpret(source);
  free(source);

  if (result == uscript::INTERPRET_COMPILE_ERROR) return 65;
  if (result == uscript::INTERPRET_RUNTIME_ERROR) return 70;
} 

int main(int argc, const char* argv[]) {
  if (argc == 1) {
    repl();
  }
  else if (argc == 2) {
    return runFile(argv[1]);
  }
  else {
    std::cout << "Usage: uscript [path]\n";
    return 64;
  }
  return 0;
}


