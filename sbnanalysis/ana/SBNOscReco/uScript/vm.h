#ifndef uscript_vm_h
#define uscript_vm_h

#include <map>
#include <vector>
#include <stdarg.h>
#include <stdio.h> 

#include "chunk.h"
#include "value.h"

namespace uscript {

enum InterpretResult {
  INTERPRET_OK,
  INTERPRET_COMPILE_ERROR,
  INTERPRET_RUNTIME_ERROR
};

class VM {
  const Chunk *chunk;
  unsigned ip;
  std::vector<Value> stack;
  std::map<const char *, Value> globals;

  uint8_t ReadInstruction();
  Value ReadConstant();
  Value Peek(unsigned distance=0) { return stack[stack.size() - distance - 1]; }
  Value Pop() { Value v = stack[stack.size() - 1]; stack.pop_back(); return v; }
  void Push(Value v) { stack.push_back(v); }
  void Reset();

  void RuntimeError(const char *format, ...);
  
  bool CallValue(Value callee, int argCount);

  bool AccessValue(Value instance, const char *name, Value *result);
  bool GetTField(ObjTInstance instance, const char *name, Value *ret);
  bool CallTMethod(ObjTMethod *method, int argCount);
  void DoAddGlobal(const char *classname, const char *name, uint8_t *data);

public:
  VM();
  InterpretResult Interpret(const char* source);
  InterpretResult Interpret(Chunk *chunk);
  void SetChunk(const Chunk *_chunk);
  InterpretResult Run(Value *ret=NULL);

  template <typename TObj>
  void inline AddGlobal(const char *name, const TObj *object) {
    DoAddGlobal(std::string(type_name<TObj>()).c_str(), name, (uint8_t*)object);
  }
};

} // end namespace
#endif
