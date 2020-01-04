#ifndef uscript_vm_h
#define uscript_vm_h

#include <map>
#include <vector>
#include <stdarg.h>
#include <stdio.h> 

#include "chunk.h"
#include "value.h"
#include "object.h"

namespace uscript {

enum InterpretResult {
  INTERPRET_OK,
  INTERPRET_COMPILE_ERROR,
  INTERPRET_RUNTIME_ERROR
};

class VM {
  Chunk chunk;
  unsigned ip;
  std::vector<Value> stack;
  std::vector<Obj *> objects;
  std::map<std::string, Value> globals;

  uint8_t ReadInstruction();
  Value ReadConstant();
  Value Peek(unsigned distance=0) { return stack[stack.size() - distance - 1]; }
  Value Pop() { Value v = stack[stack.size() - 1]; stack.pop_back(); return v; }
  void Push(Value v) { stack.push_back(v); }
  void Reset();

  void RuntimeError(const char *format, ...);
  
  bool CallValue(Value callee, int argCount);

  bool AccessValue(Value instance, const std::string &name, Value *result);
  bool GetTField(ObjTInstance *instance, const std::string &name, Value *ret);
  bool CallTMethod(ObjTMethod *method, int argCount);
  void DoRegister(const char *classname, const char *name, uint8_t *data);

public:
  VM();
  InterpretResult Interpret(const char* source);
  InterpretResult Interpret(Chunk chunk);
  void SetChunk(Chunk _chunk);
  InterpretResult Run(Value *ret=NULL);

  template <typename TObj>
  void inline Register(const char *name, const TObj *object) {
    DoRegister(std::string(type_name<TObj>()).c_str(), name, (uint8_t*)object);
  }
};

} // end namespace
#endif
