#include <iostream>
#include <cassert>

#include "vm.h"
#include "common.h"
#include "compile.h"

uscript::InterpretResult uscript::VM::Interpret(Chunk chunk) {
  ip = 0;
  SetChunk(std::move(chunk));
  return Run();
}

uscript::InterpretResult uscript::VM::Interpret(const char* source) {
  ip = 0;
  Chunk chunk;
  Compiler compiler;
  if (!compiler.Compile(source, &chunk)) {
    return uscript::INTERPRET_COMPILE_ERROR;
  }

  SetChunk(std::move(chunk));

  return Run();
}


uscript::VM::VM():
  ip(0)
  {}

void uscript::VM::SetChunk(Chunk _chunk) {
  chunk = _chunk;
  ip = 0;
}

void uscript::VM::Reset() {
  stack.clear();
}

uint8_t uscript::VM::ReadInstruction() {
  return chunk.code[ip++];
}

uscript::Value uscript::VM::ReadConstant() {
  return chunk.constants[ReadInstruction()];
}

static bool valuesEqual(uscript::Value a, uscript::Value b) {
  if (a.val != b.val) return false;

  switch (a.val) {
    case uscript::VAL_BOOL: return AS_BOOL(a) == AS_BOOL(b);
    case uscript::VAL_NIL:  return true;
    case uscript::VAL_NUMBER: return AS_NUMBER(a) == AS_NUMBER(b);
    case uscript::VAL_INTEGER: return AS_INTEGER(a) == AS_INTEGER(b);
    case uscript::VAL_OBJ: {
      return AS_STRING(a) == AS_STRING(b);
    }
  }

}

uscript::InterpretResult uscript::VM::Run(Value *ret) {
#define BINARY_OP(op) \
  do { \
    if (IS_NUMBER(Peek(0)) && IS_NUMBER(Peek(1))) { \
      double b = AS_NUMBER(Pop()); \
      double a = AS_NUMBER(Pop()); \
      Push(NUMBER_VAL(a op b)); \
    } \
    else if (IS_INTEGER(Peek(0)) && IS_INTEGER(Peek(1))) { \
      int b = AS_INTEGER(Pop()); \
      int a = AS_INTEGER(Pop()); \
      Push(INTEGER_VAL(a op b)); \
    } \
    else { \
      RuntimeError("Operands must be both integers or both numbers."); \
      return uscript::INTERPRET_RUNTIME_ERROR; \
    } \
  } while (false)

#define COMP_OP(op) \
  do { \
    if (IS_NUMBER(Peek(0)) && IS_NUMBER(Peek(1))) { \
      double b = AS_NUMBER(Pop()); \
      double a = AS_NUMBER(Pop()); \
      Push(BOOL_VAL(a op b)); \
    } \
    else if (IS_INTEGER(Peek(0)) && IS_INTEGER(Peek(1))) { \
      int b = AS_INTEGER(Pop()); \
      int a = AS_INTEGER(Pop()); \
      Push(BOOL_VAL(a op b)); \
    } \
    else { \
      RuntimeError("Operands must be both integers or both numbers."); \
      return uscript::INTERPRET_RUNTIME_ERROR; \
    } \
  } while (false)

  while (1) {

#define READ_STRING() AS_STRING(ReadConstant())
#define READ_SHORT() (ip += 2, (uint16_t)((chunk.code[ip-2] << 8) | chunk.code[ip-1])) 

#ifdef DEBUG_TRACE_EXECUTION
    std::cout << "    ";
    for (const uscript::Value &val: stack) {
      std::cout << "[ "; 
      val.Print();
      std::cout << " ]";
    }
    std::cout << std::endl;
    chunk.DisassembleInstruction(ip);
#endif
    uint8_t instruction = ReadInstruction();
    switch (instruction) {
      case uscript::OP_GET_PROPERTY: {
        Value instance = Pop();
        const std::string &name = READ_STRING();
        Value result;
        if (!AccessValue(instance, name, &result)) {
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        Push(result);
        break;
      }
      case uscript::OP_CALL: {
        int argCount = ReadInstruction();
        if (!CallValue(Peek(argCount), argCount)) {
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        break;
      }
      case uscript::OP_LOOP: {
        uint16_t offset = READ_SHORT();
        ip -= offset;
        break;
      }
      case uscript::OP_JUMP: {
        uint16_t offset = READ_SHORT();
        ip += offset;
        break;
      }
      case uscript::OP_JUMP_IF_FALSE: {
        uint16_t offset = READ_SHORT();
        if (!Peek(0)) ip += offset;
        break;
      }
      case uscript::OP_RETURN: {
        // if (ret) *ret = Pop();
        return uscript::INTERPRET_OK;
      }
      case uscript::OP_PRINT: {
        Pop().Print();
        std::cout << std::endl;
        break;
      }
      case uscript::OP_ADD: {
        if (IS_STRING(Peek(0)) && IS_STRING(Peek(1))) {
          std::string b = AS_STRING(Pop());
          std::string a = AS_STRING(Pop());
          Obj *obj = new ObjString(a+b);
          objects.push_back(obj);
          Push(OBJ_VAL(obj));
        }
        else  BINARY_OP(+);
        break;
      }
      case uscript::OP_SUBTRACT: BINARY_OP(-); break;
      case uscript::OP_MULTIPLY: BINARY_OP(*); break;
      case uscript::OP_DIVIDE:   BINARY_OP(/); break;
      case uscript::OP_CONSTANT: {
        Value constant = ReadConstant();
        Push(constant);
        break;
      }
      case uscript::OP_NOT: Push(BOOL_VAL(!Pop())); break;
      case uscript::OP_NIL: Push(NIL_VAL); break;
      case uscript::OP_TRUE: Push(BOOL_VAL(true)); break;
      case uscript::OP_FALSE: Push(BOOL_VAL(false)); break;
      case uscript::OP_POP: Pop(); break;
      case uscript::OP_GET_LOCAL: {
        uint8_t slot = ReadInstruction();
        Push(stack[slot]);
        break;
      }
      case uscript::OP_SET_LOCAL: {
        uint8_t slot = ReadInstruction();
        stack[slot] = Peek();
        break;
      }
      case uscript::OP_GET_GLOBAL: {
        const std::string &name = READ_STRING();
        if (globals.count(name)) {
          Push(globals.at(name));
        }
        else {
          RuntimeError("Undefined variable '%s'", name.c_str());
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        break;
      }
      case uscript::OP_SET_GLOBAL: {
        const std::string &name = READ_STRING();
        if (!globals.count(name)) {
          RuntimeError("Undefined variable '%s'", name.c_str());
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        globals[name] = Peek();
        break;
      }
      case uscript::OP_DEFINE_GLOBAL: {
        globals[READ_STRING()] = Pop();
        break;
      }
      case uscript::OP_EQUAL: {
        uscript::Value b = Pop();
        uscript::Value a = Pop();
        Push(BOOL_VAL(valuesEqual(a, b)));
        break;
      }
      case uscript::OP_GREATER: COMP_OP(>); break; 
      case uscript::OP_LESS:    COMP_OP(<); break;
      case uscript::OP_NEGATE: {
        if (IS_NUMBER(Peek())) {
          Push(NUMBER_VAL(-AS_NUMBER(Pop())));
        }
        else if (IS_INTEGER(Peek())) {
          Push(INTEGER_VAL(-AS_INTEGER(Pop())));
        }
        else {
          RuntimeError("Operand must be a number.");
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        break;
      }
    }
  }
#undef READ_STRING
#undef READ_SHORT
#undef BINARY_OP
#undef COMP_OP
}

bool uscript::VM::CallValue(Value callee, int argCount) {
  if (IS_OBJ(callee)) {
    switch (OBJ_TYPE(callee)) {
      case OBJ_TMETHOD:
        return CallTMethod(AS_TMETHOD(callee), argCount);
      default:
        break;
    }
  }
  RuntimeError("Cannot call on non-TMethod.");
  return false;
}

bool uscript::VM::CallTMethod(uscript::ObjTMethod *method, int argCount) {
  Value result;
  // bool success = method->Call(argCount, &stack[stack.size() - argCount], &result);
  bool success = true;
  // clean the stack
  stack.resize(stack.size() - argCount - 1);
  Push(result);
  return success;
}

bool uscript::VM::AccessValue(Value instance, const std::string &name, Value *result) {
  std::cout << "Accessing value\n";
  if (IS_OBJ(instance)) {
    switch (OBJ_TYPE(instance)) {
      case uscript::OBJ_TINSTANCE:
        return GetTField(AS_TINSTANCE(instance), name, result);
      default:
        break;
    }
  }
  RuntimeError("Cannot access on non-TInstance.");
  return false;
   
}

bool uscript::VM::GetTField(uscript::ObjTInstance *instance, const std::string &name, uscript::Value *ret) {
  if (!instance->data) return false;

  uscript::TClassInfo &classinfo = chunk.tclasslist.classes[instance->tclassIndex];
  if (classinfo.fields.count(name)) {
    uscript::TField field = classinfo.fields.at(name);
    switch(field.type) {
      case uscript::FIELD_BOOL:
        *ret = BOOL_VAL((bool)instance->data[field.offset]);
        break;
      case uscript::FIELD_INT:
        *ret = INTEGER_VAL(*((int*)(instance->data+field.offset)));
        break;
      case uscript::FIELD_UNSIGNED:
        *ret = INTEGER_VAL((int)*((unsigned*)(instance->data+field.offset)));
        break;
      case uscript::FIELD_FLOAT:
        *ret = NUMBER_VAL((double)*((float*)(instance->data+field.offset)));
        break;
      case uscript::FIELD_DOUBLE:
        *ret = NUMBER_VAL(*((double*)(instance->data+field.offset)));
        break;
      case uscript::FIELD_TINSTANCE: {
        ObjTInstance *inst = new ObjTInstance;
        inst->type = uscript::OBJ_TINSTANCE;
        inst->data = instance->data + field.offset;
        inst->tclassIndex = field.tclassIndex;
        *ret = OBJ_VAL(inst);
        break;
      }
    }
    return true;
  }

  return false;

}

void uscript::VM::DoRegister(const char *classname, const char *name, uint8_t *data) {

  int index = chunk.tglobals.at(name);
  const uscript::TClassInfo &classinfo = chunk.tclasslist.classes[index];
  assert(classinfo.name == classname);
  // make a global
  ObjTInstance *obj = new ObjTInstance; 
  obj->type = uscript::OBJ_TINSTANCE;
  obj->data = data;
  obj->tclassIndex = index; 
  
  globals[name] = OBJ_VAL(obj);
}

void uscript::VM::RuntimeError(const char *format, ...) { 
  va_list args;
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
  fputs("\n", stderr);

  Reset();
}

