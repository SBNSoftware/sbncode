#include <iostream>
#include <cassert>

#include "vm.h"
#include "common.h"
#include "compile.h"

#ifdef USCRIPT_TIMING
#include <chrono>
using namespace std::chrono;
#endif

uscript::InterpretResult uscript::VM::Interpret(Chunk *chunk) {
  ip = 0;
  SetChunk(chunk);
  return Run();
}

uscript::InterpretResult uscript::VM::Interpret(const char* source) {
  ip = 0;
  Chunk chunk;
  if (!uscript::Compiler::Compile(source, &chunk)) {
    return uscript::INTERPRET_COMPILE_ERROR;
  }

  SetChunk(&chunk);

  return Run();
}


uscript::VM::VM():
  ip(0)
  {}

void uscript::VM::SetChunk(const Chunk *_chunk) {
  chunk = _chunk;
  ip = 0;
}

void uscript::VM::Reset() {
  stack.clear();
  ip = 0;
}

uint8_t uscript::VM::ReadInstruction() {
  return chunk->code[ip++];
}

uscript::Value uscript::VM::ReadConstant() {
  return chunk->constants[ReadInstruction()];
}

static bool valuesEqual(uscript::Value a, uscript::Value b) {
  if (a.val != b.val) return false;

  switch (a.val) {
    case uscript::VAL_BOOL: return AS_BOOL(a) == AS_BOOL(b);
    case uscript::VAL_NIL:  return true;
    case uscript::VAL_NUMBER: return AS_NUMBER(a) == AS_NUMBER(b);
    case uscript::VAL_INTEGER: return AS_INTEGER(a) == AS_INTEGER(b);
    case uscript::VAL_OBJ_STRING: {
      return AS_CSTRING(a) == AS_CSTRING(b);
    }
    default:
      return false;
  }

}

uscript::InterpretResult uscript::VM::Run(Value *ret) {
#ifdef USCRIPT_TIMING
#define IMPL_USCRIPT_TIMING \
  auto stop = high_resolution_clock::now(); \
  __uscript_global_vmtime.fetch_add(duration_cast<nanoseconds>(stop - start).count());
#else
#define IMPL_USCRIPT_TIMING
#endif

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
      IMPL_USCRIPT_TIMING \
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
      IMPL_USCRIPT_TIMING \
      return uscript::INTERPRET_RUNTIME_ERROR; \
    } \
  } while (false)

#define READ_STRING() AS_CSTRING(ReadConstant())
#define READ_SHORT() (ip += 2, (uint16_t)((chunk->code[ip-2] << 8) | chunk->code[ip-1])) 

#ifdef USCRIPT_TIMING
  auto start = high_resolution_clock::now();
#endif

  while (1) {


#ifdef DEBUG_TRACE_EXECUTION
    std::cout << "    ";
    for (const uscript::Value &val: stack) {
      std::cout << "[ "; 
      val.Print();
      std::cout << " ]";
    }
    std::cout << std::endl;
    chunk->DisassembleInstruction(ip);
#endif
    uint8_t instruction = ReadInstruction();
    switch (instruction) {
      case uscript::OP_GET_PROPERTY: {
        Value instance = Pop();
        const char *name = READ_STRING();
        Value result;
        if (!AccessValue(instance, name, &result)) {
          IMPL_USCRIPT_TIMING
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        Push(result);
        break;
      }
      case uscript::OP_INDEX: {
        Value index = Pop();
        if (!IS_INTEGER(index)) {
          RuntimeError("Cannot index with non-integer.");
          IMPL_USCRIPT_TIMING
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        int ind = AS_INTEGER(index);
        Value callee = Pop();
        if (!IndexValue(callee, ind)) {
          IMPL_USCRIPT_TIMING
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        break;
      }
      case uscript::OP_CALL: {
        int argCount = ReadInstruction();
        if (!CallValue(Peek(argCount), argCount)) {
          IMPL_USCRIPT_TIMING
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
        if (ret) *ret = Pop();
        Reset(); // reset the vm state for next time
        IMPL_USCRIPT_TIMING
        return uscript::INTERPRET_OK;
      }
      case uscript::OP_PRINT: {
        Pop().Print();
        std::cout << std::endl;
        break;
      }
      case uscript::OP_LENGTH: {
        Value val = Pop();
        if (IS_TINSTANCE(val)) {
          ObjTInstance inst = AS_TINSTANCE(val);
          Value length = INTEGER_VAL(inst.data.len);
          Push(length);
          break;
        }
        else {
          RuntimeError("Cannot access length of non-tinstance.");   
          IMPL_USCRIPT_TIMING
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
      }
      case uscript::OP_FIELDS: {
        Value val = Pop();
        if (IS_TINSTANCE(val)) {
          ObjTInstance inst = AS_TINSTANCE(val);
          if (inst.data.info) {
             TClassInfo *info = inst.data.info;
             for (auto const &field: info->fields) {
               std::cout << field.first << std::endl;
             }
          }
          
          int length = inst.data.len;
          if (length >= 0) {
            std::cout << "Indexable with length: " << length << std::endl;
          }
          break;
        }
        else {
          RuntimeError("Cannot access fields of non-tinstance.");   
          IMPL_USCRIPT_TIMING
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
      }
      case uscript::OP_ADD: BINARY_OP(+); break;
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
        const char *name = READ_STRING();
        if (globals.count(name)) {
          Push(globals.at(name));
        }
        else {
          RuntimeError("Undefined variable '%s'", name);
          IMPL_USCRIPT_TIMING
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        break;
      }
      case uscript::OP_SET_GLOBAL: {
        const char *name = READ_STRING();
        if (!globals.count(name)) {
          RuntimeError("Undefined variable '%s'", name);
          IMPL_USCRIPT_TIMING
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
          IMPL_USCRIPT_TIMING
          return uscript::INTERPRET_RUNTIME_ERROR;
        }
        break;
      }
    }
  }
#undef IMPL_USCRIPT_TIMING
#undef READ_STRING
#undef READ_SHORT
#undef BINARY_OP
#undef COMP_OP
}

bool uscript::VM::IndexValue(Value callee, int index) {
  if (IS_TINSTANCE(callee)) {
    uscript::ObjTInstance inst = AS_TINSTANCE(callee);
    if (index >= 0  && index < inst.data.len) {
      // indexed value itself has no length and is not a vector
      TData field_data = inst.data;
      field_data.len = -1; 
      Value ret = GetTValue(inst.loc + index * inst.data.Size(), field_data);
      Push(ret);
      return true;
    } 
    else {
      RuntimeError("Cannot index with value (%i) into instance of size (%i)", index, inst.data.len);
      return false;
    }
  }
  return false;
}

bool uscript::VM::CallValue(Value callee, int argCount) {
  RuntimeError("Functions are not implemented.");
  return false;
}

uscript::Value uscript::VM::GetTValue(uint8_t *loc, uscript::TData data) {
  // case where this is a vector
  if (data.info != NULL && data.info->is_vec) {
    // go to the location of the start of the vector
    // FIXME: this is probably implementation defined -- needs to be fixed

    // get the contents of the vector
    uint8_t **vec = (uint8_t**)loc;
    uint8_t *start = vec[0];
    uint8_t *end = vec[1];

    // the data is now what the vector was pointing to
    data = data.info->vec_data;

    // update the size
    data.len = (end - start) / data.Size();

    // update the loc
    loc = start;
  }
  // all non-zero size Tdata's are instances
  if (data.len >= 0 || data.type == uscript::FIELD_TINSTANCE) {
    ObjTInstance inst;
    inst.loc = loc;
    inst.data = data;
    return TINSTANCE_VAL(inst);
  }
  Value ret;
  switch(data.type) {
    case uscript::FIELD_BOOL:
      ret = BOOL_VAL((bool)*loc);
      break;
    case uscript::FIELD_INT:
      ret = INTEGER_VAL(*((int*)loc));
      break;
    case uscript::FIELD_ENUM:
      ret = INTEGER_VAL((int)*((uint32_t*)loc));
      break;
    case uscript::FIELD_UNSIGNED:
      ret = INTEGER_VAL((int)*((unsigned*)loc));
      break;
    case uscript::FIELD_FLOAT:
      ret = NUMBER_VAL((double)*((float*)loc));
      break;
    case uscript::FIELD_DOUBLE:
      ret = NUMBER_VAL(*((double*)loc));
      break;
    case uscript::FIELD_TINSTANCE: 
      break; // unreachable
  }
  return ret;
}

bool uscript::VM::AccessValue(Value instance, const char *name, Value *result) {
  if (IS_TINSTANCE(instance)) {
    bool success = GetTField(AS_TINSTANCE(instance), name, result);
    if (!success) {
      uscript::ObjTInstance &inst = AS_TINSTANCE(instance);
      uscript::TClassInfo *classinfo = inst.data.info;
      RuntimeError("Could not find value %s.", name);
      //RuntimeError("Could not find value %s in class %s.", name, classinfo->name.c_str());
    }
    return success;
  }
  RuntimeError("Cannot access on non-TInstance.");
  return false;
   
}

bool uscript::VM::GetTField(uscript::ObjTInstance instance, const char *name, uscript::Value *ret) {
  if (!instance.loc) return false; // shouldn't happen
  if (!instance.data.info) return false; // possible for list of built-in type (e.g. vector<int> or double[4])

  // if this is a list, then you have to index it before accessing a value
  if (instance.data.len >= 0) return false;

  uscript::TClassInfo *classinfo = instance.data.info;
  if (classinfo->fields.count(name)) {
    uscript::TField field = classinfo->fields.at(name);
    *ret = GetTValue(instance.loc + field.offset, field.data); 
    return true;
  }

  return false;

}

void uscript::VM::DoAddGlobal(const char *classname, const char *name, uint8_t *loc) {
  uscript::TClassInfo *classinfo = uscript::Compiler::GetClassInfo(classname);
  ObjTInstance inst;
  inst.loc = loc;
  inst.data.info = classinfo;
  inst.data.info = classinfo;
  inst.data.len = -1;
  inst.data.type = uscript::FIELD_TINSTANCE;
  // intern the name
  name = uscript::Compiler::Intern(name);
  globals[name] = TINSTANCE_VAL(inst);
}

void uscript::VM::RuntimeError(const char *format, ...) { 
  va_list args;
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
  fprintf(stderr, "\nIn source: %s\n", chunk->source.c_str());

  Reset();
}

template<>
void uscript::VM::AddGlobal<int>(const char *name, const int *obj) {
  name = uscript::Compiler::Intern(name);
  globals[name] = INTEGER_VAL(*obj);
}

template<>
void uscript::VM::AddGlobal<unsigned>(const char *name, const unsigned *obj) {
  name = uscript::Compiler::Intern(name);
  globals[name] = INTEGER_VAL((int)*obj);
}

template<>
void uscript::VM::AddGlobal<float>(const char *name, const float *obj) {
  name = uscript::Compiler::Intern(name);
  globals[name] = NUMBER_VAL(*obj);
}

template<>
void uscript::VM::AddGlobal<double>(const char *name, const double *obj) {
  name = uscript::Compiler::Intern(name);
  globals[name] = NUMBER_VAL(*obj);
}

template<>
void uscript::VM::AddGlobal<bool>(const char *name, const bool *obj) {
  name = uscript::Compiler::Intern(name);
  globals[name] = BOOL_VAL(*obj);
}

void uscript::VM::AddGlobal(const char *name) {
  name = uscript::Compiler::Intern(name);
  globals[name] = NIL_VAL;
}

