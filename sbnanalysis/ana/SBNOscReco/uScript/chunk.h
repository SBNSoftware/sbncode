#ifndef uscript_chunk_h
#define uscript_chunk_h

#include <vector>
#include <string>
#include <map>

#include "value.h"
#include "common.h"

namespace uscript {

enum OpCode {
  OP_RETURN,
  OP_CONSTANT,
  OP_NIL,
  OP_TRUE,
  OP_FALSE,
  OP_POP,
  OP_GET_LOCAL,
  OP_SET_LOCAL,
  OP_SET_GLOBAL,
  OP_GET_GLOBAL,
  OP_DEFINE_GLOBAL,
  OP_EQUAL,
  OP_GREATER,
  OP_LESS,
  OP_ADD,     
  OP_SUBTRACT,
  OP_MULTIPLY,
  OP_DIVIDE, 
  OP_NOT,
  OP_NEGATE,
  OP_PRINT,
  OP_FIELDS,
  OP_LENGTH,
  OP_JUMP_IF_FALSE,
  OP_JUMP,
  OP_LOOP,
  OP_CALL,
  OP_INDEX,
  OP_GET_PROPERTY
};

class Chunk {
public:
  unsigned AddConstant(Value constant);
  void Write(uint8_t instruction);
  void Disassemble(const std::string &name) const;

  std::string source;
  std::vector<uint8_t> code;
  std::vector<Value> constants;
  unsigned DisassembleInstruction(unsigned index) const;
};

} // end namespace
#endif
