#include <iostream>
#include <cassert>

#include "chunk.h"

void uscript::Chunk::Write(uint8_t instruction) {
  code.push_back(instruction);
}

unsigned uscript::Chunk::AddConstant(uscript::Value constant) {
  constants.push_back(constant);
  return constants.size() - 1;
}

void uscript::Chunk::Disassemble(const std::string &name) const {
  std::cout << "== " << name << " ==" << std::endl; 
  for (unsigned index = 0; index < code.size();) {
    index = DisassembleInstruction(index);
  }
}

unsigned simpleInstruction(const char *name, unsigned index) {
  std::cout << name << std::endl;
  return index + 1;
}

unsigned jumpInstruction(const char *name, int sign, unsigned index, uint8_t jumpa, uint8_t jumpb) {
  uint16_t jump = ((uint16_t)jumpa) << 8;
  jump |= jumpb;

  std::cout << name << " " << index << " " << (index + 3 + sign * jump) << std::endl;
  return index + 3;
}

unsigned constantInstruction(const char *name, unsigned index, uint8_t constant, const uscript::Value &value) {
  std::cout << name << " " << (unsigned)constant << " "; 
  value.Print();
  std::cout << std::endl;
  return index + 2;
} 

unsigned byteInstruction(const char *name, unsigned index, uint8_t slot) {
  std::cout << name << " " << (unsigned)slot << std::endl;
  return index + 2;
}

unsigned uscript::Chunk::DisassembleInstruction(unsigned index) const {
  std::cout << index << " ";
  uint8_t instruction = code[index];
  switch (instruction) {
    case uscript::OP_RETURN:
      return simpleInstruction("OP_RETURN", index);
    case uscript::OP_PRINT:
      return simpleInstruction("OP_PRINT", index);
    case uscript::OP_FIELDS:
      return simpleInstruction("OP_FIELDS", index);
    case uscript::OP_LENGTH:
      return simpleInstruction("OP_LENGTH", index);
    case uscript::OP_CONSTANT:
      return constantInstruction("OP_CONSTANT", index, code[index+1], constants[code[index+1]]);
    case uscript::OP_ADD:
      return simpleInstruction("OP_ADD", index);
    case uscript::OP_SUBTRACT:
      return simpleInstruction("OP_SUBTRACT", index);
    case uscript::OP_MULTIPLY:
      return simpleInstruction("OP_MULTIPLY", index);
    case uscript::OP_DIVIDE:
      return simpleInstruction("OP_DIVIDE", index);
    case uscript::OP_NEGATE:
      return simpleInstruction("OP_NEGATE", index);
    case uscript::OP_NIL:
      return simpleInstruction("OP_NIL", index);
    case uscript::OP_NOT:
      return simpleInstruction("OP_NOT", index);
    case uscript::OP_TRUE:
      return simpleInstruction("OP_TRUE", index);
    case uscript::OP_FALSE:
      return simpleInstruction("OP_FALSE", index);
    case uscript::OP_POP:
      return simpleInstruction("OP_POP", index);
    case uscript::OP_LOOP:
      return jumpInstruction("OP_LOOP", -1, index, code[index+1], code[index+2]);
    case uscript::OP_JUMP:
      return jumpInstruction("OP_JUMP", 1, index, code[index+1], code[index+2]);
    case uscript::OP_JUMP_IF_FALSE:
      return jumpInstruction("OP_JUMP_IF_FALSE", 1, index, code[index+1], code[index+2]);
    case uscript::OP_GET_PROPERTY:
      return constantInstruction("OP_GET_PROPERTY", index, code[index+1], constants[code[index+1]]);
    case uscript::OP_CALL:
      return byteInstruction("OP_CALL", index, code[index+1]);
    case uscript::OP_INDEX:
      return simpleInstruction("OP_INDEX", index);
    case uscript::OP_SET_LOCAL:
      return byteInstruction("OP_SET_LOCAL", index, code[index+1]);
    case uscript::OP_GET_LOCAL:
      return byteInstruction("OP_GET_LOCAL", index, code[index+1]);
    case uscript::OP_DEFINE_GLOBAL:
      return constantInstruction("OP_DEFINE_GLOBAL", index, code[index+1], constants[code[index+1]]);
    case uscript::OP_GET_GLOBAL:
      return constantInstruction("OP_GET_GLOBAL", index, code[index+1], constants[code[index+1]]);
    case uscript::OP_SET_GLOBAL:
      return constantInstruction("OP_SET_GLOBAL", index, code[index+1], constants[code[index+1]]);
    case uscript::OP_EQUAL:
      return simpleInstruction("OP_EQUAL", index);
    case uscript::OP_GREATER:
      return simpleInstruction("OP_GREATER", index);
    case uscript::OP_LESS:
      return simpleInstruction("OP_LESS", index);
    default:
      std::cout << "Unknown Instruction: " << instruction << std::endl;
      return index + 1;
  }
}

