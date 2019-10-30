#include "DynamicSelector.h"
#include <string>

#include "TClassTable.h"
#include "TClass.h"
#include "TDictionary.h"
#include "TDictAttributeMap.h"
#include "TProtoClass.h"

#include <iostream>
#include <typeinfo>

#include "../Data/RecoTrack.h"

void numu::PrintClasses() {
  DictFuncPtr_t reco_track_func = gClassTable->GetDict("numu::RecoTrack");
  std::cout << "track dict func: " << reco_track_func << std::endl;
  TClass *reco_track_class = reco_track_func();

  TDataMember *reco_track_length = reco_track_class->GetDataMember("length");
  std::cout << "data mem: " << reco_track_length << std::endl;
  std::cout << "data type name: " << reco_track_length->GetTypeName() << std::endl; 
  std::cout << "data is basic: " << reco_track_length->IsBasic() << std::endl;
  std::cout << "data is enum: " << reco_track_length->IsEnum() << std::endl;
  int offset = reco_track_class->GetDataMemberOffset("length");
  std::cout << "data off: " << offset << std::endl;

  TDataMember *reco_track_start = reco_track_class->GetDataMember("start");
  std::cout << "data mem: " << reco_track_start << std::endl;
  std::cout << "data type name: " << reco_track_start->GetTypeName() << std::endl; 
  std::cout << "data is basic: " << reco_track_start->IsBasic() << std::endl;
  std::cout << "data is enum: " << reco_track_start->IsEnum() << std::endl;
  offset = reco_track_class->GetDataMemberOffset("start");
  std::cout << "data off: " << offset << std::endl;

  TDataMember *reco_track_start_fx = reco_track_class->GetDataMember("start.fX");
  std::cout << "data mem: " << reco_track_start_fx << std::endl;
  std::cout << "data type name: " << reco_track_start_fx->GetTypeName() << std::endl; 
  std::cout << "data is basic: " << reco_track_start_fx->IsBasic() << std::endl;
  std::cout << "data is enum: " << reco_track_start_fx->IsEnum() << std::endl;
  offset = reco_track_class->GetDataMemberOffset("start.fX");
  std::cout << "data off: " << offset << std::endl;

  TDataMember *reco_track_wenter = reco_track_class->GetDataMember("wall_enter");
  std::cout << "data mem: " << reco_track_wenter << std::endl;
  std::cout << "data type name: " << reco_track_wenter->GetTypeName() << std::endl; 
  std::cout << "data is basic: " << reco_track_wenter->IsBasic() << std::endl;
  std::cout << "data is enum: " << reco_track_wenter->IsEnum() << std::endl;
  offset = reco_track_class->GetDataMemberOffset("wall_enter");
  std::cout << "data off: " << offset << std::endl;

  TDataMember *reco_track_garbo = reco_track_class->GetDataMember("garbo");
  std::cout << "data mem: " << reco_track_garbo << std::endl;
}

struct LiteralValue {
  bool is_float;
  bool is_valid;
  int data_int;
  float data_num;
};

enum ROOTValueClass {
  RVTrack, RVTrueTrack, RVEvent
};

enum ROOTValueType {
  RVBool, RVFloat, RVInteger, RVUnsigned,
};

struct ROOTValue {
  ROOTValueClass base_class;
  int offset;
  int size;
  ROOTValueType type;
};

struct Value {
  bool is_literal;
  LiteralValue literal;
  ROOTValue root;
};

enum BoolOp {
 BOand, BOor
};

enum Comp {
  Ceq, Cleq, Cgeq, Clt, Cgt
};


enum TokenType {
  TVar,
  TLiteral,
  TComp,
  TBoolOp
};

struct Token {
  TokenType type;
  std::string str;
};

struct Expression {
  Token lhs;
  Token comp;
  Token rhs;
  Token boolop;
};

enum TokenizerState {
  StartExpression,
  GetComparison,
  EndExpression,
  NextExpression
};


std::vector<Expression> FailTokenize(const std::string &cutstr, unsigned ind) {
  std::cerr << "Failed to parse cutstr (" << cutstr << ") at index (" << ind << ")\n";
  return {};
}

bool StartsVar(char c) {
  return isalpha(c);
}

int ConsumeVar(const std::string &str, unsigned index, Token &token) {
  unsigned delta = 0;
  char c = str[index];
  while ((isalpha(c) || isdigit(c) || c == '_' || c == '.') && index + delta < str.size()) {
    delta ++;
    c = str[index + delta];
  }
  token.str = str.substr(index, delta);
  token.type = TVar;
  return delta;
}

bool StartsLiteral(char c) {
  return isdigit(c) || c == '-';
}

int ConsumeLiteral(const std::string &str, unsigned index, Token &token) {
  unsigned delta = 0;
  char c = str[index];
  while ((isdigit(c) || c == '.' || c == '-') && index + delta < str.size()) {
    delta ++;
    c = str[index + delta];
  }
  token.str = str.substr(index, delta);
  token.type = TLiteral;
  return delta;
}

int ConsumeComparison(const std::string &str, unsigned index, Token &token) {
  std::string one = str.substr(index, 1);
  std::string two = index+1 < str.size() ? str.substr(index, 2) : std::string();

  if (two == "==" || two == "<=" || two == ">=") {
    token.str = two;
    token.type = TComp;
    return 2;
  }
  if (one == "<" || one == ">") {
    token.str = one;
    token.type = TComp;
    return 1;
  }
  return -1;
}

int ConsumeBoolOp(const std::string &str, unsigned index, Token &token) {
  std::string two = index+1 < str.size() ? str.substr(index, 2) : std::string();
  if (two == "||" || two == "&&") {
    token.str = two;
    token.type = TBoolOp;
    return 2;
  }
  return -1;
}

std::vector<Expression> Tokenize(const std::string &cutstr) {
  unsigned len = cutstr.size();
  TokenizerState state = StartExpression;

  std::vector<Expression> ret;
  unsigned expression_index = 0;

  unsigned i = 0;
  while (i < len) {
    char c = cutstr[i];
    int delta;
    switch (state) {
      case StartExpression: 
        // make a new expression
        ret.emplace_back();
        expression_index = ret.size() - 1;

        if (StartsLiteral(c)) {
          i += ConsumeLiteral(cutstr, i, ret[expression_index].lhs);
        }
        else if (StartsVar(c)) {
          i += ConsumeVar(cutstr, i, ret[expression_index].lhs);
        }
        else {
          return FailTokenize(cutstr, i);
        }
        state = GetComparison;
        break;
      case EndExpression: // expression -- must be literal is var
        if (StartsLiteral(c)) {
          i += ConsumeLiteral(cutstr, i, ret[expression_index].rhs);
        }
        else if (StartsVar(c)) {
          i += ConsumeVar(cutstr, i, ret[expression_index].rhs);
        }
        else {
          return FailTokenize(cutstr, i);
        }
        state = NextExpression;
        break;
      case GetComparison: // get the comparison operation
        delta = ConsumeComparison(cutstr, i, ret[expression_index].comp);
        if (delta < 0) return FailTokenize(cutstr, i);
        i += delta;
        state = EndExpression; 
        break;
      case NextExpression: // get the boolean operation
        delta = ConsumeBoolOp(cutstr, i, ret[expression_index].boolop);
        if (delta < 0) return FailTokenize(cutstr, i);
        i += delta;
        state = StartExpression; 
        break;
    }
  }

  if (state != NextExpression) {
    return FailTokenize(cutstr, i);
  }

  return ret;
}

BoolOp Token2BoolOp(const Token &token) {
  assert(token.type == TBoolOp);
  if (token.str == "||") return BOor;
  else if (token.str == "&&") return BOand;
  assert(false);
}

Comp Token2Comp(const Token &token) {
  if (token.str == "==") return Ceq;
  else if (token.str == "<=") return Cleq;
  else if (token.str == ">=") return Cgeq;
  else if (token.str == "<") return Clt;
  else if (token.str == ">") return Cgt;
  assert(false);
}

LiteralValue MakeLiteral(const Value &value, const numu::RecoTrack &track, const numu::RecoEvent &event) { 
  LiteralValue ret;
  if (value.is_literal) {
    ret = value.literal;
  }
  else {
    ret.is_valid = true;
    const uint8_t *dataptr = NULL;
    switch (value.root.base_class) {
      case RVTrack:
        dataptr = ((const uint8_t *)&track) + value.root.offset;
        break;
      case RVTrueTrack:
        if (track.match.has_match) {
          dataptr = ((const uint8_t *)&event.true_tracks.at(track.match.mcparticle_id)) + value.root.offset;
        }
        break;
      case RVEvent:
        dataptr = ((const uint8_t *)&event) + value.root.offset;
        break;
    }

    if (dataptr == NULL) {
      ret.is_valid = false;
      return ret;
    }

    switch (value.root.type) {
      case RVBool:
        ret.is_float = false;
        ret.data_int = *dataptr & 1;
        break;
      case RVFloat:
        ret.is_float = true;
        ret.data_num = *(const float *)dataptr;
        break;
      case RVInteger:
        ret.is_float = false;
        ret.data_int = *(const int *)dataptr;
        break;
      case RVUnsigned: 
        assert(value.root.size == 4);
        ret.is_float = false; 
        ret.data_int = *(const unsigned *)dataptr;
        break;
    }
  }
  return ret;
}

Value Token2Val(const Token &token) {
  Value ret;
  assert(token.type == TLiteral || token.type == TVar);  
  // Compile a literal
  if (token.type == TLiteral) {
    size_t ind = 0;
    int data_int = std::stoi(token.str, &ind); 

    ret.is_literal = true;
    ret.literal.is_valid = true;

    if (ind != token.str.size()) {
      float data_num = std::stof(token.str, &ind);
      ret.literal.data_num = data_num;
      ret.literal.is_float = true;
    }
    else {
      ret.literal.data_int = data_int;
      ret.literal.is_float = false;
    }
  }
  // Compile a ROOT object
  else {
    ret.is_literal = false;
    TDataMember *data = NULL;
    if (strncmp(token.str.c_str(), "track.", 6) == 0) {
      ret.root.base_class = RVTrack; 
      TClass *reco_class = gClassTable->GetDict("numu::RecoTrack")(); 
      data = reco_class->GetDataMember(token.str.substr(6).c_str());
    }
    else if (strncmp(token.str.c_str(), "event.", 6) == 0) {
      ret.root.base_class = RVEvent;
      TClass *reco_class = gClassTable->GetDict("numu::RecoEvent")(); 
      data = reco_class->GetDataMember(token.str.substr(6).c_str());
    }
    else if (strncmp(token.str.c_str(), "true_track.", 11) == 0) {
      ret.root.base_class = RVTrueTrack;
      TClass *reco_class = gClassTable->GetDict("numu::RecoTrack")(); 
      data = reco_class->GetDataMember(token.str.substr(11).c_str());
    }
    else assert(false);

    ret.root.offset = data->GetOffset();
    ret.root.size = data->GetUnitSize();

    std::string type = data->GetTypeName();
    bool is_enum = data->IsEnum();
    if (type == "float") {
      ret.root.type = RVFloat;
    }
    else if (type == "unsigned" || is_enum) {
      ret.root.type = RVUnsigned;
    }
    else if (type == "int") {
      ret.root.type = RVInteger;
    }
    else if (type == "bool") {
      ret.root.type = RVBool;
    }
    else assert(false);
  }

  return ret;
}

template <typename T>
bool DoComp(Comp comp, const T &lhs, const T &rhs) {
  switch (comp) {
    case Ceq:   
      return lhs == rhs;
    case Cleq:
      return lhs <= rhs;
    case Cgeq:
      return lhs >= rhs;
    case Clt:
      return lhs < rhs;
    case Cgt:
      return lhs > rhs;
  }
  return false;
}

std::function<bool (const numu::RecoTrack &, const numu::RecoEvent &)> FailCompile() {
  std::cerr << "Failed to compile\n";
  return [](const numu::RecoTrack &, const numu::RecoEvent &) { return false; };
}

std::function<bool (const numu::RecoTrack &, const numu::RecoEvent &)> numu::Compile(const std::string &cutstr) {
  if (cutstr.size() == 0) {
    return [](const numu::RecoTrack &track, const numu::RecoEvent &event) { return true; };
  }

  std::vector<std::function<bool (const numu::RecoTrack &, const numu::RecoEvent &)>> closures;
  std::vector<BoolOp> boolops;

  std::vector<Expression> expressions = Tokenize(cutstr);
  if (!expressions.size()) return FailCompile();

  // turn each expression into its own closure
  for (unsigned i = 0; i < expressions.size(); i++) {
    const Expression &exp = expressions[i];
    Value lhs = Token2Val(exp.lhs);
    Value rhs = Token2Val(exp.rhs);
    Comp comp = Token2Comp(exp.comp);

    closures.push_back([lhs, rhs, comp](const numu::RecoTrack &track, const numu::RecoEvent &event) {
      LiteralValue lit_lhs = MakeLiteral(lhs, track, event); 
      LiteralValue lit_rhs = MakeLiteral(rhs, track, event); 
      if (!lit_lhs.is_valid || !lit_rhs.is_valid) return false; 
      if (lit_lhs.is_float) {
        assert(lit_rhs.is_float);
        return DoComp(comp, lit_lhs.data_num, lit_rhs.data_num);
      }
      assert(!lit_rhs.is_float);
      return DoComp(comp, lit_lhs.data_int, lit_rhs.data_int);
    });
 
    if (i+1 < expressions.size()) { 
      BoolOp op = Token2BoolOp(exp.boolop);
      boolops.push_back(op);
    }
  }
  return [closures, boolops](const numu::RecoTrack &track, const numu::RecoEvent &event) {
    for (unsigned i = 0; i < closures.size(); i++) {
      bool this_ret = closures[i](track, event);
      if (i == boolops.size()) return this_ret;
      if (this_ret && boolops[i] == BOor) return true;
      if (!this_ret && boolops[i] == BOand) return false;
    }
  };
}



