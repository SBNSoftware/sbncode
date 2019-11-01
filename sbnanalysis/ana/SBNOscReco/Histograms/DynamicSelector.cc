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

  TDataMember *reco_track_match = reco_track_class->GetDataMember("match");
  std::cout << "data mem: " << reco_track_match << std::endl;
  std::cout << "data type name: " << reco_track_match->GetTypeName() << std::endl;
  std::cout << "data is basic: " << reco_track_match->IsBasic() << std::endl;
  std::cout << "data is enum: " << reco_track_match->IsEnum() << std::endl;
  std::cout << "data off: " << reco_track_class->GetDataMemberOffset("match");
  std::cout << "data size: " << reco_track_match->GetUnitSize() << std::endl; 

  TDataMember *reco_track_start_fx = reco_track_class->GetDataMember("start.fX");
  std::cout << "data mem: " << reco_track_start_fx << std::endl;
  offset = reco_track_class->GetDataMemberOffset("start.fX");
  std::cout << "data off: " << offset << std::endl;

  TDataMember *reco_track_wenter = reco_track_class->GetDataMember("wall_enter");
  std::cout << "data mem: " << reco_track_wenter << std::endl;
  std::cout << "data type name: " << reco_track_wenter->GetTypeName() << std::endl; 
  std::cout << "data is basic: " << reco_track_wenter->IsBasic() << std::endl;
  std::cout << "data is enum: " << reco_track_wenter->IsEnum() << std::endl;
  offset = reco_track_class->GetDataMemberOffset("wall_enter");
  std::cout << "data off: " << offset << std::endl;
  std::cout << "data size: " << reco_track_wenter->GetUnitSize() << std::endl; 

  TDataMember *reco_track_garbo = reco_track_class->GetDataMember("garbo");
  std::cout << "data mem: " << reco_track_garbo << std::endl;
}

struct Value {
  bool is_literal;
  numu::LiteralValue literal;
  numu::ROOTValue root;
};

enum BoolOp {
 BOand, BOor
};

enum Comp {
  Ceq, Cneq, Cleq, Cgeq, Clt, Cgt
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
  int start_depth;
  int end_depth;
};

struct CompiledExpression {
  numu::TrackSelector closure;
  BoolOp boolop;
  int start_depth;
  int end_depth;
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

  if (two == "==" || two == "<=" || two == ">=" || two == "!=") {
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

bool IsOpenParenthesis(char c) {
  return c == '(';
}

bool IsCloseParenthesis(char c) {
  return c == ')';
}

int ConsumeParenthesis(const std::string &str, unsigned index, int *depth) {
  int delta = 0;
  int ret = 0;
  if (str[index] == ')') {
    delta = 1;
    ret = 1;
  }
  else if (str[index] == '(') {
    delta = -1;
    ret = 1;
  }
  if (depth != NULL) *depth += delta;
  return ret;
}

std::vector<Expression> Tokenize(const std::string &cutstr) {
  unsigned len = cutstr.size();
  TokenizerState state = StartExpression;

  std::vector<Expression> ret;

  // zero-sized string -- empty expression list
  if (len == 0) return ret;

  unsigned expression_index = 0;

  unsigned i = 0;
  int depth = 0;
  while (i < len) {
    char c = cutstr[i];
    int delta;
    switch (state) {
      case StartExpression: 
        if (StartsLiteral(c)) {
          // make a new expression
          ret.emplace_back();
          expression_index = ret.size() - 1;
          i += ConsumeLiteral(cutstr, i, ret[expression_index].lhs);
          state = GetComparison;
          ret[expression_index].start_depth = depth;
        }
        else if (StartsVar(c)) {
          // make a new expression
          ret.emplace_back();
          expression_index = ret.size() - 1;
          i += ConsumeVar(cutstr, i, ret[expression_index].lhs);
          state = GetComparison;
          ret[expression_index].start_depth = depth;
        }
        else if (IsOpenParenthesis(c)) {
          i += ConsumeParenthesis(cutstr, i, &depth);
          if (depth > 0) return FailTokenize(cutstr, i-1);
        }
        else {
          return FailTokenize(cutstr, i);
        }
        break;
      case EndExpression: // expression -- must be literal or var
        if (StartsLiteral(c)) {
          i += ConsumeLiteral(cutstr, i, ret[expression_index].rhs);
          ret[expression_index].end_depth = depth;
          state = NextExpression;
        }
        else if (StartsVar(c)) {
          i += ConsumeVar(cutstr, i, ret[expression_index].rhs);
          ret[expression_index].end_depth = depth;
          state = NextExpression;
        }
        else {
          return FailTokenize(cutstr, i);
        }
        break;
      case GetComparison: // get the comparison operation
        delta = ConsumeComparison(cutstr, i, ret[expression_index].comp);
        if (delta < 0) return FailTokenize(cutstr, i);
        i += delta;
        state = EndExpression; 
        break;
      case NextExpression: // get the boolean operation
        // parenthesis closing the previous expression
        if (IsCloseParenthesis(c)) {
          i += ConsumeParenthesis(cutstr, i, &depth);
          if (depth > 0) return FailTokenize(cutstr, i-1);
          ret[expression_index].end_depth = depth; // update its end-depth
        }
        else {
          delta = ConsumeBoolOp(cutstr, i, ret[expression_index].boolop);
          if (delta < 0) return FailTokenize(cutstr, i);
          i += delta;
          state = StartExpression; 
        }
        break;
    }
  }

  if (depth != 0) {
    return FailTokenize(cutstr, i);
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
  else if (token.str == "!=") return Cneq;
  else if (token.str == "<=") return Cleq;
  else if (token.str == ">=") return Cgeq;
  else if (token.str == "<") return Clt;
  else if (token.str == ">") return Cgt;
  assert(false);
}

numu::LiteralValue numu::MakeROOTLiteral(const numu::ROOTValue &value, const numu::RecoTrack &track, const numu::RecoEvent &event) {
  numu::LiteralValue ret;
  ret.is_valid = true;
  const uint8_t *dataptr = NULL;
  switch (value.base_class) {
    case RVTrack:
      dataptr = ((const uint8_t *)&track) + value.offset;
      break;
    case RVTrueTrack:
      if (track.match.has_match) {
        dataptr = ((const uint8_t *)&event.true_tracks.at(track.match.mcparticle_id)) + value.offset;
      }
      break;
    case RVEvent:
      dataptr = ((const uint8_t *)&event) + value.offset;
      break;
  }

  if (dataptr == NULL) {
    ret.is_valid = false;
    return ret;
  }

  switch (value.type) {
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
      assert(value.size == 4);
      ret.is_float = false; 
      ret.data_int = *(const unsigned *)dataptr;
      break;
  }
  return ret;
}

numu::LiteralValue MakeLiteral(const Value &value, const numu::RecoTrack &track, const numu::RecoEvent &event) { 
  if (value.is_literal) {
    return value.literal;
  }
  else {
    return MakeROOTLiteral(value.root, track, event);
  }
}

numu::ROOTValue numu::MakeROOTValue(const std::string &name) {
  numu::ROOTValue ret;

  TDataMember *data = NULL;
  TClass *reco_class = NULL;
  std::vector<std::string> split_string {""};
  for (char s: name) {
    if (s == '.') {
      split_string.push_back("");
    }
    else {
      split_string[split_string.size()-1] += s;
    }
  } 
  if (split_string[0] == "track") {
    ret.base_class = RVTrack; 
    reco_class = gClassTable->GetDict("numu::RecoTrack")(); 
  }
  else if (split_string[0] == "event") {
    ret.base_class = RVEvent;
    reco_class = gClassTable->GetDict("numu::RecoEvent")(); 
  }
  else if (split_string[0] == "true_track") {
    ret.base_class = RVTrueTrack;
    reco_class = gClassTable->GetDict("numu::RecoTrack")(); 
  }
  else assert(false);

  unsigned offset = 0;
  for (unsigned i = 1; i < split_string.size(); i++) {
    data = reco_class->GetDataMember(split_string[i].c_str());
    // end of loop, should be in basic type
    if (i+1 == split_string.size()) assert(data->IsBasic() || data->IsEnum());
    // otherwise, should not be in basic type
    else assert(!(data->IsBasic() || data->IsEnum()));

    // update the offset
    offset += data->GetOffset();
    // not end of loop, get the next class to recurse to
    if (i+1 !=  split_string.size()) {
      reco_class = gClassTable->GetDict(data->GetTypeName())();
    }
  }

  ret.offset = offset;
  ret.size = data->GetUnitSize();

  std::string type = data->GetTypeName();
  bool is_enum = data->IsEnum();
  if (type == "float") {
    ret.type = RVFloat;
  }
  else if (type == "unsigned" || is_enum) {
    ret.type = RVUnsigned;
  }
  else if (type == "int") {
    ret.type = RVInteger;
  }
  else if (type == "bool") {
    ret.type = RVBool;
  }
  else assert(false);

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
    ret.root = numu::MakeROOTValue(token.str);
  }

  return ret;
}

template <typename T>
bool DoComp(Comp comp, const T &lhs, const T &rhs) {
  switch (comp) {
    case Ceq:   
      return lhs == rhs;
    case Cneq:
      return lhs != rhs;
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
  std::vector<Expression> expressions = Tokenize(cutstr);
  if (!expressions.size()) return FailCompile();

  std::vector<CompiledExpression> cexpressions;

  // Compile each expression and build AST 
  for (unsigned i = 0; i < expressions.size(); i++) {
    cexpressions.emplace_back();
    const Expression &exp = expressions[i];
    Value lhs = Token2Val(exp.lhs);
    Value rhs = Token2Val(exp.rhs);
    Comp comp = Token2Comp(exp.comp);

    cexpressions[i].start_depth = exp.start_depth;
    cexpressions[i].end_depth = exp.end_depth;

    cexpressions[i].closure = [lhs, rhs, comp](const numu::RecoTrack &track, const numu::RecoEvent &event) {
      numu::LiteralValue lit_lhs = MakeLiteral(lhs, track, event); 
      numu::LiteralValue lit_rhs = MakeLiteral(rhs, track, event); 
      if (!lit_lhs.is_valid || !lit_rhs.is_valid) return false; 
      if (lit_lhs.is_float) {
        assert(lit_rhs.is_float);
        return DoComp(comp, lit_lhs.data_num, lit_rhs.data_num);
      }
      assert(!lit_rhs.is_float);
      return DoComp(comp, lit_lhs.data_int, lit_rhs.data_int);
    };
 
    if (i+1 < expressions.size()) { 
      cexpressions[i].boolop = Token2BoolOp(exp.boolop);
    }
  }

  // combine into a single closure
  return [cexpressions](const numu::RecoTrack &track, const numu::RecoEvent &event) {
    // empty expression -- always true 
    if (cexpressions.size() == 0) return true;

    unsigned i = 0;
    bool ret = cexpressions[i].closure(track, event);
    while (i+1 < cexpressions.size()) {
      // try to short circuit
      if ((ret && cexpressions[i].boolop == BOor) || (!ret && cexpressions[i].boolop == BOand)) {
        int check_depth = cexpressions[i].end_depth;  
        do {
          i += 1;
        } while (i < cexpressions.size() && check_depth >= cexpressions[i].end_depth);
      }
      // otherwise, evaluate the next expression
      else {
        i += 1;
        ret = cexpressions[i].closure(track, event);
      }
    }
    // return the last value after short-circuiting
    return ret;
  };
}

template <typename T>
void Increment(std::vector<unsigned> &indices, const std::vector<std::vector<T>> &lists) {
  for (unsigned i = 0; i < indices.size(); i++) {
    indices[i] += 1;
    if (indices[i] == lists[i].size()) {
      indices[i] = 0;
    }
    else {
      break;
    }
  }
}

bool NotZero(const std::vector<unsigned> &indices) {
  for (unsigned i = 0; i < indices.size(); i++) {
    if (indices[i] != 0) return true;
  }
  return false;
}

std::vector<std::string> numu::MultiplyNames(const std::vector<std::vector<std::string>> &strings) {
  std::vector<std::string> ret;
  std::vector<unsigned> indices (strings.size(), 0);
  do {
    std::string this_str;
    for (unsigned i = 0; i < indices.size(); i++) {
      this_str += strings[i][indices[i]] + "_";
    }
    ret.push_back(this_str);
    Increment(indices, strings);
  } while (NotZero(indices));
  return ret;
}

std::vector<numu::TrackSelector> numu::MultiplySelectors(const std::vector<std::vector<std::string>> &track_selector_strings) {
  std::vector<numu::TrackSelector> ret;
  std::vector<std::vector<numu::TrackSelector>> track_selectors;
  for (unsigned i = 0; i < track_selector_strings.size(); i++) {
    track_selectors.emplace_back();
    for (unsigned j = 0; j < track_selector_strings[i].size(); j++) {
      track_selectors[i].push_back(numu::Compile(track_selector_strings[i][j]));
    }
  }
  std::vector<unsigned> indices (track_selectors.size(), 0);
  do {
    std::vector<numu::TrackSelector> selectors;
    for (unsigned i = 0; i < indices.size(); i++) {
      selectors.push_back(track_selectors[i][indices[i]]);
    }
    ret.push_back(
      [selectors](const numu::RecoTrack &track, const numu::RecoEvent &event) {
        for (const numu::TrackSelector &selector: selectors) {
          if (!selector(track, event)) return false;
        }
        return true;
      });
    Increment(indices, track_selectors);
  } while (NotZero(indices));
  return ret;
}
