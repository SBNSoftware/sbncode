#ifndef OPT0FINDER_EREXCEPTION_CXX
#define OPT0FINDER_EREXCEPTION_CXX

#include "OpT0FinderException.h"

namespace flashmatch {

  OpT0FinderException::OpT0FinderException(const std::string& msg)
    : std::exception()
  {
      _msg = "\033[93m EXCEPTION \033[00m\033[95m";
      _msg += msg;
      _msg += "\033[00m\n";
  }

  const char* OpT0FinderException::what() const throw()
  { return _msg.c_str(); }

}
#endif
