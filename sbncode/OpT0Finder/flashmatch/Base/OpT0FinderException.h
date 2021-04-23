/**
 * \file OpT0FinderException.h
 *
 * \ingroup OpT0Finder
 * 
 * \brief Class def header for exception classes in OpT0Finder package
 *
 * @author kazuhiro
 */

/** \addtogroup OpT0Finder

    @{*/
#ifndef OPT0FINDER_EXCEPTION_H
#define OPT0FINDER_EXCEPTION_H

#include <iostream>
#include <exception>

namespace flashmatch {
  /**
     \class OpT0FinderException
     Generic (base) exception class for OpT0Finder package
  */
  class OpT0FinderException : public std::exception{

  public:
    /// Default ctor w/ error message input (optional)
    OpT0FinderException(const std::string& msg="");

    virtual ~OpT0FinderException() throw(){};

    virtual const char* what() const throw();

  private:
    /// Error message
    std::string _msg;
  };

}
#endif
/** @} */ // end of doxygen group 

