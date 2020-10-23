/**
 * \file GeoAlgoException.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class GeoAlgoException
 *
 * @author kazu 
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOALGOEXCEPTION_H
#define BASICTOOL_GEOALGOEXCEPTION_H

#include <iostream>
#include <exception>

namespace geoalgo {
/**
   \class GeoAlgoException
   User defined class GeoAlgoException ... these comments are used to generate
   doxygen documentation!
 */
    
 class GeoAlgoException : public std::exception{

  public:
        
    GeoAlgoException(std::string msg="") : std::exception()
    { 
      _msg  = "\n\033[93m<<EXCEPTION>>\033[00m\033[95m ";
      _msg += msg;
      _msg += "\033[00m\n";
    }
    
    virtual ~GeoAlgoException() throw(){};
    virtual const char* what() const throw() 
    {return _msg.c_str(); }

  private:
    
    std::string _msg;
 };

}

#endif
/** @} */ // end of doxygen group 

