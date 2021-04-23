/**
 * \file LoggerFeature.h
 *
 * \ingroup Base
 * 
 * \brief Class definition file of LoggerFeature
 *
 * @author Kazu - Nevis 2015
 */

/** \addtogroup Base

    @{*/

#ifndef __OPT0FINDERLOGGERFEATURE_H__
#define __OPT0FINDERLOGGERFEATURE_H__

#include <vector>
#include "OpT0FinderLogger.h"

namespace flashmatch {
    
  /**
    \class LoggerFeature
    Framework base class equipped with a logger class
  */
  class LoggerFeature {
    
  public:
    
    /// Default constructor
    LoggerFeature(const std::string logger_name="LoggerFeature")
      : _logger(nullptr)
    { _logger = &(::flashmatch::logger::get(logger_name)); }
    
    /// Default copy constructor
    LoggerFeature(const LoggerFeature &original) : _logger(original._logger) {}
    
    /// Default destructor
    virtual ~LoggerFeature(){};
    
    /// Logger getter
    inline const flashmatch::logger& logger() const
    { return *_logger; }
    
    /// Verbosity level
    void set_verbosity(::flashmatch::msg::Level_t level)
    { _logger->set(level); }

    /// Name getter, defined in a logger instance attribute
    const std::string& name() const
    { return logger().name(); }
    
  private:
    
    flashmatch::logger *_logger;   ///< logger
    
  };
}
#endif

/** @} */ // end of doxygen group
