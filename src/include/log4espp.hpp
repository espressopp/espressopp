/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/****************************************************************************
*                                                                           *
*  Author      : Thomas Brandes, SCAI, FhG                                  *
*  Copyright   : SCAI, FhG, St. Augustin, Germany                           *
*                MPI-P, Mainz, Germany                                      *
*  Date        : Sep 08                                                     *
*  Last Update : Sep 08                                                     *
*                                                                           *
*  This file is part of the Espresso++ software                             *
*                                                                           *
*  Module      : log4espp.hpp                                               *
*                                                                           *
*  Function    : Macros for logging with log4cxx|log4cpp|default            *
*                                                                           *
****************************************************************************/

#ifndef _LOG4ESPP_HPP
#define _LOG4ESPP_HPP

/************************************************************************
*                                                                       *
*  Compile time guards for LOGGING                                      *
*                                                                       *
*    make sure that the desired logging levels are enabled              *
*                                                                       *
*    LOG4ESPP_DEBUG_ENABLED  :  LOG4ESPP_DEBUG is done                  *
*    LOG4ESPP_INFO_ENABLED   :  LOG4ESPP_INFO is done                   *
*                                                                       *
*  The compile time guards itself can be set by these macros:           *
*                                                                       *
*    LOG4ESPP_LEVEL_TRACE  - compile all                                *
*    LOG4ESPP_LEVEL_DEBUG  - compile debug and higher                   *
*    LOG4ESPP_LEVEL_INFO   - compile info and higher                    *
*    LOG4ESPP_LEVEL_WARN   - compile warn and higher                    *
*    LOG4ESPP_LEVEL_ERROR  - compile error and higher                   *
*    LOG4ESPP_LEVEL_FATAL  - compile fatal only                         *
*                                                                       *
*  Please note: These guards are only for compile time, so logging      *
*               can still be switched off at runtime.                   *
*                                                                       *
************************************************************************/

#if defined(LOG4ESPP_LEVEL_TRACE)

#define LOG4ESPP_TRACE_ENABLED
#define LOG4ESPP_DEBUG_ENABLED
#define LOG4ESPP_INFO_ENABLED

#define LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_FATAL_ENABLED

#elif defined(LOG4ESPP_LEVEL_DEBUG)

#undef LOG4ESPP_TRACE_ENABLED
#define LOG4ESPP_DEBUG_ENABLED
#define LOG4ESPP_INFO_ENABLED
#define LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_FATAL_ENABLED

#elif defined(LOG4ESPP_LEVEL_INFO)

#undef LOG4ESPP_TRACE_ENABLED
#undef LOG4ESPP_DEBUG_ENABLED
#define LOG4ESPP_INFO_ENABLED
#define LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_FATAL_ENABLED

#elif defined(LOG4ESPP_LEVEL_WARN)

#undef LOG4ESPP_TRACE_ENABLED
#undef LOG4ESPP_DEBUG_ENABLED
#undef LOG4ESPP_INFO_ENABLED
#define LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_FATAL_ENABLED

#elif defined(LOG4ESPP_LEVEL_ERROR)

#undef LOG4ESPP_TRACE_ENABLED
#undef LOG4ESPP_DEBUG_ENABLED
#undef LOG4ESPP_INFO_ENABLED
#undef LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_FATAL_ENABLED

#elif defined(LOG4ESPP_LEVEL_FATAL)

#undef LOG4ESPP_TRACE_ENABLED
#undef LOG4ESPP_DEBUG_ENABLED
#undef LOG4ESPP_INFO_ENABLED
#undef LOG4ESPP_WARN_ENABLED
#undef LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_FATAL_ENABLED

#elif defined(LOG4ESPP_LEVEL_OFF)

#undef LOG4ESPP_TRACE_ENABLED
#undef LOG4ESPP_DEBUG_ENABLED
#undef LOG4ESPP_INFO_ENABLED
#undef LOG4ESPP_WARN_ENABLED
#undef LOG4ESPP_ERROR_ENABLED
#undef LOG4ESPP_FATAL_ENABLED

#else

#undef LOG4ESPP_TRACE_ENABLED
#define LOG4ESPP_DEBUG_ENABLED
#define LOG4ESPP_INFO_ENABLED
#define LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_FATAL_ENABLED

#endif

/************************************************************************
*                                                                       *
*   LOG4ESPP   <======    LOG4CPP                                       *
*                                                                       *
************************************************************************/

#if defined(HAVE_LOG4CPP) and defined(LOG4ESPP_USE_LOG4CPP)

#include <cstdio>
#include "log4cpp/Portability.hh"
#ifdef CMAKE_HEADERS
#include <unistd.h>
#endif
#include <iostream>
#include <sstream>
#include "log4cpp/Category.hh"
#include "log4cpp/Appender.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/OstreamAppender.hh"
#ifdef LOG4CPP_HAVE_SYSLOG
#include "log4cpp/SyslogAppender.hh"
#endif
#include "log4cpp/Layout.hh"
#include "log4cpp/BasicLayout.hh"
#include "log4cpp/Priority.hh"
#include "log4cpp/NDC.hh"
#include <log4cpp/SimpleConfigurator.hh>
#include <log4cpp/BasicConfigurator.hh>

  /*******************************************************
  *   LOG4ESPP_CONFIGURE                                 *
  *******************************************************/

#define LOG4ESPP_CONFIGURE() { char *logFile; \
   logFile = getenv("LOG4ESPP"); \
   if (logFile != NULL) {\
      printf ("LOG4ESPP: LOG4CPP property configuration with file %s\n", logFile); \
      log4cpp::SimpleConfigurator::configure(logFile); \
     } \
     else { \
      printf ("LOG4ESPP: LOG4CPP basic configuration\n"); \
      log4cpp::BasicConfigurator::configure(); \
     } \
   }

#define LOG4ESPP_ROOTLOGGER(aLogger) log4cpp::Category& aLogger = log4cpp::Category::getRoot()

#define LOG4ESPP_DECL_LOGGER(aLogger) log4cpp::Category& aLogger
#define LOG4ESPP_LOGGER(aLogger,name) log4cpp::Category& aLogger = log4cpp::Category::getInstance(std::string(name))

  /*******************************************************
  *   LOG4ESPP_XXXXX_ON                                  *
  *******************************************************/

#define LOG4ESPP_TRACE_ON(logger) (logger.isPriorityEnabled(log4cpp::Priority::DEBUG))
#define LOG4ESPP_DEBUG_ON(logger) (logger.isPriorityEnabled(log4cpp::Priority::DEBUG))
#define LOG4ESPP_INFO_ON(logger) (logger.isPriorityEnabled(log4cpp::Priority::INFO))
#define LOG4ESPP_WARN_ON(logger) (logger.isPriorityEnabled(log4cpp::Priority::WARN))
#define LOG4ESPP_ERROR_ON(logger) (logger.isPriorityEnabled(log4cpp::Priority::ERROR))
#define LOG4ESPP_FATAL_ON(logger) (logger.isPriorityEnabled(log4cpp::Priority::FATAL))

  /*******************************************************
  *   LOG4ESPP_TRACE                                     *
  *******************************************************/

#ifdef LOG4ESPP_TRACE_ENABLED
#define LOG4ESPP_TRACE(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::DEBUG)) \
                 { std::ostringstream omsg; omsg << msg; logger.debug(omsg.str()); }}
#else
#define LOG4ESPP_TRACE(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_DEBUG                                     *
  *******************************************************/

#ifdef LOG4ESPP_DEBUG_ENABLED
#define LOG4ESPP_DEBUG(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::DEBUG)) \
                 { std::ostringstream omsg; omsg << msg; logger.debug(omsg.str()); }}
#else
#define LOG4ESPP_DEBUG(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_INFO                                     *
  *******************************************************/

#ifdef LOG4ESPP_INFO_ENABLED
#define LOG4ESPP_INFO(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::INFO)) \
                 { std::ostringstream omsg; omsg << msg; logger.info(omsg.str()); }}
#else
#define LOG4ESPP_INFO(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_WARN                                     *
  *******************************************************/

#ifdef LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_WARN(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::WARN)) \
                 { std::ostringstream omsg; omsg << msg; logger.warn(omsg.str()); }}
#else
#define LOG4ESPP_WARN(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_ERROR                                     *
  *******************************************************/

#ifdef LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_ERROR(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::ERROR)) \
                 { std::ostringstream omsg; omsg << msg; logger.error(omsg.str()); }}
#else
#define LOG4ESPP_ERROR(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_FATAL                                     *
  *******************************************************/

#ifdef LOG4ESPP_FATAL_ENABLED
#define LOG4ESPP_FATAL(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::FATAL)) \
                 { std::ostringstream omsg; omsg << msg; logger.fatal(omsg.str()); }}
#else
#define LOG4ESPP_FATAL(logger,msg)
#endif

#define LOG4ESPP_PUSH(string) log4cpp::NDC::push(string)
#define LOG4ESPP_POP() log4cpp::NDC::pop()

/************************************************************************
*                                                                       *
*   LOG4ESPP   <======    LOG4CXX                                       *
*                                                                       *
************************************************************************/

#elif defined(HAVE_LOG4CXX) and defined(LOG4ESPP_USE_LOG4CXX)

#include <log4cxx/logstring.h>
// #include <stdlib.h>
#include <log4cxx/logger.h>
#include <log4cxx/defaultconfigurator.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/helpers/exception.h>
#include <log4cxx/logmanager.h>
#include <log4cxx/ndc.h>
// #include <locale.h>

using namespace log4cxx;
using namespace log4cxx::helpers;

  /*******************************************************
  *   LOG4ESPP_CONFIGURE                                 *
  *******************************************************/

#define LOG4ESPP_CONFIGURE() { char *logFile; \
   logFile = getenv("LOG4ESPP"); \
   if (logFile != NULL) {\
      PropertyConfigurator::configure(logFile); \
     } \
     else { \
      /* DefaultConfigurator::configure(LogManager::getLoggerRepository()); */ \
      BasicConfigurator::configure(); \
      LoggerPtr rootLogger = Logger::getRootLogger(); \
      rootLogger->setLevel(log4cxx::Level::getWarn()); \
     } \
   }

  /*******************************************************
  *   LOG4ESPP_LOGGER                                    *
  *******************************************************/

#define LOG4ESPP_ROOTLOGGER(aLogger) log4cxx::LoggerPtr aLogger = Logger::getRootLogger()
#define LOG4ESPP_LOGGER(aLogger,name) log4cxx::LoggerPtr aLogger = log4cxx::Logger::getLogger(name)
#define LOG4ESPP_DECL_LOGGER(aLogger) log4cxx::LoggerPtr aLogger
#define LOG4ESPP_SET_LOGGER(aLogger,name) aLogger = log4cxx::Logger::getLogger(name)

  /*******************************************************
  *   LOG4ESPP_XXXXX_ON                                  *
  *******************************************************/

#define LOG4ESPP_TRACE_ON(logger)  (LOG4CXX_UNLIKELY(logger->isTraceEnabled()))
#define LOG4ESPP_DEBUG_ON(logger) (LOG4CXX_UNLIKELY(logger->isDebugEnabled()))
#define LOG4ESPP_INFO_ON(logger)  (LOG4CXX_UNLIKELY(logger->isInfoEnabled()))
#define LOG4ESPP_WARN_ON(logger)  (LOG4CXX_UNLIKELY(logger->isWarnEnabled()))
#define LOG4ESPP_ERROR_ON(logger)  (LOG4CXX_UNLIKELY(logger->isErrorEnabled()))
#define LOG4ESPP_FATAL_ON(logger)  (LOG4CXX_UNLIKELY(logger->isFatalEnabled()))

  /*******************************************************
  *   LOG4ESPP_DEBUG                                     *
  *******************************************************/

#ifdef LOG4ESPP_DEBUG_ENABLED
#define LOG4ESPP_DEBUG(logger,msg) LOG4CXX_DEBUG(logger,msg)
#else
#define LOG4ESPP_DEBUG(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_TRACE                                     *
  *******************************************************/

#ifdef LOG4ESPP_TRACE_ENABLED
#define LOG4ESPP_TRACE(logger,msg) LOG4CXX_TRACE(logger,msg)
#else
#define LOG4ESPP_TRACE(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_WARN                                     *
  *******************************************************/

#ifdef LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_WARN(logger,msg) LOG4CXX_WARN(logger,msg)
#else
#define LOG4ESPP_WARN(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_INFO                                     *
  *******************************************************/

#ifdef LOG4ESPP_INFO_ENABLED
#define LOG4ESPP_INFO(logger,msg) LOG4CXX_INFO(logger,msg)
#else
#define LOG4ESPP_INFO(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_ERROR                                     *
  *******************************************************/

#ifdef LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_ERROR(logger,msg) LOG4CXX_ERROR(logger,msg)
#else
#define LOG4ESPP_ERROR(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_FATAL                                     *
  *******************************************************/

#ifdef LOG4ESPP_FATAL_ENABLED
#define LOG4ESPP_FATAL(logger,msg) LOG4CXX_FATAL(logger,msg)
#else
#define LOG4ESPP_FATAL(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_SET_XXXX                                  *
  *******************************************************/

#define LOG4ESPP_SET_TRACE(logger) logger->setLevel(log4cxx::Level::getTrace())
#define LOG4ESPP_SET_DEBUG(logger) logger->setLevel(log4cxx::Level::getDebug())
#define LOG4ESPP_SET_INFO(logger) logger->setLevel(log4cxx::Level::getInfo())
#define LOG4ESPP_SET_WARN(logger) logger->setLevel(log4cxx::Level::getWarn())
#define LOG4ESPP_SET_ERROR(logger) logger->setLevel(log4cxx::Level::getError())
#define LOG4ESPP_SET_FATAL(logger) logger->setLevel(log4cxx::Level::getFatal())

#define LOG4ESPP_PUSH(string) NDC::push(string)
#define LOG4ESPP_POP()        NDC::pop()

/************************************************************************
*                                                                       *
*   LOG4ESPP   <======   NO LOGGER  ()                                  *
*                                                                       *
************************************************************************/

#elif defined(LOG4ESPP_NONE)

class LogClass {
};

#define LOG4ESPP_CONFIGURE()

  /*******************************************************
  *   LOG4ESPP_DECL_LOGGER(logger)                       *
  *   LOG4ESPP_ROOTLOGGER(logger)                        *
  *   LOG4ESPP_LOGGER(logger,name)                       *
  *******************************************************/

#define LOG4ESPP_ROOTLOGGER(aLogger) LogClass aLogger
#define LOG4ESPP_LOGGER(aLogger,name) LogClass aLogger
#define LOG4ESPP_DECL_LOGGER(aLogger) LogClass aLogger

#define LOG4ESPP_TRACE_ON(logger) (0)
#define LOG4ESPP_DEBUG_ON(logger) (0)
#define LOG4ESPP_INFO_ON(logger)  (0)
#define LOG4ESPP_WARN_ON(logger)  (0)
#define LOG4ESPP_ERROR_ON(logger) (0)
#define LOG4ESPP_FATAL_ON(logger) (0)

#define LOG4ESPP_TRACE(logger,msg)
#define LOG4ESPP_DEBUG(logger,msg)
#define LOG4ESPP_INFO(logger,msg)
#define LOG4ESPP_WARN(logger,msg)
#define LOG4ESPP_ERROR(logger,msg)
#define LOG4ESPP_FATAL(logger,msg)

/************************************************************************
*                                                                       *
*   LOG4ESPP   <======   Generic logger (additional lib required)       *
*                                                                       *
************************************************************************/

#else

#include <sstream>
#include "esutil/Logger.hpp"

  /*******************************************************
  *   LOG4ESPP_CONFIGURE                                 *
  *******************************************************/

#define LOG4ESPP_CONFIGURE() { log4espp::Logger::configure(); }

#define LOG4ESPP_ROOTLOGGER(aLogger) log4espp::Logger& aLogger = log4espp::Logger::getRoot()

#define LOG4ESPP_DECL_LOGGER(aLogger) log4espp::Logger& aLogger
#define LOG4ESPP_LOGGER(aLogger,name) log4espp::Logger& aLogger = log4espp::Logger::getInstance(std::string(name))

  /*******************************************************
  *   LOG4ESPP_XXXXX_ON                                  *
  *******************************************************/

#define LOG4ESPP_TRACE_ON(logger) (logger.isTraceEnabled())
#define LOG4ESPP_DEBUG_ON(logger) (logger.isDebugEnabled())
#define LOG4ESPP_INFO_ON(logger) (logger.isInfoEnabled())
#define LOG4ESPP_WARN_ON(logger) (logger.isWarnEnabled())
#define LOG4ESPP_ERROR_ON(logger) (logger.isErrorEnabled())
#define LOG4ESPP_FATAL_ON(logger) (logger.isFatalEnabled())

  /*******************************************************
  *   LOG4ESPP_TRACE                                     *
  *******************************************************/

#ifdef LOG4ESPP_TRACE_ENABLED
#define LOG4ESPP_TRACE(logger,msg) { if (logger.isTraceEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.trace(LOG4ESPP_LOCATION, omsg.str()); } }
#else
#define LOG4ESPP_TRACE(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_DEBUG                                     *
  *******************************************************/

#ifdef LOG4ESPP_DEBUG_ENABLED
#define LOG4ESPP_DEBUG(logger,msg) { if (logger.isDebugEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.debug(LOG4ESPP_LOCATION, omsg.str()); } }
#else
#define LOG4ESPP_DEBUG(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_INFO                                     *
  *******************************************************/

#ifdef LOG4ESPP_INFO_ENABLED
#define LOG4ESPP_INFO(logger,msg) { if (logger.isInfoEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.info(LOG4ESPP_LOCATION, omsg.str()); } }
#else
#define LOG4ESPP_INFO(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_WARN                                     *
  *******************************************************/

#ifdef LOG4ESPP_WARN_ENABLED
#define LOG4ESPP_WARN(logger,msg) { if (logger.isWarnEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.warn(LOG4ESPP_LOCATION, omsg.str()); } }
#else
#define LOG4ESPP_WARN(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_ERROR                                     *
  *******************************************************/

#ifdef LOG4ESPP_ERROR_ENABLED
#define LOG4ESPP_ERROR(logger,msg) { if (logger.isErrorEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.error(LOG4ESPP_LOCATION, omsg.str()); } }
#else
#define LOG4ESPP_ERROR(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_FATAL                                     *
  *******************************************************/

#ifdef LOG4ESPP_FATAL_ENABLED
#define LOG4ESPP_FATAL(logger,msg) { if (logger.isFatalEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.fatal(LOG4ESPP_LOCATION, omsg.str()); } }
#else
#define LOG4ESPP_FATAL(logger,msg)
#endif

  /*******************************************************
  *   LOG4ESPP_SET_XXXX                                  *
  *******************************************************/

#define LOG4ESPP_SET_TRACE(logger) { logger->setLevel(TRACE, true); }
#define LOG4ESPP_SET_DEBUG(logger) { logger->setLevel(DEBUG, true); }
#define LOG4ESPP_SET_INFO(logger) { logger->setLevel(INFO, true); }
#define LOG4ESPP_SET_WARN(logger) { logger->setLevel(WARN, true); }
#define LOG4ESPP_SET_ERROR(logger) { logger->setLevel(ERROR, true); }
#define LOG4ESPP_SET_FATAL(logger) { logger->setLevel(FATAL, true); }

  /*******************************************************
  *   LOG4ESPP_PUSH / LOG4ESPP_POP                       *
  *******************************************************/

#define LOG4ESPP_PUSH(string) Logger::push(string)
#define LOG4ESPP_POP()        Logger::pop()

#endif

#endif
