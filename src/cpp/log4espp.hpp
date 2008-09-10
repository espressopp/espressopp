
#include "acconfig.hpp"

#if defined(HAVE_LOG4CPP) and defined(LOG4ESPP_USE_LOG4CPP)

#include <stdio.h>
#include "log4cpp/Portability.hh"
#ifdef LOG4CPP_HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <iostream>
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

#define LOG4ESPP_CONFIGURE() log4cpp::SimpleConfigurator::configure("log4cpp.init")
 
#define LOG4ESPP_ROOTLOGGER(aLogger) log4cpp::Category& aLogger = log4cpp::Category::getRoot()

#define LOG4ESPP_DEF_LOGGER(aLogger) log4cpp::Category& aLogger;
#define LOG4ESPP_LOGGER(aLogger,name) log4cpp::Category& aLogger = log4cpp::Category::getInstance(std::string(name))

#if defined(LOG4ESSP_LEVEL_FATAL)
#define LOG4ESPP_DEBUG(logger,msg) 
#define LOG4ESPP_INFO(logger,msg) 
#define LOG4ESPP_WARN(logger,msg) 
#define LOG4ESPP_ERROR(logger,msg)
#define LOG4ESPP_FATAL(logger,msg)
#else
#define LOG4ESPP_DEBUG(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::DEBUG)) \
                 { ostringstream omsg; omsg << msg; logger.debug(omsg.str()); }}
#define LOG4ESPP_INFO(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::INFO)) \
                 { ostringstream omsg; omsg << msg; logger.info(omsg.str()); }}
#define LOG4ESPP_WARN(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::WARN)) \
                 { ostringstream omsg; omsg << msg; logger.info(omsg.str()); }}
#define LOG4ESPP_ERROR(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::ERROR)) \
                 { ostringstream omsg; omsg << msg; logger.info(omsg.str()); }}
#define LOG4ESPP_FATAL(logger,msg) { if (logger.isPriorityEnabled(log4cpp::Priority::FATAL)) \
                 { ostringstream omsg; omsg << msg; logger.info(omsg.str()); }}
#endif

#define LOG4ESPP_PUSH(string) log4cpp::NDC::push(string)
#define LOG4ESPP_POP() log4cpp::NDC::pop()

#elif defined(HAVE_LOG4CXX) and defined(LOG4ESPP_USE_LOG4CXX)

#include <log4cxx/logstring.h>
#include <stdlib.h>
#include <log4cxx/logger.h>
#include <log4cxx/defaultconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/helpers/exception.h>
#include <log4cxx/logmanager.h>
#include <log4cxx/ndc.h>
#include <locale.h>

using namespace log4cxx;
using namespace log4cxx::helpers;

#define LOG4ESPP_CONFIGURE() DefaultConfigurator::configure(LogManager::getLoggerRepository())
#define LOG4ESPP_ROOTLOGGER(aLogger) log4cxx::LoggerPtr aLogger = Logger::getRootLogger()
#define LOG4ESPP_LOGGER(aLogger,name) log4cxx::LoggerPtr aLogger = log4cxx::Logger::getLogger(name)
#define LOG4ESPP_DEF_LOGGER(aLogger) log4cxx::LoggerPtr aLogger

#if defined(LOG4ESSP_LEVEL_FATAL)
#define LOG4ESPP_DEBUG(logger,msg) 
#define LOG4ESPP_INFO(logger,msg)  
#define LOG4ESPP_WARN(logger,msg)  
#define LOG4ESPP_ERROR(logger,msg) 
#define LOG4ESPP_FATAL(logger,msg) LOG4CXX_FATAL(logger,msg)
#elif defined(LOG4ESPP_LEVEL_ERROR)
#define LOG4ESPP_DEBUG(logger,msg) 
#define LOG4ESPP_INFO(logger,msg)  
#define LOG4ESPP_WARN(logger,msg)  
#define LOG4ESPP_ERROR(logger,msg) LOG4CXX_ERROR(logger,msg)
#define LOG4ESPP_FATAL(logger,msg) LOG4CXX_FATAL(logger,msg)
#elif defined(LOG4ESPP_LEVEL_WARN)
#define LOG4ESPP_DEBUG(logger,msg)
#define LOG4ESPP_INFO(logger,msg) 
#define LOG4ESPP_WARN(logger,msg)  LOG4CXX_WARN(logger,msg)
#define LOG4ESPP_ERROR(logger,msg) LOG4CXX_ERROR(logger,msg)
#define LOG4ESPP_FATAL(logger,msg) LOG4CXX_FATAL(logger,msg)
#elif defined(LOG4ESPP_LEVEL_INFO)
#define LOG4ESPP_DEBUG(logger,msg)
#define LOG4ESPP_INFO(logger,msg)  LOG4CXX_INFO(logger,msg)
#define LOG4ESPP_WARN(logger,msg)  LOG4CXX_WARN(logger,msg)
#define LOG4ESPP_ERROR(logger,msg) LOG4CXX_ERROR(logger,msg)
#define LOG4ESPP_FATAL(logger,msg) LOG4CXX_FATAL(logger,msg)
#else
#define LOG4ESPP_DEBUG(logger,msg) LOG4CXX_DEBUG(logger,msg)
#define LOG4ESPP_INFO(logger,msg)  LOG4CXX_INFO(logger,msg)
#define LOG4ESPP_WARN(logger,msg)  LOG4CXX_WARN(logger,msg)
#define LOG4ESPP_ERROR(logger,msg) LOG4CXX_ERROR(logger,msg)
#define LOG4ESPP_FATAL(logger,msg) LOG4CXX_FATAL(logger,msg)
#endif

#define LOG4ESPP_PUSH(string) NDC::push(string)
#define LOG4ESPP_POP()        NDC::pop()

#else 

#include <iostream>

#define LOG4ESPP_CONFIGURE()
#define LOG4ESPP_ROOTLOGGER(aLogger) 
#define LOG4ESPP_LOGGER(aLogger,name) 

#define LOG4ESPP_DEBUG(logger,msg) std::cout << "DEBUG: " << msg << "\n"
#define LOG4ESPP_INFO(logger,msg)  std::cout << "INFO: " << msg << "\n"
#define LOG4ESPP_WARN(logger,msg)  std::cout << "WARN: " << msg << "\n"
#define LOG4ESPP_ERROR(logger,msg) std::cout << "ERROR: " << msg << "\n"
#define LOG4ESPP_FATAL(logger,msg) std::cout << "FATAL: " << msg << "\n"

#endif
