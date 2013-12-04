/*
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

#ifndef _ESUTIL_LOGGER_HPP
#define _ESUTIL_LOGGER_HPP
/** \file Logger.hpp    Class for logging.
*/

#include <string>
#include <vector>

#if !defined(LOG4ESPP_LOCATION)
#if defined(_MSC_VER)
#if _MSC_VER >= 1300
      #define __LOG4ESPP_FUNC__ __FUNCSIG__
#endif
#else
#if defined(__GNUC__) and defined(LOG4ESPP_LONGNAMES)
      #define __LOG4ESPP_FUNC__ __PRETTY_FUNCTION__
#else
      #define __LOG4ESPP_FUNC__ __FUNCTION__
#endif
#endif

#if !defined(__LOG4ESPP_FUNC__)
#define __LOG4ESPP_FUNC__ ""
#endif

#define LOG4ESPP_LOCATION log4espp::Location(__FILE__, __LOG4ESPP_FUNC__, __LINE__)
#endif

/** Namespace for logging in Espresso */

namespace log4espp {

  /**************************************************************************
  *  class Location                                                         *
  **************************************************************************/

  /** Location is a class containing file, line and function info; it
      specifies a location in a source file, e.g. for logging statement  */

  class Location {

   public:

     const char* filename;   //!< Name of the source file
     const char* funcname;   //!< Name of the function 
     int line;               //!< Line number of location in source file

     /** Constructor of a location */

     Location(const char* _filename, const char* _funcname, int _line) {
       filename = _filename;
       funcname = _funcname;
       line     = _line;
     }

  };

  /**************************************************************************
  *  class Logger                                                           *
  **************************************************************************/

  /** Logger is a an abstract class for hierarchical organization of logging objects 

  */

  class Logger {

  public:

     typedef enum {TRACE  = 10, //!< even more detailed than DEBUG
                   DEBUG  = 20, //!< designates fine-grained informational events
                   INFO   = 30, //!< informational messages highlighting progress
                   WARN   = 40, //!< for potentially harmful situations
                   ERROR  = 50, //!< for errors that might still allow the application to continue
                   FATAL  = 60  //!< servere errors that will presumably lead to aborts
     } Level;

  protected:

    std::string name;  //!< name of the logger
    bool  setFlag;     //!< This flag indicates that level has been set explicitly.

    Level myLevel;     //!< Current level of this logger.

    class Logger* parent;  //!< points back to the parent logger, NULL for root

    std::vector<class Logger*> sons;   //!< sub loggers of this logger

    static std::vector<std::string> context;  //!< stack for different context

  public:

    /** General constructor for a logger 

        \param _name is the name of this logger (not full name)
        \param _parent is pointer to the parent logger
    */

    Logger(std::string _name, class Logger* _parent) {
        name     = _name;
        parent   = _parent;
        setFlag  = false;
        myLevel  = WARN;
        if (parent != NULL) {
           parent->sons.push_back(this);
           myLevel = parent->getLevel();
        }
    }

    /** Virtual estructor needed due to virtual functions */

    virtual ~Logger() {}

    /** Get reference to the root logger. */

    static Logger& getRoot();

    /** Get reference to a logger with a given name */

    static Logger& getInstance(std::string name);

    /** Add some content for the following log statements */

    static void push(std::string item) { context.push_back(item); }

    /** Reverse to the last push */

    static void pop() { context.pop_back(); }

    /** Asks for the full name of the logger, e.g. "X.Y.Z" */

    std::string getFullName() {
      if (parent == NULL) return name;
      std::string fullname = parent->getFullName();
      if (fullname == "") return name;
      return fullname + "." + name;
    }

    /** Configuration of the logging system. */

    static void configure();

    /** Check if logging statements of level TRACE are enabled. */

    bool isTraceEnabled() { return myLevel <= TRACE; }

    /** Check if logging statements of level DEBUG are enabled. */

    bool isDebugEnabled() { return myLevel <= DEBUG; }

    /** Check if logging statements of level INFO are enabled. */

    bool isInfoEnabled()  { return myLevel <= INFO; }

    /** Check if logging statements of level WARN are enabled. */

    bool isWarnEnabled()  { return myLevel <= WARN; }

    /** Check if logging statements of level ERROR are enabled. */

    bool isErrorEnabled() { return myLevel <= ERROR; }

    /** Check if logging statements of level FATAL are enabled. */

    bool isFatalEnabled() { return myLevel <= FATAL; }

    /** Getter routine for the logging level of this object. */

    Level getLevel() { return myLevel; }

    /** Getter routine for the effective logging level of 
        this object. If level has not been set explicitly it
        will ask the ancestors for the level 
    */

    Level getEffectiveLevel() 
    { if (setFlag || parent == NULL) return myLevel; 
      return parent->getEffectiveLevel();
    }

    /** Getter routine for the logging level of this object. 

        \return the logging level as a string.
    */

    static const char* getLevelString(int level) { 

       if (level == TRACE) return "TRACE";
       if (level == DEBUG) return "DEBUG";
       if (level == INFO)  return "INFO";
       if (level == WARN)  return "WARN";
       if (level == ERROR) return "ERROR";
       if (level == FATAL) return "FATAL";
       return "UNKNOWN";
    }

    /** Setter routine for the logging level of this object. 

        \param level is the new logging level for this object
        \param force true specifies that this is an explicit set

        This routine will set implicitly the levels recursively 
        for the descendants whose level has not been set explicitly.
    */

    void setLevel(Level level, bool force = true)
    { if (!force && setFlag) return;
      myLevel = level;
      setFlag = force;
      for (size_t i = 0; i < sons.size(); i++) {
         sons[i]->setLevel(level, false);
      }
    }

    /** Logging output for level TRACE. This routine should
        only be called if TraceEnabled() returns true.

        \param loc is the file location of the logging statement
        \param is the message to be printed

        Each derived class has to implement this routine. This
        abstract class does not handle output of logging at all.
    */

    virtual void trace(Location loc, const std::string& msg) = 0;


    /** Logging output for level DEBUG. This routine should
        only be called if isDebugEnabled() returns true.

        \param loc is the file location of the logging statement
        \param is the message to be printed

        Each derived class has to implement this routine. This
        abstract class does not handle output of logging at all.
    */

    virtual void debug(Location loc, const std::string& msg) = 0;


    /** Logging output for level INFO. This routine should
        only be called if isInfoEnabled() returns true.

        \param loc is the file location of the logging statement
        \param is the message to be printed

        Each derived class has to implement this routine. This
        abstract class does not handle output of logging at all.
    */

    virtual void info (Location loc, const std::string& msg) = 0;

    /** Logging output for level WARN. This routine should
        only be called if isWarnEnabled() returns true.

        \param loc is the file location of the logging statement
        \param is the message to be printed

        Each derived class has to implement this routine. This
        abstract class does not handle output of logging at all.
    */

    virtual void warn (Location loc, const std::string& msg) = 0;

    /** Logging output for level ERROR. This routine should
        only be called if isErrorEnabled() returns true.

        \param loc is the file location of the logging statement
        \param is the message to be printed

        Each derived class has to implement this routine. This
        abstract class does not handle output of logging at all.
    */

    virtual void error(Location loc, const std::string& msg) = 0;

    /** Logging output for level FATAL. This routine should
        only be called if isFatalEnabled() returns true.

        \param loc is the file location of the logging statement
        \param is the message to be printed

        Each derived class has to implement this routine. This
        abstract class does not handle output of logging at all.
    */

    virtual void fatal(Location loc, const std::string& msg) = 0;

    // static method to create instance of a derived logger

    template<class DerivedLogger>
    static Logger& createInstance(std::string name);

  };

  /** This routine  creates a logger instance for a given name.
      The template argument specifies the class of the logger.
  */

  template<class DerivedLogger>
  Logger& Logger::createInstance(std::string name)

  { // tokenize the name with delimiter "."
  
    std::vector<std::string> tokens;
    
    // Skip delimiters at beginning.
    std::string::size_type lastPos = name.find_first_not_of(".", 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = name.find_first_of(".", lastPos);
  
    while (std::string::npos != pos || std::string::npos != lastPos) {

      // Found a token, add it to the vector.
      tokens.push_back(name.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = name.find_first_not_of(".", pos);
      // Find next "non-delimiter"
      pos = name.find_first_of(".", lastPos);
    }
  
    // now find the logger in the hierarchy tree

    Logger* instance = &getRoot();

    for (size_t i = 0; i < tokens.size(); i++) {

        // find a son for the next token

        Logger* son = NULL;
        for (size_t s = 0; s < instance->sons.size(); s++) {
            Logger* candidate = instance->sons[s];
            if (candidate->name == tokens[i]) {
                son = candidate;
                break;
            }
        }
 
        // create a new son if not found

        if (son == NULL) {
           // logger not available, so create it
           son = new DerivedLogger(tokens[i], instance);
        }

        // got to the next deeper level

        instance = son;
    }

    return *instance;
  }

}

#endif
