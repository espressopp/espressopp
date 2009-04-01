#ifndef GenLogger_HPP
#define GenLogger_HPP

/** \file GenLogger.hpp    Python implementation for logging.

<b>Responsible:</b>
<a href="mailto:brandes@scai.fraunhofer.de">Thomas Brandes</a>

*/

#include <esutil/Logger.hpp>

namespace log4espp {

  /**************************************************************************
  *                                                                         *
  *                                                                         *
  *                                                                         *
  **************************************************************************/

  /** GenLogger is a very generic C++ implementation of the abstract Logger class */

  class GenLogger : public Logger {

  private:

    /** Generic routine for logging output that can be used for all levels 

        \param level is the string representation of the level
        \param loc is the file location of the logging statement
        \param msg is the message output of the logging statement
    */

    void log(const char* level, Location& loc, const std::string& msg);

  public:

   /** Constructor for o a generic logger.

       \param name is the name of the logger at this level
       \param parent is a pointer to the ancestor logger (NULL for root)

       \sa Logger::Logger
   */

   GenLogger(std::string name, class Logger* parent);

   ~GenLogger() {}

   /** Implementation of Logger::trace */

   virtual void trace(Location loc, const std::string& msg);

   /** Implementation of Logger::debug */

   virtual void debug(Location loc, const std::string& msg);

   /** Implementation of Logger::info */

   virtual void info (Location loc, const std::string& msg);

   /** Implementation of Logger::warn */

   virtual void warn (Location loc, const std::string& msg);

   /** Implementation of Logger::error */

   virtual void error(Location loc, const std::string& msg);

   /** Implementation of Logger::fatal */

   virtual void fatal(Location loc, const std::string& msg);

   /** Traversing this logger and all output loggers, can be used for DEBUG of Logger */

   void traverse();

  };

}

#endif
