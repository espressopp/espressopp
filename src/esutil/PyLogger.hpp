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

#ifndef PyLogger_HPP
#define PyLogger_HPP

/** \file PyLogger.hpp    Python implementation for logging.
*/

#include "python.hpp"
#include <esutil/Logger.hpp>

namespace log4espp {

  /**************************************************************************
  *                                                                         *
  *                                                                         *
  *                                                                         *
  **************************************************************************/

  /** PyLogger is an implementaton of the abstract class Logger that uses 
      corresponding Python Logger objects for the output of the logging 
      statements.

      One major efficiency aspect is that an object of this class will not
      ask for the logging level of its Python equivalent; therefore it has 
      to be ensured that changing the logging level of a Python logger
      will notity an object of this C++ class.

  */

  class PyLogger : public Logger {

  private:

   boost::python::object pyLogger;   //!< Pointer to the Python instance of this logger.

   /** This routine updates the logging level by the given level of the Python logger

       \param pyLevel is the Python coding of the level

       The Python coding of the level will be translated to Logger::Level.
   */

   void setPythonLevel(int pyLevel);

   /** Commonly used routine for logging output that can be used for all levels

       \param level is the string representation of the level
       \param loc is the file location of the logging statement
       \param msg is the message output of the logging statement

       The logging output will be printed via the Python Logger object. 
   */

   void log(int level, Location& loc, const std::string& msg);

  public:

   typedef boost::python::object object;

   /** Initialization routine to be called with initialization */

   static void initLogging();

   static void registerPython();

   /** Constructor for a Logger that will use a corresponding Python instance
       for the configuration and for the output of the logging messages  
   */

   PyLogger(std::string, class Logger* parent);

   ~PyLogger() {}

   /** This routine sets the Python instance of the logger and/or 
       updates the logging level of the C++ class.

       \param pyLogger is the Python Logger object.

       This routine must also be called to update the logging level
       of this object.
   */

   void setPythonLogger(boost::python::object pyLogger);

   /** This routine sets the Python loggers for this logger and the
       descendants.

       \param ParentName is the name of the parent of this logger.

   */

   void setPythonLoggers(std::string& parentName);

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

  };

}

#endif
