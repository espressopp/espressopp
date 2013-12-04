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
