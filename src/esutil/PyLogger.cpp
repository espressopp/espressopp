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

/** \file PyLogger.cpp    Implementation of logging via Python

<b>Responsible:</b>
<a href="mailto:brandes@scai.fraunhofer.de">Thomas Brandes</a>

*/


#include "PyLogger.hpp"
#include <iostream>

/********************************************************************
*  Implementation of trace/debug/info/warn/error/fatal              *
********************************************************************/

#undef DEBUGGING

/********************************************************************
*  Macro to ask a Python object for a certain attribute             *
*   ->  avoids throw exception used by Boost Python getattr         *
********************************************************************/

#define hasattr(pyObj, name) PyObject_HasAttrString(pyObj.ptr(), name)

using namespace log4espp;

using namespace boost::python;
using namespace std;

/********************************************************************
*   Global variable for root of all loggers                         *
********************************************************************/

static PyLogger* rootLogger = NULL;

/********************************************************************
*   Global variables for certain Python objects                     *
********************************************************************/

static object pythonLoggerClass = object();

static bool pyInit = false;

static object pyNone = object();

static object pyLogging = object();

/* Python values for logging levels:
   actual values will be set during initialization
*/

static int pyNOTSET = 0;
static int pyDEBUG  = 0;
static int pyTRACE  = 0;
static int pyINFO   = 0;
static int pyWARN   = 0;
static int pyERROR  = 0;
static int pyFATAL  = 0;

/********************************************************************
*   PyLogger :: constructor                                         *
********************************************************************/

PyLogger::PyLogger(string _name, Logger* _parent) : Logger(_name, _parent)

{
#ifdef DEBUGGING
   printf("PyLogger %s created\n", getFullName().c_str());
#endif

   // Python logger will only be set if Python has already been initialized

   if (pyInit) {
      setPythonLogger(pyLogging.attr("getLogger")(getFullName()));
   } else {
      pyLogger = object();
   }
}

/********************************************************************
*  Logger::configure()                                              *
********************************************************************/

void Logger::configure()  {

   // nothing to do as configuration is done in Python
}

/********************************************************************
*   Logger :: getRoot                                               *
********************************************************************/

Logger& Logger::getRoot() {

  if (rootLogger == NULL) {

     rootLogger = new PyLogger("", NULL);
  }

  return *rootLogger;
}

/********************************************************************
*  Logger :: getInstance                                            *
********************************************************************/

Logger& Logger::getInstance (string name) {
  return Logger::createInstance<PyLogger>(name);
}

/********************************************************************
*  Helper routine for logging via Python                            *
********************************************************************/

void PyLogger::log(int level, Location& loc, const string& msg) 

{ // return if no python logger is availabe 

  if (pyLogger == object()) return;

  object name   = pyLogger.attr("name");

  object record = pyLogger.attr("makeRecord")(name, level, 
                                              loc.filename, loc.line,
                                              msg, pyNone, pyNone);

  record.attr("funcName") = loc.funcname;

  pyLogger.attr("handle")(record);

}

/********************************************************************
*  Implementation of trace/debug/info/warn/error/fatal              *
********************************************************************/
 
/* Example of an inefficient realization for asking the 
   logging level of the Python logger.
  
   bool PyLogger::isTraceEnabled() 

   { if (pyLogger == NULL) return false;
   
     // Solution 1: take the level atrribute (must not be NOTSET)

     int lev = extract<int>(pyLogger.attr("level"));  

     // Solution 2: ask for the effective level in hierarchy

     int lev = extract<int>(pyLogger.attr("getEffectiveLevel")());

     return lev <= pyTRACE;

   }

*/

void PyLogger::trace(Location loc, const string& msg) {

  log(pyTRACE, loc, msg);
}

void PyLogger::debug(Location loc, const string& msg) {

  log(pyDEBUG, loc, msg);
}

void PyLogger::info(Location loc, const string& msg) {

  log(pyINFO, loc, msg);
}

void PyLogger::warn(Location loc, const string& msg) {

  log(pyWARN, loc, msg);
}

void PyLogger::error(Location loc, const string& msg) {

  log(pyERROR, loc, msg);
}

void PyLogger::fatal(Location loc, const string& msg) {
 
  log(pyFATAL, loc, msg);
}

void PyLogger::setPythonLevel(int pyLevel) 

{
  if (pyLevel == pyNOTSET) {
      setFlag = false;
  } else if (pyLevel == pyTRACE) {
      setLevel(TRACE, true);
  } else if (pyLevel == pyDEBUG) {
      setLevel(DEBUG, true);
  } else if (pyLevel == pyINFO) {
      setLevel(INFO, true);
  } else if (pyLevel == pyWARN) {
      setLevel(WARN, true);
  } else if (pyLevel == pyERROR) {
      setLevel(ERROR, true);
  } else if (pyLevel == pyFATAL) {
      setLevel(FATAL, true);
  } else {
      printf("ERROR: setPythonLevel for %s: %d is unknown log level of Python\n", 
              getFullName().c_str(), pyLevel);
  }

#ifdef DEBUG
  printf("setPythonLevel: level of logger %s is now %d, setflag = %d\n",
          getFullName().c_str(), myLevel, setFlag);
#endif

}

/********************************************************************
*  PyLogger :: setPythonLogger (python_logger_object)               *
********************************************************************/

void PyLogger::setPythonLogger(object _pyLogger)

{
  if (pyLogger == object()) {
      pyLogger = _pyLogger;
  } else if (pyLogger != _pyLogger) {
      printf("ATTENTION: Python Logger object for %s has changed\n", getFullName().c_str());
  } 

  // get the level of this logger and set it correctly

  int level = extract<int>(pyLogger.attr("level"));

  setPythonLevel(level);
}

/********************************************************************
*  PyLogger :: setPythonLoggers    for this logger and all sons     *
********************************************************************/

void PyLogger::setPythonLoggers(string& parentName) 

{  
   // set also the Python loggers for the sons

   string fullName;

   if (parentName == "")
      fullName  = name;
   else
      fullName  = parentName + "." + name;
 
#ifdef DEBUGGING
   printf("setPytthonLoggers for %s\n", fullName.c_str());  
#endif 

   object pyLogger = pyLogging.attr("getLogger") (fullName);

   // set the logger; will inherit the level 

   setPythonLogger(pyLogger);

   for (size_t i = 0; i < sons.size(); i++) {
       PyLogger* son = (PyLogger*)(sons[i]);
       son->setPythonLoggers(fullName);
   }
}

/********************************************************************
*                                                                   *
********************************************************************/

static void pyInitialize() {

  if (!pyInit) {

#ifdef DEBUGGING
     printf("PyLogger:: pyInitialize\n");  
#endif 
     // ToDo: make sure that Python itself has been initialized
     //       even it will probably never be called this way

     // do some more initiliaztion

     pyLogging = import("logging");

#ifdef DEBUGGING
     printf("PyLogger:: import logging done\n");  
#endif 

     // ToDo: make sure that Python itself has been initialized
     pyNOTSET =  extract<int>(pyLogging.attr("NOTSET"));
     pyDEBUG  =  extract<int>(pyLogging.attr("DEBUG"));
     pyTRACE  =  extract<int>(pyLogging.attr("TRACE"));
     pyINFO   =  extract<int>(pyLogging.attr("INFO"));
     pyWARN   =  extract<int>(pyLogging.attr("WARN"));
     pyERROR  =  extract<int>(pyLogging.attr("ERROR"));
     pyFATAL  =  extract<int>(pyLogging.attr("FATAL"));

#ifdef DEBUGGING
     printf("PyLogger:: got the logging levels\n");  
#endif 

     pyInit = true;
  }
}

/********************************************************************
*                                                                   *
********************************************************************/

/** This routine can be used to update the logging level via Python */

void loggerUpdate(object pythonLogger) {

  // get a C++ Logger for this Python Logger

  string name = extract<string>(pythonLogger.attr("name"));

  if (name == "root") name = "";

  Logger& logger = Logger::getInstance(name);

  // Dynamic cast: is safe as all logger objects are PyLogger

  PyLogger *pyLogger = dynamic_cast<PyLogger*>(&logger);

  pyLogger->setPythonLogger(pythonLogger);
}

/** Initialization routine called from Python

    \param loggerClass is the name of the class to which all loggers of
           will belong; this class guarantees that log level changes will
           be availabe here.
*/

void PyLogger::initLogging()

{
#ifdef DEBUGGING
  printf("initLogging for PyLogger\n");  
#endif 

  pyInitialize();

  // make a root logger if not available yet

  if (rootLogger == NULL) {
     rootLogger = new PyLogger("", NULL);
  }

  // now traverse all existing C++ loggers and update them

  string name = "";
  rootLogger->setPythonLoggers(name);
}

/********************************************************************
*   exportLogging()                                                 *
********************************************************************/

/** This routine exports the bindings for python */

void PyLogger::registerPython() {

  def("setLogger", &PyLogger::initLogging);
  def("setLogger", &loggerUpdate);
}
