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

/*
 This file contains the main routine for the embedded python
 interpreter version of ESPResSo.
*/
#include "acconfig.hpp"
#include "espressopp_common.hpp"
#include "python.hpp"

#include <boost/python/numpy.hpp>

static void finalize();

BOOST_PYTHON_MODULE(_espressopp)
{
    initMPIEnv();

    boost::python::numpy::initialize();

    // register all classes with python
    espressopp::registerPython();

    Py_AtExit(&finalize);
}

void finalize() { finalizeMPIEnv(); }
