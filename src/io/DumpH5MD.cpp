/*
  Copyright (c) 2015
      Jakub Krajniak (jkrajniak at gmail.com)
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

#include <fstream>
#include <iomanip>
#include <string>
#include "DumpH5MD.hpp"
#include "storage/Storage.hpp"

#include "bc/BC.hpp"

#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"

namespace espresso {
namespace io {

void DumpH5MD::dump() {
  shared_ptr<System> system = getSystem();
  ConfigurationsExt conf(system);
  conf.setUnfolded(unfolded);
  conf.gather();
}

void DumpH5MD::registerPython() {
  using namespace espresso::python;  //NOLINT

  class_<DumpH5MD, bases<ParticleAccess>, boost::noncopyable >
  ("io_DumpH5MD", init< shared_ptr< System >,
                       shared_ptr< integrator::MDIntegrator >,
                       std::string,
                       std::string,
                       bool,
                       std::string,
                       std::string>())
    .add_property("unfolded", &DumpH5MD::getUnfolded,
                              &DumpH5MD::setUnfolded)
    .def("dump", &DumpH5MD::dump);
}
}  // end namespace io
}  // end namespace espresso
