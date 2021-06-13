/*
  Copyright (C) 2017,2021
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

#include "python.hpp"
#include "types.hpp"
#include "ExtAnalyze.hpp"
#include "SystemAccess.hpp"
//#include "analysis/AnalysisBase.hpp"
#include "ParticleAccess.hpp"

namespace espressopp
{
// using namespace analysis;
namespace integrator
{
LOG4ESPP_LOGGER(ExtAnalyze::theLogger, "ExtAnalyze");

// ExtAnalyze::ExtAnalyze(std::shared_ptr< AnalysisBase > _analysis, int _interval) :
// Extension(_analysis->getSystem()), interval(_interval)
ExtAnalyze::ExtAnalyze(std::shared_ptr<ParticleAccess> _particle_access, int _interval)
    : Extension(_particle_access->getSystem()), interval(_interval)
{
    LOG4ESPP_INFO(theLogger, "Analyze observable in integrator");
    // analysis     = _analysis;
    particle_access = _particle_access;
    type = Extension::ExtAnalysis;
}

void ExtAnalyze::disconnect() { _aftIntV.disconnect(); }

void ExtAnalyze::connect()
{
    // connection to end of integrator
    _aftIntV = integrator->aftIntV.connect(std::bind(&ExtAnalyze::perform_action, this));
}

// void ExtAnalyze::performMeasurement() {
void ExtAnalyze::perform_action()
{
    LOG4ESPP_INFO(theLogger, "performing measurement in integrator");
    if ((integrator->getStep() - 1) % interval == 0)
    {
        particle_access->perform_action();
    }
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void ExtAnalyze::registerPython()
{
    using namespace espressopp::python;
    class_<ExtAnalyze, std::shared_ptr<ExtAnalyze>, bases<Extension> >(
        "integrator_ExtAnalyze", init<std::shared_ptr<ParticleAccess>, int>())
        .def("connect", &ExtAnalyze::connect)
        .def("disconnect", &ExtAnalyze::disconnect);
}
}  // namespace integrator
}  // namespace espressopp
