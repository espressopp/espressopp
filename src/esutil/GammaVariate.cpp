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

#include "python.hpp"
#include "GammaVariate.hpp"
#include "esutil/RNG.hpp"

using namespace boost;

namespace espressopp
{
namespace esutil
{
GammaVariate::GammaVariate(std::shared_ptr<RNG> _rng, const int alpha, const real beta)
    : Super(*(_rng->getBoostRNG()), DistType(alpha, beta)), rng(_rng)
{
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void GammaVariate::registerPython()
{
    using namespace espressopp::python;

    real (GammaVariate::*pyCall)() = &GammaVariate::operator();

    class_<GammaVariate>("esutil_GammaVariate", init<std::shared_ptr<RNG> >())
        .def("__call__", pyCall);
}
}  // namespace esutil
}  // namespace espressopp
