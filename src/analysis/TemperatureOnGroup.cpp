/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)

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
#include "TemperatureOnGroup.hpp"

namespace espressopp
{
namespace analysis
{
void TemperatureOnGroup::registerPython()
{
    using namespace espressopp::python;
    class_<TemperatureOnGroup, bases<Observable> >(
        "analysis_TemperatureOnGroup",
        init<std::shared_ptr<System>, std::shared_ptr<ParticleGroup> >());
}

}  // end namespace analysis
}  // end namespace espressopp
