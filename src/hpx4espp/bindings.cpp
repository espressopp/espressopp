/*
  Copyright (C) 2020-2022
      Max Planck Institute for Polymer Research & JGU Mainz

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

#if HPX4ESPP_ENABLED

#include "hpx4espp/include/hpx_version.hpp"

#if (HPX4ESPP_HPX_VERSION_FULL >= 10701)
#include <hpx/modules/runtime_configuration.hpp>
#else
#include <hpx/runtime/components/component_factory_base.hpp>
#endif

HPX_REGISTER_COMPONENT_MODULE()

#include "hpx4espp/HPXRuntime.hpp"
#include "hpx4espp/SystemHPX.hpp"
#include "hpx4espp/VerletList.hpp"

#include "hpx4espp/esutil/bindings.hpp"
#include "hpx4espp/interaction/bindings.hpp"
#include "hpx4espp/integrator/bindings.hpp"
#include "hpx4espp/storage/bindings.hpp"

#endif  // HPX4ESPP_ENABLED

#include "python.hpp"
#include "bindings.hpp"
#include "utils/HPXUtils.hpp"

namespace espressopp
{
namespace hpx4espp
{
void registerPython()
{
    hpx4espp::utils::HPXUtils::registerPython();

#if HPX4ESPP_ENABLED
    hpx4espp::HPXRuntime::registerPython();
    hpx4espp::SystemHPX::registerPython();
    hpx4espp::VerletList::registerPython();

    hpx4espp::integrator::registerPython();
    hpx4espp::interaction::registerPython();
    hpx4espp::storage::registerPython();
    hpx4espp::esutil::registerPython();
#endif  // HPX4ESPP_ENABLED
}
}  // namespace hpx4espp
}  // namespace espressopp
