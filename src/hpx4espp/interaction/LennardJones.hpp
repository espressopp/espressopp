/*
  Copyright (C) 2019-2022
      Max Planck Institute for Polymer Research & JGU Mainz
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

// ESPP_CLASS
#ifndef _HPX4ESPP_INTERACTION_LENNARDJONES_HPP
#define _HPX4ESPP_INTERACTION_LENNARDJONES_HPP

#include "vec/interaction/LennardJones.hpp"

namespace espressopp
{
namespace hpx4espp
{
namespace interaction
{
class LennardJones : public vec::interaction::LennardJones
{
public:
    typedef vec::interaction::LennardJones base;

    template <class... Ts>
    LennardJones(Ts... args) : base(args...)
    {
    }

    static void registerPython();
};

using vec::interaction::LennardJones_pickle;

}  // namespace interaction
}  // namespace hpx4espp
}  // namespace espressopp

#endif
