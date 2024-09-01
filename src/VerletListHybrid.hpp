/*
  Copyright (C) 2016
      Jakub Krajniak <jkrajniak at gmail.com>

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
#ifndef _VERLETLISTHYBRID_HPP
#define _VERLETLISTHYBRID_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"
#include "VerletList.hpp"

namespace espressopp
{
class VerletListHybridAT : public VerletList
{
public:
    VerletListHybridAT(std::shared_ptr<System> system, real cut, bool rebuildVL)
        : VerletList(system, cut, rebuildVL)
    {
    }

    /** Register this class so it can be used from Python. */
    static void registerPython();

protected:
    void checkPair(Particle &pt1, Particle &pt2);

    static LOG4ESPP_DECL_LOGGER(theLogger);
};

class VerletListHybridCG : public VerletList
{
public:
    VerletListHybridCG(std::shared_ptr<System> system, real cut, bool rebuildVL)
        : VerletList(system, cut, rebuildVL)
    {
    }

    /** Register this class so it can be used from Python. */
    static void registerPython();

protected:
    void checkPair(Particle &pt1, Particle &pt2);

    static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // namespace espressopp

#endif
