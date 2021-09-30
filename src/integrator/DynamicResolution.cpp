/*
  Copyright (C) 2015-2016,2021
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

#include "DynamicResolution.hpp"
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include "python.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp
{
namespace integrator
{
using namespace espressopp::iterator;  // NOLINT

LOG4ESPP_LOGGER(DynamicResolution::theLogger, "DynamicResolution");

DynamicResolution::DynamicResolution(shared_ptr<System> _system,
                                     shared_ptr<FixedVSList> _vslist,
                                     real _rate)
    : Extension(_system), rate_(_rate), vs_list(_vslist)
{
    LOG4ESPP_INFO(theLogger, "construct DynamicResolution");
    type = Extension::Adress;
}

DynamicResolution::~DynamicResolution()
{
    LOG4ESPP_INFO(theLogger, "~DynamicResolution");
    disconnect();
}

void DynamicResolution::connect()
{
    LOG4ESPP_INFO(theLogger, "connect");
    _aftIntV = integrator->aftIntV.connect(boost::bind(&DynamicResolution::updateWeights, this),
                                           boost::signals2::at_back);
}

void DynamicResolution::disconnect()
{
    LOG4ESPP_INFO(theLogger, "disconnect");
    _aftIntV.disconnect();
}

void DynamicResolution::set_active(bool active) { active_ = active; }

void DynamicResolution::updateWeights()
{
    LOG4ESPP_INFO(theLogger, "updateWeights");
    if (!active_) return;

    System &system = getSystemRef();

    // Update weights of the particles.
    FixedVSList::GlobalTuples vs = vs_list->globalTuples;
    FixedVSList::GlobalTuples::iterator it = vs.begin();
    for (; it != vs.end(); ++it)
    {
        Particle *vp = system.storage->lookupLocalParticle(it->first);
        if (vp)
        {
            real res = vp->lambda() + rate_;
            if (res > (1.0 + rate_)) continue;
            if (res > 1.0) res = 1.0;
            vp->lambda() = res;

            // Update weights for all underlying particles.
            for (FixedVSList::tuple::iterator itp = it->second.begin(); itp != it->second.end();
                 ++itp)
            {
                Particle *at = system.storage->lookupLocalParticle(*itp);
                if (at)
                {
                    at->lambda() = res;
                }
                else
                {
                    std::cout << " AT particle (" << *itp << ") of VP " << vp->id() << "-"
                              << vp->ghost() << " not found in tuples ";
                    std::cout << " (" << vp->position() << ")" << std::endl;
                    exit(1);
                }
            }
        }
    }
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void DynamicResolution::registerPython()
{
    using namespace espressopp::python;  // NOLINT
    class_<DynamicResolution, shared_ptr<DynamicResolution>, bases<Extension> >(
        "integrator_DynamicResolution", init<shared_ptr<System>, shared_ptr<FixedVSList>, real>())
        .add_property("active", &DynamicResolution::active, &DynamicResolution::set_active)
        .add_property("rate", &DynamicResolution::rate, &DynamicResolution::set_rate)
        .def("update_weights", &DynamicResolution::updateWeights)
        .def("connect", &DynamicResolution::connect)
        .def("disconnect", &DynamicResolution::disconnect);
}
}  // end namespace integrator
}  // end namespace espressopp
