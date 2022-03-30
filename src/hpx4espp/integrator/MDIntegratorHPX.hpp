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

#ifndef HPX4ESPP_INTEGRATOR_MDINTEGRATORHPX_HPP
#define HPX4ESPP_INTEGRATOR_MDINTEGRATORHPX_HPP

#include "hpx4espp/include/types.hpp"
#include "hpx4espp/SystemHPX.hpp"

#include "esconfig.hpp"
#include "types.hpp"
#include "integrator/MDIntegrator.hpp"
#include "Extension.hpp"
#include "esutil/ESPPIterator.hpp"
#include "log4espp.hpp"
#include <boost/signals2.hpp>
#include <iostream>

namespace espressopp
{
namespace hpx4espp
{
namespace integrator
{
struct ExtensionList : public std::vector<shared_ptr<hpx4espp::integrator::Extension> >
{
    typedef espressopp::esutil::ESPPIterator<std::vector<hpx4espp::integrator::Extension> >
        Iterator;
};

class MDIntegratorHPX : public espressopp::integrator::MDIntegrator
{
public:
    typedef espressopp::integrator::MDIntegrator baseClass;

    MDIntegratorHPX(shared_ptr<SystemHPX> system, shared_ptr<storage::StorageHPX> storageHPX)
        : storageHPX(storageHPX), systemHPX(system), baseClass(system)
    {
    }

    void addExtension(shared_ptr<hpx4espp::integrator::Extension> extension);

    shared_ptr<hpx4espp::integrator::Extension> getExtension(int k);

    int getNumberOfExtensions();

    inline shared_ptr<storage::StorageHPX> getStorageHPX() { return storageHPX; }

    // Boost signals for extensions
    // TODO: These extensions should take futures as arguments as they will be used in
    // continuation-style programming of integration steps. e.g.

    // signals to extend the integrator
    boost::signals2::signal<void()> runInit;  // initialization of run()
    boost::signals2::signal<void()> recalc1;  // inside recalc, before updateForces()
    // boost::signals2::signal<void ()> aftCalcSlow; // after calculation of slow forces
    // updateForces(true) in VerlocityVerletRESPA
    boost::signals2::signal<void()> recalc2;  // inside recalc, after  updateForces()
    // boost::signals2::signal<void ()> befIntP; // before integrate1()
    // boost::signals2::signal<void (real&)> inIntP; // inside end of integrate1()
    // boost::signals2::signal<void ()> aftIntP; // after  integrate1()
    // boost::signals2::signal<void ()> aftInitF; // after initForces()
    // boost::signals2::signal<void ()> aftCalcFLocal; // after calcForces in local cells (before
    // collectGhostForces)
    boost::signals2::signal<void()> aftCalcF;  // after calcForces()
    // boost::signals2::signal<void ()> befIntV; // before integrate2()
    // boost::signals2::signal<void ()> aftIntV; // after  integrate2()
    // boost::signals2::signal<void ()> aftIntSlow; // after integrateSlow() in VerlocityVerletRESPA

    static void registerPython();

protected:
    ExtensionList exList;

    shared_ptr<SystemHPX> systemHPX;

    shared_ptr<storage::StorageHPX> storageHPX;

private:
    static LOG4ESPP_DECL_LOGGER(logger);
};

}  // namespace integrator
}  // namespace hpx4espp
}  // namespace espressopp

#endif  // HPX4ESPP_INTEGRATOR_MDINTEGRATORHPX_HPP
