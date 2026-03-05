/*
  Copyright (C) 2022
      Data Center, Johannes Gutenberg University Mainz

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
#include <boost/signals2.hpp>
#include "CoulombScafacos.hpp"
#include "CellListAllParticlesInteractionTemplate.hpp"

int ifTuned;
fcs_int p3m_grid, p3m_cao;
fcs_float p3m_rcut, p3m_alpha;

namespace espressopp
{
namespace interaction
{
typedef class CellListAllParticlesInteractionTemplate<CoulombScafacos> CellListCoulombScafacos;

CoulombScafacos::CoulombScafacos(
    shared_ptr<System> _system, real _prefactor, real _tolerance, int _ntotal, const char* _method)
{
    system = _system;
    prefactor = _prefactor;
    tolerance = _tolerance;
    ntotal = _ntotal;
    method = _method;
    // strcpy(method,_method);

    I = Tensor(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    sum = NULL;
    totsum = NULL;

    preset();
    // getParticleNumber(); // geting the number of particles for the current node // it's done in
    // preset

    // This function calculates the square of all particle charges. It should be called ones,
    // if the total number of particles doesn't change.
    count_charges(system->storage->getRealCells());
}

CoulombScafacos::~CoulombScafacos()
{
    delete[] sum;
    delete[] totsum;
    sum = NULL;
    totsum = NULL;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void CoulombScafacos::registerPython()
{
    using namespace espressopp::python;

    class_<CoulombScafacos, bases<Potential> >(
        "interaction_CoulombScafacos", init<shared_ptr<System>, real, real, int, const char*>())
        .add_property("prefactor", &CoulombScafacos::getPrefactor, &CoulombScafacos::setPrefactor)
        .add_property("tolerance", &CoulombScafacos::getTolerance, &CoulombScafacos::setTolerance)
        .add_property("ntotal", &CoulombScafacos::getNTotal, &CoulombScafacos::setNTotal)
        .add_property("method", &CoulombScafacos::getMethod, &CoulombScafacos::setMethod);

    class_<CellListCoulombScafacos, bases<Interaction> >(
        "interaction_CellListCoulombScafacos",
        init<shared_ptr<storage::Storage>, shared_ptr<CoulombScafacos> >())
        .def("getPotential", &CellListCoulombScafacos::getPotential);
}

}  // namespace interaction
}  // namespace espressopp
