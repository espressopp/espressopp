#include "python.hpp"
#include <cmath>
#include "PressureTensor.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "VerletList.hpp"
#include "storage/NodeGrid.hpp"

#include "Tensor.hpp"

using namespace espresso;
using namespace iterator;
using namespace interaction;

namespace espresso {
  namespace analysis {
  Tensor PressureTensor::computeTensor() const {

    System& system = getSystemRef();

    // determine volume of the box
    Real3D Li = system.bc->getBoxL();
    real V = Li[0] * Li[1] * Li[2];

    // compute the kinetic contribution (2/3 \sum 1/2mv^2)
    Tensor vvLocal(0.0);
    Tensor vv;

    CellList realCells = system.storage->getRealCells();
    for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
      real mass = cit->mass();
      Real3D& vel = cit->velocity();
      vvLocal += mass * Tensor(vel, vel);
    }
    
    boost::mpi::all_reduce(*mpiWorld, vvLocal, vv, std::plus<Tensor>());

    // compute the short-range nonbonded contribution
    Tensor wijLocal(0.0);
    Tensor wij;

    const InteractionList& srIL = system.shortRangeInteractions;
    for (size_t j = 0; j < srIL.size(); j++) {
      srIL[j]->computeVirialTensor(wijLocal);
    }

    boost::mpi::all_reduce(*mpiWorld, wijLocal, wij, std::plus<Tensor>());
    return (vv + wij) / V;
  }

  Tensor PressureTensor::computeTensor1(real xmin, real xmax,
          real ymin, real ymax, real zmin, real zmax) const {

    System& system = getSystemRef();
    const bc::BC& bc = *system.bc;  // boundary conditions

    // determine the local volume size
    real V = (xmax-xmin) * (ymax-ymin) * (zmax-zmin);

    // compute the kinetic contribution (2/3 \sum 1/2mv^2)
    Tensor vvlocal(0.0);
    CellList realCells = system.storage->getRealCells();
    for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
      Real3D pos = cit->position();
      if( pos[0]>xmin && pos[0]<xmax && 
              pos[1]>ymin && pos[1]<ymax && 
              pos[2]>zmin && pos[2]<zmax ){
        real mass = cit->mass();
        Real3D& vel = cit->velocity();
        vvlocal += mass * Tensor(vel, vel);
      }
    }
    Tensor vv(0.0);
    boost::mpi::all_reduce(*mpiWorld, vvlocal, vv, std::plus<Tensor>());

    // compute the short-range nonbonded contribution
    Tensor w(0.0);
    const InteractionList& srIL = system.shortRangeInteractions;
    for (size_t j = 0; j < srIL.size(); j++) {
      srIL[j]->computeVirialTensor(w, xmin, xmax, ymin, ymax, zmin, zmax);
    }

    return ( vv + w ) / V;
  }

    real PressureTensor::compute() const {
      return -1.0;
    }

    // TODO it is fast solution. one should think about the overloading
    
    using namespace boost::python;

    void PressureTensor::registerPython() {
      using namespace espresso::python;
      class_<PressureTensor, bases< Observable > >
        ("analysis_PressureTensor", init< shared_ptr< System > >())
        .def("compute1", &PressureTensor::computeTensor)
        .def("compute2", &PressureTensor::computeTensor1)
      ;
    }
  }
}
