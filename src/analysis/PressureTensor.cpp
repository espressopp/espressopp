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

  /*
  python::list PressureTensor::computeTensorLocal() const {

    System& system = getSystemRef();

    // determine local box coordinates
    real xmin, xmax, ymin, ymax, zmin, zmax;
    python::list boxVolume;

    xmin = system.storage->getLocalBoxXMin();
    xmax = system.storage->getLocalBoxXMax();
    ymin = system.storage->getLocalBoxYMin();
    ymax = system.storage->getLocalBoxYMax();
    zmin = system.storage->getLocalBoxZMin();
    zmax = system.storage->getLocalBoxZMax();

    boxVolume.append(xmin);
    boxVolume.append(xmax);
    boxVolume.append(ymin);
    boxVolume.append(ymax);
    boxVolume.append(zmin);
    boxVolume.append(zmax);

    // determine the local volume size
    real V = (xmax-xmin) * (ymax-ymin) * (zmax-zmin);

    // compute the kinetic contribution (2/3 \sum 1/2mv^2)
    Tensor vvLocal(0.0);
    CellList realCells = system.storage->getRealCells();
    for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
      real mass = cit->mass();
      Real3D& vel = cit->velocity();
      vvLocal += mass * Tensor(vel, vel);
    }

    // compute the short-range nonbonded contribution
    Tensor wijLocal(0.0);
    const InteractionList& srIL = system.shortRangeInteractions;
    for (size_t j = 0; j < srIL.size(); j++) {
      srIL[j]->computeVirialTensor(wijLocal);
    }

    python::list res;

    res.append(boxVolume);
    res.append((vvLocal + wijLocal) / V);

    return res;
  }
   */

    // TODO: this dummy routine is still needed as we have not yet TensorObservable
    real PressureTensor::compute() const {
      return -1.0;
    }

    using namespace boost::python;

    void PressureTensor::registerPython() {
      using namespace espresso::python;
      class_<PressureTensor, bases< Observable > >
        ("analysis_PressureTensor", init< shared_ptr< System > >())
        .def("compute", &PressureTensor::computeTensor) 
        //.def("computeLocal", &PressureTensor::computeTensorLocal)
      ;
    }
  }
}
