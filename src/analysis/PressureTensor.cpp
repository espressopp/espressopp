#include "python.hpp"
#include <cmath>
#include "PressureTensor.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "VerletList.hpp"

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
      //boost::mpi::reduce(*mpiWorld, vvLocal.get(), 6, vv.get(), std::plus<real>(), 0); //reduce is uncorrect
      boost::mpi::all_reduce(*mpiWorld, vvLocal.get(), 6, vv.get(), std::plus<real>());

      // compute the short-range nonbonded contribution
      Tensor wijLocal(0.0);
      Tensor wij;

      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        srIL[j]->computeVirialTensor(wijLocal);
      }

      //boost::mpi::reduce(*mpiWorld, wijLocal.get(), 6, wij.get(), std::plus<real>(), 0); //reduce is uncorrect
      boost::mpi::all_reduce(*mpiWorld, wijLocal.get(), 6, wij.get(), std::plus<real>());
      return (vv + wij) / V;
    }

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
      ;
    }
  }
}
