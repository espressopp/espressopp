#include "python.hpp"
#include <cmath>
#include "PressureTensor.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "VerletList.hpp"

using namespace espresso;
using namespace iterator;
using namespace interaction;

namespace espresso {
  namespace analysis {
    void PressureTensor::computeTensor(real wij[6]) const {

      System& system = getSystemRef();
  
      // determine number of local particles and total particles
      int N;
      int Nsum;
      N = system.storage->getNRealParticles();
      boost::mpi::reduce(*mpiWorld, N, Nsum, std::plus<int>(), 0);

      // determine volume of the box
      Real3D Li = system.bc->getBoxL();
      real V = Li[0] * Li[1] * Li[2];

      // compute the kinetic contribution (2/3 \sum 1/2mv^2)
      real e_kinetic;
      real p_kinetic;
      real vvLocal[6], vv[6];
      for (size_t j = 0; j < 6; j++) vvLocal[j] = 0.0;

      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real mass = cit->p.mass;
        vvLocal[0] += mass * cit->m.v[0] * cit->m.v[0];
        vvLocal[1] += mass * cit->m.v[1] * cit->m.v[1];
        vvLocal[2] += mass * cit->m.v[2] * cit->m.v[2];
        vvLocal[3] += mass * cit->m.v[0] * cit->m.v[1];
        vvLocal[4] += mass * cit->m.v[0] * cit->m.v[2];
        vvLocal[5] += mass * cit->m.v[1] * cit->m.v[2];
      }
      boost::mpi::reduce(*mpiWorld, vvLocal, 6, vv, std::plus<real>(), 0);

      // compute the short-range nonbonded contribution
      real wijLocal[6];
      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        srIL[j]->computeVirialTensor(wijLocal);
      }
      boost::mpi::reduce(*mpiWorld, wijLocal, 6, wij, std::plus<real>(), 0);
      for (size_t k = 0; k < 6; k++) wij[k] = (vv[k] + wij[k]) / (3.0 * V);
    }

    real PressureTensor::compute() const {
      real wij[6];
      computeTensor(wij);
      return -1.0;
    }

    using namespace boost::python;

    static object wrapCompute(class PressureTensor* obj) {
      real wij[6];
      obj->computeTensor(wij);
      return make_tuple(wij[0], wij[1], wij[2], wij[3], wij[4], wij[5]);
    }

    void PressureTensor::registerPython() {
      using namespace espresso::python;
      class_<PressureTensor, bases< Observable > >
        ("analysis_PressureTensor", init< shared_ptr< System > >())
        .def("compute", &wrapCompute) 
      ;
    }
  }
}
