#include "python.hpp"
#include <cmath>
#include "PressureTensor.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "VerletList.hpp"
#include "mpi.hpp"

using namespace espresso;
using namespace iterator;
using namespace interaction;

namespace espresso {
  namespace analysis {
    real PressureTensor::compute() const {

      System& system = getSystemRef();
  
      // determine number of local particles and total particles
      int N;
      int Nsum;
      N = system.storage->getNRealParticles();
      boost::mpi::reduce(*mpiWorld, N, Nsum, std::plus<int>(), 0);

      // determine volume of the box
      Real3D Li = system.bc->getBoxL();
      real V = Li[0] * Li[1] * Li[2];

      // compute the kinetic contriubtion (2/3 \sum 1/2mv^2)
      real e_kinetic;
      real p_kinetic;
      real vxvx = 0.0;
      real vyvy = 0.0;
      real vzvz = 0.0;
      real vxvy = 0.0;
      real vxvz = 0.0;
      real vyvz = 0.0;
      real vxvx_sum;

      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        vxvx += cit->m.v[0] * cit->m.v[0];
        vyvy += cit->m.v[1] * cit->m.v[1];
        vzvz += cit->m.v[2] * cit->m.v[2];
        vxvy += cit->m.v[0] * cit->m.v[1];
        vxvz += cit->m.v[0] * cit->m.v[2];
        vyvz += cit->m.v[1] * cit->m.v[2];
      }

      boost::mpi::reduce(*mpiWorld, vxvx, vxvx_sum, std::plus<real>(), 0);
      e_kinetic = 0.5 * vxvx_sum;
      p_kinetic = 2.0 * e_kinetic / (3.0 * V);

      // compute the short-range nonbonded contribution
      // loop over interaction types
      real rij_dot_Fij = 0.0;
      real wij[6];
      const InteractionList& srIL = system.shortRangeInteractions;
      for (size_t j = 0; j < srIL.size(); j++) {
        srIL[j]->computeVirialTensor(wij);
        rij_dot_Fij += wij[3];
      }
      real p_nonbonded = rij_dot_Fij / (3.0 * V);

      return (p_kinetic + p_nonbonded);
    }

    void PressureTensor::registerPython() {
      using namespace espresso::python;
      class_<PressureTensor, bases< Observable > >
        ("analysis_PressureTensor", init< shared_ptr< System > >())
      ;
    }
  }
}
