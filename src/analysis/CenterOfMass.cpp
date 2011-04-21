#include "python.hpp"
#include <cmath>
#include "CenterOfMass.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"

using namespace espresso;
using namespace iterator;

namespace espresso {
  namespace analysis {
    Real3D CenterOfMass::computeVector() const {

      System& system = getSystemRef();
  
      // compute the center-of-mass of the system (r_CM = \sum m_i r_i / \sum m_i)
      real xcom = 0.0;
      real ycom = 0.0;
      real zcom = 0.0;
      real mass = 0.0;
      real xcom_sum = 0.0;
      real ycom_sum = 0.0;
      real zcom_sum = 0.0;
      real mass_sum = 0.0;

      CellList realCells = system.storage->getRealCells();
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        const Particle& p = *cit;
        xcom += p.mass() * p.position()[0];
        ycom += p.mass() * p.position()[1];
        zcom += p.mass() * p.position()[2];
        mass += p.mass();
      }
      boost::mpi::reduce(*mpiWorld, xcom, xcom_sum, std::plus<real>(), 0);
      boost::mpi::reduce(*mpiWorld, ycom, ycom_sum, std::plus<real>(), 0);
      boost::mpi::reduce(*mpiWorld, zcom, zcom_sum, std::plus<real>(), 0);
      boost::mpi::reduce(*mpiWorld, mass, mass_sum, std::plus<real>(), 0);

      //Real3D force(0.0, 0.0, 0.0);
      return Real3D(xcom / mass_sum, ycom / mass_sum, zcom / mass_sum);
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real CenterOfMass::compute() const {
      return -1.0;
    }

    void CenterOfMass::registerPython() {
      using namespace espresso::python;
      class_<CenterOfMass, bases< Observable > >
        ("analysis_CenterOfMass", init< shared_ptr< System > >())
        .def("compute", &CenterOfMass::computeVector)
      ;
    }
  }
}
