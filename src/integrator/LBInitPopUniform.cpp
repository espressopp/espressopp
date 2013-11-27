#include "python.hpp"
#include "LBInitPopUniform.hpp"

namespace espresso {
  namespace integrator {
//    LOG4ESPP_LOGGER(LBInitPopUniform::theLogger, "LBInitPopUniform");
    LBInitPopUniform::LBInitPopUniform(shared_ptr<System> system,
                                       shared_ptr< LatticeBoltzmann > latticeboltzmann)
    : LBInit(system, latticeboltzmann) {
    }

    void LBInitPopUniform::createDenVel (real _rho0, Real3D _u0) {
      printf("Creating an initial configuration with uniform density %f and \nvelocity %f, %f, %f\n",
          _rho0, _u0.getItem(0), _u0.getItem(1), _u0.getItem(2));
      printf("-------------------------------------\n");

      real invCs2 = 1. / latticeboltzmann->getCs2();

      real invCs4 = invCs2*invCs2;
      real scalp, value;
      Int3D _Ni = latticeboltzmann->getNi();

      // set initial velocity of the populations from Maxwell's distribution
      for (int i = 0; i < _Ni.getItem(0); i++) {
      // test the damping of a sin-like initial velocities:
        real trace = _u0*_u0*invCs2;
        for (int j = 0; j < _Ni.getItem(1); j++) {
          for (int k = 0; k < _Ni.getItem(2); k++) {
            for (int l = 0; l < latticeboltzmann->getNumVels(); l++) {
              scalp = _u0 * latticeboltzmann->getCi(l);
              value = 0.5 * latticeboltzmann->getEqWeight(l) * _rho0 * (2. + 2. * scalp * invCs2 + scalp * scalp * invCs4 - trace);
              latticeboltzmann->setLBFluid(Int3D(i,j,k),l,value);
              latticeboltzmann->setGhostFluid(Int3D(i,j,k),l,0.0);
            }
          }
        }
      }
    }

    /* do nothing with external forces */
    void LBInitPopUniform::setForce (Real3D _force) {}
    void LBInitPopUniform::addForce (Real3D _force) {}

    void LBInitPopUniform::registerPython() {
      using namespace espresso::python;

      class_<LBInitPopUniform, bases< LBInit > >
          ("integrator_LBInit_PopUniform", init< shared_ptr< System >,
                                             shared_ptr< LatticeBoltzmann > >())
          .def("createDenVel", &LBInitPopUniform::createDenVel)
      ;
    }
  }
}
