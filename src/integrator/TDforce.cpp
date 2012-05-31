#include "python.hpp"
#include "TDforce.hpp"
#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/InterpolationLinear.hpp"
#include "interaction/InterpolationAkima.hpp"
#include "interaction/InterpolationCubic.hpp"


namespace espresso {
  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(TDforce::theLogger, "TDforce");

    TDforce::TDforce(shared_ptr<System> system)
    :Extension(system) {

        type = Extension::TDforce;

        center = (0.0,0.0,0.0);

        LOG4ESPP_INFO(theLogger, "TDforce constructed");
    }

    TDforce::~TDforce() {}



    void TDforce::connect(){
        _applyForce = integrator->aftCalcF.connect(
            boost::bind(&TDforce::applyForce, this));
    }

    void TDforce::disconnect(){
        _applyForce.disconnect();
    }


    void TDforce::addForce(int itype, const char* _filename, int type) {
        boost::mpi::communicator world;
        filename = _filename;
        Table table;

        if (itype == 1) { // create a new InterpolationLinear
            table = make_shared <interaction::InterpolationLinear> ();
            table->read(world, _filename);
        }

        else if (itype == 2) { // create a new InterpolationAkima
            table = make_shared <interaction::InterpolationAkima> ();
            table->read(world, _filename);
        }

        else if (itype == 3) { // create a new InterpolationCubic
            table = make_shared <interaction::InterpolationCubic> ();
            table->read(world, _filename);
        }

        forces.insert(std::make_pair(type,table));
    }


    void TDforce::applyForce() {
          LOG4ESPP_DEBUG(theLogger, "apply TD force");

          System& system = getSystemRef();

          // iterate over CG particles
          CellList cells = system.storage->getRealCells();
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {

              Table table = forces.find(cit->getType())->second;

              if (table) {
                  // calculate distance from reference point
                  Real3D dist3D = cit->getPos() - center;
                  real dist = sqrt(dist3D.sqr());

                  // read fforce from table
                  real fforce = table->getForce(dist);
                  fforce /= dist;

                  // substract td force
                  cit->force() -= (dist3D * fforce);

                  /*
                  // use this if you need 1-dir force only!
                  real d1 = cit->getPos()[0] - center[0];
                  real force = table->getForce(d1);
                  cit->force()[0] -= force;
                  */
              }
          }
    }


    void TDforce::setCenter(real x, real y, real z){
            center = Real3D(x, y, z);
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void TDforce::registerPython() {

      using namespace espresso::python;

      void (TDforce::*pySetCenter)(real x, real y, real z)
                        = &TDforce::setCenter;

      void (TDforce::*pyAddForce)(int itype, const char* filename, int type)
                        = &TDforce::addForce;

      class_<TDforce, shared_ptr<TDforce>, bases<Extension> >
        ("integrator_TDforce", init< shared_ptr<System> >())
        .add_property("filename", &TDforce::getFilename)
        .def("connect", &TDforce::connect)
        .def("disconnect", &TDforce::disconnect)
        .def("setCenter", pySetCenter)
        .def("addForce", pyAddForce)
        ;
    }

  }
}
