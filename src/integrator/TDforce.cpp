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
    : SystemAccess(system) {
        //setFilename(itype, filename);

        center = (0.0,0.0,0.0);

        LOG4ESPP_INFO(theLogger, "TDforce constructed");
    }

    TDforce::~TDforce() {}

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
          LOG4ESPP_DEBUG(theLogger, "thermodynamically forcelize");

          System& system = getSystemRef();

          // iterate over CG particles
          CellList cells = system.storage->getRealCells();
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
              //frictionThermo(*cit);

              /*
              std::map<int, Table>::iterator it;
              it = forces.find(cit->getType());
              Table table = *it;
              */
              Table table = forces.find(cit->getType())->second;

              if (table) {
                  // calculate distance from reference point
                  Real3D dist3D = cit->getPos() - center;
                  real dist = sqrt(dist3D.sqr());

                  // read force from table
                  real ffactor = table->getForce(dist);

                  std::cout << "ffactor of particle " << cit->getId()
                          << " is " << ffactor << "\n";

                  // substract td force
                  cit->force() -= ffactor;
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

      class_<TDforce, shared_ptr<TDforce> >

        ("integrator_TDforce", init< shared_ptr<System> >())
        .add_property("filename", &TDforce::getFilename)
        .def("setCenter", pySetCenter)
        .def("addForce", pyAddForce)
        ;
    }

  }
}
