#include "python.hpp"
#include "FreeEnergyCompensation.hpp"
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

    LOG4ESPP_LOGGER(FreeEnergyCompensation::theLogger, "FreeEnergyCompensation");

    FreeEnergyCompensation::FreeEnergyCompensation(shared_ptr<System> system)
    :Extension(system) {

        type = Extension::FreeEnergyCompensation;

        center = (0.0,0.0,0.0);

        LOG4ESPP_INFO(theLogger, "FreeEnergyCompensation constructed");
    }

    FreeEnergyCompensation::~FreeEnergyCompensation() {}



    void FreeEnergyCompensation::connect(){
        _applyForce = integrator->aftCalcF.connect(
            boost::bind(&FreeEnergyCompensation::applyForce, this));
    }

    void FreeEnergyCompensation::disconnect(){
        _applyForce.disconnect();
    }
    
    
    
    
    


    void FreeEnergyCompensation::addForce(int itype, const char* _filename, int type) {
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


    void FreeEnergyCompensation::applyForce() {
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


    void FreeEnergyCompensation::setCenter(real x, real y, real z){
            center = Real3D(x, y, z);
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void FreeEnergyCompensation::registerPython() {

      using namespace espresso::python;

      void (FreeEnergyCompensation::*pySetCenter)(real x, real y, real z)
                        = &FreeEnergyCompensation::setCenter;

      void (FreeEnergyCompensation::*pyAddForce)(int itype, const char* filename, int type)
                        = &FreeEnergyCompensation::addForce;

      class_<FreeEnergyCompensation, shared_ptr<FreeEnergyCompensation>, bases<Extension> >
        ("integrator_FreeEnergyCompensation", init< shared_ptr<System> >())
        .add_property("filename", &FreeEnergyCompensation::getFilename)
        .def("connect", &FreeEnergyCompensation::connect)
        .def("disconnect", &FreeEnergyCompensation::disconnect)
        .def("setCenter", pySetCenter)
        .def("addForce", pyAddForce)
        ;
    }

  }
}
