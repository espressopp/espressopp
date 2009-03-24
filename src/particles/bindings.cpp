#include "bindings.hpp"
#include "types.hpp"

#include <boost/python.hpp>

#include <particles/Storage.hpp>
#include <particles/Computer.hpp>

using namespace boost::python;
using namespace espresso::base;

namespace espresso {
 namespace particles {
    /** Function object that prints data of a single particle. 
     */
    class ParticleWriter : public Computer {
    private:
      PropertyReference<Real3D>
      pos;  //<! reference to the position vector of all particles.
      ConstPropertyReference<ParticleId>
      id;   //<! reference to the identification vector of all particles.

    public:
      /** Construct a writer for a particle storage.
          \param particleStorage is needed to get access the property vectors of all particles.
      */
      ParticleWriter(Storage &particleStorage, PropertyId position) :
        pos(particleStorage.getPropertyReference<Real3D>(position)),
        id(particleStorage.getIDProperty())
      {}

      /** Function that is applied to a read-only particle.
          \sa espresso::particlestorage::Particlecomputer::operator()
      */
      virtual void operator()(const ParticleReference pref) {
        printf("Particle : id = %ld, pos = (%f,%f,%f)\n",
               size_t(id[pref]), pos[pref][0], pos[pref][1], pos[pref][2]);
      }
    };
 
    /** Function object that prints data of a single particle. 
     */

    class PythonComputer : public Computer {

    private:

       object pyCompute;
       int counter;

    public:

      /** Construct a Computer and make a reference back to Python.  */

      PythonComputer()

      {  pyCompute = object();
         counter = 0;
      }

      void setCallback(object pyCompute) {

         this->pyCompute = pyCompute;
      }

      /** The operator() has to be implemented in C++ but calls Python routine.  */

      virtual void operator()(const ParticleReference pref) {

         pyCompute.attr("each")(counter++);
      }
    };
 
    void registerPython() {

      Storage::registerPython();

      // Computer has to be exported otherwise Python will not accept
      // something as argument for foreach in Storage

      // noncopyable must be specified as it is an abstract class

      class_<Computer, boost::noncopyable> ("particles_Computer", no_init)
      ;

      class_<ParticleWriter, bases<Computer> >
         ("particles_ParticleWriter", init<Storage&, PropertyId>())
      ;
     
      class_<PythonComputer, bases<Computer> >
         ("particles_PythonComputer", init<>())
      .def("setCallback", &PythonComputer::setCallback)
      ;
    }

  }
}

