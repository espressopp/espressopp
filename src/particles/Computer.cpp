#include "particles/Computer.hpp"
#include "storage/Storage.hpp"
#include "storage/PropertyHandle.hpp"
#include "Property.hpp"
#include "python.hpp"

using namespace espresso;
using namespace espresso::storage;
using namespace espresso::particles;
using namespace espresso::python;

namespace espresso {
  namespace particles {
    class PythonComputer
      : public Computer,
	public wrapper< Computer > {
      
    public: 
      void prepare(Storage::SelfPtr storage) {
	// get the particleId property
	particleId = storage->getIdPropertyHandle();
    
	// call Python prepare
	if (override prepare = get_override("prepare"))
	  prepare(storage);
      }

      /** The operator() calling the Python objects function each */
      bool apply(const ParticleHandle pref) {
	return get_override("apply")(particleId[pref]);
      }
  
      void finalize() {
	// call Python finalize
	if (override finalize = get_override("finalize"))
	  finalize();
	else Computer::finalize();
      }

    private:
      PropertyHandle< ParticleId > particleId;
    };
  }
}

void Computer::registerPython() {
  espresso::python::class_
    < espresso::particles::PythonComputer, boost::noncopyable >
    ("particles_PythonComputer")
    .def("prepare", pure_virtual(&Computer::prepare))
    .def("finalize", &Computer::finalize)
    .def("apply", pure_virtual(&Computer::apply))
    ;

}
