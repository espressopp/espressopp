#include "particles/Computer.hpp"
#include "particles/Storage.hpp"
#include "particles/PropertyHandle.hpp"
#include "Property.hpp"
#include "python.hpp"

using namespace espresso;
using namespace espresso::particles;
using namespace espresso::python;

namespace espresso {
  namespace particles {
    class PythonComputer
      : public Computer,
	public wrapper< Computer > {
      
    public: 
      virtual void prepare(Storage::SelfPtr storage) {
	// get the particleId property
	particleId = storage->getIdPropertyHandle();
    
	// call Python prepare
	if (override prepare = get_override("prepare"))
	  prepare(storage);
	else Computer::prepare(storage);
      }

      /** The operator() calling the Python objects function each */
      virtual void apply(const ParticleHandle pref) {
	get_override("apply")(particleId[pref]);
      }
  
      virtual void finalize() {
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
    .def("prepare", &Computer::prepare)
    .def("finalize", &Computer::finalize)
    .def("apply", pure_virtual(&Computer::apply))
    ;

}
