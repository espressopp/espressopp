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
      PythonComputer(Storage::SelfPtr _storage) 
	: storage(_storage) {}
      
      PythonComputer(Property< Real3D >::SelfPtr prop) 
      { storage = prop->getStorage(); }
      
      PythonComputer(Property< real >::SelfPtr prop) 
      { storage = prop->getStorage(); }

      PythonComputer(Property< int >::SelfPtr prop) 
      { storage = prop->getStorage(); }

      virtual void prepare() {
	// get the particleId property
	particleId = storage->getIdPropertyHandle();
    
	// call Python prepare
	if (override prepare = get_override("prepare"))
	  prepare();
	else Computer::prepare();
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
      Storage::SelfPtr storage;
      ConstPropertyHandle< ParticleId > particleId;
    };
  }
}

void Computer::registerPython() {
  espresso::python::class_
    < espresso::particles::PythonComputer, boost::noncopyable >
    ("particles_PythonComputer", init< Storage::SelfPtr >())
    .def(init< Property< Real3D >::SelfPtr >())
    .def(init< Property< real >::SelfPtr >())
    .def(init< Property< int >::SelfPtr >())
    .def("prepare", &PythonComputer::prepare)
    .def("finalize", &PythonComputer::finalize)
    .def("apply", &PythonComputer::apply)
    ;

}
