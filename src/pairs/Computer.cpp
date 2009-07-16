#include "Particle.hpp"
#include "python.hpp"
#include "pairs/Computer.hpp"
#include "particles/Storage.hpp"

using namespace espresso;
using namespace espresso::pairs;
using namespace espresso::python;

using espresso::particles::Storage;
using espresso::particles::ParticleHandle;
using espresso::particles::ConstPropertyHandle;

// We need to put the wrapper class into the namespace, even though it
// is used only in this file. Otherwise, we get a strange boost runtim
// error.
namespace espresso {
  namespace pairs {
    class PythonComputer
      : public Computer,
	public wrapper< Computer > {
  
    public: 
      PythonComputer(Storage::SelfPtr _storage1, 
		     Storage::SelfPtr _storage2)
	: storage1(_storage1), storage2(_storage2) {}

      PythonComputer(Storage::SelfPtr _storage)
	: storage1(_storage), storage2(_storage) {}
  
      virtual void prepare() {
	// get the particleId properties
	idProperty1 = storage1->getIdPropertyHandle();
	idProperty2 = storage2->getIdPropertyHandle();

	// call Python prepare
	if (override prepare = get_override("prepare"))
	  prepare();
	else Computer::prepare();
      }

      virtual void apply(const Real3D dist,
			 const ParticleHandle p1,
			 const ParticleHandle p2) {
	get_override("apply")(dist, idProperty1[p1], idProperty2[p2]);
      }

      virtual void finalize() {
	// call Python finalize
	if (override finalize = get_override("finalize"))
	  finalize();
	else Computer::finalize();
      }

    private:
      Storage::SelfPtr storage1;
      Storage::SelfPtr storage2;
      ConstPropertyHandle< ParticleId > idProperty1;
      ConstPropertyHandle< ParticleId > idProperty2;
    };
  }
}

void Computer::registerPython() {
  espresso::python::class_
    < espresso::pairs::PythonComputer, boost::noncopyable >
    ("pairs_PythonComputer", 
     init< Storage::SelfPtr, optional< Storage::SelfPtr > >())
    .def("prepare", &PythonComputer::prepare)
    .def("finalize", &PythonComputer::finalize)
    .def("apply", &PythonComputer::apply)
    ;
}

