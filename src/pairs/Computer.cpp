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
      virtual void prepare(Storage::SelfPtr storage1,
			   Storage::SelfPtr storage2) {
	// get the particleId properties
	id1 = storage1->getIdPropertyHandle();
	id2 = storage2->getIdPropertyHandle();

	// call Python prepare
	if (override prepare = get_override("prepare"))
	  prepare(storage1, storage2);
	else Computer::prepare(storage1, storage2);
      }

      virtual void apply(const Real3D dist,
			 const ParticleHandle p1,
			 const ParticleHandle p2) {
	get_override("apply")(dist, id1[p1], id2[p2]);
      }

      virtual void finalize() {
	// call Python finalize
	if (override finalize = get_override("finalize"))
	  finalize();
	else Computer::finalize();
      }

    private:
      ConstPropertyHandle< ParticleId > id1;
      ConstPropertyHandle< ParticleId > id2;
    };
  }
}

void Computer::registerPython() {
  espresso::python::class_
    < espresso::pairs::PythonComputer, boost::noncopyable >
    ("pairs_PythonComputer")
    .def("prepare", &PythonComputer::prepare)
    .def("finalize", &PythonComputer::finalize)
    .def("apply", &PythonComputer::apply)
    ;
}

