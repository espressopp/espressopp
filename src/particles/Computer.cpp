#include "Computer.hpp"
#include "particles/PropertyHandle.hpp"
#include "python.hpp"

using namespace espresso;
using namespace espresso::particles;
using namespace espresso::python;

class PythonComputer
  : public Computer,
    public wrapper< Computer > {
  
public: // parts invisible to Python
  
  virtual void prepare(Storage::SelfPtr storage) {
    // get the particleId property
    particleId = storage->getIdPropertyHandle();
    
    // call Python prepare
    if (override prepare = get_override("prepare"))
      prepare();
    else Computer::prepare(storage);
  }

  /** The operator() calling the Python objects function each */
  virtual void operator()(const ParticleHandle pref) {
    get_override("__apply__")(particleId[pref]);
  }
  
  virtual void finalize() {
    // call Python finalize
    if (override finalize = get_override("finalize"))
      finalize();
    else Computer::finalize();
  }

private:
  ConstPropertyHandle<ParticleId> particleId;
};

void Computer::registerPython() {
  using namespace espresso::python;
  class_< PythonComputer, boost::noncopyable >
    ("particles_PythonComputer", init<>())
    .def("prepare", &PythonComputer::prepare)
    .def("finalize", &PythonComputer::finalize)
    ;

}
