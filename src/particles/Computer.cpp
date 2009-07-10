#include "Computer.hpp"
#include "Storage.hpp"
#include "particles/PropertyHandle.hpp"
#include "python.hpp"

using namespace espresso;
using namespace espresso::particles;
using namespace espresso::python;

class PythonComputer
  : public Computer,
    public wrapper< Computer > {
  
public: // parts invisible to Python
  PythonComputer(Storage::SelfPtr _storage) : storage(_storage) {}

  virtual void prepare() {
    // get the particleId property
    particleId = storage->getIdPropertyHandle();
    
    // call Python prepare
    if (override prepare = get_override("prepare"))
      prepare();
    else Computer::prepare();
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
  Storage::SelfPtr storage;
  ConstPropertyHandle< ParticleId > particleId;
};

void Computer::registerPython() {
  using namespace espresso::python;
  class_< PythonComputer, boost::noncopyable >
    ("particles_PythonComputer", init< Storage::SelfPtr >())
    .def("prepare", &PythonComputer::prepare)
    .def("finalize", &PythonComputer::finalize)
    .def("__apply__", &PythonComputer::operator())
    ;

}
