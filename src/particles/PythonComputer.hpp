#ifndef _PARTICLES_PYTHONCOMPUTER_HPP
#define _PARTICLES_PYTHONCOMPUTER_HPP

#include <particles/Computer.hpp>
#include <boost/python.hpp>

namespace espresso {
  namespace particles {
    /** Function object that prints data of a single particle. 
     */
    class PythonComputer:
      public Computer,
      public boost::python::wrapper<PythonComputer> {

    public: // parts visible to Python
      PythonComputer(PStorage _storage)
	: storage(_storage) {}

    public: // parts invisible to Python

      /** The operator() calling the Python objects function each */
      virtual void operator()(const ParticleHandle pref);

      static void registerPython();

    private:
      PStorage storage;
    };
  }
}

#endif
