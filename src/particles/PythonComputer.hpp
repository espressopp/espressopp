#ifndef _PARTICLES_PYTHONCOMPUTER_HPP
#define _PARTICLES_PYTHONCOMPUTER_HPP

#include <particles/Computer.hpp>
#include <particles/PropertyHandle.hpp>
#include <boost/python.hpp>

namespace espresso {
  namespace particles {
    /** Function object that prints data of a single particle. 
     */
    class PythonComputer:
      public Computer,
      public boost::python::wrapper<Computer> {

    public: // parts invisible to Python

      /** The operator() calling the Python objects function each */
      virtual void operator()(const ParticleHandle pref);

      virtual void prepare(const Storage *);

      static void registerPython();

    private:
      ConstPropertyHandle<ParticleId> particleId;
    };
  }
}

#endif
