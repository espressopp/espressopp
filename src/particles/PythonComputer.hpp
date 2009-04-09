#ifndef _PARTICLES_PYTHONCOMPUTER_HPP
#define _PARTICLES_PYTHONCOMPUTER_HPP

#include <particles/Computer.hpp>
#include <boost/python.hpp>

namespace espresso {
  namespace particles {
    /** Function object that prints data of a single particle. 
     */
    class PythonComputer: public Computer {
    private:
      boost::python::object pyCompute;

    public: // parts visible to Python

      /** Construct a Computer and make a reference back to Python.  */
      PythonComputer(): pyCompute(boost::python::object()) {}
      
      void setCallback(boost::python::object _pyCompute) { pyCompute = _pyCompute; }
      const boost::python::object getCallback() const { return pyCompute; }

    public: // parts invisible to Python

      /** The operator() calling the Python objects function each */
      virtual void operator()(const ParticleHandle pref);

      static void registerPython();
    };
  }
}

#endif
