#ifndef _PAIRS_PYTHONCOMPUTER_HPP
#define _PAIRS_PYTHONCOMPUTER_HPP

#include <Real3D.hpp>
#include <pairs/Computer.hpp>
#include <pairs/Set.hpp>
#include <particles/ParticleHandle.hpp>
#include <particles/Set.hpp>
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

namespace espresso {
  namespace pairs {
    /** Function object that prints data of a single particle. 
     */
    class PythonComputer:

      public Computer,

      public boost::python::wrapper<PythonComputer> {

    public: // parts visible to Python
      PythonComputer(boost::shared_ptr<Set> _set)
	: set(_set) {}

    public: // parts invisible to Python

      virtual void operator()(const Real3D &dist,
                              const particles::ParticleHandle p1,
                              const particles::ParticleHandle p2);

      static void registerPython();

    private:
      boost::shared_ptr<Set> set;
    };
  }
}

#endif
