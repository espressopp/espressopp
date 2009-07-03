#ifndef _PAIRS_PYTHONCOMPUTER_HPP
#define _PAIRS_PYTHONCOMPUTER_HPP

#include <Real3D.hpp>
#include <pairs/Computer.hpp>
#include <particles/Storage.hpp>
#include <particles/ParticleHandle.hpp>
#include <boost/python.hpp>

namespace espresso {
  namespace pairs {
    /** Function object that prints data of a single particle. 
     */
    class PythonComputer:
      public Computer,
      public boost::python::wrapper< PythonComputer > {

    public: // parts visible to Python
      PythonComputer(particles::Storage::SelfPtr _storage)
	: storage(_storage) {}
      
    public: // parts invisible to Python
      
      virtual void operator()(const Real3D &dist,
                              const particles::ParticleHandle p1,
                              const particles::ParticleHandle p2);
      
      static void registerPython();
      
    private:
      particles::Storage::SelfPtr storage;
    };
  }
}

#endif
