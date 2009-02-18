//method to compute the potential energy

#include "LennardJones.hpp"
#include <mpi.hpp>
#include <python.hpp>
#include <iostream>

using namespace espresso::interaction;

/* ---------------------------------------------------------------------- */

LOG4ESPP_LOGGER(LennardJones::theLogger, "interaction.LennardJones");

/* ---------------------------------------------------------------------- */

namespace espresso {
  namespace interaction {
    LennardJones::LennardJones() {
      epsilon = 1.0;
      sigma = 1.0;
      cutoff = 2.5;
      cutoffSqr = 2.5*2.5;
    }

    LennardJones::~LennardJones() {}

    real LennardJones::computeEnergy (const Real3D &dist,
				      const const_reference p1,
				      const const_reference p2) const {
      return computeEnergy(dist);
    }

    real LennardJones::computeEnergy (const Real3D &dist) const {
      return computeEnergySqr(dist.sqr());
    }

    real LennardJones::computeEnergy(const real dist) const {
      return computeEnergySqr(dist*dist);
    }
    
    real LennardJones::computeEnergySqr (const real distSqr) const {
      if (distSqr < cutoffSqr) {
	real frac2 = sigma*sigma / distSqr;
	real frac6 = frac2 * frac2 * frac2;
	real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
	return energy;
      } else return 0.0;
    }
    
    Real3D LennardJones::computeForce (const Real3D &dist,
					       const const_reference p1,
					       const const_reference p2) const {
	return computeForce(dist);
    }

    Real3D LennardJones::computeForce (const Real3D &dist) const {
      Real3D f = 0.0;
      real   frac2;
      real   frac6;
      
      real distSqr = dist.sqr();
      
      if (distSqr < cutoffSqr) {
	frac2 = sigma / distSqr;
	frac6 = frac2 * frac2 * frac2;
	real ffactor = 48.0 * epsilon * (frac6*frac6 - 0.5 * frac6) * frac2;
	
	LOG4ESPP_DEBUG(theLogger, "computeForce, distSqr = " << distSqr <<
		       ", ffactor = " << ffactor);
	
	f = dist * ffactor;
      } 
      
      return f;
    }
    
    real LennardJones::getCutoff() const { return cutoff; }
    real LennardJones::getCutoffSqr() const { return cutoffSqr; }
    void LennardJones::setCutoff(real _cutoff) { 
#ifdef HAVE_MPI
      pmiObject.invoke<&LennardJones::setCutoffWorker>();
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _cutoff, pmi::getControllerMPIRank());
#endif
      setCutoffLocal(_cutoff);
    }

#ifdef HAVE_MPI
    void LennardJones::setCutoffWorker() {
      real _cutoff;
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _cutoff, pmi::getControllerMPIRank());
      setCutoffLocal(_cutoff);
    }
#endif

    void LennardJones::setCutoffLocal(real _cutoff) {
      cutoff = _cutoff; 
      cutoffSqr = cutoff * cutoff;
    }
    
    real LennardJones::getEpsilon() const { return epsilon; }
    void LennardJones::setEpsilon(real _epsilon) { 
#ifdef HAVE_MPI
      pmiObject.invoke<&LennardJones::setEpsilonWorker>();
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _epsilon, pmi::getControllerMPIRank());
#endif
      setEpsilonLocal(_epsilon);
    }

#ifdef HAVE_MPI
    void LennardJones::setEpsilonWorker() {
      real _epsilon;
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _epsilon, pmi::getControllerMPIRank());
      setEpsilonLocal(_epsilon);
    }
#endif    

    void LennardJones::setEpsilonLocal(real _epsilon) {
      epsilon = _epsilon;
    }

    real LennardJones::getSigma() const { return sigma; }
    void LennardJones::setSigma(real _sigma) { 
#ifdef HAVE_MPI
      pmiObject.invoke<&LennardJones::setSigmaWorker>();
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _sigma, pmi::getControllerMPIRank());
#endif
      setSigmaLocal(_sigma);
    }
    
#ifdef HAVE_MPI
    void LennardJones::setSigmaWorker() {
      real _sigma;
      boost::mpi::communicator world;
      boost::mpi::broadcast(world, _sigma, pmi::getControllerMPIRank());
      setSigmaLocal(_sigma);
    }
#endif    

    void LennardJones::setSigmaLocal(real _sigma) {
      sigma = _sigma;
    }


  }


  //////////////////////////////////////////////////
  // REGISTRATION WITH PMI
  //////////////////////////////////////////////////
  PMI_REGISTER_CLASS("espresso::interaction::LennardJones", espresso::interaction::LennardJones);
  PMI_REGISTER_METHOD("setCutoffWorker", espresso::interaction::LennardJones, setCutoffWorker);
  PMI_REGISTER_METHOD("setEpsilonWorker", espresso::interaction::LennardJones, setEpsilonWorker);
  PMI_REGISTER_METHOD("setSigmaWorker", espresso::interaction::LennardJones, setSigmaWorker);


#ifdef HAVE_PYTHON  
  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
//   Real3D LennardJones::computeForceOverload(const Real3D& dist) const {
//     return LennardJones::computeForce(dist);
//   }

  void 
  LennardJones::registerPython() {
    using namespace boost::python;

    // create thin wrappers around overloaded member functions
    Real3D (LennardJones::*computeForceOverload)(const Real3D&) const =
      &LennardJones::computeForce;
    real (LennardJones::*computeEnergyOverload1)(const Real3D &) const =
      &LennardJones::computeEnergy;
    real (LennardJones::*computeEnergyOverload2)(const real) const =
      &LennardJones::computeEnergy;
    
    class_<LennardJones>("interaction_LennardJones", init<>())
      .def("getCutoff", &LennardJones::getCutoff)
      .def("setCutoff", &LennardJones::setCutoff)
      .def("getEpsilon", &LennardJones::getEpsilon)
      .def("setEpsilon", &LennardJones::setEpsilon)
      .def("getSigma", &LennardJones::getSigma)
      .def("setSigma", &LennardJones::setSigma)
      .def("computeForce", computeForceOverload)
      .def("computeEnergy", computeEnergyOverload1)
      .def("computeEnergy", computeEnergyOverload2);
  }
#endif


}
