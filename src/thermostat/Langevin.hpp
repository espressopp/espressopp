#ifndef _LANGEVIN
#define _LANGEVIN

#include <boost/signals2.hpp>
#include <boost/random.hpp>

#include "logging.hpp"
#include "types.hpp"
#include "Thermostat.hpp"
#include "integrator/VelocityVerlet.hpp"

using namespace boost;
using namespace boost::random;

namespace espresso {
  namespace thermostat {

    class Langevin: public Thermostat {

    private:
  
      static LOG4ESPP_DECL_LOGGER(theLogger);

      real gamma;

      boost::shared_ptr<Property<Real3D> > position;
      boost::shared_ptr<Property<Real3D> > velocity;
      boost::shared_ptr<Property<Real3D> > force;
      boost::shared_ptr<integrator::VelocityVerlet> integrator;

      boost::signals2::connection stepA;
      boost::signals2::connection stepB;

      // random number generator for normal distribution

      boost::minstd_rand0 linearCongruential;
      boost::normal_distribution<double> normalDist;
      variate_generator<boost::minstd_rand0&, 
                        boost::normal_distribution<double> > gauss;

    public:

      /** Method to register the Python bindings. */
      static void registerPython();

      /** Constructor with the particle set, temperature, and friction coefficient. */
      Langevin(boost::shared_ptr<particles::Set> _particles,
               real _temperature,
               real _gamma,
	       boost::shared_ptr<Property<Real3D> > _position,
               boost::shared_ptr<Property<Real3D> > _velocity,
               boost::shared_ptr<Property<Real3D> > _force);

      /** The thermostat will be connected to a VelocityVerlet integrator so that it
          can be called by this integrator. 
      */

      void connect(boost::shared_ptr<integrator::VelocityVerlet> integrator);

      void disconnect();

      /** Method to get gamma. */
      real getGamma() const;

      /** Method to set gamma. */
      void setGamma(real _gamma);

      /** Method of the thermostat to modify r and v before force computation. */
      virtual void thermalizeA();

      /** Method of the thermostat to modify r and v before force computation. */
      virtual void thermalizeB(int itimestep);

      virtual ~Langevin();

    };

  }
}

#endif
