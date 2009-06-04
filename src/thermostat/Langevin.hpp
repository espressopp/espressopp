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

    /** A Langevin thermostat models the embedding of particles into a
       viscouls liquid of constant temperature using kicks and drag.

       The Langevin thermostat is connected to a Velocity-Verlet 
       integrator (\sa espresso::integrator::VelocityVerlet).

    */

    class Langevin: public Thermostat {

    private:
  
      static LOG4ESPP_DECL_LOGGER(theLogger);

      real gamma;

      boost::shared_ptr<Property<Real3D> > position;
      boost::shared_ptr<Property<Real3D> > velocity;
      boost::shared_ptr<Property<Real3D> > force;

      // Boolean variables that is true for an existing connection to
      // an integrator. At any time only one connection is allowed.

      bool connected;

      // Connections to the integrator after step A and step B 

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

          \param thisSharedPtr must be the shared pointer to this object
          \param integrator is the Velocity-Verlet integrator to which this thermostat will be connected.

          This thermostat will not keep a reference to the integrator but to its connections.

      */

      void connect(boost::shared_ptr<Langevin> thisSharedPtr,
                   boost::shared_ptr<integrator::VelocityVerlet> integrator);

      /** Disconnect this thermostat object from its integrator. */
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
