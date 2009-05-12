#ifndef _LANGEVIN
#define _LANGEVIN

#include "types.hpp"
#include "Thermostat.hpp"

namespace espresso {
  namespace thermostat {

    class Langevin: public Thermostat {
    private:

      real gamma;

      boost::shared_ptr<Property<Real3D> > position;
      boost::shared_ptr<Property<Real3D> > velocity;
      boost::shared_ptr<Property<Real3D> > force;

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

      /** Method to get gamma. */
      real getGamma() const;

      /** Method to set gamma. */
      void setGamma(real _gamma);

      /** Method of the thermostat to modify r and v before force computation. */
      virtual void thermalizeA();

      /** Method of the thermostat to modify r and v before force computation. */
      virtual void thermalizeB();

      virtual ~Langevin();

    };

  }
}

#endif
