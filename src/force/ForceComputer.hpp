#ifndef _FORCE_FORCECOMPUTER_HPP
#define _FORCE_FORCECOMPUTER_HPP

#include <boost/enable_shared_from_this.hpp>
#include <boost/signals2.hpp>

#include "types.hpp"
#include "integrator/MDIntegrator.hpp"
#include "potential/Potential.hpp"
#include "pairs/Set.hpp"

namespace espresso {
  namespace force {

    /** A ForceComputer updates the forces by applying an potential to a set of particle pairs.

        A ForceComputer is constructed by an potential and a set of particle pairs. It must be
        connected to an integrator so that forces will be updated during integration.

    */

    class ForceComputer 
      : public boost::enable_shared_from_this< ForceComputer >
    {
    public:
      typedef boost::shared_ptr< ForceComputer > SelfPtr;

    private:

      potential::Potential::SelfPtr potential;
      pairs::Set::SelfPtr pairs;

      // variable that holds the connection to the integrator
      // At this time the ForceComputer can only connect to one integrator.

      boost::signals2::connection forceCalc;

      static LOG4ESPP_DECL_LOGGER(theLogger);

      // This method will be connect to and called by the integrator

      void updateForces(const integrator::MDIntegrator& integrator);

    public:

      /** Constructor of ForceComputer just takes the potential and the pair set */

      ForceComputer(potential::Potential::SelfPtr potential,
                    pairs::Set::SelfPtr pairs);

      /** This method connects this ForceComputer via its shared pointer to the integrator. */

      void connect(integrator::MDIntegrator::SelfPtr integrator);

      /** This method disconnects the ForceComputer from the integrator. This is also done
          automatically when the integrator is deleted.
      */

      void disconnect();

      static void registerPython();

      virtual ~ForceComputer();

   };
  }
}

#endif
