#ifndef _INTERACTION_FENE_HPP
#define _INTERACTION_FENE_HPP

#include <pmi.hpp>
#include <logging.hpp>
#include <interaction/Interaction.hpp>

namespace espresso {
  namespace interaction {

    /** This class provides routines to computer forces and energies
	based on the FENE potential.

	\f[ V(r) = -\frac{1}{2} \Delta r_\text{max}^2 K \log \left[ 1 - \left(
	\frac{r-r_0}{\Delta r_\text{max}} \right)^2 \right]
        \f]

    */

    class FENE: public Interaction {

    private:

      IF_MPI(pmi::ParallelClass<FENE> pmiObject;)

      real K;       
      real r0;
      real rMax;

      typedef espresso::particleset::ParticleSet::const_reference const_reference;

      static LOG4ESPP_DECL_LOGGER(theLogger);

    public:
 
      IF_PYTHON(static void registerPython();)

      /** Default constructor. Member variables are accessed by setter and getter. */

      FENE(); 
       
      virtual ~FENE();

      virtual real computeEnergy (const Real3D &dist,
				  const const_reference p1,
				  const const_reference p2) const ;

      virtual real computeEnergy (const Real3D &dist) const;
      
      virtual real computeEnergy (const real dist) const;

      virtual real computeEnergySqr (real distSqr) const;
      
      virtual Real3D computeForce (const Real3D &dist,
				   const const_reference p1,
				   const const_reference p2) const;

      virtual Real3D computeForce (const Real3D &dist) const;
      
      virtual void setK(real _K);
      virtual void setr0(real _r0);
      virtual void setrMax(real _rMax);

#ifdef HAVE_MPI
      // PMI Worker routines:
      virtual void setKWorker();
      virtual void setr0Worker();
      virtual void setrMaxWorker();
#endif
      
      /* FENE should probably derived from a two-body interaction
	 without a cutoff to avoid the following */
      virtual real getCutoff() const;
      virtual real getCutoffSqr() const;
      
      // Local functions

    private:

      virtual void setKLocal(real _K);
      virtual void setr0Local(real _r0);
      virtual void setrMaxLocal(real _rMax);

    };
  }
}

#endif
