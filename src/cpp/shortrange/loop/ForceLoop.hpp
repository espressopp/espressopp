#ifndef _SHORTRANGE_LOOP_FORCELOOP_HPP
#define _SHORTRANGE_LOOP_FORCELOOP_HPP

namespace espresso {
  namespace shortrange {
    //////////////////////////////////////////////////
    // Interactions
    //////////////////////////////////////////////////
    class Interaction {
    public:
      real computeEnergy(Real3D &dist, ParticleRef p1, ParticleRef p2) const;
      Real3D computeForce(Real3D &dist, ParticleRef p1, ParticleRef p2) const;
    };

    class LennardJonesInteraction : public Interaction {
    public:
      real computeEnergy(Real3D &dist, ParticleRef p1, ParticleRef p2) const;
      Real3D computeForce(Real3D &dist, ParticleRef p1, ParticleRef p2) const;
    };

    //////////////////////////////////////////////////
    // Computers
    //////////////////////////////////////////////////
    class Computer {
    public:
      virtual void compute(const Real3D dist, 
			   const ParticleRef p1, 
			   const ParticleRef p2) = 0;
    };

    template < class Interaction >
    class ForceComputer : public Computer {
      const ArrayPropertyRef<real,3> force;
      const Interaction &interaction;

      Real3D pressure;
      bool computesPressure;

    public:
      virtual void compute(const Real3D dist, 
			   const ParticleRef p1, 
			   const ParticleRef p2) {
	Real3D f = interaction.computeForce(dist, p1, p2);
	force[p1] += f;
	force[p2] -= f;
	if (computesPressure) pressure += f*dist;
      }
    };

    template < class Interaction >
    class EnergyComputer : public Computer {
      real totalEnergy;
      const PropertyRef< real > energy;
      const Interaction &interaction;

    public:
      virtual void compute(const Real3D dist, 
			   const ParticleRef p1, 
			   const ParticleRef p2) {
	real e = interaction.computeEnergy(dist, p1, p2);
	energy[p1] += 0.5*e;
	energy[p2] += 0.5*e;
	totalEnergy += e;
      }
    };

    class ComputerList : public Computer {
      vector < Computer* > compList;
    public:
      virtual void compute(const Real3D dist, 
			   const ParticleRef p1, 
			   const ParticleRef p2) {
	for (vector< Computer* >::iterator it = compList.begin();
	     it != compList.end(); it++)
	  it->compute(dist, p1, p2);
      }
    };

  }
}
#endif


// the main loop:
// ParticlePairCollection
// foreach(pc, comp.compute);

// Typical pair collection:
// typedef Collection< Pair < ParticleRef > > PairCollection;

// Interface:
// template < class valueType >
// class CollectionBaseIterator {
//   CollectionBaseIterator operator++(int);
//   bool operator!=(CollectionBaseIterator &other);
// };
//
// template < class valueType >
// class CollectionBase {
//   typedef CollectionBaseIterator iterator;
//   iterator begin() const;
//   iterator end() const;
// };
//
//
//
// A pair is an STL pair.
//
// A ParticleRef should be cheaply copyable.
//
