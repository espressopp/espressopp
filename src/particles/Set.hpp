#ifndef _PARTICLES_SET_HPP
#define _PARTICLES_SET_HPP

#include "types.hpp"
#include "particles/Storage.hpp"

namespace espresso {
  namespace particles {
    /** MOCK particle set. Provides a view onto a set of particles
        from a ParticleStorage
    */
    class Set {
    protected:
      /// the storage our particles are stored in
      Storage *theStorage;

    public:
      /** base constructor
          
          @param _store pointer to the Storage the
          particles in this set come from
      */
      Set(Storage *_store = 0): theStorage(_store) {}
      virtual ~Set() {}

      /** for a particle of the Storage of this class,
          check whether it belongs to this set
      */
      virtual bool isMember(ParticleHandle pref) const = 0;

      /** apply computer to all particles of this set
       */
      virtual void foreach(Computer &computer) = 0;
      ///
      virtual void foreach(ConstComputer &computer) const = 0;

      Storage* getStorage() { return theStorage; }
    };
  }
}

#endif
