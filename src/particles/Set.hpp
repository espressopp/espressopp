#ifndef _PARTICLES_SET_HPP
#define _PARTICLES_SET_HPP

#include "types.hpp"
#include "Storage.hpp"
#include "Computer.hpp"

namespace espresso {
  namespace particles {
    /** MOCK particle set. Provides a view onto a set of particles
        from a ParticleStorage
    */
    class Set {
    protected:
      /// the storage our particles are stored in
      Storage::SelfPtr theStorage;

    public:
      typedef shared_ptr<Set> SelfPtr;

      /** base constructor
          
          @param _store pointer to the Storage the
          particles in this set come from
      */
      Set(Storage::SelfPtr _store = Storage::SelfPtr()) 
	: theStorage(_store) {}

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

      Storage::SelfPtr getStorage() const { return theStorage; }

    public:

      /** Abstract class needs also registration in Python */

      static void registerPython();

    };
  }
}

#endif
