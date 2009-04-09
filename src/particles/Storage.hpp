#ifndef _PARTICLES_STORAGE_HPP
#define _PARTICLES_STORAGE_HPP

#include "types.hpp"
#include "esutil/TupleVector.hpp"
#include "Particle.hpp"
#include "ParticleHandle.hpp"
#include "PropertyHandle.hpp"

namespace espresso {
  // forward declarations
  template<typename> class Property;
  template<typename> class ArrayProperty;

  namespace particles {
    // forward declarations
    class Computer;
    class ConstComputer;

    /**
       MOCK implementation of particlestorage.
       Currently only provides a subset of the expected interface,
       and without guarantee :-).
    */
    class Storage {
    public:

      Storage();

      virtual ~Storage() {}

      /// @name access to particles
      //@{

      /** add a particle.
          @return the id of the created particle.
      */
      virtual ParticleId addParticle();

      /** delete a particle
          @param id the id of the particle to delete
      */
      virtual void deleteParticle(ParticleId id);

      /** get a persistent particle from a temporary reference
          @param ref temporary reference to a particle
          @return the particle
      */
      ParticleId getParticleId(ConstParticleHandle ref) const {
        return getIdPropertyHandle()[ref];
      }

      /** get a reference to a particle using its ID. This is
          potentially slow.
          @param id the id of the particle to fetch
      */
      virtual ParticleHandle getParticleHandle(ParticleId id);

      /** get a particle to a particle using its ID. This is
          potentially slow.
          @param id the id of the particle to fetch
      */
      virtual ConstParticleHandle getParticleHandle(ParticleId id) const;

      /// loop over all particles
      virtual void foreach(Computer &);
      /// loop over all particles
      virtual void foreach(ConstComputer &) const;

      //@}

      /// get a short lifetime reference to the property representing the particle ID
      const ConstPropertyHandle<ParticleId> getIdPropertyHandle() const {
        return particles.getProperty<ParticleId>(particleIdProperty);
      }

      /// make this class available at Python
      static void registerPython();

    protected:
      /// here the particle data is stored
      esutil::TupleVector particles;

    private:
      typedef esutil::TupleVector::PropertyId PropertyId;

      /// since Property/ArrayProperty are the ones that operates properties
      template<typename> friend class espresso::Property;
      template<typename> friend class espresso::ArrayProperty;

      /// @name access to particle properties. Used by the Property classes
      //@{

      /** add a property
          @param dim dimensionality of the property (e.g. 3 for a 3D vector)
          @return the ID of the property for use with getProperty or eraseProperty
          @tparam T type of the property
      */
      template<typename T>
      PropertyId addProperty(size_t dim = 1) { return particles.addProperty<T>(dim); }

      /** delete a property
          @param n ID of the property to delete as obtained from addProperty
      */
      void deleteProperty(PropertyId id);

      /** get a short lifetime reference to a property by its ID
          @throw std::out_of_range if one tries to obtain a handle to the ID property
      */
      template<typename T>
      PropertyHandle<T> getPropertyHandle(PropertyId id) {
        // no non-const reference to the ID
        if (id == particleIdProperty) {
          throw std::out_of_range("id is not writable");
        }
        return particles.getProperty<T>(id);
      }
      /// get a short lifetime reference to a property by its ID
      template<typename T>
      ConstPropertyHandle<T> getConstPropertyHandle(PropertyId id) const {
        return particles.getProperty<T>(id);
      }

      /** get a short lifetime reference to a property by its ID
          @throw std::out_of_range if one tries to obtain a handle to the ID property
      */
      template<typename T>
      ArrayPropertyHandle<T> getArrayPropertyHandle(PropertyId id) {
        // no non-const reference to the ID
        if (id == particleIdProperty) {
          throw std::out_of_range("id is not writable");
        }
        return particles.getArrayProperty<T>(id);
      }
      /** get a short lifetime reference to a property by its ID
          @throw std::range_error if the given and array dimensions mismatch
      */
      template<typename T>
      ConstArrayPropertyHandle<T> getConstArrayPropertyHandle(PropertyId id) const {
        return particles.getArrayProperty<T>(id);
      }

      //@}

      /// unique ID counter for the particles
      size_t uniqueID;
      /// ID of the particle ID property
      PropertyId particleIdProperty;

      /// private and does not exist, do not try to use
      Storage(const Storage &);
    };

  }
}

#endif
