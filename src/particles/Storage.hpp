#ifndef _PARTICLES_STORAGE_HPP
#define _PARTICLES_STORAGE_HPP

#include "types.hpp"
#include "esutil/TupleVector.hpp"
#include "Particle.hpp"
#include "Property.hpp"

namespace espresso {
  namespace particles {
    typedef esutil::real real;

    // forward declaration of Computer and ConstComputer
    class Computer;
    class ConstComputer;

    /**
       MOCK implementation of particlestorage.
       Currently only provides a subset of the expected interface,
       and without guarantee :-). For the intended use, refer to
       @ref Storage.cpp, which contains a routine intended_use,
       which shows how the class is intended to be used.
    */
    class Storage {
    public:

      Storage();

      virtual ~Storage() {}

      /// @name access to particles
      //@{

      /** add a particle
          @return a temporary reference to the created particle. This reference only identifies
          the particle as long as the Storage is not modified. For keeping a persistent
          reference to the particle, use @ref getParticleID.
      */
      virtual ParticleReference addParticle();

      /** delete a particle
          @param id the id of the particle to delete
      */
      virtual void deleteParticle(ParticleId id);

      /** get a persistent particle from a temporary reference
          @param ref temporary reference to a particle
          @return the particle
      */
      ParticleId getParticleId(ConstParticleReference ref) const {
        return getIDProperty()[ref];
      }

      /** get a reference to a particle using its ID. This is
          potentially slow.
          @param id the id of the particle to fetch
      */
      virtual ParticleReference getParticleReference(ParticleId id);

      /** get a particle to a particle using its ID. This is
          potentially slow.
          @param id the id of the particle to fetch
      */
      virtual ConstParticleReference getParticleReference(ParticleId id) const;

      /// loop over all particles
      virtual void foreach(Computer &);
      /// loop over all particles
      virtual void foreach(ConstComputer &) const;

      //@}

      /// @name access to particle properties
      //@{

      /** get a short lifetime reference to a property by its ID
          @throw std::out_of_range if one tries to obtain a handle to the ID property
      */
      template<typename T>
      PropertyReference<T> getPropertyReference(PropertyId id) {
        // no non-const reference to the ID
        if (id == particleIDProperty) {
          throw std::out_of_range("id is not writable");
        }
        return particles.getProperty<T>(id);
      }
      /// get a short lifetime reference to a property by its ID
      template<typename T>
      ConstPropertyReference<T> getPropertyReference(PropertyId id) const {
        return particles.getProperty<T>(id);
      }

      /** get a short lifetime reference to a property by its ID
          @throw std::out_of_range if one tries to obtain a handle to the ID property
      */
      template<typename T>
      ArrayPropertyReference<T> getArrayPropertyReference(PropertyId id) {
        // no non-const reference to the ID
        if (id == particleIDProperty) {
          throw std::out_of_range("id is not writable");
        }
        return particles.getArrayProperty<T>(id);
      }
      /** get a short lifetime reference to a property by its ID
          @throw std::range_error if the given and array dimensions mismatch
      */
      template<typename T>
      ConstArrayPropertyReference<T> getArrayPropertyReference(PropertyId id) const {
        return particles.getArrayProperty<T>(id);
      }

      /// get a short lifetime reference to the property representing the particle ID
      const ConstPropertyReference<ParticleId> getIDProperty() const {
        return particles.getProperty<ParticleId>(particleIDProperty);
      }

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
      void eraseProperty(PropertyId id);

      //@}

      /// @name convenience functions
      //@{

      /** fill the storage with a regular 3D lattice of particles
          @param size the size of the box to fill
          @param N the number of particle per side
          @param positions where to write the particle positions. If omitted,
          a new property is created
          @return the ID of the position property
      */
      PropertyId fillWithLattice(real size, size_t N, PropertyId positions);

      //@}

      /** make this class available at Python. */

      static void registerPython();

    protected:
      /// here the particle data is stored
      esutil::TupleVector particles;

    private:
      /// unique ID counter for the particles
      size_t uniqueID;
      /// ID of the particle ID property
      PropertyId particleIDProperty;

      /// private and does not exist, do not try to use
      Storage(const Storage &);
    };

  }
}

#endif
