#ifndef _PARTICLES_STORAGE_HPP
#define _PARTICLES_STORAGE_HPP

#include <exception>

#include "types.hpp"
#include "esutil/TupleVector.hpp"
#include "Particle.hpp"
#include "particles/ParticleHandle.hpp"
#include "particles/PropertyHandle.hpp"
#include "particles/Set.hpp"

namespace espresso {
  // forward declarations
  class PropertyBase;
  template<typename> class Property;
  template<typename> class ArrayProperty;

  namespace particles {
    typedef PropertyHandle< ParticleId > IdPropertyHandle;

    // forward declaration
    class StorageMismatch : public std::exception 
    {};

    class Storage : public Set, 
		    public enable_shared_from_this< Storage >,
		    boost::noncopyable
    {
    public:
      typedef shared_ptr< Storage > SelfPtr;

      Storage();

      virtual ~Storage();

      /// @name access to particles
      //@{

      /** add a particle with a given id. Note that this is a basic
	  and slow implementation, but derived classes can and should
	  implement this more efficiently.
          @param id the id of the particle to create.
	  @throw std::out_of_range if the particle already exists
      */
      virtual ParticleHandle addParticle(ParticleId id);
      void _addParticle(ParticleId id);

      /** delete a particle
          @param id the id of the particle to delete
	  @throw std::out_of_range if the particle does not exist
      */
      virtual void deleteParticle(ParticleId id);

      /** get a persistent particle id from a temporary reference
          @param ref temporary reference to a particle
          @return the particle
      */
      ParticleId getParticleId(ParticleHandle ref) {
        return getIdPropertyHandle()[ref];
      }

      /** get a reference to a particle using its ID. This is
          a slow implementation, but derived classes can and should
	  implement this more efficiently.
          @param id the id of the particle to fetch
	  @return a handle to the particle, or ParticleHandle() if
	  the particle does not exist
      */
      virtual ParticleHandle getParticleHandle(ParticleId id);


      //@}

      /// get a short lifetime reference to the property representing the particle ID
      IdPropertyHandle getIdPropertyHandle() {
        return particles.getProperty< ParticleId >(particleIdProperty);
      }

      // The Set interface
      
      /** for a particle of the ParticleStorage of this class,
	  check whether it belongs to this set
      */
      virtual bool contains(ParticleHandle);
      virtual bool contains(ParticleId);


      /** apply computer to all particles of this set
       */
      virtual void foreach(ApplyFunction function);
      // make the other variants of foreach available
      using Set::foreach;

      virtual SelfPtr getStorage();

      void checkProperty(shared_ptr< PropertyBase > prop);

      /// make this class available at Python
      static void registerPython();

    protected:

      /// here the particle data is stored
      esutil::TupleVector particles;

    private:
      typedef esutil::TupleVector::PropertyId PropertyId;

      /// since Property/ArrayProperty are the ones that operates properties
      template< typename > friend class espresso::Property;
      template< typename > friend class espresso::ArrayProperty;

      /// @name access to particle properties. Used by the Property classes
      //@{

      /** add a property
          @param dim dimensionality of the property (e.g. 3 for a 3D vector)
          @return the ID of the property for use with getProperty or eraseProperty
          @tparam T type of the property
      */
      template< typename T >
      PropertyId addProperty(size_t dim = 1) { return particles.addProperty<T>(dim); }

      /** delete a property
          @param n ID of the property to delete as obtained from addProperty
      */
      void deleteProperty(PropertyId id);

      /// get a short lifetime reference to a property by its ID
      template<typename T>
      PropertyHandle<T> getPropertyHandle(PropertyId id) {
        return particles.getProperty<T>(id);
      }

      /// get a short lifetime reference to a property by its ID
      template<typename T>
      ArrayPropertyHandle<T> getArrayPropertyHandle(PropertyId id) {
        return particles.getArrayProperty<T>(id);
      }

      //@}

      /// ID of the particle ID property
      PropertyId particleIdProperty;
    };


  }
}

#endif
