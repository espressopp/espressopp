#ifndef _STORAGE_STORAGE_HPP
#define _STORAGE_STORAGE_HPP

#include <exception>

#include "types.hpp"
#include "esutil/TupleVector.hpp"
#include "Particle.hpp"
#include "storage/ParticleHandle.hpp"
#include "storage/PropertyHandle.hpp"
#include "particles/Set.hpp"
#include "pairs/Set.hpp"
#include "bc/BC.hpp"

namespace espresso {
  // forward declarations
  class PropertyBase;
  template<typename> class Property;
  template<typename> class ArrayProperty;

  namespace storage {
    typedef PropertyHandle< ParticleId > IdPropertyHandle;
    typedef PropertyHandle< Real3D > PosPropertyHandle;
    typedef esutil::TupleVector::PropertyId PropertyId;

    class Storage : public particles::Set, 
		    public pairs::Set, 
		    public enable_shared_from_this< Storage >,
		    boost::noncopyable
    {
    public:
      typedef shared_ptr< Storage > SelfPtr;

      /** constructor, specifies the boundary conditions to apply.
	  Derived Storage classes restrict the possible boundary
	  condition classes according to their needs. */
      Storage(bc::BC::SelfPtr);

      virtual ~Storage();

      /// @name access to particles
      //@{

      /** add a particle with a given id. Note that this is a basic
	  and slow implementation, but derived classes can and should
	  implement this more efficiently.
          @param id the id of the particle to create.
	  @throw std::out_of_range if the particle already exists
      */
      virtual ParticleHandle addParticle(ParticleId id) = 0;
      /** This is mostly identical to addParticle(). It is only
	  required to allow exporting the function to Python, as Python
	  does know how to handle the return value of addParticle(). */
      void _addParticle(ParticleId id);

      /** delete a particle
          @param id the id of the particle to delete
	  @throw std::out_of_range if the particle does not exist
      */
      virtual void deleteParticle(ParticleId id) = 0;

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
      virtual ParticleHandle getParticleHandle(ParticleId id) = 0;

      /// get a short lifetime reference to the property representing the particle ID
      IdPropertyHandle getIdPropertyHandle();
      /// get a short lifetime reference to the property representing the particle positions
      PosPropertyHandle getPositionPropertyHandle();

      //@}

      bc::BC::SelfPtr getBoundaryConditions() const { return bc; }

      /// @name the particles::Set interface
      //@{
      
      /** for a particle of the ParticleStorage of this class,
	  check whether it belongs to this set
      */
      bool contains(ParticleHandle);
      bool contains(ParticleId);

      SelfPtr getStorage();
      //@}

      /// @name the pairs::Set interface
      //@{

      SelfPtr getLeftStorage();
      SelfPtr getRightStorage();

      //@}

      /** call this to specify the property representing the particle
	  ids.  You cannot call this function twice.  This has to be
	  done in Python so that we get a proper wrapper.

	  If no property is specified, this constructs the property internally.
	  Use this in simply C++-code, where access to the handles is sufficient.
      */
      virtual void setIdProperty(boost::shared_ptr< Property< ParticleId > > =
				 boost::shared_ptr< Property< ParticleId > >());

      /** call this to specify the property representing the particle
	  position.  This has to be done in Python so that we get a
	  proper wrapper.

	  If no property is specified, this constructs the property internally.
	  Use this in simply C++-code, where access to the handles is sufficient.
      */
      virtual void setPositionProperty(boost::shared_ptr< Property< Real3D > > =
				       boost::shared_ptr< Property< Real3D > >());

      /// make this class available at Python
      static void registerPython();

    protected:
      bool foreachApply(particles::Computer &computer) = 0;
      /// the default implementation loops over all pairs in the container
      bool foreachPairApply(pairs::Computer &computer);

      /** get reference to the underlying TupleVector. The Storage
	  uses this solely for adding/removing properties. */
      virtual esutil::TupleVector &getTupleVector() = 0;

      /// boundary condition in charge
      bc::BC::SelfPtr bc;

    private:
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
      PropertyId addProperty(size_t dim = 1) { return getTupleVector().addProperty<T>(dim); }

      /** delete a property
          @param n ID of the property to delete as obtained from addProperty
      */
      void deleteProperty(PropertyId id);

      /// get a short lifetime reference to a property by its ID
      template< typename T >
      PropertyHandle< T > getPropertyHandle(PropertyId id) {
        return getTupleVector().getProperty< T >(id);
      }

      /// get a short lifetime reference to a property by its ID
      template< typename T >
      ArrayPropertyHandle< T > getArrayPropertyHandle(PropertyId id) {
        return getTupleVector().getArrayProperty< T >(id);
      }

      //@}

      /// ID of the particle ID property
      PropertyId particleIdProperty;
      /// ID of the particle ID property
      PropertyId particlePosProperty;
    };

  }
}

#endif
