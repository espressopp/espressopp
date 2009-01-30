#ifndef _PARTICLESTORAGE_PARTICLESTORAGE_HPP
#define _PARTICLESTORAGE_PARTICLESTORAGE_HPP

// only included for the mock implementation. Remove for final
#include <vector>

#include "types.hpp"
#include "util/virtual_functional.hpp"
#include "util/TupleVector.hpp"

namespace espresso {
    namespace particlestorage {

	/**
	   MOCK implementation of particlestorage.
	   Currently only provides a subset of the expected interface,
	   and without guarantee :-). For the intended use, refer to
	   @ref ParticleStorage.cpp, which contains a routine intended_use,
	   which shows how the class is intended to be used.
	*/
	class ParticleStorage {
	public:
	    typedef util::TupleVector::reference reference;
	    typedef util::TupleVector::const_reference const_reference;

	    template<typename T>
	    struct PropertyTraits {
		typedef util::TupleVector::PropertyReference<T> Reference;
		typedef util::TupleVector::ConstPropertyReference<T> ConstReference;
		typedef util::TupleVector::ArrayPropertyReference<T> ArrayReference;
		typedef util::TupleVector::ConstArrayPropertyReference<T> ConstArrayReference;
	    };

	private:
	    util::TupleVector particles;

	    /// unique ID counter for the particles
	    size_t uniqueID;
	    /// ID of the particle ID property
	    size_t particleIDProperty;
	public:

	    ParticleStorage();

	    virtual ~ParticleStorage() {}

	    /// @name access to particles
	    //@{

	    /** add a particle
		@return the ID of the added particle, basically just for deleting it again...
	     */
	    virtual size_t addParticle();

	    /** delete a particle
		@param id the id of the particle to delete
	    */
	    virtual void deleteParticle(size_t id);

	    /** get a reference to a particle using its ID. This is
		potentially slow.
		@param id the id of the particle to fetch
	     */
	    virtual reference getParticleByID(size_t id);

	    /** get a particle to a particle using its ID. This is
		potentially slow.
		@param id the id of the particle to fetch
	     */
	    virtual const_reference getParticleByID(size_t id) const {
		return const_reference(const_cast<ParticleStorage *>(this)->getParticleByID(id));
	    }

	    /// loop over all particles
	    virtual void foreach(util::VirtualUnaryFunction<reference, void> &);
	    /// loop over all particles
	    virtual void foreach(util::VirtualUnaryFunction<const_reference, void> &) const;

	    //@}

	    /// @name access to particle properties
	    //@{

	    /// get a short lifetime reference to a property by its ID
	    template<typename T>
	    typename PropertyTraits<T>::Reference getProperty(size_t id) {
		// no non-const reference to the ID
		if (id == particleIDProperty) {
		    throw std::out_of_range("id is not writable");
		}
		return particles.getProperty<T>(id);
	    }
	    /// get a short lifetime reference to a property by its ID
	    template<typename T>
	    typename PropertyTraits<T>::ConstReference getProperty(size_t id) const {
		return particles.getProperty<T>(id);
	    }
	    /// get a short lifetime reference to a property by its ID
	    template<typename T>
	    typename PropertyTraits<T>::ArrayReference getArrayProperty(size_t id) {
		// no non-const reference to the ID
		if (id == particleIDProperty) {
		    throw std::out_of_range("id is not writable");
		}
		return particles.getArrayProperty<T>(id);
	    }
	    /// get a short lifetime reference to a property by its ID
	    template<typename T>
	    typename PropertyTraits<T>::ConstArrayReference getArrayProperty(size_t id) const {
		return particles.getArrayProperty<T>(id);
	    }

	    /// get a short lifetime reference to the property representing the particle ID
	    const PropertyTraits<size_t>::ConstReference getIDProperty() const {
		return particles.getProperty<size_t>(particleIDProperty);
	    }

	    /** add a property
		@param dim dimensionality of the property (e.g. 3 for a 3D vector)
		@return the ID of the property for use with getProperty or eraseProperty
		@tparam T type of the property
	    */
	    template<typename T>
	    size_t addProperty(size_t dim = 1) { return particles.addProperty<T>(dim); }

	    /** delete a property
		@param n ID of the property to delete as obtained from addProperty
	    */
	    void eraseProperty(size_t id) {
		// no non-const reference to the ID
		if (id == particleIDProperty) {
		    throw std::out_of_range("id cannot be erased");
		}
		particles.eraseProperty(id);
	    }

	    //@}
	};

	typedef util::VirtualUnaryFunction<ParticleStorage::reference, void> ParticleComputer;
	typedef util::VirtualUnaryFunction<ParticleStorage::const_reference, void> ConstParticleComputer;
    }
}

#endif
