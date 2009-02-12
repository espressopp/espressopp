#ifndef _PARTICLESTORAGE_PARTICLESTORAGE_HPP
#define _PARTICLESTORAGE_PARTICLESTORAGE_HPP

// only included for the mock implementation. Remove for final
#include <vector>

#include "types.hpp"
#include "esutil/virtual_functional.hpp"
#include "esutil/TupleVector.hpp"

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
	    typedef esutil::TupleVector::reference reference;
	    typedef esutil::TupleVector::const_reference const_reference;

	    template<typename T>
	    struct PropertyTraits {
		typedef esutil::TupleVector::PropertyReference<T> Reference;
		typedef esutil::TupleVector::ConstPropertyReference<T> ConstReference;
            };

            template<typename T>
	    struct ArrayPropertyTraits {
                typedef esutil::TupleVector::ArrayPropertyReference<T> Reference;
                typedef esutil::TupleVector::ConstArrayPropertyReference<T> ConstReference;
            };

	private:
	    esutil::TupleVector particles;

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
		@return a temporary reference to the created particle. This reference only identifies
		the particle as long as the ParticleStorage is not modified. For keeping a persistent
		reference to the particle, use @ref getParticleID.
	     */
	    virtual reference addParticle();

	    /** delete a particle
		@param id the id of the particle to delete
	    */
	    virtual void deleteParticle(size_t id);

	    /** get the persistent ID of a particle
		@param ref temporary reference to a particle
		@return the persistent ID of the particle
	     */
	    size_t getParticleID(const_reference ref) const {
		return getIDProperty()[ref];
	    }

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
	    virtual void foreach(esutil::VirtualUnaryFunction<reference, void> &);
	    /// loop over all particles
	    virtual void foreach(esutil::VirtualUnaryFunction<const_reference, void> &) const;

	    //@}

	    /// @name access to particle properties
	    //@{

	    /** get a short lifetime reference to a property by its ID
                @throw std::out_of_range if one tries to obtain a handle to the ID property
            */
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

	    /** get a short lifetime reference to a property by its ID
                @throw std::out_of_range if one tries to obtain a handle to the ID property
            */
	    template<typename T>
	    typename ArrayPropertyTraits<T>::Reference getVarArrayProperty(size_t id) {
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
	    typename ArrayPropertyTraits<T>::ConstReference getArrayProperty(size_t id) const {
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

            /// @name convenience functions
            //@{

            /** fill the storage with a regular 3D lattice of particles
                @param size the size of the box to fill
                @param N the number of particle per side
                @param positions where to write the particle positions. If omitted,
                a new property is created
                @return the ID of the position property
             */
            size_t fillWithLattice(real size, size_t N, size_t positions);

            //@}

        private:
            /// private and does not exist, do not try to use
            ParticleStorage(const ParticleStorage &);
	};

	typedef esutil::VirtualUnaryFunction<ParticleStorage::reference, void> ParticleComputer;
	typedef esutil::VirtualUnaryFunction<ParticleStorage::const_reference, void> ConstParticleComputer;
    }
}

#endif
