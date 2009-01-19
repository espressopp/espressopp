#ifndef _PARTICLESTORAGE_PARTICLESTORAGE_HPP
#define _PARTICLESTORAGE_PARTICLESTORAGE_HPP

// only included for the mock implementation. Remove for final
#include <vector>

#include "types.hpp"

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
	    /// just to have some particle data
	    std::vector<size_t> id;
	    /// just to have some particle data
	    std::vector<real>   position;
	    /// just to have some particle data
	    std::vector<real>   force;

	public:
	    /** a reference to an element, which then can be used to read and
		write properties of the element using @ref ParticleStorage::PropertyReference.
	
		Actually, this is nothing but the index of the element. A
		reference can only be obtained from element access to the
		ParticleSet/Storage classes.
	    */
	    class ReferenceBase {
	    protected:
		/// index
		size_t index;

		/// constructor
		ReferenceBase(size_t _index): index(_index) {}
	    public:
		/// convert to index
		operator size_t() const { return index; }
	    };

	    /** see @ref ParticleStorage::ReferenceBase; this class marks the referenced
		element as non-const. Later, this will be fetched from the MultiVector,
		and there not require several friends
	    */
	    class reference: public ReferenceBase {
		friend class ParticleStorage;
		friend class ParticleSet;
		/// constructable only through @ref ParticleStorage/Set
		reference(size_t _index): ReferenceBase(_index) {}
	    };

	    /** see @ref MultiVector::ReferenceBase; this class marks the referenced
		element as const */
	    class const_reference: public ReferenceBase {
		friend class ParticleStorage;
		friend class ParticleSet;
		/// constructable only through @ref MultiVector
		const_reference(size_t _index): ReferenceBase(_index) {}
	    public:
		/// non-const->const conversion
		const_reference(const reference &_ref): ReferenceBase(_ref.index) {}
	    };

	    /** reference to a scalar property. This is a base implementation, from
		which the const and non-const versions inherit.
	    */
	    template<class T, class VectorClass>
	    class PropertyReferenceBase {
	    protected:
		VectorClass &vec;
		/// constructor
		PropertyReferenceBase(VectorClass &_vec): vec(_vec) {}
	    public:
		/// dereference
		T &operator[](reference n) { return vec[n]; }
		/// dereference
		const T &operator[](const_reference n) const { return vec[n]; }
	    };

	    /// see @ref ParticleStorage::PropertyReferenceBase
	    template<class T>
	    class PropertyReference: public PropertyReferenceBase< T, std::vector<T> > {
		friend class ParticleStorage;
		PropertyReference(std::vector<T> &_vec)
		    : PropertyReferenceBase< T, std::vector<T> >(_vec) {}
	    };
	    
	    /// see @ref ParticleStorage::PropertyReferenceBase; this is for a constant property
	    template<class T>
	    class ConstPropertyReference: public PropertyReferenceBase<const T, const std::vector<T> > {
		friend class ParticleStorage;
		ConstPropertyReference(const std::vector<T> &_vec)
		    : PropertyReferenceBase< const T, const std::vector<T> >(_vec) {}
	    public:
		/// non-const->const conversion
		ConstPropertyReference(const PropertyReference<T> &_ref)
		    : PropertyReferenceBase<const T, const std::vector<T> >(_ref.vec) {}
	    };

	    /** reference to an array property. In contrast to
		@ref ParticleStorage::PropertyReference, this returns a pointer, since
		the element itself is in fact an array whose size
		is only known at runtime.
	    */
	    template<class T, class VectorClass>
	    class ArrayPropertyReferenceBase {
	    protected:
		VectorClass &vec;
		/// constructor
		ArrayPropertyReferenceBase(VectorClass &_vec): vec(_vec) {}

	    public:
		/// dereference
		T *operator[](reference n) { return &vec[3*n]; }
		/// dereference
		const T *operator[](const_reference n) const { return &vec[3*n]; }
	    };

	    /** see @ref ParticleStorage::ArrayPropertyReferenceBase
		@tparam T the data type of the property (or more precisely,
		the datatype to interpret the property as)
	    */
	    template<class T>
	    class ArrayPropertyReference: public ArrayPropertyReferenceBase< T, std::vector<T> > {
		friend class ParticleStorage;

		ArrayPropertyReference(std::vector<T> &_vec)
		    : ArrayPropertyReferenceBase< T, std::vector<T> >(_vec) {};
	    };

	    /** see @ref ParticleStorage::ArrayPropertyReferenceBase;
		this is for a constant property
		@tparam T the data type of the property (or more precisely,
		the datatype to interpret the property as)
	    */
	    template<class T>
	    class ConstArrayPropertyReference
		: public ArrayPropertyReferenceBase< const T, const std::vector<T> > {
		friend class ParticleStorage;
		
		ConstArrayPropertyReference(const std::vector<T> &_vec)
		    : ArrayPropertyReferenceBase< const T, const std::vector<T> >(_vec) {};
	    public:
		/// non-const->const conversion
		ConstArrayPropertyReference(const ArrayPropertyReference<T> &_ref)
		    : ArrayPropertyReferenceBase<const T, const std::vector<T> >(_ref.vec) {}
	    };

	    virtual ~ParticleStorage() {}

	    /** get reference to a particle according to its position. Later, this will not be possible.
		Particles may be obtained through their global id or through ParticleSet.
	    */
	    reference getParticle(size_t id) { return reference(id); }
	    /** get reference to a particle according to its position. Later, this will not be possible.
		Particles may be obtained through their global id or through ParticleSet.
	    */
	    const_reference getParticle(size_t id) const { return const_reference(id); }

	    /// @name mockaccess to particle properties
	    //@{

	    /// this will be later replaced by a generic interface to obtain a property
	    PropertyReference<size_t> getIDProperty() { return PropertyReference<size_t>(id); }
	    /// this will be later replaced by a generic interface to obtain a property
	    ConstPropertyReference<size_t> getIDProperty() const { return ConstPropertyReference<size_t>(id); }
	    
	    /// this will be later replaced by a generic interface to obtain a property
	    ArrayPropertyReference<real> getPosProperty() { return ArrayPropertyReference<real>(position); }
	    /// this will be later replaced by a generic interface to obtain a property
	    ConstArrayPropertyReference<real> getPosProperty() const { return ConstArrayPropertyReference<real>(position); }

	    /// this will be later replaced by a generic interface to obtain a property
	    ArrayPropertyReference<real> getForceProperty() { return ArrayPropertyReference<real>(force); }
	    /// this will be later replaced by a generic interface to obtain a property
	    ConstArrayPropertyReference<real> getForceProperty() const { return ConstArrayPropertyReference<real>(force); }
	    
	    //@}

	    /// loop
	    virtual void foreach(class ParticleComputer &);
	    virtual void foreach(const class ParticleComputer &) const;
	};
    }
}

#endif
