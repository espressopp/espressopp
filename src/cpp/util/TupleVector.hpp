#ifndef _TUPLEVECTOR_HPP
#define _TUPLEVECTOR_HPP

#include <vector>
#include <stdexcept>
#include <boost/iterator/iterator_facade.hpp>

namespace util {

    /** a vector of elements with dynamically configurable properties.
	Each property is stored in a vector of its own, making loops
	over single properties very fast. However, as a consequence,
	there is no value type, and some vector operations are slightly
	different (e.g., elements are uninitialized)
    */
    class TupleVector {

	/** Handler for a property. This really just stores the data;
	    (re)allocation and freeing needs to be done from outside */
	struct Property {
	    /// size of one element in characters
	    size_t size;
	    /// number of data type values in one element
	    size_t dimension;
	    /// memory segment storing the elements
	    void *data;

	    /** construct an empty property
		@param _size byte size of an element
		@param _dimension number of values in an element
		(i.e. _size/_dimension = sizeof(datatype))
		@param _data memory segment address
	    */
	    Property(size_t _size, size_t _dimension, void *_data)
		: size(_size), dimension(_dimension), data(_data) {}
	};

	/// handlers for all properties
	std::vector <Property> property;

	/// size in elements of each of the stored properties
	size_t eSize;

	/// reserved size in elements of each of the stored properties
	size_t maxESize;

    public:
	/// @name standard vector members
	//@{

	///
	typedef size_t size_type;

	///
	typedef ptrdiff_t difference_type;
  
	/** a reference to an element, which then can be used to read and
	    write properties of the element using TupleVector::PropertyReference.
	
	    Actually, this is nothing but the index of the element. A
	    reference can only be obtained from element access to the
	    TupleVector class or its iterators.
	*/
	class ReferenceBase {
	protected:
	    /// index of the reference particle
	    size_type index;

	    /// constructor
	    ReferenceBase(size_type _index): index(_index) {}
	public:
	    /// convert to index
	    operator size_type() const { return index; }
	};

	/** see @ref TupleVector::ReferenceBase; this class marks the referenced
	    element as non-const */
	class reference: public ReferenceBase {
	    friend class TupleVector;
	    /// constructable only through @ref TupleVector
	    reference(size_type _index): ReferenceBase(_index) {}
	};

	/** see @ref TupleVector::ReferenceBase; this class marks the referenced
	    element as const */
	class const_reference: public ReferenceBase {
	    friend class TupleVector;
	    /// constructable only through @ref TupleVector
	    const_reference(size_type _index): ReferenceBase(_index) {}
	public:
	    /// non-const->const conversion
	    const_reference(const reference &_ref): ReferenceBase(_ref.index) {}
	};

	/** an iterator over elements. This is basically a reference to
	    an element with a mutable index

	    @tparam ReferenceType
	*/
	template<class ReferenceType, class CRTP>
	class IteratorBase: public boost::iterator_facade<
	    CRTP, 
	    ReferenceType,
	    boost::random_access_traversal_tag,
	    ReferenceType,
	    difference_type>
	{
	protected:
	    /// index of particle that is referenced
	    size_type index;

	    /// regular constructor
	    IteratorBase(size_type _index): index(_index) {}
	public:
	    /// reference type. Fixes up the iterator_facade for a reference only class
	    typedef ReferenceType reference;

	    /// default constructor, undefined iterator
	    IteratorBase() {}

	private:
	    friend class boost::iterator_core_access;

	    bool equal(const IteratorBase<ReferenceType, CRTP> &other) const {
		return index == other.index;
	    }

	    void increment()          { index++; }
	    void decrement()          { index--; }
	    void advance(size_type n) { index += n; }
	    reference dereference() const { return reference(index); }
	    difference_type distance_to(const IteratorBase<ReferenceType, CRTP> &other) const{
		return other.index - index;
	    }
	};

	/// see @ref TupleVector::IteratorBase
	class iterator: public IteratorBase<reference, iterator> {
	    friend class TupleVector;
	    /// constructable only through @ref TupleVector
	    iterator(size_type _index): IteratorBase<reference, iterator>(_index) {}
	public:
	    /// default constructor, uninitialized iterator
	    iterator() {}
	};
	/// see @ref TupleVector::IteratorBase
	class const_iterator: public IteratorBase<const_reference, const_iterator> {
	    friend class TupleVector;
	    /// constructable only through @ref TupleVector
	    const_iterator(size_type _index): IteratorBase<const_reference, const_iterator>(_index) {}
	public:
	    /// default constructor
	    const_iterator() {}
	    /// non-const->const conversion
	    const_iterator(const iterator &_it)
		: IteratorBase<const_reference, const_iterator>(_it.index) {}
	};

	//@}

	/** default constructor. Constructs an empty TupleVector, but possibly
	    allocating already some particles. When properties are added, they
	    will offer space for the specified number of particles.
	*/
	explicit TupleVector(size_type _n = 0)
	    : eSize(0), maxESize(0) { resize(_n); };

	/** partial copy constructor. Creates a new TupleVector with the same
	    properties as the previous one, but contains only n uninitialized
	    elements.
	    @param vec vector to copy
	    @param n number of elements to reserve
	*/
	explicit TupleVector(const TupleVector &vec, size_type n);
	///
	~TupleVector();

	/// @name standard vector members
	//@{

	///
	iterator begin() { return iterator(0); }
	///
	const_iterator begin() const { return const_iterator(0); }

	///
	iterator end()   { return iterator(eSize); }
	///
	const_iterator end() const { return const_iterator(eSize); };

	/// get reference to nth element
	reference operator[](size_type n) { return reference(n); }
	/// get const reference to nth element
	const_reference operator[](size_type n) const {
	    return const_reference(n);
	}

	/// get checked reference to nth element
	reference at(size_type n)  {
	    if (n >= eSize) {
		throw std::out_of_range("TupleVector::at");
	    }
	    return reference(n);
	}
	/// get checked const reference to nth element
	const_reference at(size_type n) const {
	    if (n >= eSize) {
		throw std::out_of_range("TupleVector::at");
	    }
	    return const_reference(n);
	}

	/// insert an uninitialized element at position pos
	iterator insert(iterator pos);
	/// insert a copy of element e at position pos
	iterator insert(iterator pos, reference e)
	    { return insert(pos, const_reference(e)); }
	/// insert a copy of element e at position pos
	iterator insert(iterator pos, const_reference e);
	/// insert n uninitialized elements at position pos
	void insert(iterator pos, size_type n);

	/// delete element at position pos
	iterator erase(iterator pos);
	/// delete elements from sequence
	iterator erase(iterator start, iterator end);

	/// remove all elements, but keep set of properties
	void clear() { resize(0); }

	/// get number of elements
	size_type size() const { return eSize; }
	/// @return true iff the TupleVector contains elements
	bool empty() const { return eSize == 0; }
	/// @return the maximally possible size
	size_type max_size() const { return size_t(-1); }
	/// @return the current capacity, see STL vector documentation
	size_type capacity() const { return maxESize; }

	/** change the number of elements
	    @param newESize new number of elements
	*/
	void resize(size_type newESize) {
	    eSize = newESize;
	    reserve(newESize);
	}

	/** reserve space. This function can be used to preallocate memory.
	    For any resize up th minCapacity, no further reallocation takes
	    place. In contrast to STL implementations, this reserve function
	    also might reduce the currrent capacity, if the requested capacity
	    is much smaller than than the current one. Reserving a capacity
	    smaller than the current size of the vector has no effect; the
	    capacity is always larger than the size of the TupleVector.

	    @param minCapacity number of elements to reserve at least
	*/
	void reserve(size_type minCapacity);

	//@}

	/** @name full element related members
	    here are additional routines which allow to work
	    with a whole element. These functions are in a
	    way standard vector elements, however, for a standard
	    vector, these functions are provided by the reference
	    member only. Since in this class, the reference does
	    not hold the full information, the functions need to be
	    provided by the container.
	*/
	//@{
	/** copy all data of a particle
	    @param dst where to write to
	    @param src where to get the data from
	*/
	void assign(TupleVector::reference dst,
		    TupleVector::const_reference src);
	//@}

	/** @name property related members

	    dealing with the properties dimension, extending the standard
	    container interface. For efficiency, type checking is left to
	    the class using TupleVector. For both adding and retrieving a
	    property, the type is defined by a template parameter.
	*/
	//@{

	/** reference to a scalar property. In fact this is just an array
	    that is dereferenced using a @ref TupleVector::reference
	*/
	template<class T, class Property>
	class PropertyReferenceBase {
	protected:
	    /// property that is referenced
	    Property &property;
	    /// data shortcut, already type-casted
	    T *data;

	    /// constructor
	    PropertyReferenceBase(Property &_property)
		: property(_property),
		  data(static_cast<T *>(_property.data)) {};

	public:
	    /// dereference
	    T &operator[](reference n) { return data[n]; }
	    /// dereference
	    const T &operator[](const_reference n) const { return data[n]; }
	};

	/// see @ref TupleVector::PropertyReferenceBase
	template<class T>
	class PropertyReference: public PropertyReferenceBase<T, Property> {
	    friend class TupleVector;
	    PropertyReference(Property &_property)
		: PropertyReferenceBase<T, Property>(_property) {};
	};

	/// see @ref TupleVector::PropertyReferenceBase; this is for a constant property
	template<class T>
	class ConstPropertyReference
	    : public PropertyReferenceBase<const T, const Property> {
	    friend class TupleVector;
	    ConstPropertyReference(const Property &_property)
		: PropertyReferenceBase<const T, const Property>(_property) {};
	public:
	    /// non-const->const conversion
	    ConstPropertyReference(const PropertyReference<T> &_ref)
		: PropertyReferenceBase<const T, const Property>(_ref.property) {}
	};

	/** reference to an array property. In contrast to
	    @ref TupleVector::PropertyReference, this returns a pointer, since
	    the element itself is in fact an array whose size
	    is only known at runtime.
	*/
	template<class T, class Property>
	class ArrayPropertyReferenceBase {
	protected:
	    /// reference to property
	    Property &property;
	    /// data short cut from property
	    T *data;
	    /// number of values per element, short cut from property
	    int dimension;

	    /// constructor
	    ArrayPropertyReferenceBase(Property &_property)
		: property(_property),
		  data(static_cast<T *>(_property.data)),
		  dimension(_property.dimension) {};

	public:
	    /// dereference
	    T *operator[](reference n) { return data + n*dimension; }
	    /// dereference
	    const T *operator[](const_reference n) const { return data + n*dimension; }
	};

	/** see @ref TupleVector::ArrayPropertyReferenceBase
	    @tparam T the data type of the property (or more precisely,
	    the datatype to interpret the property as)
	*/
	template<class T>
	class ArrayPropertyReference: public ArrayPropertyReferenceBase<T, Property> {
	    friend class TupleVector;

	    ArrayPropertyReference(Property &_property)
		: ArrayPropertyReferenceBase<T, Property>(_property) {};
	};

	/** see @ref TupleVector::ArrayPropertyReferenceBase;
	    this is for a constant property
	    @tparam T the data type of the property (or more precisely,
	    the datatype to interpret the property as)
	*/
	template<class T>
	class ConstArrayPropertyReference
	    : public ArrayPropertyReferenceBase<const T, const Property> {
	    friend class TupleVector;

	    ConstArrayPropertyReference(const Property &_property)
		: ArrayPropertyReferenceBase<const T, const Property>(_property) {};
	public:
	    /// non-const->const conversion
	    ConstArrayPropertyReference(const ArrayPropertyReference<T> &_ref)
		: ArrayPropertyReferenceBase<const T, const Property>(_ref.property) {}
	};

	/** add a property
	    @param dim dimensionality of the stored entries
	    @tparam T type of the property
	    @return ID of the property
	*/
	template<class T>
	size_t addProperty(size_t dim = 1) {
	    return addProperty(sizeof(T)*dim, dim);
	}

	/** delete a property
	    @param n ID of the property to delete as obtained from addProperty
	*/
	void eraseProperty(size_t n);

	/** get a non-array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T>
	PropertyReference<T> getProperty(size_t n) {
	    return PropertyReference<T>(property[n]);
	}

	/** get a constant non-array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T>
	ConstPropertyReference<T> getProperty(size_t n) const {
	    return ConstPropertyReference<T>(property[n]);
	}
  
	/** get an array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T>
	ArrayPropertyReference<T> getArrayProperty(size_t n) {
	    return ArrayPropertyReference<T>(property[n]);
	}

	/** get an array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T>
	ConstArrayPropertyReference<T> getArrayProperty(size_t n) const {
	    return ConstArrayPropertyReference<T>(property[n]);
	}

	/// get number of properties
	size_type getNumProperties() const { return property.size(); }
  
	//@}

    private:

	/// not accessible and not implemented
	TupleVector(const TupleVector &);
	/// not accessible and not implemented
	void operator=(const TupleVector &);

	/** move around some particles, with potentially overlapping
	    areas. Slower than @ref TupleVector::memcpy.
	    @param dst index of first particle of destination region
	    @param src index of first particle of source region
	    @param size number of particles to move
	*/
	void memmove(size_type dst, size_type src, size_type size);

	/** move around some particles. This is faster, but
	    areas cannot overlap.
	    @param dst index of first particle of destination region
	    @param src index of first particle of source region
	    @param size number of particles to move
	*/
	void memcpy(size_type dst, size_type src, size_type size);

	/// actual function adding a property, with resolved type size
	size_t addProperty(size_t size, size_t _dimension);
    };
}
#endif
