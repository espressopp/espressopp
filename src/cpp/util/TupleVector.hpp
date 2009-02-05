#ifndef _TUPLEVECTOR_HPP
#define _TUPLEVECTOR_HPP

#include <vector>
#include <stdexcept>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/foreach.hpp>
#include <functional>

namespace util {

    /** a vector of elements with dynamically configurable properties.
	Each property is stored in a vector of its own, making loops
	over single properties very fast. However, as a consequence,
	there is no value type, and some vector operations are slightly
	different (e.g., elements are uninitialized)
    */
    class TupleVector {
    public:
	/** a reference to an element, which then can be used to read and
	    write properties of the element using TupleVector::PropertyReference.
	
	    Actually, this is nothing but the index of the element. A
	    reference can only be obtained from element access to the
	    TupleVector class or its iterators.
	*/
	class ReferenceBase {
	protected:
	    /// index of the reference particle
	    size_t index;

	    /// constructor
	    ReferenceBase(size_t _index): index(_index) {}
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
	    ptrdiff_t>
	{
	protected:
	    /// index of particle that is referenced
	    size_t index;

	    /// regular constructor
	    IteratorBase(size_t _index): index(_index) {}
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

	    void increment()       { index++; }
	    void decrement()       { index--; }
	    void advance(size_t n) { index += n; }
	    reference dereference() const { return reference(index); }
	    ptrdiff_t distance_to(const IteratorBase<ReferenceType, CRTP> &other) const{
		return other.index - index;
	    }
	};

	/// @name standard vector members
	//@{

	///
	typedef size_t size_type;

	///
	typedef ptrdiff_t difference_type;
  
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

    private:
	/** Handler for a property. This really just stores the data;
	    (re)allocation and freeing needs to be done from outside */
	struct Property {
	    /// unique property identifier
	    size_t id;
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
	    Property(size_t _id, size_t _size, size_t _dimension, void *_data)
		: id(_id), size(_size), dimension(_dimension), data(_data) {}
	};

	/// handlers for all properties
	std::vector <Property> property;

	/// size in elements of each of the stored properties
	size_t eSize;

	/// reserved size in elements of each of the stored properties
	size_t maxESize;

	size_t uniqueID;

	/// capacity is always a multiple of granularity
	int granularity;

	/** capacity can shrink if difference between current size and capacity
	    is at least this
	*/
	int shrinkThreshold;

    public:
	/** default constructor. Constructs an empty TupleVector, but possibly
	    allocating already some particles. When properties are added, they
	    will offer space for the specified number of particles.
	*/
	explicit TupleVector(size_type _n = 0);

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
	    that is dereferenced using a @ref TupleVector::reference.
	    This reference is only guaranteed to be valid between any two
	    resizings of the TupleVector.
	*/
	template<class T, class Property>
	class PropertyReferenceBase {
	protected:
	    /// data shortcut, already type-casted
	    T *data;

	    /// constructor
	    PropertyReferenceBase(Property &_property)
		: data(static_cast<T *>(_property.data)) {};
	    /// constructor
	    PropertyReferenceBase(T *_data)
		: data(_data) {};

	public:
	    /// dereference
	    T &operator[](reference n) { return data[n.index]; }
	    /// dereference
	    const T &operator[](const_reference n) const { return data[n.index]; }
	};

	/** reference to an array property. In contrast to
	    @ref TupleVector::PropertyReference, this returns a pointer, since
	    the element itself is in fact an array of constant size.

            @tparam T class of the data type this reference
            @tparam dimension number of values of type T per element
	*/
	template<class T, size_t dimension, class Property>
	class ArrayPropertyReferenceBase {
	protected:
	    /// data short cut from property
	    T *data;

	    /// constructor
	    ArrayPropertyReferenceBase(Property &_property)
		: data(static_cast<T *>(_property.data)) {};
	    /// constructor
	    ArrayPropertyReferenceBase(T *_data): data(_data) {};

	public:
	    /// dereference
	    T *operator[](reference n) { return data + n.index*dimension; }
	    /// dereference
	    const T *operator[](const_reference n) const { return data + n.index*dimension; }
	};

	/** reference to an array property. In contrast to
	    @ref TupleVector::PropertyReference, this returns a pointer, since
	    the element itself is in fact an array whose size
	    is only known at runtime.
	*/
	template<class T, class Property>
	class VarArrayPropertyReferenceBase {
	protected:
	    /// data short cut from property
	    T *data;
	    /// number of values per element, short cut from property
	    int dimension;

	    /// constructor
	    VarArrayPropertyReferenceBase(Property &_property)
		: data(static_cast<T *>(_property.data)),
		  dimension(_property.dimension) {};
	    /// constructor
	    VarArrayPropertyReferenceBase(T *_data, size_t _dimension)
		: data(_data), dimension(_dimension) {};

	public:
	    /// dereference
	    T *operator[](reference n) { return data + n.index*dimension; }
	    /// dereference
	    const T *operator[](const_reference n) const { return data + n.index*dimension; }
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
		: PropertyReferenceBase<const T, const Property>(_ref.data) {}
	};

	/** see @ref TupleVector::ArrayPropertyReferenceBase
	    @tparam T the data type of the property (or more precisely,
	    the datatype to interpret the property as)
	*/
	template<class T, size_t dimension>
	class ArrayPropertyReference: public ArrayPropertyReferenceBase<T, dimension, Property> {
	    friend class TupleVector;

	    ArrayPropertyReference(Property &_property)
		: ArrayPropertyReferenceBase<T, dimension, Property>(_property) {};
	};

	/** see @ref TupleVector::ArrayPropertyReferenceBase;
	    this is for a constant property
	    @tparam T the data type of the property (or more precisely,
	    the datatype to interpret the property as)
	*/
	template<class T, size_t dimension>
	class ConstArrayPropertyReference
	    : public ArrayPropertyReferenceBase<const T, dimension, const Property> {
	    friend class TupleVector;

	    ConstArrayPropertyReference(const Property &_property)
		: ArrayPropertyReferenceBase<const T, dimension, const Property>(_property) {};
	public:
	    /// non-const->const conversion
	    ConstArrayPropertyReference(const ArrayPropertyReference<T, dimension> &_ref)
		: ArrayPropertyReferenceBase<const T, dimension, const Property>(_ref.data) {}
	};

	/** see @ref TupleVector::ArrayPropertyReferenceBase
	    @tparam T the data type of the property (or more precisely,
	    the datatype to interpret the property as)
	*/
	template<class T>
	class VarArrayPropertyReference
            : public VarArrayPropertyReferenceBase<T, Property>
        {
	    friend class TupleVector;
            
	    VarArrayPropertyReference(Property &_property)
		: VarArrayPropertyReferenceBase<T, Property>(_property) {};
	};

	/** see @ref TupleVector::ArrayPropertyReferenceBase;
	    this is for a constant property
	    @tparam T the data type of the property (or more precisely,
	    the datatype to interpret the property as)
	*/
	template<class T>
	class ConstVarArrayPropertyReference
	    : public VarArrayPropertyReferenceBase<const T, const Property>
        {
	    friend class TupleVector;

	    ConstVarArrayPropertyReference(const Property &_property)
		: VarArrayPropertyReferenceBase<const T, const Property>(_property) {};
	public:
	    /// non-const->const conversion
	    ConstVarArrayPropertyReference(const VarArrayPropertyReference<T> &_ref)
		: VarArrayPropertyReferenceBase<const T, const Property>(_ref.data, _ref.dimension) {}
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
	    return PropertyReference<T>(const_cast<Property &>(getPropertyData(n)));
	}

	/** get a constant non-array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T>
	ConstPropertyReference<T> getProperty(size_t n) const {
	    return ConstPropertyReference<T>(getPropertyData(n));
	}
  
	/** get an array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T, size_t dimension>
	ArrayPropertyReference<T, dimension> getArrayProperty(size_t n) {
            const Property &prop = getPropertyData(n);
            if (prop.dimension != dimension) {
                throw std::range_error("getArrayProperty: property size mismatch");
            }
	    return ArrayPropertyReference<T, dimension>(const_cast<Property &>(prop));
	}

	/** get an array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T, size_t dimension>
	ConstArrayPropertyReference<T, dimension> getArrayProperty(size_t n) const {
            const Property &prop = getPropertyData(n);
            if (prop.dimension != dimension) {
                throw std::range_error("getArrayProperty: property size mismatch");
            }
	    return ConstArrayPropertyReference<T, dimension>(prop);
	}

	/** get an array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T>
	VarArrayPropertyReference<T> getVarArrayProperty(size_t n) {
	    return VarArrayPropertyReference<T>(const_cast<Property &>(getPropertyData(n)));
	}

	/** get an array property
	    @param n ID of the property
	    @tparam T type of the property
	*/
	template<class T>
	ConstVarArrayPropertyReference<T> getVarArrayProperty(size_t n) const {
	    return ConstVarArrayPropertyReference<T>(getPropertyData(n));
	}

	/// get number of properties
	size_type getNumProperties() const { return property.size(); }
  
	//@}

        /// the capacity of the vector is always resized to a multiple of the granularity
        void setGranularity(size_t _granularity)   { granularity = _granularity; }
        /** if the difference between the current container size and its capacity is
            bigger than the shrinkThreshold, it gets resized to the current size according
            to the current granularity */
        void setShrinkThreshold(size_t _threshold) { shrinkThreshold = _threshold; }

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
	/// actual function searching for a property
	const Property &getPropertyData(size_t n) const;

	/// predicate for checking the ID of a property
	struct PredicateMatchPropertyID
	    : public std::unary_function<const Property &, bool>{
	    size_t searchID;
	    PredicateMatchPropertyID(size_t _searchID): searchID(_searchID) {}
	    bool operator()(const Property &ref) { return ref.id == searchID; }
	};
    };
}

/** Teach boost foreach that TupleVector is not copyable.
    This function has to be at global scope.
*/
inline boost::mpl::true_ *
boost_foreach_is_noncopyable(util::TupleVector *&, boost::foreach::tag)
{ return 0; }

#endif
