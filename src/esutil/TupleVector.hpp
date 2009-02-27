#ifndef _ESUTIL_TUPLEVECTOR_HPP
#define _ESUTIL_TUPLEVECTOR_HPP

#include <vector>
#include <stdexcept>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/foreach.hpp>
#include <functional>
#include "Identifier.hpp"

namespace esutil {

  /** a vector of elements with dynamically configurable properties.
      Each property is stored in a vector of its own, making loops
      over single properties very fast. However, as a consequence,
      there is no value type, and some vector operations are slightly
      different (e.g., elements are uninitialized)
  */
  class TupleVector {
    // forward declarations of member classes, for making friends :-)
  private:
    class ReferenceBase;
    template<class, class> class IteratorBase;
    template<class, class> class PropertyReferenceBase;
    template<class, class> class ArrayPropertyReferenceBase;
    class ReferenceIndexAccess;

    struct Property;

  public:
    class reference;
    class const_reference;

    class iterator;
    class const_iterator;

    template<class> class PropertyReference;
    template<class> class ConstPropertyReference;

    template<class> class ArrayPropertyReference;
    template<class> class ConstArrayPropertyReference;

    class PropertyId: public Identifier {
    public:
      /// invalid id, valid ones can only be obtained by TupleVector
      PropertyId() {}

    private:
      friend class TupleVector;

      /// constructable only by TupleVector
      PropertyId(size_t _v): Identifier(_v) {}
      /// TupleVector can get the stored number
      operator size_t() const { return v; }
      /// TupleVector can increment to obtain the next id
      PropertyId &operator++() { ++v; return *this; }
    };

  private:
    /** Handler for a property. This really just stores the data;
        (re)allocation and freeing needs to be done from outside */
    struct Property {
      /// unique property identifier
      PropertyId id;
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
      Property(PropertyId _id, size_t _size, size_t _dimension, void *_data)
        : id(_id), size(_size), dimension(_dimension), data(_data) {}
    };

    /** a reference to a constant element, which then can be used to read
        properties of the element using
        @ref esutil::TupleVector::PropertyReference.
	
        Actually, this is nothing but the index of the element. A
        reference can only be obtained from element access to the
        TupleVector class or its iterators.
    */
    class ReferenceBase {
    protected:
      /// constructable only through TupleVector
      ReferenceBase() {}

      /// constructable only through TupleVector
      ReferenceBase(size_t _index): index(_index) {}

      /// index of the referenced particle
      size_t index;

    private:
      // class that grants access to the index
      friend class ReferenceIndexAccess;

      /// not possible
      void operator=(const ReferenceBase &_ref);
    };

    class ReferenceIndexAccess {
    public:
      static size_t getIndex(const ReferenceBase &ref) { return ref.index; }
    };

    /** an iterator over elements. This is basically a reference to
        an element with a mutable index

        @tparam ReferenceType type we return as reference
        (reference or const_reference)
        @tparam CRTP curiously recurrent template pattern
    */
    template<class ReferenceType, class CRTP>
    class IteratorBase: public boost::iterator_facade<
      CRTP,
      ReferenceType,
      boost::random_access_traversal_tag,
      ReferenceType,
      ptrdiff_t>
    {
      template <class,class> friend class IteratorBase;

      friend class boost::iterator_core_access;
      // classes that need to access the index
      friend class TupleVector;

    public:
      /// reference type. Fixes up the iterator_facade for a reference only class
      typedef ReferenceType reference;

      /// default constructor, undefined iterator
      IteratorBase() {}

    protected:
      /// regular constructor
      IteratorBase(size_t _index): index(_index) {}

      /// index of particle that is referenced
      size_t index;

    private:
      template<class OtherReferenceType, class OtherIteratorType>
      bool equal(const IteratorBase<OtherReferenceType, OtherIteratorType> &other) const {
        return index == other.index;
      }

      void increment()       { index++; }
      void decrement()       { index--; }
      void advance(size_t n) { index += n; }
      ReferenceType dereference() const { return ReferenceType(index); }
      ptrdiff_t distance_to(const IteratorBase<ReferenceType, CRTP> &other) const {
        return other.index - index;
      }
    };

    /** reference to a scalar property. In fact this is just an array
        that is dereferenced using a @ref esutil::TupleVector::reference.
        This reference is only guaranteed to be valid between any two
        resizings of the TupleVector.
    */
    template<class T, class Property>
    class PropertyReferenceBase: private ReferenceIndexAccess {
      /// not possible
      void operator=(const PropertyReferenceBase &);

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
      T &operator[](reference n) const { return data[getIndex(n)]; }
      /// dereference
      const T &operator[](const_reference n) const { return data[getIndex(n)]; }
    };

    /** reference to an array property. In contrast to
        @ref TupleVector::PropertyReference, this returns a pointer, since
        the element itself is in fact an array whose size
        is only known at runtime.
    */
    template<class T, class Property>
    class ArrayPropertyReferenceBase: private ReferenceIndexAccess {
      /// not possible
      void operator=(const ArrayPropertyReferenceBase &);

    protected:
      /// data short cut from property
      T *data;
      /// number of values per element, short cut from property
      int dimension;

      /// constructor
      ArrayPropertyReferenceBase(Property &_property)
        : data(static_cast<T *>(_property.data)),
          dimension(_property.dimension) {};
      /// constructor
      ArrayPropertyReferenceBase(T *_data, size_t _dimension)
        : data(_data), dimension(_dimension) {};

    public:
      /// dereference
      T *operator[](reference n) const { return data + getIndex(n)*dimension; }
      /// dereference
      const T *operator[](const_reference n) const { return data + getIndex(n)*dimension; }
    };

  public:
    /// @name standard vector members
    //@{

    ///
    typedef size_t size_type;

    ///
    typedef ptrdiff_t difference_type;
  
    /** a reference to an element, which then can be used to read and
        write properties of the element using TupleVector::PropertyReference.
	
        Actually, this is nothing but a wrapper for an iterator. A
        reference can only be obtained from element access to the
        TupleVector class or its iterators. And from the
        implementation point of view, it doesn't differ from a
        pointer.
    */
    class reference: public ReferenceBase {
      // classes that can construct a reference
      friend class TupleVector;
      template<class, class> friend class IteratorBase;
      // for const->nonconst conversion
      friend class const_reference;

    public:
      /// get pointer to the referenced particle
      iterator operator&() const { return iterator(index); }

    private:
      /// constructable only through TupleVector
      reference(size_type _index): ReferenceBase(_index) {}

      /// not possible
      void operator=(const reference &);
    };

    /** a reference to a constant element, which then can be used to read
        properties of the element using
        @ref esutil::TupleVector::PropertyReference.
	
        Actually, this is nothing but a wrapper for an iterator. A
        reference can only be obtained from element access to the
        TupleVector class or its iterators.
    */
    class const_reference: public ReferenceBase {
      // classes that can construct a const_reference
      friend class TupleVector;
      template<class, class> friend class IteratorBase;

    public:
      /// non-const->const conversion
      const_reference(const reference &_ref): ReferenceBase(_ref) {}
      /// get pointer to the referenced particle
      const_iterator operator&() const { return const_iterator(index); }

    private:

      /// constructable only through @ref esutil::TupleVector
      const_reference(size_type _index): ReferenceBase(_index) {}

      /// not possible
      void operator=(const const_reference &);
    };

    /** a random iterator over elements. This is basically a reference to
        an element with a mutable index
    */
    class iterator: public IteratorBase<reference, iterator> {
      // class that can generate an iterator
      friend class TupleVector;
      friend class TupleVector::reference;
      // for const->nonconst conversion
      friend class const_iterator;

    public:
      /// default constructor, uninitialized iterator
      iterator() {}

    private:
      /// constructable only through TupleVector
      iterator(size_type _index)
        : IteratorBase<reference, iterator>(_index) {}
    };

    /** a random iterator over constant elements. This is basically a
        reference to an element with a mutable index
    */
    class const_iterator: public IteratorBase<const_reference, const_iterator> {
      // class that can generate an iterator
      friend class TupleVector;
      friend class const_reference;

    public:
      /// default constructor
      const_iterator() {}

      /// non-const->const conversion
      const_iterator(const iterator _it)
        : IteratorBase<const_reference, const_iterator>(_it.index) {}

    private:
      /// constructable only through TupleVector
      const_iterator(size_type _index)
        : IteratorBase<const_reference, const_iterator>(_index) {}
    };

    typedef iterator pointer;
    typedef const_iterator const_pointer;

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
            
        @tparam T the data type of the property (or more precisely,
        the datatype to interpret the property as)
    */
    template<class T>
    class PropertyReference: public PropertyReferenceBase<T, Property> {
      friend class TupleVector;
      friend class ConstPropertyReference<T>;

      /// not possible
      void operator=(const PropertyReference &);

      PropertyReference(Property &_property)
        : PropertyReferenceBase<T, Property>(_property) {};
    };

    /** reference to a constant scalar property. In fact this is just
        an array that is dereferenced using a @ref TupleVector::reference.
        This reference is only guaranteed to be valid between any two
        resizings of the TupleVector.
            
        @tparam T the data type of the property (or more precisely,
        the datatype to interpret the property as)
    */
    template<class T>
    class ConstPropertyReference
      : public PropertyReferenceBase<const T, const Property> {
      friend class TupleVector;

      /// not possible
      void operator=(const ConstPropertyReference &);

      ConstPropertyReference(const Property &_property)
        : PropertyReferenceBase<const T, const Property>(_property) {};
    public:
      /// non-const->const conversion
      ConstPropertyReference(const PropertyReference<T> &_ref)
        : PropertyReferenceBase<const T, const Property>(_ref.data) {}
    };

    /** reference to an vector-like property, which consists of
        elements of a fixed number of replications of some data
        type (ie, a classical C-array). In fact this is just an
        array that is dereferenced using a
        @ref TupleVector::reference. This reference is only guaranteed
        to be valid between any two resizings of the TupleVector. Unlike
        the @ref PropertyReference, this returns a pointer, in analogy
        to a standard C-array.
            
        @tparam T the data type of the property (or more precisely,
        the datatype to interpret the property as)
    */
    template<class T>
    class ArrayPropertyReference
      : public ArrayPropertyReferenceBase<T, Property>
    {
      friend class TupleVector;
      friend class ConstArrayPropertyReference<T>;
            
      /// not possible
      void operator=(const ArrayPropertyReference &);

      ArrayPropertyReference(Property &_property)
        : ArrayPropertyReferenceBase<T, Property>(_property) {};
    };

    /** reference to an constant vector-like property, which consists of
        elements of a fixed number of replications of some data
        type (ie, a classical C-array). In fact this is just an
        array that is dereferenced using a
        @ref TupleVector::reference. This reference is only guaranteed
        to be valid between any two resizings of the TupleVector. Unlike
        the @ref PropertyReference, this returns a pointer, in analogy
        to a standard C-array.
            
        @tparam T the data type of the property (or more precisely,
        the datatype to interpret the property as)
    */
    template<class T>
    class ConstArrayPropertyReference
      : public ArrayPropertyReferenceBase<const T, const Property>
    {
      friend class TupleVector;

      /// not possible
      void operator=(const ConstArrayPropertyReference &);

      ConstArrayPropertyReference(const Property &_property)
        : ArrayPropertyReferenceBase<const T, const Property>(_property) {};
    public:
      /// non-const->const conversion
      ConstArrayPropertyReference(const ArrayPropertyReference<T> &_ref)
        : ArrayPropertyReferenceBase<const T, const Property>(_ref.data, _ref.dimension) {}
    };

    //@}

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
        member only. Typically, an element can for example be copied
        using a construct such as *a = *b. Since in this class,
        the reference does not hold the full information,
        these functions need to be provided by the container.
    */
    //@{
    /** copy all data of a particle
        @param dst where to write to
        @param src where to get the data from
    */
    void copy(TupleVector::reference dst,
              TupleVector::const_reference src);
    /** copy all data of a set of particles. This is typically faster
        than looping over the set of particles.
        @param begin start of sequence to copy
        @param end   end of sequence to copy
        @param dst   start of where to copy the data to
    */
    void copy(TupleVector::const_iterator begin,
              TupleVector::const_iterator end,
              TupleVector::iterator dst);
    //@}

    /** @name property related members

        dealing with the properties dimension, extending the standard
        container interface. For efficiency, type checking is left to
        the class using TupleVector. For both adding and retrieving a
        property, the type is defined by a template parameter.
    */
    //@{

    /** add a property
        @param dim dimensionality of the stored entries
        @tparam T type of the property
        @return ID of the property
    */
    template<class T>
    PropertyId addProperty(size_t dim = 1) {
      return addProperty(sizeof(T)*dim, dim);
    }

    /** delete a property
        @param n ID of the property to delete as obtained from addProperty
    */
    void eraseProperty(PropertyId n);

    /** get a non-array property
        @param n ID of the property
        @tparam T type of the property
    */
    template<class T>
    PropertyReference<T> getProperty(PropertyId n) {
      return PropertyReference<T>(const_cast<Property &>(getPropertyData(n)));
    }

    /** get a constant non-array property
        @param n ID of the property
        @tparam T type of the property
    */
    template<class T>
    ConstPropertyReference<T> getProperty(PropertyId n) const {
      return ConstPropertyReference<T>(getPropertyData(n));
    }

    /** get an array property
        @param n ID of the property
        @tparam T type of the property
    */
    template<class T>
    ArrayPropertyReference<T> getArrayProperty(PropertyId n) {
      return ArrayPropertyReference<T>(const_cast<Property &>(getPropertyData(n)));
    }

    /** get an array property
        @param n ID of the property
        @tparam T type of the property
    */
    template<class T>
    ConstArrayPropertyReference<T> getArrayProperty(PropertyId n) const {
      return ConstArrayPropertyReference<T>(getPropertyData(n));
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
    /// handlers for all properties
    std::vector <Property> property;

    /// size in elements of each of the stored properties
    size_t eSize;

    /// reserved size in elements of each of the stored properties
    size_t maxESize;

    PropertyId uniqueID;

    /// capacity is always a multiple of granularity
    int granularity;

    /** capacity can shrink if difference between current size and capacity
        is at least this
    */
    int shrinkThreshold;

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
    PropertyId addProperty(size_t size, size_t _dimension);
    /// actual function searching for a property
    const Property &getPropertyData(PropertyId n) const;

    /// predicate for checking the ID of a property
    struct PredicateMatchPropertyID
      : public std::unary_function<const Property &, bool>{
      PropertyId searchID;
      PredicateMatchPropertyID(PropertyId _searchID): searchID(_searchID) {}
      bool operator()(const Property &ref) { return ref.id == searchID; }
    };
  };
}

/** Teach boost foreach that TupleVector is not copyable.
    This function has to be at global scope.
*/
inline boost::mpl::true_ *
boost_foreach_is_noncopyable(esutil::TupleVector *&, boost::foreach::tag)
{ return 0; }

#endif
