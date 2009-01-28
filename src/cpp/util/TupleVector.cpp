#include "TupleVector.hpp"

#include <cstring>
#include <boost/foreach.hpp>

using namespace util;

/// capacity is always a multiple of granularity
static const int granularity = 8;

/** capacity can shrink if difference between current size and capacity
    is at least this
*/
static const int shrinkThreshold = 4*granularity;

TupleVector::~TupleVector()
{
    // deallocate
    BOOST_FOREACH(Property &prop, property) {
	free(prop.data);
    }
}

TupleVector::TupleVector(const TupleVector &vec, size_type newESize)
    : eSize(0), maxESize(0)
{
    // set number of particles for allocation below
    resize(newESize);
    // to speed up copying/allocating
    property.reserve(vec.property.capacity());
    // copy the properties, but allocate new property data arrays
    BOOST_FOREACH(const Property &otherProp, vec.property) {
	addProperty(otherProp.size, otherProp.dimension);
    }
}

size_t TupleVector::addProperty(size_t _size, size_t _dimension)
{
    property.push_back(Property(_size, _dimension, malloc(maxESize)));
    return property.size() - 1;
}

void TupleVector::eraseProperty(size_t id)
{
    std::vector<Property>::iterator it = property.begin() + id;
    free(it->data);
    property.erase(it);
}

void TupleVector::reserve(size_t minCapacity)
{
    // cannot shrink below actual size
    if (minCapacity < eSize) {
	minCapacity = eSize;
    }
    // round up for granularity
    minCapacity = granularity*((minCapacity + granularity - 1)/granularity);
    size_t newCapacity = maxESize;

    // do we need increase or can we free?
    if ((minCapacity + shrinkThreshold <= maxESize) ||
	(minCapacity > maxESize)) {
	newCapacity = minCapacity;
    }

    // finally, reallocate if necessary
    if (newCapacity != maxESize) {
	BOOST_FOREACH(Property &prop, property) {
	    prop.data = realloc(prop.data, newCapacity*prop.size);
	}
    }
    maxESize = newCapacity;
}

void TupleVector::memmove(size_type dst, size_type src, size_type size)
{
    BOOST_FOREACH(Property &prop, property) {
	/* here we have only memory addresses (void *) and element
	   sizes (size_t) to do arithmetics, temporarily we cast to
	   char *.  Ugly, but according to C/C++ standard. The
	   rationale is that sizeof always returns the object size in
	   multiples of sizeof(char).
	*/
	char *data = static_cast<char *>(prop.data);
	std::memmove(data + dst*prop.size,
		     data + src*prop.size, size*prop.size);
    }
}

void TupleVector::memcpy(size_type dst, size_type src, size_type size)
{
    BOOST_FOREACH(Property &prop, property) {
	/* here we have only memory addresses (void *) and element
	   sizes (size_t) to do arithmetics, temporarily we cast to
	   char *.  Ugly, but according to C/C++ standard. The
	   rationale is that sizeof always returns the object size in
	   multiples of sizeof(char).
	*/
	char *data = static_cast<char *>(prop.data);
	std::memcpy(data + dst*prop.size, data + src*prop.size, size*prop.size);
    }
}

TupleVector::iterator TupleVector::insert(TupleVector::iterator pos)
{
    resize(eSize + 1);
    memmove(pos.index + 1, pos.index, eSize - (pos.index + 1));
    return pos;
}

TupleVector::iterator TupleVector::insert(TupleVector::iterator pos,
					  TupleVector::const_reference e)
{
    iterator it = insert(pos);
    assign(*it, e);
    return it;
}

void TupleVector::insert(TupleVector::iterator pos, size_type n)
{
    resize(eSize + n);
    if (n) {
	memmove(pos.index + n, pos.index, eSize - (pos.index + n));
    }
}

TupleVector::iterator TupleVector::erase(TupleVector::iterator pos)
{
    resize(eSize - 1);
    if (pos.index != eSize)
	memmove(pos.index, pos.index + 1, eSize - pos.index);
    return pos;
}

TupleVector::iterator TupleVector::erase(TupleVector::iterator start,
					 TupleVector::iterator end)
{
    resize(eSize - (end - start));
    memmove(start.index, end.index, eSize - start.index);
    return start;
}

void TupleVector::assign(TupleVector::reference dst,
			 TupleVector::const_reference src)
{
    if (dst.index != src.index)
	memmove(dst.index, src.index, 1);
}

