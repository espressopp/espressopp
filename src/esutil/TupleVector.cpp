#include "TupleVector.hpp"

#include <cstring>
#include <boost/foreach.hpp>
#include <algorithm>

using namespace espresso::esutil;

static const int defaultGranularity = 1024/sizeof(int);
static const int defaultShrinkThreshold = 4*defaultGranularity;

TupleVector::~TupleVector()
{
    // deallocate
    BOOST_FOREACH(Property &prop, property) {
	free(prop.data);
    }
}

TupleVector::TupleVector(size_t newESize)
    : eSize(0), maxESize(0), uniqueID(0),
      granularity(defaultGranularity), shrinkThreshold(defaultShrinkThreshold)
{
    resize(newESize);
}

TupleVector::TupleVector(const TupleVector &vec, size_type newESize)
    : eSize(0), maxESize(0), uniqueID(0),
      granularity(vec.granularity), shrinkThreshold(vec.shrinkThreshold)
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

TupleVector::PropertyId TupleVector::addProperty(size_t _size, size_t _dimension)
{
  property.push_back(Property(++uniqueID, _size, _dimension, malloc(maxESize)));
  return uniqueID;
}

void TupleVector::eraseProperty(TupleVector::PropertyId id)
{
    std::vector<Property>::iterator it =
	std::find_if(property.begin(), property.end(), PredicateMatchPropertyID(id));
    if (it == property.end()) {
	throw std::out_of_range("TupleVector::eraseProperty: property does not exist");
    }
    free(it->data);
    property.erase(it);
}

const TupleVector::Property &TupleVector::getPropertyData(PropertyId id) const {
    std::vector<Property>::const_iterator it =
	std::find_if(property.begin(), property.end(), PredicateMatchPropertyID(id));
    if (it == property.end()) {
	throw std::out_of_range("TupleVector::eraseProperty: property does not exist");
    }
    return *it;
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
        maxESize = newCapacity;
    }
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
    copy(*it, e);
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

void TupleVector::copy(TupleVector::reference dst,
                       TupleVector::const_reference src)
{
  if (dst.index != src.index)
    memcpy(dst.index, src.index, 1);
}

void TupleVector::copy(TupleVector::const_iterator begin,
                       TupleVector::const_iterator end,
                       TupleVector::iterator dst)
{
    memcpy(dst.index, begin.index, end - begin);
}

