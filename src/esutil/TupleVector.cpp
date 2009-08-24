#include "TupleVector.hpp"

#include <cstring>
#include <boost/foreach.hpp>
#include <algorithm>

#include <iostream>

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
  void *mem = maxESize > 0 ? malloc(maxESize*_size) : 0;
  property.push_back(Property(++uniqueID, _size, _dimension, mem));
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
    throw std::out_of_range("TupleVector::getPropertyData: property does not exist");
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
      if (newCapacity > 0) {
        prop.data = realloc(prop.data, newCapacity*prop.size);
      }
      else {
        free(prop.data);
        prop.data = 0;
      }
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

TupleVector::thin_iterator TupleVector::insert(TupleVector::thin_iterator pos)
{
  resize(eSize + 1);
  memmove(pos.index + 1, pos.index, eSize - (pos.index + 1));
  return pos;
}

TupleVector::thin_iterator TupleVector::insert(TupleVector::thin_iterator pos,
                                               TupleVector::const_reference e)
{
  thin_iterator it = insert(pos);
  memmove(it.index, e.index, 1);
  return it;
}

void TupleVector::insert(TupleVector::thin_iterator pos, size_type n)
{
  resize(eSize + n);
  if (n) {
    memmove(pos.index + n, pos.index, eSize - (pos.index + n));
  }
}

TupleVector::thin_iterator TupleVector::erase(TupleVector::thin_iterator pos)
{
  if (pos.index != eSize - 1)
    memmove(pos.index, pos.index + 1, eSize - pos.index - 1);
  resize(eSize - 1);
  return pos;
}

TupleVector::thin_iterator TupleVector::erase(TupleVector::thin_iterator start,
                                              TupleVector::thin_iterator end)
{
  memmove(start.index, end.index, eSize - end.index);
  resize(eSize - (end - start));
  return start;
}

TupleVector::reference &
TupleVector::reference::operator=(TupleVector::const_reference src)
{
  if (index != src.index)
    vector->memcpy(index, src.index, 1);
 
  return *this;
}

TupleVector::iterator
TupleVector::iterator::copy(TupleVector::const_iterator begin,
                            TupleVector::const_iterator end)
{
  size_type len = end - begin;
  vector->memmove(index, begin.index, len);
  return TupleVector::iterator(vector, index + len);
}

