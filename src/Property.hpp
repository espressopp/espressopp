#ifndef _PROPERTY_HPP
#define _PROPERTY_HPP

#include <vector>
#include <algorithm>
#include "particles/PropertyHandle.hpp"
#include "particles/ParticleHandle.hpp"
#include "particles/Storage.hpp"
#include "Particle.hpp"

namespace espresso {
  /** Exception thrown when a property and a storage don't match. */
  class StorageMismatch : public std::exception 
  {};

  class PropertyBase {
  public:
    typedef shared_ptr< PropertyBase > SelfPtr;

    PropertyBase(const particles::Storage::SelfPtr _storage);

    const particles::Storage::SelfPtr getStorage() const;

    void checkFitsTo(const particles::Set::SelfPtr set) const;

  protected:
    // what storage does the property belong to
    const shared_ptr< particles::Storage > storage;
  };

  /** Scalar particle property. This does several things:
      <ul>
      <li> Lifetime management of the property. Deletion
      of the Property object removes the property from
      its storage.
      <li> Slow, but direct access to the particle properties
      without the need for handles.
      </ul>

      @tparam T type stored in the property
  */
  // TODO: Should be noncopyable!
  //  class Property :  boost::noncopyable {
  template < typename T >
  class Property : public PropertyBase {
    typedef particles::PropertyHandle< T > Handle;

  public: // visible in Python
    typedef shared_ptr< Property< T > > SelfPtr;

    Property(particles::Storage::SelfPtr _storage)
      : PropertyBase(_storage), 
	id(storage->template addProperty<T>())
    {}

    ~Property() {
      storage->deleteProperty(id);
    }

    T getItem(ParticleId part) {
      return (*this).at(part);
    }

    void setItem(ParticleId part, const T &v) {
      (*this).at(part) = v;
    }

  public: // invisible in Python
    Handle getHandle(particles::Set::SelfPtr set) {
      checkFitsTo(set);
      return getHandle();
    }

    T &operator[](particles::ParticleHandle particle) {
      return getHandle()[particle];
    }

    T &operator[](ParticleId id) {
      particles::ParticleHandle handle =
        storage->getParticleHandle(id);
      return (*this)[handle];
    }

    /// checked access
    T &at(ParticleId id) {
      particles::ParticleHandle handle =
        storage->getParticleHandle(id);
      if (handle == particles::ParticleHandle()) {
	throw std::out_of_range("Property::at");
      }
      return (*this)[handle];
    }

  private:
    Handle getHandle() {
      return storage->template getPropertyHandle< T >(id);
    }

    esutil::TupleVector::PropertyId id;
  };


  /** Array particle property. This does several things:
      <ul>
      <li> lifetime management of the property. deletion
      of the Property object removes the property from
      its storage.
      <li> slow, but direct access to the particle properties
      without the need for handles.
      </ul>

      @tparam T type stored in the property
  */
  template < typename T >
  class ArrayProperty : public PropertyBase {
    typedef particles::ArrayPropertyHandle<T> Handle;

  public: // visible in Python
    typedef shared_ptr< ArrayProperty< T > > SelfPtr;

    std::vector<T> getItem(ParticleId part) {
      Handle ref = getHandle();
      const T
        *start = (*this).at(part),
        *end   = start + ref.getDimension();
      return std::vector<T>(start, end);
    }

    void setItem(ParticleId part, const std::vector<T> &v) {
      Handle ref = getHandle();
      if (v.size() != ref.getDimension())
        throw std::range_error("ArrayProperty::setItem: incorrect dimension");
      std::copy(v.begin(), v.end(), (*this).at(part));
    }

  public: // invisible in Python
    ArrayProperty(particles::Storage::SelfPtr _storage,
                  size_t dimension)
      : PropertyBase(_storage), 
	id(_storage->template addProperty<T>(dimension))
    {}

    ~ArrayProperty() {
      storage->deleteProperty(id);
    }

    Handle getHandle(particles::Storage::SelfPtr _storage) {
      checkFitsTo(_storage);
      return getHandle();
    }

    T *operator[](particles::ParticleHandle handle) {
      return getHandle()[handle];
    }

    T *operator[](ParticleId part) {
      particles::ParticleHandle handle =
        storage->getParticleHandle(part);
      return (*this)[handle];
    }

    /// checked access
    T *at(ParticleId part) {
      particles::ParticleHandle handle =
        storage->getParticleHandle(part);
      if (handle == particles::ParticleHandle()) {
	throw std::out_of_range("ArrayProperty::at");
      }
      return (*this)[handle];
    }

  private:
    Handle getHandle() {
      return storage->template getArrayPropertyHandle<T>(id);
    }

    esutil::TupleVector::PropertyId id;
  };

  void registerPythonProperty();

}

#endif
