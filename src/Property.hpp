#ifndef _PROPERTY_HPP
#define _PROPERTY_HPP

#include <vector>
#include <algorithm>
#include "particles/Storage.hpp"
#include "Particle.hpp"

namespace espresso {
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
  class Property {
    typedef particles::PropertyHandle< T > Handle;
    typedef particles::ConstPropertyHandle< T > ConstHandle;
  public: // visible in Python
    typedef shared_ptr< Property< T > > SelfPtr;

    Property(particles::Storage::SelfPtr _storage)
      : storage(_storage) {
      id = storage->template addProperty<T>();
    }
    ~Property() {
      storage->deleteProperty(id);
    }

    T getItem(ParticleId part) const {
      return (*this).at(part);
    }
    void setItem(ParticleId part, const T &v) {
      (*this).at(part) = v;
    }

  public: // invisible in Python
    ConstHandle getConstHandle() const {
      return storage->template getConstPropertyHandle< T >(id);
    }

    Handle getHandle() {
      return storage->template getPropertyHandle< T >(id);
    }

    operator ConstHandle() const { return getConstHandle(); }
    operator Handle() { return getHandle(); }

    T operator[](particles::ConstParticleHandle particle) const {
      return ConstHandle(*this)[particle];
    }

    T &operator[](particles::ParticleHandle particle) {
      return Handle(*this)[particle];
    }

    T operator[](ParticleId id) const {
      particles::ConstParticleHandle handle =
        storage->getParticleHandle(id);
      return (*this)[handle];
    }

    T &operator[](ParticleId id) {
      particles::ParticleHandle handle =
        storage->getParticleHandle(id);
      return (*this)[handle];
    }

    /// checked access
    T at(ParticleId id) const {
      particles::ConstParticleHandle handle =
        storage->getParticleHandle(id);
      if (handle == particles::ConstParticleHandle()) {
	throw std::out_of_range("Property::at");
      }
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

    particles::Storage::SelfPtr getStorage() { return storage; }

    particles::IdPropertyHandle getIdHandle() { 
      return storage->getIdPropertyHandle(); 
    }

  private:
    particles::Storage::SelfPtr storage;
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
  class ArrayProperty  {
    typedef particles::ArrayPropertyHandle<T> Handle;
    typedef particles::ConstArrayPropertyHandle<T> ConstHandle;

  public: // visible in Python
    typedef shared_ptr< ArrayProperty< T > > SelfPtr;

    std::vector<T> getItem(ParticleId part) const {
      ConstHandle ref = *this;
      const T
        *start = (*this).at(part),
        *end   = start + ref.getDimension();
      return std::vector<T>(start, end);
    }
    void setItem(ParticleId part, const std::vector<T> &v) {
      Handle ref = *this;
      if (v.size() != ref.getDimension())
        throw std::range_error("ArrayProperty::setItem: incorrect dimension");
      std::copy(v.begin(), v.end(), (*this).at(part));
    }

  public: // invisible in Python
    ArrayProperty(particles::Storage::SelfPtr _storage,
                  size_t dimension)
      : storage(_storage) {
      id = _storage->template addProperty<T>(dimension);
    }
    ~ArrayProperty() {
      storage->deleteProperty(id);
    }

    operator ConstHandle() const {
      return storage->template getConstArrayPropertyHandle<T>(id);
    }

    operator Handle() {
      return storage->template getArrayPropertyHandle<T>(id);
    }

    const T *operator[](particles::ConstParticleHandle handle) const {
      return ConstHandle(*this)[handle];
    }

    T *operator[](particles::ParticleHandle handle) {
      return Handle(*this)[handle];
    }

    const T *operator[](ParticleId part) const {
      particles::ConstParticleHandle handle =
        storage->getParticleHandle(part);
      return (*this)[handle];
    }

    T *operator[](ParticleId part) {
      particles::ParticleHandle handle =
        storage->getParticleHandle(part);
      return (*this)[handle];
    }

    /// checked access
    const T *at(ParticleId part) const {
      particles::ConstParticleHandle handle =
        storage->getParticleHandle(part);
      if (handle == particles::ConstParticleHandle()) {
	throw std::out_of_range("ArrayProperty::at");
      }
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

    particles::Storage::SelfPtr getStorage() { return storage; }

    particles::IdPropertyHandle getIdPropertyHandle() { 
      return storage->getIdPropertyHandle(); 
    }

  private:
    particles::Storage::SelfPtr storage;
    esutil::TupleVector::PropertyId id;
  };

  void registerPythonProperties();
}

#endif
