#ifndef _PROPERTY_HPP
#define _PROPERTY_HPP

#include <boost/shared_ptr.hpp>
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
  template <typename T>
  class Property {
    typedef particles::PropertyHandle<T> RefType;
    typedef particles::ConstPropertyHandle<T> ConstRefType;
  public: // visible in Python
    Property(particles::PStorage _storage)
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

    operator ConstRefType() const {
      return storage->template getConstPropertyHandle<T>(id);
    }

    operator RefType() {
      return storage->template getPropertyHandle<T>(id);
    }

    T operator[](particles::ConstParticleHandle handle) const {
      return ConstRefType(*this)[handle];
    }

    T &operator[](particles::ParticleHandle handle) {
      return RefType(*this)[handle];
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

  private:
    particles::PStorage storage;
    esutil::TupleVector::PropertyId id;
  };

  typedef Property<real> RealProperty;
  typedef Property<Real3D> Real3DProperty;
  typedef Property<int> IntegerProperty;

  typedef boost::shared_ptr< RealProperty > PRealProperty;
  typedef boost::shared_ptr< Real3DProperty > PReal3DProperty;
  typedef boost::shared_ptr< IntegerProperty > PIntegerProperty;

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
  template <typename T>
  class ArrayProperty {
    typedef particles::ArrayPropertyHandle<T> RefType;
    typedef particles::ConstArrayPropertyHandle<T> ConstRefType;

  public: // visible in Python
    std::vector<T> getItem(ParticleId part) const {
      ConstRefType ref = *this;
      const T
        *start = (*this).at(part),
        *end   = start + ref.getDimension();
      return std::vector<T>(start, end);
    }
    void setItem(ParticleId part, const std::vector<T> &v) {
      RefType ref = *this;
      if (v.size() != ref.getDimension())
        throw std::range_error("ArrayProperty::setItem: incorrect dimension");
      std::copy(v.begin(), v.end(), (*this).at(part));
    }

  public: // invisible in Python
    ArrayProperty(particles::PStorage _storage,
                  size_t dimension)
      : storage(_storage) {
      id = _storage->template addProperty<T>(dimension);
    }
    ~ArrayProperty() {
      storage->deleteProperty(id);
    }

    operator ConstRefType() const {
      return storage->template getConstArrayPropertyHandle<T>(id);
    }

    operator RefType() {
      return storage->template getArrayPropertyHandle<T>(id);
    }

    const T *operator[](particles::ConstParticleHandle handle) const {
      return ConstRefType(*this)[handle];
    }

    T *operator[](particles::ParticleHandle handle) {
      return RefType(*this)[handle];
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

  private:
    particles::PStorage storage;
    esutil::TupleVector::PropertyId id;
  };

  typedef ArrayProperty<int> IntegerArrayProperty;
  typedef ArrayProperty<real> RealArrayProperty;

  typedef boost::shared_ptr< RealArrayProperty > PRealArrayProperty;
  typedef boost::shared_ptr< IntegerArrayProperty > PIntegerArrayProperty;

  void registerPythonProperties();
}

#endif
