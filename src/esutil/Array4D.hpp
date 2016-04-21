/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _ESUTIL_ARRAY4D_H
#define _ESUTIL_ARRAY4D_H
#include <vector>
#include <stdexcept>
#include "Array2D.hpp"

namespace espressopp {
namespace esutil {
/** \brief A four-dimensional array (i.e. a matrix).
 * Particulary we need it for 4-body potential tables

\param T The type of the object stored in the array.

\par Example:
\code
// create the 10x10 id matrix
Array4D< double > I(10, 10, ,10, 10.0, 0.0);
for (int i = 0; i < 10; i++) I(i, i, i, i) = 1.0;
\endcode
*/
template<class T, OutlierMode outlier_mode = exception>
class Array4D: private std::vector<T> {
  typedef Array4D<T> Self;
  /** \brief The super class of the array (i.e. vector<T>). */
  typedef std::vector<T> Super;

 public:
  /** \brief The type of object, \c T, stored in the Array1D. */
  typedef T value_type;
  /** \brief Reference to to a \c T object. */
  typedef T &reference;
  /** \brief Const reference to to a \c T object. */
  typedef const T &const_reference;
  /** \brief An unsigned integer type. */
  typedef typename Super::size_type size_type;
  /** \brief Iterator used to iterate through an Array1D. */
  typedef typename Super::iterator iterator;
  /** \brief Const-iterator used to iterate through an Array1D. */
  typedef typename Super::const_iterator const_iterator;

  /** \brief Default constructor. */
  Array4D() { clear(); }

  /** \brief Creates a 4-dimensional \p n x \p m x \p l x \p k array,
  initialising the values with copies of \a t. */
  Array4D(size_type n, size_type m, size_type l, size_type k, const T &t = T()) {
    init(n, m, l, k, t);
  }

  /** \brief Erases all elements. */
  void clear() {
    this->n = 0;
    this->m = 0;
    this->p = 0;
    this->q = 0;
    Super::clear();
  }

  /** \brief Creates a 3-dimensional \p n x \p m x \p p x \p q array,
  initialising the values with copies of \a t. */
  void init(size_type n, size_type m, size_type l, size_type k, const T &t = T()) {
    Super::clear();
    resize(n, m, l, k, t);
  }

  /** \brief Changes the size of the Array4D so that the size
  becomes \p n x \p m x \p p x \p q array.

  \attention The elements of the Array4D are not initialised
  consistently, so the values are undefined. The \p t parameter
  exists only for the purpose to provide a default \a T value
  when new elements are created and needs to be specified only
  when \a T has no default constructor.

  */
  void resize(size_type n, size_type m, size_type l, size_type k, const T &t = T()) {
    this->n = n;
    this->m = m;
    this->p = l;
    this->q = k;
    Super::resize(n * m * l * k, t);
  }

  /** \brief Returns the \a n -size of the array. */
  size_type size_n() const { return n; }

  /** \brief Returns the \a m -size  of the array. */
  size_type size_m() const { return m; }

  /** \brief Returns the \a p -size  of the array. */
  size_type size_p() const { return p; }

  /** \brief Returns the \a q -size of the array. */
  size_type size_q() const { return q; }

  /** \brief \c true if the \c Array4D's size is 0 x 0 x 0 x 0. */
  bool empty() const { return n == 0 && m == 0 && p == 0 && q == 0; }

  /** \brief Returns the element at position \p i x \p j x \p k.
  \attention No range checking. */
  reference operator()(size_type i, size_type j, size_type k, size_type l) {
    return Super::operator[](i + j * n + k * n * m + l * n * m * p);
  }

  /** \brief Returns the element at position \p i x \p j x \p k.
  \attention No range checking. */
  const_reference operator()(size_type i, size_type j, size_type k, size_type l) const {
    return Super::operator[](i + j * n + k * n * m + l * n * m * p);
  }

  /** \brief Returns the element at position \p pos.

  \p pos must be a vector where the expressions \c pos[0], \c pos[1]
  and \c pos[2] return the indices.
  \attention No range checking. */
  template<class Vector>
  reference operator()(const Vector &pos) {
    return operator()(pos[0], pos[1], pos[2], pos[3]);
  }

  /** \brief Returns the element at position \p pos.

  \p pos must be a vector where the expressions \c pos[0], \c pos[1] and
  \c pos[2] return the indices.
  \attention No range checking. */
  template<class Vector>
  const_reference operator()(const Vector &pos) const {
    return operator()(pos[0], pos[1], pos[2], pos[3]);
  }

  /** \brief Returns the element at position \p i x \p j x \p k.
  \exception out_of_range \p i is out of range of the Array4D */
  reference at(size_type i, size_type j, size_type k, size_type l) {
    if (i >= n) throw std::out_of_range("Index i out of bounds.");
    if (j >= m) throw std::out_of_range("Index j out of bounds.");
    if (k >= p) throw std::out_of_range("Index k out of bounds.");
    if (l >= q) throw std::out_of_range("Index l out of bounds.");
    return operator()(i, j, k, l);
  }

  /** \brief Returns the element at position \p i x \p j x \p k.
  \exception out_of_range \p i is out of range of the Array4D */
  const_reference at(size_type i, size_type j, size_type k, size_type l) const {
    if (i >= n) throw std::out_of_range("Index i out of bounds.");
    if (j >= m) throw std::out_of_range("Index j out of bounds.");
    if (k >= p) throw std::out_of_range("Index k out of bounds.");
    if (l >= q) throw std::out_of_range("Index l out of bounds.");
    return operator()(i, j, k, l);
  }

  /** \brief Returns the element at position \p pos.

  \p pos must be a vector where the expressions \c pos[0], \c pos[1]
  and \c pos[2] return the first and second element.

  \exception out_of_range \p i is out of range of the Array4D */
  template<class Vector>
  reference at(const Vector &pos) { return at(pos[0], pos[1], pos[2], pos[3]); }

  /** \brief Returns the element at position \p pos.

  \p pos must be a vector where the expressions \c pos[0], \c pos[1]
  and \c pos[2] return the first and second element.

  \exception out_of_range \p i is out of range of the Array4D */
  template<class Vector>
  const_reference at(const Vector &pos) const { return at(pos[0], pos[1], pos[2], pos[3]); }

  /** \brief Returns an iterator pointing to the beginning of the
  Array4D. */
  iterator begin() { return Super::begin(); }

  /** \brief Returns a const_iterator pointing to the beginning of the
  Array4D. */
  const_iterator begin() const { return Super::begin(); }

  /** \brief Returns an iterator pointing to the end of the
  Array4D. */
  iterator end() { return Super::end(); }

  /** \brief Returns a const_iterator pointing to the end of the
  Array4D. */
  const_iterator end() const { return Super::end(); }

  /** \brief Assigns t2 to each of the elements. */
  template<class T2>
  void operator=(const T2 &t2) {
    for (iterator it = begin(); it != end(); it++) *it = t2;
  }

  /** \brief Adds t2 to each of the elements.  */
  template<class T2>
  void operator+=(const T2 &t2) {
    for (iterator it = begin(); it != end(); it++) *it += t2;
  }

  /** \brief Adds t2 to each of the elements.  */
  template<class T2>
  void operator-=(const T2 &t2) {
    for (iterator it = begin(); it != end(); it++) *it -= t2;
  }

  /** \brief Multiplies each of the elements by t2. */
  template<class T2>
  void operator*=(const T2 &t2) { for (iterator it = begin(); it != end(); it++) *it *= t2; }

  /** \brief Divides each of the elements by t2. */
  template<class T2>
  void operator/=(const T2 &t2) { for (iterator it = begin(); it != end(); it++) *it /= t2; }

  void operator/=(const double &x) {
    double x_rec = 1.0 / x;
    for (iterator it = begin(); it != end(); it++) *it *= x_rec;
  }

  /** \brief Equality operator. */
  bool operator==(const Array4D<T, exception> &a) {
    return static_cast< std::vector<T> >(*this) ==
        static_cast< std::vector<T> >(a);
  }

  /** \brief Inequality operator. */
  bool operator!=(const Array4D<T, exception> &a) {
    return static_cast< std::vector<T> >(*this) !=
        static_cast< std::vector<T> >(a);
  }

 private:
  /** \brief The \p n size of the \p n x \p m x \p p x \p q array */
  size_type n;
  /** \brief The \p m size of the \p n x \p m x \p p x \p q array */
  size_type m;
  /** \brief The \p p size of the \p n x \p m x \p p x \p q array */
  size_type p;
  /** \brief The \p p size of the \p n x \p m x \p p x \p q array */
  size_type q;
};

template<class T>
class Array4D<T, enlarge>: public Array4D<T, exception> {
  /** \brief The super class of the array (i.e. vector<T>). */
  typedef Array4D<T, exception> Super;

 public:
  /** \brief The type of object, \c T, stored in the Array1D. */
  typedef T value_type;
  /** \brief Reference to to a \c T object. */
  typedef T &reference;
  /** \brief Const reference to to a \c T object. */
  typedef const T &const_reference;
  /** \brief An unsigned integer type. */
  typedef typename Super::size_type size_type;
  /** \brief Iterator used to iterate through an Array1D. */
  typedef typename Super::iterator iterator;
  /** \brief Const-iterator used to iterate through an Array1D. */
  typedef typename Super::const_iterator const_iterator;

  using Super::size_n;
  using Super::size_m;
  using Super::size_p;
  using Super::size_q;

  /** \brief Default constructor. */
  explicit Array4D(const T &prototype = T()) {
    init(0, 0, 0, 0, prototype);
  }

  /** \brief Creates a one-dimensional array of size \a n with \a n
  copies of \a prototype.
  Also sets the \p prototype. */
  Array4D(size_type n, size_type m, size_type p, size_type q, const T &prototype = T()) {
    init(n, m, p, q, prototype);
  }

  /** \brief Creates a one-dimensional array of size \a n with \a n
  copies of \a prototype.
  Also sets the \p prototype. */
  void init(size_type n, size_type m, size_type p, size_type q, const T &prototype = T()) {
    set_prototype(prototype);
    Super::init(n, m, p, q, prototype);
  }

  /** \brief Sets the \p prototype used when the Array4D<T,enlarge>
  is automatically resized. */
  void set_prototype(const T &prototype = T()) {
    this->prototype = prototype;
  }

  /** \brief Returns the element at position \p i x \p j x \p p x \p q.

  When the element doesn't exist, the Array4D is automatically
  extended so that the element exists. The new elements are
  initialised with copies of the prototype.
  */


  reference at(size_type i, size_type j, size_type k, size_type l) {
    if (i >= size_n() || j >= size_m() || k >= size_p() || l >= size_q()) {
      Super old = *this;
      size_type new_n = size_n();
      size_type new_m = size_m();
      size_type new_p = size_p();
      size_type new_q = size_q();
      if (i >= new_n) new_n = i + 1;
      if (j >= new_m) new_m = j + 1;
      if (k >= new_p) new_p = k + 1;
      if (l >= new_q) new_q = l + 1;
      Super::init(new_n, new_m, new_p, new_q, prototype);

      // copy the old Array4D
      for (size_type ii = 0; ii < old.size_n(); ii++)
        for (size_type jj = 0; jj < old.size_m(); jj++)
          for (size_type kk = 0; kk < old.size_p(); kk++)
            for (size_type ll = 0; ll < old.size_q(); ll++)
              (*this)(ii, jj, kk, ll) = old(ii, jj, kk, ll);

      // init the new elements with the prototype
      for (size_type ii = old.size_n(); ii < size_n(); ii++)
        for (size_type jj = 0; jj < size_m(); jj++)
          for (size_type kk = 0; kk < size_p(); kk++)
            for (size_type ll = 0; ll < size_q(); ll++)
              (*this)(ii, jj, kk, ll) = prototype;

      for (size_type ii = 0; ii < old.size_n(); ii++)
        for (size_type jj = old.size_m(); jj < size_m(); jj++)
          for (size_type kk = 0; kk < size_p(); kk++)
            for (size_type ll = 0; ll < size_q(); ll++)
              (*this)(ii, jj, kk, ll) = prototype;

      for (size_type ii = 0; ii < old.size_n(); ii++)
        for (size_type jj = 0; jj < old.size_m(); jj++)
          for (size_type kk = old.size_p(); kk < size_p(); kk++)
            for (size_type ll = 0; ll < size_q(); ll++)
              (*this)(ii, jj, kk, ll) = prototype;

      for (size_type ii = 0; ii < old.size_n(); ii++)
        for (size_type jj = 0; jj < old.size_m(); jj++)
          for (size_type kk = 0; kk < old.size_p(); kk++)
            for (size_type ll = old.size_q(); ll < size_q(); ll++)
              (*this)(ii, jj, kk, ll) = prototype;
    }
    return Super::at(i, j, k, l);
  }


  /** \brief Returns the element at position \p pos.

  \p pos must be a vector where the expressions \c pos[0], \c pos[1]
  and \c pos[2] return the first and second element.

  When the element doesn't exist, the Array4D is automatically
  extended so that the element exists. The new elements are
  initialised with copies of the prototype.
  */
  template<class Vector>
  reference at(const Vector &pos) { return at(pos[0], pos[1], pos[2], pos[3]); }

 private:
  value_type prototype;
};

}  // namespace esutil
}  // namespace espressopp
#endif
