/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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

#ifndef _ESUTIL_ARRAY3D_H
#define _ESUTIL_ARRAY3D_H
#include <vector>
#include <stdexcept>
#include "Array2D.hpp"

namespace espressopp {
  namespace esutil {
    /** \brief A three-dimensional array (i.e. a matrix).
     * Prticulary we need it for 3-body potential tables

	\param T The type of the object stored in the array.

	\par Example:
	\code
	// create the 10x10 id matrix
	Array3D< double > I(10, 10, ,10, 0.0);
	for (int i = 0; i < 10; i++) I(i, i, i) = 1.0;
	\endcode
    */
    template < class T, OutlierMode outlier_mode = exception >
    class Array3D : private std::vector< T > {
      typedef Array3D< T > Self; 
      /** \brief The super class of the array (i.e. vector<T>). */
      typedef std::vector< T > Super;

    public:
      /** \brief The type of object, \c T, stored in the Array1D. */
      typedef T value_type;
      /** \brief Reference to to a \c T object. */
      typedef T& reference;
      /** \brief Const reference to to a \c T object. */
      typedef const T& const_reference;
      /** \brief An unsigned integer type. */
      typedef typename Super::size_type size_type;
      /** \brief Iterator used to iterate through an Array1D. */
      typedef typename Super::iterator iterator;
      /** \brief Const-iterator used to iterate through an Array1D. */
      typedef typename Super::const_iterator const_iterator;

      /** \brief Default constructor. */
      Array3D() { clear(); }

      /** \brief Creates a 3-dimensional \p n x \p m x \p l array,
	  initialising the values with copies of \a t. */
      Array3D(size_type n, size_type m, size_type l, const T &t = T()) 
      { init(n, m, l, t); }

      /** \brief Erases all elements. */
      void clear() { 
        this->n = 0;
        this->m = 0;
        this->l = 0;
        Super::clear();
      }

      /** \brief Creates a 3-dimensional \p n x \p m x \p l array,
	  initialising the values with copies of \a t. */
      void init(size_type n, size_type m, size_type l, const T &t = T()) { 
        Super::clear(); 
        resize(n, m, l, t);
      }

      /** \brief Changes the size of the Array3D so that the size
	  becomes \p n x \p m x \p l.

	  \attention The elements of the Array3D are not initialised
	  consistently, so the values are undefined. The \p t parameter
	  exists only for the purpose to provide a default \a T value
	  when new elements are created and needs to be specified only
	  when \a T has no default constructor.

      */
      void resize(size_type n, size_type m, size_type l, const T &t = T()) {
        this->n = n;
        this->m = m;
        this->l = l;
        Super::resize(n*m*l, t); 
      }

      /** \brief Returns the \a n -size of the array. */
      size_type size_n() const { return n; }
      /** \brief Returns the \a n -size of the array. */
      size_type size_x() const { return n; }

      /** \brief Returns the \a m -size  of the array. */
      size_type size_m() const { return m; }
      /** \brief Returns the \a m -size  of the array. */
      size_type size_y() const { return m; }

      /** \brief Returns the \a m -size  of the array. */
      size_type size_l() const { return l; }
      /** \brief Returns the \a m -size  of the array. */
      size_type size_z() const { return l; }

      /** \brief \c true if the \c Array3D's size is 0 x 0 x 0. */
      bool empty() const { return n == 0 && m == 0 && l == 0; }

      /** \brief Returns the element at position \p i x \p j x \p k.
	  \attention No range checking. */
      reference operator() (size_type i, size_type j, size_type k)
      { return Super::operator[] (i + j*n + k*n*m); }

      /** \brief Returns the element at position \p i x \p j x \p k.
	  \attention No range checking. */
      const_reference operator() (size_type i, size_type j, size_type k) const
      { return Super::operator[] (i + j*n + k*n*m); }

      /** \brief Returns the element at position \p pos. 

	  \p pos must be a vector where the expressions \c pos[0], \c pos[1]
      and \c pos[2] return the indices.
	  \attention No range checking. */ 
      template< class Vector >
      reference operator() (const Vector& pos)
      { return operator() (pos[0], pos[1], pos[2]); }

      /** \brief Returns the element at position \p pos. 

	  \p pos must be a vector where the expressions \c pos[0], \c pos[1] and
      \c pos[2] return the indices.
	  \attention No range checking. */ 
      template< class Vector >
      const_reference operator() (const Vector& pos) const
      { return operator() (pos[0], pos[1], pos[2]); }

      /** \brief Returns the element at position \p i x \p j x \p k. 
	  \exception out_of_range \p i is out of range of the Array3D */
      reference at(size_type i, size_type j,  size_type k) {
        if (i >= n) throw std::out_of_range("Index i out of bounds.");
        if (j >= m) throw std::out_of_range("Index j out of bounds.");
        if (k >= l) throw std::out_of_range("Index k out of bounds.");
        return operator() (i,j,k);
      }

      /** \brief Returns the element at position \p i x \p j x \p k. 
	  \exception out_of_range \p i is out of range of the Array3D */
      const_reference at(size_type i, size_type j, size_type k) const {
        if (i >= n) throw std::out_of_range("Index i out of bounds.");
        if (j >= m) throw std::out_of_range("Index j out of bounds.");
        if (k >= l) throw std::out_of_range("Index k out of bounds.");
        return operator() (i,j,k); 
      }
    
      /** \brief Returns the element at position \p pos. 

	  \p pos must be a vector where the expressions \c pos[0], \c pos[1]
      and \c pos[2] return the first and second element.

	  \exception out_of_range \p i is out of range of the Array3D */ 
      template< class Vector >
      reference at(const Vector& pos)
      { return at(pos[0], pos[1], pos[2]); }

      /** \brief Returns the element at position \p pos. 

	  \p pos must be a vector where the expressions \c pos[0], \c pos[1]
      and \c pos[2] return the first and second element.

	  \exception out_of_range \p i is out of range of the Array3D */ 
      template< class Vector >
      const_reference at(const Vector& pos) const
      { return at(pos[0], pos[1], pos[2]); }

      /** \brief Returns an iterator pointing to the beginning of the
	  Array3D. */
      iterator begin() { return Super::begin(); }

      /** \brief Returns a const_iterator pointing to the beginning of the
	  Array3D. */
      const_iterator begin() const { return Super::begin(); }

      /** \brief Returns an iterator pointing to the end of the
	  Array3D. */
      iterator end() { return Super::end(); }

      /** \brief Returns a const_iterator pointing to the end of the
	  Array3D. */
      const_iterator end() const { return Super::end(); }

      /** \brief Assigns t2 to each of the elements. */
      template < class T2 >
      void operator=(const T2 &t2) {
        for (iterator it = begin(); it != end(); it++) *it = t2;
      }

      /** \brief Adds t2 to each of the elements.  */
      template < class T2 >
      void operator+=(const T2 &t2) {
        for (iterator it = begin(); it != end(); it++) *it += t2;
      } 

      /** \brief Adds t2 to each of the elements.  */
      template < class T2 >
      void operator-=(const T2 &t2) {
        for (iterator it = begin(); it != end(); it++) *it -= t2;
      } 

      /** \brief Multiplies each of the elements by t2. */
      template < class T2 >
      void operator*=(const T2 &t2) 
      { for (iterator it = begin(); it != end(); it++) *it *= t2; } 

      /** \brief Divides each of the elements by t2. */
      template < class T2 >
      void operator/=(const T2 &t2) 
      { for (iterator it = begin(); it != end(); it++) *it /= t2; } 

      void operator/=(const double &x) {
        double x_rec = 1.0 / x;
        for (iterator it = begin(); it != end(); it++) *it *= x_rec; 
      } 

      /** \brief Equality operator. */
      bool operator==(const Array3D< T, exception > &a) {
        return static_cast< std::vector<T> >(*this) == 
                static_cast< std::vector<T> >(a);
      }

      /** \brief Inequality operator. */
      bool operator!=(const Array3D< T,exception > &a) {
        return static_cast< std::vector<T> >(*this) != 
                static_cast< std::vector<T> >(a);
      }

    private:
      /** \brief The \p n size of the \p n x \p m x \p l array */
      size_type n;
      /** \brief The \p m size of the \p n x \p m x \p l array */
      size_type m;
      /** \brief The \p l size of the \p n x \p m x \p l array */
      size_type l;
    };

    template < class T >
    class Array3D< T, enlarge > : public Array3D< T, exception > {
      /** \brief The super class of the array (i.e. vector<T>). */
      typedef Array3D< T, exception > Super;

    public:
      /** \brief The type of object, \c T, stored in the Array1D. */
      typedef T value_type;
      /** \brief Reference to to a \c T object. */
      typedef T& reference;
      /** \brief Const reference to to a \c T object. */
      typedef const T& const_reference;
      /** \brief An unsigned integer type. */
      typedef typename Super::size_type size_type;
      /** \brief Iterator used to iterate through an Array1D. */
      typedef typename Super::iterator iterator;
      /** \brief Const-iterator used to iterate through an Array1D. */
      typedef typename Super::const_iterator const_iterator;

      using Super::size_n;
      using Super::size_m;
      using Super::size_l;

      /** \brief Default constructor. */
      Array3D(const T &prototype = T()) { init(0, 0, 0, prototype); }

      /** \brief Creates a one-dimensional array of size \a n with \a n
	  copies of \a prototype.
	  Also sets the \p prototype. */
      Array3D(size_type n, size_type m, size_type l, const T &prototype = T()) 
      { init(n, m, l, prototype); }

      /** \brief Creates a one-dimensional array of size \a n with \a n
	  copies of \a prototype.
	  Also sets the \p prototype. */
      void init(size_type n, size_type m, size_type l, const T &prototype = T())
      { set_prototype(prototype); Super::init(n, m, l, prototype); }

      /** \brief Sets the \p prototype used when the Array3D<T,enlarge>
	  is automatically resized. */
      void set_prototype(const T &prototype = T()) 
      { this->prototype = prototype; }

      /** \brief Returns the element at position \p i x \p j x \p k. 

	  When the element doesn't exist, the Array3D is automatically
	  extended so that the element exists. The new elements are
	  initialised with copies of the prototype.
      */
      
      
      reference at(size_type i, size_type j, size_type k) { 
        if (i >= size_n() || j >= size_m() || k >= size_l()) {
          Super old = *this;
          size_type new_n = size_n();
          size_type new_m = size_m();
          size_type new_l = size_l();
          if (i >= new_n) new_n = i+1;
          if (j >= new_m) new_m = j+1;
          if (k >= new_l) new_l = k+1;
          Super::init(new_n, new_m, new_l, prototype);

          // copy the old Array3D
          for (size_type ii = 0; ii < old.size_n(); ii++)
            for (size_type jj = 0; jj < old.size_m(); jj++)
              for (size_type kk = 0; kk < old.size_l(); kk++)
                (*this)(ii,jj,kk) = old(ii,jj,kk);

          // init the new elements with the prototype
          for (size_type ii = old.size_n(); ii < size_n();  ii++)
            for (size_type jj = 0; jj < size_m(); jj++)
              for (size_type kk = 0; kk < size_l(); kk++)
                (*this)(ii,jj,kk) = prototype;

          for (size_type ii = 0; ii < old.size_n();  ii++)
            for (size_type jj = old.size_m(); jj < size_m(); jj++)
              for (size_type kk = 0; kk < size_l(); kk++)
                (*this)(ii,jj,kk) = prototype;

          for (size_type ii = 0; ii < old.size_n();  ii++)
            for (size_type jj = 0; jj < size_m(); jj++)
              for (size_type kk = old.size_l(); kk < size_l(); kk++)
                (*this)(ii,jj,kk) = prototype;
        }
        return Super::at(i,j,k); 
      }


      /** \brief Returns the element at position \p pos. 

	  \p pos must be a vector where the expressions \c pos[0], \c pos[1]
      and \c pos[2] return the first and second element.
    
	  When the element doesn't exist, the Array3D is automatically
	  extended so that the element exists. The new elements are
	  initialised with copies of the prototype. 
      */
      template< class Vector >
      reference at(const Vector& pos)
      { return at(pos[0], pos[1], pos[2]); }

    private:
      value_type prototype;
    };
  }
}  
#endif
