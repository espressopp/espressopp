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

// ESPP_CLASS
#ifndef _ITERATOR_CELLLISTITERATOR_HPP
#define _ITERATOR_CELLLISTITERATOR_HPP

#include "Cell.hpp"

namespace espressopp {
  namespace iterator {
    /**
       Iterates all Particles in a list of cells. This is a Python-like,
       self-contained iterator: isValid() tells whether there are more
       particles to come.
    */
    class CellListIterator {
    public:
      CellListIterator();
      CellListIterator(CellList &cl);
      CellListIterator &operator++();

      bool isValid() const;
      bool isDone() const;

      Particle &operator*() const;
      Particle *operator->() const;

      Cell &getCurrentCell() const;

    private:
      void findNonemptyCell();

      CellList::Iterator cit;
      ParticleList::Iterator pit;
    };


    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    inline 
    CellListIterator::
    CellListIterator() 
    {}
    
    inline 
    CellListIterator::
    CellListIterator(CellList &cl) : cit(cl), pit() {
      if (cit.isDone()) return;
      pit = ParticleList::Iterator((*cit)->particles);
      findNonemptyCell();
    }

    inline CellListIterator&
    CellListIterator::
    operator++() {
      ++pit;
      findNonemptyCell();
      return *this;
    }

    inline bool 
    CellListIterator::
    isValid() const 
    { return cit.isValid(); }

    inline bool 
    CellListIterator::
    isDone() const 
    { return !isValid(); }
    
    inline Particle &
    CellListIterator::
    operator*() const 
    { return *pit; }

    inline Particle*
    CellListIterator::
    operator->() const 
    { return &**this; }
    
    inline Cell&
    CellListIterator::
    getCurrentCell() const
    { return **cit; }
   
    inline void 
    CellListIterator::
    findNonemptyCell() {
      while (pit.isDone()) {
        ++cit;
        if (cit.isDone()) break;
        pit = ParticleList::Iterator((*cit)->particles);
      }
    }
  }
}
#endif
