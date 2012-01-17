// ESPP_CLASS
#ifndef _ITERATOR_CELLLISTITERATOR_HPP
#define _ITERATOR_CELLLISTITERATOR_HPP

#include "Cell.hpp"

namespace espresso {
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
