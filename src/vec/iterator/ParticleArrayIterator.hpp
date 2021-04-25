/*
  Copyright (C) 2021
      Max Planck Institute for Polymer Research & JGU Mainz

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
#ifndef VEC_ITERATOR_PARTICLEARRAYITERATOR_HPP
#define VEC_ITERATOR_PARTICLEARRAYITERATOR_HPP

#include "vec/ParticleArray.hpp"

namespace espressopp { namespace vec {
  namespace iterator {
    /**
       Iterates all Particles in a list of cells. This is a Python-like,
       self-contained iterator: isValid() tells whether there are more
       particles to come.
    */
    class ParticleArrayIterator
    {
    public:
      ParticleArrayIterator(ParticleArray &pa, bool realCellsOnly=false);
      ParticleArrayIterator &operator++();

      bool isValid() const;
      bool isDone() const;

      inline size_t index() const {
        return pit;
      }

      inline size_t cellId() const {
        if(realCellsOnly)
          return cit;
        else
          return particles.realCells()[cit];
      }

    private:
      void findNonemptyCell();
      void renewParticleIterator();
      ParticleArray &particles;
      const bool realCellsOnly;

      size_t pit;
      size_t pend;
      size_t cit;
      size_t cend;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    inline
    ParticleArrayIterator::
    ParticleArrayIterator(ParticleArray &pa, bool realCellsOnly)
      : particles(pa),
        cit(0),
        cend(realCellsOnly ? particles.realCells().size() : particles.numCells()),
        realCellsOnly(realCellsOnly)
    {
      if(cit==cend) return;
      renewParticleIterator();
      findNonemptyCell();
    }

    inline ParticleArrayIterator&
    ParticleArrayIterator::
    operator++()
    {
      ++pit;
      findNonemptyCell();
      return *this;
    }

    inline bool
    ParticleArrayIterator::
    isValid() const
    {
      return (cit < cend);
    }

    inline bool
    ParticleArrayIterator::
    isDone() const
    {
      return !isValid();
    }

    inline void
    ParticleArrayIterator::
    findNonemptyCell()
    {
      while (pit==pend) {
        ++cit;
        if(cit==cend) break;
        renewParticleIterator();
      }
    }

    inline void
    ParticleArrayIterator::
    renewParticleIterator()
    {
      const size_t cell = realCellsOnly ? particles.realCells()[cit] : cit;
      pit  = particles.cellRange()[cell];
      pend = pit + particles.sizes()[cell];
    }
  }
}}
#endif
