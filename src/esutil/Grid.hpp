#ifndef _ESUTIL_GRID_HPP
#define _ESUTIL_GRID_HPP

#include "Box.hpp"

namespace espresso {
  namespace esutil {
    /** Represents a box subdivided into an equally spaced grid.
     */
    template<class _VectorClass>
    class Grid {
    public:
      typedef _VectorClass VectorClass;
      typedef typename VectorClass::ValueType ValueType;
      typedef SmallVector<size_t, VectorClass::dimension> CellIdentifier;
      
      Grid(const Box<VectorClass> &_domain, size_t _cells[VectorClass::dimension])
        : domain(_domain)
      {
        for(size_t i = 0; i < VectorClass::dimension; ++i) {
          setNumCells(i, _cells[i]);
        }
      }

      Grid(const Box<VectorClass> &_domain)
        : domain(_domain)
      {
        for(size_t i = 0; i < VectorClass::dimension; ++i) {
          setNumCells(i, 1);
        }
      }

      /// set the domain that the grid is spanning up
      void setDomain(const Box<VectorClass> &_domain) {
        domain = _domain;
        // recalculate cellsizes
        for(size_t i = 0; i < VectorClass::dimension; ++i) {
          setNumCells(i, cells[i]);
        }
      }

      /// set the number of cells in one dimension
      void setNumCells(size_t i, size_t s) {
        cells[i] = s;
        cellsizes[i] = domain.getExtend()[i]/s;
        cellsizes_i[i] = s/domain.getExtend()[i];
        updateMaxId();
      }

      /// get the domain that the grid is spanning up
      const Box<VectorClass> &getDomain()  const { return domain; }
      /// get the number of cells in each dimension
      const size_t    *getNumCells() const { return cells; }
      /// get the dimensions of a cell
      const ValueType *getCellSize() const { return cellsizes; }
      /// get the inverse of the cell sizes
      const ValueType *getInverseCellSize() const { return cellsizes_i; }

      /// return the cell that contains a given position, or Grid::noCell
      CellIdentifier locate(const VectorClass &vec) {
        CellIdentifier result;

        VectorClass rel = vec - domain.getLeft();
        for (size_t i = 0; i < VectorClass::dimension; ++i) {
          ValueType ind = rel[i]*cellsizes_i[i];
          if (ind < 0 || ind >= cells[i]) {
            return notACell;
          }
          result[i] = ind;
        }
        return result;
      }

      /// turn a cell identifier into a linear index between 0 and getMaxLinearId()
      size_t linearize(const CellIdentifier &id) {
        size_t res = id[0];
        for (size_t d = 1; d < VectorClass::dimension; ++d) {
          res = res*cells[d] + id[d];
        }
        return res;
      }
      
      size_t getMaxLinearId() const { return maxID; }

      /// not a cell, returned by various algorithms in case of failure
      static const CellIdentifier notACell;

    protected:
      void updateMaxId() {
        maxID = 1;
        for (size_t d = 0; d < VectorClass::dimension; ++d) {
          maxID *= cells[d];
        }
      }

      Box<VectorClass> domain;
      size_t cells[VectorClass::dimension];

      size_t maxID;

      typename VectorClass::ValueType cellsizes[VectorClass::dimension];
      typename VectorClass::ValueType cellsizes_i[VectorClass::dimension];
    };

    template<class VectorClass>
    const typename Grid<VectorClass>::CellIdentifier Grid<VectorClass>::notACell(static_cast<size_t>(-1));
  }
}

#endif
