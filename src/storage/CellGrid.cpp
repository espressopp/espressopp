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

#include <cmath>

#include "log4espp.hpp"

#include "CellGrid.hpp"
#include "Int3D.hpp"
#include "Real3D.hpp"

using namespace std;
using namespace espressopp;

LOG4ESPP_LOGGER(CellGrid::logger, "DomainDecomposition.CellGrid");

CellGridIllegal::CellGridIllegal()
  : std::runtime_error("cell grid dimensions have to be positive") {}

CellGrid::CellGrid(const Int3D& _size, const real _myLeft[3], const real _myRight[3], int _frame):
        Grid(_size[0]+2*_frame, _size[1]+2*_frame, _size[2]+2*_frame),
        frame(_frame), extraSize(2*frame)
{
  if (_size[0] <= 0 || _size[1] <= 0 || _size[2] <= 0) {
    throw CellGridIllegal();
  }

  for(int i = 0; i < 3; ++i) {
    myLeft[i] = _myLeft[i];
    myRight[i] = _myRight[i];

    cellSize[i]           = (myRight[i] - myLeft[i])/real(_size[i]);
    invCellSize[i]        = 1.0 / cellSize[i];
  }

  smallestCellDiameter = std::min(std::min(cellSize[0], cellSize[1]), cellSize[2]);
}

longint CellGrid::getNumberOfInnerCells() const
{
  longint res = getGridSize(0);
  for (int i = 1; i < 3; ++i) {
    res *= getGridSize(i);
  }
  return res;
}


longint CellGrid::mapPositionToCell(const Real3D& pos) const {
  int cpos[3];

  for(int i = 0; i < 3; ++i) {
    real lpos = pos[i] - myLeft[i];
    if (lpos <= 0) cpos[i] = 0;
    else cpos[i] = static_cast<int>(lpos*invCellSize[i]) + frame;
  }
  return mapPositionToIndex(cpos);
}


longint CellGrid::mapPositionToCellClipped(const Real3D& pos) const
{
  int cpos[3];

  for(int i = 0; i < 3; ++i) {
    real lpos = pos[i] - myLeft[i];

    cpos[i] = static_cast<int>(lpos*invCellSize[i]) + frame;

    if (cpos[i] < frame) {
      cpos[i] = frame;
    }
    else if (cpos[i] >= getFrameGridSize(i) - frame) {
      cpos[i] = getFrameGridSize(i) - frame - 1;
    }
  }
  return mapPositionToIndex(cpos);  
}

longint CellGrid::mapPositionToCellChecked(const Real3D& pos) const
{
  int cpos[3];

  for(int i = 0; i < 3; ++i) {
    real lpos = pos[i] - myLeft[i];

    cpos[i] = static_cast<int>(lpos*invCellSize[i]) + frame;
    
    /* particles outside our box. Still take them if
       VERY close or nonperiodic boundary */
    if (cpos[i] < frame) {
      if (lpos >= -ROUND_ERROR_PREC) {
	cpos[i] = frame;
      }
      else {
	return noCell;
      }
    }
    else if (cpos[i] >= getFrameGridSize(i) - frame) {
      if (pos[i] <= myRight[i] + ROUND_ERROR_PREC) {
	cpos[i] = getFrameGridSize(i) - frame - 1;
      }
      else {
	return noCell;
      }
    }
  }
  return mapPositionToIndex(cpos);  
}

bool CellGrid::mapPositionToCellCheckedAndClipped(longint &cell, const Real3D& pos) const
{
  int cpos[3];
  bool outside = false;

  for(int i = 0; i < 3; ++i) {
    real lpos = pos[i] - myLeft[i];

    cpos[i] = static_cast<int>(lpos*invCellSize[i]) + frame;
    
    /* particles outside our box. Still take them if
       VERY close or nonperiodic boundary */
    if (cpos[i] < frame) {
      cpos[i] = frame;

      if (lpos < -ROUND_ERROR_PREC) {
	outside = true;
      }
    }
    else if (cpos[i] >= getFrameGridSize(i) - frame) {
      cpos[i] = getFrameGridSize(i) - frame - 1;

      if (pos[i] > myRight[i] + ROUND_ERROR_PREC) {
	outside = true;
      }
    }
  }

  cell = mapPositionToIndex(cpos);

  LOG4ESPP_TRACE(logger, "mapping position ("
		 << pos[0] << ", "
		 << pos[1] << ", "
		 << pos[2] << ") to grid position ("
		 << cpos[0] << ", "
		 << cpos[1] << ", "
		 << cpos[2] << ") <-> "
		 << cell);

  return outside;
}
