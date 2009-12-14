#include <cmath>

#define LOG4ESPP_LEVEL_DEBUG
#include "log4espp.hpp"

#include "CellGrid.hpp"

using namespace std;
using namespace espresso;

LOG4ESPP_LOGGER(CellGrid::logger, "DomainDecomposition.CellGrid");

CellGridIllegal::CellGridIllegal()
  : std::runtime_error("cell grid dimensions have to be positive") {}

CellGrid::CellGrid(const int _size[3],
		   const real _myLeft[3],
		   const real _myRight[3],
		   int _frame)
  : Grid(_size[0]+2*_frame, _size[1]+2*_frame, _size[2]+2*_frame),
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

longint CellGrid::mapPositionToCellClipped(const real pos[3]) const
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

longint CellGrid::mapPositionToCellChecked(const real pos[3]) const
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

bool CellGrid::mapPositionToCellCheckedAndClipped(longint &cell, const real pos[3]) const
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
