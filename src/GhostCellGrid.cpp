#include "GhostCellGrid.hpp"
#include <cmath>

using namespace std;

CellGridIllegal::CellGridIllegal()
  : std::runtime_error("cell grid dimensions have to be positive") {}

GhostCellGrid::GhostCellGrid(const integer _size[3],
			     const real _myLeft[3],
			     const real _myRight[3],
			     integer _frame)
  : Grid(_size[0]+2*_frame, _size[1]+2*_frame, _size[2]+2*_frame),
    frame(_frame), extraSize(2*frame)
{
  if (_size[0] <= 0 || _size[1] <= 0 || _size[2] <= 0) {
    throw CellGridIllegal();
  }

  for(integer i = 0; i < 3; ++i) {
    myLeft[i] = _myLeft[i];
    myRight[i] = _myRight[i];

    cellSize[i]           = (myRight[i] - myLeft[i])/real(_size[i]);
    invCellSize[i]        = 1.0 / cellSize[i];
  }

  smallestCellDiameter = std::min(std::min(cellSize[0], cellSize[1]), cellSize[2]);
}

integer GhostCellGrid::getNumberOfInnerCells() const
{
  integer res = getGridSize(0);
  for (integer i = 1; i < 3; ++i) {
    res *= getGridSize(i);
  }
  return res;
}

integer GhostCellGrid::mapPositionToCellClipping(const real pos[3]) const
{
  integer cpos[3];

  for(integer i = 0; i < 3; ++i) {
    real lpos = pos[i] - myLeft[i];

    cpos[i] = static_cast<integer>(lpos*invCellSize[i]) + frame;

    if (cpos[i] < frame) {
      cpos[i] = frame;
    }
    else if (cpos[i] >= getGridSize(i) - frame) {
      cpos[i] = getGridSize(i) - frame - 1;
    }
  }
  return getLinearIndex(cpos);  
}

integer GhostCellGrid::mapPositionToCellChecked(const real pos[3]) const
{
  integer cpos[3];

  for(integer i = 0; i < 3; ++i) {
    real lpos = pos[i] - myLeft[i];

    cpos[i] = static_cast<integer>(lpos*invCellSize[i]) + frame;
    
    /* particles outside our box. Still take them if
       VERY close or nonperiodic boundary */
    if (cpos[i] < 1) {
      if (lpos > -ROUND_ERROR_PREC) {
	cpos[i] = 1;
      }
      else {
	return noCell;
      }
    }
    else if (cpos[i] >= getGridSize(i) - frame) {
      if (pos[i] < myRight[i] + ROUND_ERROR_PREC) {
	cpos[i] = getGridSize(i) - frame - 1;
      }
      else {
	return noCell;
      }
    }
  }
  return getLinearIndex(cpos);  
}
