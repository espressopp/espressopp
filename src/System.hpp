#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "types.hpp"

struct System {
  real boxL[3];
  real invBoxL[3];

  const real *getBoxL()    const { return boxL; }
  const real *getInvBoxL() const { return invBoxL; }
  real getBoxL(int i)      const { return boxL[i]; }
  real getInvBoxL(int i)   const { return invBoxL[i]; }

  /** fold a coordinate to primary simulation box.
      \param pos         the position...
      \param imageBox    and the box
      \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

      Both pos and image_box are I/O,
      i. e. a previously folded position will be folded correctly.
  */
  void foldCoordinate(real pos[3], integer imageBox[3], integer dir) {
    integer tmp = (integer)floor(pos[dir]*getInvBoxL(dir));

    imageBox[dir] += tmp;
    pos[dir] -= tmp*getBoxL(dir);    

    if(pos[dir] < 0 || pos[dir] >= getBoxL(dir)) {
      /* slow but safe */
      if (fabs(pos[dir]*getInvBoxL(dir)) >= INT_MAX/2) {
# warning ERRORHANDLING MISSING
#if 0
	char *errtext = runtime_error(128 + TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{086 particle coordinate out of range, pos = %f, image box = %d} ", pos[dir], image_box[dir]);
#endif
	imageBox[dir] = 0;
	pos[dir] = 0;
      }
    }
  }

  /** fold particle coordinates to primary simulation box.
      \param pos the position...
      \param imageBox and the box

      Both pos and image_box are I/O,
      i. e. a previously folded position will be folded correctly.
  */
  void foldPosition(real pos[3], integer imageBox[3])
  {
    for(int i = 0; i < 3; ++i)
      foldCoordinate(pos, imageBox, i);
  }

  /** unfold coordinates to physical position.
      \param pos the position...
      \param imageBox and the box
      
      Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
      afterwards.
  */
  void unfoldPosition(real pos[3], integer imageBox[3])
  {
    for(int i = 0; i < 3; ++i) {
      pos[i] = pos[i] + imageBox[i]*getBoxL(i);    
      imageBox[i] = 0;
    }
  }
};

#endif
