#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "types.hpp"

namespace espresso {
  class System {
  public:
    shared_ptr< class Storage > storage;
    shared_ptr< class BC > gen_bc;

    real boxL[3];
    real invBoxL[3];

    void setBoxL(const real _boxL[3]) {
      for (int i = 0; i < 3; ++i) {
        boxL[i]    = _boxL[i];
        invBoxL[i] = 1.0/_boxL[i];
      }
    }

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
    void foldCoordinate(real pos[3], int imageBox[3], int dir);

    /** fold particle coordinates to primary simulation box.
	\param pos the position...
	\param imageBox and the box

	Both pos and image_box are I/O,
	i. e. a previously folded position will be folded correctly.
    */
    void foldPosition(real pos[3], int imageBox[3]);

    /** unfold coordinates to physical position.
	\param pos the position...
	\param imageBox and the box
      
	Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
	afterwards.
    */
    void unfoldPosition(real pos[3], int imageBox[3]);

  private:
  };
}
#endif
