#include "System.hpp"

using namespace espresso;

void 
System::foldCoordinate(real pos[3], integer imageBox[3], integer dir) {
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

void 
System::foldPosition(real pos[3], integer imageBox[3]) {
  for(int i = 0; i < 3; ++i)
    foldCoordinate(pos, imageBox, i);
}

void 
System::unfoldPosition(real pos[3], integer imageBox[3]) {
  for(int i = 0; i < 3; ++i) {
    pos[i] = pos[i] + imageBox[i]*getBoxL(i);    
    imageBox[i] = 0;
  }
}
