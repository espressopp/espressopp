#ifndef _TYPES_HPP
#define _TYPES_HPP

// general type definition to choose between single and double precsion 

// typedef float real;  // single precision

typedef double real;  // double precision

template<class real>
real sign(real _r) {
    return  (_r > 0) ? 1 : -1;
}

/** Class for a vector of three real values */

class Real3D {
   private:
      real x, y, z;  //*< elements of the 3D vector 

   public:
      Real3D(real _x, real _y, real _z) { x = _x; y = _y; z = _z; }
      Real3D(real val = 0.0) { x = val; y = val; z = val; }
      Real3D operator+ (const Real3D &b) const {
         return Real3D(x + b.x, y + b.y, z + b.z);
      }
 
      Real3D& operator+= (const Real3D &b) {
         x += b.x;
         y += b.y;
         z += b.z;
         return *this;
      }
 
      Real3D& operator-= (const Real3D &b) {
         x -= b.x;
         y -= b.y;
         z -= b.z;
         return *this;
      }
 
      Real3D operator- (const Real3D &b) const {
         return Real3D(x - b.x, y - b.y, z - b.z);
      }
 
      Real3D operator* (const Real3D &b) const {
         return Real3D(x * b.x, y * b.y, z * b.z);
      }
 
      Real3D operator* (real s) const {
         return Real3D(x * s, y * s, z * s);
      }

      real sqr() const { return x*x + y*y + z*z; }

      real getX() const { return x; }

      real getY() const { return y; }

      real getZ() const { return z; }
};

inline Real3D operator*(real s, const Real3D &v) { return v * s; }

//const void* nullptr = static_cast<void*>(0);

#endif 
