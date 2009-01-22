#ifndef TYPES_HPP
#define TYPES_HPP

// general type definition to choose between single and double precsion 

// typedef float real;  // single precision

typedef double real;  // double precision

/** Class for a vector of three real values */

class Real3D {

   private:
 
      real x, y, z;  //*< elements of the 3D vector 

   public:

      Real3D(real _x, real _y, real _z) { x = _x; y = _y; z = _z; }

      Real3D(real val) { x = val; y = val; z = val; }

      Real3D operator+ (Real3D b) const

      {
         return Real3D(x + b.x, y + b.y, z + b.z);
      }
 
      Real3D& operator+= (Real3D b) 

      {
         x += b.x;
         y += b.y;
         z += b.z;
         return *this;
      }
 
      Real3D operator- (Real3D b) const

      {
         return Real3D(x - b.x, y - b.y, z - b.z);
      }
 
      Real3D operator* (Real3D b) const

      {
         return Real3D(x * b.x, y * b.y, z * b.z);
      }
 
      Real3D operator* (real s) const

      {
         return Real3D(this->x * s, this->y * s, this->z * s);
      }

      real sqr() const

      { return x*x + y*y + z*z;
      }

      void addTo (real* ptr) const {
        ptr[0] += x;
        ptr[1] += y;
        ptr[2] += z;
      }

      void subFrom (real* ptr) const {
        ptr[0] -= x;
        ptr[1] -= y;
        ptr[2] -= z;
      }

      real getX() const {
        return x;
      }

      real getY() const {
        return y;
      }

      real getZ() const {
        return z;
      }
};

#endif 
