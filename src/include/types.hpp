#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <cstddef>    // needed for size_t
#include "esutil/Vector3D.hpp"

// general type definition to choose between single and double precsion 

// typedef float real;  // single precision

typedef double real;  // double precision

template<class real>
real sign(real _r) {
    return  (_r > 0) ? 1 : -1;
}

/// make Real3D globally accessible for convenience
typedef esutil::Vector3D<real> Real3D;

//const void* nullptr = static_cast<void*>(0);

#endif 
