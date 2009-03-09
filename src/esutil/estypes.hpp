#ifndef _ESUTIL_ESTYPES_HPP
#define _ESUTIL_ESTYPES_HPP

// general type definition to choose between single and double precsion 

// typedef float real;  // single precision

namespace espresso {
  namespace esutil {
    typedef double real;  // double precision
    
    template<class real>
    real sign(real _r) {
      return  (_r > 0) ? 1 : -1;
    }
  }
}

#endif
