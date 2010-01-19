#ifndef _ESUTIL_MATH_HPP
#define _ESUTIL_MATH_HPP

namespace espresso {
  namespace esutil {
    real getDist(real dist[3], const real p[3], const real q[3]) {
      real dist[0] = p[0]-q[0];
      real dist[1] = p[1]-q[1];
      real dist[2] = p[2]-q[2];
      return dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
    }
  }
}

#endif
