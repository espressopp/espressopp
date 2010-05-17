#include "InterpolationTable.hpp"

namespace espresso {
  namespace interaction {

    InterpolationTable::InterpolationTable()
    {
      radius = NULL;
      energy = NULL;
      force  = NULL;

      energy2 = NULL;
      force2  = NULL;
    }

    InterpolationTable::~InterpolationTable()
    {
      // free the allocated memory

      delete [] radius;
      delete [] energy;
      delete [] force;

      delete [] energy2;
      delete [] force2;
    }

    /** Spline read-in values. */

    void InterpolationTable::spline(const double *x, const double *y, int n,
                                    double yp1, double ypn, double *y2)
    {
      double *u = new double[n];

      if (yp1 > 0.99e30) {
        y2[0] = 0.0;
        u[0]  = 0.0;
      } else {
        y2[0] = -0.5;
        u[0]  = (3.0/(x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
      }
    
      for (int i = 1; i < n-1; i++) {

        double sig = (x[i] - x[i-1]) / (x[i+1]-x[i-1]);
        double p   = sig * y2[i-1] + 2.0;
    
        y2[i] = (sig-1.0) / p;
        u[i]  = (y[i+1] - y[i]) / (x[i+1] - x[i]) -
                (y[i] - y[i-1]) / (x[i] - x[i-1]);
        u[i]  = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
      }
    
      double qn, un;
    
      if (ypn > 0.99e30) {
        qn = 0.0;
        un = 0.0;
      } else {
        qn = 0.5;
        un = (3.0 / (x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
      }
    
      y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);
        
      for (int k = n-2; k >= 0; k--) {
        y2[k] = y2[k] * y2[k+1] + u[k];
      }
    
      delete [] u;
    }

    void InterpolationTable::read(const char* file)
    {
      inner = 0.6;
      real outer = 3.0;
 
      N = 201;

      int nbins = N - 1;  // number of intervals is number of points - 1

      delta    = (outer - inner) / nbins;
      invdelta = 1.0 / delta;
      deltasq6 = delta * delta / 6.0;

      // allocate the arrays

      radius = new real[N];
      energy = new real[N];
      force  = new real[N];

      // start with a test file 

      for (int i=0; i < N; i++) {
         real r = inner + i * delta;
         real frac2 = 1.0 / (r * r);
         real frac6 = frac2 * frac2 * frac2;
         radius[i] = r;
         force[i] = 48.0 * (frac6 * frac6 - 0.5 * frac6);
         energy[i] = 4.0 * (frac6 * frac6 - frac6);
      }

      energy2 = new real[N];
      force2  = new real[N];

      // the derivatives of energy are known by the forces

      real yp1 = -force[0];
      real ypN = -force[N-1];

      spline(radius, energy, N, yp1, ypN, energy2);

      // the derivatives of forces are approximated

      yp1 = (force[1] - force[0]) / (radius[1] - radius[0]);
      ypN = (force[N-1] - force[N-2]) / (radius[N-1] - radius[N-2]);

      spline(radius, force, N, yp1, ypN, force2);
    }

    real InterpolationTable::splineInterpolation(real r, const double* fn, const double* fn2) const
    {
      int index = static_cast<int>((r - inner) * invdelta);

      // printf("r = %10.3f maps to index %d, N = %d\n", r, index, N);

      assert(index >=0 && index < N-1);

      real b = (r - radius[index]) * invdelta;
      real a = 1.0 - b;

      real f = a * fn[index] + b * fn[index+1] +
               ((a*a*a-a)*fn2[index] + (b*b*b-b)*fn2[index+1]) *
               deltasq6;

      // printf("interpoloate %f, a = %f, b = %f, fn = %f - %f, fn2 = %f - %f\n", 
      //        r, a, b, fn[index], fn[index+1], fn2[index], fn2[index+1]);

      return f;
    }

    real InterpolationTable::getEnergy(real r) const
    {
      return splineInterpolation(r, energy, energy2);
    }

    real InterpolationTable::getForce(real r) const
    {
      return splineInterpolation(r, force, force2);
    }
  }
}
