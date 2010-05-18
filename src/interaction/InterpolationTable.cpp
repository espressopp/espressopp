#include "InterpolationTable.hpp"

#define MAXLINE 1024

namespace espresso {
  namespace interaction {

    LOG4ESPP_LOGGER(InterpolationTable::theLogger, "InterpolationTable");

    InterpolationTable::InterpolationTable()
    {
      // NULLify pointers in case of errors

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

    int InterpolationTable::readFile(const char* file)
    {
      char line[MAXLINE];

      FILE *fp = fopen(file, "r");

      if (fp == NULL) {
        LOG4ESPP_ERROR(theLogger, "could not open file " << file);
        return 0;
      }

      int N = 0;

      while (1) {

        if (fgets(line, MAXLINE, fp) == NULL) break;

        real r, e, f;

        int k = sscanf(line, "%lg %lg %lg", &r, &e, &f);

        if (k < 3) continue;    // do not accept this line

        N++;   // increment counter for table entries
      }

      fclose(fp);

      LOG4ESPP_INFO(theLogger, "found " << N << " valid lines in file " << file);
  
      if (N < 2) {
        LOG4ESPP_ERROR(theLogger, "File " << file << " does not contain tabulated " <<
                       "potential, need at least 2 lines: radius energy force");
        return 0;
      }

      // allocate the arrays

      radius = new real[N];
      energy = new real[N];
      force  = new real[N];
      
      // now read again and add set the values

      N = 0;

      fp = fopen(file, "r");

      assert (fp != NULL);

      while (1) {

        if (fgets(line, MAXLINE, fp) == NULL) break;

        real r, e, f;

        int k = sscanf(line, "%lg %lg %lg", &r, &e, &f);

        if (k < 3) continue;    // do not accept this line

        radius[N] = r;
        energy[N] = e;
        force[N]  = f;

        N++;   // increment counter for table entries
      }

      LOG4ESPP_INFO(theLogger, "made " << N << " entries from file " << file);

      fclose(fp);

      return N;
    }

    void InterpolationTable::read(const char* file)
    {
      int rank = 0;

      if (rank == 0) {
 
        // only control processor reads the file

        N = readFile(file);

        // now broadcast the number of entries, take 0 as an error

        // mpi::broadcast(N);

        if (N < 1) return;

      } else {

        // mpi::broadcast(N);

        // allocate the arrays

        radius = new real[N];
        energy = new real[N];
        force  = new real[N];
      }
      
      // mpi::broadcast(N);

      // make some checks 

      real dist = radius[1] - radius[0];

      assert(dist > 0.0);

      for (int i = 2; i < N; i++) {
         real r = radius[i] - radius[i-1];
         if (abs(r - dist) > 0.0001) {
            LOG4ESPP_ERROR(theLogger, "distance " << r << " not same as " << dist);
         }
      }

      int nbins = N - 1;  // number of intervals is number of points - 1

      inner = radius[0];
      real outer = radius[N-1];
 
      delta    = (outer - inner) / nbins;

      LOG4ESPP_INFO(theLogger, "tab file has range " << inner << " - "
                              << outer << ", delta = " << delta);

      invdelta = 1.0 / delta;
      deltasq6 = delta * delta / 6.0;

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

      if (index < 0 || index >= N) {
         LOG4ESPP_ERROR(theLogger, "distance " << r << " out of range "
                        << inner << " - " << inner + (N - 1) * delta);
         printf("invdelta = %f\n", invdelta);
         exit(-1);
      }

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
