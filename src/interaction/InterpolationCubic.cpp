/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "InterpolationCubic.hpp"

#define MAXLINE 1024

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(InterpolationCubic::theLogger, "InterpolationCubic");

    InterpolationCubic::InterpolationCubic() {
      // NULLify pointers in case of errors
      radius = NULL;
      energy = NULL;
      force  = NULL;
      
      energy2 = NULL;
      force2  = NULL;
      
      allocated = false;
    }

    InterpolationCubic::~InterpolationCubic() {
      LOG4ESPP_INFO(theLogger, "~InterpolationCubic");
    }
    
    
    
    // public functions

    // this is called by the created InterpolationCubic object
    void InterpolationCubic::readRaw(mpi::communicator comm, const char* file) {
      int root = 0;  // control processor

      if (comm.rank() == root) N = readFile(file, true);  // dummy read

      mpi::broadcast(comm, N, root);

      if (N < 2) {
         throw std::runtime_error("illegal file for tabulated potential");
      }

      radius = new real[N];
      energy = new real[N];
      force  = new real[N];

      if (comm.rank() == root) readFile(file, false); // read entries

      mpi::broadcast(comm, radius, N, root);
      mpi::broadcast(comm, energy, N, root);
      mpi::broadcast(comm, force, N, root);

      delta = radius[1] - radius[0];

      if (delta <= 0.0) {
         LOG4ESPP_ERROR(theLogger, "illegal distance for entries")
      }

      // check of delta is the same between all entries
      for (int i = 2; i < N; i++) {
         real r = radius[i] - radius[i-1];
         if (fabs(r - delta) > 0.0001) {
            LOG4ESPP_ERROR(theLogger, "delta " << r << " not same as " << delta);
         }
      }

      int nbins = N-1;  // number of intervals is number of points - 1
      inner = radius[0];
      outer = radius[nbins]; 
      delta = (outer - inner) / nbins;

      LOG4ESPP_INFO(theLogger, "tab file has range " << inner << " - "
                              << outer << ", delta = " << delta);

      invdelta = 1.0 / delta;
      deltasq6 = delta * delta / 6.0;

      energy2 = new real[N];
      force2  = new real[N];

      // the derivatives of energy are known by the forces
      real yp1 = -force[0];
      real ypN = -force[N-1];
      
      // the derivatives of energies are approximated (if one uses tables not for energy/force tables)
      //real yp1 = (energy[1] - energy[0]) / (radius[1] - radius[0]);
      //real ypN = (energy[N-1] -energy [N-2]) / (radius[N-1] - radius[N-2]);

      spline(radius, energy, N, yp1, ypN, energy2);

      // the derivatives of forces are approximated
      yp1 = (force[1] - force[0]) / (radius[1] - radius[0]);
      ypN = (force[N-1] - force[N-2]) / (radius[N-1] - radius[N-2]);

      spline(radius, force, N, yp1, ypN, force2);
    }// read



    // private functions


    // read file, store values in r, e, f, and return number of read lines
    // (called by read())
    int InterpolationCubic::readFile(const char* file, bool dummy) {
     
      char line[MAXLINE];
      real r, e, f;
      FILE *fp = fopen(file, "r");
     
      if (fp == NULL) {
        LOG4ESPP_ERROR(theLogger, "could not open file " << file);
        return 0;
      }
     
      int N = 0;
     
      while (1) {
        if (fgets(line, MAXLINE, fp) == NULL) break;
     
        int k = sscanf(line, "%lg %lg %lg", &r, &e, &f);
        if (k < 3) continue;    // do not accept this line
     
        if (!dummy) {
          radius[N] = r;
          energy[N] = e;
           force[N] = f;
        }
     
        N++;   // increment counter for table entries
      }
     
      fclose(fp);
     
      if (dummy) {
        LOG4ESPP_INFO(theLogger, "found " << N << " table entries file " << file);
      } else {
        LOG4ESPP_INFO(theLogger, "read " << N << " table entries in file " << file);
      }
     
      return N;
    }// readfile



    /** Spline read-in values. */
    void InterpolationCubic::spline(const real *x, const real *y, int n,
                                    real yp1, real ypn, real *y2) {
      real *u = new real[n];
     
      if (yp1 > 0.99e30) {
        y2[0] = 0.0;
        u[0]  = 0.0;
      }
      else {
        y2[0] = -0.5;
        u[0]  = (3.0/(x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
      }
     
      for (int i = 1; i < n-1; i++) {
     
        real sig = (x[i] - x[i-1]) / (x[i+1]-x[i-1]);
        real p   = sig * y2[i-1] + 2.0;
     
        y2[i] = (sig-1.0) / p;
        u[i]  = (y[i+1] - y[i]) / (x[i+1] - x[i]) -
                (y[i] - y[i-1]) / (x[i] - x[i-1]);
        u[i]  = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
      }
     
      real qn, un;
     
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
    }// spline



  }//ns interaction
}//ns espressopp
