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

#include "InterpolationLinear.hpp"

#define MAXLINE 1024

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(InterpolationLinear::theLogger, "InterpolationLinear");

    InterpolationLinear::InterpolationLinear() {
      // NULLify pointers in case of errors
      radius = NULL;
      energy = NULL;
      force  = NULL;
    }

    InterpolationLinear::~InterpolationLinear() {
      LOG4ESPP_INFO(theLogger, "~InterpolcationLinear");
    }
    
    
    
    // public functions

    // this is called by the created InterpolationLinear object
    void InterpolationLinear::readRaw(mpi::communicator comm, const char* file) {
        int root = 0;  // control processor
        
        // read number of lines
        if (comm.rank() == root) N = readFile(file, true);  // dummy read
        mpi::broadcast(comm, N, root);
        if (N < 2) throw std::runtime_error("illegal file for tabulated potential");
        
        radius = new real[N];
        energy = new real[N];
        force  = new real[N];
        ae = new real[N];
        be = new real[N];
        af = new real[N];
        bf = new real[N];
        
        // read entries
        if (comm.rank() == root) readFile(file, false);
        mpi::broadcast(comm, radius, N, root);
        mpi::broadcast(comm, energy, N, root);
        mpi::broadcast(comm, force, N, root);
        
        delta = radius[1] - radius[0];
        if (delta <= 0.0) LOG4ESPP_ERROR(theLogger, "illegal distance for entries")
        
        // check if delta is the same between all entries
        for (int i = 2; i < N; i++) {
            real r = radius[i] - radius[i-1];
            if (fabs(r - delta) > 0.0001) LOG4ESPP_ERROR(theLogger, "delta " << r << " not same as " << delta);
        }
        
        inner = radius[0];
        outer = radius[N-1];
        delta  = (outer - inner) / (N-1);
        invdelta = 1.0 / delta;
        
        LOG4ESPP_INFO(theLogger, "tab file has range " << inner << " - "
                                << outer << ", delta = " << delta);
        
        spline(radius, energy, N, ae, be);
        spline(radius, force,  N, af, bf);
      
    }// readRaw


    // private functions


    // read file, store values in r, e, f, and return number of read lines (called by readRaw())
    int InterpolationLinear::readFile(const char* file, bool dummy) {
     
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
    }// readFile


    // Linear spline method adapted from The VOTCA project (http://www.votca.org)
    void InterpolationLinear::spline(const real* x, const real* y, int N,
                                    real* a, real* b) {
        
        // calculate a,b for all intervals 0..(N-2), where interval
        // [x(i),x(i+1)] shall have number i (this means that the last interval
        // has number N-2)
        for (int i = 0; i < N-1; i++) {
            a[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
            b[i] = y[i]-a[i]*x[i];
        }
        
        
    }// spline



  }//ns interaction
}//ns espressopp
