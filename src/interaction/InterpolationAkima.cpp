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

#include "InterpolationAkima.hpp"

#define MAXLINE 1024

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(InterpolationAkima::theLogger, "InterpolationAkima");

    InterpolationAkima::InterpolationAkima() {
      // NULLify pointers in case of errors
      radius = NULL;
      energy = NULL;
      force  = NULL;
    }

    InterpolationAkima::~InterpolationAkima() {
      LOG4ESPP_INFO(theLogger, "~InterpolcationAkima");
    }
    
    
    
    // public functions

    // this is called by the created InterpolationAkima object
    void InterpolationAkima::readRaw(mpi::communicator comm, const char* file) {
        int root = 0;  // control processor
        
        // read number of lines
        if (comm.rank() == root) N = readFile(file, true);  // dummy read
        mpi::broadcast(comm, N, root);
        if (N < 2) throw std::runtime_error("illegal file for tabulated potential");
        
        radius = new real[N];
        energy = new real[N];
        force  = new real[N];
		p0e = new real[N];
		p1e = new real[N];
		p2e = new real[N];
		p3e = new real[N];
		p0f = new real[N];
        p1f = new real[N];
        p2f = new real[N];
        p3f = new real[N];
        
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
        
        spline(radius, energy, N, p0e, p1e, p2e, p3e);
        spline(radius, force,  N, p0f, p1f, p2f, p3f);
      
    }// readRaw


    // private functions


    // read file, store values in r, e, f, and return number of read lines (called by readRaw())
    int InterpolationAkima::readFile(const char* file, bool dummy) {
     
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


    // Akima spline method adapted from The VOTCA project (http://www.votca.org)
    void InterpolationAkima::spline(const real* x, const real* y, int N,
                                    real* p0, real* p1, real* p2, real* p3) {
        
        real m1, m2, m3, m4, m5;
        real temp, g0, g1, g2, x1, x2, y1, y2, x4, x5, y4, y5;
        real* t = new real[N];
        
        // Akima method:
        //  estimation of two more points on each side using a degree two polyomial.
        // Resulting slopes t[0), t[1), t[N-2), t[N-1) are directly calculated

        // left side: t[0), t[1)
        temp = (x[1]-x[0])/(x[2]-x[0]);
        temp = temp*temp;
        g0 = y[0];
        g1 = ( (y[1]-y[0]) - temp*(y[2]-y[0]) ) / ( (x[1]-x[0]) - temp*(x[2]-x[0]) );
        g2 = ( (y[2]-y[0]) - g1*(x[2]-x[0]) ) / ( (x[2]-x[0])*(x[2]-x[0]) );
        x1 = x[0] - (x[2]-x[0]);
        x2 = x[1] - (x[2]-x[0]);
        y1 = g0 + g1*(x1-x[0]) + g2*(x1-x[0])*(x1-x[0]);
        y2 = g0 + g1*(x2-x[0]) + g2*(x2-x[0])*(x2-x[0]);
        m1 = (y2-y1)/(x2-x1);
        m2 = (y[0]-y2)/(x[0]-x2);
        m3 = (y[1]-y[0])/(x[1]-x[0]);
        m4 = (y[2]-y[1])/(x[2]-x[1]);
        t[0] = getSlope(m1,m2,m3,m4);
        m5 = (y[3]-y[2])/(x[3]-x[2]);
        t[1] = getSlope(m2,m3,m4,m5);

        // right side: t[N-2), t[N-1)
        temp = (x[N-2]-x[N-1])/(x[N-3]-x[N-1]);
        temp = temp*temp;
        g0 = y[N-1];
        g1 = ( (y[N-2]-y[N-1]) - temp*(y[N-3]-y[N-1]) ) / ( (x[N-2]-x[N-1]) - temp*(x[N-3]-x[N-1]) );
        g2 = ( (y[N-3]-y[N-1]) - g1*(x[N-3]-x[N-1]) ) / ( (x[N-3]-x[N-1])*(x[N-3]-x[N-1]) );
        x4 = x[N-2] + (x[N-1]-x[N-3]);
        x5 = x[N-1] + (x[N-1]-x[N-3]);
        y4 = g0 + g1*(x4-x[N-1]) + g2*(x4-x[N-1])*(x4-x[N-1]);
        y5 = g0 + g1*(x5-x[N-1]) + g2*(x5-x[N-1])*(x5-x[N-1]);
        m1 = (y[N-3]-y[N-4])/(x[N-3]-x[N-4]);
        m2 = (y[N-2]-y[N-3])/(x[N-2]-x[N-3]);
        m3 = (y[N-1]-y[N-2])/(x[N-1]-x[N-2]);
        m4 = (y4-y[N-1])/(x4-x[N-1]);
        m5 = (y5-y4)/(x5-x4);
        t[N-2] = getSlope(m1,m2,m3,m4);
        t[N-1] = getSlope(m2,m3,m4,m5);
        
        // calculate t's for all inner points [2,N-3]
        for (int i = 2; i < N-2; i++) {
            m1 = (y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
            m2 = (y[i]-y[i-1])/(x[i]-x[i-1]);
            m3 = (y[i+1]-y[i])/(x[i+1]-x[i]);
            m4 = (y[i+2]-y[i+1])/(x[i+2]-x[i+1]);
            t[i] = getSlope(m1,m2,m3,m4);
        }
        
        // calculate p0,p1,p2,p3 for all intervals 0..(N-2), where interval
        // [x[i),x[i+1)] shall have number i (this means that the last interval
        // has number N-2)
        for (int i = 0; i < N-1; i++) {
            p0[i] = y[i];
            p1[i] = t[i];
            p2[i] = ( 3.0*(y[i+1]-y[i])/(x[i+1]-x[i]) - 2.0*t[i] - t[i+1] ) / (x[i+1]-x[i]);
            p3[i] = ( t[i] + t[i+1] - 2.0*(y[i+1]-y[i])/(x[i+1]-x[i]) ) / ( (x[i+1]-x[i])*(x[i+1]-x[i]) );
        }
     
        delete [] t;
    }// spline



  }//ns interaction
}//ns espressopp
