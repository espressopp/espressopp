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

// ESPP_CLASS
#ifndef _INTERACTION_INTERPOLATION_TABLE_HPP
#define _INTERACTION_INTERPOLATION_TABLE_HPP

#include <stdio.h>
#include "types.hpp"
#include "logging.hpp"
#include "mpi.hpp"

namespace espressopp {
  namespace interaction {

    /** This class defines an interpolation table that can be used for
        tabbed potentials (pair, bond, angle).

        The interpolation tables are created by fitting cubic splines
        to the file values and interpolating energy and force values at
        each of N distances. During a simulation, these tables are used to
        interpolate energy and force values as needed. 

        The interpolation is done via a spline style. The cubic spline 
        coefficients are computed and stored at each of the N values in 
        the table. The distance is used to find the appropriate set of 
        coefficients which are used to evaluate a cubic polynomial which 
        computes the energy or force.

        An interpolation table can be used to define a mapping from a
        value (e.g. distance r, but also angles) to an energy and a force
        value. The table itself is read from a file and provides energy and
        force values for equidistant values in a given range.

    */

    class InterpolationTable {

    public:
      InterpolationTable();
      ~InterpolationTable();

      /** Read in the radius, energy, force values; creates spline tables */
      void read(mpi::communicator comm, const char* file);

      real getEnergy(real r) const;
      real getForce(real r) const;

    protected:
      /** Logger */
      static LOG4ESPP_DECL_LOGGER(theLogger);

    private:
      /** make copy constructor private because copy is not allowed. */
      InterpolationTable(const InterpolationTable&) {}

      /** Reading values from file, control processor only; returns
          number of valid entries, error if value is less than 2 

          if dummy is true, values will not be stored in arrays r, e, f
      */
      int readFile(const char* file, bool dummy);

      /** Spline read-in values. */
      void spline(const real* x, const real* y, int n,
                  real yp1, real ypn, real* y2);

      /** Spline interpolation */
      real splineInterpolation(real r, const real* fn, const real* fn2) const;

      int N;  // number of read values

      real inner;
      real delta;
      real invdelta;
      real deltasq6;

      bool allocated;

      real *radius;
      real *energy;
      real *force;

      real *energy2;  // used for spline interpolation
      real *force2;   // used for spline interpolation
    };
  }
}

#endif
