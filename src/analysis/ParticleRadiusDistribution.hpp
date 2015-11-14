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
#ifndef _ANALYSIS_PARTICLERADIUSDISTRIBUTION_HPP
#define _ANALYSIS_PARTICLERADIUSDISTRIBUTION_HPP

#include "types.hpp"
#include "AnalysisBase.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include <vector>

namespace espressopp {
  namespace analysis {
    using namespace iterator;
    /** Class to compute the particle radius distribution. */

    typedef std::vector< real > result;

    class ParticleRadiusDistribution : public AnalysisBaseTemplate< result > {
    public:
      static void registerPython();

      ParticleRadiusDistribution(shared_ptr< System > system) : AnalysisBaseTemplate< result >(system) {
    	  // default binwidth for particle radius histogram
    	  binwidth=0.1;
    	  radave=0.0;
      }
      virtual ~ParticleRadiusDistribution() {}

      result computeRaw() {
    	// get total number of particles from the system by getting
    	// local number on each CPU and then adding all up (over mpi)
        int myN, systemN;
        System& system = getSystemRef();
        mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());

        result localbin;
        result bin;
        real localaverage=0.0;

        CellList realCells = system.storage->getRealCells();
        for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
          real r = cit->radius();
          if (r<0) {
        	std::cout << "WARNING: radius of particle is negative, setting it to zero !";
        	r = 0.0;
          }
          int idx = floor( (cit->radius() + 0.5*binwidth) / binwidth );
          // increment corresponding count in bin
          ++localbin[idx];
          // add up radii (to calculate average)
          localaverage += cit->radius();
        }

        // mpi::all_reduce(*getSystem()->comm, localaverage, radave, std::plus<real>());
        // mpi::all_reduce(*getSystem()->comm, localbin, bin, std::plus<real>());

        radave /= systemN;
        
        return bin;
      }


      python::list compute() {
        python::list ret;
        result res = computeRaw();
        for (int i=0; i<res.size(); i++) {
          ret.append(res[i]);
        }
        return ret;
      }

      python::list getAverageValue() {
        python::list ret;
        real res;
        res = 0.0; //nMeasurements>0 ? newAverage : 0;
        ret.append(res);
        res = 0.0; //nMeasurements>0 ? newVariance : 0;
        ret.append(res); //sqrt(res/(nMeasurements-1)));
        return ret;
      }

      void resetAverage() {
        //newAverage   = 0;
        //lastAverage  = 0;
        //newVariance  = 0;
        //lastVariance = 0;
      }

      void updateAverage(result res) {
    	/*
    	if (nMeasurements > 0) {
    	  if (nMeasurements == 1) {
              newAverage     = res;
              lastAverage    = newAverage;
          } else {
              newAverage   = lastAverage  + (res - lastAverage) / nMeasurements;
              newVariance  = lastVariance + (res - lastAverage) * (res - newAverage);
              lastAverage  = newAverage;
              lastVariance = newVariance;
          }
    	}
        return;
        */
      }
    private:
      // bin width of the histogram for the radii
      real binwidth;
      // average radius
      real radave;
    };
  }
}

#endif
