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
#ifndef _ANALYSIS_TEMPERATURE_HPP
#define _ANALYSIS_TEMPERATURE_HPP

#include "types.hpp"
#include "AnalysisBase.hpp"
//#include "storage/DomainDecomposition.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espresso {
  namespace analysis {
    using namespace iterator;
    /** Class to compute the temperature. */
    class Temperature : public AnalysisBaseTemplate< real > {
    public:
      static void registerPython();

      Temperature(shared_ptr< System > system) : AnalysisBaseTemplate< real >(system) {}
      virtual ~Temperature() {}

      real computeRaw() {
      int myN, systemN;
      real sumT = 0.0;
      real v2sum = 0.0;
      System& system = getSystemRef();
        
      if (system.storage->getFixedTuples()){  // AdResS - hence, need to distinguish between CG and AT particles.     
          int count = 0;
          shared_ptr<FixedTupleListAdress> fixedtupleList=system.storage->getFixedTuples();
          CellList realCells = system.storage->getRealCells();
          for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {  // Iterate over all (CG) particles.              
            Particle &vp = *cit;
            FixedTupleListAdress::iterator it2;
            it2 = fixedtupleList->find(&vp);
            
            if (it2 != fixedtupleList->end()) {  // Are there atomistic particles for given CG particle? If yes, use those for calculation.
                  std::vector<Particle*> atList;
                  atList = it2->second;
                  for (std::vector<Particle*>::iterator it3 = atList.begin();
                                       it3 != atList.end(); ++it3) {
                      Particle &at = **it3;
                      Real3D vel = at.velocity();
                      v2sum += at.mass() * (vel * vel);
                      count += 1;
                  }  
            }
            
            else{   // If not, use CG particle itself for calculation.
                  Real3D vel = cit->velocity();
                  v2sum += cit->mass() * (vel * vel);
                  count += 1;
            }
            
          }
          
          myN = count;
          //std::cout << "COUNT " << count << std::endl;
      }
      
      else{  // No AdResS - just iterate over all particles          
          CellList realCells = system.storage->getRealCells();
          for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            Real3D vel = cit->velocity();
            v2sum += cit->mass() * (vel * vel);
          }
          
          myN = system.storage->getNRealParticles();
      }
      
      mpi::all_reduce(*getSystem()->comm, v2sum, sumT, std::plus<real>());
      mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());
        
      //std::cout << "Count in Temperature calculation: " << systemN << std::endl;
      //std::cout << "sumT in Temperature calculation: " << sumT << std::endl;
      return sumT / (3.0 * systemN); //NO CONSTRAINTS, SETTLE BELOW 
      //return sumT / (2.0 * systemN);  // SETTLE CONSTAINTS
      }


      python::list compute() {
        python::list ret;
        real res = computeRaw();
        ret.append(res);
        return ret;
      }

      python::list getAverageValue() {
        python::list ret;
        real res;
        res = nMeasurements>0 ? newAverage : 0;
        ret.append(res);
        res = nMeasurements>0 ? newVariance : 0;
        ret.append(sqrt(res/(nMeasurements-1)));
        return ret;
      }

      void resetAverage() {
        newAverage   = 0;
        lastAverage  = 0;
        newVariance  = 0;
        lastVariance = 0;
      }

      void updateAverage(real res) {
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
      }
    };
  }
}

#endif
