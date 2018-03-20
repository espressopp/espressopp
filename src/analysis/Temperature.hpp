/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
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
#include "Observable.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include <boost/unordered_set.hpp>

namespace espressopp {
  namespace analysis {
    using namespace iterator;
    /** Class to compute the temperature. */
    class Temperature : public Observable {
    private:
     boost::unordered_set<longint> valid_type_ids;
     bool has_types;

    public:
      static void registerPython();

      Temperature(shared_ptr< System > system): Observable(system), eKin_(0.0) {
        has_types = false;
        result_type = real_scalar;
      }
      virtual ~Temperature() {}

      real compute_real() const {
      int myN, systemN;
      real sumT = 0.0;
      real v2sum = 0.0;
      System& system = getSystemRef();
      int count = 0;

      if (system.storage->getFixedTuples()){  // AdResS - hence, need to distinguish between CG and AT particles.
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
                      if (!has_types || valid_type_ids.count(at.type())) {
                          Real3D vel = at.velocity();
                          v2sum += at.mass() * (vel * vel);
                          count += 1;
                      }
                  }
            } else {   // If not, use CG particle itself for calculation.
                  if (!has_types || valid_type_ids.count(cit->type())) {
                      Real3D vel = cit->velocity();
                      v2sum += cit->mass() * (vel * vel);
                      count += 1;
                  }
            }

          }

          myN = count;
      } else {  // No AdResS - just iterate over all particles
          CellList realCells = system.storage->getRealCells();
          for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            Real3D vel = cit->velocity();
            if (!has_types || valid_type_ids.count(cit->type())) {
                v2sum += cit->mass() * (vel * vel);
                count += 1;
            }
          }

          myN = count;
      }

      mpi::all_reduce(*getSystem()->comm, v2sum, sumT, std::plus<real>());
      mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());

      eKin_ = 0.5*sumT;
      return sumT / (3.0 * systemN);
      }

    real getEkin() const { return eKin_; }

    private:
     void addType(longint type_id) {
       valid_type_ids.insert(type_id);
       has_types = true;
     }
     bool removeType(longint type_id) {
      bool ret_val =  valid_type_ids.erase(type_id);
      has_types = valid_type_ids.size() > 0;
      return ret_val;
     }
     mutable real eKin_;
    };
  }
}

#endif
