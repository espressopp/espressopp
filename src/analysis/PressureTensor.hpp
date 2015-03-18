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
#ifndef _ANALYSIS_PRESSURE_TENSOR_HPP
#define _ANALYSIS_PRESSURE_TENSOR_HPP

#include "types.hpp"
#include "AnalysisBase.hpp"
#include "Tensor.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"

namespace espressopp {
  namespace analysis {

    using namespace iterator;
    using namespace interaction;

    /** Class to compute the pressure tensor. */
    class PressureTensor : public AnalysisBaseTemplate< Tensor > {
    public:
      static void registerPython();

      PressureTensor(shared_ptr< System > system) : AnalysisBaseTemplate< Tensor >(system) {}
      virtual ~PressureTensor() {}

      Tensor computeRaw() {
        System& system = getSystemRef();

        // determine volume of the box
        Real3D Li = system.bc->getBoxL();
        real V = Li[0] * Li[1] * Li[2];

        // compute the kinetic contribution (2/3 \sum 1/2mv^2)
        Tensor vvLocal(0.0);
        Tensor vv;

        CellList realCells = system.storage->getRealCells();
        for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
          real mass = cit->mass();
          Real3D& vel = cit->velocity();
          vvLocal += mass * Tensor(vel, vel);
        }

        boost::mpi::all_reduce(*mpiWorld, (double*)&vvLocal,6, (double*)&vv, std::plus<double>());

        // compute the short-range nonbonded contribution
        Tensor wij(0.0);
        // virial is already reduced
        const InteractionList& srIL = system.shortRangeInteractions;
        for (size_t j = 0; j < srIL.size(); j++) {
          srIL[j]->computeVirialTensor(wij);
        }

        return (vv + wij) / V;
      }
      

      python::list compute() {
        python::list ret;
        Tensor res = computeRaw();
        for (int i=0; i<6; i++) {
          ret.append(res[i]);
        }
        return ret;
      }

      python::list getAverageValue() {
        python::list ret;
        Tensor res;
        res = nMeasurements>0 ? newAverage : Tensor(0);
        for (int i=0; i<6; i++) {
          ret.append(res[i]);
        }
        res = nMeasurements>0 ? newVariance : Tensor(0);
        for (int i=0; i<6; i++) {
          ret.append(sqrt(res[i]/(nMeasurements-1)));
        }
        return ret;
      }

      void resetAverage() {
        newAverage   = Tensor(0);
        lastAverage  = Tensor(0);
        newVariance  = Tensor(0);
        lastVariance = Tensor(0);
      }

      void updateAverage(Tensor res) {
        // compare e. g. Knuth TAOCP vol 2, 3rd edition, page 232
    	if (nMeasurements > 0) {
    	  if (nMeasurements == 1) {
              newAverage     = res;
              lastAverage    = newAverage;
          } else {
              newAverage  = lastAverage  + (res - lastAverage) / nMeasurements;
              for (int i=0; i<6; i++) {
                newVariance[i] = lastVariance[i] + (res[i] - lastAverage[i]) * (res[i] - newAverage[i]);
              }
              lastAverage = newAverage;
              lastVariance = newVariance;
          }
    	}
        return;
      }
      
    };
  }
}

#endif
