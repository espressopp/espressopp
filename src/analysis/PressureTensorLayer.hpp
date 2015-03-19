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
#ifndef _ANALYSIS_PRESSURE_TENSOR_LAYER_HPP
#define _ANALYSIS_PRESSURE_TENSOR_LAYER_HPP

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
    class PressureTensorLayer : public AnalysisBaseTemplate< Tensor > {
    private:
      real z0;
      real dz;
      int direction; // 0 - x, 1 - y, 2 - z // TODO currently only z direction
                     //because all interaction templates have to be changed
    public:
      static void registerPython();

      PressureTensorLayer(shared_ptr< System > system, real _z0, real _dz) : AnalysisBaseTemplate< Tensor >(system), z0(_z0), dz(_dz){}
      virtual ~PressureTensorLayer() {}
      /*
       * calculate a pressure for z0 layer
       * particles in [z-dz, z+dz] slab will be taken into account
       */
      Tensor computeRaw() {
        System& system = getSystemRef();
        Real3D Li = system.bc->getBoxL();

        // determine the local volume size for kinetic part
        real A = Li[0] * Li[1];
        real V = A * (2 * dz);

        // compute the kinetic contribution (2/3 \sum 1/2mv^2)
        Tensor vvlocal(0.0);
        CellList realCells = system.storage->getRealCells();

        real zmin = z0 - dz;
        real zmax = z0 + dz;

        // boundary condition for orthorhombic box
        real zminBC = zmin;
        real zmaxBC = zmax;
        bool condition2 = false;
        if(zminBC<0.0){
          zminBC += Li[2];
          condition2 = true;
        }
        else if(zmaxBC>=Li[2]){
          zmaxBC -= Li[2];
          condition2 = true;
        }

        for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
          Real3D pos = cit->position();

          bool condition;
          if(condition2){
            condition = (pos[2]>=zminBC || pos[2]<zmaxBC);
          }
          else{
            condition = (pos[2]>=zminBC && pos[2]<zmaxBC);
          }

          if( condition ){
            real mass = cit->mass();
            Real3D& vel = cit->velocity();
            vvlocal += mass * Tensor(vel, vel);
          }
        }
        Tensor vv(0.0);
        boost::mpi::all_reduce(*mpiWorld, (double*)&vvlocal,6, (double*)&vv, std::plus<double>());

        vv = vv / V;

        // compute the short-range nonbonded contribution
        Tensor w(0.0);
        const InteractionList& srIL = system.shortRangeInteractions;
        for (size_t j = 0; j < srIL.size(); j++) {
          srIL[j]->computeVirialTensor(w, z0);
        }

        w = w / A;

        return ( vv + w );
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
      
      void setH0(real _h0){z0=_h0;}
      real getH0(){return z0;}
      
      void setDH(real _dh){dz=_dh;}
      real getDH(){return dz;}
    };
  }
}

#endif
