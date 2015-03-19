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
#ifndef _ANALYSIS_PRESSURE_TENSOR_MULTI_LAYER_HPP
#define _ANALYSIS_PRESSURE_TENSOR_MULTI_LAYER_HPP

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
    using namespace std;

    /** Class to compute the pressure tensor. */
    class PressureTensorMultiLayer : public AnalysisBaseTemplate< vector<Tensor> > {
    private:
      int n;
      real dz;
      int direction; // 0 - x, 1 - y, 2 - z // TODO currently only z direction
                     //because all interaction templates have to be changed
    public:
      static void registerPython();

      PressureTensorMultiLayer(shared_ptr< System > system, int _n, real _dz) : AnalysisBaseTemplate< vector<Tensor> >(system), n(_n), dz(_dz){}
      virtual ~PressureTensorMultiLayer() {}
      /*
       * calculate the pressure in 'n' layers along Z axis
       * the first layer has coordinate Lz/n the last - Lz.
       */
      vector<Tensor> computeRaw() {

        System& system = getSystemRef();
        system.storage->decompose();

        Real3D Li = system.bc->getBoxL();
        // determine the local volume size
        real A = Li[0] * Li[1];
        real V = A * (2 * dz);

        // n * lZ is always Lz
        real lZ = Li[2] / (double)n;

        Tensor *vvlocal = new Tensor[n];
        for(int i=0; i<n; i++) vvlocal[i] = Tensor(0.0);
        // compute the kinetic contribution (2/3 \sum 1/2mv^2)
        CellList realCells = system.storage->getRealCells();
        for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
          Real3D pos = cit->position();
          real mass = cit->mass();
          Real3D& vel = cit->velocity();
          Tensor vvt = mass * Tensor(vel, vel);

          real zminBC = pos[2] - dz;
          real zmaxBC = pos[2] + dz;

          // boundary condition for orthorhombic box
          bool boundary = false;
          if(zminBC<0.0){
            zminBC += Li[2];
            boundary = true;
          }
          else if(zmaxBC>=Li[2]){
            zmaxBC -= Li[2];
            boundary = true;
          }

          int minpos = (int)( zminBC/lZ );
          int maxpos = (int)( zmaxBC/lZ );

          if(boundary){

            for(int i = 0; i<=maxpos; i++){
              vvlocal[i] += vvt;
            }
            for(int i = minpos+1; i<n; i++){
              vvlocal[i] += vvt;
            }
          }
          else{
            for(int i = minpos+1; i<=maxpos; i++){
              vvlocal[i] += vvt;
            }
          }

        }

        Tensor *vv = new Tensor[n];
        mpi::all_reduce(*mpiWorld, (double*)&vvlocal, n, (double*)&vv, std::plus<double>());

        // compute the short-range nonbonded contribution
        Tensor *w = new Tensor[n];
        const InteractionList& srIL = system.shortRangeInteractions;
        for (size_t j = 0; j < srIL.size(); j++) {
          // srIL[j]->computeVirialTensor(w, n);
        }

        vector<Tensor> pijarr;
        for(int i=0; i<n;i++){
          vv[i] = vv[i]/V;
          w[i] = w[i] / A;

          pijarr.push_back(vv[i] + w[i]);
        }

        delete[] w;
        delete[] vv;
        delete[] vvlocal;

//        vector<Tensor> pijarr;
        return pijarr;
      }
      

      python::list compute() {
        python::list ret;
        vector<Tensor> res = computeRaw();
        for (int i=0; i<n; i++) {
          ret.append(res[i]);
        }
        return ret;
      }

      python::list getAverageValue() {
        python::list ret;
        vector<Tensor> res;
        res = nMeasurements>0 ? newAverage : vector<Tensor>(n, Tensor(0.0));
        for (int i=0; i<n; i++) {
          ret.append( res[i] );
        }
        res = nMeasurements>0 ? newVariance : vector<Tensor>(n, Tensor(0.0));
        for (int i=0; i<n; i++) {
          Tensor aux(0.0);
          for(int j=0;j<6;j++)
            aux[j] = sqrt(res[i][j]/(nMeasurements-1));
          
          ret.append(aux);
        }
        return ret;
      }

      void resetAverage() {
        newAverage   = vector<Tensor>(n, Tensor(0.0));
        lastAverage  = vector<Tensor>(n, Tensor(0.0));
        newVariance  = vector<Tensor>(n, Tensor(0.0));
        lastVariance = vector<Tensor>(n, Tensor(0.0));
      }

      void updateAverage(vector<Tensor> res) {
        // compare e. g. Knuth TAOCP vol 2, 3rd edition, page 232
    	if (nMeasurements > 0) {
    	  if (nMeasurements == 1) {
            newAverage     = res;
            lastAverage    = newAverage;
          } 
          else {
            for (int i=0; i<n; i++) {
              newAverage[i]  = lastAverage[i]  + (res[i] - lastAverage[i]) / nMeasurements;
              for (int j=0; j<6; j++) {
                newVariance[i][j] = lastVariance[i][j] + (res[i][j] - lastAverage[i][j]) * (res[i][j] - newAverage[i][j]);
              }
              lastAverage[i] = newAverage[i];
              lastVariance[i] = newVariance[i];
            }
          }
    	}
      }
      
      void setN(int _n){n=_n;}
      int getN(){return n;}
      
      void setDH(real _dh){dz=_dh;}
      real getDH(){return dz;}
    };
  }
}

#endif
