/*
  Copyright (C) 2012,2013,2014,2015,2016,2017,2018
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

#include "python.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
//#include "Configuration.hpp"
#include "RDFatomistic.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"

#include <boost/serialization/map.hpp>

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

using namespace espressopp;
using namespace iterator;
using namespace std;

namespace espressopp {
  namespace analysis {

    // TODO currently works correctly if L_y is shortest side length (it's fine if L_x and/or L_z are equal to L_y)
    // rdfN is a level of discretisation of rdf (how many elements it contains)
    python::list RDFatomistic::computeArray(int rdfN) const {

      System& system = getSystemRef();
      esutil::Error err(system.comm);
      Real3D Li = system.bc->getBoxL();
      Real3D Li_half = Li / 2.;

      int nprocs = system.comm->size();
      int myrank = system.comm->rank();

      real histogram [rdfN];
      for(int i=0; i<rdfN; i++) histogram[i]=0;

      real dr = Li_half[1] / (real)rdfN; // If you work with nonuniform Lx, Ly, Lz, you
      // should use for Li_half[XXX] the shortest side length. Currently using L_y

      int num_part = 0;
      map< size_t, data > config;
      for (int rank_i=0; rank_i<nprocs; rank_i++) {
        map< size_t, data > conf;
        if (rank_i == myrank) {
          shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
          CellList realCells = system.storage->getRealCells();
          for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            Particle &vp = *cit;
            FixedTupleListAdress::iterator it2;
            it2 = fixedtupleList->find(&vp);

            if (it2 != fixedtupleList->end()) {  // Are there atomistic particles for given CG particle? If yes, use those for calculation.
              std::vector<Particle*> atList;
              atList = it2->second;
              for (std::vector<Particle*>::iterator it3 = atList.begin();
                   it3 != atList.end(); ++it3) {
                Particle &at = **it3;
                int id = at.id();
                data tmp = { at.position(), at.type(), vp.id(), vp.lambda() };
                conf[id] = tmp;
              }
            }

            else {  // If not, use CG particle itself for calculation.
              std::cout << "WARNING! In RDFatomistic, no atomistic AdResS particle found!\n";
              exit(1);
            }

          }
        }

        boost::mpi::broadcast(*system.comm, conf, rank_i);

        // for simplicity we will number the particles from 0
        for (map<size_t,data>::iterator itr=conf.begin(); itr != conf.end(); ++itr) {
          data p = itr->second;
          config[num_part] = p;
          num_part ++;
        }
      }

      int num_pairs = 0;
      // now all CPUs have all particle coords and num_part is the total number of particles

      // use all cpus
      // TODO it could be a problem if   n_nodes > num_part
      int numi = num_part / nprocs + 1;
      int mini = myrank * numi;
      int maxi = mini + numi;

      if(maxi>num_part) maxi = num_part;

      int perc=0, perc1=0;
      for(int i = mini; i<maxi; i++) {
        Real3D coordP1 = config[i].pos;
        int typeP1 = config[i].type;
        int molP1 = config[i].molecule;
        real resolutionP1 = config[i].resolution;


        if (spanbased == true) {
          if((coordP1[0] < Li_half[0]+span) && (coordP1[0] > Li_half[0]-span) ) {
            for(int j = i+1; j<num_part; j++) {
              Real3D coordP2 = config[j].pos;
              int typeP2 = config[j].type;
              int molP2 = config[j].molecule;

              if (molP1 != molP2) {

                if( ( (typeP1 == target1) && (typeP2 == target2) ) || ( (typeP1 == target2) && (typeP2 == target1) ) ) {
                  Real3D distVector = (0.0, 0.0, 0.0);
                  system.bc->getMinimumImageVector(distVector, coordP1, coordP2);

                  int bin = (int)( distVector.abs() / dr );
                  if( bin < rdfN) {
                    histogram[bin] += 1.0;
                  }
                  num_pairs += 1;

                }
              }
            }
          }
        }
        else {
          if( resolutionP1 == 1.0 ) {
            for(int j = i+1; j<num_part; j++) {
              Real3D coordP2 = config[j].pos;
              int typeP2 = config[j].type;
              int molP2 = config[j].molecule;
              real resolutionP2 = config[j].resolution;

              if (molP1 != molP2) {

                if( ( (typeP1 == target1) && (typeP2 == target2) ) || ( (typeP1 == target2) && (typeP2 == target1) ) ) {

                  num_pairs += 1;

                  if( resolutionP2 == 1.0 ) {
                    Real3D distVector = (0.0, 0.0, 0.0);
                    system.bc->getMinimumImageVector(distVector, coordP1, coordP2);

                    int bin = (int)( distVector.abs() / dr );
                    if( bin < rdfN) {
                      histogram[bin] += 1.0;
                    }
                  }

                }
              }
            }
          }
        }
      }

      real totHistogram[rdfN];
      boost::mpi::all_reduce(*system.comm, histogram, rdfN, totHistogram, plus<real>());

      int tot_num_pairs = 0;
      boost::mpi::all_reduce(*system.comm, num_pairs, tot_num_pairs, plus<int>());
      int nconfigs = 1;
      real rho = (real)1.0 / (Li[0]*Li[1]*Li[2]);
      real factor = 4.0 * M_PIl * dr * rho * (real)nconfigs * (real)tot_num_pairs;

      for(int i=0; i < rdfN; i++) {
        real radius = (i + 0.5) * dr;
        totHistogram[i] /= factor * (radius*radius + dr*dr / 12.0);
      }

      python::list pyli;
      for(int i=0; i < rdfN; i++) {
        pyli.append( totHistogram[i] );
      }
      return pyli;
    }

    python::list RDFatomistic::computeArrayPathIntegral(int rdfN) const {

      System& system = getSystemRef();
      esutil::Error err(system.comm);
      Real3D Li = system.bc->getBoxL();
      Real3D Li_half = Li / 2.;

      real Xmax = Li_half[0] + span;
      real Xmin = Li_half[0] - span;

      int nprocs = system.comm->size();
      int myrank = system.comm->rank();

      real histogram [rdfN];
      for(int i=0; i<rdfN; i++) histogram[i]=0;

      real dr = Li_half[1] / (real)rdfN; // If you work with nonuniform Lx, Ly, Lz, you
      // should use for Li_half[XXX] the shortest side length. Currently using L_y

      int num_part = 0;
      map< size_t, dataPathIntegral > config;
      for (int rank_i=0; rank_i<nprocs; rank_i++) {
        map< size_t, dataPathIntegral > conf;
        if (rank_i == myrank) {
          shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
          CellList realCells = system.storage->getRealCells();
          for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

            Particle &vp = *cit;
            FixedTupleListAdress::iterator it2;
            it2 = fixedtupleList->find(&vp);

            if (it2 != fixedtupleList->end()) {
              std::vector<Particle*> atList;
              atList = it2->second;
              for (std::vector<Particle*>::iterator it3 = atList.begin();
                   it3 != atList.end(); ++it3) {
                Particle &at = **it3;

                if ( (at.position()[0] < Xmax) && (at.position()[0] > Xmin) ) {
                  int id = at.id();
                  dataPathIntegral tmp = { at.position(), at.type(), vp.id(), at.pib() };
                  conf[id] = tmp;
                }
              }
            }

            else {  // If not, use CG particle itself for calculation.
              std::cout << "WARNING! In RDFpathintegral, no atomistic AdResS particle found!\n";
              exit(1);
            }

          }
        }

        boost::mpi::broadcast(*system.comm, conf, rank_i);

        // for simplicity we will number the particles from 0
        for (map<size_t, dataPathIntegral>::iterator itr=conf.begin(); itr != conf.end(); ++itr) {
          dataPathIntegral p = itr->second;
          config[num_part] = p;
          num_part ++;
        }
      }

      int num_pairs = 0;

      int numi = num_part / nprocs + 1;
      int mini = myrank * numi;
      int maxi = mini + numi;

      if(maxi>num_part) maxi = num_part;

      int perc=0, perc1=0;
      for(int i = mini; i<maxi; i++) {
        Real3D coordP1 = config[i].pos;
        int typeP1 = config[i].type;
        int molP1 = config[i].molecule;
        int pibP1 = config[i].pib;

        if( (coordP1[0] < Xmax) && (coordP1[0] > Xmin) ) {
          for(int j = i+1; j<num_part; j++) {
            Real3D coordP2 = config[j].pos;
            int typeP2 = config[j].type;
            int molP2 = config[j].molecule;
            int pibP2 = config[j].pib;

            if ( (molP1 != molP2) && (pibP1 == pibP2) && (coordP2[0] < Xmax) && (coordP2[0] > Xmin) ) {

              if( ( (typeP1 == target1) && (typeP2 == target2) ) || ( (typeP1 == target2) && (typeP2 == target1) ) ) {

                Real3D distVector = coordP1 - coordP2;

                // minimize the distance in simulation box
                for(int ii=0; ii<3; ii++) {
                  if( distVector[ii] < -Li_half[ii] ) distVector[ii] += Li[ii];
                  if( distVector[ii] >  Li_half[ii] ) distVector[ii] -= Li[ii];
                }

                real absdist = distVector.abs();

                int bin = (int)( absdist / dr );
                if( bin < rdfN) {

                  // Normalization (see R. Potestio et al., Phys. Rev. Lett. 111, 060601 (2013), Supp. Info.)
                  real Xplus = Xmax - coordP1[0];
                  real Xminus = coordP1[0] - Xmin;

                  real h_norm = 0.0;

                  if (Xplus <= absdist) {
                    h_norm += absdist - Xplus;
                  }
                  if (Xminus <= absdist) {
                    h_norm += absdist - Xminus;
                  }

                  real v_norm = 1.0/(2.0*M_PIl*dr*absdist*(2.0*absdist - h_norm));

                  histogram[bin] += v_norm;
                }
                num_pairs += 1;

              }
            }
          }
        }

      }

      real totHistogram[rdfN];
      boost::mpi::all_reduce(*system.comm, histogram, rdfN, totHistogram, plus<real>());

      int tot_num_pairs = 0;
      boost::mpi::all_reduce(*system.comm, num_pairs, tot_num_pairs, plus<int>());

      // Normalization (see R. Potestio et al., Phys. Rev. Lett. 111, 060601 (2013), Supp. Info.)
      real subvol = 2.0*Li[1]*Li[2]*span;
      if(span > 0.5*Li[0]) {
        subvol = Li[1]*Li[2]*Li[0];
      }
      real rho = (real)1.0 / (Li[0]*Li[1]*Li[2]);
      for(int i=0; i < rdfN; i++) {
        totHistogram[i] /= ( (real)tot_num_pairs / subvol );
      }
      // Normalization

      python::list pyli;
      for(int i=0; i < rdfN; i++) {
        pyli.append( totHistogram[i] );
      }
      return pyli;
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real RDFatomistic::compute() const {
      return -1.0;
    }

    void RDFatomistic::registerPython() {
      using namespace espressopp::python;
      class_<RDFatomistic, bases< Observable > >
      ("analysis_RDFatomistic", init< shared_ptr< System >, int, int, real, bool >())
      .def("compute", &RDFatomistic::computeArray)
      .def("computePathIntegral", &RDFatomistic::computeArrayPathIntegral)
      ;
    }
  }
}
