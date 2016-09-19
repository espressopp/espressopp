/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research

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
#include "Configuration.hpp"
#include "AdressDensity.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"
#include "math.h"

#include <boost/serialization/map.hpp>

using namespace espressopp;
using namespace iterator;
using namespace std;

namespace espressopp {
  namespace analysis {

    // splitN is a level of discretisation of density profile (how many elements it contains)
    python::list AdressDensity::computeArray(int splitN) const {

      System& system = getSystemRef();
      const bc::BC& bc = *getSystemRef().bc;
      esutil::Error err(system.comm);
      real Li = system.bc->getBoxL()[1]/2.0;

      Real3D center(0.0,0.0,0.0);
      if (verletList->getAdrCenterSet()) { //adress region centre is fixed in space
        center=**(verletList->getAdrPositions().begin());
      }

      real * histogram = 0;
      histogram = new real[splitN];
      for(int i=0;i<splitN;i++) histogram[i]=0.0;

      real dr = Li / (real)splitN;
      int num_part = 0;

      Real3D pos(0.0,0.0,0.0);
      int id = 0;
      int bin = 0;
      if(system.storage->getFixedTuples()){
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
                          pos = at.position();
                          id = at.id();

                          if (!(verletList->getAdrCenterSet())) {
                            std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();   // get positions
                            Real3D pa = **it2;
                            Real3D dist3D;
                            bc.getMinimumImageVectorBox(dist3D,pos, pa); // calculate vector between particle and first center
                            real distmin = sqrt(dist3D.sqr()); // calculate absolute distance

                            ++it2;
                            for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                              pa = **it2;
                              bc.getMinimumImageVectorBox(dist3D,pos, pa); // calculate vector between particle and other centers
                              if (sqrt(dist3D.sqr()) < distmin) distmin = sqrt(dist3D.sqr());
                            }


                            if((distmin < Li) && (exclusions.count(id) == 0)) {
                              bin = floor (distmin/dr);

                              histogram[bin] += 1.0;
                              num_part += 1;
                            }
                          }
                          else{
                            Real3D dist3D;
                            bc.getMinimumImageVectorBox(dist3D,pos,center); //pos - center
                            real dist = sqrt(dist3D.sqr());
                            if((dist < Li) && (exclusions.count(id) == 0)) {
                              bin = floor (dist/dr);

                              histogram[bin] += 1.0;
                              num_part += 1;
                            }
                          }
                      }
                }

                else{   // If not, use CG particle itself for calculation.
                      pos = cit->position();
                      id = cit->id();
                      if (!(verletList->getAdrCenterSet())) {
                        std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();   // get positions
                        Real3D pa = **it2;
                        Real3D dist3D;
                        bc.getMinimumImageVectorBox(dist3D,pos, pa); // calculate vector between particle and first center
                        real distmin = sqrt(dist3D.sqr()); // calculate absolute distance

                        ++it2;
                        for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                          pa = **it2;
                          bc.getMinimumImageVectorBox(dist3D,pos, pa); // calculate vector between particle and other centers
                          if (sqrt(dist3D.sqr()) < distmin) distmin = sqrt(dist3D.sqr());
                        }

                        if((distmin < Li) && (exclusions.count(id) == 0)) {
                          bin = floor (distmin/dr);

                          histogram[bin] += 1.0;
                          num_part += 1;
                        }
                      }
                      else{
                        Real3D dist3D;
                        bc.getMinimumImageVectorBox(dist3D,pos,center); //pos - center
                        real dist = sqrt(dist3D.sqr());
                        if((dist < Li) && (exclusions.count(id) == 0)) {
                          bin = floor (dist/dr);

                          histogram[bin] += 1.0;
                          num_part += 1;
                        }
                      }
                }

            }

      }
      else{
            CellList realCells = system.storage->getRealCells();
            for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
              pos = cit->position();
              id = cit->id();
              if (!(verletList->getAdrCenterSet())) {
                std::vector<Real3D*>::iterator it2 = verletList->getAdrPositions().begin();   // get positions
                Real3D pa = **it2;
                Real3D dist3D;
                bc.getMinimumImageVectorBox(dist3D,pos, pa); // calculate vector between particle and first center
                real distmin = sqrt(dist3D.sqr()); // calculate absolute distance

                ++it2;
                for (; it2 != verletList->getAdrPositions().end(); ++it2) {
                  pa = **it2;
                  bc.getMinimumImageVectorBox(dist3D,pos, pa); // calculate vector between particle and other centers
                  if (sqrt(dist3D.sqr()) < distmin) distmin = sqrt(dist3D.sqr());
                }

                if((distmin < Li) && (exclusions.count(id) == 0)) {
                  bin = floor (distmin/dr);

                  histogram[bin] += 1.0;
                  num_part += 1;
                }
              }
              else{
                Real3D dist3D;
                bc.getMinimumImageVectorBox(dist3D,pos,center); //pos - center
                real dist = sqrt(dist3D.sqr());
                if((dist < Li) && (exclusions.count(id) == 0)) {
                  bin = floor (dist/dr);

                  histogram[bin] += 1.0;
                  num_part += 1;
                }
              }
            }
      }

      int total_num_part = 0;
      real * totHistogram = 0;
      totHistogram = new real[splitN];
      for(int i=0;i<splitN;i++) totHistogram[i]=0.0;

      boost::mpi::all_reduce(*mpiWorld, histogram, splitN, totHistogram, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, num_part, total_num_part, std::plus<int>());

      // normalizing
      int nconfigs = 1;
      real rho = (real)total_num_part * pow((dr / Li), 3.0);

      for(int i=0; i < splitN; i++){
        totHistogram[i] /= rho * (3*i*i + 3*i + 1);
      }

      python::list pyli;
      for(int i=0; i < splitN; i++){
        pyli.append( totHistogram[i] );
      }
      delete[] totHistogram;
      totHistogram = 0;
      delete[] histogram;
      histogram = 0;

      return pyli;
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real AdressDensity::compute() const {
      return -1.0;
    }

    using namespace boost::python;

    void AdressDensity::registerPython() {
      using namespace espressopp::python;
      class_<AdressDensity, bases< Observable > >
        ("analysis_AdressDensity", init< shared_ptr< System >, shared_ptr<VerletListAdress> >())
        .def("addExclpid", &AdressDensity::addExclpid)
        .def("compute", &AdressDensity::computeArray)
      ;
    }
  }
}
