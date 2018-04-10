/*
  Copyright (C) 2017,2018
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
#include "RadGyrXProfilePI.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"
#include "math.h"

#include <boost/serialization/map.hpp>

using namespace espressopp;
using namespace iterator;
using namespace std;

namespace espressopp {
  namespace analysis {

    python::list RadGyrXProfilePI::computeArray(int splitN, int ntrotter, int ptype) const {

      System& system = getSystemRef();
      esutil::Error err(system.comm);
      real Li = system.bc->getBoxL()[0];

      real * histogram = 0;
      histogram = new real[splitN];
      for(int i=0; i<splitN; i++) histogram[i]=0.0;
      int * countsarray = 0;
      countsarray = new int[splitN];
      for(int i=0; i<splitN; i++) countsarray[i]=0;

      real dr = Li / (real)splitN;
      python::list pyli;

      real pos = 0.0;
      int bin = 0;
      if(system.storage->getFixedTuples()) {
        shared_ptr<FixedTupleListAdress> fixedtupleList=system.storage->getFixedTuples();
        CellList realCells = system.storage->getRealCells();

        for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
          Particle &vp = *cit;
          FixedTupleListAdress::iterator it2;
          it2 = fixedtupleList->find(&vp);

          if(vp.type() == ptype) {
            pos = vp.position()[0];
            if (pos < 0.0)
            {
              bin = floor ((pos+Li)/dr);
            }
            else if(pos > Li)
            {
              bin = floor ((pos-Li)/dr);
            }
            else
            {
              bin = floor (pos/dr);
            }

            real radgyr = 0.0;

            if (it2 != fixedtupleList->end()) {
              std::vector<Particle*> atList;
              atList = it2->second;
              for (std::vector<Particle*>::iterator it3 = atList.begin();
                   it3 != atList.end(); ++it3) {
                Particle &at = **it3;
                radgyr += (vp.position()-at.position()).sqr();
              }
            }
            else {
              std::cout << "No atomistic particles found in RadGyrXProfilePI::computeArray function.";
              exit(1);
            }

            radgyr = sqrt(radgyr/ntrotter);

            histogram[bin] += radgyr;
            countsarray[bin] += 1;
          }
        }

      }
      else {
        std::cout << "No tuple list found in RadGyrXProfilePI::computeArray function.";
        exit(1);
      }

      int total_num_part = 0;

      real * totHistogram = 0;
      totHistogram = new real[splitN];
      for(int i=0; i<splitN; i++) totHistogram[i]=0.0;

      int * totCountsarray = 0;
      totCountsarray = new int[splitN];
      for(int i=0; i<splitN; i++) totCountsarray[i]=0;

      boost::mpi::all_reduce(*mpiWorld, histogram, splitN, totHistogram, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, countsarray, splitN, totCountsarray, std::plus<int>());

      for(int i=0; i < splitN; i++) {
        if(totCountsarray[i] != 0) totHistogram[i] /= totCountsarray[i];
      }

      for(int i=0; i < splitN; i++) {
        pyli.append( totHistogram[i] );
      }

      delete[] totCountsarray;
      totCountsarray = 0;
      delete[] totHistogram;
      totHistogram = 0;
      delete[] countsarray;
      countsarray = 0;
      delete[] histogram;
      histogram = 0;

      return pyli;
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real RadGyrXProfilePI::compute() const {
      return -1.0;
    }

    using namespace boost::python;

    void RadGyrXProfilePI::registerPython() {
      using namespace espressopp::python;
      class_<RadGyrXProfilePI, bases< Observable > >
      ("analysis_RadGyrXProfilePI", init< shared_ptr< System > >())
      .def("compute", &RadGyrXProfilePI::computeArray)
      ;
    }

  }
}
