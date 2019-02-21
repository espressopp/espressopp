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

#define BOOST_TEST_MODULE CellListAllPairsIterator

#define LOG4ESPP_LEVEL_TRACE
#include "ut.hpp"
#include "log4espp.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include <vector>

using namespace espressopp;
using namespace iterator;
using namespace std;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    // uncomment for debugging purposes
    // log4espp::Logger::getRoot().setLevel(log4espp::Logger::TRACE);
  }
};

BOOST_GLOBAL_FIXTURE(LoggingFixture);

BOOST_AUTO_TEST_CASE(defaultConstructor) {
  CellListAllPairsIterator it;
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());
}

BOOST_AUTO_TEST_CASE(emptyList) {
  CellList cl;
  CellListAllPairsIterator it(cl);
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());
}

BOOST_AUTO_TEST_CASE(emptyCellsNoNeighbors) {
  // create empty cells with no neighbors 
  const int NCELL = 3;

  vector < Cell > cell(NCELL);
  CellList cl;

  // create cells and fill in neighbors
  for (int i = 0; i < NCELL; ++i)
    cl.push_back(&cell[i]);
  
  CellListAllPairsIterator it(cl);
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());
}

BOOST_AUTO_TEST_CASE(emptyCellsEmptyNeighbors) {
  // create empty cells with empty neighbors 
  const int NCELL = 3;
  const int NNCELL = 5;

  vector < Cell > cell(NCELL);
  vector < Cell > neighbor(NNCELL);
  CellList cl;

  // create cells and fill in neighbors
  for (int i = 0; i < NCELL; ++i) {
    for (int j = 0; j < NNCELL; ++j)
      cell[i].neighborCells.push_back
	(NeighborCellInfo(neighbor[i], true));
    cl.push_back(&cell[i]);
  }

  CellListAllPairsIterator it(cl);
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());
}

BOOST_AUTO_TEST_CASE(emptyCellsFullNeighbors) {
  // create empty cells and neighbors with particles
  const int NCELL = 3;
  const int NNCELL = 5;
  const int NNP = 7;

  Particle p;

  vector < Cell > cell(NCELL);
  vector < Cell > neighbor(NNCELL);
  CellList cl;

  // create neighbor cells
  for (int i = 0; i < NNCELL; i++) {
    for (int j = 0; j < NNP; j++)
      neighbor[i].particles.push_back(p);
  }

  // create cells
  for (int i = 0; i < NCELL; ++i) {
    for (int j = 0; j < NNCELL; ++j)
      cell[i].neighborCells.push_back
	(NeighborCellInfo(neighbor[i], true));
    cl.push_back(&cell[i]);
  }

  CellListAllPairsIterator it(cl);
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());
}

BOOST_AUTO_TEST_CASE(fullCellsNoNeighbors) {
  // populate cells with particles, but add no neighbors
  const int NCELL = 5;
  const int NP = 11;

  Particle p;
  vector < Cell > cells(NCELL);
  CellList cl;

  vector< int > occupancy(NP*NCELL*NP*NCELL, 0);

  for (int i = 0; i < NCELL; ++i) {
    for (int j = 0; j < NP; ++j) {
      p.id() = i*NP+j;
      cells[i].particles.push_back(p);
    }
    cl.push_back(&cells[i]);
  }

  // run the loop
  int steps = 0;
  CellListAllPairsIterator it(cl);
  for (;it.isValid(); ++it) {
    ++steps;
    int pid1 = it->first->id();
    int pid2 = it->second->id();
    ++occupancy[pid1*NP*NCELL+pid2];
    ++occupancy[pid2*NP*NCELL+pid1];
    BOOST_TEST_CHECKPOINT("step " << steps << ", pair: (" << it->first->id() << ", "
		     << it->second->id() << ")");
  }
  BOOST_CHECK(it.isDone());

  // and check the results
  int expected = NCELL*(NP*(NP-1))/2;
  BOOST_REQUIRE_EQUAL(steps, expected);

  for (int pid1 = 0; pid1 < NP*NCELL; ++pid1) {
    for (int pid2 = 0; pid2 < NP*NCELL; ++pid2) {
      int cid1 = pid1 / NP;
      int cid2 = pid2 / NP;
      int occ = occupancy[pid1*NP*NCELL + pid2];
      if (pid1 == pid2) {
	// pairs with identical ids should not be touched
	BOOST_CHECK_MESSAGE(occ==0,
			    "occupancy(" << pid1 << ", " << pid2 << ") == " 
			    << occ << ", should be 0");
      } else if (cid1 == cid2) {
	// pairs in the same cell should be touched once
	BOOST_CHECK_MESSAGE(occ==1,
			    "occupancy(" << pid1 << ", " << pid2 << ") (same cell) == " 
			    << occ << ", should be 1");
      } else {
	// pairs in different cells should not be touched
	BOOST_CHECK_MESSAGE(occ==0,
			    "occupancy(" << pid1 << ", " << pid2 << ") (diff. cells) == " 
			    << occ << ", should be 0");
      }
    };
  }
}

BOOST_AUTO_TEST_CASE(fullCellsFullNeighbors) {
  const int NP = 7;
  const int NCELL = 11;
  // size of the neighborhood
  const int NN = 2;

  Particle p;
  vector< Cell > cell(NCELL);
  CellList cl;

  vector< int > occupancy(NP*NCELL*NP*NCELL, 0);

  // create cells
  for (int cid = 0; cid < NCELL; ++cid) {
    // put the cell into the list
    cl.push_back(&cell[cid]);

    // create particles in the cell
    for (int pid = 0; pid < NP; ++pid) {
      p.id() = cid*NP+pid;
      cell[cid].particles.push_back(p);
    }

    // set up neighborhood
    for (int j = 1; j <= NN; j++) {
      int nbcell;

      nbcell = (cid+j) % NCELL;
      cell[cid].neighborCells.push_back
	(NeighborCellInfo(cell[nbcell], true));

      nbcell = cid-j;
      if (nbcell < 0) nbcell += NCELL;
      cell[cid].neighborCells.push_back
	(NeighborCellInfo(cell[nbcell], false));
    }
  }

  // now run the loop
  int steps = 0;
  CellListAllPairsIterator it(cl);
  for (;it.isValid(); ++it) {
    ++steps;
    int pid1 = it->first->id();
    int pid2 = it->second->id();
    ++occupancy[pid1*NP*NCELL+pid2];
    ++occupancy[pid2*NP*NCELL+pid1];
    BOOST_TEST_CHECKPOINT("step " << steps << ", pair: (" << it->first->id() << ", "
		     << it->second->id() << ")");
  }
  BOOST_CHECK(it.isDone());

  // and check the results
  int expected = 
    NCELL*
    (((NP*(NP-1))/2) // sum of self pairs
     + NP*NN*NP); // sum of neighbor pairs
  BOOST_CHECK_EQUAL(steps, expected);

  for (int pid1 = 0; pid1 < NP*NCELL; ++pid1) {
    for (int pid2 = 0; pid2 < NP*NCELL; ++pid2) {
      int cid1 = pid1 / NP;
      int cid2 = pid2 / NP;
      int dcid = abs(cid1-cid2);
      dcid = min(dcid, NCELL-dcid);
      int occ = occupancy[pid1*NP*NCELL + pid2];
      if (pid1 == pid2) {
	// pairs with identical ids should not be touched
	BOOST_CHECK_MESSAGE(occ==0,
			    "occupancy(" << pid1 << ", " << pid2 << ") == " 
			    << occ << ", should be 0");
      } else if (dcid <= NN) {
	// pairs in the same neighborhood should be touched once
	BOOST_CHECK_MESSAGE(occ==1,
			    "occupancy(" << pid1 << ", " << pid2 
			    << ", cid1==" << cid1 << ", cid2==" << cid2 
			    << ") == " 
			    << occ << ", should be 1");
      } else {
	// pairs in other cells should not be touched
	BOOST_CHECK_MESSAGE(occ==0,
			    "occupancy(" << pid1 << ", " << pid2 
			    << ", cid1==" << cid1 << ", cid2==" << cid2 
			    << ") == " 
			    << occ << ", should be 0");
      }

    };
  }
}
