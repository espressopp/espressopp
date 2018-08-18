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

#define BOOST_TEST_MODULE NeighborCellListIterator

#include "ut.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "../NeighborCellListIterator.hpp"
#include <vector>
#include <iostream>

using namespace espressopp;
using namespace iterator;
using namespace std;

BOOST_AUTO_TEST_CASE(DefaultConstructor) {
  NeighborCellListIterator cit;
  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
}

BOOST_AUTO_TEST_CASE(EmptyList) {
  NeighborCellList ncl;
  NeighborCellListIterator cit(ncl, true);
  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
  
  NeighborCellListIterator cit2(ncl, false);
  BOOST_CHECK(!cit2.isValid());
  BOOST_CHECK(cit2.isDone());
}

BOOST_AUTO_TEST_CASE(ListOfEmptyCells) {
  const int NCELL = 10;

  Cell cell;
  LocalCellList cl;
  for (int i = 0; i < NCELL; i++)
    cl.push_back(cell);
  BOOST_CHECKPOINT("after CellList init");

  NeighborCellList ncl;
  for (int i = 0; i < NCELL; i++)
    ncl.push_back(NeighborCellInfo(cl[i], (i%2 == 0)));
  BOOST_CHECKPOINT("after NeighborCellList init");

  NeighborCellListIterator cit(ncl, true);
  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());

  NeighborCellListIterator cit2(ncl, false);
  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
}


BOOST_AUTO_TEST_CASE(FullCellList) {
  const int NCELL = 13;
  const int NP = 7;

  Particle p;
  Cell cell;
  LocalCellList cells;
  // create cells and particles
  for (int i = 0; i < NCELL; ++i) {
    cells.push_back(cell);
    for (int j = 0; j < NP; ++j) {
      p.p.id = i*NP + j;
      cells.back().particles.push_back(p);
    }
  }

  NeighborCellList ncl;
  for (int i = 0; i < NCELL; ++i)
    ncl.push_back(NeighborCellInfo(cells[i], (i%2==0)));

  vector< int > occ(NP*NCELL, 0);
  NeighborCellListIterator cit(ncl, false);

  for (int i = 0; i < NCELL*NP; ++i) {
    BOOST_CHECK(cit.isValid());
    BOOST_CHECK(!cit.isDone());
    ++occ[cit->p.id];
    ++cit;
  }

  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());

  // check that every particle was visited once
  for (int i = 0; i < NCELL*NP; ++i) {
    BOOST_CHECK_EQUAL(occ[i], 1);
  }
}

BOOST_AUTO_TEST_CASE(FullCellListAllPairs) {
  const int NCELL = 13;
  const int NP = 7;

  Particle p;
  Cell cell;
  LocalCellList cells;
  // create cells and particles
  for (int i = 0; i < NCELL; ++i) {
    cells.push_back(cell);
    for (int j = 0; j < NP; ++j) {
      p.p.id = i*NP + j;
      cells.back().particles.push_back(p);
    }
  }

  NeighborCellList ncl;
  for (int i = 0; i < NCELL; ++i)
    ncl.push_back(NeighborCellInfo(cells[i], (i%2==0)));

  vector< int > occ(NP*NCELL, 0);
  NeighborCellListIterator cit(ncl, true);

  for (int i = 0; i < (NCELL/2+NCELL%2)*NP; ++i) {
    BOOST_CHECK(cit.isValid());
    BOOST_CHECK(!cit.isDone());
    ++occ[cit->p.id];
    ++cit;
  }

  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());

  // check that half of the particles were visited once, the others not
  int pid = 0;
  for (int cid = 0; cid < NCELL; ++cid) {
    for (int i = 0; i < NP; ++i) {
      BOOST_CHECK_EQUAL(occ[pid], 1-(cid%2));
      ++pid;
    }
  }
}

