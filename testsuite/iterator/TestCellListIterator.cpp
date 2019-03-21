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

#define BOOST_TEST_MODULE CellListIterator

#include "ut.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "iterator/CellListIterator.hpp"
#include <vector>
#include <iostream>

using namespace espressopp;
using namespace iterator;
using namespace std;

BOOST_AUTO_TEST_CASE(DefaultConstructor) {
  CellListIterator cit;
  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
}

BOOST_AUTO_TEST_CASE(EmptyList) {
  CellList cl;
  CellListIterator cit(cl);
  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
}

BOOST_AUTO_TEST_CASE(ListOfEmptyCells) {
  const int NCELL = 10;

  Cell cell;
  LocalCellList cells;

  CellList cl;
  for (int i = 0; i < NCELL; i++)
    cells.push_back(cell);

  for (int i = 0; i < NCELL; i++)
    cl.push_back(&cells[i]);

  CellListIterator cit(cl);
  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
}

BOOST_AUTO_TEST_CASE(FullCellList) {
  const int NCELL = 5;
  const int NP = 7;

  Particle p;
  Cell cell;
  LocalCellList cells;

  CellList cl;
  for (int i = 0; i < NCELL; ++i) {
    cells.push_back(cell);
    for (int j = 0; j < NP; ++j) {
      p.id() = i*NP + j;
      cells.back().particles.push_back(p);
    }
  }

  for (int i = 0; i < NCELL; ++i)
    cl.push_back(&cells[i]);

  vector< int > occ(NP*NCELL, 0);
  CellListIterator cit(cl);

  for (int i = 0; i < NCELL*NP; ++i) {
    BOOST_CHECK(cit.isValid());
    BOOST_CHECK(!cit.isDone());
    ++occ[cit->id()];
    ++cit;
  }

  // check that every particle was visited once
  for (int i = 0; i < NCELL*NP; ++i) {
    BOOST_CHECK_EQUAL(occ[i], 1);
  }

  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
}

BOOST_AUTO_TEST_CASE(CellListWithHole) {
  const int NCELL = 11;
  const int NP = 7;

  Particle p;
  Cell cell;
  LocalCellList cells;

  // create cells and fill them with particles
  for (int i = 0; i < NCELL; ++i) {
    cells.push_back(cell);
    // create holes
    if (i == 0 || i == 3)
      cells.push_back(cell);
    for (int j = 0; j < NP; ++j) {
      p.id() = i*NP + j;
      cells.back().particles.push_back(p);
    }
  }
  cells.push_back(cell);

  // create cell list
  CellList cl;
  for (int i = 0; i < cells.size(); ++i)
    cl.push_back(&cells[i]);

  vector< int > occ(NP*NCELL, 0);
  CellListIterator cit(cl);

  for (int i = 0; i < NCELL*NP; ++i) {
    BOOST_REQUIRE_MESSAGE(cit.isValid(), "Failed in iteration " << i);
    BOOST_REQUIRE(!cit.isDone());
    ++occ[cit->id()];
    ++cit;
  }

  // check that every particle was visited once
  for (int i = 0; i < NCELL*NP; ++i) {
    BOOST_CHECK_EQUAL(occ[i], 1);
  }

  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
}
