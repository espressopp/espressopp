#define BOOST_TEST_MODULE CellListIterator

#include "ut.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "CellListIterator.hpp"
#include <vector>
#include <iostream>

using namespace espresso;
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
  for (int i = 0; i < NCELL; i++) {
    cells.push_back(cell);
    cl.push_back(&cells[i]);
  }

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
      p.p.identity = i*NP + j;
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
    ++occ[cit->p.identity];
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
      p.p.identity = i*NP + j;
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
    ++occ[cit->p.identity];
    ++cit;
  }

  // check that every particle was visited once
  for (int i = 0; i < NCELL*NP; ++i) {
    BOOST_CHECK_EQUAL(occ[i], 1);
  }

  BOOST_CHECK(!cit.isValid());
  BOOST_CHECK(cit.isDone());
}
