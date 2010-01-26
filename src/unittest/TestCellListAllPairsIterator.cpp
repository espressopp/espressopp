#define BOOST_TEST_MODULE CellListAllPairsIterator

#include "ut.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "CellListAllPairsIterator.hpp"
#include <vector>

using namespace espresso;
using namespace std;

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

BOOST_AUTO_TEST_CASE(noNeighbors) {
  const int NP = 5;
  const int NCELL = 5;

  Particle p;
  vector < Cell > cells(NCELL);
  CellList cl;

  for (int i = 0; i < NCELL; ++i) {
    for (int j = 0; j < NP; ++j)
      cells[i].particles.push_back(p);
    cl.push_back(&cells[i]);
  }

  CellListAllPairsIterator it(cl);
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());
}

BOOST_AUTO_TEST_CASE(noParticles) {
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
	(NeighborCellInfo(neighbor[i], false));
    cl.push_back(&cell[i]);
  }

  CellListAllPairsIterator it(cl);
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());
}

BOOST_AUTO_TEST_CASE(selfPairs) {
  const int NP = 29;

  Particle p;
  Cell cell;
  CellList cl;

  cell.neighborCells.push_back(&cell);
  cl.push_back(&cell);
  for (int i = 0; i < NP; ++i) {
    p.p.id = i;
    cell.particles.push_back(p);
  }

  CellListAllPairsIterator it(cl);
  for (int i = 0; i < (NP*(NP-1))/2; i++) {
    ++it;
    BOOST_CHECK(it.isValid());
    BOOST_CHECK(!it.isDone());
    BOOST_CHECK(it->first.p.id != it->second.p.id);
  }
    
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());
}

BOOST_AUTO_TEST_CASE(fullLoop1D) {
  const int NP = 7;
  const int NCELL = 11;
  // size of the neighborhood
  const int NN = 2;

  Particle p;
  vector< Cell > cell(NCELL);
  CellList cl;

  vector< int > occupancy(NP*NCELL*NP*NCELL, 0);

  // create cells
  for (int i = 0; i < NCELL; ++i) {
    // create particles
    for (int j = 0; j < NP; ++j) {
      p.p.id = i*NP+j;
      cell[i].particles.push_back(p);
    }
    // set up neighborhood
    for (int j = -NN; j <= NN; j++) {
      int nbcell = (i+j) % NCELL;
      if (nbcell < 0) nbcell += NCELL;
      cell[i].neighborCells.push_back(&cell[nbcell]);
    }
    cl.push_back(&cell[i]);
  }

  CellListAllPairsIterator it(cl);
  for (int i = 0; i < NCELL*(NP*(NCELL-1) + (NP*(NP-1))/2); ++i) {
    ++it;
    BOOST_CHECK(it.isValid());
    BOOST_CHECK(!it.isDone());
    int id1 = it->first.p.id;
    int id2 = it->second.p.id;
    BOOST_CHECK_NE(id1, id2);
    if (id1 > id2) { swap(id1, id2); }
    ++occupancy[id1*NP*NCELL + id2];
  }
    
  BOOST_CHECK(!it.isValid());
  BOOST_CHECK(it.isDone());

  for (int i = 0; i < NP*NCELL; ++i) {
    for (int j = 0; j < i; ++j) {
      int occ = occupancy[i*NP*NCELL + j];
      int icell = i / NP;
      int jcell = j / NP;
      if (abs(icell - jcell) <= NN) {
	// check that each pair in the neighborhood was looped once
	BOOST_CHECK_EQUAL(occ, 1);
      } else {
	// check that each pair not in the neighborhood is not looped
	BOOST_CHECK_EQUAL(occ, 0);
      }
    }
  }
}
