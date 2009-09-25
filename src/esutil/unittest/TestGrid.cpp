#include <acconfig.hpp>
#define BOOST_TEST_MODULE grid
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cmath>
#include "../Grid.hpp"
#include "Real3D.hpp"

using namespace espresso;
using namespace esutil;

typedef SmallVector<float, 2> Real2D;

BOOST_AUTO_TEST_CASE(grid_2d) {
  typedef Grid<Real2D>::CellIdentifier CellIdentifier;
  Real2D cornerA = (float []) {2,4};
  Real2D cornerB = (float []) {5,8};
  size_t cells[] = { 2, 3};
  Real2D cellsize   = (float []) { 3./2., 4./3. };
  Real2D cellsize_i = (float []) { 2./3., 3./4. };

  Grid<Real2D> grid(Box<Real2D>(cornerA, cornerB), cells);

  BOOST_CHECK_EQUAL(CellIdentifier(grid.getNumCells()),
                    CellIdentifier(cells));
  BOOST_CHECK_SMALL((Real2D(grid.getCellSize()) - Real2D(cellsize)).sqr(), 1e-10f);
  BOOST_CHECK_SMALL((Real2D(grid.getInverseCellSize()) - Real2D(cellsize_i)).sqr(), 1e-10f);

  CellIdentifier id;
  id = grid.locate(cornerA);
  BOOST_CHECK_EQUAL(id, CellIdentifier((size_t []) { 0, 0 }));
  BOOST_CHECK_EQUAL(grid.linearize(id), size_t(0));

  id = grid.locate(cornerB - Real2D(0.001));
  BOOST_CHECK_EQUAL(id, CellIdentifier((size_t []) { 1, 2 }));
  BOOST_CHECK_EQUAL(grid.linearize(id), size_t(5));

  id = grid.locate(Real2D((float []) { 2.9, 6.0 }));
  BOOST_CHECK_EQUAL(id, CellIdentifier((size_t []) { 0, 1 }));
  BOOST_CHECK_EQUAL(grid.linearize(id), size_t(1));

  id = grid.locate(Real2D(0.0));
  BOOST_CHECK_EQUAL(id, Grid<Real2D>::notACell);

  id = grid.locate(Real2D(20.0));
  BOOST_CHECK_EQUAL(id, Grid<Real2D>::notACell);

  id = grid.locate(Real2D((float []) { 20.0, 5.0 }));
  BOOST_CHECK_EQUAL(id, Grid<Real2D>::notACell);
}

BOOST_AUTO_TEST_CASE(grid_2d_reshape) {
  typedef Grid<Real2D>::CellIdentifier CellIdentifier;

  Grid<Real2D> grid(Box<Real2D>((float []) {2,4}, (float []) {5,8}), (size_t []) { 2, 3});
  {
    size_t cells[] = { 2, 3};
    Real2D cellsize   = (float []) { 3./2., 4./3. };
    Real2D cellsize_i = (float []) { 2./3., 3./4. };

    BOOST_CHECK_EQUAL(CellIdentifier(grid.getNumCells()),
                      CellIdentifier(cells));
    BOOST_CHECK_SMALL((Real2D(grid.getCellSize()) - Real2D(cellsize)).sqr(), 1e-10f);
    BOOST_CHECK_SMALL((Real2D(grid.getInverseCellSize()) - Real2D(cellsize_i)).sqr(), 1e-10f);
  }

  {
    grid.setNumCells(1, 10);
    size_t cells[] = { 2, 10};
    Real2D cellsize   = (float []) { 3./2., 4./10. };
    Real2D cellsize_i = (float []) { 2./3., 10./4. };

    BOOST_CHECK_EQUAL(CellIdentifier(grid.getNumCells()),
                      CellIdentifier(cells));
    BOOST_CHECK_SMALL((Real2D(grid.getCellSize()) - Real2D(cellsize)).sqr(), 1e-10f);
    BOOST_CHECK_SMALL((Real2D(grid.getInverseCellSize()) - Real2D(cellsize_i)).sqr(), 1e-10f);
  }

  {
    grid.setDomain(Box<Real2D>((float []) {0, 1.1}, (float []) {-2, 5.1}));
    size_t cells[] = { 2, 10};
    Real2D cellsize   = (float []) { 1., 4./10. };
    Real2D cellsize_i = (float []) { 1., 10./4. };

    BOOST_CHECK_EQUAL(CellIdentifier(grid.getNumCells()),
                      CellIdentifier(cells));
    BOOST_CHECK_SMALL((Real2D(grid.getCellSize()) - Real2D(cellsize)).sqr(), 1e-10f);
    BOOST_CHECK_SMALL((Real2D(grid.getInverseCellSize()) - Real2D(cellsize_i)).sqr(), 1e-10f);
  }
}

BOOST_AUTO_TEST_CASE(grid_3d) {
  typedef Grid<Real3D>::CellIdentifier CellIdentifier;

  Real3D cornerA(2,4,5);
  Real3D cornerB(5,8,6);

  Grid<Real3D> grid(Real3DBox(cornerA, cornerB), (size_t []) {6, 11, 160});

  {
    size_t cells[] = { 6, 11, 160};
    Real3D cellsize   (3./6., 4./11., 1./160.);
    Real3D cellsize_i (6./3., 11./4., 160./1.);

    BOOST_CHECK_EQUAL(CellIdentifier(grid.getNumCells()),
                      CellIdentifier(cells));
    BOOST_CHECK_SMALL((Real3D(grid.getCellSize()) -  cellsize).sqr(), 1e-10);
    BOOST_CHECK_SMALL((Real3D(grid.getInverseCellSize()) - cellsize_i).sqr(), 1e-10);

    CellIdentifier id;
    id = grid.locate(cornerA);
    BOOST_CHECK_EQUAL(id, CellIdentifier((size_t []) { 0, 0, 0 }));
    BOOST_CHECK_EQUAL(grid.linearize(id), size_t(0));

    id = grid.locate(cornerB - Real3D(0.00001));
    BOOST_CHECK_EQUAL(id, CellIdentifier((size_t []) { 5, 10, 159 }));
    BOOST_CHECK_EQUAL(grid.linearize(id), size_t(6*11*160 - 1));
  }

  {
    grid.setNumCells(1, 10);
    grid.setDomain(Real3DBox(Real3D(0, 1.1, 2), Real3D(-2, 5.1, 6)));

    size_t cells[] = { 6, 10, 160};
    Real3D cellsize   (2./6., 4./10., 4./160.);
    Real3D cellsize_i (6./2., 10./4., 160./4.);

    BOOST_CHECK_EQUAL(CellIdentifier(grid.getNumCells()),
                      CellIdentifier(cells));
    BOOST_CHECK_SMALL((Real3D(grid.getCellSize()) -  cellsize).sqr(), 1e-10);
    BOOST_CHECK_SMALL((Real3D(grid.getInverseCellSize()) - cellsize_i).sqr(), 1e-10);
  }
}
