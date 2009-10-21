#include <acconfig.hpp>
#define BOOST_TEST_MODULE Grid
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../Grid.hpp"

BOOST_AUTO_TEST_CASE(testGrid) {
  Grid testGrid(1, 2, 3);
  BOOST_REQUIRE_EQUAL(testGrid.getNumberOfCells(), integer(6));

  BOOST_CHECK_EQUAL(testGrid.getGridSize(0), integer(1));
  BOOST_CHECK_EQUAL(testGrid.getGridSize(1), integer(2));
  BOOST_CHECK_EQUAL(testGrid.getGridSize(2), integer(3));

  for (integer i = 0; i < 6; ++i) {
    integer x, y, z;
    testGrid.getGridPosition(i, x, y, z);
    BOOST_CHECK_EQUAL(testGrid.getLinearIndex(x, y, z), i);
  }
}
