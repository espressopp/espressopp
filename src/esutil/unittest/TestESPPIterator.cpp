#define BOOST_TEST_MODULE ESPPIterator

#include "ut.hpp"

#include <vector>
#include "esutil/ESPPIterator.hpp"

using namespace espresso::esutil;

BOOST_AUTO_TEST_CASE(Container) {
  const int N = 100;

  std::vector< int > v;
  for (int i = 0; i < N; i++)
    v.push_back(i);

  ESPPIterator< std::vector< int > > esppit(v);
  BOOST_CHECK_EQUAL(*esppit, 0);
  for (int i = 1; i < N; i++) {
    ++esppit;
    BOOST_CHECK_EQUAL(*esppit, i);
    BOOST_CHECK(esppit.isValid());
    BOOST_CHECK(!esppit.isDone());
  }
  
  ++esppit;
  BOOST_CHECK(!esppit.isValid());
  BOOST_CHECK(esppit.isDone());
}

BOOST_AUTO_TEST_CASE(emptyContainer) {
  std::vector< int > v;
  ESPPIterator< std::vector< int > > esppit(v);
  BOOST_CHECK(!esppit.isValid());
  BOOST_CHECK(esppit.isDone());
}
