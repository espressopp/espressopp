/*
 * Array4D.cpp
 * Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
 *
 * Distributed under terms of the GNU GPLv3 license.
 */


#define BOOST_TEST_MODULE Array4Dtest
#include <boost/test/included/unit_test.hpp>
#include "esutil/Array4D.hpp"

using namespace espressopp;  // NOLINT

BOOST_AUTO_TEST_SUITE(Array4Dtest)

class T {
 public:
  T() {
    a_ = -1;
  }

  explicit T(int a) {
    a_ = a;
  }
  int a_;
};

struct Fixture {
  Fixture() {
    a4 = esutil::Array4D<T, esutil::enlarge>(0, 0, 0, 0, T());
  }
  ~Fixture() { }

  esutil::Array4D<T, esutil::enlarge> a4;
};

BOOST_FIXTURE_TEST_CASE(test1, Fixture) {
  BOOST_CHECK_EQUAL(a4.size_n(), 0);
  BOOST_CHECK_EQUAL(a4.size_m(), 0);
  BOOST_CHECK_EQUAL(a4.size_p(), 0);
  BOOST_CHECK_EQUAL(a4.size_q(), 0);
  BOOST_CHECK(a4.empty());
}

BOOST_FIXTURE_TEST_CASE(test2, Fixture) {
  T &t = a4.at(0, 0, 0, 0);
  BOOST_CHECK_EQUAL(t.a_, -1);
}

BOOST_FIXTURE_TEST_CASE(test3, Fixture) {
  a4.at(2, 0, 0, 0) = T(2);
  a4.at(2, 3, 0, 0) = T(3);
  a4.at(2, 3, 4, 0) = T(4);
  a4.at(2, 3, 4, 5) = T(5);
  a4.at(1, 1, 1, 1) = T(6);

  T &t2 = a4.at(2, 0, 0, 0);
  BOOST_CHECK_EQUAL(t2.a_, 2);
  T &t3 = a4.at(2, 3, 0, 0);
  BOOST_CHECK_EQUAL(t3.a_, 3);
  T &t4 = a4.at(2, 3, 4, 0);
  BOOST_CHECK_EQUAL(t4.a_, 4);
  T &t5 = a4.at(2, 3, 4, 5);
  BOOST_CHECK_EQUAL(t5.a_, 5);
  T &t6 = a4.at(1, 1, 1, 1);
  BOOST_CHECK_EQUAL(t6.a_, 6);

  BOOST_CHECK_EQUAL(a4.size_n(), 3);
  BOOST_CHECK_EQUAL(a4.size_m(), 4);
  BOOST_CHECK_EQUAL(a4.size_p(), 5);
  BOOST_CHECK_EQUAL(a4.size_q(), 6);

  for (int i = 0; i < a4.size_n(); i++) {
    for (int j = 0; j < a4.size_m(); j++) {
      for (int p = 0; p < a4.size_p(); p++) {
        for (int q = 0; q < a4.size_q(); q++) {
          T &tt = a4.at(i, j, p, q);
          if (i == 2 && j == 0 && p == 0 && q == 0)
            BOOST_CHECK_EQUAL(tt.a_, 2);
          else if (i == 2 && j == 3 && p == 0 && q == 0)
            BOOST_CHECK_EQUAL(tt.a_, 3);
          else if (i == 2 && j == 3 && p == 4 && q == 0)
            BOOST_CHECK_EQUAL(tt.a_, 4);
          else if (i == 2 && j == 3 && p == 4 && q == 5)
            BOOST_CHECK_EQUAL(tt.a_, 5);
          else if (i == 1 && j == 1 && p == 1 && q == 1)
            BOOST_CHECK_EQUAL(tt.a_, 6);
          else
            BOOST_CHECK_EQUAL(tt.a_, -1);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
