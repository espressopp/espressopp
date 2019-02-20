/*
  Copyright (C) 2019
      Jakub Krajniak <jkrajniak at gmail.com>
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz

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

#define BOOST_TEST_MODULE RealNDTest

#include "RealND.hpp"
#include <boost/test/unit_test.hpp>

using namespace espressopp;  // NOLINT

BOOST_AUTO_TEST_SUITE(RealNDTest)

// testing RealND();
BOOST_AUTO_TEST_CASE(test1) {
  RealND vec = RealND();
  BOOST_REQUIRE_EQUAL(vec.getDimension(), 0);
}

// testing RealND(int d);
BOOST_AUTO_TEST_CASE(test2) {
        RealND vec = RealND(3);
        BOOST_REQUIRE_EQUAL(vec.getDimension(), 3);
}

// testing RealND::RealND(int _dim, real v)
BOOST_AUTO_TEST_CASE(test3) {
        int size = 3;
        real value = 6.7;
        RealND vec = RealND(size, value);
        BOOST_REQUIRE_EQUAL(vec.getDimension(), size);
        for (int i = 0; i < vec.getDimension(); i++) {
            BOOST_REQUIRE_EQUAL(vec[i], value);
        }
}

// testing RealND::RealND(const RealND &v), pass by value
BOOST_AUTO_TEST_CASE(test4) {
        int size = 3;
        real value = 6.7;
        RealND vec = RealND(size, value);
        RealND new_vec = RealND(vec);
        BOOST_REQUIRE_EQUAL(new_vec.getDimension(), size);
        for (int i = 0; i < new_vec.getDimension(); i++) {
            BOOST_REQUIRE_EQUAL(vec[i], value);
        }
}

// testing RealND::RealND(const RealND &v), pass a const&
BOOST_AUTO_TEST_CASE(test5) {
        int size = 3;
        real value = 6.7;
        RealND vec = RealND(size, value);
        const RealND& g = vec;
        RealND new_vec = RealND(g);
        BOOST_REQUIRE_EQUAL(new_vec.getDimension(), size);
        for (int i = 0; i < new_vec.getDimension(); i++) {
            BOOST_REQUIRE_EQUAL(new_vec[i], value);
        }
}

// RealND::RealND(const int _dim, const std::vector<real>& v)
// test case: _dim != v.size() --> expected to throw std::runtime_error
BOOST_AUTO_TEST_CASE(test6) {
        int size = 3;
        real value = 6.7;
        std::vector<real> t(size, value);
        const std::vector<real>& g = t;
        int fake_size = 5; // size != of size. Testing exception throw in these cases
        BOOST_REQUIRE_THROW(RealND(fake_size, g), std::runtime_error);
}

// RealND::RealND(const int _dim, const std::vector<real>& v)
// _dim == v.size() --> SHOULD not throw!
BOOST_AUTO_TEST_CASE(test7) {
        int size = 3;
        real value = 6.7;
        std::vector<real> t(size, value);
        const std::vector<real>& g = t;
        int fake_size = size;
        BOOST_REQUIRE_NO_THROW(RealND(fake_size, g));
}

// RealND::RealND(const int _dim, const std::vector<real>& v)
// pass v via const ref
BOOST_AUTO_TEST_CASE(test8) {
        int size = 3;
        real value = 6.7;
        std::vector<real> t(size, value);
        const std::vector<real>& g = t;
        RealND new_vec = RealND(size, g);

        BOOST_REQUIRE_EQUAL(new_vec.getDimension(), size);
        for (int i = 0; i < new_vec.getDimension(); i++) {
            BOOST_REQUIRE_EQUAL(new_vec[i], value);
        }
}

// RealND::RealND(const int _dim, const std::vector<real>& v)
// same as test8, but passing by value the std::vector
BOOST_AUTO_TEST_CASE(test9) {
        int size = 3;
        real value = 6.7;
        std::vector<real> t(size, value);
        RealND new_vec = RealND(size, t);

        BOOST_REQUIRE_EQUAL(new_vec.getDimension(), size);
        for (int i = 0; i < new_vec.getDimension(); i++) {
            BOOST_REQUIRE_EQUAL(new_vec[i], value);
        }
}

BOOST_AUTO_TEST_SUITE_END()
