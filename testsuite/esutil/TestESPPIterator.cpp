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

#define BOOST_TEST_MODULE ESPPIterator

#include "ut.hpp"

#include <vector>
#include "esutil/ESPPIterator.hpp"

using namespace espressopp::esutil;

BOOST_AUTO_TEST_CASE(DefaultConstructor) {
  ESPPIterator< std::vector< int > > esppit;
  BOOST_CHECK(!esppit.isValid());
  BOOST_CHECK(esppit.isDone());
}

BOOST_AUTO_TEST_CASE(Container) {
  const int N = 100;

  std::vector< int > v;
  for (int i = 0; i < N; i++)
    v.push_back(i);

  ESPPIterator< std::vector< int > > esppit(v);
  BOOST_CHECK_EQUAL((*esppit), 0);
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

BOOST_AUTO_TEST_CASE(FullRange) {
  std::vector< int > v;
  for (int i = 0; i < 10; i++)
    v.push_back(i);

  int i = 0;
  for (ESPPIterator< std::vector< int > > esppit(v.begin(), v.end());
       !esppit.isDone(); ++esppit) {
    BOOST_CHECK_EQUAL(*esppit, i);
    ++i;
  }
}
  
BOOST_AUTO_TEST_CASE(SubRange) {
  std::vector< int > v;
  for (int i = 0; i < 10; i++)
    v.push_back(i);

  std::vector< int >::iterator it = v.begin();
  ++it;
  ++it;
  ++it;
  
  std::vector< int >::iterator end = v.begin();
  ++end;
  ++end;
  ++end;
  ++end;
  ++end;
  ++end;

  int i = 3;
  for (ESPPIterator< std::vector< int > > esppit(it, end);
       !esppit.isDone(); ++esppit) {
    BOOST_CHECK_EQUAL(*esppit, i);
    ++i;
  }
  BOOST_CHECK_EQUAL(i, 6);
}

BOOST_AUTO_TEST_CASE(CopyConstructorEmpty) {
  ESPPIterator< std::vector< int > > esppit;
  ESPPIterator< std::vector< int > > esppit2(esppit);
  BOOST_CHECK(!esppit2.isValid());
  BOOST_CHECK(esppit2.isDone());
}

BOOST_AUTO_TEST_CASE(CopyConstructor) {
  std::vector< int > v;
  for (int i = 0; i < 10; ++i) v.push_back(i);

  ESPPIterator< std::vector< int > > esppit(v);
  BOOST_CHECK_EQUAL((*esppit), 0);
  ++esppit;
  BOOST_CHECK_EQUAL((*esppit), 1);

  // check that the copy points to the same value
  ESPPIterator< std::vector< int > > esppit2(esppit);
  BOOST_CHECK_EQUAL(*esppit2, 1);
  
  // check that if you advance the copy, the original is not affected
  ++esppit2;
  BOOST_CHECK(!esppit2.isDone());
  BOOST_CHECK_EQUAL(*esppit, 1);
  BOOST_CHECK_EQUAL(*esppit2, 2);

  // check that if you advance the original, the copy is not affected
  ++esppit;
  BOOST_CHECK(!esppit2.isDone());
  BOOST_CHECK_EQUAL(*esppit, 2);
  BOOST_CHECK_EQUAL(*esppit2, 2);

  // check that you can fully use the copy
  for (int i = 3; i < 10; ++i) {
    ++esppit2;
    BOOST_CHECK(!esppit2.isDone());
    BOOST_CHECK_EQUAL(*esppit2, i);
  }
  ++esppit2;
  BOOST_CHECK(esppit2.isDone());
  BOOST_CHECK(!esppit.isDone());

}

