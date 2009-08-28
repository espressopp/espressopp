#define BOOST_TEST_MODULE BlockVector
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <iostream>

using namespace std;

#include "../BlockVector.hpp"
#include "../TupleVector.hpp"

typedef espresso::esutil::BlockVector< vector<int> > BlockVector;

/* make this a template to allow block to be both a reference and not.
   Unfortunately, C++, does not allow non-const references to temporaries... */
template<class BlockClass>
void fill_block(BlockClass block, int start, int step, size_t size) {
  block.resize(size);
  int i = start;
  BOOST_FOREACH(int &val, block) {
    val = i;
    i += step;
  }
}

/* implement as a macro so that boost::test does not output the line here, but
   the one where the problem really occurs */
#define check_block(block, start, step, sz)                     \
  {                                                             \
    BlockVector::Block check_block_b = block;                   \
    int check_block_i = start;                                  \
    BOOST_CHECK_EQUAL(check_block_b.size(), size_t(sz));        \
    BOOST_FOREACH(int &val, check_block_b) {                    \
      BOOST_CHECK_EQUAL(val, check_block_i);                    \
      check_block_i += step;                                    \
    }                                                           \
  }

BOOST_AUTO_TEST_CASE(construct_test)
{
  // check that empty vector is empty
  BlockVector blockv;
  BOOST_CHECK_EQUAL(blockv.size(), size_t(0));
  BOOST_CHECK(blockv.empty());

  // check that blocks are empty from the start
  blockv = BlockVector(3);

  BOOST_CHECK_EQUAL(blockv.size(), size_t(3));

  BOOST_FOREACH(BlockVector::Block b, blockv) {
    BOOST_CHECK(b.size() == 0);
  }
  blockv.clear();
  BOOST_CHECK(blockv.empty());  
}

BOOST_AUTO_TEST_CASE(front_back_test)
{
  BlockVector blockv(3);
  fill_block(blockv[0], 23,  7, 2);
  fill_block(blockv[1], 42,  2, 3);
  fill_block(blockv[2], 29,  6, 4);

  check_block(blockv.front(), 23,  7, 2);
  check_block(blockv.back(), 29,  6, 4);
}

BOOST_AUTO_TEST_CASE(add_remove_block_test)
{
  BlockVector blockv(3);

  // setup blocks
  int i = 0;
  BOOST_FOREACH(BlockVector::Block block, blockv) {
    fill_block(block, 10*i, 1, 10*(i*i*i+1));
    i++;
  }
  blockv.resize(20);

  // check original blocks
  i = 0;
  int cnt = 0;
  BOOST_FOREACH(BlockVector::Block block, blockv) {
    if (cnt++ < 3) {
      check_block(block, 10*i, 1, 10*(i*i*i+1));
      i++;
    }
    else {
      BOOST_CHECK_EQUAL(block.size(), size_t(0));
      BOOST_CHECK_EQUAL(block.capacity(), size_t(blockv.target_gap_size()));
    }
  }

  // deletion of the last, large block enforces reallocation
  blockv.resize(2);

  // check remaining blocks
  i = 0;
  BOOST_FOREACH(BlockVector::Block block, blockv) {
    check_block(block, 10*i, 1, 10*(i*i*i+1));
    i++;
  }

  // check that we can clean up the mess
  blockv.clear();
  BOOST_CHECK(blockv.empty());  
  BOOST_CHECK_EQUAL(blockv.data_size(), size_t(0));
}

BOOST_AUTO_TEST_CASE(insert_erase_block_test)
{
  BlockVector blockv(3);

  fill_block(blockv[0], 23,  7, 2);
  fill_block(blockv[1], 42,  2, 3);
  fill_block(blockv[2], 29,  6, 4);

  // insert an empty block
  BlockVector::iterator it = blockv.insert(blockv.begin() + 1);
  // copy blockv[0] before that block
  blockv.insert(it, blockv[0]);
  // insert empty blocks at the end
  blockv.insert(blockv.end(), 2);
  
  BOOST_REQUIRE_EQUAL(blockv.size(), size_t(7));

  check_block(blockv[0], 23,  7, 2);
  check_block(blockv[1], 23,  7, 2);
  BOOST_CHECK_EQUAL(blockv[2].size(), size_t(0));
  check_block(blockv[3], 42,  2, 3);
  check_block(blockv[4], 29,  6, 4);
  BOOST_CHECK_EQUAL(blockv[5].size(), size_t(0));
  BOOST_CHECK_EQUAL(blockv[6].size(), size_t(0));

  // now restore
  it = blockv.erase(blockv.begin() + 1);

  BOOST_REQUIRE_EQUAL(blockv.size(), size_t(6));

  check_block(blockv[0], 23,  7, 2);
  BOOST_CHECK_EQUAL(blockv[1].size(), size_t(0));
  check_block(blockv[2], 42,  2, 3);
  check_block(blockv[3], 29,  6, 4);
  BOOST_CHECK_EQUAL(blockv[4].size(), size_t(0));
  BOOST_CHECK_EQUAL(blockv[5].size(), size_t(0));

  blockv.erase(blockv.begin() + 4, blockv.end());

  BOOST_REQUIRE_EQUAL(blockv.size(), size_t(4));

  check_block(blockv[0], 23,  7, 2);
  BOOST_CHECK_EQUAL(blockv[1].size(), size_t(0));
  check_block(blockv[2], 42,  2, 3);
  check_block(blockv[3], 29,  6, 4);
}

BOOST_AUTO_TEST_CASE(block_resize_test)
{
  BlockVector blockv(5);

  // init buffer blocks, in random order
  fill_block(blockv[0], 23,  7, 2);
  fill_block(blockv[1], 42,  2, 3);
  fill_block(blockv[4], 29,  6, 4);
  fill_block(blockv[3], 31, -1, 7);

  // since we only manipulate one block, we can keep its reference
  // - do that for testing
  BlockVector::Block block  = blockv[2];

  fill_block<BlockVector::Block &>(block,  0,  1, 5);
  
  // check content to detect possible overlaps
  check_block(blockv[0], 23,  7, 2);
  check_block(blockv[1], 42,  2, 3);
  check_block(block,  0,  1, 5);
  check_block(blockv[3], 31, -1, 7);
  check_block(blockv[4], 29,  6, 4);

  // now resize the middle block significantly, forcing a reallocation
  block.resize(1000);

  {
    // fill up the resized block
    for (size_t i = 5; i < 1000; ++i) { block[i] = i + 11; }
  }

  // check that the setup worked, and blocks do not e.g. overlap
  BOOST_CHECK_EQUAL(blockv[0].size(), size_t(2));
  BOOST_CHECK_EQUAL(blockv[1].size(), size_t(3));
  BOOST_CHECK_EQUAL(blockv[2].size(), size_t(1000));
  // double check copy and original
  BOOST_CHECK_EQUAL(block.size(), size_t(1000));
  BOOST_CHECK_EQUAL(blockv[3].size(), size_t(7));
  BOOST_CHECK_EQUAL(blockv[4].size(), size_t(4));

  /* check that the capacities match.
     All gaps after the resized blocks should be default size */
  BOOST_CHECK_EQUAL(blockv[2].capacity(), size_t(1000 + blockv.target_gap_size()));
  // double check copy and original
  BOOST_CHECK_EQUAL(block.capacity(), size_t(1000 + blockv.target_gap_size()));
  BOOST_CHECK_EQUAL(blockv[3].capacity(), size_t(7 + blockv.target_gap_size()));
  BOOST_CHECK_EQUAL(blockv[4].capacity(), size_t(4 + blockv.target_gap_size()));

  // now, check blocks again
  {
    // special one has to be check manually
    for (size_t i = 0; i < 5; ++i)    { BOOST_CHECK_EQUAL(block[i], int(i)); }
    for (size_t i = 5; i < 1000; ++i) { BOOST_CHECK_EQUAL(block[i], int(i + 11)); }
  }
  check_block(blockv[0], 23,  7, 2);
  check_block(blockv[1], 42,  2, 3);
  check_block(blockv[3], 31, -1, 7);
  check_block(blockv[4], 29,  6, 4);

  BOOST_CHECK_EQUAL(blockv.data_size(), size_t(1016));
  BOOST_CHECK_EQUAL(blockv.raw().size(), size_t(1116));
}

BOOST_AUTO_TEST_CASE(block_copy_and_assign_test)
{
  BlockVector blockv(3);

  fill_block(blockv[0], 23,  7, 2);
  fill_block(blockv[2], 42,  2, 3);
  fill_block(blockv[1], 29,  6, 4);

  blockv[1].assign(20, 3);

  blockv[2] = blockv[0];

  check_block(blockv[0], 23,  7, 2);
  check_block(blockv[1],  3,  0, 20);
  check_block(blockv[2], 23,  7, 2);
}

BOOST_AUTO_TEST_CASE(block_push_pop_test)
{
  BlockVector blockv(3);

  fill_block(blockv[0], 15,  3, 2);
  fill_block(blockv[1],  2,  1, 4);
  fill_block(blockv[2], 42,  2, 3);

  blockv[1].push_back(6);

  check_block(blockv[0], 15,  3, 2);
  check_block(blockv[1],  2,  1, 5);
  check_block(blockv[2], 42,  2, 3);

  blockv[1].pop_back();

  check_block(blockv[0], 15,  3, 2);
  check_block(blockv[1],  2,  1, 4);
  check_block(blockv[2], 42,  2, 3);
}

BOOST_AUTO_TEST_CASE(block_insert_test)
{
  BlockVector blockv(3);

  fill_block(blockv[0], 15,  3, 2);
  fill_block(blockv[1],  2,  1, 4);
  fill_block(blockv[2],  8,  1, 10);

  /* to be able to check easily, we simply continue
     filling the middle block with the pattern it
     already has, but get the data from other sources */
  // manual inserts
  blockv[1].insert(blockv[1].begin() + 4, 1, 6);
  blockv[1].insert(blockv[1].end(), 1, 7);
  // copy from the last block
  blockv[1].insert(blockv[1].end(), blockv[2].begin(), blockv[2].begin() + 2);

  check_block(blockv[0], 15,  3, 2);
  check_block(blockv[1],  2,  1, 8);
  check_block(blockv[2],  8,  1, 10);

  blockv[1].pop_back();

  check_block(blockv[0], 15,  3, 2);
  check_block(blockv[1],  2,  1, 7);
  check_block(blockv[2],  8,  1, 10);
}

BOOST_AUTO_TEST_CASE(block_erase_test)
{
  BlockVector blockv(3);

  fill_block(blockv[0], 15,  3, 2);
  fill_block(blockv[1],  2,  1, 10);
  fill_block(blockv[2],  8,  1, 4);

  blockv[1].erase(blockv[1].begin());
  blockv[1].erase(blockv[1].end() - 2, blockv[1].end());

  check_block(blockv[0], 15,  3, 2);
  check_block(blockv[1],  3,  1, 7);
  check_block(blockv[2],  8,  1, 4);
}

/* this test instantiates the tricky parts of the fat/thin reference/iterator business of
   TupleVector. In the first place, this checks that the BlockVector<TupleVector> can be
   instantiated at all.
*/
BOOST_AUTO_TEST_CASE(tuplevector_interop_text)
{
  typedef espresso::esutil::BlockVector< espresso::esutil::TupleVector > BlockVector;

  BlockVector blockv(3);
  blockv[0].resize(10);
  /* this only works if the correct specialization for TupleVector is found, i.e.
     the make_thick mechanism is in use.
   */
  blockv[1].assign(blockv[0].begin(), blockv[0].begin() + 3);
}
