#define PARALLEL_TEST_MODULE FixedPairList
#include "ut.hpp"

#include "FixedPairList.hpp"
#include "storage/Storage.hpp"

using namespace espresso;


// a storage that contains two fixed particles
class MockStorage1 : public Storage {
  

};

// add a bond
// * on the wrong proc
// * twice on different procs
// * twice
// * in the other order
// * 
// delete a bond
// * on the wrong proc
// * 
// move a bond to another proc


BOOST_FIXTURE_TEST_CASE() {
  
}
