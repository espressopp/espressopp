#include "esutil/BlockVector.hpp"
#include "esutil/Timer.hpp"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <boost/foreach.hpp>

using namespace std;
using namespace espresso;
using namespace espresso::esutil;

typedef vector<long long> VectorClass;
typedef BlockVector<VectorClass> BVClass;

const int nblocks  = 10000;
const int maxblock = 100;

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif

void bm_fill(std::vector<VectorClass> &vec) NOINLINE;
void bm_fill(std::vector<VectorClass> &vec)
{
  srand48(42);
  vec.resize(nblocks);
  BOOST_FOREACH(VectorClass &v, vec) {
    int size = lrand48() % maxblock;
    for (int i = 0; i < size; ++i) {
      v.push_back(lrand48());
    }
  }
}

void bm_fill(BVClass &vec) NOINLINE;
void bm_fill(BVClass &vec)
{
  srand48(42);
  vec.resize(nblocks);
  BOOST_FOREACH(BVClass::Block v, vec) {
    int size = lrand48() % maxblock;
    for (int i = 0; i < size; ++i) {
      v.push_back(lrand48());
    }
  }
}

void bm_shuffle(std::vector<VectorClass> &vec) NOINLINE;
void bm_shuffle(std::vector<VectorClass> &vec)
{
  srand48(42);
  for (size_t c = 0; c < 10000000; ++c) {
    size_t b1 = lrand48() % nblocks;
    if (vec[b1].empty()) continue;
    size_t p1 = lrand48() % vec[b1].size();
    size_t b2 = lrand48() % nblocks;
    int val = vec[b1][p1];
    vec[b2].push_back(val);
    vec[b1].erase(vec[b1].begin() + p1);      
  }
}

void bm_shuffle(BVClass &vec) NOINLINE;
void bm_shuffle(BVClass &vec) {
  srand48(42);
  for (size_t c = 0; c < 10000000; ++c) {
    size_t b1 = lrand48() % nblocks;
    if (vec[b1].empty()) continue;
    size_t p1 = lrand48() % vec[b1].size();
    size_t b2 = lrand48() % nblocks;
    int val = vec[b1][p1];
    vec[b2].push_back(val);
    vec[b1].erase(vec[b1].begin() + p1);      
  }
}

BVClass::Block bm_test_pushback(BVClass &vec) NOINLINE;
BVClass::Block bm_test_pushback(BVClass &vec) {
  return vec[2345];
}

VectorClass::iterator bm_test_pushback(std::vector<VectorClass> &vec) NOINLINE;
VectorClass::iterator bm_test_pushback(std::vector<VectorClass> &vec) {
  return vec[2345].begin();
}

int main()
{
  vector<VectorClass> stdvec;
  BVClass blkvec;

  WallTimer timer;

  /********************************************
   * fill BlockVector of vector
   ********************************************/
  {
    timer.reset();
    bm_fill(blkvec);
    float time = timer.getElapsedTime();
    size_t s = 0, c = 0;
    BOOST_FOREACH(BVClass::Block v, blkvec) { s += v.size(); c += v.capacity(); }
    cout << "filling BlockVector with " << s << " elements took " << time << " seconds" << endl;
    cout << "capacity is " << c << " elements" << endl;
    cout << "gap size " << blkvec.target_gap_size() << endl;
  }

  /********************************************
   * fill STL vector of vector
   ********************************************/
  {
    timer.reset();
    bm_fill(stdvec);
    float time = timer.getElapsedTime();
    size_t s = 0, c = 0;
    BOOST_FOREACH(VectorClass &v, stdvec) { s += v.size(); c += v.capacity(); }
    cout << "filling std::vector with " << s << " elements took " << time << " seconds" << endl;
    cout << "capacity is " << c << " elements" << endl;
  }

  /********************************************
   * random shuffling between neighbor blocks - vector
   ********************************************/
  {
    timer.reset();
    bm_shuffle(stdvec);
    float time = timer.getElapsedTime();
    size_t s = 0, c = 0;
    BOOST_FOREACH(VectorClass &v, stdvec) { s += v.size(); c += v.capacity(); }
    cout << "shuffling std::vector took " << time << " seconds" << endl;    
    cout << "final usage is still " << s << " out of a capacity of " << c << " elements" << endl;
  }

  /********************************************
   * random shuffling between neighbor blocks - BlockVector
   ********************************************/
  {
    timer.reset();
    bm_shuffle(blkvec);
    float time = timer.getElapsedTime();
    size_t s = 0, c = 0;
    BOOST_FOREACH(BVClass::Block v, blkvec) { s += v.size(); c += v.capacity(); }
    cout << "shuffling BlockVector took " << time << " seconds" << endl;    
    cout << "final usage is still " << s << " out of a capacity of " << c << " elements" << endl;
    cout << "gap size " << blkvec.target_gap_size() << endl;
  }

  return 0;
}
