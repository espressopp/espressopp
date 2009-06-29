// This benchmarks different versions of a loop via an iterator,
// where the iterator methods are virtual or non-virtual

#include <vector>
#include <iostream>
#include <sys/times.h>

#include "iterators_virtual.H"
#include "iterators_nonvirtual.H"
#include "iterators_block.H"
#include "iterators_fobject.H"

using namespace std;

// Test parameters
const int N = 10000;
const int NUM_TESTS = 10000;

/* whether to use a constant size loop in the vector class
   or manually unroll it. This makes a difference for gcc,
   even if unrolling loops */
const bool manual_unroll = false;

//////////////////////////////////////////////////
// DEFINE A 3D REAL DATA TYPE FOR THE TESTS
//////////////////////////////////////////////////
class Vector3D {
  double data[3];

public:
  Vector3D(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }

  Vector3D() { 
    data[0] = 0.0;
    data[1] = 0.0;
    data[2] = 0.0;
  }

  Vector3D(double v) {
    data[0] = v;
    data[1] = v;
    data[2] = v;
  }

  double operator[](unsigned int i) const {
    return data[i];
  }

  double &operator[](unsigned int i) {
    return data[i];
  }

  Vector3D& operator+=(const Vector3D &v) {
    if (manual_unroll) {
      (*this)[0] += v[0];
      (*this)[1] += v[1];
      (*this)[2] += v[2];
    }
    else {
      for (int i = 0; i < 3; ++i)
	(*this)[i] += v[i];
    }
    return *this;
  }
};

ostream &operator<<(ostream& out, 
		    const Vector3D& v) {
  out << v[0] << '\t'
      << v[1] << '\t'
      << v[2];
  return out;
}

typedef Vector3D TestType;


//////////////////////////////////////////////////
// DEFINE THE GLOBAL VARIABLES AND FUNCTIONS
//////////////////////////////////////////////////
struct tms start_t, stop_t;
const unsigned clocks_per_sec = sysconf(_SC_CLK_TCK);
vector<TestType> v;
TestType global_sum;

void benchmark_start() {
  times(&start_t);
}

#ifdef __GNUC__
void benchmark_stop(TestType &sum)
  __attribute__((noinline));
#endif
void benchmark_stop(TestType &sum) {
  times(&stop_t);
  cout << "user time = " 
       << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
       << "s\tsystem time = "
       << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
       << "s" << endl;
  global_sum += sum;
}


//////////////////////////////////////////////////
// DEFINE THE VARIOUS BENCHMARKS
//////////////////////////////////////////////////
// LOOPING A VIRTUAL ITERATOR VIA THE BASE 
#ifdef __GNUC__
TestType sum_virtual_base (VirtualIterator<TestType> *it)
  __attribute__((noinline));
#endif    

TestType sum_virtual_base (VirtualIterator<TestType> *it)
{
  TestType sum;
  for (it->reset(); !it->isDone(); it->next())
    sum += it->getCurrent();
  return sum;
}

void benchmark_virtual_base() {
  cout << "Benchmarking virtual via base..." << endl;
  TestType sum;

  benchmark_start();
  VirtualIterator<TestType> *it = 
    new VirtualVectorIterator<TestType>(v);
  for (size_t i=0; i < NUM_TESTS; i++) {
    sum += sum_virtual_base(it);
  }
  delete it;
  benchmark_stop(sum);
}

//////////////////////////////////////////////////
// LOOPING A VIRTUAL ITERATOR VIA A POINTER TO THE DERIVED CLASS
#ifdef __GNUC__
TestType sum_virtual_derived (VirtualVectorIterator<TestType> *it)
  __attribute__((noinline));
#endif    

TestType sum_virtual_derived (VirtualVectorIterator<TestType> *it)
{
  TestType sum;
  for (it->reset(); !it->isDone(); it->next())
    sum += it->getCurrent();
  return sum;
}

void benchmark_virtual_derived() {
  cout << "Benchmarking virtual via derived class pointer ..." << endl;
  TestType sum;

  benchmark_start();
  VirtualVectorIterator<TestType> *it 
    = new VirtualVectorIterator<TestType>(v);
  for (size_t i=0; i < NUM_TESTS; i++) {
    sum += sum_virtual_derived(it);
  }
  delete it;
  benchmark_stop(sum);
}

//////////////////////////////////////////////////
// LOOPING A VIRTUAL ITERATOR VIA A REFERENCE TO THE DERIVED CLASS
#ifdef __GNUC__
TestType sum_virtual_derived_ref (VirtualVectorIterator<TestType> &it)
  __attribute__((noinline));
#endif

TestType sum_virtual_derived_ref (VirtualVectorIterator<TestType> &it)
{
  TestType sum;
  for (it.reset(); !it.isDone(); it.next())
    sum += it.getCurrent();
  return sum;
}

void benchmark_virtual_derived_ref() {
  cout << "Benchmarking virtual derived class reference ..." << endl;
  TestType sum;

  benchmark_start();
  VirtualVectorIterator<TestType> it(v);
  for (size_t i=0; i < NUM_TESTS; i++) {
    sum += sum_virtual_derived_ref(it);
  }
  benchmark_stop(sum);
}

//////////////////////////////////////////////////
// LOOPING A VIRTUAL ITERATOR LOCALLY
#ifdef __GNUC__
TestType sum_virtual_local (VirtualVectorIterator<TestType> &_it)
  __attribute__((noinline));
#endif

TestType sum_virtual_local (VirtualVectorIterator<TestType> &_it)
{
  TestType sum;
  VirtualVectorIterator<TestType> it(_it);

  for (it.reset(); !it.isDone(); it.next())
    sum += it.getCurrent();
  return sum;
}

void benchmark_virtual_local() {
  cout << "Benchmarking local virtual iterator..." << endl;
  TestType sum;

  benchmark_start();
  VirtualVectorIterator<TestType> it(v);
  for (size_t i=0; i < NUM_TESTS; i++) {
    sum += sum_virtual_local(it);
  }
  benchmark_stop(sum);
}

//////////////////////////////////////////////////
// LOOPING OF A NONVIRTUAL ITERATOR
#ifdef __GNUC__
TestType sum_direct (vector<TestType> &v)
  __attribute__((noinline));
#endif

TestType sum_direct (vector<TestType> &v)
{
  TestType sum;
  std::vector<TestType>::iterator it;
  for(it = v.begin(); it != v.end(); it++) {
    sum += (*it);
  }
  return sum;
}

void benchmark_direct() {
  cout << "Benchmarking non-virtual iterator..." << endl;
  TestType sum;

  benchmark_start();
  for (size_t i=0; i < NUM_TESTS; i++) {
    sum += sum_direct(v);
  }
  benchmark_stop(sum);
}


//////////////////////////////////////////////////
// LOOPING OF A NONVIRTUAL ITERATOR, NOT USING VECTOR3D IN SUM
#ifdef __GNUC__
Vector3D sum_direct_noVector3D_sum (vector<Vector3D> &v)
  __attribute__((noinline));
#endif

Vector3D sum_direct_noVector3D_sum (vector<Vector3D> &v)
{
  double sx = 0, sy = 0, sz = 0;
  std::vector<Vector3D>::const_iterator it;
  for(it = v.begin(); it != v.end(); it++) {
    sx += (*it)[0];
    sy += (*it)[1];
    sz += (*it)[2];
  }
  return Vector3D(sx, sy, sz);
}

void benchmark_direct_noVector3D_sum() {
  cout << "Benchmarking non-virtual iterator without Vector3D sum..." << endl;
  Vector3D sum;

  benchmark_start();
  for (size_t i=0; i < NUM_TESTS; i++) {
    sum += sum_direct_noVector3D_sum(v);
  }
  benchmark_stop(sum);
}


//////////////////////////////////////////////////
// LOOPING OF A NONVIRTUAL ITERATOR, NOT USING VECTOR3D AT ALL
#ifdef __GNUC__
Vector3D sum_direct_double3 (vector<Vector3D> &v)
  __attribute__((noinline));
#endif

Vector3D sum_direct_double3 (vector<Vector3D> &v)
{
  double sum[3] = {0, 0, 0};
  std::vector<Vector3D>::iterator it;
  for(it = v.begin(); it != v.end(); it++) {
      sum[0] += (*it)[0];
      sum[1] += (*it)[1];
      sum[2] += (*it)[2];
  }
  return Vector3D(sum[0], sum[1], sum[2]);
}

void benchmark_direct_double3() {
  // benchmark using a double[3] in inner loop 
  cout << "Benchmarking double[3]..." << endl;
  Vector3D sum;

  benchmark_start();
  for (size_t i=0; i < NUM_TESTS; i++) {
    sum += sum_direct_double3(v);
  }
  benchmark_stop(sum);
}

//////////////////////////////////////////////////
// LOOPING OF A NONVIRTUAL ITERATOR WITH AN INLINED LOOP
void benchmark_direct_inline() {
  cout << "Benchmarking inlined non-virtual iterator..." << endl;

  benchmark_start();
  Vector3D sum;

  benchmark_start();
  {
    Vector3D lsum;
    std::vector<Vector3D>::const_iterator it;
    for (size_t i=0; i < NUM_TESTS; i++) {
      for(it = v.begin(); it != v.end(); it++) {
        lsum += *it;
      }
    }
    sum = lsum;
  }
  benchmark_stop(sum);
}

//////////////////////////////////////////////////
// LOOPING USING A FUNCTION OBJECT
class SumUp: public Computer<TestType> {
  TestType res;

public:
  virtual void operator()(const TestType &v) {
    res += v;
  }
  operator TestType() const { return res; }
};

#ifdef __GNUC__
TestType sum_fobject (VirtualForeach<TestType> &)
  __attribute__((noinline));
#endif

TestType sum_fobject (VirtualForeach<TestType> &f)
{
  SumUp sum;

  f.foreach(sum);
 
  return sum;
}

void benchmark_fobject() {
  // benchmark using a double[3] in inner loop 
  cout << "Benchmarking using a function object..." << endl;
  TestType sum;
  VirtualVectorForeach<TestType> f(v);

  benchmark_start();
  for (size_t i=0; i < NUM_TESTS; i++) {
    sum += sum_fobject(f);
  }
  cout << "sum = " << sum << endl;
  benchmark_stop(sum);
}

//////////////////////////////////////////////////
// LOOPING OF A BLOCKED VIRTUAL ITERATOR
template < int blocksize, int used >
void benchmark_block() {
  // benchmark a blocked virtual iterator
  cout << "Benchmarking blocked with virtual filler at size " 
       << blocksize << " using " << used << "..." << endl;
  TestType sum;

  benchmark_start();
  times(&start_t);
  for (size_t i=0; i < NUM_TESTS; i++) {
    VirtualBlockFiller<TestType, blocksize> *filler = 
      new VirtualVectorBlockFiller<TestType, blocksize>(v, used);
    for (BlockedIterator<TestType, blocksize> it(filler);
	 !it.isDone(); it.next())
      sum += it.getCurrent();
    delete filler;
  }
  cout << "expected sum = " << sum << endl;
  benchmark_stop(sum);
}

//////////////////////////////////////////////////
// MAIN PROGRAM
//////////////////////////////////////////////////
int main() {

  cout << "Creating random vector with " << N << " elements..." << endl;
	
  for (size_t i=0; i < N; i++)
    v.push_back(TestType(drand48(),drand48(),drand48()));

  int i;
  TestType sum;

  benchmark_virtual_base();
  benchmark_virtual_derived();
  benchmark_virtual_derived_ref();
  benchmark_virtual_local();

  benchmark_direct();
  benchmark_direct_noVector3D_sum();
  benchmark_direct_double3();
  benchmark_direct_inline();

  benchmark_fobject();

  benchmark_block<64, 64>();
  benchmark_block<8, 8>();

  cout << global_sum << endl;

}
