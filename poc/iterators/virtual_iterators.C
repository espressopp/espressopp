// This benchmarks different versions of a loop via an iterator,
// where the iterator methods are virtual or non-virtual

#include <vector>
#include <iostream>
#include <sys/times.h>

using namespace std;

const int N = 10000;
const int NUM_TESTS = 10000;

struct tms start_t, stop_t;
const unsigned clocks_per_sec = sysconf(_SC_CLK_TCK);

// Iterator base class with virtual methods
template <class T>
class VirtualIterator {
public:
  T t;

  virtual bool isDone() { return true; };
  virtual void Next() {};
  virtual const T &GetCurrent() { return t; };
};

// Vector Iterator, descendant of the Iterator base class
// uses virtual methods
template <class T>
class VirtualVectorIterator 
  : public VirtualIterator<T> {
public:
  int i;
  const vector <T> &v;

  // constructor
  VirtualVectorIterator(vector <T> &v) 
    : v(v)
  { i = 0; }

  virtual bool isDone() {
    return i >= v.size();
  }

  virtual void Next() {
    i++;
  }

  virtual const T &GetCurrent() {
    return v[i];
  }
};

// Vector Iterator that does not use virtual methods
template <class T>
class VectorIterator {
public:
  int i;
  const vector <T> &v;

  // constructor
  VectorIterator(vector <T> &v) 
    : v(v)
  { i = 0; }

  bool isDone() {
    return i >= v.size();
  }

  void Next() {
    i++;
  }

  const T &GetCurrent() {
    return v[i];
  }
};

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
    for (unsigned int i=0; i < 3; i++) 
      (*this)[i] += v[i];
    return *this;
  }
};

ostream &operator<<(ostream& out, 
		    const Vector3D& v) {
  out << v[0] << '\t'
      << v[1] << '\t'
      << v[2];
}

typedef Vector3D TestType;

int main() {
  cout << "Creating random vector with " << N << " elements..." << endl;
	
  vector<TestType> v;
  for (int i=0; i < N; i++)
    v.push_back(Vector3D(drand48(),drand48(),drand48()));

  int i;
  TestType sum;

  // benchmark an iterator that uses virtual methods
  // and is called via the base class
  VirtualIterator<TestType> *it = NULL;
  cout << "Testing iterators..." << endl;
  sum = Vector3D();
  times(&start_t);
  for (i=0; i < NUM_TESTS; i++) {
    for (it = new VirtualVectorIterator<TestType>(v); 
	 !it->isDone(); it->Next())
      sum += it->GetCurrent();
    delete it;
  }
  times(&stop_t);
  cout << sum << endl;
  cout << "user time = " 
       << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
       << "s\tsystem time = "
       << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
       << "s" << endl;

  // second repetition of the benchmark
  VirtualIterator<TestType> *it2 = NULL;
  sum = Vector3D();
  times(&start_t);
  for (i=0; i < NUM_TESTS; i++) {
    for (it2 = new VirtualVectorIterator<TestType>(v); 
	 !it2->isDone(); it2->Next())
      sum += it2->GetCurrent();
    delete it2;
  }
  times(&stop_t);
  cout << sum << endl;
  cout << "user time = " 
       << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
       << "s\tsystem time = "
       << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
       << "s" << endl;

  // benchmark an iterator that has virtual methods
  // but that is called directly
  VirtualVectorIterator<TestType> *vit = NULL;
  sum = Vector3D();
  times(&start_t);
  for (i=0; i < NUM_TESTS; i++) {
    for (vit = new VirtualVectorIterator<TestType>(v); 
	 !vit->isDone(); vit->Next())
      sum += vit->GetCurrent();
    delete vit;
  }
  times(&stop_t);
  cout << sum << endl;
  cout << "user time = " 
       << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
       << "s\tsystem time = "
       << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
       << "s" << endl;

  // benchmark a non-virtual iterator
  VectorIterator<TestType> *vit2 = NULL;
  sum = Vector3D();
  times(&start_t);
  for (i=0; i < NUM_TESTS; i++) {
    for (vit2 = new VectorIterator<TestType>(v); 
	 !vit2->isDone(); vit2->Next())
      sum += vit2->GetCurrent();
    delete vit2;
  }
  times(&stop_t);
  cout << sum << endl;
  cout << "user time = " 
       << static_cast<double>(stop_t.tms_utime - start_t.tms_utime) / clocks_per_sec
       << "s\tsystem time = "
       << static_cast<double>(stop_t.tms_stime - start_t.tms_stime) / clocks_per_sec
       << "s" << endl;


}
