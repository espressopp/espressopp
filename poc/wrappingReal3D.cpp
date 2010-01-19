#include <iostream>
#include <stdexcept>
#include <esutil/Timer.hpp>

using namespace std;
using namespace espresso;
using namespace esutil;

typedef double real;

template<int size>
class Real {
public:
  Real(real _v[size]): v(_v) {}
  operator real*() const { return v; }
private:
  real *v;
};

class WrapReal3 {
public:
  real operator()(Real<3> a, Real<3> b) __attribute__((noinline));
};

real WrapReal3::operator()(Real<3> a, Real<3> b) {
  real res = 0;
  for (int i = 0; i < 3; ++i) {
    res += a[i]*b[i];
  }
  return res;
}

class WrapRealP {
public:
  real operator()(real *a, real *b) __attribute__((noinline));
};

real WrapRealP::operator()(real* a, real* b) {
  real res = 0;
  for (int i = 0; i < 3; ++i) {
    res += a[i]*b[i];
  }
  return res;
};

template<class WrapFunc>
real runTest(int rounds, int size)
{
  WrapFunc wrap;
  real data[size+3];
  real res = 0;
  for (int i = 0; i < rounds; ++i) {
    for (int j = 0; j < size; ++j) {
      res += wrap(data + j, data + (42*j % size));
    }
  }
  return res;
}

real testReal3(int rounds, int size) __attribute__((noinline));
real testReal3(int rounds, int size) {
  return runTest< WrapReal3 >(rounds, size);
}

real testRealP(int rounds, int size) __attribute__((noinline));
real testRealP(int rounds, int size) {
  return runTest< WrapRealP >(rounds, size);
}

int main()
{
  const int rounds = 1000;
  const int size = 10000;
  WallTimer timer;

  timer.reset();
  testReal3(rounds, size);
  cout << "wrapped " << timer.getElapsedTime() << endl;

  timer.reset();
  testRealP(rounds, size);
  cout << "generic " << timer.getElapsedTime() << endl;

  return 0;
}
