#include <iostream>
#include <stdexcept>
#include <esutil/Timer.hpp>

using namespace std;
using namespace espresso;
using namespace esutil;

template<int size>
class Real {
  Real(real _v[size]): v(_v) {}
  operator real*() const { return v; }
private:
  real *v;
};

class Wrap<typename in> {
public:
  real operator()(in a, in b) {
    real res = 0;
    for (int i = 0; i < 3; ++i) {
      res += a[i]*b[i];  
    }
    return res;
  }
};

template<class WrapFunc>
void runTestWrapped(int rounds, int size)
{
  real data[(size+1)*3];
  real res = 0;
  for (int i = 0; i < rounds; ++i) {
    for (int j = 0; j < size; ++j) {
      res += WrapFunc(data + j, data + (3*j % size));
    }
  }
  return res;
}

int main()
{
  const int rounds = 1000;
  const int size = 1000000;
  WallTimer timer;

  timer.reset();
  runTest<WrapFunc< Real<3> > >(rounds, size);
  cout << "wrapped " << timer.getElapsedTime() << endl;

  timer.reset();
  runTest<WrapFunc< real *  > >(rounds, size);
  cout << "generic " << timer.getElapsedTime() << endl;

  return 0;
}
