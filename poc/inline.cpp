#include <iostream>

#include <esutil/Timer.hpp>

using namespace std;
using namespace espresso;
using namespace esutil;

template< typename type >
class Test {
  type *base;
  int size;

public:
  Test(type *_base, int _size): base(_base), size(_size) {}

  int getSizeInline() const { return size; }

  type *getInline() const { return base; }
  
  int getSize() const __attribute__((noinline));
  type *get() const __attribute__((noinline));

};

template< typename type >
type *Test< type >::get() const
{
  return base;
}

template< typename type >
int Test< type >::getSize() const
{
  return size;
}

template<typename type>
void runTestInline(Test<type> &t, int rounds) __attribute__((noinline));
template<typename type>
void runTestInline(Test<type> &test, int rounds)
{
  for (int r = 1; r <= rounds; ++r)
    for (int i = 0; i < test.getSizeInline(); ++i) {
      test.getInline()[i] = i + r;
    }
}

template<typename type>
void runTest(Test<type> &t, int rounds) __attribute__((noinline));
template<typename type>
void runTest(Test<type> &test, int rounds)
{
  for (int r = 1; r <= rounds; ++r)
    for (int i = 0; i < test.getSize(); ++i) {
      test.get()[i] = i + r;
    }
}

int main()
{
  const int size = 1000000;
  const int rounds = 100;

  //typedef double type;
  typedef int type;
  //typedef float type;
  Test< type > test(new type[size], size);

  WallTimer timer;

  timer.reset();
  runTestInline(test, rounds);
  cout << "inline1 " << timer.getElapsedTime() << endl;

  timer.reset();
  runTestInline(test, rounds);
  cout << "inline2 "  << timer.getElapsedTime() << endl;

  timer.reset();
  runTest(test, rounds);
  cout << "call1 " << timer.getElapsedTime() << endl;

  timer.reset();
  runTest(test, rounds);
  cout << "call2 " << timer.getElapsedTime() << endl;

  return 0;
}
