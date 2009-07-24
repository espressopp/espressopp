#include <iostream>
#include <esutil/Timer.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

using namespace std;
using namespace espresso;
using namespace esutil;

class Computer {
public:
  virtual ~Computer() {}
  virtual void apply(float v) = 0;
};

class TheComputer: public Computer {
public:
  TheComputer(): sum(0) {}

  void nonVirtualApply(float v) { sum += v; };

  virtual void apply(float v) { nonVirtualApply(v); };

  float get() { return sum; }

private:
  float sum;
};

typedef boost::function< void (float) > ApplyFunction;

float runTestDirect(int rounds) __attribute__ ((noinline));
float runTestDirect(int rounds) {
  TheComputer test;
  for (int i = 0; i < rounds; ++i) {
    test.apply(i);
  }
  return test.get();
}

float runTestVirtual(int rounds) __attribute__ ((noinline));
float runTestVirtual(int rounds) {
  Computer *test = new TheComputer;
  for (int i = 0; i < rounds; ++i) {
    test->apply(i);
  }
  return static_cast<TheComputer *>(test)->get();
}

float runTestFunction(int rounds) __attribute__ ((noinline));
float runTestFunction(int rounds) {
  TheComputer test;
  ApplyFunction apply = boost::bind(&TheComputer::nonVirtualApply, &test, _1);
  for (int i = 0; i < rounds; ++i) {
    apply(i);
  }
  return test.get();
}

float runTestFunctionVirtual(int rounds) __attribute__ ((noinline));
float runTestFunctionVirtual(int rounds) {
  TheComputer test;
  ApplyFunction apply = boost::bind(&Computer::apply, &test, _1);
  for (int i = 0; i < rounds; ++i) {
    apply(i);
  }
  return test.get();
}

void doFunctionIndirect(int rounds, ApplyFunction apply) __attribute__ ((noinline));
void doFunctionIndirect(int rounds, ApplyFunction apply) {
  for (int i = 0; i < rounds; ++i) {
    apply(i);
  }
}

float runTestFunctionIndirect(int rounds) __attribute__ ((noinline));
float runTestFunctionIndirect(int rounds) {
  TheComputer test;
  doFunctionIndirect(rounds, boost::bind(&TheComputer::nonVirtualApply, &test, _1));
  return test.get();
}

int main()
{
  const int rounds = 100000000;

  WallTimer timer;
  float res;

  timer.reset();
  res = runTestDirect(rounds);
  cout << "direct   " << res << " per operation " << 1e9*timer.getElapsedTime()/rounds << "ns" << endl;

  timer.reset();
  res = runTestVirtual(rounds);
  cout << "virtual  " << res << " per operation " << 1e9*timer.getElapsedTime()/rounds << "ns" << endl;

  timer.reset();
  res = runTestFunction(rounds);
  cout << "function " << res << " per operation " << 1e9*timer.getElapsedTime()/rounds << "ns" << endl;

  timer.reset();
  res = runTestFunctionVirtual(rounds);
  cout << "funcvirt " << res << " per operation " << 1e9*timer.getElapsedTime()/rounds << "ns" << endl;

  timer.reset();
  res = runTestFunctionIndirect(rounds);
  cout << "funcindr " << res << " per operation " << 1e9*timer.getElapsedTime()/rounds << "ns" << endl;

  return 0;
}
