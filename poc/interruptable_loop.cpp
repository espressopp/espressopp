/** Tests the best way to allow loop to be interrupted, either via an
    exception, or via a bool return value of apply.
*/

#include <iostream>
#include <esutil/Timer.hpp>

using namespace std;
using namespace espresso::esutil;

//////////////////////////////////////////////////
// Return whether to interrupt loop or not
//////////////////////////////////////////////////
struct BoolInterruptComputer {
  virtual const char* name() = 0;
  virtual bool apply(int i) = 0;
  virtual void setI(int _i) = 0;
};

struct BoolHalfSquare : BoolInterruptComputer {
  int i;
  double sum;

  BoolHalfSquare() { sum = 0.0; }

  const char* name() { return "BoolHalfSquare"; }

  void setI(int _i) { i = _i; }

  bool apply(int j) {
    if (i == j) return false;
    //    cout << i << " " << j << endl;
    sum += j*0.1;
    return true;
  }
};

struct BoolFullSquare : BoolInterruptComputer {
  double sum;

  BoolFullSquare() { sum = 0.0; }

  const char* name() { return "BoolFullSquare"; }

  void setI(int _i) {}

  bool apply(int j) {
    sum += j*0.1;
    return true;
  }
};


#ifdef __GNUC__
void foreach(int size, BoolInterruptComputer &computer) 
  __attribute__((noinline));
#endif

void foreach(int size, BoolInterruptComputer &computer) {
  for (int i = 0; i < size; i++) {
    computer.setI(i);
    for (int j = 0; j < size; j++) {
      if (!computer.apply(j)) break;
    }
  }
}

//////////////////////////////////////////////////
// Throw exception to interrupt loop
//////////////////////////////////////////////////
struct InterruptLoop {};
struct ExceptionInterruptComputer {
  virtual const char* name() = 0;
  virtual void apply(int i) = 0;
  virtual void setI(int _i) = 0;
};

struct ExceptionHalfSquare : public ExceptionInterruptComputer {
  int i;
  double sum;
  ExceptionHalfSquare() { sum = 0.0; }

  const char* name() { return "ExceptionHalfSquare"; }

  void setI(int _i) { i = _i; }

  virtual void apply(int j) {
    if (i == j) throw InterruptLoop();
    //    cout << i << " " << j << endl;
    sum += j*0.1;
  }
};

struct ExceptionFullSquare : ExceptionInterruptComputer {
  double sum;

  ExceptionFullSquare() { sum = 0.0; }

  const char* name() { return "ExceptionFullSquare"; }

  void setI(int _i) {}

  void apply(int j) {
    sum += 0.1*j;
  }
};

#ifdef __GNUC__
void foreach(int size, ExceptionInterruptComputer &computer)
  __attribute__((noinline));
#endif

void foreach(int size, ExceptionInterruptComputer &computer) {
  for (int i = 0; i < size; i++) {
    computer.setI(i);
    for (int j = 0; j < size; j++) {
      try {
	computer.apply(j);
      } catch (InterruptLoop &er) {
	break;
      }
    }
  }
}


template < class Computer, int num_tests, int square_size >
void foreachTest() {
  WallTimer timer;
  Computer computer;

  timer.reset();
  for (int i = 0; i < num_tests; i++)
    foreach(square_size, computer);
  cout << "  " << timer << " for " << num_tests << " tests of a " 
       << square_size << "x" << square_size << " square (" 
       << computer.name() << ")" << endl;
}

int main() {
  cout << "Tesing the interrupted loop (half square):" << endl;

  foreachTest< BoolHalfSquare, 100000, 100 >();
  foreachTest< ExceptionHalfSquare, 100000, 100 >();

  foreachTest< BoolHalfSquare, 10, 10000 >();
  foreachTest< ExceptionHalfSquare, 10, 10000 >();

  cout << "Tesing the uninterrupted loop (full square):" << endl;
  foreachTest< BoolFullSquare, 100000, 100 >();
  foreachTest< ExceptionFullSquare, 100000, 100 >();

  foreachTest< BoolFullSquare, 10, 10000 >();
  foreachTest< ExceptionFullSquare, 10, 10000 >();

  return 0;
}
