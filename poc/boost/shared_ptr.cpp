#include "types.hpp"

#include <iostream>

using namespace std;
using namespace espresso;

class A {
public:
  A() { cout << "A constructor" << endl; }

  ~A() { cout << "A destructor" << endl; }
};

int main() {
  shared_ptr< A > a_sp = make_shared< A >();

  A &aref = *a_sp;

  shared_ptr< A > a_sp2 = shared_ptr< A >(&aref);

  assert(a_sp == a_sp2);
  
  // this causes a segfault! 
  // Don't mix shared_ptr and std c pointers!
}
