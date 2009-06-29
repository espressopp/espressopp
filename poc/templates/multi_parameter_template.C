#include <iostream>
#include <cmath>

using namespace std;

#include "classes.C"


// macro to instantiate classes
#define INSTANTIATE_CLASS(n1,n2,n3)			\
  {							\
    XYZ<Class_##n1, Class_##n2, Class_##n3> x;		\
    sum += x.get();					\
  }

template <class A, class B, class C>
class XYZ {
public:
  double get() {
    A a;
    B b;
    C c;
    register double x;
    for (int i = 0; i < 10; i++) {
      x += a.get();
      if (x - int(x) > 0.5) {
	x += b.get();
	x += c.get();
      }
    }
    return x;
  }
};


int main() {
  double sum = 0.0;

#include "instantiations.C"

  cout << "sum=" << sum << endl;
}
