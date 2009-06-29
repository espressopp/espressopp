#include <boost/python.hpp>
#include <iostream>

using namespace std;
using namespace boost::python;

class A {

public:
  virtual void f() = 0;

  virtual void g() {
    cout << "This is g, and I'm calling my f!" << endl;
    f();
  }

};

// wrapper for class A
struct PyA : A,  wrapper< A > {
  void f() {
    this->get_override("f")();
  }
};


// dervied class
class B : public A {
public:
  void f() {
    cout << "I'm B::f()!" << endl;
  }
};

// wrapper for class B
struct PyB : B, wrapper< B > {
  void f() {
    if (override f = this->get_override("f")) {
      f(); 
      return;
    }
    B::f();
  }

  void default_f() {
    this->B::f();
  }
};

void caller(A &a) {
  a.g();
}

BOOST_PYTHON_MODULE(_inheritance)
{
  class_< PyA, boost::noncopyable >("A")
    .def("f", pure_virtual(&A::f))
    .def("g", &A::g)
    ;

  // simple export
  class_< B, bases< A > >("Bsimple")
    .def("f", &B::f)
    .def("g", &B::g)
    ;

  // wrapped export
  class_< PyB, bases< A >, boost::noncopyable >("B")
    .def("f", &B::f, &PyB::default_f)
    .def("g", &B::g)
    ;

  def("caller", &caller);
}


