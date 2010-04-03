#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
#include <iostream>
#include <string>

using namespace std;
using namespace boost;
using namespace boost::python;

class A {
public:
  A() {
    s = "Hello World!";
  }
  A(const A &other) {
    s = other.s;
  }
  string s;
};

class Container {
  A a;
public:
  Container() {
    cout << "Created container" << endl;
  }
  
  void setA(shared_ptr< A > _a) {
    a = *_a;
  }

  A& getA() {
    return a;
  }

  ~Container() {
    cout << "Destroyed container" << endl;
  }
};

BOOST_PYTHON_MODULE(_call_policies)
{
  class_< A >("A")
    .def_readwrite("s", &A::s);

  class_< Container, boost::shared_ptr< Container > >("Container")
    .def("setA", &Container::setA)
    .def("getA", &Container::getA, 
	 return_internal_reference<>())
    .def("getADanger", &Container::getA, 
	 return_value_policy<reference_existing_object>())
    ;
}
