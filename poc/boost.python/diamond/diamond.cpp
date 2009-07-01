#include <boost/python.hpp>
#include <iostream>

using namespace std;
using namespace boost::python;

class A 
{
public:
  A() {
    cout << "C++ A constructor" << endl;
  }
  
  virtual void doIt2() {
    cout << "C++ A doIt2" << endl;
  }

  virtual void doIt() {
    cout << "C++ A" << endl;
  }

  virtual ~A() {
    cout << "C++ A destructor" << endl;
  }
};

class B 
  : public A 
{
public:
  B() {
    cout << "C++ B constructor" << endl;
  }

  virtual void doIt() {
    this->A::doIt();
    cout << "C++ B" << endl;
  }
};

BOOST_PYTHON_MODULE(_diamond)
{
  class_< A, boost::noncopyable >("A")
    .def("doIt", &A::doIt)
    .def("doIt2", &A::doIt2)
    ;

  class_< B, bases< A >, boost::noncopyable >("B")
    .def("doIt", &B::doIt)
    ;
}
