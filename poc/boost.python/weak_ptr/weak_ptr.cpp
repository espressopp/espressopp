#include <boost/python.hpp>
#include <boost/shared_ptr.hpp> 
#include <boost/make_shared.hpp> 
#include <boost/weak_ptr.hpp> 
#include <iostream>

using namespace boost::python;
using namespace boost;
using namespace std;

class A {
public:
  A(int _v) : v(_v) {}
  int v;
};

weak_ptr< A > stored_a;

// store a shared_ptr
void setWeakPtr(shared_ptr< A > a) {
  stored_a = a;
}

void outputPtr() {
  cout << stored_a.lock()->v << endl;
}

BOOST_PYTHON_MODULE(_weak_ptr)
{    
  class_< A, shared_ptr< A > >("A", init<int>());
  def("setWeakPtr", &setWeakPtr);
  def("outputPtr", &outputPtr);
}
