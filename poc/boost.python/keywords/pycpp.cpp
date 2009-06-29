#include <iostream>
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/python.hpp>
#include <boost/python.hpp>

using namespace std;

// Keywords for arguments have to be prepared in the following way:
//
//     BOOST_PARAMETER_KEYWORD(<namespace_tag>,<keyword name>)
//
// The namespace_tag is an arbitrary keyword, used to group the keyword names
// (it is not necessary, that the keyword names are identical to element names).
//
BOOST_PARAMETER_KEYWORD(pstag, id)
BOOST_PARAMETER_KEYWORD(pstag, posx)
BOOST_PARAMETER_KEYWORD(pstag, posy)
BOOST_PARAMETER_KEYWORD(pstag, posz)

char const* info() {
  return "Proof of concept C++-Python interface version 0.02";
}

class Particle
{
  string id;
  double xpos, ypos, zpos;
public:
  Particle() {id="EMPTY"; xpos=ypos=zpos=0.0;}
  Particle(string s, double x, double y, double z) : id(s),xpos(x),ypos(y),zpos(z) {};
  string getid() {return id;}
  void setid(string s){ this->id=s;}
  double getxpos() {return xpos;}
  double getypos() {return ypos;}
  double getzpos() {return zpos;}
  void setxpos(double r) {xpos=r;}
  void setypos(double r) {ypos=r;}
  void setzpos(double r) {zpos=r;}
  void view() {cout << "(id=\'" << id << "\', posx=" << xpos << ", posy=" << ypos << ", posz=" << zpos << ")" << endl; }

  // The following Makro creates a member function "set" that has 4 arguments
  // one required one (id) and 3 optional ones (posx, posy, posz) with default values (0.0).
  // In the C++ context, this function can be called with keyword arguments:
  //
  //     set(id='something', posx=1, posy=2, posz=3.0)
  //
  // If all keywords are given, the order of the arguments is arbitrary.
  // If arguments are left out, the default values are used:
  //
  //     set(posz=1.3, id='nothing')
  //
  // is a valid call to this function. 
  //
  // Together with the makro BOOST_PYTHON_MODULE and the
  // structure definition "set_fwd" (see below) everything that is said above
  // is also valid for the Python context.
  
  BOOST_PARAMETER_MEMBER_FUNCTION(
				  (void), set, pstag,
				  (required (id, (string)))
                                  (optional (posx, (double), 0.0)
					    (posy, (double), 0.0)
					    (posz, (double), 0.0)
				   ))
  {
    setid(id);
    setxpos(posx);
    setypos(posy);
    setzpos(posz);
  }
};

// The following structure definition is needed to map the C++ member function "set"
// into the Python context. We need 4 definitions here in order to support the
// variable number of arguments

struct set_fwd
{
  template <class A0, class A1, class A2, class A3>
    void operator()
      (boost::type<void>, Particle& self, A0 const& a0, A1 const& a1, A2 const&a2, A3 const& a3) {self.set(a0,a1,a2,a3);}
  template <class A0, class A1, class A2>
    void operator()
      (boost::type<void>, Particle& self, A0 const& a0, A1 const& a1, A2 const&a2) {self.set(a0,a1,a2);}
  template <class A0, class A1>
    void operator()
      (boost::type<void>, Particle& self, A0 const& a0, A1 const& a1) {self.set(a0,a1);}
  template <class A0>
    void operator()
      (boost::type<void>, Particle& self, A0 const& a0) {self.set(a0);}
};

// The following makro does the actual coupling between C++ and Python
// All the member functions of the C++ Object that should be visible in
// Python have to be defined here:
//
//           .def("member function name", reference to member function)
//
// The ".add_property" function allows one to use data elements of the C++ in
// python in a transparent way. In python (assume p is an instance of the Particle class):
//
//                  p.x = 2.5
// is the same as
//                  p.setxpos(2.5)
//
// The names (the constant strings between quotes) that are used in the definitions
// and thus are used in python to acces the methods do not necessarily have to be
// the C++ names. They can be different.

BOOST_PYTHON_MODULE(pycpp)
{

  using namespace boost::python;
  namespace py = boost::parameter::python;
  namespace mpl = boost::mpl;

  def("info", info);
  class_<Particle>("Particle", init<string,double,double,double>())
    .def(init<>())
    .def("getid", &Particle::getid)
    .def("setid", &Particle::setid)
    .def("view", &Particle::view)
    .add_property("id",&Particle::getid, &Particle::setid)
    .add_property("x", &Particle::getxpos, &Particle::setxpos)
    .add_property("y", &Particle::getypos, &Particle::setypos)
    .add_property("z", &Particle::getzpos, &Particle::setzpos)
    .def("set", py::function<
                      set_fwd
                      , mpl::vector<
                           void
	                   , pstag::id(string)
                           , pstag::posx**(double)
	                   , pstag::posy**(double)
	                   , pstag::posz**(double)
                        >
                >()
    );
};

int main(int argc, char* argv[]) {

  // Due to the BOOST_PARAMETER_MEMBER_FUNCTION makro (see above)
  // keyword arguments can also be used in C++ in the same way as
  // they would be used in Python:
  // 
  // Particle q;
  //
  // cout << info() << endl;
  // q.view();
  // q.set(id="hello", posy=1.3);
  // q.view();

  // And now lets call a Python-Script from C++:

  if (argc > 1) {
    
    using namespace boost::python;

    Py_Initialize();
    PySys_SetArgv(argc, argv);

    try 
    {
      object main_module = import("__main__");
      object main_namespace = main_module.attr("__dict__");
      object ignored = exec_file(argv[1], main_namespace, main_namespace);
    }
    catch(error_already_set const &)
    {
      PyErr_Print();
    }
  }
  return 0;
}
