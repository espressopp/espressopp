#include <boost/python.hpp>
#include <sstream>
#include <iostream>
#include <boost/optional.hpp>
#include "python_optional.hpp"

using namespace std;
using boost::optional;

class Particle
{
  string id;
  double x, y, z;
public:
  Particle() { 
    id = "EMPTY";
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }
  string getId() const { return id; }
  double getX() const {return x;}
  double getY() const {return y;}
  double getZ() const {return z;}
  void set(string _id, 
	   double _x, 
	   double _y, 
	   double _z) {
    id = _id;
    x = _x;
    y = _y;
    z = _z;
  }
  const string toStr() const {
    stringstream ss;
    ss << "id=\'" << id 
       << "\', x=" << x 
       << ", y=" << y 
       << ", z=" << z;
    return ss.str();
  }


  // testing boost::optional feature
  void set_optional(optional<string> _id, 
		    optional<double> _x, 
		    optional<double> _y, 
		    optional<double> _z) { 
    if (_id) id = *_id;
    if (_x) x = *_x;
    if (_y) y = *_y;
    if (_z) z = *_z;
  }

  static void register_class() {
    using namespace boost::python;

    // testing boost::optional
    python_optional<double>();
    python_optional<string>();

    class_<Particle>("_Particle", init<>())
      .def("__str__", &Particle::toStr)
      .def("set", &Particle::set)
      .add_property("id",&Particle::getId)
      .add_property("x", &Particle::getX)
      .add_property("y", &Particle::getY)
      .add_property("z", &Particle::getZ)
      // testing boost::optional
      .def("set_optional", &Particle::set_optional);
  }
};

BOOST_PYTHON_MODULE(world)
{
  Particle::register_class();
};
