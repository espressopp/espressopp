#include <python.hpp>
#include <boost/python/implicit.hpp>
using namespace boost::python;

class Real3D;

class Real3DPtr {
  double *dataptr;

public:
  friend class Real3D;

  Real3DPtr(double *v) : dataptr(v) {}

  Real3DPtr(Real3DPtr &v) {
    dataptr = v.dataptr;
  }

  Real3DPtr(Real3D &v);

  double &operator[](const int i) {
    return dataptr[i];
  };

  double &at(const int i) {
    if (i < 0 || i > 2)
      throw std::out_of_range("Real3D::at");
    return (*this)[i];
  }
};


class Real3D {
  double data[3];
public:
  friend class Real3DPtr;

  Real3D() {
    for (int i = 0; i < 3; i++)
      data[i] = 0.0;
  }

  Real3D(Real3DPtr &v) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }

  Real3D(const double x, const double y, const double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }

  double &operator[](const int i) {
    return data[i];
  };

  const double &operator[](const int i) const {
    return data[i];
  }

  double &at(const int i) {
    if (i < 0 || i > 2)
      throw std::out_of_range("Real3D::at");
    return (*this)[i];
  }

  const double &at(const int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Real3D::at");
    return (*this)[i];
  }

  void setItem(const int i, double v) 
  { this->at(i) = v; }

  double getItem(const int i) const { return this->at(i); }

  
};

struct real3D_pickle_suite : pickle_suite
{
  static
  espresso::python::tuple
  getinitargs(const Real3D v)
  { return make_tuple(v[0], v[1], v[2]); }
};

Real3DPtr::Real3DPtr(Real3D &v)
  : dataptr(v.data) {}

void fold(Real3DPtr v) {
  v[0] = 42.0;
  v[1] = 52.0;
  v[2] = 62.0;
}

BOOST_PYTHON_MODULE(_real3D)
{
  class_<Real3D>("Real3D", init<>())
    .def(init< double, double, double >())
    .def("__getitem__", &Real3D::getItem)
    .def("__setitem__", &Real3D::setItem)
    .def_pickle(real3D_pickle_suite())
    ;

  implicitly_convertible<Real3D, Real3DPtr>();
  implicitly_convertible<Real3DPtr, Real3D>();

  def("fold", &fold);
}
