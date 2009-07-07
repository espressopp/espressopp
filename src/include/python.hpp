#ifndef _PYTHON_HPP
#define _PYTHON_HPP

#include "acconfig.hpp"
#include <boost/python.hpp>
#include "types.hpp"

namespace espresso {
  namespace python {
    using namespace boost::python;

    template <    
      class W // class being wrapped
      , class X1 = ::boost::python::detail::not_specified
      , class X2 = ::boost::python::detail::not_specified
      >
    class class_
      : public boost::python::class_< W, shared_ptr< W >, X1, X2 > 
    {
      typedef typename boost::python::class_< W, shared_ptr< W >, X1, X2 > base;
      
    public:
      // Construct with the class name, with or without docstring, and default __init__() function
      class_(char const* name, char const* doc = 0) 
	: base(name, doc) {}
      
      // Construct with class name, no docstring, and an uncallable __init__ function
      class_(char const* name, boost::python::no_init_t) 
	: base(name, boost::python::no_init) {}
      
      // Construct with class name, docstring, and an uncallable __init__ function
      class_(char const* name, char const* doc, boost::python::no_init_t)
	: base(name, doc, boost::python::no_init) {}
      
      // Construct with class name and init<> function
      template <class DerivedT>
      inline class_(char const* name, boost::python::init_base<DerivedT> const& i)
	: base(name, i) {}
      
      // Construct with class name, docstring and init<> function
      template <class DerivedT>
      inline class_(char const* name, char const* doc, boost::python::init_base<DerivedT> const& i)
	: base(name, doc, i) {}
    };
  }
}

#endif
