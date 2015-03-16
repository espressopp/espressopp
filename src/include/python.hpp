/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef _PYTHON_HPP
#define _PYTHON_HPP

#include "acconfig.hpp"
#include <boost/python.hpp>
#include "types.hpp"

namespace espressopp {
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
      class_(char const* name,  boost::python::no_init_t) 
	: base(name, boost::python::no_init) {}
      
      // Construct with class name, docstring, and an uncallable __init__ function
      class_(char const* name, char const* doc,  boost::python::no_init_t)
	: base(name, doc, no_init) {}
      
      // Construct with class name and init<> function
      template <class DerivedT>
      inline class_(char const* name, boost::python::init_base< DerivedT > const& i)
	: base(name, i) {}
      
      // Construct with class name, docstring and init<> function
      template <class DerivedT>
      inline class_(char const* name, char const* doc, boost::python::init_base<DerivedT> const& i)
	: base(name, doc, i) {}
    };
  }
}

#endif
