/*
  Copyright (C) 2018-2019
      The VOTCA Development Team
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

//! Macro to detect strictly gcc.
//! \details __GNUC__ and __GNUG__ were intended to indicate the GNU compilers.
//! However, they're also defined by Clang/LLVM and Intel compilers to indicate
//! compatibility. This macro can be used to detect strictly gcc and not clang
//! or icc.
#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#define STRICT_GNUC
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

#if (defined STRICT_GNUC) && GCC_VERSION > 50000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wplacement-new"
#endif
#include <boost/python.hpp>
#if (defined STRICT_GNUC) && GCC_VERSION > 50000
#pragma GCC diagnostic pop
#endif

#include "acconfig.hpp"
#include "types.hpp"

namespace espressopp
{
namespace python
{
using namespace boost::python;

template <class W  // class being wrapped
          ,
          class X1 = ::boost::python::detail::not_specified,
          class X2 = ::boost::python::detail::not_specified>
class class_ : public boost::python::class_<W, std::shared_ptr<W>, X1, X2>
{
    typedef typename boost::python::class_<W, std::shared_ptr<W>, X1, X2> base;

public:
    // Construct with the class name, with or without docstring, and default __init__() function
    class_(char const* name, char const* doc = 0) : base(name, doc) {}

    // Construct with class name, no docstring, and an uncallable __init__ function
    class_(char const* name, boost::python::no_init_t) : base(name, boost::python::no_init) {}

    // Construct with class name, docstring, and an uncallable __init__ function
    class_(char const* name, char const* doc, boost::python::no_init_t) : base(name, doc, no_init)
    {
    }

    // Construct with class name and init<> function
    template <class DerivedT>
    inline class_(char const* name, boost::python::init_base<DerivedT> const& i) : base(name, i)
    {
    }

    // Construct with class name, docstring and init<> function
    template <class DerivedT>
    inline class_(char const* name, char const* doc, boost::python::init_base<DerivedT> const& i)
        : base(name, doc, i)
    {
    }
};
}  // namespace python
}  // namespace espressopp

#endif
