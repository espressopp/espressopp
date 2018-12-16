/*
  Copyright (C) 2018
      Jakub Krajniak (jkrajniak at gmail.com)
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

// ESPP_CLASS
#ifndef _ANALYSIS_OBSERVABLE_HPP
#define _ANALYSIS_OBSERVABLE_HPP
#include "python.hpp"
#include "types.hpp"
#include "SystemAccess.hpp"
#include <vector>

namespace espressopp {
  namespace analysis {
    /** All quantities to be measured derive from this abstract base class. */
    class Observable : public SystemAccess {
    public:
      /** Observable can be int, real, scalar or vector */
      enum result_types {old_format=-1, none=0, real_scalar=1, int_scalar=2, real_vector=3, int_vector=4};
      enum ObservableTypes {
        POTENTIAL_ENERGY,
        KINETIC_ENERGY,
        OTHER
      };
      Observable(shared_ptr< System > system) : SystemAccess(system) {
        result_type=old_format;
        observable_type = OTHER;
      };
      virtual ~Observable() {};

    public:
      /** for compatibilty with old compute function only, will be deleted soon */
      virtual real compute() const { return 0.0; };
      /** returns observable of type real, used for Python and on C++ level*/
      virtual real compute_real() const { return 0.0;};
      /** returns observable of type int, used for Python and on C++ level*/
      virtual int compute_int() const { return 0; };
      /** computes vector of real values (e.g. pressure tensor, ...), used on C++ level only */
      virtual void compute_real_vector(){ return; };
      /** computes vector of integer values, used on C++ level only */
      virtual void compute_int_vector(){ return; };

      /** returns python list of real values (e.g. pressure tensor, ...), used on Python level*/
      virtual python::list compute_real_vector_python();
      /** returns python list of integer values, used on Python level*/
      virtual python::list compute_int_vector_python();

      /** returns the result type of the observable */
      // TODO at the moment it returns int instead, because it was causing an error
      // trying to convert result_types to python
      //result_types getResultType() { return result_type; };
      int getResultType() { return result_type; };

      static void registerPython();

     protected:
      result_types result_type;
      ObservableTypes observable_type;
      std::vector< real > result_real_vector;
      std::vector< int > result_int_vector;
      longint result_vector_size;

      static LOG4ESPP_DECL_LOGGER(logger);

    };
  }
}

#endif
