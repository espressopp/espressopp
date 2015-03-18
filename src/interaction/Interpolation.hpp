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

// ESPPP_CLASS
#ifndef _INTERACTION_INTERPOLATION_HPP
#define _INTERACTION_INTERPOLATION_HPP

#include <stdio.h>
#include "types.hpp"
#include "logging.hpp"
#include "mpi.hpp"

namespace espressopp {
    namespace interaction {
        
        class Interpolation {
            public:
                virtual real getEnergy(real r) const = 0;
                virtual real getForce(real r) const = 0;
                virtual void read(mpi::communicator comm, const char* file) = 0;
        };//class Interpolation
        
        
        template <class Derived>
        class InterpolationTemplate: public Interpolation {
            public:
                InterpolationTemplate();
                virtual real getEnergy(real r) const;
                virtual real getForce(real r) const;
                virtual void read(mpi::communicator comm, const char* file);
            
            protected:
                Derived* derived_this(){
                    return static_cast <Derived*> (this);
                }
                const Derived* derived_this() const {
                    return static_cast <const Derived*> (this);
                }
                
        };//template class InterpolationTemplate
            
            
        /*********************************************/
        /* INLINE IMPLEMENTATION                     */
        /*********************************************/
        template <class Derived>
        inline
        InterpolationTemplate <Derived>::InterpolationTemplate(){}
        
        template <class Derived>
        inline real
        InterpolationTemplate <Derived>::
        getEnergy(real r) const {
            return derived_this()->getEnergyRaw(r);
        }
        
        template <class Derived>
        inline real
        InterpolationTemplate <Derived>::
        getForce(real r) const {
            return derived_this()->getForceRaw(r);
        }
        
        template <class Derived>
        inline void
        InterpolationTemplate <Derived>::
        read(mpi::communicator comm, const char* file) {
            derived_this()->readRaw(comm, file);
        }
        
        
    }//ns interaction
} //ns espressopp


#endif
