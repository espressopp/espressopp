/*
  Copyright (C) 2012-2015
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

#ifndef _PARTICLEACCESS_HPP
#define	_PARTICLEACCESS_HPP
#include "python.hpp"
#include "SystemAccess.hpp"

namespace espressopp {
  class ParticleAccess : public SystemAccess {
  public:
    ParticleAccess(shared_ptr< System > system) : SystemAccess(system) {}
    virtual ~ParticleAccess() {}

    virtual void perform_action() = 0;
    
    static void registerPython();
  };

}
#endif	/* PARTICLEACCESS_HPP */
