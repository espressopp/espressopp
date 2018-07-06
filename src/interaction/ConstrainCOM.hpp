/*
  Copyright (C) 2017
  Max Planck Institute for Polymer Research
  
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
#ifndef _INTERACTION_CONSTRAINCOM_HPP
#define _INTERACTION_CONSTRAINCOM_HPP

#include <cmath>
#include <map>
#include <boost/signals2.hpp>

#include "mpi.hpp"
#include "Potential.hpp"
#include "FixedLocalTupleComListInteractionTemplate.hpp"
#include "esutil/Error.hpp"

using namespace std;

namespace espressopp {
    namespace interaction {
	/** This class provides methods to compute forces and energies of the
	 *  ConstrainCOM part.
	 *  The code is based on G. Zhang's work. Reference in literature
	 *  G. Zhang, L. A. Moreira, T. Stuehn, K. CH. Daoulas, K. Kremer,
	 *  Macro Lett., 2014, 3, 198
	 */
	
	class ConstrainCOM : public PotentialTemplate< ConstrainCOM > {
	private:
	    real k_com;
	    
	public:
	    static void registerPython();
	    
	    ConstrainCOM(real _k_com);
	    
	    ~ConstrainCOM();
	    
	    
/////////////////////////////////////////////////////////////////////////////////////////
	    //  setters and getters
	    void setK_com(real _k_com) {
		k_com = _k_com;
	    }
	    real getK_com() const { return k_com; }
/////////////////////////////////////////////////////////////////////////////////////////
	    
	    real _computeEnergy(Real3D diff) {
		
		real energy = 0.0;
		energy += diff.sqr();

		return k_com*energy;
	    }
	    
	    Real3D _computeForce(Real3D diff,
				 real Total_mass) {
		
		return 2.*k_com*diff/Total_mass;
		
	    }
	    
	    real _computeEnergySqrRaw(real distSqr) const {
		LOG4ESPP_INFO(theLogger, "There is no sense to call this function for constrain of COM");
		return 0.0;
	    }
	    
	    bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
		LOG4ESPP_INFO(theLogger, "There is no sense to call this function for constrain of COM");
		return false;
	    }
	};
    }
}

#endif
