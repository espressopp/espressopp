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
#ifndef _INTERACTION_CONSTRAINRG_HPP
#define _INTERACTION_CONSTRAINRG_HPP

#include <cmath>
#include <map>
#include <boost/signals2.hpp>

#include "mpi.hpp"
#include "Potential.hpp"
#include "FixedLocalTupleRgListInteractionTemplate.hpp"
#include "esutil/Error.hpp"

using namespace std;

namespace espressopp {
    namespace interaction {
	/** This class provides methods to compute forces and energies of the
	 *  ConstrainRG part.
	 *  The code is based on G. Zhang's work. Reference in literature
	 *  G. Zhang, L. A. Moreira, T. Stuehn, K. CH. Daoulas, K. Kremer,
	 *  Macro Lett., 2014, 3, 198
	 */
	
	class ConstrainRG : public PotentialTemplate< ConstrainRG > {
	private:
	    real k_rg;
	    
	public:
	    static void registerPython();
	    
	    ConstrainRG(real _k_rg);
	    
	    ~ConstrainRG();
	    
	    
/////////////////////////////////////////////////////////////////////////////////////////
	    //  setters and getters
	    void setK_rg(real _k_rg) {
		k_rg = _k_rg;
	    }
	    real getK_rg() const { return k_rg; }
/////////////////////////////////////////////////////////////////////////////////////////

	    // This potential term is based on the reference:
            // Equilibration of high molecular weight polymer melts: A hierarchical strategy,
	    // Macro Lett., 2014, 3, 198
	    real _computeEnergy(real diff_rg) {
		
		real energy = diff_rg*diff_rg;

		return k_rg*energy;
	    }
	    
	    // This force term is based on the reference:
            // Equilibration of high molecular weight polymer melts: A hierarchical strategy,
	    // Macro Lett., 2014, 3, 198
	    Real3D _computeForce(Real3D diff,
				 real diff_rg,
				 long unsigned int N_Constrain) {
		
		return 4.*k_rg*diff*diff_rg/N_Constrain;
	    }
	    
	    real _computeEnergySqrRaw(real distSqr) const {
		LOG4ESPP_INFO(theLogger, "There is no sense to call this function for constrain of RG");
		return 0.0;
	    }
	    
	    bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
		LOG4ESPP_INFO(theLogger, "There is no sense to call this function for constrain of RG");
		return false;
	    }
	};
    }
}

#endif
