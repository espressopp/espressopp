/*
  Copyright (C) 2014
      Pierre de Buyl
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

#include "python.hpp"
#include "SlabBC.hpp"
#include <cmath>
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {
  namespace bc {
    /* Constructor */
    SlabBC::
    SlabBC(shared_ptr< esutil::RNG > _rng, const Real3D& _boxL) : BC(_rng) {
      setBoxL(_boxL);
      slabDir=0;
    }

    /* Setter method for the box length */
    void SlabBC::setBoxL(const Real3D& _boxL) {
      boxL = _boxL;
      for (int i = 0; i < 3; i++) {
	    invBoxL[i] = 1.0/boxL[i];
	    boxL2[i] = 0.5*boxL[i];
      }
      onBoxDimensionsChanged();
    }
    void SlabBC::scaleVolume(real s) {
  	  boxL *= s;
  	  boxL2 *= s;
  	  invBoxL /= s;
  	  onBoxDimensionsChanged();
    }
    void SlabBC::scaleVolume(Real3D s) {
  	  boxL[0] *= s[0];
  	  boxL[1] *= s[1];
  	  boxL[2] *= s[2];
  	  boxL2[0] *= s[0];
  	  boxL2[1] *= s[1];
  	  boxL2[2] *= s[2];
  	  invBoxL[0] /= s[0];
  	  invBoxL[1] /= s[1];
  	  invBoxL[2] /= s[2];
  	  onBoxDimensionsChanged();
    }

    /* Returns the minimum image vector between two positions */
    void
    SlabBC::
    getMinimumImageVector(Real3D& dist,
			  const Real3D& pos1,
			  const Real3D& pos2) const {
      dist = pos1;
      dist -= pos2;

      for (int i=0; i<3; i++) {
	if (i!=slabDir) {
	  dist[i] -= round(dist[i] * invBoxL[i]) * boxL[i];
	}
      }
    }

    /* Returns the minimum image vector between two positions */
    void
    SlabBC::
    getMinimumImageVectorBox(Real3D& dist,
                             const Real3D& pos1,
                             const Real3D& pos2) const {
      dist = pos1;
      dist -= pos2;

      for (int i=0; i<3; i++) {
	if (i!=slabDir) {
	  if (dist[i] < -boxL2[i]) dist[i] += boxL[i];
	  else if (dist[i] > boxL2[i]) dist[i] -= boxL[i];
	}
      }
    }

    /* Fold back a nearby position into box */

    void
    SlabBC::getMinimumDistance(Real3D& dist) const {

      for (int i=0; i<3; i++) {
	if (i!=slabDir) {
	  if (dist[i] < -0.5 * boxL[i]) dist[i] += boxL[i];
	  else if (dist[i] > 0.5 * boxL[i]) dist[i] -= boxL[i];
	}
      }
    }

    /* Returns the minimum image vector between two positions */
    void
    SlabBC::
    getMinimumImageVectorX(real dist[3],
                          const real pos1[3],
                          const real pos2[3]) const {

      for (int i=0; i<3; i++) {
	if (i!=slabDir) {
	  dist[i] = pos1[i];
	  dist[i] -= pos2[i];
	  dist[i] -= round(dist[i] * invBoxL[i]) * boxL[i];
	}
      }
    }

    /* Fold an individual coordinate in the specified direction */
    void
    SlabBC::
    foldCoordinate(Real3D& pos, Int3D& imageBox, int dir) const {
      int tmp = static_cast<int>(floor(pos[dir]*invBoxL[dir]));

      imageBox[dir] += tmp;
      pos[dir] -= tmp*boxL[dir];

      if(pos[dir] < 0 || pos[dir] >= boxL[dir]) {
        /* slow but safe */
        if (fabs(pos[dir]*invBoxL[dir]) >= INT_MAX/2) {
# warning ERRORHANDLING MISSING
#if 0
// errortext: particle coordinate out of range
#endif
          imageBox[dir] = 0;
          pos[dir] = 0;
        }
      }

    }

    /* Unfold an individual coordinate in the specified direction */
    void
    SlabBC::
    unfoldCoordinate(Real3D& pos, Int3D& imageBox, int dir) const {
      pos[dir] += imageBox[dir]*boxL[dir];
      imageBox[dir] = 0;
    }

    /* Get random position in the central image box */
    void
    SlabBC::
    getRandomPos(Real3D& res) const {
      for(int k = 0; k < 3; k++)
        res[k] = boxL[k];

      res[0] *= (*rng)();
      res[1] *= (*rng)();
      res[2] *= (*rng)();
    }

    void
    SlabBC::
    registerPython() {
      using namespace espressopp::python;
      class_<SlabBC, bases< BC >, boost::noncopyable >
	("bc_SlabBC", init< shared_ptr< esutil::RNG >, Real3D& >())
	.add_property("boxL", &SlabBC::getBoxL, &SlabBC::setBoxL)
      ;
    }
  }
}
