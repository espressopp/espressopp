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

// ESPP_CLASS
#ifndef _BC_ORTHORHOMBICBC_HPP
#define _BC_ORTHORHOMBICBC_HPP

#include "BC.hpp"
#include "Real3D.hpp"

namespace espressopp {
  namespace bc {
    class OrthorhombicBC : public BC {
    private:
      Real3D boxL;
      Real3D boxL2;  // half box length
      Real3D invBoxL;

    public:
      /** Virtual destructor for boundary conditions. */
      virtual
      ~OrthorhombicBC() {}

      /** Constructor */
      OrthorhombicBC(shared_ptr< esutil::RNG > _rng, 
		     const Real3D& _boxL);

      /** Method to set the length of the side of the cubic simulation cell */
      virtual void
      setBoxL(const Real3D& _boxL);

      /** Getters for box dimensions */
      virtual Real3D getBoxL() const { return boxL; }

      /** Scale the Volume of the box by s^(1/3) ??? (Box-Length is scaled by s) */
      virtual void scaleVolume(real s);
      /** Scale the Volume of the box anisotropic case*/
      virtual void scaleVolume(Real3D s);

      /** Computes the minimum image distance vector between two
          positions. This routine must be implemented by derived
          classes (once the code stabilizes).

          \param dist is the distance vector (pos2 - pos1)
          \param pos1, pos2 are the particle positions 
      */
      virtual void
      getMinimumImageVector(Real3D& dist,
                            const Real3D& pos1,
                            const Real3D& pos2) const;

      virtual void
      getMinimumImageVectorBox(Real3D& dist,
                               const Real3D& pos1,
                               const Real3D& pos2) const;

      virtual void
      getMinimumImageVectorX(real dist[3],
                            const real pos1[3],
                            const real pos2[3]) const;

      /** Compute the minimum image distance where the distance
          is given by two positions in the box.
      */

      virtual void
      getMinimumDistance(Real3D& dist) const;

      /** fold a coordinate to the primary simulation box.
	  \param pos         the position...
	  \param imageBox    and the box
	  \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      virtual void 
      foldCoordinate(Real3D& pos, Int3D& imageBox, int dir) const;

      /** unfold a coordinate to physical position.
	  \param pos the position...
	  \param imageBox and the box
	
	  Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
	  afterwards.
      */
      virtual void 
      unfoldCoordinate(Real3D& pos, Int3D& imageBox, int dir) const;

      /** Get a random position within the central simulation box. The
          positions are assigned with each coordinate on [0, boxL]. */
      virtual void
      getRandomPos(Real3D& res) const;

      static void registerPython();
    };
  }
}

#endif
