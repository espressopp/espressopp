/*
  Copyright (C) 2017
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
#ifndef _BC_BC_HPP
#define _BC_BC_HPP

#include <boost/python/tuple.hpp>
#include <boost/signals2.hpp>
#include "types.hpp"
#include "log4espp.hpp"
#include "SystemAccess.hpp"

namespace espressopp {
  namespace bc {
    /** Abstract base class for boundary conditions. */
    class BC {
    public:
      BC(shared_ptr< esutil::RNG > _rng) : rng(_rng) {}

      /** Virtual destructor for boundary conditions. */
      virtual ~BC() {}

      /** Getter for box dimensions */
      virtual Real3D getBoxL() const = 0;
      virtual void scaleVolume(real s) = 0;
      virtual void scaleVolume(Real3D s) = 0;
      /** Getter for the RNG. */
      virtual shared_ptr< esutil::RNG > getRng() { return rng; }
      /** Setter for RNG. */
      virtual void setRng(shared_ptr< esutil::RNG > _rng) { rng = _rng; }

      /** Computes the minimum image distance vector between two
	  positions.

	  \param dist is the distance vector (pos2 - pos1)
	  \param pos1, pos2 are the particle positions 
      */
      virtual void
      getMinimumImageVector(Real3D& dist,
			    const Real3D& pos1,
			    const Real3D& pos2) const = 0;

      /** Computes the minimum image distance vector between two
          positions where both positions are in the box.

          \param dist is the distance vector (pos2 - pos1)
          \param pos1, pos2 are the particle positions 

          Be careful: this routine is faster than getMinimumImageVector
          but will only deliver corrrect result if the absolute distance
          in each dimension is less than the box size.
      */
      virtual void
      getMinimumImageVectorBox(Real3D& dist,
                               const Real3D& pos1,
                               const Real3D& pos2) const = 0;

      virtual void
      getMinimumImageVectorX(real dist[3],
			    const real pos1[3],
			    const real pos2[3]) const = 0;
      virtual Real3D 
      getMinimumImageVector(const Real3D& pos1,
			    const Real3D& pos2) const;

      /** Compute the minimum image distance where the distance
          is given by two positions in the box.
      */
      virtual void
      getMinimumDistance(Real3D& dist) const = 0;

      /** fold a coordinate to the primary simulation box.
	  \param pos         the position
	  \param imageBox    and the box
	  \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      virtual void 
      foldCoordinate(Real3D& pos, Int3D& imageBox, int dir) const = 0;

      virtual void 
      unfoldCoordinate(Real3D& pos, Int3D& imageBox, int dir) const = 0;

      /** Fold the coordinates to the primary simulation box.
	  \param pos the position...
	  \param imageBox and the box

	  Both pos and image_box are I/O,
	  i. e. a previously folded position will be folded correctly.
      */
      virtual void 
      foldPosition(Real3D& pos, Int3D& imageBox) const;

      virtual void
      foldPosition(Real3D& pos) const;

      /** Get the folded position as a Python tuple. */
      virtual boost::python::tuple
      getFoldedPosition(const Real3D& pos, 
			const Int3D& imageBox) const;

      virtual boost::python::tuple 
      getFoldedPosition(const Real3D& pos) const;

      /** unfold coordinates to physical position.
	  \param pos the position...
	  \param imageBox and the box
	
	  Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
	  afterwards.
      */
      virtual void 
      unfoldPosition(Real3D& pos, Int3D& imageBox) const;

      virtual Real3D
      getUnfoldedPosition(const Real3D& pos, const Int3D& imageBox) const;

      /** Get a random position within the central simulation box. The
	  positions are assigned with each coordinate on [0, boxL]. */
      virtual void
      getRandomPos(Real3D& res) const = 0;

      virtual Real3D
      getRandomPos() const;

      // This signal is called whenever the box dimensions did change
      // (e.g. when scaling the Volume or when setting boxL)
      boost::signals2::signal<void ()> onBoxDimensionsChanged;

      static void registerPython();

     protected:
      shared_ptr< esutil::RNG > rng;

      static LOG4ESPP_DECL_LOGGER(logger);

    }; 
  }
}

#endif
