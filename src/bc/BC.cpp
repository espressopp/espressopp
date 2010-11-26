#include <python.hpp>
#include <boost/python/tuple.hpp>

#include "BC.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "esutil/RNG.hpp"

namespace espresso {
  namespace bc {

    LOG4ESPP_LOGGER(BC::logger, "BC");

    Real3D
    BC::getMinimumImageVector(const Real3D& pos1,
			      const Real3D& pos2) const {
      Real3D res;
      getMinimumImageVector(res, pos1, pos2);
      return res;
    }

    Real3D
    BC::getRandomPos() const {
      Real3D res;
      getRandomPos(res);
      return res;
    }

    void 
    BC::foldPosition(Real3D& pos, Int3D& imageBox) const {
      for (int i = 0; i < 3; ++i)
	foldCoordinate(pos, imageBox, i);
    }

    void
    BC::foldPosition(Real3D& pos) const {
      Int3D imageBox(0, 0, 0);
      foldPosition(pos, imageBox);
    }

    void 
    BC::unfoldPosition(Real3D& pos, Int3D& imageBox) const {
      for (int i = 0; i < 3; ++i)
	unfoldCoordinate(pos, imageBox, i);
    }

    boost::python::tuple 
    BC::getFoldedPosition(const Real3D& pos, 
			  const Int3D& imageBox) const {
      Real3D foldedPos = pos;
      Int3D foldedImageBox = imageBox;
      foldPosition(foldedPos, foldedImageBox);
      return boost::python::make_tuple(foldedPos, foldedImageBox);
    }

    boost::python::tuple 
    BC::getFoldedPosition(const Real3D& pos) const {
      Real3D foldedPos = pos;
      Int3D foldedImageBox(0, 0, 0);
      foldPosition(foldedPos, foldedImageBox);
      return boost::python::make_tuple(foldedPos, foldedImageBox);
    }

    Real3D
    BC::getUnfoldedPosition(const Real3D& pos, 
			    const Int3D& imageBox) const {
      Real3D unfoldedPos = pos;
      Int3D unfoldedImageBox = imageBox;
      unfoldPosition(unfoldedPos, unfoldedImageBox);
      return unfoldedPos;
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    BC::registerPython() {
      using namespace espresso::python;
    
      // also register the abstract class BC to make virtual functions available
      // be careful: boost::noncopyable must be used for abstract classes with pure routines
      // no_init must be used as the abstract class BC has no constructor

      Real3D (BC::*pygetMinimumImageVector)(const Real3D& pos1,
					    const Real3D& pos2) const 
	= &BC::getMinimumImageVector;
      
      boost::python::tuple (BC::*pygetFoldedPosition1)
	(const Real3D& pos, const Int3D& imageBox) const
	= &BC::getFoldedPosition;

      boost::python::tuple (BC::*pygetFoldedPosition2)
	(const Real3D& pos) const
	= &BC::getFoldedPosition;

      Real3D (BC::*pygetUnfoldedPosition)
	(const Real3D& pos, const Int3D& imageBox) const 
	= &BC::getUnfoldedPosition;

      Real3D (BC::*pygetRandomPos)() const = &BC::getRandomPos;
    
      class_< BC, boost::noncopyable >("bc_BC", no_init)
	.add_property("boxL", &BC::getBoxL)
	.add_property("rng", &BC::getRng, &BC::setRng)
	.def("getMinimumImageVector", pygetMinimumImageVector)
	.def("getFoldedPosition", pygetFoldedPosition1)
	.def("getFoldedPosition", pygetFoldedPosition2)
	.def("getUnfoldedPosition", pygetUnfoldedPosition)
	.def("getRandomPos", pygetRandomPos)
	;
    }
  }
}
