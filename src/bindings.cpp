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

/* python.hpp has to be included before any system headers according
   to Python API to avoid redefining _POSIX_C_SOURCE by Python.h */
#include <python.hpp>
#include "bindings.hpp"

#include <Particle.hpp>
#include <ParticleGroup.hpp>
#include <System.hpp>
#include <VerletList.hpp>
#include <VerletListAdress.hpp>
#include <VerletListTriple.hpp>
#include <FixedSingleList.hpp>
#include <FixedPairList.hpp>
#include <FixedPairDistList.hpp>
#include <FixedPairListAdress.hpp>
#include <FixedTripleList.hpp>
#include <FixedTripleAngleList.hpp>
#include <FixedTripleListAdress.hpp>
#include <FixedQuadrupleList.hpp>
#include <FixedQuadrupleListAdress.hpp>
#include <FixedQuadrupleAngleList.hpp>
#include <FixedTupleList.hpp>
#include <FixedTupleListAdress.hpp>
#include <FixedLocalTupleList.hpp>
#include <Real3D.hpp>
#include <Quaternion.hpp>
#include <Tensor.hpp>
#include <Int3D.hpp>
#include <Version.hpp>
#include <ParticleAccess.hpp>
#include <RealND.hpp>

#include <esutil/PyLogger.hpp>
#include <esutil/bindings.hpp>
#include <bc/bindings.hpp>
#include <storage/bindings.hpp>
#include <integrator/bindings.hpp>
#include <interaction/bindings.hpp>
#include <analysis/bindings.hpp>
#include <io/bindings.hpp>

void espressopp::registerPython() {
  espressopp::Particle::registerPython();
  espressopp::ParticleGroup::registerPython();
  espressopp::System::registerPython();
  espressopp::VerletList::registerPython();
  espressopp::VerletListAdress::registerPython();
  espressopp::VerletListTriple::registerPython();
  espressopp::FixedSingleList::registerPython();
  espressopp::FixedPairList::registerPython();
  espressopp::FixedPairDistList::registerPython();
  espressopp::FixedPairListAdress::registerPython();
  espressopp::FixedTripleList::registerPython();
  espressopp::FixedTripleAngleList::registerPython();
  espressopp::FixedTripleListAdress::registerPython();
  espressopp::FixedQuadrupleList::registerPython();
  espressopp::FixedQuadrupleListAdress::registerPython();
  espressopp::FixedQuadrupleAngleList::registerPython();
  espressopp::FixedTupleList::registerPython();
  espressopp::FixedTupleListAdress::registerPython();
  espressopp::FixedLocalTupleList::registerPython();
  espressopp::Real3D::registerPython();
  espressopp::Quaternion::registerPython();
  espressopp::RealND::registerPython();
  espressopp::Tensor::registerPython();
  espressopp::Int3D::registerPython();
  espressopp::Version::registerPython();
  espressopp::ParticleAccess::registerPython();
  espressopp::RealNDs::registerPython();

  espressopp::esutil::registerPython();
  espressopp::bc::registerPython();
  espressopp::storage::registerPython();
  espressopp::integrator::registerPython();
  espressopp::interaction::registerPython();
  espressopp::analysis::registerPython();
  espressopp::io::registerPython();

  log4espp::PyLogger::registerPython();
}
