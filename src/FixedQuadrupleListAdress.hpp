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
#ifndef _FIXEDQUADRUPLELISTADRESS_HPP
#define _FIXEDQUADRUPLELISTADRESS_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Triple.hpp"

#include "Particle.hpp"
#include "FixedTupleListAdress.hpp"
#include "FixedQuadrupleList.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
//#include "FixedListComm.hpp"

namespace espressopp {
class FixedQuadrupleListAdress : public FixedQuadrupleList {
 public:
  FixedQuadrupleListAdress(shared_ptr< storage::Storage > _storage,
                           shared_ptr<FixedTupleListAdress> _fixedtupleList);
  ~FixedQuadrupleListAdress();

  /** Add the given particle quadruple to the list on this processor if the
	  particle with the lower id belongs to this processor.  Note that
	  this routine does not check whether the quadruple is inserted on
	  another processor as well.
	
	  \return whether the quadruple was inserted on this processor.
  */
  bool add(longint pid1, longint pid2, longint pid3, longint pid4);
  void beforeSendParticles(ParticleList& pl, class OutBuffer &buf);
  void beforeSendATParticles(std::vector<longint>& atpl, class OutBuffer& buf);
  void onParticlesChanged();

  static void registerPython();

 protected:
  // Fixedtuple list connects to this and triggers beforeSendATParticles()
  boost::signals2::connection sigBeforeSendAT, sigAfterRecvAT;

 private:
  shared_ptr<FixedTupleListAdress> fixedtupleList;
  using QuadrupleList::add;
  static LOG4ESPP_DECL_LOGGER(theLogger);
};
}
#endif
