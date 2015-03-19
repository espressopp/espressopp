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
#ifndef _STORAGE_DOMAINDECOMPOSITIONNONBLOCKING_HPP
#define _STORAGE_DOMAINDECOMPOSITIONNONBLOCKING_HPP
#include "DomainDecomposition.hpp"
// #include "types.hpp"

namespace espressopp {
  namespace storage {
    class DomainDecompositionNonBlocking: public DomainDecomposition {
    public:
      DomainDecompositionNonBlocking(shared_ptr< System > system,
              const Int3D& _nodeGrid,
			  const Int3D& _cellGrid);
      virtual ~DomainDecompositionNonBlocking() {}
      static void registerPython();
    protected:
      virtual void decomposeRealParticles();
      virtual void doGhostCommunication(bool sizesFirst, bool realToGhosts, const int dataElements = 0);
      mpi::request isendParticles(OutBuffer &data, ParticleList &list, longint node);
      mpi::request irecvParticles_initiate(InBuffer &data, longint node);
      void irecvParticles_finish(InBuffer &data, ParticleList &list);
    private:
      InBuffer inBufferL;
      InBuffer inBufferR;
      OutBuffer outBufferL;
      OutBuffer outBufferR;
      InBuffer inBufferG;   // used for ghost communication
      OutBuffer outBufferG;  // used for ghost communication
    };
  }
}
#endif
