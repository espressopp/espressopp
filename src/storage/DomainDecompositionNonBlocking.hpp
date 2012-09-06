// ESPP_CLASS
#ifndef _STORAGE_DOMAINDECOMPOSITIONNONBLOCKING_HPP
#define _STORAGE_DOMAINDECOMPOSITIONNONBLOCKING_HPP
#include "DomainDecomposition.hpp"
// #include "types.hpp"

namespace espresso {
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
