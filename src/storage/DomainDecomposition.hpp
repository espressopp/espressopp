#ifndef DOMAIN_DECOMPOSITION_HPP
#define DOMAIN_DECOMPOSITION_HPP

#include <boost/unordered_map.hpp>
#include "Storage.hpp"
#include "SkinHandler.hpp"
#include "esutil/BlockVector.hpp"

namespace espresso {
  namespace storage {
    class DomainDecomposition: public Storage, SkinHandler {
    public:
      DomainDecomposition(bc::BC::SelfPtr,
                          const Real3D &, const Real3D &,
                          size_t gridX, size_t gridY, size_t gridZ,
                          real);
      virtual ~DomainDecomposition();

      virtual ParticleHandle getParticleHandle(ParticleId id);
      virtual ParticleHandle addParticle(ParticleId);
      virtual void deleteParticle(ParticleId);

      /// set the property represeting the last particle positions; the corresponding Storage function
      virtual void setLastPositionProperty(boost::shared_ptr< Property< Real3D > > =
					   boost::shared_ptr< Property< Real3D > >());

      /// get the property represeting the last particle positions; the corresponding Storage function
      virtual PropertyHandle< Real3D > getLastPositionPropertyHandle();

      /// set the grid
      virtual void setGrid(size_t gridX, size_t gridY, size_t gridZ);

      /// set the corners of the domain
      virtual void setDomain(const Real3D &, const Real3D &);

      /// set the skin
      virtual void setSkin(real);

      /// get the current skin; overrides the corresponding Storage function
      virtual real getSkin();

      /** prepare domain decomposition; this has to be called if it is unsure whether the skin criterion
          still holds or there might be particles in the insert buffer, i.e. whenever coming from python.
      */
      virtual void prepare();

      /// checks the skin criterion, and if necessary, resorts
      virtual void positionPropertyModified();

      /// force sorting of the particles into the correct cells
      virtual void resortParticles();

      /// make this class available at Python
      static void registerPython();

      typedef boost::unordered_map<ParticleId, ParticleHandle> IdMap;
      typedef esutil::BlockVector<esutil::TupleVector> BlockVector;

    protected:
      virtual bool foreachApply(particles::Computer &);
      virtual esutil::TupleVector &getTupleVector();

      /// send the signal that handles might have changed and update our own
      virtual void updateAndEmitHandlesChanged();

      /// force sorting of the particles into the correct cells
      void updateLocations();

    private:
      /// cache with the particle handles for each particle
      IdMap location;

      /// here the particle data is stored
      BlockVector particles;

      Real3DGrid grid;
    };
  }
}

#endif
