#ifndef STORAGE_SKIN_HANDLER_HPP
#define STORAGE_SKIN_HANDLER_HPP

#include "esutil/MultiSignalConnections.hpp"
#include "Property.hpp"
#include "storage/PropertyHandle.hpp"

namespace espresso {
  namespace particles {
    class Set;
  }
  namespace storage {
    /** abstract class checking a skin condition. */
    class SkinHandler {
      typedef boost::signals2::signal1<void, const SkinHandler &> Signal;
    public:
      virtual ~SkinHandler() {}

      /** call this to specify the property representing the particle
	  position at the last skin event.  This has to be done in
	  Python so that we get a proper wrapper.

	  If no property is specified, this constructs the property
	  internally.  Use this in simply C++-code, where access to
	  the handles is sufficient.
      */
      virtual void setLastPositionProperty(boost::shared_ptr< Property< Real3D > > =
					   boost::shared_ptr< Property< Real3D > >()) = 0;
      virtual PropertyHandle< Real3D > getLastPositionPropertyHandle() = 0;

      virtual void setSkin(real) = 0;
      real getSkin() const { return skin; }
      real getSkinSqr() const { return skinSqr; }

      /** Signal emitted if the skin was changed. Note that the skin
	  criterion is checked anyways and if necessary, skinExceed
	  emitted.
       */
      Signal skinChanged;
      /** Signal emitted if at least one particle moved further
	  than the skin.
      */
      Signal skinExceeded;
      /// connection manager
      esutil::MultiSignalConnections connections;

    protected:
      // function checking the skin criterion, and emitting skinExceeded if applicable
      void checkSkin(particles::Set &);

      real skin, skinSqr;
    };

    /// SkinHandler living outside a storage
    class ExtraSkinHandler: public SkinHandler {
      ExtraSkinHandler(shared_ptr< particles::Set >);
      virtual void setLastPositionProperty(boost::shared_ptr< Property< Real3D > > =
					   boost::shared_ptr< Property< Real3D > >());
      virtual PropertyHandle< Real3D > getLastPositionPropertyHandle();
      virtual void setSkin(real);
    private:
      shared_ptr< particles::Set > set;
      shared_ptr< Property< Real3D > > lastPositionProperty;
    };
  }
}
#endif
