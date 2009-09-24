#include "SkinHandler.hpp"
#include "particles/Computer.hpp"
#include "particles/Set.hpp"

using namespace espresso;
using namespace storage;

namespace {
  struct SkinComputer: public particles::Computer {
    PosPropertyHandle position, lastPosition;

    void prepare(PosPropertyHandle _position,
		 PosPropertyHandle _lastPosition) {
      position     = _position;
      lastPosition =_lastPosition;
    }

    void prepare(shared_ptr<Storage> store) {}
  };

  struct SkinCheck: public SkinComputer {
    real maxSqr;

    bool apply(ParticleHandle p) {
      real dist = (position[p] - lastPosition[p]).sqr();
      if (dist > maxSqr) maxSqr = dist;
      return true;
    }
  };

  struct SkinReset: public SkinComputer {
    bool apply(ParticleHandle p) {
      lastPosition[p] = position[p];
      return true;
    }
  };
}

void SkinHandler::checkSkin(particles::Set &set)
{
  Storage *store = set.getStorage().get();
  SkinCheck checker;
  checker.prepare(store->getPositionPropertyHandle(), getLastPositionPropertyHandle()); 
  checker.maxSqr = 0;
  set.foreach(checker);

  if (checker.maxSqr > skin*skin) {
    SkinReset reset;
    reset.prepare(store->getPositionPropertyHandle(), getLastPositionPropertyHandle());     
    skinExceeded(*this);
  }
}

ExtraSkinHandler::ExtraSkinHandler(shared_ptr< particles::Set > _set)
  : set(_set)
{
  Storage::SelfPtr store = set->getStorage();
  storagePositionsChanged = store->positionsChanged.connect(boost::bind(&ExtraSkinHandler::checkSkin, this));
}

void ExtraSkinHandler::setLastPositionProperty(boost::shared_ptr< Property< Real3D > > prop)
{
  if (prop)
    lastPositionProperty = prop;
  else
    lastPositionProperty = make_shared< Property< Real3D > >(set->getStorage());
}

PropertyHandle< Real3D > ExtraSkinHandler::getLastPositionPropertyHandle()
{
  return lastPositionProperty->getHandle(set);
}

void ExtraSkinHandler::setSkin(real _skin)
{
  skin = _skin;
  skinChanged(*this);
  checkSkin();
}

void ExtraSkinHandler::checkSkin()
{ SkinHandler::checkSkin(*set); }

