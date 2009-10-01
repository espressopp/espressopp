#include "DomainDecomposition.hpp"
#include "python.hpp"

using namespace espresso;
using namespace storage;

/// special cell for storing new particles that do not have a position yet
const size_t INSERT_BUFFER = 0;
/// where the the real domain cells start
const size_t CELLS_START   = 1;

DomainDecomposition::DomainDecomposition(bc::BC::SelfPtr _bc,
                                         const Real3D &cornerA, const Real3D &cornerB,
                                         size_t gridX, size_t gridY, size_t gridZ,
                                         real _skin)
  : Storage(_bc), SkinHandler(skin), grid(Real3DBox(cornerA, cornerB))
{
  // initially, we just have the insert buffer
  particles.resize(CELLS_START + grid.getMaxLinearId());
}

DomainDecomposition::~DomainDecomposition() {}

void DomainDecomposition::setDomain(const Real3D &cA, const Real3D &cB)
{
  grid.setDomain(Real3DBox(cA, cB));
}

void DomainDecomposition::setGrid(size_t gridX, size_t gridY, size_t gridZ)
{
  grid.setNumCells(0, gridX);
  grid.setNumCells(1, gridY);
  grid.setNumCells(2, gridZ);
  particles.resize(CELLS_START + grid.getMaxLinearId());
}

esutil::TupleVector &
DomainDecomposition::getTupleVector()
{
  return particles.raw();
}

void DomainDecomposition::prepare() {
  if (particles[INSERT_BUFFER].size() != 0) {
    resortParticles();
  }
  else if (positionsModified) {
    positionPropertyModified();
  }
  positionsModified = false;
}

void DomainDecomposition::positionPropertyModified() {
  if (checkSkin(*this)) {
    resortParticles();
  }
}

ParticleHandle DomainDecomposition::addParticle(ParticleId id) {
  BlockVector::Block b = particles[INSERT_BUFFER];
  b.resize(b.size() + 1);
  updateAndEmitHandlesChanged();
  return b.end() - 1;
}

void DomainDecomposition::deleteParticle(ParticleId deleteID) {
  ParticleHandle pos = getParticleHandle(deleteID);
  if (pos == ParticleHandle()) {
    throw std::out_of_range("DomainDecomposition::deleteParticle: particle does not exist");
  }
  // find cell of the particle to delete
  BlockVector::Block::thin_iterator it = pos;
  size_t blockId;
  for(blockId = INSERT_BUFFER; blockId < particles.size(); ++blockId) {
    BlockVector::Block block = particles[blockId];
    if (block.end() - it < 0) {
      --blockId;
      break;
    }
  }
  if (blockId < INSERT_BUFFER || blockId >= particles.size()) {
    throw std::runtime_error("DomainDecomposition::deleteParticle: FATAL: particle to delete exists, but is not in any cell");
  }
  
  updateAndEmitHandlesChanged();
}

ParticleHandle DomainDecomposition::getParticleHandle(ParticleId id) {
  IdMap::iterator it = location.find(id);
  if (it == location.end()) {
    return ParticleHandle();
  } else {
    return it->second;
  }
}

void DomainDecomposition::resortParticles() {
  PosPropertyHandle pos = getPositionPropertyHandle();

  // we cannot use a simple foreach loops here, because
  // a) we need the block number
  // b) we probably remove elements
  
  for(size_t blockId = INSERT_BUFFER; blockId < particles.size(); ++blockId) {
    BlockVector::Block block = particles[blockId];

    for (BlockVector::Block::thin_iterator it = block.begin(); it != block.end(); ++it) {
      Real3DGrid::CellIdentifier curCell = grid.locate(pos[*it]);
      if (curCell == Real3DGrid::notACell) {
        throw std::runtime_error("domain decomposition cannot handle particles outside the domain so far");
      }
      size_t curBlockId = CELLS_START + grid.linearize(curCell);
      if (curBlockId != blockId) {
        // buffer position, since it can be gone by inserting into a neighbor buffer
        size_t pos = it - block.begin();
        particles[curBlockId].push_back(*it);
        it = block.erase(block.begin() + pos);
        // update position handle, can be invalid due to resizing
        pos = getPositionPropertyHandle();
      }
    }
  }
  updateAndEmitHandlesChanged();
}

void DomainDecomposition::updateAndEmitHandlesChanged() {
  updateLocations();
  handlesChanged();
}

namespace {
  struct LocationUpdate: public particles::Computer {
    DomainDecomposition::IdMap &location;
    IdPropertyHandle id;

    LocationUpdate(DomainDecomposition::IdMap &_location,
                   IdPropertyHandle _id)
      : location(_location), id(_id) {}
    virtual void prepare(Storage::SelfPtr) {}
    virtual bool apply(ParticleHandle part) {
      location[id[part]] = part;
      return true;
    }
  };
}

void DomainDecomposition::updateLocations() {
  location.clear();
  LocationUpdate lu(location, getIdPropertyHandle());
  foreachApply(lu);
}

bool
DomainDecomposition::foreachApply(particles::Computer &computer) {
  for(size_t blockId = INSERT_BUFFER; blockId < particles.size(); ++blockId) {
    BlockVector::Block b = particles[blockId];
    BOOST_FOREACH(ParticleHandle particle, b) {
      if (!computer.apply(particle)) 
        return false;
    }
  }
  return true;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void DomainDecomposition::registerPython() {
  using namespace espresso::python;

  class_< DomainDecomposition, boost::noncopyable, bases < Storage > >
    ("storage_DomainDecomposition", init< bc::BC::SelfPtr,
     const espresso::Real3D&, const espresso::Real3D&,
     size_t, size_t, size_t, real >())
    .def("getSkin", &DomainDecomposition::getSkin)
    ;
}

