/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  Copyright (C) 2019
      Max Planck Computing and Data Facility

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

#include "python.hpp"

//#include <algorithm>

#include "log4espp.hpp"

#include "System.hpp"
#include "Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListIterator.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "Particle.hpp"
#include "Buffer.hpp"
#include "esutil/Error.hpp"

#include <iostream>
#include <boost/unordered/unordered_map.hpp>
#include <boost/python/numpy.hpp>

#include "checks.hpp"

using namespace std;

using namespace boost;
using namespace espressopp::iterator;
namespace np = boost::python::numpy;

namespace espressopp
{
namespace storage
{
LOG4ESPP_LOGGER(Storage::logger, "Storage");

const int STORAGE_COMM_TAG = 0xaa;

const int Storage::dataOfUpdateGhosts = 0;
const int Storage::dataOfExchangeGhosts = DATA_PROPERTIES;

Storage::Storage(std::shared_ptr<System> system, int halfCellInt)
    : SystemAccess(system),
      halfCellInt(halfCellInt),
      inBuffer(*system->comm),
      outBuffer(*system->comm)
{
    // logger.setLevel(log4espp::Logger::TRACE);
    LOG4ESPP_INFO(logger, "Created new storage object for a system, has buffers");
}

Storage::~Storage() {}

longint Storage::getNRealParticles() const
{
    longint cnt = 0;
    for (CellList::const_iterator it = realCells.begin(), end = realCells.end(); it != end; ++it)
    {
        longint size = (*it)->particles.size();
        if (size)
        {
            LOG4ESPP_TRACE(logger, "cell " << ((*it) - getFirstCell()) << " size " << size);
        }
        cnt += size;
    }
    return cnt;
}

longint Storage::getNLocalParticles() const
{
    longint cnt = 0;
    for (CellList::const_iterator it = localCells.begin(), end = localCells.end(); it != end; ++it)
    {
        longint size = (*it)->particles.size();
        if (size)
        {
            LOG4ESPP_TRACE(logger, "cell " << ((*it) - getFirstCell()) << " size " << size);
        }
        cnt += size;
    }
    return cnt;
}
longint Storage::getNGhostParticles() const
{
    longint cnt = 0;
    for (CellList::const_iterator it = ghostCells.begin(), end = ghostCells.end(); it != end; ++it)
    {
        longint size = (*it)->particles.size();
        if (size)
        {
            LOG4ESPP_TRACE(logger, "cell " << ((*it) - getFirstCell()) << " size " << size);
        }
        cnt += size;
    }
    return cnt;
}
longint Storage::getNAdressParticles() const { return localAdrATParticles.size(); }

python::list Storage::getRealParticleIDs()
{
    python::list pids;
    for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
    {
        pids.append(cit->getId());
    }
    return pids;
}

// TODO find out why python crashes if inlined
// inline
void Storage::removeFromLocalParticles(Particle *p, bool weak)
{
    /* no pointer left, can happen for ghosts when the real particle
       e has already been removed */
    if (localParticles.find(p->id()) == localParticles.end())
    {
        return;
    }

    if (!weak || localParticles[p->id()] == p)
    {
        LOG4ESPP_TRACE(logger, "removing local pointer for particle id=" << p->id() << " @ " << p);
        localParticles.erase(p->id());
    }
    else
    {
        LOG4ESPP_TRACE(logger, "NOT removing local pointer for particle id="
                                   << p->id() << " @ " << p << " since pointer is @ "
                                   << localParticles[p->id()]);
    }
}

/** Scale coordinates of all real particles by factor s */
void Storage::scaleVolume(real s)
{
    for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
    {
        Real3D pos = cit->getPos();
        pos *= s;
        cit->setPos(pos);
    }
}
/** Scale coordinates of all real particles by factor s. Anisotropic case */
void Storage::scaleVolume(Real3D s)
{
    for (CellListIterator cit(realCells); !cit.isDone(); ++cit)
    {
        Real3D pos = cit->getPos();
        pos[0] *= s[0];
        pos[1] *= s[1];
        pos[2] *= s[2];
        cit->setPos(pos);
    }
}

void Storage::removeAdrATParticle(longint id)
{
    if (localAdrATParticles.find(id) == localAdrATParticles.end())
    {
        std::cout << "not removing AT particle " << id << ", since not found \n";
        return;
    }

    // remove from ParticleList (vector)
    Particle *dbegin = &AdrATParticles.front();  // see whether the array was moved
    Particle *p = lookupAdrATParticle(id);
    int i = p - &AdrATParticles[0];
    int newSize = AdrATParticles.size() - 1;
    if (i != newSize)
    {  // if we are not removing the particle in last place
        AdrATParticles[i] = AdrATParticles.back();
    }
    AdrATParticles.resize(newSize);

    // remove from particle map
    localAdrATParticles.erase(id);

    // update particle map with particle list
    if (dbegin != &AdrATParticles.front())
    {
        updateLocalParticles(AdrATParticles, true);
    }
    else if (i != newSize)
    {
        Particle *np = &(AdrATParticles[i]);
        updateInLocalAdrATParticles(np);
    }
}

// TODO find out why python crashes if inlined
// inline
void Storage::updateInLocalParticles(Particle *p, bool weak)
{
    if (!weak || localParticles.find(p->id()) == localParticles.end())
    {
        LOG4ESPP_TRACE(logger, "updating local pointer for particle id=" << p->id() << " @ " << p);

        localParticles[p->id()] = p;

        /*
        // AdResS testing TODO
        Particle* oldp = localParticles[p->id()];
        Particle* newp = p;

        localParticles[p->id()] = p;

        //std::cout << "old *p " << oldp << "\n";
        //std::cout << "new *p " << newp << "\n";

        // TODO reorganize
        if (oldp && fixedtupleList) {
            FixedTupleList::iterator it;
            it = fixedtupleList->find(oldp);
            if (it != fixedtupleList->end()) {
                std::vector<Particle*> tmp;
                tmp = it->second;
                fixedtupleList->insert(std::make_pair(newp, tmp));
                fixedtupleList->erase(it);
            }
            else {
                //std::cout << "updateInLocalParticles: Particle not found in tuples!\n";
            }
        }
        */
    }
    else
    {
        LOG4ESPP_TRACE(logger, "NOT updating local pointer for particle id="
                                   << p->id() << " @ " << p << " has already pointer @ "
                                   << localParticles[p->id()]);
    }
}

inline void Storage::updateInLocalAdrATParticles(Particle *p) { localAdrATParticles[p->id()] = p; }

void Storage::updateLocalParticles(ParticleList &list, bool adress)
{
    if (adress)
    {
        for (ParticleList::Iterator it(list); it.isValid(); ++it)
        {
            updateInLocalAdrATParticles(&(*it));
        }
    }
    else
    {
        for (ParticleList::Iterator it(list); it.isValid(); ++it)
        {
            updateInLocalParticles(&(*it));
        }
    }
}

void Storage::resizeCells(longint nCells)
{
    cells.resize(nCells);
    localCells.reserve(nCells);
    for (LocalCellList::iterator it = cells.begin(), end = cells.end(); it != end; ++it)
    {
        localCells.push_back(&(*it));
    }
}

Particle *Storage::addParticle(longint id, const Real3D &p, bool checkIfRealParticle)
{
    if (checkIfRealParticle && !checkIsRealParticle(id, p))
    {
        return nullptr;
    }

    Cell *cell;

    Particle n;
    n.init();
    n.id() = id;
    n.position() = p;
    n.image() = Int3D(0);
    getSystem()->bc->foldPosition(n.position(), n.image());
    cell = mapPositionToCellClipped(n.position());

    // std::cout << "add particle: " << n.id() << " (" << n.position() << ")\n";

    appendIndexedParticle(cell->particles, n);

    LOG4ESPP_TRACE(logger, "got particle id =" << id << " @ " << p << " ; put it into cell "
                                               << cell - getFirstCell());
    LOG4ESPP_TRACE(logger, "folded it to " << n.position()[0] << " " << n.position()[1] << " "
                                           << n.position()[2]);
    LOG4ESPP_TRACE(logger, "cell size is now " << cell->particles.size());

    return &cell->particles.back();
}

int Storage::removeParticle(longint id)
{
    Particle *p = lookupRealParticle(id);
    if (p)
    {
        Cell *cell = mapPositionToCellChecked(p->position());

        removeFromLocalParticles(p);

        int rem_part_ind = p - &cell->particles[0];
        int newSize = cell->particles.size() - 1;
        if (rem_part_ind != newSize)
        {
            cell->particles[rem_part_ind] = cell->particles.back();
        }
        cell->particles.resize(newSize);

        updateLocalParticles(cell->particles);

        onParticlesChanged();
        Particle *p1 = lookupRealParticle(id);
        if (p1)
        {
            // it should not be printed out
            std::cout << "Part still exst. pid: " << id << "  rpid: " << p1->id() << std::endl;
        }
        return 1;
    }
    else
    {
        return 0;
    }
    // TODO particle should be removed from different particle groups and lists too
}

void Storage::removeAllParticles()
{
    localParticles.clear();
    for (CellList::iterator it = localCells.begin(), end = localCells.end(); it != end; ++it)
    {
        (*it)->particles.clear();
    }
    onParticlesChanged();
}

Particle *Storage::addAdrATParticle(longint id, const Real3D &p, const Real3D &_vpp)
{
    if (!checkIsRealParticle(id, _vpp))
    {
        return static_cast<Particle *>(0);
    }

    Particle n;
    n.init();
    n.id() = id;
    n.position() = p;
    n.image() = Int3D(0);

    // std::cout << "add ATparticle: " << n.id() << " (" << n.position() << ")\n";

    // fold AT particles for same amount as VP
    Real3D vpp_old = _vpp;
    Real3D vpp_new = _vpp;

    getSystem()->bc->foldPosition(vpp_new);

    if (vpp_old != vpp_new)
    {
        // std::cout << "VP old pos (" << vpp_old << ") ";
        // std::cout << "new pos (" << vpp_new << ")\n";

        Real3D moved = vpp_old - vpp_new;
        n.position() = n.position() - moved;

        // std::cout << " Moved AT particle to " << n.position() << "\n";
    }

    int last_capacity = AdrATParticles.capacity();

    AdrATParticles.push_back(n);
    Particle *local = &AdrATParticles.back();

    if (last_capacity != int_c(AdrATParticles.capacity()))
    {
        updateLocalParticles(AdrATParticles, true);
    }
    else
    {
        updateInLocalAdrATParticles(local);
    }

    return local;
}

// this is called from fixedtuplelist only!
Particle *Storage::addAdrATParticleFTPL(Particle n)
{
    int last_capacity = AdrATParticles.capacity();

    AdrATParticles.push_back(n);
    Particle *local = &AdrATParticles.back();

    if (last_capacity != int_c(AdrATParticles.capacity()))
    {
        updateLocalParticles(AdrATParticles, true);
    }
    else
    {
        updateInLocalAdrATParticles(local);
    }

    return local;
}

/*Particle* Storage::addParticle(longint id, const Real3D& p, int type) {
    Particle* pt = addParticle(id, p);
    pt->setType(type);
    return pt;
}*/

Particle *Storage::appendUnindexedParticle(ParticleList &l, Particle &part)
{
    l.push_back(part);
    return &l.back();
}

// Particle *Storage::appendUnindexedAdrParticle(ParticleListAdr &l, Particle &part)
Particle *Storage::appendUnindexedAdrParticle(ParticleList &l, Particle &part)
{
    l.push_back(part);
    return &l.back();
}

void Storage::appendParticleListToGhosts(ParticleList &l) { AdrATParticlesG.push_back(l); }

Particle *Storage::appendIndexedParticle(ParticleList &l, Particle &part)
{
    int last_capacity = l.capacity();

    l.push_back(part);
    Particle *p = &l.back();

    if (last_capacity != int_c(l.capacity()))
    {
        updateLocalParticles(l);
    }
    else
    {
        updateInLocalParticles(p);
    }

    return p;
}

Particle *Storage::moveIndexedParticle(ParticleList &dl, ParticleList &sl, int i)
{
    int dlast_capacity = dl.capacity();
    int slast_capacity = sl.capacity();

    dl.push_back(sl[i]);
    int newSize = sl.size() - 1;
    if (i != newSize)
    {
        sl[i] = sl.back();
    }
    sl.resize(newSize);

    Particle *dst = &dl.back();

    // fix up destination list
    if (dlast_capacity != int_c(dl.capacity()))
    {
        updateLocalParticles(dl);
    }
    else
    {
        updateInLocalParticles(dst);
    }

    // fix up resorted source list; due to moving, the last particle
    // might have been moved to the position of the actually moved one
    if (slast_capacity != int_c(sl.capacity()))
    {
        updateLocalParticles(sl);
    }
    else if (i != newSize)
    {
        Particle *src = &(sl[i]);
        updateInLocalParticles(src);
    }

    return dst;
}

void Storage::fetchParticles(Storage &old)
{
    LOG4ESPP_DEBUG(logger, "number of received cells = " << old.getRealCells().size());

    for (CellListIterator it(old.getRealCells()); it.isValid(); ++it)
    {
        Particle &part = *it;
        Cell *nc = mapPositionToCellClipped(part.position());
        appendUnindexedParticle(nc->particles, part);
    }

    // update localParticles
    for (CellList::Iterator it(realCells); it.isValid(); ++it)
    {
        updateLocalParticles((*it)->particles);
    }
}

void Storage::sendParticles(ParticleList &list, longint node)
{
    LOG4ESPP_DEBUG(logger, "send " << list.size() << " particles to " << node);

    // pack for transport

    OutBuffer &data = outBuffer;

    data.reset();
    int size = list.size();
    data.write(size);
    for (ParticleList::Iterator it(list); it.isValid(); ++it)
    {
        removeFromLocalParticles(&(*it));
        data.write(*it);
    }

    beforeSendParticles(list, data);  // this also takes care of AdResS AT Particles

    list.clear();

    // ... and send
    data.send(node, STORAGE_COMM_TAG);

    LOG4ESPP_DEBUG(logger, "done");
}

void Storage::recvParticles(ParticleList &list, longint node)
{
    LOG4ESPP_DEBUG(logger, "recv from " << node);

    InBuffer &data = inBuffer;  // reuse storage buffer

    data.recv(node, STORAGE_COMM_TAG);

    // ... and unpack
    int size;
    data.read(size);
    int curSize = list.size();
    LOG4ESPP_DEBUG(logger, "got " << size << " particles, have " << curSize);

    if (size > 0)
    {
        list.resize(curSize + size);

        for (int i = 0; i < size; ++i)
        {
            Particle *p = &list[curSize + i];
            data.read(*p);
            updateInLocalParticles(p);
        }

        afterRecvParticles(list, data);  // this also takes care of AdResS AT Particles
    }

    LOG4ESPP_DEBUG(logger, "done");
}

void Storage::invalidateGhosts()
{
    for (CellListIterator it(getGhostCells()); it.isValid(); ++it)
    {
        /* remove only ghosts from the hash if the localParticles hash
          actually points to the ghost.  If there are local ghost cells
          to implement pbc, the real particle will be the one accessible
          via localParticles.
        */
        removeFromLocalParticles(&(*it), true);
    }
}

void Storage::decompose()
{
    invalidateGhosts();
    decomposeRealParticles();
    exchangeGhosts();
    onParticlesChanged();
}

void Storage::packPositionsEtc(OutBuffer &buf, Cell &_reals, int extradata, const Real3D &shift)
{
    ParticleList &reals = _reals.particles;

    LOG4ESPP_DEBUG(logger, "pack data from reals in " << (&_reals - getFirstCell()));
    LOG4ESPP_DEBUG(logger, "also packing " << ((extradata & DATA_PROPERTIES) ? "properties " : "")
                                           << ((extradata & DATA_MOMENTUM) ? "momentum " : "")
                                           << ((extradata & DATA_LOCAL) ? "local " : ""));
    LOG4ESPP_DEBUG(logger,
                   "positions are shifted by " << shift[0] << "," << shift[1] << "," << shift[2]);

    for (ParticleList::iterator src = reals.begin(), end = reals.end(); src != end; ++src)
    {
        buf.write(*src, extradata, shift);
    }
}

void Storage::unpackPositionsEtc(Cell &_ghosts, InBuffer &buf, int extradata)
{
    ParticleList &ghosts = _ghosts.particles;

    LOG4ESPP_DEBUG(logger, "unpack data to ghosts in " << (&_ghosts - getFirstCell()));
    LOG4ESPP_DEBUG(logger, "also unpacking " << ((extradata & DATA_PROPERTIES) ? "properties " : "")
                                             << ((extradata & DATA_MOMENTUM) ? "momentum " : "")
                                             << ((extradata & DATA_LOCAL) ? "local " : ""));

    for (ParticleList::iterator dst = ghosts.begin(), end = ghosts.end(); dst != end; ++dst)
    {
        buf.read(*dst, extradata);

        if (extradata & DATA_PROPERTIES)
        {
            updateInLocalParticles(&(*dst), true);
        }

        dst->ghost() = 1;
    }
}

void Storage::copyRealsToGhosts(Cell &_reals, Cell &_ghosts, int extradata, const Real3D &shift)
{
    ParticleList &reals = _reals.particles;
    ParticleList &ghosts = _ghosts.particles;

    LOG4ESPP_DEBUG(logger, "copy data from reals in " << (&_reals - getFirstCell())
                                                      << " to ghosts in "
                                                      << (&_ghosts - getFirstCell()));
    LOG4ESPP_DEBUG(logger, "also copying " << ((extradata & DATA_PROPERTIES) ? "properties " : "")
                                           << ((extradata & DATA_MOMENTUM) ? "momentum " : "")
                                           << ((extradata & DATA_LOCAL) ? "local " : ""));
    LOG4ESPP_DEBUG(logger,
                   "positions are shifted by " << shift[0] << "," << shift[1] << "," << shift[2]);

    ghosts.resize(reals.size());

    for (ParticleList::iterator src = reals.begin(), end = reals.end(), dst = ghosts.begin();
         src != end; ++src, ++dst)
    {
        dst->copyAsGhost(*src, extradata, shift);
    }
}

void Storage::packForces(OutBuffer &buf, Cell &_ghosts)
{
    LOG4ESPP_DEBUG(logger, "pack ghost forces to buffer from cell " << (&_ghosts - getFirstCell()));

    ParticleList &ghosts = _ghosts.particles;

    for (ParticleList::iterator src = ghosts.begin(), end = ghosts.end(); src != end; ++src)
    {
        buf.write(src->particleForce());

        LOG4ESPP_TRACE(logger, "from particle " << src->id() << ": packing force " << src->force());
    }
}

void Storage::unpackForces(Cell &_reals, InBuffer &buf)
{
    LOG4ESPP_DEBUG(logger, "add forces from buffer to cell " << (&_reals - getFirstCell()));

    ParticleList &reals = _reals.particles;

    for (ParticleList::iterator dst = reals.begin(), end = reals.end(); dst != end; ++dst)
    {
        ParticleForce f;
        buf.read(f);
        LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": unpacking force " << f.f);
        dst->particleForce() = f;
    }
}

void Storage::unpackAndAddForces(Cell &_reals, InBuffer &buf)
{
    LOG4ESPP_DEBUG(logger, "add forces from buffer to cell " << (&_reals - getFirstCell()));

    ParticleList &reals = _reals.particles;

    for (ParticleList::iterator dst = reals.begin(), end = reals.end(); dst != end; ++dst)
    {
        ParticleForce f;
        buf.read(f);
        LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": unpacking force " << f.f
                                               << " and adding to " << dst->force());
        dst->particleForce() += f;
    }
}

void Storage::addGhostForcesToReals(Cell &_ghosts, Cell &_reals)
{
    LOG4ESPP_DEBUG(logger, "add forces from ghosts in cell " << (&_ghosts - getFirstCell())
                                                             << " to reals in cell "
                                                             << (&_reals - getFirstCell()));

    ParticleList &reals = _reals.particles;
    ParticleList &ghosts = _ghosts.particles;

    for (ParticleList::iterator dst = reals.begin(), end = reals.end(), src = ghosts.begin();
         dst != end; ++dst, ++src)
    {
        LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": adding force " << src->force()
                                               << " to " << dst->force());

        dst->particleForce() += src->particleForce();
    }
}

void Storage::clearSavedPositions()
{
    savedRealPositions.clear();
    savedImages.clear();
}
void Storage::savePosition(size_t id)
{
    if (getSystemRef().comm->size() == 1)
    {
        Particle *p = lookupRealParticle(id);
        if (p)
        {
            savedRealPositions[id] = p->position();
            savedImages[id] = p->image();
        }
    }
    else
    {
        esutil::Error err(getSystem()->comm);
        stringstream msg;
        msg << " At the moment it works only for one CPU. One can not store old positions"
               " for several CPUs";
        err.setException(msg.str());
    }
}

void Storage::restorePositions()
{
    int count = 0;
    int totCount = 0;
    if (!savedRealPositions.empty())
    {
        for (map<size_t, Real3D>::iterator itr = savedRealPositions.begin();
             itr != savedRealPositions.end(); ++itr)
        {
            size_t id = itr->first;
            Particle *p = lookupRealParticle(id);
            if (p)
            {
                p->position() = itr->second;
            }
        }
        for (map<size_t, Int3D>::iterator itr = savedImages.begin(); itr != savedImages.end();
             ++itr)
        {
            size_t id = itr->first;
            Particle *p = lookupRealParticle(id);
            if (p)
            {
                p->image() = itr->second;
            }
        }
        count++;
    }

    mpi::all_reduce(*getSystem()->comm, count, totCount, std::plus<int>());

    if (totCount == 0)
    {
        esutil::Error err(getSystem()->comm);
        stringstream msg;

        msg << " There is nothing to restore. Check whether you saved positions";
        err.setException(msg.str());
    }
}

///////////////////////////////////////////////////////////////////////////
/// Faster adding of particles by storing data in numpy arrays
void partArrayCheck(np::ndarray const &part_arr)
{
    using namespace espressopp::python;
    if (!(part_arr.get_dtype() == numpy::dtype::get_builtin<real>()))
        throw std::runtime_error("Invalid dtype of part_arr");
    if (!(part_arr.get_flags() & numpy::ndarray::C_CONTIGUOUS))
        throw std::runtime_error("part_arr must have C_CONTIGUOUS flag");
}

void idxArrayCheck(np::ndarray const &idx_arr)
{
    using namespace espressopp::python;
    if (!(idx_arr.get_dtype() == numpy::dtype::get_builtin<int>()))
        throw std::runtime_error("Invalid dtype of idx_arr");
    if (idx_arr.shape(0) != 31)
        throw std::runtime_error("Invalid idx_arr shape. axis=0 must have shape=31");
}

void addParticlesCheck(np::ndarray const &part_arr, np::ndarray const &idx_arr)
{
    partArrayCheck(part_arr);
    idxArrayCheck(idx_arr);
}

void addParticlesFromArray(class Storage *obj,
                           np::ndarray const &part_arr,
                           np::ndarray const &idx_arr)
{
    addParticlesCheck(part_arr, idx_arr);

    const real *part = (real *)(part_arr.get_data());
    const int *idx = (int *)(idx_arr.get_data());
    const int npart = part_arr.shape(0);
    const int nidx = part_arr.shape(1);

    obj->addParticlesFromArrayImpl(part, idx, npart, nidx);
}

#define DEFINE_INDEX_VARS(PRE, IDX)        \
    const int PRE##_id = IDX[0];           \
    const int PRE##_posx = IDX[1];         \
    const int PRE##_posy = IDX[2];         \
    const int PRE##_posz = IDX[3];         \
    const int PRE##_modeposx = IDX[4];     \
    const int PRE##_modeposy = IDX[5];     \
    const int PRE##_modeposz = IDX[6];     \
    const int PRE##_vx = IDX[7];           \
    const int PRE##_vy = IDX[8];           \
    const int PRE##_vz = IDX[9];           \
    const int PRE##_modemomx = IDX[10];    \
    const int PRE##_modemomy = IDX[11];    \
    const int PRE##_modemomz = IDX[12];    \
    const int PRE##_fx = IDX[13];          \
    const int PRE##_fy = IDX[14];          \
    const int PRE##_fz = IDX[15];          \
    const int PRE##_fmx = IDX[16];         \
    const int PRE##_fmy = IDX[17];         \
    const int PRE##_fmz = IDX[18];         \
    const int PRE##_q = IDX[19];           \
    const int PRE##_radius = IDX[20];      \
    const int PRE##_fradius = IDX[21];     \
    const int PRE##_vradius = IDX[22];     \
    const int PRE##_type = IDX[23];        \
    const int PRE##_mass = IDX[24];        \
    const int PRE##_varmass = IDX[25];     \
    const int PRE##_adrAT = IDX[26];       \
    const int PRE##_lambda_adr = IDX[27];  \
    const int PRE##_lambda_adrd = IDX[28]; \
    const int PRE##_state = IDX[29];       \
    const int PRE##_pib = IDX[30];         \
    /*  */

void Storage::addParticlesFromArrayImpl(const real *part,
                                        const int *idx,
                                        const int npart,
                                        const int nidx)
{
    boost::mpi::communicator &comm = *(getSystem()->comm);
    int nidx_ = 0;
    for (int i = 0; i < 31; i++) nidx_ += (idx[i] >= 0);
    if (nidx_ != nidx) LOG4ESPP_ERROR(logger, "size mismatch in expected number of particles");

    DEFINE_INDEX_VARS(index, idx);

    {
        bool particleExists = false;
        for (int p = 0; p < npart; p++)
        {
            const longint id = part[p * nidx + index_id];
            Particle *lp = this->lookupRealParticle(id);
            if (lp != 0)
            {
                particleExists = true;
                std::cout << "[" << comm.rank() << "] WARNING: Particle " << id << " already exists"
                          << std::endl;
            }
        }
        const bool anyParticleExists =
            boost::mpi::all_reduce(comm, particleExists, std::logical_or<bool>());
        if (anyParticleExists)
        {
            std::cout
                << "[" << comm.rank()
                << "] WARNING: Some particles already exist. The list of particles was not added."
                << std::endl;
            return;
        }
    }

    if (index_adrAT >= 0)
    {
        LOG4ESPP_ERROR(logger, "Not implemented for adress");
    }

    for (int p = 0; p < npart; p++)
    {
        const int offset = p * nidx;
        const longint id = part[offset + index_id];
        Real3D pos(part[offset + index_posx], part[offset + index_posy], part[offset + index_posz]);
        Particle *sp = this->addParticle(id, pos);

        if (sp == NULL) continue;

        if ((index_vx >= 0) && (index_vy >= 0) && (index_vz >= 0))
            sp->velocity() =
                Real3D(part[offset + index_vx], part[offset + index_vy], part[offset + index_vz]);

        if ((index_modeposx >= 0) && (index_modeposy >= 0) && (index_modeposz >= 0))
            sp->modepos() = Real3D(part[offset + index_modeposx], part[offset + index_modeposy],
                                   part[offset + index_modeposz]);

        if ((index_modemomx >= 0) && (index_modemomy >= 0) && (index_modemomz >= 0))
            sp->modemom() = Real3D(part[offset + index_modemomx], part[offset + index_modemomy],
                                   part[offset + index_modemomz]);

        if ((index_fx >= 0) && (index_fy >= 0) && (index_fz >= 0))
            sp->force() =
                Real3D(part[offset + index_fx], part[offset + index_fy], part[offset + index_fz]);

        if ((index_fmx >= 0) && (index_fmy >= 0) && (index_fmz >= 0))
            sp->forcem() = Real3D(part[offset + index_fmx], part[offset + index_fmy],
                                  part[offset + index_fmz]);

        if (index_q >= 0) sp->q() = part[offset + index_q];

        if (index_radius >= 0) sp->radius() = part[offset + index_radius];

        if (index_fradius >= 0) sp->fradius() = part[offset + index_fradius];

        if (index_vradius >= 0) sp->vradius() = part[offset + index_vradius];

        if (index_type >= 0) sp->type() = int(part[offset + index_type]);

        if (index_pib >= 0) sp->pib() = int(part[offset + index_pib]);

        if (index_mass >= 0) sp->mass() = part[offset + index_mass];

        if (index_varmass >= 0) sp->varmass() = part[offset + index_varmass];

        if (index_lambda_adr >= 0) sp->lambda() = part[offset + index_lambda_adr];

        if (index_lambda_adrd >= 0) sp->lambdaDeriv() = part[offset + index_lambda_adrd];

        if (index_state >= 0) sp->state() = int(part[offset + index_state]);
    }
}

void addParticlesFromArrayReplicate(class Storage *obj,
                                    python::numpy::ndarray const &part_arr,
                                    python::numpy::ndarray const &idx_arr,
                                    const real Lx,
                                    const real Ly,
                                    const real Lz,
                                    const int xdim,
                                    const int ydim,
                                    const int zdim,
                                    const int pid_start)
{
    addParticlesCheck(part_arr, idx_arr);

    const real *part = (real *)(part_arr.get_data());
    const int *idx = (int *)(idx_arr.get_data());
    const int npart = part_arr.shape(0);
    const int nidx = part_arr.shape(1);

    obj->addParticlesFromArrayReplImpl(part, idx, npart, nidx, Lx, Ly, Lz, xdim, ydim, zdim,
                                       pid_start);
}

void Storage::addParticlesFromArrayReplImpl(const real *part,
                                            const int *idx,
                                            const int npart,
                                            const int nidx,
                                            const real Lx,
                                            const real Ly,
                                            const real Lz,
                                            const int xdim,
                                            const int ydim,
                                            const int zdim,
                                            const int pid_start)
{
    int nidx_ = 0;
    for (int i = 0; i < 31; i++) nidx_ += (idx[i] >= 0);
    if (nidx_ != nidx) LOG4ESPP_ERROR(logger, "size mismatch in expected number of particles");

    const int index_id = idx[0];
    const int index_posx = idx[1];
    const int index_posy = idx[2];
    const int index_posz = idx[3];

    if (pid_start < 0) throw std::runtime_error("Assertion failed: pid_start>=0");
    if (xdim < 1) throw std::runtime_error("Assertion failed: xdim>=1");
    if (ydim < 1) throw std::runtime_error("Assertion failed: ydim>=1");
    if (zdim < 1) throw std::runtime_error("Assertion failed: zdim>=1");

    typedef std::tuple<size_t, size_t, Real3D> repl_t;
    std::vector<repl_t> replicated;
    size_t newId = 0;
    for (int i = 0; i < xdim; i++)
    {
        for (int j = 0; j < ydim; j++)
        {
            for (int k = 0; k < zdim; k++)
            {
                for (int p = 0; p < npart; p++)
                {
                    const int offset = p * nidx;
                    const Real3D newPos(part[offset + index_posx] + i * Lx,
                                        part[offset + index_posy] + j * Ly,
                                        part[offset + index_posz] + k * Lz);
                    if (checkIsRealParticle(0, newPos))
                    {
                        const size_t oldId = p;
                        replicated.push_back({newId, oldId, newPos});
                    }
                    newId++;
                }
            }
        }
    }

    std::vector<real> partRepl(replicated.size() * nidx);
    for (size_t p = 0; p < replicated.size(); p++)
    {
        /// copy all entries
        const auto &r = replicated[p];
        const auto newId = std::get<0>(r);
        const auto oldId = std::get<1>(r);
        const auto &newPos = std::get<2>(r);

        /// copy all entries for this row
        const size_t newOffset = p * nidx;
        const size_t oldOffset = oldId * nidx;
        for (int ii = 0; ii < nidx; ii++)
        {
            partRepl[newOffset + ii] = part[oldOffset + ii];
        }

        /// replace id
        partRepl[newOffset + index_id] = newId;

        /// replace position
        partRepl[newOffset + index_posx] = newPos[0];
        partRepl[newOffset + index_posy] = newPos[1];
        partRepl[newOffset + index_posz] = newPos[2];
    }

    /// call:
    this->addParticlesFromArrayImpl(partRepl.data(), idx, replicated.size(), nidx);
}

np::ndarray getParticlesArrayImpl(class Storage *obj, np::ndarray const &idx_arr)
{
    idxArrayCheck(idx_arr);
    const int *idx = (int *)(idx_arr.get_data());
    int nidx = 0;
    for (int i = 0; i < 31; i++)
    {
        nidx += (idx[i] >= 0);
    }

    for (int i = 0; i < 31; i++)
    {
        ASSERT_LESS(idx[i], nidx)
    }

    DEFINE_INDEX_VARS(index, idx);

    int npart = 0;
    for (auto const &rc : obj->getRealCells())
    {
        npart += rc->particles.size();
    }

    auto part_arr = np::zeros(python::make_tuple(npart, nidx), np::dtype::get_builtin<real>());
    partArrayCheck(part_arr);
    real *part = (real *)(part_arr.get_data());

    if (index_adrAT >= 0)
    {
        std::cerr << "getParticlesArray not implemented for adress" << std::endl;
    }

    {
        int p = 0;
        for (CellListIterator cit(obj->getRealCells()); !cit.isDone(); ++cit, ++p)
        {
            const int offset = p * nidx;
            if (index_id >= 0) part[offset + index_id] = cit->id();
            if (index_posx >= 0) part[offset + index_posx] = cit->position()[0];
            if (index_posy >= 0) part[offset + index_posy] = cit->position()[1];
            if (index_posz >= 0) part[offset + index_posz] = cit->position()[2];
            if (index_modeposx >= 0) part[offset + index_modeposx] = cit->modepos()[0];
            if (index_modeposy >= 0) part[offset + index_modeposy] = cit->modepos()[1];
            if (index_modeposz >= 0) part[offset + index_modeposz] = cit->modepos()[2];
            if (index_vx >= 0) part[offset + index_vx] = cit->velocity()[0];
            if (index_vy >= 0) part[offset + index_vy] = cit->velocity()[1];
            if (index_vz >= 0) part[offset + index_vz] = cit->velocity()[2];
            if (index_modemomx >= 0) part[offset + index_modemomx] = cit->modemom()[0];
            if (index_modemomy >= 0) part[offset + index_modemomy] = cit->modemom()[1];
            if (index_modemomz >= 0) part[offset + index_modemomz] = cit->modemom()[2];
            if (index_fx >= 0) part[offset + index_fx] = cit->force()[0];
            if (index_fy >= 0) part[offset + index_fy] = cit->force()[1];
            if (index_fz >= 0) part[offset + index_fz] = cit->force()[2];
            if (index_fmx >= 0) part[offset + index_fmx] = cit->forcem()[0];
            if (index_fmy >= 0) part[offset + index_fmy] = cit->forcem()[1];
            if (index_fmz >= 0) part[offset + index_fmz] = cit->forcem()[2];
            if (index_q >= 0) part[offset + index_q] = cit->q();
            if (index_radius >= 0) part[offset + index_radius] = cit->radius();
            if (index_fradius >= 0) part[offset + index_fradius] = cit->fradius();
            if (index_vradius >= 0) part[offset + index_vradius] = cit->vradius();
            if (index_type >= 0) part[offset + index_type] = cit->type();
            if (index_mass >= 0) part[offset + index_mass] = cit->mass();
            if (index_varmass >= 0) part[offset + index_varmass] = cit->varmass();
            if (index_lambda_adr >= 0) part[offset + index_lambda_adr] = cit->lambda();
            if (index_lambda_adrd >= 0) part[offset + index_lambda_adrd] = cit->lambdaDeriv();
            if (index_state >= 0) part[offset + index_state] = cit->state();
            if (index_pib >= 0) part[offset + index_pib] = cit->pib();
        }
        CHECK_EQUAL(p, npart);
    }

    return part_arr;
}

///////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void Storage::registerPython()
{
    using namespace espressopp::python;
    class_<Storage, boost::noncopyable>("storage_Storage", no_init)
        .def("clearSavedPositions", &Storage::clearSavedPositions)
        .def("savePosition", &Storage::savePosition)
        .def("restorePositions", &Storage::restorePositions)
        .def("addParticle", &Storage::addParticle, return_value_policy<reference_existing_object>())
        .def("removeParticle", &Storage::removeParticle)
        .def("removeAllParticles", &Storage::removeAllParticles)
        .def("addAdrATParticle", &Storage::addAdrATParticle,
             return_value_policy<reference_existing_object>())
        .def("setFixedTuplesAdress", &Storage::setFixedTuplesAdress)
        //.def("addParticle", &Storage::addParticle, return_value_policy< reference_existing_object
        //>())
        .def("lookupLocalParticle", &Storage::lookupLocalParticle,
             return_value_policy<reference_existing_object>())
        .def("lookupRealParticle", &Storage::lookupRealParticle,
             return_value_policy<reference_existing_object>())
        .def("decompose", &Storage::decompose)
        .def("getRealParticleIDs", &Storage::getRealParticleIDs)
        .add_property("system", &Storage::getSystem)
        .def("addParticlesFromArray", &addParticlesFromArray)
        .def("addParticlesFromArrayReplicate", &addParticlesFromArrayReplicate)
        .def("getParticlesArrayImpl", &getParticlesArrayImpl);
}
}  // namespace storage
}  // namespace espressopp
