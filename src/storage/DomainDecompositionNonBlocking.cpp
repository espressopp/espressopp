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


//#include "python.hpp"
//#include <algorithm>
//#include <cmath>
//#include <sstream>
//#include "log4espp.hpp"
#include "System.hpp"
#include "Real3D.hpp"
#include "DomainDecomposition.hpp"
#include "DomainDecompositionNonBlocking.hpp"
#include "bc/BC.hpp"
#include "Int3D.hpp"
#include "Buffer.hpp"
#include "iterator/CellListIterator.hpp"

using namespace boost;
using namespace std;

namespace espressopp { 
  namespace storage {


  const int DD_COMM_TAG = 0xab;

  DomainDecompositionNonBlocking::
  DomainDecompositionNonBlocking(shared_ptr< System > _system,
          const Int3D& _nodeGrid,
          const Int3D& _cellGrid)
    : DomainDecomposition(_system, _nodeGrid, _cellGrid, 1 /*in DD nonblocking we do not allow halfcell at the moment*/),
      inBufferL(*_system->comm),
      inBufferR(*_system->comm),
      outBufferL(*_system->comm),
      outBufferR(*_system->comm),
      inBufferG(*_system->comm),
      outBufferG(*_system->comm) {}

  void DomainDecompositionNonBlocking::decomposeRealParticles() {

    //std::cout << getSystem()->comm->rank() << ": " << " decomposeRealParticles\n";

    LOG4ESPP_DEBUG(logger, "starting, expected comm buffer size " << exchangeBufferSize);

    // allocate send/recv buffers. We use the size as we need maximally so far, to avoid reallocation
    // TODO: This might be a problem when all particles are created on a single node initially!
    ParticleList sendBufL;
    sendBufL.reserve(exchangeBufferSize);
    ParticleList sendBufR;
    sendBufR.reserve(exchangeBufferSize);
    ParticleList recvBufL;
    recvBufL.reserve(exchangeBufferSize);
    ParticleList recvBufR;
    recvBufR.reserve(exchangeBufferSize);

    bool allFinished;
    do {
      bool finished = true;

      for (int coord = 0; coord < 3; ++coord) {
        LOG4ESPP_DEBUG(logger, "starting with direction " << coord);

        if (nodeGrid.getGridSize(coord) > 1) {
          for (std::vector<Cell*>::iterator it = realCells.begin(),
            end = realCells.end(); it != end; ++it) {

            Cell &cell = **it;

            // do not use an iterator here, since we need to take out particles during the loop
            for (size_t p = 0; p < cell.particles.size(); ++p) {
              Particle &part = cell.particles[p];
              const Real3D& pos = part.position();

              // check whether the particle is now "left" of the local domain
              if (pos[coord] - cellGrid.getMyLeft(coord) < -ROUND_ERROR_PREC) {
                LOG4ESPP_TRACE(logger, "send particle left " << part.id());
                moveIndexedParticle(sendBufL, cell.particles, p);
                // redo same particle since we took one out here, so it's a new one
                --p;
              }
              // check whether the particle is now "right" of the local domain
              else if (pos[coord] - cellGrid.getMyRight(coord) >= ROUND_ERROR_PREC) {
                LOG4ESPP_TRACE(logger, "send particle right " << part.id());
                moveIndexedParticle(sendBufR, cell.particles, p);
                --p;
              }
              // Sort particles in cells of this node during last direction
              else if (coord == 2) {
                const Real3D& pos = part.position();
                Cell *sortCell = mapPositionToCellChecked(pos);
                if (sortCell != &cell) {
                  if (sortCell == 0) {
                    // particle is not in the local domain
                    LOG4ESPP_DEBUG(logger, "take another loop: particle " << part.id()
                            << " @ " << pos <<
                            " is not inside node domain after neighbor exchange");
                    // isnan function is C99 only, x != x is only true if x == nan
                    if (pos[0] != pos[0] || pos[1] != pos[1] || pos[2] != pos[2]) {
                      // TODO: error handling
                      LOG4ESPP_ERROR(logger, "particle " << part.id() <<
                              " has moved to outer space (one or more coordinates are nan)");
                    } else {
                      // particle stays where it is, and will be sorted in the next round
                      finished = false;
                    }
                  } else {
                    // particle is in the local domain
                    moveIndexedParticle(sortCell->particles, cell.particles, p);
                    --p;
                  }
                }
              }

            }
          }

          mpi::request reqs[4];

          // Exchange particles, odd-even rule
          if (nodeGrid.getNodePosition(coord) % 2 == 0) {
            reqs[0]=isendParticles(outBufferL, sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
            reqs[1]=irecvParticles_initiate(  inBufferR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
            reqs[2]=isendParticles(outBufferR, sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
            reqs[3]=irecvParticles_initiate(  inBufferL, nodeGrid.getNodeNeighborIndex(2 * coord));
          } else {
            reqs[0]=irecvParticles_initiate(  inBufferR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
        	reqs[1]=isendParticles(outBufferL, sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
        	reqs[2]=irecvParticles_initiate(  inBufferL, nodeGrid.getNodeNeighborIndex(2 * coord));
        	reqs[3]=isendParticles(outBufferR, sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
          }

          mpi::wait_all(reqs, reqs + 4);

          if (nodeGrid.getNodePosition(coord) % 2 == 0) {
            irecvParticles_finish( inBufferR, recvBufR);
            irecvParticles_finish( inBufferL, recvBufL);
          } else {
            irecvParticles_finish( inBufferR, recvBufR);
        	irecvParticles_finish( inBufferL, recvBufL);
          }

          // sort received particles to cells
          if (appendParticles(recvBufL, 2 * coord) && coord == 2) finished = false;
          if (appendParticles(recvBufR, 2 * coord + 1) && coord == 2) finished = false;

          // reset send/recv buffers
          sendBufL.resize(0);
          sendBufR.resize(0);
          recvBufL.resize(0);
          recvBufR.resize(0);


        } else {
          /* Single node direction case (no communication)
            Fold particles that have left the box */
          for (std::vector< Cell* >::iterator it = realCells.begin(),
            end = realCells.end(); it != end; ++it) {
            Cell &cell = **it;
            // do not use an iterator here, since we have need to take out particles during the loop
            for (size_t p = 0; p < cell.particles.size(); ++p) {
                Particle &part = cell.particles[p];
                getSystem()->bc->foldCoordinate(part.position(), part.image(), coord);
                LOG4ESPP_TRACE(logger, "folded coordinate " << coord << " of particle " << part.id());

                if (coord == 2) {
                    Cell *sortCell = mapPositionToCellChecked(part.position());

                    if (sortCell != &cell) {
                        if (sortCell == 0) {
                            LOG4ESPP_DEBUG(logger, "take another loop: particle " << part.id()
                                    << " @ " << part.position()
                                    << " is not inside node domain after neighbor exchange");
                            const Real3D& pos = part.position();
                            // isnan function is C99 only, x != x is only true if x == nan
                            if (pos[0] != pos[0] || pos[1] != pos[1] || pos[2] != pos[2]) {
                                LOG4ESPP_ERROR(logger, "particle " << part.id() <<
                                        " has moved to outer space (one or more coordinates are nan)");
                            } else {
                                // particle stays where it is, and will be sorted in the next round
                                finished = false;
                            }
                        } else {
                            moveIndexedParticle(sortCell->particles, cell.particles, p);
                            --p;
                        }
                    }
                }
            }
          }
        }

        LOG4ESPP_DEBUG(logger, "done with direction " << coord);
      }

      // Communicate wether particle exchange is finished
      mpi::all_reduce(*getSystem()->comm, finished, allFinished, std::logical_and<bool>());
    } while (!allFinished);

    exchangeBufferSize = std::max(exchangeBufferSize,
                  std::max(sendBufL.capacity(),
                       std::max(sendBufR.capacity(),
                            std::max(recvBufL.capacity(),
                                 recvBufR.capacity()))));

    LOG4ESPP_DEBUG(logger, "finished exchanging particles, new send/recv buffer size " << exchangeBufferSize);

    LOG4ESPP_DEBUG(logger, "done");
  }

  void DomainDecompositionNonBlocking::
  doGhostCommunication(bool sizesFirst, bool realToGhosts, int extradata) {
    LOG4ESPP_DEBUG(logger, "do ghost communication " << (sizesFirst ? "with sizes " : "")
           << (realToGhosts ? "reals to ghosts " : "ghosts to reals ") << extradata);

    /* direction loop: x, y, z.
   Here we could in principle build in a one sided ghost
   communication, simply by taking the lr loop only over one
   value. */
    for (int _coord = 0; _coord < 3; ++_coord) {
      /* inverted processing order for ghost force communication,
        since the corner ghosts have to be collected via several
        nodes. We now add back the corner ghost forces first again
        to ghost forces, which only eventually go back to the real
        particle.
      */
      int coord = realToGhosts ? _coord : (2 - _coord);
      real curCoordBoxL = getSystem()->bc->getBoxL()[coord];

      // lr loop: left right
      for (int lr = 0; lr < 2; ++lr) {
        int dir         = 2 * coord + lr;
        int oppositeDir = 2 * coord + (1 - lr);

        Real3D shift(0, 0, 0);

        shift[coord] = nodeGrid.getBoundary(dir) * curCoordBoxL;

        LOG4ESPP_DEBUG(logger, "direction " << dir);

        if (nodeGrid.getGridSize(coord) == 1) {
          LOG4ESPP_DEBUG(logger, "local communication");

          // copy operation, we have to receive as many cells as we send
          if (commCells[dir].ghosts.size() != commCells[dir].reals.size()) {
            throw std::runtime_error("DomainDecomposition::doGhostCommunication: send/recv cell structure mismatch during local copy");
          }

          for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
            if (realToGhosts) {
              copyRealsToGhosts(*commCells[dir].reals[i], *commCells[dir].ghosts[i], extradata, shift);
            } else {
              addGhostForcesToReals(*commCells[dir].ghosts[i], *commCells[dir].reals[i]);
            }
          }
        }
        else {
          // exchange size information, if necessary
          if (sizesFirst) {
            LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes");

            // prepare buffers
            std::vector<longint> sendSizes, recvSizes;
            sendSizes.reserve(commCells[dir].reals.size());
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              sendSizes.push_back(commCells[dir].reals[i]->particles.size());
            }
            recvSizes.resize(commCells[dir].ghosts.size());

            mpi::request reqs[2];

            // exchange sizes, odd-even rule
            if (nodeGrid.getNodePosition(coord) % 2 == 0) {
              LOG4ESPP_DEBUG(logger, "sending to node " << nodeGrid.getNodeNeighborIndex(dir)
                        << ", then receiving from node " << nodeGrid.getNodeNeighborIndex(oppositeDir));
              reqs[0]=getSystem()->comm->isend(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG, &(sendSizes[0]), sendSizes.size());
              reqs[1]=getSystem()->comm->irecv(nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG, &(recvSizes[0]), recvSizes.size());
            }
            else {
              LOG4ESPP_DEBUG(logger, "receiving from node " << nodeGrid.getNodeNeighborIndex(oppositeDir)
                        << ", then sending to node " << nodeGrid.getNodeNeighborIndex(dir));
              reqs[0]=getSystem()->comm->irecv(nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG, &(recvSizes[0]), recvSizes.size());
              reqs[1]=getSystem()->comm->isend(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG, &(sendSizes[0]), sendSizes.size());
            }

            mpi::wait_all(reqs, reqs + 2);

            // resize according to received information
            for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
              commCells[dir].ghosts[i]->particles.resize(recvSizes[i]);
            }
            LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes done");
          }

          // prepare send and receive buffers
          longint receiver, sender;
          outBufferG.reset();
          if (realToGhosts) {
            receiver = nodeGrid.getNodeNeighborIndex(dir);
            sender = nodeGrid.getNodeNeighborIndex(oppositeDir);
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              packPositionsEtc(outBufferG, *commCells[dir].reals[i], extradata, shift);
            }
          }
          else {
            receiver = nodeGrid.getNodeNeighborIndex(oppositeDir);
            sender = nodeGrid.getNodeNeighborIndex(dir);
            for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
              packForces(outBufferG, *commCells[dir].ghosts[i]);
            }
          }

          mpi::request reqs[2];

          // exchange particles, odd-even rule
          if (nodeGrid.getNodePosition(coord) % 2 == 0) {
            reqs[0]=outBufferG.isend(receiver, DD_COMM_TAG);
            reqs[1]=inBufferG.irecv(sender, DD_COMM_TAG);
          } else {
            reqs[0]=inBufferG.irecv(sender, DD_COMM_TAG);
            reqs[1]=outBufferG.isend(receiver, DD_COMM_TAG);
          }

          mpi::wait_all(reqs, reqs + 2);

          // unpack received data
          if (realToGhosts) {
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              unpackPositionsEtc(*commCells[dir].ghosts[i], inBufferG, extradata);
            }
          }
          else {
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              unpackAndAddForces(*commCells[dir].reals[i], inBufferG);
            }
          }
        }
      }
    }
    LOG4ESPP_DEBUG(logger, "ghost communication finished");
  }

  mpi::request DomainDecompositionNonBlocking::isendParticles(OutBuffer &data, ParticleList &list, longint node)
  {
    LOG4ESPP_DEBUG(logger, "initiate non blocking isend " << list.size() << " particles to " << node);
    data.reset();
    int size = list.size();
    data.write(size);
    for (ParticleList::Iterator it(list); it.isValid(); ++it) {
        removeFromLocalParticles(&(*it));
        data.write(*it);
    }
    beforeSendParticles(list, data); // this also takes care of AdResS AT Particles
    list.clear();

    // ... and send
    return data.isend(node, DD_COMM_TAG);
  }

  mpi::request DomainDecompositionNonBlocking::irecvParticles_initiate(InBuffer& data, longint node)
  {
    LOG4ESPP_DEBUG(logger, "initiate non blocking irecv on " << node);
    return data.irecv(node, DD_COMM_TAG);
  }

  void DomainDecompositionNonBlocking::irecvParticles_finish(InBuffer &data, ParticleList &list)
  {
    LOG4ESPP_DEBUG(logger, "finish non blocking irecv");
    // ... and unpack
    int size;
    data.read(size);
    int curSize = list.size();
    LOG4ESPP_DEBUG(logger, "got " << size << " particles, have " << curSize);

    if (size > 0) {
      list.resize(curSize + size);

      for (int i = 0; i < size; ++i) {
        Particle *p = &list[curSize + i];
        data.read(*p);
        updateInLocalParticles(p);
      }

      afterRecvParticles(list, data); // this also takes care of AdResS AT Particles
    }
    LOG4ESPP_DEBUG(logger, "done");
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void DomainDecompositionNonBlocking::registerPython() {
    using namespace espressopp::python;
    class_< DomainDecompositionNonBlocking, bases< DomainDecomposition >, boost::noncopyable >
    ("storage_DomainDecompositionNonBlocking", init< shared_ptr< System >, const Int3D&, const Int3D& >())
    ;
  }


  }
}
