/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  Copyright (C) 2019
      Max Planck Computing and Data Facility
  Copyright (C) 2022
      Data Center, Johannes Gutenberg University Mainz

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

#include <algorithm>
#include <cmath>
#include <sstream>
#include "log4espp.hpp"
#include "System.hpp"

#include "Real3D.hpp"
#include "DomainDecomposition.hpp"
#include "bc/BC.hpp"
#include "Int3D.hpp"
#include "Buffer.hpp"

#include "iterator/CellListIterator.hpp"
#include "esutil/Error.hpp"

#include "boost/serialization/vector.hpp"

using namespace boost;
using namespace std;

namespace espressopp
{
namespace storage
{
const int DD_COMM_TAG = 0xab;

LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

std::string formatMismatchMessage(const Int3D& gridRequested, int nodesAvailable)
{
    std::ostringstream out;
    out << "requested node grid (" << gridRequested
        << ") does not match number of nodes in the communicator (" << nodesAvailable << ")";
    return out.str();
}

NodeGridMismatch::NodeGridMismatch(const Int3D& gridRequested, int nodesAvailable)
    : std::invalid_argument(formatMismatchMessage(gridRequested, nodesAvailable))
{
}

DomainDecomposition::DomainDecomposition(std::shared_ptr<System> _system,
                                         const Int3D& _nodeGrid,
                                         const Int3D& _cellGrid,
                                         int _halfCellInt)
    : Storage(_system, _halfCellInt), exchangeBufferSize(0)
{
    LOG4ESPP_INFO(logger, "node grid = " << _nodeGrid[0] << "x" << _nodeGrid[1] << "x"
                                         << _nodeGrid[2] << " cell grid = " << _cellGrid[0] << "x"
                                         << _cellGrid[1] << "x" << _cellGrid[2]);

    _system->NGridSize = _nodeGrid;
    createCellGrid(_nodeGrid, _cellGrid);
    initCellInteractions();
    prepareGhostCommunication();
    LOG4ESPP_DEBUG(logger, "done");
}

void DomainDecomposition::createCellGrid(const Int3D& _nodeGrid, const Int3D& _cellGrid)
{
    real myLeft[3];
    real myRight[3];

    nodeGrid = NodeGrid(_nodeGrid, getSystem()->comm->rank(), getSystem()->bc->getBoxL());

    if (nodeGrid.getNumberOfCells() != getSystem()->comm->size())
    {
        throw NodeGridMismatch(_nodeGrid, getSystem()->comm->size());
    }

    LOG4ESPP_INFO(logger, "my node grid position: " << nodeGrid.getNodePosition(0) << " "
                                                    << nodeGrid.getNodePosition(1) << " "
                                                    << nodeGrid.getNodePosition(2) << " -> "
                                                    << getSystem()->comm->rank());

    LOG4ESPP_DEBUG(logger, "my neighbors: " << nodeGrid.getNodeNeighborIndex(0) << "<->"
                                            << nodeGrid.getNodeNeighborIndex(1) << ", "
                                            << nodeGrid.getNodeNeighborIndex(2) << "<->"
                                            << nodeGrid.getNodeNeighborIndex(3) << ", "
                                            << nodeGrid.getNodeNeighborIndex(4) << "<->"
                                            << nodeGrid.getNodeNeighborIndex(5));

    for (int i = 0; i < 3; ++i)
    {
        myLeft[i] = nodeGrid.getMyLeft(i);
        myRight[i] = nodeGrid.getMyRight(i);
    }

    cellGrid = CellGrid(_cellGrid, myLeft, myRight, halfCellInt);

    LOG4ESPP_INFO(logger, "local box " << myLeft[0] << "-" << myRight[0] << ", " << myLeft[1] << "-"
                                       << myRight[1] << ", " << myLeft[2] << "-" << myRight[2]);

    longint nLocalCells = 1;
    longint nRealCells = 1;
    for (int i = 0; i < 3; ++i)
    {
        nRealCells *= cellGrid.getGridSize(i);
        nLocalCells *= cellGrid.getFrameGridSize(i);
    }

    resizeCells(nLocalCells);

    realCells.reserve(nRealCells);
    ghostCells.reserve(nLocalCells - nRealCells);

    markCells();

    LOG4ESPP_DEBUG(logger, "total # cells=" << nLocalCells << ", frame cell grid = ("
                                            << cellGrid.getFrameGridSize(0) << ", "
                                            << cellGrid.getFrameGridSize(1) << ", "
                                            << cellGrid.getFrameGridSize(2) << ")");
}

void DomainDecomposition::markCells()
{
    realCells.resize(0);
    ghostCells.resize(0);

    for (int o = 0; o < cellGrid.getFrameGridSize(2); ++o)
    {
        for (int n = 0; n < cellGrid.getFrameGridSize(1); ++n)
        {
            for (int m = 0; m < cellGrid.getFrameGridSize(0); ++m)
            {
                Cell* cur = &cells[cellGrid.mapPositionToIndex(m, n, o)];
                if (cellGrid.isInnerCell(m, n, o))
                {
                    LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is inner cell (" << m
                                                   << ", " << n << ", " << o << ")");
                    realCells.push_back(cur);
                }
                else
                {
                    LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is ghost cell (" << m
                                                   << ", " << n << ", " << o << ")");
                    ghostCells.push_back(cur);
                }
            }
        }
    }
}

// TODO one should take care of rc and system size
/** scale position coordinates of all real particles by factor s */
void DomainDecomposition::scaleVolume(real s, bool particleCoordinates)
{
    if (particleCoordinates) Storage::scaleVolume(s);

    real maxCut = getSystem()->maxCutoff;
    real skinL = getSystem()->getSkin();
    real cs = maxCut + skinL;
    if (cs > s * cellGrid.getSmallestCellDiameter())
    {
        Real3D Li = getSystem()->bc->getBoxL();  // getting the system size
        real minL = min(Li[0], min(Li[1], Li[2]));
        if (cs > minL)
        {
            esutil::Error err(getSystemRef().comm);
            stringstream msg;
            msg << "Error. The current system size " << minL << " smaller then cutoff+skin " << cs;
            err.setException(msg.str());
        }
        else
        {
            cellAdjust();
        }
    }
    else
    {
        cellGrid.scaleVolume(s);
        nodeGrid.scaleVolume(s);
    }
}
// anisotropic version
void DomainDecomposition::scaleVolume(Real3D s, bool particleCoordinates)
{
    if (particleCoordinates) Storage::scaleVolume(s);

    real maxCut = getSystem()->maxCutoff;
    real skinL = getSystem()->getSkin();
    real cs = maxCut + skinL;
    real cellD = cellGrid.getSmallestCellDiameter();

    real r0 = s[0] * cellD;
    real r1 = s[1] * cellD;
    real r2 = s[2] * cellD;

    if (cs > min(min(r0, r1), r2))
    {
        Real3D Li = getSystem()->bc->getBoxL();  // getting the system size
        real minL = min(Li[0], min(Li[1], Li[2]));
        if (cs > minL)
        {
            esutil::Error err(getSystemRef().comm);
            stringstream msg;
            msg << "Error. The current system size " << minL << " smaller then cutoff+skin " << cs;
            err.setException(msg.str());
        }
        else
            cellAdjust();
    }
    else
    {
        cellGrid.scaleVolume(s);
        nodeGrid.scaleVolume(s);
    }
}

Int3D DomainDecomposition::getInt3DCellGrid()
{
    return Int3D(cellGrid.getGridSize(0), cellGrid.getGridSize(1), cellGrid.getGridSize(2));
}
Int3D DomainDecomposition::getInt3DNodeGrid()
{
    return Int3D(nodeGrid.getGridSize(0), nodeGrid.getGridSize(1), nodeGrid.getGridSize(2));
}

void DomainDecomposition::cellAdjust()
{
    // create an appropriate cell grid
    Real3D box_sizeL = getSystem()->bc->getBoxL();
    real skinL = getSystem()->getSkin();
    real maxCutoffL = getSystem()->maxCutoff;

    // nodeGrid is already defined
    Int3D _nodeGrid(nodeGrid.getGridSize());
    // new cellGrid

    // skin_mod
    real rc_skin =
        (maxCutoffL * 2.0 > maxCutoffL + skinL ? (maxCutoffL * 2.0 + 0.01) : (maxCutoffL + skinL));
    int ix = (int)(box_sizeL[0] / (rc_skin * _nodeGrid[0]));
    int iy = (int)(box_sizeL[1] / (rc_skin * _nodeGrid[1]));
    int iz = (int)(box_sizeL[2] / (rc_skin * _nodeGrid[2]));
    Int3D _newCellGrid(ix, iy, iz);

    // if (getSystem()->comm->rank() == 0)
    //  std::cout << " Corrected DOMDEC [" << getInt3DNodeGrid() << "](" << _newCellGrid << ") \n";

    // save all particles to temporary vector
    std::vector<ParticleList> tmp_pl;
    size_t _N = realCells.size();
    tmp_pl.reserve(_N);
    for (CellList::Iterator it(realCells); it.isValid(); ++it)
    {
        tmp_pl.push_back((*it)->particles);
    }

    // reset all cells info
    invalidateGhosts();
    cells.clear();
    localCells.clear();
    realCells.clear();
    ghostCells.clear();
    for (int i = 0; i < 6; i++)
    {
        commCells[i].reals.clear();
        commCells[i].ghosts.clear();
    }

    // creating new grids
    createCellGrid(_nodeGrid, _newCellGrid);
    initCellInteractions();
    prepareGhostCommunication();

    // pushing the particles back to the empty cell ("do we have to check particles?")
    for (int i = 0; i < int_c(tmp_pl.size()); i++)
    {
        for (size_t p = 0; p < tmp_pl[i].size(); ++p)
        {
            Particle& part = tmp_pl[i][p];
            const Real3D& pos = part.position();
            Cell* sortCell = mapPositionToCellClipped(pos);
            appendUnindexedParticle(sortCell->particles, part);
        }
    }

    for (CellList::Iterator it(realCells); it.isValid(); ++it)
    {
        updateLocalParticles((*it)->particles);
    }

    exchangeGhosts();

    /// modify cell structure first before resorting
    /// particles and rebuilding neighbor lists
    onCellAdjust();

    onParticlesChanged();
}

void DomainDecomposition::initCellInteractions()
{
    LOG4ESPP_DEBUG(logger, "setting up neighbors for " << cells.size() << " cells");

    for (int o = cellGrid.getInnerCellsBegin(2); o < cellGrid.getInnerCellsEnd(2); ++o)
    {
        for (int n = cellGrid.getInnerCellsBegin(1); n < cellGrid.getInnerCellsEnd(1); ++n)
        {
            for (int m = cellGrid.getInnerCellsBegin(0); m < cellGrid.getInnerCellsEnd(0); ++m)
            {
                longint cellIdx = cellGrid.mapPositionToIndex(m, n, o);
                Cell* cell = &cells[cellIdx];

                LOG4ESPP_TRACE(logger, "setting up neighbors for cell "
                                           << cell - getFirstCell() << " @ " << m << " " << n << " "
                                           << o);

                cell->neighborCells.reserve(
                    (2 * halfCellInt + 1) * (2 * halfCellInt + 1) * (2 * halfCellInt + 1) - 1);

                // loop all neighbor cells
                for (int p = o - halfCellInt; p <= o + halfCellInt; ++p)
                {
                    for (int q = n - halfCellInt; q <= n + halfCellInt; ++q)
                    {
                        for (int r = m - halfCellInt; r <= m + halfCellInt; ++r)
                        {
                            if (p != o || q != n || r != m)
                            {
                                longint cell2Idx = cellGrid.mapPositionToIndex(r, q, p);
                                Cell* cell2 = &cells[cell2Idx];
                                cell->neighborCells.push_back(
                                    NeighborCellInfo(cell2, (cell2Idx < cellIdx)));

                                LOG4ESPP_TRACE(
                                    logger, "neighbor cell "
                                                << cell2 - getFirstCell() << " @ " << r << " " << q
                                                << " " << p
                                                << ((cell2Idx < cellIdx) ? " is" : " is not")
                                                << " taken");
                            }
                        }
                    }
                }
            }
        }
    }

    LOG4ESPP_DEBUG(logger, "done");
}

void DomainDecomposition::remapNeighbourCells(int cell_shift)
{
    // if (rename("FLAG_P","FLAG_P")==0 && getSystem()->comm->rank()==getSystem()->irank)
    // std::cout<<"SHIFT> "<<getSystem()->ghostShift<<" \n";
    // cell_shift=1: right shift for top ghost layer; cell_shift=-1: left shift
    LOG4ESPP_DEBUG(logger, (cell_shift < 0 ? "Left " : "Right ")
                               << "translation for ghost cells above the top layer\n");
    LOG4ESPP_DEBUG(logger, (cell_shift < 0 ? "Right " : "Left ")
                               << "translation for ghost cells under the bottom layer\n");
    int o, o_opp;
    real Lx = getSystem()->bc->getBoxL()[0];
    Cell *cell1, *cell2;
    longint cell_idx1, cell_idx2;
    int xgrid = getInt3DCellGrid()[0];  //*getInt3DNodeGrid()[0];

    if (cellGrid.getInnerCellsEnd(2) - cellGrid.getInnerCellsBegin(2) <= 0)
        throw std::runtime_error("remapNeighbourCells error: not working with <=1 cells on Z-dir");

    if (cell_shift == 0)
        throw std::runtime_error(
            "remapNeighbourCells error: The value of shift should not be zero! (Have you set up a "
            "second shear flow on a reversed direction?)\n");

    if (halfCellInt > 1)
        throw std::runtime_error(
            "remapNeighbourCells error: Currently not available for halfCellInt>1 \n");

    int x_begin = cellGrid.getInnerCellsBegin(0) - halfCellInt;
    int x_end = cellGrid.getInnerCellsEnd(0) - 1 + halfCellInt;
    int y_begin = cellGrid.getInnerCellsBegin(1) - halfCellInt;
    int y_end = cellGrid.getInnerCellsEnd(1) - 1 + halfCellInt;
    int z_begin = cellGrid.getInnerCellsBegin(2) - halfCellInt;
    int z_end = cellGrid.getInnerCellsEnd(2) - 1 + halfCellInt;
    bool ifreach;
    int incell_shift = (cell_shift % xgrid + xgrid) % xgrid;

    if (nodeGrid.getNodePosition(2) == 0 ||
        nodeGrid.getNodePosition(2) == getInt3DNodeGrid()[2] - 1)
    {
        // if (nodeGrid.getNodePosition(2)==0){
        // Dealing with the bottom layer of ghosts
        o = cellGrid.getInnerCellsBegin(2) - halfCellInt;
        for (int n = y_begin; n <= y_end; n++)
        {
            for (int m = x_begin; m < x_end; m++)
            {
                cell_idx1 = cellGrid.mapPositionToIndex(m, n, o);
                cell_idx2 = cellGrid.mapPositionToIndex(m + 1, n, o);
                cell1 = &cells[cell_idx1];
                cell2 = &cells[cell_idx2];
                if (getInt3DNodeGrid()[2] > 1 && m == x_begin)
                    for (ParticleList::Iterator it(cell1->particles); it.isValid(); ++it)
                        removeFromLocalParticles(&(*it));
                cell1->particles.clear();
                cell1->particles = cell2->particles;
                // for(Particle& p: cell2->particles) cell1->particles.push_back(p);
                if (getInt3DNodeGrid()[2] > 1) updateLocalParticles(cell1->particles);
            }
        }
        for (int n = y_begin; n <= y_end; n++)
        {
            cell_idx1 = cellGrid.mapPositionToIndex(x_begin + 1, n, o);
            cell_idx2 = cellGrid.mapPositionToIndex(x_end, n, o);
            cell1 = &cells[cell_idx1];
            cell2 = &cells[cell_idx2];
            // here it is cell2 copied from cell1
            if (getInt3DNodeGrid()[2] > 1)
                for (ParticleList::Iterator it(cell2->particles); it.isValid(); ++it)
                    removeFromLocalParticles(&(*it));
            cell2->particles.clear();
            if (getInt3DNodeGrid()[0] * getInt3DNodeGrid()[2] == 1)
            {
                cell2->particles = cell1->particles;
                // for(Particle& p: cell1->particles) cell2->particles.push_back(p);

                for (Particle& p : cell2->particles) p.position()[0] += Lx;
            }
            else if (getInt3DNodeGrid()[0] == 1 && getInt3DNodeGrid()[2] > 1)
                throw std::runtime_error("SORRY! I HAVE TO TK'ABOUT IT \n");
        }
        //}
        // onParticlesChanged();
        // rebuild commCell for Z+

        commCells[5].reals.clear();
        commCells[5].ghosts.clear();
        commCells[5].reals.reserve((x_end - x_begin + 1) * (y_end - y_begin + 1));
        commCells[5].ghosts.reserve((x_end - x_begin + 1) * (y_end - y_begin + 1));
        o_opp = z_end - 1;
        if (incell_shift == 0)
        {
            for (int m = x_begin; m <= x_end; m++)
            {
                for (int n = y_begin; n <= y_end; n++)
                {
                    cell_idx1 = cellGrid.mapPositionToIndex(m, n, o_opp);
                    cell1 = &cells[cell_idx1];
                    cell_idx2 = cellGrid.mapPositionToIndex(m, n, o);
                    cell2 = &cells[cell_idx2];
                    commCells[5].reals.push_back(cell1);
                    commCells[5].ghosts.push_back(cell2);
                }
            }
        }
        else
        {
            ifreach = false;
            for (int m = x_begin + 1; m <= x_end - 1; m++)
            {
                int m_shifted = (m - 1) % xgrid + 1 - incell_shift;
                if (m_shifted == x_begin && !ifreach) ifreach = true;
                if (ifreach)
                    for (int n = y_begin; n <= y_end; n++)
                    {
                        cell_idx1 = cellGrid.mapPositionToIndex(m, n, o_opp);
                        cell1 = &cells[cell_idx1];
                        cell_idx2 = cellGrid.mapPositionToIndex(m_shifted, n, o);
                        cell2 = &cells[cell_idx2];
                        commCells[5].reals.push_back(cell1);
                        commCells[5].ghosts.push_back(cell2);
                    }
            }
            ifreach = false;
            for (int m = x_begin + 1; m <= x_end - 1 && !ifreach; m++)
            {
                int m_shifted = (m - 1) % xgrid + 1 - (incell_shift - xgrid);
                if (!ifreach)
                    for (int n = y_begin; n <= y_end; n++)
                    {
                        cell_idx1 = cellGrid.mapPositionToIndex(m, n, o_opp);
                        cell1 = &cells[cell_idx1];
                        cell_idx2 = cellGrid.mapPositionToIndex(m_shifted, n, o);
                        cell2 = &cells[cell_idx2];
                        commCells[5].reals.push_back(cell1);
                        commCells[5].ghosts.push_back(cell2);
                    }
                if (m_shifted == x_end) ifreach = true;
            }
        }

        // if (nodeGrid.getNodePosition(2)==getInt3DNodeGrid()[2]-1){
        // Dealing with the top layer of ghosts
        o = cellGrid.getInnerCellsEnd(2) - 1 + halfCellInt;
        for (int n = y_end; n >= y_begin; n--)
        {
            for (int m = x_end; m > x_begin; m--)
            {
                cell_idx1 = cellGrid.mapPositionToIndex(m, n, o);
                cell_idx2 = cellGrid.mapPositionToIndex(m - 1, n, o);
                cell1 = &cells[cell_idx1];
                cell2 = &cells[cell_idx2];
                if (getInt3DNodeGrid()[2] > 1 && m == x_end)
                    for (ParticleList::Iterator it(cell1->particles); it.isValid(); ++it)
                        removeFromLocalParticles(&(*it));
                cell1->particles.clear();
                cell1->particles = cell2->particles;
                // for(Particle& p: cell2->particles) cell1->particles.push_back(p);
                if (getInt3DNodeGrid()[2] > 1) updateLocalParticles(cell1->particles);
            }
        }
        for (int n = y_begin; n <= y_end; n++)
        {
            cell_idx1 = cellGrid.mapPositionToIndex(x_begin, n, o);
            cell_idx2 = cellGrid.mapPositionToIndex(x_end - 1, n, o);
            cell1 = &cells[cell_idx1];
            cell2 = &cells[cell_idx2];
            if (getInt3DNodeGrid()[2] > 1)
                for (ParticleList::Iterator it(cell1->particles); it.isValid(); ++it)
                    removeFromLocalParticles(&(*it));
            cell1->particles.clear();
            if (getInt3DNodeGrid()[0] * getInt3DNodeGrid()[2] == 1)
            {
                cell1->particles = cell2->particles;
                // for(Particle& p: cell2->particles) cell1->particles.push_back(p);
                for (Particle& p : cell1->particles) p.position()[0] -= Lx;
            }
            else if (getInt3DNodeGrid()[0] == 1 && getInt3DNodeGrid()[2] > 1)
                throw std::runtime_error("SORRY! I HAVE TO TK'ABOUT IT \n");
        }
        //}

        // onParticlesChanged();
        // rebuild commCell for Z-
        commCells[4].reals.clear();
        commCells[4].ghosts.clear();
        commCells[4].reals.reserve((x_end - x_begin + 1) * (y_end - y_begin + 1));
        commCells[4].ghosts.reserve((x_end - x_begin + 1) * (y_end - y_begin + 1));
        o_opp = z_begin + 1;
        if (incell_shift == 0)
        {
            for (int m = x_begin; m <= x_end; m++)
            {
                for (int n = y_begin; n <= y_end; n++)
                {
                    cell_idx1 = cellGrid.mapPositionToIndex(m, n, o_opp);
                    cell1 = &cells[cell_idx1];
                    cell_idx2 = cellGrid.mapPositionToIndex(m, n, o);
                    cell2 = &cells[cell_idx2];
                    commCells[4].reals.push_back(cell1);
                    commCells[4].ghosts.push_back(cell2);
                }
            }
        }
        else
        {
            ifreach = false;
            for (int m = x_begin + 1; m <= x_end - 1 && !ifreach; m++)
            {
                int m_shifted = (m - 1) % xgrid + 1 + incell_shift;
                if (!ifreach)
                    for (int n = y_begin; n <= y_end; n++)
                    {
                        cell_idx1 = cellGrid.mapPositionToIndex(m, n, o_opp);
                        cell1 = &cells[cell_idx1];
                        cell_idx2 = cellGrid.mapPositionToIndex(m_shifted, n, o);
                        cell2 = &cells[cell_idx2];
                        commCells[4].reals.push_back(cell1);
                        commCells[4].ghosts.push_back(cell2);
                    }
                if (m_shifted == x_end) ifreach = true;
            }
            ifreach = false;
            for (int m = x_begin + 1; m <= x_end - 1; m++)
            {
                int m_shifted = (m - 1) % xgrid + 1 + (incell_shift - xgrid);
                if (m_shifted == x_begin && !ifreach) ifreach = true;
                if (ifreach)
                    for (int n = y_begin; n <= y_end; n++)
                    {
                        cell_idx1 = cellGrid.mapPositionToIndex(m, n, o_opp);
                        cell1 = &cells[cell_idx1];
                        cell_idx2 = cellGrid.mapPositionToIndex(m_shifted, n, o);
                        cell2 = &cells[cell_idx2];
                        commCells[4].reals.push_back(cell1);
                        commCells[4].ghosts.push_back(cell2);
                    }
            }
        }
    }
    LOG4ESPP_DEBUG(logger, "done");
}

Cell* DomainDecomposition::mapPositionToCell(const Real3D& pos)
{
    return &cells[cellGrid.mapPositionToCell(pos)];
}

Cell* DomainDecomposition::mapPositionToCellClipped(const Real3D& pos)
{
    return &cells[cellGrid.mapPositionToCellClipped(pos)];
}

Cell* DomainDecomposition::mapPositionToCellChecked(const Real3D& pos)
{
    longint c = cellGrid.mapPositionToCellChecked(pos);
    if (c == CellGrid::noCell)
    {
        return 0;
    }
    else
    {
        return &cells[c];
    }
}

longint DomainDecomposition::mapPositionToNodeClipped(const Real3D& pos)
{
    return nodeGrid.mapPositionToNodeClipped(pos);
}

bool DomainDecomposition::checkIsRealParticle(longint id, const Real3D& pos)
{
    return getSystem()->comm->rank() == mapPositionToNodeClipped(pos);
}

bool DomainDecomposition::appendParticles(ParticleList& l, int dir)
{
    bool outlier = false;

    LOG4ESPP_DEBUG(logger, "got " << l.size() << " particles");

    for (ParticleList::iterator it = l.begin(), end = l.end(); it != end; ++it)
    {
        Real3D& pos = it->position();

        if (nodeGrid.getBoundary(dir) != 0)
        {
            getSystem()->bc->foldCoordinate(pos, it->image(), nodeGrid.convertDirToCoord(dir));
            LOG4ESPP_TRACE(logger, "folded coordinate " << nodeGrid.convertDirToCoord(dir)
                                                        << " of particle " << it->id());
        }

        longint cell;
        if (cellGrid.mapPositionToCellCheckedAndClipped(cell, pos))
        {
            LOG4ESPP_TRACE(logger,
                           "particle " << it->id() << " @ " << pos << " is not inside node domain");
            outlier = true;
        }

        LOG4ESPP_TRACE(logger, "append part " << it->id() << " to cell " << cell);

        appendIndexedParticle(cells[cell].particles, *it);
    }
    return outlier;
}

void DomainDecomposition::decomposeRealParticles()
{
    // std::cout << getSystem()->comm->rank() << ": " << " decomposeRealParticles\n";

    LOG4ESPP_DEBUG(logger, "starting, expected comm buffer size " << exchangeBufferSize);

    // allocate send/recv buffers. We use the size as we need maximally so far, to avoid
    // reallocation
    // TODO: This might be a problem when all particles are created on a single node initially!
    ParticleList sendBufL;
    sendBufL.reserve(exchangeBufferSize);
    ParticleList sendBufR;
    sendBufR.reserve(exchangeBufferSize);
    ParticleList recvBufL;
    recvBufL.reserve(exchangeBufferSize);
    ParticleList recvBufR;
    recvBufR.reserve(exchangeBufferSize);

    ParticleList sendBuf2;
    ParticleList recvBuf2;

    bool allFinished;
    real offs = getSystem()->shearOffset;
    real Lz = getSystem()->bc->getBoxL()[2];
    real Lx = getSystem()->bc->getBoxL()[0];

    int node1, node2, ptmp1, ptmp2;
    int allCellGrid = getInt3DNodeGrid()[0] * getInt3DCellGrid()[0];
    real mid_1, mid_2;

    do
    {
        bool finished = true;

        for (int coord = 0; coord < 3; ++coord)
        {
            LOG4ESPP_DEBUG(logger, "starting with direction " << coord);

            if (nodeGrid.getGridSize(coord) > 1)
            {
                if (offs > .0 && coord == 2 &&
                    nodeGrid.getNodePosition(2) % (getInt3DNodeGrid()[2] - 1) == 0)
                {
                    sendBuf2.reserve(exchangeBufferSize);
                    recvBuf2.reserve(exchangeBufferSize);

                    if (nodeGrid.getNodePosition(2) == 0)
                    {
                        ptmp1 = (static_cast<int>(floor(((nodeGrid.getNodePosition(0) + .0) /
                                                             (getInt3DNodeGrid()[0] + .0) * Lx +
                                                         offs) /
                                                        Lx * getInt3DNodeGrid()[0])) %
                                     getInt3DNodeGrid()[0] +
                                 getInt3DNodeGrid()[0]) %
                                getInt3DNodeGrid()[0];
                        node1 = (getInt3DNodeGrid()[2] - nodeGrid.getNodePosition(2) - 1) *
                                    getInt3DNodeGrid()[0] * getInt3DNodeGrid()[1] +
                                nodeGrid.getNodePosition(1) * getInt3DNodeGrid()[0] + ptmp1;
                        ptmp2 = (ptmp1 + 1) % getInt3DNodeGrid()[0];
                        node2 = node1 + ptmp2 - ptmp1;
                        mid_1 = Lx * (ptmp1 + 0.5) / (getInt3DNodeGrid()[0] + .0);
                        mid_2 = Lx * (ptmp2 + 0.5) / (getInt3DNodeGrid()[0] + .0);
                    }
                    else if (nodeGrid.getNodePosition(2) == getInt3DNodeGrid()[2] - 1)
                    {
                        ptmp1 = (static_cast<int>(floor(((nodeGrid.getNodePosition(0) + 1.0) /
                                                             (getInt3DNodeGrid()[0] + .0) * Lx -
                                                         offs) /
                                                        Lx * getInt3DNodeGrid()[0])) %
                                     getInt3DNodeGrid()[0] +
                                 getInt3DNodeGrid()[0]) %
                                getInt3DNodeGrid()[0];
                        node1 = (getInt3DNodeGrid()[2] - nodeGrid.getNodePosition(2) - 1) *
                                    getInt3DNodeGrid()[0] * getInt3DNodeGrid()[1] +
                                nodeGrid.getNodePosition(1) * getInt3DNodeGrid()[0] + ptmp1;
                        ptmp2 = (ptmp1 - 1 + getInt3DNodeGrid()[0]) % getInt3DNodeGrid()[0];
                        node2 = node1 + ptmp2 - ptmp1;
                        mid_1 = Lx * (ptmp1 + 0.5) / (getInt3DNodeGrid()[0] + .0);
                        mid_2 = Lx * (ptmp2 + 0.5) / (getInt3DNodeGrid()[0] + .0);
                    }

                    for (std::vector<Cell*>::iterator it = realCells.begin(), end = realCells.end();
                         it != end; ++it)
                    {
                        Cell& cell = **it;

                        // do not use an iterator here, since we need to take out particles during
                        // the loop
                        for (size_t p = 0; p < cell.particles.size(); ++p)
                        {
                            Particle& part = cell.particles[p];
                            const Real3D& pos = part.position();
                            // check whether the particle is now "left" of the local domain

                            if (pos[coord] - cellGrid.getMyLeft(coord) < -ROUND_ERROR_PREC)
                            {
                                LOG4ESPP_TRACE(logger, "send particle left " << part.id());

                                if (nodeGrid.getNodePosition(2) == 0)
                                {
                                    real xtmp = part.position()[0] + offs;
                                    int itmp = static_cast<int>(floor(xtmp / Lx));
                                    part.position()[0] = xtmp - (itmp + .0) * Lx;

                                    if (abs(part.position()[0] - mid_1) <
                                        abs(part.position()[0] - mid_2))
                                        moveIndexedParticle(sendBufL, cell.particles, p);
                                    else
                                        moveIndexedParticle(sendBuf2, cell.particles, p);
                                }
                                else
                                    moveIndexedParticle(sendBufL, cell.particles, p);

                                // redo same particle since we took one out here, so it's a new one
                                --p;
                            }
                            // check whether the particle is now "right" of the local domain
                            else if (pos[coord] - cellGrid.getMyRight(coord) >= ROUND_ERROR_PREC)
                            {
                                LOG4ESPP_TRACE(logger, "send particle right " << part.id());

                                if (nodeGrid.getNodePosition(2) == getInt3DNodeGrid()[2] - 1)
                                {
                                    real xtmp = part.position()[0] - offs;
                                    int itmp = static_cast<int>(floor(xtmp / Lx));
                                    part.position()[0] = xtmp - (itmp + .0) * Lx;

                                    if (abs(part.position()[0] - mid_1) <
                                        abs(part.position()[0] - mid_2))
                                        moveIndexedParticle(sendBufR, cell.particles, p);
                                    else
                                        moveIndexedParticle(sendBuf2, cell.particles, p);
                                }
                                else
                                    moveIndexedParticle(sendBufR, cell.particles, p);

                                --p;
                            }
                            // Sort particles in cells of this node during last direction
                            else if (coord == 2)
                            {
                                const Real3D& pos = part.position();
                                Cell* sortCell = mapPositionToCellChecked(pos);
                                if (sortCell != &cell)
                                {
                                    if (sortCell == 0)
                                    {
                                        // particle is not in the local domain
                                        LOG4ESPP_DEBUG(logger, "take another loop: particle "
                                                                   << part.id() << " @ " << pos
                                                                   << " is not inside node domain "
                                                                      "after neighbor exchange");
                                        // isnan function is C99 only, x != x is only true if x ==
                                        // nan
                                        if (pos[0] != pos[0] || pos[1] != pos[1] ||
                                            pos[2] != pos[2])
                                        {
                                            // TODO: error handling
                                            LOG4ESPP_ERROR(logger,
                                                           "particle "
                                                               << part.id()
                                                               << " has moved to outer space (one "
                                                                  "or more coordinates are nan)");
                                        }
                                        else
                                        {
                                            // particle stays where it is, and will be sorted in the
                                            // next round
                                            finished = false;
                                        }
                                    }
                                    else
                                    {
                                        // particle is in the local domain
                                        moveIndexedParticle(sortCell->particles, cell.particles, p);
                                        --p;
                                    }
                                }
                            }
                        }
                    }

                    // Exchange particles, odd-even rule
                    if (nodeGrid.getNodePosition(2) == 0)
                    {
                        if (nodeGrid.getNodePosition(coord) % 2 == 0)
                        {
                            sendParticles(sendBufL, node1);
                            recvParticles(recvBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
                            sendParticles(sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
                            recvParticles(recvBufL, node1);
                            sendParticles(sendBuf2, node2);
                            recvParticles(recvBuf2, node2);
                        }
                        else
                        {
                            recvParticles(recvBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
                            sendParticles(sendBufL, node1);
                            recvParticles(recvBufL, node1);
                            sendParticles(sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
                            sendParticles(sendBuf2, node2);
                            recvParticles(recvBuf2, node2);
                        }
                    }
                    else if (nodeGrid.getNodePosition(2) == getInt3DNodeGrid()[2] - 1)
                    {
                        if (nodeGrid.getNodePosition(coord) % 2 == 0)
                        {
                            sendParticles(sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
                            recvParticles(recvBufR, node1);
                            sendParticles(sendBufR, node1);
                            recvParticles(recvBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
                            sendParticles(sendBuf2, node2);
                            recvParticles(recvBuf2, node2);
                        }
                        else
                        {
                            recvParticles(recvBufR, node1);
                            sendParticles(sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
                            recvParticles(recvBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
                            sendParticles(sendBufR, node1);
                            recvParticles(recvBuf2, node2);
                            sendParticles(sendBuf2, node2);
                        }
                    }

                    // sort received particles to cells
                    if (appendParticles(recvBufL, 2 * coord) && coord == 2) finished = false;
                    if (appendParticles(recvBufR, 2 * coord + 1) && coord == 2) finished = false;

                    if (nodeGrid.getNodePosition(2) == 0)
                    {
                        if (appendParticles(recvBuf2, 2 * coord) && coord == 2) finished = false;
                    }
                    else if (nodeGrid.getNodePosition(2) == getInt3DNodeGrid()[2] - 1)
                    {
                        if (appendParticles(recvBuf2, 2 * coord + 1) && coord == 2)
                            finished = false;
                    }

                    // reset send/recv buffers
                    sendBufL.resize(0);
                    sendBufR.resize(0);
                    recvBufL.resize(0);
                    recvBufR.resize(0);
                    sendBuf2.resize(0);
                    recvBuf2.resize(0);
                }
                else
                {
                    for (std::vector<Cell*>::iterator it = realCells.begin(), end = realCells.end();
                         it != end; ++it)
                    {
                        Cell& cell = **it;

                        // do not use an iterator here, since we need to take out particles during
                        // the loop
                        for (size_t p = 0; p < cell.particles.size(); ++p)
                        {
                            Particle& part = cell.particles[p];
                            const Real3D& pos = part.position();

                            // check whether the particle is now "left" of the local domain
                            if (pos[coord] - cellGrid.getMyLeft(coord) < -ROUND_ERROR_PREC)
                            {
                                LOG4ESPP_TRACE(logger, "send particle left " << part.id());
                                moveIndexedParticle(sendBufL, cell.particles, p);
                                // redo same particle since we took one out here, so it's a new one
                                --p;
                            }
                            // check whether the particle is now "right" of the local domain
                            else if (pos[coord] - cellGrid.getMyRight(coord) >= ROUND_ERROR_PREC)
                            {
                                LOG4ESPP_TRACE(logger, "send particle right " << part.id());
                                moveIndexedParticle(sendBufR, cell.particles, p);
                                --p;
                            }
                            // Sort particles in cells of this node during last direction
                            else if (coord == 2)
                            {
                                const Real3D& pos = part.position();
                                Cell* sortCell = mapPositionToCellChecked(pos);
                                if (sortCell != &cell)
                                {
                                    if (sortCell == 0)
                                    {
                                        // particle is not in the local domain
                                        LOG4ESPP_DEBUG(logger, "take another loop: particle "
                                                                   << part.id() << " @ " << pos
                                                                   << " is not inside node domain "
                                                                      "after neighbor exchange");
                                        // isnan function is C99 only, x != x is only true if x ==
                                        // nan
                                        if (pos[0] != pos[0] || pos[1] != pos[1] ||
                                            pos[2] != pos[2])
                                        {
                                            // TODO: error handling
                                            LOG4ESPP_ERROR(logger,
                                                           "particle "
                                                               << part.id()
                                                               << " has moved to outer sPace (one "
                                                                  "or more coordinates are nan)");
                                        }
                                        else
                                        {
                                            // particle stays where it is, and will be sorted in the
                                            // next round
                                            finished = false;
                                        }
                                    }
                                    else
                                    {
                                        // particle is in the local domain
                                        moveIndexedParticle(sortCell->particles, cell.particles, p);
                                        --p;
                                    }
                                }
                            }
                        }
                    }
                    // Exchange particles, odd-even rule
                    if (nodeGrid.getNodePosition(coord) % 2 == 0)
                    {
                        sendParticles(sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
                        recvParticles(recvBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
                        sendParticles(sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
                        recvParticles(recvBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
                    }
                    else
                    {
                        recvParticles(recvBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
                        sendParticles(sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
                        recvParticles(recvBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
                        sendParticles(sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
                    }

                    // sort received particles to cells
                    if (appendParticles(recvBufL, 2 * coord) && coord == 2) finished = false;
                    if (appendParticles(recvBufR, 2 * coord + 1) && coord == 2) finished = false;

                    // reset send/recv buffers
                    sendBufL.resize(0);
                    sendBufR.resize(0);
                    recvBufL.resize(0);
                    recvBufR.resize(0);
                }
            }
            else
            {
                /* Single node direction case (no communication)
                    Fold particles that have left the box */

                if (offs > .0)
                {
                    for (std::vector<Cell*>::iterator it = realCells.begin(), end = realCells.end();
                         it != end; ++it)
                    {
                        Cell& cell = **it;
                        // do not use an iterator here, since we have need to take out particles
                        // during the loop
                        for (size_t p = 0; p < cell.particles.size(); ++p)
                        {
                            Particle& part = cell.particles[p];
                            real ztmp = part.position()[2];
                            getSystem()->bc->foldCoordinate(part.position(), part.image(), coord);
                            if (coord == 2 && ztmp != part.position()[2])
                            {
                                real xtmp =
                                    part.position()[0] + offs * (part.position()[2] - ztmp) / Lz;
                                int itmp = static_cast<int>(floor(xtmp / Lx));
                                part.position()[0] = xtmp - (itmp + .0) * Lx;
                            }
                            LOG4ESPP_TRACE(logger, "folded coordinate " << coord << " of particle "
                                                                        << part.id());

                            if (coord == 2)
                            {
                                Cell* sortCell = mapPositionToCellChecked(part.position());

                                if (sortCell != &cell)
                                {
                                    if (sortCell == 0)
                                    {
                                        LOG4ESPP_DEBUG(logger, "take another loop: particle "
                                                                   << part.id() << " @ "
                                                                   << part.position()
                                                                   << " is not inside node domain "
                                                                      "after neighbor exchange");
                                        const Real3D& pos = part.position();
                                        // isnan function is C99 only, x != x is only true if x ==
                                        // nan
                                        if (pos[0] != pos[0] || pos[1] != pos[1] ||
                                            pos[2] != pos[2])
                                        {
                                            LOG4ESPP_ERROR(logger,
                                                           "particle "
                                                               << part.id()
                                                               << " has moved to outer spAce (one "
                                                                  "or more coordinates are nan)");
                                        }
                                        else
                                        {
                                            // particle stays where it is, and will be sorted in the
                                            // next round
                                            finished = false;
                                        }
                                    }
                                    else
                                    {
                                        moveIndexedParticle(sortCell->particles, cell.particles, p);
                                        --p;
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (std::vector<Cell*>::iterator it = realCells.begin(), end = realCells.end();
                         it != end; ++it)
                    {
                        Cell& cell = **it;
                        // do not use an iterator here, since we have need to take out particles
                        // during the loop
                        for (size_t p = 0; p < cell.particles.size(); ++p)
                        {
                            Particle& part = cell.particles[p];
                            getSystem()->bc->foldCoordinate(part.position(), part.image(), coord);
                            LOG4ESPP_TRACE(logger, "folded coordinate " << coord << " of particle "
                                                                        << part.id());

                            if (coord == 2)
                            {
                                Cell* sortCell = mapPositionToCellChecked(part.position());

                                if (sortCell != &cell)
                                {
                                    if (sortCell == 0)
                                    {
                                        LOG4ESPP_DEBUG(logger, "take another loop: particle "
                                                                   << part.id() << " @ "
                                                                   << part.position()
                                                                   << " is not inside node domain "
                                                                      "after neighbor exchange");
                                        const Real3D& pos = part.position();
                                        // isnan function is C99 only, x != x is only true if x ==
                                        // nan
                                        if (pos[0] != pos[0] || pos[1] != pos[1] ||
                                            pos[2] != pos[2])
                                        {
                                            LOG4ESPP_ERROR(logger,
                                                           "particle "
                                                               << part.id()
                                                               << " has moved to outer spaCe (one "
                                                                  "or more coordinates are nan)");
                                        }
                                        else
                                        {
                                            // particle stays where it is, and will be sorted in the
                                            // next round
                                            finished = false;
                                        }
                                    }
                                    else
                                    {
                                        moveIndexedParticle(sortCell->particles, cell.particles, p);
                                        --p;
                                    }
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

    exchangeBufferSize = std::max(
        exchangeBufferSize,
        std::max(sendBufL.capacity(),
                 std::max(sendBufR.capacity(),
                          std::max(recvBufL.capacity(),
                                   std::max(recvBufR.capacity(),
                                            std::max(sendBuf2.capacity(), recvBuf2.capacity()))))));

    LOG4ESPP_DEBUG(
        logger, "finished exchanging particles, new send/recv buffer size " << exchangeBufferSize);

    LOG4ESPP_DEBUG(logger, "done");
}

void DomainDecomposition::exchangeGhosts()
{
    LOG4ESPP_DEBUG(logger, "exchangeGhosts -> ghost communication sizes first, real->ghost");
    doGhostCommunication(true, true, dataOfExchangeGhosts);
}

void DomainDecomposition::updateGhosts()
{
    LOG4ESPP_DEBUG(logger, "updateGhosts -> ghost communication no sizes, real->ghost");
    doGhostCommunication(false, true, dataOfUpdateGhosts);
}

void DomainDecomposition::updateGhostsV()
{
    LOG4ESPP_DEBUG(logger, "updateGhostsV -> ghost communication no sizes, real->ghost velocities");
    doGhostCommunication(false, true, 2);  // 2 is the bitflag for particle momentum
}

void DomainDecomposition::collectGhostForces()
{
    LOG4ESPP_DEBUG(logger, "collectGhosts -> ghost communication no sizes, ghost->real");
    doGhostCommunication(false, false);
}

void DomainDecomposition::fillCells(std::vector<Cell*>& cv,
                                    const int leftBoundary[3],
                                    const int rightBoundary[3])
{
    LOG4ESPP_DEBUG(logger, "filling: " << leftBoundary[0] << "-" << (rightBoundary[0] - 1) << " "
                                       << leftBoundary[1] << "-" << (rightBoundary[1] - 1) << " "
                                       << leftBoundary[2] << "-" << (rightBoundary[2] - 1));

    longint total = 1;
    for (int i = 0; i < 3; ++i)
    {
        if (leftBoundary[i] < 0 || leftBoundary[i] > cellGrid.getFrameGridSize(i) ||
            rightBoundary[i] < 0 || rightBoundary[i] > cellGrid.getFrameGridSize(i) ||
            leftBoundary[i] >= rightBoundary[i])
        {
            throw std::runtime_error(
                "DomainDecomposition::fillCells: wrong cell grid specified internally");
        }
        total *= (rightBoundary[i] - leftBoundary[i]);
    }
    cv.reserve(total);

    for (int o = leftBoundary[0]; o < rightBoundary[0]; ++o)
    {
        for (int n = leftBoundary[1]; n < rightBoundary[1]; ++n)
        {
            for (int m = leftBoundary[2]; m < rightBoundary[2]; ++m)
            {
                int i = cellGrid.mapPositionToIndex(o, n, m);
                LOG4ESPP_TRACE(logger, "add cell " << i);
                cv.push_back(&cells[i]);
            }
        }
    }

    LOG4ESPP_DEBUG(logger, "expected " << total << " cells, filled with " << cv.size());
}

void DomainDecomposition::prepareGhostCommunication()
{
    // direction loop: x, y, z
    for (int coord = 0; coord < 3; ++coord)
    {
        // boundaries of area to send
        int leftBoundary[3], rightBoundary[3];
        /* boundaries perpendicular directions are the same for left/right send.
        We also send the ghost frame that we have already, so the data amount
        increase with each cycle.

        For a direction that was done already, i.e. is smaller than dir,
        we take the full ghost frame, otherwise only the inner frame.  */
        for (int offset = 1; offset <= 2; ++offset)
        {
            int otherCoord = (coord + offset) % 3;
            if (otherCoord < coord)
            {
                leftBoundary[otherCoord] = 0;
                rightBoundary[otherCoord] = cellGrid.getFrameGridSize(otherCoord);
            }
            else
            {
                leftBoundary[otherCoord] = cellGrid.getInnerCellsBegin(otherCoord);
                rightBoundary[otherCoord] = cellGrid.getInnerCellsEnd(otherCoord);
            }
        }

        //  lr loop: left right - loop
        for (int lr = 0; lr < 2; ++lr)
        {
            int dir = 2 * coord + lr;

            /* participating real particles from this node */
            LOG4ESPP_DEBUG(logger, "direction " << dir << " reals");

            if (lr == 0)
            {
                leftBoundary[coord] = cellGrid.getInnerCellsBegin(coord);
                rightBoundary[coord] =
                    cellGrid.getInnerCellsBegin(coord) + cellGrid.getFrameWidth();
            }
            else
            {
                leftBoundary[coord] = cellGrid.getInnerCellsEnd(coord) - cellGrid.getFrameWidth();
                rightBoundary[coord] = cellGrid.getInnerCellsEnd(coord);
            }
            fillCells(commCells[dir].reals, leftBoundary, rightBoundary);

            /* participating ghosts from this node */
            LOG4ESPP_DEBUG(logger, "direction " << dir << " ghosts");

            if (lr == 0)
            {
                leftBoundary[coord] = cellGrid.getInnerCellsEnd(coord);
                rightBoundary[coord] = cellGrid.getInnerCellsEnd(coord) + cellGrid.getFrameWidth();
            }
            else
            {
                leftBoundary[coord] = cellGrid.getInnerCellsBegin(coord) - cellGrid.getFrameWidth();
                rightBoundary[coord] = cellGrid.getInnerCellsBegin(coord);
            }
            fillCells(commCells[dir].ghosts, leftBoundary, rightBoundary);
        }
    }
}

void DomainDecomposition::doGhostCommunication(bool sizesFirst, bool realToGhosts, int extradata)
{
    LOG4ESPP_DEBUG(logger, "do ghost communication "
                               << (sizesFirst ? "with sizes " : "")
                               << (realToGhosts ? "reals to ghosts " : "ghosts to reals ")
                               << extradata);

    /* direction loop: x, y, z.
 Here we could in principle build in a one sided ghost
 communication, simply by taking the lr loop only over one
 value. */

    real offs = getSystem()->shearOffset;

    for (int _coord = 0; _coord < 3; ++_coord)
    {
        /* inverted processing order for ghost force communication,
            since the corner ghosts have to be collected via several
            nodes. We now add back the corner ghost forces first again
            to ghost forces, which only eventually go back to the real
            particle.
        */
        int coord = realToGhosts ? _coord : (2 - _coord);
        real curCoordBoxL = getSystem()->bc->getBoxL()[coord];

        // lr loop: left right
        for (int lr = 0; lr < 2; ++lr)
        {
            int dir = 2 * coord + lr;
            int oppositeDir = 2 * coord + (1 - lr);

            Real3D shift(0, 0, 0);

            shift[coord] = nodeGrid.getBoundary(dir) * curCoordBoxL;

            LOG4ESPP_DEBUG(logger, "direction " << dir);

            if (nodeGrid.getGridSize(coord) == 1)
            {
                LOG4ESPP_DEBUG(logger, "local communication");

                // copy operation, we have to receive as many cells as we send
                if (commCells[dir].ghosts.size() != commCells[dir].reals.size())
                {
                    throw std::runtime_error(
                        "DomainDecomposition::doGhostCommunication: send/recv cell structure "
                        "mismatch during local copy");
                }

                if (offs > .0 && coord == 2)
                {
                    for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i)
                    {
                        if (realToGhosts)
                        {
                            copyRealsToGhosts_LEBC(*commCells[dir].reals[i],
                                                   *commCells[dir].ghosts[i], extradata, shift);
                        }
                        else
                        {
                            addGhostForcesToReals(*commCells[dir].ghosts[i],
                                                  *commCells[dir].reals[i]);
                        }
                    }
                }
                else
                {
                    for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i)
                    {
                        if (realToGhosts)
                        {
                            copyRealsToGhosts(*commCells[dir].reals[i], *commCells[dir].ghosts[i],
                                              extradata, shift);
                        }
                        else
                        {
                            addGhostForcesToReals(*commCells[dir].ghosts[i],
                                                  *commCells[dir].reals[i]);
                        }
                    }
                }
            }
            else
            {
                // exchange size information, if necessary

                // This is used for cell communication in parallel computing with LEBC
                int allCellGrid = getInt3DNodeGrid()[0] * getInt3DCellGrid()[0];
                if (offs > .0 && coord == 2 &&
                    (nodeGrid.getNodePosition(2) == 0 ||
                     nodeGrid.getNodePosition(2) == getInt3DNodeGrid()[2] - 1)
                    //&& !(sizesFirst && realToGhosts)
                )
                {
                    if (commCells_bkp[dir - 4].reals.size() == 0)
                    {
                        commCells_bkp[dir - 4].reals.reserve(commCells[dir].reals.size());
                        commCells_bkp[dir - 4].ghosts.reserve(commCells[dir].ghosts.size());
                        for (int i = 0; i < commCells[dir].reals.size(); i++)
                        {
                            commCells_bkp[dir - 4].reals.push_back(commCells[dir].reals[i]);
                            commCells_bkp[dir - 4].ghosts.push_back(commCells[dir].ghosts[i]);
                        }
                    }

                    int incell_shift = ((getSystem()->ghostShift) % getInt3DCellGrid()[0] +
                                        getInt3DCellGrid()[0]) %
                                       getInt3DCellGrid()[0];
                    int new_dir = (nodeGrid.getNodePosition(2) > 0 ? -3 : -6) + dir;
                    int sz01 =
                        (getInt3DCellGrid()[1] + 2) * (getInt3DCellGrid()[0] + 1 - incell_shift);
                    int sz02 = (getInt3DCellGrid()[1] + 2) * (incell_shift + 1);
                    real Lx = getSystem()->bc->getBoxL()[0];

                    // define node1 & node2 for sendBuffer or recvBuffer
                    int node1, node2, ptmp1, ptmp2, shift_tmp;
                    if ((new_dir > 0 ? new_dir : -new_dir) == 2)
                    {
                        if (dir == 4)
                        {
                            shift_tmp = getSystem()->ghostShift;
                            node1 = (getInt3DNodeGrid()[2] - 1) * getInt3DNodeGrid()[0] *
                                    getInt3DNodeGrid()[1];
                        }
                        else
                        {
                            shift_tmp = -getSystem()->ghostShift;
                            node1 = 0;
                        }
                    }
                    else if ((new_dir > 0 ? new_dir : -new_dir) == 1)
                    {
                        if (oppositeDir == 4)
                        {
                            shift_tmp = getSystem()->ghostShift;
                            node1 = (getInt3DNodeGrid()[2] - 1) * getInt3DNodeGrid()[0] *
                                    getInt3DNodeGrid()[1];
                        }
                        else
                        {
                            shift_tmp = -getSystem()->ghostShift;
                            node1 = 0;
                        }
                    }
                    else
                        throw std::runtime_error(
                            "doGhostCommunication error: new_dir should be equal to one of "
                            "[-2,-1,1,2] \n");
                    if (shift_tmp >= 0)
                    {
                        ptmp1 = ((nodeGrid.getNodePosition(0) + shift_tmp / getInt3DCellGrid()[0]) %
                                     getInt3DNodeGrid()[0] +
                                 getInt3DNodeGrid()[0]) %
                                getInt3DNodeGrid()[0];
                        node1 += nodeGrid.getNodePosition(1) * getInt3DNodeGrid()[0] + ptmp1;
                        ptmp2 = (ptmp1 + 1) % getInt3DNodeGrid()[0];
                        node2 = node1 + ptmp2 - ptmp1;
                    }
                    else
                    {
                        ptmp1 = ((nodeGrid.getNodePosition(0) + shift_tmp / getInt3DCellGrid()[0]) %
                                     getInt3DNodeGrid()[0] +
                                 getInt3DNodeGrid()[0]) %
                                getInt3DNodeGrid()[0];
                        node1 += nodeGrid.getNodePosition(1) * getInt3DNodeGrid()[0] + ptmp1;
                        ptmp2 = ((ptmp1 - 1) % getInt3DNodeGrid()[0] + getInt3DNodeGrid()[0]) %
                                getInt3DNodeGrid()[0];
                        node2 = node1 + ptmp2 - ptmp1;
                    }

                    int rsize[3], gsize[3], cnode[2];
                    if (incell_shift == 0)
                    {
                        gsize[0] = 0;
                        rsize[0] = 0;
                        cnode[0] = node1;
                        gsize[1] = commCells[dir].ghosts.size();
                        rsize[1] = commCells[dir].reals.size();
                        gsize[2] = commCells[dir].ghosts.size();
                        rsize[2] = commCells[dir].reals.size();
                    }
                    else
                    {
                        gsize[0] = 0;
                        rsize[0] = 0;
                        cnode[0] = node1;
                        cnode[1] = node2;
                        gsize[1] = sz01;
                        rsize[1] = sz01;
                        gsize[2] = commCells[dir].ghosts.size();
                        rsize[2] = commCells[dir].reals.size();
                    }

                    for (int k = 0; k == 0 || (k == 1 && incell_shift > 0); k++)
                    {
                        if (sizesFirst)
                        {
                            LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes");
                            std::vector<longint> sendSizes, recvSizes;

                            if ((new_dir > 0 ? new_dir : -new_dir) == 2)
                            {
                                // prepare buffers
                                sendSizes.reserve(rsize[k + 1] - rsize[k]);
                                for (int i = rsize[k], end = rsize[k + 1]; i < end; ++i)
                                {
                                    sendSizes.push_back(commCells[dir].reals[i]->particles.size());
                                }
                                if (k == 0) recvSizes.resize(commCells_bkp[dir - 4].ghosts.size());
                            }
                            else if ((new_dir > 0 ? new_dir : -new_dir) == 1)
                            {
                                // prepare buffers
                                if (k == 0)
                                {
                                    sendSizes.reserve(commCells_bkp[dir - 4].reals.size());
                                    for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i)
                                    {
                                        sendSizes.push_back(
                                            commCells_bkp[dir - 4].reals[i]->particles.size());
                                    }
                                }
                                recvSizes.resize(gsize[k + 1] - gsize[k]);
                            }
                            else
                                throw std::runtime_error(
                                    "doGhostCommunication error: new_dir should be equal to one of "
                                    "[-2,-1,1,2] \n");

                            if ((new_dir > 0 ? new_dir : -new_dir) == 2)
                            {
                                if (nodeGrid.getNodePosition(coord) % 2 == 0)
                                {
                                    LOG4ESPP_DEBUG(logger, "sending to node " << cnode[k]);
                                    getSystem()->comm->send(cnode[k], DD_COMM_TAG, &(sendSizes[0]),
                                                            sendSizes.size());
                                    if (k == 0)
                                    {
                                        LOG4ESPP_DEBUG(logger, "receiving from node "
                                                                   << nodeGrid.getNodeNeighborIndex(
                                                                          oppositeDir));
                                        getSystem()->comm->recv(
                                            nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG,
                                            &(recvSizes[0]), recvSizes.size());
                                    }
                                }
                                else
                                {
                                    if (k == 0)
                                    {
                                        LOG4ESPP_DEBUG(logger, "receiving from node "
                                                                   << nodeGrid.getNodeNeighborIndex(
                                                                          oppositeDir));
                                        getSystem()->comm->recv(
                                            nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG,
                                            &(recvSizes[0]), recvSizes.size());
                                    }
                                    LOG4ESPP_DEBUG(logger, "sending to node " << cnode[k]);
                                    getSystem()->comm->send(cnode[k], DD_COMM_TAG, &(sendSizes[0]),
                                                            sendSizes.size());
                                }
                            }
                            else if ((new_dir > 0 ? new_dir : -new_dir) == 1)
                            {
                                // exchange sizes, odd-even rule
                                if (nodeGrid.getNodePosition(coord) % 2 == 0)
                                {
                                    LOG4ESPP_DEBUG(logger, "receiving from node " << cnode[k]);
                                    if (k == 0)
                                    {
                                        LOG4ESPP_DEBUG(logger,
                                                       "sending to node "
                                                           << nodeGrid.getNodeNeighborIndex(dir));
                                        getSystem()->comm->send(nodeGrid.getNodeNeighborIndex(dir),
                                                                DD_COMM_TAG, &(sendSizes[0]),
                                                                sendSizes.size());
                                    }
                                    getSystem()->comm->recv(cnode[k], DD_COMM_TAG, &(recvSizes[0]),
                                                            recvSizes.size());
                                }
                                else
                                {
                                    LOG4ESPP_DEBUG(logger, "receiving from node " << cnode[k]);
                                    getSystem()->comm->recv(cnode[k], DD_COMM_TAG, &(recvSizes[0]),
                                                            recvSizes.size());
                                    if (k == 0)
                                    {
                                        LOG4ESPP_DEBUG(logger,
                                                       "sending to node "
                                                           << nodeGrid.getNodeNeighborIndex(dir));
                                        getSystem()->comm->send(nodeGrid.getNodeNeighborIndex(dir),
                                                                DD_COMM_TAG, &(sendSizes[0]),
                                                                sendSizes.size());
                                    }
                                }
                            }
                            ////MPI wait
                            // getSystem()->comm->barrier();

                            // resize according to received information
                            if ((new_dir > 0 ? new_dir : -new_dir) == 2)
                            {
                                if (k == 0)
                                    for (int i = 0, end = commCells_bkp[dir - 4].ghosts.size();
                                         i < end; ++i)
                                    {
                                        commCells_bkp[dir - 4].ghosts[i]->particles.resize(
                                            recvSizes[i]);
                                    }
                                LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes done");
                            }
                            else if ((new_dir > 0 ? new_dir : -new_dir) == 1)
                            {
                                for (int i = gsize[k], end = gsize[k + 1]; i < end; ++i)
                                {
                                    commCells[dir].ghosts[i]->particles.resize(
                                        recvSizes[i - gsize[k]]);
                                }
                                LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes done");
                            }
                            sendSizes.clear();
                            recvSizes.clear();
                        }

                        // prepare send and receive buffers
                        longint receiver, sender;
                        outBuffer.reset();
                        // inBuffer.reset();

                        if ((new_dir > 0 ? new_dir : -new_dir) == 2)
                        {
                            if (realToGhosts)
                            {
                                receiver = cnode[k];
                                if (k == 0) sender = nodeGrid.getNodeNeighborIndex(oppositeDir);
                                // absolute x-positions (counts in cells) of current real&ghost
                                // cells
                                int cell_x1, cell_x2, itmp;
                                real cx_flag = (shift[2] > .0 ? 1.0 : -1.0);
                                for (int i = rsize[k], end = rsize[k + 1]; i < end; ++i)
                                {
                                    cell_x1 = nodeGrid.getNodePosition(0) * getInt3DCellGrid()[0] +
                                              (commCells[dir].reals[i] - getFirstCell()) %
                                                  (getInt3DCellGrid()[0] + 2) -
                                              1;
                                    cell_x2 =
                                        cnode[k] % getInt3DNodeGrid()[0] * getInt3DCellGrid()[0] +
                                        (commCells[dir].ghosts[i] - getFirstCell()) %
                                            (getInt3DCellGrid()[0] + 2) -
                                        1;
                                    itmp = static_cast<int>(
                                        floor(((cell_x1 - cell_x2) * cellGrid.getCellSize(0) +
                                               cx_flag * offs) /
                                                  Lx +
                                              0.5));
                                    packPositionsEtc_LEBC(outBuffer, *commCells[dir].reals[i],
                                                          extradata, shift,
                                                          cx_flag * offs - (itmp + .0) * Lx);
                                }

                                // exchange particles, odd-even rule
                                if (nodeGrid.getNodePosition(coord) % 2 == 0)
                                {
                                    outBuffer.send(receiver, DD_COMM_TAG);
                                    if (k == 0) inBuffer.recv(sender, DD_COMM_TAG);
                                }
                                else
                                {
                                    if (k == 0) inBuffer.recv(sender, DD_COMM_TAG);
                                    outBuffer.send(receiver, DD_COMM_TAG);
                                }
                            }
                            else
                            {
                                if (k == 0) receiver = nodeGrid.getNodeNeighborIndex(oppositeDir);
                                sender = cnode[k];

                                if (k == 0)
                                    for (int i = 0, end = commCells_bkp[dir - 4].ghosts.size();
                                         i < end; ++i)
                                    {
                                        packForces(outBuffer, *commCells_bkp[dir - 4].ghosts[i]);
                                    }

                                // exchange particles, odd-even rule
                                if (nodeGrid.getNodePosition(coord) % 2 == 0)
                                {
                                    if (k == 0) outBuffer.send(receiver, DD_COMM_TAG);
                                    inBuffer.recv(sender, DD_COMM_TAG);
                                }
                                else
                                {
                                    inBuffer.recv(sender, DD_COMM_TAG);
                                    if (k == 0) outBuffer.send(receiver, DD_COMM_TAG);
                                }
                            }
                        }
                        else if ((new_dir > 0 ? new_dir : -new_dir) == 1)
                        {
                            if (realToGhosts)
                            {
                                if (k == 0) receiver = nodeGrid.getNodeNeighborIndex(dir);
                                sender = cnode[k];

                                if (k == 0)
                                    for (int i = 0, end = commCells_bkp[dir - 4].reals.size();
                                         i < end; ++i)
                                    {
                                        packPositionsEtc(outBuffer,
                                                         *commCells_bkp[dir - 4].reals[i],
                                                         extradata, shift);
                                    }
                                // exchange particles, odd-even rule
                                if (nodeGrid.getNodePosition(coord) % 2 == 0)
                                {
                                    if (k == 0) outBuffer.send(receiver, DD_COMM_TAG);
                                    inBuffer.recv(sender, DD_COMM_TAG);
                                }
                                else
                                {
                                    inBuffer.recv(sender, DD_COMM_TAG);
                                    if (k == 0) outBuffer.send(receiver, DD_COMM_TAG);
                                }
                            }
                            else
                            {
                                receiver = cnode[k];
                                if (k == 0) sender = nodeGrid.getNodeNeighborIndex(dir);

                                for (int i = gsize[k], end = gsize[k + 1]; i < end; ++i)
                                {
                                    packForces(outBuffer, *commCells[dir].ghosts[i]);
                                }

                                // exchange particles, odd-even rule
                                if (nodeGrid.getNodePosition(coord) % 2 == 0)
                                {
                                    outBuffer.send(receiver, DD_COMM_TAG);
                                    if (k == 0) inBuffer.recv(sender, DD_COMM_TAG);
                                }
                                else
                                {
                                    if (k == 0) inBuffer.recv(sender, DD_COMM_TAG);
                                    outBuffer.send(receiver, DD_COMM_TAG);
                                }
                            }
                        }
                        else
                            throw std::runtime_error(
                                "doGhostCommunication error: new_dir should be equal to one of "
                                "[-2,-1,1,2] \n");

                        // unpack received data
                        if ((new_dir > 0 ? new_dir : -new_dir) == 2)
                        {
                            if (realToGhosts)
                            {
                                // unpack received data
                                if (k == 0)
                                    for (int i = 0, end = commCells_bkp[dir - 4].reals.size();
                                         i < end; ++i)
                                    {
                                        unpackPositionsEtc(*commCells_bkp[dir - 4].ghosts[i],
                                                           inBuffer, extradata);
                                    }
                            }
                            else
                            {
                                for (int i = rsize[k], end = rsize[k + 1]; i < end; ++i)
                                {
                                    unpackAndAddForces(*commCells[dir].reals[i], inBuffer);
                                }
                            }
                        }
                        else if ((new_dir > 0 ? new_dir : -new_dir) == 1)
                        {
                            if (realToGhosts)
                            {
                                // unpack received data
                                for (int i = rsize[k], end = rsize[k + 1]; i < end; ++i)
                                {
                                    unpackPositionsEtc(*commCells[dir].ghosts[i], inBuffer,
                                                       extradata);
                                }
                            }
                            else
                            {
                                if (k == 0)
                                    for (int i = 0, end = commCells_bkp[dir - 4].reals.size();
                                         i < end; ++i)
                                    {
                                        unpackAndAddForces(*commCells_bkp[dir - 4].reals[i],
                                                           inBuffer);
                                    }
                            }
                        }
                    }
                }
                else
                {
                    // The standard doGhostCommunication
                    if (sizesFirst)
                    {
                        LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes");

                        // prepare buffers
                        std::vector<longint> sendSizes, recvSizes;
                        sendSizes.reserve(commCells[dir].reals.size());
                        for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i)
                        {
                            sendSizes.push_back(commCells[dir].reals[i]->particles.size());
                        }
                        recvSizes.resize(commCells[dir].ghosts.size());

                        // exchange sizes, odd-even rule
                        if (nodeGrid.getNodePosition(coord) % 2 == 0)
                        {
                            LOG4ESPP_DEBUG(logger,
                                           "sending to node "
                                               << nodeGrid.getNodeNeighborIndex(dir)
                                               << ", then receiving from node "
                                               << nodeGrid.getNodeNeighborIndex(oppositeDir));

                            getSystem()->comm->send(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG,
                                                    &(sendSizes[0]), sendSizes.size());
                            getSystem()->comm->recv(nodeGrid.getNodeNeighborIndex(oppositeDir),
                                                    DD_COMM_TAG, &(recvSizes[0]), recvSizes.size());
                        }
                        else
                        {
                            LOG4ESPP_DEBUG(logger, "receiving from node "
                                                       << nodeGrid.getNodeNeighborIndex(oppositeDir)
                                                       << ", then sending to node "
                                                       << nodeGrid.getNodeNeighborIndex(dir));
                            getSystem()->comm->recv(nodeGrid.getNodeNeighborIndex(oppositeDir),
                                                    DD_COMM_TAG, &(recvSizes[0]), recvSizes.size());
                            getSystem()->comm->send(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG,
                                                    &(sendSizes[0]), sendSizes.size());
                        }

                        // resize according to received information
                        for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i)
                        {
                            commCells[dir].ghosts[i]->particles.resize(recvSizes[i]);
                        }
                        LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes done");
                    }

                    // prepare send and receive buffers
                    longint receiver, sender;
                    outBuffer.reset();
                    if (realToGhosts)
                    {
                        receiver = nodeGrid.getNodeNeighborIndex(dir);
                        sender = nodeGrid.getNodeNeighborIndex(oppositeDir);

                        for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i)
                        {
                            packPositionsEtc(outBuffer, *commCells[dir].reals[i], extradata, shift);
                        }
                    }
                    else
                    {
                        receiver = nodeGrid.getNodeNeighborIndex(oppositeDir);
                        sender = nodeGrid.getNodeNeighborIndex(dir);
                        for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i)
                        {
                            packForces(outBuffer, *commCells[dir].ghosts[i]);
                        }
                    }

                    // exchange particles, odd-even rule
                    if (nodeGrid.getNodePosition(coord) % 2 == 0)
                    {
                        outBuffer.send(receiver, DD_COMM_TAG);
                        inBuffer.recv(sender, DD_COMM_TAG);
                    }
                    else
                    {
                        inBuffer.recv(sender, DD_COMM_TAG);
                        outBuffer.send(receiver, DD_COMM_TAG);
                    }

                    // unpack received data
                    if (realToGhosts)
                    {
                        for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i)
                        {
                            unpackPositionsEtc(*commCells[dir].ghosts[i], inBuffer, extradata);
                        }
                    }
                    else
                    {
                        for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i)
                        {
                            unpackAndAddForces(*commCells[dir].reals[i], inBuffer);
                        }
                    }
                }
            }
        }
    }
    LOG4ESPP_DEBUG(logger, "ghost communication finished");
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void DomainDecomposition::registerPython()
{
    using namespace espressopp::python;
    class_<DomainDecomposition, bases<Storage>, boost::noncopyable>(
        "storage_DomainDecomposition",
        init<std::shared_ptr<System>, const Int3D&, const Int3D&, int>())
        .def("mapPositionToNodeClipped", &DomainDecomposition::mapPositionToNodeClipped)
        .def("getCellGrid", &DomainDecomposition::getInt3DCellGrid)
        .def("getNodeGrid", &DomainDecomposition::getInt3DNodeGrid)
        .def("cellAdjust", &DomainDecomposition::cellAdjust);
}

}  // namespace storage
}  // namespace espressopp
