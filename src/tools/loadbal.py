#  Copyright (C) 2017,2018(1H)
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2019
#      Max Planck Computing and Data Facility
#
#  This file is part of ESPResSo++ -> Powered by HeSpaDDA algorithm developed by horacio.v.g@gmail.com
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
**************************************************
loadbal - HeSpaDDA load balancing python functions
**************************************************


*  `qbicity(box_size,rc,skin)`:
    
    It's a function to check the system size cubicity, with a tolerance given by the rc+skin
    `dLnorm` - (Lx,Ly,Lz)/Lmax 

*  `changeIndex(dN,ima,imi)`:

    It's a function that sorts the nodeGrid, according to the index of the maximum size of the 
    box (ima) and its corresponding minimum (imi)
    
*  `nodeGridSizeCheck(node_gridX,node_gridY,node_gridZ)`:

    It's a function that verifies if it is worthy to take care of special DD for inhomogeneous 	  
    systems, otherwise the homogenous DD is triggered and returned as a xyz flags
    nodeGridSizeCheck(node_gridX,node_gridY,node_gridZ)

*  `halfDecomp(adrCenter1D,rc_skin,eh_size,halfCores1D,cellsX,ratioMS,idealGas)`:

    It's a function that decomposes half of the box in one Dimension(1D) and it is assuming symmetry
    in the initial simulation box. Even with irregularities in the DD it will find the half-Box-DD 

*  `addHsymmetry(halfNeilListX,eh_size,rc_skin,node_gridX,cellsX,ratioMS,idealGas)`:
    
    Its a function that uses the previously decomposed half-box (from halfDecomp) and unfolds it to
    match the whole Neighbors List. If it does not match the whole neighbor, due to any asymmetry in
    the number of cores or the number of cells. It will be redistributed region by region by  
    considering the whole simulation box size.

*  `adaptNeiList(neiListxin)`:
    
    It's a function which adapts the number of cells that go into each core into the data structure of
    left and right cell lists for example for 4 cores and 8 cells [3,4,5,8] to [0,3,4,5,8]

*  `reDistCellsHom(node_gridX,sizeX,rc_skin)`:
    
    It's a function which distributes the cells into nodes as if they where homogeneous. It also
    applies to inhomogeneous system whenever there are less than 2 cores per direction: X, Y or Z.

*  `reDistCells(halfNeilListX,cellsX,eh_size,rc_skin,node_gridX,ratioMS,idealGas)`:
    
    It's a function which is matching proportion of cells to the cores on a dual resolution region
    basis

*  `redistDeltaRandomly(wholeNeiListX,deltaCells,totNodesEH=0,biased=0)`:
    
    It's a function which distributes the remaining DELTA cells into nodes semi-randomly. By default
    the biase applies to the CG-region and it assumes ther cannot be more than 3 extra cells to
    redistribute, because this is the total number of regions in the simulation box `|CG|EH|CG|` (by
    default the cg biased is left this could be updated in the dyn load balancing case!

*  `findNodesMS(node_gridX,totCellsEH,totCellsCG,ratioMS,idealGas)`:
    
    It's a function which normalizes the number of cells to go to the EH and CG regions and find the
    ideal corresponding number of Nodes EH and CG

"""

from random import randint

__author__ = 'Dr. Horacio V Guzman'
__email__ = 'horacio.v.g at gmail dot com'
__version__ = '1.0'
__all__ = [
    'qbicity', 'changeIndex', 'nodeGridSizeCheck',
    'halfDecomp', 'addHsymmetry', 'adaptNeiList',
    'reDistCellsHom', 'reDistCells', 'redistDeltaRandomly',
    'findNodesMS'
]

# This function verifies if the simuation box is of cubic dimensions


def qbicity(box_size, rc, skin, cellDomTol=2.0):
    rc_skin = rc + skin
    box_aux = [box_size[0], box_size[1], box_size[2]]
    indMax = box_aux.index(max(box_aux))
    indMin = box_aux.index(min(box_aux))
    if (box_size[indMax] - box_size[indMin]) < (cellDomTol * rc_skin):
        qFlag = bool(1)
    else:
        qFlag = bool(0)
    return qFlag

# This function sorts the nodeGrid, according to the index of the maximum size of the box (indMax) and its corresponding minimum (indMin)


def changeIndex(dN, indMax, indMin):
    aux = [0, 0, 0]
    ndN = [0, 0, 0]
    dIndMax = dN.index(max(dN))
    dIndMin = dN.index(min(dN))
    ndN[indMax] = dN[dIndMax]
    ndN[indMin] = dN[dIndMin]
    listInd = range(3)
    if indMax > indMin:
        listInd.pop(indMax)
        listInd.pop(indMin)
    else:
        listInd.pop(indMin)
        listInd.pop(indMax)
    aux = dN[:]
    if dIndMax > dIndMin:
        aux.pop(dIndMax)
        aux.pop(dIndMin)
    else:
        aux.pop(dIndMin)
        aux.pop(dIndMax)
    ndN[listInd[0]] = aux[0]
    return ndN

# This function verifies if it is worthy to take care of special DD for inhomogeneous systems, otherwise the homogenous DD is triggered and returned as a xyz flags


def nodeGridSizeCheck(node_gridX, node_gridY, node_gridZ):
    flx = 0
    fly = 0
    flz = 0
    if node_gridX < 3:
        flx = node_gridX
    else:
        flx = 0
    if node_gridY < 3:
        fly = node_gridY
    else:
        fly = 0
    if node_gridZ < 3:
        flz = node_gridZ
    else:
        flz = 0
    return flx, fly, flz

# This function decomposes half of the box in one Dimension(1D) and it is assuming symmetry in the initial simulation box. Even with irregularities in the box dimension and cell-size relation it will find the half-Box-DD


def halfDecomp(adrCenter1D, rc_skin, eh_size, halfCores1D, cellsX, ratioMS, sizeX, idealGas, halfCellInt = 1):
    # this value is only in case the Ideal Gas will in reality improve any calculation or communication (i.e. Improve notoriously the sims parallelization, which is not the case yet)
    pLoadIG = 1
    cellSizes = []
    usedCores = 0
    totCellsEH = halfCellInt * round(2. * eh_size / rc_skin - 0.5)
    totCellsCG = cellsX - totCellsEH
    totNodesCG, totNodesEH = findNodesMS(halfCores1D * 2, totCellsEH, totCellsCG, ratioMS, sizeX, eh_size, idealGas)
    for i in xrange(halfCores1D):

        if idealGas:
            if i == 0:
                # 2do: SuperCell stuff
                cellSizes.append(halfCellInt * round((adrCenter1D - eh_size) / rc_skin - 0.5) + pLoadIG)
                usedCores = 1  # For Ideal Gas purposes only 1 core covers the low-resolution region
            else:
                [cellSizes.append(halfCellInt * round((eh_size) / rc_skin / (halfCores1D - usedCores) - 0.5)) for i in xrange(usedCores, halfCores1D)]  # 2do: SuperCell stuff
                deltaCells = halfCellInt * round((eh_size) / rc_skin - 0.5) - halfCellInt * round((eh_size) / rc_skin / (halfCores1D - usedCores) - 0.5) * (halfCores1D - usedCores)  # 2do: SuperCell stuff
                # Both applies a benefit to the usedCores=1 in the sense that will have one cell less loaded to CG-region. But also a penalty to the usedCores=1 that will have to manage any cells non distributed before
                cellSizes[usedCores] = cellSizes[usedCores] + deltaCells - pLoadIG
                return cellSizes

        else:  # Applies to all other systems besides the Ideal Gas
            if totNodesCG / 2. >= 2.:
                [cellSizes.append(halfCellInt * round((adrCenter1D - eh_size) / rc_skin / totNodesCG / 2. - 0.5)) for j in xrange(int(totNodesCG / 2.))]  # 2do: SuperCell stuff
            else:
                # 2do: SuperCell stuff
                cellSizes.append(halfCellInt * round((adrCenter1D - eh_size) / rc_skin - 0.5))
            if totNodesEH / 2. >= 2.:
                [cellSizes.append(halfCellInt * round((eh_size) / rc_skin / (totNodesEH / 2.) - 0.5)) for i in xrange(int(totNodesEH / 2.))]  # 2do: SuperCell stuff
            else:
                # 2do: SuperCell stuff
                cellSizes.append(halfCellInt * round((eh_size) / rc_skin - 0.5))
            return cellSizes
    return cellSizes

# This function uses the previously decomposed half-box (from halfDecomp) and unfolds it to match the whole Neighbors List. If it does not match the whole neighbor, due to any asymmetry in the number of cores or the number of cells. It will be redistributed region by region by considering the whole simulation box size.


def addHsymmetry(halfNeilListX, eh_size, rc_skin, node_gridX, cellsX, ratioMS, sizeX, idealGas, halfCellInt = 1):
    wholeNeilListX = []
    aux = halfNeilListX[:]
    # unfolds halfDecomp and attempts to match the whole neighbor list.
    aux.reverse()
    # here half neighbor list X turns into whole neighbor list
    halfNeilListX.extend(aux)
    aux2 = 0
    if len(halfNeilListX) < node_gridX:
        if (node_gridX - len(halfNeilListX)) == 1:  # Verifies if number of cores are even
            aux2 = halfNeilListX[len(halfNeilListX) - 1]
            halfNeilListX.append(round(halfNeilListX[len(halfNeilListX) - 1] / 2. - 0.5))
            halfNeilListX[len(halfNeilListX) - 2] = aux2 - halfNeilListX[len(halfNeilListX) - 1]
            if sum(halfNeilListX) != cellsX:
                wholeNeilListX = reDistCells(halfNeilListX, cellsX, eh_size, rc_skin, node_gridX, ratioMS, sizeX, idealGas, halfCellInt)
            else:
                if any([v == 0 for v in halfNeilListX]):  # Recently added 138, 139 and 140
                    wholeNeilListX = reDistCells(halfNeilListX, cellsX, eh_size, rc_skin, node_gridX, ratioMS, sizeX, idealGas, halfCellInt)
                else:
                    print "HeSpaDDA message: addHsymmetry all tests passed although not all cells were used"
                    wholeNeilListX = halfNeilListX[:]
        else:
            print "HeSpaDDA message: The distributed cores are not matching the available ones (++ reDistCells())"
            # in the original implementation of halfCell the halfCellInt was ommited in the following call
            # I think, however, that it was just forgotten, so I put it here like to all other calls to reDistCells
            wholeNeilListX = reDistCells(halfNeilListX, cellsX, eh_size, rc_skin, node_gridX, ratioMS, sizeX, idealGas, halfCellInt) 
            # To be determined if additional reDitsCells should be called!
    elif len(halfNeilListX) == node_gridX and aux2 == 0:
        if sum(halfNeilListX) != cellsX:
            print "HeSpaDDA message: The distributed cells are not matching the available ones"
            wholeNeilListX = reDistCells(halfNeilListX, cellsX, eh_size, rc_skin, node_gridX, ratioMS, sizeX, idealGas, halfCellInt)
        else:
            # Recently added 152, 153 and 154
            if any([v == 0 for v in halfNeilListX]):
                wholeNeilListX = reDistCells(halfNeilListX, cellsX, eh_size, rc_skin, node_gridX, ratioMS, sizeX, idealGas, halfCellInt)
            else:
                print "HeSpaDDA message: addHsymmetry all tests passed although not all cells are used"
                wholeNeilListX = halfNeilListX[:]
    else:
        print "HeSpaDDA message: The distributed cores are not matching the available ones", halfNeilListX
        halfNeilListX[len(halfNeilListX) - 2] = halfNeilListX[len(halfNeilListX) - 1] + halfNeilListX[len(halfNeilListX) - 2]
        halfNeilListX.pop(len(halfNeilListX) - 1)
        print "HeSpaDDA message: During DD a core has been reduced"
        wholeNeilListX = reDistCells(halfNeilListX, cellsX, eh_size, rc_skin, node_gridX, ratioMS, sizeX, idealGas, halfCellInt)
    return wholeNeilListX

# This function adapts the number of cells that go into each core into the data structure of left and right cell lists for example for 4 core and 8 cells [3,4,5,8] to [0,3,4,5,8]


def adaptNeiList(neiListxin):
    neiListx = []
    neiListx.append(0)
    [neiListx.append(neiListxin[i] + neiListx[i]) for i in xrange(len(neiListxin) - 1)]
    neiListx.append(neiListxin[len(neiListxin) - 1] + neiListx[len(neiListx) - 1])
    print "HeSpaDDA message: Your Cells Neighbor Lists is:", neiListx
    return neiListx

# This function distributes the cells into nodes as if they where homogeneous. It also applies to inhomogeneous system whenever there are less than 2 cores per direction: X, Y or Z.


def reDistCellsHom(node_gridX, sizeX, rc_skin, halfCellInt = 1):
    wholeNeiListX = []
    cellsX = halfCellInt * int(round(sizeX / rc_skin - 0.5))
    if node_gridX % 2 == 0 and cellsX % 2 == 0:
        [wholeNeiListX.append(cellsX / node_gridX) for i in xrange(node_gridX)]
    elif node_gridX % 2 != 0 and cellsX % 2 != 0:
        [wholeNeiListX.append(round((cellsX) / node_gridX - 0.5)) for i in xrange(node_gridX)]
        if int(cellsX - sum(wholeNeiListX)) != 0:
            # passing Delta as cellsX-sum(wholeNeiListX)
            wholeNeiListX = redistDeltaRandomly(wholeNeiListX, cellsX - sum(wholeNeiListX), 0)
        else:
            print "HeSpaDDA message: PASS appears...here, take a look at this value Px/Cx", round((cellsX) / node_gridX - 0.5)
            pass
    else:
        if node_gridX % 2 == 0 and cellsX % 2 != 0:
            [wholeNeiListX.append((cellsX - 1) / node_gridX) for i in xrange(node_gridX)]
            # Punishing the last one
            wholeNeiListX[node_gridX - 1] = wholeNeiListX[node_gridX - 1] + 1
        elif cellsX % 2 == 0 and node_gridX % 2 != 0:
            [wholeNeiListX.append(round((cellsX) / node_gridX - 0.5)) for i in xrange(node_gridX)]
            wholeNeiListX = redistDeltaRandomly(wholeNeiListX, cellsX - sum(wholeNeiListX), 0)
    return wholeNeiListX

# This This function is matching proportion of cells to the cores on a dual resolution region basis


def reDistCells(halfNeilListX, cellsX, eh_size, rc_skin, node_gridX, ratioMS, sizeX, idealGas, halfCellInt = 1):
    #global preFactCen
    preFactCen = 1.0
    print "HeSpaDDA message: Cells redistribution will improve whenever the Cells1D are at least twice as big as Nodes1D!"
    wholeNeiListX = []
    totCellsEH = halfCellInt * round(2. * eh_size / rc_skin - 0.5)
    totCellsCG = cellsX - totCellsEH
    totNodesCG, totNodesEH = findNodesMS(node_gridX, totCellsEH, totCellsCG, ratioMS, sizeX, eh_size, idealGas)
    print "HeSpaDDA message: Cores in Both LR and HR, are:", totNodesCG, totNodesEH
    if idealGas:	  # This represents the Ideal Gas (IG)!!! (OJO)
        wholeNeiListX_EH = []
        wholeNeiListX_CG = []
        wholeNeiListX = []
        # High Resolution region
        if totNodesEH % 2 == 0 and totCellsEH % 2 == 0:
            [wholeNeiListX_EH.append(totCellsEH / totNodesEH) for i in xrange(totNodesEH)]
            print "HeSpaDDA message IG: HR region: P and C are EVEN, given by:"
            print wholeNeiListX_EH

        elif totNodesEH % 2 != 0 and totCellsEH % 2 != 0:
            [wholeNeiListX_EH.append(round(totCellsEH / totNodesEH - 0.5)) for i in xrange(totNodesEH)]
            if int(totCellsEH - sum(wholeNeiListX_EH)) != 0:
                wholeNeiListX_EH[0:totNodesEH] = redistDeltaRandomly(wholeNeiListX_EH[0:totNodesEH], totCellsEH - sum(wholeNeiListX_EH[0:totNodesEH]), 0)
            else:
                print "HeSpaDDA message IG: HR region: PASS appears...here, take a look at this value Px/Cx", round(totCellsEH / totNodesEH - 0.5)
                pass
        else:
            if totNodesEH % 2 == 0 and totCellsEH % 2 != 0:
                [wholeNeiListX_EH.append((totCellsEH - 1) / totNodesEH) for i in xrange(totNodesEH)]
                wholeNeiListX_EH[totNodesEH - 1] = wholeNeiListX_EH[totNodesEH - 1] + 1
                print "HeSpaDDA message IG: HR region: P and noC are EVEN"
            elif totCellsEH % 2 == 0 and totNodesEH % 2 != 0:
                [wholeNeiListX_EH.append(round((totCellsEH) / totNodesEH - 0.5)) for i in xrange(totNodesEH)]
                # passing Delta cells to be redistributed semi-randomly (after prioritizying the EH-region, additional cells should go to the CG-region).
                wholeNeiListX_EH[0:totNodesEH] = redistDeltaRandomly(wholeNeiListX_EH[0:totNodesEH], totCellsEH - sum(wholeNeiListX_EH[0:totNodesEH]), 0)
                print "HeSpaDDA message IG: HR region: noP and C are EVEN"
        #@@@ Low Resolution region
        if totNodesCG % 2 == 0 and totCellsCG % 2 == 0:
            [wholeNeiListX_CG.append(totCellsCG / totNodesCG) for i in xrange(totNodesCG)]
            print "HeSpaDDA message IG: LR region: P and C are EVEN, given by:"
            print wholeNeiListX_CG
        elif totNodesCG % 2 != 0 and totCellsCG % 2 != 0:
            [wholeNeiListX_CG.append(round(totCellsCG / totNodesCG - 0.5)) for i in xrange(totNodesCG)]
            if int(totCellsCG - sum(wholeNeiListX_CG)) != 0:
                wholeNeiListX_CG[0:totNodesCG] = redistDeltaRandomly(wholeNeiListX_CG[0:totNodesCG], totCellsCG - sum(wholeNeiListX_CG[0:totNodesCG]), 0)
            else:
                print "HeSpaDDA message IG: LR region: PASS appears...here, take a look at this value Px/Cx", round(totCellsCG / totNodesCG - 0.5)
                pass
        else:
            if totNodesCG % 2 == 0 and totCellsCG % 2 != 0:
                [wholeNeiListX_CG.append((totCellsCG - 1) / totNodesCG) for i in xrange(totNodesCG)]
                wholeNeiListX_CG[totNodesCG - 1] = wholeNeiListX_CG[totNodesCG - 1] + 1
                print "HeSpaDDA message IG: LR region: P and noC are EVEN"
                print wholeNeiListX_CG
            elif totCellsCG % 2 == 0 and totNodesCG % 2 == 0:
                [wholeNeiListX_CG.append(round((totCellsCG) / totNodesCG - 0.5)) for i in xrange(totNodesCG)]
                # passing Delta cells to be redistributed semi-randomly (after prioritizying the EH-region, additional cells may come to the CG-region).
                wholeNeiListX_CG[0:totNodesCG] = redistDeltaRandomly(wholeNeiListX_CG[0:totNodesCG], totCellsCG - sum(wholeNeiListX_CG[0:totNodesCG]), 0)
                print "HeSpaDDA message IG: LR region: noP and C are EVEN"

        # Index of the middle LR region begin of HR
        indCG1 = int((len(wholeNeiListX_CG)) / 2)
        print "HeSpaDDA message indexing: The CG first subregion index is:", indCG1
        # Index of the start of the second LR region end of HR
        indEH1 = indCG1 + int(totNodesEH)
        # Ensembling the array of Cells Neighbors list
        wholeNeiListX.extend(wholeNeiListX_CG[0:indCG1])
        wholeNeiListX.extend(wholeNeiListX_EH)
        wholeNeiListX.extend(wholeNeiListX_CG[indCG1:len(wholeNeiListX_CG)])

    #@ not Ideal Gas
    else:
        if (totNodesCG + totNodesEH) > node_gridX:
            # minimum number of cores(nodes=cores) for the CG region version 1.0
            if totNodesCG > 2:
                # Assign exceeding number of cores to LR
                totNodesCG = totNodesCG + \
                    ((totNodesCG + totNodesEH) - node_gridX)
                totNodesEH = node_gridX - totNodesCG
            else:
                totNodesCG = 2  # At least use 2 core for the CG region
                totNodesEH = node_gridX - totNodesCG
        else:
            print "HeSpaDDA message indexing: Nodes CG and Nodes EH are respectively,", totNodesCG, totNodesEH
        if node_gridX % 2 == 0 and cellsX % 2 == 0:
            wholeNeiListX_EH = [0] * (node_gridX)
            wholeNeiListX_CG = [0] * (node_gridX)
            wholeNeiListX = [0] * (node_gridX)
            if totNodesCG % 2 == 0:
                indCG1 = int(totNodesCG / 2)
                indEH1 = indCG1 + int(totNodesEH)
                au1 = range(int(indCG1))
                # au1 contains the index of CG cells in the whole neighbor list
                au1.extend(range(indEH1, indEH1 + indCG1))
                # internal parameter which in case of dynamic LB, if more weight of DD goes to the center by defaults (EH -region Flag =1) and hence DD focus on any CG-region (default value decomposes cells to CG).	preFactCen=ratioMS per default(OLD)
                centralFlagEH = 1
                # new stuff TO BE CHECKED
                if centralFlagEH > 0:  # H's WL as well    #NHEW
                    for i in xrange(int(ratioMS), int(cellsX + 1)):
                        tempWNL = [0] * (node_gridX)
                        ratioMS2t = round(1. * (cellsX / (1. * pow(i, 1. / 3.))) - 0.5)
                        print "HeSpaDDA message indexing: the Ratio MS2Cellslot 'cells weight' in the CG region is:", ratioMS2t
                        for j in au1:  # This loop goes over the CG-regions
                            tempWNL[j] = round(ratioMS2t * totCellsCG / totNodesCG - 0.5)
                        totCellsEHtemp = cellsX - sum(tempWNL)
                        if totCellsEHtemp < totNodesEH and totCellsEHtemp > 0:
                            print "HeSpaDDA message indexing: Error with pass Factor MS-2-Cells, no worries HeSpaDDA will find another for you!", totCellsEHtemp
                        else:
                            # i was not the cubic root... yet
                            preFactCen = pow(i, 1. / 3.)
                            # if int(totCellsEHtemp)-1==int(totNodesEH):
                            #	totCellsCG=totCellsCG+1
                            print "HeSpaDDA message indexing: The rescaling preFactor for the Cells distribution is...an optimized value, like this:", preFactCen
                            break
                            break  # the first value found for preFactCen takes you out of the for
                # Now redistributing the cells to nodes with a proper dimension factor as a f(dimensions,ratioMS, coresEH, coresCG)
                for i in au1:
                    numRegBox = 3.  # 3. it is a region based parameter |CG|EH|CG| =3 fixed param, number of regions in the box
                    if cellsX > numRegBox * pow(ratioMS, 1. / 3.) and totNodesEH < cellsX - (pow(ratioMS, 1. / 3.) * totCellsCG):
                        wholeNeiListX_CG[i] = round(pow(ratioMS, 1. / 3.) * totCellsCG / totNodesCG - 0.5)
                        print "HeSpaDDA message LR: cells dist No IG if cells fit in 3 subregions..."
                    else:
                        ratioMS2 = round(pow(1. * (cellsX / (1. * preFactCen)), 1. / 3.) - 0.5)
                        volRatioX = (sizeX - 2. * eh_size) / sizeX  # Check eh_size or 2*eh_size
                        wholeNeiListX_CG[i] = round(ratioMS2 * volRatioX * totCellsCG / totNodesCG - 0.5)
                        print "HeSpaDDA message LR: cells dist No IG if cells fit volRatio used..."
                totCellsEH = cellsX - sum(wholeNeiListX_CG)
                print "HeSpaDDA message LR: wholeNeiListX_CG is, ", wholeNeiListX_CG
                # Now the redist in the EH-region occurs | NHEW-> We still need to check if the cells are enough for the EH cores
                if totNodesEH % 2 == 0 and totCellsEH >= totNodesEH:
                    for i in xrange(indCG1, indEH1):
                        wholeNeiListX_EH[i] = round(1.0 * totCellsEH / totNodesEH - 0.5)
                elif totNodesEH % 2 != 0 and totCellsEH % 2 == 0:
                    for i in xrange(indCG1, indEH1):
                        wholeNeiListX_EH[i] = round(1.0 * (totCellsEH - 1) / totNodesEH - 0.5)
                    # Punishing the last Node with an additional cell
                    wholeNeiListX[indEH1 - 1] = wholeNeiListX_EH[indEH1 - 1] + 1
                # print "Whole lists are as:",wholeNeiListX_EH,wholeNeiListX_CG
                for k in xrange(node_gridX):
                    wholeNeiListX[k] = wholeNeiListX_EH[k] + wholeNeiListX_CG[k]  # Superposing both Arrays

            else:  # TO BE IMPROVED (not fullfilling all nodes)!!! not EVEN number of nodes!
                indCG1 = int(totNodesCG / 2.0)  # gives 1
                indEH1 = indCG1 + int(totNodesEH)  # gives 6
                au2 = range(int(indCG1))  # 1
                # no assuming odd totNodesCG, before (indEH1,indEH1+indCG1))
                au2.extend(range(indEH1, indEH1 + indCG1 + 1))
                print "HeSpaDDA message: Cells CG wrong", totCellsCG
                if int(totCellsCG) % int(totNodesCG) == 0:  # NHEW
                    for i in au2:
                        wholeNeiListX_CG[i] = round(1.0 * (totCellsCG) / totNodesCG - 0.5)
                else:
                    for i in au2:
                        wholeNeiListX_CG[i] = round(1.0 * (totCellsCG - 1) / totNodesCG - 0.5)
                    # Punishing the last Node with an additional cell  NHEW (got rid of -1 in the index)
                    wholeNeiListX_CG[indEH1 + indCG1] = wholeNeiListX_CG[indEH1 + indCG1] + 1
                if totNodesEH % 2 == 0:
                    for i in xrange(indCG1, indEH1):
                        wholeNeiListX_EH[i] = round(1.0 * totCellsEH / totNodesEH - 0.5)
                else:
                    for i in xrange(indCG1, indEH1):
                        wholeNeiListX_EH[i] = round(1.0 * (totCellsEH - 1) / totNodesEH - 0.5)
                    # Punishing the last Node with an additional cell
                    wholeNeiListX_EH[indEH1 - 1] = wholeNeiListX_EH[indEH1 - 1] + 1
                for k in xrange(node_gridX):
                    wholeNeiListX[k] = wholeNeiListX_EH[k] + wholeNeiListX_CG[k]
        else:  # TO BE IMPROVED
            if node_gridX % 2 == 0:
                [wholeNeiListX.append(round((cellsX - 1) / node_gridX - 0.5)) for i in xrange(node_gridX)]   # Homogeneously distributed
                # Punishing the last Node with an additional cell
                wholeNeiListX[node_gridX - 1] = wholeNeiListX[node_gridX - 1] + 1
            elif cellsX % 2 == 0:
                [wholeNeiListX.append(round((cellsX) / node_gridX - 0.5)) for i in xrange(node_gridX)]
                wholeNeiListX = redistDeltaRandomly(wholeNeiListX, cellsX - sum(wholeNeiListX), totNodesEH, cellsX - sum(wholeNeiListX) - 1)
    # print "My Redist WholeNeiList is TODOs !:",wholeNeiListX
    return wholeNeiListX

# This function distributes the remaining DELTA cells into nodes as semi randomly. By default the biase applies to the CG-region and it assumes ther cannot be more than 3 extra cells to redistribute, because this is the total number of regions in the simulation box |CG|EH|CG| (by default the cg biased is left this could be updated in the dyn load balancing case!


def redistDeltaRandomly(wholeNeiListX, deltaCells, totNodesEH=0, biased=0):
    flagBiased = 0
    wholeNeiListXcopy = wholeNeiListX[:]
    index = len(wholeNeiListX) - 1
    indexOut = [0] * int(deltaCells)
    print "HeSpaDDA message: This are the deltaCells", deltaCells
    if deltaCells > 0.5:
        indexOut[-1] = 3  # initialization value for the index of the nodes that will get more cells, so that the random number generator is never punishing the same node with more cells
    else:
        indexOut = [0]
    return wholeNeiListXcopy
    if totNodesEH == 0:
        for p in xrange(0, int(deltaCells)):
            aux2 = randint(0, index)
            while aux2 == indexOut[p - 1]:
                aux2 = randint(0, index)
            indexOut[p] = aux2
        for i in indexOut:
            wholeNeiListXcopy[i] = wholeNeiListX[i] + 1
    else:
        for p in xrange(0, int(deltaCells)):
            index = len(wholeNeiListX) - 1
            if biased > 0 and biased < 3: 	            # Left biased!
                # Left CG region | * |  |  |
                aux2 = randint(0, index - totNodesEH - 1)
                nIndMin = 0
                if biased > 1 or flagBiased == 1:
                    nIndMax = index - totNodesEH
                    flagBiased = 1
                else:
                    nIndMax = index - totNodesEH - 1
                biased = biased - 1
            else:
                # Right CG region |  |  | * |
                aux2 = randint(totNodesEH + 1, index)
                nIndMin = totNodesEH + 1
                nIndMax = index
            while aux2 == indexOut[p - 1]:
                aux2 = randint(nIndMin, nIndMax)
            indexOut[p] = aux2
        for i in indexOut:
            wholeNeiListXcopy[i] = wholeNeiListX[i] + 1
    return wholeNeiListXcopy

# This function normalizes the number of cells to go to the EH and CG regions and find the ideal corresponding number of Nodes EH and CG


def findNodesMS(node_gridX, totCellsEH, totCellsCG, ratioMS, sizeX, eh_size, idealGas, procsWEH=1.):
    fRatioEH = pow(ratioMS, 1. / 3.) * (2.0 * eh_size / (1.0 * (sizeX) + 2.0 * eh_size * (pow(ratioMS, 1. / 3.) - 1.)))  # Seems to be wo Bu!
    # pow(ratioMS,1./3.)*(1.0*totCellsEH/(1.0*(totCellsCG+totCellsEH)+totCellsEH*(pow(ratioMS,1./3.)-1.)))
    # fRatioCG=(1./1.)*(1.0*totCellsCG/(1.0*(totCellsCG+totCellsEH)))
    if idealGas:
        if (node_gridX - 2.) <= totCellsEH and node_gridX > totCellsEH:
            # 2do: Could be tuned for every case! SuperCell!
            totNodesEH = round(totCellsEH / procsWEH - 0.5)
            totNodesCG = node_gridX - totNodesEH
        elif node_gridX > totCellsEH + 2:
            totNodesEH = totCellsEH  # 2do: Could be tuned for every case! SuperCell!
            totNodesCG = node_gridX - totNodesEH
        else:
            totNodesEH = node_gridX - 2  # 2do: Could be tuned for every case! SuperCell!
            totNodesCG = node_gridX - totNodesEH
            if totNodesEH < 1 and node_gridX > 0:
                print "HeSpaDDA message: You are using the minimum amount of cores!!!"
                totNodesEH = 1
                totNodesCG = 1
            else:
                print "HeSpaDDA message: Are you sure you need to use a Domain Decomposition? Verify that you are not trying to run this simulation on a single core"

    else:  # Applies to all other systems besides the Ideal Gas
        if node_gridX <= (totCellsEH + totCellsCG):
            totNodesEH = round(fRatioEH * node_gridX)
            print "HeSpaDDA message: According to the theory of HV Guzman article P_{HR} is :", totNodesEH
            totNodesCG = node_gridX - totNodesEH
            if (totNodesEH + totNodesCG) != node_gridX:
                # If there are more nodes than cells in EH=> redistribute nodes to EH and CG
                if totNodesEH > (totCellsEH):
                    diffNodesCells = totNodesEH - totCellsEH
                    # Add some cores to the CG nodes
                    if diffNodesCells <= (totCellsCG - totNodesCG):
                        # more weight in terms of cores in the LR region
                        totNodesCG = totNodesCG + diffNodesCells
                        totNodesEH = totNodesEH - diffNodesCells
                    else:
                        print "HeSpaDDA message: You seem to have more Cores than Cells! Hint(H): reduce the Nr. of Cores or use cherrypickTotalProcs function!"
                # If there are more nodes than cells in LR=> redistribute nodes to EH and CG
                elif totNodesCG > (totCellsCG):
                    diffNodesCells = totNodesCG - totCellsCG
                    if diffNodesCells <= (totCellsEH - totNodesEH):
                        totNodesCG = totNodesCG - diffNodesCells
                        # more weight in terms of cores in the HR region
                        totNodesEH = totNodesEH + diffNodesCells
                        if totNodesEH > totCellsEH:
                            print "HeSpaDDA message: Reduce the number of Processors to be used or try with cherrypickTotalProcs function!"
                    else:
                        print "HeSpaDDA message: You seem to have more Cores than Cells! Hint(H): reduce the Nr. of Cores"
            else:  # Everything seems to be fine, now look the size of NodesCG
                if totNodesCG < 2:  # Verify if the CG Nodes could built at least 2 CG regions, one left and one right according to its geometry
                    if totNodesCG == 1:
                        totNodesEH = totNodesEH - 1
                        totNodesCG = 2
                    else:
                        totNodesEH = totNodesEH - 2
                        totNodesCG = 2
                else:
                    pass
        else:
            print "HeSpaDDA message: You seem to have more Cores than Cells! Hint(H): reduce the Nr. of Cores"
    return totNodesCG, totNodesEH
