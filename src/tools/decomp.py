#  Copyright (C) 2012,2013,2017,2018(1H)
#      Max Planck Institute for Polymer Research
#  This file is part of ESPResSo++ > Powered by HeSpaDDA algorithm developed and conceived by horacio.v.g@gmail.com
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
***********************************************
decomp - Domain Decomposition python functions
***********************************************


*  `nodeGrid(n,box_size,rc,skin,eh_size=0,ratioMS=0,idealGas=0,slabMSDims=[0,0,0])`:

    It determines how the processors are distributed and how the cells are arranged. 
    The algorithm is dimensional sensitive for both homogeneous and inhomogeneous setups. 
    On top of such functionality it presents specific features for region divided heterogenous
    setups (e.g. for AdResS) [see H.V. Guzman et. al, Phys. Rev. E, 96, 053311 (2017)]
    Link to the paper: https://doi.org/10.1103/PhysRevE.96.053311

    `box_size` - 3D vector cointainig the size of the simulation box box_size [L_x,L_y,L_z]
    `rc` - cutoff radius of interaction
    `skin` - skin size for the verlet list calculation
    `n` - total number of processes
    `eh_size` - 1D length of the high-resolution region (e.g. for AdResS Atomistic/Explicit+Hybrid regions)
    `ratioMS` - spatial mapping ratio between high-resolution region and the low-resolution one
		(e.g. for AdResS mapping the atomistic water molecule to the CG-model; leads to a ratioMS=3)
    `idealGas` - this is a Flag for treating the low-resolution region as an Ideal Gas, when TRUE (no interactions included,
		 thus none force computations load)
    `slabMSDims` - 3D vector cointainig flags describing the type of axis, if heterogeneous value is 1, else 0

*  `nodeGridSimple(n)`:
    It determines how the processors are distributed and how the cells are arranged. Note: Use it exclusively for Lattice-Boltzmann simulations, or non-parallelized tests.
    `n` - number of processes 

*  `cherrypickTotalProcs(box_size,rc,skin,MnN,CpN,percTol=0.2,eh_size=0,ratioMS=0,idealGas=0,slabMSDims=[0,0,0])`:  
    
    To be used for heterogenous simulations where the spatial heterogeinity is known on an a-priori manner,
    where this function returns `n` as the total number of processes to be used for the best decomposition of the
    system as a function of a tolerance ratio which depending on the giving range [0, MnN*CpN] of processors availability
    different combinations of P_x,P_y and P_z can be found and hence several values of `n` this 'n' can become an array.
    Most of the parameters have been described for nodeGrid(...), except:
    
    `MnN` - M number of Nodes to be available (e.g. 128 processes/cores, 16 cores per Node gives a total of 8 Nodes)
    `CpN` - C number of Cores available in each Node (e.g. 16 Cores per Node or 20 Cores per Node)
    `percTol` - Axis base tolerance percentage to the ideal distribution of P-processors per axes P_x,P_y,P_z

*  `neiListHom(node_grid,box,rc,skin)`:

    The new domain decomposition divides the subdomains in a neighborlist of corse in a grid of 3 arrays [N_x,N_y,N_z]. In
    this case, the neighbor list is homogeneous (non a-priori load imbalance). Most of the parameters have been described 
    above, except:
    
    `node_grid` - M number of Nodes to be available (e.g. 128 processes/cores, 16 cores per Node gives a total of 8 Nodes)

*  `neiListAdress(node_grid, cell_grid,rc,skin,eh_size,adrCenter,ratioMS,idealGasFlag=True,sphereAdr=False,slabMSDims=[1,0,0])`:

    The new domain decomposition divides the subdomains in a neighborlist of corse in a grid of 3 arrays [N_x,N_y,N_z]. In
    this case, the neighbor list is homogeneous (non a-priori load imbalance). Most of the parameters have been described 
    above, except:
    
    `cell_grid` - Based on the homogenous allocation of cells per subdomain, a referential value
    `adrCenter` - Box center of the heterogeneous simulation box [adrCenter_x,adrCenter_y,adrCenter_z], commonly the middle of
                  high-resolution region.
    `idealGasFlag` - This is a Flag for treating the low-resolution region as an Ideal Gas (no interactions included,
		 thus less load)
    `sphereAdr` - Geometry of the high-resolution region, if TRUE spherical, otherwise Slab-like
    
*  `cellGrid(box_size, node_grid, rc, skin)`:

    It returns an appropriate grid of cells.
    
*  `tuneSkin(system, integrator, minSkin=0.01, maxSkin=1.2, precision=0.001)`:

    It tunes the skin size for the current system
    
*  `printTimeVsSkin(system, integrator, minSkin=0.01, maxSkin=1.5, skinStep = 0.01)`:
    
    It prints time of running versus skin size in the range [minSkin, maxSkin] with
    the step skinStep
"""


import sys
import espressopp

from espressopp import Int3D
from espressopp.Exceptions import Error

import math
import time
import numbers
from loadbal import qbicity, changeIndex, halfDecomp, addHsymmetry, adaptNeiList, reDistCellsHom, nodeGridSizeCheck

__author__ = 'Dr. Horacio V Guzman'
__email__ = 'horacio.v.g at gmail dot com'
__version__ = '1.0'
__all__ = [
    'nodeGrid', 'cellGrid', 'cherrypickTotalProcs',
    'neiListHom', 'neiListAdress'
]

# WARNING! New arguments are needed! At least...box_size,rc,skin
	
def nodeGrid(n=None, box_size=0, rc=None, skin=None, eh_size=0, ratioMS=0, idealGas=0, slabMSDims=[0, 0, 0,]):
    print "################################################# Warning #####################################################"
    print "This Domain Decomposition algorithm requires minimally, the following arguments nodeGrid(n, box_size, rc, skin)"
    print "If you prefer to use the simple Domain Decomp. algorithm which may affect the performance of your MD simulation"
    print "then go for the function nodeGridSimple(n), which is also the default one if you give only one argument 'n'.   "
    print "Important: In case you are aiming to perform Lattice Boltzmann simulations you no need to worry about the per- "
    print "formance."
    print "For further details look into ESPResSo++ documentation or H.V. Guzman et. al, Phys. Rev. E, 96, 053311 (2017). "
    print "###############################################################################################################"
    if isinstance(n, numbers.Number) and isinstance(box_size, numbers.Number):					
        return nodeGridSimple(n)
    else:
	    ijkmax = 3 * n * n + 1
	    boxList = [box_size[0], box_size[1], box_size[2]]
	    ratioEH2CG = [(2. * eh_size) / (abs(box_size[0] - (2. * eh_size))), (2. * eh_size)/(abs(box_size[1] - (2. * eh_size))), (2. * eh_size) / (abs(box_size[2] - (2. * eh_size)))]
	    ratioEHCG = min(ratioEH2CG)
	    if not idealGas:  # Non ideal Gas Simulation
		# This condition checks if the sys is AdResS, if afirmative it resizes box accordingly
		if ratioMS > 1:  # Spatially heterogeneous simulation
		    # print ratioEHCG
		    # This condition checks if the sys is slab based or cubic and resizes the box accordingly
		    if sum(slabMSDims) > 0.1:
		        boxList = [slabMSDims[0] * (box_size[0] + (2. * eh_size) * (pow(ratioMS, 0.3333) - 1.)), slabMSDims[1] * (box_size[1] + (2. * eh_size) * (pow(ratioMS, 0.3333) - 1.)), slabMSDims[2] * (box_size[2] + (2. * eh_size) * (pow(ratioMS, 0.3333) - 1.))]
		        boxList = [(1 - slabMSDims[0]) * box_size[0] * (pow(ratioMS, 0.3333)) + boxList[0], (1 - slabMSDims[1]) * box_size[1] * (pow(ratioMS, 0.3333)) + boxList[1], (1 - slabMSDims[2]) * box_size[2] * (pow(ratioMS, 0.3333)) + boxList[2]]  # Flag NOR

		    else:
		        boxList = [box_size[0] + (2. * eh_size) * (ratioMS - 1.), box_size[1] + (
		            2. * eh_size) * (ratioMS - 1.), box_size[2] + (2. * eh_size) * (ratioMS - 1.)]
		else:
		    print "HeSpaDDA message: Non AdResS DD...entering Homogeneous DD...3, 2, 1 "
	    # Ideal Gas HeSpaDDA (idealGas==1)
	    else:
		if ratioMS > 1:
		    # This condition checks if the sys is slab based or cubic and resizes the box accordingly
		    if sum(slabMSDims) > 0:
		        if n <= round((2. * eh_size) / (rc + skin) - 0.5 + 2.) and min(boxList) * ratioEHCG / (rc + skin) >= 1.:
		            boxList = [slabMSDims[0] * ((2. * eh_size) + 2. * (rc + skin)), slabMSDims[1] * ((2. * eh_size) + 2. * (rc + skin)), slabMSDims[2] * ((2. * eh_size) + 2. * (rc + skin))]
		            boxList = [(1 - slabMSDims[0]) * box_size[0] + boxList[0], (1 - slabMSDims[1])* box_size[1] + boxList[1], (1 - slabMSDims[2]) * box_size[2] + boxList[2]]
		            print "HeSpaDDA message: this option ONE of the ideal gas size box your low resolution region is small"			
		        elif n > round((2. * eh_size) / (rc + skin) - 0.5 + 2.) and min(boxList) / (rc + skin) >= 4. and round((2. * eh_size) / (rc + skin) - 0.5 + 2.) > min(boxList) / (rc + skin):
		            boxList = [slabMSDims[0] * ((2. * eh_size) + 2. * (rc + skin)), slabMSDims[1] * ((2. * eh_size) + 2. * (rc + skin)), slabMSDims[2] * ((2. * eh_size) + 2. * (rc + skin))]
		            boxList = [(1 - slabMSDims[0]) * box_size[0] * (pow(ratioMS, 0.3333)) + boxList[0], (1 - slabMSDims[1]) * box_size[1] * (pow(ratioMS, 0.3333)) + boxList[1], (1 - slabMSDims[2]) * box_size[2] * (pow(ratioMS, 0.3333)) + boxList[2]]
		            # boxList=[box_size[0],box_size[1],box_size[2]]	# UJ is here
		            print "HeSpaDDA message: this option TWO of the ideal gas size box your low resolution region is big"
		        elif n > round((2. * eh_size) / (rc + skin) - 0.5 + 2.) and min(boxList) / (rc + skin) < 4.:
		            boxList = [slabMSDims[0] * ((2. * eh_size) + 2. * (rc + skin)), slabMSDims[1] * ((2. * eh_size) + 2. * (rc + skin)), slabMSDims[2] * ((2. * eh_size) + 2. * (rc + skin))]
		            boxList = [(1 - slabMSDims[0]) * box_size[0] + boxList[0], (1 - slabMSDims[1]) * box_size[1] + boxList[1], (1 - slabMSDims[2]) * box_size[2] + boxList[2]]
		            # boxList=[box_size[0],box_size[1],box_size[2]]	# UJ is here
		            print "HeSpaDDA message: this option THREE of the ideal gas size box your low resolution region is just sufficiently big"
		        else:
		            print "HeSpaDDA message: No more space for distributing Cores...look if you could use some HalfCells, or reduce the nr. of processors used and try again"
		    else:
		        boxList = [(2. * eh_size) + 2. * (rc + skin), (2. * eh_size) + 2. * (rc + skin), (2. * eh_size) + 2. * (rc + skin)]  # IDEA: multiply here by the cells refina
		else:
		    print "HeSpaDDA message: Non AdResS DD!"

	    # Here starts the new Dimensional aware HDD algorithm which. Dependencies: loadbal
	    LoN_Avgmin = sum(boxList)
	    # ima and imi sort values according to the box dimensions
	    ima = boxList.index(max(boxList))
	    imi = boxList.index(min(boxList))
	    dN = [1, 1, 1]
	    fdN = [0, 0, 0]
	    for i in xrange(1, n + 1):
		for j in xrange(i, n + 1):
		    for k in xrange(j, n + 1):
		        if (i * j * k == n) and (i * i + j * j + k * k < ijkmax):
		            dN[0] = k
		            dN[1] = j
		            dN[2] = i
		            ijkmax = i * i + j * j + k * k
		            # This checks if the system's box is cubic if not add weighted averages
		            # NOTE: It could be further optimized on a system by system basis as f(deltaCellSize beteween cores) this corresponds to the general version.
		            if qbicity(box_size, rc, skin) == False:
		                ndN = changeIndex(dN, ima, imi)[:]
		                LoN_norm = [boxList[0] / ndN[0], boxList[1] / ndN[1], boxList[2] / ndN[2]]
		                LoN_Avg = sum(LoN_norm) / 3.0
		                if LoN_Avg <= LoN_Avgmin:
		                    LoN_Avgmin = LoN_Avg
		                    fdN = ndN[:]
		                    ijkmax = fdN[0] * fdN[0] + fdN[1] * fdN[1] + fdN[2] * fdN[2]
		                    print fdN
		                else:
		                    ijkmax = fdN[0] * fdN[0] + fdN[1] * fdN[1] + fdN[2] * fdN[2]
		                    print 'HeSpaDDA message: No update of dN req ...'
		            else:
		                print 'Cubicity check passed -> powered by HeSpaDDA'
		                fdN = [k, j, i]
	    if abs(box_size[1] - box_size[2]) < (2 * (rc + skin)):
		if fdN[2] > fdN[1]:
		    aux = fdN[2]
		    fdN[2] = fdN[1]
		    fdN[1] = aux
		    # print 'ordered fdN:',fdN[0],fdN[1],fdN[2]
		else:
		    print 'HeSpaDDA message: Size Lenghts are eq. while ordering axis with preference on X, Y and Z!'
	    else:
		print 'HeSpaDDA message: Size Lenghts are different in Y and Z!'
	    return Int3D(fdN[0], fdN[1], fdN[2])


def cellGrid(box_size, node_grid, rc, skin):
    rc_skin = rc + skin
    if rc_skin == 0:
        raise Error("interaction range (cutoff + skin) must be larger than 0")
    if (node_grid[0] <= 0 or node_grid[1] <= 0 or node_grid[2] <= 0):
        raise Error("invalid node grid %s" % str(node_grid))
    ix = box_size[0] / (rc_skin * node_grid[0])

    if ix < 1:
        raise Error("local box size in direction 0 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime, perhaps you could also try with Halfcells." % (ix, rc_skin))
    iy = box_size[1] / (rc_skin * node_grid[1])
    if iy < 1:
        raise Error("local box size in direction 1 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime, perhaps you could also try with Halfcells." % (iy, rc_skin))
    iz = box_size[2] / (rc_skin * node_grid[2])
    if iz < 1:
        raise Error("local box size in direction 2 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime, perhaps you could also try with Halfcells." % (iz, rc_skin))

    return Int3D(ix, iy, iz)


def nodeGridSimple(n):
    ijkmax = 3 * n * n + 1
    d1 = 1
    d2 = 1
    d3 = 1
    for i in xrange(1, n + 1):
        for j in xrange(i, n + 1):
            for k in xrange(j, n + 1):
                if (i * j * k == n) and (i * i + j * j + k * k < ijkmax):
                    d1 = k
                    d2 = j
                    d3 = i
                    ijkmax = i * i + j * j + k * k
    return Int3D(d1, d2, d3)


def cherrypickTotalProcs(box_size, rc, skin, MnN, CpN, percTol=0.2, eh_size=0, ratioMS=0, idealGas=0, slabMSDims=[0, 0, 0]):
    indMax = max(box_size)
    indMin = min(box_size)
    L_l = box_size[indMax]
    L_s = box_size[indMin]
    L_lHR = 2. * eh_size
    L_lLR = L_l - L_lHR
    SymFactor = int((L_l + L_lHR * (ratioMS**(1. / 3.) - 1.)) / (L_s * (ratioMS**(1. / 3.) - 1.)))
    nX = [0] * MnN
    nY = [0] * MnN
    nZ = [0] * MnN
    procArray = range(1 * CpN, (MnN + 1) * CpN, CpN)
    print "HeSpaDDA message: Your search array in terms of total number of processors is:", procArray
    for i in xrange(1, MnN + 1):
        nX[i - 1], nY[i - 1], nZ[i - 1] = nodeGrid(procArray[i - 1], box_size, rc, skin, eh_size, ratioMS, idealGas, slabMSDims)
    print "HeSpaDDA message: For your Information, we are tackling the following number of processors"
    print "HeSpaDDA message: Processors in the x-axis are:", nX, "\nProcessors in the y-axis are:", nY, "\nProcessors in the z-axis are:", nZ
    pickedP = []
    for i in xrange(0, MnN):
        if (indMax == 0 and indMin == 1) or (indMax == 0 and indMin == 2):
            if nX[i] > (SymFactor - percTol * SymFactor) * nY[i] and nX[i] < (SymFactor + percTol * SymFactor) * nY[i]:
                pickedP.append(i)
        elif (indMax == 1 and indMin == 0) or (indMax == 1 and indMin == 2):
            if nY[i] > (SymFactor - percTol * SymFactor) * nX[i] and nY[i] < (SymFactor + percTol * SymFactor) * nX[i]:
                pickedP.append(i)
        elif (indMax == 2 and indMin == 0) or (indMax == 2 and indMin == 1):
            if nZ[i] > (SymFactor - percTol * SymFactor) * nX[i] and nZ[i] < (SymFactor + percTol * SymFactor) * nX[i]:
                pickedP.append(i)
    print "HeSpaDDA message: There are a couple of total number of processors that satisfy your requirements: ", [(v+1)*CpN for v in pickedP]
    return [(v+1)*CpN for v in pickedP]

# WARNING! This is a new function to find values for the new neighborList DataStruct for Hom-Sys


def neiListHom(node_grid, box, rc, skin):
    # data structure Initialization
    print "HeSpaDDA message: Current homogeneous NodeGrid (X,Y,Z)", node_grid
    rc_skin = rc + skin
    neiListxin = []
    neiListyin = []
    neiListzin = []
    cursor = [box[0], box[1], box[2]]
    # define Neighbor vectors
    neiListx = [0] * (node_grid[0] + 1)
    neiListy = [0] * (node_grid[1] + 1)
    neiListz = [0] * (node_grid[2] + 1)
    # It makes use of reDistCellsHom and adaptNeiList form loadbal to find homogenous DD
    neiListxin = reDistCellsHom(node_grid[0], cursor[0], rc_skin)
    neiListx = adaptNeiList(neiListxin)
    neiListyin = reDistCellsHom(node_grid[1], cursor[1], rc_skin)
    neiListy = adaptNeiList(neiListyin)
    neiListzin = reDistCellsHom(node_grid[2], cursor[2], rc_skin)
    neiListz = adaptNeiList(neiListzin)
    return map(int, neiListx), map(int, neiListy), map(int, neiListz)

# WARNING! This is a new function to find values for the new neighborList DataStruct for inHom-Sys

def neiListAdress(node_grid, cell_grid, rc, skin, eh_size, adrCenter, ratioMS, idealGasFlag=True, sphereAdr=False, slabMSDims=[1, 0, 0]):
    # dataStructure Initialization
    print "HeSpaDDA message: Current heterogeneous NodeGrid (X,Y,Z)", node_grid
    rc_skin = rc + skin
    # define Neighbor vectors
    neiListx = [0] * (node_grid[0] + 1)
    neiListy = [0] * (node_grid[1] + 1)
    neiListz = [0] * (node_grid[2] + 1)
    # Check adress regions
    adrCenter = [adrCenter[0], adrCenter[1], adrCenter[2]]
    # Midle point DD
    cursor = [adrCenter[0] * 2, adrCenter[1] * 2, adrCenter[2] * 2]
    # cg_sizeR=adrCenter[0]+eh_size 	# for further devs
    # cg_sizeL=adrCenter[0]-eh_size	# for further devs
    # round(cg_sizeL/rc_skin-0.5)+round((cursor[0]-cg_sizeR)/ rc_skin-0.5)+round((cg_sizeR-cg_sizeL)/rc_skin-0.5) # old implementation rounded # Cells too much!
    cellsX = round(cursor[0] / rc_skin - 0.5)
    cellsY = round(cursor[1] / rc_skin - 0.5)
    cellsZ = round(cursor[2] / rc_skin - 0.5)
    print "HeSpaDDA message: Current heterogeneous CellGrid (X,Y,Z)", cellsX, cellsY, cellsZ
    # This condition checks if the Sys is a Slab or a Sphere. It should work for any middle based Sys
    if not sphereAdr:
        if slabMSDims[0] == 1:
            halfneilListx = halfDecomp(adrCenter[0], rc_skin, eh_size, int(round(node_grid[0] / 2. - 0.5)), cellsX, ratioMS, cursor[0], idealGasFlag)
            print "HeSpaDDA message: My halfneilListx called a half decomposition of the first half of the sim box as (algorithm S.1)...,", halfneilListx
            # Next instruction is doubling(by unfolding halfDecomp) the halfspace based DD
            neiListxin = addHsymmetry(halfneilListx, eh_size, rc_skin, node_grid[0], cellsX, ratioMS, cursor[0], idealGasFlag)
            print "HeSpaDDA message: My neiListxin if it make sense your heterogenous system can be splitted in 2 equal parts (algorithm S.2)...,", neiListxin
        else:
            neiListxin = reDistCellsHom(node_grid[0], cursor[0], rc_skin)
        # Contains cores Neighbor List in full format for X
        neiListx = adaptNeiList(neiListxin)
        if slabMSDims[1] == 1:
            halfneilListy = halfDecomp(adrCenter[1], rc_skin, eh_size, int(round(node_grid[1] / 2. - 0.5)), cellsY, ratioMS, cursor[1], idealGasFlag)
            print "HeSpaDDA message: My halfneilListx called a half decomposition of the first half of the sim box as (algorithm S.1)...,", halfneilListy
            # Next instruction is doubling(by unfolding halfDecomp) the halfspace based DD
            neiListyin = addHsymmetry(halfneilListy, eh_size, rc_skin, node_grid[1], cellsY, ratioMS, cursor[1], idealGasFlag)
            print "HeSpaDDA message: My neiListxin if it make sense your heterogenous system can be splitted in 2 equal parts (algorithm S.2)...,", neiListyin
        else:
            neiListyin = reDistCellsHom(node_grid[1], cursor[1], rc_skin)
        # Contains the homogeneously decomp cores Neighbor List in full format for Y
        neiListy = adaptNeiList(neiListyin)
        if slabMSDims[2] == 1:
            halfneilListz = halfDecomp(adrCenter[2], rc_skin, eh_size, int(round(node_grid[2] / 2. - 0.5)), cellsZ, ratioMS, cursor[2], idealGasFlag)
            print "HeSpaDDA message: My halfneilListx called a half decomposition of the first half of the sim box as (algorithm S.1)...,", halfneilListz
            # Next instruction is doubling(by unfolding halfDecomp) the halfspace based DD
            neiListzin = addHsymmetry(halfneilListz, eh_size, rc_skin, node_grid[2], cellsZ, ratioMS, cursor[2], idealGasFlag)
            print "HeSpaDDA message: My neiListxin if it make sense your heterogenous system can be splitted in 2 equal parts (algorithm S.2)...,", neiListzin
        else:
            neiListzin = reDistCellsHom(node_grid[2], cursor[2], rc_skin)
        # Contains the homogeneously decomp cores Neighbor List in full format for Z
        neiListz = adaptNeiList(neiListzin)
        # NOTE that additional DD options for slabs can be furhter added
    elif sphereAdr:
        flx, fly, flz = nodeGridSizeCheck(node_grid[0], node_grid[1], node_grid[2])
        print "HeSpaDDA message: Is it worthy to use an advanced DD? Flags for it in each dir X,Y,Z (0 means OK for heterogenoeus DD)...,", flx, fly, flz
        # on X-axis
        if flx == 0:
            halfneilListx = halfDecomp(adrCenter[0], rc_skin, eh_size, int(round(node_grid[0] / 2. - 0.5)), cellsX, ratioMS, cursor[0], idealGasFlag)
            print "HeSpaDDA message: My halfneilListx called a half decomposition of the first half of the sim box as (algorithm S.1)...,", halfneilListx
            neiListxin = addHsymmetry(halfneilListx, eh_size, rc_skin, node_grid[0], cellsX, ratioMS, cursor[0], idealGasFlag)
            print "HeSpaDDA message: My neiListxin if it make sense your heterogenous system can be splitted in 2 equal parts (algorithm S.2)...,", neiListxin
        elif flx > 1:
            neiListxin = reDistCellsHom(node_grid[0], cursor[0], rc_skin)
        neiListx = adaptNeiList(neiListxin)
        # on Y-axis
        if fly == 0:
            halfneilListy = halfDecomp(adrCenter[1], rc_skin, eh_size, int(round(node_grid[1] / 2. - 0.5)), cellsY, ratioMS, cursor[1], idealGasFlag)
            print "HeSpaDDA message: My halfneilListx called a half decomposition of the first half of the sim box as (algorithm S.1)...,", halfneilListy
            neiListyin = addHsymmetry(halfneilListy, eh_size, rc_skin, node_grid[1], cellsY, ratioMS, cursor[1], idealGasFlag)
            print "HeSpaDDA message: My neiListxin if it make sense your heterogenous system can be splitted in 2 equal parts (algorithm S.2)...,", neiListyin
        elif fly > 1:
            neiListyin = reDistCellsHom(node_grid[1], cursor[1], rc_skin)
        neiListy = adaptNeiList(neiListyin)
        # on Z-axis
        if flz == 0:
            halfneilListz = halfDecomp(adrCenter[2], rc_skin, eh_size, int(round(node_grid[2] / 2. - 0.5)), cellsZ, ratioMS, cursor[2], idealGasFlag)
            print "HeSpaDDA message: My halfneilListx called a half decomposition of the first half of the sim box as (algorithm S.1)...,", halfneilListz
            neiListzin = addHsymmetry(halfneilListz, eh_size, rc_skin, node_grid[2], cellsZ, ratioMS, cursor[2], idealGasFlag)
            print "HeSpaDDA message: My neiListxin if it make sense your heterogenous system can be splitted in 2 equal parts (algorithm S.2)...,", neiListzin
        elif flz > 1:
            neiListzin = reDistCellsHom(node_grid[2], cursor[2], rc_skin)
        neiListz = adaptNeiList(neiListzin)
    print "HeSpaDDA message: neiListX:", map(int, neiListx), "\n neiListY", map(int, neiListy), "\n neiListZ", map(int, neiListz)
    return Int3D(cellsX, cellsY, cellsZ), map(int, neiListx), map(int, neiListy), map(int, neiListz)


def tuneSkin(system, integrator, minSkin=0.01, maxSkin=1.5, precision=0.001, printInfo=True):
    if printInfo:
        print 'HeSpaDDA message: The tuning is started. It can take some time depending on your system.'

    fi = (1.0 + math.sqrt(5.0)) / 2.0  # golden ratio

    npart = espressopp.analysis.NPart(system).compute()

    # this is an empirical formula in order to get the appropriate number of steps
    nsteps = int(espressopp.MPI.COMM_WORLD.size * 1000000.0 / float(npart))

    if printInfo:
        print 'CellGrid before tuning: ', system.storage.getCellGrid()
        sys.stdout.write('\nSteps     = %d\n' % nsteps)
        sys.stdout.write('Precision = %g\n' % precision)
        sys.stdout.write('It runs till deltaSkin<precision\n')

    if printInfo:
        prnt_format1 = '\n%9s %10s %10s %10s %14s\n'
        sys.stdout.write(prnt_format1 % ('time1: ', ' time2: ', ' skin1: ', ' skin2: ', ' deltaSkin: '))

    while (maxSkin - minSkin >= precision):
        skin1 = maxSkin - (maxSkin - minSkin) / fi
        skin2 = minSkin + (maxSkin - minSkin) / fi

        system.skin = skin1
        system.storage.cellAdjust()
        start_time = time.time()
        integrator.run(nsteps)
        end_time = time.time()
        time1 = end_time - start_time

        system.skin = skin2
        system.storage.cellAdjust()
        start_time = time.time()
        integrator.run(nsteps)
        end_time = time.time()
        time2 = end_time - start_time

        if(time1 > time2):
            minSkin = skin1
        else:
            maxSkin = skin2

        if printInfo:
            prnt_format2 = '%7.3f %10.3f %11.4f %10.4f %12.6f\n'
            sys.stdout.write(prnt_format2 % (time1, time2, minSkin, maxSkin, (maxSkin - minSkin)))

            sys.stdout.write('\nNew skin: %g\n' % system.skin)
            sys.stdout.write('\nNew cell grid: %s\n' % system.storage.getCellGrid())

    system.skin = (maxSkin + minSkin) / 2.0
    system.storage.cellAdjust()

    return (maxSkin + minSkin) / 2.0


def printTimeVsSkin(system, integrator, minSkin=0.01, maxSkin=1.5, skinStep=0.005):
    npart = espressopp.analysis.NPart(system).compute()
    # this is an empirical formula in order to get the appropriate number of steps
    nsteps = int(espressopp.MPI.COMM_WORLD.size * 10000000.0 / float(npart))

    print '      Calculations is started. It will print out the dependece of time of \n\
      running of %d steps on the skin size into the file \'timeVSskin.dat\'.\n\
      The range of skin sizes is [%g, %g], skin step is %g. It can take some \n\
      time depending on your system.' % (nsteps, minSkin, maxSkin, skinStep)

    curSkin = minSkin

    fmt2 = ' %8.4f %8.4f\n'  # format for writing the data
    nameFile = 'timeVSskin.dat'
    resFile = open(nameFile, 'w')

    count = 0
    while (curSkin < maxSkin):
        system.skin = curSkin
        system.storage.cellAdjust()
        start_time = time.time()
        integrator.run(nsteps)
        end_time = time.time()
        time1 = end_time - start_time

        resFile.write(fmt2 % (system.skin, time1))

        count = count + 1
        if (count == 20):
            print 'skin: ', system.skin
            count = 0

        curSkin = curSkin + skinStep

    resFile.close()

    return
