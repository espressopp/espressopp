#!/usr/bin/env python2 
#  Copyright (C) 2016-2017(H)
#      Max Planck Institute for Polymer Research
#
#  This file is part of ESPResSo++.
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

import math

####################################################################
# Script for the preparation of a crystal of tetrahedral molecules #
####################################################################

# Bond length between atoms in tetrahedral molecules
bond = 0.95
bond = bond/(math.sqrt(2.0))

# Helper functions for the insertion of the atoms of a single molecule
def tetramolx(Mx):
	x = [Mx-bond/2.0, Mx+bond/2.0, Mx+bond/2.0, Mx-bond/2.0]
	return x

def tetramoly(My):
	y = [My-bond/2.0, My+bond/2.0, My-bond/2.0, My+bond/2.0]
	return y

def tetramolz(Mz):
	z = [Mz-bond/2.0, Mz-bond/2.0, Mz+bond/2.0, Mz+bond/2.0]
	return z

# We assume a box of the following size (Bx=36.845, By=15.0, Bz=15.0) and want to insert exacty 798 molecules. We use cells
# of the size (Cx=1.5, Cy=1.5, Cz=1.5) and insert molecules in their centers. We have exactly 2400 (24 in x direction, 10 in
# y and z direction respectively) cells but only insert in 35% of the cells. This gives use 840 filled cells. However we skip
# a few more insertions and later delete a few molecules to get exactly 798 molecules. We chose the box size and molecule number
# according to the H-AdResS paper (Potestio et. al., Phys. Rev. Let. 110, 108301 (2013)) to reproduce the results.

def makebonds(N):
	bonds = []
	for i in range(N):
		if i%4==0:
			bonds.append((i,i+1))
			bonds.append((i,i+2))
			bonds.append((i,i+3))
			bonds.append((i+1,i+2))
			bonds.append((i+1,i+3))
			bonds.append((i+2,i+3))	

	# Return bondlist
	return bonds

def crystal():
	x, y, z = [], [], []
	for i in range(24): # x
		for j in range(10): # y
			for k in range(10): # z
				if i == 7: # Skip a few insertions to get 805 molecules after insertion
					continue
				if j % 2 == 0:
					if k % 3 == 0:
						x += tetramolx(0.75+i*1.5)
						y += tetramoly(0.75+j*1.5)
						z += tetramolz(0.75+k*1.5)
				else:
					if k % 3 == 2:
						x += tetramolx(0.75+i*1.5)
						y += tetramoly(0.75+j*1.5)
						z += tetramolz(0.75+k*1.5)

        # Delete a few more molecules to get exactly 798 molecules					  
	for i in range(28):
		del x[-1]
		del y[-1]
		del z[-1]
	
	# Return coordinate lists
	return x,y,z





