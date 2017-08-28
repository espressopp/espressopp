#  Copyright (C) 2012,2013,2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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


from math import sqrt, pi, cos, sin
import random
from espressopp import Real3D
from espressopp.Exceptions import Error

def polymerRW(pid, startpos, numberOfMonomers, bondlength, return_angles=False, return_dihedrals=False, mindist=None, rng=None):
	"""
	Initializes polymers through random walk 
	"""
	
	x         = startpos[0]
	y         = startpos[1]
	z         = startpos[2]
	positions = [ Real3D(x, y, z) ]
	bonds     = []
	avecostheta = 0.0
	if return_angles == True:
	   angles    = []
	if return_dihedrals == True:
	   dihedrals    = []
	for i in xrange(numberOfMonomers-1):
	  if mindist and i > 0:
		while True:
		  if rng==None:
		    nextZ = (2.0*random.uniform(0,1)-1.0)*bondlength;
  		    phi   = 2.0*pi*random.uniform(0,1);
		  else:
		    nextZ = (2.0*rng()-1.0)*bondlength;
  		    phi   = 2.0*pi*rng();
		  rr    = sqrt(bondlength*bondlength-nextZ*nextZ);
		  nextX = rr*cos(phi);
		  nextY = rr*sin(phi);				
		 
		  ax    = positions[i][0] - positions[i-1][0]
		  ay    = positions[i][1] - positions[i-1][1]
		  az    = positions[i][2] - positions[i-1][2]
		  la    = sqrt(ax*ax + ay*ay + az*az)
		  

		  bx    = - nextX
		  by    = - nextY
		  bz    = - nextZ
		  lb    = sqrt(bx*bx + by*by + bz*bz)

	  
		  cx    = ax - bx
		  cy    = ay - by
		  cz    = az - bz
		  lc    = sqrt(cx*cx + cy*cy + cz*cz)
  
		  if lc > mindist:
		  	avecostheta += - (ax*bx + ay*by + az*bz) / (la * lb)
		  	#print "cos theta:", (ax*bx + ay*by + az*bz) / (la * lb)
		  	break
		  
		  	
	  else:
		if rng==None:
		  nextZ = (2.0*random.uniform(0,1)-1.0)*bondlength;
  		  phi   = 2.0*pi*random.uniform(0,1);
		else:
		  nextZ = (2.0*rng()-1.0)*bondlength;
  		  phi   = 2.0*pi*rng();
		rr    = sqrt(bondlength*bondlength-nextZ*nextZ);
		nextX = rr*cos(phi);
		nextY = rr*sin(phi);				

 	  x += nextX
	  y += nextY
	  z += nextZ
	  # update monomer list:
	  positions.append(Real3D(x, y, z))
	  # update bond list:
	  bonds.append((pid+i,pid+i+1))

	  if return_angles == True:
	    if i < numberOfMonomers-2:
		  angles.append((pid+i, pid+i+1, pid+i+2))

	  if return_dihedrals == True:
	    if i < numberOfMonomers-3:
		  dihedrals.append((pid+i, pid+i+1, pid+i+2, pid+i+3))

		  
	if mindist:	  
	  avecostheta /= (numberOfMonomers-2)
    
	if return_angles == True:

		if return_dihedrals == True:

			if mindist:	
				return positions, bonds, angles, dihedrals , avecostheta
			else:
				return positions, bonds, angles, dihedrals

		else:

			if mindist:	
				return positions, bonds, angles, avecostheta
			else:
				return positions, bonds, angles

	else:

		if return_dihedrals == True:

			if mindist:
				return positions, bonds, dihedrals, avecostheta
			else:	
				return positions, bonds, dihedrals
	
		else: 

			if mindist:
				return positions, bonds, avecostheta
			else:	
				return positions, bonds


"""
def polymerSAW(pid, startpos, numberOfMonomers, bondlength, excludedVolumeRadius, partlist, maxtries=100):

	for i in xrange(maxtries):
		overlap = False
		x=[startpos[0]]
		y=[startpos[1]]
		z=[startpos[2]]
		bonds = []
		for j in xrange(numberOfMonomers):
			if overlap:
				break

			if j > 0:
				# compute next monomer:
				nextZ     = (2.0*random.uniform(0,1)-1.0)*bondlength;
				rr     = math.sqrt(bondlength*bondlength-nextZ*nextZ);
				phi    = 2.0*math.pi*random.uniform(0,1);
				nextX      = rr*math.cos(phi);
				nextY      = rr*math.sin(phi);
			
				# update monomer list:
				x += [x[-1]+nextX]
				y += [y[-1]+nextY]
				z += [z[-1]+nextZ]
				# update bond list:
				bonds += [(pid+j,pid+j+1)]
				for k in xrange(j):
					# compute distance to existing particles of the same polymer:
					dist = math.sqrt((x[j]-x[k])**2 + (y[j]-y[k])**2 + (z[j]-z[k])**2)
					if dist < radius:
						overlap = True
						break
			if overlap:
				break	
				
			for part in partlist:
				# compute distance to existing particles:
				radius = max([excludedVolumeRadius,part[1]])
				dist = math.sqrt((x[j]-part[0][0])**2 + (y[j]-part[0][1])**2 + (z[j]-part[0][2])**2)
				if dist < radius:
					if j == 0:
						raise Error("Startpos is in excluded volume of existing particles.")
					else:
						overlap = True
						break

		
	if overlap:	
		raise Error("Number of tries exceeded maxtries.")
	else:
		return bonds, x, y, z  
	
def starPolymerRW(pid, startpos, bondlength, numberOfArms, numberOfMonomers):
	bonds = []
	x=[startpos[0]]
	y=[startpos[1]]
	z=[startpos[2]]
	for i in xrange(numberOfArms):
		addBonds, addX, addY, addZ = polymerRW(pid+i*(numberOfMonomers-1), startpos, numberOfMonomers, bondlength)
		addBonds = [(pid,addBonds[0][1])] + addBonds[1:-1] + [addBonds[-1]]
		bonds += addBonds
		x += addX[1:-1] + [addX[-1]]
		y += addY[1:-1] + [addY[-1]]
		z += addZ[1:-1] + [addZ[-1]]
	return bonds, x, y, z
		
					
def starPolymerSAW(pid, startpos, bondlength, numberOfArms, numberOfMonomers, excludedVolumeRadius, partlist, maxtries=100):
	bonds = []
	x=[startpos[0]]
	y=[startpos[1]]
	z=[startpos[2]]
	for i in xrange(numberOfArms):
		addBonds, addX, addY, addZ = polymerSAW(pid+i*(numberOfMonomers-1), startpos, numberOfMonomers, bondlength, excludedVolumeRadius, partlist, maxtries=100)
		addBonds = [(pid,addBonds[0][1])] + addBonds[1:-1] + [addBonds[-1]]
		bonds += addBonds
		x += addX[1:-1] + [addX[-1]]
		y += addY[1:-1] + [addY[-1]]
		z += addZ[1:-1] + [addZ[-1]]
		for j in xrange(len(x)-1):
			partlist += [ [(x[j+1],y[j+1],z[j+1]), excludedVolumeRadius] ]
	return bonds, x, y, z
	
	
"""
