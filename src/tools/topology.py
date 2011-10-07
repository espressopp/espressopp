import math
import random
from espresso import Real3D
from espresso.Exceptions import Error

def polymerRW(pid, startpos, numberOfMonomers, bondlength):
	x         = startpos[0]
	y         = startpos[1]
	z         = startpos[2]
	positions = [ Real3D(x, y, z) ]
	bonds     = []
	for i in range(numberOfMonomers-1):
		nextZ = (2.0*random.uniform(0,1)-1.0)*bondlength;
		rr    = math.sqrt(bondlength*bondlength-nextZ*nextZ);
		phi   = 2.0*math.pi*random.uniform(0,1);
		nextX = rr*math.cos(phi);
		nextY = rr*math.sin(phi);				

		x += nextX
		y += nextY
		z += nextZ
		# update monomer list:
		positions.append(Real3D(x, y, z))
		# update bond list:
		bonds.append((pid+i,pid+i+1))
		
	return positions, bonds
	

"""
def polymerSAW(pid, startpos, numberOfMonomers, bondlength, excludedVolumeRadius, partlist, maxtries=100):

	for i in range(maxtries):
		overlap = False
		x=[startpos[0]]
		y=[startpos[1]]
		z=[startpos[2]]
		bonds = []
		for j in range(numberOfMonomers):
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
				for k in range(j):
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
	for i in range(numberOfArms):
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
	for i in range(numberOfArms):
		addBonds, addX, addY, addZ = polymerSAW(pid+i*(numberOfMonomers-1), startpos, numberOfMonomers, bondlength, excludedVolumeRadius, partlist, maxtries=100)
		addBonds = [(pid,addBonds[0][1])] + addBonds[1:-1] + [addBonds[-1]]
		bonds += addBonds
		x += addX[1:-1] + [addX[-1]]
		y += addY[1:-1] + [addY[-1]]
		z += addZ[1:-1] + [addZ[-1]]
		for j in range(len(x)-1):
			partlist += [ [(x[j+1],y[j+1],z[j+1]), excludedVolumeRadius] ]
	return bonds, x, y, z
	
	
"""