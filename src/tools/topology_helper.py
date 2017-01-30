#  Copyright (C) 2015
#      Jakub Krajniak (jkrajniak at gmail.com)
#  Copyright (C) 2012,2013
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


# Some helper classes usefull when parsing the gromacs topology

"""
***************
topology_helper
***************
"""

import espressopp
import math
import gromacs
import os

class FileBuffer():
    def __init__(self):
        self.linecount=0
        self.lines=[]
        self.pos=0
    def appendline(self, line):
        self.lines.append(line)
    def readline(self):
        try:
            line=self.lines[self.pos]
        except:
            return ''
        self.pos+=1
        return line
    def readlastline(self):
        try:
            line=self.lines[self.pos-1]
        except:
            return ''
        return line
    def seek(self, p):
	self.pos=p
    def tell(self):
	return self.pos
        

def FillFileBuffer(fname, filebuffer):
    f=open(fname, 'r')
    for line in f:
        if "include" in line and not line[0]==';':
            name=(line.split()[1]).strip('\"')
            try: 
                FillFileBuffer(name, filebuffer)
            except IOError:
                #need to use relative path
                name = os.path.join(os.path.dirname(fname), name)
                FillFileBuffer(name, filebuffer)
        else:
            l=line.rstrip('\n')
            if l:
                filebuffer.appendline(l)
            
    f.close
    return


def FindType(proposedtype, typelist):
    list=[typeid for (typeid,atype) in typelist.iteritems() if atype==proposedtype ]
    if len(list)>1:
        print "Error: duplicate type definitons", proposedtype.parameters
        exit()
    elif len(list)==0:
        return None
    return list[0]
    

class InteractionType:
    def __init__(self, parameters):
        self.parameters=parameters
    def __eq__(self,other):
        # interaction types are defined to be equal if all parameters are equal
        for k, v in self.parameters.iteritems():
            if k not in other.parameters: return False
            if other.parameters[k]!=v: return False
        return True
    def createEspressoInteraction(self, system, fpl):
        print "WARNING: could not set up interaction for", self.parameters, ": Espresso potential not implemented"
        return None
    def automaticExclusion(self):
        #overwrite in derrived class if the particular interaction is automatically excluded
        return False

class HarmonicBondedInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        # interaction specific stuff here
        # spring constant kb is half the gromacs spring constant
        pot = espressopp.interaction.Harmonic(self.parameters['kb']/2.0, self.parameters['b0'])
        interb = espressopp.interaction.FixedPairListHarmonic(system, fpl, pot)
        return interb
    def automaticExclusion(self):
        return True
    
class MorseBondedInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        # interaction specific stuff here
        pot = espressopp.interaction.Morse(self.parameters['D'], self.parameters['beta'], self.parameters['rmin'])
        interb = espressopp.interaction.FixedPairListMorse(system, fpl, pot)
        return interb
    def automaticExclusion(self):
        return True
    
class FENEBondedInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        # interaction specific stuff here
        # spring constant kb is half the gromacs spring constant
        pot = espressopp.interaction.Fene(self.parameters['kb']/2.0, self.parameters['b0'])
        interb = espressopp.interaction.FixedPairListFene(system, fpl, pot)
        return interb
    def automaticExclusion(self):
        return True
    
class HarmonicAngleInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        # interaction specific stuff here
        # spring constant kb is half the gromacs spring constant. Also convert deg to rad
        pot = espressopp.interaction.AngularHarmonic(self.parameters['k']/2.0, self.parameters['theta']*2*math.pi/360)
        interb = espressopp.interaction.FixedTripleListAngularHarmonic(system, fpl, pot)
        return interb     
    

class TabulatedBondInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        spline = 3
        fg = "table_b"+str(self.parameters['tablenr'])+".xvg"
        fe = fg.split(".")[0]+".tab" # name of espressopp file
        gromacs.convertTable(fg, fe)
        potTab = espressopp.interaction.Tabulated(itype=spline, filename=fe)
        interb = espressopp.interaction.FixedPairListTabulated(system, fpl, potTab)
        return interb
    def automaticExclusion(self):
        return self.parameters['excl']
    
class TabulatedAngleInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        spline = 3
        fg = "table_a"+str(self.parameters['tablenr'])+".xvg"
        fe = fg.split(".")[0]+".tab" # name of espressopp file
        gromacs.convertTable(fg, fe)
        potTab = espressopp.interaction.TabulatedAngular(itype=spline, filename=fe)
        interb = espressopp.interaction.FixedTripleListTabulatedAngular(system, fpl, potTab)
        return interb  

class TabulatedDihedralInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        spline = 3
        fg = "table_d"+str(self.parameters['tablenr'])+".xvg"
        fe = fg.split(".")[0]+".tab" # name of espressopp file
        gromacs.convertTable(fg, fe)
        potTab = espressopp.interaction.TabulatedDihedral(itype=spline, filename=fe)
        interb = espressopp.interaction.FixedQuadrupleListTabulatedDihedral(system, fpl, potTab)
        return interb       

class HarmonicNCosDihedralInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        # DihedralHarmonicNCos coded such that k = gromacs spring constant. Convert degrees to rad 
        pot = espressopp.interaction.DihedralHarmonicNCos(self.parameters['K'], self.parameters['phi0']*2*math.pi/360, self.parameters['multiplicity'])
        interb = espressopp.interaction.FixedQuadrupleListDihedralHarmonicNCos(system, fpl, pot)
        return interb
    def automaticExclusion(self):
        return True
    
class RyckaertBellemansDihedralInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        print('RyckaertBellemans: {}'.format(self.parameters))
        pot = espressopp.interaction.DihedralRB(**self.parameters)
        return espressopp.interaction.FixedQuadrupleListDihedralRB(system, fpl, pot)

class HarmonicDihedralInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl):
        #print('RyckaertBellemans: {}'.format(self.parameters))
        pot = espressopp.interaction.DihedralHarmonic(self.parameters['K'], self.parameters['phi0']*2*math.pi/360)
        return espressopp.interaction.FixedQuadrupleListDihedralHarmonic(system, fpl, pot)

def ParseBondTypeParam(line):
    tmp = line.split() 
    btype= tmp[2]
    # TODO: handle exclusions automatically
    if btype == "8":
        p=TabulatedBondInteractionType({"tablenr":int(tmp[3]),"k":float(tmp[4]), 'excl':True})
    elif btype == "9":
        p=TabulatedBondInteractionType({"tablenr":int(tmp[3]), "k":float(tmp[4]), 'excl':False})
    elif btype == "1":
        p=HarmonicBondedInteractionType({"b0":float(tmp[3]), "kb":float(tmp[4])})
    elif btype == "3":
        p=MorseBondedInteractionType({"b0":float(tmp[3]), "D":float(tmp[4]), "beta":float(tmp[5])})
    elif btype == "7":
        p=FENEBondedInteractionType({"b0":float(tmp[3]), "kb":float(tmp[4])})
    elif btype == "9":
        p=TabulatedBondInteractionType({"tablenr":int(tmp[3]), "k":float(tmp[4])})
    else:
        print "Unsupported bond type", tmp[2], "in line:"
        print line
        exit()
    return p     

def ParseAngleTypeParam(line):
    tmp = line.split() 
    type= int(tmp[3])
    if type == 1:
        p=HarmonicAngleInteractionType({"theta":float(tmp[4]), "k":float(tmp[5])})
    elif type == 8:
        p=TabulatedAngleInteractionType({"tablenr":int(tmp[4]),"k":float(tmp[5])})
    else:
        print "Unsupported angle type", type, "in line:"
        print line
        exit()
    return p    

def ParseDihedralTypeParam(line):
    tmp = line.split() 
    type= int(tmp[4])
    if type == 8:
        p=TabulatedDihedralInteractionType({"tablenr":int(tmp[5]), "k":float(tmp[6])})
    elif type == 3:
        tmp[5:11] = map(float, tmp[5:11])
        p = RyckaertBellemansDihedralInteractionType(
            {'K0': tmp[5], 'K1': tmp[6], 'K2': tmp[7], 'K3': tmp[8], 'K4': tmp[9], 'K5': tmp[10]}
        )
    elif (type == 1) or (type == 9): 
        p=HarmonicNCosDihedralInteractionType({"K":float(tmp[6]), "phi0":float(tmp[5]), "multiplicity":int(tmp[7])})
    else:
        print "Unsupported dihedral type", type, "in line:"
        print line
        exit()
    return p    

def ParseImproperTypeParam(line):
    tmp = line.split()
    type= int(tmp[4])
    if type == 4:
        p=HarmonicNCosDihedralInteractionType({"K":float(tmp[6]), "phi0":float(tmp[5]), "multiplicity":int(tmp[7])})
    elif type == 2:
        p=HarmonicDihedralInteractionType({"K":float(tmp[6]), "phi0":float(tmp[5])})
    else:
        print "Unsupported improper type", type, "in line:"
        print line
        exit()
    return p

# Usefull code for generating the regular exclusions

class Node():
    def __init__(self, id):
	self.id=id
	self.neighbours=[]
    def addNeighbour(self, nb):
	self.neighbours.append(nb)

def FindNodeById(id, nodes):
    list=[n for n in nodes if n.id==id ]
    if len(list)>1:
        print "Error: duplicate nodes", id
        exit()
    elif len(list)==0:
        return None
    return list[0]

def FindNNextNeighbours(startnode, numberNeighbours, neighbours, forbiddenNodes):
    if numberNeighbours==0:
	return neighbours
	
    #avoid going back the same path
    forbiddenNodes.append(startnode)
    
    # Loop over next neighbours and add them to the neighbours list
    # Recursively call the function with numberNeighbours-1
    for n in startnode.neighbours:
	if not n in forbiddenNodes:
	    if n not in neighbours: neighbours.append(n) # avoid double counting in rings
	    FindNNextNeighbours(n, numberNeighbours-1, neighbours, forbiddenNodes)
            
 
def GenerateRegularExclusions(bonds, nrexcl, exclusions):
    nodes=[]
    # make a Node object for each atom involved in bonds
    for b in bonds:
        bids=b[0:2]
        for i in bids:
            if FindNodeById(i, nodes)==None:
               n=Node(i)
               nodes.append(n)

    # find the next neighbours for each node and append them   
    for b in bonds:
        permutations=[(b[0], b[1]), (b[1], b[0])]
        for p in permutations:
            n=FindNodeById(p[0], nodes)
            nn=FindNodeById(p[1], nodes)
            n.addNeighbour(nn)   
    
    # for each atom, call the FindNNextNeighbours function, which recursively
    # seraches for nrexcl next neighbours
    for n in nodes:
        neighbours=[]
        FindNNextNeighbours(n, nrexcl, neighbours, forbiddenNodes=[])
        for nb in neighbours:
            # check if the permutation is already in the exclusion list
            # this may be slow, but to do it in every MD step is even slower...
            # TODO: find a clever algorithm which does avoid permuations from the start
            if not (n.id, nb.id) in exclusions:
                if not (nb.id, n.id) in exclusions:
                    exclusions.append((n.id, nb.id))
   
    return exclusions
