#  Copyright (C) 2020(H)
#      Jozef Stefan Institute
#      Max Planck Institute for Polymer Research
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

r"""

*******************************
PDB - read and write pdb format
*******************************

.. function:: espressopp.tools.pdbwrite(filename, system, molsize=4, append=False, typenames=None)

  Writes a file in PDB format

  :param filename: output file name
  :type filename: string
  :param system: espressopp system
  :type system: espressopp System object
  :param molsize: if molsize>0, the molecule count is increased every molsize particles (default 4)
  :type molsize: int
  :param append: if True, append to filename, other over-write filename (default False)
  :type append: bool
  :param typenames: dictionary used for mapping from espressopp's integer particle types to the particle type strings written in a pdb file
  :type typenames: dict, key=int, value=string

.. function:: espressopp.tools.pdbread(filename,natoms,header)

  Reads one frame of a pdb format file

  :param filename: input file name
  :type filename: string
  :param natoms: number of atoms in pdf file
  :type natoms: int
  :param header: number of header lines to skip at start of file
  :type header: int

  Returns: index,atomname,locator,resname,resid,resseq,x,y,z,alpha,beta,segid,element (lists of type int,str,str,str,str,int,float,float,float,float,float,str,str)

"""
import espressopp
from math import sqrt
from espressopp import Real3D

def pdbwrite(filename, system, molsize=4, append=False, typenames=None):
    #typenames: a map of typeid to typename to be written, e.g. for water typenames={0:'H', 1:'O'}
    if append:
        file = open(filename, 'a')
        s = "\n"
    else:
        file = open(filename,'w')
        s = "REMARK generated by ESPResSo++\n"
    file.write(s)
    maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
    pid    = 0
    addToPid = 0 # if pid begins from 0, then addToPid should be +1
    mol    = 0
    molcnt = 0
    name='FE' # default name, overwritten when typenames map is given
    #following http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    #crystal header
    st = "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n"%('CRYST1',system.bc.boxL[0],system.bc.boxL[1],system.bc.boxL[2],90.00,90.00,90,'P 1',1) #boxes are orthorhombic for now
    file.write(st)
    while pid <= maxParticleID:
        if system.storage.particleExists(pid):
            particle = system.storage.getParticle(pid)
            if(pid==0):
                addToPid = 1
            xpos   = particle.pos[0]
            ypos   = particle.pos[1]
            zpos   = particle.pos[2]
            type   = particle.type

            if typenames:
                name=typenames[type]
            #st = "ATOM %6d  FE  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d\n"%(pid, mol, xpos, ypos, zpos, type)
            #following http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
            st = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f     T%04d%2s%2s\n"%('ATOM  ',(pid+addToPid)%100000,name,'','UNX','F',mol%10000,'',xpos,ypos,zpos,0,0,mol%10000,'','')#the additional 'T' in the string is needed to be recognized as string,%10000 to obey the fixed-width format
            file.write(st)
            pid    += 1
            molcnt += 1
            if molcnt == molsize:
                mol   += 1
                molcnt = 0
        else:
            pid   += 1

    file.write('END\n')
    file.close()

def fastwritepdb(filename, system, molsize=1000, append=False, folded=True):
    if append:
        file = open(filename, 'a')
        s = "\n"
    else:
        file = open(filename,'w')
        s = "REMARK generated by ESPResSo++\n"
    file.write(s)
    mol    = 0
    molcnt = 0
    configurations = espressopp.analysis.Configurations(system)
    configurations.gather()
    configuration = configurations[0]
    box_x = system.bc.boxL[0]
    box_y = system.bc.boxL[1]
    box_z = system.bc.boxL[2]
    #following http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    #crystal header
    st = "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n"%('CRYST1',box_x, box_y, box_z, 90.00, 90.00, 90, 'P 1', 1) #boxes are orthorhombic for now
    file.write(st)
    for pid in configuration:
        if folded:
            pos       = espressopp.Real3D(configuration[pid][0], configuration[pid][1], configuration[pid][2])
            foldedpos = system.bc.getFoldedPosition(pos)
            xpos      = foldedpos[0][0]
            ypos      = foldedpos[0][1]
            zpos      = foldedpos[0][2]
        else:
            xpos = configuration[pid][0]
            ypos = configuration[pid][1]
            zpos = configuration[pid][2]
        st = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f     T%04d%2s%2s\n"%('ATOM  ',pid % 100000,'FE','','UNX','F',mol % 10000,'',xpos,ypos,zpos,0,0,mol %1000,'','')
        file.write(st)
        molcnt += 1
        if molcnt == molsize:
            mol   += 1
            molcnt = 0
    file.write('END\n')
    file.close()

def pqrwrite(filename, system, molsize=4, append=False):
    if append:
        file = open(filename, 'a')
        s = "\n"
    else:
        file = open(filename,'w')
        s = "REMARK generated by ESPResSo++\n"
    file.write(s)
    maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
    pid    = 0
    addToPid = 0 # if pid begins from 0, then addToPid should be +1
    mol    = 0
    molcnt = 0
    #following http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    #crystal header
    st = "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n"%('CRYST1',system.bc.boxL[0],system.bc.boxL[1],system.bc.boxL[2],90.00,90.00,90,'P 1',1) #boxes are orthorhombic for now
    file.write(st)
    while pid <= maxParticleID:
        if system.storage.particleExists(pid):
            particle = system.storage.getParticle(pid)
            if(pid==0):
                addToPid = 1
            xpos   = particle.pos[0]
            ypos   = particle.pos[1]
            zpos   = particle.pos[2]
            type   = particle.type
            q      = particle.q
            radius = particle.radius
            #st = "ATOM %6d  FE  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d\n"%(pid, mol, xpos, ypos, zpos, type)
            #following http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
            st = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%8.3f%8.3f\n"%('ATOM  ',(pid+addToPid)%100000,'FE','','UNX','F',mol%10000,'',xpos,ypos,zpos,q,radius)

            # ATOM      1  N   ALA     1      46.457  12.189  21.556  0.1414 1.8240

            file.write(st)
            pid    += 1
            molcnt += 1
            if molcnt == molsize:
                mol   += 1
                molcnt = 0
        else:
            pid   += 1

    file.write('END\n')
    file.close()

def pdbread(filename,natoms,header):

    index=[]
    atomname=[]
    resname=[]
    locator=[]
    resid=[]
    resseq=[]
    x=[]
    y=[]
    z=[]
    alpha=[]
    beta=[]
    segid=[]
    element=[]

    file = open(filename,'r')
    for i in range(header):
        line = file.readline()
    for i in range(natoms):
        line = file.readline()
        findex=int(line[7:11])
        fname=line[12:16]
        flocator=line[16:17]
        fresname=line[17:20]
        fresid=line[21:22]
        fresseq=int(line[23:26])
        fx=float(line[30:38])
        fy=float(line[38:46])
        fz=float(line[46:54])
        falpha=float(line[55:60])
        fbeta=float(line[60:66])
        fsegid=line[72:76]
        felement=line[76:78]

        index.append(findex)
        atomname.append(fname)
        resname.append(fresname)
        locator.append(flocator)
        resid.append(fresid)
        resseq.append(fresseq)
        x.append(fx)
        y.append(fy)
        z.append(fz)
        alpha.append(falpha)
        beta.append(fbeta)
        segid.append(fsegid)
        element.append(felement)

    return index,atomname,locator,resname,resid,resseq,x,y,z,alpha,beta,segid,element
