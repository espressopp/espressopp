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

"""
************************************
espresso_old - read espressomd files
************************************

This Python module allows one to use ESPResSo data files as the
input to an ESPResSo++ simulation.

"""

import math

def read(file):
    """ Read ESPResSo data files.

    Keyword argument:
    file -- contains simulation variables, data of all particles, and information about bonds.
    (angles and dihedrals are currently not read)
    """

    # read file
    if file != "":
      Lx, Ly, Lz = 0, 0, 0 
      props = []
      x, y, z = [], [], []
      type, q = [], []
      vx, vy, vz = [], [], []
      fx, fy, fz = [], [], []
      bondpairs = []

      variable = False
      particles = False
      bonds = False

      f = open(file)
      for line in f:
        if line.strip() == "":
            continue # skip empty lines

        # read variables
        if line[1:9] == 'variable':
          variable = True
          continue # goto next line

        if variable == True:
          if line.strip() == "}":  # reached the end of variable section
            variable = False
            continue

          line = line.replace('{','').replace('}','')
          tmp = line.split()
          if tmp[0] == "box_l":
            Lx, Ly, Lz = map(float, [tmp[1], tmp[2], tmp[3]])
          continue # goto next line


        # read particle properties
        if line[1:10] == 'particles':
          line = line.replace('}','')
          tmp = line[12:len(line)].split()
          particles = True
          for prop in tmp:
            if prop == "id":
              props.append(prop)
            if prop == "pos":
              props.append(prop)
            if prop == "type":
              props.append(prop)
            if prop == "q":
              props.append(prop)
            if prop == "v":
              props.append(prop)
            if prop == "f":
              props.append(prop)
          continue # goto next line

        # read particle data
        if particles == True:
          if line.strip() == "}":  # reached the end of particles section
            particles = False
            continue

          line = line.replace('{','').replace('}','')
          tmp = line.split()
          index = -1
          for prop in props:
            index = index+1
            if prop == "id":
              continue
            if prop == "pos":
              x.append(float(tmp[index]))
              y.append(float(tmp[index+1]))
              z.append(float(tmp[index+2]))
              index = index+2
            if prop == "type":
              type.append(int(tmp[index]))
            if prop == "q":
              q.append(float(tmp[index]))
            if prop == "v":
              vx.append(float(tmp[index]))
              vy.append(float(tmp[index+1]))
              vz.append(float(tmp[index+2]))
              index = index+2
            if prop == "f":
              fx.append(float(tmp[index]))
              fy.append(float(tmp[index+1]))
              fz.append(float(tmp[index+2]))
              index = index+2


        # read bond information
        if line[1:6] == 'bonds':
          bonds = True
          continue # goto next line

        if bonds == True:
          if line.strip() == "}":  # reached the end of bonds section
            bonds = False
            continue

          line = line.replace('{','').replace('}','')
          tmp = line.split()

          if len(tmp) > 2:
            first = int(tmp[0])
            for idx in xrange(2, len(tmp), 2):
              bondpairs.append((first, int(tmp[idx])))
              
      f.close()

    return Lx, Ly, Lz, x, y, z, type, q, vx, vy, vz, fx, fy, fz, bondpairs   

    
