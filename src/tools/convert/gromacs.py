# -*- coding: iso-8859-1 -*-
import math


"""This Python module allows one to use GROMACS data files as the
   input to an ESPResSo++ simulation. It can also convert GROMACS
   tabulated potentials file into ESPResSo++ format"""

# read GORMACS data files
def read(gro_file, top_file, itp_file):

    # gro_file contains number of particles, positions, velocities and box size
    # top_file contains topology information
    # itp_file contains topology information

    # read gro file
    if gro_file != "":
        f = open(gro_file)
        f.readline() # skip comment line
        num_particles = int(f.readline())
        
        # store coordinates
        x = []
        y = []
        z = []
        for i in range(num_particles):
            s = f.readline()[20:44]
            rx = float(s[0:8])
            ry = float(s[8:16])
            rz = float(s[16:24])
            x.append(rx)
            y.append(ry)
            z.append(rz)
        
        # store box size
        Lx, Ly, Lz = map(float, f.readline().split()) # read last line, convert to float
        f.close()


    # read top and itp file
    if top_file != "":
        f = open(top_file)
        # find and store number of chains
        line = ''
        while not 'molecules' in line:
            line = f.readline()
        f.readline() # skip comment line
        num_chains = int(f.readline().split()[1])
        monomers = num_particles / num_chains
            
        f.close()
        
        # read itp file (TODO: there could be more itp files...)
        if itp_file == "": # itp contents can be in top file
            itp_file = top_file
        
        f = open(itp_file)
        # find and store bonds
        line = ''
        while not 'bonds' in line:
            line = f.readline()
        bonds = []
        line = f.readline()
        line = f.readline()
        while(len(line) != 1):
            pid1, pid2 = map(int, line.split()[0:2])
            bonds.append((pid1, pid2))
            line = f.readline()
        
        # find and store angles
        line = ''
        while not 'angles' in line:
            line = f.readline()
        angles = []
        f.readline() # skip comment line
        line = f.readline()
        tmp = line.split()[0:3]
        while(len(tmp) == 3):
            pid1, pid2, pid3 = map(int, tmp)
            angles.append((pid1, pid2, pid3))
            line = f.readline()
            tmp = line.split()[0:3]
        f.close()


    # extend bonds to num_chains - 1 chains
    bonds_per_chain = len(bonds)
    for i in range(num_chains - 1):
        for j in range(bonds_per_chain):
            pid1 = bonds[j][0]
            pid2 = bonds[j][1]
            bonds.append((pid1 + (i + 1) * monomers, pid2 + (i + 1) * monomers))

    # extend angles to num_chains - 1 chains
    angles_per_chain = len(angles)
    for i in range(num_chains - 1):
        for j in range(angles_per_chain):
            pid1 = angles[j][0]
            pid2 = angles[j][1]
            pid3 = angles[j][2]
            angles.append((pid1 + (i + 1) * monomers, \
                            pid2 + (i + 1) * monomers, \
                            pid3 + (i + 1) * monomers))


    return bonds, angles, x, y, z, Lx, Ly, Lz




# Convert GROMACS tabulated file into
#  ESPResSo++ tabulated file (new file is created).
# First column can be either distance or angle.
# Sigma and epsilon are optional, depending on whether you want
#  to convert units or not.
# For non-bonded files, c6 and c12 can be provided. Electrostatics
#  are not taken into account (f and fd columns).
def convertTable(gro_in_file, esp_out_file, sigma=1.0, epsilon=1.0, c6=1.0, c12=1.0):

    # determine file type
    angle, bonded = False, False
    if gro_in_file[6] == "a" or gro_in_file[6] == "d":
        angle  = True
        bonded = True
    if gro_in_file[6] == "b":
        bonded = True

    fin = open(gro_in_file, 'r')
    fout= open(esp_out_file, 'w')

    if bonded: # bonded has 3 columns
        for line in fin:
            if line[0] == "#": # skip comment lines
                continue
            
            columns = line.split()
            r = float(columns[0])
            f = float(columns[1]) # energy
            fd= float(columns[2]) # force
            
            # convert units
            if angle: # degrees to radians
                r = math.radians(r)
            else:
                r = r / sigma
            e = f / epsilon
            f = fd*sigma / epsilon
            
            if (not angle and r != 0) or (angle and r <= 3.1415 and r > 0): # skip 0 and above 3.14 if angle
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))
    
    else: # non-bonded has 7 columns
        for line in fin:
            if line[0] == "#": # skip comment lines
                continue
            
            columns = line.split()
            r = float(columns[0])
            #f = float(columns[1]) # electrostatics, not implemented yet
            #fd= float(columns[2]) # not implemented yet
            g = float(columns[3]) # dispersion
            gd= float(columns[4])
            h = float(columns[5]) # repulsion
            hd= float(columns[6])
            
            e = c6*g + c12*h
            f = c6*gd+ c12*hd
            
            # convert units
            r = r / sigma
            e = e / epsilon
            f = f*sigma / epsilon
            
            if r != 0: # skip 0
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))
    
    fin.close()
    fout.close()









