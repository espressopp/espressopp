<<<<<<< local
# -*- coding: utf-8 -*-
import math


"""This Python module allows one to use GROMACS data files as the
   input to an ESPResSo++ simulation. It can also convert GROMACS
   tabulated potentials file into ESPResSo++ format"""

# read GORMACS data files
def read(gro_file, top_file):

    # gro_file contains number of particles, positions, velocities and box size
    # top_file contains topology information
    # itp_files contains a list of included topology information

    # read gro file
    if gro_file != "":
        f = open(gro_file)
        f.readline() # skip comment line
        num_particles = int(f.readline())
        
        # store coordinates
        x, y, z = [], [], []
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


    # read top and itp files
    if top_file != "":
        f = open(top_file)
        
        print "Reading top file: "+top_file
        line = ''
        itp_files = []
        a = 0
        readtypes = False
        atomtypes = {}
        molecule_numbers = []
        readmolecules = False
        for line in f:
            if 'atomtypes' in line: # map atom types (espresso uses ints)
                readtypes = True
                print "Reading atomtypes"
                continue
            
            if readtypes:
                if line[0] == ";":  # skip comment line
                    continue
                if line.strip() == "": # end of atomtypes section
                    readtypes = False
                    continue
                #print " "+line.strip('\n')
                tmp = line.split()[0]
                if tmp not in atomtypes:
                    atomtypes.update({tmp:a})
                    a += 1
            
            
            if 'include' in line: # add included topology files
                print " adding included top file: "+ line.strip('\n')
                itp_files.append((line.split()[1]).strip('\"')) # strip " and add filename
            
            if 'molecules' in line: # store number of chains
                readmolecules = True
                print "Reading number of molecules: "
                continue
            
            if readmolecules:
                if line[0] == ";":  # skip comment line
                    continue
                if line.strip() == "": # end of molecules section
                    readmolecules = False
                    continue
                print " "+line.strip('\n')
                molecule_numbers.append(int(line.split()[1]))
            
            ## find and store number of chains
            #line = ''
            #while not 'molecules' in line:
                #line = f.readline()
            #f.readline() # skip comment line
            #num_chains = int(f.readline().split()[1])
            #monomers = num_particles / num_chains
        f.close()
        
        if len(itp_files) == 0: # itp contents can be in top file
            itp_files = [top_file]
        
        # read itp files
        molecule = 0
        types, bonds, angles, dihedrals = [], [], [], []
        print "Reading included topology:"
        for itp_file in itp_files:
            print " "+itp_file
            types_tmp, bonds_tmp, angles_tmp, dihedrals_tmp = [], [], [], []
            molecule += 1
            num_chains = molecule_numbers[molecule-1] # read number of molecules from array
            monomers = num_particles / (num_chains * len(molecule_numbers))
            f = open(itp_file)
            
            # find and store atom types
            line = ''
            while not 'atoms' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            line = f.readline()
            while(len(line) != 1):
                #print " current line: "+line.strip("\n")
                if line[0] == ";":   # skip comment lines
                    #print "skipping line: "+line.strip("\n")
                    line = f.readline()
                    continue
                types_tmp.append(atomtypes[line.split()[1]]) # map str type to int type
                line = f.readline()
            
            # extend types to num_chains - 1 chains
            types_per_chain = len(types_tmp)
            for i in range(num_chains):
                types.extend(types_tmp)
            
            
            # find and store bonds
            line = ''
            while not 'bonds' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            #bonds = []
            line = f.readline()
            while(len(line) != 1):
                #print " current line: "+line.strip("\n")
                if line[0] == ";":   # skip comment lines
                    #print "skipping line: "+line.strip("\n")
                    line = f.readline()
                    continue
                pid1, pid2 = map(int, line.split()[0:2])
                bonds_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                   pid2 + ((monomers * num_chains) * (molecule-1))))
                line = f.readline()
            #print "pid + (%1d * (%1d))" % (monomers, molecule-1)
            
            # extend bonds to num_chains - 1 chains
            bonds_per_chain = len(bonds_tmp)
            for i in range(num_chains):
                for j in range(bonds_per_chain):
                    pid1 = bonds_tmp[j][0]
                    pid2 = bonds_tmp[j][1]
                    bonds.append((pid1 + (i * monomers), \
                                   pid2 + (i * monomers)))
            #print "bonds len: "+str(len(bonds))
            #print bonds_tmp[0]
            #print bonds_tmp[bonds_per_chain-1]
            #print bonds[0]
            #print bonds[len(bonds)-1]
            #exit()
            
            # find and store angles
            line = ''
            while not 'angles' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            #angles = []
            f.readline() # skip comment line
            line = f.readline()
            tmp = line.split()[0:3]
            while(len(tmp) == 3):
                pid1, pid2, pid3 = map(int, tmp)
                angles_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                    pid2 + ((monomers * num_chains) * (molecule-1)), \
                                     pid3 + ((monomers * num_chains) * (molecule-1))))
                line = f.readline()
                tmp = line.split()[0:3]
            
            # extend angles to num_chains - 1 chains
            angles_per_chain = len(angles_tmp)
            for i in range(num_chains):
                for j in range(angles_per_chain):
                    pid1 = angles_tmp[j][0]
                    pid2 = angles_tmp[j][1]
                    pid3 = angles_tmp[j][2]
                    angles.append((pid1 + (i * monomers), \
                                    pid2 + (i * monomers), \
                                     pid3 + (i * monomers)))
            
            
            # find and store dihedrals
            line = ''
            while not 'dihedrals' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            #dihedrals = []
            f.readline() # skip comment line
            line = f.readline()
            tmp = line.split()[0:4]
            while(len(tmp) == 4):
                pid1, pid2, pid3, pid4 = map(int, tmp)
                dihedrals_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                       pid2 + ((monomers * num_chains) * (molecule-1)), \
                                        pid3 + ((monomers * num_chains) * (molecule-1)), \
                                         pid4 + ((monomers * num_chains) * (molecule-1))))
                line = f.readline()
                tmp = line.split()[0:4]
            
            # extend dihedrals to num_chains - 1 chains
            dihedrals_per_chain = len(dihedrals_tmp)
            for i in range(num_chains):
                for j in range(dihedrals_per_chain):
                    pid1 = dihedrals_tmp[j][0]
                    pid2 = dihedrals_tmp[j][1]
                    pid3 = dihedrals_tmp[j][2]
                    pid4 = dihedrals_tmp[j][3]
                    dihedrals.append((pid1 + (i * monomers), \
                                       pid2 + (i * monomers), \
                                        pid3 + (i * monomers), \
                                         pid4 + (i * monomers)))
            
            f.close()
        #exit()

    #molecule = 0
    #for num_chains in molecule_numbers:
        #molecule += 1
        #monomers = num_particles / num_chains
        
        ## extend bonds to num_chains - 1 chains
        #bonds_per_chain = len(bonds)
        #for i in range(num_chains - 1):
            #for j in range(bonds_per_chain):
                #pid1 = bonds[j][0]
                #pid2 = bonds[j][1]
                #bonds.append((pid1 + (i + 1) * monomers * molecule, \
                              #pid2 + (i + 1) * monomers * molecule))
        
        ## extend angles to num_chains - 1 chains
        #angles_per_chain = len(angles)
        #for i in range(num_chains - 1):
            #for j in range(angles_per_chain):
                #pid1 = angles[j][0]
                #pid2 = angles[j][1]
                #pid3 = angles[j][2]
                #angles.append((pid1 + (i + 1) * monomers * molecule, \
                               #pid2 + (i + 1) * monomers * molecule, \
                               #pid3 + (i + 1) * monomers * molecule))
        
        ## extend dihedrals to num_chains - 1 chains
        #dihedrals_per_chain = len(dihedrals)
        #for i in range(num_chains - 1):
            #for j in range(dihedrals_per_chain):
                #pid1 = dihedrals[j][0]
                #pid2 = dihedrals[j][1]
                #pid3 = dihedrals[j][2]
                #pid4 = dihedrals[j][3]
                #dihedrals.append((pid1 + (i + 1) * monomers * molecule, \
                                  #pid2 + (i + 1) * monomers * molecule, \
                                  #pid3 + (i + 1) * monomers * molecule, \
                                  #pid4 + (i + 1) * monomers * molecule))

    #if len(angles) != 0 and len(dihedrals) == 0:
        #return types, bonds, angles, x, y, z, Lx, Ly, Lz
    #elif len(dihedrals) != 0:
        #return types, bonds, angles, dihedrals, x, y, z, Lx, Ly, Lz
    #else:
        #return types, bonds, x, y, z, Lx, Ly, Lz

    params = []
    if len(types) != 0:
        params.append(types)
    if len(bonds) != 0:
        params.append(bonds)
    if len(angles) != 0:
        params.append(angles)
    if len(dihedrals) != 0:
        params.append(dihedrals)
    params.extend([x, y, z, Lx, Ly, Lz])
    return tuple(params)






# Convert GROMACS tabulated file into
#  ESPResSo++ tabulated file (new file is created).
# First column can be either distance or angle.
# Sigma and epsilon are optional, depending on whether you want
#  to convert units or not.
# For non-bonded files, c6 and c12 can be provided. Electrostatics
#  are not taken into account (f and fd columns).
def convertTable(gro_in_file, esp_out_file, sigma=1.0, epsilon=1.0, c6=1.0, c12=1.0):

    # determine file type
    bonded, angle, dihedral = False, False, False
    if gro_in_file[6] == "b":
        bonded = True
    if gro_in_file[6] == "a":
        angle  = True
        bonded = True
    if gro_in_file[6] == "d":
        dihedral = True
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
            if angle or dihedral: # degrees to radians
                r = math.radians(r)
            else:
                r = r / sigma
            e = f / epsilon
            f = fd*sigma / epsilon
            
            if (not angle and not dihedral and r != 0) or \
                 (angle and r <= 3.1415 and r > 0) or \
                  (dihedral and r >= -3.1415 and r <= 3.1415):
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









=======
# -*- coding: utf-8 -*-
import math


"""This Python module allows one to use GROMACS data files as the
   input to an ESPResSo++ simulation. It can also convert GROMACS
   tabulated potentials file into ESPResSo++ format"""

# read GORMACS data files
def read(gro_file, top_file):

    # gro_file contains number of particles, positions, velocities and box size
    # top_file contains topology information
    # itp_files contains a list of included topology information

    # read gro file
    if gro_file != "":
        f = open(gro_file)
        f.readline() # skip comment line
        num_particles = int(f.readline())
        
        # store coordinates
        x, y, z = [], [], []
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


    # read top and itp files
    if top_file != "":
        f = open(top_file)
        
        print "Reading top file: "+top_file
        line = ''
        itp_files = []
        a = 0
        readtypes = False
        atomtypes = {}
        molecule_numbers = []
        readmolecules = False
        for line in f:
            if 'atomtypes' in line: # map atom types (espresso uses ints)
                readtypes = True
                print "Reading atomtypes"
                continue
            
            if readtypes:
                if line[0] == ";":  # skip comment line
                    continue
                if line.strip() == "": # end of atomtypes section
                    readtypes = False
                    continue
                #print " "+line.strip('\n')
                tmp = line.split()[0]
                if tmp not in atomtypes:
                    atomtypes.update({tmp:a})
                    a += 1
            
            
            if 'include' in line: # add included topology files
                print " adding included top file: "+ line.strip('\n')
                itp_files.append((line.split()[1]).strip('\"')) # strip " and add filename
            
            if 'molecules' in line: # store number of chains
                readmolecules = True
                print "Reading number of molecules: "
                continue
            
            if readmolecules:
                if line[0] == ";":  # skip comment line
                    continue
                if line.strip() == "": # end of molecules section
                    readmolecules = False
                    continue
                print " "+line.strip('\n')
                molecule_numbers.append(int(line.split()[1]))
            
            ## find and store number of chains
            #line = ''
            #while not 'molecules' in line:
                #line = f.readline()
            #f.readline() # skip comment line
            #num_chains = int(f.readline().split()[1])
            #monomers = num_particles / num_chains
        f.close()
        
        if len(itp_files) == 0: # itp contents can be in top file
            itp_files = [top_file]
        
        # read itp files
        molecule = 0
        types, bonds, angles, dihedrals = [], [], [], []
        print "Reading included topology:"
        for itp_file in itp_files:
            print " "+itp_file
            types_tmp, bonds_tmp, angles_tmp, dihedrals_tmp = [], [], [], []
            molecule += 1
            num_chains = molecule_numbers[molecule-1] # read number of molecules from array
            monomers = num_particles / (num_chains * len(molecule_numbers))
            f = open(itp_file)
            
            # find and store atom types
            line = ''
            while not 'atoms' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            line = f.readline()
            while(len(line) != 1):
                #print " current line: "+line.strip("\n")
                if line[0] == ";":   # skip comment lines
                    #print "skipping line: "+line.strip("\n")
                    line = f.readline()
                    continue
                types_tmp.append(atomtypes[line.split()[1]]) # map str type to int type
                line = f.readline()
            
            # extend types to num_chains - 1 chains
            types_per_chain = len(types_tmp)
            for i in range(num_chains):
                types.extend(types_tmp)
            
            
            # find and store bonds
            line = ''
            while not 'bonds' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            #bonds = []
            line = f.readline()
            while(len(line) != 1):
                #print " current line: "+line.strip("\n")
                if line[0] == ";":   # skip comment lines
                    #print "skipping line: "+line.strip("\n")
                    line = f.readline()
                    continue
                pid1, pid2 = map(int, line.split()[0:2])
                bonds_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                   pid2 + ((monomers * num_chains) * (molecule-1))))
                line = f.readline()
            #print "pid + (%1d * (%1d))" % (monomers, molecule-1)
            
            # extend bonds to num_chains - 1 chains
            bonds_per_chain = len(bonds_tmp)
            for i in range(num_chains):
                for j in range(bonds_per_chain):
                    pid1 = bonds_tmp[j][0]
                    pid2 = bonds_tmp[j][1]
                    bonds.append((pid1 + (i * monomers), \
                                   pid2 + (i * monomers)))
            #print "bonds len: "+str(len(bonds))
            #print bonds_tmp[0]
            #print bonds_tmp[bonds_per_chain-1]
            #print bonds[0]
            #print bonds[len(bonds)-1]
            #exit()
            
            # find and store angles
            line = ''
            while not 'angles' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            #angles = []
            f.readline() # skip comment line
            line = f.readline()
            tmp = line.split()[0:3]
            while(len(tmp) == 3):
                pid1, pid2, pid3 = map(int, tmp)
                angles_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                    pid2 + ((monomers * num_chains) * (molecule-1)), \
                                     pid3 + ((monomers * num_chains) * (molecule-1))))
                line = f.readline()
                tmp = line.split()[0:3]
            
            # extend angles to num_chains - 1 chains
            angles_per_chain = len(angles_tmp)
            for i in range(num_chains):
                for j in range(angles_per_chain):
                    pid1 = angles_tmp[j][0]
                    pid2 = angles_tmp[j][1]
                    pid3 = angles_tmp[j][2]
                    angles.append((pid1 + (i * monomers), \
                                    pid2 + (i * monomers), \
                                     pid3 + (i * monomers)))
            
            
            # find and store dihedrals
            line = ''
            while not 'dihedrals' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            #dihedrals = []
            f.readline() # skip comment line
            line = f.readline()
            tmp = line.split()[0:4]
            while(len(tmp) == 4):
                pid1, pid2, pid3, pid4 = map(int, tmp)
                dihedrals_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                       pid2 + ((monomers * num_chains) * (molecule-1)), \
                                        pid3 + ((monomers * num_chains) * (molecule-1)), \
                                         pid4 + ((monomers * num_chains) * (molecule-1))))
                line = f.readline()
                tmp = line.split()[0:4]
            
            # extend dihedrals to num_chains - 1 chains
            dihedrals_per_chain = len(dihedrals_tmp)
            for i in range(num_chains):
                for j in range(dihedrals_per_chain):
                    pid1 = dihedrals_tmp[j][0]
                    pid2 = dihedrals_tmp[j][1]
                    pid3 = dihedrals_tmp[j][2]
                    pid4 = dihedrals_tmp[j][3]
                    dihedrals.append((pid1 + (i * monomers), \
                                       pid2 + (i * monomers), \
                                        pid3 + (i * monomers), \
                                         pid4 + (i * monomers)))
            
            f.close()
        #exit()

    #molecule = 0
    #for num_chains in molecule_numbers:
        #molecule += 1
        #monomers = num_particles / num_chains
        
        ## extend bonds to num_chains - 1 chains
        #bonds_per_chain = len(bonds)
        #for i in range(num_chains - 1):
            #for j in range(bonds_per_chain):
                #pid1 = bonds[j][0]
                #pid2 = bonds[j][1]
                #bonds.append((pid1 + (i + 1) * monomers * molecule, \
                              #pid2 + (i + 1) * monomers * molecule))
        
        ## extend angles to num_chains - 1 chains
        #angles_per_chain = len(angles)
        #for i in range(num_chains - 1):
            #for j in range(angles_per_chain):
                #pid1 = angles[j][0]
                #pid2 = angles[j][1]
                #pid3 = angles[j][2]
                #angles.append((pid1 + (i + 1) * monomers * molecule, \
                               #pid2 + (i + 1) * monomers * molecule, \
                               #pid3 + (i + 1) * monomers * molecule))
        
        ## extend dihedrals to num_chains - 1 chains
        #dihedrals_per_chain = len(dihedrals)
        #for i in range(num_chains - 1):
            #for j in range(dihedrals_per_chain):
                #pid1 = dihedrals[j][0]
                #pid2 = dihedrals[j][1]
                #pid3 = dihedrals[j][2]
                #pid4 = dihedrals[j][3]
                #dihedrals.append((pid1 + (i + 1) * monomers * molecule, \
                                  #pid2 + (i + 1) * monomers * molecule, \
                                  #pid3 + (i + 1) * monomers * molecule, \
                                  #pid4 + (i + 1) * monomers * molecule))

    #if len(angles) != 0 and len(dihedrals) == 0:
        #return types, bonds, angles, x, y, z, Lx, Ly, Lz
    #elif len(dihedrals) != 0:
        #return types, bonds, angles, dihedrals, x, y, z, Lx, Ly, Lz
    #else:
        #return types, bonds, x, y, z, Lx, Ly, Lz

    params = []
    if len(types) != 0:
        params.append(types)
    if len(bonds) != 0:
        params.append(bonds)
    if len(angles) != 0:
        params.append(angles)
    if len(dihedrals) != 0:
        params.append(dihedrals)
    params.extend([x, y, z, Lx, Ly, Lz])
    return tuple(params)






# Convert GROMACS tabulated file into
#  ESPResSo++ tabulated file (new file is created).
# First column can be either distance or angle.
# Sigma and epsilon are optional, depending on whether you want
#  to convert units or not.
# For non-bonded files, c6 and c12 can be provided. Electrostatics
#  are not taken into account (f and fd columns).
def convertTable(gro_in_file, esp_out_file, sigma=1.0, epsilon=1.0, c6=1.0, c12=1.0):

    # determine file type
    bonded, angle, dihedral = False, False, False
    if gro_in_file[6] == "b":
        bonded = True
    if gro_in_file[6] == "a":
        angle  = True
        bonded = True
    if gro_in_file[6] == "d":
        dihedral = True
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
            if angle or dihedral: # degrees to radians
                r = math.radians(r)
            else:
                r = r / sigma
            e = f / epsilon
            f = fd*sigma / epsilon
            
            if (not angle and not dihedral and r != 0) or \
                 (angle and r <= 3.1415 and r > 0) or \
                  (dihedral and r >= -3.1415 and r <= 3.1415):
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









>>>>>>> other
