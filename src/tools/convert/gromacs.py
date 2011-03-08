# -*- coding: utf-8 -*-
import math
from operator import itemgetter # for sorting a dict


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
            s = f.readline()[20:44] # atom coordinates
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
        readattypes, readbdtypes, readantypes, readdhtypes  = False, False, False, False
        atomtypes, bondtypes, angletypes, dihedraltypes = {}, {}, {}, {}
        molecule_numbers = []
        readmolecules = False
        for line in f:
            if 'atomtypes' in line: # map atom types (espresso++ uses ints)
                readattypes = True
                print "Reading atomtypes: ",
                continue
            
            if readattypes:
                if line[0] == ";":  # skip comment line
                    continue
                if line.strip() == "": # end of atomtypes section
                    readattypes = False
                    print sorted(atomtypes.items(), key=itemgetter(1)) # prints gromacs type and esp++ type
                    continue
                #print " "+line.strip('\n')
                tmp = line.split()[0]
                if tmp not in atomtypes:
                    atomtypes.update({tmp:a}) # atomtypes is used when reading the "atoms" section
                    a += 1
            
            if 'bondtypes' in line:
                readbdtypes = True
                continue
            
            if readbdtypes:
                if line[0] == ";":  # skip comment line
                    continue
                if line.strip() == "": # end of bondtypes section
                    readbdtypes = False
                    continue
                tmp = line.split()
                if tmp[2] == "8":
                    i, j, p = atomtypes[tmp[0]], atomtypes[tmp[1]], tmp[3]
                    if i in bondtypes:
                        bondtypes[i].update({j:p})
                    else:
                        bondtypes.update({i:{j:p}})
            
            if 'angletypes' in line:
                readantypes = True
                continue
            
            if readantypes:
                if line[0] == ";":  # skip comment line
                    continue
                if line.strip() == "": # end of angletypes section
                    readantypes = False
                    continue
                tmp = line.split()
                if tmp[3] == "8":
                    i, j, k, p = atomtypes[tmp[0]], atomtypes[tmp[1]], atomtypes[tmp[2]], tmp[4]
                    if i in angletypes:
                        if j in angletypes[i]:
                            angletypes[i][j].update({k:p})
                        else:
                            angletypes[i].update({j:{k:p}})
                    else:
                        angletypes.update({i:{j:{k:p}}})
            
            if 'dihedraltypes' in line:
                readdhtypes = True
                continue
            
            if readdhtypes:
                if line[0] == ";":  # skip comment line
                    continue
                if line.strip() == "": # end of angletypes section
                    readdhtypes = False
                    continue
                tmp = line.split()
                if tmp[4] == "8":
                    i, j, k, l, p = atomtypes[tmp[0]], atomtypes[tmp[1]], atomtypes[tmp[2]], atomtypes[tmp[3]], tmp[5]
                    if i in dihedraltypes:
                        if j in dihedraltypes[i]:
                            if k in dihedraltypes[i][j]:
                                dihedraltypes[i][j][k].update({l:p})
                            else:
                                dihedraltypes[i][j].update({k:{l:p}})
                        else:
                            dihedraltypes[i].update({j:{k:{l:p}}})
                    else:
                        dihedraltypes.update({i:{j:{k:{l:p}}}})
            
            if 'include' in line: # add included topology files
                #print " adding included top file: "+ line.strip('\n')
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
        
        #print bondtypes
        #print angletypes
        #print dihedraltypes
        #exit()
        
        if len(itp_files) == 0: # itp contents can be in top file
            itp_files = [top_file]
        
        # read itp files
        molecule = 0
        types = []
        bonds, angles, dihedrals = {}, {}, {}
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
            
            #print len(types), types
            #exit()
            
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
                tmp = line.split()
                pid1, pid2, fun = map(int, tmp[0:3])
                if (fun < 8):
                    pot = fun
                elif (fun > 7) and (len(tmp) == 3): # look for the potential in the bondtypes dict
                    t1, t2 = types[pid1-1], types[pid2-1]
                    if t1 > t2: # interactions in the other way
                        t1, t2 = t2, t1
                    pot = bondtypes[t1][t2]
                elif (fun > 7) and (len(tmp) > 3):
                    pot = tmp[3]
                else:
                    print "error while reading bond potential, line:"
                    print line
                    exit()
                bonds_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                   pid2 + ((monomers * num_chains) * (molecule-1)), pot))
                line = f.readline()
                
            #print bonds_tmp
            
            # extend bonds to num_chains - 1 chains
            bonds_per_chain = len(bonds_tmp)
            for i in range(num_chains):
                for j in range(bonds_per_chain):
                    pid1, pid2, pot = bonds_tmp[j][0:3]
                    #pid2 = bonds_tmp[j][1]
                    #bonds.append((pid1 + (i * monomers), \
                                   #pid2 + (i * monomers)))
                    if pot in bonds:
                        bonds[pot].append((pid1 + (i * monomers), \
                                           pid2 + (i * monomers)))
                    else:
                        bonds.update({pot:[(pid1 + (i * monomers), \
                                            pid2 + (i * monomers))]})

            #print "bonds len: "+str(len(bonds))
            #print bonds_tmp[0]
            #print bonds_tmp[bonds_per_chain-1]
            #print bonds[0]
            #print bonds[len(bonds)-1]
            #print bonds
            #exit()
            
            # find and store angles
            line = ''
            while not 'angles' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            #angles = []
            #f.readline() # skip comment line
            line = f.readline()
            while(len(line) != 1):
                if line[0] == ";": # skip comment lines
                    line = f.readline()
                    continue
                tmp = line.split()
                pid1, pid2, pid3, fun = map(int, tmp[0:4])
                if (fun < 8):
                    pot = fun
                elif (fun > 7) and (len(tmp) == 4): # look for the potential in the angletypes dict
                    t1, t2, t3 = types[pid1-1], types[pid2-1], types[pid3-1]
                    if t1 not in angletypes: # interactions in the other way
                        t1, t3 = t3, t1
                    pot = angletypes[t1][t2][t3]
                elif (fun > 7) and (len(tmp) > 4):
                    pot = tmp[4]
                else:
                    print "error while reading angle potential, line:"
                    print line
                    exit()
                angles_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                    pid2 + ((monomers * num_chains) * (molecule-1)), \
                                     pid3 + ((monomers * num_chains) * (molecule-1)), pot))
                line = f.readline()
                
            #print angles_tmp
            #exit()
            
            # extend angles to num_chains - 1 chains
            angles_per_chain = len(angles_tmp)
            for i in range(num_chains):
                for j in range(angles_per_chain):
                    pid1, pid2, pid3, pot = angles_tmp[j][0:4]
                    #pid2 = angles_tmp[j][1]
                    #pid3 = angles_tmp[j][2]
                    #angles.append((pid1 + (i * monomers), \
                                    #pid2 + (i * monomers), \
                                     #pid3 + (i * monomers)))
                    if pot in angles:
                        angles[pot].append((pid1 + (i * monomers), \
                                            pid2 + (i * monomers), \
                                            pid3 + (i * monomers)))
                    else:
                        angles.update({pot:[(pid1 + (i * monomers), \
                                             pid2 + (i * monomers), \
                                             pid3 + (i * monomers))]})
            
            
            
            
            # find and store dihedrals
            line = ''
            while not 'dihedrals' in line:
                line = f.readline()
                if not line: break # break out of while if EOF
            #dihedrals = []
            line = f.readline()
            while(len(line) != 0):
                if line[0] == ";": # skip comment lines
                    line = f.readline()
                    continue
                tmp = line.split()
                pid1, pid2, pid3, pid4, fun = map(int, tmp[0:5])
                if (fun < 8):
                    pot = fun
                elif (fun > 7) and (len(tmp) == 5): # look for the potential in the dihedraltypes dict
                    t1, t2, t3, t4 = types[pid1-1], types[pid2-1], types[pid3-1], types[pid4-1] # get types of particles
                    if t1 not in dihedraltypes: # interactions in the other way
                        t1, t2, t3, t4 = t4, t1, t2, t3
                    pot = dihedraltypes[t1][t2][t3][t4]
                elif (fun > 7) and (len(tmp) > 5):
                    pot = tmp[5]
                else:
                    print "error while reading dihedral potential, line:"
                    print line
                    exit()
                dihedrals_tmp.append((pid1 + ((monomers * num_chains) * (molecule-1)), \
                                       pid2 + ((monomers * num_chains) * (molecule-1)), \
                                        pid3 + ((monomers * num_chains) * (molecule-1)), \
                                         pid4 + ((monomers * num_chains) * (molecule-1)), pot ))
                line = f.readline()
                
            #print dihedrals_tmp
            #exit()
            
            
            # extend dihedrals to num_chains - 1 chains
            dihedrals_per_chain = len(dihedrals_tmp)
            for i in range(num_chains):
                for j in range(dihedrals_per_chain):
                    pid1, pid2, pid3, pid4, pot = dihedrals_tmp[j][0:5]
                    #pid2 = dihedrals_tmp[j][1]
                    #pid3 = dihedrals_tmp[j][2]
                    #pid4 = dihedrals_tmp[j][3]
                    if pot in dihedrals:
                        dihedrals[pot].append((pid1 + (i * monomers), \
                                               pid2 + (i * monomers), \
                                               pid3 + (i * monomers), \
                                               pid4 + (i * monomers)))
                    else:
                        dihedrals.update({pot:[(pid1 + (i * monomers), \
                                                pid2 + (i * monomers), \
                                                pid3 + (i * monomers), \
                                                pid4 + (i * monomers))]})
            
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



# Set interactions for all given particle types.
#  "potentials" is a dictionary where key is a string composed
#   of two particle types and value is a potential.
#    example: {"A_A":potAA, "A_B":potAB, "B_B":potBB}
#  "particleTypes" is a dictionary where key is the particle type, and
#  value is a list of particles of that type.
#    example: {"A":["A1m", "A2m"],"B":["B1u","B2u"]}
#  "system" is the system to which the interaction will be added
#  "interaction" is the interaction to which to add the potentials
# Return value is a system with all interactions added.
def setInteractions(potentials, particleTypes, system, interaction):
    
    allparticles = []
    for k, v in particleTypes.iteritems():
        for i in v:
            allparticles.append((i,k)) # create tuples: (particle, type)
    
    for i in range(len(allparticles)):
        for j in range(i, len(allparticles)):
            type1 = allparticles[i][1]
            type2 = allparticles[j][1]
            key = type1+"_"+type2
            interaction.setPotential(i, j, potentials[key])

    system.addInteraction(interaction)
    return system



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
