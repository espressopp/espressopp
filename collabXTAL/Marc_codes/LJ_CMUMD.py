#!/usr/bin/python

########################################################################
# specification of the main simulation parameters                      #
########################################################################

# === To Set =================================================================================================================================================================================
MODUS               = "auto"		# simulation modus					["start", "cont" or "auto"]

L_LATT_PARAM_SIGMA  = [60,27,27]	# box edge lengths in units of latt. parameter		[Lx will be rescaled]
NPART_TOT           = 224000		# total number of particles in the system
NUMDENS             = 0.5400		# number density
MOLFRAC             = 0.8000		# molar fraction solute

SIGMAWW             = 1.0000		# sigma parameter solvent-solvent
EPSILONWW           = 1.0000		# epsilon parameter solvent-solvent
CUTOFFWW            = 2.5000		# cutoff parameter solvent-solvent			[in units of sigma, only need for full LJ interaction]
WWTYPEWW            = "lj"		# interaction type solvent-solvent 			["lj" or "wca"]
SIGMASW             = "mix"		# sigma parameter solute-solvent 			["mix" or value]
EPSILONSW           = "mix"		# epsilon parameter solute-solvent 			["mix" or value]
CUTOFFSW            = 2.5000		# cutoff parameter solute-solvent			[in units of sigma, only need for full LJ interaction]
WWTYPESW            = "lj"		# interaction type solute-solvent 			["lj" or "wca"]
SIGMASS             = 1.1765		# sigma parameter solute-solute
EPSILONSS           = 1.6000		# epsilon parameter solute-solute
CUTOFFSS            = 2.5000		# cutoff parameter solute-solute			[in units of sigma, only need for full LJ interaction]
WWTYPESS            = "lj"		# interaction type solute-solute 			["lj" or "wca"]
MASSW               = 1.0000		# mass of solvent particles
MASSS               = 1.0000		# mass of solute particles

ENSEMBLE_EQ         = "npt"		# thermodynamic ensemble for equilibration		["nvt" or "npt"]
ENSEMBLE_SIM        = "npt"		# thermodynamic ensemble for simulation			["nvt" or "npt"]
DAMPING_EQ          = 1.0000		# damping constant for thermostat (eq)
DAMPING_SIM         = 1.0000		# damping constant for thermostat (sim)
TEMPERATURE         = 0.8000		# temperature						
PRESSURE            = 0.0024		# pressure Pxx						["auto" or value, only need for constant-P ensembles]
PISTON_MASS_EQ      = 1.0000		# effective piston mass for barostat (eq)		[only need for constant-P ensembles]
PISTON_MASS_SIM     = 5.0000		# effective piston mass for barostat (sim)              [only need for constant-P ensembles]

DT_WARMUP           = 0.0050		# warmup time step
DT_SIM              = 0.0050		# md time step
STRUCTURE           = "100"		# structure of crystalline substrate			["100" or "111"]
NN_DIST_XTAL        = 1.3206		# distance to nearest neighbour in substrate		[in units of sigma]
K_XTAL              = 1000.0		# spring constant for substrate bonds			

CMUMD_MODE          = "solute"		# switch to apply force in membrane			["solute", "solvent", "both" or "off"]
SIZETR              = 3.0000		# size of hybrid region
SIZEFR              = 3.0000		# size of force applying region
SIZECHECKR          = 20.000		# size of region for adaptation checker
FORCE_DIST          = 7.5000		# distance between center force region and substrate
EXP_MOLF            = 0.8000		# target mole fraction
FORCE_CONST_W       = 100.00		# force constant for solvent particles in membrane
FORCE_CONST_S       = 1000.0		# force constant for solute particles in membrane
SHAPE_PARAM         = 1.6000		# shape parameter for function G(x,XF)

TOTAL_WARMUP_STEPS  = 2500		# total number of warmup steps
WARMUP_ISTEPS       = 50		# number of warmup steps per loop
TOTAL_EQUIL_STEPS   = 10000		# total number of equilibration steps
EQUIL_ISTEPS        = 10		# number of equilibration steps per loop
TOTAL_SIMWARM_STEPS = 10000		# total number of sim preceding warmup steps
SIMWARM_ISTEPS      = 1000		# number of sim preceding warmup steps per loop
TOTAL_SIM_STEPS     = 60000		# total number of simulation steps			
SIM_ISTEPS          = 10		# number of simulation steps per loop
DATA_EVERY          = 100		# md steps between data outputs
DATA_UNITS          = "lj"		# units of data output					["lj" or "si"]
CONFIG_EVERY_EQ     = 10000		# md steps between configuration outputs (equilibration)
CONFIG_EVERY_SIM    = 4000		# md steps between configuration outputs (simulation)
RESTART_EVERY       = 50000		# md steps between restart file writings

EPSILON_START       = 0.1		# epsilon parameter at the beginning of warmup
SKIN                = 0.2		# skin for cell system

SEED                = "aut"		# seed for random number generator			["def"ault, "aut"o or value]
NODESETUP           = "aut"		# type of node grid setup				["aut"o or "man"ual]
# ============================================================================================================================================================================================


# === import modules =========================================================================================================================================================================
import sys
import os
import espresso
from   espresso import Real3D, Int3D, Tensor
import math
from   math import log, sqrt
import time
import gzip
# ============================================================================================================================================================================================

# === derived units and parameters ===========================================================================================================================================================
# Constants
SQRT2              = sqrt(2.0)
SQRT3              = sqrt(3.0)
PI                 = 3.1415926535897932384626433832
LJ_MIN  	   = pow(2.0, 1.0/6.0)

# SI units         = [ps, nm, K, kJ/mol, bar, g/mol, nm^3]
SI_UNITS           = {"t": 2.1655, "l": 0.3405, "T": 119.8, "E": 0.99607, "P": 418.98, "m": 40.288, "V": 0.039477655}

# units list       = [time, length, temperature, energy, pressure, mass, volume]
if DATA_UNITS== "lj":
  UNITS            = {"t": 1.0, "l": 1.0, "T": 1.0, "E": 1.0, "P": 1.0, "m": 1.0, "V": 1.0}
elif DATA_UNITS == "si":
  UNITS            = SI_UNITS
else:
  print "# ERROR [parameter_setting]: Unit type '%s' is not known!" % DATA_UNITS
  sys.exit()


ENSEMBLES          = {"eq": ENSEMBLE_EQ, "sim": ENSEMBLE_SIM}


if WWTYPEWW == "lj":
  CUTOFFWW        *= SIGMAWW
  SHIFTYWW         = 0.0
elif WWTYPEWW == "wca":
  CUTOFFWW         = LJ_MIN*SIGMAWW
  SHIFTYWW         = 'auto'
else:
  print "# ERROR [parameter_setting]: W-W interaction type '%s' is not known!" % WWTYPEWW
  sys.exit()
WARMUP_CUTOFFWW    = LJ_MIN*SIGMAWW

if SIGMASW == "mix":
  SIGMASW          = 0.5*(SIGMASS+SIGMAWW)
  SWPARAMSETUPS    = "mix"
else:
  SWPARAMSETUPS    = "man"
if EPSILONSW == "mix":
  EPSILONSW        = sqrt(EPSILONSS*EPSILONWW)
  SWPARAMSETUPE    = "mix"
else:PISTON_MASS_SIM
  SWPARAMSETUPE    = "man"
if WWTYPESW == "lj":
  CUTOFFSW        *= SIGMASW
  SHIFTYSW         = 0.0
elif WWTYPESW == "wca":
  CUTOFFSW         = LJ_MIN*SIGMASW
  SHIFTYSW         = 'auto'
else:
  print "# ERROR [parameter_setting]: S-W interaction type '%s' is not known!" % WWTYPESW
  sys.exit()
WARMUP_CUTOFFSW    = LJ_MIN*SIGMASW

if WWTYPESS == "lj":
  CUTOFFSS        *= SIGMASS
  SHIFTYSS         = 0.0
elif WWTYPESS == "wca":
  CUTOFFSS         = LJ_MIN*SIGMASS
  SHIFTYSS         = 'auto'
else:
  print "# ERROR [parameter_setting]: S-S interaction type '%s' is not known!" % WWTYPESS
  sys.exit()
WARMUP_CUTOFFSS    = LJ_MIN*SIGMASS


TEMPERATURE_SSUNITS= TEMPERATURE/EPSILONSS
PRESSURE_SSUNITS   = PRESSURE*pow(SIGMASS,3.0)/EPSILONSS

  
if STRUCTURE == "100":
  NPART_XTAL_SOLL       = 16*L_LATT_PARAM_SIGMA[1]*L_LATT_PARAM_SIGMA[2]
elif STRUCTURE == "111":
  NPART_TOT             = int(NPART_TOT*SQRT3)
  L_LATT_PARAM_SIGMA[1] = int(L_LATT_PARAM_SIGMA[1]*SQRT3)
  NPART_XTAL_SOLL       = 18*L_LATT_PARAM_SIGMA[1]*L_LATT_PARAM_SIGMA[2]
else:
  print "# ERROR [parameter_setting]: Structure type '%s' is not known!" % STRUCTURE
  sys.exit()

V_S                = PI*pow(SIGMASS,3)/6.0
V_W                = PI*pow(SIGMAWW,3)/6.0  
  
NPART_LIQ          = NPART_TOT
NPART_SOLUTE       = int(MOLFRAC*NPART_LIQ)
NPART_SOLVENT      = NPART_LIQ-NPART_SOLUTE
NPART_XTAL         = 0
V_TOT              = (NPART_SOLVENT+NPART_SOLUTE)/NUMDENS
NUMDENS0           = NUMDENS*MOLFRAC
NUMDENS1           = NUMDENS*(1.0-MOLFRAC)

NN_DIST_XTAL      *= SIGMAWW
if STRUCTURE == "100":
  LATTICE_PARAMETER= SQRT2*NN_DIST_XTAL
  LY               = L_LATT_PARAM_SIGMA[1]*LATTICE_PARAMETER
  LZ               = L_LATT_PARAM_SIGMA[2]*LATTICE_PARAMETER
  LX               = V_TOT/(LY*LZ)
  L                = Real3D(LX,LY,LZ)
  CUTSS            = (1.75*LATTICE_PARAMETER)+SIGMASS
  CUTSW            = (1.75*LATTICE_PARAMETER)+SIGMASW
  PREFAC           = [4 if (L_LATT_PARAM_SIGMA[l]>=7) else (3 if (L_LATT_PARAM_SIGMA[l]>=5) else L_LATT_PARAM_SIGMA[l]) for l in xrange(2)]
  BINSIZES         = [LATTICE_PARAMETER,PREFAC[0]*LATTICE_PARAMETER,PREFAC[1]*LATTICE_PARAMETER]
elif STRUCTURE == "111":
  LATTICE_PARAMETER= NN_DIST_XTAL
  LY               = L_LATT_PARAM_SIGMA[1]*LATTICE_PARAMETER
  LZ               = SQRT3*L_LATT_PARAM_SIGMA[2]*LATTICE_PARAMETER
  LX               = V_TOT/(LY*LZ)
  L                = Real3D(LX,LY,LZ)
  C_PARAM          = 0.5*SQRT3*LATTICE_PARAMETER/SQRT2
  CUTSS            = (4.0*C_PARAM)+SIGMASS
  CUTSW            = (4.0*C_PARAM)+SIGMASW
  PREFAC           = [4 if (L_LATT_PARAM_SIGMA[l]>=7) else (3 if (L_LATT_PARAM_SIGMA[l]>=5) else L_LATT_PARAM_SIGMA[l]) for l in xrange(2)]
  BINSIZES         = [2.0*C_PARAM,PREFAC[0]*LATTICE_PARAMETER,PREFAC[1]*SQRT3*LATTICE_PARAMETER]
else:
  print "# ERROR [parameter_setting]: Structure type '%s' is not known!" % STRUCTURE
  sys.exit()
HALFL              = Real3D(0.5*L[0],0.5*L[1],0.5*L[2])
V_TOT              = L[0]*L[1]*L[2]
BOX                = (L[0],L[1],L[2])

  
CMUMD_SWITCH       = {"solute": 0, "solvent": 1, "both": 2, "off": -1}
if CMUMD_MODE not in CMUMD_SWITCH:
  print "# ERROR [parameter_setting]: CMUMD mode '%s' is not known!" % CMUMD_MODE
  sys.exit()
if CMUMD_MODE == "off":
  SIZETR           = 0.0
  SIZEFR           = 0.0
  FORCE_DIST       = 0.0
MINSIZERES         = max([SIGMASS,SIGMASW,SIGMAWW]) if CMUMD_MODE == "off" else 0.5*SIZECHECKR*max([SIGMASS,SIGMASW,SIGMAWW])
INV_SHAPE_PARAM    = 1.0/SHAPE_PARAM
SIZECR             = FORCE_DIST-SIZETR-(0.5*SIZEFR)
#while CMUMD_EVERY%SIM_ISTEPS != 0:
  #CMUMD_EVERY     += 1
#CMUMD_EVERY_SIMSTEP= int(CMUMD_EVERY/SIM_ISTEPS)
#DT_CMUMD           = float(CMUMD_EVERY)*DT_SIM


WARMUP_NLOOPS         = TOTAL_WARMUP_STEPS/WARMUP_ISTEPS
EQUIL_NLOOPS          = TOTAL_EQUIL_STEPS/EQUIL_ISTEPS
SIMWARM_NLOOPS        = TOTAL_SIMWARM_STEPS/SIMWARM_ISTEPS
SIM_NLOOPS            = TOTAL_SIM_STEPS/SIM_ISTEPS
while (DATA_EVERY%EQUIL_ISTEPS!=0) or (DATA_EVERY%SIM_ISTEPS!=0):
  DATA_EVERY         += 1
DATA_EVERY_EQSTEP     = int(DATA_EVERY/EQUIL_ISTEPS)
DATA_EVERY_SIMSTEP    = int(DATA_EVERY/SIM_ISTEPS)
DT_DATA               = float(DATA_EVERY)*DT_SIM
PRINT_DATA_HEADER     = DATA_EVERY*50
while CONFIG_EVERY_EQ%EQUIL_ISTEPS!=0:
  CONFIG_EVERY_EQ    += 1
CONFIG_EVERY_EQSTEP   = int(CONFIG_EVERY_EQ/EQUIL_ISTEPS)
while CONFIG_EVERY_SIM%SIM_ISTEPS!=0:
  CONFIG_EVERY_SIM   += 1
CONFIG_EVERY_SIMSTEP  = int(CONFIG_EVERY_SIM/SIM_ISTEPS)
while (RESTART_EVERY%EQUIL_ISTEPS!=0) or (RESTART_EVERY%SIM_ISTEPS!=0):
  RESTART_EVERY      += 1
RESTART_EVERY_EQSTEP  = int(RESTART_EVERY/EQUIL_ISTEPS)PISTON_MASS_SIM
RESTART_EVERY_SIMSTEP = int(RESTART_EVERY/SIM_ISTEPS)
DT_CONFIG_EQ          = float(CONFIG_EVERY_EQ)*DT_SIM
DT_CONFIG_SIM         = float(CONFIG_EVERY_SIM)*DT_SIM
DT_RESTART            = float(RESTART_EVERY)*DT_SIM


WARMUP_GAMMA       = 10.0
EPSILON_DELTASS    = (EPSILONSS-EPSILON_START)/WARMUP_NLOOPS
EPSILON_DELTASW    = (EPSILONSW-EPSILON_START)/WARMUP_NLOOPS
EPSILON_DELTAWW    = (EPSILONWW-EPSILON_START)/WARMUP_NLOOPS
CAPRADIUSSS        = 0.6*SIGMASS
CAPRADIUSSW        = 0.6*SIGMASW
CAPRADIUSWW        = 0.6*SIGMAWW


CTIMELJ            = 0
CTIMESI            = 0
MODUS             += ",UNSET"


########################################################################
# function definitions                                                 #
########################################################################

def write_logfile():
  global SEPARATOR
  
  SEPARATOR              = "# ######################################################################################################################################################"
  
  rescaledData           = {"Lx": "Lx:  %.4f"%BOX[0], "V_tot": "%.4f"%V_TOT, "dtWarm": "%.4f"%DT_WARMUP, "dtSim": "%.4f"%DT_SIM, "temp": "%.4f"%TEMPERATURE, "damp_eq": "%.4f"%DAMPING_EQ, "damp_sim": "%.4f"%DAMPING_SIM, 
			    "press": "%.4f"%PRESSURE, "pm_baroEq": "%.4f"%PISTON_MASS_EQ, "pm_baroSim": "%.4f"%PISTON_MASS_SIM, "nnDist": "%.3f"%NN_DIST_XTAL, "lParam": "%.3f"%LATTICE_PARAMETER, 
			    "kxtal": "%.4f"%K_XTAL, "epsSS": "eps:  %.4f"%EPSILONSS, "epsSW": "eps:  %.4f"%EPSILONSW, "epsWW": "eps:  %.4f"%EPSILONWW, "mS": "%.4f"%MASSS, "mW": "%.4f"%MASSW, 
			    "sizeTR": "%.4f"%SIZETR, "sizeCR": "%.4f"%SIZECR, "sizeFR": "%.4f"%SIZEFR, "forceD": "%.4f"%FORCE_DIST, "sizeCheckR": "%.4f"%SIZECHECKR, "forceKS": "%.4f"%FORCE_CONST_S, "forceKW": "%.4f"%FORCE_CONST_W, "shape": "%.4f"%SHAPE_PARAM, 
			    "warmStep": "%.4f"%(TOTAL_WARMUP_STEPS*DT_WARMUP), "eqStep": "%.4f"%(TOTAL_EQUIL_STEPS*DT_SIM), "simwarmStep": "%.4f"%(TOTAL_SIMWARM_STEPS*DT_SIM), "simStep": "%.4f"%(TOTAL_SIM_STEPS*DT_SIM), 
			    "dtDat": "%.4f"%DT_DATA, "dtConfEq": "%.4f"%DT_CONFIG_EQ, "dtConfSim": "%.4f"%DT_CONFIG_SIM, "dtRest": "%.4f"%DT_RESTART, "warmCut": "%.3f"%WARM_CUTOFF, "simCut": "%.3f"%SIM_CUTOFF, 
			    "cellCut": "%.3f"%CELL_CUTOFF, "skin": "%.3f"%SKIN, "NDens": "%.4f"%NUMDENS}

  sep                    = len(max(rescaledData.values(),key=len))+1
  for key in rescaledData:
    for j in xrange(sep-len(rescaledData[key])):
      rescaledData[key] += " "
 
  print SEPARATOR
  print "#", espresso.Version().info(), "| LJ_CMUMD.py Script Version of", time.strftime("%b %d %Y, %H:%M:%S",time.gmtime())
  print SEPARATOR

  if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
    print "# ensemble for equilibration         = ", ENSEMBLES["eq"]
  print "# ensemble for simulation            = ", ENSEMBLES["sim"]

  if ("start" in MODUS) or (MODUS=="auto,UNSET"):
    if "p" in ENSEMBLES["eq"]:
      print "# initial box edge lengths           =  %s( %8.2f nm    )  Ly:  %.4f ( %.2f nm )  Lz:  %.4f ( %.2f nm )" % (rescaledData["Lx"],BOX[0]*SI_UNITS["l"],BOX[1],BOX[1]*SI_UNITS["l"],BOX[2],BOX[2]*SI_UNITS["l"])   
      print "# initial volume                     =  %s( %8.2f nm^3  )" % (rescaledData["V_tot"],V_TOT*SI_UNITS["V"])
      print "# initital number density            =  %s( %8.2f nm^-3 )" % (rescaledData["NDens"],NUMDENS/SI_UNITS["V"])
    else:
      print "# box edge lengths                   =  %s( %8.2f nm    )  Ly:  %.4f ( %.2f nm )  Lz:  %.4f ( %.2f nm )" % (rescaledData["Lx"],BOX[0]*SI_UNITS["l"],BOX[1],BOX[1]*SI_UNITS["l"],BOX[2],BOX[2]*SI_UNITS["l"])   
      print "# volume                             =  %s( %8.2f nm^3  )" % (rescaledData["V_tot"],V_TOT*SI_UNITS["V"])
      print "# number density                     =  %s( %8.2f nm^-3 )" % (rescaledData["NDens"],NUMDENS/SI_UNITS["V"])
    print "# mole fraction solute               =  %.4f" % MOLFRAC
    print "# init. total number of particles    = ", NPART_TOT
    print "# init. number of solvent particles  = ", NPART_SOLVENT
    print "# init. number of solute particles   = ", NPART_SOLUTE
  if "EQUILIBRATION" in MODUS:
    if "p" in ENSEMBLES["eq"]:
      print "# initial box edge lengths           =  %s( %8.2f nm    )  Ly:  %.4f ( %.2f nm )  Lz:  %.4f ( %.2f nm )" % (rescaledData["Lx"],BOX[0]*SI_UNITS["l"],BOX[1],BOX[1]*SI_UNITS["l"],BOX[2],BOX[2]*SI_UNITS["l"])   
      print "# initial volume                     =  %s( %8.2f nm^3  )" % (rescaledData["V_tot"],V_TOT*SI_UNITS["V"])
      print "# initital number density            =  %s( %8.2f nm^-3 )" % (rescaledData["NDens"],NUMDENS/SI_UNITS["V"])
    else:
      print "# box edge lengths                   =  %s( %8.2f nm    )  Ly:  %.4f ( %.2f nm )  Lz:  %.4f ( %.2f nm )" % (rescaledData["Lx"],BOX[0]*SI_UNITS["l"],BOX[1],BOX[1]*SI_UNITS["l"],BOX[2],BOX[2]*SI_UNITS["l"])   
      print "# volume                             =  %s( %8.2f nm^3  )" % (rescaledData["V_tot"],V_TOT*SI_UNITS["V"])
      print "# number density                     =  %s( %8.2f nm^-3 )" % (rescaledData["NDens"],NUMDENS/SI_UNITS["V"])
    print "# mole fraction solute               =  %.4f" % (float(NPART_SOLUTE)/float(NPART_LIQ))
    print "# init. total number of particles    = ", NPART_TOT
    print "# init. number of solvent particles  = ", NPART_SOLVENT
    print "# init. number of solute particles   = ", NPART_SOLUTE
  if "SIMULATION" in MODUS:
    print "# box edge lengths                   =  %s( %8.2f nm    )  Ly:  %.4f ( %.2f nm )  Lz:  %.4f ( %.2f nm )" % (rescaledData["Lx"],BOX[0]*SI_UNITS["l"],BOX[1],BOX[1]*SI_UNITS["l"],BOX[2],BOX[2]*SI_UNITS["l"])   
    print "# volume                             =  %s( %8.2f nm^3  )" % (rescaledData["V_tot"],V_TOT*SI_UNITS["V"])
    print "# number density                     =  %s( %8.2f nm^-3 )" % (rescaledData["NDens"],NUMDENS/SI_UNITS["V"])
    print "# mole fraction solute               =  %.4f" % (float(NPART_SOLUTE)/float(NPART_LIQ))
    print "# total number of particles          = ", NPART_TOT
    print "# particles in substrate             = ", NPART_XTAL
    print "# particles in solution              = ", NPART_LIQ
    print "# number of solvent particles        = ", NPART_SOLVENT
    print "# number of solute particles         = ", NPART_SOLUTE
  if ("start" in MODUS) or (MODUS=="auto,UNSET"):PISTON_MASS_SIM
    print "# dt                  (warmup)       =  %s( %8.2f fs     )" % (rescaledData["dtWarm"],DT_WARMUP*SI_UNITS["t"]*1000)  
  print "# dt                  (sim)          =  %s( %8.2f fs     )" % (rescaledData["dtSim"],DT_SIM*SI_UNITS["t"]*1000)
  print "# temperature                        =  %s( %8.2f K      )" % (rescaledData["temp"],TEMPERATURE*SI_UNITS["T"])
  if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
    print "# damping constant    (eq)           =  %s( %8.2f 1/ps   )" % (rescaledData["damp_eq"],DAMPING_EQ/SI_UNITS["t"])
  print "# damping constant    (sim)          =  %s( %8.2f 1/ps   )" % (rescaledData["damp_sim"],DAMPING_SIM/SI_UNITS["t"])
  if ("p" in ENSEMBLES["eq"]) or ("p" in ENSEMBLES["sim"]):
    print "# pressure                           =  %s( %8.2f bar    )" % (rescaledData["press"],PRESSURE*SI_UNITS["P"])
    if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
      print "# piston mass         (eq)           =  %s( %8.2f ps.bar )" % (rescaledData["pm_baroEq"],PISTON_MASS_EQ*SI_UNITS["t"]*SI_UNITS["P"])
    print "# piston mass         (sim)          =  %s( %8.2f ps.bar )" % (rescaledData["pm_baroSim"],PISTON_MASS_SIM*SI_UNITS["t"]*SI_UNITS["P"])
  if ("p" in ENSEMBLES["eq"]) or ("p" in ENSEMBLES["sim"]):
    if (PRESSURE_SSUNITS>0.0) and (PRESSURE_SSUNITS<1.0):
      print "# saturation index                   =  %.4f" % SATURATIONINDEX
  print SEPARATOR
  print "# xtal structure                     = ", STRUCTURE[0]+","+STRUCTURE[1]+","+STRUCTURE[2]
  print "# distance to nearest neighbour      =  %s( %8.2f sigma = %.3f nm )" % (rescaledData["nnDist"],NN_DIST_XTAL/SIGMAWW,NN_DIST_XTAL*SI_UNITS["l"])
  if STRUCTURE == "100":
    print "# lattice constant                   =  %s( %8.2f sigma = %.3f nm )" % (rescaledData["lParam"],SQRT2*NN_DIST_XTAL/SIGMAWW,LATTICE_PARAMETER*SI_UNITS["l"])
  if STRUCTURE == "111":
    print "# lattice constant                   =  %s( %8.2f sigma = %.3f nm )" % (rescaledData["lParam"],NN_DIST_XTAL/SIGMAWW,LATTICE_PARAMETER*SI_UNITS["l"])
  print "# spring constant between xtal parts =  %s( %8.2f kJ/mol.nm^2 = %8.2f K.kB/nm^2 )" % (rescaledData["kxtal"],K_XTAL*SI_UNITS["E"]/SI_UNITS["l"]**2,K_XTAL*SI_UNITS["T"]/SI_UNITS["l"]**2)
  print "# number of bonds in xtal            = ", NBONDS
  print SEPARATOR
  print "# WW interaction                     =  %s( %8.2f K.kB = %4.2f kJ/mol )  sig:  %.4f ( %.3f nm )  cut:  %.4f ( %.3f nm )  type:  %s" % (rescaledData["epsWW"],EPSILONWW*SI_UNITS["T"],EPSILONWW*SI_UNITS["E"],
																	  SIGMAWW,SIGMAWW*SI_UNITS["l"],CUTOFFWW,CUTOFFWW*SI_UNITS["l"],WWTYPEWW)
  print "# SW interaction      ( %s | %s )  =  %s( %8.2f K.kB = %4.2f kJ/mol )  sig:  %.4f ( %.3f nm )  cut:  %.4f ( %.3f nm )  type:  %s" % (SWPARAMSETUPE,SWPARAMSETUPS,rescaledData["epsSW"],EPSILONSW*SI_UNITS["T"],
																	EPSILONSW*SI_UNITS["E"],SIGMASW,SIGMASW*SI_UNITS["l"],CUTOFFSW,CUTOFFSW*SI_UNITS["l"],
																	WWTYPESW)
  print "# SS interaction                     =  %s( %8.2f K.kB = %4.2f kJ/mol )  sig:  %.4f ( %.3f nm )  cut:  %.4f ( %.3f nm )  type:  %s" % (rescaledData["epsSS"],EPSILONSS*SI_UNITS["T"],EPSILONSS*SI_UNITS["E"],
																	  SIGMASS,SIGMASS*SI_UNITS["l"],CUTOFFSS,CUTOFFSS*SI_UNITS["l"],WWTYPESS)
  print "# particle mass solvent              =  %s( %8.2f g/mol )" % (rescaledData["mW"],MASSW*SI_UNITS["m"])
  print "# particle mass solute               =  %s( %8.2f g/mol )" % (rescaledData["mS"],MASSS*SI_UNITS["m"])
  print SEPARATOR
  print "# CMuMD mode                         = ", CMUMD_MODE
  if CMUMD_MODE != "off":
    print "# size of transition region          =  %s( %6.2f nm )" % (rescaledData["sizeTR"],SIZETR*SI_UNITS["l"])
    print "# size of control region             =  %s( %6.2f nm )" % (rescaledData["sizeCR"],SIZECR*SI_UNITS["l"])
    print "# size of force region               =  %s( %6.2f nm )" % (rescaledData["sizeFR"],SIZEFR*SI_UNITS["l"])
    print "# distance of force region           =  %s( %6.2f nm )" % (rescaledData["forceD"],FORCE_DIST*SI_UNITS["l"])
    print "# size of adapt. checker region      =  %s( %6.2f nm )" % (rescaledData["sizeCheckR"],SIZECHECKR*SI_UNITS["l"])
    print "# target mole fraction solute        =  %.4f" % EXP_MOLF
    print "# membrane force constant solvent    =  %s( %6.2f kJ.nm^3/mol = %8.2f K.kB.nm^3 )" % (rescaledData["forceKW"],FORCE_CONST_W*SI_UNITS["E"]*SI_UNITS["V"],FORCE_CONST_W*SI_UNITS["T"]*SI_UNITS["V"])
    print "# membrane force constant solute     =  %s( %6.2f kJ.nm^3/mol = %8.2f K.kB.nm^3 )" % (rescaledData["forceKS"],FORCE_CONST_S*SI_UNITS["E"]*SI_UNITS["V"],FORCE_CONST_S*SI_UNITS["T"]*SI_UNITS["V"])
    print "# membrane shape parameter           =  %s( %6.2f nm )" % (rescaledData["shape"],SHAPE_PARAM*SI_UNITS["l"])
  print SEPARATOR
  if ("start" in MODUS) or (MODUS=="auto,UNSET"):
    print "# warmup nloops                      = ", WARMUP_NLOOPS
    print "# warmup isteps                      = ", WARMUP_ISTEPS
    print "# total warmup time                  =  %s( %7d dt = %8.2f ps )" % (rescaledData["warmStep"],TOTAL_WARMUP_STEPS,float(TOTAL_WARMUP_STEPS)*DT_WARMUP*SI_UNITS["t"]) 
  if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
    print "# equilibration nloops               = ", EQUIL_NLOOPS
    print "# equilibration isteps               = ", EQUIL_ISTEPS
    print "# total equilibration time           =  %s( %7d dt = %8.2f ps )" % (rescaledData["eqStep"],TOTAL_EQUIL_STEPS,float(TOTAL_EQUIL_STEPS)*DT_SIM*SI_UNITS["t"])
    print "# sim.-preceding warmup nloops       = ", SIMWARM_NLOOPS
    print "# sim.-preceding warmup isteps       = ", SIMWARM_ISTEPS
    print "# total sim.-preceding warmup time   =  %s( %7d dt = %8.2f ps )" % (rescaledData["simwarmStep"],TOTAL_SIMWARM_STEPS,float(TOTAL_SIMWARM_STEPS)*DT_SIM*SI_UNITS["t"]) 
  print "# simulation nloops                  = ", SIM_NLOOPS
  print "# simulation isteps                  = ", SIM_ISTEPS
  print "# total simulation time              =  %s( %7d dt = %8.2f ps )" % (rescaledData["simStep"],TOTAL_SIM_STEPS,float(TOTAL_SIM_STEPS)*DT_SIM*SI_UNITS["t"]) 
  print "# time between data output           =  %s( %7d dt = %8.2f ps )" % (rescaledData["dtDat"],DATA_EVERY,DT_DATA*SI_UNITS["t"])
  if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
    print "# time between configurations (eq)   =  %s( %7d dt = %8.2f ps )" % (rescaledData["dtConfEq"],CONFIG_EVERY_EQ,DT_CONFIG_EQ*SI_UNITS["t"])
  print "# time between configurations (sim)  =  %s( %7d dt = %8.2f ps )" % (rescaledData["dtConfSim"],CONFIG_EVERY_SIM,DT_CONFIG_SIM*SI_UNITS["t"])
  print "# time between restart files         =  %s( %7d dt = %8.2f ps )" % (rescaledData["dtRest"],RESTART_EVERY,DT_RESTART*SI_UNITS["t"])
  #print "# time between CMuMD application     =  %s( %7d dt = %8.2f ps )" % (rescaledData["dtCMuMD"],CMUMD_EVERY,DT_CMUMD*SI_UNITS["t"])
  print "# units of data output               = ", DATA_UNITS
  print SEPARATOR
  print "# NCPUs                              = ", NCPUS
  print "# nodeGrid            (%s)          =  (%d, %d, %d)" % (NODESETUP,NODEGRID[0],NODEGRID[1],NODEGRID[2])
  print "# cellGrid                           = ", CELLGRID
  if ("start" in MODUS) or (MODUS=="auto,UNSET"):
    print "# largest cutoff      (warmup)       =  %s( %6.2f nm )" % (rescaledData["warmCut"],WARM_CUTOFF*SI_UNITS["l"])
  print "# largest cutoff      (sim)          =  %s( %6.2f nm )" % (rescaledData["simCut"],SIM_CUTOFF*SI_UNITS["l"])
  print "# cell cutoff                        =  %s( %6.2f nm )" % (rescaledData["cellCut"],CELL_CUTOFF*SI_UNITS["l"])
  print "# skin                               =  %s( %6.2f nm )" % (rescaledData["skin"],SKIN*SI_UNITS["l"])
  print "# seed for RNG        (%s)          =  %d" % (SEEDSETUP,SEED)
  print SEPARATOR


def init_xtal():
  global XTALBONDS_LIST
  global NPART_XTAL
  global NBONDS
  
  pid                   = 0
  
  ##################
  ##################
  ##################  genrate Xtal surface of given lattice parameter 
  ##################
  ##################
   
  #print "# setting up xtal surface ..."
  
  ### fcc 100 face
  if STRUCTURE == "100":
    
    num                 = Int3D(4,L_LATT_PARAM_SIGMA[1],L_LATT_PARAM_SIGMA[2])

    lattice_vec1        = Real3D(0,0,LATTICE_PARAMETER)    
    lattice_vec2        = Real3D(0,0.5*LATTICE_PARAMETER,0.5*LATTICE_PARAMETER)
    lattice_vec3        = Real3D(0.5*LATTICE_PARAMETER,0.5*LATTICE_PARAMETER,0)
    
    for k in xrange(-num[0],num[0]):
      for j in xrange(0,2*num[1]):
	for i in xrange(0,num[2]):
	  c_Pos         = (i*lattice_vec1)+(j*lattice_vec2)+(k*lattice_vec3)
	  pid          += 1
	  pos           = Real3D(c_Pos[0]+(0.25*LATTICE_PARAMETER),c_Pos[1],c_Pos[2])
	  v             = Real3D(0,0,0)
	  part          = [pid, pos, 0, v, MASSS]
	  XTALPART_LIST.append(part)

  ### hcp face	  
  if STRUCTURE == "111":
    
    num                 = Int3D(4,L_LATT_PARAM_SIGMA[1],L_LATT_PARAM_SIGMA[2])
    
    lattice_vec1        = Real3D(0,0,SQRT3*LATTICE_PARAMETER)
    lattice_vec2        = Real3D(0,0.5*LATTICE_PARAMETER,SQRT3*0.5*LATTICE_PARAMETER)
    lattice_vec3        = Real3D(C_PARAM,0.5*LATTICE_PARAMETER,0.5*LATTICE_PARAMETER/SQRT3)    
    
    for k in xrange(-num[0],num[0]+1):
      for j in xrange(0,2*num[1]):
	for i in xrange(0,num[2]):
	  c_Pos         = (i*lattice_vec1)+(j*lattice_vec2)+(k*lattice_vec3)
	  pid          += 1
	  pos           = Real3D(c_Pos[0],c_Pos[1],c_Pos[2])
	  v             = Real3D(0,0,0)
	  part          = [pid, pos, 0, v, MASSS]
	  XTALPART_LIST.append(part)
		
  NPART_XTAL            = pid
  
  if pid != NPART_XTAL_SOLL:
    print str(pid)+"  "+str(NPART_XTAL_SOLL)
    print "# ERROR [init_xtal]: Different number of particles in substrate than expected!"
    sys.exit()
    
  system.storage.addParticles(XTALPART_LIST, 'id', 'pos', 'type', 'v', 'mass')
  system.storage.decompose()
  
  #### DEBUGGING #############################################################
  #espresso.tools.writexyz("debug_xtal.xyz", system, velocities = False, unfolded = False, append = False)
  #sys.exit()
  #### DEBUGGING #############################################################
  
  ##################
  ##################
  ##################  creating bonds between neighbors in Xtal
  ##################
  ##################

  #print "# initializing bonds in xtal ..."
    
  xtalVerletlist      = espresso.VerletList(system, NN_DIST_XTAL+0.01)
  for z in xrange(len(xtalVerletlist.getAllPairs())):
    if len(xtalVerletlist.getAllPairs()[z]):
      XTALBONDS_LIST += xtalVerletlist.getAllPairs()[z]
  NBONDS              = len(XTALBONDS_LIST)
  del xtalVerletlist
  
  system.storage.removeAllParticles()
  

def setup_system():
  pid            = 0
  pList          = []

  ##################
  ##################
  ##################  setup the liquid as AT system
  ##################
  ##################

  #print "# inserting liquid particles ..."
  
  for k in xrange(NPART_SOLUTE):
    pos          = system.bc.getRandomPos()
    pid         += 1
    v            = Real3D(0,0,0)
    part         = [pid, pos, 0, v, MASSS]
    pList.append(part)
  for k in xrange(NPART_SOLVENT):
    pos          = system.bc.getRandomPos()
    pid         += 1
    v            = Real3D(0,0,0)
    part         = [pid, pos, 1, v, MASSW]
    pList.append(part)
  
  system.storage.addParticles(pList, 'id', 'pos', 'type', 'v', 'mass')
  
  if pid != NPART_TOT:
    print "# ERROR [setup_system]: pid != total number of parts!"
    sys.exit()


def reset_system():
  global NPART_SOLUTE
  global NPART_SOLVENT
  global NPART_TOT
  global NPART_LIQ  
  
  halfLX                = 0.5*system.bc.boxL[0]
  pList                 = []

  for pid in xrange(NPART_XTAL):
    cpid                = XTALPART_LIST[pid][0]
    cpos                = XTALPART_LIST[pid][1]
    ctype               = XTALPART_LIST[pid][2]
    cvel                = XTALPART_LIST[pid][3]
    cmass               = XTALPART_LIST[pid][4]

    cpos[0]            += halfLX
        
    pList.append([cpid, cpos, ctype, cvel, cmass])
  
  ##################
  ##################
  ##################  setup the liquid as AT system
  ##################
  ##################

  #print "# inserting liquid particles ..."
  
  pid                   = NPART_XTAL
  
  NPART_SOLUTE          = 0
  NPART_SOLVENT         = 0
  minDistSS2            = CUTSS*CUTSS
  minDistSW2            = CUTSW*CUTSW
  for k in xrange(NPART_TOT):
    cpart               = system.storage.getParticle(k+1)
    cpos                = cpart.pos
    cDist2              = (cpos[0]-halfLX)**2
    if cpart.type == 0:
      if cDist2 >= minDistSS2:
        pid            += 1
	part            = [pid, cpart.pos, cpart.type, cpart.v, cpart.mass]
      	NPART_SOLUTE   += 1
	pList.append(part)
	FREEPARTS_LIST.append(pid)
    if cpart.type == 1:
      if cDist2 >= minDistSW2:
        pid            += 1
	part            = [pid, cpart.pos, cpart.type, cpart.v, cpart.mass]
	NPART_SOLVENT  += 1
	pList.append(part)
	FREEPARTS_LIST.append(pid)
	
  system.storage.removeAllParticles()
  
  NPART_TOT             = len(pList)
  system.storage.addParticles(pList, 'id', 'pos', 'type', 'v', 'mass')
  
  for cid in (cp[0] for cp in pList):
    parti               = system.storage.getParticle(cid)
    partiPOS            = system.bc.getFoldedPosition(parti.pos, parti.imageBox)
    system.storage.modifyParticle(cid, "pos", partiPOS[0])
  
  NPART_LIQ             = NPART_TOT-NPART_XTAL
  if NPART_TOT != (NPART_XTAL+NPART_SOLUTE+NPART_SOLVENT):
    print "# ERROR [reset_system]: Npart_tot != total number of parts!"
    sys.exit()
    
  if len(FREEPARTS_LIST) != NPART_LIQ:
    print "# ERROR [reset_system]: length of freeParts_list != total number of parts in solution!"
    sys.exit()
 
   
def read_in_restart():
  global NPART_XTAL
  global NPART_SOLUTE
  global NPART_SOLVENT
  global NPART_LIQ
  global NPART_TOT
  global WALL_INSTANCE
  #global STRUCTURE_L
  #global STRUCTURE_R
  #global BINS_INPLANE
  global MODUS
  global CTIMELJ
  global CTIMESI
  global V_TOT
  global BOX
  global NUMDENS
  global NUMDENS0
  global NUMDENS1
  global NBONDS
  
  NPART_XTAL                    = 0
  NPART_SOLUTE                  = 0
  NPART_SOLVENT                 = 0
  WALL_INSTANCE                 = Real3D(0.0,0.0,0.0)
  #STRUCTURE_L                   = []
  #STRUCTURE_R                   = []
  #BINS_INPLANE                  = [0,0]
  
  if os.path.isfile("restart.in.gz"):
    infile                      = gzip.open("restart.in.gz",'r')
    for line in infile:
      line                      = line.strip('\n')
      split                     = [value for value in line.split()]
      
      if split[0] == "#":
	mode                    = split[1]+split[2]
	if mode == "TIMESTEP:":
	  MODUS                 = "cont,%s"%split[5][1:-1]
	continue
	
      if mode == "TIMESTEP:":
	cTimeLJ                 = split[0]
	if split[0] != "first":
	  cTimeLJ               = float(split[0])
	CTIMELJ                 = cTimeLJ
	cTimeSI                 = split[1]
	if split[1] != "first":
	  cTimeSI               = float(split[1])
	CTIMESI                 = cTimeSI  
      elif mode == "PARTNUM:":
	NPART_TOT               = int(split[0])
      elif mode == "BOXDIM:":
	cL                      = Real3D(float(split[0]),float(split[1]),float(split[2]))
	V_TOT                   = cL[0]*cL[1]*cL[2]
	BOX                     = (cL[0],cL[1],cL[2])
      elif mode == "NUMDENS:":
	NUMDENS0                = float(split[0])
	NUMDENS1                = float(split[1])
	NUMDENS                 = NUMDENS0+NUMDENS1
      elif mode == "WALLINSTANCE:":
      	WALL_INSTANCE[0]        = float(split[0])
	WALL_INSTANCE[1]        = float(split[1])
	WALL_INSTANCE[2]        = float(split[2])
	#if split[0][-1] == "L":
	  #for bz in xrange(1,len(split)):
	    #STRUCTURE_L.append(float(split[bz]))
	#if split[0][-1] == "R":
	  #for bz in xrange(1,len(split)):
	    #STRUCTURE_R.append(float(split[bz]))
      elif mode == "ATOMPROPS:":
	if "EQUILIBRATION" in MODUS:
	  part                  = [int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),int(split[1]), 
			    Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
          ALLPARTS_LIST.append(part)
	  if split[1] == "0":
	    NPART_SOLUTE       += 1
	  if split[1] == "1":
	    NPART_SOLVENT      += 1
	elif "SIMULATION" in MODUS:
	  if split[1] == "-1":
	    part                = [int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),0, 
			    Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
	    ALLPARTS_LIST.append(part)
	    XTALPART_LIST.append(part)
	    NPART_XTAL         += 1
	  elif split[1] == "0":
	    part                = [int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),int(split[1]), 
			    Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
	    ALLPARTS_LIST.append(part)
	    FREEPARTS_LIST.append(int(split[0]))
	    NPART_SOLUTE       += 1
	  elif split[1] == "1":
	    part                = [int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),int(split[1]), 
			    Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
	    ALLPARTS_LIST.append(part)
	    FREEPARTS_LIST.append(int(split[0]))
	    NPART_SOLVENT      += 1
	  else:
	    print "# ERROR [readin_restart]: Particle type '%s' is not known!" % str(split[1])
	    sys.exit()
	else:
	  print "# ERROR [readin_restart]: Mode '%s' is not known!"%MODUS
	  sys.exit() 
      elif mode == "XTALBONDS:":
	XTALBONDS_LIST.append((int(split[0]),int(split[1])))
      else:
	print "# ERROR [readin_restart]: Restart file seems to be corrupted!"
	sys.exit()

    NBONDS                      = len(XTALBONDS_LIST)
    if (NPART_XTAL+NPART_SOLUTE+NPART_SOLVENT) != NPART_TOT:
      print "# ERROR [readin_restart]: Total number of particles != number of xtal particles + number of solute particles + number of solvent particles!"
      sys.exit()    
    NPART_LIQ                   = NPART_TOT-NPART_XTAL
    if NPART_LIQ != (NPART_SOLUTE+NPART_SOLVENT):
      print "# ERROR [readin_restart]: Number of liquid particles != number of solute particles + number of solvent particles!"
      sys.exit()
  else:
    print "# ERROR [readin_restart]: Restart file does not exist!"
    sys.exit()

    
def write_out_restart(myTime,cmode):
  global NUMDENS0
  global NUMDENS1
  
  if myTime == "first":
    restartFile  = gzip.open("first_Sim_Restart.in.gz",'w')
    restartFile.write("# TIME STEP: timeLJ timeSI (%s)\n"%cmode)
    restartFile.write("  first first\n")
  else:
    restartFile  = gzip.open("restart.in.gz",'w')
    restartFile.write("# TIME STEP: timeLJ timeSI (%s)\n"%cmode)
    restartFile.write("  %4.6f %4.6f\n"%(myTime+CTIMELJ,(myTime+CTIMELJ)*SI_UNITS["t"]))
  restartFile.write("# PART NUM: N\n")
  restartFile.write("  %d\n"%NPART_TOT)
  restartFile.write("# BOX DIM: lx ly lz\n")
  restartFile.write("  %4.6f %4.6f %4.6f\n"%(system.bc.boxL[0],system.bc.boxL[1],system.bc.boxL[2]))
  restartFile.write("# NUM DENS: solute solvent\n")
  if cmode != "SIMULATION":
    NUMDENS0     = float(NPART_SOLUTE)/(system.bc.boxL[0]*system.bc.boxL[1]*system.bc.boxL[2])
    NUMDENS1     = float(NPART_SOLVENT)/(system.bc.boxL[0]*system.bc.boxL[1]*system.bc.boxL[2])
  restartFile.write("  %4.6f %4.6f\n"%(NUMDENS0,NUMDENS1))
  if cmode == "SIMULATION":
    #global BINS_INPLANE
    #BINS_INPLANE = constMuMD.binsInPlane
    #global STRUCTURE_L
    #STRUCTURE_L  = constMuMD.structureL
    #global STRUCTURE_R
    #STRUCTURE_R  = constMuMD.structureR
    #restartFile.write("# WALL INSTANCE: center=%4.6f lborder=%4.6f rborder=%4.6f binsY=%d binsZ=%d\n"%(constMuMD.substrate[0],constMuMD.substrate[1],constMuMD.substrate[2],BINS_INPLANE[0],BINS_INPLANE[1]))
    restartFile.write("# WALL INSTANCE: center lborder rborder\n")
    restartFile.write("  %4.6f %4.6f %4.6f\n"%(constMuMD.substrate[0],constMuMD.substrate[1],constMuMD.substrate[2]))
    #for by in xrange(BINS_INPLANE[0]):
      #wline    = "  %4d,L"%by
      #for bz in xrange(BINS_INPLANE[1]):
	#myIndex= (by*BINS_INPLANE[1])+bz
	#wline += " %4.6f"%STRUCTURE_L[myIndex]
      #restartFile.write("%s\n"%wline)
    #for by in xrange(BINS_INPLANE[0]):
      #wline    = "  %4d,R"%by
      #for bz in xrange(BINS_INPLANE[1]):
	#myIndex= (by*BINS_INPLANE[1])+bz
	#wline += " %4.6f"%STRUCTURE_R[myIndex]
      #restartFile.write("%s\n"%wline)
  restartFile.write("# ATOM PROPS: pid type x y z vx vy vz m\n")
  
  snapshot     = [system.storage.getParticle(pID) for pID in xrange(1,NPART_TOT+1)]
  if cmode == "EQUILIBRATION":
    for cPart in snapshot[:]:
      restartFile.write("  %4d %d %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f\n"
			%(cPart.id,cPart.type,cPart.pos[0],cPart.pos[1],cPart.pos[2],cPart.v[0],cPart.v[1],cPart.v[2],cPart.mass))    
  if cmode == "SIMULATION":
    for cPart in snapshot[:NPART_XTAL]:
      restartFile.write("  %4d %d %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f\n"
			%(cPart.id,-1,cPart.pos[0],cPart.pos[1],cPart.pos[2],cPart.v[0],cPart.v[1],cPart.v[2],cPart.mass))
    for cPart in snapshot[NPART_XTAL:]:
      restartFile.write("  %4d %d %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f\n"
			%(cPart.id,cPart.type,cPart.pos[0],cPart.pos[1],cPart.pos[2],cPart.v[0],cPart.v[1],cPart.v[2],cPart.mass))
    
    restartFile.write("# XTAL BONDS: id_1 id_2\n")
    for cbid in XTALBONDS_LIST:
      restartFile.write("  %d %d\n"%(cbid[0],cbid[1]))
  
  restartFile.close()
  
  del snapshot
  

def jumpToCurrentTime():
  if (CTIMELJ>0.0) and (os.path.isfile(TRJFILE)):
    #Goto current time of trajectory file in order to attach
    str_Npart_tot     = "%d"%NPART_TOT
    lines             = gzip.open(TRJFILE).readlines()
    #linecounter       = 0
    #for lineindex, s in enumerate(lines):
      #lookupLine      = s.strip('\n')
      #lookupSplit     = [value for value in lookupLine.split()]
      #if (len(lookupSplit)==1) and (lookupSplit[0]==str_Npart_tot):
	#if linecounter == int((CTIMELJ/DT_CONFIG)+0.5):
	  #break
	#else:
	  #linecounter  += 1
    if "EQUILIBRATION" in MODUS:
      lineindex       = int(CTIMELJ/DT_CONFIG_EQ)*(NPART_TOT+2)
    if "SIMULATION" in MODUS:
      lineindex       = int(CTIMELJ/DT_CONFIG_SIM)*(NPART_TOT+2)
    linesToWrite      = lines[:lineindex]
    tfile             = gzip.open(TRJFILE,'w')
    for lines in linesToWrite:
      tfile.write("%s"%lines)
    tfile.close()
  else:
    tfile             = gzip.open(TRJFILE,'w')
    tfile.close()
    
  if (CTIMELJ>0.0) and (os.path.isfile(DATAFILENAME)):
    # Goto current time of datafile in order to attach
    str_cTimeLJ       = "%.4f"%CTIMELJ
    str_cTimeSI       = "%.4f"%CTIMESI
    lines             = open(DATAFILENAME).readlines()
    for lineindex, s in enumerate(lines):
      lookupLine      = s.strip('\n')
      lookupTime      = (lookupLine.split())[0]
      if (lookupTime==str_cTimeLJ) or (lookupTime==str_cTimeSI):
	if int(float(lookupTime)/DT_SIM)%PRINT_DATA_HEADER == 0:
	  lineindex  -= 1
	break
    linesToWrite      = lines[:lineindex]
    dfile             = open(DATAFILENAME,'w')
    for lines in linesToWrite:
      dfile.write("%s"%lines)
    dfile.write("# continued run from timestep %.4f (%.4f ps) of previous %s run\n"%(CTIMELJ,CTIMESI,MODUS[MODUS.index(",")+1:]))
    dfile.close()
  else:
    dfile             = open(DATAFILENAME,'w')
    dfile.write("# continued run from timestep %.4f (%.4f ps) of previous %s run\n"%(CTIMELJ,CTIMESI,MODUS[MODUS.index(",")+1:]))
    dfile.close()
 

def setupNodeGridManually():
  global NODEGRID
  
  Edges                  = [system.bc.boxL[0],system.bc.boxL[1],system.bc.boxL[2]]
  minEdge                = min(Edges)
  maxEdge                = max(Edges)
  i_minmaxEdges          = [Edges.index(minEdge),Edges.index(maxEdge)]
  for i in xrange(3):
    if i not in i_minmaxEdges:
      rem_i              = i
      break
    
  rc_skin                = CELL_CUTOFF+SKIN
  if rc_skin <= 0:
    print "# ERROR [manual_grid_setup]: interaction range (cutoff + skin) must be larger than 0"
    sys.exit()
    
  ijkmax                 = 3*NCPUS*NCPUS+1
  d1                     = 1
  d2                     = 1
  d3                     = 1
  for i in xrange(1,NCPUS+1):
    for j in xrange(i,NCPUS+1):
      for k in xrange(j,NCPUS+1):
	if (i*j*k==NCPUS) and (i<=int(minEdge/rc_skin)) and (j<=int(Edges[rem_i]/rc_skin)) and (i*i+j*j+k*k<ijkmax):
	  d1             = i
	  d2             = j
	  d3             = k
	  ijkmax         = i*i+j*j+k*k
  
  ind                    = i_minmaxEdges[0]
  NODEGRID[ind]          = d1
  NODEGRID[rem_i]        = d2
  ind                    = i_minmaxEdges[1]
  NODEGRID[ind]          = d3

  if (NODEGRID[0]<=0) or (NODEGRID[1]<=0) or (NODEGRID[2]<=0):
    print "# ERROR [manual_grid_setup]: It was not possible to set up the node grid manually!"
    sys.exit() 
  
  if NODEGRID[0]*NODEGRID[1]*NODEGRID[2] != NCPUS:
    print "# ERROR [manual_grid_setup]: Number of node grid elements does not match total number of CPUs!"
    sys.exit()  
  
  
def write_out_data(ens,stps):
  global NUMDENS
  global OLD_LX
  
  ELJ          = 0.0
  EBond        = 0.0
  P            = Tensor(0.0)
  ETotal       = 0.0

  datafile     = open(DATAFILENAME,'a')
  if STEP%int(PRINT_DATA_HEADER/stps) == 0:
    header     = "#     Time          T        Pxx        Pyy        Pzz"
    if "p" in ENSEMBLES[ens]:
      header  += "         Lx        Rho"
    if ens == "sim":
      if CMUMD_MODE != "off":
	header+= "   Rho_Xtal    Rho_Liq     Rho_CR    Rho_Res    MolF_CR"
      header  += "  Size_ResL  Size_ResR      EBond"
    header    += "        ELJ     ETotal\n"
    datafile.write("%s"%header)

  ELJ          = interaction.computeEnergy()/NPART_TOT
  Temp         = TEMP.compute()
  if ens == "sim":
    EBond      = xtal_interaction.computeEnergy()/NBONDS
  P            = PRESS.compute()
  ETotal       = EBond+ELJ

  outdata      = (((STEP*stps*DT_SIM)+CTIMELJ)*UNITS["t"],Temp*UNITS["T"],P[0]*UNITS["P"],P[1]*UNITS["P"],P[2]*UNITS["P"])
  fmt          = "%10.4f %10.4f %10.4f %10.4f %10.4f"
  V            = system.bc.boxL[0]*system.bc.boxL[1]*system.bc.boxL[2]
  
  if "p" in ENSEMBLES[ens]:
    NUMDENS    = NPART_TOT/V
    outdata   += (system.bc.boxL[0]*UNITS["l"],NUMDENS/UNITS["V"])  
    fmt       += " %10.4f %10.4f"
  if ens == "sim":
    d_Lx       = system.bc.boxL[0]-OLD_LX
    if CMUMD_MODE != "off":
      densdata = constMuMD.getDensityData()
      outdata += (densdata[0]/UNITS["V"],densdata[1]/UNITS["V"],densdata[2]/UNITS["V"],densdata[3]/UNITS["V"],densdata[4])
      fmt     += " %10.4f %10.4f %10.4f %10.4f %10.4f"
      OLD_LX   = system.bc.boxL[0]
    outdata   += ((constMuMD.sizeRes[0]-0.5*d_Lx)*UNITS["l"],(constMuMD.sizeRes[1]-0.5*d_Lx)*UNITS["l"],EBond*UNITS["E"])
    fmt       += " %10.4f %10.4f %10.4f"
  outdata     += (ELJ*UNITS["E"],ETotal*UNITS["E"])
  fmt         += " %10.4f %10.4f\n"

  datafile.write(fmt % outdata)
  datafile.close()


def calcSaturationIndex():
  global SATURATIONINDEX
  
  # melting line from Agrawal and Kofke for P=0.0...1.0
  fitparams              = [0.6866,0.0893,-0.0463,0.1927,-0.3618,0.2903,-0.0817]
  mline                  = fitparams[0]
  for i in xrange(1,7):
    mline               += fitparams[i]*pow(PRESSURE_SSUNITS,i)
  
  SATURATIONINDEX        = log(mline/TEMPERATURE_SSUNITS)
  
  
def BarostatChecker():
  if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
    if "p" in ENSEMBLES["eq"]:
      if PISTON_MASS_EQ < DT_SIM*PRESSURE:
	print "# WARNING [parameter_setting]: Barostat may not work as the piston mass of equilibration (%.4f) is too small!" % PISTON_MASS_EQ
	print SEPARATOR

  if "p" in ENSEMBLES["sim"]:
    PCHECKER         = [PRESSURE]
    if CMUMD_SWITCH != "off":
      PCHECKER.append(0.2*max([FORCE_CONST_S,FORCE_CONST_W]))
    if PISTON_MASS_SIM < DT_SIM*max(PCHECKER):
      print "# WARNING [parameter_setting]: Barostat may not work as the piston mass of simulation (%.4f) is too small!" % PISTON_MASS_SIM
      print SEPARATOR
      
  
def ToDo(routine):
  print "# The routine named '%s' needs to be worked over!" % routine

  while True:
    print "# Continue ('y' or 'n') ?"
    choice = raw_input("")
    if choice == "n":
      sys.exit()
    elif choice == "y":
      break
    else:
      continue
 
 
##########################################################################################################################################
##########################################################################################################################################
##																        ##
##				  AA            CCCCCCCC  TTTTTTTTTTTTTTT  III       OOOOOOO       NNNN      NN                         ##
##				 AAAA         CCC         TTTTTTTTTTTTTTT  III     OOO     OOO     NNNNN     NN                         ##
##				AA  AA       CC                 TTT        III    OO         OO    NN NNN    NN                         ##
##			       AA    AA     CC                  TTT        III   OO           OO   NN  NNN   NN                         ##
##			      AAAAAAAAAA    CC                  TTT        III   OO           OO   NN   NNN  NN                         ##
##                           AA        AA    CC                 TTT        III    OO         OO    NN    NNN NN                         ##
##                          AA          AA    CCC               TTT        III     OOO     OOO     NN     NNNNN                         ##
##                         AA            AA     CCCCCCCC        TTT        III       OOOOOOO       NN      NNNN                         ##
##																        ##
##########################################################################################################################################
##########################################################################################################################################

  
########################################################################
# 1. setup of the system, random number generator and parallelisation  #
########################################################################

if os.path.isfile("endofsim.txt"):
  print "############################################################################"
  print "## Simulation will not be performed as the file <endofsim.txt> was found! ##"
  print "############################################################################"
  sys.exit()

ALLPARTS_LIST     = []
XTALPART_LIST     = []
XTALBONDS_LIST    = []
FREEPARTS_LIST    = []

if ("cont" in MODUS) or (("auto" in MODUS) and os.path.isfile("restart.in.gz")):
  read_in_restart()
 
WARM_CUTOFF       = max([WARMUP_CUTOFFSS,WARMUP_CUTOFFSW,WARMUP_CUTOFFWW])
SIM_CUTOFF        = max([CUTOFFSS,CUTOFFSW,CUTOFFWW])
CELL_CUTOFF       = 1.3*max([SIM_CUTOFF,WARM_CUTOFF])
while (CELL_CUTOFF/WARM_CUTOFF) < 2.0:
  CELL_CUTOFF    *= 1.3
RDF_CUTOFF        = 1.2*NN_DIST_XTAL
if SEED == "aut":
  SEED		  = long(time.strftime("%Y%m%d%H%M%S"))
  SEEDSETUP       = "aut"
elif SEED == "def":
  SEED            = None
  SEEDSETUP       = "def"
else:
  SEED		  = long(SEED)
  SEEDSETUP       = "man"

system            = espresso.System()

if SEEDSETUP == "def":
  system.rng      = espresso.esutil.RNG()
else:
  system.rng      = espresso.esutil.RNG(SEED)

system.bc         = espresso.bc.OrthorhombicBC(system.rng,BOX)
system.skin       = SKIN

NCPUS             = espresso.MPI.COMM_WORLD.size
NODEGRID          = [0,0,0]
if NODESETUP == "man":
  setupNodeGridManually()
elif NODESETUP == "aut":
  NODEGRID        = espresso.tools.decomp.nodeGrid(NCPUS)
else:
  print "# ERROR [simulation_setup]: Type '%s' for node grid setup is not known!"%NODESETUP
  sys.exit() 
CELLGRID          = espresso.tools.decomp.cellGrid(BOX,NODEGRID,CELL_CUTOFF,SKIN)

system.storage    = espresso.storage.DomainDecomposition(system,NODEGRID,CELLGRID)

if ("cont" in MODUS) or ("auto" in MODUS):
  if "SIMULATION" in MODUS:
    system.storage.addParticles(ALLPARTS_LIST,'id', 'pos', 'type', 'v', 'mass')
  if "EQUILIBRATION" in MODUS:
    init_xtal()
    system.storage.addParticles(ALLPARTS_LIST,'id', 'pos', 'type', 'v', 'mass') 
if ("start" in MODUS) or (MODUS=="auto,UNSET"):
  init_xtal()
  setup_system()
  
system.storage.decompose()

if ("p" in ENSEMBLES["eq"]) or ("p" in ENSEMBLES["sim"]):
  if (PRESSURE_SSUNITS>0.0) and (PRESSURE_SSUNITS<1.0):
    calcSaturationIndex()
  
write_logfile()

if ("p" in ENSEMBLES["eq"]) or ("p" in ENSEMBLES["sim"]):
  BarostatChecker()

########################################################################
# 2. do the warmup loop                                                #
########################################################################

if ("cont" in MODUS) or (("auto" in MODUS) and os.path.isfile("restart.in.gz")):
  if CTIMELJ != "first":
    print "# this is a continued simulation run from timestep %.4f (%.4f ps) of previous %s run"%(CTIMELJ,CTIMESI,MODUS[MODUS.index(",")+1:])
  else:
    print "# this is a continued simulation run from a completed EQUILIBRATION run"
    CTIMELJ              = 0.0
    CTIMESI              = 0.0
  print SEPARATOR
if ("start" in MODUS) or (MODUS=="auto,UNSET"):
  integrator             = espresso.integrator.VelocityVerlet(system)
  integrator.dt          = DT_WARMUP
  
  thermostat             = espresso.integrator.LangevinThermostat(system)
  thermostat.gamma       = 10.0
  thermostat.temperature = float(TEMPERATURE)
  integrator.addExtension(thermostat)
  
  # setting up non-bonded interaction potential for warmup
  verletlist         = espresso.VerletList(system,WARM_CUTOFF)
  LJpotSS            = espresso.interaction.LennardJonesCapped(epsilon=EPSILON_START, sigma=SIGMASS, cutoff=WARMUP_CUTOFFSS, caprad=CAPRADIUSSS, shift='auto')
  LJpotSW            = espresso.interaction.LennardJonesCapped(epsilon=EPSILON_START, sigma=SIGMASW, cutoff=WARMUP_CUTOFFSW, caprad=CAPRADIUSSW, shift='auto')
  LJpotWW            = espresso.interaction.LennardJonesCapped(epsilon=EPSILON_START, sigma=SIGMAWW, cutoff=WARMUP_CUTOFFWW, caprad=CAPRADIUSWW, shift='auto')
  interaction        = espresso.interaction.VerletListLennardJonesCapped(verletlist)
  interaction.setPotential(type1=0, type2=0, potential=LJpotSS)
  interaction.setPotential(type1=0, type2=1, potential=LJpotSW)
  interaction.setPotential(type1=1, type2=1, potential=LJpotWW)

  # running the warmup loop
  system.addInteraction(interaction)

  start_time         = time.time()
  print "# starting warmup in nvt ensemble ..."

  for step in xrange(WARMUP_NLOOPS):
    
    if step%10 == 0:
      espresso.tools.info(system,integrator,per_atom=True)
    
    integrator.run(WARMUP_ISTEPS)
    
    LJpotSS.epsilon += EPSILON_DELTASS
    LJpotSW.epsilon += EPSILON_DELTASW
    LJpotWW.epsilon += EPSILON_DELTAWW
    interaction.setPotential(type1=0, type2=0, potential=LJpotSS)
    interaction.setPotential(type1=0, type2=1, potential=LJpotSW)
    interaction.setPotential(type1=1, type2=1, potential=LJpotWW)
  
  espresso.tools.info(system,integrator,per_atom=True)

  end_time     = time.time()
  duration     = end_time-start_time
  hours        = int(duration/3600.0)
  rest         = duration%3600
  minutes      = int(rest/60.0)
  seconds      = int(rest%60.0)
  print "# warmup duration: ", hours, "h ", minutes, "m ", seconds, "s"
  print "# warmup finished"
  print SEPARATOR
  
  system.removeInteraction(0)
  verletlist.disconnect()


########################################################################
# 3. setting up interaction potential for equilibration                #
########################################################################

if "EQUILIBRATION" in MODUS:
  integrator             = espresso.integrator.VelocityVerlet(system)
  integrator.dt          = DT_SIM
  
if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
  verletlist  = espresso.VerletList(system,SIM_CUTOFF)
  interaction = espresso.interaction.VerletListLennardJones(verletlist)
  potentialSS = interaction.setPotential(type1=0, type2=0,
                                       potential=espresso.interaction.LennardJones(
                                       epsilon=EPSILONSS, sigma=SIGMASS, cutoff=CUTOFFSS, shift=SHIFTYSS))
  potentialSW = interaction.setPotential(type1=0, type2=1,
                                       potential=espresso.interaction.LennardJones(
                                       epsilon=EPSILONSW, sigma=SIGMASW, cutoff=CUTOFFSW, shift=SHIFTYSW))
  potentialWW = interaction.setPotential(type1=1, type2=1,
                                       potential=espresso.interaction.LennardJones(
                                       epsilon=EPSILONWW, sigma=SIGMAWW, cutoff=CUTOFFWW, shift=SHIFTYWW))
				       
  # adding interactions
  system.addInteraction(interaction)

  # reset cell system
  system.storage.cellAdjust()


########################################################################
# 4. setting up ensemble for equilibration                             #
########################################################################

if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
  outline                  = "# "
  
  if ("start" in MODUS) or (MODUS=="auto,UNSET"):
    outline                 += "starting equilibration in %s ensemble"%ENSEMBLES["eq"]
    thermostat.gamma         = DAMPING_EQ

  if "EQUILIBRATION" in MODUS:
    outline                 += "continue equilibration in %s ensemble"%ENSEMBLES["eq"]
    thermostat               = espresso.integrator.LangevinThermostat(system)
    thermostat.gamma         = DAMPING_EQ
    thermostat.temperature   = float(TEMPERATURE)
    integrator.addExtension(thermostat)
    
  if "p" in ENSEMBLES["eq"]:
    barostat                 = espresso.integrator.BerendsenBarostat(system)
    barostat.tau             = PISTON_MASS_EQ
    barostat.pressure        = float(PRESSURE)
    barostat.fixed           = Int3D(1,0,0)
    integrator.addExtension(barostat)
  
  outline                 +=  " ..."
  print outline
  
  
########################################################################
# 5. do the equilibration loop                                         #
########################################################################

if ("start" in MODUS) or (MODUS=="auto,UNSET") or ("EQUILIBRATION" in MODUS):
  # initialization of measurements
  TEMP                = espresso.analysis.Temperature(system)
  PRESS               = espresso.analysis.PressureTensor(system)
 
  integrator.resetTimers()
  integrator.step     = 0 

  TRJFILE             = "config_series_eq.xyz.gz"
  DATAFILENAME        = "output_eq.dat"
  if ("start" in MODUS) or (MODUS=="auto,UNSET"):
    try:
      os.remove(TRJFILE)
    except OSError:
      pass
    try:
      os.remove(DATAFILENAME)
    except OSError:
      pass

  if "EQUILIBRATION" in MODUS:
    jumpToCurrentTime()

  start_time          = time.time()
  for STEP in xrange(EQUIL_NLOOPS):
    
    if STEP%CONFIG_EVERY_EQSTEP == 0:
      espresso.tools.writezippedxyz(TRJFILE,system,velocities = False,unfolded = False,append = True)

    if STEP%RESTART_EVERY_EQSTEP == 0:
      write_out_restart(STEP*EQUIL_ISTEPS*DT_SIM,"EQUILIBRATION")

    if STEP%DATA_EVERY_EQSTEP == 0:
      write_out_data("eq",EQUIL_ISTEPS)

    integrator.run(EQUIL_ISTEPS)
    
  espresso.tools.writezippedxyz(TRJFILE,system,velocities = False,unfolded = False,append = True)
  write_out_restart((STEP+1)*EQUIL_ISTEPS*DT_SIM,"EQUILIBRATION")
  
  end_time     = time.time()
  duration     = end_time-start_time
  hours        = int(duration/3600.0)
  rest         = duration%3600
  minutes      = int(rest/60.0)
  seconds      = int(rest%60.0)
  print "# equilibration duration: ", hours, "h ", minutes, "m ", seconds, "s"
  print "# equilibration finished"
  print SEPARATOR

  
########################################################################
# 6. setting up interaction potential for simulation                   #
########################################################################

if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
  reset_system()
  system.storage.decompose()
  print "# Particle numbers are set to:"
  print "#  total number of particles         = ", NPART_TOT
  print "#  particles in substrate            = ", NPART_XTAL
  print "#  particles in solution             = ", NPART_LIQ
  print "#  number of solute particles        = ", NPART_SOLUTE
  print "#  number of solvent particles       = ", NPART_SOLVENT
  print SEPARATOR
  
if "SIMULATION" in MODUS:
  integrator          = espresso.integrator.VelocityVerlet(system)
  integrator.dt       = DT_SIM
      
  verletlist          = espresso.VerletList(system,SIM_CUTOFF)
  interaction         = espresso.interaction.VerletListLennardJones(verletlist)
  potentialSS         = interaction.setPotential(type1=0, type2=0,
                                       potential=espresso.interaction.LennardJones(
                                       epsilon=EPSILONSS, sigma=SIGMASS, cutoff=CUTOFFSS, shift=SHIFTYSS))
  potentialSW         = interaction.setPotential(type1=0, type2=1,
                                       potential=espresso.interaction.LennardJones(
                                       epsilon=EPSILONSW, sigma=SIGMASW, cutoff=CUTOFFSW, shift=SHIFTYSW))
  potentialWW         = interaction.setPotential(type1=1, type2=1,
                                       potential=espresso.interaction.LennardJones(
                                       epsilon=EPSILONWW, sigma=SIGMAWW, cutoff=CUTOFFWW, shift=SHIFTYWW))
  # adding interactions
  system.addInteraction(interaction)

pair_bonds_xtal       = espresso.FixedPairList(system.storage)
pair_bonds_xtal.addBonds(XTALBONDS_LIST)
xtal_interaction      = espresso.interaction.FixedPairListHarmonic(system, pair_bonds_xtal, potential=espresso.interaction.Harmonic(K=K_XTAL,r0=NN_DIST_XTAL))
system.addInteraction(xtal_interaction)


########################################################################
# 7. running equilibration before simulation                           #
########################################################################

if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
  print "# starting simulation in %s ensemble with preceding equilibration in npt ensemble ..."%ENSEMBLES["sim"]

  thermostat.gamma           = DAMPING_SIM
    
  if "p" in ENSEMBLES["eq"]:
    barostat.tau             = PISTON_MASS_SIM
    
  if "v" in ENSEMBLES["eq"]:
    barostat                 = espresso.integrator.BerendsenBarostat(system)
    barostat.tau             = PISTON_MASS_SIM
    barostat.pressure        = float(PRESSURE)
    barostat.fixed           = Int3D(1,0,0)
    integrator.addExtension(barostat)
  
  integrator.resetTimers()
  integrator.step  = 0

  start_time       = time.time()
  for step in xrange(SIMWARM_NLOOPS):
    espresso.tools.info(system,integrator,per_atom=True)
    integrator.run(SIMWARM_ISTEPS)
  espresso.tools.info(system,integrator,per_atom=True)
  
  end_time     = time.time()
  duration     = end_time-start_time
  hours        = int(duration/3600.0)
  rest         = duration%3600
  minutes      = int(rest/60.0)
  seconds      = int(rest%60.0)
  print "# preceding equilibration duration: ", hours, "h ", minutes, "m ", seconds, "s"
  print "# preceding equilibration finished"
  print SEPARATOR


########################################################################
# 8. setting up ensemble for simulation                                #
########################################################################

if "SIMULATION" in MODUS:
  outline                    = "# continue simulation in %s ensemble ..."%ENSEMBLES["sim"]
  
  thermostat                 = espresso.integrator.LangevinThermostat(system)
  thermostat.temperature     = TEMPERATURE
  thermostat.gamma           = DAMPING_SIM
  integrator.addExtension(thermostat)
  
  if "p" in ENSEMBLES["sim"]:
    barostat                 = espresso.integrator.BerendsenBarostat(system)
    barostat.tau             = PISTON_MASS_SIM
    barostat.pressure        = float(PRESSURE)
    barostat.fixed           = Int3D(1,0,0)
    integrator.addExtension(barostat)
else:
  thermostat.gamma           = DAMPING_SIM
      
  if "v" in ENSEMBLES["sim"]:
    barostat.disconnect()


########################################################################
# 9. running the simulation loop                                       #
########################################################################

# initialization of measurements
TEMP                = espresso.analysis.Temperature(system)
TEMP.N              = NPART_XTAL
PRESS               = espresso.analysis.PressureTensor(system)
OLD_LX              = system.bc.boxL[0]

integrator.resetTimers()
integrator.step     = 0

BINS_INPLANE        = [int((system.bc.boxL[1]/BINSIZES[1])+0.5),int((system.bc.boxL[2]/BINSIZES[2])+0.5)]
CHECKFAC            = 0.5
TARGET_DENS         = [NUMDENS0*EXP_MOLF/MOLFRAC,0.0 if (MOLFRAC==1.0) else NUMDENS1*(1.0-EXP_MOLF)/(1.0-MOLFRAC)]
if CMUMD_MODE != "off":
  print "# CMuMD Setup:"
  print "#  number of bins within layer       =  (%d,%d)"%(BINS_INPLANE[0],BINS_INPLANE[1])
  print "#  size of single bin                =  (%.2f,%.2f,%.2f)"%(BINSIZES[0],BINSIZES[1],BINSIZES[2])
  print "#  BOP cutoff                        =  %.4f"%RDF_CUTOFF
  print "#  prefactor for checker method      =  %.4f"%CHECKFAC
  print SEPARATOR

if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
  WALL_INSTANCE     = Real3D(0.5*system.bc.boxL[0],0.5*system.bc.boxL[0]-max([CUTSS,CUTSW]),0.5*system.bc.boxL[0]+max([CUTSS,CUTSW]))
  print "# start the simulation ..."
if "SIMULATION" in MODUS:
  print outline

constMuMD              = espresso.integrator.ConstMuMD(system)
constMuMD.modus        = CMUMD_SWITCH[CMUMD_MODE]
constMuMD.sizeTR       = SIZETR
constMuMD.sizeCR       = SIZECR
constMuMD.sizeFR       = SIZEFR
constMuMD.sizeCheckR   = SIZECHECKR
#constMuMD.expDens      = TARGET_DENS
constMuMD.expMolF      = EXP_MOLF
constMuMD.dr           = BINSIZES
constMuMD.binsInPlane  = BINS_INPLANE
constMuMD.nThresh      = NPART_XTAL
constMuMD.cutoff       = RDF_CUTOFF
constMuMD.shapeParam   = INV_SHAPE_PARAM
constMuMD.fConst       = [FORCE_CONST_S,FORCE_CONST_W]
constMuMD.substrate    = WALL_INSTANCE
constMuMD.Ntot         = NPART_XTAL+NPART_SOLUTE
constMuMD.ratio        = CHECKFAC
integrator.addExtension(constMuMD)

 
TRJFILE             = "config_series_sim.xyz.gz"
DATAFILENAME        = "output_sim.dat"
if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
  CTIMELJ           = 0.0
  CTIMESI           = 0.0
  try:
    os.remove(TRJFILE)
  except OSError:
    pass
  try:
    os.remove(DATAFILENAME)
  except OSError:
    pass

if "SIMULATION" in MODUS:
  jumpToCurrentTime()
  
if ("start" in MODUS) or ("EQUILIBRATION" in MODUS) or (MODUS=="auto,UNSET"):
  write_out_restart("first","SIMULATION")

start_time          = time.time()
endreason           = "at end of simulation time"

integrator.run(0)
for STEP in xrange(SIM_NLOOPS):
  
  if (CMUMD_MODE != "off") and (STEP == 0):
    constMuMD.assignParts()

  if STEP%CONFIG_EVERY_SIMSTEP == 0:
    espresso.tools.writezippedxyz(TRJFILE,system,velocities = False,unfolded = False,append = True)

  if STEP%RESTART_EVERY_SIMSTEP == 0:
    write_out_restart(STEP*SIM_ISTEPS*DT_SIM,"SIMULATION")
    
  #if STEP%CMUMD_EVERY_SIMSTEP == 0:
    #constMuMD.adaptRegions()

  if STEP%DATA_EVERY_SIMSTEP == 0:
    write_out_data("sim",SIM_ISTEPS)

  if (constMuMD.sizeRes[0]<=MINSIZERES) or (constMuMD.sizeRes[1]<=MINSIZERES):
    endreason     = "as the reservoir became smaller than minSizeRes"
    endFile       = open("endofsim.txt",'w')
    endFile.close()
    break
 
  integrator.run(SIM_ISTEPS)
  
espresso.tools.writezippedxyz(TRJFILE,system,velocities = False,unfolded = False,append = True)
write_out_restart((STEP+1)*SIM_ISTEPS*DT_SIM,"SIMULATION")

end_time     = time.time()
duration     = end_time-start_time
hours        = int(duration/3600.0)
rest         = duration%3600
minutes      = int(rest/60.0)
seconds      = int(rest%60.0)
print "# simulation duration: ", hours, "h ", minutes, "m ", seconds, "s"
print "# simulation finished",endreason
print SEPARATOR
  

