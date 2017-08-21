This file lists paths on /data/isilon/fogarty to simulations and analysis for all of the following works:

1. List of location of useful files (input tex and figure files for accepted articles, posters, presentation slides, useful graphics)

2. Source code for various analysis tools

3. Multiscale lysozyme ligand-binding free energy calculations

4. Using force-based adaptive resolution simulations to calculate solvation free energies of amino acid sidechain analogues, Elio Fiorentini, Kurt Kremer, Raffaello Potestio, Aoife C Fogarty, submitted to JCP, 2017

5. A Spatially Adaptive Multiple Time-Stepping Algorithm and its Application in Adaptive Resolution Simulations, w/ Karsten Kreis, 2017

6. A multi-resolution model to capture both global fluctuations of an enzyme and molecular recognition in the ligand-binding site, AC Fogarty, R Potestio, K Kremer, Proteins: Structure, Function, and Bioinformatics 84 (12), 1902-1913, 2016

7. Adaptive resolution simulations with self-adjusting high-resolution regions, K Kreis, R Potestio, K Kremer, AC Fogarty, J Chem Theory Comput, 12 (8), 4067-4081, 2016

8. Advantages and challenges in coupling an ideal gas to atomistic models in adaptive resolution simulations, Karsten Kreis, Aoife C Fogarty, Kurt Kremer, Raffaello Potestio, The European Physical Journal Special Topics 224 (12), 2289-2304, 2015

9. Adaptive resolution simulation of a biomolecule and its hydration shell: Structural and dynamical properties, AC Fogarty, R Potestio, K Kremer, J Chem Phys 142 (19), 195101, 2015 

##########################################################################

1. List of location of useful files (input tex and figure files for accepted articles, posters, presentation slides, useful graphics)

SUBPATH=/data/isilon/fogarty/write-ups 

* Files for posters presented: $SUBPATH/posters

* Slides for talks given: $SUBPATH/presentations

* Latex and figure files for accepted articles: $SUBPATH/accepted-articles

* Useful image files for posters and talks: $SUBPATH/useful-images-posters

* Conference abstrats: $SUBPATH/abstracts

* Other documents (JURECA proposal, JURECA status report, ERC open access presentation, review of literature about protein dynamics, coding tutorial about git, initial literature survey about lysozyme):
$SUBPATH/other-writeups

* Latex and figures for the pdf 'instructions for ENM and ENM-atomistic parametrisation':
$SUBPATH/lysozyme-project/parametrising-ENM

##########################################################################

2. Source code to perform various analysis

Paths are relative to SUBPATH=/data/isilon/fogarty/code

* NMR S2 order parameters from ENM simulations
  $SUBPATH/nmr-orderparameter

* NMR S2 order parameters from fully atomistic simulations
  $SUBPATH/analyse-proteins/nmr-orderparameter

* Molecular fluctuation profiles (spherical geometry)
  $SUBPATH/analyse-proteins/concentric-molfluct

* Density profiles (spherical geometry)
  $SUBPATH/analyse-proteins/concentric-sphr

* RDFs between a certain atom, and all atoms of a certain type
  $SUBPATH/analyse-proteins/rdf

* Distances between pairs of atoms:
  $SUBPATH/analyseStructure

* Electrostatic potential
  $SUBPATH/misc-scripts/elec_potential_cylinder.py

* Electric field
  $SUBPATH/analyse-proteins/electric-field

* Orientation time correlation function
  $SUBPATH/analyse-proteins/orcorr

* Bulk solvent density
  $SUBPATH/misc-scripts/box-density-reversed.py

* Tetrahedral order parameter
  $SUBPATH/analyse-proteins/tetr-orderparam

* C++/OpenGL code to produce 3D visualisation of multiscale protein

  git clone https://github.com/acfogarty/multiscale-Visualiser 
  or
  $SUBPATH/3d-multiscale-visualiser/multiscale-vis

* Thermodynamic force 
  $SUBPATH/thermodynamic-force

* For the following analysis I used gromacs tools: PCA, protein-protein H-bonds, RMSF, radius of gyration

##########################################################################

3. Multiscale lysozyme ligand-binding free energy calculations

Detailed write-up and figures:
https://www.overleaf.com/7783660jvbzhrskntvr

Paths are relative to SUBPATH=/data/isilon/fogarty/lysozyme/free-e

* Fully atomistic equilibration
  $SUBPATH/aa-md-substr/equil/2-nvt

* Re-parametrisation of water-protein excluded volume interaction 
  $SUBPATH/enm-aa/1-parametrise-sigma

* Fully atomistic free energy calculations, protein-ligand complex

  decoupling
  $SUBPATH/aa-md-substr/espp/5-TI-JURECA-decoupling 

  annihilation
  $SUBPATH/aa-md-substr/espp/3-TI-JURECA-annihilation

* Multiscale free energy calculations, protein-ligand complex

  complex, decoupling (where * is the number of atomistic protein residues)
  $SUBPATH/enm-aa/5-TI-JURECA-decoupling/aa-*

  complex, annihilation (8 atomistic residues)
  $SUBPATH/enm-aa/2-TI-JURECA-annihilation

* Fully atomistic free energy calculations, ligand in water

  annihilation
  $SUBPATH/ligand/1-TI-JURECA-annihilation

##########################################################################

4. Using force-based adaptive resolution simulations to calculate solvation free energies of amino acid sidechain analogues
w/ Elio Fiorentini, Kurt Kremer, Raffaello Potestio, Aoife C Fogarty
submitted to J Chem Phys, 2017

Article text and figures: 
https://www.overleaf.com/7688157hybtyfbcczfv

Files submitted to JCP:
/data/isilon/fogarty/write-ups/aacid-solvation/article/files-submitted-JCP

Zip containing article text and figures for the version including argument that systematic differences in free energies might be explained by the dependence of the free energy on the density (figure and text about that deleted before submission)
/data/isilon/fogarty/write-ups/aacid-solvation/article/article-including-argument-fn-density.zip

Detailed write-up and figures of all work not in article:
/data/isilon/fogarty/write-ups/aacid-solvation/writeup-labbook
https://www.overleaf.com/7688009gqrjfskmmjfm

Paths are relative to SUBPATH=/data/isilon/fogarty/aacid-solvation

* Obtaining thermodynamic force:
  $SUBPATH/tip3p/thd-force/thd-at-*p*-hy-1p2
  $SUBPATH/tip3p/thd-force-ideal-gas/thd-at-*p*-hy-1p2

* Free energy calculations (Fig 2):

  $SUBPATH/methanol/TI-espp/adres/ex-*p* (AdResS)
  $SUBPATH/methanol/TI-espp/at/box* (fully atomistic)
  $SUBPATH/methylindole/TI-espp/adres/ex-*p* (AdResS)
  $SUBPATH/methylindole/TI-espp/at/box* (fully atomistic) 

* Fully atomistic FE calculations at different densities (Fig 3):
  $SUBPATH/methanol/TI-gromacs/vary-density

* List of all free energy values
  $SUBPATH/plots/free-energies-methanol.dat
  $SUBPATH/plots/free-energies-methylindole.dat

* SVG file for Figure 1(a):
  $SUBPATH/plots/adres_illustration.svg

* Python scripts (other figures with xmgrace):

  to produce Figure 1(b):
  $SUBPATH/plots/plot-atomistic-systems.py

  to produce Figure 2:
  $SUBPATH/plots/free-energies-atIBIideal-3label.py
  $SUBPATH/plots/free-energies-atIBI-methylindole.py

  to produce Figure 3:
  $SUBPATH/plots/free-energies-at-density.py

* Calculation of molecular fluctuations (Fig 5):
  $SUBPATH/analysis/mol-fluctuations

* Calculation of density profiles:
  $SUBPATH/analysis/densities

6-ns production simulations, RDF analysis (Fig 4), and any other FE calculations were performed by Elio

##########################################################################

5. A Spatially Adaptive Multiple Time-Stepping Algorithm and its Application in Adaptive Resolution Simulations
w/ Karsten Kreis
2017
Paths are relative to SUBPATH=/data/isilon/fogarty/multi-timestepping

* Gold small system: $SUBPATH/gold-surface
  Gold big system: $SUBPATH/gold-surface-94

* Alanine dipeptide small system: $SUBPATH/ala-dipeptide
  Alanine dipeptide big system: $SUBPATH/ala-dipeptide-20-1

* Thermodynamic force for slab geometry: $SUBPATH/gold-surface/TF

##########################################################################

6. A multiresolution model to capture both global fluctuations of an enzyme and molecular recognition in the ligand-binding site
AC Fogarty, R Potestio, K Kremer
Proteins: Structure, Function, and Bioinformatics 84 (12), 1902-1913, 2016

(Note: substrate is another word for ligand)

All paths relative to SUBPATH=/data/isilon/fogarty/lysozyme

* Simulation trajectories are in the following locations:

  ** Fully atomistic, no substrate **
  
  $SUBPATH/aa-md/lysozyme-tutorial/amber99sb/from-substrate-config/3-nvt
  
  ** Fully atomistic, with substrate **
 
  $SUBPATH/aa-md-substr/8-nvt-prod 
  
  ** Multiscale, no substrate **

  (two independent trajectories from different starting configurations)
  $SUBPATH/enm-aa/ipcorr/r0-aamd/ca-beads/no-substrate/run-aa8/conf1
  $SUBPATH/enm-aa/ipcorr/r0-aamd/ca-beads/no-substrate/run-aa8/conf2
  
  ** Multiscale, with substrate **

  (two independent trajectories from different starting configurations)
  $SUBPATH/enm-aa/ipcorr/r0-aamd/ca-beads/NAG-substrate/run-aa8/conf1
  $SUBPATH/enm-aa/ipcorr/r0-aamd/ca-beads/NAG-substrate/run-aa8/conf2

* In each of the above locations are some of the following sub-directories, containing the analysis of each trajectory:

  analysis-struct: C-alpha-C-alpha distances across the active site (as in Figure 4)
  analysis-hbcount: counting number of H-bonds between water and certain protein active site atoms (as in Figure 5)
  analysis-struct-hbonding: ligand-protein distances in the active site (as in Figure 6)
  analysis-elecpot: electrostatic potential (as in Figure 7(b)
  analysis-elecfield: electric field (multiscale case) (as in Figure 7(a))
  analysis-other: water density in each region; electric field (atomistic case) (as in Figure 7(a))
  analysis-nmr: NMR S2 order parameters (as in Figure S1)
  analysis-rmsf-gromacs: root mean square fluctuations (as in Figure 3)

  Further analysis sub-directories not used in article:

  analysis-hinge: hinge angle
  analysis-contactmap: residue contact maps
  analysis-pca-gromacs: Principal Component Analysis

* ENM and multiscale parametrisation 
  
  ENM equilibrium distances from fully atomistic simulations:
  aa-md/lysozyme-tutorial/amber99sb/from-substrate-config/3-nvt/get-r0 (for a c-beta+c-alpha ENM)
  aa-md/lysozyme-tutorial/amber99sb/from-substrate-config/3-nvt/get-r0-calphaonly (convert C-beta to C-alpha ENM)

  ENM spring constant: (Figure S1) 
  $SUBPATH/enm/mod-ENM/test-enm-aamdr0/ipcorr-IBI/map-kn-knb-exclvol-allharmonic

  Reference data for parametrising ENM spring constant: (Figure S1)
  $SUBPATH/enm/mod-ENM/test-enm-aamdr0/ipcorr-IBI/map-kn-knb-exclvol-allharmonic/reference-data.dat

  protein-water excluded volume interaction:
  $SUBPATH/enm-aa/pcorr/r0-aamd/ca-cb-beads/no-substrate (everything below pcorr is zipped and tarred)

  Choosing which residues should be atomistic:
  $SUBPATH/aa-md-substr/8-nvt-prod/choose-neighbouring-residues

* Files for producing figures 3-7, S1, S3 including list of paths to raw data:
  $SUBPATH/article-figures

* C++/OpenGL code to produce Figure 1(a) and Figure S2 (3D multiscale protein image):

  git clone https://github.com/acfogarty/multiscale-Visualiser 
  or
  /data/isilon/fogarty/code/3d-multiscale-visualiser/multiscale-vis

* A day-by-day lab-book-style record of the project, probably not interesting for anyone except me: /data/isilon/fogarty/write-ups/other-writeups/lys-labbook

##########################################################################

7. Adaptive resolution simulations with self-adjusting high-resolution regions
K Kreis, R Potestio, K Kremer, AC Fogarty
Journal of Chemical Theory and Computation 12 (8), 4067-4081, 2016

All paths relative to SUBPATH=/data/isilon/fogarty/deformable-region

** Study of peptide folding (radius of gyration, ramachandran plots, Fig 14, 15, 16) **

Calculation of radius of gyration and helicity during folding (Fig 14, 15), ramachandran plots of folded peptide (Fig 15) and tetrahedral order parameter of folded peptide (Fig 16).

 * Atomistic peptide trajectories and analysis:
   $SUBPATH/polyala/atomistic/polyala9/amber99

   Sub-directories starting 1-, 2-, 3- are preparation, equilibration of the unfolded configurations

   Sub-directories called 4-*-nvt-prod are folding simulations using Gromacs (not included in article)

   Folding simulations used in article, and trajectory analysis, are in sub-directories called 4-*-espp-nvt-prod

 * AdResS peptide analysis:
   $SUBPATH/polyala/adres/karsten

   Analysis of production folding trajectories used in article is in sub-directories called CONF*-NORMALLANGEVIN. Karsten Kreis performed the MD simulations. Path to each of Karsten's trajectories can be found in each sub-directory CONF*-NORMALLANGEVIN.

** Static peptide **

  Atomistic and AdResS RDFs and orientational TCFs (Fig 13b,c). Analysis at $SUBPATH/polyala/analyse-karstens-staticprotein (MD simulations performed by Karsten Kreis)

** Bulk water, U-shaped region **

  Atomistic and AdResS tetrahedral order parameter (Fig 11b,c). Analysis at $SUBPATH/bulk-water/analyse-karstens-shape3 (MD simulations performed by Karsten Kreis)

** Additional analysis code **

  2D density grid: /data/isilon/fogarty/code/analyse-proteins/for-karsten-selfadjusting-region/grid-density
  Tetrahedral order parameter: /data/isilon/fogarty/code/analyse-proteins/tetr-orderparam
  Orientational time correlation function: /data/isilon/fogarty/code/analyse-proteins/orcorr
  
All other simulations and analysis by Karsten Kreis.

##########################################################################

8. Advantages and challenges in coupling an ideal gas to atomistic models in adaptive resolution simulations
Karsten Kreis, Aoife C Fogarty, Kurt Kremer, Raffaello Potestio
The European Physical Journal Special Topics 224 (12), 2289-2304, 2015

All paths relative to SUBPATH=/data/isilon/fogarty/spce/for-ideal-gas

All raw data files plotted in the paper are under:
$SUBPATH/spce/for-ideal-gas/data-files-for-karsten

* Force-AdResS simulations: (Fig 2)
  $SUBPATH/thd-force/1-espp-longer (without thermodynamic force)
  $SUBPATH/thd-force-first/10-espp-a-longer (with thermodynamic force)

* Diffusion profile calculations: (Fig 6, 7)
  $SUBPATH/diffusion-profiles/diffusion-profile-atomistic
  $SUBPATH/diffusion-profiles/diffusion-profile-hadres

* Molecular fluctuations calculations (Fig 8):
  $SUBPATH/molfluctuations/atomistic
  $SUBPATH/molfluctuations/hadres

Other figures (H-AdResS part) by Karsten Kreis

##########################################################################

9. Adaptive resolution simulation of a biomolecule and its hydration shell: Structural and dynamical properties
AC Fogarty, R Potestio, K Kremer
The Journal of chemical physics 142 (19), 195101, 2015

All paths relative to SUBPATH=/data/isilon/fogarty/ubiquitin-aa-hydrshell

One atomistic trajectory, four AdResS trajectories with different AT region sizes and one "protein in vacuum" trajectory are located at:

Fully atomistic $SUBPATH/aa-md-espp/long-runs/divergence-settlev
AdResS 3.0 $SUBPATH/adres/force-adres/long-runs/run-22-settlev-300-rf
AdResS 2.5 $SUBPATH/adres/force-adres/vary-hydrshell/ex-25
AdResS 2.0 $SUBPATH/adres/force-adres/vary-hydrshell/ex-2
AdResS 1.5 $SUBPATH/adres/force-adres/vary-hydrshell/ex-15
Vacuum $SUBPATH/aa-cg-md

Each directory contains the following sub-directories with the following analysis:

analysis-nmr: NMR order parameter (Fig 2)
analysis-gromacs: RMSF, HB count (Fig 2, Fig 7)
analysis-hydrshell: hydration shell orientational TCF (Fig 5, 8a,c)
analysis-hydrshell-aa: hydration shell (lambda=1) orientational TCF (Fig 8b,d)
analysis-other: tetrahedral order parameter, bulk orientational TCF, density profiles (Fig 3, 4, 5, 6, 9) 

System illustration (as in Fig 1) $SUBPATH/illustration

