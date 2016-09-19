Thermodynamic integration
=========================

.. |espp| replace:: ESPResSo++

Theoretical explanation
-----------------------

Thermodynamic integration (TI) is a method used to calculate the free energy difference between two states A and B. For the theoretical background, see e.g. http://www.alchemistry.org. In this tutorial, we show how to perform TI calculations with ESPResSo++. We calculate the free energy of solvation of methanol in water. The complete python script is available in the ESPResSo+ source code under ``examples/thd_integration_solvation`` 

To do TI, we define states A and B, with potentials :math:`U^A` and :math:`U^B`. We then construct a pathway of intermediate states between A and B by defining a parameter :math:`\lambda` that takes values between 0 and 1 and writing the system potential :math:`U` as a function of :math:`\lambda`, :math:`U^A` and :math:`U^B`. The free energy difference between the states A and B is then given by

.. math::
  \Delta A = \int^1_0 \left<\frac{dU(\lambda)}{d\lambda}\right>_{\lambda}\mathrm{d}\lambda

In practise, we discretise :math:`\lambda` and perform a series of MD simulations with different :math:`\lambda` values between 0 and 1, sampling :math:`\frac{dU(\lambda)}{d\lambda}` in each simulation.

To calculate the solvation free energy of methanol in water, we use a box of water containing one methanol molecule. We simulate desolvation via two separate TI calculations. (Note that the procedure described here is decoupling, and solute-solute interactions will be treated differently if you're doing annihilation instead of decoupling, see Note 1.) 

----

**Step 1:** free energy change for switching off the Coulombic interactions

*State A:* methanol has full non-bonded (Coulomb and Lennard Jones) interactions with the solvent

*State B:* methanol has only Lennard Jones interactions with the solvent

----

**Step 2:** free energy change for switching off the Lennard Jones interactions

*State A:* methanol has only Lennard Jones interactions with the solvent

*State B:* methanol has no interaction with the solvent

----

Step 1 can be done using a linear function of :math:`\lambda`:

.. math::
  U(\lambda_C) = (1-\lambda_C)U_C^A + U_{unaffected}

where :math:`U_C^A` is the solute-solvent Coulombic interaction in state A. In |espp| the charges used for state A are the particle charges contained in the particle property ``charge``. The charges in state B are zero, so :math:`U_C^B(q)` does not appear in the expression. (The case where A and B both have non-zero charges is not implemented in ESPResSo++). The term :math:`U_{unaffected}` is all other parts of the potential that don't change with :math:`\lambda_C` including all bonded interactions, any solute-solute Coulombic interactions, solvent-solvent Coulombic interactions and all Lennard-Jones interactions. The parameter :math:`\lambda_C` goes from 0 to 1 in Step 1.

Step 2 must be done using a softcore potential because of the singularity in the Lennard-Jones potential at :math:`r_{ij} = 0`.

.. math::
  U(\lambda_L) = \sum_{i,j} U_L(r_{ij},\lambda_L) + U_{unaffected}

  U_L(r_{ij},\lambda_L) = (1-\lambda_L)U_H^A(r_A) + \lambda_L U_H^B(r_B)

  r_A=(\alpha\sigma^6_A\lambda^p+r_{ij}^6)^{1/6}

  r_B=(\alpha\sigma^6_B(1-\lambda)^p+r_{ij}^6)^{1/6}

The terms :math:`U_H^A(r_A)` and :math:`U_H^B(r_B)` are the normal Lennard-Jones 12-6 hardcore potentials:

.. math::
  U_H^A(r_A) = 4.0\epsilon_A(\frac{\sigma_A}{r_A}^{12} - \frac{\sigma_A}{r_A}^6)

The sum :math:`\sum_{i,j} U_L(r_{ij},\lambda_L)` is over all solute-solvent interactions. The term :math:`U_{unaffected}` is all other parts of the potential that don't change with :math:`\lambda_L` including any solute-solute Lennard-Jones interactions and solvent-solvent Lennard-Jones interactions, which are treated using standard hardcore Lennard-Jones. (In this particular example of methanol, there are no solute-solute Lennard-Jones interactions). Finally :math:`\alpha` and :math:`p` are adjustable parameters of the softcore potential.

The |espp| C++ code allows for different values of :math:`\epsilon_A`, :math:`\epsilon_B`, :math:`\sigma_A` and :math:`\sigma_B` for every pair of atomtypes interacting via this potential. In this example, we will set :math:`\epsilon_B` to 0 (we are switching off the Lennard-Jones interaction). The parameter :math:`\lambda_L` goes from 0 to 1 in Step 2.


|espp| code
-----------

We must perform many separate simulations, each with a different :math:`\lambda` value. It is convenient to define a list of :math:`\lambda` values in the python script and use an index to access a different element of the list in each separate simulation. The script for the first simulation contains these lines:

.. code-block:: python

  # Parameters for Thermodynamic Integration
  stateBIndices = [1,2,3,4,5,6] #indices of the methanol atoms
  lambdaVectorCoul = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 
                      0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.000, 
                      1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                      1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 
                      1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 
                      1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 
                      1.000, 1.000, 1.000]
  lambdaVectorVdwl = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.025, 
                      0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 
                      0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 
                      0.500, 0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 
                      0.725, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 
                      0.950, 0.975, 1.000]
  lambdaIndex = 0
  lambdaTICoul = lambdaVectorCoul[lambdaIndex]
  lambdaTIVdwl = lambdaVectorVdwl[lambdaIndex]

The list ``lambdaVectorCoul`` contains the values of :math:`\lambda_C` and the list ``lambdaVectorVdwl`` contains the values of :math:`\lambda_L`. The total number of simulations to do Step 1 and Step 2 will be ``len(lambdaVectorCoul)`` or ``len(lambdaVectorVdwl)``. We must make a copy of the python script for each simulation, changing each time the value of ``lambdaIndex``.

Next we set up the Coulombic interactions, assuming we already have created a ``system`` and a ``verletlist``. The electrostatics method used is generalised reaction field.

.. code-block:: python

  #atTypes - list of all atomtypes (integers) used in the pairs interacting via this potential 
  #epsilon1,epsilon2,kappa - reaction field parameters
  #annihilate=False means decoupling is used (see Note 1)
  #ftpl - a FixedTupleListAdResS object (see AdResS tutorial)
  #for non-AdResS simulations, simply set adress=False, and the parameter ftpl is not needed
  qq_adres_interaction = gromacs.setCoulombInteractionsTI(system, verletlist, nbCutoff, 
                                                  atTypes, epsilon1=1, epsilon2=80, 
                                                  kappa=0, lambdaTI=lambdaTICoul, 
                                                  pidlist=stateBIndices, 
                                                  annihilate=False, adress=True, ftpl=ftpl)

Now we set up the softcore Lennard Jones interaction.

.. code-block:: python

  #atomtypeparameters - dictionary of format {atomtype: {'eps': epsilon, 'sig': sigma}} 
  #                     where atomtype is integer and epsilon and sigma are real
  #defaults - dictionary containing a key 'combinationrule' with value 1 if the contents 
  #           of atomtypeparameters need to be converted from c6,c12 format to 
  #           epsilon,sigma format; can also be an empty dictionary if no conversion needed
  #sigmaSC, alphaSC, powerSC - parameters of the softcore potential
  alphaSC = 0.5
  powerSC = 1.0
  epsilonB = 0.0
  sigmaSC = 0.3
  lj_adres_interaction = gromacs.setLennardJonesInteractionsTI(system, defaults, 
                                         atomtypeparameters, verletlist, nbCutoff, 
                                         epsilonB=epsilonB, sigmaSC=sigmaSC, alphaSC=alphaSC, 
                                         powerSC=powerSC, lambdaTI=lambdaTIVdwl, 
                                         pidlist=stateBIndices, annihilate=False, 
                                         adress=True, ftpl=ftpl)

We open an output file. In the first line we write the values of :math:`\lambda_C` and :math:`\lambda_L` for this simulation. 

.. code-block:: python

  dhdlF = open("dhdl.xvg","a")
  dhdlF.write("#(coul-lambda, vdw-lambda) = ("+str(lambdaTICoul)+", "+str(lambdaTIVdwl)+")\n")

During the MD run, every x number of MD steps, we return to the python level and calculate the derivatives of the energies with respect to :math:`\lambda`.

.. code-block:: python

  dhdlCoul = qq_adres_interaction.computeEnergyDeriv()
  dhdlVdwl = lj_adres_interaction.computeEnergyDeriv()
  dhdlF.write(str(time)+" "+str(dhdlCoul)+" "+str(dhdlVdwl)+"\n")

After all simulations, we can now average :math:`\frac{dU(\lambda)}{d\lambda}` for each value of :math:`\lambda_C` or :math:`\lambda_L`, integrate over :math:`\lambda_C` and :math:`\lambda_L`, add the values :math:`\Delta A_C` and :math:`\Delta A_L`, and take the negative (because the procedure described here is desolvation and we want the free energy of solvation).

Some notes
----------

1. This example given here uses decoupling (solute-solvent interactions are a function of :math:`\lambda`, solute-solute interactions are not affected by changes in :math:`\lambda`). In |espp| it is also possible to do annihilation, where both solute-solvent and solute-solute interactions are a function of :math:`\lambda`, by setting ``annihilate=True`` when creating the non-bonded interactions.

2. The procedure described here is desolvation. To get the free energy of solvation, we take the negative of the value obtained after integration.

3. The example Python code snippets here use the helper functions ``gromacs.setLennardJonesInteractionsTI`` and ``gromacs.setCoulombInteractionsTI`` contained in ``$ESPRESSOHOME/src/tools/convert/gromacs.py``, but this is not necessary. You can do TI with |espp| without the Gromacs parser by directly calling ``espresso.interaction.LennardJonesSoftcoreTI`` and ``espresso.interaction.ReactionFieldGeneralizedTI``. See the documentation of these two classes. 

.. 2. For convenience, there is a module ``WCASoftcoreTI`` for doing TI in WCA systems, though in principal one could also use ``LJSoftcoreTI`` to simulate the WCA potential by choosing the appropriate cutoff and post-processing the output.
