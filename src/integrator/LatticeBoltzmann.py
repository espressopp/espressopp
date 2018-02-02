#  Copyright (C) 2012-2016,2017(H)
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008-2011
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
   
.. |espp| replace:: ESPResSo++

The LatticeBoltzmann (LB) class controls fluid hydrodynamics and allows for hybrid
LB/MD simulations. It is implemented as :py:class:`espressopp.integrator.Extension` in |espp|.

.. py:class:: espressopp.integrator.LatticeBoltzmann(system, nodeGrid, \
                                                    a = 1., tau = 1., \
                                                    numDims = 3, numVels = 19)

    :param shared_ptr system: system object
    :param Int3D nodeGrid: arrangement of CPUs in space
    :param real a: lattice spacing (in lattice units).
    :param real tau: time discretization (in lattice units)
    :param int numDims: dimensionality of the LB model
    :param int numVels: number of velocity vectors in the LB model
    :returns: lb object

    The LB-fluid in |espp| is aiming at simulations of complex soft matter systems. 
    They consist of MD particles (colloids, composite nanoparticles or polymer chains) 
    that are solved in the LB-fluid preserving hydrodynamic interactions.
    
    The default lattice model is D3Q19 (``numDims = 3``, ``numVels = 19``) and both lattice
    spacing ``a`` and timestep ``tau`` are set to 1. If some other lattice model is needed
    feel free to modify the code: adding 3D ones is straightforward, 
    for 2D cases one has to make more thouroughs changes.
    
    The parameters of the LB-fluid are expected in Lennard-Jones (LJ) units. This 
    strategy helps users with MD-background think of LB-fluid in term of LJ liquid. One
    only has to specify its properties such as liquid density,
    :math:`\\rho`, temperature, :math:`T`, and viscosity, :math:`\\eta`. 
    
    .. note::

        Standard LJ fluid can be characterized by
        :math:`\\rho \\sim 1 [\\sigma^{-3}]`,
        :math:`T \\sim 1 [\\epsilon]`, and 
        :math:`\\eta \\sim 5 [\\epsilon \\tau / \\sigma^3]`


    Example

    >>> L = 20
    >>>
    >>> # create cubic box
    >>> box = (L, L, L)
    >>> # The rc+skin= lattice_size
    >>> rc = 0.9
    >>>	skin= 0.1 
    >>>
    >>> # initialize empty default system with the created cubic box.
    >>> system, integrator = espressopp.standard_system.Default(box)
    >>>
    >>> # nodeGrid is determined based on the number of CPUs used for simulation among others
    >>> nodeGrid=espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc,skin)
    >>>
    >>> # initialize lb object. The dimensions of the lattice are obtained from the
    >>> # system's box dimensions employing lattice spacing 1.
    >>> lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)

    **Methods**

    .. py:method:: getLBMom(node, moment)
    
        Get hydrodynamic moment from a specific node
        
        :param Int3D node: node index
        :param int moment: hydrodynamic moment to get
        
        Use 0 to get density :math:`\\rho` and 1-3 for \
        mass flux components :math:`j_x`, :math:`j_y` and :math:`j_z`, correspondingly.
    
    .. py:method:: setLBMom(node, moment, value)
    
        Set hydrodynamic moment for a specific node
        
        :param Int3D node: node index
        :param int moment: hydrodynamic moment to set
        :param real value: value to set
    
    .. py:method:: saveLBConf()
    
        Dumps LB configuration with separate files for coupling forces, \
        LB-fluid moments and populations (the last one is a bit overkill). \
        The dump files are written for every CPU separately and \
        are put in the *dump* folder

    .. py:method:: keepLBDump()
    
        Sets a flag to keep previously dumped LB configuration. \
        Normally the previous dump is deleted after a new one is made.
        
        Example
        
        >>> # set bulk viscosity
        >>> for k in range (10):
        >>>     integrator.run(50000)
        >>>
        >>>     # output LB configuration
        >>>     lb.keepLBDump()         # flag to keep previously saved LB state
        >>>     lb.saveLBConf()         # saves current state of the LB fluid

    **Properties**

    .. py:data:: Int3D nodeGrid
    
        Array of CPUs in space
        
        Example
        
        >>> # it is advised to set nodeGrid by internal ESPResSo++ function
        >>> # based on the number of CPUs
        >>> nodeGrid=espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc,skin)


    .. py:data:: real a = 1.
    
        Lattice spacing (*lattice units*)
    
    .. py:data:: real tau = 1.
    
        Lattice time step (*lattice units*)
    
    .. py:data:: int numDims = 3
    
        Number of dimensions of the LB model (*D3Q19*)
    
    .. py:data:: int numVels = 19
    
        Number of velocity vectors of the LB model (*D3Q19*)

    .. py:data:: real visc_b
    
        Bulk viscosity (*LJ units*), affects :py:data:`gamma_b`.
        
        Example
        
        >>> # set bulk viscosity
        >>> lb.visc_b = 5.
        
    .. py:data:: real visc_s
    
        Shear viscosity (*LJ units*), affects :py:data:`gamma_s`.
        
        Example
        
        >>> # set shear viscosity
        >>> lb.visc_s = 5.
        
    .. py:data:: real gamma_b = 0.
    
        Bulk gamma (for experienced LB users)
        
    .. py:data:: real gamma_s = 0.

        Shear gamma (for experienced LB users)

    .. py:data:: real gamma_odd = 0.
    
        Odd gamma (for experienced LB users)
            
    .. py:data:: real gamma_even = 0.
    
        Even gamma (for experienced LB users)
            
    .. py:data:: real lbTemp = 0.
    
        Temperature of the LB fluid (*LJ units*)
        
        Example
        
        >>> L = 20
        >>> T = 1.
        >>> N = 200
        >>>
        >>> # create cubic box
        >>> box = (L, L, L)
	>>> rc=0.9
	>>> skin=0.1
        >>> 
        >>> # initialize Lennard Jones system with the created cubic box and given temperature.
        >>> system, integrator = espressopp.standard_system.LennardJones(N, box, temperature=T)
        >>>
        >>> # nodeGrid is determined based on the number of CPUs used for simulation
        >>> nodeGrid=espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc,skin)
        >>>
        >>> # initialize lb object. The dimensions of the lattice are obtained from the
        >>> # system's box dimensions employing lattice spacing 1.
        >>> lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
        >>>
        >>> # set LB temperature to T
        >>> lb.lbTemp = T

    .. py:data:: real fricCoeff=5.
    
        Friction coefficient of the coupling (*LJ units*)
        
        Example
        
        >>> # set friction coefficient of the coupling
        >>> lb.fricCoeff = 20.
    
    .. py:data:: int nSteps = 1
    
        Timescale contrast (ratio) between LB and MD
        
        Example
        
        >>> # set time step contrast between LB and MD
        >>> lb.nSteps = 10
        
    .. py:data:: int profStep = 10000
    
        Frequency of time profiling
    
        Example
    
        >>> # set profiling frequency
        >>> lb.profStep = 5000
    
    .. py:data:: Int3D getMyNi
            
        Number of real and halo nodes for the CPU
        
"""

from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp import Real3D
from espressopp import Int3D
from espressopp.integrator.Extension import *
from _espressopp import integrator_LatticeBoltzmann

class LatticeBoltzmannLocal(ExtensionLocal, integrator_LatticeBoltzmann):
    def __init__(self, system, nodeGrid,
                 a = 1., tau = 1., numDims = 3, numVels = 19):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LatticeBoltzmann,
                    system, nodeGrid, a, tau, numDims, numVels)

if pmi.isController :
    class LatticeBoltzmann(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
                            cls = 'espressopp.integrator.LatticeBoltzmannLocal',
                            pmiproperty = ['nodeGrid', 'a', 'tau', 'numDims', 'numVels', 'visc_b', 'visc_s', 'gamma_b', 'gamma_s', 'gamma_odd', 'gamma_even', 'lbTemp', 'fricCoeff', 'nSteps', 'profStep', 'getMyNi'],
                            pmicall = ["getLBMom","setLBMom","saveLBConf","keepLBDump"]
                            )
