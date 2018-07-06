#  Copyright (C) 2017
#      Jakub Krajniak (jkrajniak at gmail.com)
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


import random

def gaussian(T, N, particle_mass=None, zero_momentum=True, seed=7654321, kb=1.0):
    """Generates velocities with temperature T according to a Maxwell-Boltzmann distribution.

    Args:
        T: The desired temperature expre.
        N: The number of particles.
        particle_mass: The list of particle mass if not then every particle has mass 1.0
        zero_momentum: Remove the center-of-mass motion.
        seed: The seed for the random number generator.
        kb: The Boltzmann constant.

    Returns:
        The tuple with lists of x,y,z components of the velocity.
    """

    random.seed(seed)

    boltzmann_factor = kb*T
    E_kin = 0.0

    vx = []
    vy = []
    vz = []

    for i in range(N):
        if particle_mass:
            mass = particle_mass[i]
        else:
            mass = 1.0
        sd = (boltzmann_factor/mass)**0.5

        v_x = sd*random.gauss(0.0, 1.0)
        E_kin += 0.5*mass*v_x*v_x
        vx.append(v_x)

        v_y = sd*random.gauss(0.0, 1.0)
        E_kin += 0.5*mass*v_y*v_y
        vy.append(v_y)

        v_z = sd*random.gauss(0.0, 1.0)
        E_kin += 0.5*mass*v_z*v_z
        vz.append(v_z)

    # Scale the temperature
    local_temp = (2.0*E_kin)/(3*N*kb)
    if local_temp > 0:
        temp_scale = (T/local_temp)**0.5
        for i in range(N):
            vx[i] *= temp_scale
            vy[i] *= temp_scale
            vz[i] *= temp_scale

    if zero_momentum:
        # remove net momentum
        sumvx = 0.0
        sumvy = 0.0
        sumvz = 0.0
        for vx_, vy_, vz_ in zip(vx, vy, vz):
            sumvx += vx_
            sumvy += vy_
            sumvz += vz_
        sumvx = sumvx / N
        sumvy = sumvy / N
        sumvz = sumvz / N
        for i in range(N):
            vx[i] = vx[i] - sumvx
            vy[i] = vy[i] - sumvy
            vz[i] = vz[i] - sumvz
    return vx, vy, vz
