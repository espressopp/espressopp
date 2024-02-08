#  Copyright (C) 2020-2022
#      Max Planck Institute for Polymer Research & JGU Mainz
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

from espressopp.esutil import pmiimport
pmiimport('espressopp.hpx4espp')

from _espressopp import hpx4espp_enabled, hpx4espp_have_networking

def enabled():
    return hpx4espp_enabled()

def have_networking():
    return hpx4espp_have_networking()

if enabled():
    from espressopp.hpx4espp.HPXRuntime import *
    from espressopp.hpx4espp.SystemHPX import *
    from espressopp.hpx4espp.VerletList import *
    from espressopp.hpx4espp import esutil, integrator, interaction, storage, standard_system
