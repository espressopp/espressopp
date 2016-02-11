#  Copyright (C) 2016
#      Max Planck Institute for Polymer Research & Johannes
#      Gutenberg-Universitaet Mainz
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

import espressopp
import difflib

system, integrator = espressopp.standard_system.LennardJones(0,(10,10,10))
system.storage.addParticles([[1,0,espressopp.Real3D(5,5,5)]],'id','type','pos')
system.storage.addParticles([[2,0,espressopp.Real3D(5,6,5)]],'id','type','pos')
system.storage.addParticles([[3,1,espressopp.Real3D(6,5,5)]],'id','type','pos')
system.storage.addParticles([[4,1,espressopp.Real3D(6,6,5)]],'id','type','pos')
system.storage.addParticles([[5,23,espressopp.Real3D(4,1,2)]],'id','type','pos')
system.storage.addParticles([[6,14,espressopp.Real3D(7,7,7)]],'id','type','pos')
system.storage.addParticles([[7,6,espressopp.Real3D(9,9,9)]],'id','type','pos')
system.storage.addParticles([[8,8,espressopp.Real3D(8,8,8)]],'id','type','pos')


jack = espressopp.io.DumpXYZ(system, integrator, filename='test_dumpXYZ_type_not_hardcoded.xyz', unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
jack.dump()

print "Finished writing file, run the script to check for diff!"

