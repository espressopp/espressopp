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


import sys

def show(alltimers, precision=1):
  """
  Python functions to print timings from C++.
  """
  
  fmt1 = '%.' + str(precision) + 'f\n'
  fmt2 = '%.' + str(precision) + 'f (%.'+ str(precision) + 'f)\n'
  t=[]
  nprocs = len(alltimers)
  for ntimer in xrange(10):
    t.append(0.0)
    for k in xrange(nprocs):
      t[ntimer] += alltimers[k][ntimer]
    t[ntimer] /= nprocs          
          
  sys.stdout.write('Run    time (%) = ' + fmt1 % t[0])
  sys.stdout.write('Pair   time (%) = ' + fmt2 % (t[1], 100*t[1]/t[0]))
  sys.stdout.write('FENE   time (%) = ' + fmt2 % (t[2], 100*t[2]/t[0]))
  sys.stdout.write('Angle  time (%) = ' + fmt2 % (t[3], 100*t[3]/t[0]))
  sys.stdout.write('Comm1  time (%) = ' + fmt2 % (t[4], 100*t[4]/t[0]))
  sys.stdout.write('Comm2  time (%) = ' + fmt2 % (t[5], 100*t[5]/t[0]))
  sys.stdout.write('Int1   time (%) = ' + fmt2 % (t[6], 100*t[6]/t[0]))
  sys.stdout.write('Int2   time (%) = ' + fmt2 % (t[7], 100*t[7]/t[0]))
  sys.stdout.write('Resort time (%) = ' + fmt2 % (t[8], 100*t[8]/t[0]))
  sys.stdout.write('Other  time (%) = ' + fmt2 % (t[9], 100*t[9]/t[0]))
  sys.stdout.write('\n')
