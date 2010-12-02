"""Python functions to print timings from C++."""

import sys

def show(t, precision=1):
  fmt1 = '%.' + str(precision) + 'f\n'
  fmt2 = '%.' + str(precision) + 'f (%.'+ str(precision) + 'f)\n'
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
