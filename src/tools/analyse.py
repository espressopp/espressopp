import sys
import espresso

def info(system, integrator):
  NPart  = espresso.analysis.NPart(system).compute()
  T      = espresso.analysis.Temperature(system).compute()
  P      = espresso.analysis.Pressure(system).compute()
  Pij    = espresso.analysis.PressureTensor(system).compute()
  step   = integrator.step
  Ek     = (3.0/2.0) * T
  Epot   = []
  Etotal = 0.0
  tot    = '%5d %8.6f %10.6f %10.6f %12.8f' % (step, T, P, Pij[3], Ek)
  tt     = ''
  for k in range(system.getNumberOfInteractions()):
    e       = system.getInteraction(k).computeEnergy()/NPart
    Etotal += e
    tot    += ' %12.8f' % e
    tt     += '      e%i     ' % k
  tot += ' %12.8f' % Etotal
  tt  += '    etotal   '
  tot += ' %12.8f\n' % system.bc.boxL[0]
  tt  += '    boxL     \n'
  if step == 0:
    sys.stdout.write(' step     T          P        Pxy        ekinetic ' + tt)
  sys.stdout.write(tot)

def final_info(system, integrator, vl, start_time, end_time):
  NPart  = espresso.analysis.NPart(system).compute()
  espresso.tools.timers.show(integrator.getTimers(), precision=3)
  sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
  sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(NPart)))
  sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
  sys.stdout.write('Integration steps = %d\n' % integrator.step)
  sys.stdout.write('CPUs = %i CPU time per CPU = %.5f\n' % (espresso.MPI.COMM_WORLD.size, end_time - start_time))

