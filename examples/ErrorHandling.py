"""
This is an example of how error handling with exceptions works:
If you want to have reduced but better readable ESPResSo++ error messages
when your python script crashes, put a try: except: block around your script.
"""
# To see the error message start this example with more than 8 CPUs
import espresso

try:
  ncpus = espresso.MPI.COMM_WORLD.size
  box   = (5, 5, 5)
  rc    = 1.5
  skin  = 0.5
  nG    = espresso.tools.decomp.nodeGrid(ncpus)
  cG    = espresso.tools.decomp.cellGrid(box, nG, rc, skin)

  print "ncpus= %i box= %s rc+skin= %4.2f node_grid= %s cell_grid= %s" % (ncpus, box, rc+skin, nG, cG)

except espresso.Error as e:
  print e
