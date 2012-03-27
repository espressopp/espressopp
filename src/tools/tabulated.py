import espresso
from espresso import Real3D

def writeTabFile(pot, name, N, low=0.0, high=2.5, body=2):
    """
    writeTabFile can be used to create a table for any potential
    Parameters are:
    * pot     : this is any espresso.interaction potential
    * name    : filename
    * N       : number of line to write
    * low     : lowest r (default is 0.0)
    * high    : highest r (default is 2.5)
    
    This function has not been tested for 3 and 4 body interactions
    """
    outfile = open(name, "w")
    delta = (high - low) / (N - 1)

    for i in range(N):
        r = low + i * delta
        energy = pot.computeEnergy(r)
        if body == 2:# this is for 2-body potentials
            force = pot.computeForce(Real3D(r, 0.0, 0.0))[0]
            #force /= r
        else: # this is for 3- and 4-body potentials
            force = pot.computeForce(r)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))

    outfile.close()
