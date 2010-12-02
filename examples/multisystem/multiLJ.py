import espresso

comm=espresso.pmi.Communicator([1,3])
espresso.pmi.activate(comm)
system1=espresso.System(comm)
#system1.rng=espresso.esutil.RNG()
espresso.pmi.deactivate(comm)
