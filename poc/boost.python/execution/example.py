# this is how to change the search path for module
# you can use this to include the Espresso binary module
import sys
sys.path.append('cpp')

print
print "Call this proof of concept with parameters, e.g."
print "      python example.py a b c"

# The embedded version will call a function in the C++-module that
# removes them before the script even starts. This demonstrates that
# initialization of e.g. MPI before python is possible using embedded
# python.
print "Here are the arguments as python got them: " + str(sys.argv)

# this is in the bin directory and sets up the paths to
# the python script modules
import espresso

# finally, get the POC
import hello

# the C++-module modifies argv again when imported - it adds a "Hi,
# World!". This demonstrates how to initialize MPI when the module is
# loaded. However, this only works if MPI accepts a copy of the real
# arguments. Nevertheless, MPI could also here manipulate the
# arguments.
print "Here is what the C++-module does to the arguments: " + str(sys.argv)

# and here, finally, a C++-member function called from python
h=hello.HelloWorld()

h.printMessage()
