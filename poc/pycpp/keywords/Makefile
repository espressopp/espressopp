#######################################################################
# This Makefile creates a binary file "pycpp"                         #
# which can be executed directly and a library file "pycpp.so"        #
# which can be imported in Python as an extension with "import pycpp" #
#######################################################################

# where to find the include files for boost and python
# CPPFLAGS=\
# 	-I /people/thnfs/homes/lenzo/software/include/python2.5 \
# 	-I /people/thnfs/homes/lenzo/software/include/boost-1_35

# # where to find the boost and python libraries, plus their names
# LDFLAGS=-L /people/thnfs/homes/lenzo/software/lib \
# 	-lboost_python-gcc42-mt -lpython2.5

# qwghlm:
CPPFLAGS=\
	-I /home/olenz/software/include/boost-1_35 \
	-I /usr/include/python2.5
LDFLAGS=\
	-L /home/olenz/software/lib/ \
	-lboost_python-gcc42-mt \
	-lpython2.5

# there is no need to change anything below this line

CXX=g++
CC=gcc
CXXFLAGS=-fPIC -O3
CFLAGS=-fPIC -O3

all: pycpp.so pycpp
clean:
	-rm -f pycpp.so
	-rm -f pycpp.o
	-rm -f pycpp

pycpp.so: pycpp.o

%.so: %.o
	$(CXX) -shared $(CXXFLAGS) $(LDFLAGS) -o $@ $<
