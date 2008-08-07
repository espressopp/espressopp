# To be able to run the examples, put 
# 

# MPIP:
# where to find the include files for boost and python
#CPPFLAGS=\
#	-I /people/thnfs/homes/lenzo/software/include/python2.5 \
#	-I /people/thnfs/homes/lenzo/software/include/boost-1_35
#
## where to find the boost and python libraries, plus their names
#LDFLAGS=-L /people/thnfs/homes/lenzo/software/lib \
#	-lboost_python-gcc42-mt -lpython2.5

# qwghlm:
CPPFLAGS=\
	-I /home/olenz/software/include/boost-1_35 \
	-I /usr/include/python2.5

LDFLAGS=\
	-L /home/olenz/software/lib \
	-lboost_python-gcc42-mt-1_35 \
	-lpython2.5

#-Wl,-rpath /home/olenz/software/lib/ \

# there is no need to change anything below this line
CXX=g++
CC=gcc
CXXFLAGS=-fPIC -g -Wall
CFLAGS=-fPIC -g -Wall

all: world.so
clean:
	-rm -f *.o
	-rm -f *.so

%.so: %.o
	$(CXX) -shared $(CXXFLAGS) $(LDFLAGS) -o $@ $^
