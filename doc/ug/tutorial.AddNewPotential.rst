AddNewPotential
===============

.. |espp| replace:: ESPResSo++

The aim of the tutorial is to implement a new interaction potential in |espp|. We start with the Gromos fourth-power bond-stretching potential, because its functional form is simple and its implementation is somewhat similar to other potentials already implemented in |espp|. Everything you learn in this tutorial will then be relevent for implementing any other more complicated potential.

Make sure you have a working, compiled version of |espp| before starting the tutorial.

For those who are not so familiar with C++ or interfacing python and C++, you will find some helpful notes in the appendix.


Steps for adding a new interaction potential
--------------------------------------------

1. Choose the potential and derive the force.

2. Choose the appropriate interaction template from those in ``$ESPRESSOHOME/src/interaction``, e.g. ``VerletListInteractionTemplate.hpp``, ``FixedTripleListInteractionTemplate.hpp``

3. Create the .cpp, .hpp and .py files for your potential, place them in ``$ESPRESSOHOME/src/interaction`` and modify ``$ESPRESSOHOME/src/interaction/bindings.cpp`` and ``$ESPRESSOHOME/src/interaction/__init__.py``

4. Compile.

These steps are described in more detail below for our tutorial example potential.

Today's tutorial exercise
-------------------------

Step 1
......

The potential we are implementing today is a two-body bonded potential with the form 

.. math::
  V(r_{ij})=\frac{1}{4}k_{ij}(r^2_{ij}-r^2_{0,ij})^2

\noindent
where :math:`r_{ij}` is the distance between particles *i* and *j*. The potential has two input parameters :math:`r_{0}` and :math:`k`.

Derive the force.

Step 2
......

This is a 2-body interaction between a predefined (fixed) list of atom pairs. What is the appropriate interaction template to use? Choose one in ``$ESPRESSOHOME/src/interaction``

Open the interaction template file. (When you close the file later, close it without saving, or else later on your compile time will be very long, because of the number of dependencies on the interaction template!) Identify the functions ``addForces()`` and ``computeEnergy()``. Many interaction templates also contain functions such as ``computeVirial()``, ``computeVirialX()`` (for calculating the virial in slabs along the x-direction) etc.

Find the function calls:

.. code-block:: c

   potential->_computeForce(force, dist)

in ``addForces()`` and 

.. code-block:: c

   potential->_computeEnergy(r21)

in ``computeEnergy()``.

An interaction template can be combined with many different potentials (e.g. harmonic potential, Lennard Jones potential, etc.) Each potential will have its own C++ class containing functions to compute the energy and forces for that particular potential (see e.g. Harmonic.cpp/hpp, LennardJones.cpp/hpp) In turn, each potential can be combined with many different interaction templates.

You don't need to modify the interaction template file today. (Close it without saving!)

Step 3
......

In this step we create the .cpp, .hpp and .py files for our potential. Let's call the potential FourthPower. The FourthPower.py file will contain the end-user python interface, and in the FourthPower.cpp and FourthPower.hpp files we will create a C++ class for our potential. We will also write a wrapper which will allow the user to call the C++ code from the python interface.

3(a) Interfacing potential class and interaction template
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''

In many cases, it's not necessary to understand the contents of this section in order to implement a new potential. If you like, you can skip directly to Section `3(b) Creating the new potential class`_.

Now we need to understand how the interaction template will interface with our new class. This is done via a class template, e.g. in ``Potential.hpp``, ``AnglePotential.hpp``, ``DihedralPotential.hpp`` etc.

Still in ``$ESPRESSOHOME/src/interaction``, open the file ``Potential.hpp``. (When you close the file later, close it without saving, or else later on your compile time will be very long, because of the number of dependencies on the file!)

Find the functions ``_computeForce(Real3D& force, const Real3D& dist)`` and ``_computeEnergy(real dist)`` which you identified in the interaction template. Note that ``_computeForce(Real3D& force, const Real3D& dist)`` calls the function ``_computeForceRaw(force, dist, distSqr)`` and ``_computeEnergy(real dist)`` calls ``_computeEnergySqr(dist*dist)`` which calls ``_computeEnergySqrRaw(distSqr)``. The functions ``_computeForceRaw()`` and ``_computeEnergySqrRaw()`` are the new functions we need to write for our new potential. They will be member methods of our new C++ class FourthPower.

You don't need to modify anything in ``Potential.hpp`` today. (Close it without saving!)

3(b) Creating the new potential class
'''''''''''''''''''''''''''''''''''''

An easy way to implement the new C++ class is to identify a previously implemented potential which somewhat resembles your new potential, e.g. here we could take the Harmonic potential, which is also a 2-body potential, and which has also been interfaced with the FixedPairListInteractionTemplate.

Still in ``$ESPRESSOHOME/src/interaction``, copy the files ``Harmonic.py``, ``Harmonic.cpp`` and ``Harmonic.hpp`` to new files ``FourthPower.py``, ``FourthPower.cpp`` and ``FourthPower.hpp``. In the new files, find and replace all occurences of 'Harmonic' with 'FourthPower', and 'HARMONIC' with 'FOURTHPOWER'.

First modify ``FourthPower.hpp``.

Note the ``#include`` statement for ``FixedPairListInteractionTemplate.hpp`` and ``Potential.hpp``, the files you examined in `Step 2`_ and Step `3(a) Interfacing potential class and interaction template`_.

The Harmonic potential had parameters called ``K`` and ``r0``. You can reuse these for the FourthPower potential, along with the setters and getters ``setK``, ``getK``, ``setR0`` and ``getR0``. For better efficiency, you could also create a new variable which contains the square of ``r0``.

Now we need functions ``_computeForceRaw()`` and ``_computeEnergySqrRaw()``, as explained in Step `3(a) Interfacing potential class and interaction template`_. Modify these functions to use the functional form of the fourth power potential as derived in `Step 1`_. Note that ``Real3D dist``, which contains the vector between the two particles, has been defined as :math:`r_{p1} - r_{p2}` (see ``addForces()`` in ``FixedPairListInteractionTemplate.hpp``).

Next open ``Harmonic.py`` and ``FourthPower.py``.

Here is an example of an end-user's python script to add an interaction using the harmonic potential:

.. code-block:: python

   harmonicbondslist = espresso.FixedPairList(system.storage)
   harmonicbondslist.addBonds(bond_list) #bond_list is a list of tuples [(particleindex_i,particleindex_j),...]
   harmonic_potential = espresso.interaction.Harmonic(K=10.0, r0=1.0, cutoff = 5.0, shift = 0.0)
   harmonic_interaction = espresso.interaction.FixedPairListHarmonic(system, harmonicbondslist, potential=harmonic_potential)
   system.addInteraction(harmonic_interaction)

Compare this to the contents of ``Harmonic.py`` to understand the python source code.

Our new potential FourthPower can be called by the end-user in a similar way. Since the Harmonic and FourthPower potentials have similar input parameters (``K``, ``r0``) and both use the FixedPairListInteractionTemplate, you don't need to make any further modifications to the file ``FourthPower.py``, besides replacing 'Harmonic' with 'FourthPower'.

Next open ``FourthPower.cpp``.

Here you will find the C++/python interface, in the function ``registerPython()``. If you want to understand this function, you will find details in `Exposing a C++ class or struct to python using boost`_. You don't need to make any further modifications to this file, besides replacing 'Harmonic' with 'FourthPower'.

3(c) Including the new class in espressopp
''''''''''''''''''''''''''''''''''''''''''

Finally, update the files ``$ESPRESSOHOME/src/interaction/bindings.cpp`` and ``$ESPRESSOHOME/src/interaction/__init__.py`` (for example by copying and modifying all the lines referring to the Harmonic potential so that they now refer to the FourthPower potential). You need to make three modifications: to include the new .hpp file, to call the new registerPython() wrapper, and to import everything in the new python module.

Step 4
......

Move to the directory ``$ESPRESSOHOME``. Update the makefiles and compile using the commands:

.. code-block:: c

   cmake .
   make

.. Step 5
.. ......

.. Now test your code using the sample python script and input configuration ``FourthPowerSystem`` supplied with this pdf before the tutorial (or on a USB stick during the tutorial). Remember to run

.. .. code-block:: c

..   source $ESPRESSOHOME/ESPRC

.. in your working directory if necessary. Analyse the bond-length fluctuations using the script provided and compare them to the reference data provided.

Advanced exercise
-----------------

For an interaction potential of your choosing, follow the steps above to implement it, e.g. a non-bonded two-body interaction, probably using ``VerletListInteractionTemplate`` and based on the ``LennardJones`` potential, or a bonded three-body interaction, probably using ``FixedTripleListInteractionTemplate.hpp`` and based on the ``AngularHarmonic`` potential.

You will probably have to write setters and getters for the parameters in your potential in your .hpp file, and make the corresponding modifications to the function ``registerPython()`` in the .cpp file and the python user interface in the .py file.

Appendices
==========

Exposing a C++ class or struct to python using boost
----------------------------------------------------

(See http://www.boost.org/doc/libs/1_56_0/libs/python/doc/tutorial/doc/html/python/exposing.html)

Say we have a C++ struct called World:

.. code-block:: c

   struct World
   {
       World(std::string msg): msg(msg) {} 		// constructor
       void set(std::string msg) { this->msg = msg; }	// function set
       std::string greet() { return msg; }		// function greet
       std::string msg;					// member variable
   };

Now we write the C++ class wrapper for struct World to expose the constructor and the functions greet and set to python:

.. code-block:: c

   {
       class_<World>("World", init<std::string>())
           .def("greet", &World::greet)
           .def("set", &World::set)
       ;
   }

If there are additional constructors we can also expose them using ``def()``, e.g. for an additional constructor which takes two doubles:

.. code-block:: c

   class_<World>("World", init<std::string>())
       .def(init<double, double>())
       .def("greet", &World::greet)
       .def("set", &World::set)
   ;

We can also expose the data members of the C++ class or struct and the associated access (getter and setter) functions using ``add_property()``, e.g. for the variable myValue with access functions ``getMyValue`` and ``setMyValue``:

.. code-block:: c

    .add_property("myValue",&World::getMyValue,&World::setMyValue)

C++ classes and structs may be derived from other classes. Say we have the C++ struct myDerivedStruct which is derived from the struct myBaseStruct:

.. code-block:: c

   struct myBaseStruct { virtual ~myBaseStruct(); };
   struct myDerivedStruct : myBaseStruct {};

We can wrap the base class myBaseStruct as explained above:

.. code-block:: c

   <Base>("Base")
      /*...*/
      ;

Now when we want to wrap the class myDerivedStruct, we tell boost that it is derived from the base class myBaseStruct:

.. code-block:: c

   class_<myDerivedStruct, bases<myBaseStruct> >("myDerivedStruct")
       /*...*/
       ;

C++ templates
-------------

See http://www.cplusplus.com/doc/oldtutorial/templates/

typedef
-------

typedef declaration allows you to create an alias that can be used anywhere in place of a (possibly complex) type name

.. code-block:: c

   typedef DataType AliasName;

Python notes
------------

Syntax for classes in python
............................

(See also https://docs.python.org/2/tutorial/classes.html)

Here is a python class called DerivedClassName which is derived from two other base classes (BaseClassName1 and BaseClassName1), is initialised with two variables x and y which have default values 1 and 2, and contains a function myFunction.

.. code-block:: python

   class DerivedClassName(BaseClassName1, BaseClassName2):
       """docstring"""		#a way of providing some documentation for the class
       def __init__(self,x=1,y=2): #takes two variable which have default values 1 and 2
           self.x = x
           self.y = y
       def myFunction(self):
           return self.x * self.y

PMI
...

PMI = parallel method invocation. For more details see the file ``$ESPRESSOHOME/src/pmi.py``

.. `isController` is True when used on the controller (MPI root), False otherwise

.. `isWorker` = not isController






