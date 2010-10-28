Frequently Asked Questions
==========================

.. |espp| replace:: ESPResSo++

.. rubric:: Do I need to learn Python when using |espp|?

The short answer is "no". Most of the example scripts are self-explaining
and it is easy to adapt them for your purposes. In principle you can 
use |espp| also like other MD simulation software that is driven
by some kind of configuration file.

The long answer is "yes". If you want to take advantage of all the
nice features of |espp| you need a certain understanding of how the
Python interpreter works. 

But don't be afraid of learning Python:

 - Python is easy to learn
 - Looking at the |espp| simulation scripts gives you a very fast understanding of how Python works.
 - Writing programs in Python is much easier than writing programs in C++
 - Python programs are easier to understand than Tcl or Perl programs.

And here are some arguments why it is worth while:

 - A lot of Python programs are available that you can use in your applications
 - Python gives you a very flexible way of running MD simulations with |espp|

.. rubric:: Do you support other script languages, e.g. Tcl/Tk?

Certainly it might be possible to support also other script languages to access
the |espp| software, but we restricted ourselves to one script language to 
make it easier to exchange scripts between the user community.

.. rubric:: Can Tcl scripts converted to Python automatically?

The recommendation is - don't do it!  Instead, a Tcl interpreter can  be loaded via
Python and given the job to do. That is similar to what Tkinter does; Tkinter
is a wrapper to use the Tk toolkit from Python.

.. rubric:: Why should I use Python if C++ programs are much faster?

Our experience has also shown that Python programs are a factor of 30 up to 50 slower
than equivalent C++ programs.

But the use of Python is only intended to set up and steer the simulation while the 
simulation system itself is efficient C++ code.

.. rubric:: Can I run |espp| on parallel machines?

Yes. The parallel version uses MPI and is therefore as portable as MPI is. And
typical MD simulations scale rather well.

.. rubric:: Do I need to write parallel scripts for parallel machines?

No. The Python scripts are executed only by the first processor that will broadcast
the |espp| commands to all other processors automatically via a so called
PMI interface (Parallel Method Invocation). For you, it will look a serial script. But
the particles of the simulation are distributed among the available processors and the
issued commands for |espp| will be executed by each processor.

.. rubric:: How efficient is |espp|?

Efficiency got a very high priority but a less one than the extendability of the system. So you
should expect a good performance but might be that |espp| is less efficient than some other
simulation programs that are around. 

But if you see that |espp| is more than a factor of 2 slower than other simulation systems
you have detected a performance bug.

.. rubric:: Do I need the source code distribution or can I use binary version?

Currently, we only provide a source code distribution. But this might change in the future.
One major problem of providing binary versions is due to the fact that there are so many
different Python versions around and the binary versions of the libraries must not be
mixed.

.. rubric:: Which build systems are used for |espp|?

Compilation and installation of |espp| is due to the many shared libraries and loadable
modules rather complex and so we use a build system to make our job and maitainability 
easier.

Currently we support building the system with autotools.

.. rubric:: What means extendability?

Each software is in a certain sense extendable by just adding some functionality somewhere
in the code. But we understand extendability in the following sense:

 * You can add functionality without changing existent interfaces. For the object-oriented
   approach this means in practice: take an available base class and define a new derived 
   class with your needed functionality. 


