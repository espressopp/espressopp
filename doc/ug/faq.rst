Frequently Asked Questions
==========================

.. |espp| replace:: ESPResSo++

.. rubric:: Do I need to learn Python when using |espp|?

The short answer is "no". Most of the example scripts are
self-explanatory and can be adapted for your purposes by simple
changes. You can also use |espp| like other MD simulation software,
that is driven by some kind of configuration file.

The long answer is "yes". If you want to take advantage of all
features of |espp| you need some knowledge of how the
Python interpreter works. 

But don't be afraid of learning Python:

 - Python is easy to learn
 - The |espp| example simulation scripts gives you a very fast
   insight of how Python works.
 - Writing programs in Python is much easier than writing programs in C++
 - Python programs are easier to read than Tcl or Perl programs.

And here are some arguments why it is worth while:

 - There are many Python programs you can use in your applications
 - Python gives you a flexible way of running MD simulations with |espp|

.. rubric:: Do you support other script languages, e.g. Tcl/Tk?

No.  We choose the support only Python as |espp| scripting language.
This enables |espp| users to read and adapt scripts written by other
|espp| users.


.. rubric:: Can Tcl scripts converted to Python automatically?

The recommendation is - don't do it!  Instead, a Tcl interpreter can
be loaded via Python and given the job to do. That is similar to what
Tkinter does; Tkinter is a wrapper to use the Tk toolkit from Python.

.. rubric:: Why should I use Python if C++ programs are much faster?

Python is the driver of your simulation which will still run in the
|espp| C++ engine.

Python programs are about 30 to 50 times slower than the same programs
written in C++.

That is why we use Python to set up and control simulations while the 
simulation system itself is written in efficient C++ code.

.. rubric:: Can I run |espp| on parallel machines?

Yes. The parallel version uses MPI and is therefore as portable as MPI is. 
Typical MD simulations scale rather well.

.. rubric:: Do I need to write parallel scripts for parallel machines?

No. The Python scripts are executed only by the first processor which
will broadcast the |espp| commands to the other processors
automatically using the PMI interface (Parallel Method
Invocation). For you, it will look a serial script. But the particles
of the simulation are distributed among the available processors and
the commands issued for |espp| will be executed by each processor.

.. rubric:: How efficient is |espp|?

Efficiency is a high priority though less than the extendability of the system. You
should expect a good performance but might be that |espp| is less efficient than other
simulation programs that are around. 

If you experience that |espp| is more than 2 times slower than other simulation systems
you have found a performance bug.

.. rubric:: Do I need the source code distribution or can I use binary version?

Currently, we only provide a source code distribution. This might
change in the future.  Our major problem to provide binary versions is
that there are many different Python versions and the binary versions
of the |espp| and python libraries must not be of mixed versions.

.. rubric:: Which build systems are used for |espp|?

Compilation and installation of |espp| is due to the many shared
libraries and loadable modules rather complex and so we use a build
system to make our job and maintainability easier.

Currently we support building the system with cmake.

.. rubric:: What means extendability?

Each software is in a certain sense extendable by adding some
functionality somewhere in the code. But we understand extendability
in the following sense:

 * You can add functionality without changing existent interfaces. For the object-oriented
   approach this means in practice: take an available base class and define a new derived 
   class with your needed functionality. 


