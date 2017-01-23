ESPResSo++
==========

ESPResSo++ is an extensible, flexible, fast and parallel simulation software for
soft matter research. It is a highly versatile software package for the
scientific simulation and analysis of coarse-grained atomistic or bead-spring
models as they are used in soft matter research. ESPResSo and ESPResSo++ have
common roots and share parts of the developer/user community. However their
development is independent and they are different software packages. ESPResSo++
is free, open-source software published under the GNU General Public License
(GPL).

QUICKSTART:
===========

To get a copy of the developer version (most recent version) of ESPResSo++, you can use git or docker. Using docker will give you a binary release (nothing to compile, but performance may not be optimal). If you use git clone or download a tarball, you will have to compile ESPResSo+ yourself, which might lead to better performance.

Using [docker](https://www.docker.com):
```
# docker pull espressopp/espressopp
# docker run -it espressopp/espressopp /bin/bash
```

Using git:
```
# git clone https://github.com/espressopp/espressopp.git
```

Alternatively, you can download a tarball or zip file of [previous release versions](https://github.com/espressopp/espressopp/releases) of ESPResSo++.

QUICKINSTALL:
=============

```
# cd espressopp
# cmake -DEXTERNAL_BOOST=OFF -DEXTERNAL_MPI4PY=OFF .
# make -j2
# source ESPRC
```

After building go to the `examples` directory and have a look at the python scripts.

DOCUMENTATION:
==============

Documentation for the developer version of ESPResSo++ is at:

http://espressopp.github.io

Documentation for release versions from v1.9.4.1 onward is at:

http://espressopp.github.io/vXXX

where XXX is the version number, e.g.: 

http://espressopp.github.io/v1.9.4.1

ISSUES:
=======

Report bugs on the [github issues site](https://github.com/espressopp/espressopp/issues)


[![Build Status](https://travis-ci.org/espressopp/espressopp.svg?branch=master)](https://travis-ci.org/espressopp/espressopp)
[![Code Climate](https://codeclimate.com/github/espressopp/espressopp/badges/gpa.svg)](https://codeclimate.com/github/espressopp/espressopp)
[![codecov.io](https://codecov.io/github/espressopp/espressopp/coverage.svg?branch=master)](https://codecov.io/github/espressopp/espressopp?branch=master)
[![Analysis Status](https://scan.coverity.com/projects/8143/badge.svg?flat=1)](https://scan.coverity.com/projects/espressopp-espressopp)

