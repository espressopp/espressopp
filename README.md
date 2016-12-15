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

Using [docker](https://www.docker.com):
```
# docker pull espressopp/espressopp
# docker run -it espressopp/espressopp /bin/bash
```


QUICKINSTALL:
=============

```
# cd espressopp
# cmake -DEXTERNAL_BOOST=OFF -DEXTERNAL_MPI4PY=OFF .
# make -j -l$(nproc)-2
# source ESPRC
```

After building go to the `examples` directory and have a look at the python scripts.

ISSUES:
=======

Report bugs on the [github issues site](https://github.com/espressopp/espressopp/issues)


[![Build Status](https://travis-ci.org/espressopp/espressopp.svg?branch=master)](https://travis-ci.org/espressopp/espressopp)
[![Code Climate](https://codeclimate.com/github/espressopp/espressopp/badges/gpa.svg)](https://codeclimate.com/github/espressopp/espressopp)
[![codecov.io](https://codecov.io/github/espressopp/espressopp/coverage.svg?branch=master)](https://codecov.io/github/espressopp/espressopp?branch=master)
[![Analysis Status](https://scan.coverity.com/projects/8143/badge.svg?flat=1)](https://scan.coverity.com/projects/espressopp-espressopp)

