#ifndef _PMI_HPP
#define _PMI_HPP

#include "acconfig.hpp"

#ifdef HAVE_MPI
#include "pmi/pmi.hpp"
#else
#define PMI_REGISTER_CLASS(a,b)
#define PMI_REGISTER_METHOD(a,b,c)
#endif

#endif
