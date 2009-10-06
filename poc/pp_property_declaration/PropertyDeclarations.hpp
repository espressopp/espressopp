#if defined(PROPERTY) || defined(LOCALPROPERTY)

#ifndef PROPERTY
#define PROPERTY(name, type)
#endif
#ifndef LOCALPROPERTY
#define LOCALPROPERTY(name, type)
#endif

#define PROPERTIES_HPP_DECLARATIONS

#include "Properties.hpp"

#undef PROPERTY_HPP_DECLARATIONS

#undef PROPERTY
#undef LOCALPROPERTY

#endif
