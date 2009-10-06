#ifdef PROPERTIES_HPP_DECLARATIONS

/*************************************
 * Property declarations
 *************************************/

PROPERTY(type,     PropertyType)
PROPERTY(pos,      PropertyPosition)
LOCALPROPERTY(img, PropertyImg)

#endif

#ifndef PROPERTIES_HPP_DEFINITIONS
#define PROPERTIES_HPP_DEFINITIONS

/*************************************
 * Property definitions
 *************************************/

struct PropertyType {
  typedef int Type;
  static void init(Type &v) { v = 1; };
};

struct PropertyPosition {
  typedef float Type[3];
  static void init(Type &v) { v[0] = 0; v[1] = 1; v[2] = 2; };
};

struct PropertyImg {
  typedef int Type[3];
  static void init(Type &v) { v[0] = 3; v[1] = 4; v[2] = 5; };
};

#endif
