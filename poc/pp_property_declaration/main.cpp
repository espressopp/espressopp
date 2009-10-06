#include "Properties.hpp"

/************************************************
 * build particle struct
 ************************************************/

struct Particle {
  struct Basic {

#define PROPERTY(name, descriptor) descriptor::Type name;
#include "PropertyDeclarations.hpp"

  } basic;

  struct Local {

#define LOCALPROPERTY(name, descriptor) descriptor::Type name;
#include "PropertyDeclarations.hpp"

  } local;
};

/************************************************
 * set initial value
 ************************************************/


int main()
{
  Particle testParticle;

#define PROPERTY(name, descriptor)      descriptor::init(testParticle.basic.name);
#define LOCALPROPERTY(name, descriptor) descriptor::init(testParticle.local.name);
#include "PropertyDeclarations.hpp"

  return 0;
}
