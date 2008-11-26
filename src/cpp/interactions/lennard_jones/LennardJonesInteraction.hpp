#include "Interaction.hpp"

//define LennardJonesInteraction class
class LennardJonesInteraction:public Interaction {
private:
  double sigma;
  double epsilon;
public:
  LennardJonesInteraction() {};
  LennardJonesInteraction(double _rc) {rc=_rc; rcsq=pow(_rc,2);};
  double computeLennardJonesEnergy(double _rsq);
};
