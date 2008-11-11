#include "Interaction.hpp"

//define LennardJonesInteraction class
class LennardJonesInteraction:public Interaction {
private:
  double sigma;
  double epsilon;
public:
  LennardJonesInteraction() {};
  LennardJonesInteraction(double rcut) {rc=rcut; rc2=pow(rcut,2);};
  double computeLennardJonesEnergy(double dist2);
};
