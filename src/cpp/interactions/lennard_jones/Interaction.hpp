//base class

#include "types.hpp"

class Interaction {
protected:
  real rc;
  real rc2;
public:
  Interaction() {}
  Interaction(real rcut) {rc=rcut; rc2=pow(rcut,2);};
  virtual ~Interaction() {}
  virtual real computeEnergy(real dist2) const = 0;
  real getCutoff() const {return rc;}
  void setCutoff(real rcut) {rc=rcut; rc2=pow(rcut,2);}
};
