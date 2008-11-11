//base class
class Interaction {
protected:
  double rc;
  double rc2;
public:
  Interaction() {};
  Interaction(double rcut) {rc=rcut; rc2=pow(rcut,2);};
  double getCutoff() {return rc;};
  void setCutoff(double rcut) {rc=rcut; rc2=pow(rcut,2);};
};
