//base class
class Interaction {
protected:
  double rc;
  double rcsq;
public:
  Interaction() {};
  Interaction(double _rc) {rc=_rc; rcsq=pow(_rc,2);};
  double getCutoff() {return rc;};
  double getCutoffSq() {return rcsq;};
  void setCutoff(double _rc) {rc=_rc; rcsq=pow(_rc,2);};
};
