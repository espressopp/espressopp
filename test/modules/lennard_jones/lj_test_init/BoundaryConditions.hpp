int sign(double _r) {
    if(_r > 0)
      return 1;
    else
      return -1;
  };

class BoundaryConditions {};

class Cube:public BoundaryConditions {
  private:
  double L;
  double half_L;

  public:
  Cube() {};
  void setSide(double _L) {L=_L; half_L=_L/2;};
  double getSide() {return L;};
  double getHalfSide() {return half_L;};
  double getMinimumImageDistance(Particle _Pi, Particle _Pj) {
    double xij;
    double yij;
    double zij;
    double xijabs;
    double yijabs;
    double zijabs;

    xij = _Pi.getx() - _Pj.getx();
    yij = _Pi.gety() - _Pj.gety();
    zij = _Pi.getz() - _Pj.getz();

    xijabs = fabs(xij);
    yijabs = fabs(yij);
    zijabs = fabs(zij);

    if(xijabs > half_L) xij = (xijabs - L) * sign(xij);
    if(yijabs > half_L) yij = (yijabs - L) * sign(yij);
    if(zijabs > half_L) zij = (zijabs - L) * sign(zij);
    
    return (xij * xij + yij * yij + zij * zij);
  }
};
