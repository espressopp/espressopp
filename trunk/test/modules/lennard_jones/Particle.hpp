//a particle is defined by a point in R3
class Particle {
private:
  double x;
  double y;
  double z;
public:
  Particle() {};
  Particle(double _x, double _y, double _z);
  std::string toString();
  double getx() {return x;};
  double gety() {return y;};
  double getz() {return z;};
  void setx(double _x) {x = _x;};
  void sety(double _y) {y = _y;};
  void setz(double _z) {z = _z;};
};
