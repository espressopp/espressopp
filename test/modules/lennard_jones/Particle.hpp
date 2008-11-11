//a particle is defined by a point in R3
class Particle {
private:
  double x;
  double y;
  double z;
public:
  Particle() {};
  Particle(double xx, double yy, double zz);
  std::string toString();
  double getx() {return x;};
  double gety() {return y;};
  double getz() {return z;};
};
