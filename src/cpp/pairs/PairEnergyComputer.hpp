class EnergyComputer : public ParticlePairComputer {
  real totalEnergy;
  const PropertyRef< real > energy;
  const Interaction &interaction;

public:
  virtual void compute(const Real3D dist, 
		       const ParticleRef p1, 
		       const ParticleRef p2) {
    real e = interaction.computeEnergy(dist, p1, p2);
    energy[p1] += 0.5*e;
    energy[p2] += 0.5*e;
    totalEnergy += e;
  }
};
