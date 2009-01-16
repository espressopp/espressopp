    class ComputerList : public Computer {
      vector < Computer* > compList;
    public:
      virtual void compute(const Real3D dist, 
			   const ParticleRef p1, 
			   const ParticleRef p2) {
	for (vector< Computer* >::iterator it = compList.begin();
	     it != compList.end(); it++)
	  it->compute(dist, p1, p2);
      }
    };
