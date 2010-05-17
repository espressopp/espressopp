// ESPP_CLASS
#ifndef _INTERACTION_TABULATED_HPP
#define _INTERACTION_TABULATED_HPP

#include "Potential.hpp"
#include "InterpolationTable.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	a tabulated potential.

        The potential and forces must be provided in a file.
    */
    class Tabulated : public PotentialTemplate< Tabulated > {
    private:

      std::string filename;
      InterpolationTable table;

    public:

      static void registerPython();

      Tabulated() {
	setShift(0.0);
	setCutoff(infinity);
      }

      Tabulated(const char* filename, real cutoff) {
        setFilename(filename);
	setShift(0.0);
	setCutoff(cutoff);
      }

      // Setter and getter for filename

      void setFilename(const char* _filename);

      const char* getFilename() const { return filename.c_str(); }

      real _computeEnergySqrRaw(real distSqr) const {
        // make an interpolation
	return table.getEnergy(sqrt(distSqr));
      }

      bool _computeForceRaw(Real3DRef force,
			    ConstReal3DRef dist,
			    real distSqr) const {
	real ffactor = table.getForce(sqrt(distSqr));
	force = dist * ffactor;
	return true;
      }
    };
  }
}

#endif
