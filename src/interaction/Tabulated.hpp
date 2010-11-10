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

        Be careful: default and copy constructor of this class are used.
    */

    class Tabulated : public PotentialTemplate< Tabulated > {

    private:

      std::string filename;

      shared_ptr< InterpolationTable > table;

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

      /** Setter for the filename will read in the table. */

      void setFilename(const char* _filename);

      /** Getter for the filename. */

      const char* getFilename() const { return filename.c_str(); }

      real _computeEnergySqrRaw(real distSqr) const {
        // make an interpolation
	return table->getEnergy(sqrt(distSqr));
      }

      bool _computeForceRaw(real force[3],
                            const real dist[3],
                            real distSqr) const {

	real ffactor = table->getForce(sqrt(distSqr));
        force[0] = dist[0] * ffactor;
        force[1] = dist[1] * ffactor;
        force[2] = dist[2] * ffactor;
        return true;
      }

    };
  }
}

#endif
