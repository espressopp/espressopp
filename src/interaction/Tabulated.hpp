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

    class Tabulated: public PotentialTemplate <Tabulated> {

        private:
            std::string filename;
            shared_ptr <InterpolationTable> table;

        public:
            static void registerPython();
         
            Tabulated() {
                setShift(0.0);
                setCutoff(infinity);
            }
         
            // used for fixedpairlist (2-body bonded interaction)
            Tabulated(const char* filename){
                setFilename(filename);
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
         
            bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
                real ffactor = table->getForce(sqrt(distSqr));
                force = dist * ffactor;
                return true;
            }

    };//class
  }
}

#endif
