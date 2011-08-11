// ESPP_CLASS
#ifndef _INTERACTION_TABULATED_HPP
#define _INTERACTION_TABULATED_HPP

//#include <stdexcept>
#include "Potential.hpp"
#include "Interpolation.hpp"

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
            shared_ptr <Interpolation> table;

        public:
            static void registerPython();
         
            Tabulated() {
                setShift(0.0);
                setCutoff(infinity);
                //std::cout << "using default tabulated potential ...\n";
            }
         
            // used for fixedpairlist (2-body bonded interaction)
            Tabulated(int itype, const char* filename){
                setFilename(itype, filename);
                setShift(0.0);
                setCutoff(infinity);
                std::cout << "using tabulated potential " << filename << "\n";
            }
         
            Tabulated(int itype, const char* filename, real cutoff) {
                setFilename(itype, filename);
                setShift(0.0);
                setCutoff(cutoff);
                std::cout << "using tabulated potential " << filename << "\n";
            }
         
            /** Setter for the filename will read in the table. */
            void setFilename(int itype, const char* _filename);
         
            /** Getter for the filename. */
            const char* getFilename() const { return filename.c_str(); }
         
            real _computeEnergySqrRaw(real distSqr) const {
                // make an interpolation
                if (table) 
                    return table->getEnergy(sqrt(distSqr));
                else {
                    throw std::runtime_error("Tabulated potential table not available.");
                    //return 0.0;
                }
            }
         
            bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
                real ffactor;
                if (table) 
                   ffactor = table->getForce(sqrt(distSqr));
                else {
                    throw std::runtime_error("Tabulated potential table not available.");
                    //return 0.0;
                }
                force = dist * ffactor;
                return true;
            }

    };//class
  }
}

#endif
