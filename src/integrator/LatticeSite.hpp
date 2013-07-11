// ESPP_CLASS
#ifndef _INTEGRATOR_LATTICEMODEL_HPP
#define _INTEGRATOR_LATTICEMODEL_HPP

namespace espresso {
  namespace integrator {
    class LBSite {
      /**
      * \brief Description of the properties of the Site class
      *
      * This is a Site class for a Lattice Boltzmann method. It includes
      * declaration of the populations, local and equilibrium moments (MRT model).
      * The functions to calculate them are defined. Also, back-transformation
      * functions are included.
      *
      * In short, all the quantities and operations taking place at the lattice site
      * are defined/listed/declared here.
      *
      * This is a first try. Further extensions might be needed.
      */
      public:

        LBSite (int _numVels, real _a, real _tau);
        ~LBSite ();

        /* SET AND GET DECLARATION */
        void setF_i (int _i, real _f);	// set f_i population to _f
        real getF_i (int _i);		// get f_i population

        void setM_i (int _i, real _m);	// set m_i moment to _m
        real getM_i (int _i);		// get m_i moment

        void setMeq_i (int _i, real _meq);	// set meq_i moment to _meq
        real getMeq_i (int _i);		// get meq_i moment

        void setInvB (int _i, real _b);  // set invLov_b value to _b
        real getInvB (int _i);   // get invLoc_b value

        void setEqWLoc (int _i, real _w);  // set eqWeightLoc value to _w
        real getEqWLoc (int _i);   // get eqWeightLoc value

        void setGammaB (real _gamma_b); // set gamma for bulk
        real getGammaB ();  // get gamma for bulk

        void setGammaS (real _gamma_s); // set gamma for shear
        real getGammaS ();  // get gamma for shear

        void setGammaOdd (real _gamma_odd); // set gamma odd
        real getGammaOdd ();  // get gamma odd

        void setGammaEven (real _gamma_even); // set gamma even
        real getGammaEven ();  // get gamma even

        void setALocal (real _a); // set aLocal
        real getALocal ();    // get aLocal

        void setTauLocal (real _tau); // set tauLocal
        real getTauLocal ();  // get tauLocal
        /* END OF SET AND GET DECLARATION */

        void scaleF_i (int _i, real _value);  // scale f_i population by _value
        void scaleM_i (int _i, real _value);  // scale m_i moment by _value

        /* FUNCTIONS DECLARATION */
        void calcLocalMoments ();	// calculate local moments
        void calcEqMoments ();		// calculate equilibrium moments
        void relaxMoments (int _numVels);		// relax loc. moments towards eq.values
	      void btranMomToPop (int _numVels);		// back-transform moments to populations

      private:
        std::vector<real> f;        // populations on a lattice site
        std::vector<real> m;        // moments on a site
        std::vector<real> meq;      // eq. moments on a site
        static real aLocal;        // local variable for lattice spacing
        static real tauLocal;      // local variable for lattice time
        static real gamma_b, gamma_s, gamma_odd, gamma_even;  // gammas
        static std::vector<real> invLoc_b;    // local inverse coefficients b_i
        static std::vector<real> eqWeightLoc; // local eq. weights
    };

    class GhostLattice {
      public:
        GhostLattice (int _numVels);
        ~GhostLattice ();

        void setPop_i (int _i, real _pop);  // set f_i population to _f
        real getPop_i (int _i);             // get f_i population
      private:
        std::vector<real> pop;
    };
  }
}

#endif
