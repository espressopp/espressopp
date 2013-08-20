// ESPP_CLASS
#ifndef _ANALYSIS_ORDERPARAMETER_HPP
#define _ANALYSIS_ORDERPARAMETER_HPP

#include "types.hpp"
#include "AnalysisBase.hpp"
#include "RealND.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "Cell.hpp"

#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"

#include "boost/serialization/vector.hpp"
#include "boost/serialization/complex.hpp"

#include <boost/math/special_functions/spherical_harmonic.hpp>

using namespace std;
using namespace boost;

typedef complex<double> dcomplex;

// the following constant is not defined everywhere (e.g. not in Mac OS X)
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif


/*
 * 
 * Currently code is not parallel!!!
 * 
 */

namespace espresso {
  namespace analysis {
    using namespace iterator;
    class OrderParticleProps{
    protected:
      vector<dcomplex> qlm; //  depends on angular_momentum =2*angular_momentum+1
      real d;
      real qlmSumSqrt;
      int nnns;     // number of near neighbors
      
      real nbond;
      
      int ang_m;
      
      //boost::unordered_multiset< longint > nns;        // set of near neighbor id's
      vector<longint> nns;
      
      bool solid;
      
    public:
      OrderParticleProps(int am){
        qlm = vector<dcomplex>(2*am+1, dcomplex(0.0, 0.0) );
        d = 0;
        qlmSumSqrt = 0;
        nnns = 0;
        
        ang_m = am;
        
        solid = false;
      }
      ~OrderParticleProps(){}
      
      void insertNN(longint i){
        nnns++; // increase the number of near neighbors
        //nns.insert( i );
        nns.push_back( i );
      }
      longint getNN(int i){
        return nns[i];
      }
      
      int getNumNN(){ return nnns; }
      
      // index defined as -l, -l+1, ..., 0, ..., l-1, l for convenience 
      void setQlm(int indx, dcomplex v){
        int hh = indx + ang_m;
        if(hh<0 || hh>=2*ang_m+1) cout<<"OUT OF RANGE!!"<<endl;
        qlm[ indx + ang_m ] = v;
      }
      dcomplex getQlm(int indx){
        int hh = indx + ang_m;
        if(hh<0 || hh>=2*ang_m+1) cout<<"OUT OF RANGE!!"<<endl;
        return qlm[ indx + ang_m ];
      }

      void calculateSumQlm(){
        for(vector<dcomplex>::iterator it = qlm.begin(); it!=qlm.end(); ++it){
          qlmSumSqrt += norm( *it );
        }
        qlmSumSqrt = sqrt(qlmSumSqrt);
      }
      real getSumQlm(){
        return qlmSumSqrt;
      }
      
      void setD(real v){ d = v;}
      real getD(){ return d;}
      
      void setSolid(bool v){ solid = v;}
      bool getSolid(){ return solid;}
    };
    
    /** Class to compute order parameter. */
    class OrderParameter : public AnalysisBaseTemplate< RealND > {
    private:
      real cutoff;     // cut off in order to define pairs
      real cutoff_sq;  // cutoff^2
      real threshold;  // for local order parameter
      
      int num_solid;   // number of solid particles
      
      vector<OrderParticleProps> opp;   // additional properties
      boost::unordered_multimap <longint, OrderParticleProps> opp_map;
      
      boost::unordered_multimap <longint, longint> pairs;
      
      int angular_momentum;   // angular momentum
    public:
      static void registerPython();

      // important - verlet list should be separate from any other vl in system
      // could be done in python as well as disconnection
      OrderParameter(shared_ptr< System > system, 
                     real _cutoff,
                     int _angular_momentum,
                     real _threshold
                     ) :
                     AnalysisBaseTemplate< RealND >(system),
                     cutoff(_cutoff),
                     angular_momentum(_angular_momentum),
                     threshold(_threshold){
        cutoff_sq = cutoff * cutoff;
      }
      virtual ~OrderParameter() {
      }
      
      dcomplex SphHarm(int l_, int m_, Real3D r_);
      
      int getAngularMomentum(){ return angular_momentum; }
      void setAngularMomentum(int v){ angular_momentum = v; }
      real getCutoff(){ return cutoff; }
      void setCutoff(real v){
        cutoff = v;
        cutoff_sq = cutoff * cutoff;
      }
      real getThreshold(){ return threshold; }
      void setThreshold(real v){ threshold = v; }
      
      RealND computeRaw() {
        
        // number of local particles
        int localN = getSystem()->storage->getNRealParticles();
        
        RealND ret(localN, 0.0);
        
        // additional properties for each particle
        opp = vector<OrderParticleProps>(localN, OrderParticleProps(angular_momentum) );
        // iterate over local particles
        CellList cells = getSystem()->storage->getRealCells();
        vector<OrderParticleProps>::iterator opp_it = opp.begin();
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          opp_map.insert( make_pair(p.id(), *opp_it) );
          ++opp_it;
        }

        for (CellListAllPairsIterator it(cells); it.isValid(); ++it) {
          Real3D r = it->first->position() - it->second->position();
          real dist_sq = r.sqr();
          // setup NN list within cutoff
          if (dist_sq <= cutoff_sq){
            pairs.insert( make_pair( it->first->id(), it->second->id() ) );
            
            OrderParticleProps *opp_i1 = &(opp_map.find( it->first->id() ))->second;
            OrderParticleProps *opp_i2 = &(opp_map.find( it->second->id() ))->second;
            opp_i1->insertNN( it->second->id() );
            opp_i2->insertNN( it->first->id() );
            
            dcomplex tmpVar = SphHarm(angular_momentum, 0, r);
            dcomplex tmpVar1 = SphHarm(angular_momentum, 0, (-1)*r);
            opp_i1->setQlm( 0, (opp_i1-> getQlm(0) + tmpVar) );
            opp_i2->setQlm( 0, (opp_i2-> getQlm(0) + tmpVar1) );
            
            //
            for (int m = 1; m <= angular_momentum; m++) {
                tmpVar = SphHarm(angular_momentum, m, r);
                tmpVar1 = SphHarm(angular_momentum, m, (-1)*r);
                
                opp_i1->setQlm( m, (opp_i1->getQlm(m) + tmpVar) );
                opp_i2->setQlm( m, (opp_i2->getQlm(m) + tmpVar1) );
                
                dcomplex conj_tmpVar = pow(-1.0, (real)m) * conj( tmpVar );
                dcomplex conj_tmpVar1 = pow(-1.0, (real)m) * conj( tmpVar1 );
                opp_i1->setQlm( -m, (opp_i1->getQlm(-m) + conj_tmpVar) );
                opp_i2->setQlm( -m, (opp_i2->getQlm(-m) + conj_tmpVar1) );
            }
            
          }
        }
        
        // loop over particles and calc SumQlm
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          OrderParticleProps *opp_i = &(opp_map.find( p.id() ))->second;
          opp_i->calculateSumQlm();
        }
        
        // loop over pairs
        for(boost::unordered_multimap<longint, longint>::iterator pit = pairs.begin(); pit != pairs.end(); ++pit){
          longint first = (*pit).first;
          longint second = (*pit).second;
          OrderParticleProps *opp_i1 = &(opp_map.find( first ))->second;
          OrderParticleProps *opp_i2 = &(opp_map.find( second ))->second;
          
          for (int m = -angular_momentum; m <= angular_momentum; m++) {
            real d1 = (opp_i1->getQlm(m) * conj( opp_i2->getQlm(m) )).real() / opp_i2->getSumQlm();
            opp_i1->setD( opp_i1->getD() + d1 );
            real d2 = (opp_i2->getQlm(m) * conj( opp_i1->getQlm(m) )).real() / opp_i1->getSumQlm();
            opp_i2->setD( opp_i2->getD() + d2 );
          }
        }
        
        // loop over particles and normalize d
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps *opp_i = &(opp_map.find( pid ))->second;
          
          //catch particles without neighbors 
          if ( opp_i->getNumNN() == 0){
            cerr << "Bead: "<< pid << " has no neighbors - setting OP to zero" << endl;
            opp_i->setD( 0.0 );
            continue;
          }
          
          opp_i->setD( opp_i->getD() / opp_i->getSumQlm() / opp_i->getNumNN() );
          
          //check for NAN
          if ( opp_i->getD() != opp_i->getD() ) {
              cerr << "Current bead: " << p.id() << endl;
              cerr << "qlmSumSqrt: " << opp_i->getSumQlm() << endl;
              cerr << "near neibs: " << opp_i->getNumNN() << endl;
              throw std::runtime_error("result: NAN - this should never happen..." );
          }
          
          // assigning to the list of solid particles
          if( opp_i->getD() >= threshold ){
            num_solid++;
            opp_i->setSolid( true );
          }
          
          ret[pid] = opp_i->getD();
        }
        
        return ret;
      }
      
      python::list compute() {
        python::list ret;
        
        RealND res = computeRaw();
        
        CellList cells = getSystem()->storage->getRealCells();
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps *opp_i = &(opp_map.find( pid ))->second;
          
          python::tuple pt = python::make_tuple( pid, opp_i->getD(), opp_i->getNumNN());
          
          //ret.insert( pid, pt );
          ret.append( pt );
        }
        
        ret.sort();
        
        return ret;
      }

      // functions below should be revised. it doesn't work good for RealND
      
      python::list getAverageValue() {
        python::list ret;
        /*
        RealND res;
        res = nMeasurements>0 ? newAverage : 0;
        ret.append(res);
        res = nMeasurements>0 ? newVariance : 0;
        ret.append(sqrt(res/(nMeasurements-1)));
        */
        return ret;
      }
      
      void resetAverage() {
        newAverage   = 0;
        lastAverage  = 0;
        newVariance  = 0;
        lastVariance = 0;
      }

      void updateAverage(RealND res) {
        /*
    	if (nMeasurements > 0) {
    	  if (nMeasurements == 1) {
              newAverage     = res;
              lastAverage    = newAverage;
          } else {
              newAverage   = lastAverage  + (res - lastAverage) / nMeasurements;
              newVariance  = lastVariance + (res - lastAverage) * (res - newAverage);
              lastAverage  = newAverage;
              lastVariance = newVariance;
          }
    	}
        */
        return;
      }
    };


    
  }
}

#endif
