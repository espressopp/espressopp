// ESPP_CLASS
#ifndef _ANALYSIS_ORDERPARAMETER_HPP
#define _ANALYSIS_ORDERPARAMETER_HPP

#include "mpi.hpp"
#include "types.hpp"
#include "AnalysisBase.hpp"
#include "RealND.hpp"
#include "storage/Storage.hpp"
#include "esutil/Error.hpp"

#include "iterator/CellListIterator.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "Cell.hpp"

#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"

#include "boost/serialization/vector.hpp"
#include "boost/serialization/complex.hpp"

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <algorithm>

using namespace std;
using namespace boost;

typedef complex<double> dcomplex;

// the following constant is not defined everywhere (e.g. not in Mac OS X)
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif


namespace espresso {
  namespace analysis {
    
    using namespace iterator;
    
    // auxiliary class for storing additional properties of each bead for order analysis
    class OrderParticleProps{
    private:
      real d;
      real qlmSumSqrt;
      int nnns;     // number of near neighbors
      
      int ang_m;
      int particle_id;
      
      bool is_solid;   // if true the particle belongs to solid phase
      bool is_surface; // if true the particle belongs to the surface
        
      vector<int> nns;
      vector<dcomplex> qlm; //  depends on angular_momentum =2*angular_momentum+1
      
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
        ar & d;
        ar & qlmSumSqrt;
        ar & nnns;
        ar & ang_m;
        ar & particle_id;
        ar & nns;
        ar & qlm;
        ar & is_solid;
      }
      
    public:
      OrderParticleProps() : d(0),
                             qlmSumSqrt(0),
                             nnns(0),
                             ang_m(0),
                             particle_id(-1) {
        is_solid = false;
        is_surface = false;
      }
      
      OrderParticleProps(int am) : d(0),
                                   qlmSumSqrt(0),
                                   nnns(0),
                                   ang_m(am),
                                   particle_id(-1) {
        qlm = vector<dcomplex>(2*am+1, dcomplex(0.0, 0.0) );
        is_solid = false;
        is_surface = false;
      }
      OrderParticleProps(int am, int pid) : d(0),
                                            qlmSumSqrt(0),
                                            nnns(0),
                                            ang_m(am),
                                            particle_id(pid) {
        qlm = vector<dcomplex>(2*am+1, dcomplex(0.0, 0.0) );
        is_solid = false;
        is_surface = false;
      }
      ~OrderParticleProps(){}
      
      void insertNN(int i){
        nnns++; // increase the number of near neighbors
        nns.push_back( i );
      }
      int getNN(int i){ return nns[i]; }
      
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
      vector<dcomplex> getQlmVector(){ return qlm; }
      void addQlmVector(vector<dcomplex> v){
        if( v.size() != qlm.size() )
          cout<<"Vectors have not the same size. Local: "<< qlm.size() << "  added  "<< v.size() <<endl;
        for(int i=0; i<qlm.size(); i++){
          qlm[ i ] += v[ i ];
        }
      }
      void setQlm(vector<dcomplex> v){
        if( v.size() != qlm.size() )
          cout<<"Vectors have not the same size. Local: "<< qlm.size() << "  new  "<< v.size() <<endl;
        for(int i=0; i<qlm.size(); i++){
          qlm[ i ] = v[ i ];
        }
      }

      void calculateSumQlm(){
        qlmSumSqrt = 0;
        for(vector<dcomplex>::iterator it = qlm.begin(); it!=qlm.end(); ++it){
          qlmSumSqrt += norm( *it );
        }
        qlmSumSqrt = sqrt(qlmSumSqrt);
      }
      real getSumQlm(){ return qlmSumSqrt; }
      
      void setD(real v){ d = v;}
      real getD(){ return d;}
      
      void setPID(int v){ particle_id = v;}
      int getPID(){ return particle_id;}
      
      void setSolid(bool v){ is_solid = v;}
      bool getSolid(){ return is_solid;}
      void setSurface(bool v){ is_surface = v;}
      bool getSurface(){ return is_surface;}
    };
    
    
    
    /** compute order parameter. */
    class OrderParameter : public AnalysisBaseTemplate< RealND > {
    private:
      real cutoff;     // cut off in order to define pairs
      real cutoff_sq;  // cutoff^2
      int angular_momentum;   // angular momentum
      
      vector<OrderParticleProps> opp;   // additional properties
      boost::unordered_multimap <int, OrderParticleProps> opp_map;
      
      boost::unordered_multimap <int, int> pairs;

      /*
       * Cluster analysis.
       */
      bool do_cl_an;  // if true, then cluster analysis will be performed after calculation of d
      bool incl_surface;  // if true, then the surface particle will be included as well
      /*
       * -1 <= d <= 1, thus d_min and d_max should maintain the same property.
       * 
       *  if d_min < d_max, then the solid particle will be defined
       *  in range d_min <= d <= d_max
       * 
       *  if d_min > d_max, then 
       *     -1.0 <= d <= d_min 
       *    d_max <= d <= 1.0
       */
      real d_min, d_max;
      
      
    public:
      static void registerPython();

      OrderParameter(shared_ptr< System > system, 
                     real _cutoff,
                     int _angular_momentum) :
                        AnalysisBaseTemplate< RealND >(system),
                        cutoff(_cutoff),
                        angular_momentum(_angular_momentum){
        cutoff_sq = cutoff * cutoff;
        
        do_cl_an = true;
        incl_surface = true;
        
        d_min = -1.0;
        d_max = -0.75;
      }
      virtual ~OrderParameter() {
      }
      
      dcomplex SphHarm(int l_, int m_, Real3D r_);
      
      void setAngularMomentum(int v){ angular_momentum = v; }
      int getAngularMomentum(){ return angular_momentum; }
      void setCutoff(real v){
        cutoff = v;
        cutoff_sq = cutoff * cutoff;
      }
      real getCutoff(){ return cutoff; }

      // **************************** cluster analysis
      void setDo_cl_an(bool v){ do_cl_an = v; }
      bool getDo_cl_an(){ return do_cl_an; }

      void setIncl_surface(bool v){ incl_surface = v; }
      bool getIncl_surface(){ return incl_surface; }
      
      void setD_min(int v){ d_min = v; }
      int getD_min(){ return d_min; }
      void setD_max(int v){ d_max = v; }
      int getD_max(){ return d_max; }
      
      /*
       * It is efficient only when communications are optimized.
       */
      RealND computeRaw() {
        
        shared_ptr< storage::Storage > stor = getSystem()->storage;
        shared_ptr< mpi::communicator > cmm = getSystem()->comm;
        int this_node = getSystem() -> comm -> rank();
        
        // ------------------------------------------------------------------------------
        // iterate over local particles, create a map of additional properties
        CellList cells_loc = stor->getLocalCells();
        for(CellListIterator cit(cells_loc); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          opp_map.insert( make_pair(p.id(), OrderParticleProps(angular_momentum, p.id()) ) );
        }
        
        // ------------------------------------------------------------------------------
        // create pairs
        CellList cells_real = stor->getRealCells();
        for (CellListAllPairsIterator it(cells_real); it.isValid(); ++it) {
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
        
        // ------------------------------------------------------------------------------
        // here communicate, send all ghost info 
        //   TODO not the best way all to all communication.
        vector <OrderParticleProps> sendGhostInfo;
        for(boost::unordered_multimap<int, OrderParticleProps>::iterator opm = opp_map.begin(); opm!=opp_map.end(); ++opm){
          int id = (*opm).first;
          if( !getSystem()->storage->lookupRealParticle(id) ){
            OrderParticleProps &op = (*opm).second;
            if( !op.getNumNN()==0 ) sendGhostInfo.push_back( (*opm).second );
          }
        }
        
        int maxSize, vecSize  = sendGhostInfo.size();
        mpi::all_reduce( *cmm, vecSize, maxSize, mpi::maximum<int>() );
        while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );

        vector< OrderParticleProps > totID;
        boost::mpi::all_gather( *getSystem()->comm, &sendGhostInfo[0], maxSize, totID);
        
        // TODO use set instead of vector
        vector<int> realHere;
        for(vector<OrderParticleProps>::iterator it = totID.begin(); it!=totID.end(); ++it){
          OrderParticleProps &gop = *it;
          if( gop.getPID()!=-1 && stor->lookupRealParticle( gop.getPID() ) ){
            if( find( realHere.begin(), realHere.end(), gop.getPID() ) == realHere.end() )
              realHere.push_back( gop.getPID() );
            
            OrderParticleProps *opp_i = &(opp_map.find( gop.getPID() ))->second;
            opp_i->addQlmVector( gop.getQlmVector() );
            int numPadd = gop.getNumNN();
            for(int  i = 0; i< numPadd; i++){
              opp_i->insertNN( gop.getNN(i) );
            }
          }
        }
        
        // ------------------------------------------------------------------------------
        // loop over all real particles and calculate SumQlm
        for(CellListIterator cit(cells_real); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          OrderParticleProps *opp_i = &( opp_map.find( p.id() ) )->second;
          opp_i->calculateSumQlm();
        }
        
        // ------------------------------------------------------------------------------
        // communicate back info to ghost particles
        sendGhostInfo.clear();
        for(vector<int>::iterator it = realHere.begin(); it!=realHere.end(); ++it){
          sendGhostInfo.push_back( opp_map.find( *it )->second );
        }
        
        maxSize, vecSize  = sendGhostInfo.size();
        mpi::all_reduce( *cmm, vecSize, maxSize, mpi::maximum<int>() );
        while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );
        
        totID.clear();
        boost::mpi::all_gather( *cmm, &sendGhostInfo[0], sendGhostInfo.size(), totID);
        
        for(vector<OrderParticleProps>::iterator it = totID.begin(); it!=totID.end(); ++it){
          OrderParticleProps &gop = *it;
          if( gop.getPID()!=-1 && stor->lookupGhostParticle( gop.getPID() ) ){
            (opp_map.find( gop.getPID() ))->second = gop;
          }
        }
        
        // ------------------------------------------------------------------------------
        // loop over pairs
        for(boost::unordered_multimap<int, int>::iterator pit = pairs.begin(); pit != pairs.end(); ++pit){
          int first = (*pit).first;
          int second = (*pit).second;
          OrderParticleProps *opp_i1 = &(opp_map.find( first ))->second;
          OrderParticleProps *opp_i2 = &(opp_map.find( second ))->second;
          
          // checking 0 sum Qlm
          if(opp_i2->getSumQlm() == 0){
            cerr<<" Bead2: "<< second << " in pair has 0 sum: "<< opp_i2->getSumQlm()<<endl;
          }
          if(opp_i1->getSumQlm() == 0){
            cerr<<" Bead1: "<< first << " in pair has 0 sum: "<< opp_i1->getSumQlm()<<endl;
          }
          
          for (int m = -angular_momentum; m <= angular_momentum; m++) {
            real d1 = (opp_i1->getQlm(m) * conj( opp_i2->getQlm(m) )).real() / opp_i2->getSumQlm();
            opp_i1->setD( opp_i1->getD() + d1 );
            real d2 = (opp_i2->getQlm(m) * conj( opp_i1->getQlm(m) )).real() / opp_i1->getSumQlm();
            opp_i2->setD( opp_i2->getD() + d2 );
          }
        }

        // ------------------------------------------------------------------------------
        // communicate D information
        sendGhostInfo.clear();
        for(boost::unordered_multimap<int, OrderParticleProps>::iterator opm = opp_map.begin(); opm!=opp_map.end(); ++opm){
          int id = (*opm).first;
          if( getSystem()->storage->lookupGhostParticle(id) ){
            OrderParticleProps &op = (*opm).second;
            if( !op.getNumNN()==0 ) sendGhostInfo.push_back( (*opm).second );
          }
        }
        
        maxSize, vecSize  = sendGhostInfo.size();
        mpi::all_reduce( *cmm, vecSize, maxSize, mpi::maximum<int>() );
        while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );
        totID.clear();
        boost::mpi::all_gather( *getSystem()->comm, &sendGhostInfo[0], maxSize, totID);
        
        for(vector<OrderParticleProps>::iterator it = totID.begin(); it!=totID.end(); ++it){
          OrderParticleProps &gop = *it;
          if( gop.getPID()!=-1 && getSystem()->storage->lookupRealParticle( gop.getPID() ) ){
            OrderParticleProps *opp_i = &(opp_map.find( gop.getPID() ))->second;
            opp_i->setD(  opp_i->getD() + gop.getD() );
          }
        }
        
        // ------------------------------------------------------------------------------
        // loop over particles and normalize d
        for(CellListIterator cit(cells_real); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps *opp_i = &(opp_map.find( pid ))->second;
          
          //catch particles without neighbors 
          if ( opp_i->getNumNN() == 0 ){
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
           
          //ret[pid] = opp_i->getD();
        }
        
        RealND ret(1, 0.0);
        return ret;
      }
      
      void define_solid(){
        CellList cells = getSystem()->storage->getRealCells();
        
        int count = 0;
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps *opp_i = &(opp_map.find( pid ))->second;
          if( opp_i->getD() >= d_min && opp_i->getD() <= d_max ){
            opp_i->setSolid( true );
            count++;
          }
        }
        cout<< " Number of solid particles without surface: "<< count<< endl;
        
        if(incl_surface){
          int count1 = 0;
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
            Particle& p = *cit;
            int pid = p.id();
            OrderParticleProps *opp_i = &(opp_map.find( pid ))->second;
            if( opp_i->getSolid() && !opp_i->getSurface() ){
              int nnn = opp_i->getNumNN();
              for(int i=0; i<nnn; i++){
                OrderParticleProps *opp_i_surf = &(opp_map.find( opp_i->getNN(i) ))->second;
                if( !opp_i_surf->getSolid() ){
                  //opp_i_surf->setSolid(true);
                  opp_i_surf->setSurface(true);
                  count1++;
                }
              }
            }
            
          }
          cout<< " count1: "<< count1 << endl;
          
          // -----------------------------------------------------------------------
          // send ghost info
          shared_ptr< mpi::communicator > cmm = getSystem()->comm;
          
          vector <OrderParticleProps> sendGhostInfo;
          for(boost::unordered_multimap<int, OrderParticleProps>::iterator opm = opp_map.begin(); opm!=opp_map.end(); ++opm){
            int id = (*opm).first;
            if( !getSystem()->storage->lookupRealParticle(id) ){
              OrderParticleProps &op = (*opm).second;
              if( !op.getNumNN()==0 ) sendGhostInfo.push_back( (*opm).second );
            }
          }

          int maxSize, vecSize  = sendGhostInfo.size();
          mpi::all_reduce( *cmm, vecSize, maxSize, mpi::maximum<int>() );
          while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );

          vector< OrderParticleProps > totID;
          boost::mpi::all_gather( *getSystem()->comm, &sendGhostInfo[0], maxSize, totID);

          for(vector<OrderParticleProps>::iterator it = totID.begin(); it!=totID.end(); ++it){
            OrderParticleProps &gop = *it;
            if( gop.getPID()!=-1 && getSystem()->storage->lookupRealParticle( gop.getPID() ) ){
              OrderParticleProps *opp_i_loc = &(opp_map.find( gop.getPID() ))->second;
              opp_i_loc->setSolid( gop.getSolid() );
            }
          }
          //--------------------------------------------------------------------------
        }
        cout<< " Number of solid particles: "<< count<< endl;
        
        
      }
      
      void cluster_analysis(){
        
      }
      
      python::list compute() {
        python::list ret;
        
        RealND res = computeRaw();
        
        if( do_cl_an ){
          define_solid();
          
          cluster_analysis();
        }
        
        
        CellList cells = getSystem()->storage->getRealCells();
        int this_node = getSystem() -> comm -> rank();
        
        vector <OrderParticleProps> sendGhostInfo;
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps opp_i = (opp_map.find( pid ))->second;
          sendGhostInfo.push_back( opp_i );
        }
        
        vector<int> sendSizes;
        int ss  = sendGhostInfo.size();
        boost::mpi::all_gather( *getSystem()->comm, ss, sendSizes);
        int maxSize=0;
        for(vector<int>::iterator it = sendSizes.begin(); it!=sendSizes.end(); ++it){
          maxSize = max(maxSize, *it);
        }
        while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );

        int numProc = getSystem()->comm->size();
        vector<OrderParticleProps> totID = vector<OrderParticleProps>( numProc * maxSize, OrderParticleProps() );
        boost::mpi::gather( *getSystem()->comm, &sendGhostInfo[0], maxSize, totID, 0);

        if( getSystem()->comm->rank()==0 ){
          for(vector<OrderParticleProps>::iterator it=totID.begin(); it!=totID.end(); ++it) {
            OrderParticleProps &gop = *it;
            if(gop.getPID()>=0){
              python::list nn;
              
              for(int i=0; i<gop.getNumNN(); i++){
                nn.append( gop.getNN(i) );
              }
              
              //cout<< "Pid: "<< gop.getPID() << "  solid:  "<< gop.getSolid() << endl;
              python::tuple pt = python::make_tuple( gop.getPID(), gop.getD(),
                      gop.getNumNN(), gop.getSolid(), gop.getSurface(), nn);
              ret.append( pt );
            }
          }

          ret.sort();
        }
        else{
          ret.append(0);
        }
        
        return ret;
      }

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
