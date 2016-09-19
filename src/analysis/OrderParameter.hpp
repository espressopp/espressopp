/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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

// #include <boost/math/special_functions/spherical_harmonic.hpp>

#include <algorithm>

using namespace std;
using namespace boost;

typedef complex<espressopp::real> dcomplex;

// the following constant is not defined everywhere (e.g. not in Mac OS X)
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif


// TODO create the total map: cpu <-> pid

namespace espressopp {
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

      int label; // variable for cluster analysis

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
        ar & is_surface;
        ar & label;
      }

    public:
      OrderParticleProps() : d(0),
                             qlmSumSqrt(0),
                             nnns(0),
                             ang_m(0),
                             particle_id(-1),
                             label(-1){
        is_solid = false;
        is_surface = false;
      }

      OrderParticleProps(int am) : d(0),
                                   qlmSumSqrt(0),
                                   nnns(0),
                                   ang_m(am),
                                   particle_id(-1),
                                   label(-1){
        qlm = vector<dcomplex>(2*am+1, dcomplex(0.0, 0.0) );
        is_solid = false;
        is_surface = false;
      }
      OrderParticleProps(int am, int pid) : d(0),
                                            qlmSumSqrt(0),
                                            nnns(0),
                                            ang_m(am),
                                            particle_id(pid),
                                            label(-1){
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

      void setLabel(int v){ label = v;}
      int getLabel(){ return label;}
    };



    /** compute order parameter. */
    class OrderParameter : public AnalysisBaseTemplate< RealND > {
    private:
      real cutoff;     // cut off in order to define pairs
      real cutoff_sq;  // cutoff^2
      int angular_momentum;   // angular momentum

      //vector<OrderParticleProps> opp;   // additional properties
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

      int max_cluster, number_of_clust; // a size of a biggest cluster and number of clusters

    public:
      static void registerPython();

      OrderParameter(shared_ptr< System > system,
                     real _cutoff,
                     int _angular_momentum,
                     bool _do_cl_an,
                     bool _incl_surface,
                     real _d_min,
                     real _d_max) :
                        AnalysisBaseTemplate< RealND >(system),
                        cutoff(_cutoff),
                        angular_momentum(_angular_momentum),
                        do_cl_an(_do_cl_an),
                        incl_surface(_incl_surface),
                        d_min(_d_min),
                        d_max(_d_max){
        cutoff_sq = cutoff * cutoff;
        max_cluster = 0;
        number_of_clust = 0;
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

      void setMax_Cl(int v){ max_cluster = v; }
      int getMax_Cl(){ return max_cluster; }
      void setNum_of_Cl(int v){ number_of_clust = v; }
      int getNum_of_Cl(){ return number_of_clust; }
      /*
       * It is efficient only when communications are optimized.
       */
      RealND computeRaw() {

        opp_map.clear();
        pairs.clear();

        shared_ptr< storage::Storage > stor = getSystem()->storage;
        shared_ptr< mpi::communicator > cmm = getSystem()->comm;
        int this_node = cmm -> rank();

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
          if( !stor->lookupRealParticle(id) ){
            OrderParticleProps &op = (*opm).second;
            if( !(op.getNumNN()==0) ) sendGhostInfo.push_back( (*opm).second );
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
        //for(CellListIterator cit(cells_loc); !cit.isDone(); ++cit) {
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
          if( stor->lookupGhostParticle(id) ){
            OrderParticleProps &op = (*opm).second;
            if( !(op.getNumNN()==0) ) sendGhostInfo.push_back( (*opm).second );
          }
        }

        maxSize, vecSize  = sendGhostInfo.size();
        mpi::all_reduce( *cmm, vecSize, maxSize, mpi::maximum<int>() );
        while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );
        totID.clear();
        boost::mpi::all_gather( *getSystem()->comm, &sendGhostInfo[0], maxSize, totID);

        for(vector<OrderParticleProps>::iterator it = totID.begin(); it!=totID.end(); ++it){
          OrderParticleProps &gop = *it;
          if( gop.getPID()!=-1 && stor->lookupRealParticle( gop.getPID() ) ){
          //if( gop.getPID()!=-1 && stor->lookupLocalParticle( gop.getPID() ) ){
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

        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps *opp_i = &(opp_map.find( pid ))->second;
          if( opp_i->getD() >= d_min && opp_i->getD() <= d_max ){
            opp_i->setSolid( true );
          }
        }

        if(incl_surface){
          // -----------------------------------------------------------------------
          // send ghost info
          shared_ptr< mpi::communicator > cmm = getSystem()->comm;
          vector <OrderParticleProps> sendGhostInfo;
          for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
            Particle& p = *cit;
            int pid = p.id();
            OrderParticleProps *opp_i = &(opp_map.find( pid ))->second;
            sendGhostInfo.push_back( *opp_i );
          }

          int maxSize, vecSize  = sendGhostInfo.size();
          mpi::all_reduce( *cmm, vecSize, maxSize, mpi::maximum<int>() );
          while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );
          vector< OrderParticleProps > totID;
          mpi::all_gather( *cmm, &sendGhostInfo[0], maxSize, totID);

          for(vector<OrderParticleProps>::iterator it = totID.begin(); it!=totID.end(); ++it){
            OrderParticleProps &gop = *it;

            if( gop.getPID()!=-1 ){
              if( gop.getSolid() && !gop.getSurface() ){
                int nnn = gop.getNumNN();
                for(int i=0; i<nnn; i++){

                    if( getSystem()->storage->lookupRealParticle( gop.getNN(i) ) ){
                      OrderParticleProps *opp_i_surf = &(opp_map.find( gop.getNN(i) ))->second;
                      if( !opp_i_surf->getSolid() ){
                        opp_i_surf->setSurface(true);
                      }
                    }
                }
              }
            }

          }
          //--------------------------------------------------------------------------
        }

      }










      void cluster_analysis(){
        wrong_label_local0.clear();
        outside_id.clear();
        wrong_label_local.clear();

        local_cluster_search();

        /*
        CellList cells = getSystem()->storage->getRealCells();
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps &opp_i = (opp_map.find( pid ))->second;
          if ( opp_i.getLabel() > -1){
            cout<< "cpu:  " << getSystem()->comm->rank() << "  pid: "<< pid<< "  label:"<< opp_i.getLabel() << endl;
          }
        }
        */

        global_cluster_analysis();

        biggest_cluster();
      }

      void global_cluster_analysis(){

        // ***************** communicate outside_id
        int locSize = outside_id.size(); int maxSize=0;
        mpi::all_reduce(*getSystem()->comm, locSize, maxSize, mpi::maximum<int>());

        communicate_label cccc; cccc.pid = -1; cccc.label = -1; cccc.cpu = -1;
        while(outside_id.size()<maxSize) outside_id.push_back( cccc );

        int numProc = getSystem()->comm->size();
        vector<communicate_label> tot_label_info = vector<communicate_label>( numProc * maxSize, cccc );

        mpi::all_gather( *getSystem()->comm, &outside_id[0], maxSize, tot_label_info);
        // ***************** communicate outside_id. END

        for(vector<communicate_label>::iterator it = tot_label_info.begin(); it!=tot_label_info.end(); ++it){
          communicate_label &cur_lab = *it;

          if( cur_lab.cpu != -1 && getSystem()->storage->lookupRealParticle( cur_lab.pid ) ){
              OrderParticleProps &opp = (opp_map.find( cur_lab.pid ))->second;
              int labll = opp.getLabel();

              if ( (opp.getSolid() || opp.getSurface() ) ){
                int min_lab=0, max_lab=0;
                if ( labll == -1 && opp.getSurface()){
                  opp.setLabel( cur_lab.label );
                }
                else if(labll != -1 && opp.getSurface() ){
                  // TODO add action if meet another label
                  //min_lab = min(labll, cur_lab.label);
                  //max_lab = max(labll, cur_lab.label);
                  //relab(min_lab, max_lab);

                  //relab11(min_lab, max_lab);
                }
                else if(labll != cur_lab.label){
                  // TODO add action if meet another label
                  min_lab = min(labll, cur_lab.label);
                  max_lab = max(labll, cur_lab.label);
                  //relab(min_lab, max_lab);

                  relab11(min_lab, max_lab);
                }
                /*
                else{
                  cout<<"!!!WRONG global_cluster_analysis!!!  cpu: "<< getSystem()->comm->rank() <<
                        " particle: "<< opp.getPID()<<
                        " is "<<  opp.getSolid() << "  " << opp.getSurface() <<
                        " and has label:  "<< opp.getLabel() << endl;
                }
                */
              }
          }
        }

        vector<int> wl, cl; // wrong and correct labels

        for(unordered_map<int, int>::iterator it = wrong_label_local.begin();
                it!=wrong_label_local.end(); ++it ){
          wl.push_back( (*it).first );
          cl.push_back( (*it).second );
        }

        // ******* collect all wrong labels to 0 cpu
        locSize = wl.size(); maxSize=0;
        mpi::all_reduce(*getSystem()->comm, locSize, maxSize, mpi::maximum<int>());
        while(wl.size()<maxSize){ wl.push_back( -1 ); cl.push_back( -1 ); }

        numProc = getSystem()->comm->size();
        vector<int> totWL = vector<int>( numProc * maxSize, -1 );
        vector<int> totCL = vector<int>( numProc * maxSize, -1 );
        mpi::gather( *getSystem()->comm, &wl[0], maxSize, totWL, 0);
        mpi::gather( *getSystem()->comm, &cl[0], maxSize, totCL, 0);

        wrong_label_local0.clear();
        wl.clear(); cl.clear();

        vector<int> wl0, cl0; // wrong and correct labels


        if( getSystem()->comm->rank()==0 ){
          for(int i=0; i<totWL.size(); i++){
            if( totWL[i]>-1 && totCL[i]>-1 ){
              relab22(totCL[i], totWL[i]);
            }
          }

          for(unordered_map<int, int>::iterator it = wrong_label_local0.begin();
                  it!=wrong_label_local0.end(); ++it ){
            wl0.push_back( (*it).first );
            cl0.push_back( (*it).second );
          }
        }

        int vec_size;
        if( getSystem()->comm->rank()==0 ){
          vec_size=wl0.size();
        }
        mpi::broadcast( *getSystem()->comm, vec_size, 0);
        if( getSystem()->comm->rank()!=0 ){
          wl0 = vector<int>(vec_size, -1);
          cl0 = vector<int>(vec_size, -1);
        }

        mpi::broadcast( *getSystem()->comm, &wl0[0], wl0.size(), 0);
        mpi::broadcast( *getSystem()->comm, &cl0[0], cl0.size(), 0);

        wrong_label_local.clear();
        for(int i=0; i<wl0.size(); i++){
          wrong_label_local.insert( make_pair(wl0[i], cl0[i]) );
        }

        relabel_particles();

      }

      unordered_map< int, int > wrong_label_local0;

      void relab22(int new_lab, int old_lab){
        if( wrong_label_local0.find( old_lab ) == wrong_label_local0.end() ){
          wrong_label_local0.insert( make_pair(old_lab, new_lab) );
        }
        else if( ( wrong_label_local0.find( old_lab ) )->second != new_lab ){
          int correct_old = ( wrong_label_local0.find( old_lab ) )->second;

          int mmin = min(correct_old, new_lab);
          int mmax = max(correct_old, new_lab);

          relab22(mmin, mmax);
        }
      }


      void local_cluster_search(){
        int num_proc = getSystem()->comm->size();
        int loc_proc = getSystem()->comm->rank();

        int i_label_local = 0; // local increment
        int cur_label = loc_proc + i_label_local * num_proc; // unique label for each cpu

        outside_id.clear();
        CellList cells = getSystem()->storage->getRealCells();
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps &opp_i = (opp_map.find( pid ))->second;
          if ( opp_i.getLabel() == -1 && opp_i.getSolid() && !(opp_i.getSurface()) ){
//          if ( opp_i.getLabel() == -1 && (opp_i.getSolid() || opp_i.getSurface()) ){
            cluster_walk(pid, cur_label, opp_i);

            i_label_local++;
            cur_label=loc_proc + i_label_local * num_proc;
          }
        }

        relabel_particles();
      }

// ***************************************************************************************
      struct communicate_label{
        int pid;
        int label;
        int cpu;

        template <typename Archive>
        void serialize(Archive& ar, const unsigned int version){
          ar & pid;
          ar & label;
          ar & cpu;
        }
      };

      vector< communicate_label > outside_id;

      unordered_map< int, int > wrong_label_local;

// ***************************************************************************************

      void cluster_walk(int pid, int cur_lab, OrderParticleProps &opp){
        opp.setLabel(cur_lab);
        int num_nn = opp.getNumNN();

        for( int j=0; j< num_nn; j++) {
          int pid_neib = opp.getNN(j);

          if( !getSystem()->storage->lookupRealParticle( pid_neib ) ){
            communicate_label cc;
            cc.pid = pid_neib;
            cc.label = cur_lab;
            cc.cpu = getSystem()->comm->rank();
            outside_id.push_back( cc );
          }
          else{
            OrderParticleProps &opp_neib = (opp_map.find( pid_neib ))->second;
            int lab_neib = opp_neib.getLabel();

            if ( lab_neib != cur_lab && pid_neib!=pid &&
                    opp_neib.getSolid() && !(opp_neib.getSurface()) ){
/*
            if ( lab_neib != cur_lab && pid_neib!=pid &&
                    (opp_neib.getSolid() || opp_neib.getSurface()) ){
 */

              int min_lab=0, max_lab=0;
              if ( lab_neib == -1){
                min_lab = cur_lab;
              }
              else{
                min_lab = min(lab_neib, cur_lab);
                max_lab = max(lab_neib, cur_lab);
                relab11(min_lab, max_lab);
              }

              cluster_walk(pid_neib, min_lab, opp_neib);
            }
            else if( opp_neib.getSurface() && lab_neib==-1){
              opp_neib.setLabel(cur_lab);
            }


          }

        }
      }

      void relab(int new_lab, int old_lab){
        for(boost::unordered_multimap<int, OrderParticleProps>::iterator opm = opp_map.begin(); opm!=opp_map.end(); ++opm){
          OrderParticleProps &opp = (*opm).second;
          if( opp.getLabel() == old_lab ){
            opp.setLabel(new_lab);
          }
        }
      }

      void relab11(int new_lab, int old_lab){
        if( wrong_label_local.find( old_lab ) == wrong_label_local.end() ){
          wrong_label_local.insert( make_pair(old_lab, new_lab) );
        }
        else if( ( wrong_label_local.find( old_lab ) )->second != new_lab ){
          int correct_old = ( wrong_label_local.find( old_lab ) )->second;

          int mmin = min(correct_old, new_lab);
          int mmax = max(correct_old, new_lab);

          relab11(mmin, mmax);
        }
      }

      int getCorrectLabel(int current_label){
        if( wrong_label_local.find( current_label ) == wrong_label_local.end() ){
          return current_label;
        }
        else{
          return getCorrectLabel( (wrong_label_local.find( current_label ))->second );
        }
      }

      void relabel_particles(){
        for(boost::unordered_multimap<int, OrderParticleProps>::iterator opm = opp_map.begin(); opm!=opp_map.end(); ++opm){
          OrderParticleProps &opp = (*opm).second;
          if( opp.getSolid() || opp.getSurface() ){
            if( opp.getLabel()<0 && opp.getSolid() ){
              cout<<"!!!relabel_particles!!!  cpu: "<< getSystem()->comm->rank() <<
                      " particle: "<< opp.getPID()<<
                      " is "<<  opp.getSolid() << "  " << opp.getSurface() <<
                      " and has label:  "<< opp.getLabel() << endl;
            }
            else if(opp.getLabel()>=0){
              opp.setLabel( getCorrectLabel( opp.getLabel() ) );
            }
          }
        }
      }


      void biggest_cluster(){
        vector< vector<int> > clusters;


        CellList cells = getSystem()->storage->getRealCells();
        vector <OrderParticleProps> sendGhostInfo;
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps opp_i = (opp_map.find( pid ))->second;
          sendGhostInfo.push_back( opp_i );
        }

        vector<int> sendSizes;
        int ss  = sendGhostInfo.size();
        mpi::all_gather( *getSystem()->comm, ss, sendSizes);
        int maxSize=0;
        for(vector<int>::iterator it = sendSizes.begin(); it!=sendSizes.end(); ++it){
          maxSize = max(maxSize, *it);
        }
        while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );

        int numProc = getSystem()->comm->size();
        vector<OrderParticleProps> totID = vector<OrderParticleProps>( numProc * maxSize, OrderParticleProps() );
        mpi::gather( *getSystem()->comm, &sendGhostInfo[0], maxSize, totID, 0);

        if( getSystem()->comm->rank()==0 ){

          for(vector<OrderParticleProps>::iterator it=totID.begin(); it!=totID.end(); ++it) {
            OrderParticleProps &gop = *it;
            if(gop.getPID()>=0){

              if( gop.getLabel()>=0 ){
                while( (int)clusters.size() <= (int)gop.getLabel() ) clusters.push_back( vector<int>(0, 0) );
                clusters[ gop.getLabel() ].push_back( gop.getPID() );
              }

            }
          }





        /*
        for(boost::unordered_multimap<int, OrderParticleProps>::iterator opm = opp_map.begin(); opm!=opp_map.end(); ++opm){
          OrderParticleProps &opp = (*opm).second;

          if( opp.getLabel()>=0 ){
            while( (int)clusters.size() <= (int)opp.getLabel() ){
              clusters.push_back( vector<int>(0, 0) );
            }
            clusters[ opp.getLabel() ].push_back( opp.getPID() );
          }
        }
        */

          int max_pos = -1;
          int max_size = 0;
          int empty = 0;
          for (vector< vector<int> >::iterator it = clusters.begin() ; it != clusters.end(); ++it){
            if( max_size != max(max_size, (int)( (*it).size() ) ) ){
              max_pos = it - clusters.begin();
            }
            max_size = max(max_size, (int)( (*it).size() ) );

            if( (int)( (*it).size() )==0 ){
              empty++;
              //cout<<"  cpu: "<< getSystem()->comm->rank() <<   " label: "<< (it - clusters.begin()) << "  is empty" << endl;
            }
          }

          /*
          cout<<"  cpu: "<< getSystem()->comm->rank() <<  " max pos: "<< max_pos << endl;
          int count=0;
          for (vector<int>::iterator it = clusters[max_pos].begin() ; it != clusters[max_pos].end(); ++it){
            cout<<count<<"  BIGGEST CL  cpu: "<< getSystem()->comm->rank() <<  " id: "<< *it << endl;
            count++;
          }
          */

          int num_clusters = clusters.size();
          setNum_of_Cl(num_clusters-empty);
          setMax_Cl(max_size);
        }

      }




      python::list compute() {
        python::list ret;

        RealND res = computeRaw();

        if( do_cl_an ){
          define_solid();

          cluster_analysis();
        }

        /*
        CellList cells = getSystem()->storage->getRealCells();
        vector <OrderParticleProps> sendGhostInfo;
        for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
          Particle& p = *cit;
          int pid = p.id();
          OrderParticleProps opp_i = (opp_map.find( pid ))->second;
          sendGhostInfo.push_back( opp_i );
        }

        vector<int> sendSizes;
        int ss  = sendGhostInfo.size();
        mpi::all_gather( *getSystem()->comm, ss, sendSizes);
        int maxSize=0;
        for(vector<int>::iterator it = sendSizes.begin(); it!=sendSizes.end(); ++it){
          maxSize = max(maxSize, *it);
        }
        while(sendGhostInfo.size()<maxSize) sendGhostInfo.push_back( OrderParticleProps() );

        int numProc = getSystem()->comm->size();
        vector<OrderParticleProps> totID = vector<OrderParticleProps>( numProc * maxSize, OrderParticleProps() );
        mpi::gather( *getSystem()->comm, &sendGhostInfo[0], maxSize, totID, 0);

        int myN = getSystem()->storage->getNRealParticles();
        int num_part=0;
        mpi::all_reduce(*getSystem()->comm, myN, num_part, std::plus<int>());
        */

        if( getSystem()->comm->rank()==0 ){

          /*
          for(vector<OrderParticleProps>::iterator it=totID.begin(); it!=totID.end(); ++it) {
            OrderParticleProps &gop = *it;
            if(gop.getPID()>=0){
              python::list nn;

              for(int i=0; i<gop.getNumNN(); i++){
                nn.append( gop.getNN(i) );
              }

              python::tuple pt = python::make_tuple( gop.getPID(), gop.getD(),
                      gop.getNumNN(), gop.getSolid(), gop.getSurface(), nn);
              ret.append( pt );
            }
          }

          ret.sort();

          cout << "\n   inside OP:  max_cl: "<< getMax_Cl() << "  num_cl: "<< getNum_of_Cl()<< endl;
          */

          ret.append(getMax_Cl());
          ret.append(getNum_of_Cl());
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
