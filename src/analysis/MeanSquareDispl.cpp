/*
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

#include "MeanSquareDispl.hpp"
//#include <algorithm> //for std::sort
using namespace std;
//using namespace espressopp;

namespace espressopp {
  namespace analysis {

    //using namespace iterator;
    
    /* 
     * calculates the mean square displacement of the particles/monomers in the COM of the whole system
     * 
     * calc <r^2> the output is the average mean sq. displacement over 3 directions.
     * !! Important!! For D calculation factor 1/6 is already taken into account.
     * !! all confs should contain the same number of particles
    */
    
    python::list MeanSquareDispl::compute() const{
      
      int M = getListSize(); //number of snapshots/configurations
      real* totZ; //will store the mean squared displacement
      totZ = new real[M];
      real* Z;
      Z = new real[M];

      python::list pyli;
      
      System& system = getSystemRef();
      
      //creating vector which stores particleIDs for each CPU
      vector<longint> localIDs; 
      for (map<size_t,int>::const_iterator itr=idToCpu.begin(); itr!=idToCpu.end(); ++itr) {
        size_t i = itr->first;
        int whichCPU = itr->second;
        if(system.comm->rank()==whichCPU){
          localIDs.push_back(i);
        }
      }
      
      // COM calculation
      vector<Real3D> centerOfMassList;
      for(int m=0; m<M; m++){
        Real3D posCOM = Real3D(0.0,0.0,0.0);
        real mass = 0.0;
        Real3D posCOM_sum = Real3D(0.0,0.0,0.0);
        real mass_sum = 0.0;

        for (vector<longint>::iterator itr=localIDs.begin(); itr!=localIDs.end(); ++itr) {
            size_t i = *itr;
            Real3D pos = getConf(m)->getCoordinates(i);
            posCOM += pos;
            mass += 1;          
        }

        boost::mpi::all_reduce(*mpiWorld, posCOM, posCOM_sum, std::plus<Real3D>());
        boost::mpi::all_reduce(*mpiWorld, mass, mass_sum, std::plus<real>());

        centerOfMassList.push_back( posCOM_sum / mass_sum );
      }
      
      // MSD calculation
      int perc=0, perc1=0;
      real denom = 100.0 / (real)M;
      for(int m=0; m<M; m++){
        
        totZ[m] = 0.0;
        Z[m] = 0.0;
        for(int n=0; n<M-m; n++){
          for (vector<longint>::iterator itr=localIDs.begin(); itr!=localIDs.end(); ++itr) {
            size_t i = *itr;
            
            Real3D pos1 = getConf(n + m)->getCoordinates(i); // - centerOfMassList[n+m];
            Real3D pos2 = getConf(n)->getCoordinates(i); //     - centerOfMassList[n];
            Real3D delta = pos2 - pos1;
            Z[m] += delta.sqr();
          }
        }
        if(print_progress && system.comm->rank()==0){
          perc = (int)(m*denom);
          if(perc%5==0){
            cout<<"calculation progress (mean square displacement): "<< perc << " %\r"<<flush;
          }
        }
      }
      if(system.comm->rank()==0)
        cout<<"calculation progress (mean square displacement): 100%"<<endl;
      //summation of results from different CPUs
      boost::mpi::all_reduce( *system.comm, Z, M, totZ, plus<real>() );
      
      for(int m=0; m<M; m++){
        totZ[m] /= (real)(M - m);
      }
      
      real inv_coef = 1.0 / (6.0 * num_of_part);
      
      for(int m=0; m<M; m++){
        totZ[m] *= inv_coef;
        pyli.append( totZ[m] );
      }
      
      delete [] Z;
      Z = NULL;
      delete [] totZ;
      totZ = NULL;
      
      return pyli;
    }
    
    /* 
     * calculates mean square displacement of monomers in COM of their chains
     *         
     * !! currently only works for particles numbered like 0, 1, 2,... !!
     * !! with each chain consisting particles with subsequent ids     !!
     * 
     * calc <r^2> the output is the average mean sq. displacement over 3 directions.
     * !! Important!! For D calculation factor 1/6 is already taken into account.
     * !! all confs should contain the same number of particles
    */
    python::list MeanSquareDispl::computeG2() const{
//      cout << "0 got here!\n";
//      int M = getListSize(); //number of snapshots/configurations
//      real* totZ; //will store the mean squared displacement
//      totZ = new real[M];
//      real* Z;
//      Z = new real[M];
        cout << "Sorry. Computation of G2 does not work correctly.";
        cout << "Update ESPResSO++ when it is repaired. \n";

      python::list pyli;
      
//      System& system = getSystemRef();
//      
//      //creating vector which stores particleIDs for each CPU
//      vector<longint> localIDs; //for each CPU this will store particle IDs of particles calculated by CPU
//      for (map<size_t,int>::const_iterator itr=idToCpu.begin(); itr!=idToCpu.end(); ++itr) {
//        size_t i = itr->first; //particle ID
//        int whichCPU = itr->second; //CPU number
//        printf("id %u  CPU %i \n", i, whichCPU);
//        if(system.comm->rank()==whichCPU){
//          localIDs.push_back(i);
//        }
//      }
//      sort(localIDs.begin(), localIDs.end()); //sorts entries from low to high
//      //should not be necessary as long as the above iterator 
//      //iterates in ascending order according to the keys 
//           
//      // COM calculation      
//      Real3D posCOM = Real3D(0.0,0.0,0.0);
//      real mass = 0.0;
//      int count = 0; //counts number of particles of one chain
//      vector< vector<Real3D> > local_chainCOMlist; //will store COM of conf n and chain cid as chainCOMlist[n][cid]
//      for(int m=0; m<M; m++){  
//        vector<Real3D> innerList; //will store the local chains' COM of one conf/snapshot
//
//        //loop over local particles
//        for(int entry = 0; entry < localIDs.size(); entry++){           
//            longint i = localIDs[entry]; //pid
//            Real3D pos = getConf(m)->getCoordinates(i);
//            posCOM += pos;
//            mass += 1;
//            count += 1;
//            // this is the right if request. remember that particle 0 also has a mass of 1
//            if (count == chainlength){
//                innerList.push_back( posCOM / mass);
//                posCOM = Real3D(0.0,0.0,0.0);
//                mass = 0;
//                count = 0;
//            }            
//        } //now innerList contains the local chains' COMs of snapshot m
//        local_chainCOMlist.push_back(innerList);
//      }
//      //now local_chainCOMlist contains the local chains' COMs of each snapshot      
//          
//      // MSD calculation
//      int perc=0, perc1=0;
//      real denom = 100.0 / (real)M;
//      for(int m=0; m<M; m++){
//        totZ[m] = 0.0;
//        Z[m] = 0.0;
//        for(int n=0; n<M-m; n++){
//          int part_count = 0;
//          int local_cid = 0; //local chainID. each CPU starts with chain local_cid = 0, so it is not a global id
//          //loop over local particles
//          for(int entry = 0; entry < localIDs.size(); entry++){           
//            longint i = localIDs[entry]; //pid   
//            if(part_count == chainlength){
//                ++local_cid;
//                part_count = 0;
//            } 
//            //cout << "n, m, i, local_cid, part_count   " 
//            //        << n << "\t" << m << "\t"<< i << "\t" << local_cid << "\t" << part_count << "\n";
//            Real3D pos1 = getConf(n + m)->getCoordinates(i) - local_chainCOMlist[n+m][local_cid];
//            Real3D pos2 = getConf(n)->getCoordinates(i)     - local_chainCOMlist[n][local_cid];
//            Real3D delta = pos2 - pos1;
//            Z[m] += delta.sqr();
//            part_count++;
//          }
//        }
//        if(print_progress && system.comm->rank()==0){
//          perc = (int)(m*denom);
//          if(perc%5==0){
//            cout<<"calculation progress (mean square displacement): "<< perc << " %\r"<<flush;
//          }
//        }
//      }
//     
//      if(system.comm->rank()==0)
//        cout<<"calculation progress (mean square displacement): 100%"<<endl;
//      //summation of results from different CPUs
//      boost::mpi::all_reduce( *system.comm, Z, M, totZ, plus<real>() );
//      
//      for(int m=0; m<M; m++){
//        totZ[m] /= (real)(M - m);
//      }
//      
//      real inv_coef = 1.0 / (6.0 * num_of_part);
//      
//      for(int m=0; m<M; m++){
//        totZ[m] *= inv_coef;
//        pyli.append( totZ[m] );
//      }
//      
//      delete [] Z;
//      Z = NULL;
//      delete [] totZ;
//      totZ = NULL;
      
      return pyli;
    }

    python::list MeanSquareDispl::computeG3() const{
      cout << "0 got here!\n";
      int M = getListSize(); //number of snapshots/configurations
      real* totZ; //will store the mean squared displacement
      totZ = new real[M];
      real* Z;
      Z = new real[M];

      python::list pyli;
      
      System& system = getSystemRef();
      
      /** stores particleIDs of local particles calculated by CPU */
      vector<longint> localIDs;
      // filling each CPUs localIDs-vector
      for (map<size_t,int>::const_iterator itr=idToCpu.begin(); itr!=idToCpu.end(); ++itr) {
        size_t i = itr->first; //particle ID
        int whichCPU = itr->second; //CPU number
        printf("id %lu  CPU %i \n", i, whichCPU);
        if(system.comm->rank()==whichCPU){
          localIDs.push_back(i);
        }
      }
      sort(localIDs.begin(), localIDs.end()); //sorts entries from low to high
      //should not be necessary as long as the above iterator 
      //iterates in ascending order according to the keys 
      
      /** total number of chains */
      int num_chains = num_of_part / chainlength;
      cout << "num_of_part " << num_of_part << "\n";
      cout << "chainlenght " << chainlength << "\n";
      cout << "localID.size " << localIDs.size() << "\n";
       
      // COM calculation *******************************************************      
      // explanation of 'local_chainCOMlist' and 'innerList':
      // innerList:
      // each CPU calculates the COM of certain chains: its local chains   
      // innerlist contains these COMs (but only for one snapshot)
      // local_chainCOMlist:   
      // contains the inner lists of each snapshot
      //          
      // So:
      // innerList = [(),(),(),...]
      // where '()' are the COMs of the local chains (local = on that cpu)
      // local_chainCOMlist = [ [], [], [], ...]
      // where '[]' are the inner lists of the snapshots
      //
      // local_chainCOMlist[m=0][cid=2] is the COM of local_chain 2 of snapshot 0       
      
      
      /** stores COM of conf m and chain cid as local_chainCOMlist[m][cid] */
      vector< vector<Real3D> > local_chainCOMlist; 
      /** stores COM of the whole system for each snapshot, 
       *  called 'centerOfMassList in compute() */
      vector<Real3D> sysCOMlist; 
      
      for(int m=0; m<M; m++){        
        Real3D posCOM = Real3D(0.0,0.0,0.0); /**< position of COM of chain */
        real mass = 0.0; /**< mass of chain */  
        int part_count = 0; /**< counts number of particles of one chain */ 
        int chain_count = 0; /**< counts chains of one CPU */
        
        //testout
        //printf("snapshot %i... \n", m);
        
        Real3D sysCOM = Real3D(0.0,0.0,0.0); /**< COM of particles on one CPU */
                vector<Real3D> innerList; /**< local chains' COM of one conf/snapshot */

        //loop over local particles
        for(int entry = 0; entry < localIDs.size(); entry++){           
            longint i = localIDs[entry]; /**< pid */        
            Real3D pos = getConf(m)->getCoordinates(i);
            //testout
            //printf("pid %i :", i); printReal3D(pos); printf("\n");            
            posCOM += pos;
            mass += 1;
            part_count += 1;
            // this is the right if request. remember that particle 0 also has a mass of 1
            if (part_count == chainlength){
                Real3D chainCOM = posCOM / mass;
                innerList.push_back( chainCOM );
                sysCOM += chainCOM;             
                
                //testout
                //printf("COM "); printReal3D( posCOM / mass ); printf("\n");
                chain_count += 1;
                part_count = 0;
                posCOM = Real3D(0.0,0.0,0.0);
                mass = 0;                
            }            
        } //now innerList contains the local chains' COMs of snapshot m
        //sysCOM contains the sum of the local chains' COMs (sum of innerlist)
        local_chainCOMlist.push_back(innerList);
        
        Real3D sysCOM_sum = Real3D(0.0,0.0,0.0); /**< COM of system */
        int chain_sum = 0; //stores the total number of chains
        //FM do I have the communicator??
        boost::mpi::all_reduce(*mpiWorld, sysCOM, sysCOM_sum, std::plus<Real3D>());
        boost::mpi::all_reduce(*mpiWorld, chain_count, chain_sum, std::plus<real>());
        
        sysCOMlist.push_back(sysCOM_sum / chain_sum );
        //check if chain_count matches the number of chains
        if (chain_sum != num_chains){
            printf("ERROR: snapshot %i calculated %i chains.\t", m, chain_sum);
            printf("number of chains is: %i\n", num_chains);
            printf("calculation continues... \n");
            chain_count = 0;
        }
        else {chain_count = 0;}
        
      }
      //now local_chainCOMlist contains the local chains' COMs of each snapshot
      //sysCOMlist contains the systems COM for each snapshot
       
      //FM continue writing g3 code here... TEST first!!!
      
      //output for testing
      for(int m = 0; m < M; m++){
        //local_chainCOMlist[m].size() can be different on each CPU 
        //but should be the same for each snapshot,
        //so local_chainCOMlist[m].size() == local_chainCOMlist[n].size]
        for(int i = 0; i < local_chainCOMlist[m].size(); i++){      
          printf("COM snapshot %i, chain %i : ", m, i); 
          printReal3D(local_chainCOMlist[m][i]); printf("\n");
        }
      }
          
      // MSD calculation *******************************************************
      int perc=0; //, perc1=0;
      real denom = 100.0 / (real)M;
      for(int m=0; m<M; m++){
        totZ[m] = 0.0;
        Z[m] = 0.0;
        for(int n=0; n<M-m; n++){          
          
          //loop over local chains
          int local_chains = local_chainCOMlist[0].size(); //looked up only for snapshot 0
          //int local_chains = local_chainCOMlist[n].size(); //same value as above, but looked up new for every snapshot         
          for(int cid = 0; cid < local_chains; cid++){           
            //cout << "n, m, i, local_cid, part_count   " 
            //        << n << "\t" << m << "\t"<< i << "\t" << local_cid << "\t" << part_count << "\n";
            Real3D pos1 = local_chainCOMlist[n+m][cid] - sysCOMlist[n+m];
            Real3D pos2 = local_chainCOMlist[n][cid] - sysCOMlist[n];
            Real3D delta = pos2 - pos1;
            Z[m] += delta.sqr();           
          }
        }
        if(print_progress && system.comm->rank()==0){
          perc = (int)(m*denom);
          if(perc%5==0){
            cout<<"calculation progress (mean square displacement): "<< perc << " %\r"<<flush;
          }
        }
      }
     
      if(system.comm->rank()==0)
        cout<<"calculation progress (mean square displacement): 100%"<<endl;
      //summation of results from different CPUs
      boost::mpi::all_reduce( *system.comm, Z, M, totZ, plus<real>() );      
      
      real inv_coef = 1.0 / (6.0 * num_chains); 
      
      //averaging over the chains and the number of snapshot-pairs (n)
      //and creating the python list
      //snapshot pairs: (snap[0], snap[m]), (snap[1], snap[1+m]), ... 
      //all correspond to a difference m between snaps. there are n=M-m of each
      for(int m=0; m<M; m++){      
        totZ[m] = inv_coef * totZ[m] / (real)(M - m);        
        pyli.append( totZ[m] );
      }
      
      delete [] Z;
      Z = NULL;
      delete [] totZ;
      totZ = NULL;
      
      return pyli;
    }    

    
    
    // Python wrapping
    void MeanSquareDispl::registerPython() {
      using namespace espressopp::python;

      class_<MeanSquareDispl, bases<ConfigsParticleDecomp> >
      ("analysis_MeanSquareDispl", init< shared_ptr< System > >() )
      .def(init< shared_ptr< System >, int>() )
      .def("computeG2", &MeanSquareDispl::computeG2)
      .def("computeG3", &MeanSquareDispl::computeG3)
      .add_property("print_progress", &MeanSquareDispl::getPrint_progress, &MeanSquareDispl::setPrint_progress)
      ;
    }
  }
}
