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

#include "python.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "ConstMuMD.hpp"
#include "bc/BC.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "math.h"

#include <boost/serialization/map.hpp>
#if SYSTEM == 0
#include <boost/math/special_functions/spherical_harmonic.hpp>
#endif

namespace espressopp {

  using namespace espressopp;
  using namespace iterator;
  using namespace std;

  namespace integrator {

    ConstMuMD::ConstMuMD(shared_ptr<System> system): Extension(system){

      PI                =   3.1415926535897932384626433832795029;
      SIZE3D            =   5E4;
      MAXPARTNUM        =   5E5;
      RB_WIDTH          =   0.2;

      modus             =    -1;
      sizeTR            =  10.0;
      sizeCR            =  10.0;
      sizeFR            =   5.0;
      sizeCheckR        =  10.0;
      sizeResL          =  73.0;
      sizeResR          =  73.0;
      LSizeRes          = std::vector<real>(5,-1.0);
      RSizeRes          = std::vector<real>(5,-1.0);
      dr                = Real3D(1.0,1.0,1.0);
      nthresh           =  1000;
      cutoff            =   1.6;
      shapeParam        =   0.2;
      fConst0           = 500.0;
      fConst1           =  50.0;
      //expDens0          =   0.5;
      //expDens1          =   0.5;
      expMolF0          =   0.5;
      expMolF1          =   1.0-expMolF0;
      substrate         = Real3D(100.0,98.0,102.0);
      Ntot              = 10000;
      ratio             =   0.5;
      center_CR_L       =  83.0;
      center_CR_R       = 117.0;
      center_FR_L       =  75.5;
      center_FR_R       = 124.5;
      reservoir_L       =  73.0;
      reservoir_R       = 127.0;
      numDensInXtal     =   0.6;
      numDensInLiq      = Real3D(0.5,0.5,0.5);
      numDensInCR0      = Real3D(0.5,0.5,0.5);
      numDensInCR1      = Real3D(0.5,0.5,0.5);
      numDensInCRTot    = Real3D(0.5,0.5,0.5);
      numDensInRes      = Real3D(0.5,0.5,0.5);
      molFracInCR       = Real3D(0.8,0.8,0.8);

      type              = Extension::ConstMu;

    }

    void ConstMuMD::disconnect(){
      _initialize.disconnect();
      _perform.disconnect();
    }

    void ConstMuMD::connect(){
      // connection to initialization
      _initialize       = integrator->runInit.connect( boost::bind(&ConstMuMD::initialize, this));

      // connection to initialization
      //_applyForce     = integrator->aftInitF.connect( boost::bind(&ConstMuMD::applyForce, this));
      _perform          = integrator->aftInitF.connect( boost::bind(&ConstMuMD::perform, this));
    }

    // set and get performance modus
    void ConstMuMD::setModus(int _modus){
      modus            = _modus;
    }
    int ConstMuMD::getModus(){
      return modus;
    }

    // set and get size of transition region
    void ConstMuMD::setSizeTR(real _sizeTR){
      sizeTR            = _sizeTR;
    }
    real ConstMuMD::getSizeTR(){
      return sizeTR;
    }

    // set and get size of control region
    void ConstMuMD::setSizeCR(real _sizeCR){
      sizeCR            = _sizeCR;
    }
    real ConstMuMD::getSizeCR(){
      return sizeCR;
    }

    // set and get size of force region
    void ConstMuMD::setSizeFR(real _sizeFR){
      sizeFR            = _sizeFR;
    }
    real ConstMuMD::getSizeFR(){
      return sizeFR;
    }

    // set and get size of adapt. checker region
    void ConstMuMD::setSizeCheckR(real _sizeCheckR){
      sizeCheckR        = _sizeCheckR;
    }
    real ConstMuMD::getSizeCheckR(){
      return sizeCheckR;
    }

    // set and get sizes of reservoirs
    void ConstMuMD::setSizeRes(python::list& _sizeRes){
      sizeResL          = boost::python::extract<real>(_sizeRes[0]);
      sizeResR          = boost::python::extract<real>(_sizeRes[1]);
    }
    python::list ConstMuMD::getSizeRes(){
      python::list pyli;
      pyli.append( sizeResL );
      pyli.append( sizeResR );

      return pyli;
    }

    // set and get 3-dim. bin sizes
    void ConstMuMD::setDR(python::list& _dr){
      real drX           = boost::python::extract<real>(_dr[0]);
      real drY           = boost::python::extract<real>(_dr[1]);
      real drZ           = boost::python::extract<real>(_dr[2]);
      dr                 = Real3D(drX,drY,drZ);
      dV                 = drX*drY*drZ;
    }
    python::list ConstMuMD::getDR(){
      python::list pyli;
      pyli.append( dr[0] );
      pyli.append( dr[1] );
      pyli.append( dr[2] );

      return pyli;
    }

    // set and get number of bins in plane
    void ConstMuMD::setBinsInPlane(python::list& _binsInPlane){
      binsInPlaneY      = boost::python::extract<int>(_binsInPlane[0]);
      binsInPlaneZ      = boost::python::extract<int>(_binsInPlane[1]);
      sumBins           = binsInPlaneY*binsInPlaneZ;
    }
    python::list ConstMuMD::getBinsInPlane(){
      python::list pyli;
      pyli.append( binsInPlaneY );
      pyli.append( binsInPlaneZ );

      return pyli;
    }

    // set and get number of particles in xtal substrate
    void ConstMuMD::setNThresh(int _nthresh){
      nthresh           = _nthresh;
    }
    int ConstMuMD::getNThresh(){
      return nthresh;
    }

    // set and get cutoff radius
    void ConstMuMD::setCutoff(real _cutoff){
      cutoff            = _cutoff;
    }
    real ConstMuMD::getCutoff(){
      return cutoff;
    }

    // set and get shape parameter of membrane force
    void ConstMuMD::setShapeParam(real _shapeParam){
      shapeParam        = _shapeParam;
    }
    real ConstMuMD::getShapeParam(){
      return shapeParam;
    }

    // set and get prefactor of membrane force
    void ConstMuMD::setForceConst(python::list& _fConst){
      fConst0           = boost::python::extract<real>(_fConst[0]);
      fConst1           = boost::python::extract<real>(_fConst[1]);
    }
    python::list ConstMuMD::getForceConst(){
      python::list pyli;
      pyli.append( fConst0 );
      pyli.append( fConst1 );

      return pyli;
    }

//     // set and get expected number densities
//     void ConstMuMD::setExpDens(python::list& _expDens){
//       expDens0          = boost::python::extract<real>(_expDens[0]);
//       expDens1          = boost::python::extract<real>(_expDens[1]);
//     }
//     python::list ConstMuMD::getExpDens(){
//       python::list pyli;
//       pyli.append( expDens0 );
//       pyli.append( expDens1 );
//
//       return pyli;
//     }

    // set and get expected mole fraction solute
    void ConstMuMD::setExpMolF(real _expMolF){
      expMolF0          = _expMolF;
      expMolF1          = 1.0-_expMolF;
    }
    python::list ConstMuMD::getExpMolF(){
      python::list pyli;
      pyli.append( expMolF0 );
      pyli.append( expMolF1 );

      return pyli;
    }

    // set and get wall instance
    void ConstMuMD::setWallInstance(Real3D _substrate){
      substrate         = Real3D(_substrate[0],_substrate[1],_substrate[2]);
    }
    Real3D ConstMuMD::getWallInstance(){
      return substrate;
    }

    // set and get structure of left surface
    //void ConstMuMD::setWallL(python::list& _wallL){
      //for(int i=0; i<sumBins; i++){
	//wallInstanceL[i]= boost::python::extract<real>(_wallL[i]);
      //}
    //}
    //python::list ConstMuMD::getWallL(){
      //python::list pyli;
      //for(short i=0; i<sumBins; i++){
	//pyli.append( wallInstanceL[i] );
      //}
      //return pyli;
    //}

    // set and get structure of right surface
    //void ConstMuMD::setWallR(python::list& _wallR){
      //for(int i=0; i<sumBins; i++){
	//wallInstanceR[i]= boost::python::extract<real>(_wallR[i]);
      //}
    //}
    //python::list ConstMuMD::getWallR(){
      //python::list pyli;
      //for(short i=0; i<sumBins; i++){
	//pyli.append( wallInstanceR[i] );
      //}
      //return pyli;
    //}

    // set and get total number of particles
    void ConstMuMD::setNtot(int _Ntot){
      Ntot              = _Ntot;
    }
    int ConstMuMD::getNtot(){
      return Ntot;
    }

    // set and get prefactor for threshold checker
    void ConstMuMD::setRatio(real _ratio){
      ratio            = _ratio;
    }
    real ConstMuMD::getRatio(){
      return ratio;
    }

    // initialization
    void ConstMuMD::initialize(){
      int minbins                  = 100*sumBins;
      if( SIZE3D<minbins )
	SIZE3D                    *= (int(minbins/SIZE3D)+1);

      if( MAXPARTNUM<Ntot )
	MAXPARTNUM                *= (int(Ntot/MAXPARTNUM)+1);

      System& system               = getSystemRef();
      CellList realCells           = system.storage->getRealCells();
      Real3D Li                    = system.bc->getBoxL();
      Real3D Li_half               = 0.5*Li;
      Int3D dhBins                 = Int3D(int((Li[0]/dr[0])+0.5),binsInPlaneY,binsInPlaneZ);

      //SIZE3D                       = dhBins[0]*dhBins[1]*dhBins[2];

//       int  neighbours[MAXPARTNUM];
//       int  neighbours_Tot[MAXPARTNUM];
//       real q6_Re[7*MAXPARTNUM];
//       real q6_Im[7*MAXPARTNUM];
//       real q6_Re_Tot[7*MAXPARTNUM];
//       real q6_Im_Tot[7*MAXPARTNUM];
      neighbours                   = std::vector<int>(MAXPARTNUM,0);
      neighbours_Tot               = std::vector<int>(MAXPARTNUM,0);
      q6_Re                        = std::vector<real>(7*MAXPARTNUM,0.0);
      q6_Im                        = std::vector<real>(7*MAXPARTNUM,0.0);
      q6_Re_Tot                    = std::vector<real>(7*MAXPARTNUM,0.0);
      q6_Im_Tot                    = std::vector<real>(7*MAXPARTNUM,0.0);
//       for(int i=0; i<MAXPARTNUM; i++){
// 	neighbours[i]              = 0;
// 	neighbours_Tot[i]          = 0;
// 	for(int j=0; j<=6; j++){
// 	  q6_Re[7*i+j]             = 0.0;
// 	  q6_Im[7*i+j]             = 0.0;
// 	  q6_Re_Tot[7*i+j]         = 0.0;
// 	  q6_Im_Tot[7*i+j]         = 0.0;
// 	}
//       }

      partsInCell                  = std::vector<int>(SIZE3D,0);
      totpartsInCell               = std::vector<int>(SIZE3D,0);
      neighboursInCell             = std::vector<int>(SIZE3D,0);
      totneighboursInCell          = std::vector<int>(SIZE3D,0);
      q6InCell                     = std::vector<real>(SIZE3D,0.0);
      totq6InCell                  = std::vector<real>(SIZE3D,0.0);

      for(CellListAllPairsIterator it(realCells); it.isValid(); ++it){
        if( (it->first->type()==0) && (it->second->type()==0) ){
	  Real3D distVector        = it->first->position()-it->second->position();
          // minimize the distance in simulation box
          for(short ii=0; ii<3; ii++){
            if( distVector[ii]<-Li_half[ii] )
	      distVector[ii]      += Li[ii];
            if( distVector[ii]>Li_half[ii] )
	      distVector[ii]      -= Li[ii];
          }
          real distAbs             = distVector.abs();
	  if( distAbs<cutoff ){
	    real theta             = 0.0;
	    real phi               = 0.0;
	    getAngles(distVector,distAbs,theta,phi);
	    #if SYSTEM_SWITCH == 1
	    sphHarm((real *)&q6_Re,(real *)&q6_Im,7*(it->first->id()),7*(it->second->id()),theta,phi);
	    #endif
	    #if SYSTEM_SWITCH == 0
	    for(int m=0; m<=6; m++){
	      q6_Re[7*(it->first->id())+m] += boost::math::spherical_harmonic_r(6,m,theta,phi);
	      q6_Im[7*(it->first->id())+m] += boost::math::spherical_harmonic_i(6,m,theta,phi);
	      q6_Re[7*(it->second->id())+m]+= boost::math::spherical_harmonic_r(6,m,theta,phi);
	      q6_Im[7*(it->second->id())+m]+= boost::math::spherical_harmonic_i(6,m,theta,phi);
	    }
	    #endif
	    neighbours[it->first->id()]++;
	    neighbours[it->second->id()]++;
	  }
        }
      }
      boost::mpi::all_reduce(*mpiWorld, (int*)&neighbours[0], MAXPARTNUM, (int*)&neighbours_Tot[0], std::plus<int>());
      boost::mpi::all_reduce(*mpiWorld, (real*)&q6_Re[0], 7*MAXPARTNUM, (real*)&q6_Re_Tot[0], std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, (real*)&q6_Im[0], 7*MAXPARTNUM, (real*)&q6_Im_Tot[0], std::plus<real>());

      real borderL                 = substrate[1]-sizeTR-sizeCR-sizeFR;
      real borderR                 = substrate[2]+sizeTR+sizeCR+sizeFR;
//      real partsInReg              = 0.0;
//      real totpartsInReg           = 0.0;
      real nullInReg               = 0.0;
      real totnullInReg            = 0.0;
      real neighboursInReg         = 0.0;
      real totneighboursInReg      = 0.0;
      real q6InReg                 = 0.0;
      real totq6InReg              = 0.0;
      real neighboursInXtal        = 0.0;
      real totneighboursInXtal     = 0.0;
      real q6InXtal                = 0.0;
      real totq6InXtal             = 0.0;
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit){
	Real3D loc1                = cit->position();
	if( cit->type()==0 ){
	  real q6_bet              = 0.0;
	  if( neighbours_Tot[cit->id()] > 0 ){
	    for(int m=0; m<=6; m++){
	      q6_Re_Tot[7*cit->id()+m]/= (real)neighbours_Tot[cit->id()];
	      q6_Im_Tot[7*cit->id()+m]/= (real)neighbours_Tot[cit->id()];
	    }
	  }
	  q6_bet                  += (q6_Re_Tot[7*cit->id()]*q6_Re_Tot[7*cit->id()]+q6_Im_Tot[7*cit->id()]*q6_Im_Tot[7*cit->id()]);
	  for(int f=1; f<=6; f++){
	    q6_bet                += (2.0*(q6_Re_Tot[7*cit->id()+f]*q6_Re_Tot[7*cit->id()+f]+q6_Im_Tot[7*cit->id()+f]*q6_Im_Tot[7*cit->id()+f]));
	  }
	  q6_bet                  *= 4.0*PI/13.0;
	  q6_bet                   = sqrt(q6_bet);
	  if( (loc1[0]<borderL) || (loc1[0]>borderR) ){
	    q6InReg               += q6_bet;
	    neighboursInReg       += (real)neighbours_Tot[cit->id()];
	    nullInReg             += 1.0;
	  }
	  if( cit->id()<=nthresh ){
	    q6InXtal              += q6_bet;
	    neighboursInXtal      += (real)neighbours_Tot[cit->id()];
	  }
	}
// 	if( (loc1[0]<borderL) || (loc1[0]>borderR) ){
// 	  partsInReg              += 1.0;
// 	}
      }
//      boost::mpi::all_reduce(*mpiWorld, partsInReg, totpartsInReg, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, nullInReg, totnullInReg, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, neighboursInReg, totneighboursInReg, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, q6InReg, totq6InReg, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, neighboursInXtal, totneighboursInXtal, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, q6InXtal, totq6InXtal, std::plus<real>());

//      expDens                      = totpartsInReg/((borderL+Li[0]-borderR)*Li[1]*Li[2]);
      checkerQ6                    = ratio*totq6InXtal/nthresh+(1.0-ratio)*totq6InReg/totnullInReg;
      checkerCN                    = ratio*totneighboursInXtal/nthresh+(1.0-ratio)*totneighboursInReg/totnullInReg;
      maxGrowth                    = 0.5*sizeCheckR;

    }

    void ConstMuMD::perform(){
      adaptRegions();
      if( modus>=0 ){
	assignParts();
	applyForce();
      }
    }

    // get densities for the different regions
    python::list ConstMuMD::getDensityData(){
      python::list pyli;
      pyli.append( numDensInXtal );
      pyli.append( numDensInLiq[2] );
      pyli.append( numDensInCRTot[2] );
      pyli.append( numDensInRes[2] );
      pyli.append( molFracInCR[2] );

      return pyli;
    }

    void ConstMuMD::adaptRegions(){

      // ############################################################################################
      // ## SET REGIONS                                                                            ##
      // ############################################################################################

      System& system               = getSystemRef();
      CellList realCells           = system.storage->getRealCells();
      Real3D Li                    = system.bc->getBoxL();
      Real3D Li_half               = 0.5*Li;
      Int3D dhBins                 = Int3D(int((Li[0]/dr[0])+0.5),binsInPlaneY,binsInPlaneZ);
      int binsToCheckL[2]          = {int((substrate[1]-(0.5*sizeCheckR))/dr[0]),int((substrate[1]+(0.5*sizeCheckR))/dr[0])};
      int binsToCheckR[2]          = {int((substrate[2]-(0.5*sizeCheckR))/dr[0]),int((substrate[2]+(0.5*sizeCheckR))/dr[0])};

//       int  neighbours[MAXPARTNUM];
//       int  neighbours_Tot[MAXPARTNUM];
//       real q6_Re[7*MAXPARTNUM];
//       real q6_Im[7*MAXPARTNUM];
//       real q6_Re_Tot[7*MAXPARTNUM];
//       real q6_Im_Tot[7*MAXPARTNUM];
      for(int i=0; i<MAXPARTNUM; i++){
	neighbours[i]              = 0;
	neighbours_Tot[i]          = 0;
	for(int j=0; j<=6; j++){
	  q6_Re[7*i+j]             = 0.0;
	  q6_Im[7*i+j]             = 0.0;
	  q6_Re_Tot[7*i+j]         = 0.0;
	  q6_Im_Tot[7*i+j]         = 0.0;
	}
      }

      for(CellListAllPairsIterator it(realCells); it.isValid(); ++it){
	int checkBin                 = int(it->first->position()[0]/dr[0]);
	if( ((checkBin>=binsToCheckL[0]) && (checkBin<=binsToCheckL[1])) || ((checkBin>=binsToCheckR[0]) && (checkBin<=binsToCheckR[1])) ){
	  if( (it->first->type()==0) && (it->second->type()==0) ){
	    Real3D distVector        = it->first->position()-it->second->position();
	    // minimize the distance in simulation box
	    for(short ii=0; ii<3; ii++){
	      if( distVector[ii]<-Li_half[ii] )
		distVector[ii]      += Li[ii];
	      if( distVector[ii]>Li_half[ii] )
		distVector[ii]      -= Li[ii];
	    }
	    real distAbs             = distVector.abs();
	    if( distAbs<cutoff ){
	      real theta             = 0.0;
	      real phi               = 0.0;
	      getAngles(distVector,distAbs,theta,phi);
	      #if SYSTEM_SWITCH == 1
	      sphHarm((real *)&q6_Re,(real *)&q6_Im,7*(it->first->id()),7*(it->second->id()),theta,phi);
	      #endif
	      #if SYSTEM_SWITCH == 0
	      for(int m=0; m<=6; m++){
		q6_Re[7*(it->first->id())+m] += boost::math::spherical_harmonic_r(6,m,theta,phi);
		q6_Im[7*(it->first->id())+m] += boost::math::spherical_harmonic_i(6,m,theta,phi);
		q6_Re[7*(it->second->id())+m]+= boost::math::spherical_harmonic_r(6,m,theta,phi);
		q6_Im[7*(it->second->id())+m]+= boost::math::spherical_harmonic_i(6,m,theta,phi);
	      }
	      #endif
	      neighbours[it->first->id()]++;
	      neighbours[it->second->id()]++;
	    }
	  }
	}
      }
      boost::mpi::all_reduce(*mpiWorld, (int*)&neighbours[0], MAXPARTNUM, (int*)&neighbours_Tot[0], std::plus<int>());
      boost::mpi::all_reduce(*mpiWorld, (real*)&q6_Re[0], 7*MAXPARTNUM, (real*)&q6_Re_Tot[0], std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, (real*)&q6_Im[0], 7*MAXPARTNUM, (real*)&q6_Im_Tot[0], std::plus<real>());


//       int  partsInCell[SIZE3D];
//       int  totpartsInCell[SIZE3D];
//       int  neighboursInCell[SIZE3D];
//       int  totneighboursInCell[SIZE3D];
//       real q6InCell[SIZE3D];
//       real totq6InCell[SIZE3D];
      for(int i=0; i<SIZE3D; i++){
	partsInCell[i]             = 0;
	totpartsInCell[i]          = 0;
	neighboursInCell[i]        = 0;
	totneighboursInCell[i]     = 0;
	q6InCell[i]                = 0.0;
	totq6InCell[i]             = 0.0;
      }

      real localWCenterX           = 0.0;
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit){
	if( cit->id()<=nthresh ){
	  localWCenterX           += cit->position()[0];
	}
	int checkBin               = int(cit->position()[0]/dr[0]);
	if( ((checkBin>=binsToCheckL[0]) && (checkBin<=binsToCheckL[1])) || ((checkBin>=binsToCheckR[0]) && (checkBin<=binsToCheckR[1])) ){
	  if( cit->type()==0 ){
	    Int3D bin                  = Int3D(0,0,0);
	    Real3D loc1                = cit->position();
	    for(short i=0; i<3; i++){
	      if( loc1[i]<0.0 ){
		bin[i]                 = int((loc1[i]+Li[i])/dr[i]);
	      } else if( loc1[i]>Li[i] ) {
		bin[i]                 = int((loc1[i]-Li[i])/dr[i]);
	      } else {
		bin[i]                 = int(loc1[i]/dr[i]);
	      }
	    }

	    real q6_bet                = 0.0;
	    if( neighbours_Tot[cit->id()] > 0 ){
	      for(int m=0; m<=6; m++){
		q6_Re_Tot[7*cit->id()+m]/= (real)neighbours_Tot[cit->id()];
		q6_Im_Tot[7*cit->id()+m]/= (real)neighbours_Tot[cit->id()];
	      }
	    }
	    q6_bet                    += (q6_Re_Tot[7*cit->id()]*q6_Re_Tot[7*cit->id()]+q6_Im_Tot[7*cit->id()]*q6_Im_Tot[7*cit->id()]);
	    for(int f=1; f<=6; f++){
	      q6_bet                  += (2.0*(q6_Re_Tot[7*cit->id()+f]*q6_Re_Tot[7*cit->id()+f]+q6_Im_Tot[7*cit->id()+f]*q6_Im_Tot[7*cit->id()+f]));
	    }
	    q6_bet                    *= 4.0*PI/13.0;
	    q6_bet                     = sqrt(q6_bet);
	    int myIndex                = (bin[1]*dhBins[0]*dhBins[2])+(bin[2]*dhBins[0])+bin[0];
	    q6InCell[myIndex]         += q6_bet;
	    neighboursInCell[myIndex] += neighbours_Tot[cit->id()];
	    partsInCell[myIndex]++;
	  }
	}
      }
      real totalWCenterX           = 0.0;
      boost::mpi::all_reduce(*mpiWorld, localWCenterX, totalWCenterX, std::plus<real>());
      totalWCenterX               /= nthresh;
      int binWCenterX              = int(totalWCenterX/dr[0]);

      boost::mpi::all_reduce(*mpiWorld, (int*)&partsInCell[0], SIZE3D, (int*)&totpartsInCell[0], std::plus<int>());
      boost::mpi::all_reduce(*mpiWorld, (int*)&neighboursInCell[0], SIZE3D, (int*)&totneighboursInCell[0], std::plus<int>());
      boost::mpi::all_reduce(*mpiWorld, (real*)&q6InCell[0], SIZE3D, (real*)&totq6InCell[0], std::plus<real>());


      real localWallL[sumBins];
      real localWallR[sumBins];
      for(int i=0; i<sumBins; i++){
	localWallL[i]               = totalWCenterX;
	localWallR[i]               = totalWCenterX;
      }

      int limitsL[2]            = {std::max(binsToCheckL[0],0),std::min(binsToCheckL[1],binWCenterX-1)};
      int limitsR[2]            = {std::max(binsToCheckR[0],binWCenterX+1),std::min(binsToCheckR[1],dhBins[0]-1)};
      for(int i=0; i<dhBins[1]; i++){
	for(int j=0; j<dhBins[2]; j++){
          int myWallIndex           = (i*dhBins[2])+j;
	  for(int k=limitsL[1]; k>=limitsL[0]; k--){
	    int myCellIndex         = (myWallIndex*dhBins[0])+k;
	    if( (totpartsInCell[myCellIndex]) && (((real)totq6InCell[myCellIndex]/(real)totpartsInCell[myCellIndex])>=checkerQ6) && (((real)totneighboursInCell[myCellIndex]/(real)totpartsInCell[myCellIndex])>=checkerCN) ){
	      real cpos             = k*dr[0];
	      if( cpos<localWallL[myWallIndex] )
		localWallL[myWallIndex] = cpos;
	    }
	  }
	  for(int k=limitsR[0]; k<=limitsR[1]; k++){
	    int myCellIndex         = (myWallIndex*dhBins[0])+k;
	    if( (totpartsInCell[myCellIndex]) && (((real)totq6InCell[myCellIndex]/(real)totpartsInCell[myCellIndex])>=checkerQ6) && (((real)totneighboursInCell[myCellIndex]/(real)totpartsInCell[myCellIndex])>=checkerCN) ){
	      real cpos             = k*dr[0];
	      if( cpos>localWallR[myWallIndex] )
		localWallR[myWallIndex] = cpos;
	    }
	  }
        }
      }

      real minWallL                = Li[0];
      real maxWallR                = 0.0;
      for(int i=0; i<sumBins; i++){
	//wallInstanceL[i]           = localWallL[i];
	//minWallL                   = std::min(minWallL,wallInstanceL[i]);
	//wallInstanceR[i]           = localWallR[i];
	//maxWallR                   = std::max(maxWallR,wallInstanceR[i]);
	minWallL                   = std::min(minWallL,localWallL[i]);
	maxWallR                   = std::max(maxWallR,localWallR[i]);
      }
      //real minAllWallL;
      //real maxAllWallR;
      //boost::mpi::all_reduce(*mpiWorld, minWallL, minAllWallL, boost::mpi::minimum<real>());
      //boost::mpi::all_reduce(*mpiWorld, maxWallR, maxAllWallR, boost::mpi::maximum<real>());

      substrate[0]                 = totalWCenterX;
      real d_WallL                 = fabs(minWallL-substrate[1]);
      if( (d_WallL<maxGrowth) && ((d_WallL/dr[0])>0.5) )
	substrate[1]               = minWallL;
      real d_WallR                 = fabs(maxWallR-substrate[2]);
      if( (d_WallR<maxGrowth) && ((d_WallR/dr[0])>0.5) )
	substrate[2]               = maxWallR;

      center_CR_L                  = substrate[1]-sizeTR-(0.5*sizeCR);
      center_CR_R                  = substrate[2]+sizeTR+(0.5*sizeCR);

      center_FR_L                  = substrate[1]-sizeTR-sizeCR-(0.5*sizeFR);
      center_FR_R                  = substrate[2]+sizeTR+sizeCR+(0.5*sizeFR);

      reservoir_L                  = substrate[1]-sizeTR-sizeCR-sizeFR;
      reservoir_R                  = substrate[2]+sizeTR+sizeCR+sizeFR;

      sizeResL                     = reservoir_L;
      sizeResR                     = Li[0]-reservoir_R;
      calcSizeRes(sizeResL,sizeResR);

    }

    void ConstMuMD::assignParts(){

      // ############################################################################################
      // ## ASSIGN PARTICLES                                                                       ##
      // ############################################################################################

      System& system       = getSystemRef();
      CellList realCells   = system.storage->getRealCells();
      Real3D Li            = system.bc->getBoxL();
      real Lx              = Li[0];
      real latArea         = Li[1]*Li[2];
      real edgeshift       = 4.0*RB_WIDTH;
      real c_CutXtal[2]    = {substrate[1]+edgeshift,substrate[2]-edgeshift};
      real c_CutLiq_L[2]   = {edgeshift,substrate[1]-edgeshift};
      real c_CutLiq_R[2]   = {substrate[2]+edgeshift,Lx-edgeshift};
      real c_CutCR_L[2]    = {center_CR_L-0.5*sizeCR+edgeshift,center_CR_L+0.5*sizeCR-edgeshift};
      real c_CutCR_R[2]    = {center_CR_R-0.5*sizeCR+edgeshift,center_CR_R+0.5*sizeCR-edgeshift};
      real c_CutRes_L[2]   = {edgeshift,center_FR_L-0.5*sizeFR-edgeshift};
      real c_CutRes_R[2]   = {center_FR_R+0.5*sizeFR+edgeshift,Lx-edgeshift};

      real c_npartsXtal    = 0.0;
      real c_npartsLiq_L   = 0.0;
      real c_npartsLiq_R   = 0.0;
      real c_npartsCR_L[2] = {0.0,0.0};
      real c_npartsCR_R[2] = {0.0,0.0};
      real c_npartsRes_L   = 0.0;
      real c_npartsRes_R   = 0.0;
      real c_molFCR_L      = 0.0;
      real c_molFCR_R      = 0.0;

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit){
	real xpos          = cit->position()[0];
	if( xpos<0.0 ){
	  xpos            += Lx;
	}
	if( xpos>Lx ){
	  xpos            -= Lx;
	}

	c_npartsXtal                += smoothedHeaviside(xpos,c_CutXtal[0],c_CutXtal[1],RB_WIDTH);
	if( xpos<substrate[0] ){
	  c_npartsLiq_L             += smoothedHeaviside(xpos,c_CutLiq_L[0],c_CutLiq_L[1],RB_WIDTH);
	  c_npartsCR_L[cit->type()] += smoothedHeaviside(xpos,c_CutCR_L[0],c_CutCR_L[1],RB_WIDTH);
	  if( cit->type()==0 )
	    c_molFCR_L              += smoothedHeaviside(xpos,c_CutCR_L[0],c_CutCR_L[1],RB_WIDTH);
	  c_npartsRes_L             += smoothedHeaviside(xpos,c_CutRes_L[0],c_CutRes_L[1],RB_WIDTH);
	}
	if( xpos>=substrate[0] ){
	  c_npartsLiq_R             += smoothedHeaviside(xpos,c_CutLiq_R[0],c_CutLiq_R[1],RB_WIDTH);
	  c_npartsCR_R[cit->type()] += smoothedHeaviside(xpos,c_CutCR_R[0],c_CutCR_R[1],RB_WIDTH);
	  if( cit->type()==0 )
	    c_molFCR_R              += smoothedHeaviside(xpos,c_CutCR_R[0],c_CutCR_R[1],RB_WIDTH);
	  c_npartsRes_R             += smoothedHeaviside(xpos,c_CutRes_R[0],c_CutRes_R[1],RB_WIDTH);
	}
      }


      real tot_npartsXtal    = 0.0;
      real tot_npartsLiq_L   = 0.0;
      real tot_npartsLiq_R   = 0.0;
      real tot_npartsCR_L[2] = {0.0,0.0};
      real tot_npartsCR_R[2] = {0.0,0.0};
      real tot_npartsRes_L   = 0.0;
      real tot_npartsRes_R   = 0.0;
      real tot_molFCR_L      = 0.0;
      real tot_molFCR_R      = 0.0;

      boost::mpi::all_reduce(*mpiWorld, c_npartsXtal, tot_npartsXtal, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, c_npartsLiq_L, tot_npartsLiq_L, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, c_npartsLiq_R, tot_npartsLiq_R, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, (real*)&c_npartsCR_L, 2, (real*)&tot_npartsCR_L, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, (real*)&c_npartsCR_R, 2, (real*)&tot_npartsCR_R, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, c_npartsRes_L, tot_npartsRes_L, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, c_npartsRes_R, tot_npartsRes_R, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, c_molFCR_L, tot_molFCR_L, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, c_molFCR_R, tot_molFCR_R, std::plus<real>());

      if( substrate[2]!=substrate[1] ){
	numDensInXtal    = tot_npartsXtal/((c_CutXtal[1]-c_CutXtal[0])*latArea);
      }
      numDensInLiq       = Real3D(tot_npartsLiq_L/((c_CutLiq_L[1]-c_CutLiq_L[0])*latArea),tot_npartsLiq_R/((c_CutLiq_R[1]-c_CutLiq_R[0])*latArea),1.0);
      numDensInLiq[2]    = 0.5*(numDensInLiq[0]+numDensInLiq[1]);
      numDensInCR0       = Real3D(tot_npartsCR_L[0]/((c_CutCR_L[1]-c_CutCR_L[0])*latArea),tot_npartsCR_R[0]/((c_CutCR_R[1]-c_CutCR_R[0])*latArea),1.0);
      numDensInCR0[2]    = 0.5*(numDensInCR0[0]+numDensInCR0[1]);
      numDensInCR1       = Real3D(tot_npartsCR_L[1]/((c_CutCR_L[1]-c_CutCR_L[0])*latArea),tot_npartsCR_R[1]/((c_CutCR_R[1]-c_CutCR_R[0])*latArea),1.0);
      numDensInCR1[2]    = 0.5*(numDensInCR1[0]+numDensInCR1[1]);
      numDensInCRTot     = Real3D(numDensInCR0[0]+numDensInCR1[0],numDensInCR0[1]+numDensInCR1[1],1.0);
      numDensInCRTot[2]  = 0.5*(numDensInCRTot[0]+numDensInCRTot[1]);
      numDensInRes       = Real3D(tot_npartsRes_L/((sizeResL-2.0*edgeshift)*latArea),tot_npartsRes_R/((sizeResR-2.0*edgeshift)*latArea),1.0);
      numDensInRes[2]    = 0.5*(numDensInRes[0]+numDensInRes[1]);
      if( (tot_npartsCR_L[0]+tot_npartsCR_L[1]>0.0) && (tot_npartsCR_R[0]+tot_npartsCR_R[1]>0.0) ){
	molFracInCR      = Real3D(tot_molFCR_L/(tot_npartsCR_L[0]+tot_npartsCR_L[1]),tot_molFCR_R/(tot_npartsCR_R[0]+tot_npartsCR_R[1]),1.0);
	molFracInCR[2]   = 0.5*(molFracInCR[0]+molFracInCR[1]);
      }

    }

    void ConstMuMD::applyForce(){

      //assignParts();

      // ############################################################################################
      // ## APPLY FORCE                                                                            ##
      // ############################################################################################

      System& system         = getSystemRef();
      CellList realCells     = system.storage->getRealCells();
      real Lx                = system.bc->getBoxL()[0];
      real c_CutFR           = (0.5*sizeFR)*(0.5*sizeFR);

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit){
	real xpos          = cit->position()[0];
	if( xpos<0.0 ){
	  xpos            += Lx;
	}
	if( xpos>Lx ){
	  xpos            -= Lx;
	}

	real cdist2_FR_L   = (xpos-center_FR_L)*(xpos-center_FR_L);
	real cdist2_FR_R   = (xpos-center_FR_R)*(xpos-center_FR_R);
	if( cdist2_FR_L<c_CutFR ){
	  if( cit->type()==0 ){
            if( (modus==0) || (modus==2) ){
	      real distFac     = 1.0+cosh(sqrt(cdist2_FR_L)*shapeParam);
	      //real memForce    = fConst0*(expDens0-numDensInCR0[0])*0.5*shapeParam/distFac;
	      real memForce    = fConst0*(expMolF0-molFracInCR[0])*0.5*shapeParam/distFac;
	      cit->force()[0] += memForce;
            }
	  }
	  if( cit->type()==1 ){
            if( (modus==1) || (modus==2) ){
	      real distFac     = 1.0+cosh(sqrt(cdist2_FR_L)*shapeParam);
	      //real memForce    = fConst1*(expDens1-numDensInCR1[0])*0.5*shapeParam/distFac;
	      real memForce    = fConst1*(expMolF1+molFracInCR[0]-1.0)*0.5*shapeParam/distFac;
	      cit->force()[0] += memForce;
            }
	  }
	  continue;
	}
	if( cdist2_FR_R<c_CutFR ){
	  if( cit->type()==0 ){
            if( (modus==0) || (modus==2) ){
	      real distFac     = 1.0+cosh(sqrt(cdist2_FR_R)*shapeParam);
	      //real memForce    = fConst0*(numDensInCR0[1]-expDens0)*0.5*shapeParam/distFac;
	      real memForce    = fConst0*(molFracInCR[1]-expMolF0)*0.5*shapeParam/distFac;
	      cit->force()[0] += memForce;
            }
	  }
	  if( cit->type()==1 ){
            if( (modus==1) || (modus==2) ){
	      real distFac     = 1.0+cosh(sqrt(cdist2_FR_R)*shapeParam);
	      //real memForce    = fConst1*(numDensInCR1[1]-expDens1)*0.5*shapeParam/distFac;
	      real memForce    = fConst1*(1.0-molFracInCR[1]-expMolF1)*0.5*shapeParam/distFac;
	      cit->force()[0] += memForce;
            }
	  }
	  continue;
	}
      }

    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/
    void ConstMuMD::registerPython() {

      using namespace espressopp::python;

      class_<ConstMuMD, shared_ptr<ConstMuMD>, bases<Extension> >

        ("integrator_ConstMuMD", init< shared_ptr<System> >())

        .add_property("modus", &ConstMuMD::getModus, &ConstMuMD::setModus)
        .add_property("sizeTR", &ConstMuMD::getSizeTR, &ConstMuMD::setSizeTR)
        .add_property("sizeCR", &ConstMuMD::getSizeCR, &ConstMuMD::setSizeCR)
        .add_property("sizeFR", &ConstMuMD::getSizeFR, &ConstMuMD::setSizeFR)
        .add_property("sizeCheckR", &ConstMuMD::getSizeCheckR, &ConstMuMD::setSizeCheckR)
        .add_property("sizeRes", &ConstMuMD::getSizeRes, &ConstMuMD::setSizeRes)
        .add_property("dr", &ConstMuMD::getDR, &ConstMuMD::setDR)
        .add_property("binsInPlane", &ConstMuMD::getBinsInPlane, &ConstMuMD::setBinsInPlane)
        .add_property("nThresh", &ConstMuMD::getNThresh, &ConstMuMD::setNThresh)
        .add_property("cutoff", &ConstMuMD::getCutoff, &ConstMuMD::setCutoff)
        .add_property("shapeParam", &ConstMuMD::getShapeParam, &ConstMuMD::setShapeParam)
        .add_property("fConst", &ConstMuMD::getForceConst, &ConstMuMD::setForceConst)
        // 	.add_property("expDens", &ConstMuMD::getExpDens, &ConstMuMD::setExpDens)
        .add_property("expMolF", &ConstMuMD::getExpMolF, &ConstMuMD::setExpMolF)
        .add_property("substrate", &ConstMuMD::getWallInstance, &ConstMuMD::setWallInstance)
        .add_property("Ntot", &ConstMuMD::getNtot, &ConstMuMD::setNtot)
        .add_property("ratio", &ConstMuMD::getRatio, &ConstMuMD::setRatio)
        //.add_property("structureL", &ConstMuMD::getWallL, &ConstMuMD::setWallL)
        //.add_property("structureR", &ConstMuMD::getWallR, &ConstMuMD::setWallR)

        .def("connect", &ConstMuMD::connect)
        .def("disconnect", &ConstMuMD::disconnect)
        .def("initialize", &ConstMuMD::initialize)
        .def("perform", &ConstMuMD::perform)
        .def("adaptRegions", &ConstMuMD::adaptRegions)
        .def("assignParts", &ConstMuMD::assignParts)
        .def("applyForce", &ConstMuMD::applyForce)
        .def("getDensityData", &ConstMuMD::getDensityData)
      ;
    }
  }
}
