// ESPP_CLASS
#ifndef _ITERATOR_CELLLISTALLTRIPLESITERATOR_HPP
#define _ITERATOR_CELLLISTALLTRIPLESITERATOR_HPP

#include "Cell.hpp"
#include "log4espp.hpp"
#include "iterator/CellListIterator.hpp"
//#include "iterator/NeighborCellListIterator.hpp"

#include "esutil/Error.hpp"

using namespace std;

namespace espresso {
  namespace iterator {
    class CellListAllTriplesIterator {
    public:
      CellListAllTriplesIterator();
      CellListAllTriplesIterator(CellList &cl);
      
      CellListAllTriplesIterator &operator++();
      
      bool isValid() const;
      bool isDone() const;
      
      const ParticleTriple &operator*() const;
      const ParticleTriple *operator->() const;
      
    private:
      static LOG4ESPP_DECL_LOGGER(theLogger);

      ParticleTriple current;

      /* In general, we have 4 possibilities (?5):
       * 0) cit1 = cit2 = cit3:  inSelfLoop = 0
       * 1) cit1 != cit2; cit1 = cit3 and cit1 = cit2; cit1 != cit3: inSelfLoop = 1
       *      (!may be not the same)
       * 2) cit1 != cit2; cit2 = cit3:  inSelfLoop = 2
       * 3) cit1 != cit2 != cit3; cit1 != cit3: inSelfLoop = 3
       * 4) next cit1
       */
      short int inSelfLoop; // can be 0, 1, 2, 3
      
      bool inSelfLoop1;
      //bool inSelfLoop2;

      // current cell1
      CellList::Iterator cit1;
      // current particle1
      ParticleList::Iterator pit1;

      // we need two neighbors because three particles may be in different cells
      // current neighbor cell
      NeighborCellList::Iterator cit2;
      // current particle2
      ParticleList::Iterator pit2;

      // current neighbor cell
      NeighborCellList::Iterator cit3;
      // current particle3
      ParticleList::Iterator pit3;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    inline 
    CellListAllTriplesIterator::
    CellListAllTriplesIterator()
    {}

    inline 
    CellListAllTriplesIterator::
    CellListAllTriplesIterator(CellList &cl) {

      cit1 = CellList::Iterator(cl);
      if (cit1.isDone()) return;
      inSelfLoop = 0;
      pit1 = ParticleList::Iterator((*cit1)->particles);
      while (pit1.isDone()) {
        ++cit1;
        if (cit1.isDone()) return;
        pit1 = ParticleList::Iterator((*cit1)->particles);
      }
      
      pit2 = pit1;
      ++pit2;
      if(pit2.isDone()){
        // means that in cit1 only 1 particle
        inSelfLoop = 2;
        cit2 = NeighborCellList::Iterator((*cit1)->neighborCells);
        pit2 = ParticleList::Iterator(cit2->cell->particles);
      }
      while (pit2.isDone()) {
        ++cit2;
        if (cit2.isDone()){
          // it is clear that there is no more particles in system
          while(cit1.isValid()) ++cit1;
          return;
        }
        pit2 = ParticleList::Iterator(cit2->cell->particles);
      }
              
      pit3 = pit2;
      this->operator++();
    }

    inline CellListAllTriplesIterator &
    CellListAllTriplesIterator::
    operator++() {
//pit1 = ParticleList::Iterator((*cit1)->particles);
//cit2 = NeighborCellList::Iterator((*cit1)->neighborCells);
//pit2 = ParticleList::Iterator(cit2->cell->particles);
//std::cout << "p1: "<< pit1->id() << " p2isd: "<<pit2.isValid()<<std::endl;
//if(pit2.isDone()) cout << "Something strange pit2 isDone by inSelfLoop="<< inSelfLoop <<endl;
//exit(0);
      bool shouldNotReturn = true;
      while(shouldNotReturn){
        
        shouldNotReturn = false;
        ++pit3;
        if(pit3.isValid() && pit3->id()==pit1->id()){
          ++pit3;
        }
        if(inSelfLoop==2 && pit3.isValid() && pit3->id()==pit2->id()){
          ++pit3;
        }
        while (pit3.isDone()) {
          switch ( inSelfLoop ) {
            case 0:
              ++pit2;
              if(pit2.isValid() && pit2->id()==pit1->id()){
                ++pit2;
              }
              break;
            case 1:
              ++cit3;
              if(cit3.isDone()){
                ++pit2;
                if(pit2.isValid() && pit2->id()==pit1->id()){
                  ++pit2;
                }
                cit3 = NeighborCellList::Iterator((*cit1)->neighborCells);
              }
              break;
            case 2:
              ++pit2;
              break;
            case 3:
              ++cit3;
              if(cit3->cell == cit2->cell){
                ++cit3;
              }
              if(cit3.isDone()){
                ++pit2;
                // **
                //cit3 = NeighborCellList::Iterator((*cit1)->neighborCells);
                cit3 = cit2;
                ++cit3;
                if(cit3.isDone()){
                  cout<<"cit3 is done. But it should not happen"<<endl;
                }
              }
              
              break;
            default:
              cout<<"It is default. pit3 is done. BEGIN"<<endl;
          }
          
          while (pit2.isDone()) {
            
            switch ( inSelfLoop ) {
              case 0:
                ++pit1;
                break;
              case 1:
                ++pit1;
                break;
              case 2:
                ++cit2;
                if(cit2.isDone()){
                  ++pit1;
                  cit2 = NeighborCellList::Iterator((*cit1)->neighborCells);
                }
                break;
              case 3:
                ++cit2;
                if(cit2.isDone()){
                  ++pit1;
                  cit2 = NeighborCellList::Iterator((*cit1)->neighborCells);
                }
                cit3 = cit2;
                ++cit3;
                if(cit3.isDone()){
                  ++pit1;
                  cit2 = NeighborCellList::Iterator((*cit1)->neighborCells);
                  cit3 = cit2;
                  ++cit3;
                }
                
                // **
                //cit3 = NeighborCellList::Iterator((*cit1)->neighborCells);
                //if(cit3->cell == cit2->cell){
                //  ++cit3;
                //}
                break;
              default:
                cout<<"It is default. pit3 is done. BEGIN"<<endl;
            }
            
            while (pit1.isDone()) {

              switch ( inSelfLoop ) {
                case 0:
                  inSelfLoop = 1;
                  pit1 = ParticleList::Iterator((*cit1)->particles);
                  cit3 = NeighborCellList::Iterator((*cit1)->neighborCells);
                  //pit3 = ParticleList::Iterator(cit3->cell->particles);
                  //cout << "HHHH1: "<< inSelfLoop <<endl;
                  break;
                case 1:
                  inSelfLoop = 2;
                  pit1 = ParticleList::Iterator((*cit1)->particles);
                  cit2 = NeighborCellList::Iterator((*cit1)->neighborCells);
                  //shouldNotReturn = true;
                  break;
                case 2:
                  inSelfLoop = 3;
                  pit1 = ParticleList::Iterator((*cit1)->particles);
                  cit2 = NeighborCellList::Iterator((*cit1)->neighborCells);
                  cit3 = cit2;
                  ++cit3;
                  break;
                case 3:
                  ++cit1;
                  if(cit1.isDone())return *this;
                  inSelfLoop = 0;
                  pit1 = ParticleList::Iterator((*cit1)->particles);
                  while (pit1.isDone()) {
                    ++cit1;
                    if (cit1.isDone()) return *this;
                    pit1 = ParticleList::Iterator((*cit1)->particles);
                  }
                  shouldNotReturn = true;
                  break;
                default:
                  cout<<"It is default"<<endl;
              }
              //cout << "After "<< inSelfLoop << " pit1 "<< pit1.isValid() <<endl;

            }
            // pit 2 is done
            switch ( inSelfLoop ) {
              case 0:
                pit2 = ParticleList::Iterator((*cit1)->particles);
                if(pit2->id()==pit1->id()){
                  ++pit2;
                }
                break;
              case 1:
                pit2 = ParticleList::Iterator((*cit1)->particles);
                if(pit2->id()==pit1->id()){
                  ++pit2;
                }
                break;
              case 2:
                pit2 = ParticleList::Iterator(cit2->cell->particles);
                break;
              case 3:
                pit2 = ParticleList::Iterator(cit2->cell->particles);
                break;
              default:
                cout<<"It is default. pit2 is done"<<endl;
            }
          }

          // pit 3 is done
          switch ( inSelfLoop ) {
            case 0:
              pit3 = pit2;
              ++pit3;
              if(pit3->id()==pit1->id()){
                ++pit3;
              }
              break;
            case 1:
              pit3 = ParticleList::Iterator(cit3->cell->particles);
              break;
            case 2:
              pit3 = pit2;
              ++pit3;
              break;
            case 3:
              pit3 = ParticleList::Iterator(cit3->cell->particles);
              break;
            default:
              cout<<"It is default. pit3 is done"<<endl;
          }
        }
      
      }
      
      current.second = &*pit1;
      if(pit2->id()<pit3->id()){
        current.first  = &*pit2;
        current.third  = &*pit3;
      }
      else{
        current.first  = &*pit3;
        current.third  = &*pit2;
      }

      if(current.first->id()==0 && current.second->id()==6 &&  current.third->id()==5){
cout << "case1 part3 "<< inSelfLoop <<endl;
        
      std::cout << "current triple: (" << current.first->id() << ", " 
              << current.second->id() << ", " << current.third->id() << ")";
      std::cout << "  isGhost: (" << current.first->ghost() << ", " 
              << current.second->ghost() << ", " << current.third->ghost() << ")";
      std::cout << "  real: (" << pit2->id() << ", " 
              << pit1->id() << ", " << pit3->id() << ")"<<std::endl;
      }
      
      LOG4ESPP_TRACE(theLogger,
             "current triple: (" << current.first->p.id <<
             ", " << current.second->p.id << ", " << current.third->p.id << ")"
             );
      return *this;
    }

    inline bool 
    CellListAllTriplesIterator::
    isValid() const 
    { return cit1.isValid(); }

    inline bool 
    CellListAllTriplesIterator::
    isDone() const 
    { return !isValid(); }
      
    inline const ParticleTriple &
    CellListAllTriplesIterator::
    operator*() const 
    { return current; }
    
    inline const ParticleTriple *
    CellListAllTriplesIterator::
    operator->() const 
    { return &(**this); }

  }
}

#endif
