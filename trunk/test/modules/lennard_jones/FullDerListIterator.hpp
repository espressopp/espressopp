
#include "ListIterator.hpp"

class FullListIterator : public ListIterator {

   private:

   int*  neighbors; 
   int   neighsize;

   public:

   FullListIterator(int _N) {
 
     inum       = _N;
     ilist      = (int *) malloc (sizeof(int) * inum);
     numneigh   = (int *) malloc (sizeof(int) * inum);
     firstneigh = (int**) malloc (sizeof(int*) * inum);

     neighsize = inum * (inum - 1) / 2;
     neighbors = (int *) malloc (sizeof(int) * neighsize);

     printf ("%d is size of neighbor array\n", neighsize);

     int* ptr = neighbors;

     for (int i = 0; i < inum; i++) {

        ilist[i]      = i;
        firstneigh[i] = ptr;
        numneigh[i]   = inum - (i+1);

        int nneigh = 0; 

        for (int j = i+1; j < inum; j++) {

            *ptr++ = j;
        }
  
     }

   }

};

class OptFullListIterator : public ListIterator {

   private:

   int*  neighbors; 
   int   neighsize;

   public:

   OptFullListIterator(int _N) {
 
     inum       = _N;
     ilist      = (int *) malloc (sizeof(int) * inum);
     numneigh   = (int *) malloc (sizeof(int) * inum);
     firstneigh = (int**) malloc (sizeof(int*) * inum);

     neighsize = inum - 1,
     neighbors = (int *) malloc (sizeof(int) * neighsize);

     // fill neighbors 

     for (int j = 0; j < inum; j++)

         neighbors[j] = j+1;

     for (int i = 0; i < inum; i++) {

        ilist[i]      = i;
        firstneigh[i] = neighbors + i;
        numneigh[i]   = inum - (i+1);

     }

   }

};

