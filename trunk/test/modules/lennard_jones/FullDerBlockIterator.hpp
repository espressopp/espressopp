#include "BlockIterator.hpp"

class FullBlockIterator : public BlockIterator {

   public:

   FullBlockIterator(int _N) {

    N =  _N;
    
    reset();
   }

   void reset() {

    // set all pointers to the beginning of the list

    i = 0;
    j = i+1;

   }

   void get (particleRef &_i, int &jk, int max, particleRef jlist[]) {

    _i = i;
    jk = 0;

    for (int k = 0; k < max; k++) {
       if (j >= N) break;
       jlist[jk++] = j++;
       if (jk >= max) break;
    }

   }

   void next() {

    if (j < N) return;

    i++;
    j = i+1;

   // move pointers to the next pair in verlet list
   }

   bool done() {

    return i >= N;

   }

   private:

   int N;  // number of particles

   int i; 
   int j;
};

