
typedef int particleRef;

class FullPairIterator {

   public:

   FullPairIterator(int _N) {

    N =  _N;
    
    reset();
   }

   void reset() {

    // set all pointers to the beginning of the list

    i = 0;
    j = i+1;

   }

   int first() {

    return i;
   }

   int second() {

    return j;
   }

   void get (particleRef& _i, particleRef& _j) {

     _i = i;
     _j = j;

   }

   void next() {

    j++;

    if (j < N) return;

    i++;
    j = i+1;

   // move pointers to the next pair in verlet list
   }

   bool done() {

    return i >= N;

   }


   int N;  // number of particles

   int i; 
   int j;
};

