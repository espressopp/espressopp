
typedef int particleRef;

class ListIterator {

   public:

   int   inum;       // number of particles
   int*  ilist;      // index list of all my particles
   int*  numneigh;   // numneigh[i] nuber of neighbors for i
   int** firstneigh; // firstneigh[i] is first neighbor of i

   ListIterator() {
 
      inum = 0;
      ilist = NULL;
      numneigh = NULL;
      firstneigh = NULL;
      
   }

};

