
typedef int particleRef;

class BlockIterator {

   public:

   virtual ~BlockIterator() {}; 

   virtual void reset() = 0;
   virtual void next() = 0;
   virtual bool done() = 0;

   virtual void get (particleRef &_i, int &jk, int max, particleRef jlist[]) = 0;

};

