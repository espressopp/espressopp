
typedef int particleRef;

class PairIterator {

   public:

   virtual void reset() = 0;
   virtual particleRef first() = 0;
   virtual particleRef second() = 0;
   virtual void next() = 0;
   virtual bool done() = 0;
};

