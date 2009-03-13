
typedef int particleRef;

class PairIterator {

   public:

   virtual ~PairIterator() {}
   virtual void reset() {}
   virtual void get (particleRef&, particleRef&) {}
   virtual void next() {};
   virtual bool done() {}
};

