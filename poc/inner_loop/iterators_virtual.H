#include <vector>

// Iterator base class with virtual methods
template <class T>
class VirtualIterator {
public:
  virtual bool isDone() const = 0;
  virtual void next() = 0;
  virtual void reset() = 0;
  virtual const T &getCurrent() const = 0;
};

// Vector Iterator, descendant of the Iterator base class
// uses virtual methods
// template <class T>
// class VirtualVectorIterator 
//   : public VirtualIterator<T> {
// public:
//   size_t i;
//   size_t size;
//   const std::vector <T> &v;

//   // constructor
//   VirtualVectorIterator(std::vector <T> &_v) 
//     : v(_v)
//   { reset(); }

//   virtual void reset() { i = 0; size=v.size();}

//   virtual bool isDone() const {
//     return i >= size;
//   }

//   virtual void next() { i++; }

//   virtual const T &getCurrent() const {
//     return v[i];
//   }
// };

template <class T>
class VirtualVectorIterator 
  : public VirtualIterator<T> {
public:
  const std::vector <T> &v;
  typename std::vector<T>::const_iterator it;

  // constructor
  VirtualVectorIterator(std::vector <T> &_v) 
    : v(_v)
  { reset(); }

  virtual void reset() {
    it = v.begin();
  }

  virtual bool isDone() const {
    return it == v.end();
  }

  virtual void next() {
    it++;
  }

  virtual const T &getCurrent() const {
    return *it;
  }
};


