// blocked iteration. Here, the block size is a template parameter to easily
// try out several block sizes. In a implementation, one will probably want
// to fix this somewhere as a constant.

#include <vector>

// base class for an iterator that hands out blocks of elements at once
template<class T, int csize>
class VirtualBlockFiller {
public:
  virtual ~VirtualBlockFiller() {};
  virtual size_t fillCacheBlock(const T *cache[csize]) = 0;
  virtual void reset() = 0;
};

// Block iterator class, providing the data cache and access to it.
// Also fills up the cache through calls to the abstract filler class
template <class T, int csize>
class BlockedIterator {
public:
  size_t i, n;
  const T *cache[csize];

  VirtualBlockFiller<T, csize> *filler;

  BlockedIterator(VirtualBlockFiller<T, csize> *_filler): i(0), filler(_filler) {
    filler->reset();
    n = filler->fillCacheBlock(cache);
  };

  bool isDone() { return n == 0; };
  void next() {
    if (++i < n) { return; }
    i = 0;
    n = filler->fillCacheBlock(cache);
  };
  const T &getCurrent() { return *cache[i]; };
};

// Block Vector Filler, descendant of the Block filler base class
template <class T, int csize>
class VirtualVectorBlockFiller 
  : public VirtualBlockFiller<T, csize> {
public:
  size_t i, fill;
  const std::vector <T> &v;

  VirtualVectorBlockFiller(const std::vector <T> &v, size_t _fill = csize)
    : fill(_fill), v(v)
  { reset(); }

  virtual void reset() { i = 0; };

  virtual size_t fillCacheBlock(const T *cache[csize]) {
    size_t copy = min(v.size() - i, fill);
    // make a stack copy - using the class variable i in the loop directly
    // costs about 20% performance
    const T *cur = &(v[i]);
    for (size_t c = 0; c < copy; ++c) {
      cache[c] = cur++;
    }
    i += copy;
    return copy;
  }
};


