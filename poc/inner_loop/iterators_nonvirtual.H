#include <vector>

// Vector Iterator that does not use virtual methods

template <class T>
class VectorIterator {
public:
  size_t i;
  const std::vector <T> &v;

  // constructor
  VectorIterator(std::vector <T> &v) 
    : v(v)
  { i = 0; }

  bool isDone() {
    return i >= v.size();
  }

  void next() {
    i++;
  }

  const T &getCurrent() {
    return v[i];
  }
};
