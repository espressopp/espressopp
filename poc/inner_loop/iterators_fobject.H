#include <vector>
#include <functional>

template<class T>
class Computer: public std::unary_function<const T, T> {
public:
  virtual void operator()(const T &) = 0;
};

// class with virtual foreach
template <class T>
class VirtualForeach {
public:
  virtual void foreach(Computer<T> &) const = 0;
};

template <class T>
class VirtualVectorForeach
  : public VirtualForeach<T> {
public:
  const std::vector <T> &v;

  // constructor
  VirtualVectorForeach(std::vector <T> &_v) 
    : v(_v) {}

  virtual void foreach(Computer<T> &func) const {
    typename std::vector<T>::const_iterator it = v.begin();
    typename std::vector<T>::const_iterator end = v.end();
    for (; it != end; ++it)
      func(*it);
  }
};


