#ifndef ALLOCATOR_H_INCLUDED

#define ALLOCATOR_H_INCLUDED

#if __cplusplus >= 201103L
#include <memory>
using std::shared_ptr;
#else
#include <tr1/memory>
using std::tr1::shared_ptr;
#endif

template<class T>

class Allocator {

public:
  typedef shared_ptr< Allocator<T> > SelfPtr;


  virtual ~Allocator() {};

  virtual T* allocate() = 0;

};




#endif // ALLOCATOR_H_INCLUDED

