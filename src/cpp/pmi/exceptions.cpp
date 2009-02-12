#include "pmi/exceptions.hpp"

#include <sstream>

using namespace std;
using namespace pmi;

UserError::
UserError(const string &_message) 
  : message("user error: " + _message) {}

const char* 
UserError::what() 
  const throw() {
  return message.c_str();
}

UserError::
~UserError() throw () {};

InternalError::
InternalError(const string &_message) 
  : message("internal error: " + _message) {}

const char* 
InternalError::what() 
  const throw() {
  return message.c_str();
}

InternalError::
~InternalError() throw () {};
