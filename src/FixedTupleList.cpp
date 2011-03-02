#include "python.hpp"

#include "FixedTupleList.hpp"

//#include <vector>
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

//using namespace std;

namespace espresso {
  FixedTupleList::FixedTupleList(shared_ptr< storage::Storage > _storage)
  : FixedListComm (_storage){};





  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedTupleList::registerPython() {

    using namespace espresso::python;

    bool (FixedTupleList::*pyAdd)(pvec pids)
      = &FixedTupleList::add;

    class_< FixedTupleList, shared_ptr <FixedTupleList> >
      ("FixedTupleList", init <shared_ptr <storage::Storage> >())
      .def("add", pyAdd)
      ;
  }
}
