#ifndef __PYX_HAVE__mpi4py__MPI
#define __PYX_HAVE__mpi4py__MPI
#ifdef __cplusplus
#define __PYX_EXTERN_C extern "C"
#else
#define __PYX_EXTERN_C extern
#endif

/* "include/mpi4py/MPI.pxd":50
 * ctypedef MPI_Offset Offset
 * 
 * ctypedef public api class Status [type PyMPIStatus_Type, object PyMPIStatusObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Status ob_mpi
 *     cdef int        flags
 */

struct __pyx_obj_6mpi4py_3MPI_Status {
  PyObject_HEAD
  MPI_Status ob_mpi;
  int flags;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Status PyMPIStatusObject;

/* "include/mpi4py/MPI.pxd":54
 *     cdef int        flags
 * 
 * ctypedef public api class Datatype [type PyMPIDatatype_Type, object PyMPIDatatypeObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Datatype ob_mpi
 *     cdef int          flags
 */

struct __pyx_obj_6mpi4py_3MPI_Datatype {
  PyObject_HEAD
  MPI_Datatype ob_mpi;
  int flags;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Datatype PyMPIDatatypeObject;

/* "include/mpi4py/MPI.pxd":58
 *     cdef int          flags
 * 
 * ctypedef public api class Request [type PyMPIRequest_Type, object PyMPIRequestObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Request ob_mpi
 *     cdef int         flags
 */

struct __pyx_obj_6mpi4py_3MPI_Request {
  PyObject_HEAD
  MPI_Request ob_mpi;
  int flags;
  PyObject *ob_buf;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Request PyMPIRequestObject;

/* "include/mpi4py/MPI.pxd":63
 *     cdef object      ob_buf
 * 
 * ctypedef public api class Prequest(Request) [type PyMPIPrequest_Type, object PyMPIPrequestObject]:             # <<<<<<<<<<<<<<
 *     pass
 * 
 */

struct __pyx_obj_6mpi4py_3MPI_Prequest {
  struct __pyx_obj_6mpi4py_3MPI_Request __pyx_base;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Prequest PyMPIPrequestObject;

/* "include/mpi4py/MPI.pxd":66
 *     pass
 * 
 * ctypedef public api class Grequest(Request) [type PyMPIGrequest_Type, object PyMPIGrequestObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Request ob_grequest
 * 
 */

struct __pyx_obj_6mpi4py_3MPI_Grequest {
  struct __pyx_obj_6mpi4py_3MPI_Request __pyx_base;
  MPI_Request ob_grequest;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Grequest PyMPIGrequestObject;

/* "include/mpi4py/MPI.pxd":69
 *     cdef MPI_Request ob_grequest
 * 
 * ctypedef public api class Op [type PyMPIOp_Type, object PyMPIOpObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Op ob_mpi
 *     cdef int    flags
 */

struct __pyx_obj_6mpi4py_3MPI_Op {
  PyObject_HEAD
  MPI_Op ob_mpi;
  int flags;
  PyObject *(*ob_func)(PyObject *, PyObject *);
  int ob_usrid;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Op PyMPIOpObject;

/* "include/mpi4py/MPI.pxd":75
 *     cdef int    ob_usrid
 * 
 * ctypedef public api class Group [type PyMPIGroup_Type, object PyMPIGroupObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Group ob_mpi
 *     cdef int       flags
 */

struct __pyx_obj_6mpi4py_3MPI_Group {
  PyObject_HEAD
  MPI_Group ob_mpi;
  int flags;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Group PyMPIGroupObject;

/* "include/mpi4py/MPI.pxd":79
 *     cdef int       flags
 * 
 * ctypedef public api class Info [type PyMPIInfo_Type, object PyMPIInfoObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Info ob_mpi
 *     cdef int      flags
 */

struct __pyx_obj_6mpi4py_3MPI_Info {
  PyObject_HEAD
  MPI_Info ob_mpi;
  int flags;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Info PyMPIInfoObject;

/* "include/mpi4py/MPI.pxd":83
 *     cdef int      flags
 * 
 * ctypedef public api class Errhandler [type PyMPIErrhandler_Type, object PyMPIErrhandlerObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Errhandler ob_mpi
 *     cdef int            flags
 */

struct __pyx_obj_6mpi4py_3MPI_Errhandler {
  PyObject_HEAD
  MPI_Errhandler ob_mpi;
  int flags;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Errhandler PyMPIErrhandlerObject;

/* "include/mpi4py/MPI.pxd":87
 *     cdef int            flags
 * 
 * ctypedef public api class Comm [type PyMPIComm_Type, object PyMPICommObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Comm ob_mpi
 *     cdef int      flags
 */

struct __pyx_obj_6mpi4py_3MPI_Comm {
  PyObject_HEAD
  MPI_Comm ob_mpi;
  int flags;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Comm PyMPICommObject;

/* "include/mpi4py/MPI.pxd":91
 *     cdef int      flags
 * 
 * ctypedef public api class Intracomm(Comm) [type PyMPIIntracomm_Type, object PyMPIIntracommObject]:             # <<<<<<<<<<<<<<
 *     pass
 * 
 */

struct __pyx_obj_6mpi4py_3MPI_Intracomm {
  struct __pyx_obj_6mpi4py_3MPI_Comm __pyx_base;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Intracomm PyMPIIntracommObject;

/* "include/mpi4py/MPI.pxd":94
 *     pass
 * 
 * ctypedef public api class Cartcomm(Intracomm) [type PyMPICartcomm_Type, object PyMPICartcommObject]:             # <<<<<<<<<<<<<<
 *     pass
 * 
 */

struct __pyx_obj_6mpi4py_3MPI_Cartcomm {
  struct __pyx_obj_6mpi4py_3MPI_Intracomm __pyx_base;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Cartcomm PyMPICartcommObject;

/* "include/mpi4py/MPI.pxd":97
 *     pass
 * 
 * ctypedef public api class Graphcomm(Intracomm) [type PyMPIGraphcomm_Type, object PyMPIGraphcommObject]:             # <<<<<<<<<<<<<<
 *     pass
 * 
 */

struct __pyx_obj_6mpi4py_3MPI_Graphcomm {
  struct __pyx_obj_6mpi4py_3MPI_Intracomm __pyx_base;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Graphcomm PyMPIGraphcommObject;

/* "include/mpi4py/MPI.pxd":100
 *     pass
 * 
 * ctypedef public api class Distgraphcomm(Intracomm) [type PyMPIDistgraphcomm_Type, object PyMPIDistgraphcommObject]:             # <<<<<<<<<<<<<<
 *     pass
 * 
 */

struct __pyx_obj_6mpi4py_3MPI_Distgraphcomm {
  struct __pyx_obj_6mpi4py_3MPI_Intracomm __pyx_base;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Distgraphcomm PyMPIDistgraphcommObject;

/* "include/mpi4py/MPI.pxd":103
 *     pass
 * 
 * ctypedef public api class Intercomm(Comm) [type PyMPIIntercomm_Type, object PyMPIIntercommObject]:             # <<<<<<<<<<<<<<
 *     pass
 * 
 */

struct __pyx_obj_6mpi4py_3MPI_Intercomm {
  struct __pyx_obj_6mpi4py_3MPI_Comm __pyx_base;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Intercomm PyMPIIntercommObject;

/* "include/mpi4py/MPI.pxd":106
 *     pass
 * 
 * ctypedef public api class Win [type PyMPIWin_Type, object PyMPIWinObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_Win ob_mpi
 *     cdef int     flags
 */

struct __pyx_obj_6mpi4py_3MPI_Win {
  PyObject_HEAD
  MPI_Win ob_mpi;
  int flags;
};
typedef struct __pyx_obj_6mpi4py_3MPI_Win PyMPIWinObject;

/* "include/mpi4py/MPI.pxd":110
 *     cdef int     flags
 * 
 * ctypedef public api class File [type PyMPIFile_Type, object PyMPIFileObject]:             # <<<<<<<<<<<<<<
 *     cdef MPI_File ob_mpi
 *     cdef int      flags
 */

struct __pyx_obj_6mpi4py_3MPI_File {
  PyObject_HEAD
  MPI_File ob_mpi;
  int flags;
};
typedef struct __pyx_obj_6mpi4py_3MPI_File PyMPIFileObject;

#ifndef __PYX_HAVE_API__mpi4py__MPI

__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIStatus_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIDatatype_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIRequest_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIPrequest_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIGrequest_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIOp_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIGroup_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIInfo_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIErrhandler_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIComm_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIIntracomm_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPICartcomm_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIGraphcomm_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIDistgraphcomm_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIIntercomm_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIWin_Type;
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMPIFile_Type;

#endif

PyMODINIT_FUNC initMPI(void);

#endif
