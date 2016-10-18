from mpi4py import MPI
import mpiunittest as unittest
import arrayimpl

typemap = MPI._typedict

def mkzeros(n):
    import sys
    if not hasattr(sys, 'pypy_version_info'):
        return bytearray(n)
    return b'\0' * n

def memzero(m):
    n = len(m)
    if n == 0: return
    try:
        zero = b'\0'
        m[0] = zero
    except TypeError:
        zero = 0
        m[0] = zero
    for i in range(n):
        m[i] = zero

class BaseTestRMA(object):

    COMM = MPI.COMM_NULL
    INFO = MPI.INFO_NULL

    COUNT_MIN = 0

    def setUp(self):
        nbytes = 100*MPI.DOUBLE.size
        try:
            self.mpi_memory = MPI.Alloc_mem(nbytes)
            self.memory = self.mpi_memory
            memzero(self.memory)
        except MPI.Exception:
            from array import array
            self.mpi_memory = None
            self.memory = array('B',[0]*nbytes)
        self.WIN = MPI.Win.Create(self.memory, 1, self.INFO, self.COMM)

    def tearDown(self):
        self.WIN.Free()
        if self.mpi_memory:
            MPI.Free_mem(self.mpi_memory)

    def testPutGet(self):
        group = self.WIN.Get_group()
        size = group.Get_size()
        group.Free()
        for array in arrayimpl.ArrayTypes:
            for typecode in arrayimpl.TypeMap:
                for count in range(self.COUNT_MIN, 10):
                    for rank in range(size):
                        sbuf = array(range(count), typecode)
                        rbuf = array(-1, typecode, count+1)
                        #
                        self.WIN.Fence()
                        self.WIN.Put(sbuf.as_mpi(), rank)
                        self.WIN.Fence()
                        self.WIN.Get(rbuf.as_mpi_c(count), rank)
                        self.WIN.Fence()
                        for i in range(count):
                            self.assertEqual(sbuf[i], i)
                            self.assertNotEqual(rbuf[i], -1)
                        self.assertEqual(rbuf[-1], -1)
                        #
                        sbuf = array(range(count), typecode)
                        rbuf = array(-1, typecode, count+1)
                        target  = sbuf.itemsize
                        self.WIN.Fence()
                        self.WIN.Put(sbuf.as_mpi(), rank, target)
                        self.WIN.Fence()
                        self.WIN.Get(rbuf.as_mpi_c(count), rank, target)
                        self.WIN.Fence()
                        for i in range(count):
                            self.assertEqual(sbuf[i], i)
                            self.assertNotEqual(rbuf[i], -1)
                        self.assertEqual(rbuf[-1], -1)
                        #
                        sbuf = array(range(count), typecode)
                        rbuf = array(-1, typecode, count+1)
                        datatype = typemap[typecode]
                        target  = (sbuf.itemsize, count, datatype)
                        self.WIN.Fence()
                        self.WIN.Put(sbuf.as_mpi(), rank, target)
                        self.WIN.Fence()
                        self.WIN.Get(rbuf.as_mpi_c(count), rank, target)
                        self.WIN.Fence()
                        for i in range(count):
                            self.assertEqual(sbuf[i], i)
                            self.assertNotEqual(rbuf[i], -1)
                        self.assertEqual(rbuf[-1], -1)

    def testAccumulate(self):
        group = self.WIN.Get_group()
        size = group.Get_size()
        group.Free()
        for array in arrayimpl.ArrayTypes:
            for typecode in arrayimpl.TypeMap:
                for count in range(self.COUNT_MIN, 10):
                    for rank in range(size):
                        sbuf = array(range(count), typecode)
                        rbuf = array(-1, typecode, count+1)
                        for op in (MPI.SUM, MPI.PROD, MPI.MAX, MPI.MIN):
                            self.WIN.Fence()
                            self.WIN.Accumulate(sbuf.as_mpi(), rank, op=op)
                            self.WIN.Fence()
                            self.WIN.Get(rbuf.as_mpi_c(count), rank)
                            self.WIN.Fence()
                            for i in range(count):
                                self.assertEqual(sbuf[i], i)
                                self.assertNotEqual(rbuf[i], -1)
                            self.assertEqual(rbuf[-1], -1)

    def testGetAccumulate(self):
        group = self.WIN.Get_group()
        size = group.Get_size()
        rank = group.Get_rank()
        group.Free()
        self.WIN.Fence()
        obuf = MPI.Alloc_mem(1); memzero(obuf)
        rbuf = MPI.Alloc_mem(1); memzero(rbuf)
        try:
            try:
                self.WIN.Get_accumulate([obuf, 0, MPI.BYTE], [rbuf, 0, MPI.BYTE], rank)
            finally:
                MPI.Free_mem(obuf)
                MPI.Free_mem(rbuf)
        except NotImplementedError:
            return
        self.WIN.Fence()
        for array in arrayimpl.ArrayTypes:
            for typecode in arrayimpl.TypeMap:
                for count in range(self.COUNT_MIN, 10):
                    for rank in range(size):
                        ones = array([1]*count, typecode)
                        sbuf = array(range(count), typecode)
                        rbuf = array(-1, typecode, count+1)
                        gbuf = array(-1, typecode, count+1)
                        for op in (MPI.SUM, MPI.PROD,
                                   MPI.MAX, MPI.MIN,
                                   MPI.REPLACE, MPI.NO_OP):
                            self.WIN.Lock(rank)
                            self.WIN.Put(ones.as_mpi(), rank)
                            self.WIN.Get_accumulate(sbuf.as_mpi(),
                                                    rbuf.as_mpi_c(count),
                                                    rank, op=op)
                            self.WIN.Get(gbuf.as_mpi_c(count), rank)
                            self.WIN.Unlock(rank)
                            #
                            for i in range(count):
                                self.assertEqual(sbuf[i], i)
                                self.assertEqual(rbuf[i], 1)
                                self.assertEqual(gbuf[i], op(1, i))
                            self.assertEqual(rbuf[-1], -1)
                            self.assertEqual(gbuf[-1], -1)

    def testPutProcNull(self):
        self.WIN.Fence()
        self.WIN.Put(None, MPI.PROC_NULL, None)
        self.WIN.Fence()

    def testGetProcNull(self):
        self.WIN.Fence()
        self.WIN.Get(None, MPI.PROC_NULL, None)
        self.WIN.Fence()

    def testAccumulateProcNullReplace(self):
        self.WIN.Fence()
        zeros = mkzeros(8)
        self.WIN.Fence()
        self.WIN.Accumulate([zeros, MPI.INT], MPI.PROC_NULL, None, MPI.REPLACE)
        self.WIN.Fence()
        self.WIN.Accumulate([zeros, MPI.INT], MPI.PROC_NULL, None, MPI.REPLACE)
        self.WIN.Fence()

    def testAccumulateProcNullSum(self):
        self.WIN.Fence()
        zeros = mkzeros(8)
        self.WIN.Fence()
        self.WIN.Accumulate([zeros, MPI.INT], MPI.PROC_NULL, None, MPI.SUM)
        self.WIN.Fence()
        self.WIN.Accumulate([None, MPI.INT], MPI.PROC_NULL, None, MPI.SUM)
        self.WIN.Fence()

    def testFence(self):
        win = self.WIN
        LMODE = [0, MPI.MODE_NOSTORE, MPI.MODE_NOPUT,
                 MPI.MODE_NOSTORE|MPI.MODE_NOPUT]
        GMODE = [0, MPI.MODE_NOPRECEDE, MPI.MODE_NOSUCCEED]
        win.Fence()
        for lmode in LMODE:
            for gmode in GMODE:
                assertion =  lmode | gmode
                win.Fence(assertion)
        win.Fence()

    def testFenceAll(self):
        win = self.WIN
        assertion = 0
        modes = [0,
                 MPI.MODE_NOSTORE,
                 MPI.MODE_NOPUT,
                 MPI.MODE_NOPRECEDE,
                 MPI.MODE_NOSUCCEED]
        win.Fence()
        for mode in modes:
            win.Fence(mode)
            assertion |= mode
            win.Fence(assertion)
        win.Fence()

    def testStartComplete(self):
        self.WIN.Start(MPI.GROUP_EMPTY)
        self.WIN.Complete()

    def testPostWait(self):
        self.WIN.Post(MPI.GROUP_EMPTY)
        self.WIN.Wait()

    def testStartCompletePostWait(self):
        win = self.WIN
        wingroup = win.Get_group()
        size = wingroup.Get_size()
        rank = wingroup.Get_rank()
        if size < 2: return wingroup.Free()
        if rank == 0:
            group = wingroup.Excl([0])
            win.Start(group)
            win.Complete()
            win.Post(group)
            win.Wait()
            group.Free()
        else:
            group = wingroup.Incl([0])
            win.Post(group)
            win.Wait()
            win.Start(group)
            win.Complete()
            group.Free()
        wingroup.Free()

    def testStartCompletePostTest(self):
        comm = self.COMM
        win = self.WIN
        wingroup = win.Get_group()
        size = wingroup.Get_size()
        rank = wingroup.Get_rank()
        if size < 2: return wingroup.Free()
        if rank == 0:
            group = wingroup.Excl([0])
            win.Start(group)
            comm.Barrier()
            win.Complete()
            comm.Barrier()
            group.Free()
        else:
            group = wingroup.Incl([0])
            win.Post(group)
            flag = win.Test()
            self.assertFalse(flag)
            comm.Barrier()
            comm.Barrier()
            flag = win.Test()
            self.assertTrue(flag)
            group.Free()
        wingroup.Free()

    def testFlushSync(self):
        size = self.COMM.Get_size()
        win = self.WIN
        try:
            win.Lock_all()
            win.Sync()
            win.Unlock_all()
        except NotImplementedError:
            return
        for rank in range(size):
            win.Lock(rank)
            win.Flush(rank)
            win.Unlock(rank)
        win.Lock_all()
        win.Flush_all()
        win.Unlock_all()
        for rank in range(size):
            win.Lock(rank)
            win.Flush_local(rank)
            win.Unlock(rank)
        win.Lock_all()
        win.Flush_local_all()
        win.Unlock_all()

class TestRMASelf(BaseTestRMA, unittest.TestCase):
    COMM = MPI.COMM_SELF

class TestRMAWorld(BaseTestRMA, unittest.TestCase):
    COMM = MPI.COMM_WORLD

try:
    MPI.Win.Create(None, 1, MPI.INFO_NULL, MPI.COMM_SELF).Free()
except NotImplementedError:
    del TestRMASelf, TestRMAWorld
else:
    name, version = MPI.get_vendor()
    if name == 'Open MPI':
        if version == (1,10,0):
            del BaseTestRMA.testFlushSync
        if version == (1,8,7):
            del BaseTestRMA.testStartCompletePostTest
            del BaseTestRMA.testStartCompletePostWait
        if version == (1,8,6):
            del BaseTestRMA.testPostWait
            del BaseTestRMA.testStartComplete
            del BaseTestRMA.testStartCompletePostTest
            del BaseTestRMA.testStartCompletePostWait
        if version == (1,8,1):
            del BaseTestRMA.testFenceAll
        if version < (1,4,0):
            if MPI.Query_thread() > MPI.THREAD_SINGLE:
                del TestRMAWorld
    if name == 'HP MPI':
        BaseTestRMA.COUNT_MIN = 1

if __name__ == '__main__':
    unittest.main()
