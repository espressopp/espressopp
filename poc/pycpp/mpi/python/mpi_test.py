import cpp.mpi_test

def test(code):
    'Test some MPI'
    print "handover mpi_test -> cpp.mpi_test"
    cpp.mpi_test.test(code)
