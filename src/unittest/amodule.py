# Test module to test PMI imports

class AOS :
    pass

class A(object) :
    def __init__(self, arg=None) :
        object.__init__(self)
        self.arg = arg
        print('Got arg: ' + str(arg))
