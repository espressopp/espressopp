from _espresso import hello_HelloWorld as _HelloWorld

# map the class _HelloWorld to the class in the pseudo-namespace
#_HelloWorld = __import__(_espresso.hello_HelloWorld)

# create the wrapping python class around the C++ class
class HelloWorld(_HelloWorld):
    'Python wrapper from C++ class HelloWorld'
    def getMessage(self) :
	'Returns the messages.' 
        return _HelloWorld.getMessage(self)
