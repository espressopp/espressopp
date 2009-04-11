from espresso import pmi

from _espresso import hello_HelloWorld as HelloWorldLocal

if pmi.IS_CONTROLLER:
    pmi.exec_('from espresso.hello import HelloWorldLocal')
    pmi.exec_('joinstr=lambda a,b: \"\\n\".join((a,b))')
    class HelloWorld(object):
        def __init__(self) :
            self.local = pmi.create('HelloWorldLocal')
            return object.__init__(self)
            
        def __str__(self) :
            'Returns the messages.'
            return pmi.reduce('joinstr',
                              'HelloWorldLocal.getMessage',
                              self.local)
