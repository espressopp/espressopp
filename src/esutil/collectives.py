import _espresso
from espresso import pmi

ResultNone = _espresso.esutil_Collectives_ResultNone

def locateItem(here):
    """locate the node with here=True (e.g. indicating that data of a
    distributed storage is on the local node). This is a collective
    SPMD function.

    here is a boolean value, which should be True on at most one
    node. Returns on the controller the number of the node with
    here=True, or an KeyError exception if no node had the item,
    i.e. had here=True.
    """
    res = _espresso.esutil_Collectives_locateItem(here, pmi.CONTROLLER)
    if pmi.IS_CONTROLLER:
        if res == ResultNone:
            raise IndexError("collectives.locateItem could not find anything")
        return res
