"""
Python module oldespresso
=========================

This modules provides a reader for (old) ESPResSo input blockfiles.
"""

class Reader:
   """
   Class that reads a blockfile used for ESPResSo configuration file.

   This reader does nothing else than reading the blocks and representing
   them in an equivalent Python list structure. 
   """

   def __init__(self, filename):

       self.regionNest = []
       self.list = []

       file = open(filename)

       linenr = 0    # for error messages

       for line in file:

          linenr = linenr + 1

          items = line.split()

          for item in items:

              if item.startswith("{"):

                 region = item[1:]

                 self.enter(region)

              elif item.endswith("}"):

                 item = item.rstrip("}")

                 self.addItem(item)

                 self.leave()

              else :

                 self.addItem(item)

   def enter(self, region):

       endFlag = region.endswith("}")

       if endFlag: region = region.rstrip("}")

       # save current list on stack

       self.regionNest.append(self.list)

       self.list = []

       self.addItem(region)

       if endFlag: leave(region)

   def leave(self):

       last = self.regionNest.pop()

       last.append(self.list)

       self.list = last

   def addItem(self, item):

       if item == "": return

       self.list.append(item)

   def getValue(self, key, entries):

       if not hasattr(list, "__iter__"): raise "%s not iterable"%entries

       for entry in entries:

           if entry[0] == key: return entry[1:]

       raise "entry with key %s not found"%key

   def __getitem__(self, index):

      if isinstance(index, str):
 
         return self.getValue(index, self.list)

      else:

         result = self.list

         for key in index:
             result = self.getValue(key, result)

         return result

