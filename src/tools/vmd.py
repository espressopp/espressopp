#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


import subprocess
import time
import socket
import select
import struct
import espressopp

IMD_HANDSHAKE = 4
IMDVERSION    = 2

IMD_DISCONNECT, \
IMD_ENERGIES,   \
IMD_FCOORDS,    \
IMD_GO,         \
IMD_HANDSHAKE,  \
IMD_KILL,       \
IMD_MDCOMM,     \
IMD_PAUSE,      \
IMD_TRATE,      \
IMD_IOERROR = range(10)

def handshake(initsock):
  print "before"
  (sock, sock_port) = initsock.accept()
  print "after"
  header = struct.Struct('!II')
  msg    = header.pack(IMD_HANDSHAKE, IMDVERSION)
  sock.send(msg)
  res = []
  cnt = 100
  while len(res) == 0 and cnt > 0:
    res = select.select([sock],[],[],0)[0]
    cnt -= 1
    time.sleep(0.1)
  if len(res) == 1:
    msg = struct.unpack('!II',sock.recv(8))
    if msg[0] == IMD_GO:
      print "VMD sent IMD_GO"
      return sock
    else:
      print "unexpected answer from VMD"
      return 0
  else:
    print "VMD did not answer."
    return 0

def drain_socket(sock):
  res = select.select([sock],[],[],0)[0]
  while len(res) > 0:
    buf  = sock.recv(8)
    bufl = len(buf)
    if bufl == 0:
      break
    # msg = struct.unpack('B'*bufl,buf)
    res = select.select([sock],[],[],0)[0]
  return

def connect(system, molsize=10, pqrfile=False, vmd_path='vmd'):
  """Connects to the VMD.

    :param espressopp.system system: The system object.
    :param int molsize: The optional size of the molecule.
    :param bool pqrfile: If set to True then the pqr vmd.pqr file will be used otherwise (default)
      the vmd.pdb file will be used.
    :param str vmd_path: The path to the executable of vmd, by default it is set to 'vmd'.

    :returns: Socket to the VMD.

  """
  espressopp.tools.psfwrite("vmd.psf", system, molsize=molsize, maxdist=system.bc.boxL[0]/2)
  if pqrfile==True:
    espressopp.tools.pqrwrite("vmd.pqr", system, molsize=molsize)
  else:
    espressopp.tools.pdbwrite("vmd.pdb", system, molsize=molsize)
  initsock = socket.socket(socket.AF_INET, socket.SOCK_STREAM, 0)
  hostname = socket.gethostname()
  port     = 10000
  while port < 65000:
    try:
      initsock.bind((hostname, port))
      initsock.listen(1)
      break
    except:
      port += 1
  if port == 65000:
    print "no free port for vmd socket found."
    return initsock

  vmdfile = open("vmd.tcl","w")
  if pqrfile == True:
    vmdfile.write("mol load psf vmd.psf pqr vmd.pqr\n")
  else:
    vmdfile.write("mol load psf vmd.psf pdb vmd.pdb\n")
  vmdfile.write("logfile vmd.log\n")
  vmdfile.write("rotate stop\n")
  vmdfile.write("logfile off\n")
  # vmdfile.write("mol modstyle 0 0 CPK 0.500000 0.500000 8.000000 6.000000\n")
  vmdfile.write("mol modstyle 0 0 VDW 0.4 20\n")
  vmdfile.write("mol modcolor 0 0 SegName\n")
  vmdfile.write("color Segname {T000} 3\n")

  #vmdfile.write("mol delrep 0 top\n")
  #vmdfile.write("mol representation CPK 0.500000 0.500000 8.000000 6.000000\n")
  #vmdfile.write("mol color SegName\n")
  #vmdfile.write("mol selection {segname T000}\n")
  #vmdfile.write("mol material Opaque\n")
  #vmdfile.write("mol addrep top\n")
  #vmdfile.write("mol selupdate 0 top 0\n")
  #vmdfile.write("mol colupdate 0 top 0\n")
  #vmdfile.write("mol scaleminmax top 0 0.000000 0.000000\n")
  #vmdfile.write("mol smoothrep top 0 0\n")
  #vmdfile.write("mol drawframes top 0 {now}\n")
  #vmdfile.write("color Segname {T000} red\n")
  #vmdfile.write("color Display {Background} silver\n")

  st = "imd connect %s %i\n" % (hostname, port)
  vmdfile.write(st)
  vmdfile.write("imd transfer 1\n")
  vmdfile.write("imd keep 1\n")
  vmdfile.close()
  subprocess.Popen([vmd_path, '-e', 'vmd.tcl'])

  sock = handshake(initsock)
  if (sock != 0):
      time.sleep(0.25)
      drain_socket(sock)

  return sock

def imd_positions(system, sock, folded=True):
  maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
  count  = 0
  pid    = 0
  coords = struct.pack('')
  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
        particle = system.storage.getParticle(pid)
        if folded:
          p = particle.pos
        else:
          p = system.bc.getUnfoldedPosition(particle.pos, particle.imageBox)
        coords += struct.pack('!fff', p[0], p[1], p[2])
        count  += 1
        pid    += 1
    else:
        pid    += 1
  header = struct.pack('!II',IMD_FCOORDS,count)
  msg    = header + coords
  sock.send(msg)
