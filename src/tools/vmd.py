import subprocess
import time
import socket
import select
import struct
import espresso

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
  (sock, sock_port) = initsock.accept()
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

def connect(system):
  espresso.tools.psfwrite("vmd.psf", system)
  espresso.tools.pdbwrite("vmd.pdb", system)
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
  vmdfile.write("mol load psf vmd.psf pdb vmd.pdb\n")
  vmdfile.write("logfile vmd.log\n")
  vmdfile.write("rotate stop\n")
  vmdfile.write("logfile off\n")
  vmdfile.write("mol modstyle 0 0 CPK 1.800000 0.300000 8.000000 6.000000\n")
  vmdfile.write("mol modcolor 0 0 SegName\n")
  st = "imd connect %s %i\n" % (hostname, port)
  vmdfile.write(st)
  vmdfile.write("imd transfer 1\n")
  vmdfile.write("imd keep 1\n")
  vmdfile.close()
  subprocess.Popen(['/sw/linux/vmd-1.8.7/bin/vmd', '-e', 'vmd.tcl'])

  sock = handshake(initsock)
  if (sock != 0):
      time.sleep(0.25)
      drain_socket(sock)
  
  return sock

def imd_positions(system, sock):
  maxParticleID = int(espresso.analysis.MaxPID(system).compute())
  count  = 0
  pid    = 0
  coords = struct.pack('')
  while pid <= maxParticleID:
    particle = system.storage.getParticle(pid)
    if particle.pos:
        coords += struct.pack('!fff', particle.pos[0], particle.pos[1], particle.pos[2])
        count  += 1
        pid    += 1
    else:
        pid    += 1
  header = struct.pack('!II',IMD_FCOORDS,count)
  msg    = header + coords
  sock.send(msg)
