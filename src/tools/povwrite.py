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

"""
*****************************
povwrite - write povray files
*****************************
"""

import espressopp
from espressopp import Real3D

def camera():
  camera = "camera {\n\
    up <0, 6.0000, 0>\n\
    right <5.8417, 0, 0>\n\
    location <0.0000, 0.0000, -2.0000>\n\
    look_at <0.0000, 0.0000, -0.0000>\n\
    direction <-0.0000, -0.0000, 4.0000>\n\
  }\n"

  return camera

def lightsource1():
  lightsource1 = "light_source { \n\
    <-0.1000, 0.1000, -1.0000> \n\
    color rgb<1.000, 1.000, 1.000> \n\
    parallel \n\
    point_at <0.0, 0.0, 0.0> \n\
  }\n"


  return lightsource1

def lightsource2():
  lightsource2 = "light_source { \n\
    <1.0000, 2.0000, -0.5000>  \n\
    color rgb<1.000, 1.000, 1.000> \n\
    parallel \n\
    point_at <0.0, 0.0, 0.0> \n\
  }\n"

  return lightsource2

def background():
  background = "background { \n\
    color rgb<0.000, 0.000, 0.000> \n\
  }\n\
  fog {\n\
    distance 3.1250 \n \
    fog_type 1 \n\
    color rgb<0.000, 0.000, 0.000> \n\
  } \n"

  return background



def povwrite(system, integrator, filename, append=False, box=True):

  red = 1.0
  green = 0.5
  blue = 0.0
  framered = 1.0
  framegreen = 1.0
  frameblue = 1.0
  transparent = 0.0   

  boxsizex = system.bc.boxL[0]/2  
  boxsizey = system.bc.boxL[1]/2
  boxsizez = system.bc.boxL[2]/2
  maxParticleID = int(espressopp.analysis.MaxPID(system).compute())

  f = 0.7 #size
  pid = 0
  r = []
  xpos = []
  ypos = []
  zpos = []

  if append:
    file = open(filename, 'a')
  else:    
    file = open(filename,'w')

  if boxsizex>=boxsizey and boxsizex>=boxsizey:
    g = boxsizex

  if boxsizey>=boxsizex and boxsizey>=boxsizez:
    g = boxsizey

  if boxsizez>=boxsizex and boxsizez>=boxsizey:
    g = boxsizez

  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
      particle = system.storage.getParticle(pid)

      xpos.append(particle.pos[0])
      ypos.append(particle.pos[1])
      zpos.append(particle.pos[2])
      r.append(particle.radius)

      pid += 1

    else:
      pid += 1

  head = "// \n\
// Molecular graphics export from VMD 1.9 \n\
// http://www.ks.uiuc.edu/Research/vmd/ \n\
// Requires POV-Ray 3.5 or later \n\
// \n\
// POV 3.x input script : overlay.pov \n\
// try povray +W812 +H834 -Ioverlay.pov -Ooverlay.pov.tga +P +X +A +FT +C \n\
#if (version < 3.5) \n\
#error \"VMD POV3DisplayDevice has been compiled for POV-Ray 3.5 or above. Please upgrade POV-Ray or recompile VMD. \" \n\
#end \n\
#declare VMD_clip_on=array[3] {0, 0, 0};\n\
#declare VMD_clip=array[3];\n\
#declare VMD_scaledclip=array[3];\n\
#declare VMD_line_width=0.0020;\n\
#macro VMDC ( C1 )\n\
  texture { pigment { rgbt C1 }}\n\
#end\n\
#macro VMD_point (P1, R1, C1)\n\
  #local T = texture { finish { ambient 1.0 diffuse 0.0 phong 0.0 specular 0.0 } pigment { C1 } }\n\
  #if(VMD_clip_on[2])\n\
  intersection {\n\
    sphere {P1, R1 texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
    VMD_clip[2]\n\
  }\n\
  #else\n\
  sphere {P1, R1 texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
  #end\n\
#end\n\
#macro VMD_line (P1, P2, C1)\n\
  #local T = texture { finish { ambient 1.0 diffuse 0.0 phong 0.0 specular 0.0 } pigment { C1 } }\n\
  #if(VMD_clip_on[2])\n\
  intersection {\n\
    cylinder {P1, P2, VMD_line_width texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
    VMD_clip[2]\n\
  }\n\
  #else\n\
  cylinder {P1, P2, VMD_line_width texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
  #end\n\
#end\n\
#macro VMD_sphere (P1, R1, C1)\n\
  #local T = texture { pigment { C1 } }\n\
  #if(VMD_clip_on[2])\n\
  intersection {\n\
    sphere {P1, R1 texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
    VMD_clip[2]\n\
  }\n\
  #else\n\
  sphere {P1, R1 texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
  #end\n\
#end\n\
#macro VMD_cylinder (P1, P2, R1, C1, O1)\n\
  #local T = texture { pigment { C1 } }\n\
  #if(VMD_clip_on[2])\n\
  intersection {\n\
    cylinder {P1, P2, R1 #if(O1) open #end texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
    VMD_clip[2]\n\
  }\n\
  #else\n\
  cylinder {P1, P2, R1 #if(O1) open #end texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
  #end\n\
#end\n\
#macro VMD_cone (P1, P2, R1, C1)\n\
  #local T = texture { pigment { C1 } }\n\
  #if(VMD_clip_on[2])\n\
  intersection {\n\
    cone {P1, R1, P2, VMD_line_width texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
    VMD_clip[2]\n\
  }\n\
  #else\n\
  cone {P1, R1, P2, VMD_line_width texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
  #end\n\
#end\n\
#macro VMD_triangle (P1, P2, P3, N1, N2, N3, C1)\n\
  #local T = texture { pigment { C1 } }\n\
  smooth_triangle {P1, N1, P2, N2, P3, N3 texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
#end\n\
#macro VMD_tricolor (P1, P2, P3, N1, N2, N3, C1, C2, C3)\n\
  #local NX = P2-P1;\n\
  #local NY = P3-P1;\n\
  #local NZ = vcross(NX, NY);\n\
  #local T = texture { pigment {\n\
    average pigment_map {\n\
      [1 gradient x color_map {[0 rgb 0] [1 C2*3]}]\n\
      [1 gradient y color_map {[0 rgb 0] [1 C3*3]}]\n\
      [1 gradient z color_map {[0 rgb 0] [1 C1*3]}]\n\
    }\n\
    matrix <1.01,0,1,0,1.01,1,0,0,1,-.002,-.002,-1>\n\
    matrix <NX.x,NX.y,NX.z,NY.x,NY.y,NY.z,NZ.x,NZ.y,NZ.z,P1.x,P1.y,P1.z>\n\
  } }\n\
  smooth_triangle {P1, N1, P2, N2, P3, N3 texture {T} #if(VMD_clip_on[1]) clipped_by {VMD_clip[1]} #end no_shadow}\n\
#end\n "

  file.write(head)
  file.write(camera())
  file.write(lightsource1())
  file.write(lightsource2())
  file.write(background())

  file.write("#default { texture {  \n\
 finish { ambient 0.000 diffuse 0.650 phong 0.1 phong_size 40.000 specular 0.500 }\n\
} }\n")

  file.write("#declare VMD_line_width=0.0020;\n\
#declare VMD_line_width=0.0060; \n")

  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*f,boxsizey/g*f,boxsizez/g*-f,boxsizex/g*f,boxsizey/g*-f,boxsizez/g*-f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*f,boxsizey/g*-f,boxsizez/g*-f,boxsizex/g*-f,boxsizey/g*-f,boxsizez/g*-f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*-f,boxsizey/g*-f,boxsizez/g*-f,boxsizex/g*-f,boxsizey/g*f,boxsizez/g*-f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*-f,boxsizey/g*f,boxsizez/g*-f,boxsizex/g*f, boxsizey/g*f,boxsizez/g*-f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*f,boxsizey/g*f,boxsizez/g*-f,boxsizex/g*f,boxsizey/g* f,boxsizez/g*f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*f,boxsizey/g*-f,boxsizez/g*-f,boxsizex/g*f, boxsizey/g*-f,boxsizez/g*f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*-f,boxsizey/g*-f,boxsizez/g*-f,boxsizex/g*-f,boxsizey/g*-f,boxsizez/g*f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*-f,boxsizey/g*f,boxsizez/g*-f,boxsizex/g*-f,boxsizey/g*f,boxsizez/g*f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*f,boxsizey/g*f,boxsizez/g*f,boxsizex/g*f,boxsizey/g*-f,boxsizez/g*f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*f,boxsizey/g*-f,boxsizez/g*f,boxsizex/g*-f,boxsizey/g* -f,boxsizez/g*f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*-f,boxsizey/g*-f,boxsizez/g*f,boxsizex/g*-f,boxsizey/g*f,boxsizez/g*f, framered, framegreen, frameblue, transparent))
  file.write("VMD_line(<%f,%f,%f>,<%f,%f,%f>,rgbt<%f,%f,%f,%f>)\n"%(boxsizex/g*-f,boxsizey/g*f,boxsizez/g*f,boxsizex/g*f,boxsizey/g* f,boxsizez/g*f, framered, framegreen, frameblue, transparent))


  file.write("// MoleculeID: 0 ReprID: 0 Beginning VDW\n")
  pid = 0
  while pid < maxParticleID:
    x = (boxsizex-xpos[pid])/g*f
    y = (boxsizey-ypos[pid])/g*f
    z = (boxsizez-zpos[pid])/g*f
    radius = r[pid]/2/g*f
    file.write("VMD_sphere(<%f,%f,%f>,%f,rgbt<%f,%f,%f>)\n"%(x,y,z,radius, red, green, blue)) 
    pid += 1
