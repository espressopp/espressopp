mol load psf vmd.psf pdb vmd.pdb
logfile vmd.log
rotate stop
logfile off
mol modstyle 0 0 VDW 0.4 20
mol modcolor 0 0 SegName
color Segname {T000} 3
imd connect pckr146 10000
imd transfer 1
imd keep 1
