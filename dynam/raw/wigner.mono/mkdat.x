#!/bin/csh

set nfile = 20
set nlast = 13
foreach file (01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20)
  cp photo.dat photo$file.dat
  sed -i s/geo01/geo01_$file/ photo$file.dat
  if ($file == $nfile) sed -i s/NTRAJ=20/NTRAJ=$nlast/ photo$file.dat
# /home/gio/mopac2002/mopac2002.x photo$file >& mopac$file.log & 
end
