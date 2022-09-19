#!/bin/csh

set photopun='/home/carlotta/azometano/dynam/400traj_wigner/Stat.pun.photo.tot'

/home/mau/azometano/dynam/azmtraj << EOF > azmtraj.out
&DAT nflo=14, ntraj=500, ntime=5000, dtime=2.0, ir=8, ie=4,
     rmax=2.0, rdiss=3.2, emin=0.0, emax=99.0, badend=T,
     ratd=2.0, nhist=20, filepun="$photopun" &END
EOF

exit
namelist/DAT/ nflo, ir, ie, filepun, rdiss, rmax, angmax, rotmax, &
   ntraj, ntime, dtime, emin, emax
