#!/bin/csh

set filepun = \
   '/home/carlotta/azometano/dynam/thermal/400traj_boltz/Stat.pun.trans_equi'
set fileout = \
   '/home/carlotta/azometano/dynam/thermal/400traj_boltz/trans_equi.out'
set files1 = \
   '/home/carlotta/azometano/dynam/thermal/400traj_boltz/Stat.pun.S1'

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=50,filein="$filepun",
     fileout='Epot',nsample=10010,nskip=1,ncol=8,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=50,filein="$filepun",
     fileout='Ecin',nsample=10010,nskip=1,ncol=7,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=50,filein="$filepun",
     fileout='Etot',nsample=10010,nskip=1,ncol=6,fac=27.2113957
&END
EOF

# raccolta energie armoniche
grep 'harmonic pot. e.' $fileout | awk '{print $9}' >harm_ene.out
/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=50,filein='harm_ene.out',
     fileout='Eharm',nsample=10010,nskip=0,ncol=1,fac=27.2113957
&END
EOF

# raccolta energie iniziali
grep '      1         0.0000000000' $files1 >Stat.pun.ini 
/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=9.0,nh=30,filein='Stat.pun.ini',
     fileout='E_S0',nsample=1000,nskip=0,ncol(1)=7,ncol(2)=18,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=9.0,nh=30,filein='Stat.pun.ini',
     fileout='E_S1',nsample=1000,nskip=0,ncol=6,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=2.5,xmax=4.5,nh=20,filein='Stat.pun.ini',
     fileout='E_ecc',nsample=1000,nskip=0,ncol(1)=19,ncol(2)=-18,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=1.5,nh=10,filein='azmkine.pun',
     fileout='Ek_az',nsample=1000,nskip=1,ncol(1)=3,fac=1.0,
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=1.5,nh=10,filein='azmkine.pun',
     fileout='Ek_met1',nsample=1000,nskip=1,ncol(1)=4,fac=1.0,
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=1.5,nh=10,filein='azmkine.pun',
     fileout='Ek_met2',nsample=1000,nskip=1,ncol(1)=5,fac=1.0,
&END
EOF

paste histo.Ek_met1 histo.Ek_met2 >histo.Ek_met
