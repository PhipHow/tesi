#!/bin/csh

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=100,filein='Stat.pun.trans_equi',
     fileout='Epot',nsample=10010,nskip=1,ncol=8,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=100,filein='Stat.pun.trans_equi',
     fileout='Ecin',nsample=10010,nskip=1,ncol=7,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=100,filein='Stat.pun.trans_equi',
     fileout='Etot',nsample=10010,nskip=1,ncol=6,fac=27.2113957
&END
EOF

# raccolta energie iniziali
grep '      1         0.0000000000' Stat.pun.S1 >Stat.pun.ini 
/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=8.0,nh=100,filein='Stat.pun.ini',
     fileout='E_S0',nsample=1000,nskip=0,ncol(1)=7,ncol(2)=18,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=8.0,nh=100,filein='Stat.pun.ini',
     fileout='E_S1',nsample=1000,nskip=0,ncol=6,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=2.5,xmax=4.5,nh=20,filein='Stat.pun.ini',
     fileout='E_ecc',nsample=1000,nskip=0,ncol(1)=19,ncol(2)=-18,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=1.6,nh=8,filein='azmkine.pun',
     fileout='Ek_az',nsample=1000,nskip=1,ncol(1)=3,fac=1.0,
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=1.6,nh=8,filein='azmkine.pun',
     fileout='Ek_met1',nsample=1000,nskip=1,ncol(1)=4,fac=1.0,
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=1.6,nh=8,filein='azmkine.pun',
     fileout='Ek_met2',nsample=1000,nskip=1,ncol(1)=5,fac=1.0,
&END
EOF

paste histo.Ek_met1 histo.Ek_met2 >histo.Ek_met
