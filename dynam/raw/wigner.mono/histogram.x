#!/bin/csh

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=50,filein='Stat.pun.trans_wigner',
     fileout='Epot',nsample=10010,nskip=1,ncol=8,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=50,filein='Stat.pun.trans_wigner',
     fileout='Ecin',nsample=10010,nskip=1,ncol=7,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=50,filein='Stat.pun.trans_wigner',
     fileout='Etot',nsample=10010,nskip=1,ncol=6,fac=27.2113957
&END
EOF

# raccolta energie armoniche
grep 'harmonic pot. e.' trans_wigner.out | awk '{print $9}' >harm_ene.out
/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=5.0,nh=50,filein='harm_ene.out',
     fileout='Eharm',nsample=10010,nskip=0,ncol=1,fac=27.2113957
&END
EOF

# raccolta energie iniziali
grep '      1         0.0000000000' Stat.pun.S1 >Stat.pun.ini 
/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=8.0,nh=40,filein='Stat.pun.ini',
     fileout='E_S0',nsample=1000,nskip=0,ncol(1)=7,ncol(2)=18,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=0.0,xmax=8.0,nh=40,filein='Stat.pun.ini',
     fileout='E_S1',nsample=1000,nskip=0,ncol=6,fac=27.2113957
&END
EOF

/home/mau/azometano/dynam/histogram << EOF >histogram.out
&DAT xmin=3.4,xmax=3.6,nh=40,filein='Stat.pun.ini',
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
