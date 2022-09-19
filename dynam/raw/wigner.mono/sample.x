#!/bin/csh

/home/mau/mopac2002/util/sample << EOF >sample.out
&sam filegeo='wigner.mc_geo',
     filedyn='trans_wigner.dyn', fileinf='trans_wigner.inf',
     opsam='BIS',ngeo=8500,prob=T,refprob=0.01,hv=3.53,dhv=0.05, &end
EOF

exit
  
namelist /sam/ filegeo,filedyn,fileinf, &
  & opsam,ngeo,iseed,prob,refprob,hv,dhv,ngw_min,istat
  !
  !     ********  DATI DI INPUT *********
  !     filedyn=file .dyn di mopac
  !     fileinf=file .inf di mopac
  !     filegeo=file di geometrie sampled che produrra il programma
  !     opsam=opzioni di sampling (REG=regular, RAN=random, BIS=Bisezione)
  !     ngeo=numero di geometrie che si desiderano
  !     iseed=random number seed
  !     prob=comando logico per scegliere di includere o no la probabilita` di transizione
  !     refprob=valore di riferimento per la probabilita` di transizione
  !     hv=massimo di energia (in eV)
  !     dhv=larghezza dell'intervallo di energia
  !     ngw_min=per ogni traiettoria, minimo valore di ngw dal quale si inizia a campionare
  !     istat=se prob=F, istat e' lo stato per il quale si calcola la
  !           differenza di energia col fondamantale
  !
