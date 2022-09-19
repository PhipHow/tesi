program DISS

implicit none
integer:: i,nflo, ntraj, ir, idum, nold,ntrak,istep,ndiss1,ndiss2,nodiss,nbad,imin
double precision:: E0, E, tol,r1,r2,endtime,rdiss,time,BohrtoA,rmax,angmax,rotmax
double precision, dimension(2):: avnnc1,avnnc2,avcnnc
double precision, allocatable,dimension(:):: flo
character*64:: filepun
character*12, allocatable,dimension(:):: varname  
character*12:: var1,var2
logical:: zdiss1,zdiss2
integer, dimension(2):: ntrans,ncis,nrot,ninv,nang
namelist/DAT/nflo,ir,filepun,rdiss,rmax,angmax,rotmax

BohrtoA=0.52917726
rdiss=4.D0
rmax=2.D0
angmax=150.D0
rotmax=60.D0

read(5,nml=DAT) 
allocate(varname(nflo),flo(nflo))
open(10,file=filepun,form='formatted',status='old')
open(7,file='diss.pun',form='formatted',status='unknown')
write(7,*)'endtime, nodiss, diss1, diss2, trans-S0, trans-S1, cis-S0, cis-S1, rot-S0, rot-S1, inv-S0, inv-S1, &
        <nnc1>, <nnc2>, <cnnc>, bad'

!!!!!!!!!!!!!!!!!!Ciclo sui tempi finali!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do
read(5,*,end=9999,err=9999)endtime
rewind(10)
read(10,*) var1,var2,varname
ndiss1=0 
ndiss2=0 
nodiss=0 
nbad=0
nold=0
ntrans=0
ncis=0
nrot=0
ninv=0
nang=0
avnnc1=0.D0
avnnc2=0.D0
avcnnc=0.D0
zdiss1=.false.
zdiss2=.false.
 do
  read(10,*,end=999)ntraj,idum,flo
  if (ntraj /= nold) then
     if(time < endtime) then
       if(zdiss2) then
         ndiss2=ndiss2+1
       else
         nbad=nbad+1
       endif
     endif
     nold=ntraj
     zdiss1=.false.
     zdiss2=.false.
  else
     time=flo(1)
    r1=999.D0
    do i=0,3
     if(flo(ir+i) < r1) then
       r1=flo(ir+i)
       imin=i 
     endif
    enddo
    if(imin == 0) then
     ! r1=flo(ir)
      r2=flo(ir+1)
    elseif (imin == 1) then
     ! r1=flo(ir+1)
      r2=flo(ir)
    elseif (imin == 2) then
     ! r1=flo(ir+2)
      r2=flo(ir+3)
    elseif (imin == 3) then
     ! r1=flo(ir+3)
      r2=flo(ir+2)
    endif
    r1=r1*BohrtoA
    r2=r2*BohrtoA

    if(r2 > rdiss) then
      if(r1 > rdiss) then
        zdiss2=.true.
        zdiss1=.false.
      else 
        zdiss2=.false.
        zdiss1=.true.
      endif
    endif

    if(abs(time-endtime) < 1.E-8) then
      if(zdiss2) ndiss2=ndiss2+1 
      if(zdiss1) ndiss1=ndiss1+1
      if((.not. zdiss1) .and. (.not. zdiss2)) nodiss=nodiss+1
      if(r2 < rmax) then
        nang(nint(flo(2)))=nang(nint(flo(2)))+1
        avcnnc(nint(flo(2)))=avcnnc(nint(flo(2)))+flo(ir+6)
        if(flo(ir+4) > flo(ir+5)) then
          avnnc1(nint(flo(2)))=avnnc1(nint(flo(2)))+flo(ir+5)   !!!angolo piccolo
          avnnc2(nint(flo(2)))=avnnc2(nint(flo(2)))+flo(ir+4)   !!!angolo grande
        else
          avnnc1(nint(flo(2)))=avnnc1(nint(flo(2)))+flo(ir+4)   !!!angolo piccolo
          avnnc2(nint(flo(2)))=avnnc2(nint(flo(2)))+flo(ir+5)   !!!angolo grande
        endif
        if((flo(ir+4) > angmax) .or. (flo(ir+5) > angmax)) then
          ninv(nint(flo(2)))=ninv(nint(flo(2)))+1
        else
          if(flo(ir+6) > (180.D0-rotmax)) then
            ntrans(nint(flo(2)))=ntrans(nint(flo(2)))+1
          elseif(flo(ir+6) < rotmax) then
            ncis(nint(flo(2)))=ncis(nint(flo(2)))+1
          else
            nrot(nint(flo(2)))=nrot(nint(flo(2)))+1
          endif
        endif
      endif        
    endif 
 endif
enddo

999 if(time < endtime) nbad=nbad+1
write(6,*) 'tempo finale=',endtime
write(6,*) 'numero di traiettorie=',ntraj
ntraj=ntraj-nbad
write(6,*) 'numero di traiettorie indissociate=',nodiss,dble(100*nodiss)/dble(ntraj),' % '
write(6,*) 'numero di singole dissociazioni=',ndiss1,dble(100*ndiss1)/dble(ntraj),' % '
write(6,*) 'numero di doppie dissociazioni=',ndiss2,dble(100*ndiss2)/dble(ntraj),' % '
write(6,*) 'numero di bad traiettorie=',nbad
write(6,*) 'numero di molecole trans=',ntrans!,dble(ntrans)/dble(ntraj),' % '
write(6,*) 'numero di molecole cis=',ncis!,dble(ncis)/dble(ntraj),' % '
write(6,*) 'numero di ts rotazione=',nrot!,dble(nrot)/dble(ntraj),' % '
write(6,*) 'numero di ts inversione=',ninv!,dble(ninv)/dble(ntraj),' % '
if(nang(1) == 0) nang(1)=1
if(nang(2) == 0) nang(2)=1
write(7,'(f10.1,11(I4,f6.3),6f7.2,I6)') endtime,nodiss,dble(nodiss)/dble(ntraj),&
       ndiss1,dble(ndiss1)/dble(ntraj),ndiss2,dble(ndiss2)/dble(ntraj), &
       ntrans(1),dble(ntrans(1))/dble(ntraj),ntrans(2),dble(ntrans(2))/dble(ntraj), &
       ncis(1),dble(ncis(1))/dble(ntraj),ncis(2),dble(ncis(2))/dble(ntraj), &
       nrot(1),dble(nrot(1))/dble(ntraj),nrot(2),dble(nrot(2))/dble(ntraj), &
       ninv(1),dble(ninv(1))/dble(ntraj),ninv(2),dble(ninv(2))/dble(ntraj), &
       avnnc1(1)/dble(nang(1)),avnnc1(2)/dble(nang(2)),avnnc2(1)/dble(nang(1)),avnnc2(2)/dble(nang(2)), &
       avcnnc(1)/dble(nang(1)),avcnnc(2)/dble(nang(2)),nbad
enddo
!!!!!!!!!!!!!!!!!!!!Fine ciclo tempi finali!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
9999 stop

end program DISS
