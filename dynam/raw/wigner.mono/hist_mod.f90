module hist_mod
  implicit none
  integer, parameter :: dpr = kind(1.d0)
  private
  public :: hist, hist2, hist2var
contains
  !*************************************************************************
  subroutine hist(xx,ndat,nh,x1,x2,name,iopw)
    implicit none
    real (kind=dpr), dimension (ndat), intent (inout) :: xx
    real (kind=dpr), intent (inout)  :: x1,x2
    integer, optional, intent (in) :: iopw
    integer, intent (in)    :: ndat
    integer, intent (inout) :: nh
    character (len=*), intent(in) :: name

    character (len=80) :: fname
    character (len=*), parameter :: form1='(1x,a1,f14.8,2es20.8)', &
         & form2='(f16.8,2es20.8)'
    integer, allocatable, dimension (:) :: nbox
    integer :: i,nlow,nhigh,ibox,iop
    real (kind=dpr)  :: dh,box,dat,h,frac

    if (present(iopw)) then
       iop = iopw
    else
       iop = 0
    endif

    if (nh < 0) then
       nh=20
       x1=minval(xx)
       x2=maxval(xx)
    endif

    if (iop /= 0 .and. nh > ndat) then
       write(6,"(a,3i6)") "  Error in HIST.   iop,nh,ndat =",iop,nh,ndat
       stop 12
    endif

    allocate (nbox(nh))
    nlow=0
    nhigh=0
    nbox=0
    dh=1.d0
    if (nh.ge.1) dh=(x2-x1)/real(nh,kind=dpr)
    do i=1,ndat
       box=(xx(i)-x1)/dh+0.999999999999_dpr
       ibox=box
       if (ibox.lt.1) then
          nlow=nlow+1
       elseif (ibox.gt.nh) then
          nhigh=nhigh+1
       else
          nbox(ibox)=nbox(ibox)+1
       endif
    enddo

    dat=real(ndat,kind=dpr)*dh
    if (iop == 0) then
       fname='histo.'//trim(name)
       open(8,file=fname,form='formatted',status='unknown')
       write (6,'(a,a,/,20x,a)') ' Histogram ',name,' number  fraction'
       write (6,form1) '<',x1,real(nlow)/dh,real(nlow)/dat
    endif
    do i=1,nh
       h=x1+dh*(real(i-1,kind=dpr)+0.5_dpr)
       frac=real(nbox(i),kind=dpr)/dat
       if (iop == 0) then
          write (6,form2) h,nbox(i)/dh,frac
          write (8,form2) h,nbox(i)/dh,frac
       else
          xx(i) = frac
       endif
    enddo

    if (iop == 0) then
       write (6,form1) '>',x2,real(nhigh)/dh,real(nhigh)/dat
       close(8)
    endif
    deallocate (nbox)
  end subroutine hist
  !*************************************************************************
  subroutine hist2(xx,yy,ndat,nh,x1,x2,name)
    implicit none
    real (kind=kind(1.d0)), dimension (ndat), intent (in) :: xx,yy
    real (kind=kind(1.d0)), intent (inout)  :: x1,x2
    integer, intent (in)    :: ndat
    integer, intent (inout) :: nh
    character (len=*), intent(in) :: name
    !
    character (len=80) :: fname
    character (len=*), parameter :: form1='(1x,a1,f20.8,i10,f20.8)', &
         & form2='(f16.8,e20.9,e20.10,2e15.6)'
    integer :: i,nlow,nhigh,ibox
    real (kind=kind(1.d0))  :: dh,box,dat,hdh,h,frac
    real (kind=kind(1.d0)), allocatable, dimension (:) :: sbox,av,sd,nb

    if (nh < 0) then
       nh=20
       x1=minval(xx)
       x2=maxval(xx)
    endif

    allocate (sbox(nh))
    allocate (nb(nh),av(nh),sd(nh))
    sbox=0.d0
    nlow=0
    nhigh=0
    dh=1.d0
    nb=0
    av=0.d0
    sd=0.d0
    if (nh.ge.1) dh=(x2-x1)/real(nh)
    do i=1,ndat
       box=(xx(i)-x1)/dh+0.999999999999d0
       ibox=box
       if (ibox.lt.1) then
          nlow=nlow+1
       elseif (ibox.gt.nh) then
          nhigh=nhigh+1
       else
          sbox(ibox)=sbox(ibox)+yy(i)
          nb(ibox)=nb(ibox)+1
          sd(ibox)=sd(ibox)+yy(i)**2
       endif
    enddo

    where (nb > 0)
       av=sbox/nb
       sd=sd/nb
       sd=sqrt(sd-av**2)
    end where

    sbox=sbox/dh

    fname='histo.'//trim(name)
    open(8,file=fname,form='formatted',status='unknown')
    write (6,'(a,a,/,20x,a,15x,a,20x,a,10x,a)')  &
         & ' Histogram ',name,' number','fraction','av','sd'
    dat=real(ndat)
    hdh=dh/2.d0
    write (6,form1) '<',x1,nlow,real(nlow)/dat
    do i=1,nh
       h=x1+dh*(real(i-1)+0.5d0)
       frac=sbox(i)/dat
       write (6,form2) h,sbox(i),frac,av(i),sd(i)
       write (8,form2) h,sbox(i),frac,av(i),sd(i)
    enddo
    write (6,form1) '>',x2,nhigh,real(nhigh)/dat
    !
    close(8)
    deallocate (sbox)
    deallocate (nb,av,sd)
  end subroutine hist2
  !*************************************************************************
  subroutine hist2var(xx,yy,ndat,nhx,nhy,x1,x2,y1,y2,name)
    implicit none
    real (kind=kind(1.d0)), dimension (ndat), intent (in) :: xx,yy
    real (kind=kind(1.d0)), intent (inout)  :: x1,x2,y1,y2
    integer, intent (in)    :: ndat
    integer, intent (inout) :: nhx,nhy
    character (len=*), intent(in) :: name
    !
    character (len=80) :: fname
    character (len=*), parameter :: form1='(1x,a1,2f20.8,i10,f20.8)', &
         & form2='(f16.8,e20.9,e20.10,2e15.6)'
    integer :: i,j,nlow,nhigh,ibox,xbox,ybox
    real (kind=kind(1.d0))  :: dhx,dhy,box,dat,hx,hy,frac
    real (kind=kind(1.d0)), allocatable, dimension (:,:) :: sbox

    if (nhx < 0) then
       nhx=20
       x1=minval(xx)
       x2=maxval(xx)
    endif
    if (nhy < 0) then
       nhy=20
       y1=minval(yy)
       y2=maxval(yy)
    endif

    allocate (sbox(nhx,nhy))
    sbox=0.d0
    nlow=0
    nhigh=0
    dhx=1.d0
    dhy=1.d0
    if (nhx.ge.1) dhx=(x2-x1)/real(nhx)
    if (nhy.ge.1) dhy=(y2-y1)/real(nhy)
    do i=1,ndat
       xbox=int((xx(i)-x1)/dhx+0.999999999999d0)
       ybox=int((yy(i)-y1)/dhy+0.999999999999d0)
       if (xbox.lt.1 .or. ybox.lt.1) then
          nlow=nlow+1
       elseif (xbox.gt.nhx .or. ybox.gt.nhy) then
          nhigh=nhigh+1
       else
          sbox(xbox,ybox)=sbox(xbox,ybox)+1
       endif
    enddo

    sbox=sbox/(dhx*dhy)

    fname='histo.'//trim(name)
    open(8,file=fname,form='formatted',status='unknown')
    write (6,'(a,3f15.9)')  'dhx,dhy=',dhx,dhy
    write (6,'(a,a,/,20x,a,15x,a)')  &
         & ' Histogram ',name,' number','fraction'
    dat=real(ndat)
    write (6,form1) '<',x1,y1,nlow,real(nlow)/dat
    do i=1,nhx
       hx=x1+dhx*(real(i-1)+0.5d0)
       do j=1,nhy
          hy=y1+dhy*(real(j-1)+0.5d0)
          frac=sbox(i,j)/dat
          write (6,form2) hx,hy,sbox(i,j),frac
          write (8,form2) hx,hy,sbox(i,j),frac
       end do
       write (8,*)
    enddo
    write (6,form1) '>',x2,y2,nhigh,real(nhigh)/dat
    close(8)
    deallocate (sbox)
  end subroutine hist2var
end module hist_mod
