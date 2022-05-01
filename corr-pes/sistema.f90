program sistema 

implicit none 

integer :: i,j,n 
integer, parameter :: dpr = kind(1.d0) 
real (kind = dpr) :: r0, ra, theta0, r1, r2, ang1, ang2, cnnc, det, &
                     k1, v1, v2, pi, torad, arg1, carg1, carg2, zp, ang1or
real (kind = dpr), allocatable, dimension(:) :: zpe, par, w 
real (kind = dpr), allocatable, dimension(:,:) :: a 
integer, allocatable, dimension(:) :: iw 

character*80 :: string
 

DATA PI,TORAD/ & 
3.141592653589793D0,0.01745329251994329D0/

 
r0 = 1.6 
ra = 3.5
theta0 = 126.8176185 

read(5,*) n 
allocate(zpe(1:n), par(1:n), a(1:n,1:n), w(1:n), iw(1:n)) 

write(6,*) "MATRICE A" 

do i = 1,n 
    
    read(5,*) r1, r2, ang1, ang2, cnnc, zpe(i) 
    
    !calcolo dei coefficienti di A1,A2,A3 
    call quartic(r1,v1)
    call quartic(r2,v2) 

    a(i,1) = v1 + v2 
    a(i,2) = v1*v2 
    a(i,3) = v1**2 + v2**2 

 
    !calcolo dei coefficienti di Ainv1,Ainv2,Ainv3 
    a(i,4) = ((1.0_dpr-v1)*cos(torad*ang1))+((1.0_dpr-v2)*cos(torad*ang2)) 
    a(i,5) = ((1.0_dpr-v1)*cos(2.0_dpr*torad*ang1))+((1.0_dpr-v2)*cos(2.0_dpr*torad*ang2)) 

    !calcolo dei coefficienti di Arot1,Arot2,Arot3,Arot4
    call cubic(ang1,carg1)
    call cubic(ang2,carg2)
    
    
    a(i,6) = cos(torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
    a(i,7) = cos(2.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
    a(i,8) = cos(3.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
    if (n == 9) a(i,9) = cos(4.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
    
    write(6,'(10f12.6)') a(i,1:n) 

end do 

write(6,*)
write(6,*) "ZPE"
write(6,'(10f12.6)') zpe
write(6,*)

par = zpe 

call elgau(1, n, n, a, par, det, w, iw) 

write(6,*) "PARAMETRI"
write(6,'(4f20.12)') par(1:3) 
write(6,'(4f20.12)') par(4:5) 
write(6,'(4f20.12)') par(6:n)



!CALCOLO ZPE ALLE GEOMETRIE FORNITE

read(5,*) string
open(1,file = string, form = 'formatted', status = 'old')
read(5,*) string
open(2,file = string, form = 'formatted', status = 'unknown')
read(1,'(a80)') string
write(2,'(a80)') string

do
    read(1,'(a80)', end = 99) string
    
    if (index(string, '.') == 0) then
        write(2,'(a80)') string
    else
        read(string,*) r1, r2, ang1, ang2, cnnc

        ang1or = ang1

        if (ang1 > 180.0_dpr) then
            cnnc = 180.0_dpr
            ang1 = 360.0_dpr - ang1
        end if


        zp = 0.0_dpr

        call quartic(r1,v1)
        call quartic(r2,v2) 

        a(1,1) = v1 + v2 
        a(1,2) = v1*v2 
        a(1,3) = v1**2 + v2**2

        a(1,4) = ((1.0_dpr-v1)*cos(torad*ang1))+((1.0_dpr-v2)*cos(torad*ang2)) 
        a(1,5) = ((1.0_dpr-v1)*cos(2.0_dpr*torad*ang1))+((1.0_dpr-v2)*cos(2.0_dpr*torad*ang2))

        call cubic(ang1,carg1)
        call cubic(ang2,carg2)
    
        a(1,6) = cos(torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
        a(1,7) = cos(2.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
        a(1,8) = cos(3.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
        a(1,9) = cos(4.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)

        do j = 1,n
            zp = zp + a(1,j)*par(j)
        end do

        write(2,'(6f12.6)') r1, r2, ang1or, ang2, cnnc, zp
    end if

end do 

99 stop


!SUBROUTINES PER QUARTICA E CUBICA
contains

subroutine quartic(r,qua)
 
    implicit none
    real(kind = dpr), intent(in) :: r
    real(kind = dpr), intent(out) :: qua

    k1 = (r-r0)/(ra-r0) 
        if (k1 <= 0) then 
            qua = 0.0_dpr 
        else if (k1 >= 1) then 
            qua = 1.0_dpr 
        else if (k1>0 .and. k1<1) then 
            qua = 6.0_dpr*(k1**2) - 8.0_dpr*(k1**3) + 3.0_dpr*(k1**4) 
        end if

end subroutine quartic


subroutine cubic(ang,cub)

    implicit none
        
    real(kind = dpr), intent(in) :: ang
        real(kind = dpr), intent(out) :: cub

        arg1 = (ang - theta0)/(180.0_dpr - theta0)
        if(arg1 <= 0) then
            cub = 1.0_dpr
        else if (arg1 >=1 ) then
            cub = 0.0_dpr
        else if (arg1 > 0 .and. arg1 < 1) then
            cub = 1.0_dpr - 3.0_dpr*(arg1**2) + 2.0_dpr*(arg1**3)
        end if

end subroutine cubic

end program sistema