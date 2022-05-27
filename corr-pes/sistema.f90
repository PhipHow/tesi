program sistema 

implicit none 

integer :: i,j,n 
integer, parameter :: dpr = kind(1.d0) 
real (kind = dpr) :: r0, ra, theta0, r1, r2, ang1, ang2, cnnc, det, &
                     k1, v1, v2, pi, torad, arg1, carg1, carg2, zp, ang1or, &
                     e0, e1, e0z, e1z, f1, f2, zp1
real (kind = dpr), allocatable, dimension(:) :: zpe, par, w, b1, b2 !b1-->vettore dei B old , b2--> vettore dei B new
real (kind = dpr), allocatable, dimension(:,:) :: a 
integer, allocatable, dimension(:) :: iw 

character*80 :: string 
 

DATA PI,TORAD/ & 
3.141592653589793D0,0.01745329251994329D0/

!PARAMETRI NON LINEARI UTILIZZATI
r0 = 1.6 
ra = 3.5
theta0 = 126.8176185 

read(5,*) n 
allocate(zpe(1:n), par(1:n), a(1:n,1:n), w(1:n), iw(1:n), b1(0:9), b2(0:9)) 

!INIZIO CALCOLO MATRICE A
write(6,*) "MATRICE A" 

do i = 1,n 
    
    !LETTURA DELLA GEOMETRIA + ZPE
    read(5,*) r1, r2, ang1, ang2, cnnc, zpe(i) 
    
    !CALCOLO DEGLI ELEMENTI DELLA MATRICE A
    call coeff(n,a(i,1:n))
    
    !SCRITTURA DEGLI ELEMENTI DELLA MATRICE A
    write(6,'(10f12.6)') a(i,1:n) 

end do 
!FINE CALCOLO MATRICE A

write(6,*)
!SCRITTURA DEL VETTORE ZPE UTILIZZATO
write(6,*) "ZPE"
write(6,'(10f12.6)') zpe
write(6,*)

par = zpe 

!RISOLUZIONE DEL SISTEMA TRAMITE ELIMINAZIONE DI GAUSS
call elgau(1, n, n, a, par, det, w, iw) 

!SCRITTURA DEI PARAMETRI
write(6,*) "PARAMETRI"
write(6,'(4f20.12)') par(1:3) !A1, A2, A3
write(6,'(4f20.12)') par(4:5) !Ainv1, Ainv2
write(6,'(4f20.12)') par(6:n) !Arot1, Arot2, Arot3, Arot4

read(5,*) b1 !VENGONO LETTI I PARAMETRI B OLD
read(5,*) b2 !VENGONO LETTI I PARAMETRI B NEW


!-----INIZIO CALCOLO ZPE ALLE GEOMETRIE FORNITE-------!

!VIENE CERCATO NEL FILE DI DATI IL NOME DEL FILE DA CUI PRENDERE LE GEOMETRIE 
read(5,*) string
!IL FILE "STRING" VIENE APERTO
open(1,file = string, form = 'formatted', status = 'old')
read(5,*) string
!VIENE APERTO IL FILE SU CUI VENGONO SCRITTI I RISULTATI
open(2,file = string, form = 'formatted', status = 'unknown')
!VENGONO LETTI 80 CARATTERI DAL FILE DELLE GEOMETRIE (TUTTA LA RIGA)
read(1,'(a80)') string
!LA PRIMA RIGA LETTA VIENE SCRITTA SUL FILE DEI RISULTATI
write(2,'(a80)') string

do
    read(1,'(a80)', end = 99) string
    
    !SE NON VIENE TROVATO UN "." NEL FILE LA RIGA VIENE UNICAMENTE COPIATA
    if (index(string, '.') == 0) then
        write(2,'(a80)') string
    else
        !LETTURA PRIMA GEOMETRIA SUL FILE DI DATI
        read(string,*) r1, r2, ang1, ang2, cnnc, e0, e1

        ang1or = ang1

        if (ang1 > 180.0_dpr) then
            cnnc = 180.0_dpr
            ang1 = 360.0_dpr - ang1
        end if


        zp = 0.0_dpr

        !CALCOLATI I COEFFICIENTI PER LE GEOMETRIE PASSATE DA FILE
        call coeff(n,w)
        
        !call quartic(r1,v1)
        !call quartic(r2,v2) 

        !a(1,1) = v1 + v2 
        !a(1,2) = v1*v2 
        !a(1,3) = v1**2 + v2**2

        !a(1,4) = ((1.0_dpr-v1)*cos(torad*ang1))+((1.0_dpr-v2)*cos(torad*ang2)) 
        !a(1,5) = ((1.0_dpr-v1)*cos(2.0_dpr*torad*ang1))+((1.0_dpr-v2)*cos(2.0_dpr*torad*ang2))

        !call cubic(ang1,carg1)
        !call cubic(ang2,carg2)
    
        !a(1,6) = cos(torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
        !a(1,7) = cos(2.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
        !if (n > 7) a(1,8) = cos(3.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
        !if (n == 9) a(1,9) = cos(4.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)

        !VIENE CALCOLATA LA ZPE (ZP) PER LE VARIE GEOMETRIE
        do j = 1,n
            zp = zp + w(j)*par(j)
        end do

        f1 = b1(0) !F1 ASSUME IL VALORE DI BO (LA COSTANTE OLD)

        !CALCOLO DEL FATTORE F1 (THE OLD ONE)
        call coeff(9,w)

        do j = 1,9 
            f1 = f1 + w(j)*b1(j)
        end do

        !CALCOLO DEL FATTORE F2 (THE NEW ONE)
        f2 = b2(0) !F2 ASSUME IL VALORE DI BO (LA COSTANTE NEW)

        do j = 1,9 
            f2 = f2 + w(j)*b2(j)
        end do

        e0z = e0 + zp                   !S0 CORRETTA CON ZPE
        e1z = e0z + (e1 - e0)*(f2/f1)    !S1 CORRETTA CON ZPE
        zp1 = e1z - e1             !CORREZIONE ZPE S1
        
        !NEL FILE DEI RISULTATI VENGONO SCRITTI QUESTI VALORI
        write(2,'(9f12.6)') r1, r2, ang1or, ang2, cnnc, e0z, e1z, zp, zp1
    end if

end do 


99 stop


!SUBROUTINES PER QUARTICA E CUBICA
contains

subroutine coeff(nn,ww)
    
    implicit none
    
    integer, intent(in) :: nn
    real(kind = dpr), dimension(:), intent(out) :: ww

    !calcolo dei coefficienti di A1,A2,A3 
    call quartic(r1,v1)
    call quartic(r2,v2) 

    ww(1) = v1 + v2 
    ww(2) = v1*v2 
    ww(3) = v1**2 + v2**2 


    !calcolo dei coefficienti di Ainv1,Ainv2,Ainv3 
    ww(4) = ((1.0_dpr-v1)*cos(torad*ang1))+((1.0_dpr-v2)*cos(torad*ang2)) 
    ww(5) = ((1.0_dpr-v1)*cos(2.0_dpr*torad*ang1))+((1.0_dpr-v2)*cos(2.0_dpr*torad*ang2)) 

    !calcolo dei coefficienti di Arot1,Arot2,Arot3,Arot4
    call cubic(ang1,carg1)
    call cubic(ang2,carg2)


    ww(6) = cos(torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
    ww(7) = cos(2.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
    if (nn > 7) ww(8) = cos(3.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)
    if (nn == 9) ww(9) = cos(4.0_dpr*torad*cnnc)*carg1*carg2*(1.0_dpr-v1)*(1.0_dpr-v2)

end subroutine coeff


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