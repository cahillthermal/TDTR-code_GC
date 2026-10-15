!     Last change:  DGC   6 Dec 2006    8:43 pm
       module slabpar
        implicit none
        INTEGER :: nlayers, callrad=0
        REAL, DIMENSION(10) :: h, hcap, cond, diff
        COMPLEX, parameter :: i = (0.0, 1.0)
        real :: w0, w02, diffmax
        REAL, parameter :: pisq = 3.14159**2, ulmax = 10.0
        COMPLEX, DIMENSION(10) :: k2
       end module slabpar

       program freqd

!            calculates  thermal response of mulilayer slab in a mode-locked
!            pump probe experiment
!            modified to include full calculation of radial heat flow

       USE slabpar
       implicit none

       interface
        subroutine gausscalc(omega,tt)
         implicit none
         COMPLEX , INTENT(out) :: tt
         REAL, INTENT(IN) :: omega
        end subroutine gausscalc
       END interface

       integer :: n, nend, j, jj

       complex :: temp, tempin, tempout, tm1, tp1, phasef
       real :: freq, time, pumpf, omega
       REAL, parameter :: pi2=2*3.14159, btime = 1e-12
       REAL, PARAMETER :: efreq = 1.0e12, etime = 4e-9, probef = 80.0e6
       INTEGER, parameter :: nmax=1.01*efreq/probef
       COMPLEX, Dimension (nmax) :: tm, tp

        open(1, status='old', file='freqdv4.par')

!
!           first get the spot 1/e^2 radius and modulation frequency of the pump
!
        WRITE(*,*) ' laser repetition rate:', probef/1e6, ' MHz'
        read(1,*) w0, pumpf
        w02 = w0**2
        write(*,*) ' 1/e^2 radius:', w0,' modulation frequency:', pumpf/1e6, 'MHz'
        WRITE(*,*) ' '
!
!                 then read in the properties of each layer
!
        write(*,*) 'conductivity, heat capacity, thickness:'
        diffmax = 0.0
        do j = 1,10
         read(1,*,end=999) cond(j),hcap(j),h(j)
         diff(j) = cond(j)/hcap(j)
         diffmax = MAX(diff(j),diffmax)
         write(*,110) cond(j), hcap(j), h(j)
 110     FORMAT(3(ES12.3))
        end do
!
!                 define values of the nlayer+1 layer for zero heat flow at bottom
!
 999    cond(j) = 0.0
        diff(j) = 0.0
        nlayers = j-1
        write(*,*) 'total number of layers:', nlayers
        write(*,*) ' '


!               open a file for the data

       open(2, status='unknown', file='freqdv4.dat')

!                 loop over all frequencies starting with 0

       write(*,*) '           Hz, creating frequency response table'
       jj = 1

       freq = -probef
       do WHILE (freq .LE. efreq)
       WRITE(*,111) freq
 111   FORMAT('+', es10.2 )
        freq = freq+probef

        omega = ABS((freq-pumpf)*pi2)
        call gausscalc(omega,tm1)

        omega = (freq+pumpf)*pi2
        call gausscalc(omega,tp1)

!             add a gaussian roll-off to reduce some of the ringing

        tm(jj) = tm1*EXP(-2.0*(freq/efreq)**2)
        tp(jj) = tp1*EXP(-2.0*(freq/efreq)**2)
       jj = jj+1

       END do

       tm(1) = CONJG(tm(1))

       nend = jj-1

       WRITE(*,*) ' made ', callrad, ' calls to function radial'
       WRITE(*,*)

!         at each time, sum up the real part of the frequency response

        WRITE(*,*) 'calculating time domain signal'

!                calculate a few points at negative times

        time = -100e-12
        do WHILE (TIME .Lt. 0)
         time = time+10e-12
         tempin  = (tp(1)+tm(1))/2.0
         tempout = (tp(1)-tm(1))/2.0

        do n = 2, nend
          freq = (n-1)*probef
          phasef = EXP(i*time*freq*pi2)
          tempin = tempin+(tp(n)+tm(n))*phasef
          tempout = tempout+(tp(n)-tm(n))*phasef
        END do

        temp = REAL(tempin) + i*AIMAG(tempout)

!           rotate the phase of temp because we are delaying the pump not the probe

          temp = probef*temp*EXP(i*time*pumpf*pi2)

        write(2,100) time, temp, -REAL(temp)/AIMAG(temp)
 100    FORMAT(4(ES15.4))
        end do

        time = 0.0
        do WHILE (TIME .LE. etime)
         time = (time**.5+btime**.5/3)**2
         tempin  = (tp(1)+tm(1))/2.0
         tempout = (tp(1)-tm(1))/2.0

        do n = 2, nend
          freq = (n-1)*probef
          phasef = EXP(i*time*freq*pi2)
          tempin = tempin+(tp(n)+tm(n))*phasef
          tempout = tempout+(tp(n)-tm(n))*phasef
        END do

        temp = REAL(tempin) + i*AIMAG(tempout)

          temp = probef*temp*EXP(i*time*pumpf*pi2)
        write(2,100) time, temp, -REAL(temp)/AIMAG(temp)
        end do

       stop
       end program freqd

!---------------------------------------------------------------------------------
       subroutine gausscalc(omega,tt)

       USE slabpar
       USE nr, ONLY : qromo, midpnt, midinf, qromb
!
!              subroutine to calculate the complex temperature of a multilayer
!              slab at frequency omega, full calculation or radial heat flow
!
       implicit none
       interface
         function radialr(m)
         USE nrtype
         implicit none
         REAL(SP), DIMENSION(:), INTENT(IN) :: m
         REAL(SP), DIMENSION(size(m)) :: radialr
        END function radialr
        end interface

        interface
         function radiali(m)
         USE nrtype
         implicit none
         REAL(SP), DIMENSION(:), INTENT(IN) :: m
         REAL(SP), DIMENSION(size(m)) :: radiali
        END function radiali
        end interface

        interface
          function radial(m1)
          implicit none
          REAL, INTENT(IN) :: m1
          COMPLEX :: radial
          END function radial
        end interface


       REAL, INTENT(IN) :: omega
       COMPLEX, INTENT(out) :: tt
       INTEGER :: j
       REAL :: ttr, tti, cross, k2min

!                calculate iw/D for each layer at this frequency w

       do j=1, nlayers
        k2(j) = i*omega/diff(j)
       end do
!                check if full integration is needed

        k2min = omega/diffmax
        if (k2min*w02 .lt. 1e5 ) then
!                do the real and imaginary parts of the integration

        cross = 1.5/w0

        ttr = qromb( radialr, 0.0, cross )
        tti = qromb( radiali, 0.0, cross )

        tt = ttr+i*tti
        else
         tt = radial(0.0)*0.5/pisq/w02
        END if
        tt = tt*2.0*3.14159

       return
       END subroutine gausscalc
!
!--------------------------------------------------------------------------
!
       function radial(m)
        USE slabpar
        implicit none

        interface
         function u(j,m2)
         implicit none
         COMPLEX :: u
         integer, INTENT(IN) :: j
         REAL, INTENT(IN) :: m2
        END function u                                               
        end interface

        real, INTENT(IN) :: m
        COMPLEX :: radial
        REAL :: m2

        COMPLEX, DIMENSION(2) :: B, Btemp, Btemp2
        COMPLEX, DIMENSION(2,2) :: alpha, beta=0.0
        INTEGER :: j
        COMPLEX :: gj, gjp1, uj, ujp1

!             keep track of how many times this function is called

        callrad=callrad+1

!              initialize the B vector

       B(1) = CMPLX(0.0, 0.0)
       B(2) = CMPLX(1.0, 0.0)

       m2 = m**2

       j = nlayers

       do WHILE(j .ge. 1)

        uj = u(j, m2)

        IF( REAL(uj)*h(j) .LT. ulmax ) then

!             calculate the new B vector iteratively

           ujp1 = u(j+1, m2)

           gjp1 = cond(j+1)*ujp1
           gj = cond(j)*uj

           alpha(1,1) = gj+gjp1
           alpha(1,2) = gj-gjp1
           alpha(2,1) = alpha(1,2)
           alpha(2,2) = alpha(1,1)

           Btemp = matmul(alpha,B)

           beta(1,1) = EXP(-uj*h(j))
           beta(2,2) = 1.0/beta(1,1)

           Btemp2 = MATMUL(beta,Btemp)

           B = Btemp2/gj
        else

!              heat flow doesn't make it to the bottom of this layer so reset B vector

           B(1) = CMPLX(0.0, 0.0)
           B(2) = CMPLX(1.0, 0.0)
        END if

        j = j-1

        END do

        radial = (B(1)+B(2))/(B(2)-B(1))/u(1,m2)/cond(1)

       return
       END function radial

!-------------------------------------------------------------------------------

       function u(j,m2)
       USE slabpar
       implicit none
       INTEGER, INTENT(IN) :: j
       REAL, INTENT(IN) :: m2
       COMPLEX :: u

       u = (k2(j)+4.0*pisq*m2)**0.5
       return
       end function u

!----------------------------------------------------------------------------------

       function radiali(m)
         USE slabpar
         use nrtype
         implicit none

         interface
            function radial(m1)
            implicit none
            REAL, INTENT(IN) :: m1
            COMPLEX :: radial
            END function radial
         end interface

         REAL(SP), DIMENSION(:), INTENT(IN) :: m
         REAL(SP), DIMENSION(size(m)) :: radiali
         REAL :: msc
         INTEGER :: j

         do j = 1, SIZE(m)
           msc = m(j)
           radiali(j) = aimag(radial(msc))*EXP(-pisq*msc**2*w02)*msc
         end do

         return
        END function radiali
!
!---------------------------------------------------------------------------------
!
       function radialr(m)
         USE slabpar
         use nrtype
         implicit none

         interface
            function radial(m1)
            implicit none
            REAL, INTENT(IN) :: m1
            COMPLEX :: radial
            END function radial
         end interface

         REAL(SP), DIMENSION(:), INTENT(IN) :: m
         REAL(SP), DIMENSION(size(m)) :: radialr
         INTEGER :: j
         REAL :: msc

         do j = 1, SIZE(m)
           msc = m(j)
           radialr(j) = real(radial(msc))*EXP(-pisq*msc**2*w02)*msc
         end do

         return
        END function radialr
