CCCCCC##################################################################
C     sshctomach.f
C     A program to make Mach number distribution from SSH components.
C     
C     Ver. 1.0; Made by Takashi Yoshida in 2020, 3, 29th.
C     
C     Please prepare the input data of the range of a convective region
C     in the 3D polar coordinate.
C     One needs to write the input file name, 
C     when this program is running.
C     The requested format in the input file is as follows.
C     Line 1: the raddi of the inner and outer boundaries (rin1 rout) 
C             in the units of 10^8 cm as real(8) variables.
C     Line 2: the numbers of the radial, polar, and azimuthal meshes
C             (ie je ke) as integer variables.
C     Lines 3 - 3+ie: the radii of the meshes in units of 10^8 cm
C             (x1b(1:ie)) as real(8) variables.
C     Lines 4+ie - 4+ie+je: the polar coordinates of meshes
C             in units of radian (x2b(1:je)) as real(8) variables.
C     Lines 5+ie+je - 5+ie+je+ke: the azimuthal coordinates of meshes
C             in units of radian (x3b(1:ke)) as real(8) variables.
C
C     The output file name is "machdata_out.txt".
C
CCCCCC##################################################################

      module coordinate_data
      integer :: ie,je,ke,iin1,iou1
      real(8) :: rin1,rou1
      real(8), allocatable :: x1b(:),x2b(:),x3b(:),xmu(:)
      real(8), allocatable :: mach(:,:,:)
      end module coordinate_data


      subroutine input_coordinate_data
      use coordinate_data
      integer :: i,j,k
      character*80 :: fname

      write(*,'("Please write input file name.")')
      read(*,'(a80)') fname
      open(10,file = fname)
c      open(10,file='meshdata_25M.txt')
      read(10,*,end=91) rin1,rou1
      read(10,*) ie,je,ke,iin1,iou1
      allocate(x1b(ie),x2b(je),x3b(ke),xmu(je),mach(ie,je,ke))
      do i = 1, ie
         read(10,*) x1b(i)
      enddo
      do j = 1, je
         read(10,*) x2b(j)
      enddo
      do k = 1, ke
         read(10,*) x3b(k)
      enddo
      close(10)
      xmu(1:je) = cos(x2b(1:je))

      return

 91   write(*,'("Input file is not found.")')
      stop

      end subroutine input_coordinate_data


      subroutine output

      use coordinate_data
      implicit NONE
      integer :: i,j,k

      open(12,file='machdata_out.txt')
      write(12,'("#",a15,3a16)') "x1b","x2b","x3b","Mach"
      do i = 1, ie
         do j = 1, je
            do k = 1, ke
               write(12,'(1p4e16.8)') x1b(i),x2b(j),x3b(k),mach(i,j,k)
            enddo
         enddo
      enddo
      close(12)

CC   This part is used for recovery check.
C      open(31,file='machdata1_m1_l05e8_119new2.txt')
Cc      open(31,file='../m22B0a23d/machdata1_l05e8_000new2.txt')
Cc      open(31,file='../m27B1b3d/machdata1_l25e8_000new2.txt')
C      write(31,*) iin1,iou1
Cc      write(31,*) 1,ie
C      do i = 1, ie
C         do j = 1, je
C            write(31,'(1p16e23.15)') mach(i,j,1:ke)
C         enddo
C      enddo
C      close(31)

      return
      end subroutine output


      program sshctomach

      use coordinate_data
      implicit NONE
      real(8), parameter :: pi=acos(-1.d0)
      integer :: mode,mbound
      integer :: nend,lend,n,l,m
      integer :: i,j,k,id
!     Boundary condisions
      real(8) :: rin0,rou0
!     Spherical Bessel functions
      real(8), allocatable :: gln(:,:),xk(:,:),xn(:,:)
      real(8) :: xk0,xk1,xkp,dk
      real(8) :: x1,sj1,sy1,sjp1,syp1,sjm1,sym1,sjpm1,sypm1
      real(8) :: x2,sj2,sy2,sjp2,syp2,sjm2,sym2,sjpm2,sypm2
      real(8) :: x,sj,sy,sjp,syp
      real(8) :: tll,tllp,tll0,dtlldk,tllm11,tllm12
!     Spherical harmonics
      real(8) :: xl,xm,cylm,powlpm,powlmm
      real(8), allocatable :: plxx(:),ylmp(:,:,:),ylmm(:,:,:)
!     SSH decomposition components
      real(8), allocatable :: amrenlm(:,:,:),amimnlm(:,:,:)
!     Renormalization coefficient for Mach number
      real(8) :: renorm


C     Specify model
      write(*,'("Please set the model number")')
      write(*,'("1: model 25M, 2: model 22L, 3: model 27LA")')
      read(*,*) mode
      if(mode .eq. 1) then
         rin0 = 2.30d0
         rou0 = 1.299d1
         nend = 42
         lend = 48
      elseif(mode .eq. 2) then
         rin0 = 2.13d0
         rou0 = 1.299d1
         nend = 43
         lend = 47
      elseif(mode .eq. 3) then
         rin0 = 5.45d0
         rou0 = 6.325d1
         nend = 46
         lend = 43
      else
         write(*,'("Please input number 1, 2, or 3.")')
         stop
      endif

C     Specify the coordinate (r,theta,phi) for output data
      call input_coordinate_data

      
C     Allocation and input of a_nlm,xk_nl,xn_nl
      allocate(amrenlm(0:nend,0:lend,-lend:lend),
     &     amimnlm(0:nend,0:lend,-lend:lend),
     &     xk(0:nend,0:lend),xn(0:nend,0:lend),gln(0:nend,ie),
     &     plxx(je),ylmp(je,ke,2),ylmm(je,ke,2))
      amrenlm(0:nend,0:lend,-lend:lend) = 0d0
      amimnlm(0:nend,0:lend,-lend:lend) = 0d0
      if(mode .eq. 1) then
         open(11,file='tab1.txt')
      elseif(mode .eq. 2) then
         open(11,file='tab2.txt')
      else
         open(11,file='tab3.txt')
      endif
      do n = 1,nend+14
         read(11,*)
      enddo
      do l = 0, lend
         do m = 0, l
            read(11,'(8x,50e16.8)') amrenlm(0:nend,l,m)
            read(11,'(8x,50e16.8)') amimnlm(0:nend,l,m)
         enddo
      enddo
      close(11)
      do l = 0, lend
         do m = 1, l
            amrenlm(0:nend,l,-m) = (-1d0)**m * amrenlm(0:nend,l,m)
            amimnlm(0:nend,l,-m) = (-1d0)**(m+1) * amimnlm(0:nend,l,m)
         enddo
      enddo

      do l=0,lend
         xl = dble(l)
!     calculation of xk(n,l) using midpoint method
         dk = pi / (rou1 - rin1)
         do n = 0, nend
            if(n .ge. 2) then
               dk = 1.d-1 * (xk(n-1,l) - xk(n-2,l))
               xk0 = xk(n-1,l) + dk
            elseif(n. eq. 1) then
               dk = 1.d-1*pi / (rou1 - rin1)
               xk0 = xk(n-1,l) + dk
            else
               dk = 2.d-1*pi / (rou1 - rin1)
               xk0 = 2.d-1*pi / (rou1 - rin1)
            endif
            x1 = xk0*rin1
            x2 = xk0*rou1
            call sphbes(l,x1,sj1,sy1,sjp1,syp1)
            call sphbes(l,x2,sj2,sy2,sjp2,syp2)
            tll0 = sy1*sj2 - sy2*sj1
            do j = 1, 200
               xkp = xk0
               tllp = tll0
               xk0  = xk0 + dk
               x1 = xk0*rin1
               x2 = xk0*rou1
               call sphbes(l,x1,sj1,sy1,sjp1,syp1)
               call sphbes(l,x2,sj2,sy2,sjp2,syp2)
               tll0 = sy1*sj2 - sy2*sj1
c               write(*,'(i4,1p4e12.4)') j,xk0,tll0,tllp
               if(tll0*tllp .le. 0.d0) exit
               if(j .eq. 200) stop
            enddo
C     xk1 should be between xkp and xk0
            do j = 1, 100
               xk1 = 5.d-1*(xkp + xk0)
               x1 = xk1*rin1
               x2 = xk1*rou1
               call sphbes(l,x1,sj1,sy1,sjp1,syp1)
               call sphbes(l,x2,sj2,sy2,sjp2,syp2)
               tll = sy1*sj2 - sy2*sj1
               if(tll*tllp .ge. 0.d0) then
                  xkp = xk1
                  tllp = tll
               elseif(tll*tll0 .ge. 0.d0) then
                  xk0 = xk1
                  tll0 = tll
               else
                  write(*,'("midpoint method is wrong!")')
               endif
               if(abs(dmax1(xk0-xk1, xk1-xkp)/xk1) .lt. 1d-10) exit
            enddo
            if(j .ge. 100) then
               write(*,'("Failed iterations: n=",i2," l=",i2)') n,l
               stop
            endif
            xk(n,l) = xk1
!     Calculation of normalization coefficient of gln(r), xn(n,l)
            x1 = xk(n,l)*rin1
            x2 = xk(n,l)*rou1
            if(l .ne. 0) then
               call sphbes(l-1,x1,sjm1,sym1,sjpm1,sypm1)
               call sphbes(l-1,x2,sjm2,sym2,sjpm2,sypm2)
               tllm11 = sy1*sjm1 - sj1*sym1
               tllm12 = sy1*sjm2 - sj1*sym2
            else
               call sphbes(0,x1,sj1,sy1,sjp1,syp1)
               call sphbes(0,x2,sj2,sy2,sjp2,syp2)
               tllm11 = -(sy1*sy1 + sj1*sj1)
               tllm12 = -(sy1*sy2 + sj1*sj2)
            endif
            xn(n,l) = 1.d0/(5.d-1*(rou1*rou1*rou1*tllm12*tllm12
     &           - rin1*rin1*rin1*tllm11*tllm11))**5.d-1
         enddo
c        write(*,'(i4,1p4e12.4)') l,xk(1,l),xn(1,l),xk(nend,l),xn(nend,l)
      enddo


C     Calculation of Mach number
      do l = 0, lend
         xl = dble(l)
!     calculation of gln(r)
         do n = 0, nend
            x1 = xk(n,l)*rin1
            call sphbes(l,x1,sj1,sy1,sjp1,syp1)
            do i = 1, ie
               x = xk(n,l)*x1b(i)
               call sphbes(l,x,sj,sy,sjp,syp)
               gln(n,i) = xn(n,l)*(sy1*sj - sj1*sy)
            enddo
         enddo

!     m = 0
!     Calculation of spherical harmonics
         cylm = sqrt(2.5d-1 * (2.d0*xl+1.d0) / pi)
         do j = 1, je
            call plgndr_func(l,0,xmu(j),plxx(j))
            plxx(j) = cylm*plxx(j)
         enddo
!     Mach number for m = 0 component
         do n = 0, nend
            do i = 1, ie
               do j = 1, je
                  mach(i,j,1:ke) = mach(i,j,1:ke)
     &                 + amrenlm(n,l,0)*gln(n,i)*plxx(j)
               enddo
            enddo
         enddo

!     m != 0 Ylm
         do m = 1, l
            xm = dble(m)
            call calc_pow(l-m,powlmm)
            call calc_pow(l+m,powlpm)

            cylm = sqrt(2.5d-1 * (2.d0*xl+1.d0) / pi * powlmm / powlpm)
            do j = 1, je
               call plgndr_func(l,m,xmu(j),plxx(j))
               do k = 1, ke
                  ylmp(j,k,1) = cylm*plxx(j)*cos(xm*x3b(k)) !real part
                  ylmp(j,k,2) = cylm*plxx(j)*sin(xm*x3b(k)) !imaginary part
                  ylmm(j,k,1) = (-1.d0)**m     * ylmp(j,k,1) !real part
                  ylmm(j,k,2) = (-1.d0)**(m+1) * ylmp(j,k,2) !imaginary part
               enddo
            enddo

            do n = 0, nend
               do i = 1, ie
                  do j = 1, je
                     do k = 1, ke
                        mach(i,j,k) = mach(i,j,k) + gln(n,i)
     &                       *(amrenlm(n,l,m)*ylmp(j,k,1)
     &                       + amimnlm(n,l,m)*ylmp(j,k,2)
     &                       + amrenlm(n,l,-m)*ylmm(j,k,1)
     &                       + amimnlm(n,l,-m)*ylmm(j,k,2))
                     enddo
                  enddo
               enddo
            enddo

         enddo
      enddo

      renorm = sqrt((rou1**3-rin1**3)/(rou0**3-rin0**3))
      mach(1:ie,1:je,1:ke) = renorm * mach(1:ie,1:je,1:ke)

      call output

      deallocate(amrenlm,amimnlm,xk,xn,gln,plxx,ylmp,ylmm,mach)

      end



      subroutine plgndr_func(l,m,x,plgndr)

C     This subroutine is based on FUNCTION plgndr(l,m,n)
C     in Section 6.8, Numerical Recipes in Fortran 77, 2nd Edition, 
C     Press, W. H. et al. (1992).
      implicit NONE
      integer :: l,m
      real(8) :: x, plgndr
c     Computes the associated Legendre polynominal Plm(x).
c     Here m and l are integers satisfying 0 <= m <= l,
c     while x lies in the range -1 <= x <= 1.
      integer :: i,ll
      real(8) :: fact,pll,pmm,pmmp1,somx2

      if(m .lt. 0 .or. m .gt. l .or. abs(x) .gt. 1d0) then
         write(*,*) "bad arguments in plgndr"
         stop
      endif

      pmm = 1.d0                !Compute Pmm
      if(m .gt. 0) then
         somx2 = sqrt((1.d0-x)*(1.d0+x))
         fact = 1.d0
         do i = 1, m
            pmm = -pmm * fact * somx2
            fact = fact + 2.d0
         enddo
      endif

      if(l .eq. m) then
         plgndr = pmm
      else
         pmmp1 = x * (2*m+1) * pmm !Compute P^m_m+1
         if(l .eq. m+1) then
            plgndr = pmmp1
         else                   !Compute P^m_l, l > m+1
            do ll = m+2, l
               pll = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/dble(ll-m)
               pmm = pmmp1
               pmmp1 = pll
            enddo
            plgndr = pll
         endif
      endif

      return
      end subroutine plgndr_func


      subroutine calc_pow(m,powm)

      implicit NONE
      integer :: m
      integer :: i
      real(8) :: powm

      if(m .ge. 0) then
         powm=1.d0
         if(m .gt. 1) then
            do i = 2, m
               powm = powm * dble(i)
            enddo
         endif
      else
         write(*,'("m should be zero or positive integer.")')
         write(*,'("Error: m=",i2)') m
         stop
      endif

      return
      end subroutine calc_pow


      subroutine sphbes(n,x,sj,sy,sjp,syp)

C     This subroutine is from Section 6.7 in Numerical Recipes in
C     Fortran 77, 2nd Edition, Press, W. H. et al. (1992).
      implicit NONE
      integer :: n
      real(8) :: sj,sjp,sy,syp,x
C     USES bessjy
c     Returns spherical Bessel functions jn(x), yn(x), and their
c        derivatives j'n(x), y'n(x) for integer n.
      real(8), parameter :: RTPIO2=1.2533141d0
      real(8) :: factor,order,rj,rjp,ry,ryp

      if(n.lt.0 .or. x.le.0d0) then
         write(*,'("Bad arguments in sphbes")')
         stop
      endif

      order = n + 5d-1
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor = RTPIO2/sqrt(x)
      sj = factor*rj
      sy = factor*ry
      sjp = factor*rjp - sj/(2d0*x)
      syp = factor*ryp - sy/(2d0*x)

      return
      end subroutine sphbes



      subroutine bessjy(x,xnu,rj,ry,rjp,ryp)

C     This subroutine is from Section 6.7 in Numerical Recipes in
C     Fortran 77, 2nd Edition, Press, W. H. et al. (1992).
      implicit NONE
      integer, parameter :: MAXIT=10000
      real(8), parameter :: EPS=1d-16,FPMIN=1d-30,XMIN=2d0
     &     ,PI=3.141592653589793d0
C     USES beschb
c     Returns the Bessel functions rj = J_nu, ry = Y_nu and their
c     derivatives rjp = J'_nu, ryp = Y'_nu, for positive x and for
c     xnu = nu >= 0. The relative accuracy is within one or two
c     significant digits of EPS, except for near a zero of one of
c     the functions, where EPS controls its abusolute accuracy.
c     FPMIN is a number close to the machine's smallest floating-point
c     number. All internal arithmetic is in double precision.
c     To convert the entire routine to double presition, change the
c     REAL declararation above and decrease EPS to 10^{-16}.
c     Aslo convert the subroutine beschb.
      real(8) :: rj,rjp,ry,ryp,x,xnu

      integer :: i,isign,l,nl
      real(8) :: a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,
     &     dr,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,
     &     p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,
     &     rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2

      if(x .le. 0d0 .or. xnu .lt. 0d0) then
         write(*,'("Bad arguments in bessjy")')
         stop
      endif
      if(x .lt. XMIN) then
         nl = int(xnu + 0.5d0)
      else
         nl = max(0, int(xnu-x+1.5d0))
      endif

      xmu = xnu - nl
      xmu2 = xmu*xmu
      xi = 1d0/x
      xi2 = 2d0*xi
      w = xi2/PI
      isign = 1
      h = xnu*xi
      if(h .lt. FPMIN) h = FPMIN
      b = xi2*xnu
      d = 0d0
      c = h

      do i = 1, MAXIT
         b = b + xi2
         d = b - d
         if(abs(d) .lt. FPMIN) d = FPMIN
         c = b - 1d0/c
         if(abs(c) .lt. FPMIN) c = FPMIN
         d = 1d0/d
         del = c*d
         h = del*h
         if(d .lt. 0d0) isign = -isign
         if(abs(del-1d0) .lt. EPS) exit
      enddo
      if(i .ge. MAXIT) then
         write(*,'("x too large in bessjy; try asymptotic expansion")')
         stop
      endif

      rjl = isign*FPMIN
      rjpl = h*rjl
      rjl1 = rjl
      rjp1 = rjpl
      fact = xnu*xi
      do l = nl, 1, -1
         rjtemp = fact*rjl + rjpl
         fact = fact - xi
         rjpl = fact*rjtemp - rjl
         rjl = rjtemp
      enddo
      if(rjl .eq. 0d0) rjl = EPS
      f = rjpl/rjl

      if(x .lt. XMIN) then
         x2 = 5d-1*x
         pimu = PI*xmu
         if(abs(pimu) .lt. EPS) then
            fact = 1d0
         else
            fact = pimu/sin(pimu)
         endif
         d = -log(x2)
         e = xmu*d
         if(abs(e) .lt. EPS) then
            fact2 = 1d0
         else
            fact2 = sinh(e)/e
         endif
         call beschb(xmu,gam1,gam2,gampl,gammi)
         ff = 2d0/PI * fact * (gam1*cosh(e) + gam2*fact2*d)
         e = exp(e)
         p = e/(gampl*PI)
         q = 1d0/(e*PI*gammi)
         pimu2 = 5d-1*pimu
         if(abs(pimu2) .lt. EPS) then
            fact3 = 1d0
         else
            fact3 = sin(pimu2)/pimu2
         endif
         r = PI*pimu2*fact3*fact3
         c = 1d0
         d = -x2*x2
         sum = ff + r*q
         sum1 = p

         do i = 1, MAXIT
            ff = (i*ff+p+q)/(i*i-xmu2)
            c = c*d/i
            p = p/(i-xmu)
            q = q/(i+xmu)
            del = c*(ff+r*q)
            sum = sum+del
            del1 = c*p - i*del
            sum1 = sum1 + del1
            if(abs(del) .lt. (1d0+abs(sum))*EPS) exit
         enddo
         if(i .ge. MAXIT) then
            write(*,'("bessy series faild to converge")')
            stop
         endif

         rymu = -sum
         ry1 = -sum1*xi2
         rymup = xmu*xi*rymu - ry1
         rjmu = w/(rymup - f*rymu)
      else
         a = 0.25d0 - xmu2
         p = -0.5d0*xi
         q = 1d0
         br = 2d0*x
         bi = 2d0
         fact = a*xi/(p*p+q*q)
         cr = br + q*fact
         ci = bi + p*fact
         den = br*br + bi*bi
         dr = br/den
         di = -bi/den
         dlr = cr*dr - ci*di
         dli = cr*di + ci*dr
         temp = p*dlr - q*dli
         q = p*dli + q*dlr
         p = temp
         do i = 2, MAXIT
            a = a+2*(i-1)
            bi = bi + 2d0
            dr = a*dr + br
            di = a*di + bi
            if(abs(dr) + abs(di) .lt. FPMIN) dr = FPMIN
            fact = a/(cr*cr + ci*ci)
            cr = br + cr*fact
            ci = bi - ci*fact
            if(abs(cr) + abs(ci) .lt. FPMIN) cr = FPMIN
            den = dr*dr + di*di
            dr = dr/den
            di = -di/den
            dlr = cr*dr - ci*di
            dli = cr*di + ci*dr
            temp = p*dlr - q*dli
            q = p*dli +q*dlr
            p = temp
            if(abs(dlr-1d0)+abs(dli) .lt. EPS) exit
         enddo
         if(i .ge. MAXIT) then
            write(*,'("cf2 failed in bessjy")')
            stop
         endif
         
         gam = (p-f)/q
         rjmu = sqrt(w/((p-f)*gam + q))
         rjmu = sign(rjmu,rjl)
         rymu = rjmu*gam
         rymup = rymu*(p+q/gam)
         ry1 = xmu*xi*rymu - rymup
      endif

      fact = rjmu/rjl
      rj = rjl1*fact
      rjp = rjp1*fact
      do i = 1, nl
         rytemp = (xmu + i)*xi2*ry1 - rymu
         rymu = ry1
         ry1 = rytemp
      enddo
      ry = rymu
      ryp = xnu*xi*rymu -ry1
      
      return
      end subroutine bessjy


      subroutine beschb(x,gam1,gam2,gampl,gammi)
      
C     This subroutine is from Section 6.7 in Numerical Recipes in
C     Fortran 77, 2nd Edition, Press, W. H. et al. (1992).
      implicit NONE
      integer, parameter :: NUSE1=7,NUSE2=8
!     USES chebev
!     Evaluates gamma1 and gamma2 by Chebyshev expansion for |x|=<1/2.
!     Also returns 1/gamma(1+x) and 1/gamma(1-x).
!     If converting to double precision, set NUSE1 = 7, NUSE2 = 8.
      real(8) :: gam1,gam2,gammi,gampl,x
      real(8) :: xx,c1(7),c2(8),chebev
      save :: c1,c2
      
      data c1/-1.142022680371168d0,6.5165112670737d-3,
     &     3.087090173086d-4,-3.4706269649d-6, 6.9437664d-9,
     &     3.67795d-11,-1.356d-13/
      data c2/1.843740587300905d0,-7.68528408447867d-2,
     &     1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8,
     &     2.423096d-10,-1.702d-13,-1.49d-15/

      xx = 8d0*x*x - 1d0
      gam1 = chebev(-1d0,1d0,c1,NUSE1,xx)
      gam2 = chebev(-1d0,1d0,c2,NUSE2,xx)
      gampl = gam2 - x*gam1
      gammi = gam2 + x*gam1

      return
      end subroutine beschb


      function chebev(a,b,c,m,x)

C     This function is from Section 5.8 in Numerical Recipes in
C     Fortran 77, 2nd Edition, Press, W. H. et al. (1992).
      implicit NONE
      integer :: m
      real(8) :: chebev,a,b,x,c(m)
c     Chebyshev evaluation: All arguments are input. c(1:m) is an array
c     of Chebyshev coefficients, the first m elements of c output from
c     chebft (which must have been called with the same a and b).
c     The Chebyeshev polynominal Sigma_{k=1}^{m} c_k T_{k-1}(y)-c_1/2
c     is evaluated at a point y = [x-(b+a)/2]/[(b-a)/2], and the result
c     is returned as the function value.
      integer :: j
      real(8) :: d,dd,sv,y,y2
      
      if((x-a)*(x-b) .gt. 0d0) then
         write(*,'("x not in range in chebev")')
         stop
      endif
      d = 0d0
      dd = 0d0
      y = (2d0*x - a - b)/(b - a)
      y2 = 2d0*y
      do j = m, 2, -1
         sv = d
         d = y2*d - dd + c(j)
         dd = sv
      enddo
      chebev = y*d - dd + 5d-1*c(1)
      
      return
      end function chebev
