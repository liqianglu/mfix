!
!******************************************************************
!
!  DES input generator
!  Sample generator to generate particle positions and outputs in DES readable format
!  Please use this as a basis and perform any additional customizations as
!  needed and this module does not come with any guarantees
!
!  The particles generated by this code are greater than or equal to specified np and
!  MFIX will pick up the correct particles according to value of np specified in the code
!
!  TO DO:
!  1) Fix so that only np particle information is created
!  2) In MFIX read the entire input file and if the number of particles does not correspond to
!     the input deck, than flag an error.
!
!  Additions:
!  (Sept. 4, 2007): Added multiparticle capability
!
!  Author: Jay Boyalakuntla  (May-12-06)
!  Last Modified: S. Pannala (Sep-4-07)
!
!******************************************************************
!
!
	Program DES_Particle_Genrator

        implicit none

	integer, parameter :: maxp = 20000 ! Static allocation

	integer dim
	integer random		! 1=random distribution; 0=alligned
        integer accept, touch, ia
 	integer i, j, k, iseed
        integer n(100), np(100), np1, nptotal, nt
	integer nx, ny, nz
        integer ntypes

	real*8 xp, yp, zp, x(maxp), y(maxp), z(maxp)
	real*8 u, v, w
        real*8 xl, yl, zl
        real*8 radius(100), dia(100), dist, rad1(100)
        real*8 density(100)
        real*8 rad, maxdia, maxradius

        real*8 rnd1
        integer nrnd
        integer ntype_index(maxp)

        open(unit=10, file="Pgen.in", status='old')
        open(unit=20, file="particle_input.dat", status='replace')

	read (10,*) dim
	read (10,*) random
	read (10,*) ntypes
        read (10,*) (np(i), i=1,ntypes)
        read (10,*) (radius(i), i=1,ntypes)
        read (10,*) (density(i), i=1,ntypes)
        read (10,*) xl
        read (10,*) yl
        read (10,*) zl

	dia = 2.1*radius	! so that particles don't touch each other to begin with
        rad1 = 1.05*rad

        maxdia = maxval(dia(1:ntypes))
        maxradius = maxdia/2.0d0

        nptotal = sum(np(1:ntypes))

	nx = int(xl/maxdia)
	nz = int(zl/maxdia)
        if(dim.eq.2) nz = 1

        k = 1
        j = 1
        i = 1

        n(:) = 0

!       Specifying the initial distribution of the particles

	if(random.eq.0) then	! Ordered particle lattice
           do nt = 1, nptotal
 101	      continue
	      call random_number(rnd1)
	      nrnd = int(rnd1*ntypes) + 1
	      if(n(nrnd).lt.np(nrnd)) then
		 n(nrnd) = n(nrnd)+1
		 ntype_index(nt) = nrnd
	      else
		 go to 101
	      end if
	      j = int((nt-1)/(nx*nz)) + 1
	      k = int(((nt-(j-1)*(nx*nz)) - 1)/nx) + 1
	      i = nt-(k-1)*nx-(j-1)*(nx*nz)

	      z(nt) = 0.5*maxdia + (k-1)*maxdia
	      y(nt) = 0.5*maxdia + (j-1)*maxdia
	      x(nt) = 0.5*maxdia + (i-1)*maxdia
              if(xp.gt.xl.or.yp.gt.yl.or.zp.gt.zl) then
                 write(*,*) 'Particle location outside the domain', xp, yp, zp
                 stop
              end if
           end do
        else if(random.eq.1) then
           write(*,*) 'Random option not debugged/correct for multiparticle'
           iseed = 98765432
           call random_seed(iseed)
	   do k = 1, nz
              do i = 1, np1
 10		 continue
		 call random_particle(maxradius,xp,yp,zp,xl,yl,zl,dim)
		 x(i) = xp
		 y(i) = yp
		 z(i) = zp
		 dist = 0d0
		 do j = 1, i
		    if(j.ne.i) then
		       dist = sqrt((x(i)-x(j))**2 + (y(i)-y(j))**2 +&
		       (z(i)-z(j))**2)
		       if(dist.le.maxdia) then
			  go to 10
		       end if
		    end if
		 end do
              end do
           end do
        end if

!       setting the velocities to zero
	u = 0d0
	v = 0d0
	w = 0d0

        if(dim.eq.2) then
           do i = 1, nptotal
	      write(20,11) x(i), y(i), radius(ntype_index(i)), density(ntype_index(i)), u, v
           end do
	else
           do i = 1, nptotal
	      write(20,12) x(i), y(i), z(i), radius(ntype_index(i)), density(ntype_index(i)), u, v, w
           end do
	end if

 11     FORMAT (6(d10.4,2x))
 12     FORMAT (8(d10.4,2x))

        stop
        end

!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       Random particle
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
        subroutine random_particle(rad, xp1, yp1, zp1, xl1, yl1, zl1, dim1)

        integer i, dim1, ic
        real*8 rad, xp1, yp1, zp1, xl1, yl1, zl1
        real*8 rad1
        real*8 pxy(3)

	ic = 100000
	do i = 1, ic
	   call random_number(pxy)
	   xp1 = dble(pxy(1))*xl1
	   yp1 = dble(pxy(2))*yl1
	   zp1 = dble(pxy(3))*zl1
	   rad1 = 1.05*rad
	   if(dim1.eq.2) zp1 = rad1
	   rad1 = 1.05*rad
	   if(dim1.eq.2) then
	      if((xp1.ge.rad1).and.(xp1.le.xl1-rad1).and.(yp1.ge.rad1)&
	      .and.(yp1.le.yl1-rad1)) exit
	   else
	      if((xp1.ge.rad1).and.(xp1.le.xl1-rad1).and.(yp1.ge.rad1)&
	      .and.(yp1.le.yl1-rad1).and.(zp1.ge.rad1).and.&
	      (zp1.le.zl1-rad1)) exit
	   end if
	end do
	if(i.gt.ic) then
	   print *,'not able to place particle'
	   stop
	end if

        return
        end subroutine random_particle

!
