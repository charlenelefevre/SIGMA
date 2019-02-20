!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

	subroutine RegridDataLNK(input,grid,e1,e2,n,loglog)
	use Tools
	IMPLICIT NONE
	real*8 grid(n)
	real*8 e1(n),e2(n),x0,y01,y02,x1,y11,y12,wp,gamma
	integer i,j,n, n_line, iostatus
	external input
	real*8,allocatable :: x(:),y1(:),y2(:)
	integer n0,i0
	logical loglog
	complex*16 m,m0,m1,cdlog10

	allocate(x(5000))
	allocate(y1(5000))
	allocate(y2(5000))

	call input(x,y1,y2,n0)


	i=1
	i0=1
	x0=x(i0)
	y01=y1(i0)
	y02=y2(i0)
	wp=(1d0-y01)/x0**2
	gamma=y02/x0**3
103	if(x0.ge.grid(i)) then
		e1(i)=1d0-wp*grid(i)**2
		e2(i)=gamma*grid(i)**3
		e1(i)=y01
		e2(i)=y02
		i=i+1
		goto 103
	endif
100	i0=i0+1
	if(i0.gt.n0) goto 102
	x1=x(i0)
	y11=y1(i0)
	y12=y2(i0)
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		e1(i)=y11+(grid(i)-x1)*(y01-y11)/(x0-x1)
		e2(i)=y12+(grid(i)-x1)*(y02-y12)/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y01=y11
	y02=y12
	goto 100
102	continue

	if(loglog) then
		m0=dcmplx(e1(i-1),e2(i-1))
		if(abs(m0).gt.2d0.and..false.) then
! don't use the conducting extrapolation since it is not very accurate
			do j=i,n
				m=m0*sqrt(grid(j)/grid(i-1))
				e1(j)=real(m)
				e2(j)=dimag(m)
			enddo
		else
! use loglog extrapolation
			m0=dcmplx(e1(i-2),e2(i-2))
			m1=dcmplx(e1(i-1),e2(i-1))
			do j=i,n
				m=10d0**(cdlog10(m0)+cdlog10(m1/m0)*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
				e1(j)=real(m)
				e2(j)=dimag(m)
				e1(j)=10d0**(log10(e1(i-2))+log10(e1(i-1)/e1(i-2))*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
				e2(j)=10d0**(log10(e2(i-2))+log10(e2(i-1)/e2(i-2))*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
			enddo
		endif
	else
! use the dielectric extrapolation, this is the default
		do j=i,n
			e1(j)=e1(i-1)
			e2(j)=e2(i-1)*grid(i-1)/grid(j)
		enddo
	endif

	deallocate(x)
	deallocate(y1)
	deallocate(y2)

	return
	end


	subroutine RegridDataLNK_file(input,grid,e1,e2,n,loglog,rho)
	use Tools
	IMPLICIT NONE
	REAL (KIND=dp) :: grid(n)
	REAL (KIND=dp) :: e1(n),e2(n),x0,y01,y02,x1,y11,y12,wp,gamma
	REAL (KIND=dp) :: rho
	integer i,j,n, n_line, iostatus
	character*100 input
	real (KIND=dp),allocatable :: x(:),y1(:),y2(:)
	integer n0,i0
	logical loglog
	COMPLEX (KIND=dp) :: m,m0,m1,cdlog10

	allocate(x(5000))
	allocate(y1(5000))
	allocate(y2(5000))

	!read data from file
	input = "DATA/" // input
	input = trim(input) // ".dat"
	IF (verbose) print*, input
	open(99,file=input)
	n0=0
	read(99,fmt=*) n0, rho
	do i=2, n0+1
		read (99, fmt=* ) x(i-1), y1(i-1), y2(i-1)
	end do
	close(99)


	i=1
	i0=1
	x0=x(i0)
	y01=y1(i0)
	y02=y2(i0)
	wp=(1d0-y01)/x0**2
	gamma=y02/x0**3
103	if(x0.ge.grid(i)) then
		e1(i)=1d0-wp*grid(i)**2
		e2(i)=gamma*grid(i)**3
		e1(i)=y01
		e2(i)=y02
		i=i+1
		goto 103
	endif
100	i0=i0+1
	if(i0.gt.n0) goto 102
	x1=x(i0)
	y11=y1(i0)
	y12=y2(i0)
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		e1(i)=y11+(grid(i)-x1)*(y01-y11)/(x0-x1)
		e2(i)=y12+(grid(i)-x1)*(y02-y12)/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y01=y11
	y02=y12
	goto 100
102	continue

	if(loglog) then
		m0=dcmplx(e1(i-1),e2(i-1))
		if(abs(m0).gt.2d0.and..false.) then
! don't use the conducting extrapolation since it is not very accurate
			do j=i,n
				m=m0*sqrt(grid(j)/grid(i-1))
				e1(j)=real(m)
				e2(j)=dimag(m)
			enddo
		else
! use loglog extrapolation
			m0=dcmplx(e1(i-2),e2(i-2))
			m1=dcmplx(e1(i-1),e2(i-1))
			do j=i,n
				m=10d0**(cdlog10(m0)+cdlog10(m1/m0)*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
				e1(j)=real(m)
				e2(j)=dimag(m)
				e1(j)=10d0**(log10(e1(i-2))+log10(e1(i-1)/e1(i-2))*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
				e2(j)=10d0**(log10(e2(i-2))+log10(e2(i-1)/e2(i-2))*log10(grid(i-2)/grid(j))/log10(grid(i-2)/grid(i-1)))
			enddo
		endif
	else
! use the dielectric extrapolation, this is the default
		do j=i,n
			e1(j)=e1(i-1)
			e2(j)=e2(i-1)*grid(i-1)/grid(j)
		enddo
	endif

	deallocate(x)
	deallocate(y1)
	deallocate(y2)

	return
	end


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
