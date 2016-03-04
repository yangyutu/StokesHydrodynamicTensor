	SUBROUTINE ludcmp(a,np,indx,d)
c

c
c	See Numerical Recipes for more description on this subroutine
c

c
      INTEGER np,indx(np),n
	double precision d
	double precision, PARAMETER :: TINY = 1D-20
	INTEGER i,imax,j,k
	double precision aamax,dum,sum 
	double precision a(np, np) 

	double precision, allocatable :: vv(:)

	allocate ( vv(np) )

	n=np
	d=1d0

	do i=1,n 
	  aamax=0d0
	  do j=1,n
	    if (abs(a(i,j)).gt.aamax)aamax=abs(a(i,j))
	  enddo
	  if (aamax.eq.0d0) stop 'singular matrix in ludcmp' 
	  vv(i)=1.d0/aamax 
	enddo
	do j=1,n 
	  do i=1,j-1 
	    sum=a(i,j)
	    do k=1,i-1
	      sum=sum-a(i,k)*a(k,j)
	    enddo
	    a(i,j)=sum
	  enddo
	  aamax=0
	  do i=j,n 
	    sum=a(i,j)
	    do k=1,j-1
	      sum=sum-a(i,k)*a(k,j)
	    enddo
	    a(i,j)=sum
	    dum=vv(i)*abs(sum)
	    if (dum.ge.aamax)then 
	      imax=i
	      aamax=dum
	    endif
	  enddo
	  if (j.ne.imax)then 
	    do k=1,n 
	      dum=a(imax,k)
	      a(imax,k)=a(j,k)
		  a(j,k)=dum
		enddo
		d=-d 
		vv(imax)=vv(j) 
	  endif
	  indx(j)=imax
	  if(a(j,j).eq.0d0)a(j,j)=TINY	
	  if(j.ne.n)then 
	    dum=1d0/a(j,j)
		do i=j+1,n
		  a(i,j)=a(i,j)*dum
		enddo
	  endif
	enddo 

	deallocate ( vv )

	return
	END