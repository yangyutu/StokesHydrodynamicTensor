	SUBROUTINE ppinf(mabinf, p1, p2, r,rimage, nxyz)

        common / num_tot / n
	  common / n_3 / n3
        common / np_3 / np3
        common / num_par / np
	  common / boxlen_xy / boxlen
	  common / boxlen_z / boxlenz

	  common / radius / a
	  common / delta_fn / delta

c    *****************************************************************************
c    ** subroutine to calculate the far-field particle-particle mobility tensor **
c    ** by the multipole expansion of the Green's function                      **
c    *****************************************************************************

	  integer	n1, n2, i, j, alpha, beta, k, l, n, n3,np,np3

	  integer	nxyz(3, np),p1,p2

	  double precision		r(np3), a, mabinf(6,6), d, delta(6,6), 
     +						stklet, dipole, rij(3), boxlen, boxlenz,
     + imagestklet,	imagerij(3),imagedist,imagedipole,rimage(np3),
     + wallgreen		   

c    **	calculate particle-particle separation   **
	  
	  mabinf = 0.d0
	  n1=p1
	  n2=p2

	  do i=1, 3
	   rij(i) = r(nxyz(i, p1)) - r(nxyz(i, p2))
          imagerij(i)=r(nxyz(i,p1))-rimage(nxyz(i,p2))

	  enddo

	  d=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
         imagedist=sqrt(imagerij(1)**2+imagerij(2)**2+imagerij(3)**2)
        if(d.le. 2*a .and. n1 .ne. n2) then
        d=2*a+1
        end if

c    ** run loops over particles and particle coordinates. i,j denote particle **
c    ** indices. alpha and beta denote particle coordinates                    **

	  do i=1,2
	    if(i.eq.1)then
	      k=0
	    else
	      k=3
	    endif
	    do j=1,2
		  if(j.eq.1)then
	        l=0
		  else
		    l=3
		  endif
	      do alpha=1,3
	        do beta=1,3
		      if(i.ne.j)then
			  			  		
			  			  		
			  			  					    
			    stklet=delta(alpha,beta)/d+rij(alpha)*rij(beta)
     +                      /d**3
			    dipole=delta(alpha,beta)/d**3-3*rij(alpha)
     +										   *rij(beta)/d**5
     
c        	      imagestklet=delta(alpha,beta)/imagedist+
c     +      imagerij(alpha)*imagerij(beta)
c     +                      /imagedist**3
     
c          imagedipole=delta(alpha,beta)/imagedist**3-3*imagerij(alpha)
c     +        *imagerij(beta)/imagedist**5
     
c                       call calwalgnsub(alpha,beta,r,rimage,nxyz,p1,p2,
c     + imagerij,imagedist,wallgreen)
                

			  else
			    stklet=0d0
			    dipole=0d0
c			    wallgreen=0.0

			  endif

		  	  mabinf(k+alpha,l+beta)=delta(i,j)*delta(alpha,beta)
     +		  +3.0*a/4.0*(stklet+dipole*2.0/3.0*a**2)
c     +        + wallgreen
     
     		if( i .gt. 10) then
		write(*,*) 'h'
		end if
		    end do
		  end do
	    end do
	  end do

	  call invert (mabinf, 6)

	  return
	  end
