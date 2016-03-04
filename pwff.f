	SUBROUTINE pwff(pwinf, p1, p2, r,rimage, nxyz)

        common / num_tot / n
	  common / n_3 / n3
	  common / num_par / np
	  common / num_wall / nw
	  common / nw_3 / nw3
         common / np_3 / np3
	  common / radius / a
	  common / delta_fn / delta
	  	  
	  common / boxlen_z / boxlenz
	  common / boxlen_xy / boxlen
	  
c
c
	integer		np, nw, nw3, i, j, k, l, n, n3, pw, n1,p1,p2

	integer		nxyz(3, n)

	double precision  a, r(np3), rijd, delta(6,6), stklet, dipole, 
     +				  rij(3), pwinf(3, 3), boxlenz, boxlen, 
     +				  mw2b(nw3, nw3),rimage(np3),imagerij(3), imagedist,
     +  wallgreen

c	double precision, allocatable :: mpw(:, :)
	
c	allocate ( mpw(nw3+3, nw3+3) )
c
c	the particle-wall far-field resistances are added in these loops.
c	the procedure is as follows: a mobility tensor is calculated
c	for the free particle and all the wall particles. the mobility
c	tensor is inverted to get the resistance tensor. for calculating
c	particle velocities, only a 3x3 tensor for the single particle
c	is required
c
	pwinf = 0.0

	  do i = 1, 3

	    rij(i) = r(nxyz(i,p1)) - r(nxyz(i,np+pw))
            imagerij(i)=r(nxyz(i,p1))-rimage(nxyz(i,p2))
c		if(i.ne.3) rij(i) = rij(i) - anint(rij(i)/boxlen) * boxlen

	  enddo

	  rijd = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )

        imagedist = sqrt ( imagerij(1)**2 + imagerij(2)**2 + 
     +   imagerij(3)**2 )
        
        
	  do i = 1, 3

	    do j = 1, 3

		  stklet = delta(i, j) / rijd + rij(i) * rij(j) / rijd ** 3

		  dipole = delta(i, j) / rijd ** 3 
     +			 - 3.0 * rij(i) * rij(j) / rijd ** 5
     
       
                call calwalgnsub(i,j,r,rimage,nxyz,p1,p2,
     + imagerij,imagedist,wallgreen)
     
     
         pwinf(i, j) = pwinf(i, j)+delta(i,j)  
     + +   wallgreen
     
  
     
     
     

		 end do

	  end do

c	end do
c
c
	call invert (pwinf, 3)
c
c	
	
c
c	deallocate ( mpw )
	return
	end
