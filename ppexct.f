	SUBROUTINE ppexct(rexact, n1, n2, t, r, nxyz)

	  common / num_par / np
	  common / num_tot / n
	  common / num_wall / nw
	  common / np_3 / np3

	  common / boxlen_xy / boxlen
	  common / boxlen_z / boxlenz

	  common / radius / a
	  common / delta_fn / delta

	  common / p_p_lub / fx, fy, g1, g2, g3

c    *************************************************************************
c    **  subroutine for calculating the exact two particle mobility tensor. **
c    **  (analytical solution as given by Batchelor, 1976)                  **
c    *************************************************************************
 
	  integer	n1, n2, i, j, k, l, alpha, beta, m2, m, nw, n, n3,np3,np

	  integer	nxyz(3, np)

	  double precision a, r(np3), rij(3), delta(6,6), X(2,2), Y(2,2), 
     +			       boxlen, rijsq, rijsep, rexact(6,6), s, m1, 
     +			       fx(0:11), fy(0:11), g1, g2, g3, t, one_4ssq,
     +			       two_pow_m, boxlenz

c    ** calculate separation between the two particles. define dimensionless **
c    **	separation 'u'                                                     **
	
	rexact = 0.d0
c
	do i=1, 3

	  rij(i) = r(nxyz(i,n2)) - r(nxyz(i,n1))
	  

c		rij(i) = rij(i) - anint (rij(i)/boxlenz) * boxlenz

c	  endif

	enddo

	rijsq = rij(1)**2 + rij(2)**2 + rij(3)**2
	rijsep = sqrt ( rijsq )
	s = rijsep / a

	if(s.le.2)then

c	  if((n1.lt.nw).or.(n2.lt.nw))write(*,*)t,n1,n2,s

	  s = 2.0 + 1d-8
	  rijsep = s * a
	  rijsq = rijsep ** 2

	endif
c
c
c
			one_4ssq = 1.0 - 4.0 / s**2

			X(1,1) = g1 / one_4ssq - g2 * log(one_4ssq)
     +			   - g3 * one_4ssq * log(one_4ssq) + fx(0) - g1

			do m2 = 1, 5
			  
			  m = 2 * m2

			  two_pow_m = 1.0 / real(2.0 ** m)

			  if(m.eq.2)then
				m1=-2d0
			  else
				m1=real(m-2)
			  endif

			  X(1,1) = X(1,1) + (two_pow_m*two_pow_m*fx(m)
     +				 - g1 - 2d0/real(m) * g2 + 4d0/real(m)/m1 * g3) 
     +				 * (2d0/s) ** m

			end do

			X(2,2)=X(1,1)
c
c
			X(1,2) = -2d0/s * g1 / one_4ssq 
     +			   - g2 * log((s+2.0)/(s-2.0))
     +			   - g3 * one_4ssq * log((s+2)/(s-2)) - 4d0 * g3 / s


			do m2=1,6
			  
			  m = 2 * m2 - 1

			  two_pow_m = 1.0 / (2.0 ** m)

			  if(m.eq.2)then
				m1=-2d0
			  else
				m1=real(m-2)
			  endif

			  X(1,2) = X(1,2) - ( two_pow_m*two_pow_m*fx(m)
     +				 - g1 - 2d0/real(m) * g2 + 4d0/real(m)/m1*g3 ) 
     +				 * (2d0/s)**m

			end do

			X(2,1)=X(1,2)
c
c
			Y(1,1) = -0.167d0 * log(one_4ssq) + fy(0)

			do m2=1,5
			  
			  m = 2 * m2

			  if(m.eq.2)then
				m1=-2d0
			  else
				m1=real(m-2)
			  endif

			  Y(1,1) = Y(1,1) + (1d0/real(2**m)*1d0/real(2**m)*fy(m)
     +				 - 2d0/real(m) * 0.167) * (2d0/s) ** m

			end do

			Y(2,2) = Y(1,1)
c
c
			Y(1,2) = -0.167 * log((s+2)/(s-2))

			do m2=1,6
			  
			  m=2d0*m2-1d0

			  if(m.eq.2)then
				m1=-2d0
			  else
				m1=real(m-2)
			  endif

			  Y(1,2) = Y(1,2) - ( 1d0/real(2**m)*1d0/real(2**m)*fy(m)
     +				 - 2d0/m*0.167 ) * (2d0/s)**m

			end do

			Y(2,1)=Y(1,2)
c
c    **	begin loops over particles and coordinates. i,j are particle indices. **
c    **   alpha, beta are particle coordinates                                  **
c

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

c    **   calculate the exact particle-particle mobility tensor   **

					  rexact(k+alpha,l+beta)=rexact(k+alpha,l+beta)+
     +				  X(i,j)*rij(alpha)*rij(beta)/rijsq
     +				 +Y(i,j)*(delta(alpha,beta)
     +				 -rij(alpha)*rij(beta)/rijsq)

c    **	end loops over particle coordinates   **

					end do
				  end do

c    **	end loops over particle indices   **

				end do
			  end do  

	  return
	  end