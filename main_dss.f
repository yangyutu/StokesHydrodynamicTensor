	PROGRAM SBD

	  common / num_par / np
	  common / num_tot / n
	  common / num_wall / nw
	  common / np_3 / np3
	  common / n_3 / n3
	  common / nw_3 / nw3

	  common / boxlen_xy / boxlen
	  common / boxlen_z / boxlenz

	  common / radius / a
	  common / delta_fn / delta
	  common / brsh / brush_length

	  common / selftm_vol / ms, vol
	  common / beta_parm / zi

	  common / twopi_bxlenx / tpi_bxln
	  common / twopi_bxlenz / tpi_bxlnz

	  common / p_p_lub / fx, fy, g1, g2, g3

	  common / r_cut / rcut, rmax
	  common / rr / rr
	  common / k_cut / kcut, kmax
	  common / kk / kk


c    ******************************************************************
c    ** This program runs the actual simulation and writes           **
c    ** coordinates into an output file.                             **
c    ******************************************************************
	  
	  integer, parameter :: nstepmax = 10000000
	  
	  integer          i, j, k, l, np, nw, n3, np3, n, nw3, dum
	  
	  integer		   nstep, iprint, istart,fileindex,filestart,fileend,
     +   skip

	  integer, allocatable :: nxyz(:, :)

	  double precision, parameter :: pi = 3.1416, sqpi = 1.7725, 
     +								 twopi = 6.2831853, tol=1d-5

	  double precision, allocatable :: D(:,:), rgrnd(:,:),
     +								   rmove(:),
     +								   rinfi(:,:), mwall(:,:),
     +								   mw2b(:,:), r(:)

	  double precision t, delta(6,6), a, tempr, dt, 
     +				   boxlen, zi, ms(3,3), brush_length, dt_m,
     +				   m_2, fac1, fac2, istart_dt, t_min_isdt, vol, 
     +				   tpi_bxln, fx(0:11), fy(0:11), g1, g2, g3, 
     +				   boxlenz, tpi_bxlnz, phi, deltat, acore, 
     +				   Dssx, Dssy, Dssz, Dssxt, Dssyt, Dsszt,
     +  xmean,ymean,proheight

        character        par_in*30, wall_in*150, dssout*30, system*2,
     +  filepath*160,cnfile*160,opfile*160,dsfileout*150,filetemp*3,
     +  filein*160,opfiletemp*160,rdistfile*160,xyzanddsfile*160  

	  double precision p, rcut, kcut

	  integer rbinnum

	open(1,file='run.txt')
c	read(1,*)
c	read(1,*)nw
c	read(1,*)
c	read(1,'(a150)') filepath
	read(1,*) 
	read(1,'(a150)')cnfile
	read(1,*)
	read(1,'(a150)')opfile
c	read(1,*) 
c	read(1,'(a150)')wall_in
	read(1,*)
	read(1,*) nframe
	read(1,*)
	read(1,*)acore
	read(1,*) 
	read(1,*) filestart,fileend
	read(1,*)
	read(1,*) skip
      read(1,*)
      read(1,*) rbinnum
      read(1,*) 
      read(1,*) proheight
      read(1,*)
	read(1,*) np
        
        do fileindex=filestart,fileend

        if(fileindex .lt. 10) then
        write(filetemp,'(i1)') fileindex
        else if(fileindex .lt. 100) then
        write(filetemp,'(i2)') fileindex
        else
        write(filetemp,'(i3)') fileindex
        end if
        
      	  filein=trim(cnfile)//trim(filetemp)//'.txt'
	
	  opfiletemp=trim(opfile)//trim(filetemp)//'.txt'
	  dsfileout='opandds'//trim(filetemp)//'.txt'
	  rdistfile='rdistrids'//trim(filetemp)//'.txt'
        xyzanddsfile='xyzandds'//trim(filetemp)//'.txt'
		write(*,*) filein
		write(*,*) opfiletemp

        	  open ( 200, file = filein )
	  open ( 300, file = dsfileout )
        open (211, file = opfiletemp)
        open(400,file=rdistfile)
        open(500,file=xyzanddsfile)
        	 
	  
	  

	  n = np + nw	  	
	  np3 = np * 3
	  nw3 = nw * 3
	  n3 = np3 + nw3

	  allocate ( D(np3, np3), rgrnd(np3, np3), 
     +			 rmove(np3), 
     +			 rinfi(np3, np3), mwall(nw3, nw3), 
     +			 mw2b(nw3, nw3), r(n3), nxyz(3, n) )

	  do i = 1, n
		
	    nxyz(1, i) = 3 * (i-1) + 1
	    
		nxyz(2, i) = nxyz(1, i) + 1

		nxyz(3, i) = nxyz(2, i) + 1

	  end do  


	  a = acore



	  do i=1,6
	   do j=1,6
	     delta(i,j)=0d0
	 	 if(i.eq.j)delta(i,j)=1d0
	   end do
	  end do



c
c	compute the self term (needs to be calculated only once)
c
	fx(0)=1.d0
	fx(1)=3.d0
	fx(2)=9.d0
	fx(3)=19.d0
	fx(4)=93.d0
	fx(5)=387.d0
	fx(6)=1197.d0
	fx(7)=5331.d0
	fx(8)=19821.d0
	fx(9)=76115.d0
	fx(10)=320173.d0
	fx(11)=1178451.d0

	fy(0)=1.d0
	fy(1)=1.5d0
	fy(2)=2.25d0
	fy(3)=7.375d0
	fy(4)=29.063d0
	fy(5)=70.594d0
	fy(6)=230.391d0
	fy(7)=700.336d0
	fy(8)=2226.629d0
	fy(9)=8694.131d0
	fy(10)=32889.478d0
	fy(11)=13034.1382d0

	g1=0.25
	g2=0.225
	g3=0.0268
c
               
c
c	since the wall particles are always fixed, we can calculate
c	mobility terms of wall-wall interactions beforehand and use
c	them throughout the simulation
c
	
c
c	for the far-field multi body calculations (ewald summed)
c
c	  call ffwall(mwall, r, nxyz)
c
c	for the far-field particle-wall calculations
c
c	  call ffw2b(mw2b, r, nxyz)
c	
c    **	Open files   **
c      

c	initialize arrays
c
	  do j=1,np3
		  
		  rmove(j) = 0.d0
	    
	  end do
c
c    **   Now start the dynamics   **
c 
	  Dssxt = 0.0
	  Dssyt = 0.0
	  Dsszt = 0.0

	  do 50 l=1, nframe
	  
	      read(211,*,end=1000) dum,c6, rg,psi6

c    ******************************************************************
c    ** call subroutines for calculating diffusion coefficients,     ** 
c    ** forces and random displacement terms                         **
c    ******************************************************************
	      xmean=0.0
            ymean=0.0
	    do i = 1, np
            
		read(200,*,end=1000)dum,r(nxyz(1, i)),r(nxyz(2, i)),r(nxyz(3,i))
    						       
		r(nxyz(1, i))=r(nxyz(1, i))*acore
		r(nxyz(2, i))=r(nxyz(2, i))*acore
		r(nxyz(3, i))=r(nxyz(3, i))*acore
		  xmean=xmean+r(nxyz(1, i))
 
		  ymean=ymean+r(nxyz(2, i))
c		  r(nxyz(3, i)) = r(nxyz(3, i)) * acore
                r(nxyz(3, i)) = proheight+acore

	    end do
	    xmean=xmean/dble(np)
	       ymean=ymean/dble(np)
	      do i = 1, np
	    r(nxyz(1, i)) = r(nxyz(1, i))-xmean 
		  r(nxyz(2, i)) = r(nxyz(2, i))-ymean
		    end do 
	  if( mod(l,skip) .eq. 0) then     
		
		call grndrm(rinfi, rmove, rgrnd, t, r, nxyz)

	    do j=1,np3
            
	      do k=1,np3

		    D(j,k)=rgrnd(j,k)
		  
		  end do
c	    write(*,*) d(j,j)
		end do
c
		call invert(D, np3)
c	do i=1,np
c	write(*,*)i
c	write(*,*)D(nxyz(1,i):nxyz(3,i),nxyz(1,i):nxyz(3,i))
c	pause
c	end do

		Dssx = 0
		Dssy = 0
		Dssz = 0

		do j=1,np

		  Dssx = Dssx + D(nxyz(1,j), nxyz(1,j))
		  Dssy = Dssy + D(nxyz(2,j), nxyz(2,j))
		  Dssz = Dssz + D(nxyz(3,j), nxyz(3,j))
		  
		  write(500,'(2i5,4f14.4)' ) l,np,r(nxyz(1, j)), r(nxyz(2, j)), 
     +	  D(nxyz(1,j), nxyz(1,j)), D(nxyz(2,j), nxyz(2,j))

		end do

		Dssx = Dssx / dble(np)
		Dssy = Dssy / dble(np)
		Dssz = Dssz / dble(np)



	 write(300, '(i5, 7f14.4)') l, c6, rg,psi6,Dssx, Dssy,
     + (Dssx+Dssy)/2.0, Dssz
     
        call sortds(D,r,nxyz,rbinnum)

		

c    **	end loop over time steps   **
c
       end if

50	  continue
c
c
1000	  close(1)
	  close(200)
	  close(300)
        close(211)
        close(400)
        close(500)
	  deallocate ( D, rgrnd, rmove, rinfi, mwall, mw2b, r, nxyz )
        end do
	  stop
	  end
c
c
