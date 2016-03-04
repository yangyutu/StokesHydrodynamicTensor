        subroutine sortds(D,r,nxyz,rbinnum)
        
        common / num_par / np
	  common / np_3 / np3
	  common / num_tot / n
	  	  common / n_3 / n3
	  common / radius / a
	  
	  integer, parameter:: nmax2=200
	  integer  n, np, n3, np3, nw3, i, j, k, l
	  integer nxyz(3, n),rbin,rbinnum
	  
	  
	  double precision D(np3,np3),a,xmean,ymean,dist,Dxhist(nmax2),
     +  Dyhist(nmax2),Dxhistcount(nmax2),Dyhistcount(nmax2),r(n3)  
	  
	  Dxhist=0.0
	  Dyhist=0.0
	  Dxhistcount=0
	  Dyhistcount=0
	  
	  do i=1,np
	  
	  dist=(r(nxyz(1,i)))**2+(r(nxyz(2,i)))**2
	  
	  dist=dist**0.5
	  
	  rbin=int(dist/a)+1
	  
	  Dxhist(rbin)=Dxhist(rbin)+D(nxyz(1,i), nxyz(1,i))
	  Dyhist(rbin)=Dyhist(rbin)+D(nxyz(2,i), nxyz(2,i))
	  Dxhistcount(rbin)=Dxhistcount(rbin)+1
	  Dyhistcount(rbin)=Dyhistcount(rbin)+1
	  
	  end do
	  
	  do i=1,rbinnum
	  if(Dxhistcount(i) .eq. 0) then
	  write(400,'(f12.5,f14.4,i5,3f12.5)') dble(i-0.5),a*(i-0.5),
     +   0, 0.0,0.0,0.0
	  else
	  Dxhist(i)=Dxhist(i)/Dxhistcount(i)
	  Dyhist(i)=Dyhist(i)/Dyhistcount(i)
	  write(400,'(f12.5,f14.4,i5,3f12.5)') dble(i-0.5), a*(i-0.5), 
     +	  int(Dxhistcount(i)), Dxhist(i),
     +   Dyhist(i), 
     +  0.5*(Dxhist(i)+Dyhist(i))
     
        end if
        
        end do
        
        end