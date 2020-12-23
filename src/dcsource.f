	subroutine  gauss(sgm,dt,ts,source,nso)
	   dimension source(*)
           sum=0.0
	   source(1)=0.0
	   nso=int(2.*ts/dt)
	   do 1 it=2,nso
	      cc1=dt*it-ts
	      cc=-2.*cc1*cc1/sgm/sgm
	      cc2=ts/sgm
	      cc=(exp(cc)-exp(-2.*cc2*cc2))/1.2533/sgm
	      source(it)=cc
	      sum=sum+(source(it)+source(it-1))*0.5
1	   continue
	   sum=sum*dt
	   do it=1,nso
	       source(it)=source(it)/sum
           enddo
	return
	end

	subroutine fault(theat,delt,ambd,s,a)
	    dimension a(*)

	    fm0=2./3.1415926 
            theata=(s-theat)*.0174533
            delta=delt*.0174533
            ambda=ambd*.0174533
            a(3)=fm0*(sin(2*theata)*cos(ambda)*sin(delta)+
     $                .5*cos(2*theata)*sin(ambda)*sin(2*delta))
            a(2)=fm0*(cos(theata)*cos(ambda)*cos(delta)-
     $                sin(theata)*sin(ambda)*cos(2*delta))
            a(1)=fm0*(.5*sin(ambda)*sin(2*delta))
            a(5)=fm0*(cos(2*theata)*cos(ambda)*sin(delta)-
     $                .5*sin(2*theata)*sin(ambda)*sin(2*delta))
            a(4)=fm0*(-sin(theata)*cos(ambda)*cos(delta)-
     $                cos(theata)*sin(ambda)*cos(2*delta))
	return
	end
      subroutine convt0(nx,ny,y,d,dt,n,l2n,nx1,
     *      ncent,fscl)
      dimension y(*),d(*)
      nz=max0(nx,ny)
      call log2fd(nz,n,l2n)
      l2n=l2n+1
      n=2*n
      nx1=nx+1
      ny1=ny+1
      do 3 i=1,ny
      j2=2*i
      j1=j2-1
      d(j1)=y(i)
3     d(j2)=0.
      do 4 i=ny1,n
      j2=2*i
      j1=j2-1
      d(j1)=0.
4     d(j2)=0.
      call colb(l2n,d,-1.)
      nhalf=n/2
      ncent=nhalf+1
      fscl=dt/float(n)
      return
      end
