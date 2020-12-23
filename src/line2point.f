
C********************************************************************C
C*                                                                  *C
C*                hy2sac.f                                          *C 
C*                                                                  *C
C*       Convert Green's function from hybrid method to SAC files   *C
C*                                                                  *C
C*   Written by: Lianxing Wen                                       *C
C*               Seismological Lab. Caltech                         *C
C*                                                                  *C
C*   Last modified: Oct. 22, 1996                                   *C
C*                                                                  *C
C********************************************************************C

       parameter(nf=130,ndd=20)
       dimension distx(2000),distz(2000),val(2000),val1(ndd)
       common/aaa/a(ndd,20000),b(20000)
       integer nnn(ndd)
       character*20 iname
       character*3 kstrn1
       character*8 kstrn2
	common/ff/pp(20000)

       integer OUT

          OUT  = 7
	  pi   = 3.1415926
	  ite  = 0
	  key  = 1
	  ndis2 = 1
          idev  = 0

	  call setpar()

	  call mstpar("accfile",       "s", iname)
	  call mstpar("Moment",        "f", amom)

	  call mstpar("Conv_Source",   "d", ite)
	  call mstpar("Dt1",           "f", dt1)
	  call mstpar("Dt2",           "f", dt2)
	  call mstpar("Dt3",           "f", dt3)
	  call mstpar("beta",          "f", rhos)
	  call getpar("point_source",  "d", key)
	  call mstpar("nt_kir",        "d", n)
	  call mstpar("dt_WKM",        "f", dt)

	  call getpar("nrecv",         "d", ndis2)
	  call getpar("finaldev",      "d", idev)
	  call endpar()

c         /* read the 2D output Green's function */
          open(OUT,file=iname,status='old')

          ndis1  =1
          do j=1,n
	      read(OUT,*)(a(i,j),i=ndis1,ndis2)
          enddo

c	  close(OUT)    
c         /* end read */

	  amom=(amom*(1.e-20))/(4.*pi*rhos)


          do i = ndis1, ndis2
	      
	      ifound = 0
	      nnn(i) = n 
	      do j=n,1,-1
		  if(abs(a(i,j)) .le. 1.e-30 .and. ifound .ne. 1) then
		      nnn(i) = j
	          else		
		      ifound = 1
                  endif
              enddo
	      nnn(i) = nnn(i) -1
		  
              if(key.eq.1) then
c                 /*convolve with 1/sqrt(t) */
                  write(*,*)"Convolving with 1/sqrt(t)"
                  call step(nnn(i),1,dt,i,i)
                  nnn(i)=nnn(i)/2
              endif
          enddo
          if(key.eq.1) then
               dt=dt*2.
          else 
               amom = 1.0
          endif
    
	  if(ite.eq.1)then
c             /* call source time function */	  
	      call stime(dt1,dt2,dt3,dt,pp,nfa)
	  endif    

	  xmin1=20000.
	  xman=0.0

	  integer0=ichar('0')
          do 200 id=ndis1,ndis2
	      n = nnn(id)
c             /* get the position of the receivers */
              integer_id = integer0+id

	      if(id.lt.10)kstrn1=char(integer_id)
	      if(id.ge.10.and.id.le.99)then
	          id10=id/10+integer0
	          id1=id-(id10-integer0)*10+integer0
                  kstrn1=char(id10)//char(id1)
	      endif

              kstrn2='syn.SAC.'
              iname=kstrn2//kstrn1
              if(distx(id).gt.xman)xman=distx(id)
              if(distx(id).lt.xmin1)xmin1=distx(id)

c              write(*,*)'n=',n
              do i=1,n
    		 b(i)=a(id,i)
              enddo

	      if(key .eq. 1 .and. idev .eq. 1) then

c                 /* take derivative */
		  call diff(n,b,dt)
c
c                 modify with 1/sqrt(r)

       	          xs  = 0.
 	          cc1 = distx(id)-xs
	          cc2 = zs-distz(id)
	          roo = sqrt(cc1*cc1+cc2*cc2)
	          do i=1,n
	              cc   = 1.4142/sqrt(roo*abeta)
 	              b(i) = b(i)*cc
                  enddo
	      endif

c              write(*,*)'n1=',n
	      do i=1,n
 	          b(i) = b(i)*amom
              enddo

	      if(ite.eq.1)then
c                 /* convolve source time function */	  
                  write(*,*)"Convolving with source time function"
	          call convt(b,n,pp,nfa,b,n,dt)
	      endif    

c              write(*,*)'n2=',n
              do i=1,n    
                  a(id,i)=b(i)
c		  write(*,*)'b(i)=',i,' ',b(i)
              enddo
2	 continue

	      call WSAC1(iname,b,n,tstart,dt,NERR)
	      call SETFHV('DEPMIN',ymin,NERR)
	      call SETFHV('DEPMAX',ymax,NERR)
	      call SETFHV('DEPMEN',sum,NERR)
	      call SETFHV('B',tstart,NERR)
	      call SETFHV('DELTA',dt,NERR)
	      distance = distx(id)
	      call SETFHV('DIST',distance,NERR)
	      GCARC = distx(id)/111.
	      call SETFHV('GCARC',GCARC,NERR)
	      call SETNHV('NPTS',n,NERR)

 200	  continue

	stop
	end

      subroutine colb(nn,datai,signi)
      dimension datai(*)
      n=2**(nn+1)
      j=1
      do 5 i=1,n,2
      if(i-j)1,2,2
    1 tempr=datai(j)
      tempi=datai(j+1)
      datai(j)=datai(i)
      datai(j+1)=datai(i+1)
      datai(i)=tempr
      datai(i+1)=tempi
    2 m=n/2
    3 if(j-m)5,5,4
    4 j=j-m
      m=m/2
      if(m-2)5,3,3
    5 j=j+m
      mmax=2
    6 if(mmax-n)7,10,10
    7 istep=2*mmax
      theta=signi*6.28318531/float(mmax)
      sinth=sin(theta/2.)
      wstpr=-2.0  *sinth*sinth
      wstpi= sin(theta)
      wr=1.
      wi=0.
      do 9 m=1,mmax,2
      do 8 i=m,n,istep
      j=i+mmax
      tempr=wr*datai(j)-wi*datai(j+1)
      tempi=wr*datai(j+1)+wi*datai(j)
      datai(j)=datai(i)-tempr
      datai(j+1)=datai(i+1)-tempi
      datai(i)=datai(i)+tempr
    8 datai(i+1)=datai(i+1)+tempi
      tempr=wr
      wr=wr*wstpr-wi*wstpi+wr
    9 wi=wi*wstpr+tempr*wstpi+wi
      mmax=istep
      go to 6
   10 return
      end

      subroutine convt(x,nx,y,ny,z,nz,dt)
      dimension x(*),y(*),z(*),c(80200),d(80200)
      nz=max0(nx,ny)
      call log2fd(nz,n,l2n)
      l2n=l2n+1
      n=2*n
      nx1=nx+1
      ny1=ny+1
      do 1 i=1,nx
      j2=2*i
      j1=j2-1
      c(j1)=x(i)
1     c(j2)=0.
      do 2 i=nx1,n
      j2=2*i
      j1=j2-1
      c(j1)=0.
2     c(j2)=0.
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
      call colb(l2n,c,-1.)
      call colb(l2n,d,-1.)
      nhalf=n/2
      ncent=nhalf+1
      do 6 i=1,ncent
       j2=2*i
      j1=j2-1
      e1=c(j1)*d(j1)-c(j2)*d(j2)
      e2=c(j1)*d(j2)+c(j2)*d(j1)
      c(j1)=e1
6     c(j2)=e2
      call conj(c,ncent)
      call colb(l2n,c,1.)
      fscl=dt/float(n)
      do 7 i=1,nz
      j1=2*i-1
7     z(i)=c(j1)*fscl
      return
      end
      subroutine conj(c,ncent)
      dimension c(*)
      nhm1=ncent-2
      ncent2=ncent*2
      do 1 i=1,nhm1
      l2=2*i
      l1=l2-1
      l3=l2+1
      k1=ncent2+l1
      k2=ncent2+l2
      j1=ncent2-l3
      j2=ncent2-l2
      c(k1)=c(j1)
1     c(k2)=-c(j2)
      return
      end
      subroutine log2fd(np,n,l2n)
      n1=0
      n=1
      l2n=0
1     continue
      if((np.gt.n1).and.(np.le.n)) go to 2
      n1=n1*2
      n=n*2
      l2n=l2n+1
      go to 1
2     return
      end
      subroutine step(nfa,nfad,dp,ndis1,ndis2)
       parameter(nf=130,ndd=20)
       common/aaa/rp(ndd,20000),b(20000)
	common/ff/pp(20000)
	dimension f(20000),p(20000)

        call fa(dp,nfa,f)
        nff=nfa/2
	do 3 k=ndis1,ndis2
	l=0
	do 41 j=1,nfa
 41	pp(j)=rp(k,j)
       do 20 n=1,nff,nfad
      l=l+1
      p(l)=convs(dp,nfa-1,n-1,f)
      rp(k,l)=p(l)
 20     continue
 3       continue
	lz=l
       return
      end

	subroutine stime(dt1,dt2,dt3,dt,f,nf)
	dimension f(20000)

	    n1  = int(dt1/dt+0.1)
	    n2  = int(dt2/dt+0.1)
	    n3  = int(dt3/dt+0.1)
	    n12 = n1+n2+1
	    n   = n1+n2+n3+1

	    nf =n

	    sum=0.
	    f(1)=0.0
	    do i=2,n
	         k=i
	        f(k)=1.
	        if(i.lt.(n1+1))f(k)=(i-1.0)*(1./n1)
	        if(i.gt.n12)f(k)=1.0-(i-n12)*(1./n3)
 	        sum = sum+(f(k)+f(k-1))*0.5
            enddo

	    sum=sum*dt
	    do i=1,nf
 	        f(i)=f(i)/sum
            enddo
	return
	end

      subroutine fa(del,j2,f)
	dimension f(20000)
      do 5 j=2,j2   
      f(j)     = 1./sqrt((j-1)*del)    
 5      continue  
      f(1)     = (11.-4.*2.**.5)/sqrt(2.*del)  
      return
         end
      function convs(del,nf,n,fp)
	common/ff/fa(20000)
	dimension fp(20000)
      nn=n
      dn=del
      if(nn.lt.1) go to 2
      ndo=min0(nn,(nf-1)/2)
      ip=2
      np=2*nn
      even=fp(ip)*fa(np)
      odd=0
      if(ndo.lt.2) go to 11
      do 10 i=2,ndo
      ip=ip+1
      np=np-1
      odd=odd+fp(ip)*fa(np)
      ip=ip+1
      np=np-1
      even=even+fp(ip)*fa(np)
  10      continue
 11     continue
      ends=fp(1)*fa(2*nn+1)+fp(ip+1)*fa(np-1)
      convs=dn*(ends+4.*even+2.*odd)/3.
	go to 5
 2      convs=0.
 5	continue
      return
      end 

