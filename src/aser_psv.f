
      subroutine GRT(isourcec,iflat,ispc,mtdc,nbc,jo0,cc0,
     *               ss0,dd0, tth0,
     *               lfinal0,nen0,na0,nray0,ncoun0,
     *               dpc,nnc,xxx,tsc,nstressc,so,der,green)

c   isourcec The source type =0: isotropic source
c                            =1: dislocation source
c   iflat = 0       flat layer 1: flatting
c   nbc      The layer where the source is
c   jo0
c   cc0
c   ss0
c   dd0
c   tth0
c   dpc      The time interval
c   nnc      The number of time steps
c   xxx      The epicentral distance
c   tsc      The starting time
c   so       The radiation pattern
c   nstressc The control of output Green's function
c            For dislocation source:
c               =0   w    velocity
c               =12  Txz  stress
c               =22  Tzz  stress
c               =11  Txx  stress
c            For isotropic source:
c               <=1  The Green's function
c               > 1  The z-deravitive of Green's function
c   green    The output of the Green's function


      include 'psv.p'

      integer nbc,jo0,lfinal0,nnc,nstressc,isourcec
      integer ispc, der, mtdc
      real dpc,xxx,tsc
      real so(*), green(*)
      real cc0(*),ss0(*),dd0(*),tth0(*)
      integer nen0(*),na0(*),nray0(*),ncoun0(*)

      common/trace/tracing
      common/down/ndown, nstress
      common/pathc/po,to,k
      common/flt/phh(iflt,inz),ns,beta,nz,nup,nst
      common/str/ph(itt,inz)
      common/orst/cc(icc),ss(icc),dd(icc),tth(icc),xx
      common/rays/na(ina1,ina2),nray(ina1,ina2),nend,lmax,nsp,love,
     *nzr,iwh
      common/tfix/tn1,tn2,tn3,tn4,jn1,jn2,jn3,jn4,det,ncont
      common/term/iordr
      common/coeff/itrnm
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
      common/thy/tt(itt),pp(itt),ff(itt),l
      common/plotc/con
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
      common/cagcon/tstart,dp,nn,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
      common/soray/nsorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
      common/dscfil/far1
      common/tmpwr/iwr
      common/sourcetype/isource
      common/scase/ispecial

      dimension nen(ina2),ncoun(ina2)

      common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prnts,prntc,flat,prcon,prampl

c DLTM and mtd used in the contor
c dltm max seperation between time points
c mtd=3 controls rate of approaching this limit
c Dltp the last time point before TO must be less than
c Dtim the seperation between TO and the first time point after TO
c muse be less than dtim
	iwr=0
	nfm=0

	prnt=.false.
	prnts=.false.
	prntc=.false.
	prcon=.false.
	prampl=.false.
	tracing = -1

	do j=1,nnc
	    green(j)=0.0
        enddo

c    isource=1 dislocation source
c           =0 explosition source
        isource=isourcec
        nplnw=2
        itrnm=2
	iwh=0
	iordr=1
	ndirt=2

        nb =nbc
	dp =dpc
	nstress=nstressc
	tstart=tsc
        jo=jo0
	dt=dp
	ispecial = ispc

	nn=nnc
        depth=0.0
        do jj=1,nbc
            depth=depth+tth0(jj)
	enddo
	sdepth=depth
	rhos=dd0(nbc)

	do j=1,jo0
	    cc(j)=cc0(j)
	    ss(j)=ss0(j)
	    dd(j)=dd0(j)
	    tth(j)=tth0(j)
        enddo

        ns=nb
        ncont = 2

	flat = .true.
	if(iflat .eq. 0)flat=.false.

       lfinal=lfinal0
       kk=1
       do  nsp=1,lfinal0
	   nen(nsp)=nen0(nsp)
 	   ncoun(nsp)=1
	   do j=1,nen(nsp)
	       na(j,nsp)=na0(kk)
	       nray(j,nsp)=nray0(kk)
               kk=kk+1
	   enddo
        enddo

	do i=1,jo
	   cc(i)=cc0(i)
	   ss(i)=ss0(i)
	   dd(i)=dd0(i)
	   tth(i)=tth0(i)
        enddo

        if(tracing .gt. 0) then
            write(*,*)'nb=',nb,' dp=',dp,' nstress=',nstress
            write(*,*)'tstart=',tstart,' jo=',jo,' nn=',nnc
 	    write(*,*)"jo=",jo
 	    do j=1,jo
 	     write(*,*)cc(j), ss(j), dd(j), tth(j)
 	    enddo
 	    write(*,*)"lfinal=",lfinal
 	    do nsp=1,lfinal
 	       write(*,*)nen(nsp),(na(j,nsp),j=1,nen(nsp))
 	       write(*,*)ncoun(nsp), (nray(j,nsp),j=1,nen(nsp))
 	   enddo
        endif

	tracing = -1
      jpt = 0
      nz=8
      nconjt = 0
      lmax = jo+1
      call curay(jo)
      ll=lmax
      dtim = dp/5.
      dtim = 0.001
      dltp = dtim
c      dtim = dp/5.
      dltm = dp*5.
      dltm =5.0
c Changed temporatorly
      mtd=mtdc
      nnn=15
      ns=nz
      beta=s(nb)

      itry = 0
c     /* quick and dirty fix for small distances */
      xmin = 8
112   continue
      x=xxx
c      write(*,*)"x=",x
      x=abs(x)
      if(x .lt. 0.1) x = 0.1
      far1=0.0
      xx=x
      con=1.
c           con=sqrt(2./xx)
      hs =0.
      l=1
      nnsp=1
      do  nsp=nnsp,lfinal

	   ndown=ncoun0(nsp)
c comments, the nrec is the layer where the receivers are,
c it needs caution, when an interface is present.
           nrec=na(nen(nsp)-1,nsp)

c          comment by Lianxing Wen ->testing
           if(nray(nen(nsp)-1,nsp) .eq. 5)then
c	       ispecial = 1
           else
c	       ispecial = 0
           endif
c          /* parameters at receiver */
           tqc=cc(nrec)
           tqs=ss(nrec)
           tqd=dd(nrec)
c          write(*,*)tqc,tqs,tqd

           love=0
           if(nray(1,nsp).eq.4) love=2

           nst=nray(1,nsp)

           nsorce=nray(1,nsp)
           if(nsorce.eq.4) nsorce=3

           nup=0
           if(na(1,nsp).gt.nb) nup=2

           ndiret=0
           ncount=ncoun(nsp)
           nend=nen(nsp)-1
c	   nend=nen(nsp)

           nx=1
           do nm=1,nend
               nx=max0(nx,na(nm,nsp))
           enddo

           lmax=nx
           tq=0
           call trav(ncount)
           if (love.eq.2) nray(1,nsp)=4
      enddo

      if(abs(xxx) .lt. xmin .and. itry .eq. 0) then
          to0 = to
          xxx = xmin
          itry = 1
          nn  = nnc + 30
          goto 112
      endif

      ishift = 0
      if(itry .eq. 1) ishift = (to - to0) / dp

      ymax=0.0
      jz=5
      if(nstress.eq.0 .or. nstress.eq.12)jz=8
      do jk=1,nn
          pp(jk)=so(1)*ph(jk,jz)+so(2)*ph(jk,jz-1)+so(3)*ph(jk,jz-2)
	  if(abs(pp(jk)).gt.ymax)ymax=abs(pp(jk))
      enddo

       if(der .eq. 0) goto 3456
       ymax=0.0
       y1=pp(1)
       pp(1)=0.0
       do jk=2,nn
	    y2=pp(jk)
	    pp(jk)=(pp(jk)-y1)/dp
	    y1=y2
	     if(abs(pp(jk)).gt.ymax)ymax=abs(pp(jk))
        enddo

3456    continue
       write(*,*)'ymax=',ymax

	    do jk=1,nnc
	         green(jk)=pp(jk+ishift)
            enddo

527    continue
      return
      end
