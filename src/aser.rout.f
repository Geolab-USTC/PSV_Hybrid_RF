
c********************************************************************c
c*       The subroutines used for GRT calculation                   *c
c*                                                                  *c
c*       Modify History:                                            *c
c*            Feb 27, 1997, Lianxing Wen   Add receiver functions   *c
c*                                         for grad(div) of P wave  *c
c*                                                                  *c
c********************************************************************c


       subroutine high(k)
       include 'psv.p'
       common/trace/tracing
       common/tstar/tq1,tvl,qa(icc),qb(icc)
         common/flt/phh(iflt,inz),ns,beta,nz,nup,nst
        common/rays/na(ina1,ina2),nray(ina1,ina2),ndnd,lmax,nsp,love,
     *nzr,iwh
        dimension pr(imgc),pi(imgc)
       common/soray/nsorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
        common/fixp/ddn(icc),arn(icc)
       common/orst/cc(icc),ss(icc),dd(icc),tth(icc),xx
c     dimension space(200)
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
      common/cagcon/tstart,dp,nn,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
      common/exact/phi(iexa),td(iexa),tw(iexa),nend,nm
      common/magic/pp(imgc),ddpt(imgc),tt(imgc)
       common/pathc/po,to,kk
      common/spe/derlp(ispe),dd1,dd2,dd3,dd4,no
      common/tfix/tn1,tn2,tn3,tn4,jn1,jn2,jn3,jn4,dett,ncont
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prntc,prnts,flat,prcon,prampl
      common/dscfil/far1
      common/scase/ispecial
       complex pp
c     real*8 tt,td,to,po,ddn,arn,s2,v2,p,ttp,tg,tc,rg,rcsq,rssq,xm
c     real*8 c,s,d,th,x
c     dimension e(100)
      complex ddpt
c       real*8 drc,drs
c
c     high is called by subroutine setup
        if(tracing .gt. 0)then
	  write(*,*)' high is being called'
	endif  
c
        km=k
      kk=k
      tmx=nn*dp
c     looking for branch cut nearest po,or highest velocity in t(p)
      xm=0.
      do 17 j=1,k
      if(rlp(j).le.0.001) go to 28
       if(c(j).gt.xm) xm=c(j)
      go to 17
 28     if(s(j).gt.xm) xm=s(j)
 17     continue
       if((k.eq.1).and.(nsorce.eq.5)) xm=c(2)
       if((k.eq.1).and.(nsorce.eq.3)) xm=s(2)
      s2=s(k+1)
c      k  is max layer number used in ray spec
c      Modified by Lianxing Wen
       iref = 0
       do j=1,nnend
	   if(na(j,nsp) .eq. lmax) iref = iref + 1
       enddo

        v2=0.
      if(love.eq.2) go to 62
       do 61 j=1,k
 61      v2=amax1(v2,c(j))
	if(lmax.gt.nb .and. iref .gt. 1) v2=amax1(v2,c(k+1))
      go to 663
  62  do 662 j=1,k
       v2=amax1(v2,s(j))
  662  continue
      if(lmax.gt.nb .and. iref .gt. 1) v2=amax1(v2,s(k+1))
 663  continue

      v1 = 1.e5
      if(love .eq. 2)then
c        /* CMB only */
	  do j=k-4,k
	      v1 = amin1(v1,s(j))
          enddo
      else
	  do j=k-4,k
	      v1 = amin1(v1,c(j))
	      if(s(j) .gt. 0.1) v1 = amin1(v1,s(j))
          enddo
      endif

      tn1=.8
      tn2=.4
      tn3=.1
       tn4=.001
       jn1=30
       jn2=30
       jn3=45
       jn4=100
        jn3=30
       jn4=30
 333  format(1x,6i10)
      del=1./xm
 81      q=-1.e-6
      det=1.e+12
c      computes (po) and (to)
       call find2(q,kk,del,det,po,to)
 222   format(5x,' to= ',e18.8,5x,' po=',e18.8)
      write(*,222) to,po
c     if(far1.eq.0.0.or.far1.gt.to) far1=to
      if (nsp.eq.3) far1=to
      if(ispecial .ne. 1) then
          p=1./v2
      else
          p=1./v1
      endif

      rg=abs(p-po)
      nk=2
      if(tracing.gt.0)then
          write(*,*)'in high calling herlp'
      endif    
      call herlp(k,p,ttp,dtp)
      if(tracing.gt.0)then
          write(*,*)'out calling herlp'
      endif    
      tc=ttp
       tg=to-ttp
      tg=abs(tg)
      mo=2
c     testing location of branch cut, looking for first motion
      if(po.le.1./v2) go to 6
c	temporary
      if(tg.gt.tn1) go to 6
      jn=jn1
      if(tg.gt.tn2) go to 18
      jn=jn2
      if(tg.gt.tn3) go to 18
      jn=jn3
   18 qz = rg/(jn+1)
      do 15 j = 1,jn
      derlp(j) = qz
 15     continue
      no=jn
      if(tg.lt.tn4) go to 2
      go to 19
c    nk= power of sin function
c       kn and n are not used
 6      continue
c    nnn=no. of points between pc and  po
c    the following if statements are changes made by t.wallace and make
c    a noticeable change in greens fns. output from previous runs when
c    nnn was a constant.
      nnn=15
       if((to-tc).gt.30.) nnn=60
       if((to-tc).gt.40.) nnn=65
       if((to-tc).gt.70.) nnn=70
       if((to-tc).gt.100.) nnn=85
       if((to-tc).gt.120.) nnn=90
       if(ispecial .eq. 1) nnn = nnn*25
       call delps(nnn,rg,1,nk)
      if (.not.prnt) go to 19
      write(*,200) v2,tc,xm,derlp(1),po
 19     if(po.le.1./v2) go to 2
c     pln1 computes refracted part of response
      call pln1(po,to,k,n,tc,kn,v2)
 2     mo=no+2
      if(tg.lt.tn4) mo=2
      if(po.lt.1./v2) mo=2
c      contor computes contour in p-plane
c      ncont=1, first motion approximation only.
       go to (36,13),ncont
 13     if((nfm.gt.0).and.((ttp-to).gt..5)) go to 36
c   if nfm is greater than 0. we use the first motion approx.
c     contor is used to find the complex  p contour
      if(tracing.gt.0)then
          write(*,*)'in high calling contour'
      endif    
      call contor(tmx,m,kn,n,mo,tg)
      if(tracing.gt.0)then
          write(*,*)'in high out calling contour'
      endif    
c     write(*,200) (pp(j),ddpt(j),tt(j),j=jj,m)
      nft=0
      go to 37
c  uses first motion approximation
 36    j=0
      nft=2
      i=mo-1
      pp(i)=po
      tt(i)=to
      sf=sf2(po,k)
34     continue
	i=i+1
      j=j+1
       if(i.gt.300) write (1,100) i
  100   format(1x,'in high,dimension of pp,etc is greater than 300',i10)
      pp(i)=po
      tt(i)=tt(i-1)+(dp/2.)*j*j
      ddpt(i)=1./sqrt(2.*sf*(tt(i)-to)) *(0.,1.)
      test=tt(i)-tstart
      if((test.gt.nn*dp).and.((tt(i)-to).gt.dp)) go to 35
      go to 34
 35     continue
      m=i
 37     continue
      if(.not.prcon) go to 620
      if(k.lt.58) go to 620
      if(k.gt.60) go to 620
       write(*,5)
    5 format (1x,1h0,13x,'pp',27x,'ddpt',24x,'tt')
       jj=mo
  200    format(1x,2e20.8,2e15.4,e20.8)
      mk=m
      do 95 j=mo,m
       pr(j)=pp(j)
        pi(j)=-pp(j)*(0.,1.)
      tq=tt(j)-to
      if(tq.gt.5.) go to 96
 95     continue
 96      m=j
c    first 5 seconds of contour
        ms=mo-1
       pr(ms)=po
        pi(ms)=0.
       mm=m-mo+1
        m=mk
c  pln2- computes reflected part of response
620	continue
      call pln2(po,to,k,mo,m,kn,nft)
      nend=m
      nm=no
      if(po.lt.1./v2) nm=0
      if(tg.lt.tn4)nm=0
c     if(prnts) print99, (td(llm),phi(llm),phh(llm,1),phh(llm,2),
c    2 phh(llm,3),llm=1,nend)
   99 format(1x,1h0,14x,'td',23x,'phi'/(d25.8,4e15.4))
       if(tracing .gt. 0)then
       	  write(*,*)' high called'
       endif	  
      return
      end
      subroutine setup(k,ncount)
      include 'psv.p'
       common/trace/tracing
      common/tstar/tq1,tvl,qa(icc),qb(icc)
        common/str/ph(itt,inz)
       common/pathc/po,to,kk
         common/flt/phh(iflt,inz),ns,beta,nz,nup,nst
      common/thy/tt(itt),pp(itt),ff(itt),l
      common/exact/phi(iexa),tw(iexa),td(iexa),nend,nm
      common/yntcom/ifd,smx
      common/cagcon/tstart,dp,nn,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
      common/plotc/con
      dimension phy(itt)
      logical prnt,prnts,prntc,flat,prcon,prampl
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      if(tracing.gt.0) write(*,*)'entering setup'
c     real *8 tw
c   called by trav
c     ncount-no. of rays with same response
c    k=lmax (lmax is no. of layers + 1)
c  the following write statement was for debugging.
      itest=1
      if(l .gt.1) go to 11
      del=dp
      tt(1)=tstart
      do 10 j=2,nn
      tt(j)=tt(j-1) +del
 10    continue
         do 9 jj=1,nz
        do 12 j=1,nn
 12     ph(j,jj)=0.
 9      continue
 11     continue
       l=l+1
       nf=ncount
c   adjust moves response so that to falls on a tt(j)
c     high generates response
      call high(k)
      call adjust(nn,nfix)
        do 19 j=1,nn
        if(tt(j).gt.td(mkp)) go to 20
  19   continue
  20   nquad=j
c       write(*,468) j,tt(j),nquad,tt(nquad),mkp,td(mkp)
  468  format(1x,'j,tt(j),nquad,tt(nquad),mkp,td(mkp)=',(3(i5,e15.5)))
c   adds response to pp(from here to 99)
        do 32 i=1,nn
        if(td(1).lt.tt(i)) go to 34
   32 continue
        write(*,300)
  300 format (1x,' entire array before tt(1)')
        return
   34 continue
        do 35 ii=1,nn
        if (td(nend).lt.tt(ii)) go to 36
   35 continue
        nnnew=nn
        go to 38
   36 nnnew=ii-1
   38 continue
        do 13 jj=1,nz
        do 15 j=1,nend
 15     phi(j)=phh(j,jj)
c      if(jj.eq.3) write(*,400) jj
  400  format(1x,'in setup, td & phi for component no. ',i5)
c     if(jj.eq.3) write(*,410) (td(j),phi(j),j=1,nend)
 410   format(1x,4(2e14.6,1x))
c     if(jj.eq.3) write(*,410) (ph(j,jj),j=1,nnnew)
      n2=1
      m=nm+1
      if(nfix.lt.1) go to 37
      n2=nfix+1
      n1=nfix-1
      nbtm=n2
      if(n1.le.2)  go to 42
      if (i.gt.n1) go to 42
c      if(jj.eq.1) write(*,437) nf,con
 437  format(1x,' nf,con =',i5,e12.4)
      nlin=1
      do 31 j=i,n1
           if(j.ge.nquad) nlin=2
           y=ynterp(td,phi,tt(j),nend,nlin)
           phy(j)=y*nf*con
           ph(j,jj)=ph(j,jj)+phy(j)
 31     continue
 42    continue
      if (i.gt.nfix) go to 37
        y=phi(m)
        j=nfix
        phy(j)=y*nf*con
        ph(j,jj)=ph(j,jj)+phy(j)
   37    continue
      if (i.gt.n2) nbtm=i
      nlin=1
   	do 33 j=nbtm,nnnew
            if(j.ge.nquad) nlin=2
            y=ynterp(td,phi,tt(j),nend,nlin)
            phy(j)=y*nf*con
            ph(j,jj)=ph(j,jj)+phy(j)
 33     continue
 99     continue
 1000 format(1x,10e11.3)
 13     continue
       if(tracing .gt. 0)write(*,*)'exiting setup'
      return
 200   format(1x,3i10)
  100 format (1x,5e15.6)
   14 format (1x,1h0,15x,'tt',15x,'pp'/(2g18.6))
      end
        subroutine rad(p,rss,rds,r45)
	include 'psv.p'
         common/flt/phh(iflt,inz),ns,beta,nz,nup,nst
       common/soray/nsorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
      common/orst/c(icc),s(icc),d(icc),th(icc),x
        common/rays/na(ina1,ina2),nray(ina1,ina2),nend,lmax,nsp,love,
     *nzr,iwh
        common/testin/iflag1
        common/sourcetype/isource
        complex p,rss,rds,r45
        complex ea,eb,cr
        
        iflag1=0
        if(nsorce.eq.3) vv=s(nb)
        if(nsorce.eq.5) iflag1=1
        if(nsorce.eq.5) vv=c(nb)

        eb=cr(p,vv)
c     nup=0 up going
c     nup=2 down going
        if(nup.eq.0) el=-1.
        if(nup.eq.2) el=1.
c      see  page  eq (7)  of langston  and  helmberger
c        eb=ea for p sorce potential
        if(nsorce.eq.5) go to 20
        if(love.gt.1)  go to 50
c     omega potential
        rss=-1.*el*p*eb
        rds=eb*eb-p*p
        r45=3.*el*p*eb
        go to 30
c       ss same as hark
c       ds opp to hark
c   45 same as hark
 20     continue
c     phi potential
        rss=-p*p
        rds=+2.*p*eb*el
        r45=p*p-2.*eb*eb
         go to 30
 50     continue
c       sh pattern
          vs=vv*vv
        rss=1./vs
         rds=-el*eb/(vs*p)
         r45=0.
 30     continue
c       write(*,nrad)
c    puts in rad. field plus eta factor see (18)
c
c
c
c
        if(isource.eq.0)then
            rss=1./(2.*3.1415926*d(nb)*vv*vv)
            rds=0.0
            r45=0.0
        endif    
        
        rss=rss/eb
        rds=rds/eb
        r45=r45/eb
      return
        end
      subroutine pln2(po,to,k,mo,m,kn,nft)
      include 'psv.p'
       common/trace/tracing
      common/displ/wss(iexa,ior),wds(iexa,ior),w45(iexa,ior),
     *qss(iexa,ior),qds(iexa,ior),q45(iexa,ior),vss(iexa,ior),
     *vds(iexa,ior)
      common/term/iordr
        common/rays/na(ina1,ina2),nray(ina1,ina2),nenn,lmax,nsp,love,
     *nzr,iwh
         common/flt/phh(iflt,inz),ns,beta,nz,nup,nst
      common/magic/pp(imgc),ddpt(imgc),tt(imgc)
      common/exact/phr(iexa),ttt(iexa),tw(iexa),nend,nm
      common/cagcon/tstart,dp,nn,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
      common/coeff/itrnm
      common/orst/c(icc),s(icc),d(icc),th(icc),x
        complex rss,rds,r45
      complex px,pxx,prss,prds,pr45
      common/tmpwr/iwr
        complex rer,ret,rev
        complex t1,t2,t3,t4,t5,t6
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prnts,prntc,flat,prcon,prampl
      complex p,gc,gencof,bt,ddpt,rpp,cr
       complex pp
c     real*8 ttt,tt,po,to
      complex ea,eb,r1,r2,r3,ph
c   if you want fmt 201 written, make itest=1
      itest=0
      if(tracing .gt. 0)write(*,*)'entering pln2'
      do 5 i=mo,m
        tw(i)=tt(i)
        ttt(i)=tt(i)
        p=pp(i)
        gc=1.
       if(nnend.lt.2.or.itrnm.eq.1) go to 17
c      if(i.eq.mo) prnt=.true.
c      if(i.eq.mo) write(*,1055) i,mo,tw(i)
c1055  format(1x,'in pln2, do 10, i,mo,tw(i)=',2i5,e15.5)
c
c
c
c
c      the statement after do 10 will turn on write statements
c      in tranm and refft. (warning statements for temp. change.)
c
c
c
c
       do 10 j=1,nnend
      if(j.eq.1.and.i.eq.mo) iwr=1
 10     gc=gc*gencof(j,p)
      iwr=0
c     prnt=.false.
17      continue
c      if (i.eq.mo) write(*,1000) mo,p
1000  format(1x,' in pln2 mo=',i5,'    p= ',2e13.5)
c       t1=csqrt(p)
	t1=1.
      ph=ddpt(i)*gc*t1
        phr(i)=aimag(ph)
       px=8.*p*x
c     if (i.eq.mo) write(*,1001) p,t1,ph,ddpt(i),gc,phr(i),x,px
1001  format(1x,'p,t1,ph = ',/(6e13.5),/'ddpt(i),gc,phr(i)= ',/(6e13.5),
     +/' x,px = ',/(4e13.5))
      call recev(p,rev,rer,ret)
      call rad(p,rss,rds,r45)
       call dis(p,i,x,ph,rev,rer,ret,rss,rds,r45)
      if(prnt) write(*,111)p,ph,gc
      if(prnt) write(*,112) (vss(i,j),j=1,iordr)
5     continue
110   format(1x,2e15.4)
111   format(1x,6e12.2)
      if(prnt) write(*,111)p,ph,gc
      if(prnt) write(*,111)rev,rer,ret
      if(prnt) write(*,111)rss,rds,r45
      if(prnt) write(*,111)px,pxx
       p=po
      i=mo-1
      pp(i)=po
c   spreading factor, see app. thesis
 100    format(1x,4e15.4)
      sf=sf2(po,k)
       tw(i)=to
       ttt(i)=to
        do 19 j=1,m
 19       tw(j)=ttt(j)
       gc=1.
      if(nnend.lt.2.or.itrnm.eq.1) go to 18
c     prnt=.true.
c     write(*,1057) i,mo,tw(i)
1057  format(1x,'in pln2,do 11, i,mo,tw(i)=',2i5,e15.5)
         do 11 j=1,nnend
 11     gc=gc*gencof(j,p)
c     prnt=.false.
  18     continue
        r3=1.
c       frt=po/(2.*sf)
	frt=1./(2.*sf)
      rpp=gc*sqrt(frt)
      rpp=rpp*r3
 2    pre=real(rpp)
      pim=aimag(rpp)
    3 no=mo-2
      f1=0.
      sum=0.
      if(mo.le.3) go to 46
      tnn=to-dp
      call interp( tw,phr,m,tnn,y)
      f1=y
      sum=sum+2.*pim*(sqrt(abs(to-ttt(no))))
      if (prnt) write(*,4) sum
    4 format (5x,'sum = ',g18.6)
      j=no-1
  92   if(ttt(j).le.tnn) go to 91
      sum=sum+(phr(j+1)+phr(j))*(ttt(j+1)-ttt(j))/2.
      if (prnt) write(*,4) sum
      if(j.le.1) go to 97
      j=j-1
      go to 92
 91    sum=sum+(f1+phr(j+1))*(ttt(j+1)-tnn)/2.
 97   if (prnt) write(*,4) sum
 46     tnf=to+dp
      if(ttt(mo)-to.gt.dp) go to 43
      sum=sum+2.*pre*(sqrt(abs(ttt(mo)-to)))
      if (prnt) write(*,4) sum
      call interp( tw,phr,m,tnf,y)
      f3=y
      j=mo+1
 94   if(ttt(j).ge.tnf) go to 96
      sum=sum+(phr(j)+phr(j-1))*(tt(j)-ttt(j-1))/2.
      if (prnt) write(*,4) sum
      j=j+1
       go to 94
 96     sum=sum+(f3+phr(j-1))*(tnf-ttt(j-1))/2.
      if (prnt) write(*,4) sum
      phr(i)=(3.*sum/dp-f1-f3)/4.
      go to 44
 43    continue

      phr(mo)=pre/(dp**.5)
      f3=phr(mo)
      sum=sum+2.*pre*(dp)**.5
      phr(i)=(3.*sum/dp-f1-f3)/4.
      if (prnt) write(*,4) sum
 44     continue
        ph1=phr(i)
       ph=cmplx(0.,ph1)
c      write(*,100) ph,sum
        p=po
c        write(*,*)'before call recev rad dis'
      call recev(p,rev,rer,ret)
      call rad(p,rss,rds,r45)
       call dis(p,i,x,ph,rev,rer,ret,rss,rds,r45)
c        write(*,*)'after call recev rad dis'
c        write(*,111)p,ph,gc
c        write(*,111)rev,rer,ret
c         write(*,111)rss,rds,r45
 112    format(1x,4e15.4)
c     write(*,112) (wss(i,j),j=1,iordr)
c     write(*,112) (wds(i,j),j=1,iordr)
c     write(*,112) (w45(i,j),j=1,iordr)
c     write(*,112) (qss(i,j),j=1,iordr)
c     write(*,112) (qds(i,j),j=1,iordr)
c     write(*,112) (q45(i,j),j=1,iordr)
c     write(*,112) (vss(i,j),j=1,iordr)
c     write(*,112) (vds(i,j),j=1,iordr)
 300    format(1x,5e14.4)
 200     format(1x,2e15.4)
c    for p   and   sv
c       phh 1   vss
c       phh 2   vds
c       phh 3   qss
c       phh 4    qds
c     phh 5   q45
c      phh 6   wss
c       phh 7   wds
c    phh 8   w45
c       write(*,*)'before call seryy'
        call ser(wss,tw,m)
        call ser(wds,tw,m)
        call ser(w45,tw,m)
        call ser(qss,tw,m)
        call ser(qds,tw,m)
        call ser(q45,tw,m)
        call ser(vss,tw,m)
        call ser(vds,tw,m)
c       write(*,*)'after call seryy'
      if(itest.eq.1) write(*,201) m
  201   format(1x,'in pln2, m=',i10)
        do 125 j=1,m
       phh(j,1)=vss(j,1)
       phh(j,2)=vds(j,1)
        phh(j,3)=qss(j,1)
        phh(j,4)=qds(j,1)
       phh(j,5)=q45(j,1)
        phh(j,6)=wss(j,1)
        phh(j,7)=wds(j,1)
        phh(j,8)=w45(j,1)
  125     continue
      j=mo-1
c      write(*,113) tw(j),phh(j,1)
c      if(prampl) write(*,113)(tw(j),phh(j,1),phh(j,2),phh(j,3),
c      write(*,113)tw(j),phh(j,1),phh(j,2),phh(j,3),
c     +phh(j,4),phh(j,5),phh(j,6),phh(j,7),phh(j,8)
c     write(*,113)(tw(j),phh(j,5),phh(j,6),phh(j,7),phh(j,8),j=1,m)
c113   format(1x,e18.6,4e15.4)
  113  format(1x,e15.6,8e13.4)
      if(tracing .gt. 0)write(*,*)'exiting pln2'
      return
      end
      subroutine ser (a,b,m)
      include 'psv.p'
      dimension a(iexa,ior),b(iexa),t(ior)
      common/term/iordr
c           integration of series
c      write(*,*)'iordr=',iordr
      do 3 i=1,iordr
      t(i) = 0.
    3 continue
       if (iordr.eq.1) return
      do 1 j=2,m
      det=(b(j)-b(j-1))*.5
      k = iordr
      l = iordr-1
      do 4 i=1,l
      t(i) = t(i)+det*(a(j,k)+a(j-1,k))
      k=k-1
      a(j,k) = a(j,k)+t(i)
    4 continue
    1 continue
      return
      end
      subroutine dis(p,i,x,ph,rev,rer,ret,rss,rds,r45)
      include 'psv.p'
         common/flt/phh(iflt,inz),ns,beta,nz,nup,nst
      common/displ/wss(iexa,ior),wds(iexa,ior),w45(iexa,ior),
     *qss(iexa,ior),qds(iexa,ior),q45(iexa,ior),vss(iexa,ior),
     *vds(iexa,ior)
      common/term/iordr
       complex aknu,akder
      dimension aknu(ior,3),akder(ior,3)
      complex p,ph,rev,rer,ret,rss,rds,r45
      complex px,t1,t2,t3,t4,t5,t6
        complex prss,prds,pr45
      px = p*x
       jj=i
      do 87 i=1,iordr
       vds(1,i)=0.
       vss(1,i)=0.
       q45(1,i)=0.
       qds(1,i)=0.
       qss(1,i)=0.
       w45(1,i)=0.
       wds(1,i)=0.
       wss(1,i)=0.
 87     continue
        i=jj
      prss=ph*rss
      prds=ph*rds
      pr45=ph*r45
c        write(*,*)'ph,rss,rds,r45',ph,rss,rds,r45
      call kfunc(aknu,akder,px,iordr)
      if(nst.eq.4) go to 80
c  assuming p or sv ray
        t1=prss*rev
       t2=prds*rev
        t3=pr45*rev
        t4=prss*rer
        t5=prds*rer
        t6=pr45*rer
      do 83 j=1,iordr
      wss(i,j) = aimag(t1*aknu(j,3))
      wds(i,j) = aimag(t2*aknu(j,2))
      w45(i,j) = aimag(t3*aknu(j,1))
      qss(i,j) = aimag(t4*akder(j,3))
      qds(i,j) = aimag(t5*akder(j,2))
      q45(i,j) = aimag(t6*aknu(j,2))
   83 continue
      t1 = ret*prss/(x*.5)
      t2 = ret*prds/x
      vss(i,1) = 0.
      vds(i,1) = 0.
       if (iordr.eq.1) go to 81
      do 84 j=2,iordr
      k = j-1
      vss(i,j) = aimag(t1*aknu(k,3))
      vds(i,j) = aimag(t2*aknu(k,2))
   84 continue
      go to 81
   80 continue
c  sh  potential
      do 82 j=1,iordr
      wss(i,j) = 0.
      wds(i,j) = 0.
      w45(i,j) = 0.
   82 q45(i,j) = 0.
      t1 = prss*ret
      t2 = prds*ret
      t4 = prss*rer*(-2./x)
      t5 = prds*rer*(-1./x)
      qss(i,1) = 0.
      qds(i,1) = 0.
      do 85 j=1,iordr
      vss(i,j) = aimag(t1*akder(j,3))
      vds(i,j) = aimag(t2*akder(j,2))
      if (j.eq.1) go to 85
      k = j-1
      qss(i,j) = aimag(t4*aknu(k,3))
      qds(i,j) = aimag(t5*aknu(k,2))
   85 continue
   81 continue
      return
      end
      subroutine kfunc(t,tder,px,iordr)
      include 'psv.p'
      dimension t(ior,3),tder(ior,3)
      complex t,tder,px
      do 5 i=1,3
      t(1,i)=1.0
    5 continue
       if (iordr.eq.1) go to 6
      do 1 i=1,3
      mu = 4*((i-1)**2)
      do 2 j=2,iordr
      k=j-1
      t(j,i) = t(k,i)*(mu-(2*k-1)**2)/(px*k*8.)
    2 continue
    1 continue
    6  continue
      do 3 i=1,3
      tder(1,i)=1.0
    3 continue
       if (iordr.eq.1) return
      do 4 j=2,iordr
      k=j-1
      tder(j,2) = t(j,1)+t(k,2)/px
      tder(j,3) = t(j,2)+2*t(k,3)/px
    4 continue
      return
      end
      subroutine pln1(po,to,k,n,tc,kn,v2)
      include 'psv.p'
       common/trace/tracing
      common/spe/derlp(ispe),dd1,dd2,dd3,dd4,no
      common/orst/c(icc),s(icc),d(icc),th(icc),x
      common/cagcon/tstart,dp,nn,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
      common/coeff/itrnm
      common/exact/phi(iexa),tt(iexa),tw(iexa),nend,nm
        complex rss,rds,r45,rvz,rvr
         common/flt/phh(iflt,inz),ns,beta,nz,nup,nst
        common/rays/na(ina1,ina2),nray(ina1,ina2),nenn,lmax,nsp,love,
     *nzr,iwh
        complex  prss,prds,pr45,rev,rer,ret,a,ph
        complex   t1,t2,t3,t4,t5,t6
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prnts,prntc,flat,prcon,prampl
      complex rpr,z,tq,gencof
c     real*8 po,to,tc,v2,y,ttp,pp,tt,ww
c
c     pln1 is called by subroutine high
c
 100   format(1x,6e15.5)
c    v2  is max. p vel. in layers used in ray spec
       if(tracing.gt.0)then
           write(*,*)'pln1 is being called'
       endif    
       y=1./v2
       r3=1.
c   no is defined in subroutine high
      do 80 i=2,no
      j=i-1
       y=y+derlp(j)
         u=y
       z=cmplx(u,0.)
      p=y
      call herlp(k,y,ttp,dtp)
      tt(i)=ttp
        tq=1.
       if(nnend.lt.2.or.itrnm.eq.1) go to 18
c     if(i.eq.no) prnt=.true.
c      if(i.eq.no) write(*,1055) i,no,tt(i)
1055  format(1x,'in pln1, do 7, i,no,tt(i)=',2i5,e15.5)
       do 7 j=1,nnend
 7      tq=tq*gencof(j,z)
c     prnt=.false.
 18     continue
c       ph=tq*dtp*p**.5
       ph=tq*dtp
 2      phi(i)=aimag(ph)
      call rad(z,rss,rds,r45)
      call recev(z,rev,rer,ret)
       call dis(z,i,x,ph,rev,rer,ret,rss,rds,r45)
 80     continue
       i=no
  4   if(to-ttp.lt.dltp) go to 3
      ww=po-y
       y=y+ww/2.
         u=y
       z=cmplx(u,0.)
      p=y
      i=i+1
      no = i
      call herlp(k,y,ttp,dtp)
      tt(i)=ttp
         tq=1.
       if(nnend.lt. 2.or.itrnm.eq.1) go to 19
        do 12 j=1,nnend
 12     tq=tq*gencof(j,z)
 19     continue
c       ph=tq*dtp*p**.5
       ph=tq*dtp
        phi(i)=aimag(ph)
       px=8.*p*x
      call recev(z,rev,rer,ret)
      call rad(z,rss,rds,r45)
       call dis(z,i,x,ph,rev,rer,ret,rss,rds,r45)
      go to 4
 3    tt(1)=tc
      phi(1)=0.
*      write(*,200) tc
 200   format(1x,' tc = ',e14.4)
      return
       end
       subroutine recev(p,rev,rer,ret)
       include 'psv.p'
        complex rer,ret
       common/orst/cc(icc),ss(icc),dd(icc),tth(icc),xx
        common/rays/na(ina1,ina2),nray(ina1,ina2),nend,lmax,nsp,love,
     *nzr,iwh
       common/soray/nsorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
      common/cagcon/tstart,dp,nn,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
       common/testin/iflag1
       common/down/ndown, nstress
        common/sourcetype/isource
        complex rst,rpt
        complex rpz,rpr,rsz,rsr
       complex ea,eb,r1,r2,r3,r4,r5,rev,p,q,cr
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prntc,prnts,flat,prcon,prampl
       complex lamda, mu
c    two factor taken out
       iflag1=0
       if(love.gt. 1)  go to 13
        nde=nray(nend,nsp)
c      nzr=0 vertical comp.
c        nzr=2     radial comp.
        q=p
	mu   =tqd*tqs*tqs
	lamda=tqd*tqc*tqc-2.0*mu
c        write(*,*)"mu lamda ", mu, " ",lamda
        ea=cr(p,tqc)
        eb=cr(p,tqs)
        r1=eb**2-p*p
        r2=r1**2+4.*p*p*ea*eb
        r2=r2*(tqs**2)
        rpz=(r1*ea)/r2
        rsz=(2.*p*ea*eb)/r2
        rsr=(eb*r1)/r2
        rpr=rsz
c       z  positive downward
c        see a3,a4,a7,a9, helm    1974
         rpr=-rsz
c         see a9
        rst=(-1.*eb*r1)/(p* r2)
c      near field a6  rpt
        rpt=-1.*rpr/p
        if (nde.eq.5)  go to 14
c   iwh=0  whole space,   iwh=2  half space
c     sv  wave
      if (iwh.eq.0) go to 2
        rev=rsz
        rer=rsr
        ret=rst
c     write(*,100) rev,rer,ret,p,ea,eb,r2,r1
  100   format(1x,'in recev, rev,rer,ret,p='/8e11.3/
     *' ea,eb,r2,r1='/8e11.3)
      go to 3
    2 continue
       if(nstress.le.1)then
           rev=p/2.
           rer=eb/2.
           if(ndown.eq.-1)rer=-rer
           ret=-eb/p/2.
       endif
       if(nstress.eq.22)then
	   rev=0.0
	   rer=mu*p*eb
	   ret=1.0
	   if(ndown.eq.-1)rer=-rer
       endif
       if(nstress.eq.11)then
	   rev=0.0
	   rer=-mu*p*eb
	   ret=1.0
	   if(ndown.eq.-1)rer=-rer
       endif
       if(nstress.eq.12)then
	   rev=0.5*mu*(eb*eb-p*p)
	   rer=0.0
	   ret=1.0
       endif
       if(isource.eq.0)then
           rer=1.0
           ret=0.0
           rev=eb    
           if(ndown.eq.-1)rev=-rev
       endif        
    3 continue
        return
 14      continue
c        p   wave
       if (iwh.eq.0) go to 4
        rev=rpz
        rer=rpr
        ret=rpt
        go to 5
    4 continue
	if(nstress.le.1)then
            rev=ea/2.
	    if(ndown.eq.-1)rev=-rev
            rer=-p/2.
            ret=1./2.
        endif
	if(nstress.eq.22)then
	    rev=0.0
	    rer=0.5*(lamda*p*p+(lamda+2.*mu)*ea*ea)
	    ret=1.0
        endif
	if(nstress.eq.11)then
	    rev=0.0
	    rer=0.5*((lamda+2.*mu)*p*p+lamda*ea*ea)
	    ret=1.0
        endif
	if(nstress.eq.12)then
	    rev=-mu*p*ea
	    rer=0.0
	    ret=1.0
	    if(ndown.eq.-1)rev=-rev
        endif
	if(nstress.eq.31)then
c           /* grad (div) x and Dgrad(div)x/Dz  */
	    rer = -(0.5/tqc/tqc)*p
            rev = rer * ea 
	    if(ndown.eq.-1)rev=-rev
	    ret=1.0
        endif
	if(nstress.eq.32)then
c           /* grad (div) z and Dgrad(div)z/Dz  */
	    rer = (0.5/tqc/tqc)*ea
	    if(ndown.eq.-1)rer=-rer
	    rev = (0.5/tqc/tqc)*ea*ea
	    ret=1.0
        endif
       if(isource.eq.0)then
           rer=1.0
           ret=0.0
           rev=ea    
       endif        
    5 continue
        return
 13      continue
c     sh   wave
        rev=0.
        rer=1.
        ret=p
        return
       end
      subroutine trav(ncount)
      include 'psv.p'
      common/cagcon/tstart,dp,nn,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
      common/orst/c(icc),s(icc),d(icc),th(icc),x

        common/rays/na(ina1,ina2),nray(ina1,ina2),nend,lmax,nsp,love,
     +nzr,iwh
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
      common/tstar/tq,tvl,qa(icc),qb(icc)
c     ncount= number of rays with the same responce
c     nend+1= number of legs in either mode
c     nray code: p=5, s=3
c     ll is the number of layers involved, including the water layer
c     nend is the number of layers traversed -1, the last leg to the rec
c     assumed to be p
c     testing na and nray for lp and ls information
  100 format(1x,1x,e18.6)
      nnend=nend
         do 1 j=1,ll
      rlp(j)=0.
    1 rls(j)=0.
         do 5 i=1,nend

      nry=nray(i,nsp)
       ny=na(i,nsp)
      if(nry    .eq.5) rlp(ny   )=rlp(ny   )+1.
      if(nry.eq.3.or.nry.eq.4) rls(ny   )=rls(ny   )+1.
    5 continue
c   lmax computed in main, is deepest layer penetration of ray.
      k=lmax
 200    format(1x,6i6)
c       write(*,113) nsp
c     write(*,113) nsp
 113    format(1x,/' ray no =',i10)
 8      call setup(k,ncount)
 6      return
      end
      function gencof(i,q)
      include 'psv.p'
       common/soray/nsorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
        common/rays/na(ina1,ina2),nray(ina1,ina2),nend,lmax,nsp,love,
     +nzr,iwh
      common/cagcon/tstart,dp,nz,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
      common/orst/c(icc),s(icc),d(icc),th(icc),x
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prnts,prntc,flat,prcon
      complex q,tpp,tps,tsp,tss,rpp,rps,rsp,rss,gencof
c     nslb  source is situated at the bottom of layer nslb
      gencof=1.
c     i is the index before the interaction
c     (i+1) is the index after interaction
c     stuff model parameters
        nslb=nb
       m=na(i+1,nsp)
       k=na(i,nsp)
      mr=nray(i+1,nsp)
      kr=nray(i,nsp)
       if(m.eq.k) go to 14
      if(ndirt.eq.0) go to 33
      if(m.eq.1) go to 34
      if(k.eq.1) go to 34
 33     continue
      call tranm(q,c(k),s(k),d(k),c(m),s(m),d(m),tpp,tps,tsp,tss)
      if(m.gt.k) go to 11
      tps=-tps
      tsp=-tsp
   11 continue
      if(kr.eq.mr) go to 12
      if(kr     .eq.5) gencof=tps
      if(kr.eq.3.or.kr.eq.4) gencof=tsp
      go to 20
 12    if(kr.eq.3.or.kr.eq.4) gencof=tss
      if(kr     .eq.5) gencof=tpp
      go to 20
 14     m=k-1
      if((i.eq.1).and.(k.eq.nslb)) go to 200
       m=k+1
      if((i.eq.1).and.(nslb.lt.k)) go to 200
      if(i.ne.1) go to 300
      m=2
      go to 200
  300 ij=i
   15 ij=ij-1
      if(ij.eq.1) go to 13
       if(na(ij,nsp).eq.na(i,nsp)) go to 15
   13 nn=i-ij
      if(na(ij,nsp).le.k) nfix=1
      if(na(ij,nsp).gt.k) nfix=-1
      r=(-1.)**nn*nfix
      if(r.gt.0.) m=k-1
      if(r.le.0.) m=k+1
 200    continue
      if(m.eq.1) go to 61
      call refft(q,c(k),s(k),d(k),c(m),s(m),d(m),rpp,rps,rsp,rss)
      go to 63
 61    call refft(q,c(2),s(2),d(2),c(m),s(m),d(m),rpp,rps,rsp,rss)
 63     continue
c
c     check polarity for z downward
c
c	if(m.eq.1)write(*,*)'rpp,rps,rsp,rss',rpp,rps,rsp,rss
       if (m.gt.k) go to 21
        rps=-rps
      rsp=-rsp
   21 continue
      if(kr.eq.mr) go to 22
      if(kr     .eq.5) gencof=rps
      if(kr.eq.3.or.kr.eq.4) gencof=rsp
      go to 20
   22 if(kr.eq.3.or.kr.eq.4) gencof=rss
      if(kr     .eq.5) gencof=rpp
 34     continue
 20      continue
        qr=q
         qi=aimag(q)
        qt=qi*10000.
c       if(qr.gt.qt) go to 80
c      if(prnt) go to 80
c       if(ndirt.eq.2) return
c       if(m.ne.1.and.k.ne.1) return
      if(.not.prnt) return
 80     continue
      write(*,110) tpp,tps,tsp,tss
       write(*,111)rpp,rps,rsp,rss
      write(*,112) c(k),s(k),d(k),c(m),s(m),d(m)
 110  format(15x,'tpp',15x,'tps',15x,'tsp',15x,'tss'/(8g10.3))
 111  format(15x,'rpp',15x,'rps',15x,'rsp',15x,'rss'/(8g10.3))
112   format(18x,'stuff(k) and (m)'/(6g15.4))
        return
      end
      subroutine tranm(p,v1,s1,d1,v2,s2,d2,tpp,tps,tsp,tss)
      include 'psv.p'
        common/rays/na(ina1,ina2),nray(ina1,ina2),nend,lmax,nsp,love,
     +nzr,iwh
       common/testin/iflag1
      common/tmpwr/iwr
           complex t,ts,tpp,tsp,tps,tss,bs
           complex p,cr,e1p,e2p
           complex a,b,ap,bp
           complex e1,e2,c1,c2,c3,c4,c5,c6
      real k1,k2,k3,k4,d1,d2
       iflag1=0
	tiny = 1.e-6
c	if( abs(s2-s1) .lt. tiny .and. abs(d2-d1) .lt. tiny)then
c	    if(love .gt. 0)then
c		tss = cmplx(1.0,0.0)
c		return
c             else if( abs(v2-v1) .lt. tiny) then
c		tss = cmplx(1.0,0.0)
c		tpp = cmplx(1.0,0.0)
c		tps = cmplx(0.0,0.0)
c		tsp = cmplx(0.0,0.0)
c		return
c	     endif
c        endif
      if(abs(s2-s1) .lt. tiny) s2 = s1 + tiny
      if(abs(d2-d1) .lt. tiny) d2 = d1 + tiny
      if(abs(v2-v1) .lt. tiny) v2 = v1 + tiny
      rho1=d1
      rho2=d2
ctemp write
c     write(*,601) s1,s2,v1,v2
c     if(iwr.eq.1) write(*,601) s1,s2,v1,v2
 601     format(1x,'s1,s2,v1,v2=',4e12.4)
      k4=rho2*s2**2/(rho1*s1**2)
      b1=.5/(1.-k4)
      b2=.5*k4/(k4-1.)
      k1=b1/s1**2
      k2=b2/s2**2
      k3=k1+k2
ctemp write
c     write(*,602) k1,k2,k3,k4,b1,b2
 602   format    (1x,'k1,k2,k3,k4,b1,b2='/6e12.4)
c     write(*,603) p
 603    format(1x,'before cr calls, p=',2e13.4)
c     if(iwr.eq.1) write(*,602) k1,k2,k3,k4,b1,b2
c     if(iwr.eq.1) write(*,603) p
      e1=cr(p,v1)
      e2=cr(p,v2)
      e2p=cr(p,s2)
      e1p=cr(p,s1)
ctemp write
c     write(*,604) e1,e2,e1p,e2p
c     write(*,605) love
 605    format(1x,'love=',i5)
 604     format(1x,'e1,e2=',4e13.4,'e1p,e2=',4e13.4)
c     if(iwr.eq.1) write(*,604) e1,e2,e1p,e2p
c     if(iwr.eq.1) write(*,605) love
         if(love.gt.0) go to 10
      c1=(p**2)*(k3-p**2)**2
      c2=p**2*e1*e1p*e2p
c     write(*,607) c1,c2
c     if(iwr.eq.1) write(*,607) c1,c2
 607    format(1x,'c1,c2=',/4e13.4)
      c3=(e1*e1p)*(k2-p**2)**2
      c4=e2p*(k1-p*p)**2
c     write(*,608) c3,c4
c     if(iwr.eq.1) write(*,608) c3,c4
 608    format(1x,'c3,c4=',/4e13.4)
      c5=k1*k2*e1*e2p
      c6=k1*k2*e1p
c     write(*,606) c5,c6
 606   format(1x,'c5,c6=',/4e13.4)
c     if(iwr.eq.1) write(*,606) c5,c6
      ap=c1+c3-c5
      bp=c2+c4-c6
      bs=ap+e2*bp
c     write(*,609) ap,bp,bs
 609    format(1x,'ap,bp,bs=',/ 6e13.4)
c     if(iwr.eq.1) write(*,609) ap,bp,bs
      t=2.*k1*e1*(e2p*(k1-p**2)-e1p*(k2-p**2))
      tpp=t/bs
c     write(*,610) t,tpp
 610   format(1x,'t,tpp=',/4e13.4)
c     if (iwr.eq.1) write(*,610) t,tpp
      t=2.*k1*p*e1*(e1p*e2-(k3-p**2))
      tps=t/bs
c     write(*,611) t,tps
c     if(iwr.eq.1) write(*,611) t,tps
 611   format(1x,'t,tps=',/4e13.4)
      ts=2.*k1*p*e1p*((k3-p**2)-e1*e2p)
      tsp=ts/bs
c     write(*,612) ts,tsp
c     if(iwr.eq.1) write(*,612) ts,tsp
 612    format(1x,'ts,tsp=',/4e13.4)
      ts=-2.*k1*e1p*(e1*(k2-p**2)-e2*(k1-p*p))
      tss=ts/bs
c     write(*,613)  ts,tss
 613    format(1x,'ts,tss=',/4e13.4)
c     write(*,614)
c     if(iwr.eq.1) write(*,613)  ts,tss
c     if(iwr.eq.1) write(*,614)
 614    format(1x,'return to gencof')
       return
 10     r1=d1*s1**2
        r2=d2*s2**2
        tss=(2.*r1*e1p)/(r1*e1p+r2*e2p)
c     write(*,615) tss
c     if(iwr.eq.1) write(*,615) tss
 615   format(1x,' in tranm, tss=',2e13.4)
       return
        end
       subroutine refft(p,v1,s1,d1,v2,s2,d2,rpp,rps,rsp,rss)
       include 'psv.p'
        common/rays/na(ina1,ina2),nray(ina1,ina2),nend,lmax,nsp,love,
     +nzr,iwh
           complex p,e1,e2,e1p,e2p
           complex c1,c2,c3,c4,c5,c6,a,b,ap,bp,cr,rft
           complex rpp,rps,rss,rsp  ,bt,aps,bps,asp,bsp
           real k1,k2,k3,k4 ,d
        common/testin/iflag1
        iflag1=0
	tiny = 1.e-6
c	if( abs(s2-s1) .lt. tiny .and. abs(d2-d1) .lt. tiny)then
c	    if(love .gt. 0)then
c		rss = cmplx(0.0,0.0)
c		return
c             else if( abs(v2-v1) .lt. tiny) then
c		rss = cmplx(0.0,0.0)
c		rpp = cmplx(0.0,0.0)
c		rps = cmplx(0.0,0.0)
c		rsp = cmplx(0.0,0.0)
c		return
c	     endif
c        endif
      if(abs(s2-s1) .lt. tiny) s2 = s1 + tiny
      if(abs(d2-d1) .lt. tiny) d2 = d1 + tiny
      if(abs(v2-v1) .lt. tiny) v2 = v1 + tiny
 
      k4=s2**2*d2/(s1**2*d1)
      b1=.5/(1.0-k4)
      b2=.5*k4/(k4-1.0)
      k1=b1/s1**2
      k2=b2/s2**2
      k3=k1+k2
      e1=cr(p,v1)
      e2=cr(p,v2)
      e1p=cr(p,s1)
      e2p=cr(p,s2)
         if(love.gt.0) go to 10
      c1=(p**2)*(k3-p**2)**2
      c2=p**2*e1*e1p*e2p
      c3=(e1*e1p)*(k2-p**2)**2
      c4=e2p*(k1-p*p)**2
      c5=k1*k2*e1*e2p
      c6=k1*k2*e1p
      ap=c1+c3-c5
      bp=c2+c4-c6
      a=-c1+c3-c5
      b=-c2+c4-c6
      bt=ap+bp*e2
      rpp=(a-b*e2)/bt
      aps=2.*p*e1 *(k2-p*p)*(k3-p*p)
      bps=2.*p*e1*(k1-p*p)*e2p
      rps=(aps-bps*e2)/bt
      a=-c1 +c3 +c5
      b=-c2 +c4 +c6
       rss=(a-b*e2)/bt
      asp=2.*p*e1p*(k2-p*p)*(k3-p*p)
      bsp=2.*p*e1p*(k1-p*p)*e2p
      rsp=-(asp-bsp*e2)/bt
        return
 10     r1=d1*s1**2
        r2=d2*s2**2
        rss=(r1*e1p-r2*e2p)/(r1*e1p+r2*e2p)
       return
      end
      subroutine curay(jo)
      include 'psv.p'
      dimension drcsq(icc),drssq(icc)
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
       common/orst/cc(icc),ss(icc),dd(icc),tth(icc),xx
        common/fixp/ddn(icc),arn(icc)
      logical prnt,prnts,prntc,flat,prcon
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      common/depths/depth(icc)
      common/tstar/tq,tvl,qa(icc),qb(icc)
       dimension dt(icc)
c     real*8 ddn,arn,rcsq,rssq
c     real*8 c,s,d,th,x
c     real*8 drcsq,drssq
      x=xx
7	continue
  8   depth(1)=tth(1)/2.0
      do 10 j = 2,jo
 10    depth(j)=depth(j-1)+(tth(j)+tth(j-1))/2.
      do 5 j = 1,jo
      q = 6371.0 / (6371.0-depth(j))
      if(flat) q=1.
      arn(j)=1./q
      c(j)=(cc(j))*q
      s(j)=(ss(j))*q
      d(j)=(dd(j))*q
      th(j)=(tth(j))*q
      drcsq(j)=1./c(j)**2
      drssq(j)=1./s(j)**2
      rcsq(j) = drcsq(j)
      rssq(j) = drssq(j)
 5      continue
      dt(1)=tth(1)
      do 20 j=2,jo
      dt(j)=dt(j-1) +tth(j)
 20     continue
      do 25 j=1,jo
      if(flat) dt(j)=0.0
      ddn(j)=6371.-dt(j)/6371.
 25     continue
        do 88 j=1,jo
         cc(j)=(c(j))
         ss(j)=(s(j))
         dd(j)=(d(j))
        tth(j)=(th(j))
        ddn(j)=1.
        arn(j)=1.
 88     continue
  9   return
      end
c /* This subroutine is changed by Lianxing Wen on Nov. 20, 1995 */
c /* There are several problems for this subroutine: */
c /*     1. sometimes, it is very difficult to find the point p  */
c /*                   on the contour, which makes t real        */
c /*     2. when you reduce det criterion, sometimes, it finds   */
c /*                   the wrong t, which could be smaller than  */
c /*                   the previous one                          */ 
      subroutine time2(pr,pi,dl,q,dpt,t,kn,n,pil,det,njtm)
      include 'psv.p'
      common/pathc/po,to,k
      common/trace/tracing
      common/soray/nsorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
      dimension e(icc),f(icc)
      complex pc,dpt,tq,cr
       complex p,e,bl,t,q,psq,f
c       real*8 c,s,d,th,r
       common/testin/iflag1
c       real*8 rcsq,rssq
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),r
      common / lprint/ prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prnts,prntc,flat,prcon
      real*8 pii, prr, ej, fj, rex 
      real*8 ti,del, x1m, x4m

      if(tracing.gt.0)then
          write(*,*)'entering time2'
      endif    

      max_num = 100

      det0 = det
      ko=k
3     continue
      pii = pi
      del = dl
     
c      first ti
 	p=pr*(1.,0.)+pii*(0.,1.)

        t=p*r
       psq=p*p
       do j=1,ko
           e(j)=cr(p,c(j))
           f(j)=cr(p,s(j))
           t=t+rlp(j)*th(j)*e(j)+rls(j)*th(j)*f(j)
       enddo    

       ti = -t*(0.,1.)
       if(ti.le.0) del = -del

       i = 1
 1     continue

           tl = ti

           x1m = pii
	   pii = pii + del
	   x4m = pii
 	   p   = pr*(1.,0.)+pii*(0.,1.)
	   i  = i +1

           t=p*r
           do j=1,ko
              e(j)=cr(p,c(j))
              f(j)=cr(p,s(j))
              t=t+rlp(j)*th(j)*e(j)+rls(j)*th(j)*f(j)
           enddo    
           ti = -t*(0.,1.)
           if(ti.lt.0.) del =  pii*(i**2)
           if(ti.gt.0.) del = -pii/2.

           tii = ti
           if(abs(tii).le.det0) go to 2

          if(ti*tl.le.0.0) go to 44
	  if(i.gt.max_num) then
	      write(*,*)"Fail to find the cross point"
	      goto 40
          endif

      go to 1
      
 44    continue 
       tl = ti

 45    continue 

           i = i+1
           pii=(x1m+x4m)*.5
	   prr = pr
 	   p   = pr*(1.,0.)+pii*(0.,1.)

           ti = pii*r
           do j=1,ko
              call doublecr(prr, pii, c(j), rex, ej)
              call doublecr(prr, pii, s(j), rex, fj)
              ti = ti + rlp(j)*th(j)*ej+rls(j)*th(j)*fj
           enddo    

           if(ti*tl.le.0.)then
c               /* crossing the point */
                x1m = pii
           else
		x4m = pii
           endif

           tii = ti
           if(abs(tii).le.det0) go to 2

	   if(i.gt.max_num)then
	       write(*,*)"Fail to locate the ti=",ti, " for det0=", det0
	       goto 4
           endif

      go to 45 

40       continue
c           /* fail to find pi reduce  and find again */
	    write(*,*)"Search Fails, Reduce pi and search again"
            pi = 0.5 * pi 
	    goto 3

4       continue
c           /* fail to find pi reduce det0 and find again */
	    write(*,*)"Search Fails, Reduce det0 and search again"
            det0 = 2.0 *det0
	    goto 3


      

2       continue
       bl=r
       t = p*r
      do j=1,ko
          e(j)=cr(p,c(j))
          f(j)=cr(p,s(j))
          bl=bl-p*th(j)      *(rlp(j)/e(j)+rls(j)/f(j))
          t=t+rlp(j)*th(j)*e(j)+rls(j)*th(j)*f(j)
      enddo    

c       write(*,*)"Successfully found p=", p, " for t=",t
      q=p
      dpt=1./bl
      k2=k
      det = det0

      if(i.gt.max_num) go to 17
      
      tdpt=real(dpt)
      if(tdpt.le.0)then
           pi=pi*5.
c        I don't know why
           write(*,*) "Found tdpt <= 0 go back again "
c  comment out for certain purpose (it may need to adjust back later)
c           go to 3
      endif     
      if(tracing.gt.0)then
          write(*,*)'exiting time2'
      endif    
      return
 17     write(*,111)
 111   format(1x,' time2 failed'/)
      write(*,110) p,e(k2),t,det,r,pr
  110 format(4x,'p=',2e18.6,20x,'e(1)=',2e17.6/5x,'t=',2e18.6,
     +20x,'det =',e18.6,'r=',e18.6,'pr=',e18.6)
        return
       end
      subroutine adjust(nn,nfix)
      include 'psv.p'
      common/thy/t(itt),ph(itt),tt(itt),lz
      common/exact/phi(iexa),td(iexa),tw(iexa),nend,nm
      m=nm+1
      tr=td(m)
      i=0
 80    i=i+1
      if(i.gt.nn) go to 70
      if(t(i).gt.tr) go to 81
      go to 80
 81    dne=tr-t(i-1)
      dpl=t(i)-tr
      if(abs(dne).gt.abs(dpl)) go to 83
      delta=-dne
      nfix=i-1
      go to 85
 83     delta=dpl
      nfix=i
 85    do 84 j=1,nend
      td(j)=td(j)+delta
      tw(j)=td(j)
 84     continue
      return
 70     nfix=0
      do 89 j=1,nend
 89    tw(j)=td(j)
        return
      end
         function ts(kst,kend)
	 include 'psv.p'
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
        common/fixp/ddn(icc),arn(icc)
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
      logical prnt,prnts,prntc,flat,prcon
c     real*8 arn,ddn,x1,po,to,rcsq,rssq
c     real*8 c,s,d,th,x
      dimension t(icc)
      det=1.e+12
      do 30 j=1,kst
        rls(j)=0
 30     rlp(j)=2
      n=0
      x1=0.
      do 102 j=kst,kend
      rlp(j)=2
      rls(j)=0
      n=n+1
      px=1./c(j+1)
      t(n)=1.e+7
c      if(px.gt.po) go to 97
      tx=ptim(px,j)
      write(*,100) j,tx
      t(n)=tx
 97     continue
 102    continue
      ts=1.e+6
       do 106 j=1,n
       ts=amin1(t(j),ts)
      if(t(j).le.ts) mray=j
 106   continue
      if (prnt) write(*,1) (t(j),j=1,n)
    1 format (5x,'t(j),j=1,n'/(4g18.4))
 100   format(1x,i10,2e18.6)
      return
      end
      subroutine herlp(k,p,ttp,dtp)
      include 'psv.p'
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
c     real*8 c,s,d,th,x
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
c     real*8 p,ttp,psq,totem,bltem,bl,rcsq,rssq,e,f
c     call stout('entering help|')
      psq = p**2
      bltem = 0.0
c     write(*,101) p,k
 101    format(1x,' p=',e12.3,'k=',i10)
      totem = 0.0
      do 11 j = 1,k
c      e=sqrt(abs(rcsq(j)-psq))
      e=abs(rcsq(j)-psq)
        if(e.lt.1.e-7) e=1.e-7
c     f=sqrt(abs(rssq(j)-psq))
      f=abs(rssq(j)-psq)
      if (f.lt.1.e-7) f=1.e-7
c     write(*,103) j,f,e
 103     format(1x,'j=',i10,'   f=',e12.3,'   e=',e12.3)
c     totem=totem+e*th(j)*rlp(j)+f*th(j)*rls(j)
      rlpsq=rlp(j)*rlp(j)
      thsq=th(j)*th(j)
      rlssq=rls(j)*rls(j)
      totem=totem+sqrt(e*thsq*rlpsq)+sqrt(f*thsq*rlssq)
      bltem=bltem-sqrt(rlpsq*thsq/e)-sqrt(rlssq*thsq/f)
c      bltem=bltem-rlp(j)*th(j)/e-rls(j)*th(j)/f
  11   continue
      bl = x + p*bltem
      ttp=p*x+totem
c     write(*,100) x,p,bltem,bl
 100    format(1x,'in help,x=',e12.3,'p=',e12.3,'bltem=',e12.3,
     *'bl=',e12.3)
      dtp=1./bl
      return
      end
      subroutine delps (nnn,rg,nn,n)
      include 'psv.p'
      dimension pp(ispe)
      common/spe/derlp(ispe),dd1,dd2,dd3,dd4,no
c      real*8 rg
 100   format(1x,'rg in delps = ',e11.4)
c    set itest = 1  if you want n,nn,nnn,and rg info from delps
      itest=0
      if(itest.eq.1) write(*,100) rg
      rg=rg-1.e-08
      pi=3.141593
      if(itest.eq.1) write(*,110) nnn
 110    format(1x,'in delps, nnn=',i10)
      an=pi/(nnn*2.)
      if(itest.eq.1) write(*,111) nn
111    format(1x,'in delps, nn=',i10)
      j=nn
      a=an
      if(itest.eq.1) write(*,112) j,nn,n
112     format(1x,'in delps, j,nn,n=',3i10)
      derlp(j)=rg*(sin (a)**n)
      to=derlp(j)
      a=a+an
       k=1
      pp(1)=derlp(1)
 1      j=j+1
      if (j.gt.ispe) go to 4
      k=k+1
      if (k.gt.ispe) go to 3
      pp(k)=rg*sin(a)**n
      derlp(j)=pp(k)-pp(k-1)
      derlp(j)=abs (derlp(j))
      to=to+derlp(j)
      a=a+an
      if(to.lt.rg) go to 1
 2      no=j-1
      return
  3   write(*,101)k
      return
  4   write(*,102) j
      return
101   format(1x,'array pp in delps is too large ',i5)
102   format(1x,'array delp in delps is too large ',i5)
      end

      subroutine doublecr(pr, pi, c, rex, imx)
      real *8 pr, pi, rex, imx
      real c

      real *8 c2, u, x, r

      c2 = 1.0/c/c
      u  = c2 - (pr *pr - pi * pi)
      x  = -2.0*pr*pi
      r  = dsqrt(x*x+u*u)
      rex = dabs(r+u)/2.
      imx = dabs(r-u)/2.
      rex = dsqrt(rex)
      imx = -dsqrt(imx)
      return
      end

      complex function cr(p,c)
       common/testin/iflag1
      real*8 u,x,r,w1,w2,r1,r2
      complex p,cz
      iflag1=2
      cz=1./c**2-p*p
      if(iflag1.eq.1) write(*,100) c,p,cz
 100  format(1x,' c= ',e16.8,'p= ',2e16.8,'cz= ',2e16.8)
      u=real(cz)
      x=aimag(cz)
      r=dsqrt(x*x+u*u)
      if (iflag1.eq.1) write(*,101) u,x,r
 101  format(1x,' u= ',d16.8,'  x= ',d16.8,'  r= ',d16.8)
      w1=dabs(r+u)/2.
      w2=dabs(r-u)/2.
      r1=dsqrt(w1)
      r2=dsqrt(w2)
      if (iflag1.eq.1) write(*,102) w1,w2,r1,r2
 102  format(1x,' w1=',d16.8,'  w2= ',d16.8,/'r1=',d16.8,'r2=',d16.8)
      r11=r1
      r22=r2
      cr=r11-r22*(0.,1.)
      if (iflag1.eq.1) write(*,103) r11,r22,cr
 103  format(1x,'  r11=',e16.8,'  r22=',e16.8,'  cr=',2e16.8)
        return
      end

      function sf2(p,k)
      include 'psv.p'
c     real*8 p,psq,te,rcsq,rssq,esq,e ,f,fsq
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
c      real*8 c,s,d,th,x
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
      common/soray/msorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
      psq = p ** 2
      go to (1,2),nplnw
  2   te = 0.0
c      do 5 j = 1,k
      do 5 j = 2,k
      esq=abs(rcsq(j)-psq)
      fsq=abs(rssq(j)-psq)
      e=sqrt(esq)
      f=sqrt(fsq)
  5   te=te+th(j)*rlp(j)*rcsq(j)/(esq*e) +rls(j)*th(j)*rssq(j)/(fsq*f)
      sf2=te
c      write(*,*)'sf2=',sf2
c   computes second derivative of  t  with respect to  p
c   see a19 of paper bssa,feb. 1968
      return
  1     esq=abs(rcsq(nb)-psq)
      fsq=abs(rssq(nb)-psq)
      sf2=(p/x)*(rlp(nb)/esq+rls(nb)/fsq)
      write(*,100) nb,rlp(nb),rls(nb),esq,fsq,p,x,sf2
 100  format(1x,'insf2,nb,rlp(nb),rls(nb),esq,fsq,p,x,sf2=',/3i5,5e11.3)
      return
      end
       function ptim(p,k)
       include 'psv.p'
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
c     real*8 psq,rr,e,rcsq,rssq,c,s,d,th,x,f
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
      psq = p ** 2
      rr = 0.0
      do 5 j = 1,k
      e=sqrt(abs(rcsq(j)-psq))
      f=sqrt(abs(rssq(j)-psq))
 5     rr=rr+rlp(j)*th(j)*e +rls(j)*th(j)*f
      ptim = p*x + rr
      return
      end
      subroutine contor(tmx,m,kn,n,mo,tg)
      include 'psv.p'
       common/soray/nsorce,ntype,ndiret,tqc,tqs,tqd,jpt,nconjt,dconjt,nb
      common/cagcon/tstart,dp,nz,ndirt,nnend,dltp,nnn,dtim,mtd,dltm,nfm
     +,mkp
         common/conf/tfx
      common/magic/pp(imgc),ddpt(imgc),tt(imgc)
      common/spe/derlp(ispe),dd1,dd2,dd3,dd4,no
      common/pathc/po,to,k
      common/tfix/tn1,tn2,tn3,tn4,jn1,jn2,jn3,jn4,det1,ntntt
      complex pt,ddpt,dev
      complex pp,p,ct
      common/scase/ispecial
c     real*8 q,pi,tt,po,to,tg
c     dtim=nearest point to(to)
c       dltm=max. spacing between tt(j) points
c     mtd=rate of approaching this limit, see statement  14
c      real*8 q, pi, dl
        det=(dtim/10.)*po
        dett=det
        jn=10
        deltm=1.e-6
        
 1     q=po
       i=mo-1
       pil=1.e-12
       dett=det
       dttm=dtim
       if(dtim.le..001 )dttm=dtim
       if(tg.le.tn3) det=.01*det
       pp(i)=po
       tt(i)=to
        njtm=0
       km=105

       if(derlp(1).le.deltm) derlp(1)=deltm
 33     i=km-1
        l=0
	if(ispecial .eq. 1) jn = no
       do 10 j=1,jn
           l=l+1
           i=i+1
           q=q+derlp(j)
           pi=po*0.1
           dl=pi*.45
c           write(*,*) "calling this time2 i=",i
           call time2(q,pi,dl,p,dev,ct,kn,n,pil,det,njtm)
           timei=-ct*(0.,1.)
           dr=real(dev)
           ddpt(i)=dev
           tt(i)=ct
           pp(i)=p
           if(tt(i)-tt(i-1).le.dltm) mkp=i
	   if(ispecial .eq. 0) then
               if(tt(i)-to.le.dttm) then
c                  commented out by Lianxing Wen
*		   write(*,*)"tt(i) - to < dttm goto 32 "
c               jj=j+1
c               pt=(po-pp(i))/2.
c               derlp(jj)=real(pt)
		   go to 32
               endif

               if(dr.le.0.) go to 30
c               if(abs(timei).gt.det  ) go to 30
               if(tt(i)-to.lt.0.) then
		   write(*,*)"tt(i) - to < 0 goto 30"
		   go to 30
               endif
               jj=j+1
               pt=(po-pp(i))/2.
               derlp(jj)=real(pt)
           endif
 10     continue
       if(ispecial .eq. 0) then
 30        if(i.le.km) derlp(1)=5.*derlp(1)
           if(i.le.km) go to 33
       endif

       i=i-1
32      jj=i-km+1
       do 31 j=1,jj
         ll=mo+j-1
         nn=i-j+1
         if(mkp.eq.nn) mkp=ll
         tt(ll)=tt(nn)
         pp(ll)=pp(nn)
         ddpt(ll)=ddpt(nn)
c        write(*,*)"ll=", ll, " nn=", nn," pp", pp(ll)
 31    continue
      i=ll
        pt=pp(ll)-pp(ll-1)
        derlp(jn)=real(pt)
         j=jn
         mm=1
         tm=tmx+tstart
         dev=pp(ll)
        pi=aimag(dev)
        q=real(dev)
         det=dett
 12     q=q+mm*derlp(jn)
        tg=tt(i)-to
        i=i+1
        dev=pp(i-1)-pp(i-2)
c        write(*,*)"i-1=", i-1, " pp", pp(i-1), " pp=",pp(i-2)
       
        delpr=real(dev)
        delpi=aimag(dev)
        delpi=abs(delpi)
        pi=(delpi*mm*derlp(j))/delpr +pi
        if((nconjt.gt.0).and.(tg.gt.dconjt)) njtm=2
       if((nconjt.gt.0).and.(tg.gt.dconjt)) go to 17
        dl=pi*.45
 17     call time2(q,pi,dl,p,dev,ct,kn,n,pil,det,njtm)
c           write(*,*) "that time2 ct ", ct, " p=",p, "dev=",dev
        pp(i)=p
        ddpt(i)=dev
         tt(i)=ct
          ctt=ct-to
        if ((tt(i).gt.tm).and.(ctt.gt.dp)) go to 13
14      if(tt(i)-tt(i-1).le.dltm) mm=mtd*mm
       if(tt(i)-tt(i-1).le.dltm) mkp=i
      go to 12
 13      continue
      m=i
c     write(*,500) m,mkp,tt(mkp)
 500  format(1x,'in contor,m,mkp,tt(mkp)=',2i5,e15.5)
        return
      end
      function ynterp(xa,ya,x,l,n)
      dimension xa(*),ya(*),za(21)
      common /yntcom/ ifd,smx
      logical smx
      data sg1/0.0/

      if(sg1 .ne. 0.0) go to 12
      sg1 = 1.0
      smx = .false.
      ifd = 0
   12 sg2 = 0.0
      if(smx) go to 44
      iu = l + 1
      if(l .lt. 2) go to 43
      nin = n
      if(n .lt. 1 .or. n .gt. 20) go to 40
      np1 = n + 1
      if (np1 -l) 14,48,43
   14 if(xa(l) - xa(1))16,43,18
   16 sg2 = 1.0
      il = iu
      iu = 0
      go to 20
   18 il = 0
   20 ia = (iu + il)/2
      if(iabs(iu-il) .lt. 3) go to 26
      if(x - xa(ia))22,50,24
   22 iu = ia
      go to 20
   24 il = ia
      go to 20
   26 if(x .eq. xa(ia)) go to 50
      ifd = ia
      ia = ia - np1/2
      if(mod(nin,2) .eq. 0) go to 30
      if(x .gt. xa(ifd)) go to 28
      if(sg2 .ne. 0.0) ia = ia + 1
      go to 30
   28 if(sg2 .eq. 0.0) ia = ia + 1
   30 if(ia .lt. 1) ia = 1
      if(l - nin .lt. ia) ia = l - nin
   32 ja = ia
      do 34 j = 1,np1
      za(j) = ya(ja)
   34 ja = ja + 1
      ib = ia
      do 38 i = 1,nin
      t1 = x - xa(ib)
      do 36 j = i,nin
      ja = j + ia

1     t2 = xa(ja) - xa(ib)
c      if(t2.lt.0.0005) write(*,1001) ja,ib,xa(ja),xa(ib),t2
c1001  format(1x,'in ynterp,ja,ib,xa(ja),xa(ib),t2=',2i5,3e17.7)
      if(t2. lt. 1.e-30)then
          ja = ja+1
          goto 1
      endif    

   36 za(j+1) = za(i) + t1*(za(j+1) - za(i))/t2
   38 ib = ib + 1
      iflag = 1
      ynterp = za(nin+1)
      return
   40 write(*,42)
   42 format(1x,/' error return from ''ynterp''--n is out of range.')
   43 ifd = 1
      iflag =3
      ynterp = ya(1)
      return
   44 smx = .false.
      if(iflag - 2) 32,46,43
   46 ynterp=ya(ia)
      ifd = ia
      return
   48 ia = 1
      ifd = 1
      go to 32
   50 iflag = 2
      go to 46
      end
      subroutine interp(xp,yp,n,x,y)
      dimension xp(n),yp(n)
      real dif1,dif2,dify,dr
    1 if (x.gt.xp(n))  go to 6
      if (x .lt. xp(1))  go to 6
    2 do 10  i=1,n
      if (xp(i)-x) 10,102,3
   10 continue
    3 k=i-1
      dif1=xp(i)-xp(k)
      dif2=xp(i)-x
      ratio = dif2/dif1
      dify=abs(yp(i) - yp(k))
      dr = dify*ratio
      if (yp(i) .gt. yp(k))  go to 4
    5 y=yp(i) + dr
      return
    4 y=yp(i)-dr
      return
  102 y=yp(i)
      return
    6 y=0.0
      return
      end
      subroutine find2 (q,k,del,det,po,to)
      include 'psv.p'
       common/trace/tracing
      common/travel/rlp(icc),rls(icc),ll,fdp,map,nplnw
      common/stuff/c(icc),s(icc),d(icc),th(icc),rcsq(icc),rssq(icc),x
      common/tstar/tq,tvl,qa(icc),qb(icc)
       common/lprint/prnt,prnts,prntc,flat,prcon,prampl
      logical prnt,prnts,prntc,flat,prcon
       real*8 e(1000),f(1000),p1,p2,xx,yy,p
      tde = del
1000  format(1x,'entering find2, qa,qb=',4e12.4)
    8 p = q
	 p1=q
	 p2=del
 5      p=(p1+p2)*0.5
      psq = p ** 2
      bltem = 0.0
      do 10 j = 1,k
	 xx=dabs(1./s(j)-p)
	 yy=dabs(1./s(j)+p)
      f(j) =dsqrt(xx*yy)
	 yy=dabs(1./c(j)+p)
	 xx=dabs(1./c(j)-p)
      e(j) =dsqrt(xx*yy)
1001   format(1x,'j,psq,rssq,rcsq,e,f=',i3,5e12.4)
c     write(*,1001) j,psq,rssq(j),rcsq(j),e(j),f(j)
 10     bltem=bltem-rlp(j)*th(j)/e(j)-rls(j)*th(j)/f(j)
      bl = x + bltem*p
       if(abs (p2-p1).le.1.e-15) go to 1
    6 if (abs(bl).le.1.e-14) go to 1
	 if(bl.gt.0.0) then
	 p1=p
	 else
	 p2=p
	 end if
       go to 5
 1    po = p
      totem = 0.0
	do 11 j=1,k
 11	totem=totem+e(j)*th(j)*rlp(j)+f(j)*th(j)*rls(j)
	to=p*x+totem
	if(tracing .gt. 0)then
             write(*,17) po, to, bl
        endif     
  17  format(4x,'po= ',e18.6,10x,'to= ',e18.6,10x,'bl= ',e18.6)
      return
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
      subroutine convt(x,nx,z,d,n,l2n,nx1,ncent,fscl)
      dimension x(*),z(*),c(40200),d(*)
      nz=nx
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
      call colb(l2n,c,-1.)
      do 6 i=1,ncent
       j2=2*i
      j1=j2-1
      e1=c(j1)*d(j1)-c(j2)*d(j2)
      e2=c(j1)*d(j2)+c(j2)*d(j1)
      c(j1)=e1
6     c(j2)=e2
      call conj(c,ncent)
      call colb(l2n,c,1.)
      do 7 i=1,nz
      j1=2*i-1
7     z(i)=c(j1)*fscl
      return
      end
