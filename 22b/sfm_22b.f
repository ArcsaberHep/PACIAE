cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sfm(iiii,ijk,ss,kfa,kfb)
c	perfome the hadronization by calling 'pyexec' (string fragmentation)
c	it was written by Ben-Hao Sa on 31/07/02
c	its input messages are in 'pyjets'
c	its internal wooking block is 'sa1_h' 
c	its output message is in 'pyjets' ('sa1_h' the same)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000,mplis=40000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 161007
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)   ! 141208
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c      COMMON/PYSUBS/MSEL,NONSUB,MSUB(200),KFIN(2,-40:40),CKIN(200)
c      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   ! 221203 240407
c080104
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21   ! 051108
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
c080104
	common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
c	arraies in above statement are for hadronized and 
c	 decayed particles used in hadronic cascade processes
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio 
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 240209  
        character chaf*16,chaua*16,chaub*16   ! 141208
	dimension peo(4),rc(3)

c	ich1=0.
c	do i1=1,n
c	kf=k(i1,2)
c	ich1=ich1+pychge(kf)
c	enddo
c	write(mstu(11),*)'in hadniz before, sum of charge=',ich1/3   !
c	call pylist(1)
	rrp=1.16
	nn=0
	do i1=1,kszj
	do j1=1,5
	kn(i1,j1)=0
	pn(i1,j1)=0.
	rn(i1,j1)=0.
	enddo
	enddo
c051108
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        itime=0
        gtime=0.
        astr=0.
        do i1=1,5
        akapa(i1)=0.
        enddo
        endif
c051108
	mstj(1)=mstj1_2
c	write(9,*)'in hadniz iii,mstu,mstj(1)-(3)=',iii,mstu(21),mstj(1),
c     c	 mstj(2),mstj(3)   
	mstj(21)=0
c       particle decay is inhibited
c	produced hadron from calling 'pyexec' is arranged at the position   
c	 of parent, decayed hadrons do not have proper position so 
c	 we inhibite first the decay
	call pyexec
c141208
        if(n.eq.0)then
c       fragment that hh collision pair by pythia directly
        call pyname(kfa,chaua)
        call pyname(kfb,chaub)
c       write(22,*)'cha,chb=',chaua,chaub
        call  pyinit('cms',chaua,chaub,ss)
        call pyevnt
        nbh=0   ! 111210
        ijk=1
c240209
        n44=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=44   ! '44': prompt direct photon
        n44=n44+1
        endif
        enddo
c       move "44" from 'pyjets' to 'sgam'
        if(n44.gt.0)call remo_gam(44)
c240209
        endif
c141208
c051108
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        do i1=1,n
        if(k(i1,2).eq.92)astr=astr+1.
        enddo
c       parj(1)=parj1
c       parj(2)=parj2
c       parj(3)=parj3
c       parj(21)=parj21
        atime=dfloat(itime)
        if(atime.gt.0.)then
        akapa(1)=akapa(1)/atime
        akapa(2)=akapa(2)/atime
        akapa(3)=akapa(3)/atime
        akapa(4)=akapa(4)/atime
        akapa(5)=akapa(5)/atime
        gtime=gtime/atime
c       gtime: averaged # of gluons in a string in current event
        endif
c       write(9,*)'af call luexec and kjp22,n=',n   !
        endif
c051108
        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813
c	call pylist(1)
c	ich1=0.
c	do i1=1,n
c	kf=k(i1,2)
c	ich1=ich1+pychge(kf)
c	enddo
c	write(mstu(11),*)'after, sum of charge=',ich1/3   !
c080104
c       'sbh' to 'pyjets'
        if(nbh.ge.1)then
        do l=1,nbh
        l1=n+l
        do m=1,5
        k(l1,m)=kbh(l,m)
        p(l1,m)=pbh(l,m)
        v(l1,m)=vbh(l,m)
        enddo
        enddo
        n=n+nbh
	endif
c	write(9,*)'be decay n=',n
c	write(22,*)'be decay'
c	call pylist(1)
c080104

c	return   ! temporal

c	'pyjets' to 'sa1_h'
	nn=n
        do j1=1,n
        do j2=1,5
	kn(j1,j2)=k(j1,j2)
        pn(j1,j2)=p(j1,j2)
	rn(j1,j2)=v(j1,j2)
        enddo
	enddo

c	do i=1,nn
c	write(9,102)(rn(i,j),j=1,4)   ! sa
c	enddo
c	arrange produced particles on the surface of sphere (radius rrp)
c        and centred on parent position (produced particle 
c	 is put on its parent position originally)
	ipp=1
100	ip=0
	r1=rn(ipp,1)	
	r2=rn(ipp,2)	
	r3=rn(ipp,3)
	rc(1)=r1
	rc(2)=r2
	rc(3)=r3
	rn(ipp,4)=0.
c	in corresponding with the time set in 'posi'
c	find out the hadrons with common parent (rc)
c	note: produced particles with same position are arranged continuously
c        in pyjets
	do j2=ipp+1,n
	s1=rn(j2,1)	
	s2=rn(j2,2)	
	s3=rn(j2,3)
	if(dabs(s1-r1).le.1.e-6 .and. dabs(s2-r2).le.1.e-6 .and.
     c	 dabs(s3-r3).le.1.e-6)ip=ip+1
	enddo
	ipp1=ipp+ip
	call posi(ipp,ipp1,rc,rrp)
	ipp=ipp1+1
	if(ipp.lt.n)goto 100
c	write(9,*)'after rearrangement, nn=',nn   ! sa
c	do i=1,nn
c	write(9,102)(rn(i,j),j=1,4)   ! sa
c	enddo

c	transfer four position messages from 'sa1_h' to 'pyjets'
	do i=1,nn
	do j=1,4
	v(i,j)=rn(i,j)
	enddo
	enddo
c	write(9,*)'af position n=',n

c       decay of unstable hadrons
	call decayh(rrp)
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine decayh(rrp)
c       decay of unstable hadrons
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000,mplis=40000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 161007
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 240209
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio   ! 060813
        dimension rc(3)
c	particle decay before rescattering is set in paciae.f

c	find out a unstable hadron 
c080104	nnn=n
	jd=0
c	jd: statistics of number of unstable hadrons adjudged
c110604	200	continue ! nnn1=nnn-jd 080104
c080104	if(nnn1.le.0)goto 500
	i1o=1   ! 080104
200	do 400 i1=i1o,n    ! 110604
	kf=k(i1,2)
	md=mdcy(pycomp(kf),1)
	if(md.eq.1)then   ! if 1
c	write(9,*)'iii,nnn,jd,nnn1=',iii,nnn,jd,nnn1
c	write(9,*)'i1,kf,md=',i1,kf,md
	jd=jd+1
c	decay of unstable hadron i1
	call pydecy(i1)
c	'pyjets' is filled up simultaneously 
c	remove decaying hadron from 'pyjets'
        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813
c	write(22,*)'i1=',i1
c	call pylist(1)
c	store the position of decaying hadron
	do i=1,3
	rc(i)=rn(i1,i)
	enddo
c	write(9,*)'rc=',(rc(i),i=1,3)   ! sa
c       remove decaying hadron from 'sa1_h' as well
c	 i. e. move the hadron list ('sa1_h') 1 step downward from i1+1 
c	 to nn
	call updah(nn,i1+1,1)
c	write(9,*)'nn=',nn   ! sa 
	nn=nn-1
c	write(9,*)'nn=',nn   ! sa 
	nn1=nn   ! decayed particles are located above nn1, 240209
c240209
c       move "22" from 'pyjets' to 'sgam'
        jb1=0
700     do i3=nn1+jb1+1,n
        kf=k(i3,2)
        if(kf.ne.22)then
        jb1=jb1+1
        goto 800
        endif
        ngam=ngam+1
        do i2=1,5
        kgam(ngam,i2)=k(i3,i2)
        pgam(ngam,i2)=p(i3,i2)
        vgam(ngam,i2)=v(i3,i2)
        enddo  
        if(i3.eq.n)then      
        n=n-1
        goto 900
        endif
c       move particle list 'pyjets' one step downward from i3+1 to n
        do j=i3+1,n
        j1=j-1
        do jj=1,5
        k(j1,jj)=k(j,jj)
        p(j1,jj)=p(j,jj)
        v(j1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 700
800     enddo
900     continue
c240209

	nn=n
c	write(9,*)'nn1,nn,n=',nn1,nn,n   ! sa

c	fill produced hadrons (from decay, no gamma) into 'sa1_h'
	do i=nn1+1,nn
	do j=1,5
	kn(i,j)=k(i,j)
	pn(i,j)=p(i,j)
	rn(i,j)=v(i,j)
	enddo
        enddo

c	write(9,*)'after decay, i1,kf,n,nn,nn1,jd=',i1,kf,n,nn,nn1,jd! sa
c	do i=1,nn
c	write(9,102)(rn(i,j),j=1,4)   ! sa
c	enddo

c	arrange decaied hadrons (from nn1+1 to nn) on the surface of
c        sphere with radius rrp and centred on parent position
	call posi(nn1,nn,rc,rrp)
c	write(9,*)'after decay and rearrangement,n,nn,nn1=',n,nn,nn1! sa
c	do i=1,nn
c	write(9,102)(rn(i,j),j=1,4)   ! sa
c	enddo
c       transfer four position messages of decaied hadrons to 'pyjets'
	do i=nn1+1,nn
        do j=1,5
        v(i,j)=rn(i,j)
        enddo
        enddo
c	write(9,*)'after decay and rearrangement,v,n,nn,nn1=',n,nn,nn1! sa
c	do i=1,n
c	write(9,102)(v(i,j),j=1,4)   ! sa
c	enddo
c	write(mstu(11),*)'after decay'
c	call pylist(1)

c	if(jd.eq.2)goto 500   ! temporal
	i1o=i1   ! 080104
	goto 200   ! 300   110604
	endif   ! endif 1 
400	continue
c110604	goto 500
c110604	300	goto 200	

500	continue
102     format(4(1x,f10.4))
c	if(nout.eq.1 .or. iii.eq.1 .or. iii.eq.neve)then
c	write(mstu(11),*)'after hadronization and decay'
c	call prt_sa1_h(nn)
c	endif
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updah(j2,jc,i)
c	move the hadron list i steps downward from jc till j2
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000)
	common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
	do j=jc,j2
	do jj=1,5
	kn(j-i,jj)=kn(j,jj)
	pn(j-i,jj)=pn(j,jj)
	rn(j-i,jj)=rn(j,jj)
	enddo
	enddo
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine posi(nn1,nn2,rc,rrp)
c	arrange produced particles (from nn1+1 to nn2) on the surface of 
c	 sphere with radius rrp and centred on parent position
c	rc : the coordinate of center of the parent
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (mplis=40000)
        PARAMETER (KSZJ=40000)
	common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
	dimension rc(3),rr(3,kszj)
	iii=0
	do 100 ii=nn1+1,nn2
200	call samp1(rrp,ii,rr)
c	non-overlapping demand among baryons
	if(ii.eq.1)goto 100
	if(iabs(kn(ii,2)).lt.1000.or.iabs(kn(ii,2)).gt.10000)goto 100
	do 300 jj=1,ii-1
	if(iabs(kn(jj,2)).lt.1000.or.iabs(kn(jj,2)).gt.10000)goto 300
	dcc=dsqrt((rr(1,jj)-rr(1,ii))**2+(rr(2,jj)-rr(2,ii))**2
     c	 +(rr(3,jj)-rr(3,ii))**2)
	if(dcc.lt.0.8)then
	iii=iii+1
	if(iii.lt.10000)goto 200
c	write(9,*)'ii,jj,rrp,dcc=',ii,jj,rrp,dcc
c	write(9,*)'rii=',(rr(1,ii),i=1,3)	
c	write(9,*)'rjj=',(rr(1,jj),i=1,3)
c	write(9,*)'non-overlapping demand violating'
	goto 100   ! 121699	
c072400	stop 10000
	endif
c	0.8 is the hard core of baryon
300	continue
100	continue
	do ii=nn1+1,nn2
c	give zero to the time of hadronized and decayed particles
	rn(ii,4)=0.
	do jj=1,3
	rn(ii,jj)=rr(jj,ii)+rc(jj)
	enddo
	enddo
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine samp1(xf,i,rr)
c       give the position,on the surface of sphere with
c	 radius xf,to particle i
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000)
	dimension rr(3,kszj)
        cita=2*pyr(1)-1.
        fi=2.*3.1416*pyr(1)
        sita=dsqrt(1.-cita**2)
        rr(1,i)=xf*sita*dcos(fi)
        rr(2,i)=xf*sita*dsin(fi)
        rr(3,i)=xf*cita
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa1_h(nn1)
c       print particle list and the sum of four momentum and charge
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        dimension peo(4)
        do i=1,nn1
c	write(mstu(11),*)i,kn(i,2),(pn(i,j),j=1,4)
	write(9,*)i,kn(i,2),(pn(i,j),j=1,4)
        enddo
        call psum(pn,1,nn1,peo)
        ich1=0.
        do i1=1,nn1
        kf=kn(i1,2)
        ich1=ich1+pychge(kf)
        enddo
c	write(mstu(11),*)peo,ich1/3   !
	write(9,*)peo,ich1/3   !
        return
        end



c051108ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine strtension(ip,np)
c       calculate the effective string tension after tuning parj(1),(2),
c        (3),(21) to the rapidity density etc. for each string  
c       the string takes up items from n+1 to n+np in "pyjets" 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000,KSZ1=30)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21   ! 051108
        vfr24=3.5   ! parameter alpha
        vfr25=0.8   ! \sqrt(s_0) in GeV
        toteng=0.0
        toten=0.0
        totglukt=0.0
        pmax=0.
        ggg=0.
c 	do i=n+1,n+np
        do i=ip,ip+np-1
 	toten=toten+p(i,4)   ! s, string total energy
 	pp=dsqrt(p(i,1)**2+p(i,2)**2)
 	if(k(i,2).eq.21.and.pp.gt.vfr25)then
        toteng=toteng+dlog(pp/vfr25)   ! sum over gluons in a string
        ggg=ggg+1.
        endif	
 	if(k(i,2).eq.21.and.pp.gt.pmax)pmax=pp   ! k_{Tmax}^2
 	enddo
        if(pmax.gt.vfr25)totglukt=totglukt+dlog(pmax/vfr25)   ! numerator	
 	ss=dlog(toten/vfr25)+toteng   ! denominator
 	effk2=(1.-totglukt/ss)**(-vfr24)
c       string tension of the pure qqbar string, kapa0, is assumed to be 1 	
	akapa(1)=akapa(1)+effk2
	itime=itime+1
        gtime=gtime+ggg
	return
	end



c051108ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine strtension1(ip,np)
c       calculate the effective string tension for each string and then
c       calculate the new parj(1),parj(2) etc. in the current event.
c       when this subroutine is called in pythia the string has been boosted
c       to its own cms frame. the whole string takes up items from n+1 to n+np
c       in pyjets
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000,KSZ1=30)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)         
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21   ! 051108
        common/sa29/effk1,lcub   ! 051108
        vfr24=3.5   ! parameter alpha
        vfr25=0.8   ! \sqrt(s_0) in GeV
        toteng=0.0
        toten=0.0
        totglukt=0.0
        pmax=0.
        ggg=0.
ch      write(9,*)'in strtension1 ip,np,n,effk1=',ip,np,n,effk1   !
c        write(9,*)'in strtension1 ip,np,n,effk1=',ip,np,n,effk1   !
c 	do i=n+1,n+np
        do i=ip,ip+np-1
 	toten=toten+p(i,4)   ! s, string total energy
 	pp=dsqrt(p(i,1)**2+p(i,2)**2)
 	if(k(i,2).eq.21.and.pp.gt.vfr25)then
        toteng=toteng+dlog(pp/vfr25)   ! sum over gluons in a string
        ggg=ggg+1.
        endif	
 	if(k(i,2).eq.21.and.pp.gt.pmax)pmax=pp   ! k_{Tmax}^2
c       write(9,*)'i,toten,pp,toteng,ggg=',i,toten,pp,toteng,ggg   !
 	enddo
        if(pmax.gt.vfr25)totglukt=totglukt+dlog(pmax/vfr25)   ! numerator	
 	ss=dlog(toten/vfr25)+toteng   ! denominator
        div=totglukt/ss
 	effk2=(1.-div)**(-vfr24)
c       string tension of the pure qqbar string, kapa0, is assumed to be 1 	
c       write(9,*)'pmax,nu.,de.,div.,vfr24,effk2=',pmax,totglukt,ss,div,
c    c   vfr24,effk2   !
	akapa(1)=akapa(1)+effk2
 	parj(21)=parj21*((effk2/effk1)**(0.5))
 	parj(1)=parj1**(effk1/effk2)
 	parj(2)=parj2**(effk1/effk2)
 	parj(3)=parj3**(effk1/effk2)
	akapa(2)=akapa(2)+parj(2)
	akapa(3)=akapa(3)+parj(21)
	akapa(4)=akapa(4)+parj(1)
	akapa(5)=akapa(5)+parj(3)
	itime=itime+1
        gtime=gtime+ggg
ch       write(9,*)'2 old parj1,parj2,parj3,parj21=',parj1,parj2,parj3
ch     c   ,parj21   !
ch       write(9,*)'new=',parj(1),parj(2),parj(3),parj(21)   !
c       write(9,*)'akapk=',akapa(1),akapa(2),akapa(3),akapa(4),akapa(5)   !
ch       write(9,*)'2 out of strtension ip,np,n=',ip,np,n   !
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updad_pyj(j2,j1,i)
c       move the parton list (pyjets) i steps downward from j1 to j2
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        do j=j1,j2
        do jj=1,5
        k(j-i,jj)=k(j,jj)
        p(j-i,jj)=p(j,jj)
        v(j-i,jj)=v(j,jj)
        enddo
        enddo
c        n=n-i
        return
        end

