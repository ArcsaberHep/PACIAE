cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sfm
c	perfome the hadronization by calling 'pyexec' 
c	it was written by Ben-Hao Sa on 31/07/02
c	its input messages are in 'pyjets'
c	its internal wooking block is 'sa1_h' 
c	its output message is in 'pyjets' ('sa1_h' the same)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=80000,mplis=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 161007
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c      COMMON/PYSUBS/MSEL,NONSUB,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)   ! 010418
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   ! 221203 240407
c080104
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3,
     c   parj21,parj4,adiv,gpmax,nnc   ! 051108 070417 010518
        common/sa34/itorw,iikk,cp0,cr0,kkii   ! 050920
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
c080104
	common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
c	arraies in above statement are for hadronized and 
c	 decayed particles used in hadronic cascade processes
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 250209
c141208
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio   
c       ithroq : the total # of quarks thrown away
c       ithrob : the total # of antiquarks thrown away
c       throe : total momentum and energy of the partons thrown away
c       ithroc : total charge of the partons thrown away
c141208
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
c151021
        if(ipden.eq.2 .and. itden.eq.2)then
        mstj(1)=1
        else
        mstp(111)=1
        endif
c151021 
	call pyexec
        if(kkii.eq.2)return   ! 050920
        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813

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

c	arrange produced particles on the surface of sphere (with radius rrp
c        and centred on parent position),produced particle 
c	 is put on its parent position originally
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

c	find out the hadrons with same parent (rc)
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
c070417	call decayh(rrp)
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine decayho(rrp)
c       decay of unstable hadrons
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=80000,mplis=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 161007
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 250209
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
	nn1=nn   ! decayed particles are located above nn1, 250209
c250209
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
c250209

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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine decayh(rrp)
c       decay of unstable hadrons
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=80000,mplis=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 161007
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 250209
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio   ! 060813
        dimension rc(3)

c	find out unstable hadron 
	jd=0
c	jd: statistics of number of unstable hadrons adjudged
	i1o=1   ! 080104
	nn1=n   ! 131218
200	continue
	do 400 i1=i1o,nn1    ! 110604 141218
	kf=k(i1,2)
	md=mdcy(pycomp(kf),1)


	if(md.eq.1)then   ! 1
c	write(9,*)'i1,kf,md=',i1,kf,md
	jd=jd+1

c	decay of unstable hadron i1
	call pydecy(i1)
c	'pyjets' is filled up simultaneously 
c141218	decayed particles are located above nn1

c	store the position of decaying hadron
	do i=1,3
	rc(i)=rn(i1,i)
	enddo
c	write(9,*)'rc=',(rc(i),i=1,3)   ! sa

c	update particle list

c	remove decaying hadron from 'pyjets'
        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813
c	write(22,*)'i1,kf=',i1,k(i1,2)
c	call pylist(1)

c       remove decaying hadron from 'sa1_h' as well
c	 i. e. move the hadron list ('sa1_h') 1 step downward from i1+1 
c	 to nn
	call updah(nn,i1+1,1)
c	write(9,*)'nn=',nn   ! sa 
	nn=nn-1
c	write(9,*)'nn=',nn   ! sa 
	nn1o=nn1   ! 141218
	nn1=nn   

c       move "22" (decay gamma) from 'pyjets' to 'sgam'
	jb1=0
700	continue
	do 800 i3=nn1+jb1+1,n
        kf=k(i3,2)
        if(kf.ne.22)then   !!
	jb1=jb1+1
        goto 800

        elseif(kf.eq.22)then   !!
        ngam=ngam+1
        do i2=1,5
        kgam(ngam,i2)=k(i3,i2)
        pgam(ngam,i2)=p(i3,i2)
        vgam(ngam,i2)=v(i3,i2)
        enddo
c	update particle list
        if(i3.eq.n)then   !!!   
        n=n-1
        goto 900
        elseif(i3.lt.n)then   !!!
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
        goto 700   ! 131218
	else   !!!
	endif   !!!

	else   !!
	endif   !!
800     enddo
900	continue

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

	if(i1.eq.nn1o)goto 401 !141218, originally no this statement->dead loop 
	i1o=i1   ! 131218 
	goto 200   ! 131218
	endif   ! 1


400	enddo
401	continue   ! 141218
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
      PARAMETER (KSZJ=80000)
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
        parameter (mplis=80000)
        PARAMETER (KSZJ=80000)
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
	parameter (kszj=80000)
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
        parameter (kszj=80000)
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updad_pyj(j2,j1,i)
c       move the parton list (pyjets) i steps downward from j1 to j2
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=80000)
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
