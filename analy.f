        subroutine analy(nmin,nminf,ncha,nchaf)   
c       analyses an event based on the messages in 'pyjets' 
c       it is composed by Ben-Hao Sa on 28/12/19
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     c  iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
	common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
	common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph
     c	 ,vnlep   ! 260314
	common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003
        common/sa30/vneump,vneumt,mstptj
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        common/anly1/an(40,5,20),bn(20),anf(40,5,20),bnf(20)
        real nmin,nminf,ncha,nchaf
	dimension c(5),dinel(600),dineli(600),sthroe(4),wthroe(4)

c       ispmax: maximum # of particle KF code wanted to statistics
c	ispkf: array of particle KF code 
c       ispkf(i): KF code of i-th particle
c	kfmax: dimention of ispkf array, i.e. maximum # of KF code considered   
c	bn(i) (bnf(i)): multiplicity of i-th particle in partial (full)
c        phase space
c	an(l,i,j) (anf(l,i,j)): three dimension array in partial (full) 
c         phase space
c	 l: order # of y (pt, eta, ...) interval
c	 i: kind of distribution 
c         for pp, NA,AN and AB collisions
c	   i=1: for y
c	   i=2: for pt
c           .      . 
c           .      .
c           .      .
c         for lp and lA collisions
c          i=1 : z
c          i=2 : \nu 
c          i=3 : Q^2
c           .        .
c           .        .
c           .        .
c	 j: order # of particle in array of 'ispkf' (1 to ispmax)
c	isdmax: maximum # of distributions considered
c	asd(i): interval segmented for i-th distribution  
c	iflmax: maximum # of windows, =0 means no window at all
c	afl(j,i,1): lower bound of i-th window for j-th particle
c	afl(j,i,2): upper bound of i-th window for j-th particle
c260314  for pp,NA,AN and AB collisions
c 	  i=1: y window
c	  i=2: pt window
c         .        .
c         .        . 
c         .        .
c260314  for lp and lA collisions
c         i=1 : Q^2=-q^2 (fq2 named in program) window
c         i=2 : W^2 (w2l) window
c         i=3 : y (yyl) window
c         i=4 : P_h (pph) window
c         i=5 : z (zl) window
c         .        .
c         .        .
c         .        .

c       for lp and lA collisions
c       pincl (pscal): four momentum and mass of incident (scatterd) lepon
c       pinch: four momentum and mass of incident hadron
c        vnu: \nu; fq2: Q^2=-q^2; w2l: W^2; yyl: y; zl: z; xb: x_B; pph: P_h

        dpmax=adj1(27)   ! largest momentum allowed for particle
        adj140=adj1(40)
c        print*,'adj140=',adj140

c	initiates the variales
	nmin=0
	nminf=0
        ncha=0
        nchaf=0
	do i1=1,20
	bn(i1)=0.
	bnf(i1)=0.
        enddo   
        do i1=1,40   
	do i2=1,5
	do i3=1,20
	an(i1,i2,i3)=0.
	anf(i1,i2,i3)=0.
	enddo
	enddo
	enddo

c	analyses an event
	do 400 j=1,n
	ik=k(j,2)
        plu6=pyp(j,6)
	eta=pyp(j,19)
	w=1.
	p1=p(j,1)
	p2=p(j,2)
	p3=p(j,3)
	p4=p(j,4)
	p5=p(j,5)   ! sa 26/05/99
c       write(99,191) ik, p1, p2, p3, p4   ! 240119
	yy=pyp(j,17)
	ppt=pyp(j,10)
c281104
	ppm=dsqrt(p1*p1+p2*p2+p3*p3)
	if(ppm.le.dpmax.and.p4.le.dpmax)then
	goto 3000
	else
	goto 400   ! through away that particle
	endif
c281104
3000	continue
c060718	if((itden.eq.0.and.ipden.eq.1).or.(itden.eq.1.and.ipden.eq.0)
c060718     c	 .or.(itden.eq.1.and.ipden.eq.1))then   ! 260314
	if(ipden.lt.11)then   !  260314 060718 for pp,pA (Ap) & AB
	c(1)=yy   
c150622 if(ifram.eq.1)c(1)=eta   
	c(2)=ppt
c	.
c	.
c	.
	kkk=1
c033101
c	statistics of negative multiplicity
	if(adj140.ge.3 .and. plu6.lt.-0.9)then   ! for hadron 140414
	nminf=nminf+1
	do i=1,iflmax
	if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 700
	enddo
	nmin=nmin+1	
700	endif
	if(adj140.lt.3 .and. plu6.lt.-0.2)then   ! for parton 140414
        nminf=nminf+1   ! -plu6 230206
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 702
        enddo
        nmin=nmin+1   ! -plu6 230206
702     endif
c010220
c033101	statistics of positive multiplicity
        if(adj140.ge.3 .and. plu6.gt.0.9)then   ! for hadron 140414
        nchaf=nchaf+1   ! +plu6 
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 701
        enddo
        ncha=ncha+1   ! +plu6 
701     endif
	if(adj140.lt.3 .and. plu6.gt.0.2)then   ! for parton 140414
c070802
        nchaf=nchaf+1   ! +plu6 230206
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 703
        enddo
        ncha=ncha+1   ! +plu6 230206
703     endif
	endif   ! 260314

c       statistics of y, pt, ... distributions (for NA,AN,and BA); z, \nu,
c        ... distributions (for lp and -lA)
	do 500 kk=1,ispmax
	kf=ispkf(kk)
	if(ik.ne.kf)goto 500

c????????????????????????????????????????????????????????????????????
c       excludes the projectile and the target spectator nucleons
        if(ik.eq.2212 .or. ik.eq.2112)then
c	eep=dsqrt(win*win+0.938*0.938)
c	if(ifram.eq.1)eep=0.5*win
c	eepd=eep-0.001
c	eepu=eep+0.001
c	if(ifram.eq.0 .and. ((p4.gt.eepd .and. p4.lt.eepu) .or.
c     c   p4.le.0.940))then   ! 111899
c	goto 500
c	endif   ! 111899
c	if(ifram.eq.1 .and. (p4.gt.eepd .and. p4.lt.eepu))then
	if(ppt.le.1.e-5)goto 500 ! 310521 e-4->e-5
        endif
c????????????????????????????????????????????????????????????????????

c       analyses for pp,pA and AB
	if(ipden.lt.11)call stati_h(yy,ppt,eta,p5,ik,kk,w,bn,an,bnf,anf)  
c       analyses for lp and lA
	if(ipden.ge.11.and.ipden.le.16)
     c   call stati_l(p1,p2,p3,p4,p5,ik,kk,w,bn,an,bnf,anf)
c260314
	goto 400
500	continue
400	continue

c	statistics of multiplicity distributions,
c	 spectator nucleons are excluded
	do kkk=1,ispmax
c	ik=ispkf(kkk)
c        if(iabs(ik).eq.2212.or.iabs(ik).eq.2112.or.iabs(ik).eq.3122.or.
c     &   iabs(ik).eq.3212.or.ik.eq.3222.or.ik.eq.3112)then
c	multiplicity is located at which interval
c	if(ik.eq.2212)then
	idf=bnf(kkk)/asd(5)+1
c	the 5-th distribution is particle multiplicity distribution
	if(idf.lt.1 .or. idf.gt.20)goto 405   ! 131204
	anf(idf,5,kkk)=anf(idf,5,kkk)+1./asd(5)
c       active the window 
405	do i=1,iflmax   ! 131204
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 404
        enddo
	idd=bn(kkk)/asd(5)+1
	if(idd.lt.1 .or. idd.gt.20)goto 404   ! 131204
	an(idd,5,kkk)=an(idd,5,kkk)+1./asd(5)
404	continue
c	endif
	enddo
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine stati_h(y,pt,eta,p5,ik,kk,ww,a,b,af,bf)   ! 260314
c	on line statistics for NA,AN,and AA collisions   ! 260314
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &  iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
	dimension a(20),b(40,5,20),c(5),af(20),bf(40,5,20),id(5) ! 070419
	amass=p5   ! 010600
        amass2=amass*amass
        pt2=pt*pt
        et=dsqrt(pt2+amass2)
        do 10000 i=1,iflmax
        goto (10,20,30,40,50) i
10      c(i)=y   
c150622 if(ifram.eq.1)c(i)=eta   
	goto 10000
20	c(i)=pt
	goto 10000
30	continue 
40	continue 
50	continue
10000	continue
c	calculate the abscissa one by one
40000	do 20000 i=1,isdmax
	goto (100,200,300,400,500) i
c	y is located in which interval?
100	ii=dabs(y)/asd(i)+1
	if(ifram.eq.1 .and. y.gt.0.)ii=ii+20   ! 311019
        if(ifram.eq.1 .and. y.lt.0.)ii=20-ii+1   ! 311019
c       note: 20 here should be change together with the dimension
c        40   ! 311019
	id(i)=ii
c100	id(i)=y/asd(i)+1
c	if(id(i).le.0)id(i)=1
	goto 20000
c	pt is located in which interval?
200	id(i)=pt/asd(i)+1
c	if(id(i).le.0)id(i)=1
	goto 20000
c	eta is located in which interval?
300	ii=dabs(eta)/asd(i)+1
	if(ifram.eq.1 .and. eta.gt.0.)ii=ii+20   ! 311019
        if(ifram.eq.1 .and. eta.lt.0.)ii=20-ii+1   ! 311019
c       note: 20 here should be change together with the dimension
c        40   ! 311019
	id(i)=ii
	goto 20000
c	et is counted in which eta interval ?
400	id(i)=id(i-1)
	if(kk.eq.3)id(i)=p5/0.1+1
c	statistcs of the mass dis. of rho0
        goto 20000
500	continue
20000	continue
c	make statistics of particle yield and desired distributions
	af(kk)=af(kk)+ww
	do i=1,isdmax
        ii=id(i)
	if(ii.lt.1 .or. ii.gt.40)goto 30000   ! 010218 070419
c	if(i.eq.2)then
c	bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)/pt
c	goto 30000
c	endif
        if(i.eq.1)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)
        if(i.eq.2)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)/pt
        if(i.eq.3)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)
c010600
	if(i.eq.4 .and. kk.eq.3)then
	bf(ii,i,kk)=bf(ii,i,kk)+ww/0.1
	goto 30000
	endif 
c010600
        if(i.eq.4)bf(ii,i,kk)=bf(ii,i,kk)+et/asd(i-1)
30000	enddo
c	if(iflmax.eq.0)return
c	put the filter to be effective
	do i=1,iflmax 
	if(c(i).lt.afl(kk,i,1) .or. c(i).gt.afl(kk,i,2))return
	enddo

c	make statistics for particle yield and desired distributions
	a(kk)=a(kk)+ww
	do i=1,isdmax
	ii=id(i)
	if(ii.lt.1 .or. ii.gt.40)goto 50000   ! 010218  070419
c	if(i.eq.2)then
c	b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)/pt
c	goto 50000
c	endif
        if(i.eq.1)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
        if(i.eq.2)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)/pt
        if(i.eq.3)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
c010600
	if(i.eq.4 .and. kk.eq.3)then
	b(ii,i,kk)=b(ii,i,kk)+ww/0.1
	goto 50000
	endif 
c010600
        if(i.eq.4)b(ii,i,kk)=b(ii,i,kk)+et/asd(i-1)
50000	enddo
	return
	end



c260314cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine stati_l(p1,p2,p3,p4,p5,ik,kk,ww,a,b,af,bf)   
c       on line statistics for case of lepton incidence
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &  iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
	common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph
     c	 ,vnlep
	dimension a(20),b(40,5,20),af(20),bf(40,5,20),c(5),id(5) ! 070419
c        vnu: \nu; fq2: Q^2=-q^2; w2l: W^2; yyl: y; zl: z; xb: x_B; pph: P_h
c	calculate kinematic variable relevant to the produced hadron
	pph=p1*p1+p2*p2+p3*p3
	pph=dmax1(pph,1.d-20)
	pph=dsqrt(pph)
	zln=pinch(4)*p4-pinch(1)*p1-pinch(2)*p2-pinch(3)*p3   ! numerator of z
	zld=pinch(5)*vnu   ! denominator of z
	zld=dmax1(zld,1.d-20)
	zl=zln/zld
c	write(9,*)'ik,kk,pph,zl=',ik,kk,pph,zl
        do 10000 i=1,iflmax   ! iflmax: total number of kinematic constrants
        goto (10,20,30,40,50) i
10      c(i)=fq2   ! -q^2
	goto 10000
20      c(i)=w2l   ! W^2
        goto 10000
30      c(i)=yyl   ! y
        goto 10000
40      c(i)=pph   ! p_h
	goto 10000
50      c(i)=zl   ! z
10000   continue
c	write(9,*)'c(i)=',(c(i),i=1,5)
c       calculate the abscissa one by one
40000   do 20000 i=1,isdmax
        goto (100,200,300,400,500) i
c       z is located in which interval?
100     id(i)=zl/asd(i)+1
        goto 20000
c       \nu is located in which interval?
200     id(i)=vnu/asd(i)+1
        goto 20000
c       -q^2 is located in which interval?
300     id(i)=fq2/asd(i)+1
        goto 20000
400     continue
        goto 20000
500     continue
20000   continue
c	write(9,*)'asd(i)=',(asd(i),i=1,5)
c	write(9,*)'id(i)=',(id(i),i=1,3)
c       make statistics of particle yield and desired distributions
        af(kk)=af(kk)+ww
        do i=1,isdmax
        ii=id(i)
        if(ii.lt.1 .or. ii.gt.40)goto 30000   ! 070419
        if(i.eq.1)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)
        if(i.eq.2)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)
        if(i.eq.3)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)
30000   enddo
c       put kinematic constraints to be effective
        do i=1,iflmax
        if(c(i).lt.afl(kk,i,1) .or. c(i).gt.afl(kk,i,2))return
        enddo
c       make statistics of particle yield and desired distributions
        a(kk)=a(kk)+ww
        do i=1,isdmax
        ii=id(i)
	if(ii.lt.1 .or. ii.gt.40)goto 50000   ! 070419
        if(i.eq.1)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
        if(i.eq.2)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
        if(i.eq.3)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
50000   enddo
        return
        end
