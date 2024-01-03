	subroutine coales(ijk,neve,nnout,nap,nat,nzp,nzt)   ! 300713
c     A simple coalescence model writen by Sa Ben-Hao on 04/06/2004
c       Its input messages are in 'parlist'
c       Its working block is 'parlist' ('pyjets' either in 'decayh')
c	  Its output message is in 'pyjets' (in 'sa1_h' either)
c	ijk: the run number
c	neve: total number of runs
c     nnout: a internal printing per nnout runs
c	nap and nzp: atomic and charge number of projectile
c     nat and nzt: atomic and charge number of target
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP

	PARAMETER (KSZJ=40000)
        parameter (mplis=40000)
	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c	Those variables in above four statements are only used here and   
c	in subroutine 'decayh','findb' and 'thephi'. 
c	PYDAT1,PYDAT2,PYDAT3 and PYJETS are the subroutines in PYTHIA
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
	common/sa6_t/ithroq_t,ithrob_t,ich_t,non6_t,throe_t(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  napp,natt,nzpp,nztt,pio   !
	dimension rc(3),b(3),pstr(3),rstr(3),numb(3),jk(20)
	dimension p0(4),pf(20,4),ppt(50,2),peo(5),ppp(20,5)
	dimension ppr(kszj,4),peoh(5),pnn(kszj,5),rr(3)
	dimension spglu(4),ppglu(4) !080512, the sum of E-p of unsplitted gluon

	rrp=1.16
	nn=0
      do i1=1,kszj
      do j1=1,5
        kn(i1,j1)=0.
        pn(i1,j1)=0.
        rn(i1,j1)=0.
	enddo
	enddo

	ithroq_t=0
      ithrob_t=0
      ich_t=0
      do i=1,4
        throe_t(i)=0.
      enddo

	nout=nnout
	imc=adj1(13)
      ibc=adj1(14)
	iphas=adj1(21)

c     The probability of gluon spliting into u,d & s quark
        adj132=adj1(32)
        prosum=1.+1.+adj132
        prod=1./prosum   ! 0.4286 originally
        pros=adj132/prosum   ! 0.1428 originally
        prods=prod+pros   ! 0.5714 originally

c	Conservation of net baryon.
	netba=0
	do i1=1,nbh
	kf=kbh(i1,2)
	kfab=iabs(kf)
	if(kf.gt.0.and.(kf.gt.1000 .and. kf.lt.10000))netba=netba+1
	if(kf.lt.0.and.(kfab.gt.1000 .and. kfab.lt.10000))netba=netba-1
	enddo
c1210214	if(nap.gt.1 .and. nat.gt.1)netba=nap+nat-netba
c120214	if(nap.eq.1.and.nat.eq.1.and.ipden.eq.0)netba=nzp+nzt-netba ! pp 300713
c120214	if(nap.eq.1.and.nat.eq.1.and.ipden.eq.2.and.itden.eq.0)
c120214	c   netba=nzt-netba   ! ep 300713
	if(ipden.eq.0.and.itden.eq.0)netba=nzp+nzt-netba ! pp 300713 120214
	if(ipden.gt.2.and.itden.eq.0)
     c	 netba=nzt-netba   ! for lepton+p 300713 120214
c120214	for e+e-, need not to recalculate netba

c160110 if(ijk.eq.1)then
c       napt=nap+nat
c       nzpt=nzp+nzt
c       sbaryi=dfloat(napt)
c       schgei=3.*dfloat(nzpt)
c160110 endif

c     Throw away t quark (antiquark) for the moment.
888     continue
        do i1=1,iprl
        kf=idp(i1)
        if(iabs(kf).eq.6)then
        if(kf.eq.6)ithroq_t=ithroq_t+1
        if(kf.eq.-6)ithrob_t=ithrob_t+1
        ich_t=ich_t+pychge(kf)
        do i2=1,4
        throe_t(i2)=throe_t(i2)+pp(i2,i1)
        enddo
c       Move parton list one step downward since i1+1
        call updad(iprl,i1+1,1)
        iprl=iprl-1
        goto 888
        endif
        enddo

c     make the partons in order of g, qba and q.
      iii=0
      jjj=0
        do ii=1,3
c	1: refers to g, 2: qba, 3: q
        kf=21
        do j=iii+1,iprl
        call ord_c(jjj,j,kf,ii)
        enddo
        iii=jjj
        numb(ii)=jjj
c       numb(1),(2) and (3): the order # of last g,qba & q
        enddo
        n1=numb(1)
        n2=numb(2)
        n3=numb(3)

c100	split forcibly gulon (after 'parcas') into qqba pair
c	amu=pymass(2)
c	amd=pymass(1)
c	ams=pymass(3)
	amd=0.0099D0  
	amu=0.0056D0 
	ams=0.199D0
c	write(9,*)'amu,amd,ams=',amu,amd,ams   ! sa
	amuu=2*amu
	amdd=2*amd
	amss=2*ams
	deles=0.
	do i1=1,4  ! 080512
	spglu(i1)=0.0
	enddo  ! 080512
100	continue
	if(n1.le.0)goto 102
	do ii=1,n1   ! sa

c140604	if(n1.eq.0 .or. ii.gt.n1)goto 102
	dele=0.
	eg=pp(4,ii)
	if(eg.lt.amuu)goto 200
	if(eg.ge.amuu .and. eg.lt.amdd)then
	kf=2   ! u
	amq=amu
	endif
	if(eg.ge.amdd .and. eg.lt.amss)then
	kf=2   ! u
	amq=amu
	if(pyr(1).gt.0.5)then  
	kf=1   ! d
	amq=amd
	endif
	endif
	if(eg.gt.amss)then
	rand=pyr(1)
	kf=3   ! s
	amq=ams
	if(rand.gt.pros .and. rand.le.prods)then
	kf=1   ! d
	amq=amd
	endif
	if(rand.gt.prods)then
	kf=2   ! u
	amq=amu
	endif
	endif
	do j=1,4
	p0(j)=pp(j,ii)
	enddo

c	Fill the q & qba splited from g into parton list
c	arrange qba after n2
c	Move parton list one step forward since n2+1 upto n3
	call updaf(n3,n2)
	n2=n2+1 
	n3=n3+1
c	arrange q after n3
	n3=n3+1

c	Splited qba takes the four coordinate of g, splited q is arranged around 
c	 g within 0.5 fm randumly in each one of the three coordinates and 
c	 has same 4-th coordinate as g
	do i=1,4
	rp(i,n2)=rp(i,ii)
	enddo
	do i=1,3
	rr(i)=pyr(1)*0.5
	rp(i,n3)=rp(i,ii)+rr(i)
        if(pyr(1).gt.0.5)rp(i,n3)=rp(i,ii)-rr(i)
	enddo
	rp(4,n3)=rp(4,ii)

c	Give momentum to q & qba
	do j1=1,20
	do j2=1,5
	ppp(j1,j2)=0.
	enddo
	enddo
 
	decsuc=1
        call decmom(p0,ppp,amq,amq,decsuc)
        do j1=1,4
        pnn(1,j1)=ppp(1,j1)
        pnn(2,j1)=ppp(2,j1)
        enddo
        if(decsuc.eq.0)then   !
        do i1=1,3
	pii=pyr(1)*p0(i1)
	pnn(1,i1)=pii
	pnn(2,i1)=p0(i1)-pii
	enddo
	pn11=pnn(1,1)
        pn12=pnn(1,2)
        pn13=pnn(1,3)
	pnn(1,4)=dsqrt(amq*amq+pn11*pn11+pn12*pn12+pn13*pn13)
        pn21=pnn(2,1)
        pn22=pnn(2,2)
        pn23=pnn(2,3)
        pnn(2,4)=dsqrt(amq*amq+pn21*pn21+pn22*pn22+pn23*pn23)
	pnn(1,5)=amq
	pnn(2,5)=amq
	endif   ! 1

c	'pnn' is a internal array 
	do j=1,4
	pp(j,n2)=pnn(1,j)
	pp(j,n3)=pnn(2,j)
	enddo
	dele=eg-pnn(1,4)-pnn(2,4)
	deles=deles+dele

c	Give other characters to q & qba 
	idp(n2)=-kf
	idp(n3)=kf
	rmp(n2)=amq
	rmp(n3)=amq
	goto 400   ! 140604 100 originally

200	continue

c	Treating the case of gluon can not splite into qqba 
	if(ii.eq.n1)goto 300

c	Transfer the four momentum of that gluon to anyone of other gluons 
c	randomly.
	ran=pyr(1)
	iran=ran*(n1-ii)+1
	do i1=1,4
	pp(i1,iran)=pp(i1,iran)+pp(i1,ii)
	enddo
	goto 400
300	continue

c	Transfer the four momentum of that gluon to anyone of qba randomly.
c	  ran=pyr(1) ! 080512
c	  iran=ran*(n2-n1)+1
c        do i1=1,4
c	    pp(i1,iran)=pp(i1,iran)+pp(i1,ii)
c	  enddo ! 080512

c	record the four momentum of that gluon, and evenly share that  
c	 record with the quarks and antiquarks later. yan, ! 080512 sa
	  do i1=1,4
		spglu(i1)=spglu(i1)+pp(i1,ii)
	  enddo   ! 080512
400	  continue
	call updad(n3,ii+1,1)

c	Move parton list one step downward since ii+1
	n1=n1-1
	n2=n2-1
	n3=n3-1

c	Share the dele
c	Share the dele, abandon, yan, 080512 
c	  if(n3.gt.0)then
c		dele=dele/float(n3)
c		do i1=1,n3
c		  pp(4,i1)=pp(4,i1)+dele
c	      if(dele.lt.0.)then
c		  if(pp(4,i1).lt.0.)pp(4,i1)=pp(4,i1)-dele
c		    pabs=abs(pp(3,i1))
c		    if(pabs.ge.pp(4,i1))pp(4,i1)=pp(4,i1)-dele
c		  endif
c	    enddo
c	  endif  080512 
	goto 100

	enddo   ! sa	
102	continue
c	share the dele, yan, 080512 
c	write(99,*)'g n1,n3,dele,deles,e=',n1,n3,dele,deles
c     c	(pp(4,i1),i1=1,n3)   ! sa
	if(n3.gt.0)then
	  do i1=1,3
		ppglu(i1)=spglu(i1)/float(n3)
	  enddo
	  ppglu(4)=(ppglu(4)+deles)/float(n3)
	  do i1=1,n3
		do i2=1,4
		  pp(i2,i1)=pp(i2,i1)+ppglu(i2)
	    enddo
	  enddo
	endif ! 080512 
101	format(4(1x,f10.4))
600   format(20(1x,i3))
c     Split forcibly gluon into qqba pair, finished.   ! 080512 sa

	iprl=n3
	igens=0
	adj12=adj1(12)
	adj16=adj1(16)
	adj17=adj1(17)
	adj17=max(4.0,adj17) ! 070612, yan
	if(adj12.eq.2)goto 900   
c		write(9,*)'adj16=',adj16   ! 080512 
c     Parton production according to Field-Feynman model
	iprloo=1   ! 080512 sa
700	do 800 i1=iprloo,iprl   ! 080512 sa   
        kf0=idp(i1)
	ee=pp(4,i1)
        iflav=1
        if(kf0.lt.0)iflav=-1
c       iflav = 1 : if source parton is quark
c             =-1 : if source parton is antiquark
c     if(ee .gt. adj17) call ffm(i1,kf0,igen,iflav,n1,n2,n3) ! 080512
	rand=pyr(1)  ! 080512
c	rand=rlu(1)
      if(ee.gt.adj17.and.rand.lt.adj16) then  ! 080512
	  call ffm(i1,kf0,igen,iflav,n1,n2,n3) 
	endif  ! 080512
c       igen : times of energetic quark can excite qqbar pair from     
c        vacuum (in 'ffm'), which is controled by four momenta of energetic quark ! 080512  

800     continue
	iprlo=iprl   ! 080512 sa
	iprl=n3
	igens=igens+1 
c080512	igens: repeating times of considering deexcitation of energetic 
c080512	 parton over parton list   
c	if(igens.gt.adj16)goto 900  ! 080512
	iadj16=int(adj16) ! 080512
	if(igens.ge.iadj16)goto 900 ! 080512
c080512 sa do i1=1,iprl
c	ee=pp(4,i1)
c	if(ee .gt. adj17)goto 700
c080512 sa	enddo
	iprloo=iprlo+1   ! 080512 sa
	goto 700   ! 080512 sa
900	continue
c     Parton production according to Field-Feynman model, finished.

	if(adj12.eq.2)goto 1000 ! no need in order of qba and q. 

c     Make the partons in order of qba and q
        iii=0
        jjj=0
        do ii=2,3   ! 1,3
c       1: refers to g, 2: qba, 3: q
        kf=21
        do j=iii+1,iprl
        call ord_c(jjj,j,kf,ii)
        enddo
        iii=jjj
        numb(ii)=jjj
c       numb(1),(2) and (3): the order # of last g,qba & q
        enddo
c271004	n1=numb(1)
        n2=numb(2)
        n3=numb(3)

1000	continue
        iqba=n2

c     Order the qba according to energy from the maximal to minimal.
	call eord(n1+1,n2)

c     Order the q according to energy from the maximal to minimal.
        call eord(n2+1,n3)

c	Parton coalescence 
c160110 if(ijk.eq.1)call tabhb
        call tabhb   ! 160110
c	read the table of hadron (pseudoscalar 0- & vector 1- only) &
c	 primary baryon (spin 1/2 octet & 3/2 decuplet only)
	
	call coal(n3,iqba,ijk,rrp,iphas,netba)
c	n3: total number of partons (qba and q) 
c	iqba: total # of qba (qba is ordered before q) 
c	ithroq : the total # of quarks thrown away
c	ithrob : the total # of antiquarks thrown away
c	throe : total 4-momentum of the partons thrown away
c	ich : total charge of the partons thrown away

	iqba=ithrob   ! 110905 not active originally
	n3=ithroq+ithrob   ! 110905 not active originally
	if(iphas.eq.1 .and. n3.ge.2)then
	call coal(n3,iqba,ijk,rrp,0,0)   ! 110905
	endif
 
	ichth=ich   ! 092600

c	Transfer the data form 'sbh' to 'sa1_h'.
        if(nbh.ge.1)then
        do l=1,nbh
        l1=nn+l
        do m=1,5
        kn(l1,m)=kbh(l,m)
        pn(l1,m)=pbh(l,m)
        rn(l1,m)=vbh(l,m)
        enddo
        enddo
        nn=nn+nbh
        endif

c	Transfer the data form 'sa1_h' to 'pyjets'.
	n=nn
	do j1=1,nn
        do j2=1,5
        k(j1,j2)=kn(j1,j2)
        p(j1,j2)=pn(j1,j2)
        v(j1,j2)=rn(j1,j2)
        enddo
        enddo

c     Decay of unstable hadrons
	call decayh(rrp)

	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ord_c(ipi,j,kf,ii)
c	Make order for particle j
c	ipi: particle j should be ordered after ipi
c	kf: a control variable
c	ii: a control variable 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      parameter (mplis=40000)
      common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	dimension rr(4),p1(4)
	ik=idp(j)
      if(ii.eq.1 .and. ik.eq.kf)goto 100   ! g
	if(ii.eq.2 .and. ik*kf.lt.0)goto 100   ! qba
	if(ii.eq.3 .and. ik*kf.gt.0)goto 100   ! q
	goto 200
100	ipi=ipi+1
	idpp=idp(ipi)
	rmpp=rmp(ipi)
	do jj=1,4
	  p1(jj)=pp(jj,ipi)
	  rr(jj)=rp(jj,ipi)
	enddo
	idp(ipi)=idp(j)
	rmp(ipi)=rmp(j)
	do jj=1,4
	  pp(jj,ipi)=pp(jj,j)
	  rp(jj,ipi)=rp(jj,j)
	enddo
	idp(j)=idpp
	rmp(j)=rmpp
	do jj=1,4
	  pp(jj,j)=p1(jj)
	  rp(jj,j)=rr(jj)
	enddo
200	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ord_th(ipi,j,kf,ii)
c	Make order for particle j
c	ipi: particle j should be ordered after ipi
c	kf: a control variable
c	ii: a control variable 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (mplis=40000)
	common/throqb/iprlth,nonth,rpth(4,mplis),ppth(4,mplis),
     c	 idpth(mplis),rmpth(mplis)   ! 110905
	dimension rr(4),p1(4)
	ik=idpth(j)
        if(ii.eq.1 .and. ik.eq.kf)goto 100   ! g
	if(ii.eq.2 .and. ik*kf.lt.0)goto 100   ! qba
	if(ii.eq.3 .and. ik*kf.gt.0)goto 100   ! q
	goto 200
100	ipi=ipi+1
	idpp=idpth(ipi)
	rmpp=rmpth(ipi)
	do jj=1,4
	p1(jj)=ppth(jj,ipi)
	rr(jj)=rpth(jj,ipi)
	enddo
	idpth(ipi)=idpth(j)
	rmpth(ipi)=rmpth(j)
	do jj=1,4
	ppth(jj,ipi)=ppth(jj,j)
	rpth(jj,ipi)=rpth(jj,j)
	enddo
	idpth(j)=idpp
	rmpth(j)=rmpp
	do jj=1,4
	ppth(jj,j)=p1(jj)
	rpth(jj,j)=rr(jj)
	enddo
200	return
	end        



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pes(pei,il,ih,peo)
c	sum up momentum and energy  
c	pei : two dimension array of input momentum and energy
c	il and ih : lower and higher limits of summation
c	peo : one dimension array of output momentum,energy & sqrt(s)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	PARAMETER (KSZJ=40000)
        parameter (mplis=40000)
        dimension pei(kszj,4),peo(5)
	do i=1,5
	peo(i)=0.
	enddo
	do i=il,ih
	do j=1,4
	peo(j)=peo(j)+pei(i,j)
	enddo
	enddo
	peo(5)=peo(4)*peo(4)
	do i=1,3
	peo(5)=peo(5)-peo(i)*peo(i)
	enddo
	peo5=peo(5)
	if(peo5.lt.0.)then
	peo5=dabs(peo5)
	endif

100	format('            px           py          pz         e      
     c	 sqrt(s)')
200	format(4x,5(1x,f9.3))
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine tabhb
c	The table of primary hadron (pseudoscalar 0- & vector 1- only) &
c	 primary baryon (spin 1/2 octet & 3/2 decuplet only)
c	kqh,kfh,proh: the quark composition,the flavor code,the
c	 probability of hadron
c	kqb,kfb,prob: the quark composition,the flavor code,the
c	 probability of baryon
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
	data (kqh(i,1),i=1,80)/2*2,2*1,2*3,3*2,3*1,2*3,4,1,4,2,4,3,4,
     c	 2,5,1,5,5,54*0/	
	data (kqh(i,2),i=1,80)/-1,-3,-2,-3,-2,-1,-2,-2,-2,-1,-1,-1,-3,
     c	     -3,-1,-4,-2,-4,-3,-4,-4,-5,-2,-5,-1,-5,54*0/
	data (kfh(i,1),i=1,80)/211,321,-211,311,-321,-311,111,221,331,
     c	     111,221,331,221,331,411,-411,421,-421,431,-431,441,521,
     c	     -521,511,-511,553,54*0/
	data (kfh(i,2),i=1,80)/213,323,-213,313,-323,-313,113,223,0,
     c	     113,223,0,333,0,413,-413,423,-423,433,-433,443,2*0,513,
     c	     -513,0,54*0/
	data (proh(i,1),i=1,80)/6*1.,0.5,2*0.25,0.5,2*0.25,2*0.5,12*1,
     c	     54*0./
	data (proh(i,2),i=1,80)/6*1.,2*0.5,0.,2*0.5,0.,1.,0.,7*1,
     c	     2*0,2*1,0,54*0./
	data (kqb(i,1),i=1,80)/5*2,3*1,3,4*2,1,2,1,3,2,62*0/
	data (kqb(i,2),i=1,80)/3*2,1,3,2*1,2*3,3*1,2,1,3*3,1,62*0/
	data (kqb(i,3),i=1,80)/2,1,3,1,3,1,5*3,6*4,5,62*0/
	data (kfb(i,1),i=1,80)/0,2212,3222,2112,3322,0,3112,3312,0,
     c	     3122,3212,4122,4222,4112,4232,4132,4332,5122,62*0/
	data (kfb(i,2),i=1,80)/2224,2214,3224,2114,3324,1114,3114,3314,
     c	     3334,0,3214,4212,4222,0,4232,4132,2*0,62*0/
	data (prob(i,1),i=1,80)/0.,4*1,0.,2*1.,0.,2*0.5,7*1.,62*0./
	data (prob(i,2),i=1,80)/9*1.,0.,1.,2*1.,0.,2*1.,2*0.,62*0./
	do i1=1,imc
	kf1=kfh(i1,1)
	kf2=kfh(i1,2)
	amash(i1,1)=pymass(kf1)
	amash(i1,2)=pymass(kf2)
	enddo
	do i1=1,ibc
	kf1=kfb(i1,1)
	kf2=kfb(i1,2)
	amasb(i1,1)=pymass(kf1)
	amasb(i1,2)=pymass(kf2)
	enddo

	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updaf(n2,n1)
c     Move parton list one step forward since n1+1
c	n2 : the last order # of parton list
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (mplis=40000)
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	do j=n2,n1+1,-1
	jj=j+1
	idp(jj)=idp(j)
	rmp(jj)=rmp(j)
	do j1=1,4
	pp(j1,jj)=pp(j1,j)
	rp(j1,jj)=rp(j1,j)
	enddo
	enddo
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updad(j2,jc,i)
c	Move the parton list i steps downward from jc till j2.
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (mplis=40000)
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	do j=jc,j2
	idp(j-i)=idp(j)
	rmp(j-i)=rmp(j)
	do jj=1,4
	pp(jj,j-i)=pp(jj,j)
	rp(jj,j-i)=rp(jj,j)
	enddo
	enddo
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine eord(ni,nc)
c	Order particle set (ni to nc) according to energy
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (mplis=40000)
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	dimension rr(4),p1(4)
	do 100 i1=ni,nc
	ii=i1
	if(ii.eq.nc)goto 100
	j=ii
	alar=pp(4,ii)
	do 200 i2=ii+1,nc
	ee=pp(4,i2)
	if(ee.ge.alar)then
	j=i2
	alar=ee
	endif
200	continue
	idpp=idp(ii)
	rmpp=rmp(ii)
	do jj=1,4
	p1(jj)=pp(jj,ii)
	rr(jj)=rp(jj,ii)
	enddo
	idp(ii)=idp(j)
	rmp(ii)=rmp(j)
	do jj=1,4
	pp(jj,ii)=pp(jj,j)
	rp(jj,ii)=rp(jj,j)
	enddo
	idp(j)=idpp
	rmp(j)=rmpp
	do jj=1,4
	pp(jj,j)=p1(jj)
	rp(jj,j)=rr(jj)
	enddo
100	continue
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coal(n1,nqb,ijk,rrp,iphas,netba)
c	Parton coalescence (hadronization)
c	n1 : total # of partons (q & qba only)
c	nqb : total # of qba (qba is ordered before q)
c	ijk : the run number
c	iphas=1: with phase space adjudgment
c	     =0: without phase space adjudgment
c	netba: number of baryons to be first generated keeping
c	 baryon conservation 
c	ithroq : total # of quarks thrown away
c	ithrob : total # of antiquarks thrown away
c	throe : total momentum and energy of the partons thrown away
c	ich : total charge of the partons thrown away
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=40000,mplis=40000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa24/adj1(40),nnstop,non24,zstop
	common/throqb/iprlth,nonth,rpth(4,mplis),ppth(4,mplis),
     c	 idpth(mplis),rmpth(mplis)   ! 110905
	dimension pc(4),rc(3),iar(3),rcp(3)
	dimension psu(3),peo(5),pnn(kszj,5)
	dimension numb(3)   ! 110905
	ithroq=0
	ithrob=0
	ich=0
	do i=1,4
	throe(i)=0.
	enddo
c110905
	iprlth=0
        do j=1,mplis
        idpth(j)=0
        rmpth(j)=0.
        do i=1,4
        rpth(i,j)=0.
        ppth(i,j)=0.
        enddo
        enddo
c110905

c	Coalescence
	ibarp=0
	ibarm=0
	imes=0
c2	sumch=0.

c	Generate first 'netba' baryons (if 'netba'>0) or -'netba' antibaryons 
c	 (if 'netba'<0) keeping baryon conservation 
	if(netba.gt.0)call barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,1)
	if(netba.lt.0)then
302     continue
	do ii1=1,nqb
	call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,1)
	if(isu.eq.1 .and. ibarm.lt.-netba)goto 303
        if(isu.eq.1 .and. ibarm.eq.-netba)goto 304
        if(isu.eq.0) stop 6666
	enddo
303     goto 302
304     continue
	endif
	
c     First compose antiquark into hadron one by one.
c     A antiquark can be a component of antibaryon or meson.
	jj=0
300	continue
	jj=jj+1
	if(nqb.eq.0)goto 100   ! for baryon composing
	if(nqb.lt.3 .and. n1.eq.nqb)goto 301 ! no q, throw away nqb
c	do ii1=1,nqb
	if(nqb.eq.1)ii1=1
	if(nqb.gt.1)ii1=nqb*pyr(1)+1
        kf=idp(ii1)
        relpr=adj1(31)
        relpr=relpr/(1.+relpr)
c       probability of antibaryon
        if(kf.eq.-3)then
        relpr=adj1(31)*adj1(33)
        relpr=relpr/(1.+relpr)
c       probability of strange antibaryon
        endif
	rand=pyr(1)

        if(rand.le.relpr .and. nqb.ge.3)then   !
        call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,2)

        if(isu.eq.1 .and. nba.eq.1)then   ! one antibryon produced
	goto 300
	endif   
        if(isu.eq.0)then   !! ii1 can not -> antibryon
        call mespro(n1,nqb,ii1,iphas,nme,imes,rrp,isu)
        if(isu.eq.1 .and. nme.eq.1)then   ! one hadron produced
	goto 300   
	endif
        if(isu.eq.0)then   !!! ii1 can not -> meson either, throw away
	iprlth=iprlth+1   ! 110905
        ithrob=ithrob+1
        kk=idp(ii1)
        ich=ich+pychge(kk)
        idpth(iprlth)=kk
        rmpth(iprlth)=rmp(ii1)
        do i2=1,4
        throe(i2)=throe(i2)+pp(i2,ii1)
        ppth(i2,iprlth)=pp(i2,ii1)
        rpth(i2,iprlth)=rp(i2,ii1)
        enddo
c     Move parton list one steps downward since ii1+1
        call updad(n1,ii1+1,1)   ! changed
        n1=n1-1
        nqb=nqb-1
        goto 300
        endif   !!!
        endif   !!
        elseif(rand.gt.relpr .or. nqb.lt.3)then   !
        call mespro(n1,nqb,ii1,iphas,nme,imes,rrp,isu)
        if(isu.eq.1 .and. nme.eq.1)then   ! one hadron produced
	goto 300   
	endif
        if(isu.eq.0)then   !! ii1 can not -> meson
        call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,2)
        if(isu.eq.1 .and. nba.eq.1)then   ! one antibryon produced
	goto 300   
	endif
        if(isu.eq.0)then   !!! ii1 can not -> antibryon either, throw away
	iprlth=iprlth+1   ! 110905
        ithrob=ithrob+1
        kk=idp(ii1)
        ich=ich+pychge(kk)
        idpth(iprlth)=kk
        rmpth(iprlth)=rmp(ii1)
        do i2=1,4
        throe(i2)=throe(i2)+pp(i2,ii1)
        ppth(i2,iprlth)=pp(i2,ii1)
        rpth(i2,iprlth)=rp(i2,ii1)
        enddo
c     Move parton list one steps downward since ii1+1
        call updad(n1,ii1+1,1)   ! changed
        n1=n1-1
	nqb=nqb-1
        goto 300
        endif   !!!
        endif   !!
	else   !
        endif   !
c	enddo

c	Throw away those qba remained
301	continue
	if(nqb.eq.0)goto 100
        do i1=1,nqb
	iprlth=iprlth+1   ! 110905
        ithrob=ithrob+1
        kk=idp(i1)
        ich=ich+pychge(kk)
        idpth(iprlth)=kk
        rmpth(iprlth)=rmp(i1)
        do i2=1,4
        throe(i2)=throe(i2)+pp(i2,i1)
        ppth(i2,iprlth)=pp(i2,i1)
        rpth(i2,iprlth)=rp(i2,i1)
        enddo
        enddo

c     Move parton list nqb steps downward since nqb+1
	call updad(n1,nqb+1,nqb)   ! changed
	n1=n1-nqb
	nqb=0

100	continue   ! for baryon (composed of quark only) production

600	format(20(1x,i3))
400     if(n1-nqb.eq.0)goto 505   ! get out of 'coal'
        if(n1-nqb.lt.3)goto 500   ! throw away those quarks

c     Proceed for baryon production
	call barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,0)
	if(n1-nqb.eq.0)goto 505
500     continue   

c     Throw away those q remained
503     continue
	if(n1-nqb.eq.0)goto 505
	ithroq=ithroq+n1-nqb
	do i1=nqb+1,n1
	iprlth=iprlth+i1   ! 110905
	j=iprlth
	kk=idp(i1)
	ich=ich+pychge(kk)
	idpth(j)=idp(i1)
	rmpth(j)=rmp(i1)
	do i2=1,4
	throe(i2)=throe(i2)+pp(i2,i1)
	ppth(i2,j)=pp(i2,i1)
	rpth(i2,j)=rp(i2,i1)
	enddo
	enddo
	n1=0   ! 110905

505	continue

c	Reconstruct parton list 'parlist'
	if(iphas.eq.1)then
c110905
c     Make the partons in order of qba and q
        iii=0
        jjj=0
        do ii=2,3   ! 1,3
c       1: refers to g, 2: qba, 3: q
        kf=21
        do j=iii+1,iprlth
        call ord_th(jjj,j,kf,ii)
        enddo
        iii=jjj
        numb(ii)=jjj
c       numb(1),(2) and (3): the order # of last g,qba & q
        enddo
c271004 n1=numb(1)
        n2=numb(2)
        n3=numb(3)
c110905
	nqb=ithrob
	n1=nqb+ithroq
	do i1=1,n1
	idp(i1)=idpth(i1)
	rmp(i1)=rmpth(i1)
	do j1=1,4
	pp(j1,i1)=ppth(j1,i1)
        rp(j1,i1)=rpth(j1,i1)
        enddo
        enddo
        do i1=n1+1,mplis
	idp(i1)=0
	rmp(i1)=0.
        do j1=1,4
        pp(j1,i1)=0.
        rp(j1,i1)=0.
        enddo
        enddo
	endif
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine findm(kf1,kf2,cm,kfii,amasi,isucc,iflav)
c   Find out the primary meson from mesonic table according to kf1 & kf2
c	cm : invariant mass of kf1 & kf2
c	kfii : flavor code of the primary meson
c	amasi : mass of the primary meson
c	isucc = 1 : success
c           = 0 : fail
c	iflav = 1 : kf1>0,do not need to permute kf1 & kf2
c	     = -1 : kf1<0,need to permute kf1 & kf2
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
	dimension ikf(2)

	if(iflav.eq.1)then   ! 1
	  if1=kf1
	  if2=kf2
	  do 500 i4=1,imc
          kfi=kqh(i4,1)
          kfj=kqh(i4,2)
		amas1=amash(i4,1)
		amas1=abs(cm-amas1)
		amas2=amash(i4,2)
		amas2=abs(cm-amas2)
        if(kfi.eq.if1 .and. kfj.eq.if2)then   ! 2 success
		if(proh(i4,1).eq.0 .and. proh(i4,2).eq.0.)goto 500
		if(proh(i4,1).eq.0 .and. proh(i4,2).ne.0.)goto 506   ! vector
		if((proh(i4,1).ne.0.and.proh(i4,2).ne.0.).and.amas2.le.amas1)
     c	goto 506   ! vector

c         Proceed for pseudoscalar
          kfii=kfh(i4,1)
          amasi=amash(i4,1)
          proi=proh(i4,1)
          ran1=pyr(1)
          if(ran1.gt.proi)goto 500
          goto 504   ! success

506       kfii=kfh(i4,2)   ! vector
          amasi=amash(i4,2)
          proi=proh(i4,2)
          ran1=pyr(1)
          if(ran1.gt.proi)goto 500
	    goto 504   ! success

        endif   ! 2
500	  continue
        isucc=0   ! fail
	  return
	endif   ! 1

	ikf(1)=kf1
	ikf(2)=kf2
c	Two body permutation = arrangement (2,2)
	do 501 i1=1,2
	  if1=ikf(i1)
	  do 502 i2=1,2
		if(i2.eq.i1)goto 502
		if2=ikf(i2) 
		do 503 i4=1,imc
        kfi=kqh(i4,1)
        kfj=kqh(i4,2)
        amas1=amash(i4,1)
        amas1=abs(cm-amas1)
        amas2=amash(i4,2)
        amas2=abs(cm-amas2)
        if(kfi.eq.if1 .and. kfj.eq.if2)then   ! success
          if(proh(i4,1).eq.0 .and. proh(i4,2).eq.0.)goto 503
          if(proh(i4,1).eq.0 .and. proh(i4,2).ne.0.)goto 505   ! vector
          if((proh(i4,1).ne.0.and.proh(i4,2).ne.0.).and.amas2.le.amas1)
     c    goto 505   ! vector

c         Proceed for pseudoscalar
          kfii=kfh(i4,1)
          amasi=amash(i4,1)
          proi=proh(i4,1)
          ran1=pyr(1)
          if(ran1.gt.proi)goto 503
          goto 504   ! success

505       kfii=kfh(i4,2)   ! vector
          amasi=amash(i4,2)
          proi=proh(i4,2)
          ran1=pyr(1)
          if(ran1.gt.proi)goto 503
	    goto 504

        endif
503       continue
502	  continue
501	continue
      isucc=0   ! fail
	return
504   isucc=1   ! success
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine findb(kf0,kf1,kf2,cm,kfii,amasi,isucc,iflav)
c     Find out the primary baryon (antibaryon) from baryonic table  
c	 according to kf0,kf1,and kf2,these flavor codes are all > 0
c	cm: invariant mass of kf0,kf1 & kf2
c	kfii : flavor code of the primary baryon
c	amasi : mass of the primary baryon
c	isucc = 1 : success
c       isucc = 0 : fail
c	iflav = 1 : if composing parton is quark
c             =-1 : if composing parton is antiquark
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
	dimension ikf(3)
	ikf(1)=kf0
	ikf(2)=kf1
	ikf(3)=kf2

c	Three body permutation = arrangement (3,3)
	do 104 i1=1,3
	if1=ikf(i1)
	do 105 i2=1,3
	if(i2.eq.i1)goto 105
	if2=ikf(i2)
	do 106 i3=1,3
	if(i3.eq.i2)goto 106
	if(i3.eq.i1)goto 106
	if3=ikf(i3)
	do 107 i4=1,ibc
	kfi=kqb(i4,1)
	kfj=kqb(i4,2)
        kfk=kqb(i4,3)
        amas1=amasb(i4,1)
        amas1=dabs(cm-amas1)
        amas2=amasb(i4,2)
        amas2=dabs(cm-amas2)
        if(kfi.eq.if1 .and. kfj.eq.if2 .and. kfk.eq.if3)then ! success
	if(prob(i4,1).eq.0.and.prob(i4,2).eq.0)goto 107
        if(prob(i4,1).eq.0 .and. prob(i4,2).ne.0.)goto 108   ! 3/2
        if((prob(i4,1).ne.0.and.prob(i4,2).ne.0.).and.amas2.le.amas1)
     c   goto 108   ! 3/2 
c	Goto 108, for spin 3/2 decuplet.
c	Proceed for spin 1/2 octet.
        kfii=kfb(i4,1)
        amasi=amasb(i4,1)
	if(iflav.eq.-1)then
	kfii=-kfb(i4,1)
	endif
        proi=prob(i4,1)
        ran1=pyr(1)
        if(ran1.gt.proi)goto 107
	goto 109   ! success
108     kfii=kfb(i4,2)   ! spin 3/2 decuplet
        amasi=amasb(i4,2)
	if(iflav.eq.-1)then
	kfii=-kfb(i4,2)
	endif
        proi=prob(i4,2)
        ran1=pyr(1)
        if(ran1.gt.proi)goto 107
	goto 109   ! success
	endif
107	continue
106	continue
105	continue
104	continue
	isucc=0   ! fail
	return
109	isucc=1   ! success

	return   
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,iway)
c     To compose baryon
c     n1 : total # of partons (q & qba)
c     nqb : total # of qba (qba is ordered before q)  
c     iphas=1: with phase space adjudgment
c          =0: without phase space adjudgment
c	ibarp: statistic number of baryon
c	netba: number of baryons keeping baryon conservation 
c	iway: a switch
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000,mplis=40000)
      common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
      common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
      common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
      common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
      common/sa24/adj1(40),nnstop,non24,zstop
	dimension rcp(3)

	dpmax=adj1(27)
      drmax=adj1(28)	
	nba=0
400	continue

	if(n1-nqb.lt.3)return   ! number of quarks less than three
      do 404 i1=nqb+1,n1-2
	  kf1=idp(i1)
        do 405 i2=i1+1,n1-1
          kf2=idp(i2)
          do 406 i3=i2+1,n1
            kf3=idp(i3)

	  sume=pp(4,i1)+pp(4,i2)+pp(4,i3)
	  sump1=pp(1,i1)+pp(1,i2)+pp(1,i3)
	  sump2=pp(2,i1)+pp(2,i2)+pp(2,i3)
        sump3=pp(3,i1)+pp(3,i2)+pp(3,i3)
	  cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
	  if(cm.gt.1.e30)cm=1.e30
	  if(cm.le.0.)goto 406   ! fail 071204
        cm=sqrt(cm)

c     Find out the primary baryon from hadron table due to kf1,kf2 & kf3
        call findb(kf1,kf2,kf3,cm,kfii,amasi,isucc,1)

        if(isucc.eq.0)goto 406   ! fail, and keep on cycle, try again.

c	Proceed for success   
c	Phase space adjudgment
        if(iphas.eq.1)then
	    call phas(i1,i2,i3,isucc,3)
          if(isucc.eq.0)goto 406   ! fail
        endif

c	  Proceed for success
        ibarp=ibarp+1
	  nba=nba+1

c       Give proper variables to the primary baryon
        nnol=nn
        nn=nn+1
	  kn(nn,1)=1
        kn(nn,2)=kfii
        kn(nn,3)=0
        pn(nn,5)=amasi
	  pn(nn,1)=sump1
	  pn(nn,2)=sump2
	  pn(nn,3)=sump3
	  pnnm=sump1*sump1+sump2*sump2+sump3*sump3
	  pnnmm=amasi*amasi+pnnm
	  if(pnnmm.gt.1.e30)pnnmm=1.e30
        if(pnnmm.le.0.)pnnmm=1.e-30
	  pnnn=sqrt(pnnmm)
	  pn(nn,4)=pnnn
	  dele=sume-pnnn

c     Produced hadron is arranged among constituent partons randomly.
	  pyrx=pyr(1)
	  pyry=pyr(1)
	  pyrz=pyr(1)
	  rn(nn,1)=pyrx*rp(1,i1)+pyry*rp(1,i2)+pyrz*rp(1,i3)
	  rn(nn,2)=pyrx*rp(2,i1)+pyry*rp(2,i2)+pyrz*rp(2,i3)
	  rn(nn,3)=pyrx*rp(3,i1)+pyry*rp(3,i2)+pyrz*rp(3,i3)

c     Move parton list one step downward from i3+1 to n1
411	  call updad(n1,i3+1,1)
	  n1=n1-1

c     Move parton list one step downward from i2+1 to n1
	  call updad(n1,i2+1,1)
	  n1=n1-1

c     Move parton list one step downward from i1+1 to n1
	  call updad(n1,i1+1,1)
	  n1=n1-1

c     Share the surplus energy. 
	  if(n1+nn.gt.0)then
		dele=dele/float(n1+nn)
		if(n1.gt.0)then
		  do i4=1,n1
			pp(4,i4)=pp(4,i4)+dele
			if(dele.lt.0.)then
			  if(pp(4,i4).lt.0.)pp(4,i4)=pp(4,i4)-dele
			  pabs=abs(pp(3,i4)) ! yan, 240512
			  if(pabs.ge.pp(4,i4))pp(4,i4)=pp(4,i4)-dele
			endif
		  enddo
		endif
		if(nn.gt.0)then
		  do i4=1,nn
			pn(i4,4)=pn(i4,4)+dele
			if(dele.lt.0.)then
			  if(pn(i4,4).lt.0.)pn(i4,4)=pn(i4,4)-dele
			  pabs=abs(pp(3,i4)) ! yan, 240512
			  if(pabs.ge.pn(i4,4))pn(i4,4)=pn(i4,4)-dele
			endif
		  enddo
		endif
	  endif
	  if(iway.eq.1)then 
	    if(nba.lt.netba)goto 400 ! recycle all the partons remained, do again.
	    if(nba.eq.netba)return
	  endif
        if(iway.eq.2 .and. nba.eq.1)return
	  if(iway.eq.0)goto 400   ! 121204
c   iway=1: when the baryon number equal net baryon, return. Used in creat net baryon.
c   iway=2: it can return when there generate one baryon. 
c   iway=0: check all the probability of parton constituent baryon, then return.
406       continue   ! fail
405     continue
404   continue
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine an_barpro(n1,nqb,i1,iphas,nba,ibarm,netba,rrp,isu,
     c	iway)
c   To compose an anti-baryon
c     n1 : total # of partons (q & qba)
c     nqb : total # of qba (qba is ordered before q)
c	i1: antiquark wanted to compose an antibaryon   
c     iphas=1: with phase space adjudgment
c          =0: without phase space adjudgment
c	ibarm: statistic number of anti-baryon
c	-netba: number of anti-baryons keeping baryon conservation 
c	isu: =1 success
c	     =0 fail
c	iway: a switch
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000,mplis=40000)
      common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
      common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
      common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
      common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
      common/sa24/adj1(40),nnstop,non24,zstop
	dimension rcp(3)
	
	dpmax=adj1(27)
	drmax=adj1(28)	
	isu=1
	nba=0

	if(nqb.lt.3)goto 100   ! number of quarks less than three,return
	kf1=idp(i1)
      do 405 i2=1,nqb
	  if(i2.eq.i1)goto 405
        kf2=idp(i2)
        do 406 i3=1,nqb
		if(i3.eq.i1)goto 406
		if(i3.eq.i2)goto 406
		kf3=idp(i3)

	  sume=pp(4,i1)+pp(4,i2)+pp(4,i3)
	  sump1=pp(1,i1)+pp(1,i2)+pp(1,i3)
	  sump2=pp(2,i1)+pp(2,i2)+pp(2,i3)
        sump3=pp(3,i1)+pp(3,i2)+pp(3,i3)
	  cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
	  if(cm.gt.1.e30)cm=1.e30
	  if(cm.le.0.)goto 406   ! fail 071204
        cm=sqrt(cm)

c     Find out the primary antibaryon from hadron table according to kf1,kf2 
c	 & kf3
        call findb(-kf1,-kf2,-kf3,cm,kfii,amasi,isucc,-1)

        if(isucc.eq.0)goto 406   ! fail

c	Proceed for success   
c	Phase space adjudgment
        if(iphas.eq.1)then
		call phas(i1,i2,i3,isucc,3)
          if(isucc.eq.0)goto 406   ! fail
        endif

c	Proceed for success
	  goto 400
406     continue   ! fail

405   continue
	goto 100
400	continue
      ibarm=ibarm+1
	nba=nba+1

c     Give proper variables to the primary antibaryon
      nnol=nn
      nn=nn+1
	kn(nn,1)=1
      kn(nn,2)=kfii
      kn(nn,3)=0
      pn(nn,5)=amasi
	pn(nn,1)=sump1
	pn(nn,2)=sump2
	pn(nn,3)=sump3
	pnnm=sump1*sump1+sump2*sump2+sump3*sump3
	pnnmm=amasi*amasi+pnnm
	if(pnnmm.gt.1.e30)pnnmm=1.e30
	if(pnnmm.le.0.)pnnmm=1.e-30
	pnnn=sqrt(pnnmm)
	pn(nn,4)=pnnn
	dele=sume-pnnn

c     Arrange produced particle on the surface of sphere with radius
c     rrp and centered at the center of mass

c	Produced hadron is arranged among contituent partons randomly
	pyrx=pyr(1)
	pyry=pyr(1)
	pyrz=pyr(1)
	rn(nn,1)=pyrx*rp(1,i1)+pyry*rp(1,i2)+pyrz*rp(1,i3)
	rn(nn,2)=pyrx*rp(2,i1)+pyry*rp(2,i2)+pyrz*rp(2,i3)
	rn(nn,3)=pyrx*rp(3,i1)+pyry*rp(3,i2)+pyrz*rp(3,i3)

c     Move parton list one step downward from i3+1 to n1
411	call updad(n1,i3+1,1)
	if(i1.gt.i3)i1=i1-1
	if(i2.gt.i3)i2=i2-1
	nqb=nqb-1
	n1=n1-1

c     Move parton list one step downward from i2+1 to n1
	call updad(n1,i2+1,1)
	if(i1.gt.i2)i1=i1-1
	nqb=nqb-1
	n1=n1-1

c     Move parton list one step downward from i1+1 to n1
	call updad(n1,i1+1,1)
	nqb=nqb-1
	n1=n1-1

c	Share the surplus energy.
	if(n1+nn.gt.0)then
	  dele=dele/float(n1+nn)
	  if(n1.gt.0)then
		do i4=1,n1
		  pp(4,i4)=pp(4,i4)+dele
		  if(dele.lt.0.)then
			if(pp(4,i4).lt.0.)pp(4,i4)=pp(4,i4)-dele
			  pabs=abs(pp(3,i4)) ! yan, 240512
			if(pabs.ge.pp(4,i4))pp(4,i4)=pp(4,i4)-dele
		  endif
		enddo
	  endif
	  if(nn.gt.0)then
          do i4=1,nn
		  pn(i4,4)=pn(i4,4)+dele
		  if(dele.lt.0.)then
			if(pn(i4,4).lt.0.)pn(i4,4)=pn(i4,4)-dele
			  pabs=abs(pp(3,i4)) ! yan, 240512
			if(pabs.ge.pn(i4,4))pn(i4,4)=pn(i4,4)-dele
		  endif
		enddo
	  endif
	endif
c   iway=1: creat an antibaryon and return
c   iway=2: creat an antibaryon and a baryon, then return
	if(iway.eq.1 .and. nba.eq.1)return
      if(iway.eq.2 .and. nba.eq.1)then
c     An antibaryon generation must be followed immediately a baryon
c     generation keeping baryon conservation
	  call barpro(n1,nqb,iphas,nbaa,ibarp,netba,rrp,2)
	  return
	endif

100	isu=0
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mespro(n1,nqb,i1,iphas,nme,imes,rrp,isu)
c   Compose a meson
c     n1 : total # of partons (q & qba)
c     nqb : total # of qba (qba is ordered before q)
c	i1: antiquark wanted to compose a meson
c     iphas=1: with phase space adjudgment
c          =0: without phase space adjudgment
c	imes: statistic number of meson
c     isu: =1 success
c          =0 fail
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000,mplis=40000)
      common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
      common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
      common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
      common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
      common/sa24/adj1(40),nnstop,non24,zstop
      dimension rcp(3)

	dpmax=adj1(27)
      drmax=adj1(28)
	isu=1
	nme=0

	if(n1.eq.nqb)goto 100   ! 300105, no quark, return
	kf1=idp(i1)
	do 102 i2=nqb+1,n1
	  kf2=idp(i2)
	  sume=pp(4,i1)+pp(4,i2)
	  sump1=pp(1,i1)+pp(1,i2)
	  sump2=pp(2,i1)+pp(2,i2)
        sump3=pp(3,i1)+pp(3,i2)
	  cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
	  if(cm.gt.1.e30)cm=1.e30
	  if(cm.le.0.)goto 102   ! fail 071204
        cm=sqrt(cm)

c     Find out primary meson from hadronic table according to kf2 & kf1
        call findm(kf2,kf1,cm,kfii,amasi,isucc,1)

        if(isucc.eq.0)goto 102   ! fail
c	Proceed for success
c	Phase space adjudgment
	  if(iphas.eq.1)then
	    call phas(i1,i2,0,isucc,2)
          if(isucc.eq.0)goto 102   ! fail
        endif

c	Proceed for success
	  imes=imes+1
	  nme=nme+1

c     Give proper variables to the primary meson
        nnol=nn
        nn=nn+1
        kn(nn,1)=1
        kn(nn,2)=kfii
        kn(nn,3)=0
        pn(nn,5)=amasi
	  pn(nn,1)=sump1
	  pn(nn,2)=sump2
	  pn(nn,3)=sump3
	  pnnm=sump1*sump1+sump2*sump2+sump3*sump3
	  pnnmm=amasi*amasi+pnnm
        if(pnnmm.gt.1.e30)pnnmm=1.e30
        if(pnnmm.le.0.)pnnmm=1.e-30
	  pnnn=sqrt(pnnmm)
	  pn(nn,4)=pnnn
	  dele=sume-pnnn

c	Produced hadron is arranged among contituent partons randomly
	  pyrx=pyr(1)
	  pyry=pyr(1)
	  rn(nn,1)=pyrx*rp(1,i1)+pyry*rp(1,i2)
	  rn(nn,2)=pyrx*rp(2,i1)+pyry*rp(2,i2)
	  rn(nn,3)=pyrx*rp(3,i1)+pyry*rp(3,i2)

c     Move parton list one step downward since i2+1
111	  call updad(n1,i2+1,1)
        n1=n1-1

c     Move parton list one step downward since i1+1
        call updad(n1,i1+1,1)
        nqb=nqb-1
        n1=n1-1
	  if(n1+nn.gt.0)then
	    dele=dele/float(n1+nn)
	    if(n1.gt.0)then
		  do i4=1,n1
			pp(4,i4)=pp(4,i4)+dele
			if(dele.lt.0.)then
			  if(pp(4,i4).lt.0.)pp(4,i4)=pp(4,i4)-dele
			  pabs=abs(pp(3,i4)) ! yan, 240512
			  if(pabs.ge.pp(4,i4))pp(4,i4)=pp(4,i4)-dele
			endif
		  enddo
		endif
	    if(nn.gt.0)then
            do i4=1,nn
		    pn(i4,4)=pn(i4,4)+dele
		    if(dele.lt.0.)then
		  	  if(pn(i4,4).lt.0.)pn(i4,4)=pn(i4,4)-dele
			  pabs=abs(pp(3,i4)) ! yan, 240512
			  if(pabs.ge.pn(i4,4))pn(i4,4)=pn(i4,4)-dele
		    endif
            enddo
	    endif
	  endif
	  if(nme.eq.1)return
102	continue
100	isu=0   ! 300105
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conse_c(pp,ps,npl,np,nstep)
c     Adjust four momentum conservation by iteration,no more than
c	 5000 iterations
c       pp : four momentum of particle
c       ps : the four momentum should be conserved to
c	npl : order # of the first particle
c       np : order # of last particle
c	nstep : interval of the step
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm   ! 140604
        dimension pp(kszj,5),ps(4),ff(kszj),pxyz(3),arp(3)

        ps4=ps(4)
        do i=1,3
        pxyz(i)=0.
        enddo
        jj=0
100     es=0.
        do i=npl,np,nstep
        es=es+pp(i,4)
        enddo
        fr=es/ps4

        do i=npl,np,nstep
        amas=pp(i,5)
        amas2=amas*amas
        ppm=pp(i,4)
        ppf=ppm/fr
c090700
	den2=ppm*ppm-amas2
	if(den2.le.0.)then
	den2=1.e-15
	endif
c090700
        ff(i)=sqrt(abs(ppf*ppf-amas2)/den2)
        do j=1,3
        ppp=ff(i)*pp(i,j)
        pp(i,j)=ppp
        pxyz(j)=pxyz(j)+ppp
        enddo
        enddo
        do i=1,3
        arp(i)=abs(1.-pxyz(i)/ps(i))
        enddo
        if(abs(1.-fr).le.dep .and. arp(1).le.dep .and. arp(2).le.dep
     c   .and. arp(3).le.dep) goto 200
        do i=1,3
        pxyz(i)=pxyz(i)-ps(i)
        pxyz(i)=pxyz(i)/(float(np-npl)/float(nstep)+1)
        enddo
        do i=npl,np,nstep
        do j=1,3
        pp(i,j)=pp(i,j)-pxyz(j)
        enddo
	pp5=pp(i,5)
	pp52=pp5*pp5
        pp(i,4)=sqrt(pp52+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
        enddo
        jj=jj+1
        if(jj.lt.5000)goto 100
200     return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
	subroutine ffm(ii,kf0,igen,iflav,n1,n2,n3)
c	Quark (antiquark) generation according to Field-Feynman model
c	ii : the order # of source quark (or antiquark)
c	kf0 : flavor code of source quark (or antiquark)
c	igen : times of source quark can excite qqbar pair from
c        vacuum, which is controled by four momenta  ! 080512 sa
c	iflav = 1 : if source parton is quark (kf0>0)
c           =-1 : if source parton is antiquark (kf0<0)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      parameter (mplis=40000)
      common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
      common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
      common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa24/adj1(40),nnstop,non24,zstop
	dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c	 peo(5),pppp(50,2)

	kapa=adj1(15)
	adj17=adj1(17) 
	adj17=max(4.0,adj17) ! 070612, yan
	adj23=adj1(23)

c     The probability of gluon spliting into u,d & s
      adj132=adj1(32)
      prosum=1.+1.+adj132
      prod=1./prosum   ! 0.4286 originally
      pros=adj132/prosum   ! 0.1428 originally
      prods=prod+pros   ! 0.5714 originally

	do i1=1,3
	  rc(i1)=rp(i1,ii)
	enddo
	do i2=1,4
	  p0(i2)=pp(i2,ii)
        p00(i2)=p0(i2)
	enddo
	kf00=kf0
	n30=n3
	e0=p0(4) ! E
c	w0=e0+p0(3) ! E+p_z

	p00(3)=p0(3) ! for the negative direction p_z
	if(p00(3).lt.0.0)p0(3)=-p0(3)  ! for the negative direction p_z
	w0=e0+p0(3) ! E+p_z

	if(w0.lt.0.)return   ! stop generation

c	qqba creation
	igen=0
100	continue

	ie1=0
c	sample the flavor for qqba 
c	amd=pymass(1)   ! d
c	amu=pymass(2)   ! u
c	ams=pymass(3)   ! s
	amd=0.0099D0    !041107
	amu=0.0056D0 
	ams=0.199D0

	amdd=amd*2
	amuu=amu*2
	amss=ams*2
	eg=e0
      if(eg.lt.amuu)return   ! stop generation
      if(eg.ge.amuu .and. eg.lt.amdd)then
        kf=2   ! u
        amasi=amuu
      endif
      if(eg.ge.amdd .and. eg.lt.amss)then
        kf=2   ! u
        amasi=amuu
        if(pyr(1).gt.0.5)then
          kf=1   ! d
          amasi=amdd
        endif
      endif
      if(eg.gt.amss)then
        rand=pyr(1)
        kf=3   ! s
        amasi=amss
        if(rand.gt.pros .and. rand.le.prods)then
          kf=1   ! d
          amasi=amdd
	  endif
        if(rand.gt.prods)then
          kf=2   ! u
          amasi=amuu
        endif
      endif
	kf1=kf
	kf2=-kf
	if(iflav.eq.1)then
	  kf1=-kf
	  kf2=kf
	endif
	
104	continue   

c     Sample transverse momenta from two dimensional Gaussian distribution       
      sigm2=adj1(29)   ! 0.26
      ptmax=adj1(30)   ! 2.
      call tdgaus(sigm2,ptmax,1,pppp)
      p1(1)=pppp(1,1)
      p1(2)=pppp(1,2)
	ppt=p1(1)*p1(1)+p1(2)*p1(2)

c	Sample z from LUND or FF fragmentation function
	if(adj23.eq.1)then
	  call funcz(z1)
	else
	  prr=amasi*amasi/4.+ppt
	  call pyzdis(kf1,kf2,prr,z1)
	endif

	w1=z1*w0 ! E1+p_z1=z*(E0+p_z0)
	amt12=amasi*amasi+ppt ! m_uu^2+p_T^2  
	e1=0.5*(w1+amt12/w1)
	p1(3)=0.5*(w1-amt12/w1)

	if(e1.gt.e0.or.p1(3).lt.0.0)then ! p1(3).lt.0.0 means p_z reversed
	  ie1=ie1+1
	  if(ie1.gt.100)goto 106   ! stop generation
	  goto 104    
	endif	

	p1(4)=e1 ! E
	ee1=sqrt(amt12+p1(3)*p1(3)) ! m^2=p^2+m0^2

c	Fill the q & qba generated into parton list after n3.
c     qbar and q generated are arranged around sourve parton within 0.5 fm 
c     randumly in each one of the three coordinates and has same fourth 
c     coordinate as sourve parton
	do i=1,3
	  rr(i)=pyr(1)*0.5
	  rp(i,n3+1)=rc(i)+rr(i)
        if(pyr(1).gt.0.5)rp(i,n3+1)=rc(i)-rr(i)
        rr(i)=pyr(1)*0.5
	  rp(i,n3+2)=rc(i)+rr(i)
        if(pyr(1).gt.0.5)rp(i,n3+2)=rc(i)-rr(i)
	enddo
	rp(4,n3+1)=rp(4,ii) ! ii is the # of sourse parton
      rp(4,n3+2)=rp(4,ii) ! n+1 and n+2 is the new generated parton

c	Give momentum to q & qba generated (random three momentum method).
	amq=0.5*amasi
      do i1=1,3
	  pii=pyr(1)*p1(i1)
	  pp(i1,n3+1)=pii
	  pnn(1,i1)=pii
	  pp(i1,n3+2)=p1(i1)-pii
	  pnn(2,i1)=p1(i1)-pii
	enddo

	if(p00(3).lt.0.0)then  ! for the negative direction p_z
	  pp(3,n3+1)=-pp(3,n3+1)
	  pp(3,n3+2)=-pp(3,n3+2)
	endif

	pn11=pnn(1,1)
      pn12=pnn(1,2)
      pn13=pnn(1,3)
	pp(4,n3+1)=sqrt(amq*amq+pn11*pn11+pn12*pn12+pn13*pn13)
      pn21=pnn(2,1)
      pn22=pnn(2,2)
      pn23=pnn(2,3)
      pp(4,n3+2)=sqrt(amq*amq+pn21*pn21+pn22*pn22+pn23*pn23)

c	Give other properties to q and qba generated
	rmp(n3+1)=amq
      rmp(n3+2)=amq
	idp(n3+1)=kf1
	idp(n3+2)=kf2
	n3=n3+2

c	Give proper variables to the remnant parton
	w1c=w0-w1   ! conservation
	do i3=1,3
	  p1c(i3)=p0(i3)-p1(i3)   
	enddo
c     p0 is the original parton, p1 is the generated parton pair, 
c     p1c is the original parion after generated parton pair.

c	amt1c2=(w1c-2*p1c(3))*w1c
	p1c1=p1c(1)
	p1c2=p1c(2)
	amt1c2=amq*amq+p1c1*p1c1+p1c2*p1c2
c	e1c=0.5*(w1c+amt1c2/w1c)
	e1c=w1c-p1c(3)   
	p1c(4)=e1c

	igen=igen+1

	kf0=kf00
	if(kf0.gt.0)iflav=1
	if(kf0.lt.0)iflav=-1
	w0=w1c
	e0=e1c
	do i3=1,4
	  p0(i3)=p1c(i3)
	enddo
	if(e0.le.adj17)goto 106   ! stop generation
	if(w0.lt.0.)goto 106   ! stop generation
	rand=pyr(1)  ! 080512
c	rand=rlu(1)
	adj16=adj1(16)
c	write(9,*)'adj16=',adj16
	if(rand.gt.adj16)goto 106 ! 080512
	w00=p0(4)+p0(3)

	goto 100

106	continue
	if(igen.eq.0)return
	do i3=1,3
	  pp(i3,ii)=p0(i3)
c	rp(i3,ii)=rc(i3)
	enddo
	if(p00(3).lt.0.0)then  ! for the negative direction p_z
	  pp(3,ii)=-pp(3,ii)
	endif
	pp(4,ii)=p0(4)

	do i1=n30+1,n3
	  j1=i1-n30
	  do i2=1,4
	    pnn(j1,i2)=pp(i2,i1)
	  enddo
	  pnn(j1,5)=rmp(i1)
	enddo
	ii1=n3-n30   ! new generated parton from ii-th parton

c	Incluse remnant of ii-th parton
	ii1=ii1+1   
	do i2=1,4
	  pnn(ii1,i2)=pp(i2,ii)
	enddo
	pnn(ii1,5)=rmp(ii)

	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine funcz(zz)
c     Sample zz from fragmentation function using selecting sample method.	  
c	Distribution function: f(z)dz=[1-a+3*a*(1-z)**2]dz.
c	The largest value of function: fmax=f(0)=1-a+3*a.
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      common/sa24/adj1(40),nnstop,non24,zstop
	a=adj1(6)   ! 0.77
	b=3.*a   ! 3*0.77
	fmax=1-a+b   
100	ran1=pyr(1)
	ran2=pyr(1)
	fm=ran1*fmax
c101204	fz=1.-0.77+3*0.77*(1.-ran2)**2.
c101204	if(fm.le.fz)goto 100
	fz=1.-a+b*(1.-ran2)**2.   ! 101204
      if(fm.gt.fz)goto 100   ! 101204
	zz=ran2  	
	return 
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine phas(i1,i2,i3,isucc,j)
c	The phase space judgement.
c	j=2 for meson
c	j=3 for baryon
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(mplis=40000)
      common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
      common/sa24/adj1(40),nnstop,non24,zstop
	dimension ri1(3),ri2(3),ri3(3),pi1(3),pi2(3),pi3(3)
	delc=adj1(22)

	if(j.eq.2)goto 100   ! for meson
c	 proceed for baryon 
        ri1(1)=rp(1,i1)
        ri1(2)=rp(2,i1)
        ri1(3)=rp(3,i1)
        ri2(1)=rp(1,i2)
        ri2(2)=rp(2,i2)
        ri2(3)=rp(3,i2)
        ri3(1)=rp(1,i3)
        ri3(2)=rp(2,i3)
        ri3(3)=rp(3,i3)
        ri121=ri1(1)-ri2(1)
        ri122=ri1(2)-ri2(2)
        ri123=ri1(3)-ri2(3)
        ri131=ri1(1)-ri3(1)
        ri132=ri1(2)-ri3(2)
        ri133=ri1(3)-ri3(3)
        ri231=ri2(1)-ri3(1)
        ri232=ri2(2)-ri3(2)
        ri233=ri2(3)-ri3(3)
	delr12=sqrt(ri121*ri121+ri122*ri122+ri123*ri123)
	delr13=sqrt(ri131*ri131+ri132*ri132+ri133*ri133)
	delr23=sqrt(ri231*ri231+ri232*ri232+ri233*ri233)
        pi1(1)=pp(1,i1)
        pi1(2)=pp(2,i1)
        pi1(3)=pp(3,i1)
        pi2(1)=pp(1,i2)
        pi2(2)=pp(2,i2)
        pi2(3)=pp(3,i2)
        pi3(1)=pp(1,i3)
        pi3(2)=pp(2,i3)
        pi3(3)=pp(3,i3)
        pi121=pi1(1)-pi2(1)
        pi122=pi1(2)-pi2(2)
        pi123=pi1(3)-pi2(3)
        pi131=pi1(1)-pi3(1)
        pi132=pi1(2)-pi3(2)
        pi133=pi1(3)-pi3(3)
        pi231=pi2(1)-pi3(1)
        pi232=pi2(2)-pi3(2)
        pi233=pi2(3)-pi3(3)
	delp12=sqrt(pi121*pi121+pi122*pi122+pi123*pi123)
	delp13=sqrt(pi131*pi131+pi132*pi132+pi133*pi133)
	delp23=sqrt(pi231*pi231+pi232*pi232+pi233*pi233)
	del12=delr12*delp12
	del13=delr13*delp13
	del23=delr23*delp23

	if(del12.le.delc.and.del13.le.delc.and.del23.le.delc)then
	  isucc=1
	else
	  isucc=0
	endif
	return

100	continue   ! for meson
        ri1(1)=rp(1,i1)
        ri1(2)=rp(2,i1)
        ri1(3)=rp(3,i1)
        ri2(1)=rp(1,i2)
        ri2(2)=rp(2,i2)
        ri2(3)=rp(3,i2)
        ri121=ri1(1)-ri2(1)
        ri122=ri1(2)-ri2(2)
        ri123=ri1(3)-ri2(3)
        delr=sqrt(ri121*ri121+ri122*ri122+ri123*ri123)
        pi1(1)=pp(1,i1)
        pi1(2)=pp(2,i1)
        pi1(3)=pp(3,i1)
        pi2(1)=pp(1,i2)
        pi2(2)=pp(2,i2)
        pi2(3)=pp(3,i2)
        pi121=pi1(1)-pi2(1)
        pi122=pi1(2)-pi2(2)
        pi123=pi1(3)-pi2(3)
        delp=sqrt(pi121*pi121+pi122*pi122+pi123*pi123)
        delrp=delr*delp
	if(delrp.le.delc)then
	  isucc=1
	else
	  isucc=0
	endif
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tdgaus(v,pmax,np,pp)
c.... 2-d Gaussian distribution with width v, i.e., e^(-p2/v)dp2, 0<p2<pmax
c.... set pmax < 0 if pmax should be infinity.
c.... np : the total # of particles wanted to sample their transverse momenta
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/sa33/smadel,ecce,parecc,iparres   ! 220312 240412
      dimension pp(50,2)
      do 30 i=1,np
        p2 = 0
        if(v.le.1.e-8)return
        if(pmax.lt.1.E-9)return
        if(pmax.lt.0)then
          a = 1.
          goto 10
        endif
        aa=-pmax/v
        if(aa.lt.-70)then
          a=1.
          goto 10
        endif
        a = 1. - exp(aa)
10      p2 = -v*log(max(1.e-20,1. - a*pyr(1)))
        if(p2.LT.0.)goto 10
        ps=sqrt(p2)
        fi=2.*3.1415926*pyr(1)
c220312 randomly sample [px,py] on circle of sphere with radius ps
c220312	pp(i,1)=ps*cos(fi)
c220312	pp(i,2)=ps*sin(fi)
c220312 randomly sample [px,py] on circle of ellipsoid with half major axis 
c220312 of ps*(1+smadel) and half minor axis of ps*(1-smadel)
	pp(i,1)=ps*cos(fi)*(1+smadel)   ! 220312
	pp(i,2)=ps*sin(fi)*(1-smadel)   ! 220312
c220312	note: ps is not in the dimension list
30    continue
      return
      end
