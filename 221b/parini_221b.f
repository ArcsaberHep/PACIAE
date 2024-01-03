	subroutine parini(time_neu,time_par,time_had,parp21,parp22
     c	 ,win,psno,ijk)   ! 111010 121110
c	generate an event (parton or hadron) for a relativistic  
c	 nucleus-nucleus collision based on 'pythia' 
c	it was composed by Ben-Hao Sa on 30/July/05
c	the intermediate working arries are in common statement 'sa2' 
c	output message is in 'PYJETS'
        parameter(kszj=40000,ksz1=30)
        parameter(nsize=240000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	dimension bst(4),bbb(3),bb(3)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c	those variables in above common blocks are defined in 'jetset'
        COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
	COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)   ! 221203
c	those variables in above common block are defined in 'pythia'
	COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
        common/sa6/kfmaxi,nwhole
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
c	csen: e+p total x section in fm^2
c080104
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
	common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph 
     c	 ,vnlep   ! 260314
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 280809
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21,adiv,gpmax,nnc   !   070417
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
c080104
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5) 
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
c	iii : number of current event 
c	neve : total number of events 
c	bp : impact parameter
c       c17(i,1-3) : position of i-th nucleon (origin is set at the center of
c       target nucleus)
c       tp(i) : time of i-th nucleon counted since collision of two nuclei
c       p17(i,1-4) : four momentum of i-th nucleon 
c       ishp(i)=1 if i-th particle inside the simulated volume
c              =0 if i-th particle outside the simulated volume
c	cspsn : total cross section of J/Psi (Psi') + n
c	cspsm : total cross section of J/Psi (Psi') + meson
c	iabsb = 0 : without J/Psi (Psi') + baryon
c	      = 1 : with J/Psi (Psi') + baryon
c	iabsm = 0 : without J/Psi (Psi') + meson
c	      = 1 : with J/Psi (Psi') + meson
        common/count/isinel(600)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
c       ecsen: largest interaction distance between lepton and p   ! 060813
c060813 120214 note: incident lepton collides with nucleon in nucleus once 
c        only due to low total x-section. that collision is the one with 
c        lowest minimum approaching distance
c       sig (fm^2): cross section of pion + pion to kaon + kaon
c       edipi: largest interaction distance between two pions.
c       epin: largest interaction distance between pion and nucleon.
c       ekn: largest interaction distance between kaon and nucleon.
c       ecsnn: largest interaction distance between two nucleons.
c       t0 : average proper formation time at rest.
c       ddt : time accuracy used in parton initiation 
c	time accuracy used in parton cascade is dddt 
c       rou0 : normal nucleon density.
c       rao : enlarged factor for the radius of simulated volume.

	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c	nap,nat,nzp,nzt,pio
c 	nap and nzp (nat and nzt) are the mass and charge numbers of 
c	 projectile (target) nucleus
c       r0p=rnp     : the standard radius of projectile
c       r0t=rnt     : the standard radius of target
	common/ctllist/nctl,noinel(600),nctl0,noel
c       nctl: number of collision pairs in current collision list
c       nctl0: number of collision pairs in initial collision list
c       noinel(1): statistics of nn elas. colli.;
c	noinel(592): statistics of nn colli. calling pythia
c	noel: statistics of the blocked nn colli. #
	common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
	common/sa23/kpar,knn,kpp,knp,kep   ! 200601 060813
	common/sa30/vneump,vneumt
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 270312 240412
c	vneump (vneumt): # of participant nucleons of projectile (target) 
c	 calculated geometrically or by Glauber model 

	dimension peo(4),pi(4),pj(4),xi(4),xj(4)
	dimension inoin(kszj)
        dimension lc(nsize,5),tc(nsize),tw(nsize)

	kpar=0
        knn=0
        kpp=0
        knp=0
	kep=0   ! 060813
c060813 120214 kep uses to statistics of # of times calling pythia in 
c        case of lepton is projectile
	if(iii.eq.1)then
	if(nchan.eq.1)then
c	Non Single Difractive (NSD)  
        msel=0
        msub(11)=1
        msub(12)=1
        msub(13)=1
        msub(28)=1
        msub(53)=1
        msub(68)=1
c	msub(91)=1
c       msub(92)=1
c       msub(93)=1
	msub(94)=1
        msub(95)=1
	elseif(nchan.eq.2)then
c	qqb --> gamma^*/Z^0 used to generate Drell-Yen
        msel=0
        msub(1)=1
	elseif(nchan.eq.3)then
c	J/psi production
	msel=0
	msub(86)=1
c        msub(87)=1
c        msub(88)=1
c        msub(89)=1
        elseif(nchan.eq.4)then
c 	heavy-flavor production
	msel=0
        msub(81)=1
        msub(82)=1
        elseif(nchan.eq.5)then
c	direct photon
	msel=0
        msub(14)=1
        msub(18)=1
        msub(29)=1
        msub(114)=1
        msub(115)=1
        elseif(nchan.eq.6)then
c	soft only
        msel=0
	msub(92)=1
	msub(93)=1
	msub(94)=1
	msub(95)=1
	msub(96)=0   ! 270705
	else
	msel=1
c	pythia 
	endif
	endif

c240209 flavor code 22: hardonic decay photon
c                   44: prompt direct photon (<- pythia)
c                   55: photon from parton-parton scattering
c                       qg->q(gamma) and q(-q)->g(gamma)
c                   66: hardonic direct photon
c240209                 pi+pi->rho+(gamma) and pi+rho->pi+(gamma)
c270312	initiation of x,y,xy,x^2,y^2,and sump (statistics of the number of 
c270312	 nucleons in overlap region)   ! 131212
	sumx=0.
	sumy=0.
	sumxy=0.   ! 131212
	sumx2=0.
        sumy2=0.
	sump=0.
c270312

c	initiate the nucleus-nucleus collision system

c	creat the initial particle list (nucleon)
c230311 in position phase space
c191110
c       A+B (nucleus-nucleus)   ! 230311
        if(ipden.eq.1 .and. itden.eq.1)then   !! 230311
c	distribute projectile nucleons
	napt=nap
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.gt.27.and.napt.lt.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
        elseif(napt.eq.56)then
        alpt=0.49
        elseif(napt.eq.64)then
        alpt=0.49     
        elseif(napt.eq.108)then
        alpt=0.495
        elseif(napt.eq.184)then
        alpt=0.53
        elseif(napt.eq.197)then
        alpt=0.54
        elseif(napt.eq.207)then
        alpt=0.545
        elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	alp=alpt
	r0=r0p
	am=suppm   ! upper bound in sampling the radius of projectile nucleon
	ac=suppc   ! maximum radius for projectile
	ratps=vneump/nap   ! ratio of projectile participant nucleons to total 
c	if(iii.eq.10)write(9,*)'nap,vneump,ratps=',nap,vneump,ratps
	do i1=1,nap
	rann=pyr(1)
	if(rann.lt.ratps)then
c       sample position of projectile nucleon in overlap region of colliding 
c	 nuclei
	call arrove(i1,1,sumx,sumy,sumxy,sumx2,sumy2,sump,
     c	 alp,r0,am,ac)   ! 270312 131212 101014
	else
c	sample position of projectile nucleon according to Woods-Saxon
c	 distribution
	call woodsax_samp(i1,1,alp,r0,am,ac,1)   ! 230311
c230311 last argument here is 'iway', iway=1: particle i1 must be outside the 
C230311  overlap region of colliding nuclei, iway=0: no more requirement 
	endif
	enddo
c230311
c	distribute target nucleons
	napt=nat
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.gt.27.and.napt.lt.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
        elseif(napt.eq.56)then
        alpt=0.49
        elseif(napt.eq.64)then
        alpt=0.49     
        elseif(napt.eq.108)then
        alpt=0.495
        elseif(napt.eq.184)then
        alpt=0.53
        elseif(napt.eq.197)then
        alpt=0.54
        elseif(napt.eq.207)then
        alpt=0.545
        elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	alp=alpt
	r0=r0t
	am=suptm   ! upper bound in sampling the radius of target
	ac=suptc   ! maximum radius for target
	ratps=vneumt/nat   ! ratio of target participant nucleons to total 
	do i1=1,nat
	i2=i1+nap
	rann=pyr(1)
	if(rann.lt.ratps)then
c       sample position of target nucleon in overlap region of colliding nuclei
	call arrove(i2,0,sumx,sumy,sumxy,sumx2,sumy2,sump,
     c	 alp,r0,am,ac)   ! 270312 131212 101014
	else
c	sample position of target nucleon according to Woods-Saxon
c	 distribution
	call woodsax_samp(i2,0,alp,r0,am,ac,1)
	endif
	enddo
c       p+A or lepton+A   ! 060813 120214
        elseif((ipden.eq.0.or.ipden.gt.1) .and. itden.eq.1)then ! 060813 050214
c	write(9,*)'be. eA spatial initiation sump=',sump
        do i=1,3
        c17(1,i)=0.
        enddo
c	distribute target nucleons
	napt=nat
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.gt.27.and.napt.lt.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
        elseif(napt.eq.56)then
        alpt=0.49
        elseif(napt.eq.64)then
        alpt=0.49     
        elseif(napt.eq.108)then
        alpt=0.495
        elseif(napt.eq.184)then
        alpt=0.53
        elseif(napt.eq.197)then
        alpt=0.54
        elseif(napt.eq.207)then
        alpt=0.545
        elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	alp=alpt
	r0=r0t
	am=suptm   ! upper bound in sampling the radius of target
	ac=suptc   ! maximum radius for target
	do i1=1,nat
	i2=i1+nap
	call woodsax_samp(i2,0,alp,r0,am,ac,0)
	enddo
c	write(9,*)'af. eA spatial initiation sump=',sump
c240513
c       A+p
        elseif(ipden.eq.1 .and. itden.eq.0)then   !!
c	distribute projectile nucleons
	napt=nap
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.gt.27.and.napt.lt.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
        elseif(napt.eq.56)then
        alpt=0.49
        elseif(napt.eq.64)then
        alpt=0.49     
        elseif(napt.eq.108)then
        alpt=0.495
        elseif(napt.eq.184)then
        alpt=0.53
        elseif(napt.eq.197)then
        alpt=0.54
        elseif(napt.eq.207)then
        alpt=0.545
        elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	alp=alpt
	r0=r0p
	am=suppm   ! upper bound in sampling the radius of projectile nucleon
	ac=suppc   ! maximum radius for projectile
	do i1=1,nap
	call woodsax_samp(i1,1,alp,r0,am,ac,0)
	enddo
        do i=1,3
        c17(nap+1,i)=0.
        enddo
c240513
c	p+p or lepton+p   ! 070417
c070417	elseif((ipden.eq.0 .and. itden.eq.0) .or. 
c	c	 (itden.eq.0 .and. ipden.ge.11))then   !! 070417
c	do i=1,3
c	c17(1,i)=0.
c	c17(2,i)=0.
c	enddo        
        endif   !!
c230311
c270312
	if(sump.ne.0.)then
	asumx=sumx/sump
	sigmx2=sumx2/sump-asumx*asumx
	asumy=sumy/sump
        sigmy2=sumy2/sump-asumy*asumy
	asumxy=sumxy/sump   ! 131212
	sigmxy=asumxy-asumx*asumy   ! 131212
        sigmsu=sigmy2+sigmx2   ! change from sigmxy to sigmsu 131212
	sigmde=sigmy2-sigmx2   ! 131212
	argu=sigmde*sigmde+4*sigmxy*sigmxy   ! 131212   
c	reaction plane eccentricity of participant nucleons
c131212	if(sigmsu.gt.0.)ecce=(sigmy2-sigmx2)/sigmsu   
c131212
c       participant eccentricity of participant nucleons
	if(argu.gt.0. .and. sigmsu.gt.0.)
     c	 ecce=sqrt(argu)/sigmsu !131212
c	calculate \epsilon{2}=\sqrt(<\epsilon_{part}^2>)
cc	ecce=ecce*ecce
c	note, \epsilon{2} should be \sqrt(aecceo), aecceo is a output 
c	 in paciae_21b.f
c	calculate transverse overlap area
	argu1=sigmx2*sigmy2-sigmxy*sigmxy
	if(argu1.gt.0.)secce=3.1416*sqrt(argu1) ! overlop area 250113 
c131212
c	assuming ecce=geometric eccentricity of ellipsoid (\sqrt{(1-b^2/a^2)}) 
c	 with half major axis b=pt*(1+smadel) and half minor axis 
c	 a=pt*(1-smadel), the resulted smadel=-ecce*ecce/4 (if neglecting 
c	 the samll term of ecce*ecce*(-2*smadel+smadel*smadel)
	ecc2=ecce*ecce   ! 250113
	smadel_a=parecc*ecc2/4. ! approximated deformation parameter 250113
c250113
	delta1=(2.-ecc2+2.*(1.-ecc2)**0.5)/ecc2   
        delta2=(2.-ecc2-2.*(1.-ecc2)**0.5)/ecc2
	if(delta1.le.1.)then
	smadel=parecc*delta1  ! exact deformation parameter
	elseif(delta2.le.1.)then
        smadel=parecc*delta2  ! exact deformation parameter
	else
	endif
c	write(9,*)'ecce,smadel_a,smadel=',ecce,smadel_a,smadel
c250113
c	here a sign change is introduced because of asymmetry of initial 
c	 spatial space is opposed to the final momentum space 
c	write(9,*)'vneump,vneumt,sump,ecce,smadel=',
c     c	 vneump,vneumt,sump,ecce,smadel
	endif	
c270312
c191110
c	for A+B,p+A,A+p,lepton+A 230311 240513 060813 120214
	r0pt=r0p+r0t
c240513	if(itden.ne.0)then   ! 060813
	do i=1,nap
	c17(i,1)=c17(i,1)+bp
	enddo
c240513	endif   ! 230311 
c191110
c	if(iii.eq.10)then
c	write(9,*)'af. boost,bp=',bp
c	do i=1,nat+nap
c	write(9,*)c17(i,1),c17(i,2),c17(i,3)
c	enddo
c	endif
c191110
c	the beam direction is identified as the z axis
c	the origin of position space is set on the center of
c	target nucleus and the origin of time is set at the moment of 
c	first nn colission assumed to be 1.e-5
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	continue
c230311 in momentum phase space
	if(ifram.eq.1)then
	ep1=0.5*win
	et1=ep1
	ep2=0.5*win
	et2=ep2
	pm2=pmas(pycomp(2212),1)**2	
	pp1=sqrt(ep1*ep1-pm2)
	pt1=-sqrt(et1*et1-pm2)
	pm2=pmas(pycomp(2112),1)**2	
	pp2=sqrt(ep2*ep2-pm2)
	pt2=-sqrt(et2*et2-pm2)
c260314 set four momentum and mass for incident lepton and nucleon 
	if(ipden.ge.11.and.ipden.le.16)then   ! in cms
	pincl(1)=0.
	pincl(2)=0.
	pincl(4)=0.5d0*win
	pincl(5)=pmas(pycomp(ipden),1)
	pincl3=pincl(4)*pincl(4)-pincl(5)*pincl(5)
	pincl3=dmax1(pincl3,1.d-20)
	pincl(3)=dsqrt(pincl3)
	pinch(1)=0.
	pinch(2)=0.
	pinch(4)=0.5d0*win
	pinch(5)=pmas(pycomp(2212),1)
	pinch3=pinch(4)*pinch(4)-pinch(5)*pinch(5)
	pinch3=dmax1(pinch3,1.d-20)
	pinch(3)=dsqrt(pinch3)
	endif
c260314
	endif	
	if(ifram.eq.0)then
	pp1=win
	pt1=1.e-20
	pp2=win
	pt2=1.e-20
	pm2=pmas(pycomp(2212),1)**2	
	ep1=sqrt(pp1*pp1+pm2)
	et1=sqrt(pt1*pt1+pm2)
	pm2=pmas(pycomp(2112),1)**2
	ep2=sqrt(pp2*pp2+pm2)
	et2=sqrt(pt2*pt2+pm2)
c260314 set four momentum and mass for incident lepton and nucleon
	if(ipden.ge.11.and.ipden.le.16)then   ! in lab
	pincl(1)=0.
	pincl(2)=0.
	pincl(3)=win
	pincl(5)=pmas(pycomp(ipden),1)
	pincl4=pincl(3)*pincl(3)+pincl(5)*pincl(5)
	pincl4=dmax1(pincl4,1.d-20)
	pincl(4)=dsqrt(pincl4)
	pinch(1)=0.
	pinch(2)=0.
	pinch(3)=0.
	pinch(5)=pmas(pycomp(2212),1)
	pinch(4)=pinch(5)
	endif
c260314
	endif
c260314
c	if(ipden.ge.11.and.ipden.le.16)then
c	write(9,*)'pincl=',(pincl(i),i=1,5)
c	write(9,*)'pinch=',(pinch(i),i=1,5)
c	endif
c260314	
100	inzp=iabs(nzp)
        inzt=iabs(nzt)
	do i=1,nap
	p17(i,1)=0.
	p17(i,2)=0.
	if(i.le.inzp)then
	p17(i,3)=pp1
	p17(i,4)=ep1
	else
	p17(i,3)=pp2
	p17(i,4)=ep2
	endif
	enddo
	napt=nap+nat
	do i=nap+1,napt
	p17(i,1)=0.
	p17(i,2)=0.
	if(i.le.nap+inzt)then	
	p17(i,3)=pt1
	p17(i,4)=et1
	else
	p17(i,3)=pt2
	p17(i,4)=et2
	endif	
	enddo

	do i=1,napt
	tp(i)=0.
	enddo
c	write(9,*)'nap,nzp,nat,nzt=',nap,nzp,nat,nzt
c	write(9,*)'p17(1,3),p17(2,3)=',p17(1,3),p17(2,3)
c	write(9,*)'p17(1,4),p17(2,4)=',p17(1,4),p17(2,4)

c	calculate the velocity of the nucleus-nucleus CM in LAB or nucleon-
c	 nucleon CM system
	bst(1)=p17(1,1)*nap+p17(nap+1,1)*nat
	bst(2)=p17(1,2)*nap+p17(nap+1,2)*nat
	bst(3)=p17(1,3)*nap+p17(nap+1,3)*nat
	bst(4)=p17(1,4)*nap+p17(nap+1,4)*nat
	bst(1)=-bst(1)/bst(4)
	bst(2)=-bst(2)/bst(4)
	bst(3)=-bst(3)/bst(4)

	n=0
        naf=0   ! 080104
        nsa=0
        do i1=1,kszj
        do j1=1,5
	k(i1,j1)=0
	p(i1,j1)=0.
	v(i1,j1)=0.
        kaf(i1,j1)=0
        paf(i1,j1)=0.
        vaf(i1,j1)=0.
        ksa(i1,j1)=0
        psa(i1,j1)=0.
        vsa(i1,j1)=0.
        enddo
        ishp(i1)=0
        tau(i1)=0.
        enddo

        nctl=0
        do i=1,nsize
        do j=1,5
        lc(i,j)=0
        enddo
        tc(i)=0.
        tw(i)=0.
        enddo

        do i=1,100
        numb(i)=0
        numbs(i)=0
        enddo

c      '1 -> |nzp|' are projectile protons or lepton, '|nzp|+1 -> nap' 
c	 are projectile neutrons (if projectile is nucleus); 'nap+1 -> 
c	 nap+nzt' are targer protons, the rest are target nuctrons after 
c	 initiated 'pyjets'   ! 060813 120214
	n=napt
	do i=1,n
	k(i,1)=1
	k(i,2)=2112
	p(i,5)=pmas(pycomp(2112),1)
	if((i.le.iabs(nzp).and.ipden.lt.2).or.(i.gt.nap .and. i.le.nap+
     c	 nzt))then   ! 060813 120214
	k(i,2)=2212
	p(i,5)=pmas(pycomp(2212),1)
c060813 120214
        elseif(i.le.nap.and.(ipden.eq.11.and.nzp.eq.-1))then
        k(i,2)=11
        p(i,5)=pmas(pycomp(11),1)
        elseif(i.le.nap.and.(ipden.eq.11.and.nzp.eq.1))then
        k(i,2)=-11
        p(i,5)=pmas(pycomp(11),1)
        elseif(i.le.nap.and.(ipden.eq.12.and.nzp.eq.-1))then
        k(i,2)=12
        p(i,5)=pmas(pycomp(12),1)
        elseif(i.le.nap.and.(ipden.eq.12.and.nzp.eq.1))then
        k(i,2)=-12
        p(i,5)=pmas(pycomp(12),1)
        elseif(i.le.nap.and.(ipden.eq.13.and.nzp.eq.-1))then
        k(i,2)=13
        p(i,5)=pmas(pycomp(13),1)
        elseif(i.le.nap.and.(ipden.eq.13.and.nzp.eq.1))then
        k(i,2)=-13
        p(i,5)=pmas(pycomp(13),1)
	elseif(i.le.nap.and.(ipden.eq.14.and.nzp.eq.-1))then
        k(i,2)=14
        p(i,5)=pmas(pycomp(14),1)
        elseif(i.le.nap.and.(ipden.eq.14.and.nzp.eq.1))then
        k(i,2)=-14
        p(i,5)=pmas(pycomp(14),1)
        elseif(i.le.nap.and.(ipden.eq.15.and.nzp.eq.-1))then
        k(i,2)=15
        p(i,5)=pmas(pycomp(15),1)
        elseif(i.le.nap.and.(ipden.eq.15.and.nzp.eq.1))then
        k(i,2)=-15
        p(i,5)=pmas(pycomp(15),1)
        elseif(i.le.nap.and.(ipden.eq.16.and.nzp.eq.-1))then
        k(i,2)=16
        p(i,5)=pmas(pycomp(16),1)
        elseif(i.le.nap.and.(ipden.eq.16.and.nzp.eq.1))then
        k(i,2)=-16
        p(i,5)=pmas(pycomp(16),1)
        else
c060813 120214
	endif
	do j=1,3
	p(i,j)=p17(i,j)
	v(i,j)=c17(i,j)
	enddo
	p(i,4)=p17(i,4)
	v(i,4)=tp(i)
	enddo
c	write(9,*)'initial particle list n=',n
c	do i1=1,n
c	write(9,*)'i1,k(i1,2)=',i1,k(i1,2)
c	enddo
500	continue    ! 031103
c	v, vbh and vsa arraies are the position four vector
c	note: for v etc., we do not take care of their fifth component
c	 for array k, we take care of only first three components
c	write(9,*)'after initialion iii=',iii   !s
c	call psum(p,1,n,peo)   !!
c	write(9,*)'psum, after initializing nucleus-nucleus collision'   !!
c	write(9,*)peo   !!

c	boost PYJETS into cms of initial nucleus-nucleus collision system 
c	 from lab or initial nucleon-nucleon cms system.
c	call pyrobo(1,n,0.0,0.0,bst(1),bst(2),bst(3))
c	Lorentz contraction
	bzp3=0.
	bzp4=0.
	bzt3=0.
	bzt4=0.
	do i=1,nap
	bzp3=bzp3+p(i,3)
	bzp4=bzp4+p(i,4)
	enddo
	do i=nap+1,napt
	bzt3=bzt3+p(i,3)
	bzt4=bzt4+p(i,4)
	enddo
	bzp=bzp3/bzp4
	bzt=bzt3/bzt4
	gamp=1./sqrt(dmax1(1.d-20,(1.0d0-bzp*bzp)))
c060813 120214 no Lowrantz contraction for incident lepton
	if(ipden.ge.2)gamp=1.   ! 060813 120214
	gamt=1./sqrt(dmax1(1.d-20,(1.0d0-bzt*bzt)))
c	try no lorentz contract for target
c	gamt=1.
	do i=1,nap
	c17(i,3)=c17(i,3)/gamp
	v(i,3)=v(i,3)/gamp	
	enddo
	do i=nap+1,napt
	c17(i,3)=c17(i,3)/gamt
	v(i,3)=v(i,3)/gamt
	enddo
c	filter out those kind of particles wanted to study and make 
c	 the order of proton, neutron, ... (cf. 'filt')
	call filt
c060813 120214
c	since lepton was moved to last position after calling filt, one has to 
c	 remove it to the first position in pyjets  
	if(ipden.ge.2)call ltof(n)    
c060813 050214
c	write(9,*)'af. ltof particle list n=',n
c	do i1=1,n
c	write(9,*)'i1,k(i1,2)=',i1,k(i1,2)
c	enddo
c	write(9,*)'numbs=',(numbs(i1),i1=1,kfmax)
	nsa=n
        do i=1,n
        do j=1,5
	ksa(i,j)=k(i,j)
        psa(i,j)=p(i,j)
        vsa(i,j)=v(i,j)
        enddo
	ishp(i)=1
	enddo
        do m=1,kfmax
        numb(m)=numbs(m)
        enddo
c       note: particle list is composed of the arraies in common block
c        'sa2', the array 'ishp' in common block 'wz', the array 'tau' in 
c        common block 'sa4', and the array 'numb' in common block 'sa5'
	time=0.
	irecon=0

	call copl(time)
c       calculate the position for the center of mass of the
c	non-freeze-out system. The distance of a particle, when checking
c	is it freezing out or not, is measured with respect to this center
c       creat the initial collision list, note: be sure that the initial  
c	collision list must not be empty
	call ctlcre(lc,tc,tw)
c	write(9,*)'af. ctlcre nctl=',nctl
c	do i1=1,nctl
c	write(9,*)'lc(i1,1),lc(i1,2),tc(i1)=',lc(i1,1),lc(i1,2),tc(i1)
c	enddo

c070417 move origin of time to collision time of first nucleon-nucleon collision
c	find out colli. pair with least colli. time
	call find(icp,tcp,lc,tc,tw,0)
c	write(9,*)'af. find icp,tcp=',icp,tcp
	if(icp.eq.0)stop 'initial collision list is empty'   !
	time=tcp
c070417 perform classical Newton motion in Lab. system for all particles 
	call his(time,lc,tc,tw,istop)
c	write(9,*)'af. his nctl=',nctl
c	do i1=1,nctl
c	write(9,*)'lc(i1,1),lc(i1,2),tc(i1)=',lc(i1,1),lc(i1,2),tc(i1)
c	enddo
	do ij=1,nsa
	vsa(ij,4)=0.
	enddo
	do ij=1,nctl
	tc(ij)=tc(ij)-time+1.e-5
	enddo
c	write(9,*)'af. moving to time origin nctl=',nctl
c	do i1=1,nctl
c	write(9,*)'lc(i1,1),lc(i1,2),tc(i1)=',lc(i1,1),lc(i1,2),tc(i1)
c	enddo
	time=0.
	call copl(time) 
400	continue
	time_neu=0.d0   ! 121110
        time_par=0.d0   ! 111010
        time_had=0.d0   ! 111010
c	write(9,*)'be. scat, iii,nsa,naf,nctl=',iii,nsa,naf,nctl   ! sa
c       administrate a nucleus-nucleus collision 
        call scat(time_neu,time_par,time_had,lc,tc,tw,win,parp21,
     c	 parp22,psno,ijk,ipau,irecon,gamt)   ! 211107 111010
c	write(9,*)'af. scat, iii,nsa,ijk,nctl=',iii,nsa,ijk,nctl   ! sa
        if(ijk.eq.1)return   

800	continue
c-----------------------------------------------------------------------
c       put 'saf' after 'sa2'
        do l=1,naf
        l1=nsa+l
        do m=1,5
        ksa(l1,m)=kaf(l,m)
        psa(l1,m)=paf(l,m)
        vsa(l1,m)=vaf(l,m)
        enddo
        enddo
        nsa=nsa+naf
c       'sa2' to 'pyjets' after finish calculation
        call trans
c       change K0,K0ba to K0L and K0S
	do j=1,n
        kf=k(j,2)
        if(kf.eq.311 .or. kf.eq.-311)then
        rrlu=pyr(1)
        k(j,2)=130
        if(rrlu.gt.0.5)k(j,2)=310
        endif
        enddo
c 	P(N,5)=SQRT(MAX(-P(N,1)**2-P(N,2)**2-P(N,3)**2+P(N,4)**2,0.0))
c 	P(N-1,5)=SQRT(MAX(-P(N-1,1)**2-P(N-1,2)**2-P(N-1,3)**2
c     &	+P(N-1,4)**2,0.0))
c	call pyrobo(1,n,0.0,0.0,-bst(1),-bst(2),-bst(3))
c	boost PYJETS back to lab or nucleon-nucleon cms system.

	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sysini(win)   ! 060813
c       give the initial values to quantities needed in calculation
        parameter (KSZ1=30,kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
        common/count/isinel(600)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        common/sa6/kfmaxi,nwhole
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   ! 221203 250204
        anat=nat
        anap=nap
	param(1)=para1_1   ! 250204 200504 110106
c	write(9,*)'para1_1=',param(1)   ! 250204 200504
c	rou0=PARAM(11)
c       considering the nucleus as a sphere with radii rnt for target
c        and rnp for projectile.
c        rnt=(3.*anat/(4.*3.1415926*rou0))**(0.33333)
c        rnp=(3.*anap/(4.*3.1415926*rou0))**(0.33333)
	rp00=1.12   ! 1.05 to 1.12 070613
	rt00=1.12   ! 1.05 to 1.12 070613
c070613	if(nap.gt.16)rp00=1.16*(1-1.16*anap**(-0.666666))
c070613	if(nat.gt.16)rt00=1.16*(1-1.16*anat**(-0.666666))
	if(itden.eq.0)rnt=rt00*anat**(0.33333)   ! 310805
        if(itden.eq.1)rnt=rt00*anat**(0.33333)   ! +0.54 160511
	if(nat.eq.2 .and. nzt.eq.1)rnt=4.0   ! 2.60 2.095  1.54 2603141
c060813 120214 if(itden.eq.2)rnt=0.5
	if(ipden.eq.0)rnp=rp00*anap**(0.33333)   ! 31080
        if(ipden.eq.1)rnp=rp00*anap**(0.33333)   ! +0.54 160511
        if(ipden.ge.2)rnp=0.5   ! lepton   ! 060813 120214
	if(nap.eq.2 .and. nzp.eq.1)rnp=4.0   ! 2.60 2.095  1.54
	rou0=3./4./3.1416*anat/(rnt*rnt*rnt)   ! 310805
        r0p=rnp
        r0t=rnt
C       set initial values to some quantities
c       in the program the x-sections are given in a unit of fm^2   ! 060813.
        csnn=PARAM(1)*0.1
        cspin=PARAM(2)*0.1
        cskn=PARAM(3)*0.1
        cspipi=PARAM(4)*0.1
        cspsn=PARAM(13)*0.1
        cspsm=PARAM(14)*0.1
        csspn=PARAM(15)*0.1
        csspm=PARAM(16)*0.1
c060813 120214
	if(ipden.ge.2)then   
	if(ifram.eq.0)then
	ept=sqrt(win*win+0.938*0.938)
        rots=sqrt((ept+0.938)*(ept+0.938)-win*win)
	endif
        if(ifram.eq.1)rots=win
	call crosep(rots,csen)   ! temporary using e^-p total x-section
c       if(nzp.le.0)call crosep(rots,csen)   ! e^-p total x-section
c       if(nzp.gt.0)call crosepp(rots,csen)   ! e^+p total x-section   
	csen=csen*0.1
	endif   
c060813 120214
c       largest collision distance between two colliding particles.
        edipi=sqrt(cspipi/3.1416)
        epin=sqrt(cspin/3.1416)
        ekn=sqrt(cskn/3.1416)
        ecsnn=sqrt(csnn/3.1416)
	ecspsn=sqrt(cspsn/3.1416)
	ecspsm=sqrt(cspsm/3.1416)
	ecsspn=sqrt(csspn/3.1416)
	ecsspm=sqrt(csspm/3.1416)
	ecsen=sqrt(csen/3.1416)   ! 060813

        anp=nap**.3333
        ant=nat**.3333
        do ia=1,2
        if(ia.eq.1)napt=nap
        if(ia.eq.2)napt=nat
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.ge.27.and.napt.le.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
         elseif(napt.eq.56)then
        alpt=0.49
          elseif(napt.eq.64)then
        alpt=0.49     
           elseif(napt.eq.108)then
        alpt=0.495
           elseif(napt.eq.184)then
        alpt=0.53
           elseif(napt.eq.197)then
        alpt=0.54
           elseif(napt.eq.207)then
        alpt=0.545
            elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	if(ia.eq.1)alp=alpt
  	if(ia.eq.2)alt=alpt
   	enddo
        suppc=rp00*anp+2.*alp
        suptc=rt00*ant+2.*alt
        suppm=1./(1+exp(0.-r0p/alp))
        suptm=1./(1+exp(0.-r0t/alt))

        sig=PARAM(5)*0.1
        rcsit=PARAM(6)
	t0=PARAM(7)
c230805	t0=0.   ! 221102   
c221102	proper formation time of particle from 'pythia'
        dep=PARAM(9)
        ddt=PARAM(8)
        rao=PARAM(10)
        kfmax=KFMAXT
        do i=1,100
        kfaco(i)=KFACOT(i)
        enddo
        do j=1,600
        isinel(j)=ISINELT(j)
        enddo
        do i=1,100
        do j=1,100
        disbe(i,j)=0.
        enddo
        enddo
        do j=1,kfmax
c       something might be missing here ?
        disbe(1,j)=DISDET(j)
        disbe(2,j)=DISDET(j)
        disbe(3,j)=DISDET(j)
        disbe(4,j)=DISDET(j)
c       disbe(26,j)=DISDET(j)
c       disbe(27,j)=DISDET(j)
c       disbe(28,j)=DISDET(j)
c       disbe(29,j)=DISDET(j)
        enddo
400     do i=1,99
        do j=i+1,100
        disbe(j,i)=disbe(i,j)
        enddo
        enddo
c        kfmaxi=kfmax
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine woodsax_samp(ii,jj,alp,r0,am,ac,iway)   ! 191110 230311
c       sample position of nucleon ii in nucleus according to
c        Woods-Saxon distribution 
c       jj=0 and 1 for target and projectile, respectively
c	alp: diffusion length
c	r0: radius of nucleus
c	am: upper bound in sampling the radius
c	ac: maximum radius 
c230311 iway=1: ii must be outside overlap region of colliding nuclei
c230311 iway=0: no more requirement
	PARAMETER (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc  
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	b=bp 
c       selecting sample method
c	Woods-Sax distribution
	iii=0
100	iii=iii+1
        if(iii.eq.100000)then
        write(9,*)'difficult to arrange produced nucleon in'
        write(9,*)'subroutine woodsax,infinitive loop may happen'
	goto 200   ! set larget number of try is equal to 100000
        endif
	a1=pyr(1)
	xf=ac*(a1)**(1./3.)
	b1=pyr(1)
	deno2=1.+exp((xf-r0)/alp)
c       if(deno2.eq.0.)deno2=1.e-10
	yf=1./deno2
c	Gaussan distribution
cc	yf=exp(-xf*xf/2./r0)
	if(b1.gt.yf/am) goto 100
	call samp(xf,ii)
c       subroutine 'samp' is sampling the direction according to isotropic
c        distribution
        x=c17(ii,1)
        y=c17(ii,2)
        z=c17(ii,3)
        if(iway.eq.0)goto 200   ! 230311
c	ii must be outside overlap region of colliding nuclei
        if(jj.eq.0)then   ! ii in target (origin)
c       relative to projectile center, above x, y, and z are b-x, y, and z, 
c	 respectively
c       (b-x,y,z) is inside or not inside the sphere of projectile
        r1=sqrt((b-x)*(b-x)+y*y+z*z)
        if(r1.lt.r0p)goto 100
c        c17(ii,1)=x
c        c17(ii,2)=y
c        c17(ii,3)=z
        endif
        if(jj.eq.1)then   ! ii in projectile
c       relative to target center, they are x+b, y, and z, respectively
c       (x+b,y,z) is inside or not inside the sphere of target
        r1=sqrt((x+b)*(x+b)+y*y+z*z)
        if(r1.lt.r0t)goto 100
c        c17(ii,1)=x
c        c17(ii,2)=y
c        c17(ii,3)=z
	endif
200	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine samp(xf,i)
c       arrange i-th particle on the surface of sphere with radius xf
	PARAMETER (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	cita=2*pyr(1)-1.
	fi=2.*pio*pyr(1)
	sita=sqrt(1.-cita**2)
	c17(i,1)=xf*sita*cos(fi)
	c17(i,2)=xf*sita*sin(fi)
	c17(i,3)=xf*cita
c	if(xf.gt.20) write(*,*) 'xf=',xf
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine filt
c       filter out particles wanted to study and make them in      
c       the order of proton,neutron,pba,nba,pi+,pi-,pi0,k-,k0-,sigma0,
c       sigma-,sigma+,sigma0ba,sigma-ba,sigma+ba,lamda,lamdaba,k0,k+,
c       cascade-,cascade-ba,cascade0,cascade0ba,omega-,omega+,Delta-,
c       Delta0,Delta+,Delta++,rho+,rho-,rho0,J/Psi,Psi',x0c,x1c,x2c,
c       D,Dba,D0,D0ba,lamdac+,sigmac0,sigmac+,sigmac++,omega,k*+,K*0,
c       D*,D*ba,D*0,D*0ba (52 kind of particle altogether) 
c060813	120214 in case of lepton+A, one images lepton as a initial 
c	 projectile proton 
      PARAMETER (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
	common/sa6/kfmaxi,nwhole
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
	SAVE/PYCIDAT2/
	iii=0
	jjj=0
	do i=1,kfmax
	kf=kfaco(i)
	do j=iii+1,n
	call ord(jjj,j,kf)
	enddo
	iii=jjj
	numbs(i)=jjj
	enddo
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ord(ipi,j,kf)
c       order particles according to flavor code
c       j: the particle needed to order
c       ipi: j-th particle should order after ipi
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	parameter(kszj=40000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	dimension kk(5),pp(5),vv(5)
	ik=k(j,2)
	if(ik.eq.kf)then
	ipi=ipi+1
	do jj=1,5
	kk(jj)=k(ipi,jj)
	pp(jj)=p(ipi,jj)
	vv(jj)=v(ipi,jj)
	enddo
	do jj=1,5
	k(ipi,jj)=k(j,jj)
	p(ipi,jj)=p(j,jj)
	v(ipi,jj)=v(j,jj)
	enddo
	do jj=1,5
	k(j,jj)=kk(jj)
	p(j,jj)=pp(jj)
	v(j,jj)=vv(jj)
	enddo
	endif
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine prt_pyj(nn)
c	print particle list and sum of momentum and energy
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/PYJETS/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension peo(4)
	do i=1,nn
	write(mstu(11),*)i,ksa(i,2),(psa(i,j),j=1,4)
	enddo
	call psum(psa,1,nsa,peo)
	ich1=0.
	do i1=1,nn
	kf=ksa(i1,2)
	ich1=ich1+pychge(kf)
	enddo
        write(mstu(11),*)ich1/3,peo   ! 
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbh(nn)
c       print particle list and sum of momentum and energy
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
        do i=1,nn
        write(mstu(11),*)i,kbh(i,2),(pbh(i,j),j=1,4)
c	write(9,*)i,kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        write(mstu(11),*)ich1/3,peo   !
c	write(9,*)ich1/3,peo   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa2(nn)
c       print particle list and sum of momentum and energy
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
        do i=1,nn
        write(9,*)i,kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        write(9,*)ich1/3,peo   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbe(nn)   ! 280809
c       print particle list and sum of momentum and energy
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sbe/nbh,nonbe,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
        do i=1,nn
        write(22,*)i,kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        write(22,*)ich1/3,peo   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_saf(nn)
c       print particle list and sum of momentum and energy
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/saf/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
        do i=1,nn
        write(9,*)i,kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        write(9,*)ich1/3,peo   !
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine psum(pei,il,ih,peo)
c       calculate sum of momentum and energy
c       pei: two dimension array of input momentum and energy
c       il and ih: lower and upper limits of sum
c       peo : one dimension array of output momentum and energy  
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        dimension pei(kszj,5),peo(4)
        do i=1,4
        peo(i)=0.
        enddo
        do i=il,ih
        do j=1,4
        peo(j)=peo(j)+pei(i,j)
        enddo
        enddo
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine scat(time_neu,time_par,time_had,lc,tc,tw,win,parp21
     c	 ,parp22,psno,ijk,ipau,irecon,gamt)   ! 211107 110101
c	administrate a proton-nucleus,nucleus-nucleus,or 
c	 lepton+A collision ! 060813 120214
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      PARAMETER (kszj=40000,mplis=40000,KSZ1=30)
        parameter(nsize=240000)
	dimension b(3),bkk(3)   ! 010600
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c  afl(20,5,2)
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
	common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23 
	common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
 	common/sa16/dtt,dni(10),dpi(10),edi(10),bmin,bmax
     &   ,bar(10),abar(10),barf(10),abarf(10)   
     &   ,emin(10),eminf(10),eplu(10),epluf(10)   
        common/sa18/tdh,itnum,non18,cptl,cptu,cptl2,cptu2,snum(4,20),
     &	 v1(4,20),v2(4,20),v12(4,20),v22(4,20)
	common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph 
     c	 ,vnlep   ! 260314
	common/sa23/kpar,knn,kpp,knp,kep   ! 060813   
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 280809
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21,adiv,gpmax,nnc   !   070417
        common/sa28/nstr,nstr00,nstra(kszj),nstrv(kszj)   ! 280809
	common/sa33/smadel,ecce,secce,parecc,iparres   ! 240412
        common/sa34/iikk   ! 060617
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
	common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5) 
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 240209
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5) ! 240209
        common/ctllist/nctl,noinel(600),nctl0,noel
        common/ctllist_p/nreac(9),nrel   ! 280809
	common/show/vip(mplis),xap(mplis)
        common/schuds/schun,schudn,schudsn,sfra,cmes   !she042021 
        dimension lc(nsize,5),tc(nsize),tw(nsize),tcp0(nsize)   ! 121110
	dimension pi(4),pj(4),pii(4),pjj(4),peo(4),pint(4)	
	dimension nni(10),ndi(10),npi(10)
	dimension pkk(kszj,4),kk6(5),pp6(5),vv6(5)   ! 060803   
	dimension cc(5),nreaco(9)   ! 280809
	dimension kdiq(kszj,5),dgmas(kszj),pl(100,5)   ! 260314
        dimension skapa(5),ksin(kszj,5),psin(kszj,5),vsin(kszj,5)   ! 070417 110517 
        integer winel   	
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 120214
c       arraies in 'pyjets' are used in the processes after calling 'pythia' 
c        in nn or lepton+n collision and after proton-nucleus,nucleus-nucleus,
c        or lepton-nucleus collision 060813 120412
c       arraies in 'sa2' are used in the processes in nn or lepton+n 
c        collision
c       arraies in 'sbh' are used to store hadron after nn or lepton+n 
c        collision and after proton-nucleus,nucleus-nucleus or 
c        lepton-nucleus collision 060813 120414
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 120214
c       numb(i) is used in the scattering processes, numbs(i) is used in 
c        the process of calling 'pythia'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c	0. is the hard distance between two pions
c	0.5 is the hard distance between two nucleons
c	0. is the hard distance between pion and nucleon
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       lc(i,1) and lc(i,2) are the line # of colliding particles of i-th 
c        collision pair in particle list, respectively
c       lc(i,3) and lc(i,4) are the flavor codes of scattered particles
c        in i-th collision.
c       lc(i,5) identifies the different inelastic processes,
c       lc(i,5)=592 refers to the process calling 'pythia'
c       tc(i) is the collision time of i-th colli.
c       tw(i) is the cross section ratio of (i-th inelas.)/tot
c       common block 'sbe' stores cumulatively parton (q,qq, and g) 
c        configuration before breaking the diquarks
c       common block 'saf' stores cumulatively parton (q and g)
c        configuration after breaking the diquarks
c       idi: counts cumunatively the number of diquark (anti-diquark)
c       idio: value of idi after last nn collision
c       ndiq(j): = 0 if j is quark (antiquark)
c                = idi if j is diquark (anti-diquark)
c       note: j is line number in 'sbe' ('saf')
c       npt(idi): = line number of idi-th diquark (anti-diquark)
c        partner in 'saf'
c280809 ifcom(idi): line number of first component of idi-th diquark
c       nstr: statitics of number of strings in a hh collis.
c       nstr00: number of strings after call remo
c       nstra(i): line number of first component of i-th string
c280809 nstrv(i): line number of last component of i-th string
c	write(9,*)'in scat iparres=',iparres   ! 240412
	adj112=adj1(12)   ! 300705
	adj140=adj1(40)   ! 200905
        time=time_neu   ! 111010 121110
        ijk=0
	do i=1,10
	nni(i)=0
	ndi(i)=0
	npi(i)=0
	dni(i)=0.
	dpi(i)=0.
	edi(i)=0.
c033101
	bar(i)=0.
	abar(i)=0.
	barf(i)=0.
	abarf(i)=0.
	emin(i)=0.
	eplu(i)=0.
	eminf(i)=0.
	epluf(i)=0.
c033101
	enddo

csa**********************************************************
        do i1=1,20
        do i2=1,4
        snum(i2,i1)=0.
        v1(i2,i1)=0.
        v2(i2,i1)=0.
        v12(i2,i1)=0.
        v22(i2,i1)=0.
        enddo
        enddo

c240209
        ngam=0
        do i1=1,kszj
        do j1=1,5
        kgam(i1,j1)=0
        pgam(i1,j1)=0.
        vgam(i1,j1)=0.
        enddo
        enddo
c240209
        nctl0=nctl
c060805 mstj(1)=mstj1_1   ! 221203
c070417 loop over hadron-hadron collisions in a nucleus-nucleus collision
c       statistics of number of nucleon-nucleon collisions (nnc) in a nucleus-nucleus
c        collision, statistics of variables of itime etc. over nucleon-nucleon
c        collisions in a nucleus-nucleus collision
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        nnc=0
        sgtime=0.
        sitime=0.
        sastr=0.
        sadiv=0.
        sgpmax=0.
        do i1=1,5
        skapa(i1)=0.
        enddo
        endif
c070417
        iii=1   ! 280809, iii-th hadron-hadron collis.
c	write(9,*)'in scat nctl=',nctl
c	do i1=1,nctl
c	write(9,*)'lc(i1,1),lc(i1,2),tc(i1)=',lc(i1,1),lc(i1,2),tc(i1)
c	enddo
10 	if(iii.eq.1)goto 1000
101	call copl(time)
c	find out the binary colli. with minimum collsion time
1000	call find(icp,tcp,lc,tc,tw,1)
c	write(9,*)'scat af. find icp,tcp=',icp,tcp
	time=tcp   ! 121110
	tcp0(iii)=tcp   ! 121110
c121110	if(icp.eq.0)goto 100
c	icp=0 means the collision list is empty
c121110
	if(icp.eq.0)then
	time_neu=tcp0(iii-1)
	goto 100
	endif
c121110
c050805
	nbe=0
	idi=0
        idio=0
	do i1=1,kszj
        do j1=1,5
	kbe(i1,j1)=0
        pbe(i1,j1)=0.
        vbe(i1,j1)=0.
	enddo
        ndiq(i1)=0
        npt(i1)=0
        ifcom(i1)=0   ! 280809
	enddo
c050805
c070417	mstj(1)=mstj1_1   ! 221203 060805
	l=lc(icp,1)
	l1=lc(icp,2)

	time0=time
	kfa=ksa(l,2)
	kfb=ksa(l1,2)
        ikfa=iabs(kfa)   ! 070417
        ikfb=iabs(kfb)   ! 070417
	kfaab=iabs(kfa)   ! 060813 120214
	kfbab=iabs(kfb)   ! 060813 120214
	time=tcp
c	record this collision time
c	write(9,*)'iii,l,l1,kfa,kfb=',iii,l,l1,kfa,kfb

c?????????????????????????????????????????????????????????????????
csa	if(time.le.ttt)then
c	record the spatial and momentum coordinates
csa	if(iiii.eq.1)call txp(time0,time)
c	calculate the directed and elliptic flow
csa	call flow(time0,time)
csa	endif
csa***************************************************************

cc	tlco(l,4)=tcp
cc	tlco(l1,4)=tcp
20	continue
	ilo=0
	pi(4)=psa(l,4)
	pj(4)=psa(l1,4)
	if(pi(4).lt.1.e-20)pi(4)=1.e-20   ! 041204
	if(pj(4).lt.1.e-20)pj(4)=1.e-20   ! 041204
	do i=1,3
	pi(i)=psa(l,i)
	pj(i)=psa(l1,i)
c	if(pi(4).lt.1.e-5.or.pj(4).lt.1.e-5)then
c    	 write(*,*)'pi,pj,n,l,l1,icp=',
c     &	pi(4),pj(4),nsa,l,l1,icp,nctl,tcp,tc(icp)
c	do iop=1,nctl
c	write(*,*)lc(iop,1),lc(iop,2),tc(iop)
c	enddo
c	endif
	b(i)=(pi(i)+pj(i))/(pi(4)+pj(4))
	enddo
c200601
        pti=sqrt(pi(1)**2+pi(2)**2)
        ptj=sqrt(pj(1)**2+pj(2)**2)
c200601
c	boost to CMS frame of colliding pair
	call lorntz(ilo,b,pi,pj)
	ss=pi(4)+pj(4)
        if(ss.lt.1.e-18)ss=1.e-18
c	perform classical Newton motion   
	call his(time,lc,tc,tw,istop)
c	write(9,*)'scat af.his ss,nctl,istop=',ss,nctl,istop
c 	do iop=1,nctl
c	write(9,*)lc(iop,1),lc(iop,2),tc(iop)
c	enddo
	if(istop.eq.1)goto 100
c	istop=1 means all particles have get out of considered volume
	m1=numb(1)
        m2=numb(2)
	m3=numb(3)
        m4=numb(4)
	m7=numb(7)
c	write(9,*)'m1 - m7=',m1,m2,m3,m4,m7
	if((ipden.lt.2.and.(l.le.m2.and.l1.le.m2).and.ss.ge.parp21)
     c   .or.(ipden.gt.2.and.(((kfaab.ge.11.and.kfaab.le.16).and.l1.le.
     c   m2).or.(l.le.m2.and.(kfbab.ge.11.and.kfbab.le.16))).and.ss.ge.
     c   parp21))then   ! if 1 011210 060813 120214 
c       calculate the angular 'seta' of the momenta pi and pj
	ctai=pyangl(pi(3),sqrt(pi(1)**2+pi(2)**2))
	ctaj=pyangl(pj(3),sqrt(pj(1)**2+pj(2)**2))
	cctai=cos(ctai)
	cctaj=cos(ctaj)
	if(cctai.gt.0.)then
c       calculate the 'orentation' of the vector pi
	call codi(pi,cfi1,sfi1,ccta1,scta1)
	else
	call codi(pj,cfi1,sfi1,ccta1,scta1)
	endif

	if(kfa.eq.2212.and.kfb.eq.2212)then   
c200601
        kpp=kpp+1
        if(pti.le.1.e-4)kpar=kpar+1
        if(ptj.le.1.e-4)kpar=kpar+1
c200601
	call pyinit('cms','p','p',ss)   
c	write(9,*)'after calling pyeinit'   !s
        call pyevnt
c	write(9,*)'after calling pyevnt'   !s
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.2212.and.kfb.eq.2112)then   ! 070417
c200601
        knp=knp+1
        if(pti.le.1.e-4)kpar=kpar+1
        if(ptj.le.1.e-4)kpar=kpar+1
c200601
	if(cctai.gt.0.)then
	call pyinit('cms','p','n0',ss)   
	else
	call pyinit('cms','n0','p',ss)   
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.2212.and.kfa.eq.2112)then    ! 070417
c200601
        knp=knp+1
        if(pti.le.1.e-4)kpar=kpar+1
        if(ptj.le.1.e-4)kpar=kpar+1
c200601
	if(cctai.gt.0.)then
	call pyinit('cms','n0','p',ss)   
	else
	call pyinit('cms','p','n0',ss)   
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.2112.and.kfb.eq.2112)then    ! 070417
c200601
        knn=knn+1
        if(pti.le.1.e-4)kpar=kpar+1
        if(ptj.le.1.e-4)kpar=kpar+1
c200601
	call pyinit('cms','n0','n0',ss)   
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417
c060813
	if(kfa.eq.11.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
c	call pyinit('cms','gamma/e-','p',ss)
	call pyinit('cms','e-','p',ss)
	else
c	call pyinit('cms','p','gamma/e-',ss)
	call pyinit('cms','p','e-',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.11.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
c	call pyinit('cms','p','gamma/e-',ss)
	call pyinit('cms','p','e-',ss)
	else
c	call pyinit('cms','gamma/e-','p',ss)
	call pyinit('cms','e-','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.11.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
c	call pyinit('cms','gamma/e-','n0',ss)
	call pyinit('cms','e-','n0',ss)
	else 
c	call pyinit('cms','n0','gamma/e-',ss)
        call pyinit('cms','n0','e-',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.11.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
c	call pyinit('cms','n0','gamma/e-',ss)
        call pyinit('cms','n0','e-',ss)
	else
c	call pyinit('cms','gamma/e-','n0',ss)
	call pyinit('cms','e-','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417
c120214
	if(kfa.eq.-11.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
c	call pyinit('cms','gamma/e+','p',ss)
	call pyinit('cms','e+','p',ss)
	else
c	call pyinit('cms','p','gamma/e+',ss)
	call pyinit('cms','p','e+',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-11.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
c	call pyinit('cms','p','gamma/e+',ss)
	call pyinit('cms','p','e+',ss)
	else
c	call pyinit('cms','gamma/e+','p',ss)
	call pyinit('cms','e+','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.-11.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
c	call pyinit('cms','gamma/e+','n0',ss)
	call pyinit('cms','e+','n0',ss)
	else 
c	call pyinit('cms','n0','gamma/e+',ss)
        call pyinit('cms','n0','e+',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-11.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
c	call pyinit('cms','n0','gamma/e+',ss)
        call pyinit('cms','n0','e+',ss)
	else
c	call pyinit('cms','gamma/e+','n0',ss)
	call pyinit('cms','e+','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.12.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_e','p',ss)
	else
	call pyinit('cms','p','nu_e',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.12.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','nu_e',ss)
	else
	call pyinit('cms','nu_e','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.12.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_e','n0',ss)
	else 
        call pyinit('cms','n0','nu_e',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.12.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','nu_e',ss)
	else
	call pyinit('cms','nu_e','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.-12.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_ebar','p',ss)
	else
	call pyinit('cms','p','nu_ebar',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-12.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','nu_ebar',ss)
	else
	call pyinit('cms','nu_ebar','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.-12.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_ebar','n0',ss)
	else 
        call pyinit('cms','n0','nu_ebar',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-12.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','nu_ebar',ss)
	else
	call pyinit('cms','nu_ebar','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.13.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','mu-','p',ss)
	else
	call pyinit('cms','p','mu-',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.13.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','mu-',ss)
	else
	call pyinit('cms','mu-','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.13.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','mu-','n0',ss)
	else 
        call pyinit('cms','n0','mu-',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.13.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','mu-',ss)
	else
	call pyinit('cms','mu-','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.-13.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','mu+','p',ss)
	else
	call pyinit('cms','p','mu+',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-13.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','mu+',ss)
	else
	call pyinit('cms','mu+','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.-13.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','mu+','n0',ss)
	else 
        call pyinit('cms','n0','mu+',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-13.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','mu+',ss)
	else
	call pyinit('cms','mu+','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.14.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_mu','p',ss)
	else
	call pyinit('cms','p','nu_mu',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.14.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','nu_mu',ss)
	else
	call pyinit('cms','nu_mu','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.14.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_mu','n0',ss)
	else 
        call pyinit('cms','n0','nu_mu',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.14.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','nu_mu',ss)
	else
	call pyinit('cms','nu_mu','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.-14.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_mubar','p',ss)
	else
	call pyinit('cms','p','nu_mubar',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-14.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','nu_mubar',ss)
	else
	call pyinit('cms','nu_mubar','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.-14.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_mubar','n0',ss)
	else 
        call pyinit('cms','n0','nu_mubar',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-14.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','nu_mubar',ss)
	else
	call pyinit('cms','nu_mubar','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.15.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','tau-','p',ss)
	else
	call pyinit('cms','p','tau-',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.15.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','tau-',ss)
	else
	call pyinit('cms','tau-','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.15.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','tau-','n0',ss)
	else 
        call pyinit('cms','n0','tau-',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.15.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','tau-',ss)
	else
	call pyinit('cms','tau-','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.-15.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','tau+','p',ss)
	else
	call pyinit('cms','p','tau+',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-15.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','tau+',ss)
	else
	call pyinit('cms','tau+','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.-15.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','tau+','n0',ss)
	else 
        call pyinit('cms','n0','tau+',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-15.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','tau+',ss)
	else
	call pyinit('cms','tau+','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.16.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_tau','p',ss)
	else
	call pyinit('cms','p','nu_tau',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.16.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','nu_tau',ss)
	else
	call pyinit('cms','nu_tau','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.16.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_tau','n0',ss)
	else 
        call pyinit('cms','n0','nu_tau',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.16.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','nu_tau',ss)
	else
	call pyinit('cms','nu_tau','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfa.eq.-16.and.kfb.eq.2212)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_taubar','p',ss)
	else
	call pyinit('cms','p','nu_taubar',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-16.and.kfa.eq.2212)then   ! 070417
	kep=kep+1
        if(cctai.gt.0.)then
	call pyinit('cms','p','nu_taubar',ss)
	else
	call pyinit('cms','nu_taubar','p',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

        if(kfa.eq.-16.and.kfb.eq.2112)then   ! 070417
        kep=kep+1
	if(cctai.gt.0.)then
	call pyinit('cms','nu_taubar','n0',ss)
	else 
        call pyinit('cms','n0','nu_taubar',ss)
	endif
        call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

	if(kfb.eq.-16.and.kfa.eq.2112)then   ! 070417
	kep=kep+1
	if(cctai.gt.0.)then
        call pyinit('cms','n0','nu_taubar',ss)
	else
	call pyinit('cms','nu_taubar','n0',ss)
	endif
	call pyevnt
	goto 2222   ! 070417
	endif   ! 070417

c120214
c060813

2222	if(ipden.lt.11)call pyedit(2)   ! 060813 070417
        if(ipden.ge.11)call pyedit(1)   ! 060813

500     continue
c241003
c	call pyedit(2)   ! sa
c	call psum(p,1,n,peo)   !!!
c	write(9,*)'before, rot., sum=',n,peo   ! sa
c	write(mstu(11),*)'in scat'
c	call pylist(1)
c	do i1=1,n
c	write(9,*)'p=',(p(i1,j1),j1=1,4)
c	enddo
c241003
	do j=1,n
	do j1=1,4
	pint(j1)=p(j,j1)
	enddo
c     c	 write(9,*)'before pint=',(pint(i1),i1=1,3)   ! sa

c	if(cctai.gt.0.99)goto 1002

c	'cctai.gt.0.99' means pi (or pj) nearly on the z axis, don't need 
c	  rotation
c	perform the rotate for produced particle from calling 'pythia'
	call rosa(cfi1,sfi1,ccta1,scta1,cfis,sfis,cctas,sctas,pint)
c	if(kf.eq.443)then
c	write(9,*)'cfis,sfis,cctas,sctas=',cfis,sfis,cctas,sctas   !!!
c	endif   !!!
	do j1=1,4
	p(j,j1)=pint(j1)
	enddo
	enddo

c241003
c	do i1=1,n
c	write(9,*)'p=',(p(i1,j1),j1=1,4)
c	enddo
c241003

1002	continue
c	boost back to Lab.
c	write(9,*)'before Lorentz back,b=',b
c	call psum(p,1,n,peo)   ! sa
c	write(9,*)'sum=',peo   ! sa
	ilo=1
	do j=1,n,2
	if(j.eq.n)then
	do j1=1,4
	pi(j1)=p(j,j1)
	enddo
	call lorntz(ilo,b,pi,pi)
	do j1=1,4
	p(j,j1)=pi(j1)
	enddo
	goto 510
	endif
	do j1=1,4
	pi(j1)=p(j,j1)
	pj(j1)=p(j+1,j1)
	enddo
	call lorntz(ilo,b,pi,pj)
	do j1=1,4
	p(j,j1)=pi(j1)
	p(j+1,j1)=pj(j1)
	enddo
510	enddo
c110517 rotation and boost for particles in sbh
	do j=1,nbh
	do j1=1,4
	pint(j1)=pbh(j,j1)
	enddo

c	'cctai.gt.0.99' means pi (or pj) nearly on the z axis, don't need 
c	  rotation
c	perform the rotate for produced particle from calling 'pythia'
	call rosa(cfi1,sfi1,ccta1,scta1,cfis,sfis,cctas,sctas,pint)
c	if(kf.eq.443)then
c	write(9,*)'cfis,sfis,cctas,sctas=',cfis,sfis,cctas,sctas   !!!
c	endif   !!!
	do j1=1,4
	pbh(j,j1)=pint(j1)
	enddo
	enddo

c	boost back to Lab.
	ilo=1
	do j=1,nbh,2
	if(j.eq.nbh)then
	do j1=1,4
	pi(j1)=pbh(j,j1)
	enddo
	call lorntz(ilo,b,pi,pi)
	do j1=1,4
	pbh(j,j1)=pi(j1)
	enddo
	goto 511
	endif
	do j1=1,4
	pi(j1)=pbh(j,j1)
	pj(j1)=pbh(j+1,j1)
	enddo
	call lorntz(ilo,b,pi,pj)
	do j1=1,4
	pbh(j,j1)=pi(j1)
	pbh(j+1,j1)=pj(j1)
	enddo
511	enddo
c110517
c	write(9,*)'after Lorentz back'
c	call psum(p,1,n,peo)   ! sa
c	write(9,*)'sum=',peo   ! sa

c-------------------------------------------------------------------
c       give four position to the particles after calling pyevnt
c110517 for particles in pyjets
        call ptcre(l,l1,time,1)
c       for particles in sbh
        call ptcre(l,l1,time,2)
c110517
c       arrange particles (quark,diquark, and gluon mainly) after
c        calling pyevnt into the overlap region randomly
c061207	call ptcre(l,l1,time,gamt)
c--------------------------------------------------------------------
	icp5=592
c	592-th scattering process is referred to calling 'pythia'
        noinel(icp5)=noinel(icp5)+1
c	write(9,*)'# of calling pythia=',noinel(icp5)   ! 070802
c260314	statistics of number of leptons studied, identify scattered lepton,  
c	 and fill up pscal(5) 
	if(ipden.ge.11.and.ipden.le.16)then   !
c	identify the studied leptons
	kfl=ipden
	if(nzp.gt.0.)kfl=-ipden
	nlep=0
	do j=1,n
	ikl=k(j,2)
	if(ikl.eq.kfl)then
	nlep=nlep+1
	pl(nlep,1)=p(j,1)
	pl(nlep,2)=p(j,2)
	pl(nlep,3)=p(j,3)
	pl(nlep,4)=p(j,4)
	pl(nlep,5)=p(j,5)
	endif
	enddo
c	find the scattered lepton (with largest energy among studied leptons)
	if(nlep.gt.1)then   !!
	vnlep=vnlep+nlep
	elep=1.d0
	do j1=1,nlep
	plj14=pl(j1,4)
	if(plj14.ge.elep)then
	elep=plj14
	jj=j1
	endif
	enddo
	do j2=1,5
	pscal(j2)=pl(jj,j2)
	enddo
	elseif(nlep.eq.1)then   !!
	vnlep=vnlep+nlep
	do j2=1,5
        pscal(j2)=pl(nlep,j2)
        enddo
	else   !!
	endif   !!
c	write(9,*)'pscal=',(pscal(i1),i1=1,5)
c	calculate kinematic variables relevant to incident and scattered 
c	 lepton only, in cms
	pdotk=pinch(4)*pincl(4)-pinch(1)*pincl(1)-pinch(2)*pincl(2)   
     c   -pinch(3)*pincl(3)   ! P.k  
	q11=pincl(1)-pscal(1)
	q22=pincl(2)-pscal(2)
	q33=pincl(3)-pscal(3)
	q44=pincl(4)-pscal(4)
	q112=q11*q11
	q222=q22*q22
	q332=q33*q33
	q442=q44*q44
	pdotq=pinch(4)*q44-pinch(1)*q11-pinch(2)*q22-pinch(3)*q33   ! P.q
	vnu=pdotq/pinch(5)   ! \nu
	fq2=-(q442-q112-q222-q332)   ! Q^2=-q^2
	w2l=(pinch(4)+q44)**2-(pinch(1)+q11)**2-(pinch(2)+q22)**2-
     c	 (pinch(3)+q33)**2   ! W^2
	pdotk=dmax1(pdotk,1.d-20)
	yyl=pdotq/pdotk   ! y
	pdotq=dmax1(pdotq,1.d-20)
	xb=fq2/2./pdotq   ! x_b
c	write(9,*)'q11,q22,q33,q44=',q11,q22,q33,q44
c	write(9,*)'q112,q222,q332,q442=',q112,q222,q332,q442
c	write(9,*)'pdotk,pdotq,vnu=',pdotk,pdotq,vnu
c	write(9,*)'fq2,w2l,yyl,xb',fq2,w2l,yyl,xb
	endif   !
c260314
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
c060617
        if(iikk.eq.1)then
        kkii=kkii+1
        nbh=0
        do i1=1,kszj
        do j1=1,5
        kbh(i1,j1)=0
        pbh(i1,j1)=0.
        vbh(i1,j1)=0.
        enddo
        enddo
        endif
c060617
c110517
        if(iikk.eq.0)then   ! 060617
c	pyjets to 'sin' which is a internal array
	nsin=n
	do i1=1,n
        do i2=1,5
        ksin(i1,i2)=k(i1,i2)
        psin(i1,i2)=p(i1,i2)
        vsin(i1,i2)=v(i1,i2)
        enddo
	enddo
c	sbh to pyjets
	n=nbh
        do i1=1,n
        do i2=1,5
        k(i1,i2)=kbh(i1,i2)
        p(i1,i2)=pbh(i1,i2)
        v(i1,i2)=vbh(i1,i2)
        enddo
        enddo
        do i1=n+1,kszj
        do i2=1,5
	k(i1,i2)=0
	p(i1,i2)=0.
	v(i1,i2)=0.
	enddo
	enddo
c	'sin' to sbh
	nbh=nsin
        do i1=1,nbh
        do i2=1,5
        kbh(i1,i2)=ksin(i1,i2)
        pbh(i1,i2)=psin(i1,i2)
        vbh(i1,i2)=vsin(i1,i2)
        enddo
        enddo
        do i1=nbh+1,kszj
        do i2=1,5
        kbh(i1,i2)=0
        pbh(i1,i2)=0.
        vbh(i1,i2)=0.
        enddo
        enddo
c110517
        endif   ! 060617
	if(adj140.eq.5) goto 30000   ! 140414
c070802
	igq=0
	do j1=1,n
	kfj1=iabs(k(j1,2))
	if(kfj1.le.8.or.kfj1.eq.2101.or.kfj1.eq.3101.or.kfj1.eq.3201
     c	 .or.kfj1.eq.1103.or.kfj1.eq.2103.or.kfj1.eq.2203.or.kfj1.eq.
     c	 3103.or.kfj1.eq.3203.or.kfj1.eq.3303.or.kfj1.eq.21)igq=igq+1! 040805'
	enddo
cs	write(9,*)'af. ptcre iii,igq,mstj(1)=',iii,igq,mstj(1)   ! sa
	if(igq.eq.0)goto 9003   ! no q, diquark, and g at all 130206
	if(igq.eq.0)then   ! no q, diquark, and g at all
c	remove current nn collision pair from collision list
	do j1=icp+1,nctl
	j=j1-1
	tc(j)=tc(j1)
        tw(j)=tw(j1)
        do m=1,5
        lc(j,m)=lc(j1,m)
        enddo
	enddo
	nctl=nctl-1
	iii=iii+1   ! 060805
	goto 10
	endif
cs	write(22,*)'be. remo iii,n,nbh=',iii,n,nbh   ! sa
c     c	 ,iiii   ! sa 
cs	call pylist(1)
c       remove hadrons from 'PYJETS' to 'sbh' and truncate 'PYJETS'
c	 correspondingly
c110517	call remo
c280809
c070417 find number of strings and line number of first and last components of each
c        string as well as consider the reduction of strange quark suppression
        nstr=0

c070417
        if(iikk.eq.0 .and. (kjp22.eq.0 .or. kjp22.eq.1))then   ! 060617
	nnc=nnc+1
        itime=0
        gtime=0.
        astr=0.
        adiv=0.
        gpmax=0.
        do i1=1,5
        akapa(i1)=0.
        enddo
        vfr24=3.5   ! parameter alpha
        vfr25=0.8   ! \sqrt(s_0) in GeV
	vfr252=vfr25*vfr25
        
        jb=0
10000   do i1=jb+1,n
c	find a string
        if(k(i1,1).eq.2)then   ! i1 is 'A'
        do i2=i1+1,n
        if(k(i2,1).eq.1)then   ! i2 is 'V'
        nstr=nstr+1
        nstra(nstr)=i1   ! line number of first component of nstr-th string
        nstrv(nstr)=i2   ! line number of last component of nstr-th string

c	calculate the effective string tension and parj(1) etc. for current string
        toteng=0.0
        toten=0.0
        totglukt=0.0
        pmax=0.d0
        ggg=0.
	do i3=i1,i2
        toten=toten+p(i3,4)   ! root_s, string total energy 
        pp2=p(i3,1)**2+p(i3,2)**2   ! pt*pt 
        ppp=dsqrt(pp2)
	if(k(i3,2).eq.21)then  
        if(ppp.gt.pmax)pmax=ppp   ! k_{Tmax} 
        if(pp2.ge.vfr252)then   
        toteng=toteng+dlog(pp2/vfr252)   ! sum over gluons in a string 
        ggg=ggg+1.
        endif
	endif
        enddo
        pmax2=pmax*pmax   
        if(pmax2.ge.vfr252)totglukt=totglukt+dlog(pmax2/vfr252)   ! numerator 
        sss=dlog(toten*toten/vfr252)+toteng   ! denominator   
        div=totglukt/sss
c       div: factor related to number of gluons and hardest gluon in current string
c       pmax: transverse momentum of hardest gluon in current string
        adiv=adiv+div
        gpmax=gpmax+pmax
c       string tension of the pure qqbar string, kapa0, is assumed to be 1
        effk2=(1.-div)**(-vfr24)
        akapa(1)=akapa(1)+effk2
        itime=itime+1
        gtime=gtime+ggg
	if(kjp22.eq.0)goto 10001
        akapa(2)=akapa(2)+parj2**(1./effk2)
        akapa(3)=akapa(3)+parj21*((effk2/1.)**(0.5))
        akapa(4)=akapa(4)+parj1**(1./effk2)
        akapa(5)=akapa(5)+parj3**(1./effk2)
		
10001	jb=i2   ! 160317
        if(jb.lt.n)goto 10000
        if(jb.eq.n)goto 20000
        endif
        enddo
        endif
        enddo
20000   nstr00=nstr

        do i1=1,n
        if(k(i1,2).eq.92)astr=astr+1.
        enddo
c	average over strings in a nucleon-nucleon collision
        atime=dfloat(itime)
        if(atime.gt.0.)then
        akapa(1)=akapa(1)/atime
        akapa(2)=akapa(2)/atime
        akapa(3)=akapa(3)/atime
        akapa(4)=akapa(4)/atime
        akapa(5)=akapa(5)/atime
        gtime=gtime/atime
        adiv=adiv/atime
        gpmax=gpmax/atime
c       gtime: averaged # of gluons in a string in current nucleon-nucleon (NN) collision
        endif
c	statistics of variables of itime etc. over NN collisions in a nucleus-nucleus 
c	 collision 
        sadiv=sadiv+adiv
        sgpmax=sgpmax+gpmax
        skapa(1)=skapa(1)+akapa(1)
        skapa(2)=skapa(2)+akapa(2)
        skapa(3)=skapa(3)+akapa(3)
        skapa(4)=skapa(4)+akapa(4)
        skapa(5)=skapa(5)+akapa(5)
        sgtime=sgtime+gtime
        sastr=sastr+astr
        sitime=sitime+itime
        endif
c070417

        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813
c        if(iii.eq.12.or.iii.eq.29.or.iii.eq.151.or.iii.eq.213)then
c        write(22,*)'nstr00,nstr,nstra,nstrv=',nstr00,nstr
c        write(22,*)(nstra(i1),i1=1,nstr)
c        write(22,*)(nstrv(i1),i1=1,nstr)
c        endif
c280809
c080104
c	'pyjets' to 'sbe'. etc.
	if(n.ge.1)then   ! 1
	do i1=1,n
c050805	i3=i1+nbe
	i3=i1   ! 050805
	kf=k(i1,2)
        kfab=iabs(kf)
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then   ! 2
c     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then   ! 2
        idi=idi+1
c050805	ndiq(i1+naf)=idi
	ndiq(i1)=idi   ! 050805
	endif   ! 2
        do i2=1,5
        kbe(i3,i2)=k(i1,i2)
        pbe(i3,i2)=p(i1,i2)
        vbe(i3,i2)=v(i1,i2)
        enddo
	enddo
c050805	nbeo=nbe   ! 190204
	nbeo=0
	nbe=i3
	endif   ! 1
cs	write(9,*)'af. sbe n,nbh,nbe=',n,nbh,nbe   ! sa
c	goto 200
c	write(22,*)'af. sbe iii,nbe=',iii,nbe   ! 240412
c	call pylist(1)   ! 240412
c	call prt_sbe(nbe)   ! 240412
c080104
c       break up diquark and give four momentum and four position
c        to the broken quarks (working in 'pyjets')
	call break
        n00=n   ! 280809
c280809 n00: 'largest line number' in 'pyjets' in iii-th nucleon-nucleon colli.
c        after 'break'. partons above n00 are all produced partons from
c        inelastic collis.
c280809
c101210
c       if(iii.eq.141)then
c       write(22,*)'af break'
c       call pyedit(2)
c       call pylist(1)
c       call prt_sbh(nbh)  
c       endif
c101210
	if(adj1(18).eq.0)goto 777   ! without Pauli
c191202
c       Pauli effect (working in 'pyjets'_current and 'saf'_old)
        tpaul=1.
c       tpaul: product of the unoccupation probabilities
        do i1=1,n   ! current
        kfp=k(i1,2)
        kfp=iabs(kfp)
        ppaul=1.
        if(kfp.eq.1 .or. kfp.eq.2 .or. kfp.eq.3)then
        
	call pauli(i1,ppaul)
c       ppaul: the unoccupation probability of particle i1
        if(ppaul .lt. 0.)then   ! over occupation, should be blocked
        tpaul=0.
        goto 666
        endif
        endif
        tpaul=tpaul*ppaul
        enddo
666     if(pyr(1) .ge. tpaul)then   ! blocked "1"
        noel=noel+1   ! statistics of number of collisions blocked
c080104 
	ipau=1
c190204
c	remove "current part of 'sbe'" from 'sbe' and truncate 'ndiq' and 
c	 'npt' correspondingly
c050805	do i1=nbeo+1,nbe
        do i1=1,nbe   ! 050805
        do i2=1,5
        kbe(i1,i2)=0
        pbe(i1,i2)=0.
        vbe(i1,i2)=0.
        enddo
c050805	ndiq(i1+naf)=0
	ndiq(i1)=0   ! 050805
        enddo
c050805	do i1=idio,idi
	do i1=1,idi   ! 050805
	npt(i1)=0
        ifcom(i1)=0   ! 280809
	enddo
c050805	nbe=nbeo
	nbe=0   ! 050805
c190204
c050805	idi=idio
	idi=0   ! 050805
c080104
c	remove current nn collision pair from collision list
	do j1=icp+1,nctl
	j=j1-1
	tc(j)=tc(j1)
        tw(j)=tw(j1)
        do m=1,5
        lc(j,m)=lc(j1,m)
        enddo
	enddo
	nctl=nctl-1
	goto 10
        endif   ! "1"
777	continue   ! unblocked
c080104
c       'pyjets' to 'saf'   
c	do i1=1,n
c	naf=naf+1
c	if(naf.gt.kszj)then
c	write(9,*)'iiii,naf,kszj=',iiii,naf,kszj   ! sa
c	stop 11111
c	endif
c	do i2=1,5
c	kaf(naf,i2)=k(i1,i2)
c	paf(naf,i2)=p(i1,i2)
c	vaf(naf,i2)=v(i1,i2)
c	enddo
c	enddo
c080104
200	continue
c050805	idio=idi   ! 080104
	idio=0   ! 080104 050805
c080104
c300705
c       add CME charge separation for u d s,c   !she042021
c        print*,"n=",n
        if(cmes.eq.0)goto 902
        if((nap.eq.nat).and.(nzp.eq.nzt).and.(cmes.eq.1))
     c    call chargecme(win)
c
902     if(adj140.eq.1)goto 9002   ! run ended before parton cascade
c	'pyjets' to 'parlist'
        iprl=n
        do i=1,n
        idp(i)=k(i,2)
        rp(4,i)=0.
c       parton cascade process is assumed to be start at time 0.
        eee=p(i,4)
        pp(4,i)=eee
        do j=1,3
        rp(j,i)=v(i,j)
        ppp=p(i,j)
        pp(j,i)=ppp
        vp(j,i)=ppp/eee
        enddo
        rmp(i)=p(i,5)
        taup(i)=0.
        vip(i)=0.
        xap(i)=0.
        enddo
        do i=n+1,mplis
        do j=1,3
        rp(j,i)=0.
        pp(j,i)=0.
        vp(j,i)=0.
        enddo
        rp(4,i)=0.
        pp(4,i)=0.
        taup(i)=0.
	vip(i)=0.
        xap(i)=0.
        idp(i)=0
        rmp(i)=0.
        enddo
c201203
c	goto 890
c	write(22,*)'be. call parcas n,iprl=',n,iprl   ! 240412
c	call prt_pyj(n)   ! 240412 
	iijk=0   ! 151203
c280809
        do i1=1,9
        nreaco(i1)=nreac(i1)
        enddo
c280809
	call parcas(time_par,iii,iijk,win,nap,rnt,rnp,n00)   ! 120603 280809
c120603
	if(iijk.eq.1)then
        iiii=iiii-1   ! 10/08/98
        ijk=1  ! 10/08/98
        return   ! 10/08/98
	endif
	if(iijk.eq.2)siijk=siijk+1   ! 201203
c201203
c       'parlist' to 'pyjets'
        n=iprl
c121110	if(n.eq.0)goto 300   ! no parton at all, give up this event
        do i=1,n
c280809 k(i,1)=1
        k(i,2)=idp(i)
c280809 do j=3,5
c       k(i,j)=0
c280809 enddo
        do j=1,4
        v(i,j)=rp(j,i)
        p(i,j)=pp(j,i)
        enddo
        p(i,5)=rmp(i)
        v(i,5)=0.
        enddo
c	do i=n+1,kszj
c	do j=1,5
c	k(i,j)=0
c	p(i,j)=0.
c	v(i,j)=0.
c	enddo
c	enddo
c	write(22,*)'af. call parcas n,iprl=',n,iprl   ! 240412
c	call prt_pyj(n)   ! 240412
c101210
c       if(iii.eq.141)then   
c       write(9,*)'af parcas'
c       write(22,*)'af parcas'
c	call pyedit(2)
c       call pylist(1)
c       call prt_sbh(nbh)   
c       endif
c101210
c201203
c121110 
        n55=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=55
        n55=n55+1
        endif
        enddo
c       move "55" from 'pyjets' to 'sgam'
        if(n55.gt.0)call remo_gam(55)
c121110
	if(adj140.eq.2)goto 9002   ! 200905 run ended after parton cascade  

c280809 if(adj112.ne.0)goto 889   ! coalescence
        if(adj112.ne.0.or.(adj112.eq.0.and.(nreac(4).gt.nreaco(4).or.
     c   nreac(6).gt.nreaco(6).or.nreac(7).gt.nreaco(7))))goto 889   
c	 coalescence ! 280809 150512
c	recover parton configuration in 'sbe' 
c280809
        if(idi.gt.0)then
c       if(iii.eq.12.or.iii.eq.29.or.iii.eq.151.or.iii.eq.213)then
c       write(22,*)'be. recover n,nbe,n00,idi=',n,nbe,n00,idi
c       call prt_sbe(nbe)
c       endif
c280809
c	loop over 'sbe'
	idii=0
cs	write(9,*)'be. recover n,nbe=',n,nbe   ! sa
	do i=1,nbe
	kf=kbe(i,2)
        kfab=iabs(kf)
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then
c060805     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then
	idii=idii+1
	do j=1,5
	kdiq(idii,j)=kbe(i,j)   ! diquark k array
	enddo
	dgmas(idii)=pbe(i,5)   ! diquark mass
c	write(9,*)'i,kfab,idii,kdiq=',i,kfab,idii,(kdiq(idii,j),j=1,5)
	endif
	enddo
c	loop over 'pyjets'
c	write(9,*)'pyjets n=',n
	idij=0
	jb=0
	dele=0.
880	do 980 i=jb+1,n
	jb=jb+1
	ndiqi=ndiq(i)
	if(ndiqi.ne.0)then   ! diquark (anti-diquark)
	idij=idij+1
	j=npt(ndiqi)   ! diquark partner 
	do i1=1,5
        k(i,i1)=kdiq(idij,i1)   ! diquark k array
	enddo
c	write(9,*)'i,p=',i,(p(i,i1),i1=1,5)
c	write(9,*)'j,p=',j,(p(j,i1),i1=1,5)
	do i1=1,3   
	p(i,i1)=p(i,i1)+p(j,i1)
	enddo
	dimass=dgmas(idij)   ! diquark mass
	pi1=p(i,1)
	pi2=p(i,2)
	pi3=p(i,3)
	pi4=sqrt(pi1*pi1+pi2*pi2+pi3*pi3+dimass*dimass)
	dele=dele+p(i,4)+p(j,4)-pi4
c	write(9,*)'dele,p(i,4),p(j,4),pi4=',dele,p(i,4),p(j,4),pi4
	p(i,4)=pi4
	p(i,5)=dimass
c	write(9,*)'ndiqi,pi1,pi2,pi3,pi4,dimass=',ndiqi,pi1,pi2,
c     c	pi3,pi4,dimass
c060805
c	if(j.eq.n)goto 1800
	if(j.eq.n)then
	n=n-1
	goto 1800
	endif
c060805
	goto 1100
	endif
980	continue
	goto 1800
1100	continue
c       move particle list,'pyjets',and 'ndiq'one step downward from 
c	 j+1 to n
        do j1=j+1,n
	ndiq(j1-1)=ndiq(j1)
        do jj=1,5
        k(j1-1,jj)=k(j1,jj)
        p(j1-1,jj)=p(j1,jj)
        v(j1-1,jj)=v(j1,jj)
        enddo
        enddo
        n=n-1
c	subtract 'npt' by one from idij+1 to idi
	if(idij.lt.idi)then
	do j1=idij+1,idi
	npt(j1)=npt(j1)-1
	enddo
	endif
        goto 880
1800	continue
c060805	n=jb
cs	write(22,*)'af recover 1 jb,n=',jb,n
cs	call prt_pyj(n)   ! sa
cs	call prt_sbh(nbh)

c       share energy 'dele' into particles
c	write(9,*)'n,dele=',n,dele
	del=dele/float(n)
	do j3=1,n
        p(j3,4)=p(j3,4)+del
        if(del.lt.0.)then
        if(p(j3,4).lt.0.)p(j3,4)=p(j3,4)-del
        pabs=abs(p(j3,3))
        if(pabs.ge.p(j3,4))p(j3,4)=p(j3,4)-del
        endif
        enddo
c	write(22,*)'af. recover 2 iii,n=',iii,n   !240412
c	call pylist(1)   ! 240412
c	call prt_pyj(n)   ! 240412   
cs	call prt_sbh(nbh)

        if(iparres.eq.0 .or. (iparres.eq.1.and.(nreaco(4).eq.nreac(4))
     c   .and.(nreaco(6).eq.nreac(6)).and.(nreaco(7).eq.nreac(7))))
     c   then   !    ela. parton-parton collisions only 240412
	do i1=1,n
	do j1=1,5
	pbe(i1,j1)=p(i1,j1)
	enddo
	enddo
c	'sbe' to 'PYJETS'
	call tran_sbe 
	endif   ! 240412
c	write(22,*)'af. call tran_sbe iii,n=',iii,n   !240412
c	call pylist(1)   ! 240412
c	call prt_pyj(n)   ! 240412

c	write(9,*)'idii,idij,jb,nbe,n,j=',idii,idij,jb,nbe,n,j
c010305	if(nout.eq.1 .or. iii.eq.1 .or. iii.eq.neve)then
c300705	if(nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
c     c	 .eq.neve)then
c280809
c       if(iii.eq.12.or.iii.eq.29.or.iii.eq.151.or.iii.eq.213)then
c       write(22,*)'af. recover n,nbh,nbe,noo,idi=',n,nbh,nbe,n00,idi   ! sa
cs      write(9,*)'af. recover n,nbh=',n,nbh   ! sa
c       call pyedit(2)
c       call pylist(1)
cs      call prt_sbh(nbh)
c       endif
        endif   ! 280809
c280809   
889	continue
c	hadronization
c	write(9,*)'before hadronization iii=',iii   ! sa
c141208	if(adj112.eq.0)call sfm   ! string fragmentation 110604
        ijjk=0
c280809 if(adj112.eq.0)then   ! 240209
        if(adj112.eq.0.and.nreaco(4).eq.nreac(4).and.nreaco(6)
     c   .eq.nreac(6).and.nreaco(7).eq.nreac(7))then   ! 280809 150512
        call sfm(iii,ijjk,ss,kfa,kfb)   ! string fragmentation
c101210
c       chcd=0.
c       do i1=1,n
c       kfi1=k(i1,2)
c       chcd=chcd+pychge(kfi1)
c       enddo
c       chcd=chcd/3.
c       if(chcd.ne.chab)then
c       write(22,*)'iii,kfa,kfb,ss,chab,chcd=',iii,kfa,kfb,ss,chab,chcd
c       call pylist(1)
c       endif
c101210
        endif   ! 240209
c141208
c280809	if(adj112.ne.0)then
        if(adj112.ne.0.or.(adj112.eq.0.and.(nreac(4).gt.nreaco(4).or.
     c   nreac(6).gt.nreaco(6).or.nreac(7).gt.nreaco(7))))then ! 280809 150512
	call coales(iiii,neve,nout,nap,nat,nzp,nzt,kfa,kfb)   ! coalescence
	endif
887	format(20(1x,i3))
cs	write(9,*)'after hadniz iii,n,nbh=',iii,n,nbh   ! sa
c161203
c010305	if(nout.eq.1 .or. iii.eq.1 .or. iii.eq.neve)then
c300705	if(nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
c     c	 .eq.neve)then   
cs	write(mstu(11),*)'af hadniz n,nbh=',n,nbh   ! sa
c	call pyedit(2)
cs	call pylist(1)   ! sa
c130805	call prt_sbh(nbh)
c	write(22,*)'throe_t=',throe_t
c	write(22,*)'ithroq_t,ithrob_t,ich_t=',ithroq_t,ithrob_t,
c     c	 float(ich_t)/3.
c	write(22,*)'throe_p=',throe_p
c	write(22,*)'ithroq_p,ithrob_p,ich_p=',ithroq_p,ithrob_p,
c     c	 float(ich_p)/3.
c	write(22,*)'throe=',throe
c	write(22,*)'ithroq,ithrob,ithroc=',ithroq,ithrob,
c     c	 float(ithroc)/3.
c	endif
c121110
30000	n22=0   ! 140414
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        n22=n22+1
        endif
        enddo
c       move "22" from 'pyjets' to 'sgam'
        if(n22.gt.0)call remo_gam(22)
c	change K0S and K0L to K0 and K0ba
	do j=1,n
	kf=k(j,2)
        if(kf.eq.130 .or. kf.eq.310)then
        rrlu=pyr(1)
        k(j,2)=311
        if(rrlu.gt.0.5)k(j,2)=-311
        endif
        enddo
c130805	write(9,*)'af. sbh to pyjets n,nbh,naf=',n,nbh,naf   ! sa
c130805	call pylist(1)
c121110
	if(kjp21.eq.1)then
c       'pyjets' to 'sa1_h'
        nn=n
        do i1=1,n
        do i2=1,5
        kn(i1,i2)=k(i1,i2)
        pn(i1,i2)=p(i1,i2)
        rn(i1,i2)=v(i1,i2)
        enddo
        enddo
        do i=nn+1,kszj
        do j=1,5
        kn(i,j)=0
        pn(i,j)=0.
        rn(i,j)=0.
        enddo
        enddo
c       hadronic cascade (rescattering)
	time_had=0.d0
	call hadcas(iii,neve,nout,time_had,ijkk)
c       write(9,*)'after hadcas iii,ijkk,time_had=',
c     c  iii,ijkk,time_had   ! sa
        if(ijkk.eq.1)stop 'stop in hadcas'
c       'sa1_h' to 'pyjets'
        n=nn
        do i1=1,n
        do i2=1,5
        k(i1,i2)=kn(i1,i2)
        p(i1,i2)=pn(i1,i2)
        v(i1,i2)=rn(i1,i2)
        enddo
        enddo
        do i=n+1,kszj
        do j=1,5
        k(i,j)=0
        p(i,j)=0.
        v(i,j)=0.
        enddo
        enddo
        n66=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=66
        n66=n66+1
        endif
        enddo
c       move "66" from 'pyjets' to 'sgam'
        if(n66.gt.0)call remo_gam(66)
	endif
c121110
9002    if(adj140.lt.4)then   ! 200905 .ne.4 originally 140414
c270905
        naf0=naf
c       'pyjets' to 'saf'
        if(n.ge.1)then   ! 3
        do i1=1,n
        naf=naf+1
        if(naf.gt.kszj)then   ! 4
        write(9,*)'iiii,naf,kszj=',iiii,naf,kszj   ! sa
        stop 11111
        endif   ! 4
        do i2=1,5
        kaf(naf,i2)=k(i1,i2)
        paf(naf,i2)=p(i1,i2)
        vaf(naf,i2)=v(i1,i2)
        enddo
        enddo
        endif   ! 3
c       'sbh' to 'pyjets'
        if(nbh.eq.0)then   !
        n=0
        nup=0
	goto 9000   ! 200905
	endif   !
	do i1=1,nbh
	do m=1,5
	k(i1,m)=kbh(i1,m)
	p(i1,m)=pbh(i1,m)
	v(i1,m)=vbh(i1,m)
	enddo
	enddo
	n=nbh
9000	do i=n+1,kszj   ! 261103
	do j=1,5
	k(i,j)=0
	p(i,j)=0.
	v(i,j)=0.
	enddo
	enddo
	endif   ! 200905
9003	continue   ! 130206
c	change K0S and K0L to K0 and K0ba
	do j=1,n
	kf=k(j,2)
        if(kf.eq.130 .or. kf.eq.310)then
        rrlu=pyr(1)
        k(j,2)=311
        if(rrlu.gt.0.5)k(j,2)=-311
        endif
        enddo
        if(n.eq.0)nup=0   ! 270905
        if(n.ge.1)then   ! 270905
        call filt
c060813 120214 for pA,Ap,and AA return nucleons only
        do i=1,2   ! 221110
        nup=numbs(i)
        enddo
c	write(9,*)'n,k(i1,2)=',n,(k(i1,2),i1=1,n)
c	write(9,*)'numbs(1-2),nup=',numbs(1),numbs(2),nup
c060813	120214 for lepton+A, add lepton to nup, i.e. return lepton either
	if(ipden.ge.2)then   ! 060813 120214
	do i1=nup+1,n
	ki1ab=iabs(k(i1,2))   ! 060813 120214
	if(ki1ab.ge.11.and.ki1ab.le.16)ii=i1   ! 060813 120214
	enddo
	do jj=1,5
        kk6(jj)=k(ii,jj)
        pp6(jj)=p(ii,jj)
        vv6(jj)=v(ii,jj)
        enddo
	do j1=ii-1,nup+1,-1
        do jj=1,5
        k(j1+1,jj)=k(j1,jj)
        p(j1+1,jj)=p(j1,jj)
        v(j1+1,jj)=v(j1,jj)
        enddo
        enddo
	do jj=1,5
        k(nup+1,jj)=kk6(jj)
        p(nup+1,jj)=pp6(jj)
        v(nup+1,jj)=vv6(jj)
        enddo
	nup=nup+1
c	write(9,*)'2 n,nup,k(i1,2)=',n,nup,(k(i1,2),i1=1,n)
	endif
c060813
        naf1=n-nup
c       naf1 is the number of particles from 'pythia' etc. after filtor
c        and needs to be stored in 'saf'
c       nup is the number of particles from 'pythia' etc. after filtor
c        and needs to be updated
        if(naf1.eq.0)goto 9001
        do i=nup+1,n
        naf=naf+1
        do j=1,5
        kaf(naf,j)=k(i,j)
        paf(naf,j)=p(i,j)
        vaf(naf,j)=v(i,j)
        enddo
        enddo
        endif   ! 270905
9001    continue
c       update particle list after calling 'pythia' etc. (i. e. 'pyjets' to 
c        'sa2' and truncate 'sa2')
	call updpip(l,l1,icp,lc,tc,tw,time,nup,iii)
c080104
c	write(9,*)'af. updpip nsa,nctl=',nsa,nctl   ! sa
cs      call prt_sa2(nsa)   ! sa
cs      call prt_saf(naf)
c	do i1=1,nctl
c	write(9,*)'i1,lci,lcj,t=',i1,lc(i1,1),lc(i1,2),tc(i1)   ! sa
c	enddo
c080104
        l=lc(icp,1)
        l1=lc(icp,2)
c	write(9,*)'be. updtlp iii,nsa,nctl,l,l1=',iii,nsa,nctl,l,l1
c       update n-n collision (time) list after calling 'pythia' etc.
        call updtlp(l,l1,time,lc,tc,tw,nup,iii,kjp21)
c	write(9,*)'af. updtlp iii,nsa,nctl=',iii,nsa,nctl   ! sa
cs      call prt_sa2(nsa)
cs      call prt_saf(naf)
c	do i=1,nctl
c	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
c	enddo
	if(nctl.eq.0)goto 100   ! 021204
        goto 300   ! ss is enough to call pythia
        endif   ! if 1
c221110
	winel=0
c	if ss is not enough for calling pythia treat as els. collision then
	call coelas(l,l1,ss,pi,pj)   ! 010600
c       calculate four-momentum of two particles after elastic reaction, pi
c       and pj in CMS frame
c	write(9,*)'af. coelas iii,l,l1,ss=',iii,l,l1,ss
c	write(9,*)'pi=',(pi(i1),i1=1,4)
c	write(9,*)'pj=',(pj(i1),i1=1,4)
        call updple(l,l1,b,pi,pj,time)
c       update the particle list for elastic scattering,pi and pj have been
c       boosted back to Lab fram or cms of nucleus-nucleus collision
c	write(9,*)'af. updple nsa,l,l1=',nsa,l,l1
c	do i1=1,nsa
c	write(9,*)'i1,ksa(i1,2)=',i1,ksa(i1,2)
c	enddo
c	write(9,*)'af. updple pi=',(pi(i1),i1=1,4)
c	write(9,*)'af. updple pj=',(pj(i1),i1=1,4)
	noinel(1)=noinel(1)+1
c       update the collision list
	call updatl(l,l1,time,lc,tc,tw,winel,iii,kjp21) ! 010530
c	write(9,*)'af. updatl iii,nsa,nctl=',iii,nsa,nctl
c	do i1=1,nctl
c	write(9,*)'i1,lci,lcj,t=',i1,lc(i1,1),lc(i1,2),tc(i1)
c	enddo
	
c221110
300	iii=iii+1
c	if(iii.eq.8)return   ! temporal 240412
	if(iii.gt.100*(nctl0))then
        write(9,*)'infinite loop may have happened in'
        write(9,*)'subroutine scat iiii=',iiii
c10/08/98       stop 'infinite loop occurs'
        iiii=iiii-1   ! 10/08/98
        ijk=1  ! 10/08/98
        return   ! 10/08/98
        endif

	goto 10
100	continue
	call copl(time)
c070417
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        dnnc=dfloat(nnc)
        adiv=sadiv/dnnc
        gpmax=sgpmax/dnnc
        akapa(1)=skapa(1)/dnnc
        akapa(2)=skapa(2)/dnnc
        akapa(3)=skapa(3)/dnnc
        akapa(4)=skapa(4)/dnnc
        akapa(5)=skapa(5)/dnnc
        gtime=sgtime/dnnc
        astr=sastr/dnnc
        itime=sitime/dnnc
        endif
c070417
c121110	time_had=time   ! 111010
	return
	end

c************************************************************she042021
        subroutine chargecme(win)
c   the CME-induced charge initial charge separation by switching the 
c   py values of a fraction of the downward(upward) moving(u,d,s,c)quarks 
c   for symmetrical collision systems,i.e., Ru&Ru Zr&Zr at RHIC and LHC.
c   Here in symmetrical systems, nap=nat,nzp=nzt, and the fraction and
c   magnetic field function is A*bp-B*bp^3 type.  by shezl 2021

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      parameter(kszj=40000)
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
      common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
      common/schuds/schun,schudn,schudsn,sfra,cmes
      dimension numk(kszj)
      real(kind=8) p2u,erhic,erela,rerzcp,ruzcp

       erhic=200.                     ! RHIC energy 200 
       erela=0.45+0.55*(win/erhic)  !RHIC energy as a base
       rnzp=real(nzp)
       rnap=real(nap)
       rerz=rnzp/rnap
       ruzcp=((96./42.)*rerz)**(0.667) !isobar Zr Ru(96,42)as a base

       sfra=3.1*(2448.135*nap**(-1.667)*bp-160.810*nap**(-2.333)*bp**3.)
     c         *erela*ruzcp*0.01

c      print*,"erela,ruzcp,nap,nzp,bp,sfra",erela,ruzcp,nap,nzp,bp,sfra

       do 140  i=1,n
       if(abs(k(i,2)).eq.1.or.abs(k(i,2)).eq.2.or.abs(k(i,2)).eq.3
     c  .or.(k(i,2)).eq.4)then
        schun=schun+1
        if(pyr(1).gt.0..and.pyr(1).le.sfra)then
        numk(i)=0
        schudn=schudn+1
        do 150 ii=1,n
         if(numk(ii).eq.1) cycle
        do 160 jj=ii+1,n
         if(numk(jj).eq.1) cycle
          if((k(ii,2)+k(jj,2)).eq.0.and.(k(ii,2)*p(ii,2).lt.0).and.
     c         (k(jj,2)*p(jj,2).lt.0))then

             p2u=p(ii,2)
             p(ii,2)=p(jj,2)
             p(jj,2)=p2u
             schudsn=schudsn+1
             numk(ii)=1
             numk(jj)=1
             goto 160
          endif
160       enddo
150       enddo
          endif
          endif
140       enddo
          end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine lorntz(ilo,b,pi,pj)
c	perform Lorentz (or inverse Lorentz) transformation
c	implicit real*8 (a-h,o-z)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	dimension pi(4),pj(4),b(3),dpi(4),dpj(4)
	bb=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
	DB=DSQRT(bb)
        eps1=1d0-1d-12   ! 121108
c121108	IF(DB.GT.0.99999999D0) THEN
        IF(DB.GT.eps1) THEN   ! 121108 
	do i=1,3
c       rescale boost vector if too close to unity. 
c121108	b(i)=b(i)*(0.99999999D0/DB)
        b(i)=b(i)*(eps1/DB)   ! 121108 	
	enddo
c121108	DB=0.99999999D0
        DB=eps1   ! 121108 
	bb=DB**2
	endif
	bbb=1d0-bb
c	if(bbb.le.1.d-10)bbb=1.d-10
	gam=1d0/dsqrt(bbb)
	ga=gam*gam/(gam+1d0)
	do i=1,4
	dpi(i)=pi(i)
	dpj(i)=pj(i)
	enddo
	if(ilo.eq.1) goto 100
c	Lorentz transformation
	pib=dpi(1)*b(1)+dpi(2)*b(2)+dpi(3)*b(3)
	pjb=dpj(1)*b(1)+dpj(2)*b(2)+dpj(3)*b(3)
	do i=1,3
	pi(i)=dpi(i)+b(i)*(ga*pib-gam*dpi(4))
	pj(i)=dpj(i)+b(i)*(ga*pjb-gam*dpj(4))
	enddo
	pi(4)=gam*(dpi(4)-pib)
	pj(4)=gam*(dpj(4)-pjb)
	return
100	continue
c	inverse Lorentz transformation
	pib=dpi(1)*b(1)+dpi(2)*b(2)+dpi(3)*b(3)
	pjb=dpj(1)*b(1)+dpj(2)*b(2)+dpj(3)*b(3)
	do i=1,3
	pi(i)=dpi(i)+b(i)*(ga*pib+gam*dpi(4))
	pj(i)=dpj(i)+b(i)*(ga*pjb+gam*dpj(4))
	enddo
	pi(4)=gam*(dpi(4)+pib)
	pj(4)=gam*(dpj(4)+pjb)
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine his(t1,lc,tc,tw,istop)
c	classical Newton motion in Lab. system
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	parameter (kszj=40000)
	parameter(nsize=240000)
        common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa4/tau(kszj),tlco(kszj,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/ctllist/nctl,noinel(600),nctl0,noel
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	dimension lc(nsize,5),tc(nsize),tw(nsize)
	istop=1
	in=0
	do 200 i=1,nsa
c	if(t1.le.tau(i))goto 100
c	do move particles which have not produced
	if(ishp(i).eq.1) goto 10
c	pp4=psa(i,4)
c	do j=1,3
c	vp=psa(i,j)/pp4
c	vsa(i,j)=vsa(i,j)+vp*(t1-vsa(i,4))	
c	enddo
	in=in+1
	goto 200   ! 100 271004
10	aa=0.
c090609
c220110        pp1=psa(i,1)
c        pp2=psa(i,2)
c        ppt=pp1*pp1+pp2*pp2
c        if(ppt.lt.1.e-12)ppt=1.e-12
c        ppt=sqrt(ppt)
c220110        if(ppt.le.1.e-4)goto 100
c       unwounded nucleon do not make Newton motion
c090609
	pp4=psa(i,4)
c	due to the fast speed of bayons, we could not use a limited interaction
c	region
	r0=rao*max(rnt,rnp)
c	if(abs(k(i,2)).gt.1000)r0=1.E+10*r0
	do j=1,3
	vp=psa(i,j)/pp4
	vsa(i,j)=vsa(i,j)+vp*(t1-vsa(i,4))
	aa=aa+(vsa(i,j)-coor(j))**2
	enddo
c251004	vsa(i,4)=t1
	aa=sqrt(aa)
	if(aa.lt.r0) goto 100
c	if freeze-out occurs deduct the distance between the last collision 
c	and current collision
	do j=1,3
	vp=psa(i,j)/pp4
	vsa(i,j)=vsa(i,j)-vp*(t1-vsa(i,4))
	enddo
	ishp(i)=0
	do il=1,nctl
	if(lc(il,1).eq.i.or.lc(il,2).eq.i) tc(il)=0.
	enddo
	goto 200   ! 271004
100	continue
	vsa(i,4)=t1   ! 251004
200	continue
	if(in.eq.nsa) return
	istop=0
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ptcre1(l,l1,time)   ! 110517 
c	give four position to the particles after calling pythia  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	PARAMETER (kszj=40000,KSZ1=30)
	COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	do i=1,n
c060813	if(ipden.ne.0 .or. itden.ne.0)then
	rl=pyr(1)
	do m=1,3
c	write(*,*)'v=',v(i,m),k(i,2),time
	v(i,m)=v(i,m)+vsa(l,m)*rl+vsa(l1,m)*(1.-rl)
	enddo
c060813	endif
	if(ipden.eq.0 .and. itden.eq.0)then
	cita=2*pyr(1)-1.
        fi=2.*pio*pyr(1)
        sita=sqrt(1.-cita**2)
	v(i,1)=sita*cos(fi)
	v(i,2)=sita*sin(fi)
	v(i,3)=cita
	endif
	v(i,4)=time   ! 230805
	enddo
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ptcre(l,l1,time,ii)   ! 110517   
c	give four position to the particles after calling pythia  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	PARAMETER (kszj=40000,KSZ1=30)
	COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)   ! 110517   
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	if(ii.eq.1)then   ! for particles in pyjets 110517
	do i=1,n
	if(ipden.ne.0 .or. itden.ne.0)then
	rl=pyr(1)
	do m=1,3
c	write(*,*)'v=',v(i,m),k(i,2),time
	v(i,m)=v(i,m)+vsa(l,m)*rl+vsa(l1,m)*(1.-rl)
	enddo
	endif
	if(ipden.eq.0 .and. itden.eq.0)then
	cita=2*pyr(1)-1.
        fi=2.*pio*pyr(1)
        sita=sqrt(1.-cita**2)
	v(i,1)=sita*cos(fi)
	v(i,2)=sita*sin(fi)
	v(i,3)=cita
	endif
	v(i,4)=time   ! 230805
	enddo
	endif   ! 110517
c110517
        if(ii.eq.2)then   ! for particles in sbh
        do i=1,nbh
        if(ipden.ne.0 .or. itden.ne.0)then
        rl=pyr(1)
        do m=1,3
c       write(*,*)'v=',v(i,m),k(i,2),time
        vbh(i,m)=vbh(i,m)+vsa(l,m)*rl+vsa(l1,m)*(1.-rl)
        enddo
        endif
        if(ipden.eq.0 .and. itden.eq.0)then
        cita=2*pyr(1)-1.
        fi=2.*pio*pyr(1)
        sita=sqrt(1.-cita**2)
        vbh(i,1)=sita*cos(fi)
        vbh(i,2)=sita*sin(fi)
        vbh(i,3)=cita
        endif
        vbh(i,4)=time   ! 230805
        enddo
        endif   
c110517
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine arrove(ii,jj,sumx,sumy,sumxy,sumx2,sumy2,sump,
     c	 alp,r0,am,ac)   ! 101014   
c	 ! 191110 270312 131212
c	arrange randomly particle ii in overlap region of colliding nuclei 
c	jj=0 and 1 for target and projectile, respectively  
	PARAMETER (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc 
	common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	b=bp
        iii=0
54      iii=iii+1
        if(iii.eq.10000)then
        write(9,*)'difficult to arrange produced nucleon in'
        write(9,*)'subroutine arrove,infinitive loop may happen'
	goto 55   ! set larget number of try is equal to 10000
        endif
c       sample a point in the unit sphere 
c101014        x=1.-2.*pyr(1)
c        y=1.-2.*pyr(1)
c        z=1.-2.*pyr(1)
c        rr=x*x+y*y+z*z
c101014        if(rr.gt.1) goto 54
	if(jj.eq.0)then   ! ii in target (origin)
c101014	x=x*r0t
c	y=y*r0t
c101014	z=z*r0t
c101014
c	sample a point according to woodsax distribution
	call woodsax_samp(ii,jj,alp,r0,am,ac,0)
	x=c17(ii,1)
        y=c17(ii,2)
        z=c17(ii,3)   
c101014
c       relative to projectile center, they are b-x, y, and z, respectively 
c	(x-b,y,z) is in the sphere of projectile ?
        r1=sqrt((b-x)*(b-x)+y*y+z*z)
        if(r1.gt.r0p)goto 54
c 101014        c17(ii,1)=x
c        c17(ii,2)=y
c101014        c17(ii,3)=z
c270312
	sumx=sumx+x
	sumy=sumy+y
	sumxy=sumxy+x*y   ! 131212
        sumx2=sumx2+x*x
        sumy2=sumy2+y*y
	sump=sump+1.
c270312
	endif
	if(jj.eq.1)then   ! ii in projectile
c101014	x=x*r0p
c	y=y*r0p
c101014	z=z*r0p
c101014
c	sample a point according to woodsax distribution
        call woodsax_samp(ii,jj,alp,r0,am,ac,0)
        x=c17(ii,1)
        y=c17(ii,2)
        z=c17(ii,3)
c101014
c       relative to target center, they are x+b, y, and z, respectively 
c	(x+b,y,z) is in the sphere of target ?
        r1=sqrt((x+b)*(x+b)+y*y+z*z)
        if(r1.gt.r0t)goto 54
c101014	c17(ii,1)=x
c        c17(ii,2)=y
c101014        c17(ii,3)=z
c270312
	xb=x+b   ! 101014 chen
        sumx=sumx+xb   ! 101014 chen 
        sumy=sumy+y
	sumxy=sumxy+xb*y   ! 131212 101014
        sumx2=sumx2+xb*xb   ! 101014
        sumy2=sumy2+y*y
	sump=sump+1.
c270312
	endif
55	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rotate(cctas,sctas,cfis,sfis,pp3,pi,pj)
c	perform rotation
c       pi,pj: input, four momentum of colliding pair before scattering
c              output, four momentum of scattered particles after rotation
c       pp3: momentum modulus of pi or pj, both are equal in their cms,
c        after scattering
c       cctas,sctas,cfis,sfis: direction cosines of momentum of one of
c        scattered particle after relative to the momentum
c        of corresponding particle before scattering
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	dimension pi(4),pj(4)
c	fi1=atan2(pi(2),pi(1))
c	cta1=atan2(sqrt(pi(1)**2+pi(2)**2),pi(3))
	fi1=pyangl(pi(1),pi(2))
	cta1=pyangl(pi(3),sqrt(pi(1)**2+pi(2)**2))
	cfi1=cos(fi1)
	sfi1=sin(fi1)
	ccta1=cos(cta1)
	scta1=sin(cta1)
	pi(1)=cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis
	pi(2)=sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis
	pi(3)=ccta1*cctas-scta1*sctas*cfis
	pi(1)=pp3*pi(1)
	pi(2)=pp3*pi(2)
	pi(3)=pp3*pi(3)
	do i=1,3
	pj(i)=0.-pi(i)
	enddo
	return
	end



c********************************************************************
        subroutine tran_saf
c       'saf' to 'pyjets' 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        PARAMETER (kszj=40000,KSZ1=30)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/saf/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)   ! 080104
        do l=1,nsa
        do m=1,5
        k(l,m)=ksa(l,m)
        p(l,m)=psa(l,m)
        v(l,m)=vsa(l,m)
        enddo
        enddo
        n=nsa
	do l=n+1,kszj
	do m=1,5
	k(l,m)=0
	p(l,m)=0.
	v(l,m)=0.
	enddo
	enddo
        return
        end

c********************************************************************
        subroutine tran_sbe
c       'sbe' to 'pyjets' 
        PARAMETER (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)   ! 080104
        do l=1,nsa
        do m=1,5
        k(l,m)=ksa(l,m)
        p(l,m)=psa(l,m)
        v(l,m)=vsa(l,m)
        enddo
        enddo
        n=nsa
	do l=n+1,kszj
	do m=1,5
	k(l,m)=0
	p(l,m)=0.
	v(l,m)=0.
	enddo
	enddo
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	subroutine conse(np,pp,ps,ii,jj)
c	keep four momentum conservation
c	np : the # of particles
c	ps : four momentum to which the four momenta of particles should 
c            conserves 
c	pp : four momenta of particles
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
      	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	dimension pp(250,5),ps(4),ff(250),pxyz(3),arp(3)
	ps4=ps(4)
	do i=1,3
	pxyz(i)=0.
	enddo
	jj=0
100	es=0.
	do i=1,np
	es=es+pp(i,4)
	enddo
	fr=es/ps4
	if(abs(1.-fr) .le. dep)goto 200
	do i=1,np
	ppm=pp(i,4)/0.938
	ppf=ppm/fr
	ff(i)=sqrt(abs(ppf*ppf-1.)/(ppm*ppm-1.))
	do j=1,3
	ppp=ff(i)*pp(i,j)
	pp(i,j)=ppp
	pxyz(j)=pxyz(j)+ppp
	enddo
	enddo
	do i=1,3
	arp(i)=abs(1.-pxyz(i)/ps(i))
	pxyz(i)=pxyz(i)-ps(i)
	enddo
	if(abs(1.-fr).le.dep .and.arp(1).le.dep .and. arp(2).le.dep  
     c   .and. arp(3).le.dep) goto 200
	do i=1,3
	pxyz(i)=pxyz(i)/np
	enddo
	do i=1,np
	do j=1,3
	pp(i,j)=pp(i,j)-pxyz(j)
	enddo
	pp(i,4)=sqrt(0.880+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
c	0.880 = 0.938*0.938
	enddo
	jj=jj+1
	if(jj.eq.4000)then
	write(9,*)'infinitive loop may occur in subroutine conse(),'
        write(9,*)'which means four-momentum conservation' 
        write(9,*)'needed is hard to be achieved,check value' 
        write(9,*)'of PARAM(9)'
	return
	endif
	goto 100
200	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine codi(pis,cfi1,sfi1,ccta1,scta1)
c	calculate the 'orientation' of the vector pis
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	dimension pis(4),pi(4)
c	do i=1,4
c	pi(i)=pis(i)
c	enddo
c	if(pi(1).lt.1.d-15)pi(1)=1.d-15 
c	fi1=atan2(pi(2),pi(1))
c        cta1=atan2(sqrt(pi(1)**2+pi(2)**2),pi(3))
	fi1s=pyangl(pis(1),pis(2))
        cta1s=pyangl(pis(3),sqrt(pis(1)**2+pis(2)**2))
	fi1=fi1s
	cta1=cta1s
        cfi1=dcos(fi1)
        sfi1=dsin(fi1)
        ccta1=dcos(cta1)
        scta1=dsin(cta1)
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rosa(cfi1,sfi1,ccta1,scta1,cfis,sfis,cctas,sctas,
     c	 pis)
c       perform rotate for produced particles from 'pythia'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        dimension pis(4),pi(4)
	do i=1,4
	pi(i)=pis(i)
	enddo
	pp=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)
	pp=dsqrt(pp)   ! 110517
	call codi(pis,cfis,sfis,cctas,sctas)
        pi(1)=pp*(cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis)
        pi(2)=pp*(sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis)
        pi(3)=pp*(ccta1*cctas-scta1*sctas*cfis)
	do i=1,4
	pis(i)=pi(i)
	enddo
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine flow(tt,tt1)
c	calculate directed and elliptic flow
	parameter (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa18/tdh,itnum,non18,cptl,cptu,cptl2,cptu2,snum(4,20),
     &	 v1(4,20),v2(4,20),v12(4,20),v22(4,20)
	do i1=1,itnum
	ti=(i1-1)*tdh
	if(tt.le.ti .and. tt1.gt. ti)then
c       since the particles after collision at time0 have contribution
c        to the statistics at any moment in between this collision (
c        time0) and next collision (time)

        do j=1,nsa
        ik=ksa(j,2)
c	yy=pyp(j,17)
c        if(ik.eq.2212 .and. (yy.gt.0.9 .and. yy.le.5.))then
	if(ik.eq.2212)then
	px=psa(j,1)
	py=psa(j,2)
	px2=px*px
	py2=py*py
	pt2=px2+py2
	pt=sqrt(pt2)
	if(pt.gt.cptl .and. pt.le.cptu)then
	snum(1,i1)=snum(1,i1)+1.
	pxt=px/pt
	pxt2=(px2-py2)/pt2
	v1(1,i1)=v1(1,i1)+pxt
	v2(1,i1)=v2(1,i1)+pxt2
	v12(1,i1)=v12(1,i1)+pxt*pxt
	v22(1,i1)=v22(1,i1)+pxt2*pxt2
	endif
	endif
	if(iabs(ik).eq.211)then
	px=psa(j,1)
	py=psa(j,2)
	px2=px*px
	py2=py*py
	pt2=px2+py2
	pt=sqrt(pt2)
	if(pt.gt.cptl2 .and. pt.le.cptu2)then
	snum(2,i1)=snum(2,i1)+1.
	pxt=px/pt
	pxt2=(px2-py2)/pt2
	v1(2,i1)=v1(2,i1)+pxt
	v2(2,i1)=v2(2,i1)+pxt2
	v12(2,i1)=v12(2,i1)+pxt*pxt
	v22(2,i1)=v22(2,i1)+pxt2*pxt2
	endif
	endif
	enddo

	endif
	enddo
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo
c	move q,qbar,g,anti-diquark and diquark from 'pyjets' to 'sbh' 110517 
      PARAMETER (KSZJ=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	nbh=0
	do i1=1,kszj
        do j1=1,5
        kbh(i1,j1)=0
        pbh(i1,j1)=0.
        vbh(i1,j1)=0.
        enddo
        enddo
        jb=0

201     do i1=jb+1,n
        kf=k(i1,2)
        kfab=iabs(kf)
        if(kfab.gt.8 .and. kfab.ne.2101 .and. kfab.ne.3101
     c   .and. kfab.ne.3201 .and. kfab.ne.1103 .and. kfab.ne.2103
     c   .and. kfab.ne.2203 .and. kfab.ne.3103 .and. kfab.ne.3203
     c   .and. kfab.ne.3303 .and. kfab.ne.21)then
        jb=jb+1
        goto 202
        endif
c	write(9,*)'n,i1,jb=',n,i1,jb   ! sa
        nbh=nbh+1
        do i2=1,5
        kbh(nbh,i2)=k(i1,i2)
        pbh(nbh,i2)=p(i1,i2)
        vbh(nbh,i2)=v(i1,i2)
        enddo
        if(i1.eq.n)then
        n=n-1
        goto 203
        endif
c	move particle list one step downward from i1+1 to n
        do j=i1+1,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 201
202     enddo

203     continue
	return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break
c       break up diquark and give four momentum and four position
c        to the broken quarks
      PARAMETER (KSZJ=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 080104 280809
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        amd=0.00990   ! pymass(1)
        amu=0.00560   ! pymass(2)
        ams=0.199   ! pymass(3)
        amuu=2*amu
        amdd=2*amd
        amss=2*ams
c170205
c       the probability of gluon spliting into u,d & s
        adj132=adj1(32)
        prosum=1.+1.+adj132
        prod=1./prosum   ! 0.4286 originally
        pros=adj132/prosum   ! 0.1428 originally
        prods=prod+pros   ! 0.5714 originally
c170205
	jb=0
	ii=idio   ! 080104
100     do i1=jb+1,n
        kf=k(i1,2)
	kfab=iabs(kf)
        if(kfab.ne.2101 .and. kfab.ne.3101
     c   .and. kfab.ne.3201 .and. kfab.ne.1103 .and. kfab.ne.2103
     c   .and. kfab.ne.2203 .and. kfab.ne.3103 .and. kfab.ne.3203
     c   .and. kfab.ne.3303)then
c     c   .and. kfab.ne.3303 .and. kfab.ne.21)then
        jb=jb+1
        goto 300
        endif

	if(kf.eq.2101)then
	kf1=2
	kf2=1
	goto 200
	endif
	if(kf.eq.3101)then
        kf1=3
        kf2=1
        goto 200
	endif
	if(kf.eq.3201)then
        kf1=3
        kf2=2
        goto 200
	endif
	if(kf.eq.1103)then
        kf1=1
        kf2=1
        goto 200
	endif
	if(kf.eq.2103)then
        kf1=2
        kf2=1
        goto 200
	endif
	if(kf.eq.2203)then
        kf1=2
        kf2=2
        goto 200
	endif
	if(kf.eq.3103)then
        kf1=3
        kf2=1
        goto 200
	endif
	if(kf.eq.3203)then
        kf1=3
        kf2=2
        goto 200
	endif
	if(kf.eq.3303)then
        kf1=3
        kf2=3
        goto 200
	endif
c251103
	if(kf.eq.-2101)then
	kf1=-2
	kf2=-1
	goto 200
	endif
	if(kf.eq.-3101)then
        kf1=-3
        kf2=-1
        goto 200
	endif
	if(kf.eq.-3201)then
        kf1=-3
        kf2=-2
        goto 200
	endif
	if(kf.eq.-1103)then
        kf1=-1
        kf2=-1
        goto 200
	endif
	if(kf.eq.-2103)then
        kf1=-2
        kf2=-1
        goto 200
	endif
	if(kf.eq.-2203)then
        kf1=-2
        kf2=-2
        goto 200
	endif
	if(kf.eq.-3103)then
        kf1=-3
        kf2=-1
        goto 200
	endif
	if(kf.eq.-3203)then
        kf1=-3
        kf2=-2
        goto 200
	endif
	if(kf.eq.-3303)then
        kf1=-3
        kf2=-3
        goto 200
	endif
c251103
	if(kf.eq.21)then   ! 1
	eg=p(i1,4)
        if(eg.lt.amdd)then   ! 2  
        kf1=2
        kf2=-2
        goto 200
        elseif(eg.ge.amdd .and. eg.lt.amss)then   ! 2
        kf1=2
        kf2=-2
        if(pyr(1).gt.0.5)then
        kf1=1
        kf2=-1
        endif
        goto 200
        elseif(eg.gt.amss)then   ! 2
        kf1=3
        kf2=-3
        rand=pyr(1)
        if(rand.gt.pros .and. rand.le.prods)then
        kf1=1
        kf2=-1
        endif
        if(rand.gt.prods)then
        kf1=2
        kf2=-2
        endif
        goto 200
	else   ! 2
	goto 300
	endif   ! 2
	endif   ! 1
200	k(i1,2)=kf1
        k(n+1,2)=kf2
c221203	k(i1,1)=1
        k(n+1,1)=3   ! 280809, 1 originally
c221203	k(i1,3)=0
        k(n+1,3)=0
c221203
        k(n+1,4)=0
        k(n+1,5)=0
c221203
c080104
	ii=ii+1
c050805	npt(ii)=n+1+naf
	npt(ii)=n+1   ! 050805
        ifcom(ii)=i1   ! 280809
c080104
c       give four momentum to the breaked quarks
	call bream(i1,kf1,kf2)
c       give four coordinate to the breaked quarks
        call coord(i1)
	if(i1.eq.n)then
	n=n+1
	goto 400
	endif
        n=n+1
        goto 100

300     enddo
400	return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine bream(ii,kf1,kf2)
c       give four momentum to the broken quarks
c       ii: line number of diquark in 'pyjets'
c       kf1,kf2: flavor codes of broken quarks
      PARAMETER (KSZJ=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)   ! 260503
        am1=pymass(kf1)
        am2=pymass(kf2)
        pp(1,5)=am1
        pp(2,5)=am2
c       pp : four momentum of breaked quarks, local variable 
        do i1=1,4
        ps(i1)=p(ii,i1)
        enddo
c       ps : four momentum of diquark, local variable 
	goto 400   ! activate it for 'decay method'
c       broken quarks share out diquark four momentum randomly,
c        denoted as 'random four momentum method'
c       do i1=1,4   ! activate it for 'random four momentum method'
c       broken quarks share out diquark three momentum randomly,
c        denoted as 'random three momentum method'
401	do i1=1,3   ! activate it for 'random three momentum method'
        pi(i1)=pyr(1)*p(ii,i1)
        pp(2,i1)=ps(i1)-pi(i1)
        pp(1,i1)=pi(i1)
        enddo
c	goto 300   ! activate it for 'random four momentum method'
c250503
	pp11=pp(1,1)
	pp12=pp(1,2)
	pp13=pp(1,3)
c021005
        pp14=am1*am1+pp11*pp11+pp12*pp12+pp13*pp13
        if(pp14.le.0.)pp14=1.e-20
        pp(1,4)=sqrt(pp14)
c021005
	pp21=pp(2,1)
	pp22=pp(2,2)
	pp23=pp(2,3)
c021005
        pp24=am2*am2+pp21*pp21+pp22*pp22+pp23*pp23
        if(pp24.le.0.)pp24=1.e-20
        pp(2,4)=sqrt(pp24)
c021005
	goto 300   ! activate it for 'random three momentum method'
c250503
c260503
400	continue
	decsuc=1
	call decmom(ps,pp,am1,am2,decsuc)
	if(decsuc.eq.0)goto 401   ! return to random three momentum method
300	continue
c       adjust four momentum conservation by iteration,no more than
c        4000 iterations
c	call conser(2,pp,ps)   
c260503
c        write(9,*)'after	'
c        do i=1,2
c        write(9,*)(pp(i,j),j=1,5)
c        enddo
c260503
        do i1=1,4
        p(ii,i1)=pp(1,i1)
        enddo
        p(ii,5)=am1
        do i1=1,4
        p(n+1,i1)=pp(2,i1)
        enddo
        p(n+1,5)=am2
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine decmom(ps,pp,am1,am2,decsuc)
c	calculate four momentum of decayed particles
c	ps: four momentum of decaying particle
c	am1 and am2: mass of decayed pair
	parameter(kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4),pj(4),ps(4),pp(20,5)   
        double precision bb(3)
c       calculate the E and |p| of broken quark in rest frame of diquark
        sm2=ps(4)*ps(4)-ps(1)*ps(1)-ps(2)*ps(2)-ps(3)*ps(3)
c       one problem here is that 'sm2' may not equal to square of diquark 
c	 (gluon) rest mass,'bream' is called for spliting g especially
c030603
c1	if(sm2.lt.1.e-10)then
c1	sm2=1.e-10
c1	endif
c       write(9,*)'in decmom sm2=',sm2   ! sa
c       write(9,*)'ps=',ps   ! sa
c030603
	if(sm2.lt.0.005)then   ! 110211
	decsuc=0   ! go back to random three momentum method
	return
	endif
        sm=sqrt(sm2)   ! M (should be diquark mass)
c       pp(1,4)=(sm2-am2*am2+am1*am1)/2./sm
c       pp(2,4)=(sm2-am1*am1+am2*am2)/2./sm
        ppp=(sm2-(am1+am2)*(am1+am2))*(sm2-(am1-am2)*(am1-am2))
c161204	ppp=abs(ppp)   ! 030603 ?
	if(ppp.lt.1.e-28)ppp=1.e-28   !161204
        ppp=sqrt(ppp)/2./sm
c110211 goto 500   ! activate it for exponential cos(seta) distribution
c       the direction of broken quark is sampled isotropically in '4pi'
        coset=1.-2.*pyr(1)
        if(abs(coset).gt.1.)then
        coset=coset/abs(coset)
        endif
c021005
        siset=1.-coset*coset
        if(siset.lt.1.e-28)siset=1.e-28
c021005
        siset=sqrt(siset)   ! 021005
100     cosi1=pyr(1)
        cosi12=cosi1*cosi1
        eta2=2.*pyr(1)-1.
        eta22=eta2*eta2
        coseta=cosi12+eta22
        if(coseta.gt.1.)goto 100
        if(coseta.lt.1.e-28)coseta=1.e-28
        cofi=(cosi12-eta22)/coseta
        sifi=2.*cosi1*eta2/coseta
        goto 600
500     continue
c       cos(seta) is sampled from exponential distribution when
c        0<seta<pi/2 and its absolute value is assumed to be symmetry
c        about seta=pi/2. 'fi' is assumed to be isotropic in 2pi
        coset=log(1.+1.7183*pyr(1))
        if(pyr(1).lt.0.5)coset=-coset
c021005
        siset=1.-coset*coset
        if(siset.lt.1.e-28)siset=1.e-28
c021005
        siset=sqrt(siset)
        fi=2.*3.1416*pyr(1)
        cofi=cos(fi)
        sifi=sin(fi)
600     continue
        pi(1)=ppp*siset*cofi
        pi(2)=ppp*siset*sifi
        pi(3)=ppp*coset
c021005
        pi4=ppp*ppp+am1*am1
        if(pi4.lt.1.e-28)pi4=1.e-28
c021005
        pi(4)=sqrt(pi4)
        pj(1)=-pi(1)
        pj(2)=-pi(2)
        pj(3)=-pi(3)
c021005
        pj4=ppp*ppp+am2*am2
        if(pj4.lt.1.e-28)pj4=1.e-28
c021005
        pj(4)=sqrt(pj4)
c       write(9,*)'before rotation'   ! sa
c       write(9,*)ppp,(pi(i),i=1,4)   ! sa
c       write(9,*)(pj(i),i=1,4)   ! sa
c050603
c       rotate to the frame where diquark (gluon), ps, is described
c       calculate the direction cosines of ps
        fi1=pyangl(ps(1),ps(2))
c021005
        ps12=ps(1)*ps(1)+ps(2)*ps(2)
        if(ps12.lt.1.e-28)ps12=1.e-28
        ps12=sqrt(ps12)
c021005
        cta1=pyangl(ps(3),ps12)
        cfi1=cos(fi1)
        sfi1=sin(fi1)
        ccta1=cos(cta1)
        scta1=sin(cta1)
        sctas=siset
        cctas=coset
        sfis=sifi
        cfis=cofi
        pi(1)=cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis
        pi(2)=sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis
        pi(3)=ccta1*cctas-scta1*sctas*cfis
        pi(1)=ppp*pi(1)
        pi(2)=ppp*pi(2)
        pi(3)=ppp*pi(3)
        do i=1,3
        pj(i)=0.-pi(i)
        enddo
c021005
        pi4=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)+am1*am1
        if(pi4.lt.1.e-28)pi4=1.e-28
        pj4=pj(1)*pj(1)+pj(2)*pj(2)+pj(3)*pj(3)+am2*am2
        if(pj4.lt.1.e-28)pj4=1.e-28
        pi(4)=sqrt(pi4)
        pj(4)=sqrt(pj4)
c021005
c       write(9,*)'after rotation'   ! sa
c       write(9,*)(pi(i),i=1,4)   ! sa
c       write(9,*)(pj(i),i=1,4)   ! sa
c       boost to moving frame of diquark
        ee=ps(4)
        if(ee.lt.1.e-14)ee=1.e-14   ! 021005
        do i1=1,3
        bb(i1)=ps(i1)/ee
        enddo
c       write(9,*)'b=',(bb(i),i=1,3)
        call lorntz(1,bb,pi,pj)
c       write(9,*)'after boost back, ps=',(ps(i),i=1,4)   ! sa
c       write(9,*)(pi(i),i=1,4)   ! sa
c       write(9,*)(pj(i),i=1,4)   ! sa
c       write(9,*)(pi(i)+pj(i),i=1,4)   ! sa
c050603
        pp(1,1)=pi(1)
        pp(1,2)=pi(2)
        pp(1,3)=pi(3)
        pp(1,4)=pi(4)
        pp(2,1)=pj(1)
        pp(2,2)=pj(2)
        pp(2,3)=pj(3)
        pp(2,4)=pj(4)
c050603
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conser(np,pp,ps)
c       adjust four momentum conservation
c       np: the # of particles
c       pp: four momenta of particles have to be conserved
c       ps: the four momentum should be conserved to
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        dimension pp(20,5),ps(4),ff(20),pxyz(3),arp(3)
        ps4=ps(4)
        do i=1,3
        pxyz(i)=0.
        enddo
        jj=0
100     es=0.
        do i=1,np
        es=es+pp(i,4)
        enddo
        fr=es/ps4
        if(abs(1.-fr) .le. dep)goto 200
        do i=1,np
        amas=pp(i,5)
        ppm=pp(i,4)/amas
        ppf=ppm/fr
        ff(i)=sqrt(abs(ppf*ppf-1.)/(ppm*ppm-1.))
        do j=1,3
        ppp=ff(i)*pp(i,j)
        pp(i,j)=ppp
        pxyz(j)=pxyz(j)+ppp
        enddo
        enddo
        do i=1,3
        arp(i)=abs(1.-pxyz(i)/ps(i))
        pxyz(i)=pxyz(i)-ps(i)
        enddo
        if(abs(1.-fr).le.dep .and. arp(1).le.dep .and. arp(2).le.dep
     c   .and. arp(3).le.dep) goto 200
        do i=1,3
        pxyz(i)=pxyz(i)/np
        enddo
        do i=1,np
        do j=1,3
        pp(i,j)=pp(i,j)-pxyz(j)
        enddo
        pp5=pp(i,5)
        pp52=pp5*pp5
        pp(i,4)=sqrt(pp52+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
        enddo
        jj=jj+1
        if(jj.eq.4000)then
        write(9,*)'infinitive loop may occur in subroutine conser(),'
        write(9,*)'which means four-momentum conservation'
        write(9,*)'needed is hard to be achieved,check value'
        write(9,*)'of PARAM(9)'
        return
        endif
        goto 100
200     return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coord(ii)
c       give four position to broken quarks
c       first broken quark takes the four position of diquark
c       second broken quark is arranged around first ones within
c        0.5 fm randumly in each of three position coordinates and has same
c        fourth position coordinate as diquark
      PARAMETER (KSZJ=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension rr(3)
        do i1=1,3
c261002        rr(i1)=pyr(1)*v(ii,i1)
        rr(i1)=pyr(1)*0.5   ! 261002
        v(n+1,i1)=v(ii,i1)+rr(i1)
        if(pyr(1).gt.0.5)v(n+1,i1)=v(ii,i1)-rr(i1)
        enddo
        v(n+1,4)=v(ii,4)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine copl(tt)
c	calculate coordinate of center of mass of non-freeze-out system
c	position of a particle, checking is it freezes out or not, is 
c	 calculated with respect to this origin.
	parameter(kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        COMMON/SA2/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
     	 COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	common/sa4/tau(kszj),tlco(kszj,4)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	do ii=1,3
	coor(ii)=0.
	enddo
	samass=0.
	do 110 ii=1,n
c	if(tau(ii).gt.tt)goto 110
	if(ishp(ii).eq.0)goto 110
	kf=k(ii,2)
	amass=pmas(pycomp(kf),1)
	if(abs(kf).eq.213 .or. kf.eq.113)amass1=p(ii,5)   ! 010600 
	if((abs(kf).eq.213 .or. kf.eq.113) .and. abs(amass-amass1)
     &	 .gt.0.001)amass=amass1   ! 010600 
	samass=samass+amass
	do 100 jj=1,3
	coor(jj)=coor(jj)+amass*v(ii,jj)
100	continue
110	continue
	do ii=1,3
	coor(ii)=coor(ii)/max(0.14,samass)
	enddo
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ctlcre(lc,tc,tw)
c	create initial collision list  
	parameter(nsize=240000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
      	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
	common/ctllist/nctl,noinel(600),nctl0,noel
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c	nap,nat,nzp,nzt,pio
	dimension lc(nsize,5),tc(nsize),tw(nsize)
	time=0.
	nctl=1
	dminf=100.
	nzpab=iabs(nzp)   ! in order to consider ppbar or pbarp
	nztab=iabs(nzt)
	nzpt=nzpab+nztab
	napt=nap+nat
c	write(9,*)'ctlcre nzpab,nztab,nzpt,napt=',nzpab,nztab,nzpt,napt
	do 10 l=1,nzpab   ! projectile proton or lepton ! 060813 120213
	do l1=nzpab+1,nzpt   ! target proton
	tc(nctl)=0.
	mtc=0
	call coij(l,l1,nctl,lc,tc,tw,mtc,dminf,if,jf)
c	write(9,*)'tp af. coij l,l1,mtc,nctl=',l,l1,mtc,nctl
	if(mtc.gt.0)then
	nctl=nctl+1
	mtc=0
	endif
	enddo
	do l1=nap+nztab+1,napt   ! target neutron
        tc(nctl)=0.
	mtc=0
        call coij(l,l1,nctl,lc,tc,tw,mtc,dminf,if,jf)
c	write(9,*)'tn af. coij l,l1,mtc,nctl=',l,l1,mtc,nctl
        if(mtc.gt.0)then
	nctl=nctl+1
	mtc=0
	endif
        enddo
10      continue
        do 20 l=nzpt+1,nap+nztab   ! projectile neutron
c	write(9,*)'pn nzpt+1,nap+nztab=',nzpt+1,nap+nztab
	do l1=nzpab+1,nzpt   ! target proton	
        tc(nctl)=0.
	mtc=0
        call coij(l,l1,nctl,lc,tc,tw,mtc,dminf,if,jf)
        if(mtc.gt.0)then
	nctl=nctl+1
	mtc=0
	endif
        enddo
	do l1=nap+nztab+1,napt   ! target neutron
        tc(nctl)=0.
	mtc=0
        call coij(l,l1,nctl,lc,tc,tw,mtc,dminf,if,jf)
        if(mtc.gt.0)then
	nctl=nctl+1
	mtc=0
	endif
        enddo
20	continue
        if(mtc.eq.0)nctl=nctl-1
c	do iop=1,nctl
c	write(*,*)lc(iop,1),lc(iop,2),tc(iop)
c	enddo	
	if(nctl.eq.0)then
c 	at least one collision should occur. this collision has the smallest 
c	 'least approaching distance', that is guaranteed by the variable 
c	 'dminf'
	lc(1,1)=if
	lc(1,2)=jf
	tc(1)=0.02
	nctl=1
	endif
	do i=nctl+1,nsize
	do m=1,5
	lc(i,m)=0
	enddo
	tc(i)=0.
	tw(i)=0.
	enddo
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coij(i,j,icp,lc,tc,tw,mtc,dminf,if,jf)
c       calculate collision time & fill up lc(i,1-2) as well as tc(i)
c	 for creating the initial collsion list 
	PARAMETER (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	parameter(nsize=240000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension lc(nsize,5),tc(nsize),tw(nsize)
	dimension dr(3),db(3),dv(3),px(4),py(4),vi(4),vj(4)
     c	,pi(4),pj(4),b(3)
	ki=ksa(i,2)
	kj=ksa(j,2)
	pi(4)=psa(i,4)
	if(pi(4).lt.1.e-20)pi(4)=1.e-20
	pj(4)=psa(j,4)
	if(pj(4).lt.1.e-20)pj(4)=1.e-20
	deno6=pi(4)+pj(4)
	do k=1,3
	pi(k)=psa(i,k)
	pj(k)=psa(j,k)
	b(k)=(pi(k)+pj(k))/deno6
	enddo
	ilo=0
	call lorntz(ilo,b,pi,pj)
	do l=1,3
	px(l)=vsa(i,l)
	py(l)=vsa(j,l)
	enddo
	px(4)=0.
	py(4)=0.
	call lorntz(ilo,b,px,py)
	rb=0.
	bb=0.
	rr=0.
	rtai=0.
	do k=1,3
	vi(k)=pi(k)/pi(4)
	vj(k)=pj(k)/pj(4)
	enddo
	do k=1,3
	dr(k)=px(k)-py(k)-(vi(k)*px(4)-vj(k)*py(4))
	db(k)=vi(k)-vj(k)
	dv(k)=db(k)
	rb=rb+dr(k)*db(k)
	bb=db(k)**2+bb
	rr=rr+dr(k)*dr(k)
	enddo
	if(bb.le.1.e-10) return
	tcol=0.-rb/bb
c        if(tcol-px(4) .le. ddt)return
c        if(tcol-py(4) .le. ddt)return
c       for collision to occur,time must one step ahead
csa	if(tcol.lt.1.0e-7)return
	do iik=1,3
	dr(iik)=px(iik)-py(iik)-(vi(iik)*
     &	px(4)-vj(iik)*py(4))+tcol*db(iik)
	rtai=rtai+dr(iik)*dr(iik)
	enddo
	sg=rtai
c	sg=rr+tcol*rb
c	if(sg.lt.0)then
c	write(*,*)'sg=',sg   !
c	return
c	endif
	dmin=sqrt(sg)
	if(dmin.lt.dminf)then
	dminf=dmin
	if=i
	jf=j
	endif
c	write(9,*)'coij i,j,dmin,ecsen=',i,j,dmin,ecsen
	if(ipden.lt.2 .and. (dmin.gt.ecsnn))return   ! 060813 120214
	if(ipden.gt.2 .and. (dmin.gt.ecsen))return   ! 060813 120214
c	distance between the two particles should be smaller than ecsnn (ecsen)
c	 060813
	do ik=1,3
	px(ik)=px(ik)+vi(ik)*(tcol-px(4))
	py(ik)=py(ik)+vj(ik)*(tcol-py(4))
	enddo
c	move along Newton trajectory in CMS
	px(4)=tcol
	py(4)=tcol
c	write(*,*)'CMStcol=',tcol
	ilo=1

	call lorntz(ilo,b,px,py)

c	transform back to Lab.
	if(px(4).gt.py(4)) px(4)=py(4)
	tcol=px(4)
c041204
	drmax=rao*max(rnt,rnp)
c	write(9,*)'i,j,tcol,drmax=',i,j,tcol,drmax
	if(tcol.le.drmax)goto 180   ! 041204
        return   ! 041204
c041204
180	tc(icp)=tcol
	mtc=1
        lc(icp,1)=i
        lc(icp,2)=j
c	write(*,*)'LABtcol=',tcol
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find(icp,tcp,lc,tc,tw,ico)
c	find out the binary collision with least collision time
	parameter(nsize=240000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	common/ctllist/nctl,noinel(600),nctl0,noel
	dimension lc(nsize,5),tc(nsize),tw(nsize)
	icp=0
	tcp=20000.
	do i=1,nctl
	if(ico.eq.0)goto 100
	if(tc(i).le.1.0e-7) goto 241
100	if(tcp.lt.tc(i))  goto 241
	icp=i
	tcp=tc(i)
241	continue
	enddo
	if(nctl.eq.0)icp=0
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine tcolij(l,l1,time,icp,lc,tc,tw)
c	calculate collision time & fill up lc(i,1-2) as well as tc(i) 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	parameter (kszj=40000)
	parameter(nsize=240000)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa4/tau(kszj),tlco(kszj,4)
      	COMMON/SA2/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)

	dimension lc(nsize,5),tc(nsize),tw(nsize)
	dimension dr(3),db(3),pi(4),pj(4),vi(3),vj(3)
	dimension ri(4),rj(4),rfi(4),rfj(4),b(3)
       	 pel=p(l,4)
	pel1=p(l1,4)
	if(pel.lt.1.e-20)pel=1.e-20   ! 041204
	if(pel1.lt.1.e-20)pel1=1.e-20   ! 041204
	pi(4)=pel
	pj(4)=pel1
	do i=1,3
c	write(*,*)'p,psa=',p(l,i)
c	write(*,*)'p1,psa1=',p(l1,i)
	pi(i)=p(l,i)
	pj(i)=p(l1,i)
	b(i)=(pi(i)+pj(i))/(pi(4)+pj(4))
	enddo
	ilo=0
	call lorntz(ilo,b,pi,pj)
c	perform Lorentz transf. to CMS frame for momentum.
	bta=dsqrt(b(1)**2+b(2)**2+b(3)**2)
c	if boost is too violent,put particles on mass shell by hand.
	if(bta.gt.0.99999d+0)then
	kl=k(l,2)
	kl1=k(l1,2)
	klab=iabs(kl)
	kl1ab=iabs(kl1)
	bmi=pmas(pycomp(kl),1)
	bmj=pmas(pycomp(kl1),1)
	pi(4)=sqrt(bmi**2+pi(1)**2+pi(2)**2+pi(3)**2)
	pj(4)=sqrt(bmj**2+pj(1)**2+pj(2)**2+pj(3)**2)
	endif
	ss=pi(4)+pj(4)
c	do not pair into the collision list if the threshold is too small.
c	if(((klab.eq.2211.or.klab.eq.2112).and.
c     &	(kl1ab.eq.2112.or.kl1ab.eq.2212)).and.ss.le.parp(2))goto 10

	do i=1,4
	ri(i)=v(l,i)
	rj(i)=v(l1,i)
c	write(*,*)'v,vsa=',v(l,i)
c	write(*,*)'v1,vsa1=',v(l1,i)
	enddo
cc	ri(4)=time
cc	rj(4)=time
	call lorntz(ilo,b,ri,rj)
c	perform Lorentz transf. to CMS frame for coordinate.
	rb=0.
	bb=0.
	rr=0.
	rtai=0.
	kflag=0
	do ik=1,3
	vi(ik)=pi(ik)/pi(4)
	vj(ik)=pj(ik)/pj(4)
	enddo

	do i=1,3
	rfi(i)=v(l,i)+(tau(l)-time)*(p(l,i)/p(l,4))
	rfj(i)=v(l1,i)+(tau(l1)-time)*(p(l1,i)/p(l1,4))
	enddo
	rfi(4)=tau(l)
	rfj(4)=tau(l1)
	call lorntz(ilo,b,rfi,rfj)
c	gamli=p(l,4)/p(l,5)
c	gamlj=p(l1,4)/p(l1,5)
 	ctaui=rfi(4)
  	ctauj=rfj(4)
 	tcol=ctaui
	if(ctaui.lt.ctauj)tcol=ctauj
 	do ik=1,3
 	db(ik)=(vi(ik)-vj(ik))*tcol
	dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))+db(ik)
	rtai=rtai+dr(ik)*dr(ik)	
	enddo
  	dot=0.
	do ik=1,3
	dot=dr(ik)*pi(ik)+dot
	enddo
c	dot=-1
	if(dot.ge.0.)then
	kflag=1
	if(tcol.le.ri(4) )goto 10
	if(tcol.le.rj(4) )goto 10	
	else
	rtai=0.
	do ik=1,3
	dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))
	db(ik)=vi(ik)-vj(ik)
	rb=rb+dr(ik)*db(ik)
	bb=bb+db(ik)*db(ik)
	rr=rr+dr(ik)*dr(ik)
	enddo
	if(bb .le. 1.e-10)goto 10
	tcol=0.-rb/bb
	if(tcol.le.ri(4) )goto 10
	if(tcol.le.rj(4) )goto 10
	if(tcol-ctaui .le. 0.)goto 10
	if(tcol-ctauj .le. 0.)goto 10
c	for collision to occur,time must one step ahead
cTai
	do ik=1,3
	dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))+tcol*db(ik)
	rtai=rtai+dr(ik)*dr(ik)
	enddo
c	gamai=pi(4)/pmas(pycomp(k(l,2)),1)
c	gamaj=pj(4)/pmas(pycomp(k(l1,2)),1)

C TAIAN

c	when collision happens,particles should already be produced
c	we give a zero formation time for particles produced after
c       calling 'pythia'
	sg1=rr+tcol*rb
c		write(*,*)'sar',sg1
ctai
	endif
	sg=rtai
c		write(*,*)'tair',sg
c	if(sg1.lt.0.and.iii.le.50)then
c	write(*,*)'sar,tair=',sg1,sg
c	dmin=0.
c	tcol=-rr/rb
c	goto 20
c	endif

	dmin=sqrt(sg)
c	calculate the interaction distance between particles l & l1.
20	call intdis(l,l1,ss,rsig)
c	distance between the two particles should be smaller than rsig
	if(dmin.gt.rsig)goto 10
c	move along Newton trajectory in CMS
	do ik=1,3
	ri(ik)=ri(ik)+vi(ik)*(tcol-ri(4))
	rj(ik)=rj(ik)+vj(ik)*(tcol-rj(4))
	enddo
	ri(4)=tcol
	rj(4)=tcol
	ilo=1

	call lorntz(ilo,b,ri,rj)

c	transform back to Lab.
	tcol1=ri(4)
	tcol2=rj(4)
	if(kflag.eq.0)then
	if(tcol1-tau(l).lt.0.) goto 10
	if(tcol2-tau(l1).lt.0.) goto 10
	else
	if(tcol1-tau(l).lt.-1.E-4) goto 10
	if(tcol2-tau(l1).lt.-1.E-4) goto 10
	endif
	if(ri(4).gt.rj(4)) ri(4)=rj(4)
	tcol=ri(4)
	if(tcol.le.time)goto 10
c	collision happens in the future
c	if(ifram.eq.0)coor(3)=coor(3)+rnt
	do i=1,3
	ri(i)=v(l,i)+p(l,i)*(tcol-time)/pel-coor(i)
	rj(i)=v(l1,i)+p(l1,i)*(tcol-time)/pel1-coor(i)
	enddo
	rri=sqrt(ri(1)*ri(1)+ri(2)*ri(2)+ri(3)*ri(3))
	rrj=sqrt(rj(1)*rj(1)+rj(2)*rj(2)+rj(3)*rj(3))
c	the rnt in rao*max(rnt,rnp)+rnt is due to the fact that
c        we could not know the postion of the mass-center in the future.
	rrr=rao*max(rnt,rnp)
	if(ifram.eq.0)rrr=rao*max(rnt,rnp)+rnt
c	if(abs(k(l,2)).gt.1000.or.abs(k(l1,2)).gt.1000)rrr=1.E+10*rrr
	if(rri.gt.rrr)goto 10
	if(rrj.gt.rrr)goto 10
c	particles under consideration must be still within considered region
c	 when the collision happens
	if(tcol.le.rrr)goto 18   ! 041204
        return   ! 041204
18	tc(icp)=tcol
	lc(icp,1)=l
	lc(icp,2)=l1
c	write(9,*)'tcolij l,l1,icp,tcol=',l,l1,icp,tcol
10	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pauli(ii,ppaul)
c	calculate the unoccupation probability (ppaul) of particle ii
c	 in 'PYJETS'
	parameter(kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/saf/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	dimension rkk(kszj,4),pkk(kszj,4),rr(4),pp(4),b(3)
	kf=k(ii,2)   ! new
	xxii=v(ii,1) 
	yyii=v(ii,2) 
	zzii=v(ii,3) 
	ttii=v(ii,4)
	pxii=p(ii,1) 
	pyii=p(ii,2) 
	pzii=p(ii,3)
	eeii=p(ii,4) 
	b(1)=pxii/eeii
	b(2)=pyii/eeii
	b(3)=pzii/eeii
c	pick up the partons with same flavor as ii from 'PYJETS' and 
c	 'saf' 
	nkk=0 
c	the new produced partons, except ii, in 'PYJETS' are also to be
c	 considered
	do i=1,n   ! loop over new
	if(i.eq.ii)goto 100
	kfi1=k(i,2)
        if(kfi1.eq.kf)then
	nkk=nkk+1
	do j=1,4
	rkk(nkk,j)=v(i,j)
	pkk(nkk,j)=p(i,j)
	enddo
	endif
100	enddo
	do i1=1,nsa   ! loop over old
        kfi1=ksa(i1,2)
        if(kfi1.eq.kf)then
	nkk=nkk+1
        do j=1,4
        rkk(nkk,j)=vsa(i1,j)
        pkk(nkk,j)=psa(i1,j)
        enddo
	endif
        enddo
c	boost to the rest frame of ii
        ilo=0
        do 200 j2=1,nkk
	do j1=1,4
        rr(j1)=rkk(j2,j1)
        pp(j1)=pkk(j2,j1)
        enddo
        call lorntz(ilo,b,rr,pp)
        do j1=1,4
        rkk(j2,j1)=rr(j1)
        pkk(j2,j1)=pp(j1)
        enddo
200	enddo
	rr(1)=xxii
	rr(2)=yyii
	rr(3)=zzii
	rr(4)=ttii
	call lorntz(ilo,b,rr,rr)
	xxii=rr(1)	
	yyii=rr(2)	
	zzii=rr(3)	
	ttii=rr(4)	
c	calculate the number of partons occupied in or on the surface of 
c	 six dimension cub (around ii): (dr*dp)**3=h**3, dr*dp=h, h=1.24 
c	 GeV*fm/c, if dr=1.0 fm then dp=1.24 GeV/c
	anq=0   ! statistics of ocuupation number in (dr*dp)**3=h**3
	do i1=1,nkk   
	dxx=xxii-rkk(i1,1)
	dyy=yyii-rkk(i1,2)
	dzz=zzii-rkk(i1,3)
c	following three staments for without boost
cc	dpx=pxii-pkk(i1,1)
cc	dpy=pyii-pkk(i1,2)
cc	dpz=pzii-pkk(i1,3)
c       following three staments for with boost
	dpx=pkk(i1,1)
	dpy=pkk(i1,2)
	dpz=pkk(i1,3)
	dxx=abs(dxx)
	dyy=abs(dyy)
	dzz=abs(dzz)
	dpx=abs(dpx)
	dpy=abs(dpy)
	dpz=abs(dpz)
	if(dxx.le.0.5.and.dyy.le.0.5.and.dzz.le.0.5.and.
     c	 dpx.le.0.62.and.dpy.le.0.62.and.dpz.le.0.62)anq=anq+1.
	enddo
	proba=anq/6.
c	6=2*3, spin and colour degeneracies of quark (antiquark)
	ppaul=1.-proba
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_gam(ii)   ! 240209
c       move particles with flavor code ii (ii='22' or '44' or '55' or '66')
c        from 'pyjets' to 'sgam'
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        jb=0
201     do i1=jb+1,n
        kf=k(i1,2)
        if(kf.ne.ii)then
        jb=jb+1
        goto 202
        endif
        ngam=ngam+1
        do i2=1,5
        kgam(ngam,i2)=k(i1,i2)
        pgam(ngam,i2)=p(i1,i2)
        vgam(ngam,i2)=v(i1,i2)
        enddo
        if(i1.eq.n)then
        n=n-1
        goto 203
        endif
c       move particle list 'pyjets' one step downward from i1+1 to n
        do j=i1+1,n
        j1=j-1
        do jj=1,5
        k(j1,jj)=k(j,jj)
        p(j1,jj)=p(j,jj)
        v(j1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 201
202     enddo
203     continue
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updpip(l,l1,icp,lc,tc,tw,time,nup,iii)
c       update particle list after calling 'pythia' and 
c        truncate collision list correspondingly.
c	nup is the number of particles from 'pythia' after filtor 
c	 and needs to put in particle list (sa2)   ! 060813
c	iii : current collision number
        parameter(kszj=40000)
        parameter(nsize=240000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/ctllist/nctl,noinel(600),nctl0,noel
	common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio   ! 060813
c060813	ipyth: stord the order number in particle list (sa2) of 'pythia' 
c060813	 particles  
        dimension lc(nsize,5),tc(nsize),tw(nsize)
	dimension peo(4)
	do m=1,2000
	ipyth(m)=0
	enddo
c        n1=nwhole-n
c       put the 'pythia' particles into particle list
c        (i.e. update particle list in inelastic scattering case)
c        & turncate collision list correspondingly.
c	write(9,*)'be. updpip nsa=',nsa
c	do i1=1,nsa
c	write(9,*)'i1,ksa(i1,2)=',i1,ksa(i1,2)
c	enddo
	ll=l
	ll1=l1
	if(nup.eq.0)goto 900   ! 221110
        do 500 i=1,nup
	kf=k(i,2)
	kfab=iabs(kf)   ! 060813 210214
        do 600 j=1,kfmax
	if(kf.ne.kfaco(j))goto 600
        jj=numb(j)+1
c       update the particle list.
	do m=nsa,jj,-1
        mm=m+1
        ksa(mm,2)=ksa(m,2)
        ksa(mm,1)=1
        ksa(mm,3)=ksa(m,3)
        do m1=1,5
        psa(mm,m1)=psa(m,m1)
	vsa(mm,m1)=vsa(m,m1)
	enddo
        ishp(mm)=ishp(m)
        tau(mm)=tau(m)
	enddo
        if(ll.ge.jj)ll=ll+1
        if(ll1.ge.jj)ll1=ll1+1
c       give proper values to particle jj.
        ksa(jj,2)=kf
        ksa(jj,1)=1
        ksa(jj,3)=0
        do m=1,5
        psa(jj,m)=p(i,m)
	vsa(jj,m)=v(i,m)
	enddo
        ishp(jj)=1
        tau(jj)=time+t0*p(i,4)/p(i,5)
c	the values of 'ishp' and 'tau' for particles from 'pythia' 
c	 are given here, the proper formation time of 'pythia' particle 
c	 is assume to be equal to t0 fm/c, except nucleon and j/psi
	if(kf.eq.2212 .or. kf.eq.2112)then
	tau(jj)=time+t0*p(i,4)/p(i,5)*taup
        elseif(kf.eq.443.or.kf.eq.30443)then
        tau(jj)=time+t0*p(i,4)/p(i,5)*taujp
        endif
	do m=j,kfmax
        numb(m)=numb(m)+1
        enddo
	ipyth(i)=jj
	do m=1,2000
	ipym=ipyth(m)
	if(ipym.gt.jj)ipyth(m)=ipym+1
	enddo
c       update the values of lc(m,1-2) when they are .ge. jj
	do m=1,nctl
        lc1=lc(m,1)
        if(lc1.ge.jj)lc(m,1)=lc1+1
        lc2=lc(m,2)
        if(lc2.ge.jj)lc(m,2)=lc2+1
        enddo
	if(i.gt.2)goto 1000
c	jj is instead of ll all over the collision list for the moment and 
c	 throw away those collision pairs in 'updtlp'
        do m=1,nctl
        lc1=lc(m,1)
        lc2=lc(m,2)
        if(lc1.eq.ll)lc(m,1)=jj
        if(lc2.eq.ll)lc(m,2)=jj
        enddo
1000	goto 200
600     continue   ! will not come to here for pA, Ap, and AA due to nup
c060813
	if(ipden.ge.2 .and. (kfab.ge.11.and.kfab.le.16))then   ! 120214
c	put lepton on last position in particle list (sa2) 
	do m=1,5
	ksa(nsa+1,m)=k(i,m)
	psa(nsa+1,m)=p(i,m)
	vsa(nsa+1,m)=v(i,m)
        enddo
        ishp(nsa+1)=1
        tau(nsa+1)=time+t0*p(i,4)/p(i,5)
        ipyth(i)=nsa+1
	nsa=nsa+1
	goto 500
	endif
c060813
200	continue
	nsa=nsa+1
c	if(iii.eq.1)then   !!
c	write(9,*)'600, i,jj,ll,ll1,=',i,jj,ll,ll1   !!
c	write(9,*)'numb=',(numb(i1),i1=1,kfmax)   !!
c	n2=numb(2)   !!
c	call prt(n2)   !!  
c        do i1=1,nctl   !!
c        write(9,*)'updpip i1,lci,lcj,t=',i1,lc(i1,1),lc(i1,2),tc(i1)   !!
c        enddo   !!
c	endif   !!
	if(i.gt.2)goto 500
        ll2=ll
        ll=ll1
        ll1=ll2
500	continue
c	write(9,*)'updpip af. put in nsa=',nsa
c	do i1=1,nsa
c	write(9,*)'i1,ksa(i1,2)=',i1,ksa(i1,2)
c	enddo
900	l=ll
        l1=ll1
c       take out the scattering particles from particle list (i.e.
c        update particle list in inelastic scattering case) &
c        truncate the collision list correspondingly.
        kf1=ksa(l,2)
        kf2=ksa(l1,2)
c	write(9,*)'l,l1,kf1,kf2=',l,l1,kf1,kf2   !!
	kf=kf1
        ll=l
        do 700 i=1,2
c011210
        if(ll.eq.nsa)then   !
        do i1=1,kfmax
	if(kf.ne.kfaco(i1))goto 400
        do m=i1,kfmax
        numb(m)=numb(m)-1
        enddo
        if(i1.gt.1)then
        numba=numb(i1)
        do m=1,i1-1
        if(numb(m).eq.numba)numb(m)=numb(m)-1
        enddo
        endif
        goto 100
400     enddo
        endif   !
c011210
        do j=ll+1,nsa
        jj=j-1
        ksa(jj,2)=ksa(j,2)
        ksa(jj,1)=1
        ksa(jj,3)=ksa(j,3)
        do m=1,5
        psa(jj,m)=psa(j,m)
	vsa(jj,m)=vsa(j,m)
	enddo
        ishp(jj)=ishp(j)
	tau(jj)=tau(j)
	enddo
        do m=1,nctl
        lc1=lc(m,1)
        lc2=lc(m,2)
        if(lc1.gt.ll)lc(m,1)=lc1-1
        if(lc2.gt.ll)lc(m,2)=lc2-1
        enddo
        do 800 j=1,kfmax
        if(kf.ne.kfaco(j))goto 800
        do m=j,kfmax
        numb(m)=numb(m)-1
        enddo
        goto 100
800     continue
100     continue
        nsa=nsa-1
        if(l1.gt.ll)l1=l1-1
	do m=1,2000
	ipym=ipyth(m)
	if(ipym.gt.ll)ipyth(m)=ipym-1
	enddo
c	if(iii.eq.1)then   !!
c	write(9,*)'after taking out, nsa=',nsa   !!
c	write(9,*)'numb=',(numb(i1),i1=1,kfmax)   !!
c	n2=numb(2)   !!
c	call prt(n2)   !!
c        do i1=1,nctl   !!
c        write(9,*)'updpip i1,lci,lcj,t=',i1,lc(i1,1),lc(i1,2),tc(i1)   !!
c        enddo   !!
c	endif   !!
        if(i.eq.2)goto 700
        ll=l1
        kf=kf2
700     continue
c	write(9,*)'updpip af. take out nsa=',nsa
c	do i1=1,nsa
c	write(9,*)'i1,ksa(i1,2)=',i1,ksa(i1,2)
c	enddo
	return
        end



c************************************************************************
        subroutine updtlp(ic,jc,time,lc,tc,tw,nup,iii,kjp21)
c	update collision list after calling 'pythia'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000,KSZ1=30)
        parameter(nsize=240000)
c        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa14/ipyth(2000),idec(2000),iwide(2000)
c010530        common/sa19/kji   ! 16/09/99
        common/ctllist/nctl,noinel(600),nctl0,noel
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        dimension lc(nsize,5),tc(nsize),tw(nsize)
c	ipyth: store line number of particles from calling 'pythia' 
c       loop over old colliding pairs
        j=0
        do i=1,nctl
        i1=lc(i,1)
        j1=lc(i,2)
        if(i1.eq.ic .or. i1.eq.jc)goto 400
        if(j1.eq.ic .or. j1.eq.jc)goto 400
        if((tc(i)-time).le.ddt) goto 400
c       through away the pair whih tc<= time
        j=j+1
        tc(j)=tc(i)
        tw(j)=tw(i)
        do m=1,5
        lc(j,m)=lc(i,m)
        enddo
400     continue
	enddo
        do i=j+1,nctl+1
        tc(i)=0.0
        tw(i)=0.0
        do m=1,5
        lc(i,m)=0
        enddo
        enddo
c	write(9,*)'updtlp af. loop over old colliding pairs j=',j
c	do i1=1,j
c	write(9,*)'i1,lci,lcj,t=',i1,lc(i,1),lc(i,2),tc(i)
c	enddo

        nctl=j+1
c060813	loop over particle list for each generated particle from pythia
	if(nup.eq.0)goto 700   ! 221110
	m2=numb(2)
c	m7=numb(7)  
	do j11=1,nup
	j1=ipyth(j11)
	ksaj1=ksa(j1,2)
	ksaab=iabs(ksaj1)   ! 060813 120214
	if(j1.le.m2.or.(ipden.ge.2.and.(ksaab.ge.11.and.ksaab.le.16)))
     c	 goto 301 !221110 060813 120214
c060813	consider only the reinteraction among nucleons & nucleon with 
c120214	 lepton
	goto 300
c060813	loop over particle list
301	mm=m2   
	do i=1,mm   
c060813	120214 consider only the reinteraction of j11 with nucleons 
        if(nctl.gt.nsize)stop 30000
c010600
	do j22=1,nup
	j2=ipyth(j22)
	if(i.eq.j2)goto 600
	enddo
c010600
602     i1=i
        iflag=0
        call rsfilt(j1,i1,iflag)
c	write(9,*)'updtlp j1,i1,kj1,ki1,iflag=',
c     c	 j1,i1,ksa(j1,2),ksa(i1,2),iflag  
        if(iflag.eq.0)goto 100
        tc(nctl)=0.0
        call tcolij(i1,j1,time,nctl,lc,tc,tw)
c	write(9,*)'af. tcolij nctl,tc(nctl)=',nctl,tc(nctl)
        if(tc(nctl).gt.1.0e-7) nctl=nctl+1
100     continue
600     enddo   ! loop for i
300	enddo   ! loop for j11
700     if(tc(nctl).le.1.e-7) nctl=nctl-1
        do i=nctl+1,nsize
        do m=1,5
        lc(i,m)=0
        enddo
        tc(i)=0.
        tw(i)=0.
        enddo
        return
        end



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	subroutine coelas(ic,jc,eij,pi,pj)
c	perform elastic scattering
	parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/SA2/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c       note the name of the arraies in 'sa2' in this subroutine
	dimension pi(4),pj(4)
	iic=k(ic,2)
	jjc=k(jc,2)
	d=3.65*(eij-pmas(pycomp(iic),1)-pmas(pycomp(jjc),1))
	if(d.lt.1.e-10)return
	pt=0.2
	a=min(10.3,1./(1.12*pt)/(1.12*pt))
	d6=d**6
	b=d6*a/(1.+d6)
	if(b.lt.1.e-20)then
	b=1.e-20
	endif
	pm2=pi(1)**2+pi(2)**2+pi(3)**2
	pm=sqrt(pm2)
	t0=-4.*pm2
	if(abs(t0).lt.1.e-20)then
	cctas=1.
	goto 100
	endif
	cc=pyr(1)
	if(abs(b*t0).lt.0.0001)then
	abt=1.
c	elseif(b*t0.lt.-50.)then
c	abt=0.
	else
	abt=dexp(dmax1(-7.0D2,dble(b*t0)))
	endif
	tt1=dlog(cc+(1.-cc)*abt)
	if(abs(tt1).lt.1.e-30 .and. b.le.1.e-20)then
	cctas=1.
	goto 100
	endif
	tt=tt1/b
	if(abs(tt).lt.1.e-20)then
	cctas=1.
	goto 100
	endif
	cctas=1.-tt*2./t0
	if(abs(cctas).gt.1.)then
	cctas=sign(1.d0,cctas)   ! 250910
	endif
100	continue
	sctas=sqrt(1.-cctas**2)
	fis=2.*3.1416*pyr(1)
	cfis=cos(fis)
	sfis=sin(fis)
	call rotate(cctas,sctas,cfis,sfis,pm,pi,pj)
	return
	end


c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	subroutine updple(ic,jc,b,pi,pj,time)
c	update particle list for elastic scattering
	parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/SA2/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c       note the name of the arrays in 'sa2' in this subroutine
	common/sa4/tau(kszj),tlco(kszj,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	dimension pi(4),pj(4),b(3)
c	write(9,*)'updple be. Lorentz'
c	write(9,*)'pi=',(pi(i1),i1=1,4)
c	write(9,*)'pj=',(pj(i1),i1=1,4)
	ilo=1
c	ilo=1 for inverse Lorentz transformation
	call lorntz(ilo,b,pi,pj)
	do i=1,4
	p(ic,i)=pi(i)
	p(jc,i)=pj(i)
	enddo
	return
	end



c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	subroutine updatl(ic,jc,time,lc,tc,tw,winel,iii,kjp21)!010530
c	update collision time list for usual scattering
	parameter (kszj=40000,KSZ1=30)
	parameter(nsize=240000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/ctllist/nctl,noinel(600),nctl0,noel
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	dimension lc(nsize,5),tc(nsize),tw(nsize)
	integer winel

c	loop over old colliding pairs
	j=0
	do i=1,nctl
	i1=lc(i,1)
	j1=lc(i,2)
c	ia=(i1-ic)*(j1-jc)*(i1-jc)*(j1-ic)
c	if(ia.eq.0) goto 400
	if(i1.eq.ic .or. i1.eq.jc)goto 400
	if(j1.eq.ic .or. j1.eq.jc)goto 400
	if((tc(i)-time).le.ddt) goto 400
c	through away the pair whih tc<= time
	j=j+1
	tc(j)=tc(i)
	tw(j)=tw(i)
	do m=1,5
	lc(j,m)=lc(i,m)
	enddo
400	continue
	enddo
	do i=j+1,nctl+1
	tc(i)=0.0
	tw(i)=0.0
	do m=1,5
	lc(i,m)=0
	enddo
	enddo
c	write(9,*)'updatl af. loop over old colliding pairs j=',j
c	do i1=1,j
c	write(9,*)'i1,lci,lcj,t=',i1,lc(i1,1),lc(i1,2),tc(i1)
c	enddo

	nctl=j+1
c	loop over particle list
	m2=numb(2)
c	m7=numb(7) 
c060813	120214 consider only the reinteraction among nucleons or 
c	lepton with nucleons
	j1=ic
	do ik=1,2
	ksaj1=ksa(j1,2)
	ksaab=iabs(ksaj1)   ! 060812 120214
        if(j1.le.m2.or.(ipden.ge.2.and.(ksaab.ge.11.and.ksaab.le.16)))
     c	 goto 301   ! 060813 120214
        goto 300
301	mm=m2
	do i=1,mm
	if(j1.eq.ic .and. i.eq.jc)goto 600 
	if(j1.eq.jc .and. i.eq.ic)goto 600
c	forbiden scattered particles colliding with each other
	if(nctl.gt.nsize)then
        write(MSTU(11),*)'size of array "nsize" needs to be extended'
        write(MSTU(11),*)'error is serious,stop running'
        stop 30000
        endif
 	i1=i
	iflag=0
	call rsfilt(j1,i1,iflag)
c	write(9,*)'updatl j1,i1,kj1,ki1,iflag=',
c     c   j1,i1,ksa(j1,2),ksa(i1,2),iflag
	if(iflag.eq.0)goto 100
	tc(nctl)=0.0
	call tcolij(i1,j1,time,nctl,lc,tc,tw)
c	write(9,*)'af. tcolij nctl,tc(nctl)=',nctl,tc(nctl)
        if(tc(nctl).gt.1.0e-7) nctl=nctl+1
100	continue		
600	enddo
300	if(ik.eq.2)goto 500
	j1=jc
500	enddo
700	if(tc(nctl).le.1.e-7) nctl=nctl-1
	do i=nctl+1,nsize
	do m=1,5
	lc(i,m)=0
	enddo
	tc(i)=0.
	tw(i)=0.
	enddo
	return
	end



c******************************************************************
	subroutine rsfilt(l,l1,iflag)
c       play the role of first range filter and guarantee the collision list 
c        is composed according to the entrance channels of considered 
c        inelastic reactions
c       subroutine intdis plays the role of second range filter
c       collision pairs not interested can not filter through both of rsfilt 
c        and intdis
	parameter (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/SA2/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	m2=numb(2)
	m4=numb(4)
	kl=k(l,2)
	kl1=k(l1,2)
	klab=iabs(kl)   ! 060813 120214
	kl1ab=iabs(kl1)   ! 060813 120214
	if(l.eq.l1) goto 10
	if(ishp(l).eq.0.or.ishp(l1).eq.0) goto 10

c	constraints on the direct reactions
	if(kl.eq.211 .and. (kl1.eq.-211 .or. kl1.eq.111 .or.
     &	 abs(kl1).eq.2212 .or. abs(kl1).eq.2112 .or. kl1.eq.
     &   3112 .or. kl1.eq.-3122 .or. kl1.eq.-3222 .or. kl1
     &   .eq.-3212 .or. kl1.eq.3212 .or. kl1.eq.3122 .or. kl1.eq.3312
     &   .or. kl1.eq.-3322))goto 11
	if(kl1.eq.211 .and. (kl.eq.-211 .or. kl.eq.111 .or.
     &	 abs(kl).eq.2212 .or. abs(kl).eq.2112 .or. kl.eq.
     &   3112 .or. kl.eq.-3122 .or. kl.eq.-3222 .or. kl
     &   .eq.-3212 .or. kl.eq.3212 .or. kl.eq.3122 .or. kl.eq.3312
     &   .or. kl.eq.-3322))goto 11
	if(kl.eq.-211 .and. (kl1.eq.111 .or.
     &	 abs(kl1).eq.2212 .or. abs(kl1).eq.2112 .or. kl1.eq.
     &   -3112 .or. kl1.eq.3122 .or. kl1.eq.3222 .or. kl1
     &   .eq.3212 .or. kl1.eq.-3212 .or. kl1.eq.-3122 .or. kl1.eq.
     &   -3312 .or. kl1.eq.3322))goto 11
	if(kl1.eq.-211 .and. (kl.eq.111 .or.
     &	 abs(kl).eq.2212 .or. abs(kl).eq.2112 .or. kl.eq.
     &   -3112 .or. kl.eq.3122 .or. kl.eq.3222 .or. kl
     &   .eq.3212 .or. kl.eq.-3212 .or. kl.eq.-3122 .or. kl.eq.
     &   -3312 .or. kl.eq.3322))goto 11
	if(kl.eq.111 .and. (kl1.eq.111 .or. abs(kl1).eq.2212
     &   .or. abs(kl1).eq.2112 .or. abs(kl1).eq.3122 .or.
     &   abs(kl1).eq.3112 .or. abs(kl1).eq.3212 .or. abs(kl1).eq.3222
     &   .or. abs(kl1).eq.3312 .or. abs(kl1).eq.3322))goto 11
	if(kl1.eq.111 .and. (kl.eq.111 .or. abs(kl).eq.2212
     &   .or. abs(kl).eq.2112 .or. abs(kl).eq.3122 .or.
     &   abs(kl).eq.3112 .or. abs(kl).eq.3212 .or. abs(kl).eq.3222
     &   .or. abs(kl).eq.3312 .or. abs(kl).eq.3322))goto 11
        if(kl.eq.321 .and. (kl1.eq.-2212 .or. kl1.eq.-2112 .or. 
     &   kl1.eq.-3122 .or. kl1.eq.-3222 .or. kl1.eq.-3112
     &	 .or. kl1.eq.-3212 .or.kl1.eq.-3312 .or.kl1.eq.-3322))goto 11
        if(kl1.eq.321 .and. (kl.eq.-2212 .or. kl.eq.-2112 .or. 
     &   kl.eq.-3122 .or. kl.eq.-3222 .or. kl.eq.-3112
     &	 .or. kl.eq.-3212 .or.kl.eq.-3312 .or.kl.eq.-3322))goto 11
        if(kl.eq.-321 .and. (kl1.eq.2212 .or. kl1.eq.2112 .or. 
     &   kl1.eq.3122 .or. kl1.eq.3222 .or. kl1.eq.3112
     &	 .or. kl1.eq.3212 .or.kl1.eq.3312 .or. kl1.eq.3322))goto 11
        if(kl1.eq.-321 .and. (kl.eq.2212 .or. kl.eq.2112 .or. 
     &   kl.eq.3122 .or. kl.eq.3222 .or. kl.eq.3112
     &	 .or. kl.eq.3212 .or.kl.eq.3312 .or. kl.eq.3322))goto 11
	if(kl.eq.311 .and. (kl1.eq.-2212 .or. kl1.eq.-2112 .or. kl1.eq. 
     &   -3122 .or. kl1.eq.-3212 .or. kl1.eq.-3112 .or.kl1.eq.-3222
     &   .or. kl1.eq.-3312 .or. kl1.eq.-3322))goto 11
	if(kl1.eq.311 .and. (kl.eq.-2212 .or. kl.eq.-2112 .or. kl.eq. 
     &   -3122 .or. kl.eq.-3212 .or. kl.eq.-3112 .or.kl.eq.-3222
     &   .or. kl.eq.-3312 .or. kl.eq.-3322))goto 11
	if(kl.eq.-311 .and. (kl1.eq.2212 .or. kl1.eq.2112 .or. kl1.eq. 
     &   3122 .or. kl1.eq.3212 .or. kl1.eq.3112 .or.kl1.eq.3222
     &   .or.kl1.eq.3312 .or.kl1.eq.3322))goto 11
	if(kl1.eq.-311 .and. (kl.eq.2212 .or. kl.eq.2112 .or. kl.eq. 
     &   3122 .or. kl.eq.3212 .or. kl.eq.3112 .or.kl.eq.3222
     &   .or.kl.eq.3312 .or.kl.eq.3322))goto 11

c       constraints on the J/Psi (Psi') induced reactions
	if((kl.eq.443 .or. kl.eq.30443) .and. (kl1.eq.2212 .or. kl1.eq.
     &	 2112 .or. kl1.eq.211 .or. kl1.eq.111 .or. kl1.eq.-211 
     &	.or. kl1.eq.213 .or. kl1.eq.113 .or. kl1.eq.-213))goto 11 
	if((kl1.eq.443 .or. kl1.eq.30443) .and. (kl.eq.2212 .or. kl.eq.
     &	 2112 .or. kl.eq.211 .or. kl.eq.111 .or. kl.eq.-211 
     &	.or. kl.eq.213 .or. kl.eq.113 .or. kl.eq.-213))goto 11
 
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c NN scattering is chosen
	if(kl.eq.2112.and.(kl1.eq.2112.or.kl1.eq.2212))goto 11
	if(kl.eq.2212.and.(kl1.eq.2112.or.kl1.eq.2212))goto 11
c060813	120214 consider interaction between lepton and nucleon
	if((klab.ge.11.and.klab.le.16).and.(kl1.eq.2112.or.kl1.eq.2212))
     c	 goto 11
	if((kl.eq.2112.or.kl.eq.2212).and.(kl1ab.ge.11.and.kl1ab.le.
     c	 16))goto 11 
c060813	120214
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c	constraints on the annihilation reactions
	if(kl.eq.-3212 .and. (kl1.eq.2212 .or. kl1.eq.2112))goto 11
	if(kl.eq.-3122 .and. (kl1.eq.2212.or. kl1.eq.2112))goto 11
	if(kl1.eq.-3212 .and. (kl.eq.2212 .or. kl.eq.2112))goto 11
	if(kl1.eq.-3122 .and. (kl.eq.2212.or. kl.eq.2112))goto 11
	if(kl.eq.-2212 .and. (kl1.eq.2212 .or. kl1.eq.2112))goto 11
	if(kl1.eq.-2212 .and. (kl.eq.2212 .or. kl.eq.2112))goto 11
	if(kl.eq.-2112 .and. (kl1.eq.2212.or. kl1.eq.2112))goto 11
	if(kl1.eq.-2112 .and. (kl.eq.2212.or. kl.eq.2112))goto 11

c	constraints on the reverse reactions
        if(kl.eq.211 .and. (kl1.eq.3112 .or. abs(kl1).eq.3122 .or. kl1
     &   .eq.-3222 .or. abs(kl1).eq.3212 .or. abs(kl1).eq.3312 .or.
     &   abs(kl1).eq.3322 .or. abs(kl1).eq.3334.or.kl1.eq.1114
     &    .or.kl1.eq.2114.or.kl1.eq.2214))goto 11
        if(kl1.eq.211 .and. (kl.eq.3112 .or. abs(kl).eq.3122 .or. kl
     &   .eq.-3222 .or. abs(kl).eq.3212 .or. abs(kl).eq.3312 .or.
     &   abs(kl).eq.3322 .or. abs(kl).eq.3334.or.kl.eq.1114
     &    .or.kl.eq.2114.or.kl.eq.2214))goto 11
        if(kl.eq.-211 .and. (kl1.eq.3222 .or. abs(kl1).eq.3122 .or. kl1
     &   .eq.-3112 .or. abs(kl1).eq.3212 .or. abs(kl1).eq.3312 .or.
     &   abs(kl1).eq.3322 .or. abs(kl1).eq.3334.or.
     &   kl1.eq.2224.or.kl1.eq.2114.or.kl1.eq.2214))goto 11
        if(kl1.eq.-211 .and. (kl.eq.3222 .or. abs(kl).eq.3122 .or. kl
     &   .eq.-3112 .or. abs(kl).eq.3212 .or. abs(kl).eq.3312 .or.
     &   abs(kl).eq.3322 .or. abs(kl).eq.3334.or.
     &   kl.eq.2224.or.kl.eq.2114.or.kl.eq.2214))goto 11
        if(kl.eq.111 .and. (abs(kl1).eq.3112 .or. abs(kl1).eq.3122 .or. 
     &   abs(kl1).eq.3222 .or. abs(kl1).eq.3212 .or. abs(kl1).eq.3312
     &   .or. abs(kl1).eq.3322 .or. abs(kl1).eq.3334.or.kl1.eq.2224.
     &    or.kl1.eq.2114.or.kl1.eq.2214.or.kl1.eq.1114))goto 11
        if(kl1.eq.111 .and. (abs(kl).eq.3112 .or. abs(kl).eq.3122 .or. 
     &   abs(kl).eq.3222 .or. abs(kl).eq.3212 .or. abs(kl).eq.3312
     &   .or. abs(kl).eq.3322 .or. abs(kl).eq.3334.or.kl.eq.2224.
     &    or.kl.eq.2114.or.kl.eq.2214.or.kl.eq.1114))goto 11
        if(kl.eq.321 .and. (kl1.eq.-321 .or. kl1.eq.-311 .or. kl1.eq.
     &   3222 .or. kl1.eq.3212 .or. kl1.eq.3112 .or. kl1.eq.3122 .or. 
     &   kl1.eq.3312 .or. kl1.eq.3322 .or. kl1.eq.3334))goto 11
        if(kl1.eq.321 .and. (kl.eq.-321 .or. kl.eq.-311 .or. kl.eq.
     &   3222 .or. kl.eq.3212 .or. kl.eq.3112 .or. kl.eq.3122 .or. 
     &   kl.eq.3312 .or. kl.eq.3322 .or. kl.eq.3334))goto 11
        if(kl.eq.-321 .and. (kl1.eq.311 .or. kl1.eq.-3222 .or. kl1.eq.
     &   -3212 .or. kl1.eq.-3112 .or. kl1.eq.-3122 .or. kl1.eq.-3312 
     &  .or. kl1.eq.-3322 .or. kl1.eq.-3334))goto 11
        if(kl1.eq.-321 .and. (kl.eq.311 .or. kl.eq.-3222 .or. kl.eq.
     &   -3212 .or. kl.eq.-3112 .or. kl.eq.-3122 .or. kl.eq.-3312
     &   .or. kl.eq.-3322 .or. kl.eq.-3334))goto 11
        if(kl.eq.311 .and. (kl1.eq.-311 .or. kl1.eq.3222 .or. kl1.eq.
     &   3212 .or. kl1.eq.3112 .or. kl1.eq.3122 .or. 
     &   kl1.eq.3312 .or. kl1.eq.3322 .or. kl1.eq.3334))goto 11
        if(kl1.eq.311 .and. (kl.eq.-311 .or. kl.eq.3222 .or. kl.eq.
     &   3212 .or. kl.eq.3112 .or. kl.eq.3122 .or. 
     &   kl.eq.3312 .or. kl.eq.3322 .or. kl.eq.3334))goto 11
        if(kl.eq.-311 .and. (kl1.eq.-3222 .or. kl1.eq.-3212 .or. kl1
     &   .eq.-3112 .or. kl1.eq.-3122 .or. kl1.eq.-3312 .or.
     &   kl1.eq.-3322 .or. kl1.eq.-3334))goto 11
        if(kl1.eq.-311 .and. (kl.eq.-3222 .or. kl.eq.-3212 .or. kl
     &   .eq.-3112 .or. kl.eq.-3122 .or. kl.eq.-3312 .or.
     &   kl.eq.-3322 .or. kl.eq.-3334))goto 11
	if(kl1.eq.2212.and.(kl.eq.1114.or.kl.eq.2114.or.
     &       kl.eq.2214.or.abs(kl).eq.213.or.kl.eq.113))goto 11
	if(kl1.eq.2112.and.(kl.eq.2224.or.kl.eq.2114.or.
     &       kl.eq.2214.or.abs(kl).eq.213.or.kl.eq.113))goto 11
	if(kl.eq.2212.and.(kl1.eq.1114.or.kl1.eq.2114.or.
     &       kl1.eq.2214.or.abs(kl1).eq.213.or.kl1.eq.113))goto 11
	if(kl.eq.2112.and.(kl1.eq.2224.or.kl1.eq.2114.or.
     &       kl1.eq.2214.or.abs(kl1).eq.213.or.kl1.eq.113))goto 11
	goto 10

11	iflag=1
10	continue
	return
	end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	subroutine intdis(l,l1,ss,rsig)
c	calculate interaction distance between particles l and l1.
c	It plays also the role of second range filter
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000)
        COMMON/SA2/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
	rsig=0.
	kl=k(l,2)
	kl1=k(l1,2)
	klab=iabs(kl)   ! 060813 120214
	kl1ab=iabs(kl1)   ! 060813 120214

	if(abs(kl).eq.2212 .or. abs(kl).eq.2112)idpl=1 
c060813	projectile is nucleon

	if(abs(kl).eq.443 .or. abs(kl).eq.30443)idpl=2   ! J/Psi or psi 060813

	if(abs(kl).eq.211 .or. kl.eq.111)idpl=3   ! pion 060813
	if(abs(kl).eq.321 .or. abs(kl).eq.311)idpl=4   ! kaon 060813
	if(abs(kl).eq.3212 .or. abs(kl).eq.3112 .or. abs(kl).eq.3222
     c	 .or. abs(kl).eq.3122 .or. abs(kl).eq.3312 .or. abs(kl).eq.
     c   3322 .or. abs(kl).eq.3334)idpl=5 !\Sigma,\Lambda,\Xi, or \Omega 060813
	if(abs(kl).eq.213 .or. kl.eq.113)idpl=6   ! \rho 060813
	if(kl.eq.1114 .or. kl.eq.2114.or.kl.eq.2214 .or. kl.eq.2224)idpl=7
c060813	\Delta
	if(klab.ge.11.and.klab.le.16)idpl=8   ! lepton 060813 120214
	if(abs(kl1).eq.2212 .or. abs(kl1).eq.2112)idpl1=1

	if(abs(kl1).eq.443 .or. abs(kl1).eq.30443)idpl1=2   ! 98/03/24

	if(abs(kl1).eq.211 .or. kl1.eq.111)idpl1=3
	if(abs(kl1).eq.321 .or. abs(kl1).eq.311)idpl1=4
	if(abs(kl1).eq.3212 .or. abs(kl1).eq.3112 .or. abs(kl1)
     c	.eq.3222 .or. abs(kl1).eq.3122 .or. abs(kl1).eq.3312
     c  .or. abs(kl1).eq.3322 .or. abs(kl1).eq.3334)idpl1=5
	if(abs(kl1).eq.213 .or. kl1.eq.113)idpl1=6
	if(kl1.eq.1114 .or. kl1.eq.2114.or.kl1.eq.2214 .or. kl1.eq.2224)
     c	 idpl1=7
	if(kl1ab.ge.11.and.kl1ab.le.16)idpl1=8   ! 060813 050214

	if(idpl.eq.1 .and. idpl1.eq.1)rsig=ecsnn
	if(idpl.eq.8 .and. idpl1.eq.1)rsig=ecsen   ! 060813
	if(idpl.eq.1 .and. idpl1.eq.8)rsig=ecsen   ! 060813

	
	if(idpl.eq.2 .and. idpl1.eq.1)then
	rsig=ecspsn
	if(kl.eq.30443)rsig=ecsspn
	endif
	if(idpl.eq.1 .and. idpl1.eq.2)then
	rsig=ecspsn
        if(kl1.eq.30443)rsig=ecsspn
        endif

       	if(idpl.eq.2 .and. idpl1.eq.3)then
	rsig=ecspsm
c	write(9,*)'psi-pi kl,kl1,ss, cont. rsig=',kl,kl1,ss,rsig   ! sa
	if(kjp20.eq.0)call cspspi(kl,kl1,ss,rsig,th)
c	write(9,*)'psi-pi ss,rsig,the=',ss,rsig,th   ! sa
        if(kl.eq.30443)then
	rsig=ecsspm
	if(kjp20.eq.0)call cspppi(kl,kl1,ss,rsig,th)
	endif
        endif
	if(idpl.eq.3 .and. idpl1.eq.2)then
	rsig=ecspsm
        if(kjp20.eq.0)call cspspi(kl,kl1,ss,rsig,th)
        if(kl1.eq.30443)then
	rsig=ecsspm
	if(kjp20.eq.0)call cspppi(kl,kl1,ss,rsig,th)
	endif
        endif

	if(idpl.eq.2 .and. idpl1.eq.6)then
	rsig=ecspsm
c	write(9,*)'psi-rhp kl,kl1,ss,cont. rsig=',kl,kl1,ss,rsig   ! sa
	if(kjp20.eq.0)call cspsro(kl,kl1,ss,rsig,th)   
	if(kl1.eq.113)rsig=rsig*1.414
c	in case of rho0, cross section enlarges a factor 2 to 
c	 consider the effect of omega
c	write(9,*)'psi-rhp ss, rsig,the=',ss,rsig,th   ! sa
	if(kl.eq.30443)then
        rsig=ecsspm
	if(kjp20.eq.0)call csppro(kl,kl1,ss,rsig,th)
        if(kl1.eq.113)rsig=rsig*1.414
	endif
	endif
	if(idpl.eq.6 .and. idpl1.eq.2)then
	rsig=ecspsm
	if(kjp20.eq.0)call cspsro(kl,kl1,ss,rsig,th)   
	if(kl.eq.113)rsig=rsig*1.414
        if(kl1.eq.30443)then
        rsig=ecsspm
	if(kjp20.eq.0)call csppro(kl,kl1,ss,rsig,th)   
        if(kl.eq.113)rsig=rsig*1.414
        endif
	endif

	if(idpl.eq.3 .and. idpl1.eq.3)rsig=edipi
	if(idpl.eq.1 .and. idpl1.eq.3)rsig=epin
	if(idpl.eq.3 .and. idpl1.eq.1)rsig=epin
	if(idpl.eq.3 .and. idpl1.eq.5)rsig=epin
	if(idpl.eq.5 .and. idpl1.eq.3)rsig=epin
c	assume the total cross section of (pion)y,((pion)cascade) and  
c	 ((pion)omiga) = (pion)n
	if(idpl.eq.4 .and. (idpl1.eq.1 .or. idpl1.eq.5))rsig=ekn
	if((idpl.eq.1 .or. idpl.eq.5) .and. idpl1.eq.4)rsig=ekn
	if(idpl.eq.4 .and. idpl1.eq.4)rsig=edipi

c	assume the total cross section of ky (k cascade) and (k omiga)
c        = kn
c	assume the total cross section of kk=(pion)(pion)
	if(idpl.eq.1 .and. idpl1.eq.6)rsig=epin
	if(idpl.eq.1 .and. idpl1.eq.7)rsig=ecsnn
	if(idpl.eq.3 .and. idpl1.eq.7)rsig=epin
	if(idpl1.eq.1 .and. idpl.eq.6)rsig=epin
	if(idpl1.eq.1 .and. idpl.eq.7)rsig=ecsnn
	if(idpl1.eq.3 .and. idpl.eq.7)rsig=epin
	if(idpl.eq.1 .and. idpl1.eq.5)rsig=ecsnn
	if(idpl.eq.5 .and. idpl1.eq.1)rsig=ecsnn
c	assume the total cross section of (rho)n=(pion)n
c	assume the total cross section of n(delta)=nn
c	assume the total cross section of (pion)(delta)=(pion)n
c	assume the total cross section of ny, n(cascade) and n(omiga)=nn
	return
	end



c********************************************************************
        subroutine trans
c       'sa2' to 'pyjets' after finish calculation
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (kszj=40000,KSZ1=30)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        do l=1,nsa
        do m=1,5
        k(l,m)=ksa(l,m)
        p(l,m)=psa(l,m)
        v(l,m)=vsa(l,m)
        enddo
        enddo
        n=nsa
        return
        end



c*********************************************************************
	subroutine cspspi(kl,kl1,ss,rsig,the)
c	calculate the interaction distance of J/Psi + pion
c	the corresponding cross section has been renewed   ! 112299
	parameter (ii=37)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	dimension ab(ii),co(ii)
        data ab/0.6536,    0.6688,
     c      0.6840,    0.6992,    0.7143,    0.7295,    0.7447,
     c      0.7599,    0.7750,    0.7902,    0.8054,    0.8206,
     c      0.8357,    0.8509,    0.8661,    0.8813,    0.8965,
     c      0.9116,    0.9268,    0.9420,    0.9572,    0.9723,
     c      0.9875,    1.0027,    1.0179,    1.0331,    1.0482,
     c      1.0634,    1.0786,    1.0938,    1.1089,    1.1241,
     c      1.1393,    1.1545,    1.1696,    1.1848,    1.2000/
        data co/0.7228,    0.9916,
     c      1.1119,    1.1595,    1.1637,    1.1405,    1.0993,
     c      1.0468,    0.9873,    1.0600,    1.0663,    1.0220,
     c      0.9584,    0.8863,    0.8115,    0.7370,    0.6651,
     c      0.5968,    0.5331,    0.4741,    0.4201,    0.3711,
     c      0.3269,    0.2873,    0.2520,    0.2207,    0.1932,
     c      0.1690,    0.1479,    0.1295,    0.1137,    0.1000,
     c      0.0883,    0.0784,    0.0699,    0.0627,    0.0567/
	the=pmas(pycomp(kl),1)+pmas(pycomp(kl1),1)
	the=ss-the
c	write(9,*)'kl,kl1,ss,the=',kl,kl1,ss,the   ! sa
	if(the.lt.ab(1))then
	rsig=0.
	goto 100
	endif
	if(the.gt.ab(ii))then 
	rsig=sqrt(co(ii)*1.2*0.1/3.1416)
c	1.2 is the ratio of total cross section to absorption
c	0.1 is the transfermation factor from mb to fm^2
	goto 100
	endif
	do i=1,ii-1
	a1=ab(i)
	a2=ab(i+1)
	if(the.ge.a1 .and. the.lt.a2)then
	c1=co(i)
	c2=co(i+1)
	rsig=c1+(the-a1)/(a2-a1)*(c2-c1)
	rsig=sqrt(rsig*1.2*0.1/3.1416)
	goto 100
	endif
	enddo
100	return
	end
c**********************************************************************



c**********************************************************************
        subroutine cspsro(kl,kl1,ss,rsig,the)
c       calculate the interaction distance of J/Psi + rho
c       the corresponding cross section has been renewed   ! 112299
        parameter (ii=80)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        dimension ab(ii),co(ii)
        data ab/    0.0100,    0.0251,    0.0401,    0.0552,    0.0703,
     c      0.0853,    0.1004,    0.1154,    0.1305,    0.1456,
     c      0.1606,    0.1757,    0.1908,    0.2058,    0.2209,
     c      0.2359,    0.2510,    0.2661,    0.2811,    0.2962,
     c      0.3113,    0.3263,    0.3414,    0.3565,    0.3715,
     c      0.3866,    0.4016,    0.4167,    0.4318,    0.4468,
     c      0.4619,    0.4770,    0.4920,    0.5071,    0.5222,
     c      0.5372,    0.5523,    0.5673,    0.5824,    0.5975,
     c      0.6125,    0.6276,    0.6427,    0.6577,    0.6728,
     c      0.6878,    0.7029,    0.7180,    0.7330,    0.7481,
     c      0.7632,    0.7782,    0.7933,    0.8084,    0.8234,
     c      0.8385,    0.8535,    0.8686,    0.8837,    0.8987,
     c      0.9138,    0.9289,    0.9439,    0.9590,    0.9741,
     c      0.9891,    1.0042,    1.0192,    1.0343,    1.0494,
     c      1.0644,    1.0795,    1.0946,    1.1096,    1.1247,
     c      1.1397,    1.1548,    1.1699,    1.1849,    1.2000/
        data co/   26.0928,   25.2027,   21.2192,   18.0128,   15.3798,
     c     13.1673,   11.2818,    9.6619,    8.2637,    7.0540,
     c      7.9089,    7.7564,    7.0910,    6.2957,    5.4907,
     c      4.7272,    4.0277,    3.4011,    2.8488,    2.3682,
     c      1.9548,    1.6027,    1.3059,    1.0582,    0.8536,
     c      0.6866,    0.5521,    0.4454,    0.3622,    0.2989,
     c      0.2521,    0.2190,    0.1969,    0.1838,    0.1777,
     c      0.1771,    0.1806,    0.1870,    0.1955,    0.2053,
     c      0.2156,    0.2261,    0.2363,    0.2459,    0.2547,
     c      0.2625,    0.2693,    0.2749,    0.2793,    0.2826,
     c      0.2847,    0.2857,    0.2857,    0.2847,    0.2828,
     c      0.2801,    0.2766,    0.2724,    0.2677,    0.2624,
     c      0.2567,    0.2506,    0.2442,    0.2375,    0.2306,
     c      0.2236,    0.2165,    0.2094,    0.2022,    0.1951,
     c      0.1880,    0.1809,    0.1740,    0.1672,    0.1605,
     c      0.1539,    0.1475,    0.1413,    0.1353,    0.1294/
        the=pmas(pycomp(kl),1)+pmas(pycomp(kl1),1)
        the=ss-the
        if(the.lt.ab(1))then
        rsig=0.
        goto 100
        endif
        if(the.gt.ab(ii))then
        rsig=sqrt(co(ii)*1.2*0.1/3.1416)
c       1.2 is the ratio of total cross section to absorption
c       0.1 is the transfermation factor from mb to fm^2
        goto 100
        endif
        do i=1,ii-1
        a1=ab(i)
        a2=ab(i+1)
        if(the.ge.a1 .and. the.lt.a2)then
        c1=co(i)
        c2=co(i+1)
        rsig=c1+(the-a1)/(a2-a1)*(c2-c1)
        rsig=sqrt(rsig*1.2*0.1/3.1416)
        goto 100
        endif
        enddo
100     return
        end
c**********************************************************************



c*********************************************************************
        subroutine cspppi(kl,kl1,ss,rsig,the)
c       calculate the interaction distance of Psi' + pion
        parameter (ii=173)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        dimension ab(ii),co(ii)
	data ab/    0.0565,    0.0632,    0.0698,    0.0765,    0.0831,
     c    0.0898,    0.0964,    0.1031,    0.1097,    0.1164,
     c    0.1230,    0.1297,    0.1363,    0.1430,    0.1496,
     c    0.1563,    0.1629,    0.1696,    0.1762,    0.1828,
     c    0.1895,    0.1961,    0.2028,    0.2094,    0.2161,
     c    0.2227,    0.2294,    0.2360,    0.2427,    0.2493,
     c    0.2560,    0.2626,    0.2693,    0.2759,    0.2826,
     c    0.2892,    0.2959,    0.3025,    0.3092,    0.3158,
     c    0.3225,    0.3291,    0.3358,    0.3424,    0.3490,
     c    0.3557,    0.3623,    0.3690,    0.3756,    0.3823,
     c    0.3889,    0.3956,    0.4022,    0.4089,    0.4155,
     c    0.4222,    0.4288,    0.4355,    0.4421,    0.4488,
     c    0.4554,    0.4621,    0.4687,    0.4754,    0.4820,
     c    0.4887,    0.4953,    0.5020,    0.5086,    0.5153,
     c    0.5219,    0.5285,    0.5352,    0.5418,    0.5485,
     c    0.5551,    0.5618,    0.5684,    0.5751,    0.5817,
     c    0.5884,    0.5950,    0.6017,    0.6083,    0.6150,
     c    0.6216,    0.6283,    0.6349,    0.6416,    0.6482,
     c    0.6549,    0.6615,    0.6682,    0.6748,    0.6815,
     c    0.6881,    0.6947,    0.7014,    0.7080,    0.7147,
     c    0.7213,    0.7280,    0.7346,    0.7413,    0.7479,
     c    0.7546,    0.7612,    0.7679,    0.7745,    0.7812,
     c    0.7878,    0.7945,    0.8011,    0.8078,    0.8144,
     c    0.8211,    0.8277,    0.8344,    0.8410,    0.8477,
     c    0.8543,    0.8609,    0.8676,    0.8742,    0.8809,
     c    0.8875,    0.8942,    0.9008,    0.9075,    0.9141,
     c    0.9208,    0.9274,    0.9341,    0.9407,    0.9474,
     c    0.9540,    0.9607,    0.9673,    0.9740,    0.9806,
     c    0.9873,    0.9939,    1.0006,    1.0072,    1.0138,
     c    1.0205,    1.0272,    1.0338,    1.0404,    1.0471,
     c    1.0537,    1.0604,    1.0670,    1.0737,    1.0803,
     c    1.0870,    1.0936,    1.1003,    1.1069,    1.1136,
     c    1.1202,    1.1269,    1.1335,    1.1402,    1.1468,
     c    1.1535,    1.1601,    1.1668,    1.1734,    1.1801,
     c    1.1867,    1.1933,    1.2000/
	data co/    6.0942,    7.9450,    8.7050,    8.9397,    8.8716,
     c    8.6166,    8.2434,    7.7964,    7.3055,    6.7919,
     c    6.2703,    5.7519,    5.2444,    4.7537,    4.2839,
     c    3.8379,    3.4177,    3.0245,    2.6588,    2.3210,
     c    2.0107,    2.2416,    4.6256,    5.2157,    5.3807,
     c    5.3276,    5.1470,    4.8893,    4.5862,    4.2591,
     c    3.9227,    3.5875,    3.2610,    2.9484,    2.6534,
     c    2.3784,    2.1250,    1.8939,    1.6853,    1.4991,
     c    1.3348,    1.1915,    1.0683,    0.9642,    0.8778,
     c    0.8081,    0.7536,    0.7131,    0.6854,    0.6692,
     c    0.6632,    0.6663,    0.6773,    0.6952,    0.7189,
     c    0.7475,    0.7801,    0.8159,    0.8540,    0.8938,
     c    0.9346,    2.2641,    4.1154,    4.8313,    5.1627,
     c    5.2807,    5.2649,    5.1609,    4.9984,    4.7977,
     c    4.5733,    4.3358,    4.0929,    3.8506,    3.6130,
     c    3.3832,    3.1635,    2.9554,    2.7598,    2.5774,
     c    2.4084,    2.2527,    2.1100,    1.9800,    1.8621,
     c    1.7557,    1.6602,    1.5750,    1.4991,    1.4320,
     c    1.3730,    1.3212,    1.2761,    1.2370,    1.2031,
     c    1.1740,    1.1491,    1.1279,    1.1097,    1.0943,
     c    1.0811,    1.0698,    1.0600,    1.0514,    1.0438,
     c    1.0367,    1.0302,    1.0238,    1.0175,    1.0111,
     c    1.0045,    0.9975,    0.9901,    0.9822,    0.9737,
     c    0.9646,    0.9549,    0.9445,    0.9334,    0.9217,
     c    0.9094,    0.8964,    0.8828,    0.8686,    0.8539,
     c    0.8386,    0.8230,    0.8068,    0.7904,    0.7735,
     c    0.7564,    0.7391,    0.7216,    0.7039,    0.6861,
     c    0.6683,    0.6505,    0.6327,    0.6150,    0.5974,
     c    0.5800,    0.5627,    0.5457,    0.5290,    0.5125,
     c    0.4963,    0.4805,    0.4650,    0.4499,    0.4352,
     c    0.4209,    0.4070,    0.3936,    0.3806,    0.3681,
     c    0.3560,    0.3444,    0.3333,    0.3226,    0.3124,
     c    0.3027,    0.2935,    0.2847,    0.2764,    0.2686,
     c    0.2612,    0.2542,    0.2477,    0.2416,    0.2360,
     c    0.2307,    0.2258,    0.2213/
        the=pmas(pycomp(kl),1)+pmas(pycomp(kl1),1)
        the=ss-the
        if(the.lt.ab(1))then
        rsig=0.
        goto 100
        endif
        if(the.gt.ab(ii))then
        rsig=sqrt(co(ii)*1.2*0.1/3.1416)
c	1.2 is the ratio of total cross section to absorption
c	0.1 is the transfermation factor from mb to fm^2
        goto 100
        endif
        do i=1,ii-1
        a1=ab(i)
        a2=ab(i+1)
        if(the.ge.a1 .and. the.lt.a2)then
        c1=co(i)
        c2=co(i+1)
        rsig=c1+(the-a1)/(a2-a1)*(c2-c1)
        rsig=sqrt(rsig*1.2*0.1/3.1416)
	goto 100
        endif
        enddo
100     return
        end



c*********************************************************************
        subroutine csppro(kl,kl1,ss,rsig,the)
c       calculate the interaction distance of Psi' + rho
        parameter (ii=80)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        dimension ab(ii),co(ii)
	data ab/     0.0200,    0.0387,    0.0575,    0.0762,    0.0949,
     c    0.1137,    0.1324,    0.1511,    0.1699,    0.1886,
     c    0.2073,    0.2261,    0.2448,    0.2635,    0.2823,
     c    0.3010,    0.3197,    0.3385,    0.3572,    0.3759,
     c    0.3947,    0.4134,    0.4322,    0.4509,    0.4696,
     c    0.4884,    0.5071,    0.5258,    0.5446,    0.5633,
     c    0.5820,    0.6008,    0.6195,    0.6382,    0.6570,
     c    0.6757,    0.6944,    0.7132,    0.7319,    0.7506,
     c    0.7694,    0.7881,    0.8068,    0.8256,    0.8443,
     c    0.8630,    0.8818,    0.9005,    0.9192,    0.9380,
     c    0.9567,    0.9754,    0.9942,    1.0129,    1.0316,
     c    1.0504,    1.0691,    1.0878,    1.1066,    1.1253,
     c    1.1441,    1.1628,    1.1815,    1.2003,    1.2190,
     c    1.2377,    1.2565,    1.2752,    1.2939,    1.3127,
     c    1.3314,    1.3501,    1.3689,    1.3876,    1.4063,
     c    1.4251,    1.4438,    1.4625,    1.4813,    1.5000/
	data co/   34.9004,   23.6715,   18.1408,   14.5462,   11.9041,
     c    9.8312,    8.1433,    6.7397,    5.5588,    4.5603,
     c    3.7153,    3.0020,    2.4032,    1.9043,    1.4930,
     c    1.1583,    0.8907,    0.6814,    0.5224,    0.4067,
     c    0.3278,    0.2800,    0.2580,    0.2573,    0.2736,
     c    0.3034,    0.3436,    0.3913,    0.4441,    0.5002,
     c    0.5576,    0.6152,    0.6716,    0.7259,    0.7775,
     c    0.8256,    0.8700,    0.9104,    0.9464,    0.9782,
     c    1.0056,    1.0287,    1.0477,    1.0627,    1.0738,
     c    1.0813,    1.0853,    1.0862,    1.0841,    1.0793,
     c    1.0720,    1.0625,    1.0510,    1.0376,    1.0226,
     c    1.0063,    0.9887,    0.9700,    0.9505,    0.9302,
     c    0.9093,    0.8879,    0.8662,    0.8442,    0.8221,
     c    0.7999,    0.7777,    0.7556,    0.7336,    0.7118,
     c    0.6902,    0.6689,    0.6479,    0.6272,    0.6068,
     c    0.5867,    0.5668,    0.5472,    0.5277,    0.5082/
        the=pmas(pycomp(kl),1)+pmas(pycomp(kl1),1)
        the=ss-the
        if(the.lt.ab(1))then
        rsig=0.
        goto 100
        endif
        if(the.gt.ab(ii))then
        rsig=sqrt(co(ii)*1.2*0.1/3.1416)
c	1.2 is the ratio of total cross section to absorption
c	0.1 is the transfermation factor from mb to fm^2
        goto 100
        endif
        do i=1,ii-1
        a1=ab(i)
        a2=ab(i+1)
        if(the.ge.a1 .and. the.lt.a2)then
        c1=co(i)
        c2=co(i+1)
        rsig=c1+(the-a1)/(a2-a1)*(c2-c1)
        rsig=sqrt(rsig*1.2*0.1/3.1416)
	goto 100
        endif
        enddo
100     return
        end


c060813 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine crosep(rots,csen)
!!exponential interpolation for ep total cross section (in mbarn) 
!! calculated with herafitter by Xing-Long Li on 10/Dec./2013
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      parameter (rotsMin=2.d0) !!the lower limit of rots (GeV)
      parameter (rotsMax=1003) !!the upper limit of rots (GeV)
      parameter (Ndat=315)
!!sect: ep cross section data for sqrtS range: 2~1003GeV
!!the unit is mbarn
!!when i=1,sqrtS=2.0 GeV; sqrtS=2.0*(1.02)^(i-1) GeV otherwise
      dimension sect(Ndat)
      parameter (sect=[
     & 0.125377D-04,0.136063D-04,0.147199D-04,0.158792D-04,0.170824D-04,
     & 0.183314D-04,0.196228D-04,0.209579D-04,0.223386D-04,0.237574D-04,
     & 0.252260D-04,0.267283D-04,0.282729D-04,0.298683D-04,0.314933D-04,
     & 0.331676D-04,0.348724D-04,0.366252D-04,0.384099D-04,0.402340D-04,
     & 0.421004D-04,0.439997D-04,0.459366D-04,0.479080D-04,0.499176D-04,
     & 0.519633D-04,0.540432D-04,0.561548D-04,0.583026D-04,0.604863D-04,
     & 0.626984D-04,0.649432D-04,0.672178D-04,0.695343D-04,0.718663D-04,
     & 0.742366D-04,0.766461D-04,0.790697D-04,0.815331D-04,0.840216D-04,
     & 0.865379D-04,0.890822D-04,0.916557D-04,0.942586D-04,0.968830D-04,
     & 0.995384D-04,0.102212D-03,0.104924D-03,0.107655D-03,0.110411D-03,
     & 0.113187D-03,0.115993D-03,0.118824D-03,0.121673D-03,0.124552D-03,
     & 0.127461D-03,0.130369D-03,0.133312D-03,0.136294D-03,0.139258D-03,
     & 0.142289D-03,0.145318D-03,0.148358D-03,0.151438D-03,0.154537D-03,
     & 0.157648D-03,0.160778D-03,0.163928D-03,0.167102D-03,0.170296D-03,
     & 0.173515D-03,0.176740D-03,0.179988D-03,0.183258D-03,0.186551D-03,
     & 0.189855D-03,0.193175D-03,0.196515D-03,0.199880D-03,0.203265D-03,
     & 0.206648D-03,0.210060D-03,0.213497D-03,0.216952D-03,0.220408D-03,
     & 0.223888D-03,0.227392D-03,0.230909D-03,0.234442D-03,0.237992D-03,
     & 0.241560D-03,0.245143D-03,0.248743D-03,0.252364D-03,0.256010D-03,
     & 0.259659D-03,0.263324D-03,0.267009D-03,0.270715D-03,0.274444D-03,
     & 0.278182D-03,0.281933D-03,0.285693D-03,0.289480D-03,0.293325D-03,
     & 0.297112D-03,0.300951D-03,0.304824D-03,0.308704D-03,0.312598D-03,
     & 0.316505D-03,0.320433D-03,0.324381D-03,0.328354D-03,0.332347D-03,
     & 0.336349D-03,0.340362D-03,0.344403D-03,0.348459D-03,0.352536D-03,
     & 0.356641D-03,0.360763D-03,0.364893D-03,0.369052D-03,0.373232D-03,
     & 0.377425D-03,0.381629D-03,0.385859D-03,0.390121D-03,0.394418D-03,
     & 0.398725D-03,0.403041D-03,0.407380D-03,0.411739D-03,0.416122D-03,
     & 0.420532D-03,0.424966D-03,0.429429D-03,0.433908D-03,0.438403D-03,
     & 0.442923D-03,0.447466D-03,0.452044D-03,0.456646D-03,0.461274D-03,
     & 0.465905D-03,0.470575D-03,0.475274D-03,0.480000D-03,0.484746D-03,
     & 0.489511D-03,0.494317D-03,0.499153D-03,0.504004D-03,0.508874D-03,
     & 0.513782D-03,0.518727D-03,0.523706D-03,0.528696D-03,0.533711D-03,
     & 0.538757D-03,0.543840D-03,0.548953D-03,0.554107D-03,0.559291D-03,
     & 0.564487D-03,0.569720D-03,0.574975D-03,0.580275D-03,0.585612D-03,
     & 0.590979D-03,0.596380D-03,0.601802D-03,0.607256D-03,0.612757D-03,
     & 0.618283D-03,0.623860D-03,0.629465D-03,0.635094D-03,0.640766D-03,
     & 0.646474D-03,0.652220D-03,0.658001D-03,0.663819D-03,0.669669D-03,
     & 0.675552D-03,0.681488D-03,0.687464D-03,0.693475D-03,0.699515D-03,
     & 0.705602D-03,0.711726D-03,0.717890D-03,0.724097D-03,0.730344D-03,
     & 0.736634D-03,0.742995D-03,0.749397D-03,0.755772D-03,0.762216D-03,
     & 0.768727D-03,0.775280D-03,0.781876D-03,0.788516D-03,0.795202D-03,
     & 0.801928D-03,0.808704D-03,0.815513D-03,0.822366D-03,0.829276D-03,
     & 0.836237D-03,0.843243D-03,0.850306D-03,0.857418D-03,0.864568D-03,
     & 0.871777D-03,0.879027D-03,0.886321D-03,0.893664D-03,0.901059D-03,
     & 0.908518D-03,0.916017D-03,0.923583D-03,0.931201D-03,0.938865D-03,
     & 0.946581D-03,0.954347D-03,0.962169D-03,0.970049D-03,0.977973D-03,
     & 0.985961D-03,0.993995D-03,0.100210D-02,0.101027D-02,0.101850D-02,
     & 0.102677D-02,0.103510D-02,0.104349D-02,0.105192D-02,0.106040D-02,
     & 0.106897D-02,0.107761D-02,0.108635D-02,0.109510D-02,0.110390D-02,
     & 0.111274D-02,0.112163D-02,0.113063D-02,0.113971D-02,0.114886D-02,
     & 0.115806D-02,0.116732D-02,0.117663D-02,0.118600D-02,0.119543D-02,
     & 0.120495D-02,0.121453D-02,0.122419D-02,0.123391D-02,0.124371D-02,
     & 0.125356D-02,0.126348D-02,0.127347D-02,0.128353D-02,0.129367D-02,
     & 0.130387D-02,0.131415D-02,0.132451D-02,0.133493D-02,0.134544D-02,
     & 0.135602D-02,0.136666D-02,0.137737D-02,0.138816D-02,0.139902D-02,
     & 0.140998D-02,0.142099D-02,0.143209D-02,0.144327D-02,0.145453D-02,
     & 0.146587D-02,0.147729D-02,0.148878D-02,0.150035D-02,0.151200D-02,
     & 0.152373D-02,0.153554D-02,0.154745D-02,0.155945D-02,0.157155D-02,
     & 0.158371D-02,0.159591D-02,0.160822D-02,0.162062D-02,0.163313D-02,
     & 0.164572D-02,0.165839D-02,0.167115D-02,0.168400D-02,0.169693D-02,
     & 0.170996D-02,0.172306D-02,0.173628D-02,0.174958D-02,0.176297D-02,
     & 0.177646D-02,0.179003D-02,0.180371D-02,0.181748D-02,0.183135D-02,
     & 0.184530D-02,0.185935D-02,0.187350D-02,0.188775D-02,0.190209D-02
     & ])
!!check if rots is in range of data set,if not in range 2~1003,return -1
      if(rots.lt.rotsMin.or.rots.gt.rotsMax)then
            csen=-1.d0
            return
      endif
!!calculate csen by Interpolation
      x=dlog(rots/2.d0)/dlog(1.02d0)+1.d0
      i=floor(x)
      csen=(sect(i+1)-sect(i))*(x-i)+sect(i)

      return
      end



c060813 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ltof(ii)
c120214	move ii-th particle (lepton) in pyjets to first position
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=40000,ksz1=30)
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
	dimension kk(5),pp(5),vv(5)
        do jj=1,5
        kk(jj)=k(ii,jj)
        pp(jj)=p(ii,jj)
        vv(jj)=v(ii,jj)
        enddo
c	move particle list (pyjets) one step forward from ii-1 to 1
	do j1=ii-1,1,-1
	do jj=1,5
        k(j1+1,jj)=k(j1,jj)
        p(j1+1,jj)=p(j1,jj)
        v(j1+1,jj)=v(j1,jj)
        enddo
        enddo
	do jj=1,5
	k(1,jj)=kk(jj)
	p(1,jj)=pp(jj)
	v(1,jj)=vv(jj)
	enddo
	do i1=1,kfmax
	numbs(i1)=numbs(i1)+1   ! now first particle in pyjets is e-
	enddo
	return
	end



c******************************************************************************
	BLOCK DATA PYCIDATA
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
	COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
	common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
	SAVE /PYCIDAT1/,/PYCIDAT2/
      	DATA KFACOT/2212,2112,-2212,-2112,211,-211,111,-321,-311,
     &        3212,3112,3222,-3212,-3112,-3222,3122,-3122,311,
     &     321,3312,-3312,3322,-3322,3334,-3334,1114,2114,2214,2224,
     &	213,-213,113,443,30443,10441,20443,445,411,-411,421,-421,
     &	4122,4112,4212,4222,223,323,313,413,-413,423,-423,48*0/
      	DATA DISDET/0.5,0.5,0.5,0.5,46*0.,0.5,0.5,0.5,0.5,46*0./
      	DATA ISINELT/384*1,208*0,8*1/  ! with delta and rho
      	DATA KFMAXT/52/
      	DATA PARAM/40.,25.,21.,10.,2.0,0.85,1.0,0.02,0.1,4.0,0.16,0.04,
     &        6.0,3.0,12.,6.,4*0/   ! 060813 	
                  DATA WEIGH/600*1.0/
	data kjp20,vjp20,vjp21,vjp22,vjp23/1,0.3,4.0,1.5,8.0/

	END
C******************************************************************
C...........Main switches and parameters...........................
C\item[KFACOT] flavor order of considered particles
C  \item[DISDET] allowable minimum distance between two
C  particles,=0.5 between two necleons,=0 otherwise
C  \item[ISINELT] switch for i-th inelastic channel
C  =0 closed,=1,opened
C \item[KFMAXT](D=12) KFMAXT kinds of particles are involved in rescattering
C PARAM(1)(D=40.0mb) totle cross-section of nucleon-nucleon 
C PARAM(2)(D=25.0mb)  totle cross-section of pi-nucleon 
C PARAM(3)(D=21.0mb) totle cross-section of K-nucleon 
C PARAM(4)(D=10.0mb)  totle cross-section of pi-pi
C PARAM(5)(D=2.0mb)  cross-section of pi+pi -->K K 
C PARAM(6)(D=0.85) ratio of inelastic cross-section to totle cross-section
C PARAM(7)(D=1.0fm) formation time at rest-frame of particle
C PARAM(8)(D=0.02fm) time accuracy used in hadron cascade
C PARAM(9)(D=0.1) accuracy of four-momentum conservation
C PARAM(10)(D=4.0) size of effective rescattering region is product of 
C  PARAM(10) and radius of target, origin is set on center of target nucleus
C PARAM(11)(D=0.16fm^-3) nucleon density of nucleus
C PARAM(12)(D=0.04 GeV^2/c^2) The <Pt^2> for the Gaussian distribution of 
C	spectator, no used anymore
C PARAM(13)(D=6.0mb) totle cross-section of J/Psi + n
C PARAM(14)(D=3.0mb) totle cross-section of J/Psi + meson
C PARAM(15)(D=12.0mb) totle cross-section of Psi' + n
C PARAM(16)(D=6.0mb) totle cross-section of Psi' + meson
c	kjp20 = 0 : energy dependent cross section
c             = 1 : constant cross section 
c	vjp20 : constant cross section of strangeness production
c	vjp21 : cross section of pion + p to pion + delta
c	vjp22 : cross section of pion + p to rho + p
c	vjp23 : cross section of n + n to n + delta
C@@@@@@@@@@@@@@@@@@@@@  END  @@@@@@@@@@@@@@@@@@@@@@@@
