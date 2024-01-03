	subroutine parini(time_ini,parp21,parp22,win,psno,ijk)   ! 081010
c	generate the partonic initial state for the relativistic nucleus- 
c	 nucleus collision based on 'pythia' 
c	it was composed by Ben-Hao Sa on 04/12/03
c	the intermediate working arraies are in common statement 'sa2'
c       'saf' also consists of intermediate working arraies 
c       'saf' to 'pyjets' after call 'scat'   ! 220110 
c	output message is in 'pyjets' (partons) and 'sbh' (hadrons) 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter(kszj=40000,ksz1=30)
        parameter(nsize=240000)
	double precision bst(4),bzp,bzt,bbb(3),bb(3)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c	those variables in above common blocks are defined in 'jetset'
        COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
	COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)   ! 221203
c	those variables in above common block are defined in 'pythia'
	COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa6/kfmaxi,nwhole
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
c080104
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 220110
	common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21,adiv,gpmax,nnc   !   070417
	common/sa30/vneump,vneumt   ! 241110
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
c080104
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5) 
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
c	iii : number of current event
c	neve : total number of events 
c	bp : impact parameter
c       'sbe': store initial parton confiquration (with diquark) of a A+A
c       'saf': store parton configuration after parton re scattering 
c              (w/o diquark) 
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
c       ecsen: largest collision distance between e- and p   ! 060813
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
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 270312 240412 131212

	dimension peo(4),pi(4),pj(4),xi(4),xj(4)
	dimension inoin(kszj)
        dimension lc(nsize,5),tc(nsize),tw(nsize)

	kpar=0
        knn=0
        kpp=0
        knp=0
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

c270312 initiation of x,y,xy,x^2,y^2 and sump (statistics of the number of
c	 nucleons in overlap region)   ! 131212
        sumx=0.
        sumy=0.
	sumxy=0.   ! 131212
        sumx2=0.
        sumy2=0.
        sump=0.
c270312

c	initiate the nucleus-nucleus collision system
c241110
c       creat the initial particle list (nucleon)
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
	call arrove(i1,1,sumx,sumy,sumxy,sumx2,sumy2,sump)   ! 270312 131212
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
	call arrove(i2,0,sumx,sumy,sumxy,sumx2,sumy2,sump)   ! 270312 131212
	else
c	sample position of target nucleon according to Woods-Saxon
c	 distribution
	call woodsax_samp(i2,0,alp,r0,am,ac,1)
	endif
	enddo
c       p+A    ! 060813
        elseif(ipden.eq.0 .and. itden.eq.1)then   !! 060813
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
c       p+p
c070417	elseif(ipden.eq.0 .and. itden.eq.0)then   !!
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
c       reaction plane eccentricity of participant nucleons
c131212	if(sigmsu.gt.0.)ecce=(sigmy2-sigmx2)/sigmsu
c131212
c       participant eccentricity of participant nucleons
        if(argu.gt.0. .and. sigmsu.gt.0.)
     c   ecce=sqrt(argu)/sigmsu !131212
c       calculate \epsilon{2}=\sqrt(<\epsilon_{part}^2>)
cc      ecce=ecce*ecce
c       note, \epsilon{2} should be \sqrt(aecceo), aecceo is a output
c        in paciae_21c.f
c       calculate transverse overlap area
        argu1=sigmx2*sigmy2-sigmxy*sigmxy
        if(argu1.gt.0.)secce=3.1416*sqrt(argu1) ! overlop area 250113
c131212
c       assuming ecce=geometric eccentricity of ellipsoid (\sqrt{(1-b^2/a^2)})
c        with half major axis b=pt*(1+smadel) and half minor axis
c        a=pt*(1-smadel), the resulted smadel=-ecce*ecce/4 (if neglecting
c        the samll term of ecce*ecce*(-2*smadel+smadel*smadel)
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
c       write(9,*)'ecce,smadel_a,smadel=',ecce,smadel_a,smadel
c250113
c       here a sign change is introduced because of asymmetry of initial
c        spatial space is oppsed to the final momentum space
c       write(9,*)'vneump,vneumt,sump,ecce,smadel=',
c     c  vneump,vneumt,sump,ecce,smadel
        endif   
c270312
c191110
c       for A+B or p+A or A+p or e+A 230311 240513 060813
	r0pt=r0p+r0t
c240513	if(itden.ne.0)then   ! 060813
	do i=1,nap
	c17(i,1)=c17(i,1)+bp
	enddo
c240513	endif   ! 230311
c191110
c	if(iii.eq.10)then
c	write(9,*)'after woodnat,bp=',bp
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
	pp1=dsqrt(ep1*ep1-pm2)
	pt1=-dsqrt(et1*et1-pm2)
	pm2=pmas(pycomp(2112),1)**2	
	pp2=dsqrt(ep2*ep2-pm2)
	pt2=-dsqrt(et2*et2-pm2)
	endif	
	if(ifram.eq.0)then
	pp1=win
	pt1=1.e-20
	pp2=win
	pt2=1.e-20
	pm2=pmas(pycomp(2212),1)**2	
	ep1=dsqrt(pp1*pp1+pm2)
	et1=dsqrt(pt1*pt1+pm2)
	pm2=pmas(pycomp(2112),1)**2
	ep2=dsqrt(pp2*pp2+pm2)
	et2=dsqrt(pt2*pt2+pm2)
	endif	
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
        nbe=0   ! 080104
        naf=0   ! 080104
        nsa=0
	idi=0
	idio=0
        do i1=1,kszj
        do j1=1,5
	k(i1,j1)=0
	p(i1,j1)=0.
	v(i1,j1)=0.
        kbe(i1,j1)=0
        pbe(i1,j1)=0.
        vbe(i1,j1)=0.
        kaf(i1,j1)=0
        paf(i1,j1)=0.
        vaf(i1,j1)=0.
        ksa(i1,j1)=0
        psa(i1,j1)=0.
        vsa(i1,j1)=0.
        enddo
        ndiq(i1)=0
        npt(i1)=0
        ifcom(i1)=0   ! 220110
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

c      '1 -> nzp' are projectile protons or e-, 'nzp+1 -> nap' are projectile 
c	neutrons; 'nap+1 -> nap+nzt' are targer protons, the rest are target
c	nuctrons. initiate 'pyjets'   ! 060813
	n=napt
	do i=1,n
	k(i,1)=1
	k(i,2)=2112
	p(i,5)=pmas(pycomp(2112),1)
c190517	if((i.le.nzp.and.ipden.ne.2).or.(i.gt.nap .and. i.le.nap+nzt))
c190517	     c	 then   ! 060813
        if(i.le.nzp .or. (i.gt.nap .and. i.le.nap+nzt))then
	k(i,2)=2212
	p(i,5)=pmas(pycomp(2212),1)
c060813
c190517	elseif(i.le.nap.and.ipden.eq.2)then
c	k(i,2)=11
c	p(i,5)=pmas(pycomp(11),1)
c	else
c060813
	endif
	do j=1,3
	p(i,j)=p17(i,j)
	v(i,j)=c17(i,j)
	enddo
	p(i,4)=p17(i,4)
	v(i,4)=tp(i)
	enddo
500	continue    ! 031103
c	v, vbh and vsa arraies are the position four vector
c	note: for v etc., we do not take care of their fifth component
c	 for array k, we take care of only first three components
c	write(9,*)'after initializing iii=',iii   !s
c	call psum(p,1,n,peo)   !!
c	write(9,*)'psum, after initializing nucleus-nucleus collision'   !!
c	write(9,*)peo   !!

c	boost PYJETS into cms of initial nucleus-nucleus collision system 
c	 from lab or initial nucleon-nucleon cms system.
c	call pyrobo(1,n,0.0,0.0,bst(1),bst(2),bst(3))
c	Lorentz contract
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
	gamp=1.d0/dsqrt(dmax1(1.d-20,(1.0d0-bzp*bzp)))
c190517	if(ipden.eq.2)gamp=1.   ! 060813
	gamt=1.d0/dsqrt(dmax1(1.d-20,(1.0d0-bzt*bzt)))
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
c060813
c       since e- was moved to last position after calling filt, one has to
c        remove it to the fist position
c190517	if(ipden.eq.2)then
c	call ltof(n)   ! move last particle (e-) in pyjets to first position
c190517	endif
c060813
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
c	note: particle list is composed of the arraies in common block 
c	 'sa2', the array 'ishp' in common block 'wz', the array 'tau' in 
c        common block 'sa4', and the array 'numb' in common block 'sa5'
	time=time_ini   ! 081010
	irecon=0
	time=time_ini   ! 081010
	call copl(time)
c       calculate the position for the center of mass of the
c	non-freeze-out system. The distance of a particle, when checking
c	is it freezing out or not, is measured with respect to this center
c       creat the initial collision list, note: be sure that the initial  
c	collision list must not be empty
	call ctlcre(lc,tc,tw)

c070417 move origin of time to collision time of first nucleon-nucleon collision
c	find out colli. pair with least colli. time
	call find(icp,tcp,lc,tc,tw,0)
	if(icp.eq.0)stop 'initial collision list is empty'   !
	time=tcp
c070417 perform classical Newton motion in Lab. system for all particles 
	call his(time,lc,tc,tw,istop)
	do ij=1,nsa
	vsa(ij,4)=0.
	enddo
	do ij=1,nctl
	tc(ij)=tc(ij)-time+1.e-5
	enddo
	time=time_ini   ! 081010
	call copl(time) 
400	continue

c       administrate a nucleus-nucleus collision 
        call scat(time,lc,tc,tw,win,parp21,parp22,psno,ijk,ipau,irecon,
     c   gamt)   ! 021207
        if(ijk.eq.1)return   
        time_ini=time   ! 081010
c	write(9,*)'af scat, iii,time,time_ini=',iii,time,time_ini   !s

800	continue
c	'sbe' to 'pyjets'
c	call tran_sbe
c       'saf' to 'pyjets'
	call tran_saf 
c241110
c       if(iii.eq.5)then
c       write(22,*)'af scat'
c       call pylist(1)
c       endif
c241110
        n00=n   ! 220110
c220110 n00: 'largest line number' in 'pyjets'
c220110 partons above n00 appear after inelastic collision  
c       'sa2' to 'sbh'
        nbh=0
        if(nsa.ge.1)then
        nbh=nsa
        do i1=1,nsa
        do i2=1,5
        kbh(i1,i2)=ksa(i1,i2)
        pbh(i1,i2)=psa(i1,i2)
        vbh(i1,i2)=vsa(i1,i2)
        enddo
        enddo
        endif
	do i1=nbh+1,kszj
        do i2=1,5
        kbh(i1,i2)=0
        pbh(i1,i2)=0.
        vbh(i1,i2)=0.
        enddo
        enddo
c 	P(N,5)=SQRT(MAX(-P(N,1)**2-P(N,2)**2-P(N,3)**2+P(N,4)**2,0.0))
c 	P(N-1,5)=SQRT(MAX(-P(N-1,1)**2-P(N-1,2)**2-P(N-1,3)**2
c     &	+P(N-1,4)**2,0.0))
c	call pyboro(1,n,0.0,0.0,-bst(1),-bst(2),-bst(3))
c	boost PYJETS back to lab or nucleon-nucleon cms system.

	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sysini(win)   ! 060813
c       give the initial values to quantities needed in calculation
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (KSZ1=30,kszj=40000)
	COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/count/isinel(600)
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
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
	param(1)=para1_1   ! 250204 200504
c	write(9,*)'para1_1=',param(1)   ! 250204 200504
c	rou0=PARAM(11)
c       considering the nucleus as a sphere with radii rnt for target
c        and rnp for projectile.
c        rnt=(3.*anat/(4.*3.1415926*rou0))**(0.33333)
c        rnp=(3.*anap/(4.*3.1415926*rou0))**(0.33333)
	rp00=1.12   ! 1.05 to 1.12 070613
	rt00=1.12   ! 1.05 to 1.12 070613
c070613	if(nap.gt.16)rp00=1.16*(1-1.16*anap**(-0.666666)) !rp00=1.122 (nat=208)
c070613	if(nat.gt.16)rt00=1.16*(1-1.16*anat**(-0.666666)) ! rt00=1.12 (nat=197)
	if(itden.eq.0)rnt=rt00*anat**(0.33333)   ! 310805
        if(itden.eq.1)rnt=rt00*anat**(0.33333)   ! +0.54  160511
        if(itden.eq.2)rnt=0.5
	if(ipden.eq.0)rnp=rp00*anap**(0.33333)   ! 310805
        if(ipden.eq.1)rnp=rp00*anap**(0.33333)   ! +0.54  160511
        if(ipden.eq.2)rnp=0.5
	if(nap.eq.2 .and. nzp.eq.1)rnp=4.0   ! 2.60 2.095  1.54
	rou0=3./4./3.1416*anat/(rnt*rnt*rnt)   ! 310805
        r0p=rnp
        r0t=rnt
C       set initial values to some quantities
c       in the program the x-sections are given in a unit of fm^2   ! 060813
        csnn=PARAM(1)*0.1
        cspin=PARAM(2)*0.1
        cskn=PARAM(3)*0.1
        cspipi=PARAM(4)*0.1
        cspsn=PARAM(13)*0.1
        cspsm=PARAM(14)*0.1
        csspn=PARAM(15)*0.1
        csspm=PARAM(16)*0.1
c060813
        if(ipden.eq.2)then
        if(ifram.eq.0)then
        ept=sqrt(win*win+0.938*0.938)
        rots=sqrt((ept+0.938)*(ept+0.938)-win*win)
        endif
        if(ifram.eq.1)rots=win
        call crosep(rots,csen)   
        csen=csen*0.1
	endif   
c060813
c       largest collision distance between two colliding particles.
        edipi=dsqrt(cspipi/3.1416)
        epin=dsqrt(cspin/3.1416)
        ekn=dsqrt(cskn/3.1416)
        ecsnn=dsqrt(csnn/3.1416)
	ecspsn=dsqrt(cspsn/3.1416)
	ecspsm=dsqrt(cspsm/3.1416)
	ecsspn=dsqrt(csspn/3.1416)
	ecsspm=dsqrt(csspm/3.1416)
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
        suppm=1.d0/(1+dexp(0.d0-r0p/alp))
        suptm=1.d0/(1+dexp(0.d0-r0t/alt))

        sig=PARAM(5)*0.1
        rcsit=PARAM(6)
	t0=PARAM(7)   ! 230805
c230805t0=0.   ! 221102   
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine arrove(ii,jj,sumx,sumy,sumxy,sumx2,sumy2,sump)   
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
        write(9,*)'subroutine arrove,infinitive loop may occur'
	goto 55   ! set larget try is equal to 10000
        endif
c       sample a point in the unit sphere 
        x=1.-2.*pyr(1)
        y=1.-2.*pyr(1)
        z=1.-2.*pyr(1)
        rr=x*x+y*y+z*z
        if(rr.gt.1) goto 54
	if(jj.eq.0)then   ! ii in target (origin)
	x=x*r0t
	y=y*r0t
	z=z*r0t
c       relative to projectile center, they are b-x, y, and z, respectively 
c	adjudge does (x-b,y,z) is in the sphere of projectile
        r1=sqrt((b-x)*(b-x)+y*y+z*z)
        if(r1.gt.r0p)goto 54
        c17(ii,1)=x
        c17(ii,2)=y
        c17(ii,3)=z
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
	x=x*r0p
	y=y*r0p
	z=z*r0p
c       relative to target center, they are x+b, y, and z, respectively 
c	adjudge does (x+b,y,z) is in the sphere of target
        r1=sqrt((x+b)*(x+b)+y*y+z*z)
        if(r1.gt.r0t)goto 54
	c17(ii,1)=x
        c17(ii,2)=y
        c17(ii,3)=z
c270312
        sumx=sumx+x
        sumy=sumy+y
        sumxy=sumxy+x*y   ! 131212
        sumx2=sumx2+x*x
        sumy2=sumy2+y*y
        sump=sump+1.
c270312
	endif
55	return
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
	goto 200   ! set larget try is equal to 100000
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
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
	PARAMETER (kszj=40000,KSZ1=30)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	cita=2*pyr(1)-1.
	fi=2.*pio*pyr(1)
	sita=dsqrt(1.-cita**2)
	c17(i,1)=xf*sita*dcos(fi)
	c17(i,2)=xf*sita*dsin(fi)
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
c	D,Dba,D0,D0ba,lamdac+,sigmac0,sigmac+,sigmac++,omega,k*+,K*0,
c	D*,D*ba,D*0,D*0ba (52 kind of particle altogether)
c060813 in case of e+A, one images e as a initial projectile proton
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
      PARAMETER (kszj=40000,KSZ1=30)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
	common/sa6/kfmaxi,nwhole
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
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
c	order particles according to flavor code
c	j: the particle needed to order
c	ipi: j-th particle should order after ipi
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
	parameter(kszj=40000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	dimension pp(5),vv(5),kk(5)
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
	subroutine prt_pyj(nn,cc)
c	print particle list and sum of momentum and energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=40000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/pyjets/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
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
        cc=ich1/3.
        write(mstu(11),*)'c & p sum=',cc,peo   ! 
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbh(nn,cc)
c       print particle list and sum of momentum and energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=40000)
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
        cc=ich1/3.
        write(mstu(11),*)'c & p sum=',cc,peo   !
c	write(9,*)peo,ich1/3   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa2(nn,cc)
c       print particle list and sum of momentum and energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=40000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
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
        cc=ich1/3.
        write(22,*)'c & p sum=',cc,peo   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbe(nn,cc)   ! 220110
c       print particle list and sum of momentum and energy
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sbe/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
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
        cc=ich1/3.
        write(22,*)'c & p sum=',cc,peo   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_saf(nn,cc)   ! 220110
c       print particle list and sum of momentum and energy
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/saf/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
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
        cc=ich1/3.
        write(22,*)'c & p sum=',cc,peo   !
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine psum(pei,il,ih,peo)
c       calculate sum of momentum and energy
c       pei: two dimension array of input momentum and energy
c       il and ih: lower and upper limits of sum
c       peo : one dimension array of output momentum and energy  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=40000)
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
        subroutine scat(time,lc,tc,tw,win,parp21,parp22,psno,ijk,
     c	 ipau,irecon,gamt)   ! 021207
c	administrate a nucleus-nucleus collision and e+A collision ! 060813
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
      PARAMETER (kszj=40000,KSZ1=30)
        parameter(nsize=240000)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23 
	common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
 	common/sa16/dtt,dni(10),dpi(10),edi(10),bmin,bmax
     &   ,bar(10),abar(10),barf(10),abarf(10)   
     &   ,emin(10),eminf(10),eplu(10),epluf(10)   
        common/sa18/tdh,itnum,non18,cptl,cptu,cptl2,cptu2,snum(4,20),
     &	 v1(4,20),v2(4,20),v12(4,20),v22(4,20)
	common/sa23/kpar,knn,kpp,knp,kep   ! 060813   
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 220110
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21,adiv,gpmax,nnc   !   070417
        common/sa28/nstr,nstr00,nstra(kszj),nstrv(kszj)   ! 220110
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
	common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5) 
       	common/ctllist/nctl,noinel(600),nctl0,noel
	common/sa34/iikk   ! 060617
        dimension lc(nsize,5),tc(nsize),tw(nsize)
	dimension pi(4),pj(4),pii(4),pjj(4),peo(4),pint(4)	
	dimension nni(10),ndi(10),npi(10)
	dimension pkk(kszj,4)   
	dimension cc(5),b(3),bkk(3)    	
        dimension skapa(5),ksinn(kszj,5),psinn(kszj,5),vsinn(kszj,5)   !070417 110517
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c	arraies in 'pyjets' are used in the processes in calling 'pythia', 
c	 after nn scattering, and after nucleus-nucleus and/or e+A collision
c        060813
c	arraies in 'sa2' are used in the processes in nn collision
c       arraies in 'sbh' are used to store hadron after nn collision
c        and after nucleus-nucleus collision as well as e+A collision 060813
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c	numb(i) is used in the scattering processes, numbs(i) is used in 
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
c	common block 'sbe' stores cumulatively parton (q,qq, and g) 
c	 configuration before breaking the diquarks 
c	common block 'saf' stores cumulatively parton (q and g) 
c	 configuration after breaking the diquarks
c	idi: counts cumunatively the number of diquark (anti-diquark)
c	idio: value of idi after last nn collision
c       ndiq(j): = 0 if j is quark (antiquark) 
c                = idi if j is diquark (anti-diquark) 
c       note: j is line number in 'sbe' ('saf')
c       npt(idi): = line number of idi-th diquark (anti-diquark) 
c	 partner in 'saf'
c220110 ifcom(idi): line number of first component of idi-th diquark
c       nstr: statitics of number of strings in a hh collis.
c       nstr00: number of strings after call remo
c       nstra(i): line number of first component of i-th string
c220110 nstrv(i): line number of last component of i-th string
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

        nctl0=nctl
c070417	mstj(1)=mstj1_1   ! 221203
c070417 loop over hadron-hadron collisions in a nucleus-nucleus collision
        if(kjp22.eq.0 .or. kjp22.eq.1)then
c       nnc: statistics of number of nucleon-nucleon collisions in a nucleus-nucleus
c        collision, statistics of variables of itime etc. over nucleon-nucleon
c        collisions in a nucleus-nucleus collision
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
        iii=1   ! 220110, iii-th hadron-hadron collis.
10 	if(iii.eq.1)goto 1000
101	call copl(time)
c	find out the binary colli. with minimum collsion time
1000	call find(icp,tcp,lc,tc,tw,1)
        time=tcp   !   190517
	if(icp.eq.0)goto 100
c	icp=0 means the collision list is empty
	l=lc(icp,1)
	l1=lc(icp,2)
cm	write(9,*)'af find, iii,icp,l,l1,tcp=',iii,icp,l,l1,tcp   ! sa
	time0=time
	kfa=ksa(l,2)
	kfb=ksa(l1,2)
        ikfa=iabs(kfa)   ! 070417
        ikfb=iabs(kfb)   ! 070417
	time=tcp
c	record this collision time


c?????????????????????????????????????????????????????????????????
csa	if(time.le.ttt)then
c	record the spatial and momentum coordinates
csa	if(iiii.eq.1)call txp(time0,time)
c	calculate the time dependent direct and elliptic flow
csa	call flow_t(time0,time)
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
        pti=dsqrt(pi(1)**2+pi(2)**2)
        ptj=dsqrt(pj(1)**2+pj(2)**2)
c200601
c	boost to CMS frame of colliding pair
	call lorntz(ilo,b,pi,pj)
	ss=pi(4)+pj(4)
        if(ss.lt.1.e-18)ss=1.e-18
c	perform classical Newton motion in Lab. system  
	call his(time,lc,tc,tw,istop)
	if(istop.eq.1)goto 100
c	istop=1 means all particles have get out of considered volume
	m1=numb(1)
        m2=numb(2)
	m3=numb(3)
        m4=numb(4)
        m7=numb(7)
c	write(9,*)'m1-m4=',m1,m2,m3,m4
c060805	if((l.le.m4 .and. l1.le.m4) .and. ss.ge.parp21)then   ! if 1
c241110	if(((l.le.m2 .and. l1.le.m2).or.(kfa.eq.2212.and.kfb.eq.-2212)
c241110     c   .or.(kfb.eq.2212.and.kfa.eq.-2212)) .and. ss.ge.parp21)
c241110     c   then   ! if 1
        if((l.le.m2 .and. l1.le.m2) .and. ss.ge.parp21)then  ! if 1 011210 060813 070417
c060813	m7 to m2
c	calculate the angular 'seta' of the momenta pi and pj
	ctai=pyangl(pi(3),dsqrt(pi(1)**2+pi(2)**2))
	ctaj=pyangl(pj(3),dsqrt(pj(1)**2+pj(2)**2))
	cctai=dcos(ctai)
	cctaj=dcos(ctaj)
	if(cctai.gt.0.)then
c       calculate the 'orentation' of the vector pi
	call codi(pi,cfi1,sfi1,ccta1,scta1)
	else
	call codi(pj,cfi1,sfi1,ccta1,scta1)
	endif

	if(kfa.eq.2212.and.kfb.eq.2212)then   ! 070417
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

	if(kfb.eq.2212.and.kfa.eq.2112)then   ! 070417
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

	if(kfa.eq.2112.and.kfb.eq.2112)then   ! 070417
c200601
        knn=knn+1
        if(pti.le.1.e-4)kpar=kpar+1
        if(ptj.le.1.e-4)kpar=kpar+1
c200601
	call pyinit('cms','n0','n0',ss)   
        call pyevnt
	endif   ! 070417
c011210

2222	call pyedit(2)   

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
c	give four position to the particles after calling pyevnt
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
c110517
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
        if(iikk.eq.0)then   ! 060617 
c	pyjets to 'sinn' which is a internal array
	nsin=n
	do i1=1,n
        do i2=1,5
        ksinn(i1,i2)=k(i1,i2)
        psinn(i1,i2)=p(i1,i2)
        vsinn(i1,i2)=v(i1,i2)
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
c	'sinn' to sbh
	nbh=nsin
        do i1=1,nbh
        do i2=1,5
        kbh(i1,i2)=ksinn(i1,i2)
        pbh(i1,i2)=psinn(i1,i2)
        vbh(i1,i2)=vsinn(i1,i2)
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
c070802
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
	if(iikk.eq.0)then   ! 060617
	igq=0
	do j1=1,n
	kfj1=iabs(k(j1,2))
c140805	if(kfj1.le.21)igq=igq+1
        if(kfj1.le.8.or.kfj1.eq.2101.or.kfj1.eq.3101.or.kfj1.eq.3201
     c   .or.kfj1.eq.1103.or.kfj1.eq.2103.or.kfj1.eq.2203.or.kfj1.eq.
     c   3103.or.kfj1.eq.3203.or.kfj1.eq.3303.or.kfj1.eq.21)igq=igq+1! 140805
	enddo
	if(igq.eq.0)then   ! no q, diquark, and g at all
c130206	
c	'pyjets' to 'sbh'
c241110	nbh=n
c	do j1=1,n
c	do j2=1,5
c	kbh(j1,j2)=k(j1,j2)
c	pbh(j1,j2)=p(j1,j2)
c	vbh(j1,j2)=v(j1,j2)
c	enddo
c	enddo
c	do j1=n+1,kszj
c	do j2=1,5
c	kbh(j1,j2)=0
c	pbh(j1,j2)=0.
c	vbh(j1,j2)=0.
c	enddo
c	enddo
c241110	goto 200
c130206
c	remove current nn collision pair from collision list
	do j1=icp+1,nctl   ! active on 241110
	j=j1-1   ! 241110
	tc(j)=tc(j1)   ! 241110
        tw(j)=tw(j1)   ! 241110
        do m=1,5   ! 241110
        lc(j,m)=lc(j1,m)   ! 241110
        enddo   ! 241110   
        enddo   ! 241110
        nctl=nctl-1   ! 241110
        iii=iii+1   ! 060805 241110
        goto 10   ! 241110
	endif
c	reconstruct leading nucleon
	if(nap.ne.1.and.nat.ne.1)then
	irecon=irecon+1
c       call recons(irecon)   
c	write(22,*)'af. recons n,irecon,iii,event=',n,irecon,iii,iiii
c	call pylist(1)
	endif
c241110
c       if(iiii.eq.5)then
c       write(22,*)'be remo iii=',iii
c       call pyedit(2)
c       call pylist(1)
c       endif
c241110
c       remove hadrons from 'pyjets' to 'sbh' and truncate 'pyjets'
c	 correspondingly
c110517	call remo
	call pyedit(2)   ! 220517
c241110
c       if(iiii.eq.5 .and. (iii.ge.12.and.iii.lt.15))then
c       write(22,*)'af remo'
c       call pylist(1)
c       call prt_sbh(nbh,cc)
c       endif
c241110
c220110
c080104
c	'pyjets' to 'sbe'. etc.
	if(n.ge.1)then   ! 1
	do i1=1,n
	i3=i1+nbe
	kf=k(i1,2)
        kfab=iabs(kf)
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then   ! 2
c     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then   ! 2
        idi=idi+1
        ndiq(i1+naf)=idi
	endif   ! 2
        do i2=1,5
        kbe(i3,i2)=k(i1,i2)
        pbe(i3,i2)=p(i1,i2)
        vbe(i3,i2)=v(i1,i2)
        enddo
	enddo
	nbeo=nbe   ! 190204
	nbe=i3
	endif   ! 1
c	write(9,*)'af. fill sbe n,nbe,irecon,iii,event=',n,nbe,irecon,
c     c   iii,iiii   ! sa
c	goto 200
c080104
c       break up diquark and give four momentum and four position
c        to the broken quarks (working in 'pyjets')
	call break
c241110
c        if(iiii.eq.5 .and. (iii.ge.12.and.iii.lt.15))then
c        write(22,*)'af break'
c        call pylist(1)
c        call prt_sbh(nbh,cc)
c        endif
c241110
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
        do i1=nbeo+1,nbe
        do i2=1,5
        kbe(i1,i2)=0
        pbe(i1,i2)=0.
        vbe(i1,i2)=0.
        enddo
	ndiq(i1+naf)=0
        enddo
	do i1=idio,idi
	npt(i1)=0
        ifcom(i1)=0   ! 220110
	enddo
	nbe=nbeo
c190204
	idi=idio    
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
c       'pyjets' to 'saf'. etc.
        if(n.ge.1)then   
        do i1=1,n
        naf=naf+1
	if(naf.gt.kszj)then
	write(9,*)'iiii,naf,kszj=',iiii,naf,kszj   ! sa
	stop 11111
	endif
        do i2=1,5
        kaf(naf,i2)=k(i1,i2)
        paf(naf,i2)=p(i1,i2)
        vaf(naf,i2)=v(i1,i2)
        enddo
        enddo
        endif
c080104
200	continue
	idio=idi   ! 080104
	endif   ! 060617
c060617
	if(iikk.eq.1)then 
c       'pyjets' to 'sbh'
	nbh=n
	do j1=1,n
	do j2=1,5
	kbh(j1,j2)=k(j1,j2)
	pbh(j1,j2)=p(j1,j2)
	vbh(j1,j2)=v(j1,j2)
	enddo
	enddo
	do j1=n+1,kszj
	do j2=1,5
	kbh(j1,j2)=0
	pbh(j1,j2)=0.
	vbh(j1,j2)=0.
	enddo
	enddo
	endif
c060617
c080104 241110
c       if(iiii.eq.5 .and. iii.eq.13)then
c       write(22,*)'be updpip'
c       write(9,*)'be updpip iiii,iii=',iiii,iii
c       call prt_sa2(nsa,cc)
c       call prt_sbh(nbh,cc)
c       endif
c080104 241110
c	update particle list after  calling 'pythia' (i. e. 'sbh' to 
c	 'sa2' and truncate 'sa2')
cm	write(9,*)'be updpip l,l1,kfa,kfb,nsa,nbh,kbh(1-nbh,2)=',
cm     c	 l,l1,kfa,kfb,nsa,nbh,(kbh(i1,2),i1=1,nbh)   ! sa
c	call prt_sbh(nbh,cc)
	call updpip(l,l1,icp,lc,tc,tw,time,iii)   
c241110
c       if(iiii.eq.5 .and. (iii.ge.12.and.iii.lt.15))then
c       write(22,*)'af updpip'
c       write(9,*)'af updpip iiii,iii=',iiii,iii
c       call prt_sa2(nsa,cc)
c       call prt_saf(naf,cc)
c       endif
c241110
c011204	l=lc(icp,1)
c011204	l1=lc(icp,2)
c       update collision list after calling 'pythia'
        call updtlp(time,lc,tc,tw,iii)
cm	write(9,*)'af updtlp iii,nsa,nctl=',iii,nsa,nctl   ! sa
cm	call prt_sa2(nsa,cc)   
cm	do i=1,nctl
cm	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
cm	enddo
	if(nctl.eq.0)goto 100   ! 021204
	goto 300   ! ss is enough to call pythia	
	endif   ! if 1

c       if ss is not enough to call pythia then treatting as elastic 
c	 collision
c	write(9,*)'elastic ss=',ss
c       calculate four-momentum of two particles after elastic reaction, pi
c       and pj in CMS frame
9004	call coelas(l,l1,ss,pi,pj)   ! 220517
c       update the particle list for elastic scattering, pi and pj have been
c       boosted back to Lab fram 
	call updple(l,l1,b,pi,pj)
c	statistics of the number of elastic nn collisions
        noinel(1)=noinel(1)+1
c       update the collision list after elastic scattering
        call updatl(l,l1,time,lc,tc,tw,iii)

c	write(9,*)'ela. iiii,iii,icp=',iiii,iii,icp   ! sa
c	call prt_sa2(nsa,cc)   
c	do i=1,nctl
c	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
c	enddo
   
300	continue
c241110
c        if(iiii.eq.5)then
c        nzpt=nzp+nzt
c        sumch=0.
c        write(22,*)'iiii,iii,nsa,naf,ss,l,l1,kfa,kfa=',
c     c   iiii,iii,nsa,naf,ss,l,l1,kfa,kfb
c        write(22,*)'numb=',(numb(i1),i1=1,52)
c        call prt_sbh(nbh,cc)
c        call prt_sa2(nsa,charge)
c        sumch=sumch+charge
c        write(22,*)'sa2, charge=',charge
c        call prt_saf(naf,charge)
c        sumch=sumch+charge
c        write(22,*)'saf, charge=',charge
c        if(sumch.ne.nzpt)write(22,*)'charge not cons.,iii,sumch=',
c     c   iii,sumch
c        endif
c241110
        iii=iii+1
c	if(iii.eq.2)stop   ! temporal
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
	return
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine his(t1,lc,tc,tw,istop)
c	classical Newton motion in Lab. system
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000)
	parameter(nsize=240000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa4/tau(kszj),tlco(kszj,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/ctllist/nctl,noinel(600),nctl0,noel
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	dimension lc(nsize,5),tc(nsize),tw(nsize)
	istop=1
	in=0
	do 200 i=1,nsa
	r0=rao*dmax1(rnt,rnp)   ! 060813
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
	pp4=psa(i,4)
c	due to the fast speed of bayons, we could not use a limited interaction
c	region
c060813	r0=rao*dmax1(rnt,rnp)
c	if(iabs(k(i,2)).gt.1000)r0=1.E+10*r0
	do j=1,3
	vp=psa(i,j)/pp4
	vsa(i,j)=vsa(i,j)+vp*(t1-vsa(i,4))
	aa=aa+(vsa(i,j)-coor(j))**2
	enddo
c251004	vsa(i,4)=t1
	aa=dsqrt(aa)
	if(aa.lt.r0) goto 100
c	if freeze-out occurs deduct the distance between the last collision 
c	and now
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
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
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
        sita=dsqrt(1.-cita**2)
	v(i,1)=sita*dcos(fi)
	v(i,2)=sita*dsin(fi)
	v(i,3)=cita
	endif
	v(i,4)=time
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
	subroutine ptcre_n(l,l1,time,gamt)   ! 021207
c	arrange particles (quark,diquark, and gluon mainly) after 
c	 calling pythia into the overlap region randomly  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	PARAMETER (kszj=40000,KSZ1=30)
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	if(ipden.ne.0 .or. itden.ne.0)then
        b=bp/r0t
	do i=1,n
        iii=0
54      iii=iii+1
        if(iii.eq.10000)then
        write(22,*)'difficult to arrange produced particles in'
        write(22,*)'subroutine ptcre,infinitive loop may occur'
        endif
c       sample a point in the unit sphere of target
        x=1.-2.*pyr(1)
        y=1.-2.*pyr(1)
        z=1.-2.*pyr(1)
        rr=dsqrt(x*x+y*y+z*z)
        if(rr.gt.1) goto 54
c       x and y components of that point in the system of unit sphere of
c        projectile are x and y-b, respectively. Adjudge that does (x,y-b) is 
c        in the sphere of projectile
        r1=r0p*dsqrt(x*x+(b-y)*(b-y))
        if(r1.gt.r0p)goto 54
        xx=x*r0t
        yy=y*r0t
        zz=z*r0t/gamt
        v(i,1)=xx
        v(i,2)=yy
        v(i,3)=zz  
c       write(5,*)xx,yy,zz       
        enddo
	endif
	if(ipden.eq.0 .and. itden.eq.0)then
        do i=1,n
	cita=2*pyr(1)-1.
        fi=2.*pio*pyr(1)
        sita=dsqrt(1.-cita**2)
	v(i,1)=sita*dcos(fi)
	v(i,2)=sita*dsin(fi)
	v(i,3)=cita
        v(i,3)=v(i,3)/gamt
        enddo
	endif
        do i=1,n
	v(i,4)=time   ! 230805
	enddo
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updpip(l,l1,icp,lc,tc,tw,time,iii)   
c       update particle list 'sa2' after calling pythia  
c	 (i. e. 'sbh' to 'sa2' and truncate 'sa2')   
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter(kszj=40000)
        parameter(nsize=240000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        COMMON/SBH/N,NONBH,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 080104
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/ctllist/nctl,noinel(600),nctl0,noel
	common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio   ! 060813
c060813	ipyth: stord the order number in the particle list (sa2) of 'pythia' 
c060813	 particles
        dimension lc(nsize,5),tc(nsize),tw(nsize)  
	dimension peo(4)
        do m=1,2000
        ipyth(m)=0
        enddo
c       put 'sbh' to 'sa2'   
c241110
c        if(iiii.eq.5)then
c        write(22,*)'in updpip iiii,iii,nsa,l,l1,nbh=',
c     c   iiii,iii,nsa,l,l1,n
c       write(22,*)'kfaco=',(kfaco(i1),i1=1,52)
c        write(22,*)'numb=',(numb(i1),i1=1,52)
c        call prt_sa2(nsa,cc)
c        call prt_sbh(n,cc)
c        endif
c241110
        ll=l
        ll1=l1
        if(n.eq.0)goto 200   ! 241110
        do 500 i=1,n
	kf=k(i,2)
        do 600 j=1,kfmax
        if(kf.ne.kfaco(j))goto 600
        jj=numb(j)+1
c        if(iiii.eq.5 .and. iii.eq.13)write(22,*)'j,jj=',j,jj
c       update particle list etc.
	do m=nsa,jj,-1
        mm=m+1
c080104	ksa(mm,2)=ksa(m,2)
c080104	ksa(mm,1)=1
c080104	ksa(mm,3)=ksa(m,3)
        do m1=1,5
        ksa(mm,m1)=ksa(m,m1)   ! 080104
        psa(mm,m1)=psa(m,m1)
	vsa(mm,m1)=vsa(m,m1)
	enddo
        ishp(mm)=ishp(m)
        tau(mm)=tau(m)
	enddo
        do m=1,2000
        ipym=ipyth(m)
        if(ipym.ge.jj)ipyth(m)=ipym+1
        enddo
        if(ll.ge.jj)ll=ll+1
        if(ll1.ge.jj)ll1=ll1+1
c       update the values of lc(m,1-2) with value.ge.jj
        do m=1,nctl
        lc1=lc(m,1)
        if(lc1.ge.jj)lc(m,1)=lc1+1
        lc2=lc(m,2)
        if(lc2.ge.jj)lc(m,2)=lc2+1
        enddo
c       give proper values to particle jj.
c221203	ksa(jj,2)=kf
c221203	ksa(jj,1)=1
c221203	ksa(jj,3)=0
        do m=1,5
        ksa(jj,m)=k(i,m)   ! 221203
        psa(jj,m)=p(i,m)
	vsa(jj,m)=v(i,m)
	enddo
        ishp(jj)=1
        tau(jj)=time+t0*p(i,4)/p(i,5)
c	the values of 'ishp' and 'tau' for hadrons from 'pythia' 
c	 are given here, the proper formation time of 'pythia' particle 
c	 is assume to be equal to t0 fm/c, except nucleon and j/psi
	if(kf.eq.2212 .or. kf.eq.2112)then
	tau(jj)=time+t0*p(i,4)/p(i,5)*taup
        elseif(kf.eq.443.or.kf.eq.30443)then
        tau(jj)=time+t0*p(i,4)/p(i,5)*taujp
        endif
	ipyth(i)=jj
	do m=j,kfmax
        numb(m)=numb(m)+1
        enddo
	nsa=nsa+1
	goto 500
600     continue
c	if produced hadron is not in considered hadron list
	nsa=nsa+1
        do m=1,5
        ksa(nsa,m)=k(i,m)
        psa(nsa,m)=p(i,m)
        vsa(nsa,m)=v(i,m)
        enddo
        ishp(nsa)=0
        tau(nsa)=0.
	ipyth(i)=nsa   
500	continue
200     continue   ! 241110
c241110
c        if(iiii.eq.5)then
c        write(22,*)'updpip af remove iiii,iii,l,l1,nbh=',iiii,iii,l,l1,n
c        write(22,*)'numb=',(numb(i1),i1=1,52)
c        call prt_sa2(nsa,cc)
c        call prt_sbh(n,cc)
c        endif
c241110
        l=ll
        l1=ll1
c       remove colli. pair composed of l or l1  
	jj=0
	do 300 ii=1,nctl
	i1=lc(ii,1)
        j1=lc(ii,2)
        if(i1.eq.l .or. i1.eq.l1)goto 300
        if(j1.eq.l .or. j1.eq.l1)goto 300
	jj=jj+1
        tc(jj)=tc(ii)
        tw(jj)=tw(ii)
        do m=1,5
        lc(jj,m)=lc(ii,m)
        enddo	
300	continue
	do ii=jj+1,nctl+1
        tc(ii)=0.0
        tw(ii)=0.0
        do m=1,5
        lc(ii,m)=0
        enddo
        enddo
	nctl=jj
cm	write(9,*)'updpip af remove collis. pairs nctl=',nctl   ! sa
cm	do i=1,nctl
cm	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
cm	enddo
c	remove hadrons l and l1 from 'sa2'
	kf1=ksa(l,2)   
        kf2=ksa(l1,2)
c	write(9,*)'l,l1,kf1,kf2=',l,l1,kf1,kf2   !!
	kf=kf1
        ll=l
        do 700 i=1,2
        if(ll.eq.nsa)then   ! 
	do i1=1,kfmax
	if(kf.ne.kfaco(i1))goto 400
c241110	numbm=numb(i1)
c241110	do i2=1,i1
c241110	if(numb(i2).eq.numbm)numb(i2)=numb(i2)-1
c241110	enddo
        do m=i1,kfmax
        numb(m)=numb(m)-1
        enddo
        if(i1.gt.1)then
        numba=numb(i1)
        do m=1,i1-1
        if(numb(m).eq.numba)numb(m)=numb(m)-1
        enddo
        endif
c241110
	goto 100
400	enddo
        endif   !
        do j=ll+1,nsa
        jj=j-1
c080504	ksa(jj,2)=ksa(j,2)
c080504	ksa(jj,1)=1
c080504	ksa(jj,3)=ksa(j,3)
        do m=1,5
	ksa(jj,m)=ksa(j,m)   ! 080504
        psa(jj,m)=psa(j,m)
	vsa(jj,m)=vsa(j,m)
	enddo
        ishp(jj)=ishp(j)
	tau(jj)=tau(j)
	enddo
	if(nctl.eq.0)goto 900
        do m=1,nctl
        lc1=lc(m,1)
        lc2=lc(m,2)
        if(lc1.gt.ll)lc(m,1)=lc1-1
        if(lc2.gt.ll)lc(m,2)=lc2-1
        enddo
900	do 800 j=1,kfmax
        if(kf.ne.kfaco(j))goto 800
        do m=j,kfmax
        numb(m)=numb(m)-1
        enddo
	if(j.gt.1)then
	numba=numb(j)
	do m=1,j-1
	if(numb(m).eq.numba)numb(m)=numb(m)-1
        enddo
	endif
        goto 100
800     continue
100     continue
        nsa=nsa-1
        if(l1.gt.ll)l1=l1-1
        do m=1,2000
        ipym=ipyth(m)
        if(ipym.gt.ll)ipyth(m)=ipym-1
        enddo
        if(i.eq.2)goto 700
        ll=l1
        kf=kf2
700     continue
	return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coelas(ic,jc,eij,pi,pj)
c	perform elastic scattering
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000)
	COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c       note the name of the arraies in 'sa2' in this subroutine
	dimension pi(4),pj(4)
	iic=k(ic,2)
	jjc=k(jc,2)
	d=3.65*(eij-pmas(pycomp(iic),1)-pmas(pycomp(jjc),1))
	if(d.lt.1.e-10)return
	pt=0.2
	a=dmin1(10.3d0,1.d0/(1.12d0*pt)/(1.12d0*pt))
	d6=d**6
	b=d6*a/(1.+d6)
	if(b.lt.1.e-20)then
	b=1.e-20
	endif
	pm2=pi(1)**2+pi(2)**2+pi(3)**2
	pm=dsqrt(pm2)
	t0=-4.*pm2
	if(dabs(t0).lt.1.d-20)then
	cctas=1.
	goto 100
	endif
	cc=pyr(1)
	if(dabs(b*t0).lt.0.0001d0)then
	abt=1.
c	elseif(b*t0.lt.-50.)then
c	abt=0.
	else
	abt=dexp(dmax1(-7.0D2,dble(b*t0)))
	endif
	tt1=dlog(cc+(1.-cc)*abt)
	if(dabs(tt1).lt.1.d-30 .and. b.le.1.d-20)then
	cctas=1.
	goto 100
	endif
	tt=tt1/b
	if(dabs(tt).lt.1.d-20)then
	cctas=1.
	goto 100
	endif
	cctas=1.-tt*2./t0
	if(dabs(cctas).gt.1.d0)then
	cctas=dsign(1.d0,cctas)
	endif
100	continue
	sctas=dsqrt(1.-cctas**2)
	fis=2.*3.1416*pyr(1)
	cfis=dcos(fis)
	sfis=dsin(fis)
	call rotate(cctas,sctas,cfis,sfis,pm,pi,pj)
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rotate(cctas,sctas,cfis,sfis,pp3,pi,pj)
c	perform rotation
c       pi,pj: input, four momentum of colliding pair before scattering
c              output,four momentum of scattered particles after rotation
c       pp3: momentum modulus of pi or pj, both are equal in their cms,
c        after scattering
c       cctas,sctas,cfis,sfis: direction cosines of momentum of one of 
c        scattered particle relative to the momentum
c        of corresponding particle before scattering
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	dimension pi(4),pj(4)
c	fi1=atan2(pi(2),pi(1))
c	cta1=atan2(dsqrt(pi(1)**2+pi(2)**2),pi(3))
	fi1=pyangl(pi(1),pi(2))
	cta1=pyangl(pi(3),dsqrt(pi(1)**2+pi(2)**2))
	cfi1=dcos(fi1)
	sfi1=dsin(fi1)
	ccta1=dcos(cta1)
	scta1=dsin(cta1)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updple(ic,jc,b,pi,pj)
c	update particle list for elastic scattering 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000)
      COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c       note the name of the arrays in 'sa2'
	dimension pi(4),pj(4),b(3)
	ilo=1
c       ilo=1 for inverse Lorentz transformation
	call lorntz(ilo,b,pi,pj)
	do i=1,4
	p(ic,i)=pi(i)
	p(jc,i)=pj(i)
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
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (kszj=40000,KSZ1=30)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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
c            conserve
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
	if(dabs(1.-fr) .le. dep)goto 200
	do i=1,np
	ppm=pp(i,4)/0.938
	ppf=ppm/fr
	ff(i)=dsqrt(dabs(ppf*ppf-1.)/(ppm*ppm-1.))
	do j=1,3
	ppp=ff(i)*pp(i,j)
	pp(i,j)=ppp
	pxyz(j)=pxyz(j)+ppp
	enddo
	enddo
	do i=1,3
	arp(i)=dabs(1.-pxyz(i)/ps(i))
	pxyz(i)=pxyz(i)-ps(i)
	enddo
	if(dabs(1.-fr).le.dep .and.arp(1).le.dep .and. arp(2).le.dep  
     c   .and. arp(3).le.dep) goto 200
	do i=1,3
	pxyz(i)=pxyz(i)/np
	enddo
	do i=1,np
	do j=1,3
	pp(i,j)=pp(i,j)-pxyz(j)
	enddo
	pp(i,4)=dsqrt(0.880+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
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
c        cta1=atan2(dsqrt(pi(1)**2+pi(2)**2),pi(3))
	fi1s=pyangl(pis(1),pis(2))
        cta1s=pyangl(pis(3),dsqrt(pis(1)**2+pis(2)**2))
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
	pp=dsqrt(pp)
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
	subroutine flow_t(tt,tt1)
c	calculate direct and elliptic flow within time interval of 
c        [0,itnum*tdh] in partonic initialization stage 
c       tdh,itnum: the time step and the number of time steps (used
c        in subroutine 'flow_t')
c       cpt,cptu;cptl2,cptu2 : pt cut for g;u and ubar 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000,KSZ1=30)
        common/pyjets/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa18/tdh,itnum,non18,cptl,cptu,cptl2,cptu2,snum(4,20),
     &	 v1(4,20),v2(4,20),v12(4,20),v22(4,20)
	do i1=1,itnum
	ti=(i1-1)*tdh
c       calculate directed and elliptic flow at time of ti
	if(tt.le.ti .and. tt1.gt. ti)then
c       since the particles after collision at tt have contribution
c        to the statistics at any moment in between this collision (
c        tt) and next collision (tt1)

        do j=1,nsa
        ik=ksa(j,2)
c	yy=pyp(j,17)
c        if(ik.eq.2212 .and. (yy.gt.0.9 .and. yy.le.5.))then
	if(ik.eq.21)then
	px=psa(j,1)
	py=psa(j,2)
	px2=px*px
	py2=py*py
	pt2=px2+py2
        if(pt2.lt.1.e-20)pt2=1.e-20   ! 120607
	pt=dsqrt(pt2)
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
	if(iabs(ik).eq.2)then
	px=psa(j,1)
	py=psa(j,2)
	px2=px*px
	py2=py*py
	pt2=px2+py2
        if(pt2.lt.1.e-20)pt2=1.e-20   ! 120607
	pt=dsqrt(pt2)
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine recons(irecon)
c	reconstruct diquark-quark 'A and V' pair into nucleon 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=40000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
	dimension ps(4),rs(4),idi2(2,2),dele(kszj),pp(20,5),isuc(2)   ! 230407
	dimension iiglu(100),pk(5),vk(5),rr(3),pppp(50,2),kk(5)
c	idi2(i,1): line number of 'A' quark (or diquark) in i-th quark-
c	 diquark 'A and V' pair
c       idi2(i,2): line number of 'V' quark (or diquark) in i-th quark-
c        diquark 'A and V' pair
c	isuc(i): =1 if i-th quark-diquark 'A and V' pair can compose into
c	 nucleon, Delta(0), and Delta(+), otherwise =0
c	iiglu(i): line number of i-th gluon in 'pyjets'
        sigm2=adj1(29)   ! 0.26
        ptmax=adj1(30)   ! 2.
	do j1=1,2
	isuc(j1)=0
	do j2=1,2
        idi2(j1,j2)=0
	enddo
	enddo
	delte=0.
	do j1=1,kszj   ! 230407, 2 originally
	dele(j1)=0.
	enddo
	do j1=1,100
	iiglu(j1)=0
	enddo
        imc=adj1(13)
        ibc=adj1(14)
	adj23=adj1(23)   ! 180405
c050505	goto 600   ! 050505
c170205
c       the probability of gluon spliting into u,d & s 
	adj132=adj1(32)  
	prosum=1.+1.+adj132
	prod=1./prosum   ! 0.4286 originally
	pros=adj132/prosum   ! 0.1428 originally
	prods=prod+pros   ! 0.5714 originally
c170205
c	write(9,*)'irecon,event=',irecon,iiii
c	call pyedit(2)
c	write(22,*)'inter recons irecon,event=',irecon, iiii
c	call pylist(1)
c	count number of gluons
	jjj=0
	do j=1,n
	ik=k(j,2)
	if(ik.eq.21)then
	jjj=jjj+1
	iiglu(jjj)=j   ! line number of jjj-th gluon
	endif
	enddo
	jglu=jjj   ! number of gluons
c	write(9,*)'jglu,iiglu=',jglu,(iiglu(j),j=1,jglu)
	if(jglu.eq.0)goto 600
        if(jglu.eq.1)then   ! if 1
c	force breaking that gluon
	ii1=iiglu(1)
	ps(1)=p(ii1,1)
	ps(2)=p(ii1,2)
	ps(3)=p(ii1,3)
	ps(4)=p(ii1,4)
	rs(1)=v(ii1,1)
	rs(2)=v(ii1,2)
	rs(3)=v(ii1,3)
	rs(4)=v(ii1,4)
	eg=ps(4)
        amd=0.00990   ! pymass(1)
        amu=0.00560   ! pymass(2)
        ams=0.199   ! pymass(3)
        amuu=2*amu
        amdd=2*amd
        amss=2*ams
c090505
        if(eg.lt.amuu)then   ! 2
        delte=eg   ! thrown away that gluon
c       move particle list,'pyjets',one step downward since ii1+1
        do j=ii1+1,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 600
c090505
        elseif(eg.lt.amdd)then   ! 2  090505
        kf1=2
        kf2=-2
	am1=amu
	am2=amu
        goto 700
        elseif(eg.ge.amdd .and. eg.lt.amss)then   ! 2
        kf1=2
        kf2=-2
	am1=amu
	am2=amu
        if(pyr(1).gt.0.5)then
        kf1=1
        kf2=-1
	am1=amd
	am2=amd
        endif
        goto 700
        elseif(eg.gt.amss)then   ! 2
        kf1=3
        kf2=-3
	am1=ams
	am2=ams
        rand=pyr(1)
        if(rand.gt.pros .and. rand.le.prods)then
        kf1=1
        kf2=-1
	am1=amd
	am2=amd
        endif
        if(rand.gt.prods)then
        kf1=2
        kf2=-2
	am1=amu
	am2=amu
        endif
	goto 700   ! 090505
	else   ! 2
c090505	goto 800   ! that gluon do not have enough energy to break
        endif   ! 2
700	continue 
c	exchange that gluon with the parton ahead
	j1=ii1+1
        do jj=1,5
	kk(jj)=k(j1,jj)
	pk(jj)=p(j1,jj)
	vk(jj)=v(j1,jj)
        k(j1,jj)=k(ii1,jj)
        p(j1,jj)=p(ii1,jj)
        v(j1,jj)=v(ii1,jj)
	k(ii1,jj)=kk(jj)
	p(ii1,jj)=pk(jj)
	v(ii1,jj)=vk(jj)
	enddo
c	write(22,*)'af. exchange event=',iiii
c	call pylist(1)
c       move particle list,'pyjets',one step forward since j1+1
        do j=n,j1+1,-1
        j2=j+1
        do jj=1,5
        k(j2,jj)=k(j,jj)
        p(j2,jj)=p(j,jj)
        v(j2,jj)=v(j,jj)
        enddo
        enddo
        n=n+1
c	write(22,*)'af. move forward event=',iiii
c	call pylist(1)
	k(j1,1)=2   ! A
	k(j1,2)=kf1
	k(j1,3)=0
        k(j1+1,1)=1   ! V
        k(j1+1,2)=kf2
        k(j1+1,3)=0
c       breaked q and qbar forms a string
c       give four momentum to breaked quarks
	decsuc=1   ! c1
	call decmom(ps,pp,am1,am2,decsuc)   ! c1
c	as mass of gluon from 'pyjets' may be negative it may be better
c	 (from energy conservation point of view) not using 'decmom' but
c	 random three momentum method if square root s less than 0.1
	if(decsuc.eq.0)then   ! c1
        do i4=1,3  
        pi=pyr(1)*ps(i4)
        pp(2,i4)=ps(i4)-pi
        pp(1,i4)=pi
        enddo
        pp11=pp(1,1)
        pp12=pp(1,2)
        pp13=pp(1,3)
	ampp=am1*am1+pp11*pp11+pp12*pp12+pp13*pp13
	if(ampp.gt.1.d40)ampp=1.d40
	if(ampp.lt.1.d-40)ampp=1.d-40
        pp(1,4)=dsqrt(ampp)
        pp21=pp(2,1)
        pp22=pp(2,2)
        pp23=pp(2,3)
	ampp=am2*am2+pp21*pp21+pp22*pp22+pp23*pp23
	if(ampp.gt.1.d40)ampp=1.d40
        if(ampp.lt.1.d-40)ampp=1.d-40
        pp(2,4)=dsqrt(ampp)
	endif   ! c1
	p(j1,1)=pp(1,1)
	p(j1,2)=pp(1,2)
	p(j1,3)=pp(1,3)
	p(j1,4)=pp(1,4)
	p(j1+1,1)=pp(2,1)
	p(j1+1,2)=pp(2,2)
	p(j1+1,3)=pp(2,3)
	p(j1+1,4)=pp(2,4)
c       give four coordinate to breaked quarks
c        first breaked quark takes the four coordinate of diquark
c        second breaked quark is arranged around first ones within
c        0.5 fm randumly in each of three coordinates and has same
c        fourth coordinate as diquark
	v(j1,1)=rs(1)
	v(j1,2)=rs(2)
	v(j1,3)=rs(3)
	v(j1,4)=rs(4)
	do j4=1,3
	rr(j4)=pyr(1)*0.5
	v(j1+1,j4)=rs(j4)+rr(j4)
	if(pyr(1).gt.0.5)v(j1+1,j4)=rs(j4)-rr(j4)
	enddo
	v(j1+1,4)=rs(4)
	delte=eg-p(j1,4)-p(j1+1,4)
c	write(9,*)'ii1,kf1,kf2,eg,delte=',iiglu(1),kf1,kf2,eg,delte
c	write(9,*)'pp(1,)=',(pp(1,ii),ii=1,4)
c	write(9,*)'pp(2,)=',(pp(2,ii),ii=1,4)
c	write(9,*)'sum=',(pp(1,ii)+pp(2,ii),ii=1,4)
c	write(9,*)'ps=',(ps(ii),ii=1,4)
c	write(9,*)'v(j1,)=',(rs(ii),ii=1,4)
c	write(22,*)'af. break g event=',iiii
c	call pylist(1)
	goto 600
	endif   ! 1
c	move particle list, 'pyjets', jglu steps forward since n to 1
800	do j=n,1,-1
	j1=j+jglu
	do jj=1,5
	k(j1,jj)=k(j,jj)
        p(j1,jj)=p(j,jj)
        v(j1,jj)=v(j,jj)
	enddo
	enddo
	n=n+jglu
	do j1=1,jglu
	iiglu(j1)=iiglu(j1)+jglu
	enddo
c	move g to the beginning of 'pyjets'
	jjj=0
	do j=jglu+1,n   ! do 1
	ik=k(j,2)
	if(ik.eq.21)then
	jjj=jjj+1
	do jj=1,5
	k(jjj,jj)=k(j,jj)
	p(jjj,jj)=p(j,jj)
        v(jjj,jj)=v(j,jj)
	enddo
	endif
	enddo   ! do 1
	do j2=1,jglu   ! do 2
	j11=iiglu(j2)
	do j1=j11+1,n
	do jj=1,5
        k(j1-1,jj)=k(j1,jj)
        p(j1-1,jj)=p(j1,jj)
        v(j1-1,jj)=v(j1,jj)
        enddo
	enddo
	n=n-1
	if(j2.lt.jglu)then
	do j3=j2+1,jglu
	iiglu(j3)=iiglu(j3)-1
	enddo
	endif
	enddo   ! do 2
c090505
c        mglu=jglu-1   ! jglu is odd
c        if(mod(jglu,2).eq.0)mglu=jglu   ! jglu is even
c        do j=1,mglu,2
c        k(j,1)=2   ! A
c        k(j1+1,1)=1   ! V
c        enddo
c        if(mod(jglu,2).ne.0)then
c        k(jglu-2,1)=2   ! A
c        k(jglu-1,1)=2   ! A
c        k(jglu,1)=1   ! V
c        endif
c090505
	k(jglu,1)=1   ! V
c       in order arranging gluons into a string, note: gluon in "pyjest"
c        always has k(i,1)=2 (i=1,2,...,jglu), so k(1,1)=2   ! A
c	write(22,*)'jglu,jjj,k(jglu,1),event=',jglu,jjj,k(jglu,1),iiii
c	call pylist(1)
600	jjj=0
	iii=1    
c	find out pair of 'A and V' composed of diquark (quark) and quark (diquark) 
	do i2=iii,n   ! 2
	kf=k(i2,2)
        kfab=iabs(kf)
        if(kfab.ne.2101 .and. kfab.ne.1103 .and. kfab.ne.2103
     c   .and. kfab.ne.2203 .and. kfab.gt.2)then   ! composed of u, d only
c	if(kfab.ne.2101 .and. kfab.ne.3101
c     c   .and. kfab.ne.3201 .and. kfab.ne.1103 .and. kfab.ne.2103
c     c   .and. kfab.ne.2203 .and. kfab.ne.3103 .and. kfab.ne.3203
cc     c   .and. kfab.ne.3303 .and. kfab.gt.10 .and. kfab.ne.21)then
c     c   .and. kfab.ne.3303 .and. kfab.gt.10)then
	iii=iii+1
	goto 500
	endif
	k1=k(i2,1)
c	if(k1.eq.2.and.kf.ne.21)then   ! k1=2 means 'A'  if 1
	if(k1.eq.2)then   ! k1=2 means 'A'  if 1
	do i3=i2+1,n   ! 3
	kf4=k(i3,2)
	kf4ab=iabs(kf4)
	k2=k(i3,1)
	if(k2.eq.1.and.((kfab.le.2.and.(kf4ab.eq.2101.or.kf4ab.eq.1103
     c   .or.kf4ab.eq.2103.or.kf4ab.eq.2203)).or.(kf4ab.le.2.and.(kfab
     c   .eq.2101.or.kfab.eq.1103.or.kfab.eq.2103.or.
     c   kfab.eq.2203))))then   ! k2=1 means 'V'  230407 if 2  
c230407        if(k2.eq.1.and.((kfab.le.2.and.kf4ab.gt.1000).or.
c230407     c   (kf4ab.le.2.and.kfab.gt.1000)))then   ! k2=1 means 'V'  if 2  
c	if(k2.eq.1.and.((kfab.lt.10.and.kf4ab.gt.1000).or.
c     c	 (kf4ab.lt.10.and.kfab.gt.1000)))then   ! k2=1 means 'V'  if 2    
c	'A and V' pair is diquark and quark (or quark and diquark) pair
c	if(k2.eq.1 .and. kf+kf4.ge.20)then
	jjj=jjj+1
        idi2(jjj,1)=i2
	idi2(jjj,2)=i3
	p1x=p(i2,1)
        p1y=p(i2,2)
        p1z=p(i2,3)
	p1e=p(i2,4)
	p2x=p(i3,1)
        p2y=p(i3,2)
        p2z=p(i3,3)
	p2e=p(i3,4)
	p12x=p1x+p2x
	p12y=p1y+p2y
	p12z=p1z+p2z
	p12e=p1e+p2e
	cm2=p12e*p12e-p12x*p12x-p12y*p12y-p12z*p12z
	if(cm2.gt.1.d40)cm2=1.d40
	if(cm2.lt.1.d-40)cm2=1.d-40
	cm=dsqrt(cm2)
c	compose diquark-quark 'A-V' pair into baryon which is inclueded 
c	 in PYTHIA as beam or target or is Delta
	if(kf.gt.10)then   ! i2 is diquark
	kfbb=kf/1000
	kf1=kfbb
	kf2=(kf-kfbb*1000)/100
	kf3=kf4
	sdir=dsign(1d0,p(i2,3))
	else   ! i2 is quark
        kfbb=kf4/1000
	kf1=kfbb
        kf2=(kf4-kfbb*1000)/100
        kf3=kf
	sdir=dsign(1d0,p(i3,3))
	endif
c	write(9,*)'iii,jjj,i3,p1e,p1x,p1y=',iii,jjj,i3,p1e,p1x,p1y
c	write(9,*)'p1z,p2e,p2x,p2y,p2z=',p1z,p2e,p2x,p2y,p2z
c	write(9,*)'p12e,p12x,p12y,p12z=',p12e,p12x,p12y,p12z
c	write(9,*)'kf,kf4,k1,k2,cm,=',kf,kf4,k1,k2,cm
c	write(9,*)'jjj,idi2=',jjj,idi2(jjj,1),idi2(jjj,2)
c	write(9,*)'kf1,2,3=',kf1,kf2,kf3
	call tabhb
c       compose diquark-quark 'A-V' pair into baryon
	if(kf1.gt.0.and.kf2.gt.0.and.kf3.gt.0)then
	call findb(kf1,kf2,kf3,cm,kfii,amasi,isucc,1)
	elseif(kf1.lt.0.and.kf2.lt.0.and.kf3.lt.0)then
	call findb(-kf1,-kf2,-kf3,cm,kfii,amasi,isucc,-1)   ! 020605 Tan
	else
	endif
c	write(9,*)'n,isucc,kfii,amasi=',n,isucc,kfii,amasi
c	if that baryon is inclueded in PYTHIA as beam or target, or that 
c	 baryon is Delta(0) or Delta(+) proceed then
	kiab=iabs(kfii)
c290805	if(isucc.eq.1.and.(kiab.eq.2212.or.kiab.eq.2112.or.kiab.eq.2214)
c290805     c	 )then   ! if 3 260805
c280805	if(isucc.eq.1.and.(kiab.eq.2212.or.kiab.eq.2112))then   ! if 3 260805
	if(isucc.eq.1.and.(kiab.eq.2212.or.kiab.eq.2112.or.kfii.eq.2114
     c	 .or.kfii.eq.2214))then   ! if 3
c	if(isucc.eq.1.and.(kiab.eq.2212.or.kiab.eq.2112.or.kfii.eq.3122
c     c	 .or.kfii.eq.3112.or.kfii.eq.3212.or.kfii.eq.3222.or.kfii.eq.
c     c	 3312.or.kfii.eq.3322.or.kfii.eq.3334.or.kfii.eq.2114
c     c	 .or.kfii.eq.2214))then   ! if 3
	isuc(jjj)=1
c	put that baryon on line i2 in 'pyjets' and 
c	 give proper variables to that baryon
	if(kfii.eq.2114)then   ! tread Detla(0) as neutron 
	kfii=2112
	amasi=0.940   ! pymass(2112)
	endif
	if(kfii.eq.2214)then   ! tread Detla(+) as proton
	kfii=2212
	amasi=0.938   ! pymass(2212)
	endif
	k(i2,1)=1
	k(i2,2)=kfii
	k(i2,3)=0
	p(i2,5)=amasi
c200405	p(i2,1)=p12x
c200405	p(i2,2)=p12y
c200405	p(i2,3)=p12z
c200405
300	call tdgaus(sigm2,ptmax,1,pppp)   ! 120505
	pi21=pppp(1,1)
	p(i2,1)=pi21
	pi22=pppp(1,2)
	p(i2,2)=pi22
	ppt=pi21*pi21+pi22*pi22
c	write(9,*)'ppt=',ppt   ! sa
	if(ppt.gt.0.1)goto 300   ! 120505 080805 010905
	ppr=p12x*p12x+p12y*p12y+p12z*p12z
	ppl=ppr-ppt
	if(ppl.gt.1.d40)ppl=1.d40
	if(ppl.lt.1.d-40)ppl=1.d-40
	p(i2,3)=dsqrt(ppl)*sdir
c130705	if(pyr(1).ge.0.5)p(i2,3)=-p(i2,3)   ! 030605
	pnnmm=amasi*amasi+ppr
	if(pnnmm.gt.1.d40)pnnmm=1.d40
        if(pnnmm.le.1.d-40)pnnmm=1.d-40
	pnnn=dsqrt(pnnmm)
c200405
	p(i2,4)=pnnn
	dele(jjj)=p12e-pnnn
c	write(9,*)'kfii,amasi,pnnn,dele(jjj)=',kfii,amasi,pnnn,dele(jjj)
c	write(9,*)'ppt,pt12=',ppt,dsqrt(p12x*p12x+p12y*p12y)   ! sa
	goto 888   ! 240805
	endif   ! if 3
c	if that baryon is Delta(-) or Delta(++) let it decays and put 
c	 decayed nucleon on i2 and decayed pion on i3 in 'pyjets'  
	if(isucc.eq.1.and.kfii.eq.1114)then   ! Delt(-) ->n+pi(-) if 4
c170405	ps(1)=p12x
c170405	ps(2)=p12y
c170405	ps(3)=p12z
c170405	ps(4)=p12e
	am1=0.940   ! pymass(2112)
	am2=0.140   ! pymass(-211)
c170405
301	call tdgaus(sigm2,ptmax,1,pppp)   ! 120505
	pi21=pppp(1,1)
	p(i2,1)=pi21
	pi22=pppp(1,2)
	p(i2,2)=pi22
	pi2t=pi21*pi21+pi22*pi22
c	write(9,*)'pi2t=',pi2t   ! sa
	if(pi2t.gt.0.1)goto 301   ! 120505 080805
302	call tdgaus(sigm2,ptmax,1,pppp)   ! 120505
	pi31=pppp(1,1)
	p(i3,1)=pi31
	pi32=pppp(1,2)
	p(i3,2)=pi32
	pi3t=pi31*pi31+pi32*pi32
c	write(9,*)'pi3t=',pi3t   ! sa
	if(pi3t.gt.1.0)goto 302   ! 120505 
	ppx=p(i2,1)+p(i3,1)
	ppy=p(i2,2)+p(i3,2)
	ppt=ppx*ppx+ppy*ppy
	ppl=(p12x*p12x+p12y*p12y+p12z*p12z-ppt)
	if(ppl.gt.1.d40)ppl=1.d40
	if(ppl.lt.1.d-40)ppl=1.d-40
	pplsr=dsqrt(ppl)*sdir   ! 130705
c130705	if(pyr(1).ge.0.5)pplsr=-pplsr   ! 030605
c030605	if(adj23.eq.1)then
c	call funcz(z1)
c	else
c	prr=am1*am1/4.+pi2t   ! 120505 pp2t originally
c	call pyzdis(kf1,kf2,prr,z1)
c030605	endif
	ppl1=pplsr*pyr(1)   ! 030605 z1*pplsr
	ppl2=pplsr-ppl1
        p(i2,3)=ppl1   ! 030605 *sdir
	p(i3,3)=ppl2
	pi24=am1*am1+pi2t+ppl1*ppl1
	if(pi24.gt.1.d40)pi24=1.d40
	if(pi24.lt.1.d-40)pi24=1.d-40
	p(i2,4)=dsqrt(pi24)
	pi34=am2*am2+pi3t+ppl2*ppl2
	if(pi34.gt.1.d40)pi34=1.d40
	if(pi34.lt.1.d-40)pi34=1.d-40
	p(i3,4)=dsqrt(pi34)
c170405
	k(i2,1)=1
	k(i2,2)=2112
	k(i2,3)=0
	p(i2,5)=am1
	k(i3,1)=1
	k(i3,2)=-211
	k(i3,3)=0
	p(i3,5)=am2
	dele(jjj)=p12e-p(i2,4)-p(i3,4)
c	write(9,*)'kfii,amasi,p(i2,4),p(i3,4),dele(jjj)=',kfii,amasi,
c     c	p(i2,4),p(i3,4),dele(jjj)
c	write(9,*)'pi2t,pi3t,pt12=',pi2t,pi3t,dsqrt(p12x*p12x+p12y*p12y)   ! sa
	endif   ! if 4
	if(isucc.eq.1.and.kfii.eq.2224)then   ! Delt(++) ->p+pi(+) if 5 
c170405	ps(1)=p12x
c170405	ps(2)=p12y
c170405	ps(3)=p12z
c170405	ps(4)=p12e
	am1=0.938   ! pymass(2212)
	am2=0.140   ! pymass(211)
c170405
303	call tdgaus(sigm2,ptmax,1,pppp)   ! 120505
	pi21=pppp(1,1)
	p(i2,1)=pi21
	pi22=pppp(1,2)
	p(i2,2)=pi22
	pi2t=pi21*pi21+pi22*pi22
	if(pi2t.gt.0.1)goto 303   ! 120505 080805
c	write(9,*)'pi2t=',pi2t   ! sa
304	call tdgaus(sigm2,ptmax,1,pppp)   ! 120505
	pi31=pppp(1,1)
	p(i3,1)=pi31
	pi32=pppp(1,2)
	p(i3,2)=pi32
	pi3t=pi31*pi31+pi32*pi32
c	write(9,*)'pi3t=',pi3t  ! sa
	if(pi3t.gt.1.0)goto 304   ! 120505 
	ppx=p(i2,1)+p(i3,1)
	ppy=p(i2,2)+p(i3,2)
	ppt=ppx*ppx+ppy*ppy
	ppl=(p12x*p12x+p12y*p12y+p12z*p12z-ppt)
	if(ppl.gt.1.d40)ppl=1.d40
	if(ppl.lt.1.d-40)ppl=1.d-40
	pplsr=dsqrt(ppl)*sdir   ! 130705
c130705	if(pyr(1).ge.0.5)pplsr=-pplsr   ! 030605
c030605	if(adj23.eq.1)then
c	call funcz(z1)
c	else
c	prr=am1*am1/4.+pi2t   ! 120505 pp2t originally
c	call pyzdis(kf1,kf2,prr,z1)
c030605	endif
	ppl1=pplsr*pyr(1)   ! 030605 z1*pplsr
	ppl2=pplsr-ppl1
        p(i2,3)=ppl1   ! 030605 *sdir
	p(i3,3)=ppl2
	pi24=am1*am1+pi2t+ppl1*ppl1
	if(pi24.gt.1.d40)pi24=1.d40
	if(pi24.lt.1.d-40)pi24=1.d-40
	p(i2,4)=dsqrt(pi24)
	pi34=am2*am2+pi3t+ppl2*ppl2
	if(pi34.gt.1.d40)pi34=1.d40
	if(pi34.lt.1.d-40)pi34=1.d-40
	p(i3,4)=dsqrt(pi34)
c170405
	k(i2,1)=1
	k(i2,2)=2212
	k(i2,3)=0
	p(i2,5)=am1
	k(i3,1)=1
	k(i3,2)=211
	k(i3,3)=0
	p(i3,5)=am2
	dele(jjj)=p12e-p(i2,4)-p(i3,4)
c	write(9,*)'kfii,amasi,p(i2,4),p(i3,4),dele(jjj)=',kfii,amasi,
c     c	p(i2,4),p(i3,4),dele(jjj)
c	write(9,*)'pi2t,pi3t,pt12=',pi2t,pi3t,dsqrt(p12x*p12x+p12y*p12y)   ! sa
	endif   ! if 5
888	iii=i3+1   ! 240805
	goto 400
	endif   ! if 2
        if(k2.eq.1)then
	iii=i3+1
        goto 400
        endif   
	enddo   ! 3
	endif   ! if 1
400	continue
500	enddo   ! 2
c	write(9,*)'jjj,isuc,idi2=',jjj,isuc(1),isuc(2),idi2(1,1),
c     c	 idi2(1,2),idi2(2,1),idi2(2,2)
c	write(22,*)'be. remove di event=',iiii
c	call pylist(1)
c230407	do ij=1,jjj
	do ij=jjj,1,-1   ! 230407
	isu=isuc(ij)
	if(isu.eq.1)then
	j1=idi2(ij,2)
c	move particle list,'pyjets',one step downward since j1+1
	do j=j1+1,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
c230407	if(ij.eq.1)idi2(ij+1,2)=idi2(ij+1,2)-1
	n=n-1
	endif
	enddo
c	share 'del' energy into particles
c230407
	delee=0.
	do ij=1,jjj
	delee=delee+dele(ij)
	enddo
c230407
c230407	del=delte+dele(1)+dele(2)
	del=delte+delee   ! 230407
	if(n.gt.0)then
	del=del/dfloat(n)
c	write(9,*)'share n,del=',n,del
	do j3=1,n
	p(j3,4)=p(j3,4)+del
	if(del.lt.0.)then
	if(p(j3,4).lt.0.)p(j3,4)=p(j3,4)-del
	pabs=dabs(p(j3,3))
	if(pabs.ge.p(j3,4))p(j3,4)=p(j3,4)-del
	endif
	enddo
	endif
c	write(22,*)'out of recons'
c	call pylist(1)
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
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 080104 220110
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
        k(n+1,1)=1
c221203	k(i1,3)=0
        k(n+1,3)=0
c221203
        k(n+1,4)=0
        k(n+1,5)=0
c221203
c080104
	ii=ii+1
        npt(ii)=n+1+naf
        ifcom(ii)=i1+naf   ! 220110
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
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)   ! 260503
c261108
        dimension pamass(3)
        kfab1=iabs(kf1)
        kfab2=iabs(kf2)

        pamass(1)=0.0099D0
        pamass(2)=0.0056D0
        pamass(3)=0.199D0

        if(kf1.le.3)then
          am1=pamass(kfab1)
        else
          am1=pymass(kf1)
        endif

        if(kf2.le.3)then
          am2=pamass(kfab2)
        else
          am2=pymass(kf2)
        endif
c261108
c        am1=pymass(kf1)
c        am2=pymass(kf2)
        pp(1,5)=am1
        pp(2,5)=am2
c       pp : four momentum of broken quarks, local variable 
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
        pp(1,4)=dsqrt(pp14)
c021005
	pp21=pp(2,1)
	pp22=pp(2,2)
	pp23=pp(2,3)
c021005
        pp24=am2*am2+pp21*pp21+pp22*pp22+pp23*pp23
        if(pp24.le.0.)pp24=1.e-20
        pp(2,4)=dsqrt(pp24)
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
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=40000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)   
c       calculate the E and |p| of broken quark in rest frame of diquark
        sm2=ps(4)*ps(4)-ps(1)*ps(1)-ps(2)*ps(2)-ps(3)*ps(3)
c       one problem here is that 'sm2' may not equal to square of diquark 
c	 (gluon) rest mass,'bream' is called for spliting g especially
c030603
c1	if(sm2.lt.1.d-10)then
c1	sm2=1.d-10
c1	endif
c       write(9,*)'in decmom sm2=',sm2   ! sa
c       write(9,*)'ps=',ps   ! sa
c030603
	if(sm2.lt.0.005)then   ! 110211
	decsuc=0   ! go back to random three momentum method
	return
	endif
        sm=dsqrt(sm2)   ! M (should be diquark mass)
c       pp(1,4)=(sm2-am2*am2+am1*am1)/2./sm
c       pp(2,4)=(sm2-am1*am1+am2*am2)/2./sm
        ppp=(sm2-(am1+am2)*(am1+am2))*(sm2-(am1-am2)*(am1-am2))
c161204	ppp=dabs(ppp)   ! 030603 ?
	if(ppp.lt.1.d-28)ppp=1.d-28   !161204
        ppp=dsqrt(ppp)/2./sm
c110211 goto 500   ! activate it for exponential cos(seta) distribution
c       the direction of broken quark is sampled isotropically in '4pi'
        coset=1.-2.*pyr(1)
        if(dabs(coset).gt.1.d0)then
        coset=coset/dabs(coset)
        endif
c021005
        siset=1.-coset*coset
        if(siset.lt.1.d-28)siset=1.d-28
c021005
        siset=dsqrt(siset)   ! 021005
100     cosi1=pyr(1)
        cosi12=cosi1*cosi1
        eta2=2.*pyr(1)-1.
        eta22=eta2*eta2
        coseta=cosi12+eta22
        if(coseta.gt.1.)goto 100
        if(coseta.lt.1.d-28)coseta=1.d-28
        cofi=(cosi12-eta22)/coseta
        sifi=2.*cosi1*eta2/coseta
        goto 600
500     continue
c       cos(seta) is sampled from exponential distribution when
c        0<seta<pi/2 and its absolute value is assumed to be symmetry
c        about seta=pi/2. 'fi' is assumed to be isotropic in 2pi
        coset=dlog(1.d0+1.7183*pyr(1))
        if(pyr(1).lt.0.5d0)coset=-coset
c021005
        siset=1.-coset*coset
        if(siset.lt.1.d-28)siset=1.d-28
c021005
        siset=dsqrt(siset)
        fi=2.*3.1416*pyr(1)
        cofi=dcos(fi)
        sifi=dsin(fi)
600     continue
        pi(1)=ppp*siset*cofi
        pi(2)=ppp*siset*sifi
        pi(3)=ppp*coset
c021005
        pi4=ppp*ppp+am1*am1
        if(pi4.lt.1.d-28)pi4=1.d-28
c021005
        pi(4)=dsqrt(pi4)
        pj(1)=-pi(1)
        pj(2)=-pi(2)
        pj(3)=-pi(3)
c021005
        pj4=ppp*ppp+am2*am2
        if(pj4.lt.1.d-28)pj4=1.d-28
c021005
        pj(4)=dsqrt(pj4)
c       write(9,*)'before rotation'   ! sa
c       write(9,*)ppp,(pi(i),i=1,4)   ! sa
c       write(9,*)(pj(i),i=1,4)   ! sa
c050603
c       rotate to the frame where diquark (gluon), ps, is described
c       calculate the direction cosines of ps
        fi1=pyangl(ps(1),ps(2))
c021005
        ps12=ps(1)*ps(1)+ps(2)*ps(2)
        if(ps12.lt.1.d-28)ps12=1.d-28
        ps12=dsqrt(ps12)
c021005
        cta1=pyangl(ps(3),ps12)
        cfi1=dcos(fi1)
        sfi1=dsin(fi1)
        ccta1=dcos(cta1)
        scta1=dsin(cta1)
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
        if(pj4.lt.1.d-28)pj4=1.d-28
        pi(4)=dsqrt(pi4)
        pj(4)=dsqrt(pj4)
c021005
c       write(9,*)'after rotation'   ! sa
c       write(9,*)(pi(i),i=1,4)   ! sa
c       write(9,*)(pj(i),i=1,4)   ! sa
c       boost to moving frame of diquark
        ee=ps(4)
        if(ee.lt.1.d-14)ee=1.d-14   ! 021005
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
        if(dabs(1.-fr) .le. dep)goto 200
        do i=1,np
        amas=pp(i,5)
        ppm=pp(i,4)/amas
        ppf=ppm/fr
        ff(i)=dsqrt(dabs(ppf*ppf-1.d0)/(ppm*ppm-1.))
        do j=1,3
        ppp=ff(i)*pp(i,j)
        pp(i,j)=ppp
        pxyz(j)=pxyz(j)+ppp
        enddo
        enddo
        do i=1,3
        arp(i)=dabs(1.-pxyz(i)/ps(i))
        pxyz(i)=pxyz(i)-ps(i)
        enddo
        if(dabs(1.-fr).le.dep .and. arp(1).le.dep .and. arp(2).le.dep
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
        pp(i,4)=dsqrt(pp52+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
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
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension rr(3)
        do i1=1,3
c261002        rr(i1)=pyr(1)*v(ii,i1)
        rr(i1)=pyr(1)*0.5   ! 261002
        v(n+1,i1)=v(ii,i1)+rr(i1)
        if(pyr(1).gt.0.5d0)v(n+1,i1)=v(ii,i1)-rr(i1)
        enddo
        v(n+1,4)=v(ii,4)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine copl(tt)
c	calculate coordinate of center of mass of non-freeze-out system
c	position of a particle, checking is it freezes out or not, is 
c	 calculated with respect to this origin.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=40000)
        COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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
	if(iabs(kf).eq.213 .or. kf.eq.113)amass1=p(ii,5)   ! 010600 
	if((iabs(kf).eq.213 .or. kf.eq.113) .and. dabs(amass-amass1)
     &	 .gt.0.001d0)amass=amass1   ! 010600 
	samass=samass+amass
	do 100 jj=1,3
	coor(jj)=coor(jj)+amass*v(ii,jj)
100	continue
110	continue
	do ii=1,3
	coor(ii)=coor(ii)/dmax1(0.14d0,samass)
	enddo
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ctlcre(lc,tc,tw)
c	create initial collision list  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(nsize=240000)
      	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non10,
     c   disbe(100,100)
	common/ctllist/nctl,noinel(600),nctl0,noel
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c	nap,nat,nzp,nzt,pio
	dimension lc(nsize,5),tc(nsize),tw(nsize)
c081010	time=0.
	nctl=1
	dminf=100.
	nzpab=iabs(nzp)   ! in order to consider ppbar or pbarp
	nztab=iabs(nzt)
	nzpt=nzpab+nztab
	napt=nap+nat
	do 10 l=1,nzpab   ! projectile proton or e-   ! 060813
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
10      continue
        do 20 l=nzpt+1,nap+nztab   ! projectile neutron
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
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	PARAMETER (kszj=40000,KSZ1=30)
	parameter(nsize=240000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
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
	dmin=dsqrt(sg)
	if(dmin.lt.dminf)then
	dminf=dmin
	if=i
	jf=j
	endif
	if(ipden.ne.2 .and. dmin.gt.ecsnn)return   ! 060813
        if(ipden.eq.2 .and. dmin.gt.ecsen)return   ! 060813
c	distance between the two particles should be smaller than ecsnn (ecsen)
c        060813
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
	drmax=rao*dmax1(rnt,rnp)
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
c	find out the binary collision with minimum collision time
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


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updtlp(time,lc,tc,tw,iii)
c	update collision list after calling 'pythia' successfully
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000,KSZ1=30)
        parameter(nsize=240000)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
	common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa14/ipyth(2000),idec(2000),iwide(2000)
c010530        common/sa19/kji   ! 16/09/99
        COMMON/SBH/N,NON,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/ctllist/nctl,noinel(600),nctl0,noel
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        dimension lc(nsize,5),tc(nsize),tw(nsize)
c	ipyth: store line number of  hardons from calling 'pythia' 
c	write(9,*)'in updtlp nctl=',nctl
c	do i=1,nctl
c	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
c	enddo
c       loop over old colliding pairs
        j=0
	if(nctl.eq.0)goto 200
        do 400 i=1,nctl
        if((tc(i)-time).le.ddt) goto 400
c       through away the pair whih tc<= time
        j=j+1
        tc(j)=tc(i)
        tw(j)=tw(i)
        do m=1,5
        lc(j,m)=lc(i,m)
        enddo
400     continue
        do i=j+1,nctl+1
        tc(i)=0.0
        tw(i)=0.0
        do m=1,5
        lc(i,m)=0
        enddo
        enddo

c	write(9,*)'updtlp rmove nctl,nbh,m4=',j,n,numb(4)   ! sa
c	do i=1,j
c	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
c	enddo

200	nctl=j+1
c060813	loop over particle list for each generated particle from pythia

	m2=numb(2)   ! 060813
        m4=numb(4)
	m7=numb(7)   ! 241110
c	m9=numb(9)
c	m17=numb(17)
c	m19=numb(19)
c	m25=numb(25)
c	m29=numb(29)
c	m32=numb(32)
c	m34=numb(34)
c        m34=numb(kfmax-11)
c       subtract 11, since we do not consider the rescattering of x0c, etc
	do j11=1,n
	j1=ipyth(j11)
        kfj=ksa(j1,2)   
c	write(9,*)'n,j11,j1,kf,m4=',n,j11,j1,kfj,m4   
	if(kfj.ne.11.and.j1.gt.m2)goto 300   ! 241110 060813 m7 to m2
c060813	consider only the reinteraction among nucleons & nucleon with e-
c060813 loop over particle list
c	mm=m34   
c060813	mm=m7      ! 241110
        mm=m2   ! 130913 m7 to m2
c060813	consider only the reinteraction j11 with nucleons
	do i=1,mm   

        if(nctl.gt.nsize)then
	write(9,*)'iiii,nsize,n,nctl=',iiii,nsize,n,nctl   ! sa
	stop 22222
	endif
c010600
	do j22=1,n
	j2=ipyth(j22)
	if(i.eq.j2)goto 600
	enddo
c010600

	i1=i
	kfi=ksa(i1,2)
c	write(9,*)'j1,i1,kfj,kfi=',j1,i1,kfj,kfi   ! sa
        iflag=0
        call rsfilt(j1,i1,iflag)
        if(iflag.eq.0)goto 100
        tc(nctl)=0.0
c011204	call tcolij(i1,j1,time,nctl,lc,tc,tw)
        call tcolij(j1,i1,time,nctl,lc,tc,tw)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rsfilt(l,l1,iflag)
c	 subroutine rsfilt plays the role of first range filter 
c	 subroutine intdis plays the role of second range filter
c       collision pairs not interested can not filter through both of rsfilt 
c        and intdis
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000,KSZ1=30)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	kl=k(l,2)
	kl1=k(l1,2)
	klab=iabs(kl)
	kl1ab=iabs(kl1)
	if(l.eq.l1) goto 10
	if(ishp(l).eq.0.or.ishp(l1).eq.0) goto 10

c060813 consider nn collision 
	if(klab.eq.2212.and.(kl1ab.eq.2112.or.kl1ab.eq.2212))goto 11
	if(klab.eq.2112.and.(kl1ab.eq.2112.or.kl1ab.eq.2212))goto 11
c060813 consider ep (en) and pe (ne)
        if(kl.eq.11.and.(kl1.eq.2112.or.kl1.eq.2212))goto 11
        if((kl.eq.2112.or.kl.eq.2212).and.kl1.eq.11)goto 11
c060813
	goto 10

11	iflag=1
10	continue
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
      	COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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
	pi(4)=dsqrt(bmi**2+pi(1)**2+pi(2)**2+pi(3)**2)
	pj(4)=dsqrt(bmj**2+pj(1)**2+pj(2)**2+pj(3)**2)
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

	dmin=dsqrt(sg)
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
	rri=dsqrt(ri(1)*ri(1)+ri(2)*ri(2)+ri(3)*ri(3))
	rrj=dsqrt(rj(1)*rj(1)+rj(2)*rj(2)+rj(3)*rj(3))
c	the rnt in rao*max(rnt,rnp)+rnt is due to the fact that
c        we could not know the postion of the mass-center in the future
	rrr=rao*dmax1(rnt,rnp)
	if(ifram.eq.0)rrr=rao*dmax1(rnt,rnp)+rnt
c	if(dabs(k(l,2)).gt.1000.or.dabs(k(l1,2)).gt.1000)rrr=1.E+10*rrr
	if(rri.gt.rrr)goto 10
	if(rrj.gt.rrr)goto 10
c	particles under consideration must be still within considered region
c	 when the collision happens
	if(tcol.le.rrr)goto 18   ! 041204
        return   ! 041204
18	tc(icp)=tcol
	lc(icp,1)=l
	lc(icp,2)=l1
10	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine intdis(l,l1,ss,rsig)
c	calculate interaction distance between particles l and l1.
c	It plays also the role of second range filter
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000)
        COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
	rsig=0.
	kl=k(l,2)
	kl1=k(l1,2)

	if(iabs(kl).eq.2212 .or. iabs(kl).eq.2112)idpl=1
        if(iabs(kl1).eq.2212 .or. iabs(kl1).eq.2112)idpl1=1
	if(kl.eq.11)idpl=8   ! e- 060813
	if(kl1.eq.11)idpl1=8   ! 060813
        if(idpl.eq.1 .and. idpl1.eq.1)rsig=ecsnn
	if(idpl.eq.8 .and. idpl1.eq.1)rsig=ecsen   ! 060813
        if(idpl.eq.1 .and. idpl1.eq.8)rsig=ecsen   ! 060813
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updatl(ic,jc,time,lc,tc,tw,iii)
c	update collision time list after elastic scattering
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000,KSZ1=30)
	parameter(nsize=240000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/ctllist/nctl,noinel(600),nctl0,noel
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	dimension lc(nsize,5),tc(nsize),tw(nsize)

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

	nctl=j+1
c	loop over particle list

        m2=numb(2)   ! 130913
        m4=numb(4)
        m7=numb(7)   ! 241110
	j1=ic
	do ik=1,2

c010530
	kf=ksa(j1,2)
	if(kf.ne.11.and.j1.gt.m2)goto 300   ! 241110 060813 m7 to m2
c130913	consider only the reinteraction among nucleons & nucleon with e- 060813 
	mm=m2   ! 241110 060813 m7 to m2 
	do i=1,mm
	if(j1.eq.ic .and. i.eq.jc)goto 600 
	if(j1.eq.jc .and. i.eq.ic)goto 600
c	forbiden scattered particles colliding with each other
	if(nctl.gt.nsize)then
        write(MSTU(11),*)'size of array "nsize" needs to be extended'
        write(MSTU(11),*)'error is serious,stop running'
        stop 22222
        endif

 	i1=i
	iflag=0
	call rsfilt(j1,i1,iflag)
	if(iflag.eq.0)goto 100
	tc(nctl)=0.0
	call tcolij(i1,j1,time,nctl,lc,tc,tw)
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pauli(ii,ppaul)
c	calculate the unoccupation probability (ppaul) of particle ii
c	 in 'pyjets'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=40000)
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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
c	pick up the partons with same flavor as ii from 'pyjets' and 
c	 'saf' 
	nkk=0 
c	the new produced partons, except ii, in 'pyjets' are also to be
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
	dxx=dabs(dxx)
	dyy=dabs(dyy)
	dzz=dabs(dzz)
	dpx=dabs(dpx)
	dpy=dabs(dpy)
	dpz=dabs(dpz)
	if(dxx.le.0.5.and.dyy.le.0.5.and.dzz.le.0.5.and.
     c	 dpx.le.0.62.and.dpy.le.0.62.and.dpz.le.0.62)anq=anq+1.
	enddo
	proba=anq/6.
c	6=2*3, spin and colour degeneracies of quark (antiquark)
	ppaul=1.-proba
	return
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
c       move ii-th particle (e-) in pyjets to first position
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
c       move particle list (pyjets) one step forward from ii-1 to 1
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
	COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
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
     &        6.0,3.0,6*0/ 	
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
