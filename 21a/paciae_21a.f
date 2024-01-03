        program paciae_20a  
c	for hadron-hadron collisions and e+e annihilation 
c	it should be run together with parcas_20a.f,coales_20a.f,
c        sfm_20a.f, and p20a.f
c       paciae_20a.f: user program including generation of partonic initial 
c        state etc.
c       parcas_20a.f: perform parton rescattering, where only 2->2 processes
c        are considered and LO pQCD cross section or its regularized
c        approximation is used
c       sfm_20a.f: hadronization according to LUND string fragmentation model
c       coales_20a.f: hadronization according to Monte Carlo coalescence model
c       selecting sfm_20a.f or coales_20a.f by parameter adj1(12) 
c       hadcas_20a.f: hadronic rescattering
c       note: dimension of arraies in common block 'pyjets' should be the 
c	 same among all programs
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=40000,mplis=40000)   ! 061007
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 161007
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)   ! 061007
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)   ! 151107
	common/sa1/kjp21,non1,bp,i,neve,nout,nosc   ! 061007
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2),ifram
c061007
	common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23   ! 171108  
	common/sa24/adj1(40),nnstop,non24,zstop 
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   ! 151107  
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 151107 160110
c051108
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21
        common/sa28/nstr,nstr00,nstra(kszj),nstrv(kszj)   ! 160110
        common/sa29/effk1,lcub   ! 160110
        common/sa33/smadel,ecce,parecc,iparres   ! 220312 240412
c220312	smadel: small perpurbation of ellipse from circle
c240412 iparres: =0 consider ela. parton-parton collisions only
c240412 iparres: =1 otherwise
c       effk1,parj1,parj2,parj3,and parj21 are the tuned effective string
c        tension,parj(1),(2),(3),and (21),respectively
c       lcub: size of cubic used in energy density calculation
        dimension skapa(5),skapao(5)
c       itime: the number of strings in current event
c       astr: the number of strings in current event
c       gtime: the number of gluon in current event
c       akapa(1): sum of string tension over strings in the current event
c       akapa(2): sum of parj(2) over strings in the current event
c       akapa(3): sum of parj(21) over strings in the current event
c       akapa(4): sum of parj(1) over strings in the current event
c       akapa(5): sum of parj(3) over strings in the current event
c051108
	common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
c       iprl: current total number of partons in particle list
c       rp  : position and time coordinates of particle
c       pp  : momentum of particle
c       tp  : current time of particle
c       taup: format time of particle
c       idp : flavor of particle
c       rmp : rest mass of particle
c       ep  : total energy of particle
c       vp  : velocity of particle
	common/show/vip(mplis),xap(mplis)
	common/ctllist_p/nreac(9),nrel   ! 051207
c	nreac(i): statistics of # of successful i-th collision
c	nrel: statistics of # of collision blocked   
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
c       ithroq_p : total # of quarks thrown away
c       ithrob_p : total # of antiquarks thrown away
c       throe_p : total momentum and energy of the partons thrown away
c       ichh_p : total charge of the partons thrown away
c061007
c161007
        common/sa1_h/nn4,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)   ! 171108
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,iifram,
     &  iabsb,iabsm,non10,csspn,csspm   ! 171108
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_t/ithroq_t,ithrob_t,ich_t,non6_t,throe_t(4)
        common/wz/coor(3)   ! 131108
c       ithroq_t : total # of t quarks thrown away
c       ithrob_t : total # of t antiquarks thrown away
c       throe_t : total momentum and energy of the t quarks (t antiquarks) 
c        thrown away
c       ichh_t : total charge of t quarks (t antiquarks) thrown away
c       ithroq: total # of quarks thrown away in coales.f
c       ithrob: total # of antiquarks thrown away in coales.f
c       throe: total momentum and energy of the partons thrown away in coales.f
c       ichh: total charge of the partons thrown away in coales.f
c161007
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio   ! 141208
c220312
        common/sa18_pt/snum_pt(21,20),v1_pt(21,20),v2_pt(21,20),
     c   v12_pt(21,20),v22_pt(21,20)   
        common/sa18_eta/snum_eta(21,20),v1_eta(21,20),v2_eta(21,20),
     c   v12_eta(21,20),v22_eta(21,20)
 
        dimension ssnum_pt(21,20),sv1_pt(21,20),sv2_pt(21,20),
     c   sv12_pt(21,20),sv22_pt(21,20),sv1o_pt(21,20),sv2o_pt(21,20),
     c   sv12o_pt(21,20),sv22o_pt(21,20)   
        dimension ssnum_eta(21,20),sv1_eta(21,20),sv2_eta(21,20),
     c   sv12_eta(21,20),sv22_eta(21,20),sv1o_eta(21,20),sv2o_eta(21,20)
     c   ,sv12o_eta(21,20),sv22o_eta(21,20)   

        dimension ssv1_pt(21,20),ssv2_pt(21,20),ssv12_pt(21,20),
     c   ssv22_pt(21,20),ssv1o_pt(21,20),ssv2o_pt(21,20),
     c   ssv12o_pt(21,20),ssv22o_pt(21,20)
        dimension ssv1_eta(21,20),ssv2_eta(21,20),ssv12_eta(21,20),
     c   ssv22_eta(21,20),ssv1o_eta(21,20),ssv2o_eta(21,20),
     c   ssv12o_eta(21,20),ssv22o_eta(21,20)  
c220312
c010412
        common/fsa18_pt/fsnum_pt(21,20),fv1_pt(21,20),fv2_pt(21,20),
     c   fv12_pt(21,20),fv22_pt(21,20)
        common/fsa18_eta/fsnum_eta(21,20),fv1_eta(21,20),fv2_eta(21,20)
     c   ,fv12_eta(21,20),fv22_eta(21,20)

        dimension fssnum_pt(21,20),fsv1_pt(21,20),fsv2_pt(21,20),
     c   fsv12_pt(21,20),fsv22_pt(21,20),fsv1o_pt(21,20),
     c   fsv2o_pt(21,20),fsv12o_pt(21,20),fsv22o_pt(21,20)
        dimension fssnum_eta(21,20),fsv1_eta(21,20),fsv2_eta(21,20),
     c   fsv12_eta(21,20),fsv22_eta(21,20),fsv1o_eta(21,20),
     c   fsv2o_eta(21,20),fsv12o_eta(21,20),fsv22o_eta(21,20)

        dimension fssv1_pt(21,20),fssv2_pt(21,20),fssv12_pt(21,20),
     c   fssv22_pt(21,20),fssv1o_pt(21,20),fssv2o_pt(21,20),
     c   fssv12o_pt(21,20),fssv22o_pt(21,20)
        dimension fssv1_eta(21,20),fssv2_eta(21,20),fssv12_eta(21,20),
     c   fssv22_eta(21,20),fssv1o_eta(21,20),fssv2o_eta(21,20),
     c   fssv12o_eta(21,20),fssv22o_eta(21,20)
c010412
        dimension an(20,5,20),bn(20),san(20,5,20),sbn(20),c(5)
        dimension anf(20,5,20),bnf(20),sanf(20,5,20),sbnf(20)
	dimension sao(20,5,20),sbo(20),saof(20,5,20),sbof(20)
        dimension sthroe(4),wthroe(4)   ! 161007
        dimension kdiq(kszj,5),dgmas(kszj),snreac(9)   ! 151107 051207
        dimension nreaco(9)   ! 160110
c       ispmax : maximum # of the different kind of particles wanted
c        to statistics
c       isdmax : maximum # of distributions wanted to statistics
c       ispkf(i): flavor code of i-th particle wanted to statistics
c       asd(i): interval of the i-th variable
c       array an(l,i,j) (san(l,i,j)):
c        l: value of distribution argument (e. g. value of y or pt)
c        i: identify the distribution (e. g. i=2 is y distribution)
c        j: order # of kf code of particle wanded to statistics (1 to ispmax)
c       iflmax : maximum # of filters,it equal 0 means no filter at all
c       afl(j,i,1): lower limit of i-th filter for j-th particle
c       afl(j,i,2): upper limit of i-th filter for j-th particle
c        i=1 : y filter
c        i=2 : pt filter
c         .        .
c         .        .
c         .        .

c151107
c       common block 'sbe': store parton (q,qq,and g) configuration before 
c        breaking diquark
c       idi: counts cumunatively the number of diquark (anti-diquark)
c       idio: value of idi after last nn collision
c       ndiq(j): = 0 if j is quark (antiquark)
c                = idi if j is diquark (anti-diquark)
c       note: j is line number in 'sbe' 
c       npt(idi): = line number of idi-th diquark (anti-diquark)
c        partner in 'pyjets'
c151107
c160110 ifcom(idi): line number of first component of idi-th diquark
c       nstr: statitics of number of strings in a hh collis.
c       nstr00: number of strings after call remo
c       nstra(i): line number of first component in i-th string
c160110 nstrv(i): line number of last component in i-th string
c       not write header 
        mstp(127)=0
        mstp(122)=0
        mstu(11)=22
	open(5,file='usux.dat',status='old')
        open(9,file='rms0.out',status='unknown')
	open(22,file='main.out',status='unknown')
c       open(88,file='outayut.out',status='unknown')
c	open(99,file='check.out',status='unknown')

	read(5,*)nn,nout,ee
	read(5,*)nap,nzp,nat,nzt
	read(5,*)mstu21,mstj1_1,mstj1_2,mstj2,mstj3,mstp33   ! 051207
	read(5,*)ifram,ispmax,isdmax,iflmax,nchan,itden
	read(5,*)(ispkf(i),i=1,10)
        read(5,*)(ispkf(i),i=11,ispmax)
	read(5,*)(asd(i),i=1,isdmax)
        if(iflmax.eq.0)goto 200
        do kk=1,ispmax
        do i=1,iflmax
        read(5,*)(afl(kk,i,j),j=1,2)
        enddo
        enddo
c061007
200	read(5,*)(adj1(i),i=1,10)
        read(5,*)(adj1(i),i=11,20)
        read(5,*)(adj1(i),i=21,30)
        read(5,*)(adj1(i),i=31,40)
c061007
c051108
        read(5,*)effk1,kjp20,kjp21,kjp22,lcub   ! 171108
        read(5,*)ttaup,taujp,iabsb,iabsm   ! 171108
	read(5,*)smadel,iparres   ! 220312 240412
c171108 kjp21 = 0 (1): without (with) hadronic rescattering   
        close(5)   
c       kjp22 = 0 : constant string tension with string tension calculation 
c       kjp22 = 1 : variable effective string tension
c       kjp22 = 2 : original case (constant string tension)
        parj(2)=adj1(32)
        parj2=parj(2)
c       suppression of s quark pair production (D=0.3) relative to u (d)
        parj(1)=adj1(31)
        parj(3)=adj1(33)
        parj(21)=adj1(34)
        parj1=parj(1)
        parj3=parj(3)
        parj21=parj(21)
c051108 mstp(91)=adj1(35)   ! 220808
c051108 if(mstp(91).eq.1)parp(91)=adj1(39)   ! 220808
c051108 if(mstp(91).eq.2)parp(92)=adj1(39)   ! 220808
c051108
        adj112=adj1(12)   ! 151107
	adj140=adj1(40)   ! 061007
c	nn: number of run
c	ee: cms energy if ifram=1
c	    incident momentum if ifram=0
c       adj140=1: run terminated after partonic initiation
c	adj140=2: run terminated after parton scattering 
c       adj140=3: hadronization followes parton initiation 
c       adj140=4: pythia 
c       adj140=5: hadronization followes parton scattering 
c       itden=0: for pp,ppbar,or pbarp
c       itden=2: for e+e- 
c061007
	do i1=1,5
        ppsa(i1)=0.
        enddo
	bp=1.   ! impact parameter
c061007
	mstu(21)=mstu21   ! 051207
c	decide the coefficients in Lund string fragmentation function	
	parj(41)=adj1(6)   ! 151107
	parj(42)=adj1(7)   ! 151107
c	decide the k factor in primary parton-parton hard scattering
	mstp(33)=mstp33
        parp(31)=adj1(10)   ! 151107

        if(nchan.eq.0)then
c       inelastic (INEL)
        msel=0
        msub(11)=1
        msub(12)=1
        msub(13)=1
        msub(28)=1
        msub(53)=1
        msub(68)=1
c       msub(91)=1
        msub(92)=1
        msub(93)=1
        msub(94)=1
        msub(95)=1
	elseif(nchan.eq.1)then
c       Non Single Diffractive (NSD)
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
c       qqb --> gamma^*/Z^0 used to generate Drell-Yen
        msel=0
        msub(1)=1
        elseif(nchan.eq.3)then
c       J/psi production
        msel=0
        msub(86)=1
c        msub(87)=1
c        msub(88)=1
c        msub(89)=1
        elseif(nchan.eq.4)then
c       heavy-flavor production
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
c	msub(95)=1
	elseif(nchan.eq.6)then
c       soft only
	msel=0
	msub(92)=1
	msub(93)=1
	msub(94)=1
	msub(95)=1
	msub(96)=0
	else
	msel=1
c       pythia default
	endif
c	forbiden decay of particle, if mdcy(...)=0 
	mdcy(pycomp(211),1)=0   ! pion+
	mdcy(pycomp(-211),1)=0   ! pion-
c	mdcy(pycomp(111),1)=0   ! pion0   ! 231006
c	mdcy(pycomp(221),1)=0   ! eta
c	mdcy(pycomp(331),1)=0   ! eta'
	mdcy(pycomp(310),1)=0   ! k0_S
	mdcy(pycomp(321),1)=0   ! kaon+
	mdcy(pycomp(-321),1)=0   ! kaon-
c	mdcy(pycomp(311),1)=0   ! kaon0
c	mdcy(pycomp(-311),1)=0   ! kaon0-
c	mdcy(pycomp(323),1)=0   ! kaon*+
c	mdcy(pycomp(-323),1)=0   ! kaon*-
c	mdcy(pycomp(313),1)=0   ! kaon*0
c	mdcy(pycomp(-313),1)=0
	mdcy(pycomp(2212),1)=0   ! proton
	mdcy(pycomp(-2212),1)=0   ! antiproton
	mdcy(pycomp(2112),1)=0   ! neutron
	mdcy(pycomp(-2112),1)=0   ! antineutron
        mdcy(pycomp(333),1)=0   ! phi
        mdcy(pycomp(3122),1)=0   ! Lambda0
        mdcy(pycomp(-3122),1)=0
	mdcy(pycomp(3312),1)=0   ! Xi-
	mdcy(pycomp(-3312),1)=0
        mdcy(pycomp(3334),1)=0   ! Omega-
        mdcy(pycomp(-3334),1)=0   ! 
c        mdcy(pycomp(3212),1)=0   ! Sigma0
c        mdcy(pycomp(-3212),1)=0
c080410	mdcy(pycomp(3112),1)=0   ! Sigma-
c080410	mdcy(pycomp(3222),1)=0   ! Sigma+

        nminf=0
	npluf=0
        nmin=0
	nplu=0
        do i1=1,20
        sbn(i1)=0.
        sbnf(i1)=0.
        do i2=1,5
        do i3=1,20
        san(i1,i2,i3)=0.
        sanf(i1,i2,i3)=0.
        enddo
        enddo
        enddo
        
c051108
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        mkapa=0
        sgtime=0.
        sgtimeo=0.
        sastr=0.
        do i1=1,5
        skapa(i1)=0.
        skapao(i1)=0.
        enddo
        endif
c051108
c161007

	sthroq=0.
        sthrob=0.
        sthroc=0.
        do i1=1,4
        sthroe(i1)=0.
        enddo
c161007
c160110
        nrel=0
        nrea=0
        do i=1,9
        nreac(i)=0
        enddo
c160110
        sedi=0.   ! 131108
        spdi=0.   ! 131108
        ssdi=0.   ! 131108
        stime_par=0.   ! 171108
        stime_had=0.   ! 171108
        stimep=0.d0   ! 280910
        stimeh=0.d0   ! 280910
c131108 statistics of energy density of charged particle (parton or hadron)
	if((itden.eq.0) .and. (nzp.eq.1 .and. nzt.eq.1))
     c   call pyinit('cms','p','p',ee)   ! 110607
	if((itden.eq.0) .and. (nzp.eq.1 .and. nzt.eq.-1))
     c   call pyinit('cms','p','pbar',ee)   ! 110607
	if((itden.eq.0) .and. (nzp.eq.-1 .and. nzt.eq.1))
     c   call pyinit('cms','pbar','p',ee)   ! 110607
        write(9,*)'nn,nout,ee=',nn,nout,ee
	write(9,*)'nap,nzp,nat,nzt=',nap,nzp,nat,nzt
	write(9,*)'mstu21,mstj1_1,mstj1_2,mstj2,mstj3,mstp33=',
     c	 mstu21,mstj1_1,mstj1_2,mstj2,mstj3,mstp33   ! 051207
        write(9,*)'ifram,ispmax,isdmax,nchan,adj140,itden=',
     c   ifram,ispmax,isdmax,nchan,adj140,itden
        write(9,*)'ispkf=',(ispkf(i),i=1,10)
        write(9,*)(ispkf(i),i=11,ispmax)
        write(9,*)'asd=',(asd(i),i=1,isdmax)
        write(9,*)'effk1,kjp20,kjp21,kjp22,lcub=',effk1,kjp20,
     c   kjp21,kjp22,lcub   !171108
        write(9,*)'par1,par2,par3,par21=',parj(1),parj(2),parj(3),
     c   parj(21)   ! 171108
        write(9,*)'ttaup,taujp,iabsb,iabsm=',ttaup,taujp,iabsb,iabsm ! 171108
        do kk=1,ispmax
        do i=1,iflmax
        write(9,*)(afl(kk,i,j),j=1,2)
        enddo
        enddo
c091007
        write(9,*)'adj1=',(adj1(i),i=1,10)
        write(9,*)'adj1=',(adj1(i),i=11,20)
        write(9,*)'adj1=',(adj1(i),i=21,30)
        write(9,*)'adj1=',(adj1(i),i=31,40)
	write(9,*)'smadel,iparres=',smadel,iparres   ! 220312 240412
c091007
c061007	do 300 i=1,nn
c220312
        do i2=1,21
        do i1=1,20
        ssnum_pt(i2,i1)=0.   
        ssnum_eta(i2,i1)=0.   
        sv1_pt(i2,i1)=0.
        sv2_pt(i2,i1)=0.
        sv1o_pt(i2,i1)=0.
        sv2o_pt(i2,i1)=0.
        sv12_pt(i2,i1)=0.
        sv22_pt(i2,i1)=0.
        sv12o_pt(i2,i1)=0.
        sv22o_pt(i2,i1)=0.
        sv1_eta(i2,i1)=0.
        sv2_eta(i2,i1)=0.
        sv1o_eta(i2,i1)=0.
        sv2o_eta(i2,i1)=0.
        sv12_eta(i2,i1)=0.
        sv22_eta(i2,i1)=0.
        sv12o_eta(i2,i1)=0.
        sv22o_eta(i2,i1)=0.

        ssv1_pt(i2,i1)=0.
        ssv2_pt(i2,i1)=0.
        ssv1o_pt(i2,i1)=0.
        ssv2o_pt(i2,i1)=0.
        ssv12_pt(i2,i1)=0.
        ssv22_pt(i2,i1)=0.
        ssv12o_pt(i2,i1)=0.
        ssv22o_pt(i2,i1)=0.
        ssv1_eta(i2,i1)=0.
        ssv2_eta(i2,i1)=0.
        ssv1o_eta(i2,i1)=0.
        ssv2o_eta(i2,i1)=0.
        ssv12_eta(i2,i1)=0.
        ssv22_eta(i2,i1)=0.
        ssv12o_eta(i2,i1)=0.
        ssv22o_eta(i2,i1)=0.
        enddo
        enddo
c220312
c010412
        do i2=1,21
        do i1=1,20
        fssnum_pt(i2,i1)=0.   
        fssnum_eta(i2,i1)=0.   
        fsv1_pt(i2,i1)=0.
        fsv2_pt(i2,i1)=0.
        fsv1o_pt(i2,i1)=0.
        fsv2o_pt(i2,i1)=0.
        fsv12_pt(i2,i1)=0.
        fsv22_pt(i2,i1)=0.
        fsv12o_pt(i2,i1)=0.
        fsv22o_pt(i2,i1)=0.
        fsv1_eta(i2,i1)=0.
        fsv2_eta(i2,i1)=0.
        fsv1o_eta(i2,i1)=0.
        fsv2o_eta(i2,i1)=0.
        fsv12_eta(i2,i1)=0.
        fsv22_eta(i2,i1)=0.
        fsv12o_eta(i2,i1)=0.
        fsv22o_eta(i2,i1)=0.

        fssv1_pt(i2,i1)=0.
        fssv2_pt(i2,i1)=0.
        fssv1o_pt(i2,i1)=0.
        fssv2o_pt(i2,i1)=0.
        fssv12_pt(i2,i1)=0.
        fssv22_pt(i2,i1)=0.
        fssv12o_pt(i2,i1)=0.
        fssv22o_pt(i2,i1)=0.
        fssv1_eta(i2,i1)=0.
        fssv2_eta(i2,i1)=0.
        fssv1o_eta(i2,i1)=0.
        fssv2o_eta(i2,i1)=0.
        fssv12_eta(i2,i1)=0.
        fssv22_eta(i2,i1)=0.
        fssv12o_eta(i2,i1)=0.
        fssv22o_eta(i2,i1)=0.
        enddo
        enddo
c010412
	i=0   ! 061007
	nncoll=0   ! 061007
300	i=i+1   ! 061007
c220312
        do i1=1,21
        do i2=1,20
        snum_pt(i1,i2)=0.   
        snum_eta(i1,i2)=0.   
        v1_pt(i1,i2)=0.
        v2_pt(i1,i2)=0.
        v12_pt(i1,i2)=0.
        v22_pt(i1,i2)=0.
        v1_eta(i1,i2)=0.
        v2_eta(i1,i2)=0.
        v12_eta(i1,i2)=0.
        v22_eta(i1,i2)=0.
        enddo
        enddo
c220312
c010412
        do i1=1,21
        do i2=1,20
        fsnum_pt(i1,i2)=0.
        fsnum_eta(i1,i2)=0.
        fv1_pt(i1,i2)=0.
        fv2_pt(i1,i2)=0.
        fv12_pt(i1,i2)=0.
        fv22_pt(i1,i2)=0.
        fv1_eta(i1,i2)=0.
        fv2_eta(i1,i2)=0.
        fv12_eta(i1,i2)=0.
        fv22_eta(i1,i2)=0.
        enddo
        enddo
c010412
c160110
        do i1=1,9
        nreaco(i1)=nreac(i1)
        enddo
c160110
c	decide fragmentation scheme
	mstj(1)=mstj1_1    ! 151107
	mstj(2)=mstj2
	mstj(3)=mstj3
c060607 decay or not   
        mstj(21)=0   ! no decay 060607
	siijk=0.   ! 061007
c211006
        do i1=1,20
        bn(i1)=0.
	bnf(i1)=0.
        do i2=1,5
        do i3=1,20
        an(i1,i2,i3)=0.
	anf(i1,i2,i3)=0.
        enddo
        enddo
        enddo
c151107
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
        ifcom(i1)=0   ! 160110
        enddo
c151107
c141208
c        if(itden.eq.0 .or. itden.eq.2)then
c        non6_c=123456
c        ithroq=0
c        ithrob=0
c        ich=0
c        do i=1,4
c        throe(i)=0.
c        enddo
c        endif
c141208
c060607 call pyevnt
        if(itden.eq.0)call pyevnt   ! for pp, pbarp, ppbar; 060607 110607
        if(itden.eq.2)call pyeevt(0,ee)   ! for e+e-; 060607 110607
	call pyedit(2)
c060607
c131108        do j=1,n
c        kf=k(j,2)
c        if(kf.eq.22)then
c        k(j,2)=44   ! '44': prompt direct photon
c        endif
c131108        enddo
c080607
        if(adj140.ne.4)then   !!
c       remove hadrons from 'pyjets' to 'sbh' and truncate 'pyjets'
c        correspondingly
        call remo
c160110
c       find number of strings and line number of first and last components
c        of each string
        nstr=0
        jb=0
10000   do i1=jb+1,n
        if(k(i1,1).eq.2)then   ! i1 is 'A'
        do i2=i1+1,n
        if(k(i2,1).eq.1)then   ! i2 is 'V'
        nstr=nstr+1
        nstra(nstr)=i1   ! line number of first component of nstr-th string
        nstrv(nstr)=i2   ! line number of first component of nstr-th string
        jb=i2
        if(jb.lt.n)goto 10000
        if(jb.eq.n)goto 20000
        endif
        enddo
        endif
        enddo
20000   nstr00=nstr
        call pyedit(2)
c        if(i.eq.1409 .or. i.eq.1531 .or. i.eq.2771)then
cs      write(9,*)'af. remo iii,n,nbh=',iii,n,nbh   ! sa
c        write(22,*)'af. remo i,n,nbh=',i,n,nbh   ! sa
c        write(22,*)'nstr00,nstr,nstra,nstrv=',nstr00,nstr
c        write(22,*)(nstra(i1),i1=1,nstr)
c        write(22,*)(nstrv(i1),i1=1,nstr)
c        call pylist(1)
c        call prt_sbh(nbh)
c        endif
c160110
cc	if(nout.eq.1 .or. mod(i,nout).eq.0)then
cc	write(22,*)'after partonic initialization'   ! sa
cc	call pylist(1)   ! sa
cc	endif
	if(n.le.0)then   ! no parton at all
c       write(9,*)'nn,nncoll,n=',nn,nncoll,n   ! sa
        nncoll=nncoll+1
        if(nncoll.gt.nn)then
        write(9,*)'nncoll=',nncoll   ! sa
        stop 8888
        endif
        i=i-1
        goto 300
        endif
c151107
        if((adj140.eq.3.or.adj140.eq.5).and.adj112.eq.0.and.n.ge.1)then
c	'pyjets' to 'sbe'. etc.
	do i1=1,n
	i3=i1   ! 050805
	kf=k(i1,2)
        kfab=iabs(kf)
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then   ! 2
c     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then   ! 2
        idi=idi+1
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
        endif
c151107
c       break up diquark and give four momentum and four position 
c        to broken quarks (working in 'pyjets')
        call break
        n00=n   ! 160110
c160110 n00: 'largest line number' in 'pyjets'
c160110 partons above n00 appear after inelastic collision
        call pyedit(2)
c160110
c        if(i.eq.1409 .or. i.eq.1531 .or. i.eq.2771)then
c        write(22,*)'af. break n,nbh,nbe,n00=',n,nbh,nbe,n00   ! sa
c        call pylist(1)   ! sa
c        call prt_sbh(nbh)
c        endif
c160110
        if(adj140.eq.1)goto 304   ! run terminates after parton initiation 
c061007
c       'pyjets' to 'parlist'
        iprl=n
        do i1=1,n
        idp(i1)=k(i1,2)
        rp(4,i1)=0.
        eee=p(i1,4)
        pp(4,i1)=eee
        do j=1,3
        rp(j,i1)=v(i1,j)
        ppp=p(i1,j)
        pp(j,i1)=ppp
        vp(j,i1)=ppp/eee
        enddo
        rmp(i1)=p(i1,5)
        taup(i1)=0.
	vip(i1)=0.
	xap(i1)=0.
        enddo
        do i1=n+1,mplis
        do j=1,3
        rp(j,i1)=0.
        pp(j,i1)=0.
        vp(j,i1)=0.
        enddo
        rp(4,i1)=0.
        pp(4,i1)=0.
        taup(i1)=0.
	vip(i1)=0.
	xap(i1)=0.
        idp(i1)=0
        rmp(i1)=0.
        enddo
        if(adj140.eq.3)goto 305
	if(adj140.eq.2 .or. adj140.eq.5)then   !!! parton scattering 151107 
c       parton cascade process is assumed to be start at time 0.
	time_par=0.
	if(n.eq.1)goto 302   ! one parton only
c201203
c2	write(9,*)'before parcas'   ! sa
	adj1(28)=200.*1.   ! 1. is the radius of proton	
        iijk=0   ! 151203
        call parcas(time_par,jjj,iijk,ee,nap,1.,1.,n00)   ! 120603 160110
c160110 jjj: loop number in parton-parton collisions in a hh collision
c	component 1. and 1. above: radius of projectile and target, respectively
c       write(9,*)'after parcas i,iijk,time_par=',i,iijk,time_par! sa
c120603
        if(iijk.eq.1)then
c       write(9,*)'i,iijk=',i,iijk
        i=i-1   ! 151107
        goto 300   ! give up current event avoiding infinite collision loop
        endif
c120603
        if(iijk.eq.2)siijk=siijk+1 ! iijk=2: initial collis. list is empty
c051108 if(siijk.ne.0.)write(9,*)'empty initial collis,siijk=',siijk ! 151107
c051108
        if(siijk.ne.0.)then
c       write(9,*)'i,siijk=',i,siijk
        i=i-1
        goto 300   ! give up the event with empty initial parton collis. list
        endif
c051108 
c201203
c151107 if(iprl.eq.0)then   ! no parton at all, give up this event
c	i=i-1
c       goto 300  
c151107	endif
c       'parlist' to 'pyjets'
        n=iprl
        do i1=1,n
c160110 k(i1,1)=1
        k(i1,2)=idp(i1)
c160110 do j=3,5
c       k(i1,j)=0
c160110 enddo
        do j=1,4
        v(i1,j)=rp(j,i1)
        p(i1,j)=pp(j,i1)
        enddo
        p(i1,5)=rmp(i1)
        v(i1,5)=0.
        enddo
c160110 if(i.eq.1409 .or. i.eq.1531 .or. i.eq.2771)then
c        write(22,*)'af. parcas,n,nbh,nbe,n00=',n,nbh,nbe,n00
c        call pyedit(2)
c        call pylist(1)
c        call prt_sbh(nbh)
c160110 endif
c201203
        endif   !!! parton rescattering finished 151107
c151107
        if(adj140.eq.2)goto 304   ! run terminats after parton rescattering
305     continue   ! 030512
c030512	if(adj112.ne.0)goto 302   ! coalescence
	if(adj112.ne.0.or.(adj112.eq.0.and.(nreac(4).gt.nreaco(4).or.
     c   nreac(6).gt.nreaco(6))))goto 302   ! coalescence
c	recover parton configuration in 'sbe' 
        if(idi.gt.0)then   ! 160110
c160110
c        if(i.eq.1409 .or. i.eq.1531 .or. i.eq.2771)
c     c  write(22,*)'be. recove idi=',idi
c       if(iii.eq.12.or.iii.eq.29.or.iii.eq.151.or.iii.eq.213)then
c       write(22,*)'be. recover n,nbe,n00,idi=',n,nbe,n00,idi
c       call prt_sbe(nbe)
c       endif
c160110
c	loop over 'sbe'
	idii=0
cs	write(9,*)'be. recover n,nbe=',n,nbe   ! sa
	do ii=1,nbe
	kf=kbe(ii,2)
        kfab=iabs(kf)
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then
c060805     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then
	idii=idii+1
	do j=1,5
	kdiq(idii,j)=kbe(ii,j)   ! diquark k array
	enddo
	dgmas(idii)=pbe(ii,5)   ! diquark mass
c	write(9,*)'ii,kfab,idii,kdiq=',ii,kfab,idii,(kdiq(idii,j),j=1,5)
	endif
	enddo
c	loop over 'pyjets'
c	write(9,*)'pyjets n=',n
	idij=0
	jb=0
	dele=0.
880	do 980 ii=jb+1,n
	jb=jb+1
	ndiqi=ndiq(ii)
	if(ndiqi.ne.0)then   ! diquark (anti-diquark)
	idij=idij+1
	j=npt(ndiqi)   ! diquark partner 
	do i1=1,5
        k(ii,i1)=kdiq(idij,i1)   ! diquark k array
	enddo
c	write(9,*)'ii,p=',ii,(p(ii,i1),i1=1,5)
c	write(9,*)'j,p=',j,(p(j,i1),i1=1,5)
	do i1=1,3   
	p(ii,i1)=p(ii,i1)+p(j,i1)
	enddo
	dimass=dgmas(idij)   ! diquark mass
	pi1=p(ii,1)
	pi2=p(ii,2)
	pi3=p(ii,3)
	pi4=dsqrt(pi1*pi1+pi2*pi2+pi3*pi3+dimass*dimass)
	dele=dele+p(ii,4)+p(j,4)-pi4
c	write(9,*)'dele,p(ii,4),p(j,4),pi4=',dele,p(ii,4),p(j,4),pi4
	p(ii,4)=pi4
	p(ii,5)=dimass
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
c       move particle list, 'pyjets' and 'ndiq', one step downward from 
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
cs	call prt_luj(n)   ! sa
cs	call prt_sbh(nbh)

c       share energy 'dele' into particles
c	write(9,*)'n,dele=',n,dele
	del=dele/dfloat(n)
	do j3=1,n
        p(j3,4)=p(j3,4)+del
        if(del.lt.0.)then
        if(p(j3,4).lt.0.)p(j3,4)=p(j3,4)-del
        pabs=dabs(p(j3,3))
        if(pabs.ge.p(j3,4))p(j3,4)=p(j3,4)-del
        endif
        enddo
cs	write(22,*)'af share jb,n=',jb,n
cs	call prt_luj(n)
cs	call prt_sbh(nbh)

        if(iparres.eq.0 .or. (iparres.eq.1.and.(nreaco(4).eq.nreac(4))
     c   .and.(nreaco(6).eq.nreac(6)).and.(nreaco(7).eq.nreac(7))))
     c   then   !    ela. parton-parton collisions only 240412
	do i1=1,n
	do j1=1,5
	pbe(i1,j1)=p(i1,j1)
	enddo
	enddo
c       'sbe' to 'pyjets'
	call tran_sbe
	endif   ! 240412
        n00=n   ! 160110
c        if(i.eq.1409 .or. i.eq.1531 .or. i.eq.2771)then
c       write(22,*)'af. recover n,nbh,nbe,n00,idi=',n,nbh,nbe,n00,idi   ! sa
cs      write(9,*)'af. recover n,nbh=',n,nbh   ! sa
c       call pyedit(2)
c       call pylist(1)
c        call prt_sbh(nbh)
c       endif
        endif   ! 160110 
c151107
c151107	endif   !!!
c061007
302	if(adj140.eq.3 .or. adj140.eq.5)then   ! 151107
c       hadronization
c160110
        if(adj112.eq.0)then   ! 1
        if(nreac(4).gt.nreaco(4) .or. nreac(6).gt.nreaco(6))then
c       for inela. processes 4 and 6
c       write(9,*)'i,nthro,4,4o,6,6o=',
c     c   i,nthro,nreac(4),nreaco(4),nreac(6),nreaco(6)
        call coales(i,nn,nout,nap,nat,nzp,nzt)
        else
c       otherwise
        call sfm   ! string fragmentation 151107
        endif
        endif   ! 1
c160110
        if(adj112.ne.0)call coales(i,nn,nout,nap,nat,nzp,nzt)   ! coalescence
        endif   ! 151107
304     if(adj140.eq.1 .or. adj140.eq.2)then   ! 151107
c	'sbh' to 'pyjets'
        if(nbh.eq.0)goto 303   ! 261103
        do l=1,nbh
        l1=n+l
        do m=1,5
        k(l1,m)=kbh(l,m)
        p(l1,m)=pbh(l,m)
        v(l1,m)=vbh(l,m)
        enddo
        enddo
        n=n+nbh
303	do i1=n+1,kszj   ! 261103
        do j=1,5
        k(i1,j)=0
        k(i1,j)=0
        p(i1,j)=0.
        v(i1,j)=0.
        enddo
        enddo
        endif   ! 151107
        endif   !!
c080607
        if(adj140.eq.4)then   ! 080607
c131108
        nbh=0
        do i1=1,kszj
        do j1=1,5
        kbh(i1,j1)=0
        pbh(i1,j1)=0.
        vbh(i1,j1)=0.
        enddo
        enddo
        call sfm   ! string fragmentation
c131108
        endif   ! 080607
c171108
	if(kjp21.eq.1)then   ! 1
c	write(22,*)'be hadcas event=',i   !sa
c	call pylist(1)
c        write(22,*)'throe_t=',throe_t   ! 161007
c        write(22,*)'ithroq_t,ithrob_t,ich_t=',ithroq_t,ithrob_t,
c     c   dfloat(ich_t)/3.   ! 161007
c	write(22,*)'throe_p=',throe_p
c	write(22,*)'ithroq_p,ithrob_p,ich_p=',ithroq_p,ithrob_p,
c     c  dfloat(ich_p)/3.
c        write(22,*)'throe=',throe   ! 161007
c        write(22,*)'ithroq,ithrob,ich=',ithroq,ithrob,
c     c   dfloat(ich)/3.   ! 161007
c       change K0S, K0L to K0, K0ba
        do j=1,n
        kf=k(j,2)
        if(kf.eq.130 .or. kf.eq.310)then
        rrlu=pyr(1)
        k(j,2)=311
        if(rrlu.gt.0.5)k(j,2)=-311
        endif
        enddo
c       'pyjets' to 'sa1_h'
        nn4=n
        do i1=1,n
        do i2=1,5
        kn(i1,i2)=k(i1,i2)
        pn(i1,i2)=p(i1,i2)
        rn(i1,i2)=v(i1,i2)
        enddo
        enddo
	do ii=nn4+1,kszj
	do j=1,5
	kn(ii,j)=0
	pn(ii,j)=0.
	rn(ii,j)=0.
	enddo
	enddo
        do ii=1,100
        numb(ii)=0
        enddo
        if(i.eq.1)call sysini_h
        call filt_h
        do ii=1,kfmax
        nup=numb(ii)
        enddo
        nbh1=nn4-nup
c        write(9,*)'i,nn4,nup,nbh1,kfmax=',i,nn4,nup,nbh1,kfmax
c        write(9,*)(numb(ii),ii=1,10)
c        write(9,*)(numb(ii),ii=11,20)
c        write(9,*)(numb(ii),ii=21,30)
c        write(9,*)(numb(ii),ii=31,40)
c        write(9,*)(numb(ii),ii=41,50)
c        write(9,*)(numb(ii),ii=51,52)
c       write(9,*)(numb(ii),ii=61,100)
c       nup is the number of particles kept in 'sa1_h'
c       nbh1 is the number of particles storing in 'sbh'
c       'sa1_h' to 'sbh'
        if(nbh1.eq.0)goto 7000
        do ii=nup+1,nn4
        nbh=ii-nup
        do j=1,5
        kbh(nbh,j)=kn(ii,j)
        pbh(nbh,j)=pn(ii,j)
        vbh(nbh,j)=rn(ii,j)
        enddo
        enddo
7000    continue
	nn4=nup
	nbh=nbh1   ! 261103
	do ii=nn4+1,kszj
	do j=1,5
	kn(ii,j)=0
	pn(ii,j)=0.
	rn(ii,j)=0.
	enddo
	enddo
c241103
c	hadronic cascade (rescattering)
c	write(9,*)'be hadcas i,n,nbh=',i,n,nbh   ! sa
c	call prt_sa1_h(nn4)
c	call prt_sbh(nbh)
        time_had=0.   ! 270910
	call hadcas(i,nn,nout,time_had,ijkk)   ! 241103
c	write(9,*)'after hadcas i,ijkk,time_had,nn4,nbh=',
c     c	 i,ijkk,time_had,nn4,nbh   ! sa
	if(ijkk.eq.1)then   ! 161203   
c110603	i=i-1   ! it has been executed in 'scat' in 'hadcas' 
	goto 300   ! give up current event avoiding infinite collision loop
	endif
c241103
c	write(22,*)'af hadcas i=',i   !sa
c	write(9,*)'af hadcas i=',i   ! sa
c	call prt_sa1_h(nn4)
c	call prt_sbh(nbh)
c201203
c	'sa1_h' to 'pyjets'
        n=nn4
        do i1=1,n
        do i2=1,5
        k(i1,i2)=kn(i1,i2)
        p(i1,i2)=pn(i1,i2)
        v(i1,i2)=rn(i1,i2)
        enddo
        enddo
c201203
c       'sbh' to 'pyjets'
	if(nbh.eq.0)goto 9000   ! 261103
        do l=1,nbh
        l1=n+l
        do m=1,5
        k(l1,m)=kbh(l,m)
        p(l1,m)=pbh(l,m)
        v(l1,m)=vbh(l,m)
        enddo
        enddo
        n=n+nbh
9000	continue
	do ii=n+1,kszj   ! 261103
	do j=1,5
	k(ii,j)=0
	p(ii,j)=0.
	v(ii,j)=0.
	enddo
	enddo
	endif   ! 1 241103
c       change K0,K0ba to K0L and K0S
	do j=1,n
	kf=k(j,2)
        if(kf.eq.311 .or. kf.eq.-311)then
        rrlu=pyr(1)
        k(j,2)=130
        if(rrlu.gt.0.5)k(j,2)=310
        endif
	enddo
c	if(i.eq.7)then
c	write(mstu(11),*)'be pyexec'
c	call pyedit(2)
c	call pylist(1)
c	write(22,*)'throe_t=',throe_t
c	write(22,*)'ithroq_t,ithrob_t,ich_t=',ithroq_t,ithrob_t,
c     c	 dfloat(ich_t)/3.
c	write(22,*)'throe_p=',throe_p
c	write(22,*)'ithroq_p,ithrob_p,ich_p=',ithroq_p,ithrob_p,
c     c	 dfloat(ich_p)/3.
c	write(22,*)'throe=',throe
c	write(22,*)'ithroq,ithrob,ithroc=',ithroq,ithrob,
c     c	 dfloat(ithroc)/3.
c	call prt_sbh(nbh)
c	endif  
c       perform particle, declared unstable in the 'mdcy' array, decay
c130205	call pyexec
c171108	rrp=1.16   ! 130205 
c171108	call decayh(rrp)   ! 130205 
c181003
	continue
c	write(9,*)'paciae, nnstop,zstop=',nnstop,zstop
        if(nnstop.ne.0)then
        sstop1=zstop/dfloat(nnstop)
        sstop=sstop+sstop1
        else
        nzstop=nzstop+1
        endif
c181003
c171108
c060607
        if((mod(i,nout).eq.0) .or. i.eq.nn)then   ! 160110
        write(22,*)'event=',i
	call pylist(1)
        write(22,*)'throe_t=',throe_t   ! 161007
        write(22,*)'ithroq_t,ithrob_t,ich_t=',ithroq_t,ithrob_t,
     c   dfloat(ich_t)/3.   ! 161007
	write(22,*)'throe_p=',throe_p
	write(22,*)'ithroq_p,ithrob_p,ich_p=',ithroq_p,ithrob_p,
     c   dfloat(ich_p)/3.
        write(22,*)'throe=',throe   ! 161007
        write(22,*)'ithroq,ithrob,ich=',ithroq,ithrob,
     c   dfloat(ich)/3.   ! 161007
c061007
	endif
c       analyse an event
c131108
        edi=0.   ! 131108
        pdi=0.   ! 131108
        sdi=0.   ! 131108
        call copl
c       write(9,*)'coor=',coor
c131108
        coox=coor(1)
        cooy=coor(2)
        cooz=coor(3)
c       write(88,*)'tota levents,n=',nn
c       write(88,*)'event',i
	do 400 j=1,n   
	w=1.
	ik=k(j,2)
        ikab=iabs(ik)   ! 091007
	plu6=pyp(j,6)
        abplu6=dabs(plu6)   ! 131108
        p1=p(j,1)
        p2=p(j,2)
        p3=p(j,3)
        p4=p(j,4)
        p5=p(j,5)   
        ppt=pyp(j,10)
        yy=pyp(j,17)
        eta=pyp(j,19)
        c(1)=yy
	if(ifram.eq.1)c(1)=eta
        c(2)=ppt
c       .
c       .
c       .
c131108
        iadj=0
        r1=v(j,1)
        r2=v(j,2)
        r3=v(j,3)
        r1=r1-coox
        r2=r2-cooy
        r3=r3-cooz
        v(j,1)=r1
        v(j,2)=r2
        v(j,3)=r3
        if((adj140.eq.1.or.adj140.eq.2) .and. ikab.le.6 .and.
     c   abplu6.gt.0.2d0)then   ! for parton
        call adjins(r1,r2,r3,p1,p2,p3,iadj,sdi)
        if(iadj.eq.1)then
        edi=edi+p4
        pdi=pdi+1.
        endif
        endif
c       charged parton energy
        if((adj140.eq.3.or.adj140.eq.4.or.adj140.eq.5) .and.
     c   abplu6.gt.0.9d0)then   ! for hadron
        call adjins(r1,r2,r3,p1,p2,p3,iadj,sdi)
        if(iadj.eq.1)then
        edi=edi+p4
        pdi=pdi+1.
        endif
        endif
c       charged hadron energy
c131108
        kkk=1
c       statistics negative multiplicity
        if((adj140.eq.3.or.adj140.eq.4.or.adj140.eq.5) .and. 
     c   plu6.lt.-0.9)then   ! for hadron 151107 
	nminf=nminf-plu6
        do i1=1,iflmax
        if(c(i1).lt.afl(kkk,i1,1) .or. c(i1).gt.afl(kkk,i1,2))goto 700
        enddo
        nmin=nmin-plu6
700	endif
c091007
        if((adj140.eq.1.or.adj140.eq.2) .and. ikab.le.6 .and. 
     c   plu6.lt.-0.2)then   ! for parton 151107
        nminf=nminf+1
        do i1=1,iflmax
        if(c(i1).lt.afl(kkk,i1,1) .or. c(i1).gt.afl(kkk,i1,2))goto 702
        enddo
        nmin=nmin+1
702     endif
c091007
c	statistics of positive multiplicity
        if((adj140.eq.3.or.adj140.eq.4.or.adj140.eq.5) .and. 
     c   plu6.gt.0.9)then   ! for hadron 151107 
	npluf=npluf+plu6
        do i1=1,iflmax
        if(c(i1).lt.afl(kkk,i1,1) .or. c(i1).gt.afl(kkk,i1,2))goto 701
        enddo
	nplu=nplu+plu6
701	endif
c091007
        if((adj140.eq.1.or.adj140.eq.2) .and. ikab.le.6 .and. 
     c   plu6.gt.0.2)then   ! for parton 151107
        npluf=npluf+1
        do i1=1,iflmax
        if(c(i1).lt.afl(kkk,i1,1) .or. c(i1).gt.afl(kkk,i1,2))goto 703
        enddo
        nplu=nplu+1
703     endif
c091007
c       statistics of y, pt, ect. distributions
c       if(i.eq.1)then
c       write(9,*)'be. 500 i,j,ik,ispmax=',i,j,ik,ispmax
c       write(9,*)'ispkf=',(ispkf(l),l=1,10)
c       write(9,*)(ispkf(l),l=11,ispmax)
c       write(9,*)'i,j,ik,w=',i,j,ik,w
c       write(9,*)'yy,ppt,eta,p5=',yy,ppt,eta,p5
c       write(9,*)'be. particle multiplicity=',(bn(ll),ll=1,ispmax)
c       endif
        do 500 kk=1,ispmax
        kf=ispkf(kk)
        if(kf.ne.ik)goto 500
c250510
	if(ik.eq.2212 .or. ik.eq.2112)then
	if(ppt.le.2.d-1)goto 500
	endif
c       exclude the projectile and target spectator nucleons 
c        (exclude pT<0.01 proton and neutron for pp) 
c250510
c080610
c	if((ik.eq.1 .or. ik.eq.2) .and. (p5.ge.(0.330d0-0.001d0) .and. 
c	c  p5.lt.(0.330d0+0.001d0)))goto 500
        if((ik.eq.1 .or. ik.eq.2) .and. p5.eq.0.330d0)goto 500
c       exclude the valence quark
c080610
c       if(i.eq.1)then
c       write(9,*)'i,j,ik,kk,w=',i,j,ik,kk,w
c       write(9,*)'yy,ppt,eta,p5=',yy,ppt,eta,p5
c       write(9,*)'be. particle multiplicity=',(bn(ll),ll=1,ispmax)
c       endif
        call stati(yy,ppt,eta,p5,ik,kk,w,bn,an,bnf,anf,p3)
        call flow_f(p1,p2,ppt,yy,eta,p5,ik,kk)   ! 220312 010412
c       write(9,*)'af. particle multiplicity=',(bn(ll),ll=1,ispmax)
        goto 400
500     continue
cccccccc--ayut-start
c       if(nout.eq.1)then
c       write(88,888)j,k(j,2),(p(j,i1),i1=1,4),(v(j,i1),i1=1,4)
c       endif
888     format(i5,1x,i5,4f7.3,1x,4f7.3)
cccccccc--ayut-stop
400     continue
c       write(9,*)'event=',i
c       write(9,*)'particle multiplicity=',(bn(ll),ll=1,ispmax)
c220312
        do i1=1,21
        do i2=1,20
        ssnum_pt(i1,i2)=ssnum_pt(i1,i2)+snum_pt(i1,i2)   
        ssnum_eta(i1,i2)=ssnum_eta(i1,i2)+snum_eta(i1,i2)   
        sv1_pt(i1,i2)=sv1_pt(i1,i2)+v1_pt(i1,i2)
        sv2_pt(i1,i2)=sv2_pt(i1,i2)+v2_pt(i1,i2)
        sv12_pt(i1,i2)=sv12_pt(i1,i2)+v12_pt(i1,i2)
        sv22_pt(i1,i2)=sv22_pt(i1,i2)+v22_pt(i1,i2)
        sv1_eta(i1,i2)=sv1_eta(i1,i2)+v1_eta(i1,i2)
        sv2_eta(i1,i2)=sv2_eta(i1,i2)+v2_eta(i1,i2)
        sv12_eta(i1,i2)=sv12_eta(i1,i2)+v12_eta(i1,i2)
        sv22_eta(i1,i2)=sv22_eta(i1,i2)+v22_eta(i1,i2)

	sspt=snum_pt(i1,i2)
        if(sspt.gt.0.)then
        ssv1_pt(i1,i2)=ssv1_pt(i1,i2)+v1_pt(i1,i2)/sspt
        ssv2_pt(i1,i2)=ssv2_pt(i1,i2)+v2_pt(i1,i2)/sspt
        ssv12_pt(i1,i2)=ssv12_pt(i1,i2)+v12_pt(i1,i2)/sspt
        ssv22_pt(i1,i2)=ssv22_pt(i1,i2)+v22_pt(i1,i2)/sspt
        endif
        sseta=snum_eta(i1,i2)
        if(sseta.gt.0.)then
        ssv1_eta(i1,i2)=ssv1_eta(i1,i2)+v1_eta(i1,i2)/sseta
        ssv2_eta(i1,i2)=ssv2_eta(i1,i2)+v2_eta(i1,i2)/sseta
        ssv12_eta(i1,i2)=ssv12_eta(i1,i2)+v12_eta(i1,i2)/sseta
        ssv22_eta(i1,i2)=ssv22_eta(i1,i2)+v22_eta(i1,i2)/sseta
        endif
        enddo
        enddo
c220312
c010412
        do i1=1,21
        do i2=1,20
        fssnum_pt(i1,i2)=fssnum_pt(i1,i2)+fsnum_pt(i1,i2)   
        fssnum_eta(i1,i2)=fssnum_eta(i1,i2)+fsnum_eta(i1,i2)   
        fsv1_pt(i1,i2)=fsv1_pt(i1,i2)+fv1_pt(i1,i2)
        fsv2_pt(i1,i2)=fsv2_pt(i1,i2)+fv2_pt(i1,i2)
        fsv12_pt(i1,i2)=fsv12_pt(i1,i2)+fv12_pt(i1,i2)
        fsv22_pt(i1,i2)=fsv22_pt(i1,i2)+fv22_pt(i1,i2)
        fsv1_eta(i1,i2)=fsv1_eta(i1,i2)+fv1_eta(i1,i2)
        fsv2_eta(i1,i2)=fsv2_eta(i1,i2)+fv2_eta(i1,i2)
        fsv12_eta(i1,i2)=fsv12_eta(i1,i2)+fv12_eta(i1,i2)
        fsv22_eta(i1,i2)=fsv22_eta(i1,i2)+fv22_eta(i1,i2)

	fsspt=fsnum_pt(i1,i2)
        if(fsspt.gt.0.)then
        fssv1_pt(i1,i2)=fssv1_pt(i1,i2)+fv1_pt(i1,i2)/fsspt
        fssv2_pt(i1,i2)=fssv2_pt(i1,i2)+fv2_pt(i1,i2)/fsspt
        fssv12_pt(i1,i2)=fssv12_pt(i1,i2)+fv12_pt(i1,i2)/fsspt
        fssv22_pt(i1,i2)=fssv22_pt(i1,i2)+fv22_pt(i1,i2)/fsspt
        endif
        fsseta=fsnum_eta(i1,i2)
        if(fsseta.gt.0.)then
        fssv1_eta(i1,i2)=fssv1_eta(i1,i2)+fv1_eta(i1,i2)/fsseta
        fssv2_eta(i1,i2)=fssv2_eta(i1,i2)+fv2_eta(i1,i2)/fsseta
        fssv12_eta(i1,i2)=fssv12_eta(i1,i2)+fv12_eta(i1,i2)/fsseta
        fssv22_eta(i1,i2)=fssv22_eta(i1,i2)+fv22_eta(i1,i2)/fsseta
        endif
        enddo
        enddo
c010412
        do kk=1,ispmax
        sbn(kk)=sbn(kk)+bn(kk)
	sbnf(kk)=sbnf(kk)+bnf(kk)
        do i1=1,20
        do i2=1,isdmax
        san(i1,i2,kk)=san(i1,i2,kk)+an(i1,i2,kk)
	sanf(i1,i2,kk)=sanf(i1,i2,kk)+anf(i1,i2,kk)
        enddo
        enddo
        enddo
        sedi=sedi+edi   ! 131108
        spdi=spdi+pdi   ! 131108
        ssdi=ssdi+sdi   ! 131108
        stime_par=stime_par+time_par   ! 171108
        stime_had=stime_had+time_had   ! 171108
c        if(time_par.ge.4.d-316)stimep=stimep+1.d0   ! 280910
c        if(time_had.ge.4.d-316)stimeh=stimeh+1.d0   ! 280910
        if(time_par.gt.0.d0)stimep=stimep+1.d0   ! 280910
        if(time_had.gt.0.d0)stimeh=stimeh+1.d0   ! 280910
c       write(9,*)'i,stimep,stimeh=',i,stimep,stimeh   ! 280910
c051108
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        skapa(1)=skapa(1)+akapa(1)
        skapa(2)=skapa(2)+akapa(2)
        skapa(3)=skapa(3)+akapa(3)
        skapa(4)=skapa(4)+akapa(4)
        skapa(5)=skapa(5)+akapa(5)
        sgtime=sgtime+gtime
        sastr=sastr+astr
        if(itime.ne.0)mkapa=mkapa+1
        endif
c       mkapa: number of events with string
c051108
c161007
        sthroq=sthroq+ithroq+ithroq_p+ithroq_t
        sthrob=sthrob+ithrob+ithrob_p+ithrob_t
        sthroc=sthroc+ich+ich_p+ich_t
        do i1=1,4
        sthroe(i1)=sthroe(i1)+throe(i1)+throe_p(i1)+throe_t(i1)
        enddo
c161007
	if(mod(i,nout).eq.0 .or. i.eq.nn)then   ! 211006
	open(8,file='rms.out',status='unknown')   ! 051108
	flaa=dfloat(i)   ! 241006
        cubv=dfloat(lcub)
        cubv=cubv*cubv*cubv
        sedio=sedi/cubv/flaa   ! 131108
        spdio=spdi/cubv/flaa   ! 131108
        ssdio=ssdi/cubv/flaa   ! 131108
c280910
        if(stimep.gt.0.d0)then
        stime_paro=stime_par/stimep   ! 171108
        else
        stime_paro=0.d0
        endif 
        if(stimeh.gt.0.d0)then
        stime_hado=stime_had/stimeh   ! 171108
        else
        stime_hado=0.d0
        endif 
c280910
	dnminf=nminf/flaa
	dnpluf=npluf/flaa
        dnmin=nmin/flaa
        dnplu=nplu/flaa
        do kk=1,ispmax
        sbo(kk)=sbn(kk)/flaa   !241006
        sbof(kk)=sbnf(kk)/flaa   !241006
        do i1=1,20
        do i2=1,isdmax
        sao(i1,i2,kk)=san(i1,i2,kk)/flaa   !241006
	saof(i1,i2,kk)=sanf(i1,i2,kk)/flaa   !241006
        enddo
        enddo
        enddo
c051108
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        sastro=sastr/flaa
        if(mkapa.gt.0)then
        sgtimeo=sgtime/dfloat(mkapa)
c       sgtimeo: average number of gluons in a string over strings with gluon
        do i1=1,5
        skapao(i1)=skapa(i1)/dfloat(mkapa)
        enddo
        endif
        endif
c051108
c220312
        do i1=1,21
        do i2=1,20
        sss=ssnum_pt(i1,i2)   ! 280607
        if(sss.gt.1.d-5)then   ! 280607
        sv1op=sv1_pt(i1,i2)/sss   ! direct, average
        sv1o_pt(i1,i2)=sv1op
        sv12op=sv12_pt(i1,i2)/sss   ! direct, squared average
c121207 if(sss-1.gt.1.d-5)then   ! 280607
c121207 sv12op=(sv12op-sv1op*sv1op)/(sss-1.)
        sv12op=sv12op-sv1op*sv1op   ! 121207
        if(sv12op.le.1.d-20)sv12op=1.d-20
        sv12o_pt(i1,i2)=sqrt(sv12op)   !  direct, square root of variance
c121207 endif   ! 280607
        sv2op=sv2_pt(i1,i2)/sss   ! elliptic, average
        sv2o_pt(i1,i2)=sv2op
        sv22op=sv22_pt(i1,i2)/sss   ! elliptic, squared average
c121207 if(sss-1.gt.1.d-5)then   ! 280607
c121207 sv22op=(sv22op-sv2op*sv2op)/(sss-1.)
        sv22op=sv22op-sv2op*sv2op   ! 121207
        if(sv22op.le.1.d-20)sv22op=1.d-20
        sv22o_pt(i1,i2)=sqrt(sv22op)   ! elliptic, square root of variance
c121207 endif   ! 280607
        endif   ! 280607
        sss=ssnum_eta(i1,i2)   ! 280607
        if(sss.gt.1.d-5)then   ! 280607
        sv1oe=sv1_eta(i1,i2)/sss   ! direct, average
        sv1o_eta(i1,i2)=sv1oe
        sv12oe=sv12_eta(i1,i2)/sss   ! direct, squared average
c121207 if(sss-1.gt.1.d-5)then   ! 280607
c121207 sv12oe=(sv12oe-sv1oe*sv1oe)/(sss-1.)
        sv12oe=sv12oe-sv1oe*sv1oe   ! 121207
        if(sv12oe.le.1.d-20)sv12oe=1.d-20
        sv12o_eta(i1,i2)=sqrt(sv12oe)   !  direct, square root of variance
c121207 endif   ! 280607
        sv2oe=sv2_eta(i1,i2)/sss   ! elliptic, average
        sv2o_eta(i1,i2)=sv2oe
        sv22oe=sv22_eta(i1,i2)/sss   ! elliptic, squared average
c121207 if(sss-1.gt.1.d-5)then   ! 280607
c121207 sv22oe=(sv22oe-sv2oe*sv2oe)/(sss-1.)
        sv22oe=sv22oe-sv2oe*sv2oe   ! 121207   
        if(sv22oe.le.1.d-20)sv22oe=1.d-20
        sv22o_eta(i1,i2)=sqrt(sv22oe)   ! elliptic, square root of variance
c121207 endif   ! 280607
        endif   ! 280607
        enddo
        enddo

        do i1=1,21
        do i2=1,20
        sv1op=ssv1_pt(i1,i2)/flaa
        ssv1o_pt(i1,i2)=sv1op   ! direct, average
        sv12op=ssv12_pt(i1,i2)/flaa   ! direct, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv12op=(sv12op-sv1op*sv1op)/(flaa-1.)
        sv12op=sv12op-sv1op*sv1op   ! 121207
        if(sv12op.le.1.d-20)sv12op=1.d-20
        ssv12o_pt(i1,i2)=sqrt(sv12op)   !  direct, square root of variance
c121207 endif
        sv2op=ssv2_pt(i1,i2)/flaa   
        ssv2o_pt(i1,i2)=sv2op   ! elliptic, average
        sv22op=ssv22_pt(i1,i2)/flaa   ! elliptic, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv22op=(sv22op-sv2op*sv2op)/(flaa-1.)
        sv22op=sv22op-sv2op*sv2op   !121207
        if(sv22op.le.1.d-20)sv22op=1.d-20
        ssv22o_pt(i1,i2)=sqrt(sv22op)   ! elliptic, square root of variance
c121207 endif
        sv1op=ssv1_eta(i1,i2)/flaa
        ssv1o_eta(i1,i2)=sv1op   ! direct, average
        sv12op=ssv12_eta(i1,i2)/flaa   ! direct, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv12op=(sv12op-sv1op*sv1op)/(flaa-1.)
        sv12op=sv12op-sv1op*sv1op   ! 121207
        if(sv12op.le.1.d-20)sv12op=1.d-20
        ssv12o_eta(i1,i2)=sqrt(sv12op)   !  direct, square root of variance
c121207 endif
        sv2op=ssv2_eta(i1,i2)/flaa   
        ssv2o_eta(i1,i2)=sv2op   ! elliptic, average
        sv22op=ssv22_eta(i1,i2)/flaa   ! elliptic, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv22op=(sv22op-sv2op*sv2op)/(flaa-1.)
        sv22op=sv22op-sv2op*sv2op   ! 121207
        if(sv22op.le.1.d-20)sv22op=1.d-20
        ssv22o_eta(i1,i2)=sqrt(sv22op)   ! elliptic, square root of variance
c121207 endif
        enddo
        enddo
c220312
c010412
        do i1=1,21
        do i2=1,20
        fsss=fssnum_pt(i1,i2)   ! 280607
        if(fsss.gt.1.d-5)then   ! 280607
        fsv1op=fsv1_pt(i1,i2)/fsss   ! direct, average
        fsv1o_pt(i1,i2)=fsv1op
        fsv12op=fsv12_pt(i1,i2)/fsss   ! direct, squared average
c121207 if(sss-1.gt.1.d-5)then   ! 280607
c121207 sv12op=(sv12op-sv1op*sv1op)/(sss-1.)
        fsv12op=fsv12op-fsv1op*fsv1op   ! 121207
        if(fsv12op.le.1.d-20)fsv12op=1.d-20
        fsv12o_pt(i1,i2)=sqrt(fsv12op)   !  direct, square root of variance
c121207 endif   ! 280607
        fsv2op=fsv2_pt(i1,i2)/fsss   ! elliptic, average
        fsv2o_pt(i1,i2)=fsv2op
        fsv22op=fsv22_pt(i1,i2)/fsss   ! elliptic, squared average
c121207 if(sss-1.gt.1.d-5)then   ! 280607
c121207 sv22op=(sv22op-sv2op*sv2op)/(sss-1.)
        fsv22op=fsv22op-fsv2op*fsv2op   ! 121207
        if(fsv22op.le.1.d-20)fsv22op=1.d-20
        fsv22o_pt(i1,i2)=sqrt(fsv22op)   ! elliptic, square root of variance
c121207 endif   ! 280607
        endif   ! 280607
        fsss=fssnum_eta(i1,i2)   ! 280607
        if(fsss.gt.1.d-5)then   ! 280607
        fsv1oe=fsv1_eta(i1,i2)/fsss   ! direct, average
        fsv1o_eta(i1,i2)=fsv1oe
        fsv12oe=fsv12_eta(i1,i2)/fsss   ! direct, squared average
c121207 if(sss-1.gt.1.d-5)then   ! 280607
c121207 sv12oe=(sv12oe-sv1oe*sv1oe)/(sss-1.)
        fsv12oe=fsv12oe-fsv1oe*fsv1oe   ! 121207
        if(fsv12oe.le.1.d-20)fsv12oe=1.d-20
        fsv12o_eta(i1,i2)=sqrt(fsv12oe)   !  direct, square root of variance
c121207 endif   ! 280607
        fsv2oe=fsv2_eta(i1,i2)/fsss   ! elliptic, average
        fsv2o_eta(i1,i2)=fsv2oe
        fsv22oe=fsv22_eta(i1,i2)/fsss   ! elliptic, squared average
c121207 if(sss-1.gt.1.d-5)then   ! 280607
c121207 sv22oe=(sv22oe-sv2oe*sv2oe)/(sss-1.)
        fsv22oe=fsv22oe-fsv2oe*fsv2oe   ! 121207   
        if(fsv22oe.le.1.d-20)fsv22oe=1.d-20
        fsv22o_eta(i1,i2)=sqrt(fsv22oe)   ! elliptic, square root of variance
c121207 endif   ! 280607
        endif   ! 280607
        enddo
        enddo

        do i1=1,21
        do i2=1,20
        fsv1op=fssv1_pt(i1,i2)/flaa
        fssv1o_pt(i1,i2)=fsv1op   ! direct, average
        fsv12op=fssv12_pt(i1,i2)/flaa   ! direct, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv12op=(sv12op-sv1op*sv1op)/(flaa-1.)
        fsv12op=fsv12op-fsv1op*fsv1op   ! 121207
        if(fsv12op.le.1.d-20)fsv12op=1.d-20
        fssv12o_pt(i1,i2)=sqrt(fsv12op)   !  direct, square root of variance
c121207 endif
        fsv2op=fssv2_pt(i1,i2)/flaa   
        fssv2o_pt(i1,i2)=fsv2op   ! elliptic, average
        fsv22op=fssv22_pt(i1,i2)/flaa   ! elliptic, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv22op=(sv22op-sv2op*sv2op)/(flaa-1.)
        fsv22op=fsv22op-fsv2op*fsv2op   !121207
        if(fsv22op.le.1.d-20)fsv22op=1.d-20
        fssv22o_pt(i1,i2)=sqrt(fsv22op)   ! elliptic, square root of variance
c121207 endif
        fsv1op=fssv1_eta(i1,i2)/flaa
        fssv1o_eta(i1,i2)=fsv1op   ! direct, average
        fsv12op=fssv12_eta(i1,i2)/flaa   ! direct, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv12op=(sv12op-sv1op*sv1op)/(flaa-1.)
        fsv12op=fsv12op-fsv1op*fsv1op   ! 121207
        if(fsv12op.le.1.d-20)fsv12op=1.d-20
        fssv12o_eta(i1,i2)=sqrt(fsv12op)   !  direct, square root of variance
c121207 endif
        fsv2op=fssv2_eta(i1,i2)/flaa   
        fssv2o_eta(i1,i2)=fsv2op   ! elliptic, average
        fsv22op=fssv22_eta(i1,i2)/flaa   ! elliptic, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv22op=(sv22op-sv2op*sv2op)/(flaa-1.)
        fsv22op=fsv22op-fsv2op*fsv2op   ! 121207
        if(fsv22op.le.1.d-20)fsv22op=1.d-20
        fssv22o_eta(i1,i2)=sqrt(fsv22op)   ! elliptic, square root of variance
c121207 endif
        enddo
        enddo
c010412
c161007
        wthroq=sthroq/flaa
        wthrob=sthrob/flaa
        wthroc=sthroc/flaa
        do i1=1,4
        wthroe(i1)=sthroe(i1)/flaa
        enddo
c161007
c051207
c160110 nrea=0
        srea=0.   ! 080410
        do i1=1,9
c080410 nrea=nrea+nreac(i1)
        snreac(i1)=nreac(i1)/flaa
        srea=srea+snreac(i1)   ! 080410
        enddo
c080410 srea=dfloat(nrea)/flaa
        write(8,*)'average collision # in parton cascade=',srea
        write(8,*)'total # of scaterring processes in parton cascade='
        write(8,600)(snreac(i1),i1=1,9)
600     format(9(1x,e10.3))
c051207
c051108
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        write(8,*)'kjp22=0, par1,par2,par3,par21=',
     c   parj1,parj2,parj3,parj21
        write(8,*)'keff2,par2,par21,par1,par3=',
     c   (skapao(i1),i1=1,5)
        write(8,*)'averaged # of gluon in a string,averaged # of' 
        write(8,*)' strings in an event and averaged # of events with' 
        write(8,*)'string=',sgtimeo,sastro,mkapa
        endif
c051108
        write(8,*)'multiplicity of negative particles=',dnmin,dnminf
        write(8,*)'multiplicity of positive particles=',dnplu,dnpluf
        write(8,*)'particle multiplicity,p=',(sbo(ll),ll=1,ispmax)   !241006
	write(8,*)'particle multiplicity,f=',(sbof(ll),ll=1,ispmax)   !241006
	write(8,*)'time_par,time_had,sedio,spdio,ssdio,i=',stime_paro,
     c   stime_hado,sedio,spdio,ssdio,i   ! 131108
c161007
        write(8,*)'wthroe=',wthroe   ! 161007
        write(8,*)'wthroq,wthrob,wthroc=',wthroq,wthrob,wthroc/3.
c161007
c220312
	write(8,*)'particle-wise average, partial phase space'   ! 010412
        write(8,*)'final state v1(pt)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(sv1o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v1(pt)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(sv12o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v2(pt)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(sv2o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v2(pt)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(sv22o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v1(eta)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(sv1o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v1(eta)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(sv12o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v2(eta)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(sv2o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v2(eta)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(sv22o_eta(i1,i2),i2=1,20)
        enddo
602     format(6(1x,e10.3))

	write(8,*)'event-wise average, partial phase space'   ! 010412
        write(8,*)'final state v1(pt) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(ssv1o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v1(pt) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(ssv12o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v2(pt) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(ssv2o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v2(pt) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(ssv22o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v1(eta) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(ssv1o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v1(eta) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(ssv12o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v2(eta) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(ssv2o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v2(eta) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(ssv22o_eta(i1,i2),i2=1,20)
        enddo
c220312
c010412
	write(8,*)'particle-wise average, full phase space'
        write(8,*)'final state v1(pt)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fsv1o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v1(pt)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fsv12o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v2(pt)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fsv2o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v2(pt)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fsv22o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v1(eta)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fsv1o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v1(eta)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fsv12o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v2(eta)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fsv2o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v2(eta)'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fsv22o_eta(i1,i2),i2=1,20)
        enddo

	write(8,*)'event-wise average, full phase space'
        write(8,*)'final state v1(pt) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fssv1o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v1(pt) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fssv12o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v2(pt) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fssv2o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v2(pt) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fssv22o_pt(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v1(eta) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fssv1o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v1(eta) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fssv12o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'final state v2(eta) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fssv2o_eta(i1,i2),i2=1,20)
        enddo
        write(8,*)'square root of variance of v2(eta) e'
        do i1=1,21
        write(8,*)'i1=',i1
        write(8,602)(fssv22o_eta(i1,i2),i2=1,20)
        enddo
c010412
        do m2=1,isdmax
        write(8,*)'ID of distribution m2=',m2
        do m3=1,ispmax
        write(8,*)'distribution belong to m3=',m3
        write(8,*)(sao(m1,m2,m3),m1=1,20)   !241006
	write(8,*)(saof(m1,m2,m3),m1=1,20)   !241006
        enddo
        enddo
	endif   ! 211006
c061007	300	continue
c051108
	close(8)
c       following four statements are moved from the beginning of loop over
c        event to the end
	if(mod(i,nout).eq.0)print*,'event=',i
	open(7,file='nout.out',status='unknown')
        write(7,*)'i=',i
        close(7)
        open(5,file='usux.dat',status='old')
        read(5,*)nn,nout
        close(5)
c051108
c061007
	if(i.eq.nn)goto 301
	if(i.lt.nn)goto 300   
c061007
c       statistics of processes generated
301	call pystat(0)   ! 061007
        close(9)
	close(22)
	stop
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine stati(y,pt,eta,p5,ik,kk,ww,a,b,af,bf,p3)
c       on line statistics
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2),ifram
	dimension a(20),b(20,5,20),af(20),bf(20,5,20),c(5),id(5)
        do 10000 i=1,iflmax
        goto (10,20,30,40,50) i
10      c(i)=y
        if(ifram.eq.1)c(i)=eta
        goto 10000
20      c(i)=pt
        goto 10000
30      continue
40      continue
50      continue
10000   continue
c       calculate the abscissa one by one
40000   do 20000 i=1,isdmax
        goto (100,200,300,400,500) i
c       y is located in which interval?
100     ii=dabs(y)/asd(i)+1
        if(ifram.eq.1 .and. y.gt.0.)ii=ii+10
        if(ifram.eq.1 .and. y.lt.0.)ii=10-ii+1
c       note: 10 here should be changed together with the dimension of
c       20
        id(i)=ii
c100    id(i)=y/asd(i)+1
c       if(id(i).le.0)id(i)=1
        goto 20000
c       pt is located in which interval?
200     id(i)=pt/asd(i)+1
	ptt=asd(i)*id(i)-asd(i)/2.   ! 150307
        goto 20000
c       eta is located in which interval?
300     ii=dabs(eta)/asd(i)+1
        if(ifram.eq.1 .and. eta.gt.0.)ii=ii+10
        if(ifram.eq.1 .and. eta.lt.0.)ii=10-ii+1
c       note: 10 here should be changed together with the dimension of
c        20
        id(i)=ii
        goto 20000
400     continue
        goto 20000
500     continue
20000   continue
c       make statistics of particle yield and desired distributions
        af(kk)=af(kk)+ww
        do i=1,isdmax
        ii=id(i)
        if(ii.lt.1 .or. ii.gt.20)goto 30000
        if(i.eq.1)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)
        if(i.eq.2)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)/ptt   ! 150307
        if(i.eq.3)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)
30000   enddo
c       put the filter to be effective
        do i=1,iflmax
        if(c(i).lt.afl(kk,i,1) .or. c(i).gt.afl(kk,i,2))return
        enddo
c       make statistics of particle yield and desired distributions
        a(kk)=a(kk)+ww
        do i=1,isdmax
        ii=id(i)
	if(ii.lt.1 .or. ii.gt.20)goto 50000
        if(i.eq.1)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
        if(i.eq.2)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)/ptt   ! 150307
        if(i.eq.3)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
50000   enddo
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
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 151107 160110
	jb=0
        ii=idio   ! 151107
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

c	according to composition of diquark 190107
        iy=isign(1,kf)
        kf1=kfab/1000
        kf2=(kfab-1000*kf1)/100
        kf1=iy*kf1
        kf2=iy*kf2

c151107 first component of diquark takes the line # of diquark 
c        and put second component at n+1
	k(i1,2)=kf1
        k(n+1,2)=kf2
        k(n+1,1)=1 
        k(n+1,3)=0
        k(n+1,4)=0
        k(n+1,5)=0
        ii=ii+1   ! 151107
        npt(ii)=n+1   ! 151107
        ifcom(ii)=i1   ! 160110
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

c        am1=pymass(kf1)
c        am2=pymass(kf2)
        pp(1,5)=am1
        pp(2,5)=am2
c       pp : four momentum of broken quark, local variable 
        do i1=1,4
        ps(i1)=p(ii,i1)
        enddo
c       ps : four momentum of diquark, local variable 
	goto 400   ! activate it for 'decay method'
c       broken quarks share out diquark four momentum randomly,
c	 denoted as 'random four momentum method'
c	do i1=1,4   ! activate it for 'random four momentum method'
c	broken quarks share out diquark three momentum randomly,
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
c	am1 and am2: mass of decayed particle pair
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=40000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4),pj(4),ps(4),pp(20,5)   
        double precision bb(3)
c       calculate the E and |p| of broken quark in rest frame of diquark
        sm2=ps(4)*ps(4)-ps(1)*ps(1)-ps(2)*ps(2)-ps(3)*ps(3)
c       one problem here is that 'sm2' may not equal to square of diquark 
c	 (or gluon) rest mass, in gluon splitting especially
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
        sm=dsqrt(sm2)   ! M (should be diquark mass)
c       pp(1,4)=(sm2-am2*am2+am1*am1)/2./sm
c       pp(2,4)=(sm2-am1*am1+am2*am2)/2./sm
        ppp=(sm2-(am1+am2)*(am1+am2))*(sm2-(am1-am2)*(am1-am2))
c161204	ppp=dabs(ppp)   ! 030603 ?
	if(ppp.lt.1.e-28)ppp=1.e-28   !161204
        ppp=dsqrt(ppp)/2./sm
c110211 goto 500   ! activate it for exponential cos(seta) distribution
c       the direction of broken quark is sampled isotropically in '4pi'
        coset=1.-2.*pyr(1)
        if(dabs(coset).gt.1.)then
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
        if(coseta.lt.1.e-28)coseta=1.e-28
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
c       rotate to the frame where diquark (or gluon), ps, is described
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
        if(pi4.lt.1.d-28)pi4=1.d-28
        pj4=pj(1)*pj(1)+pj(2)*pj(2)+pj(3)*pj(3)+am2*am2
        if(pj4.lt.1.e-28)pj4=1.e-28
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
        subroutine coord(ii)
c       give four position to broken quarks
c       first broken quark takes the four posotion of diquark
c       second broken quark is arranged around first ones within
c        0.5 fm randumly in each of three coordinates and has same
c        fourth position as diquark
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
        if(pyr(1).gt.0.5)v(n+1,i1)=v(ii,i1)-rr(i1)
        enddo
        v(n+1,4)=v(ii,4)
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
        subroutine remo
c       remove hadrons from 'pyjets' to 'sbh' 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=40000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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
        if(kfab.le.8 .or. kfab.eq.2101 .or. kfab.eq.3101
     c   .or. kfab.eq.3201 .or. kfab.eq.1103 .or. kfab.eq.2103
     c   .or. kfab.eq.2203 .or. kfab.eq.3103 .or. kfab.eq.3203
     c   .or. kfab.eq.3303 .or. kfab.eq.21)then
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rotate(cctas,sctas,cfis,sfis,pp3,pi,pj)
c	perform rotation
c       pi,pj: input,four momentum of colliding pair before scattering
c              output,four momentum of scattered pair after scattering and 
c              rotation
c       pp3: momentum modulus of pi or pj, both are equal in their cms,
c        after scattering
c       cctas,sctas,cfis,sfis: direction cosines of momentum of one
c        colliding particle after scattering, relative to the momentum
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



c********************************************************************
        subroutine tran_sbe
c       'sbe' to 'pyjets' 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (kszj=40000,KSZ1=30)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nsa,nonbe,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)   ! 080104
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine psum(pei,il,ih,peo)
c       calculate sum of momentum and energy
c       pei: two dimension array for momentum and energy of particles,inputed
c       il and ih: lower and upper limits of summation
c       peo: one dimension array for momentum and energy,outputed
c        implicit real*8 (a-h,o-z)
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



c131108cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine copl
c       calculate position of charged particles center of mass  
c        (parton or hadron) 
c       distance of a particle from center of mass is used to check that 
c        is it freezes out or not
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter(kszj=40000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
         COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/wz/coor(3)
        common/sa24/adj1(40),nnstop,non24,zstop
        adj140=adj1(40)
        do ii=1,3
        coor(ii)=0.
        enddo
        samass=0.
        if(adj140.eq.1 .or. adj140.eq.2)then   ! for parton
        do 110 ii=1,n
        kf=k(ii,2)
        plu6=pyp(ii,6)
        ikab=iabs(kf)
        abplu6=dabs(plu6)
        if(ikab.le.6 .and. abplu6.gt.0.2d0)then
        amass=pmas(pycomp(kf),1)
        samass=samass+amass
        do 100 jj=1,3
        coor(jj)=coor(jj)+amass*v(ii,jj)
100     continue
        endif
110     continue
        do ii=1,3
        coor(ii)=coor(ii)/dmax1(0.14d0,samass)
        enddo
        endif   ! for parton
        if(adj140.eq.3.or.adj140.eq.4.or.adj140.eq.5)then   ! for hadron
        do 210 ii=1,n
        kf=k(ii,2)
        plu6=pyp(ii,6)
        ikab=iabs(kf)
        abplu6=dabs(plu6)
        if(abplu6.gt.0.9d0)then
        amass=pmas(pycomp(kf),1)
        samass=samass+amass
        do 200 jj=1,3
        coor(jj)=coor(jj)+amass*v(ii,jj)
200     continue
        endif
210     continue
        do ii=1,3
        coor(ii)=coor(ii)/dmax1(0.14d0,samass)
        enddo
        endif   ! for hadron
        return
        end



c131108ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine adjins(r1,r2,r3,p1,p2,p3,iadj,sdi)
c       does the point of (r1,r2,r3) inside the sphere with radius 1
c        fm and centered at coor(3)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/sa29/effk1,lcub   ! 160110
        parameter(kszj=40000)
        common/wz/coor(3)
        cubl=dfloat(lcub)/2.d0
        x1=r1-coor(1)
        x2=r2-coor(2)
        x3=r3-coor(3)
        x12=x1*x1
        x22=x2*x2
        x32=x3*x3
        rrr=dsqrt(x12+x22+x32)
        if(rrr.ge.sdi)sdi=rrr
        ax1=dabs(x1)
        ax2=dabs(x2)
        ax3=dabs(x3)
c       ap1=dabs(p1)
c       ap2=dabs(p2)
c       ap3=dabs(p3)
c       if(ax1.le.0.5 .and. ax2.le.0.5 .and. ax3.le.0.5 .and.
c       c   ap1.le.0.25 .and. ap2.le.0.25 .and. ap3.le.0.25)iadj=1
        if(ax1.le.cubl .and. ax2.le.cubl .and. ax3.le.cubl)iadj=1
        return
        end



c171108cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbh(nn)
c       print particle list and sum of momentum and energy
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
        do i=1,nn
c       write(mstu(11),*)i,kbh(i,2),(pbh(i,j),j=1,4)
	write(9,*)i,kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
c       write(mstu(11),*)peo,ich1/3   !
	write(9,*)peo,ich1/3   !
        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbe(nn)   ! 160110
c       print particle list and sum of momentum and energy
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000)
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
        write(22,*)ich1/3,peo   !
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine flow_f(px,py,pt,y,eta,p5,ik,kk)   ! 220312 010412
c	calculate pt (eta) dependent direct and elliptic flow at final
c        state (partonic or hadronic)
	parameter (kszj=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c  afl(20,5,2),ifram
        common/sa18_pt/snum_pt(21,20),v1_pt(21,20),v2_pt(21,20),
     c   v12_pt(21,20),v22_pt(21,20)   ! 280607
	common/sa18_eta/snum_eta(21,20),v1_eta(21,20),v2_eta(21,20),
     c   v12_eta(21,20),v22_eta(21,20)   ! 280607
        common/sa24/adj1(40),nnstop,non24,zstop   ! 120707
c010412
        common/fsa18_pt/fsnum_pt(21,20),fv1_pt(21,20),fv2_pt(21,20),
     c   fv12_pt(21,20),fv22_pt(21,20)   
        common/fsa18_eta/fsnum_eta(21,20),fv1_eta(21,20),fv2_eta(21,20),
     c   fv12_eta(21,20),fv22_eta(21,20)  
        yoe=y
        if(ifram.eq.1)yoe=eta
c010412
c	snum(1,i1)=snum(1,i1)+1.
c	write(9,*)'px,py,pt,p5,ik,kk=',px,py,pt,p5,ik,kk   ! sa
	amass=p5   
        amass2=amass*amass
		if(pt.lt.1.e-20)pt=1.e-20 ! 080512, yan
        pt2=pt*pt
	iik=iabs(ik)
        adj140=adj1(40)

c       pt is located at which interval ?
	idpt=pt/asd(2)+1
c       eta is located at which interval ?
        ii=abs(eta)/asd(3)+1
        if(ifram.eq.1 .and. eta.gt.0.)ii=ii+10
        if(ifram.eq.1 .and. eta.lt.0.)ii=10-ii+1
c       note: 10 here should be change together with the dimension
c        20
	idmt=ii

	px2=px*px
        py2=py*py
	pxt=px/pt   ! v1
	pxt2=(px2-py2)/pt2   ! v2
c	write(9,*)'v1,v2=',pxt,pxt2   ! sa
        if(idpt.lt.1 .or. idpt.gt.20)goto 100
c010412	
c	statistics of v1(pt) (v2(pt)) in full eta phase space
        fsnum_pt(kk,idpt)=fsnum_pt(kk,idpt)+1   ! 280607
	fv1_pt(kk,idpt)=fv1_pt(kk,idpt)+pxt
	fv2_pt(kk,idpt)=fv2_pt(kk,idpt)+pxt2
	fv12_pt(kk,idpt)=fv12_pt(kk,idpt)+pxt*pxt
	fv22_pt(kk,idpt)=fv22_pt(kk,idpt)+pxt2*pxt2
c120607 direct and elliptic flows of charged particle in hadronic final state
        if(adj140.eq.3 .or. adj140.eq.4 .or. adj140.eq.5)then   ! 120707
	if(iik.eq.211 .or. iik.eq.321 .or. iik.eq.2212)then
        fsnum_pt(21,idpt)=fsnum_pt(21,idpt)+1   ! 280607
	fv1_pt(21,idpt)=fv1_pt(21,idpt)+pxt
	fv2_pt(21,idpt)=fv2_pt(21,idpt)+pxt2
	fv12_pt(21,idpt)=fv12_pt(21,idpt)+pxt*pxt
	fv22_pt(21,idpt)=fv22_pt(21,idpt)+pxt2*pxt2
        endif
        endif   ! 120707
c120707
        if(adj140.eq.1 .or. adj140.eq.2)then
        if(iik.eq.1 .or. iik.eq.2 .or. iik.eq.3)then
        fsnum_pt(21,idpt)=fsnum_pt(21,idpt)+1   ! 280607
        fv1_pt(21,idpt)=fv1_pt(21,idpt)+pxt
        fv2_pt(21,idpt)=fv2_pt(21,idpt)+pxt2
        fv12_pt(21,idpt)=fv12_pt(21,idpt)+pxt*pxt
        fv22_pt(21,idpt)=fv22_pt(21,idpt)+pxt2*pxt2
        endif
        endif
c	statistics of v1(pt) (v2(pt)) in partial eta phase space
	if(yoe.lt.afl(kk,1,1) .or. yoe.gt.afl(kk,1,2))goto 100   
c010412
        snum_pt(kk,idpt)=snum_pt(kk,idpt)+1   ! 280607
	v1_pt(kk,idpt)=v1_pt(kk,idpt)+pxt
	v2_pt(kk,idpt)=v2_pt(kk,idpt)+pxt2
	v12_pt(kk,idpt)=v12_pt(kk,idpt)+pxt*pxt
	v22_pt(kk,idpt)=v22_pt(kk,idpt)+pxt2*pxt2
c120607 direct and elliptic flows of charged particle in hadronic final state
        if(adj140.eq.3 .or. adj140.eq.4 .or. adj140.eq.5)then   ! 120707
	if(iik.eq.211 .or. iik.eq.321 .or. iik.eq.2212)then
        snum_pt(21,idpt)=snum_pt(21,idpt)+1   ! 280607
	v1_pt(21,idpt)=v1_pt(21,idpt)+pxt
	v2_pt(21,idpt)=v2_pt(21,idpt)+pxt2
	v12_pt(21,idpt)=v12_pt(21,idpt)+pxt*pxt
	v22_pt(21,idpt)=v22_pt(21,idpt)+pxt2*pxt2
        endif
        endif   ! 120707
c120707
        if(adj140.eq.1 .or. adj140.eq.2)then
        if(iik.eq.1 .or. iik.eq.2 .or. iik.eq.3)then
        snum_pt(21,idpt)=snum_pt(21,idpt)+1   ! 280607
        v1_pt(21,idpt)=v1_pt(21,idpt)+pxt
        v2_pt(21,idpt)=v2_pt(21,idpt)+pxt2
        v12_pt(21,idpt)=v12_pt(21,idpt)+pxt*pxt
        v22_pt(21,idpt)=v22_pt(21,idpt)+pxt2*pxt2
        endif
        endif
c120707
100     if(idmt.lt.1 .or. idmt.gt.20)goto 200
c010412
c	statistics of v1(eta) (v2(eta)) in full pt phase space
	
        fsnum_eta(kk,idmt)=fsnum_eta(kk,idmt)+1   ! 280607
	fv1_eta(kk,idmt)=fv1_eta(kk,idmt)+pxt
        fv2_eta(kk,idmt)=fv2_eta(kk,idmt)+pxt2
        fv12_eta(kk,idmt)=fv12_eta(kk,idmt)+pxt*pxt
        fv22_eta(kk,idmt)=fv22_eta(kk,idmt)+pxt2*pxt2
        if(adj140.eq.3 .or. adj140.eq.4 .or. adj140.eq.5)then   ! 120707
        if(iik.eq.211 .or. iik.eq.321 .or. iik.eq.2212)then
        fsnum_eta(21,idmt)=fsnum_eta(21,idmt)+1   ! 280607
	fv1_eta(21,idmt)=fv1_eta(21,idmt)+pxt
        fv2_eta(21,idmt)=fv2_eta(21,idmt)+pxt2
        fv12_eta(21,idmt)=fv12_eta(21,idmt)+pxt*pxt
        fv22_eta(21,idmt)=fv22_eta(21,idmt)+pxt2*pxt2
	endif
        endif   ! 120707
c120707
        if(adj140.eq.1 .or. adj140.eq.2)then
        if(iik.eq.1 .or. iik.eq.2 .or. iik.eq.3)then
        fsnum_eta(21,idmt)=fsnum_eta(21,idmt)+1   ! 280607
        fv1_eta(21,idmt)=fv1_eta(21,idmt)+pxt
        fv2_eta(21,idmt)=fv2_eta(21,idmt)+pxt2
        fv12_eta(21,idmt)=fv12_eta(21,idmt)+pxt*pxt
        fv22_eta(21,idmt)=fv22_eta(21,idmt)+pxt2*pxt2
        endif
        endif
c120707
c	statistics of v1(eta) (v2(eta)) in partial pt phase space
	if(pt.lt.afl(kk,2,1) .or. pt.gt.afl(kk,2,2))goto 200   
c010412
        snum_eta(kk,idmt)=snum_eta(kk,idmt)+1   ! 280607
	v1_eta(kk,idmt)=v1_eta(kk,idmt)+pxt
        v2_eta(kk,idmt)=v2_eta(kk,idmt)+pxt2
        v12_eta(kk,idmt)=v12_eta(kk,idmt)+pxt*pxt
        v22_eta(kk,idmt)=v22_eta(kk,idmt)+pxt2*pxt2
        if(adj140.eq.3 .or. adj140.eq.4 .or. adj140.eq.5)then   ! 120707
        if(iik.eq.211 .or. iik.eq.321 .or. iik.eq.2212)then
        snum_eta(21,idmt)=snum_eta(21,idmt)+1   ! 280607
	v1_eta(21,idmt)=v1_eta(21,idmt)+pxt
        v2_eta(21,idmt)=v2_eta(21,idmt)+pxt2
        v12_eta(21,idmt)=v12_eta(21,idmt)+pxt*pxt
        v22_eta(21,idmt)=v22_eta(21,idmt)+pxt2*pxt2
	endif
        endif   ! 120707
c120707
        if(adj140.eq.1 .or. adj140.eq.2)then
        if(iik.eq.1 .or. iik.eq.2 .or. iik.eq.3)then
        snum_eta(21,idmt)=snum_eta(21,idmt)+1   ! 280607
        v1_eta(21,idmt)=v1_eta(21,idmt)+pxt
        v2_eta(21,idmt)=v2_eta(21,idmt)+pxt2
        v12_eta(21,idmt)=v12_eta(21,idmt)+pxt*pxt
        v22_eta(21,idmt)=v22_eta(21,idmt)+pxt2*pxt2
        endif
        endif
c120707
200	return
	end



c171108************************************************************************
	BLOCK DATA PYCIDATA
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
	COMMON/PYCIDAT2/KFMAXT,nont2,PARAM(20),WEIGH(600)
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



C171108************************************************************
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