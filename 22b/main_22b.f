	program main_22b   ! 121110 071213
c	user main program of partron and hadron cascade model for relativistic 
c	 pA, AA (eA) collision 
c	composed of main_22b.f, parini_22b.f, parcas._22b.f, 
c        sfm_22b.f, coales_22b.f, hadcas_22b.f, and p22b.f
c       main_22b.f: an example of user main program 
c	parini_22.f: generate a partonic initial state for a
c        nucleus-nucleus collision  
c	parcas_22b.f: perform parton rescattering, where only 2->2 processes 
c	 are considered and LO pQCD cross section or its regularized 
c	 approximation is used
c	sfm_22b.f: hadronization according to LUND string fragmentation model 
c	coales_22b.f: hadronization according to Monte Carlo coalescence model
c	select sfm_22b.f or coales_22b.f by parameter adj1(12) 
c	hadcas_22b.f: perform hadronic rescattering
c       p22b.f: pythia 6.4 with a little bit modifications
c	read paciae_guide for the details
c       note: the statistics made here is for either parton or    
c        hadron according to purpose 
	parameter (kszj=40000,mplis=40000,KSZ1=30)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
        common/pyint1/mint(400),vint(400)
        COMMON/pycidat1/KFACOT(100),DISDET(100),ISINELT(600)
	common/pycidat2/kfmaxt,nont2,param(20),weigh(600)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
	common/sa6/kfmaxi,nwhole
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c  afl(20,5,2)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &  iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
	common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
        common/sa16/dtt,dni(10),dpi(10),edi(10),bmin,bmax
     &   ,bar(10),abar(10),barf(10),abarf(10)   ! 033101
     &	 ,emin(10),eminf(10),eplu(10),epluf(10)   ! 033101
	common/sa18/tdh,itnum,non18,cptl,cptu,cptl2,cptu2,snum(4,20),
     &	 v1(4,20),v2(4,20),v12(4,20),v22(4,20)
	common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph   
     c	 ,vnlep   ! 260314
	common/sa23/kpar,knn,kpp,knp,kep   ! 200601 060813
	common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003
	common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21,adiv,gpmax,nnc   ! 020708 070417   
	common/sa29/parp78,lcub   ! 150612 yan 070417
	common/sa30/vneump,vneumt   ! 191110
        common/sa31/rmax,bbb(200),TA1(200),TA2(200),TA1A2(200),
     c  part1(200),part2(200)   ! 020511
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 220312 240412 131212
        common/sa6_c/ithroq,ithrob,ithroc,non_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa6_t/ithroq_t,ithrob_t,ich_t,non6_t,throe_t(4)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)   ! 050603
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 240209
        common/sa1_h/nn,non_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)   ! 050603
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	common/show/vip(mplis),xap(mplis)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	common/count/isinel(600)
	common/ctllist/npctl,npinel(600),npctl0,npel   ! 061103
        common/ctllist_p/nreac(9),nrel   ! 071103
	common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/throqb/iprlth,nonqb,rpth(4,mplis),ppth(4,mplis),
     c  idpth(mplis),rmpth(mplis)   ! 070905
        dimension an(20,5,20),bn(20),san(20,5,20),sbn(20),
     c	 anf(20,5,20),bnf(20),sanf(20,5,20),sbnf(20)
	dimension sao(20,5,20),sbo(20),saof(20,5,20),sbof(20)
        dimension skapa(5),skapao(5),snreac(9)   ! 020708 280809
	dimension c(5),dinel(600),dineli(600),sthroe(4),wthroe(4)
	dimension bpp(20),kdiq(kszj,5),dgmas(kszj)
	dimension acoll(20),acollp(20),acollt(20)
	dimension sbp(20) 
	dimension numbth(3),pl(100,5)  ! 070905 260314
        real nmin,nminf,ncha,nchaf!020203
c260314 pl(ii,5): four momentum and mass of ii-th lepton
c	common block sa24: adjustable variables
c	adj1(i), i=
c       1: k factor used in parton cascade
c       2: parameter \alpha_s (as in program) in parton cascade
c       3: parameter tcut in program, to avoid divergence in calculating
c          parton-parton differential cross section in parton cascade
c       4: parameter idw, number of intervals in numerical integration in parton
cascade
c       5: =1 and 0, with and without nuclear shadowing, respectively
c       6: parameter a, i.e. parj(41) in Lund string fragmentation function
c          if adj1(12)=0
c          parameter a in FF string fragmentation function if adj1(12)=1
c       7: parameter b, i.e. parj(42) in Lund string fragmentation function
c       8: mstp(82) in PYTHIA 6.4
c       9: parp(81), (D=1.4 GeV/c), effective minimum transverse momentum for
c          multiple interactions with mstp(82)=1
c       10: parp(31),k factor in pythia
c       11: time accuracy used in hadron cascade (time accuracy used in
c           parton cascade is dddt)
c       12: model of hadronization: =0 string fragmentation; =1: coalescence
c       13: dimension of meson table considered if adj1(12)=1
c       14: dimension of baryon table considered if adj1(12)=1
c       15: string tension
c       16: number of loops in deexcitation of energetic quark in coalescence
c       17: the threshold energy in deexcitation of energetic quark in
coalescence
c       18: =0 and 1 without and with partonic Pauli blocking in parton cascade,
c           respectively
c       19: time accuracy used in parton cascade (dddt in program)
c       20: =0 exact pQCD parton-parton cross section
c           =1 limited and regularized parton-parton cross section (B. Zhang)
c           =2 the same as 0 but flat scattering angle distribution is assumed
c           =3 the same as 1 but flat scattering angle distribution is assumed
c       21: =0 and 1 without and with phase space adjudgment, respectively,
c        in coalescence
c       22: critical value of the product of radii both in coordinate and
c        momentum phase space (4 is assumed)
c       23: =0 LUND fragmentation function used in subroutine 'ffm' in
coalescence
c           =1 IF fragmentation function used
c       24: the virtuality cut ('tl0') in time-like radiation in parton cascade
c       25: \Lambda_QCD in parton cascade
c       26: number of random number thrown away
c       27: largest momentum allowed for particle ('dpmax' in program)
c       28: largest position allowed for particle (drmax=para10*max(rnt,rnp) 
c	    in program), giving 1 to it in usu.dat and calculating it in the 
c	    running)
c       29: width of two dimension Gaussian distribution sampling px and py of
c           produced quark pair in deexcitation of the energetic quark in
coalescence
c       30: maximum $p_T^2$ in above two dimension Gaussian distribution
c       31: parj(1) in pythia
c       32: parj(2) in pythia
c       33: parj(3) in pythia
c       34: parj(21) in pythia
c       35: mstp(91) in pythia,parton transverse momentum (k_{\perp})
c           distribution inside hadron;
c           =1, Gaussian;
c           =2, exponential
c       36: =0 without phenomenological parton energy loss in parton cascade
c           =1 with phenomenological parton energy loss
c       37: the coefficient ('c') in phenomenological parton energy loss
c       38: pt cut in phenomenological parton energy loss
c       39: width of Gaussian k_{\perp} distribution in hadron if mstp(91)=1
c           width of exponential k_{\perp} distribution in hadron if mstp(91)=2
c       40: =1 event endded after parton initiation
c           =2 event endded after parton rescattering
c           =4 event endded after hadron rescattering
c210803
c220312 smadel: small perpurbation of ellipse from circle
c220312	parecc: a parameter converting initial spatial space eccentricity 
c220312	 to final momentum space
c240412	iparres: =0 consider ela. parton-parton collisions only
c240412 iparres: =1 otherwise
c260314 pincl (pscal): four momentum and mass of incident (scatterd) lepon
c       pinch: four momentum and mass of incident hadron
c        vnu: \nu; fq2: Q^2=-q^2; w2l: W^2; yyl: y; zl: z; xb: x_B; pph: P_h

c       para1_1: total cross section of nn, used in parton initiation
c        for nuclus-nucleus collision, this is irrelevant
c        to the inelastic cross section of nn used in PYTHIA
c       para1_2: total cross section of nn, used in hadron cascade
c       dni: nucleon number density
c       dpi: particle (nucleon, pion, kaon ...) number density
c       edi: energy density of particle (nucleon, pion, kaon ...)
c       all of above three densities are calculated in r less or equal
c        2 fm and t less or equal 10 fm/c
c       pj and ej: (pt)**2 and Et of J/psi
c       pjp and ejp: (pt)**2 and Et of (J/psi) prime
c       acoll: array, the demension of which should be larger
c        than or equal to 'nmax'
c       note: the dimension of 'bpp' must be < or = nmax
c       ipden: =0,if projectile is proton
c              =1, projectile is nucleus
c060813 120214
c              =11, projectile is e- (e+)  
c              =12, projectile is nu_e (nu_ebar)  
c              =13, projectile is mu- (mu+)  
c              =14, projectile is nu_mu (nu_mubar)
c              =15, projectile is tau- (tau+)  
c              =16, projectile is nu_tau (nu_taubar)  
c060813 120214
c       itden: =0, if target is proton
c              =1, target is nucleus
c       suppm: the upper bound in sampling the radius of projectile nucleon
c       suptm: the upper bound in sampling the radius of target nucleon
c       suppc: the maximum radius in sampling for projectile
c       suptc: the maximum radius in sampling for target
c       r0p: projectile radius
c       r0t: target radius
c       pio: 3.1416
c       bp: impact parameter
c       iii: current run number
c       coor: position of CM
c       ispmax: maximum # of kinds of particles wanted to statistics
c       ispkf(i): flavor of i-th kind of particle wanted to statistics
c       kfmax: the maximum # of particles with different flavor considered
c       kfaco(i): i-th particle flavor
c       numb(i): sum of particles up to the last one with flavor code of
c         kfaco(i) in particle list
c       an(l,i,j) (san(l,i,j)):
c        l: value of distribution argument (e. g. value of y or pt)
c        i: identify the distribution (e. g. i=2 is y distribution)
c        j: order # of kf code of particle wanded to statistics (1 to ispmax)
c       bn,sbn,sbo: record the multiplicity of particles
c       anf(l,i,j) (sanf(l,i,j)), for instance, is corresponding to an(l,i,j)
c        (san(l,i,j)) but is statistics of full phase space instead of partical
c        phase space
c       isdmax: maximum # of distributions wanted to calculate
c       asd(i): interval segmented for i-th distribution
c260314 for NA,AN and AA collisions
c        i=1: for y
c        i=2: for pt
c         .      .
c         .      .
c         .      .
c260314 for lepton-nucleus collision
c        i=1 : z
c        i=2 : \nu
c        i=3 : Q^2
c         .        .
c         .        .
c         .        .
c       iflmax: maximum # of filters,=0 means no filter at all
c       afl(j,i,1): lower limit of i-th filter for the j-th particle
c       afl(j,i,2): upper limit of i-th filter for the j-th particle
c260314 for NA,AN and AA collisions
c        i=1: y filter
c        i=2: pt filter
c         .        .
c         .        .
c         .        .
c260314 for lepton-nucleus collision
c        i=1 : Q^2=-q^2 (name in program, fq2) filter
c        i=2 : W^2 (w2l) filter
c        i=3 : y (yyl) filter
c        i=4 : P_h (pph) filter
c        i=5 : z (zl) filter
c         .        .
c         .        .
c         .        .
c       ifram: = 0 for fixed target, = 1 for collider
c       cspipi (fm^2): total cross section of pion + pion
c       sig (fm^2): cross section of pion + pion to kaon + kaon
c       cspin (fm^2): total cross section of pion + nucleon interaction
c       cskn (fm^2): total cross section of kaon + nucleon interaction
c       csnn (fm^2): total cross section of n + n interaction
c       rcsit: ratio of inelastic to total cross section
c       disbe(i,j): allowable minimum distance between two particles of
c        kfaco(i) & kfaco(j).
c       c17(i,1-3): three position of particle i
c       tp(i): time of particle i
c       ishp(i): =1 if i-th particle inside the simulated volume
c                =0 otherwise 
c       tau(i): formation time of particle i.
c       isinel(i): = 0 without i-th inelastic process
c                  = 1 with i-th inelastic process
c       nreac(i): statistics of successful i-th collision in parton cascade
c       nrel: statistics of blocked collisions in parton cascade
c       npel: statistics of bolcked nn collisions in parton initiation
c       npinel(600): statistics of successful nn collisions in parton initiation
c       noel : statistics of elastic collisions in hadron cascade
c       noinel(i): statistics the i-th inelastic channel in hadron cascade
c       nosc = 1 : pythia type output only
c              2 : also OSC1999A standard output event-by-event
c              3 : no used
c020708
c       itime is the number of strings in current event
c       astr is the number of strings in current event
c       gtime is the number of gluon in current event
c       akapa(1): sum of string tension over strings in the current event
c       akapa(2): sum of parj(2) over strings in the current event
c       akapa(3): sum of parj(21) over strings in the current event
c       akapa(4): sum of parj(1) over strings in the current event
c       akapa(5): sum of parj(3) over strings in the current event
c020708
c	?????????????????????????????????????
c	0. is the hard distance between two pions
c	0.5 is the hard distance between two nucleons
c	0. is the hard distance between pion and nucleon
c	??????????????????????????????????????
c	win: =incident momentum if iflam=0 
c            =sqrt(s) if iflam=1
c250209 flavor code 22: hardonic decay photon
c                   44: prompt direct photon (<- pythia)
c                   55: photon from parton-parton scattering
c                       qg->q(gamma) and q(-q)->g(gamma)
c                   66: hardonic direct photon
c250209                 pi+pi->rho+(gamma) and pi+rho->pi+(gamma)
c	open(1,file='test.dat',status='old')
c	read(1,100)frame
c	read(1,100)beam
c	read(1,100)target
c	read(1,*)win
c100	format(1x,a)
c	close(1)
2000    open(2,file='prx.out',status='unknown')
        open(3,file='pix.out',status='unknown')
	open(11,file='usu.dat',status='unknown')
	mstu(11)=22
	open(22,file='main.out',status='unknown')
	read(11,*)neve,nout,nosc   
	read(11,*)nap,nzp,nat,nzt
	read(11,*)ddt,dtt,bmin,bmax,nmax   ! 201208 
	read(11,*)kjp21,ifram,para7,para10,kjp20
	read(11,*)pio,ipden,itden
	read(11,*)ispmax,isdmax,iflmax
	read(11,*)(ispkf(i),i=1,10)
	read(11,*)(ispkf(i),i=11,ispmax)
	read(11,*)(asd(i),i=1,isdmax)
	if(iflmax.eq.0)goto 200
	do kk=1,ispmax
	do i=1,iflmax
	read(11,*)(afl(kk,i,j),j=1,2)
	enddo
	enddo
200     read(11,*)parp21,parp22,win   
	read(11,*)ttaup,taujp,iabsb,iabsm,nchan   ! 201208
	read(11,*)para13,para14,psno,para15,para16,ajpsi,vneum
	read(11,*)para1_1,para1_2,para2,para4   
	read(11,*)tdh,cptl,cptu,cptl2,cptu2,itnum   ! 201208
	read(11,*)mstu21,mstp81,mstj1_2,mstj2,mstj3   ! 160617
c210803
	read(11,*)(adj1(i),i=1,10)
	read(11,*)(adj1(i),i=11,20)
        read(11,*)(adj1(i),i=21,30)
        read(11,*)(adj1(i),i=31,40)
c210803
        read(11,*)kjp22,kjp23,kjp24,parp78   ! 020708 020511 150612 yan 070417
	read(11,*)parecc,iparres,smadel   ! 220312 240412 300513
	close(11)
c	tdh and itnum: time step and number of time steps used in subroutine 'flow' 
c	cptl,cptu;cptl2,cptu2 : pt cut in 'flow' for particle 1;particle 2
c       nchan=0: inelastic (INEL)
c	nchan=1: Non Single Difractive (NSD) 
c	nchan=2: qqb --> gamma^*/Z^0, used to generate Drell-Yen
c	nchan=3: J/psi production
c	nchan=4: heavy-flavor production
c	nchan=5: direct photon
c	nchan=6: soft only
c	nchan=7: pythia 

c	neve : # of events to be generate
c	nap (nzp) : # of nucleons (protons) in projectile nucleus
c	nat (nzt) : # of nucleons (protons) in target nucleus
c060813 for e^-A: formally set nap=1,nzp=-1,ipden=11,itden=1, kf=11; 
c       for e^+A: formally set nap=1,nzp=1,ipden=11,itden=1, kf=-11;
c       for nu_eA: formally set nap=1,nzp=-1,ipden=12,itden=1, kf=12;
c       for nu_ebarA: formally set nap=1,nzp=1,ipden=12,itden=1, kf=-12;        
c060813 in hadronic initiation, for instance
c	t0 : average proper formation time at rest
c	ddt : time accuracy
c	dep : the accuracy in four momentum conservation
c	rou0 : normal nuclear density
c	rao : enlarge factor in the radius of simulated volume
c	bmin,bmax : minimum and maximum impact parameters, bmin=bmax means
c	 definite impact parameter, 2*nmax: the number of 
c	 intervals segmented in [bmin,bmax]
c       kjp20: =1 constant cross sections 
c              =0 energy dependent cross sections
c       kjp21: = 0 without hadron rescattering,  
c              = 1 with hadron rescattering
c020708
c       kjp22: = 0 constant string tension with calculation for string 
c                  effective tension
c       kjp22: = 1 constant string tension with calculation of akapa(1 - 5) ! 070417
c       kjp22: = 2 constant string tension   ! 070417
c020511 kjp23: = 1 npart calculated by geometric model
c020511 kjp23: = 2 npart calculated by Glauber model
c020511 kjp24: = 1 sharp sphere in Glauber model
c020511 kjp24: = 2 Woods-Saxon in Glauber model
c020708	
c	param(1)=para1
	param(2)=para2
	param(4)=para4
	param(7)=para7
	param(8)=ddt
	param(10)=para10
	param(13)=para13
c	totle cross-section of J/Psi + n
	param(14)=para14
c	totle cross-section of J/Psi + meson
	param(15)=para15
c	totle cross-section of Psi' + n
	param(16)=para16
c	totle cross-section of Psi' + meson
        idw=adj1(4)   ! 020511
c020511 # of segments in integration   
c       mstp(81)=21   ! let PYEVNW do the whole job.
        mstp(81)=mstp81   ! 160617
c160617 mstp(81)=1 (default) executing pyevnt,=21 executing pyevnw     
	mstp(82)=adj1(8)
c        =0: soft (two-string) only; =1: both of soft and hard
        parp(81)=adj1(9)
c       effective minimum transverse momentum
	parj(2)=adj1(32)
        parj2=parj(2)   ! 020708
c	the suppression of s quark pair production (D=0.3)
	parj(1)=adj1(31)
	parj(3)=adj1(33)
	parj(21)=adj1(34)
c020708
        parj1=parj(1)
        parj3=parj(3)
        parj21=parj(21)
c020708
c	mstp(82)=0   ! without any hard interactions (D=1)
	mstp(91)=adj1(35)
	if(mstp(91).eq.1)parp(91)=adj1(39)
	if(mstp(91).eq.2)parp(92)=adj1(39)
c051108	parp(93)=adj1(30) ! upper cut-off for k_perp distribution in hadron  
c290705	parj(21)=adj1(29)
        parp(2)=parp21
c       parp21: lowest CM energy for calling 'pythia' (D=10.), for 
c	 case of nchan=6
c       parp22: lowest CM energy for calling 'pythia' (D=10.), for 
c	 case of nchan=3
	mstp(33)=1   
c	inclusion of k factor in hard cross sections for parton-parton 
c	 interactions (default=0)
	parp(31)=adj1(10)   ! D=1.5
c070417 contral the strength of colour reconnection
        parp(78)=parp78   ! 070417
c	mstj1_1: =0, no jet fragmentation at all used in parini.f
c	mstj1_2: =1, Lund string fragmentation used in sfm.f
c       no writing of header
        mstp(122)=0   ! 060813 
c	independent fragmentation
	mstu(21)=mstu21   ! 120603
c	gluon jet fragmentation scheme in IF
	mstj(2)=mstj2
c	how the particles share the momentum    
	mstj(3)=mstj3
c
c	parameters in Lund string fragmentation function  
	parj(41)=adj1(6)   ! D=0.3
	parj(42)=adj1(7)   ! D=0.58

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
c070417
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        snnc=0.
        sgtime=0.
        sgtimeo=0.
        sitime=0.
        sastr=0.
        sadiv=0.
        sgpmax=0.
        do i1=1,5
        skapa(i1)=0.
        skapao(i1)=0.
        enddo
        endif
c070417
c061103
	pinel=0.   
	pel=0.   
c061103
	sinel=0.
	sel=0.
	do i1=1,600
	dinel(i1)=0.
	enddo
	do i1=1,20
	acoll(i1)=0.
	acollp(i1)=0.
	acollt(i1)=0.
	sbp(i1)=0.
	enddo

	volum=4.*3.1416/3.*2.**3
c261002	volume of sphere with radius of 2 fm in position space 
c010220
        ncha=0
        nchaf=0
c010220
c200601
	skpar=0.
	sknn=0.
	skpp=0.
	sknp=0.
	skep=0. ! statistic of lepton-p collisions with calling pythia 060813
c200601
	sthroq=0.
	sthrob=0.
        sthroc=0.
	do i1=1,4
	sthroe(i1)=0.
	enddo
	adj12=adj1(12)
	dpmax=adj1(27)   ! 281194
	adj140=adj1(40)   ! 290505
c	give values to some important variables
	call sysini(win)   ! 060813
	adj1(28)=para10*max(rnt,rnp)
	iii=0
	open(5,file='sxp.out',status='unknown')   ! sa 26/05/99
	open(9,file='rms0.out',status='unknown')
	open(34,file='oscar.out',status='unknown')
c020511
        csnn1=csnn*10   ! csnn in fm^2 csnn1 in mb
        idw1=idw/50   ! *100
c        write(9,*)'csnn1,kjp24,idw1=',csnn1,kjp24,idw1
        if(ipden.lt.2)call overlap(nap,nat,rnp,rnt,csnn1,kjp23,kjp24,
     c	 rou0,idw1)   ! 060813 050214
c020511
	write(9,*)'nap,nzp,nat,nzt,win=',nap,nzp,nat,nzt,win
	write(9,*)'neve,nout,nosc=',neve,nout,nosc
	write(9,*)'bmin,bmax,dtt,nmax,parp78=',bmin,bmax,dtt,nmax,
     c	 parp78 c150612 yan 070417
	write(9,*)'kjp21,ifram,para7,para10,kjp20,kjp22,kjp23,kjp24=',
     c   kjp21,ifram,para7,para10,kjp20,kjp22,kjp23,kjp24   ! 020511
	write(9,*)ispmax,isdmax,iflmax
	write(9,*)(ispkf(i),i=1,ispmax)
	write(9,*)(asd(i),i=1,isdmax)
	write(9,*)'parp21,parp22,ttaup,taujp=',parp21,parp22,ttaup,taujp
	write(9,*)'iabsb,iabsm,nchan=',iabsb,iabsm,nchan
        write(9,*)'para13 - 16 =',para13,para14,para15,para16
	write(9,*)'psno,ajpsi,vneum=',psno,ajpsi,vneum
        write(9,*)'para1_1,para1_2,para2,para4=',para1_1,para1_2,para2,
     c	 para4   
	write(9,*)'tdh,itnum,cptl=',tdh,itnum,cptl
	write(9,*)'cptu,cptl2,cptu2=',cptu,cptl2,cptu2
	write(9,*)'mstu21,mstp81,mstj1_2,mstj2,mstj3=',
     c	 mstu21,mstp81,mstj1_2,mstj2,mstj3   ! 160617  
c210803
	write(9,*)'adj1=',(adj1(i),i=1,10)
	write(9,*)'adj1=',(adj1(i),i=11,20)
        write(9,*)'adj1=',(adj1(i),i=21,30)
        write(9,*)'adj1=',(adj1(i),i=31,40)
c210803
c	write(9,*)'41-45=',parj(41),parj(42),parj(43),parj(44),
c     c	 parj(45)
        write(9,*)'parecc,iparres,smadel=',parecc,iparres,smadel
c220312 240412 300513
	if(iflmax.ne.0)then
	do kk=1,ispmax
	do i=1,iflmax
	write(9,*)(afl(kk,i,j),j=1,2)
	enddo
	enddo
	endif
c	write(9,*)'adj12=',adj12   ! sa
	nmin=0
	nminf=0
	ich=0	
	time=0.
c       write(9,*)'rnt,para10,adj1(28)=',rnt,para10,adj1(28)   ! sa
	jjj=1
c280113	psno: =0 fixed impact parameter 
c291207       =1 systematic sampling method
c291207       =2 random sampling method
c	for given b (impact parameter)
	if(abs(bmin-bmax).lt.10d-4)then   ! i. e. case of psno=0. 280113
	bp=bmin
        r4=rnp
        if(rnt.gt.rnp)r4=rnt
        rr4=bp/r4
        vneu=exp(-rr4*rr4)
c	calculate the overlap region of two nuclei at given b, geometric
c020511  or Glauber method
c060605
c280113	if(nap.ne.1 .and. nat.ne.1)then   ! 020511
c020511
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
	anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
        pir=part1(ibpp)
        tir=part2(ibpp)
c       write(9,*)'bp,ibpp,part1,part2=',bp,ibpp,pir,tir
c280113	endif
c020511
	if(ipden.lt.2.and.nap.eq.1)pir=1.
	if(nat.eq.1)tir=1.
c060605
	vneump=pir   ! 111399
	vneumt=tir   ! 111399
	write(9,*)'psno,b,N_part_p,N_part_t,N_bin=',
     c	 psno,bp,vneump,vneumt,anbin*csnn   ! 190309 280113
	goto 300
	endif
	if(psno.eq.1.)then   ! 280113
c	systematic sampling method for given interval of b according to b**2 law
	nmax2=2*nmax
	bmaxn=(bmax*bmax-bmin*bmin)/(2.*nmax)
	bmin2=bmin*bmin
	i2=0
	do i1=1,nmax2,2
	i2=i2+1
	bpp(i2)=sqrt(i1*bmaxn+bmin2)
	enddo
	write(9,*)'b=',(bpp(i1),i1=1,i2)   !!
	stab=0.   ! 280113
	stb=0.
	stbp=0.   ! 111399
	stbt=0.   ! 111399
	do i1=1,i2
	bp=bpp(i1)
	r4=rnp
	if(rnt.gt.rnp)r4=rnt
	rr4=bp/r4
	acoll(i1)=exp(-rr4*rr4)
c	calculate the overlap region of two nuclei at given b 
c060605
c280113	if(nap.ne.1 .and. nat.ne.1)then   ! 020511
c020511
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
	anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
        pir=part1(ibpp)
        tir=part2(ibpp)
c280113	endif
c020511
	if(ipden.lt.2.and.nap.eq.1)pir=1.
	if(nat.eq.1)tir=1.
c060605
	stab=stab+bp   ! 280113
	acoll(i1)=anbin   ! 280113
	acollp(i1)=pir
	acollt(i1)=tir
c	write(*,*)'bp,din,vneu,tb=',bp,din,acoll(i1),tb
	stbp=stbp+acollp(i1)
	stbt=stbt+acollt(i1)
	stb=stb+acoll(i1)
	enddo
	stab=stab/float(i2)   ! 280113
	aneump=stbp/float(i2)   ! 191110
	aneumt=stbt/float(i2)   ! 191110
	vneum=stb/float(i2)
	write(9,*)'psno,ave. b=',psno,stab   ! 280113
        write(9,*)'N_bin=',(acoll(i1)*csnn,i1=1,i2)   !! 280113
	write(9,*)'(N_part)_p=',(acollp(i1),i1=1,i2)   ! 191110
	write(9,*)'(N_part)_t=',(acollt(i1),i1=1,i2)   ! 191110
	write(9,*)'ave. N_part_p,N_part_t,N_bin=',
     c	 aneump,aneumt,vneum*csnn ! 191110	280113
	endif   ! 280113
c	average b in [bmin,bmax]
        avb=2./3.*(bmin+bmax)
c	above equation is correct when bmin=0 only
	r4=rnp
	if(rnt.gt.rnp)r4=rnt
	rr4=avb/r4
	avneu=exp(-rr4*rr4)
c	calculate the overlap region of two nuclei at given b 
c060605
c280113	if(nap.ne.1 .and. nat.ne.1)then   ! 020511
c020511
        ibpp=int(avb/0.1+1.0)
        ibpp=min(ibpp,200)
	anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
        pir=part1(ibpp)
        tir=part2(ibpp)
c280113	endif
c020511
	if(pden.lt.2.and.nap.eq.1)pir=1.
	if(nat.eq.1)tir=1.
c060605
	aanbin=anbin   ! 280113
	astbp=pir
	astbt=tir
c280809
        nrel=0
        nrea=0
        do i=1,9
        nreac(i)=0
        enddo
c280809
c280113
	if(psno.eq.2.)then
	averb=0.
	psnon=0.   ! N_bin in case of psno=2
	psnop=0.   ! parojectile N_part in case of psno=2
	psnot=0.   ! target N_part in case of psno=2
	endif 
c280113
c	generate a event for nucleus-nucleou collision
	iran=adj1(26)
	if(iran.eq.0)goto 300
	do i1=1,iran
	thrr=pyr(1)
	enddo
	vnlep=0.d0 ! statistics of the number of studied leptons 260314
300	iii=iii+1
c070905
        ithroq=0
        ithrob=0
        ithroc=0
        ithroq_p=0
        ithrob_p=0
        ich_p=0
        ithroq_t=0
        ithrob_t=0
        ich_t=0
        do i=1,4
        throe(i)=0.
        throe_t(i)=0.
        throe_p(i)=0.
        enddo
	iprlth=0
        do j=1,mplis
        idpth(j)=0
        rmpth(j)=0.
        do i=1,4
        rpth(i,j)=0.
        ppth(i,j)=0.
        enddo
        enddo
c070905
	nncoll=0
c061103
	npel=0   
        do i=1,600
        npinel(i)=0
        enddo   
c071103
c071103
	noel=0
        do i=1,600
        noinel(i)=0
        enddo
c061103
	nspe=0   ! 111899
	if(abs(bmin-bmax).lt.10d-4)goto 800   ! moved two lines forward 280113
	if(psno.eq.1.)then   ! 280113
	sbp(jjj)=sbp(jjj)+1
	bp=bpp(jjj) 
c191110
	vneump=acollp(jjj)
	vneumt=acollt(jjj)
	goto 800   ! 280113
	endif   ! 280113
c191110 if(psno.eq.2)bp=sqrt(pyr(1)*(bmax*bmax-bmin2)+bmin2)   ! 291207 
        if(psno.eq.2)then
	bp=sqrt(pyr(1)*(bmax*bmax-bmin2)+bmin2)
c       calculate the overlap region of two nuclei at given bp 
c280113	if(nap.ne.1 .and. nat.ne.1)then   ! 020511
c020511
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
	anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
        pir=part1(ibpp)
        tir=part2(ibpp)
c280113	endif
c020511
        if(ipden.lt.2.and.nap.eq.1)pir=1.
        if(nat.eq.1)tir=1.
        vneump=pir   
        vneumt=tir  
	if(iii.eq.1)write(9,*)'psno,b,N_part_p,N_part_t,N_bin=',
     c	 psno,bp,vneump,vneumt,anbin*csnn   ! 280113
	endif
c191110
800     continue
c280113
	if(psno.eq.2.)then
	averb=averb+bp
	psnon=psnon+anbin
	psnop=psnop+vneump
	psnot=psnot+vneumt
	endif 
c280113
        open(11,file='usu.dat',status='unknown')
        read(11,*)neve,nout
        close(11)
c	write(9,*)'after close(11)'   !
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
	do i1=1,5
	ppsa(i1)=0.
	enddo
c	forbiden decay of particle, if mdcy(...)=0
c291207 mdcy(pycomp(111),1)=0   ! 190606
c291207 mdcy(pycomp(221),1)=0   ! 190606
c251108        mdcy(pycomp(310),1)=0   ! k0_S
c	mdcy(pycomp(333),1)=0   ! phi
c210313	mdcy(pycomp(3122),1)=0   ! Lambda0
c210313	mdcy(pycomp(-3122),1)=0
c	mdcy(pycomp(443),1)=0   ! j/psi
c	mdcy(pycomp(10441),1)=0
c	mdcy(pycomp(20443),1)=0
c	mdcy(pycomp(445),1)=0
c	mdcy(pycomp(411),1)=0
c	mdcy(pycomp(-411),1)=0
c	mdcy(pycomp(421),1)=0
c	mdcy(pycomp(-421),1)=0
c	mdcy(pycomp(4122),1)=0
c	mdcy(pycomp(4112),1)=0
c	mdcy(pycomp(4212),1)=0
c	mdcy(pycomp(4222),1)=0
c	mdcy(pycomp(30443),1)=0   ! j/psi

c	mdcy(pycomp(3212),1)=0
c	mdcy(pycomp(-3212),1)=0
c251108        mdcy(pycomp(3112),1)=0   ! Sigma-
c	mdcy(pycomp(-3112),1)=0
c251108        mdcy(pycomp(3222),1)=0   ! Sigma+
c	mdcy(pycomp(-3222),1)=0
c251108        mdcy(pycomp(3312),1)=0   ! Xi-
c251108        mdcy(pycomp(-3312),1)=0
c	mdcy(pycomp(3322),1)=0
c	mdcy(pycomp(-3322),1)=0
c251108        mdcy(pycomp(3334),1)=0   ! Omega-
c251108        mdcy(pycomp(-3334),1)=0
c       mdcy(pycomp(1114),1)=0
c       mdcy(pycomp(2114),1)=0
c       mdcy(pycomp(2214),1)=0
c       mdcy(pycomp(2224),1)=0
c       mdcy(pycomp(213),1)=0   ! rho+
c       mdcy(pycomp(-213),1)=0   ! rho-  
c	mdcy(pycomp(113),1)=0   ! rho0
c	mdcy(pycomp(223),1)=0   ! omega 041202
c       mdcy(pycomp(413),1)=0
c       mdcy(pycomp(-413),1)=0
c       mdcy(pycomp(423),1)=0
c       mdcy(pycomp(-423),1)=0
c	mdcy(pycomp(13),1)=0
c	mdcy(pycomp(-13),1)=0
	ijk=0   ! 10/08/98
c181003
	nnstop=0
        zstop=0.
c181003
c       ??????????????? oscar stander output ??????????????????????
        if(nosc.eq.2 .or. nosc.eq.3)then   ! osc
        write(34,501)  ! osc
501     format(8hOSC1999A)   ! osc
        write(34,502)   ! osc
502     format(18hfull-event-history)   ! osc
        write(34,503)   ! osc
503     format(10hPACIAE 1.0)   ! osc
        ntest=1   ! osc
        if(ifram.eq.0)then   ! osc
        ebeam=sqrt(win*win+0.938*0.938)   ! osc
        write(34,504)nap,nzp,nat,nzt,ebeam,ntest   ! osc
        endif   ! osc
        if(ifram.eq.1)then   ! osc
        ebeam=win*win/2./0.938-1.   ! osc
        write(34,505)nap,nzp,nat,nzt,ebeam,ntest   ! osc
        endif   ! osc
504     format(1h(,i3,1h,,i3,1h),1h+,1h(,i3,1h,,i3,1h),1x,3hlab,f9.3,i2)
     c     ! osc
505     format(1h(,i3,1h,,i3,1h),1h+,1h(,i3,1h,,i3,1h),1x,4haacm,f9.3,
     c   i2)  ! osc
        write(34,*)'event number=',iii   ! osc
        write(34,*)'impact parameter=',bp   ! osc
        write(34,*)'phi=0.'   ! osc
        endif   ! osc
c       ??????????????? oscar stander output ??????????????????????
c	creat an event
c	parton initiation and administration for an event generation  
	ijk=0 
cs	write(9,*)'be. parini iii,itden=',iii,itden   ! sa
	call parini(time_neu,time_par,time_had,parp21
     c	 ,parp22,win,psno,ijk)   ! 111010 121110 240513
c060605	for p+A,A+p,lepton+A, and A+A   ! 240513 060813 120214
cs	write(9,*)'af. parini ijk=',ijk ! sa
        if(ijk.eq.1)goto 300 ! 071005 to avoide infinite loop in parcas 060813 
	rrp=1.16   ! 130205
c070905
c	goto 777   ! temporary
	if(adj12.ne.0)then   ! 2
        if(ipden.lt.11)call pyedit(2)   ! 060813
	if(ipden.ge.11)call pyedit(1)   ! 060813
c	if(nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
c	c   .eq.neve)then   ! 3
c	write(mstu(11),*)'event=',iii
c	call pylist(1)
c	write(mstu(11),*)'ppsa=',(ppsa(i1),i1=1,5)
c	write(22,*)'throe_t=',throe_t
c	write(22,*)'ithroq_t,ithrob_t,ich_t=',ithroq_t,ithrob_t,
c	c   float(ich_t)/3.
c	write(22,*)'throe_p=',throe_p
c	write(22,*)'ithroq_p,ithrob_p,ich_p=',ithroq_p,ithrob_p,
c	c   float(ich_p)/3.
c	write(22,*)'throe=',throe
c	write(22,*)'ithroq,ithrob,ithroc=',ithroq,ithrob,
c	c   float(ithroc)/3.
c	call prt_sgam(ngam)   ! 240209
c       call prt_sbh(nbh)
c	endif   ! 3
c	recoalescent the q and qbar
c       reconstruct parton list ('parlist')
c       make the partons in order of qba and q
        iij=0
        jji=0   ! 291207
        do ii=2,3   ! 1,3
c       1: refers to g, 2: qba, 3: q
        kf=21
        do j=iij+1,iprlth
        call ord_th(jji,j,kf,ii)   ! 291207
        enddo
        iij=jji   ! 291207
        numbth(ii)=jji   ! 291207
c       numbth(1),(2) and (3): the line # of last g,qba and q
        enddo
c271004 n1=numbth(1)
        n2=numbth(2)
        n3=numbth(3)
c	write(9,*)'af ord_th n1,n2,n3,iprlth,ithroq,ithrob=',
c     c   n1,n2,n3,iprlth,ithroq,ithrob   ! sa
c061005	write(9,600)(idpth(i2),i2=1,n3)   ! sa
        iqba=ithrob
        n3=iqba+ithroq
        do i1=1,n3
        idp(i1)=idpth(i1)
        rmp(i1)=rmpth(i1)
        do j1=1,4
        pp(j1,i1)=ppth(j1,i1)
        rp(j1,i1)=rpth(j1,i1)
        enddo
        enddo
        do i1=n3+1,mplis
        idp(i1)=0
        rmp(i1)=0.
        do j1=1,4
        pp(j1,i1)=0.
        rp(j1,i1)=0.
        enddo
        enddo
        ithroq=0
        ithrob=0
        ithroc=0
        do i=1,4
        throe(i)=0.
        enddo
	iprlth=0
        do j=1,mplis
        idpth(j)=0
        rmpth(j)=0.
        do i=1,4
        rpth(i,j)=0.
        ppth(i,j)=0.
        enddo
        enddo
        nn=0
        do i1=1,kszj
        do j1=1,5
        kn(i1,j1)=0.
        pn(i1,j1)=0.
        rn(i1,j1)=0.
        enddo
        enddo
        call coal(n3,iqba,iii,rrp,0,0)
c       'sa1_h' to 'PYJETS'
        if(nn.ne.0)then   ! 4
        do l=1,nn
        l1=n+l
        do m=1,5
        k(l1,m)=kn(l,m)
        p(l1,m)=pn(l,m)
        v(l1,m)=rn(l,m)
        enddo
        enddo
        n=n+nn
        endif   ! 4
	endif   ! 2
c070905 
	goto 777
c       perform particle, declared unstable in 'mdcy' array, decay
c130205	call pyexec
777	continue   ! 070905
	if(adj140.ge.4)call decayh(rrp)   ! 130205 140414
	if(iii.eq.1)then
	write(9,*)'iii,neve=',iii,neve
	write(9,*)'nap,nzp,nat,nzp,bp=',nap,nzp,nat,nzt,bp
	write(9,*)'sig,t0,ddt,dep=',sig,t0,ddt,dep
	write(9,*)'rou0,rao,rnp,rnt=',rou0,rao,rnp,rnt
	write(9,*)'csnn,cspin,cskn=',csnn,cspin,cskn
	write(9,*)'cspipi,cspsn,cspsm=',cspipi,cspsn,cspsm
	write(9,*)'ifram,rcsit,kfmax,ipden,itden=',
     c	 ifram,rcsit,kfmax,ipden,itden   ! 060813
	write(9,*)(kfaco(i),i=1,kfmax)
	write(9,*)(disbe(i,i),i=1,kfmax)
	write(9,*)(disbe(1,i),i=1,8)
	write(9,*)'isinel='
	write(9,600)isinel
	endif
c050603
        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813
	if(nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
     c	 .eq.neve)then 
	write(mstu(11),*)'event=',iii
	call pylist(1)
        write(mstu(11),*)'ppsa=',(ppsa(i1),i1=1,5)
	write(22,*)'throe_t=',throe_t
	write(22,*)'ithroq_t,ithrob_t,ich_t=',ithroq_t,ithrob_t,
     c	 float(ich_t)/3.
	write(22,*)'throe_p=',throe_p
	write(22,*)'ithroq_p,ithrob_p,ich_p=',ithroq_p,ithrob_p,
     c	 float(ich_p)/3.
	write(22,*)'throe=',throe
	write(22,*)'ithroq,ithrob,ithroc=',ithroq,ithrob,
     c	 float(ithroc)/3.
        call prt_sgam(ngam)   ! 240209
	if(itden.ne.1)call prt_sbh(nbh)
	endif  
c050603
600	format(25(1x,i2)/)
	if(mod(iii,nout).eq.0)then
	print*,'event=',iii
	endif

c       ??????????????? oscar stander output ??????????????????????
        if((nosc.eq.1.and.(mod(iii,nout).eq.0)).or.nosc.eq.2)
     c	 call oscar(iii)   ! 160705   
c       ??????????????? oscar stander output ??????????????????????
c	analyse an event
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
	yy=pyp(j,17)
	ppt=pyp(j,10)
c281104
	ppm=sqrt(p1*p1+p2*p2+p3*p3)
	if(ppm.le.dpmax.and.p4.le.dpmax)then
	goto 3000
	else
	goto 400
	endif
c281104
3000	continue
        if((itden.eq.0.and.ipden.eq.1).or.(itden.eq.1.and.ipden.eq.0)
     c   .or.(itden.eq.1.and.ipden.eq.1))then   ! 260314
	c(1)=yy   
	if(ifram.eq.1)c(1)=eta   
	c(2)=ppt
c	.
c	.
c	.
	kkk=1
c033101
c	statistics of negative multiplicity
	if(adj140.ge.4. .and. plu6.lt.-0.9)then ! for hadron 140414 (.eq.)
	nminf=nminf-plu6
	do i=1,iflmax
	if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 700
	enddo
	nmin=nmin-plu6	
700	endif
	if(adj140.lt.4. .and. plu6.lt.-0.2)then ! for parton 140414 (ne.)
        nminf=nminf+1.   ! -plu6 230206
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 702
        enddo
        nmin=nmin+1.   ! -plu6 230206
702     endif
c010220
c       statistics of charged particles multiplicity
c033101	statistics positive multiplicity
c070802
c       if(abs(plu6).gt.0.9 .and. iabs(ik).ne.11)then
        if(adj140.ge.4. .and. plu6.gt.0.9)then   ! for hadron 140414
c070802
        nchaf=nchaf+plu6
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 701
        enddo
        ncha=ncha+plu6
701     endif
	if(adj140.lt.4. .and. plu6.gt.0.2)then   ! for parton 140414
c070802
        nchaf=nchaf+1.   ! +plu6 230206
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 703
        enddo
        ncha=ncha+1.   ! +plu6 230206
703     endif
	endif   ! 260314
c260314	statistics of y, pt, ect. distributions (for NA,AN and AA); z, \nu,
c        ect. distributions (for lepton-nucleus)
	do 500 kk=1,ispmax
	kf=ispkf(kk)
	if(ik.ne.kf)goto 500

c????????????????????????????????????????????????????????????????????
        if(ik.eq.2212 .or. ik.eq.2112)then
c	eep=sqrt(win*win+0.938*0.938)
c	if(ifram.eq.1)eep=0.5*win
c	eepd=eep-0.001
c	eepu=eep+0.001
c	if(ifram.eq.0 .and. ((p4.gt.eepd .and. p4.lt.eepu) .or.
c     c   p4.le.0.940))then   ! 111899
c	nspe=nspe+1   ! 111899
c	goto 500
c	endif   ! 111899
c	if(ifram.eq.1 .and. (p4.gt.eepd .and. p4.lt.eepu))then
	if(ppt.le.1.d-4)then   
	nspe=nspe+1   ! 111899
	goto 500
	endif   ! 111899 
        endif
c       exclude the projectile and the target spectator nucleons
c????????????????????????????????????????????????????????????????????
840	continue   ! 260314
c260314	case of nucleon or nucleus incidence
        if((itden.eq.0.and.ipden.eq.1).or.(itden.eq.1.and.ipden.eq.0)
     c   .or.(itden.eq.1.and.ipden.eq.1))
     c	 call stati_h(yy,ppt,eta,p5,ik,kk,w,bn,an,bnf,anf)   ! 010600
c       case of lepton incidence 
        if(ipden.ge.11.and.ipden.le.16)
     c   call stati_l(p1,p2,p3,p4,p5,ik,kk,w,bn,an,bnf,anf)
c260314
	goto 901   ! sa 26/05/99 (originally  'goto 400')
500	continue
c?????? follow three statements are increaded on 26/05/99 by sa ******
901     if(nout.ne.1)then
csa        write(5,900)ik,x1,y1,z1,t1,p1,p2,p3,p4,p5
	endif
c********************************************************************
400	continue
c	statistics of multiplicity distributions,
c	 spectator nucleons are excluded
	do kkk=1,ispmax
c	ik=ispkf(kkk)
c        if(abs(ik).eq.2212.or.abs(ik).eq.2112.or.abs(ik).eq.3122.or.
c     &   abs(ik).eq.3212.or.ik.eq.3222.or.ik.eq.3112)then
c	multiplicity is located at which interval
c	if(ik.eq.2212)then
	idf=bnf(kkk)/asd(5)+1
c	number five distribution is multiplicity distribution
	if(idf.lt.1 .or. idf.gt.20)goto 405   ! 131204
	anf(idf,5,kkk)=anf(idf,5,kkk)+1./asd(5)
c       put the filter to effect
405	do i=1,iflmax   ! 131204
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 404
        enddo
	idd=bn(kkk)/asd(5)+1
	if(idd.lt.1 .or. idd.gt.20)goto 404   ! 131204
	an(idd,5,kkk)=an(idd,5,kkk)+1./asd(5)
404	continue
c	endif
	enddo
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
c070417
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        snnc=snnc+nnc
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
	pel=pel+npel   ! 061103   
	pinel=pinel+npinel(592)   
	sel=sel+noel
	do i1=1,600
	sinel=sinel+noinel(i1)
	dinel(i1)=dinel(i1)+noinel(i1)
	enddo
c071103
        rinel=0.   ! 230210
	do i1=1,9
	rinel=rinel+nreac(i1)
	enddo
c071103
c200601
	skpar=skpar+kpar
	sknn=sknn+knn
	skpp=skpp+kpp
	sknp=sknp+knp
	skep=skep+kep   ! 060813
c200601
	sthroq=sthroq+ithroq+ithroq_p+ithroq_t
	sthrob=sthrob+ithrob+ithrob_p+ithrob_t
        sthroc=sthroc+ithroc+ich_p+ich_t
	do i1=1,4
	sthroe(i1)=sthroe(i1)+throe(i1)+throe_p(i1)+throe_t(i1)
	enddo

	open(8,file='nout.out',status='unknown')
	write(8,*)'iii=',iii
	close(8)
c	internal and final printing and controled return
	if(mod(iii,nout).eq.0 .or. iii.eq.neve)then
	open(10,file='rms.out',status='unknown')
	flaa=float(iii-ich)   ! July/20/98
	if(flaa.le.1.d-20)goto 1200
c280113
	if(psno.eq.2)then
	averbo=averb/flaa
	psnono=psnon/flaa
	psnopo=psnop/flaa
	psnoto=psnot/flaa
	endif 
c280113
	dnmino=nmin/flaa
	dnminfo=nminf/flaa
c010200
        dncha=ncha/flaa
        dnchaf=nchaf/flaa
c010220
	do kk=1,ispmax
	sbo(kk)=sbn(kk)/flaa
	sbof(kk)=sbnf(kk)/flaa
	do i1=1,20
	do i2=1,isdmax
	sao(i1,i2,kk)=san(i1,i2,kk)/flaa
	saof(i1,i2,kk)=sanf(i1,i2,kk)/flaa
	enddo
	enddo
	enddo
c070417
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        snnco=snnc/flaa
        sastro=sastr/flaa
        sgtimeo=sgtime/flaa
c       sgtimeo: average number of gluons in a string over strings with gluon
        do i1=1,5
        skapao(i1)=skapa(i1)/flaa
        enddo
        sadivo=sadiv/flaa
        sgpmaxo=sgpmax/flaa
        sitimeo=sitime/flaa
        endif
c070417
c061103
        peli=pel/flaa
        pineli=pinel/flaa
c061103
c071103
	reli=nrel/flaa   ! original =rel/flaa
	rineli=rinel/flaa
c071103
	seli=sel/flaa
	sineli=sinel/flaa
	do i1=1,600
	dineli(i1)=dinel(i1)/flaa
	enddo
c200601
	skparo=skpar/flaa
	sknno=sknn/flaa
	skppo=skpp/flaa
	sknpo=sknp/flaa
	skepo=skep/flaa   ! 060813
c200601
	wthroq=sthroq/flaa
	wthrob=sthrob/flaa
        wthroc=sthroc/flaa
	do i1=1,4
	wthroe(i1)=sthroe(i1)/flaa
	enddo
c280809
        srea=0.   ! 230210
        do i1=1,9
c230210 nrea=nrea+nreac(i1)
        snreac(i1)=nreac(i1)/flaa
        srea=srea+snreac(i1)   ! 230210
        enddo
c230210 srea=float(nrea)/flaa
c280809
1200	continue
	write(10,*)'mstp82,parp81,bp=',mstp(82),parp(81),bp   ! 291207
	write(10,*)'nn colli. # and blocked # in parton initialization=',
     c	 pineli,peli,peli+pineli
c071103
	write(10,*)'# of successful and blocked collision in parton 
     c	 cascade=',rineli,reli,reli+rineli
c071103
c280809
        write(10,*)'average collision # in parton cascade=',srea
        write(10,*)'total # of scaterring processes in parton cascade'
        write(10,*)(snreac(i1),i1=1,9)
c280809
	write(10,*)'el. and inela. colli. # and sum in hadron cascade=',
     c	 seli,sineli,seli+sineli
c200601
        write(10,*)'(Npart)mini-jet,Nnn,Npp=',skparo,sknno,skppo
        write(10,*)'Nnp,Ntot,Nep=',sknpo,sknno+skppo+sknpo,skepo   ! 060813
c200601
	if(psno.eq.2)write(10,*)'psno, ave. b,N_part and N_bin=',
     c	 psno,averbo,psnopo,psnoto,psnono*csnn   ! 280113
	if(ipden.ge.11.and.ipden.le.16)
     c   write(10,*)'event average number of lepton studied=',vnlep/flaa !260314
c070417
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        write(10,*)'kjp22=0, par1,par2,par3,par21=',
     c   parj1,parj2,parj3,parj21
        write(10,*)'keff2,par2,par21,par1,par3=',
     c   (skapao(i1),i1=1,5)
        write(10,*)'averaged # of gluon in a string,averaged # of'
        write(10,*)' strings in an event (counted with KF=92) and'
        write(10,*)'average # of NN collisions in an event=',sgtimeo,
     c   sastro,snnco   ! 070417
        write(10,*)'event averaged value of the factor related to # of'
        write(10,*)'gluons and hardest gluon in a string,event averaged'
        write(10,*)'transverse momentum of hardest gluon,event averaged'
        write(10,*)'# strings=',sadivo,sgpmaxo,sitimeo
        endif
c070417
	write(10,*)'multiplicity of negative particles=',dnmino
	write(10,*)'multiplicity of negative particles=',dnminfo
	write(10,*)'multiplicity of positive particles,partial=',dncha
	write(10,*)'multiplicity of positive particles,full=',dnchaf
	write(10,*)'throw away ithroq,ithrob,ithroc=',
     c	 wthroq,wthrob,wthroc/3.
	write(10,*)'throe=',wthroe
        write(10,*)'avb,avneu,astbp,astbt,aanbin=',
     c	 avb,avneu,astbp,astbt,aanbin   ! 280113
        write(10,*)'particle multiplicity=',(sbof(ll),ll=1,ispmax)
	write(10,*)'particle multiplicity=',(sbo(ll),ll=1,ispmax)
601	format(8(1x,f6.4))
csa****************************************************************

	do m2=1,isdmax
	write(10,*)'ID of distribution m2=',m2
	do m3=1,ispmax
	write(10,*)'distribution belong to m3=',m3
	write(10,*)(sao(m1,m2,m3),m1=1,20)
	write(10,*)(saof(m1,m2,m3),m1=1,20)
	enddo
	enddo
c260314
	if(ipden.ge.11.and.ipden.le.16)then
        do kk=1,10
	sbn1=sbn(1)
	sbn1=dmax1(sbn1,1.d-20)
	sbnf1=sbnf(1)
	sbnf1=dmax1(sbnf1,1.d-20)
        if(kk.ne.1)sbo(kk)=sbn(kk)/sbn1
	if(kk.ne.1)sbof(kk)=sbnf(kk)/sbnf1  
        do i1=1,20
        do i2=1,isdmax
	san1=san(i1,i2,1)
	san1=dmax1(san1,1.d-20)
	sanf1=sanf(i1,i2,1)
	sanf1=dmax1(sanf1,1.d-20)
        if(kk.ne.1)sao(i1,i2,kk)=san(i1,i2,kk)/san1   
	if(kk.ne.1)saof(i1,i2,kk)=sanf(i1,i2,kk)/sanf1   
        enddo
        enddo
        enddo
        write(10,*)'relative multiplicity,p=',(sbo(ll),ll=1,10)
	write(10,*)'relative multiplicity,f=',(sbof(ll),ll=1,10)   
        do m2=1,isdmax
        write(10,*)'ID of relative distribution m2=',m2
        do m3=1,10
        write(10,*)'distribution belong to m3=',m3
        write(10,*)(sao(m1,m2,m3),m1=1,20)  
	write(10,*)(saof(m1,m2,m3),m1=1,20)   
        enddo
        enddo
	endif
c260314
	write(10,*)'average frequency of the occurring of each inela.'
	write(10,*)dineli

	close(10)
	endif

1000	if(iii.lt.neve)then
	if(abs(bmin-bmax).lt.10d-4)goto 300
	if(psno.eq.2.)goto 300   ! 280113
	if(jjj.ge.10)then
c	10 : the total number of impact paremeters in systematic sampling for
c            impact  parameter
	jjj=1
	goto 300
	endif
	jjj=jjj+1
	goto 300
	endif
	
c060813	statistics of processes generated
	call pystat(0)   ! 060813
	close(2)
	close(3)
	close(5)
	close(9)
	close(22)
	close(34)
900 	format(i5,8(1x,f9.3),1x,f5.3)
c	timeb=dtime(ty)
c	write(9,*)'time consuming =',timeb
	stop
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine stati_h(y,pt,eta,p5,ik,kk,ww,a,b,af,bf)   ! 260314
c	on line statistics for NA,AN,and AA collisions   ! 260314
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c  afl(20,5,2)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &  iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
	dimension a(20),b(20,5,20),c(5),af(20),bf(20,5,20),id(5)
	amass=p5   ! 010600
        amass2=amass*amass
        pt2=pt*pt
        et=sqrt(pt2+amass2)
        do 10000 i=1,iflmax
        goto (10,20,30,40,50) i
10	c(i)=y   
	if(ifram.eq.1)c(i)=eta   
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
100	ii=abs(y)/asd(i)+1
	if(ifram.eq.1 .and. y.gt.0.)ii=ii+10
        if(ifram.eq.1 .and. y.lt.0.)ii=10-ii+1
c       note: 10 here should be change together with the dimension
c        20
	id(i)=ii
c100	id(i)=y/asd(i)+1
c	if(id(i).le.0)id(i)=1
	goto 20000
c	pt is located in which interval?
200	id(i)=pt/asd(i)+1
c	if(id(i).le.0)id(i)=1
	goto 20000
c	eta is located in which interval?
300	ii=abs(eta)/asd(i)+1
	if(ifram.eq.1 .and. eta.gt.0.)ii=ii+10
        if(ifram.eq.1 .and. eta.lt.0.)ii=10-ii+1
c       note: 10 here should be change together with the dimension
c        20
	id(i)=ii
	goto 20000
c	et is counted in which eta interval ?
400	id(i)=id(i-1)
	if(kk.eq.3)id(i)=p5/0.1+1
c	statistcs the mass dis. of rho0
        goto 20000
500	continue
20000	continue
c	make statistics for particle yield and desired distributions
	af(kk)=af(kk)+ww
	do i=1,isdmax
        ii=id(i)
c010218	if(ii.gt.20)goto 30000
	if(ii.lt.1 .or. ii.gt.20)goto 30000   ! 010218
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
c010218	if(ii.gt.20)goto 50000
	if(ii.lt.1 .or. ii.gt.20)goto 50000   ! 010218
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
c       on line statistics for lepto-nucleus collision
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c  afl(20,5,2)
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &  iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
	common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph
     c	 ,vnlep
	dimension a(20),b(20,5,20),af(20),bf(20,5,20),c(5),id(5)
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
        if(ii.lt.1 .or. ii.gt.20)goto 30000
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
	if(ii.lt.1 .or. ii.gt.20)goto 50000
        if(i.eq.1)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
        if(i.eq.2)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
        if(i.eq.3)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
50000   enddo
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine oscar(iii)
c       record history of spatial and momentum coordinates due to
c        OSC1999A
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000,KSZ1=30)
        common/PYJETS/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
c        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        write(34,*)'history orderring=',iii
c       write(34,*)'event orderring=',iii
        do j=1,nsa
        id=ksa(j,2)
        x=vsa(j,1)   ! -coor(1)
        y=vsa(j,2)   ! -coor(2)
        z=vsa(j,3)   ! -coor(3)
        t=vsa(j,4)
        px=psa(j,1)
        py=psa(j,2)
        pz=psa(j,3)
        e=psa(j,4)
        am=psa(j,5)
        write(34,100)j,id,px,py,pz,e,am,x,y,z,t
        enddo
100     format(i5,1x,i5,4f8.3,1x,f7.3,1x,4f9.3)   ! 170705 'm' f4.3 -> f7.3
        return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sgam(nn)   ! 240209
c       print particle list 'sgam' and sum of momentum and energy
        parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        dimension peo(4)
        call psum(pgam,1,ngam,peo)
        ch1=0.
        do i1=1,nn
        kf=kgam(i1,2)
        if(kf.eq.22 .or. kf.eq.44 .or. kf.eq.55 .or. kf.eq.66)goto 100
        ch1=ch1+pychge(kf)
100     enddo
        write(22,*)ch1/3.,peo   !
        do i=1,nn
        write(22,*)i,kgam(i,1),kgam(i,2),(pgam(i,j),j=1,4)
        enddo
        return
        end


                                                              
*******************************************************************************
c This code calculates the nuclear overlap functions which are needed 
c to scale pp to pA and AB.
c Units are fm.
c TAB is overlap function in 1/fm^2; divide the value by 10 to get it in mb^-1
c scal1 is sigma(A1,A2)/sigma(p,p)
c scal2 is n(A1,A2)/sigma(p,p)
c The calculation is done in steps of 0.1 fm. 
c The cutoff R and b are 10 and 20 fm, respectively.
c D.Miskowiec 1997, updated in 2001. 
*******************************************************************************
      subroutine overlap(A1,A2,rnp,rnt,sigma_NN,kjp23,denflag,density0,
     c  nshot)   ! 020511
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      integer denflag,A1,A2
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200),
     c  part1(200),part2(200)   ! 020511
      rmax=10.0
c        write(9,*)'sigma_NN,kjp23,denflag,density0,nshot=',
c     c  sigma_NN,kjp23,denflag,density0,nshot
c calculate thickness functions for A1 and A2
      do i=1,200
        bb(i)=dble(i)/10.0-0.05   !longhy:change "i" type, b range is 
c        [0.05,19.95 ]
        if(kjp23.eq.2)then
        TA1(i)=TA(A1,bb(i),denflag,density0,nshot)
        TA2(i)=TA(A2,bb(i),denflag,density0,nshot)
c        if(i.eq.36)write(9,*)'kjp23,b,density0,ta1,ta2=',kjp23,bb(i),
c     c  density0,ta1(i),ta2(i)
        endif
      enddo
c Calculate overlap function for A1+A2 collision using densities
        if(kjp23.eq.2)then   ! 190511
      do i=1,200
        bbb=bb(i)
        TA1A2(i)=TAB(A1,A2,bbb,denflag,density0,nshot)
      enddo
        endif
c Calculate npart using thickness functions
        do i=1,200
        bbb=bb(i)
        if(kjp23.eq.2)then
        call PART(A1,A2,bbb,sigma_NN,denflag,density0,nshot,
     c  part1(i),part2(i))
c        if(i.eq.36)write(9,*)'i,bbb,pir,tir=',i,bbb,part1(i),part2(i)
        endif
        if(kjp23.eq.1)then
        call irpt1(bbb,1,a1,a2,rou,rnp,rnt,1,0.d0,pir,tir,1.d0,
     c   1.d0,1.d0,1.d0)
        part1(i)=pir
        part2(i)=tir
c        if(i.eq.36)write(9,*)'i,bbb,pir,tir=',i,bbb,pir,tir
        endif
        enddo
c Print data
c-----------
      if(kjp23.eq.2)write(9,901)'i','b','TA','TB','TAB','Apart','Bpart'
        if(kjp23.eq.1)write(9,902)'i','b','Apart','Bpart'
      do i=1,200
        if(kjp23.eq.2)write(9,905)i,bb(i),TA1(i),TA2(i),TA1A2(i),
     c  part1(i),part2(i)
        if(kjp23.eq.1)write(9,906)i,bb(i),part1(i),part2(i)
      enddo
 901  format(2a6,5a10)
 902    format(2a6,2a10)
 905  format(i5,2x,f6.2,5(f10.3))
 906  format(i5,2x,f6.2,2(f10.3,1x))
      return
      end



*******************************************************************************
      function TA(A,b,denflag,density0,nshot)   ! 020511
c       one dimension trapezoid integral
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP      
      integer denflag,A
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200),
     c  part1(200),part2(200)   ! 020511
c      write(9,*)'ta a,b,denflag,density0,nshot=',
c     c  a,b,denflag,density0,nshot
      TA=0.0
cc      do j=1,nshot
cc        z=pyr(1)*rmax
cc        TA=TA+density(A,b,0.0d0,z,denflag,density0)
cc      enddo
cc      ta=TA/nshot*2*rmax   ! 2: left and right symmetry, sa
        nshot2=2*nshot
        delta=rmax/nshot
        do j=0,nshot2
        if(j.le.nshot)then
        z=-dfloat(nshot-j)*delta
        else
        z=dfloat(j-nshot)*delta
        endif
        density1=density(A,b,0.0d0,z,denflag,density0)
        if(j.eq.0 .or. j.eq.nshot2)density1=0.5*density1
        TA=TA+density1
        enddo
        TA=TA*delta
c         write(9,*)'out of ta ta=',ta
      return
      end



*******************************************************************************
c direct calculation of TAB (not via TA*TB)
      function TAB(A1,A2,b,denflag,density0,nshot)
c       four dimensions trapezoid integral
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      integer denflag,a1,a2
      common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200),
     c  part1(200),part2(200)
      TAB=0.0
cc      do j=1,nshot
cc        x=(pyr(1)-0.5)*2*rmax
cc        y=(pyr(1)-0.5)*2*rmax
cc        za=(pyr(1)-0.5)*2*rmax
cc        zb=(pyr(1)-0.5)*2*rmax
        nshot2=2*nshot
        delta=rmax/nshot
        do 100 i=0,nshot2
        if(i.le.nshot)then
        x=-dfloat(nshot-i)*delta
        else
        x=dfloat(i-nshot)*delta
        endif
        taby=0.0
        do 200 j=0,nshot2
        if(j.le.nshot)then
        y=-dfloat(nshot-j)*delta
        else
        y=dfloat(j-nshot)*delta
        endif
        tabz=0.0
        do 300 k=0,nshot2
        if(k.le.nshot)then
        z=-dfloat(nshot-k)*delta
        else
        z=dfloat(k-nshot)*delta
        endif
        tabz1=0.0
        do 400 k1=0,nshot2
        if(k1.le.nshot)then
        z1=-dfloat(nshot-k1)*delta
        else
        z1=dfloat(k1-nshot)*delta
        endif
        o=density(A1,x+b,y,z,denflag,density0)*
     c    density(A2,x,y,z1,denflag,density0)
cc        TAB=TAB+o
cc      enddo
cc      TAB=TAB/nshot*(2*rmax)**4
        if(k1.eq.0 .or. k1.eq.nshot2)o=0.5*o
        tabz1=tabz1+o
400     enddo
        if(k.eq.0 .or. k.eq.nshot2)tabz1=0.5*tabz1
        tabz=tabz+tabz1
300     enddo
        if(j.eq.0 .or. j.eq.nshot2)tabz=0.5*tabz
        taby=taby+tabz
200     enddo
        if(i.eq.0 .or. i.eq.nshot2)taby=0.5*taby
        tab=tab+taby
100     enddo
        tab=tab*delta**4
      return
      end



*******************************************************************************
c denflag=1 - sharp sphere
c denflag=2 - Woods-Saxon
      function density(A,x,y,z,denflag,density0)   ! 020511
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP      
      integer denflag,A
      r=sqrt(x*x+y*y+z*z)      
      if (denflag.eq.1) then 
        RA=(A/4.0*3.0/3.1415/density0)**0.3333333
        density=0.0
        if (r.le.RA) density=density0
      elseif (denflag.eq.2) then 
        RA=1.12*A**0.333333-0.86/A**0.333333   ! 020511 070613 recovered
c       RA=1.19*A**0.333333-1.61/A**0.333333
c       density0=3./4.*A/3.1416/RA**3/(1+3.1416**2*0.54**2/RA**2)
        DR=0.54
        density=density0/(1+exp((r-RA)/dr))
      else 
        write(*,*)'wrong density profile flag'
        stop
      endif
      return
      end



*******************************************************************************
c npart calculation via TA and TB
c part1 and part2 are the respective numbers of participants from A1 and A2
c sigma_NN is in millibarn and TA1 and TA2 are in fm^-2 
c sigma_NN*TA/10 is mean number of NN collisions.
c Vector b points from A to B.
      subroutine PART(A1,A2,b,sigma_NN,denflag,density0,nshot,art1,art2)! 020511
c       two dimensions trapezoid integral
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP      
      integer denflag,A1,A2
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200),
     c  part1(200),part2(200)   ! 020511
c      write(9,*)'part A1,A2,b,sigma_NN=',A1,A2,b,sigma_NN        
      art1=0.0
      art2=0.0
cc      do j=1,nshot
cc        x=(pyr(1)-0.5)*2*rmax
cc        y=(pyr(1)-0.5)*2*rmax
        nshot2=2*nshot
        delta=rmax/nshot
        do 100 i=0,nshot2
        if(i.le.nshot)then
        x=-dfloat(nshot-i)*delta
        else
        x=dfloat(i-nshot)*delta
        endif
        art1y=0.0
        art2y=0.0
        do 200 j=0,nshot2
        if(j.le.nshot)then
        y=-dfloat(nshot-j)*delta
        else
        y=dfloat(j-nshot)*delta
        endif
        b1=sqrt(x**2+y**2)
        b2=sqrt((x-b)**2+y**2) 
        ib1=int(b1*10+1.0)
        ib2=int(b2*10+1.0)
        ib1=min(ib1,200)
        ib2=min(ib2,200)
        if(j.eq.0 .or. j.eq.nshot2)then
        art1y=art1y+0.5*TA1(ib1)*(1-exp(-sigma_NN*TA2(ib2)/10))
        art2y=art2y+0.5*TA2(ib2)*(1-exp(-sigma_NN*TA1(ib1)/10))
        else
        art1y=art1y+TA1(ib1)*(1-exp(-sigma_NN*TA2(ib2)/10))
        art2y=art2y+TA2(ib2)*(1-exp(-sigma_NN*TA1(ib1)/10))
        endif
200     enddo
        if(i.eq.0 .or. i.eq.nshot2)then
        art1y=0.5*art1y
        art2y=0.5*art2y
        endif
        art1=art1+art1y
        art2=art2+art2y
100     enddo
c        art1=art1+TA1*(1-exp(-sigma_NN*TA2/10))
c        art2=art2+TA2*(1-exp(-sigma_NN*TA1/10))
c       art1=art1+TA1(ib1)*(1-(1-sigma_NN*TA2(ib2)/10/A2)**A2)
c       art2=art2+TA2(ib2)*(1-(1-sigma_NN*TA1(ib1)/10/A1)**A1)
cc      enddo
c      write(9,*)'part1,part2=',part1,part2
cc      art1=art1/nshot*(2*rmax)**2
cc      art2=art2/nshot*(2*rmax)**2
      art1=art1*delta**2
      art2=art2*delta**2
c      write(9,*)'out of part part1,part2=',part1,part2
      return
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shanul(x,q2,rag,raq)   ! 181213
c       calculate nuclear ratio R^A_i=f_{i/A}(x,Q2)/f_i(x,Q2) according
c        to Xin-Nian Wang's paper (PL, B527(2002)85), multiply it to
c        the parton distribution function in pythia, resulted parton
c        distribution function is including nuclear shadowing effect
c       it was proved in Eur. Phys. J. C9(1999)61 that nuclear ratio does
c        not depend strongly on the choice for the parton distribution
c        function in nucleon f_i(x,Q2)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        sq=0.1
        sg=0.26
        x35=x**0.35
        xsqr=sqrt(x)
        x2=x*x
        x3=x2*x
        xx=x3-1.2*x2+0.21*x
        a=nat   ! nap originally 181213
c       what is the definition of "a" for asymmetry reaction system ?
        a13=a**0.3333
        aa=(a13-1)**0.6
        coa=log(a)
        coa16=coa**0.16666
        bbq=1.-3.5*xsqr
        bbg=1.-1.5*x35
        eq=exp(-x2/0.01)
        eg=exp(-x2/0.004)
c       raq=a*(1.+1.19*coa16*xx-sq*aa*bbq*eq)
c       rag=a*(1.+1.19*coa16*xx-sg*aa*bbg*eg)
	raq=1.+1.19*coa16*xx-sq*aa*bbq*eq
        rag=1.+1.19*coa16*xx-sg*aa*bbg*eg
c	write(9,*)'rag,raq=',rag,raq   ! 181213	
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE IRPT1(BB1,IB1,IAP1,IAT1,ROU1,RP1,RT1,ISETA1,CSETA1,
     $ PIR1,TIR1,AP1,BP1,AT1,BT1)
C     'IRPT' IS A PROGRAM FOR COMPUTING THREE DIMENSIONS INTEGRAL BY GAUSSIAN 
C      TYPE QUADRATURE RULE.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      REAL*4 BB,ROU,RP,RT,CSETA,PIR,TIR,AP,BP,AT,BT
        BB=BB1
        IB=IB1
        IAP=IAP1
        IAT=IAT1
        ROU=ROU1
        RP=RP1
        RT=RT1
        ISETA=ISETA1
        CSETA=CSETA1
        PIR=PIR1
        TIR=TIR1
        AP=AP1
        BP=BP1
        AT=AT1
        BT=BT1
        CALL IRPT(BB,IB,IAP,IAT,ROU,RP,RT,ISETA,CSETA,PIR,TIR,
     $ AP,BP,AT,BT)
        PIR1=PIR
        TIR1=TIR
      RETURN  
      END 



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE IRPT(BB,IB,IAP,IAT,ROU,RP,RT,ISETA,CSETA,PIR,TIR,
     $ AP,BP,AT,BT)
C     'IRPT' IS A PROGRAM FOR COMPUTING THREE DIMENSIONS INTEGRAL BY GAUSSIAN 
C      TYPE QUADRATURE RULE.
C     HERE IT'S USED TO CALCULATE THE INTERSECT REGION OF PROJECTILE AND TARGET.
C     IT IS WRITTEN BY SA BEN-HAO ON NOV. 6,1988.
C     IAP AND IAT = MASS NUMBER OF PROJECTILE AND TARGET .
C     RP AND RT = RADII OF PROJECTILE AND TARGET (IN SPHERE CASE).
C     AP,BP AND AT,BT = THE LENGTH OF SINGLE, DOUBLE AXES OF PROJECTILE
C      AND TARGET (IN SPHEROID CASE).
C     AP=BP IF PROJECTILE IS A SPHERE,THE SAME FOR TARGET.
C     IB = NUMBER OF IMPACT PARAMETERS NEEDED TO CALCULATE IN SINGLE RUN.
C     BB = CURRENT IMPACT PARAMETER
C     ISETA = NUMBER OF SETA (THE COSIN OF ANGLE OF SINGLE AXIS OF TARGET 
C      RELATIVE TO THE INCIDENT AXIS) NEEDED TO CALCULATE IN SINGLE RUN.
C      FOR THE CASE OF SPHERE TARGET ISETA SHOULE BE EQUAL TO 1.
C     CSETA = 0. FOR THE CASE OF PROJECTILE AND TARGET ARE SPHERES. 
C     CSETA = 1. FOR THE CASE OF TARGET IS SPHEROID BUT PROJECTILE IS SPHERE.
C     CSETA = 2. FOR THE CASE OF TARGET AND PROJECTILE ARE SPHEROIDS.
C     IRP OR PIR = THE NUMBER OF NUCLEONS IN INTERSECT REGION FROM PROJECTILE.
C     IRT OR TIR = THE NUMBER OF NUCLEONS IN INTERSECT REGION FROM TARGET.
C     IRPT = IRP*IRT. 
C     APIR OR AIRP = THE AVERAGE IRP OVER SETA. 
C     ATIR OR AIRT = THE AVERAGE IRT OVER SETA. 
C     RLARGE = THE LARGER ONE BETWEEN AP AND AT.
C     D1,D2 AND D3 = LOWER LIMITS IN X,Y AND Z DIRECTIONS.
C     U1,U2 AND U3 = UPPER LIMITS IN X,Y AND Z DIRECTIONS.
C     N1,N2,N3 = THE NUMBER OF NODES IN X,Y AND Z DIRECTIONS. 
C     X11,X22 AND X33(200) = ARRAIES OF NODES IN X,Y AND Z DIRECTIONS.
C     X1,X2 AND X3(200) = ARRAIES OF TRANSFORMED NODES IN X,Y AND Z DIRECTIONS. 
C     W1,W2 AND W3(200) = ARRAIES OF WEIGHTS IN X,Y AND Z DIRECTIONS. 
      COMMON/COM1/SINS,COSS
      COMMON/AA18/DNI(10),DPI(10),EDI(10),BMIN,BMAX
      DIMENSION ENDPTS(2),BSCR(200) 
      DIMENSION X1(200),X2(200),X3(200),SUM(2),W1(200),W2(200),W3(200)
      DIMENSION X11(200),X22(200),X33(200)
      N1=100
      N2=N1
      N3=N1
      IF(ISETA.EQ.1)GOTO 210
      RLARGE=AP
      IF(AT.GT.AP)RLARGE=AT
      GOTO 211  
210   RLARGE=RP
      IF(RT.GT.RP)RLARGE=RT
211   D1=-RLARGE
      D2=D1
      D3=D1
      U1=RLARGE
      U2=U1
      U3=U1
CC      COSS=SETA
CC      SINS=SQRT(1-SETA*SETA)
C     THE SYMMETRY AROUND 3.1416/2 IS TAKEN INTO ACCOUNT.
      SSS=1./FLOAT(ISETA) 
      BBB=(BMAX-BMIN)/FLOAT(IB)
      CALL GAUSSQ(1,N1,0.,0.,0,ENDPTS,BSCR,X11,W1)
      CALL GAUSSQ(1,N2,0.,0.,0,ENDPTS,BSCR,X22,W2)
      CALL GAUSSQ(1,N3,0.,0.,0,ENDPTS,BSCR,X33,W3)
C      WRITE(9,*)'B=',BB 
C      WRITE(9,*)'AP,BP,AT,BT=',AP,BP,AT,BT
C      WRITE(9,*)'BMIN,BMAX,IB,ISETA,CSETA=',BMIN,BMAX,IB,ISETA,CSETA
C      WRITE(9,*)'N1,N2,N3=',N1,N2,N3
C      WRITE(9,*)'D1,D2,D3=',D1,D2,D3
C      WRITE(9,*)'U1,U2,U3=',U1,U2,U3
C      WRITE(9,110)(X11(I),I=1,N1)
C      WRITE(9,111)(W1(I),I=1,N1)
C 110  FORMAT(1X,'X1=',10E12.4/)
C 111  FORMAT(1X,'W1=',10E12.4/)
C     CALCULATE THE TRANSFORMED NODES.
C        F((U+D)/2+(U-D)/2*X).
C     X IS NODE.
      DO 101 I=1,N1 
 101  X1(I)=0.5*(U1+D1+(U1-D1)*X11(I))
      DO 102 I=1,N2 
 102  X2(I)=0.5*(U2+D2+(U2-D2)*X22(I))
      DO 103 I=1,N3 
 103  X3(I)=0.5*(U3+D3+(U3-D3)*X33(I))
C     CALCULATE THE SUMS (INTEGRALS). 
C       N
C      SUM  W  * F(X ) .
C      I=1   I      I
      IF(IB.EQ.1)GOTO 201
      IIB=1
202   BB=BBB*IIB+BMIN 
C      WRITE(9,*)'B=',BB
201   PIR=0.
      TIR=0.
      IF(ISETA.EQ.1)GOTO 203
      JS=1
208   SETA=SSS*JS 
      COSS=SETA
      SINS=SQRT(1-SETA*SETA)
203   ROUP=3./4./3.1416*IAP/(RP*RP*RP)
      ROUT=3./4./3.1416*IAT/(RT*RT*RT)
      DO 100 II=1,2 
	ROU=ROUP
	IF(II.EQ.2)ROU=ROUT
      SUM(II)=0.
      DO 105 I1=1,N1
      XX1=X1(I1)
      SUM2=0.
      DO 106 I2=1,N2
      XX2=X2(I2)
      SUM3=0.
      DO 107 I3=1,N3
      XX3=X3(I3)
 107  SUM3=SUM3+W3(I3)*
     $ FUNC(XX1,XX2,XX3,II,IAP,IAT,RP,RT,ROU,BB,AP,BP,AT,BT,CSETA)
C
C
C
 106  SUM2=SUM2+W2(I2)*SUM3
 105  SUM(II)=SUM(II)+W1(I1)*SUM2 
      SUM(II)=SUM(II)*(U1-D1)/2.*(U2-D2)/2.*(U3-D3)/2.*ROU
100   CONTINUE
      SUMM1=SUM(1)
      SUMM2=SUM(2)
      PIR=PIR+SUMM1 
      TIR=TIR+SUMM2 
C      WRITE(9,*)'PIR,TIR,PTIR=',SUMM1,SUMM2,SUMM1*SUMM2
      IF(ISETA.EQ.1)GOTO 204
      IF(JS.GE.ISETA)GOTO 209 
      JS=JS+1
      GOTO 208
209   APIR=PIR/ISETA
      ATIR=TIR/ISETA
      APTIR=APIR*ATIR 
C      WRITE(9,205)
C205   FORMAT(/12X,'B',13X,'AIRP',13X,'AIRT',15X,'AIRPT')
C      WRITE(9,206)BB,APIR,ATIR,APTIR
C206   FORMAT(/3X,E12.4,4X,E12.4,4X,E12.4,10X,E12.4) 
204   IF(IB.EQ.1)GOTO 207 
      IF(IIB.GE.IB)GOTO 207 
      IIB=IIB+1
      GOTO 202
207   RETURN
      END 




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SUBROUTINE FOR THE INTEGRANT. 
      FUNCTION FUNC(X1,X2,X3,II,IAP,IAT,RP,RT,ROU,BB,AP,BP,AT,BT,CSETA)
      COMMON/COM1/SINS,COSS
      X12=X1*X1 
      X22=X2*X2
      X32=X3*X3
      BY=BB-X2
      BY2=BY*BY
      IF(CSETA.NE.0.)GOTO 400
      IF(II.EQ.2)GOTO 200
C	CALCULATE THE VOLUME OF INTERSECTION FOR PROJECTILE
      ARG1=SQRT(X12+BY2+X32)
      F1=0.
      IF(RP.GE.ARG1)F1=1.
      ARG2=SQRT(X12+X22)
      F2=0.
      IF(RT.GE.ARG2)F2=1.
      GOTO 300
C	CALCULATE THE VOLUME OF INTERSECTION FOR TARGET
200   ARG1=SQRT(X12+X22+X32)
      F1=0.
      IF(RT.GE.ARG1)F1=1.
      ARG2=SQRT(X12+BY2)
      F2=0.
      IF(RP.GE.ARG2)F2=1.
300   FUNC=F1*F2
      GOTO 500
400   AT2=AT*AT 
      BT2=BT*BT
      IF(II.EQ.2)GOTO 600
      ARG1=SQRT(X12+BY2+X32)
      F1=0.
      IF(RP.GE.ARG1)F1=1.
      ARG2=SQRT(X12+X22)
      RR=AT*SINS
      RR=AMAX1(BT,RR)
C     WRITE(6,*)'RR=',RR
      F2=0.
      IF(RR.GE.ARG2)F2=1.
      GOTO 700
600   YS=X2*COSS-X3*SINS
      ZS=X2*SINS+X3*COSS
      ARG1=(X12+YS*YS)/BT2+ZS*ZS/AT2
      F1=0.
      IF(1.GE.ARG1)F1=1.
      ARG2=SQRT(X12+BY2)
      F2=0.
      IF(RP.GE.ARG2)F2=1.
700   FUNC=F1*F2
500   RETURN
      END




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                                                                               
C        TITLE:  GAUSSQ                                                         
C                                                                               
C                                1/20/75                                        
C                                                                               
C           THIS SET OF ROUTINES COMPUTES THE NODES T(J) AND WEIGHTS            
C        W(J) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED              
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE                  
C                                                                               
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX                          
C                                                                               
C                              N                                                
C        BY                   SUM W  F(T )                                      
C                             J=1  J    J                                       
C                                                                               
C        (NOTE W(X) AND W(J) HAVE NO CONNECTION WITH EACH OTHER.)               
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT                   
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE                              
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY        
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT                  
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.                    
C                                                                               
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF               
C        ORTHOGONAL POLYNOMIALS.  THE NODES T(J) ARE JUST THE ZEROES            
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.                                  
C                                                                               
C     INPUT PARAMETERS (ALL REAL NUMBERS ARE IN DOUBLE PRECISION)               
C                                                                               
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF                 
C                 QUADRATURE RULE:                                              
C                                                                               
C        KIND = 1:  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)                    
C        KIND = 2:  CHEBYSHEV QUADRATURE OF THE FIRST KIND                      
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)                          
C        KIND = 3:  CHEBYSHEV QUADRATURE OF THE SECOND KIND                     
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)                             
C        KIND = 4:  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON                     
C                   (-INFINITY, +INFINITY)                                      
C        KIND = 5:  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**            
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.                       
C                   NOTE: KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.              
C        KIND = 6:  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*            
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1                   
C                                                                               
C        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE             
C        ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-          
C        BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE--        
C                 (OTHERWISE USE 0.D0)                                          
C        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-           
C                 POINT (OR BOTH OF THE INTERVAL IS REQUIRED TO BE A            
C                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO             
C                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED                
C                 ENDPOINTS (1 OR 2).                                           
C        ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF               
C                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.                        
C        B        REAL SCRATCH ARRAY OF LENGTH N                                
C                                                                               
C     OUTPUT PARAMETERS (BOTH DOUBLE PRECISION ARRAYS OF LENGTH N)              
C                                                                               
C        T        WILL CONTAIN THE DESIRED NODES.                               
C        W        WILL CONTAIN THE DESIRED WEIGHTS W(J).                        
C                                                                               
C     SUBROUTINES REQUIRED                                                      
C                                                                               
C        SLVE, CLASS, AND IMTQL2 ARE PROVIDED.  UNDERFLOW MAY SOMETIMES         
C        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE              
C        TURNED OFF.                                                            
C                                                                               
C     ACCURACY                                                                  
C                                                                               
C        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,          
C        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP           
C        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,             
C        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT         
C        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR           
C        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME           
C        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,             
C        COMPLETELY HARMLESS.                                                   
C                                                                               
C     METHOD                                                                    
C                                                                               
C           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION              
C        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE                
C        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE                     
C        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH                 
C        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF           
C        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,                
C        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A         
C        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.         
C        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF           
C        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.                         
C                                                                               
C     REFERENCES                                                                
C                                                                               
C        1.  GOLUB, G. H., AND WELSCH, J. H., "CALCULATION OF GAUSSIAN          
C            QUADRATURE RULES," MATHEMATICS OF COMPUTATION 23 (APRIL,           
C            1969), PP. 221-230.                                                
C        2.  GOLUB, G. H., "SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,"          
C            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).             
C        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE-        
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C                                                                               
      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)            
C                                                                               
      REAL B(N), T(N), W(N), ENDPTS(2), MUZERO, T1,                             
     X GAM, SLVE, DSQRT, ALPHA, BETA                                            
C                                                                               
      CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)                           
C                                                                               
C           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.              
C           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY               
C           B THE OFF-DIAGONAL ELEMENTS.                                        
C           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2                  
C           SUBMATRIX.                                                          
C                                                                               
      IF (KPTS.EQ.0)  GO TO 100                                                 
      IF (KPTS.EQ.2)  GO TO  50                                                 
C                                                                               
C           IF KPTS=1, ONLY T(N) MUST BE CHANGED                                
C                                                                               
      T(N) = SLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)                     
      GO TO 100                                                                 
C                                                                               
C           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED                       
C                                                                               
   50 GAM = SLVE(ENDPTS(1), N, T, B)                                            
      T1 = ((ENDPTS(1) - ENDPTS(2))/(SLVE(ENDPTS(2), N, T, B) - GAM))           
      B(N-1) = SQRT(T1)                                                         
      T(N) = ENDPTS(1) + GAM*T1                                                 
C                                                                               
C           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1        
C           AND THUS THE VALUE OF B(N) IS ARBITRARY.                            
C           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL            
C           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.                       
C           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING            
C                                                                               
  100 W(1) = 1.0E0                                                              
      DO 105 I = 2, N                                                           
  105    W(I) = 0.0E0                                                           
C                                                                               
      CALL IMTQL2 (N, T, B, W, IERR)                                            
      DO 110 I = 1, N                                                           
  110    W(I) = MUZERO * W(I) * W(I)                                            
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      FUNCTION SLVE(SHIFT, N, A, B)                                             
C                                                                               
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE                    
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION                    
C                                                                               
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,                              
C                                                                               
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN                    
C       THE N-TH POSITION.                                                      
C                                                                               
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL                   
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION               
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER           
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.            
C                                                                               
C                                                                               
      REAL A(N), B(N), ALPHA                                                    
C                                                                               
      ALPHA = A(1) - SHIFT                                                      
      NM1 = N - 1                                                               
      DO 10 I = 2, NM1                                                          
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA                                 
      SLVE = 1.0E0/ALPHA                                                        
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)                      
C                                                                               
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE          
C        RECURRENCE RELATION                                                    
C                                                                               
C             B P (X) = (X - A ) P   (X) - B   P   (X)                          
C              J J            J   J-1       J-1 J-2                             
C                                                                               
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS,         
C        AND THE ZERO-TH MOMENT                                                 
C                                                                               
C             MUZERO = INTEGRAL W(X) DX                                         
C                                                                               
C        OF THE GIVEN POLYNOMIAL'S WEIGHT FUNCTION W(X).  SINCE THE             
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS             
C        GUARANTEED TO BE SYMMETRIC.                                            
C                                                                               
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND             
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR            
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS               
C        REQUIRE THE GAMMA FUNCTION.                                            
C                                                                               
C     ..................................................................        
C                                                                               
      REAL A(N), B(N), MUZERO, ALPHA, BETA                                      
      REAL ABI, A2B2, GAMMA, PI, DSQRT, AB                                     
      DATA PI / 3.141592653589793E0/                                            
C                                                                               
      NM1 = N - 1                                                               
      GO TO (10, 20, 30, 40, 50, 60), KIND                                      
C                                                                               
C              KIND = 1:  LEGENDRE POLYNOMIALS P(X)                             
C              ON (-1, +1), W(X) = 1.                                           
C                                                                               
   10 MUZERO = 2.0E0                                                            
      DO 11 I = 1, NM1                                                          
         A(I) = 0.0E0                                                           
         ABI = I                                                                
   11    B(I) = ABI/SQRT(4*ABI*ABI - 1.0E0)                                     
      A(N) = 0.0E0                                                              
      RETURN                                                                    
C                                                                               
C              KIND = 2:  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)          
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)                            
C                                                                               
   20 MUZERO = PI                                                               
      DO 21 I = 1, NM1                                                          
         A(I) = 0.0E0                                                           
   21    B(I) = 0.5E0                                                           
      B(1) = SQRT(0.5E0)                                                        
      A(N) = 0.0E0                                                              
      RETURN                                                                    
C                                                                               
C              KIND = 3:  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)         
C              ON (-1, +1), W(X) = SQRT(1 - X*X)                                
C                                                                               
   30 MUZERO = PI/2.0E0                                                         
      DO 31 I = 1, NM1                                                          
         A(I) = 0.0E0                                                           
   31    B(I) = 0.5E0                                                           
      A(N) = 0.0E0                                                              
      RETURN                                                                    
C                                                                               
C              KIND = 4:  HERMITE POLYNOMIALS H(X) ON (-INFINITY,               
C              +INFINITY), W(X) = EXP(-X**2)                                    
C                                                                               
   40 MUZERO = SQRT(PI)                                                         
      DO 41 I = 1, NM1                                                          
         A(I) = 0.0E0                                                           
   41    B(I) = SQRT(I/2.0E0)                                                   
      A(N) = 0.0E0                                                              
      RETURN                                                                    
C                                                                               
C              KIND = 5:  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON               
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND           
C              BETA GREATER THAN -1                                             
C                                                                               
   50 AB = ALPHA + BETA                                                         
      ABI = 2.0E0 + AB                                                          
      MUZERO = 2.0E0 ** (AB + 1.0E0) * GAMMA(ALPHA + 1.0E0) * GAMMA(          
     X BETA + 1.0E0) / GAMMA(ABI)                                              
      A(1) = (BETA - ALPHA)/ABI                                                 
      B(1) = SQRT(4.0E0*(1.0E0 + ALPHA)*(1.0E0 + BETA)/((ABI + 1.0E0)*          
     1  ABI*ABI))                                                               
      A2B2 = BETA*BETA - ALPHA*ALPHA                                            
      DO 51 I = 2, NM1                                                          
         ABI = 2.0E0*I + AB                                                     
         A(I) = A2B2/((ABI - 2.0E0)*ABI)                                        
   51    B(I) = SQRT (4.0E0*I*(I + ALPHA)*(I + BETA)*(I + AB)/                  
     1   ((ABI*ABI - 1)*ABI*ABI))                                               
      ABI = 2.0E0*N + AB                                                        
      A(N) = A2B2/((ABI - 2.0E0)*ABI)                                           
      RETURN                                                                    
C                                                                               
C              KIND = 6:  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON                   
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER         
C              THAN -1.                                                         
C                                                                               
   60 MUZERO = GAMMA(ALPHA + 1.0E0)                                            
      DO 61 I = 1, NM1                                                          
         A(I) = 2.0E0*I - 1.0E0 + ALPHA                                         
   61    B(I) = SQRT(I*(I + ALPHA))                                             
      A(N) = 2.0E0*N - 1 + ALPHA                                                
      RETURN                                                                    
      END                                                                       
C                                                                               
C     ------------------------------------------------------------------        
C                                                                               
C     TITLE:  IMTQL2                                                            
C                                                                               
C                                                                               
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,           
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,                     
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.                      
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).           
C     THIS IS A MODIFIED VERSION OF THE 'EISPACK' ROUTINE IMTQL2.               
C                                                                               
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE         
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL         
C     METHOD.                                                                   
C                                                                               
C     ON INPUT:                                                                 
C                                                                               
C        N IS THE ORDER OF THE MATRIX;                                          
C                                                                               
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;                  
C                                                                               
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX                
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;                      
C                                                                               
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.                       
C                                                                               
C      ON OUTPUT:                                                               
C                                                                               
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN                  
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT                  
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;                             
C                                                                               
C        E HAS BEEN DESTROYED;                                                  
C                                                                               
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS        
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS            
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED         
C          EIGENVALUES;                                                         
C                                                                               
C        IERR IS SET TO                                                         
C          ZERO       FOR NORMAL RETURN,                                        
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN                       
C                     DETERMINED AFTER 30 ITERATIONS.                           
C                                                                               
      SUBROUTINE IMTQL2(N, D, E, Z, IERR)                                       
C                                                                               
      REAL D(N), E(N), Z(N), B, C, F, G, P, R, S, MACHEP                        
      REAL DSQRT, DABS, DSIGN                                                   
C                                                                               
      MACHEP=2.2E-17                                                            
      IERR=0                                                                    
      IF (N.EQ.1)GOTO 1001                                                      
      E(N) = 0.0E0                                                              
      DO 240 L = 1, N                                                           
         J = 0                                                                  
  105    DO 110 M = L, N                                                        
      IF(M.EQ.N)GOTO 120                                                        
      IF(ABS(E(M)).LE.MACHEP*(ABS(D(M))+ABS(D(M+1))))GOTO 120                   
110   CONTINUE                                                                  
C                                                                               
120   P=D(L)                                                                    
      IF(M.EQ.L)GOTO 240                                                        
      IF(J.EQ.30)GOTO 1000                                                      
      J=J+1                                                                     
C                                                                               
      G=(D(L+1)-P)/(2.0E0*E(L))                                                 
      R=SQRT(G*G+1.0E0)                                                         
      G=D(M)-P+E(L)/(G+SIGN(R,G))                                               
      S=1.0E0                                                                   
      C=1.0E0                                                                   
      P=0.0E0                                                                   
      MML=M-L                                                                   
C                                                                               
C                                                                               
      DO 200 II=1,MML                                                           
            I = M - II                                                          
            F = S * E(I)                                                        
            B = C * E(I)                                                        
            IF (ABS(F) .LT. ABS(G)) GO TO 150                                   
            C = G / F                                                           
            R = SQRT(C*C+1.0E0)                                                 
      E(I+1)=F*R                                                                
      S=1.0E0/R                                                                 
      C=C*S                                                                     
      GOTO 160                                                                  
 150  S=F/G                                                                     
      R=SQRT(S*S+1.0E0)                                                         
      E(I+1)=G*R                                                                
      C=1.0E0/R                                                                 
      S=S*C                                                                     
 160  G=D(I+1)-P                                                                
            R = (D(I) - G) * S + 2.0E0 * C * B                                  
            P = S * R                                                           
            D(I+1) = G + P                                                      
            G = C * R - B                                                       
C                                                                               
            F = Z(I+1)                                                          
      Z(I+1)=S*Z(I)+C*F                                                         
 200  Z(I)=C*Z(I)-S*F                                                           
C                                                                               
      D(L)=D(L)-P                                                               
      E(L)=G                                                                    
      E(M)=0.0E0                                                                
         GO TO 105                                                              
  240 CONTINUE                                                                  
C                                                                               
C                                                                               
      DO 300 II = 2, N                                                          
      I=II-1                                                                    
         K = I                                                                  
         P = D(I)                                                               
C                                                                               
         DO 260 J = II, N                                                       
            IF (D(J) .GE. P) GO TO 260                                          
            K = J                                                               
            P = D(J)                                                            
  260    CONTINUE                                                               
C                                                                               
         IF (K .EQ. I) GO TO 300                                                
         D(K) = D(I)                                                            
         D(I) = P                                                               
         P = Z(I)                                                               
         Z(I) = Z(K)                                                            
         Z(K) = P                                                               
  300 CONTINUE                                                                  
C                                                                               
      GO TO 1001                                                                
C                EIGENVALUE           AFTER 30 ITERATIONS                       
 1000 IERR = L                                                                  
 1001 RETURN                                                                    
C                                                                               
      END                                                                       
      FUNCTION GAMMA(X)                                                        
      C=1.0E0                                                                   
      Y=X                                                                       
      IF(Y-2.0E0)1,2,2                                                          
 2    IF(Y-3.0E0)3,3,4                                                          
 4    IF(Y-10.0E0)5,5,6                                                         
 5    Y=Y-1.0E0                                                                 
      C=C*Y                                                                     
      IF(Y-3.0E0)3,3,5                                                          
 1    C=Y*C                                                                     
      Y=Y+1.0E0                                                                 
      IF(Y-2.0E0)1,9,9                                                          
    9 C=1.0E0/C                                                                 
    3 B=Y-2.0E0                                                                 
      G=(((((((-.5113262726698E-6*B+.51063592072582E-5)*B                       
     1  -.248410053848712E-4)*B+.815530498066373E-4)*B                          
     2  -.2064476319159326E-3)*B+.4677678114964956E-3)*B                        
     3  -.9083465574200521E-3)*B+.002099759035077063E0)*B                       
      G=(((((((G-.002851501243034649E0)*B+.0111538196719067E0)*B                
     1  -.2669510287555266E-3)*B+.07424900794340127E0)*B                        
     2  +.08157691940138868E0)*B+.4118403304219815E0)*B                         
     3  +.4227843350985181E0)*B+.9999999999999999E0                             
      VALUE=G*C                                                                 
      GO TO 7                                                                   
    6 D=1.0E0/Y                                                                 
      C=D*D                                                                     
      G=(((((((((-1.392432216905901E0*C+.1796443723688306E0)*C                  
     1  -.02955065359477124E0)*C+.00641025641025641E0)*C                        
     2  -.0019175269175269175E0)*C+.8417508417508418E-3)*C                      
     3  -.5952380952380952E-3)*C+.79365079365079365E-3)*C                       
     4  -.002777777777777778E0)*C+.08333333333333333E0)*D                       
     5  +.9189385332046727E0+(Y-.5E0)*LOG(Y)-Y                                  
      VALUE=EXP(G)                                                              
    7 GAMMA=VALUE                                                              
      RETURN                                                                    
      END 
