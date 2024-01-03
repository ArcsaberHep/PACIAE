        program main_23   
c       administrates MC simulation for relativistic pp(lp), pA (Ap), 
c        AB (lA), & e+e- collisions 
c	the program composes of main_23.f,parini_23.f,parcas_23.f,sfm_23.f, 
c        coales_23.f,hadcas_23.f,p_23.f, and analy.f
c       main_23.f: administrates the MC simulation
c       parini_23.f: generates a partonic initial state for colliding system
c       parcas_23.f: performs parton rescattering, where 2->2 processes
c        are considered only and LO pQCD cross section or its regularized
c        approximation is used
c       sfm_23.f: hadronization with LUND string fragmentation model
c       coales_23.f: hadronization with Monte Carlo coalescence model
c        selects sfm_23.f or coales_23.f by parameter adj1(12) 
c       hadcas_23.f: performs hadronic rescattering
c       p_23.f: pythia 64.28 with a little bit modifications
c       analy.f: an example of event analysis subroutine, user is free
c        to replace it with your own one

c250420 note: because of history reason, all of the Fortran programs are 
c        written with mode of regardless either capital lette or small 
c        letter. typing ": set ic + enter key" before searching.

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000,mplis=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
        common/pyint1/mint(400),vint(400)
        COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
	common/pycidat2/kfmaxt,nont2,param(20),weigh(600)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
	common/sa6/kfmaxi,nwhole
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
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
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 220110
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3,
     c   parj21,parj4,adiv,gpmax,nnc   ! 020708 070417 010518
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620   
        common/sa29/parp78,lcub   ! 150612 yan 070417
	common/sa30/vneump,vneumt,mstptj   ! 241110 100821 230722
        common/sa31/rmax,bbb(200),TA1(200),TA2(200),TA1A2(200),
     c  part1(200),part2(200),binn(200)   ! 020511 020718
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 220312 240412 131212
        common/sa34/itorw,iikk,cp0,cr0,kkii   ! 010418 010518 050920a
        common/sa35/ncpart,ncpar(kszj)   ! 280722
        common/sa6_c/ithroq,ithrob,ithroc,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)   ! 210921  
	common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5) ! 010518
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)   ! 050603
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)   ! 050603
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	common/show/vip(mplis),xap(mplis)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	common/count/isinel(600)
        common/ctllist/npctl,npinel(600),npctl0,npctlm ! 061103 180121 230121 
        common/ctllist_p/nreac(9),nrel   ! 071103
	common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 250209
        common/throqb/iprlth,nonqb,rpth(4,mplis),ppth(4,mplis),
     c  idpth(mplis),rmpth(mplis)   ! 010718
        common/anly1/an(40,5,20),bn(20),anf(40,5,20),bnf(20)   ! 281219
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5) ! 280620
        common/ancoal/icoal1,icoal2,xkappa1,xkappa2 ! yan, 200820
        dimension san(40,5,20),sbn(20),sanf(40,5,20),sbnf(20)   ! 070419
	dimension sao(40,5,20),sbo(20),saof(40,5,20),sbof(20)   ! 070419
        dimension skapa(6),skapao(6),snreac(9)   ! 020708 220110 010518
	dimension c(5),dinel(600),dineli(600),sthroe(4),wthroe(4)
        dimension einel(600),eineli(600)   ! 140820
	dimension bpp(20),kdiq(kszj,5),dgmas(kszj)
	dimension acoll(20),acollp(20),acollt(20)
	dimension sbp(20),numbth(3)   ! 010718	 
        dimension nreaco(9),pl(100,5)   ! 220110 260314
	dimension ksin(kszj,5),psin(kszj,5),vsin(kszj,5)   ! 010518 230618
        real nmin,nminf,ncha,nchaf   ! 020203
c010418 itorw: =1 executing pyevnt, =2 executing pyevnw
c260314 pl(ii,5): four momentum and mass of ii-th lepton
c	adj1(i), i=
c	1: k factor used in parton rescattering,k=0: no parton rescattering
c	2: parameter \alpha_s ('as' in program) in parton rescattering
c	3: parameter (\mu)^2 ('tcut' in program): to avoid divergence in 
c        calculating parton-parton differential cross section in parcas_23.f 
c	4: parameter idw, # of integration intervals in parcas_23.
c	5: =0, w/o nuclear shadowing, 
c          =1, Wang's nuclear shadowing (PL, B527(2002)85),
c          =2, EPS09 nuclear shadowing (JHEP 07(2012)073)
c	6: parameter 'a', i.e. parj(41) in Lund string fragmentation function
c	   if adj1(12)=0
c	   parameter 'a' in Field-Feynman string fragmentation function 
c           if adj1(12)=1
c       7: parameter 'b', i.e. parj(42) in Lund string fragmentation function
c       8: i.e. mstp(82) in PYTHIA64
c       9: i.e. parp(81) (D=1.9 GeV/c), effective minimum transverse momentum  
c	   of multiple interactions if mstp(82)=1                      
c	10: parp(31),k factor in pythia64 
c	11: time accuracy used in hadas_23.f  
c	12: model of hadronization, =0 string fragmentation; =1: coalescence
c	13: dimension of meson table considered if adj1(12)=1
c       14: dimension of baryon table considered if adj1(12)=1
c	15: default string tension  
c	16: number of loops in deexcitation of energetic quark in coales_23.f
c	17: threshold energy in deexcitation of energetic quark in coales_23.f
c	18: =0, rest partons hadronize by string fragmentation
c           =1, rest partons hadronize by coalescence
c	19: time accuracy used in parcas_23.f ('dddt' in program)
c	20: =0 exact pQCD parton-parton cross section
c	    =1 limited and regularized parton-parton cross section (B. Zhang,
c              Comput. Phys. Commun. 109(1998)193)
c	    =2 the same as 0 but flat scattering angle distribution is assumed
c           =3 the same as 1 but flat scattering angle distribution is assumed
c	21: =0 and 1, without and with phase space adjudgment in coales_23.f
c	22: critical value of the product of radii both in coordinate and 
c	 momentum phase space (4 is assumed) used in coales_23.f 
c	23: =0 LUND fragmentation function is used in subroutine 'ffm' in 
c              coales_23.f
c	    =1 Field-Feymman fragmentation function is used 
c       24: =tl0,the virtuality cut in time-like radiation in parcas_23.f, 
c            4*tl0 is assumed
c	25: \Lambda_QCD in parcas_23.f
c	26: number of random number thrown away
c       27: largest momentum allowed for particle ('dpmax')
c       28: largest position allowed for particle (drmax=para10*max(rnt,rnp),  
c           which is 1 in usu.dat and is recalculated in the running)  
c	29: width of two dimension Gaussian distribution sampling px and py of
c           produced quark pair in deexcitation of the energetic parton   
c           in coales_23.f
c	30: maximum $p_T^2$ in above two dimension Gaussian distribution 
c	31: parj(1) in pythia64
c	32: parj(2) in pythia64
c	33: parj(3) in pythia64
c	34: parj(21) in pythia64
c	35: mstp(91) in pythia64, selects parton transverse momentum 
c           (k_{\perp}) distribution inside hadron; 
c	    =1, Gaussian; 
c           =2, exponential
c	36: =0 without phenomenological parton energy loss in parcas_23.f
c	    =1 otherwise
c	37: the coefficient ('c') in phenomenological parton energy loss
c	38: pt cut in phenomenological parton energy loss 
c	39: width of Gaussian k_{\perp} distribution in hadron if mstp(91)=1
c           width of exponential k_{\perp} distribution in hadron if mstp(91)=2
c	40: =1 simulation ends after partonic initiation
c	    =2 simulation ends after partonic rescattering
c	    =3 hadronization by coalesce model follows partonic initiation 
c	    =4 simulation ends after hadron rescattering

c       smadel: small perpurbation of ellipse relative to circle
c       parecc: a parameter converting initial spatial space eccentricity 
c        to final momentum space ellipticity
c       iparres: =0 considers elastic parton-parton collisions only in 
c                   parcas_23.f
c	         =1 with inelastic parton-parton collisions as well

c       pincl (pscal): four momentum and mass of incident (scatterd) lepon
c       pinch: four momentum and mass of incident hadron
c        vnu: \nu; fq2: Q^2=-q^2; w2l: W^2; yyl: y; zl: z; xb: x_B; pph: P_h

c       para1_1: nn total cross section, used in parini_23.f 
c       para1_2: nn total cross section  used in hadcas_23.f
c	dni: nucleon number density

c	pj and ej: (pt)**2 and transverse energy of J/psi 
c	pjp and ejp: (pt)**2 and transverse energy of (J/psi) prime

c	acoll: a array, demension of which should be larger
c	 than or equal to 'nmax' (interval # of impact parameter)
c       the dimension of 'bpp' must be small than or equal to 'nmax'

c       ipden: =0, if projectile is proton
c              =2, for e+e- collisions  ! 180921 yan
c              =1, if projectile is nucleus
c	       =11, projectile is e- (e+)  
c	       =12, projectile is nu_e (nu_ebar)  
c	       =13, projectile is mu- (mu+)  
c	       =14, projectile is nu_mu (nu_mubar)
c	       =15, projectile is tau- (tau+)  
c	       =16, projectile is nu_tau (nu_taubar)  
c       itden: =0, if target is proton
c              =1, if target is nucleus
c              =2, for e+e- collisions  ! 180921 yan

c       suppm: the upper bound in sampling the radius of projectile nucleon
c       suptm: the upper bound in sampling the radius of target nucleon
c       suppc: the maximum radius in sample for projectile
c       suptc: the maximum radius in sample for target

c       r0p: projectile radius 
c       r0t: target radius 
c       pio: 3.1416
c	bp: impact parameter (=0 for pp,lp and lA)
c	iii: current run number
c	coor: CM position of collision system
c	ispmax: maximum # of particle KF code wanted to statistics
c	ispkf(i): KF code of i-th particle wanted to statistics
c	kfmax: the maximum # of particle KF code considered 
c	kfaco(i): KF code of i-th particle
c       numb(i): order number of last particle among particles with same flavor 
c        code of kfaco(i) in particle list

c	ifram: = 0 for fixed target, = 1 for collider 
c	cspipi (fm^2): total cross section of pion + pion
c	sig (fm^2): cross section of pion + pion to kaon + kaon
c	cspin (fm^2): total cross section of pion + nucleon interaction
c	cskn (fm^2): total cross section of kaon + nucleon interaction
c	csnn (fm^2): total cross section of n + n interaction
c	rcsit: ratio of inelastic to total cross section
c	disbe(i,j): allowable minimum distance between two particles of
c        kfaco(i) & kfaco(j).
c	c17(i,1-3): three position of particle i
c	tp(i): time of particle i
c	ishp(i): =1 if i-th particle inside the simulated volume
c	         =0 otherwise 
c	tau(i): formation time of particle i.
c	isinel(i): = 0 without i-th inelastic process in hadcas_23.f
c	           = 1 with i-th inelastic process in hadcas_23.f

c       nrel: statistics of blocked parton-parton scattering process in
c        parcas_23.f
c       nreac(i): statistics of successful i-th parton-parton scattering 
c        process in parcas_23.f
c       npinel(592): # of nn collision calling pythia' in parini_23.f 
c       npinel(593): # of nn collision not calling pythia' in parini_23.f
c	noel : statistics of elastic collisions in hadcas_23.f
c       noinel(i): statistics the i-th inelastic collisions in hadcas_23.f
c	nosc = 1 : pythia type output only
c	       2 : OSCA1999 standard output (event-by-event) as well
c	       3 : no used 
c020708
c       itime: number of strings in current event
c       astr: no use 
c       gtime: number of gluons in current event 
c       akapa(1): parj(1)   
c       akapa(2): parj(2) 
c       akapa(3): parj(3) 
c       akapa(4): parj(4) 
c       akapa(5): parj(21) 
c       akapa(6): effective string tension
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
c240219             77: photons from hadronization

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
c	open(98,file='databs_g.dat',status='unknown')   ! 260219
c260219	databs_g.dat: record gammas E-by-E
c       open(99,file='databs_h.dat',status='unknown')   ! 240119 260219
c260219 databs_h.dat: record hadrons E-by-E
c120214 if projectile is lepton, set nap =1 

c       reads input variables for event generation
	read(11,*)neve,nout,nosc   
	read(11,*)nap,nzp,nat,nzt
	read(11,*)ddt,dtt,bmin,bmax,nmax   ! 241108
	read(11,*)kjp21,ifram,para7,para10,kjp20
	read(11,*)pio,ipden,itden
c       reads input variables for event analyses
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
c       reads input variables for event generation
200     read(11,*)parp21,parp22,win   
	read(11,*)ttaup,taujp,iabsb,iabsm,nchan   ! 241108
	read(11,*)para13,para14,psno,para15,para16,ajpsi,vneum
	read(11,*)para1_1,para1_2,para2,para4   
	read(11,*)tdh,cptl,cptu,cptl2,cptu2,itnum   ! 241108
	read(11,*)mstu21,mstj1_1,mstj1_2,mstj2,mstj3,itorw   ! 160617 010418 
	read(11,*)(adj1(i),i=1,10)
	read(11,*)(adj1(i),i=11,20)
        read(11,*)(adj1(i),i=21,30)
        read(11,*)(adj1(i),i=31,40)
        read(11,*)kjp22,kjp23,kjp24,parp78,mstptj   !  100821 230722
	read(11,*)parecc,iparres,smadel,dparj4,cp0,cr0,seco   ! 120219 260219
	close(11)

c	bmin,bmax : minimum and maximum impact parameters, bmin=bmax means
c	 definite impact parameter, 2*nmax: the number of 
c	 intervals segmented in [bmin,bmax]	
c       nchan=0: inelastic (INEL)
c	nchan=1: Non Single Difractive (NSD) 
c	nchan=2: qqb --> gamma^*/Z^0, used to generate Drell-Yen
c	nchan=3: J/psi production
c	nchan=4: heavy-flavor production
c	nchan=5: direct photon
c	nchan=6: soft only
c	nchan=7: W+/- production
c       nchan=8: pythia default (msel=1)
c       nchan=9: Z0 production  
c       setting of nchan=0,1,3,4,5,7,8 and 9 is ready

c	neve : # of events to be generate
c	nap (nzp) : # of nucleons (protons) in projectile nucleus
c	nat (nzt) : # of nucleons (protons) in target nucleus
c       for e^-A: formally sets nap=1,nzp=-1,ipden=11,itden=1, kf=11; 
c       for e^+A: formally sets nap=1,nzp=1,ipden=11,itden=1, kf=-11;
c       for nu_eA: formally sets nap=1,nzp=-1,ipden=12,itden=1, kf=12;
c       for nu_ebarA: formally sets nap=1,nzp=1,ipden=12,itden=1, kf=-12; 
c210921 here formally seting nzp=-1 is precently, when one uses it later, 
c        one should make nzp=iabs(nzp)
c            .
c            .
c            .
c060813

c	t0 : average proper formation time at rest
c	ddt : time accuracy
c	dep : the accuracy in four momentum conservation
c	rou0 : normal nuclear density
c	rao : enlarge factor in the radius of simulated volume
c	kjp20: =1 constant cross sections 
c	       =0 energy dependent cross sections
c	kjp21: = 0 without hadron rescattering, 
c	       = 1 with hadron rescattering
c
c       kjp22 = 1 : variable single string tension and parj(1) etc. 
c       kjp22 = 2 : variable multiple string tension and parj(1) etc.
c	kjp22 = 3 : variable (single+multiple) string tension and parj(1) etc.
c	kjp22 = 4 : default string tension and parj(1) etc.
c       kjp23: = 1 npart is calculated by geometric model
c       kjp23: = 2 npart is calculated by Glauber model
c       kjp24: = 1 sharp sphere in Glauber model
c       kjp24: = 2 Woods-Saxon in Glauber model
c

c	pathn: collision numer suffered by projectile nucleon in target nucleus
c	tdh and itnum: time step and number of time steps used in subroutine 
c	 'flow_t'
c	cptl,cptu (cptl2,cptu2) : pt cut in 'flow_t' for particle 1 (2)

c	param(1)=para1   ! given in 'sysini' in parini_23.f
	param(2)=para2   ! read from usu.dat
	param(4)=para4
	param(7)=para7
	param(8)=ddt
	param(10)=para10
c	totle cross-section of J/Psi + n
	param(13)=para13
c	totle cross-section of J/Psi + meson
	param(14)=para14
c	totle cross-section of Psi' + n
	param(15)=para15
c	totle cross-section of Psi' + meson
	param(16)=para16
c020511 # of segments used in integration
        idw=adj1(4)   
	mstp(82)=adj1(8)
	parp(81)=adj1(9)
	parj(2)=adj1(32)
        parj2=parj(2)   
	parj(1)=adj1(31)
	parj(3)=adj1(33)
	parj(4)=dparj4   
	parj(21)=adj1(34)
c020708
        parj1=parj(1)
        parj3=parj(3)
	parj4=parj(4)   
        parj21=parj(21)
c020708

c270219 recalculate parj(1) with popcorn mechanism correction
        parj10=parj1
        wx0=parj3 
        wy0=parj4
        wrho0=parj2
        wnumer=1.+2.*wx0*wrho0+9.*wy0+6.*wx0*wy0*wrho0+
     c   3.*wx0*wx0*wy0*wrho0*wrho0
        wdenom=2.+wrho0
        ww0=wnumer/wdenom
        wweff=ww0
        parj(1)=seco*wweff*(parj10/seco/ww0)

        akapa(1)=parj(1)
        akapa(2)=parj(2)
        akapa(3)=parj(3)
        akapa(4)=parj(4)
        akapa(5)=parj(21)
        akapa(6)=1.   ! default effective string tension

	mstp(91)=adj1(35)
	if(mstp(91).eq.1)parp(91)=adj1(39)
	if(mstp(91).eq.2)parp(92)=adj1(39)
c	parp(93)=adj1(30) ! upper cut-off for k_perp distribution in hadron  
c	parj(21)=adj1(29)
        parp(2)=parp21
c       parp21: lowest CM energy for calling 'pythia' (D=10.), for 
c	 the case of nchan not equal to 3   
c       parp22: lowest CM energy for calling 'pythia' (D=10.), for 
c	 the case of nchan=3

	mstp(33)=1      
c	inclusion of k factor multiplying the differential cross sections for 
c       hard parton-parton process
	parp(31)=adj1(10)   ! D=1.5

c070417 contral the strength of colour reconnection
	parp(78)=parp78   

c	mstj1_1: =6, with inelastic processes 4, 6, and 7 (parcas_23.f)
c                =7, with inelastic process 7 only (parcas_23.f)
c	mstj1_2: =0, w/o final state time-like parton shower if iparres=1
c                =1, w/ final state time-like parton shower if iparres=1
c230722 mstptj: input value of mstp(111) (mstj(1)) for pp, pA (Ap), and AA 
c               (for e+e-)
c       mstptj=0: simulation pass through partonic initial state, partonic
c        rescattering, fragmentation, and hadronic rescttering states
c       mstptj=1: PYTHIA-like simulation
c       gluon jet fragmentation scheme in IF
	mstj(2)=mstj2
c	how the particles share the momentum in IF
	mstj(3)=mstj3

c       check on possible errors during execution
	mstu(21)=mstu21   ! 120603
c
c	parameters in Lund string fragmentation function  
	parj(41)=adj1(6)   ! D=0.3
	parj(42)=adj1(7)   ! D=0.58

c	initiation of event averaged variales
        dnmin=0.
        dnminf=0.
        dncha=0.
        dnchaf=0.
	do i1=1,20   
	sbn(i1)=0.
	sbnf(i1)=0.
        enddo   
        do i1=1,40   
	do i2=1,5
	do i3=1,20
	san(i1,i2,i3)=0.
	sanf(i1,i2,i3)=0.
	enddo
	enddo
	enddo

        stime_ini=0.
        stime_par=0.
        stime_had=0.
        snnc=0.
        sgtime=0.
        sgtimeo=0.
        sitime=0.
        sadiv=0.
        sgpmax=0.
        do i1=1,6   
        skapa(i1)=0.
        skapao(i1)=0.
        enddo
        segam=0.   
        segam1=0.   
        segam2=0.   
        segam3=0.   
	spathn=0.   
	sevbin=0.   
        rinel=0.
        rel=0.
	sinel=0.
	sel=0.
c170121
        swoun=0.
        snpctl0=0.
        snpctlm=0.   ! 180121
c170121       
        snpar=0.   ! 280722
	do i1=1,600
        dinel(i1)=0.
        einel(i1)=0.
        enddo
	do i1=1,20
	acoll(i1)=0.
	acollp(i1)=0.
	acollt(i1)=0.
	sbp(i1)=0.
	enddo

	volum=4.*3.1416/3.*2.**3
c261002	volume of sphere with radius of 2 fm in position phase space 
c200601
	skpar=0.
	sknn=0.
	skpp=0.
	sknp=0.
        skep=0. ! statistic of ep and en collisions with calling pythia 060813
c200601
	sthroq=0.
	sthrob=0.
        sthroc=0.
	do i1=1,4
	sthroe(i1)=0.
	enddo
	adj12=adj1(12)
        adj18=adj1(18)   
	adj140=adj1(40)  

c	gives values to some important variables
	call sysini(win)
	open(5,file='sxp.out',status='unknown')   ! sa 26/05/99
	open(9,file='rms0.out',status='unknown')
	open(34,file='oscar.out',status='unknown')


c020511	calculates the nuclear overlap function etc.
c       if((ipden.eq.0 .and. itden.eq.0) .or. ipden.ge.11)goto 80001!180921 yan
        if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.
     c   itden.eq.2) .or. ipden.ge.11)goto 80001   ! 180921 yan        
        csnn1=csnn*10   ! csnn in fm^2 csnn1 in mb
        idw1=idw/50   ! *100
c        write(9,*)'csnn1,kjp24,idw1=',csnn1,kjp24,idw1
        if(ipden.lt.2)call overlap(nap,nat,rnp,rnt,csnn1,kjp23,kjp24,
     c	 rou0,idw1)   ! 060813 120214 changed from .ne. to .lt.
c020511
80001	continue   ! 310518
	adj1(28)=para10*dmax1(rnt,rnp)
	iii=0


	write(9,*)'win,nap,nzp,nat,nzt=',win,nap,nzp,nat,nzt
	write(9,*)'neve,nout,nosc=',neve,nout,nosc
	write(9,*)'dtt,bmin,bmax,nmax,parp78,mstptj=',dtt,bmin,bmax,
     c   nmax,parp78,mstptj   ! 150612 yan 070417 100821 230722
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
	write(9,*)'tdh,cptl,itnum=',tdh,cptl,itnum
	write(9,*)'cptu,cptl2,cptu2=',cptu,cptl2,cptu2
	write(9,*)'mstu21,mstj1_1,mstj1_2,mstj2,mstj3,itorw=',
     c	 mstu21,mstj1_1,mstj1_2,mstj2,mstj3,itorw   ! 160617 010418    
c210803
	write(9,*)'adj1=',(adj1(i),i=1,10)
	write(9,*)'adj1=',(adj1(i),i=11,20)
        write(9,*)'adj1=',(adj1(i),i=21,30)
        write(9,*)'adj1=',(adj1(i),i=31,40)
c210803
	write(9,*)'parecc,iparres,smadel,dparj4,cp0,cr0,seco=',
     c	 parecc,iparres,smadel,dparj4,cp0,cr0,seco   ! 120219 260219
	if(iflmax.ne.0)then
	do kk=1,ispmax
	do i=1,iflmax
	write(9,*)(afl(kk,i,j),j=1,2)
	enddo
	enddo
	endif
c	write(9,*)'adj12=',adj12   ! sa
	ich=0	
	time=0.
c       write(9,*)'rnt,para10,adj1(28)=',rnt,para10,adj1(28)   ! sa
	jjj=1


c       psno: =0 fixed impact parameter
c       psno: =1 systematic sampling method
c       psno: =2 random sampling method

c	calculate the impact parameter etc. 
c	for given b (impact parameter)
c260718	if(dabs(bmin-bmax).lt.10d-4)then   ! i. e. case of psno=0. 280113
	if(psno.eq.0)then   ! 260718
	bp=bmin
        r4=rnp
        if(rnt.gt.rnp)r4=rnt
        rr4=bp/r4
        vneu=dexp(-rr4*rr4)
c	calculates the overlap region of two nuclei at given b by  
c020511	 interpolation method
	if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.   ! 180921 yan
     c   itden.eq.2) .or. ipden.ge.11)goto 80002   ! 020718, 180921 yan
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
c291118
c180219	if(ipden.eq.1 .and. itden.eq.1)then   ! A+B
	anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
c180219	elseif(ipden.eq.0 .and. itden.eq.1)then   ! p+A
c	anbin=ta2(ibpp)   ! overlap function of B nucleus (1/fm^2) 
c	elseif(ipden.eq.1 .and. itden.eq.0)then   ! A+p
c	anbin=ta1(ibpp)   ! overlap function of A nucleus (1/fm^2) 020718
c	else
c180219	endif
c291118
        pir=part1(ibpp)
        tir=part2(ibpp)
	evbin=anbin*csnn   !  Nbin in current event 020718
	pirr=pir   ! 140219
	tirr=tir   ! 180219
c       write(9,*)'bp,ibpp,part1,part2,evbin=',bp,ibpp,pir,tir,evbin   ! 020718
c280113	endif
c020511
80002	if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.   ! 180921 yan
     c   itden.eq.2) .or. ipden.ge.11)then   ! 020718, 180921 yan
	pir=1.
	tir=1.
	evbin=1.
	pirr=1.   ! 140219
	tirr=1.
	endif
c060605
        vneump=pir   ! 111399
        vneumt=tir   ! 111399
        write(9,*)'psno,b,N_part_p,N_part_t,N_bin=',
     c   psno,bp,vneump,vneumt,evbin   ! 190309 280113
	goto 80003   ! 010518 changed from 300 to 80003
	endif

	if(psno.eq.1.)then   ! 280113
c	systematic sampling method for given interval of b according to b**2 law
	nmax2=2*nmax
	bmaxn=(bmax*bmax-bmin*bmin)/(2.*nmax)
	bmin2=bmin*bmin
	i2=0
	do i1=1,nmax2,2
	i2=i2+1
	bpp(i2)=dsqrt(i1*bmaxn+bmin2)
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
	acoll(i1)=dexp(-rr4*rr4)
c	calculate the overlap region of two nuclei at given b 
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
c291118
c180219	if(ipden.eq.1 .and. itden.eq.1)then   ! A+B
	anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
c180219	elseif(ipden.eq.0 .and. itden.eq.1)then   ! p+A
c	anbin=ta2(ibpp)   ! overlap function of B nucleus (1/fm^2) 
c	elseif(ipden.eq.1 .and. itden.eq.0)then   ! A+p
c	anbin=ta1(ibpp)   ! overlap function of A nucleus (1/fm^2) 020718
c	else
c180219	endif
c291118
        pir=part1(ibpp)
        tir=part2(ibpp)
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
	aneump=stbp/dfloat(i2)   ! 241110
	aneumt=stbt/dfloat(i2)   ! 241110
	vneum=stb/dfloat(i2)
	write(9,*)'psno,ave. b=',psno,stab   ! 280113
        write(9,*)'N_bin=',(acoll(i1)*csnn,i1=1,i2)   !! 280113
	write(9,*)'(N_part)_p=',(acollp(i1),i1=1,i2)   ! 241110
        write(9,*)'(N_part)_t=',(acollt(i1),i1=1,i2)   ! 241110
	write(9,*)'ave. N_part_p,N_part_t,N_bin=',
     c	 aneump,aneumt,vneum*csnn ! 191110      280113

c	average b in [bmin,bmax]
        avb=2./3.*(bmin+bmax)
c	above equation is correct when bmin=0 only
	r4=rnp
	if(rnt.gt.rnp)r4=rnt
	rr4=avb/r4
	avneu=dexp(-rr4*rr4)
c	calculate the overlap region of two nuclei at given b 
        ibpp=int(avb/0.1+1.0)
        ibpp=min(ibpp,200)
c291118
c180219	if(ipden.eq.1 .and. itden.eq.1)then   ! A+B
	anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
c180219	elseif(ipden.eq.0 .and. itden.eq.1)then   ! p+A
c	anbin=ta2(ibpp)   ! overlap function of B nucleus (1/fm^2) 
c	elseif(ipden.eq.1 .and. itden.eq.0)then   ! A+p
c	anbin=ta1(ibpp)   ! overlap function of A nucleus (1/fm^2) 020718
c	else
c180219	endif
c291118
        pir=part1(ibpp)
        tir=part2(ibpp)
	aanbin=anbin   ! 280113
	astbp=pir
	astbt=tir
c220110
	endif   ! 280113

80003	continue   ! 310518
        nrel=0
        nrea=0
        do i=1,9
        nreac(i)=0
        enddo
c220110
c280113
        if(psno.eq.2.)then
        averb=0.
        psnon=0.   ! N_bin in case of psno=2
        psnop=0.   ! parojectile N_part in case of psno=2
        psnot=0.   ! target N_part in case of psno=2
        endif
	nncoll=0
	vnlep=0.d0 ! statistics of the number of studied leptons 260314
c280113
	iran=adj1(26)
	if(iran.eq.0)goto 300
	do i1=1,iran
	thrr=pyr(1)
	enddo


c	loop over event 
300	iii=iii+1
c	print*,'loop over event, iii=',iii   ! sa
c220110
        do i1=1,9
        nreaco(i1)=nreac(i1)
        enddo
c220110
c180520
c260620 write(9,*)'begin iii=',iii
c       write(9,100)(nreac(i1),i1=1,9)
c260620 100     format(9(1x,i3))   
c180520
c250209
        ngam=0
        nsin=0   
        do i1=1,kszj
        do j1=1,5
        kgam(i1,j1)=0
        pgam(i1,j1)=0.
        vgam(i1,j1)=0.
        ksin(i1,j1)=0   
        psin(i1,j1)=0.   
        vsin(i1,j1)=0.  
c280620 
        ktrs(i1,j1)=0   
        ptrs(i1,j1)=0.   
        vtrs(i1,j1)=0.  
c280620 
        enddo
        enddo
	siijk=0   ! 201203
c061103
	noel=0
        do i=1,600
        noinel(i)=0
        npinel(i)=0   ! 140820
        enddo
c061103
c010507
        non6_c=123456   ! 141208
	if(psno.eq.0)goto 800   ! 260718  
 
	if(psno.eq.1)then   ! 280113
	sbp(jjj)=sbp(jjj)+1
	bp=bpp(jjj) 
	vneump=acollp(jjj)   ! 241110
	vneumt=acollt(jjj)   ! 241110
	evbin=acoll(jjj)*csnn   ! 020718
	pirr=vneump   ! 140219
	tirr=vneumt   ! 180219
        goto 800   ! 280113
        endif   ! 280113  

	if(psno.eq.2)then   ! calcaulate impact parameter etc. for psno=2
        bmin2=bmin*bmin   ! 190620
        bp=sqrt(pyr(1)*(bmax*bmax-bmin2)+bmin2)
c       calculate the overlap region of two nuclei at given bp 
c020511
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
c291118
c180219	if(ipden.eq.1 .and. itden.eq.1)then   ! A+B
	anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
c180219	elseif(ipden.eq.0 .and. itden.eq.1)then   ! p+A
c	anbin=ta2(ibpp)   ! overlap function of B nucleus (1/fm^2) 
c	elseif(ipden.eq.1 .and. itden.eq.0)then   ! A+p
c	anbin=ta1(ibpp)   ! overlap function of A nucleus (1/fm^2) 020718
c	else
c	endif
c291118
        pir=part1(ibpp)
        tir=part2(ibpp)
        evbin=anbin*csnn   ! 020718
	pirr=pir   ! 140219
        tirr=tir   ! 180219
c280113        endif
c020511
        vneump=pir
        vneumt=tir
        if(iii.eq.1)write(9,*)'psno,b,N_part_p,N_part_t,N_bin=',
     c   psno,bp,vneump,vneumt,evbin   ! 280113 020718
        endif
c241110
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
	do i1=1,5
	ppsa(i1)=0.
	enddo

c	forbiden decay of particle (sets mdcy(...)=0)
c	mdcy(pycomp(111),1)=0
        mdcy(pycomp(310),1)=0   ! k0_S
        mdcy(pycomp(333),1)=0   ! phi
        mdcy(pycomp(3122),1)=0   ! Lambda
        mdcy(pycomp(-3122),1)=0
	mdcy(pycomp(443),1)=0   ! j/psi
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
        mdcy(pycomp(3112),1)=0   ! Sigma-
c	mdcy(pycomp(-3112),1)=0
        mdcy(pycomp(3222),1)=0   ! Sigma+
c	mdcy(pycomp(-3222),1)=0
        mdcy(pycomp(3312),1)=0   ! Xi-
        mdcy(pycomp(-3312),1)=0
c	mdcy(pycomp(3322),1)=0
c	mdcy(pycomp(-3322),1)=0
        mdcy(pycomp(3334),1)=0   ! Omega-
        mdcy(pycomp(-3334),1)=0
        mdcy(pycomp(23),1)=0   ! Z0
        mdcy(pycomp(24),1)=0   ! W+
        mdcy(pycomp(-24),1)=0   ! W-
c       mdcy(pycomp(1114),1)=0
c       mdcy(pycomp(2114),1)=0
c       mdcy(pycomp(2214),1)=0
c       mdcy(pycomp(2224),1)=0
c       mdcy(pycomp(213),1)=0   ! rho+
c       mdcy(pycomp(-213),1)=0   ! rho-  
c       mdcy(pycomp(113),1)=0   ! rho0
c       mdcy(pycomp(223),1)=0   ! omega 041202
c       mdcy(pycomp(413),1)=0
c       mdcy(pycomp(-413),1)=0
c       mdcy(pycomp(423),1)=0
c       mdcy(pycomp(-423),1)=0
c	mdcy(pycomp(13),1)=0
c	mdcy(pycomp(-13),1)=0
c060620
c       default pythia or user own selection for subprocesses
        if(nchan.eq.0)then
c       INEL
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
        endif
        if(nchan.eq.1)then
c       Non Single Difractive (NSD)
        msel=0
        msub(11)=1
        msub(12)=1
        msub(13)=1
        msub(28)=1
        msub(53)=1
        msub(68)=1
c       msub(91)=1
c       msub(92)=1
c       msub(93)=1
        msub(94)=1
        msub(95)=1
        endif
c090921
        if(nchan.eq.3)then
c       J/Psi production
        msel=0
        msub(86)=1
        msub(106)=1
        endif        
c090921        
        if(nchan.eq.7)then
c       W+/- production (nchan=7,isub=2,16,20,23,25,31)
        msel=0
        msub(2)=1
        msub(16)=1
        msub(20)=1
        msub(23)=1
        msub(25)=1
        msub(31)=1
        endif
        if(nchan.eq.9)then
c       Z0 production (nchan=9,isub=1,15,19,22,23,30)
        msel=0
        msub(1)=1
        msub(15)=1
        msub(19)=1
        msub(22)=1
        msub(23)=1
        msub(30)=1
        endif
c060620
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
        ebeam=dsqrt(win*win+0.938*0.938)   ! osc
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


c	creats pp (pA,Ap,AB,lp, lA, & e+e-) collision events


c	partonic initiation for the collision system   
	ijk=0
        time_ini=0.d0   ! 081010
c120620 mstp(111)=0   ! 050620 
        npctlm=0.   ! 180121
csa     print*,'be. parini iii,n=',iii,n   ! sa
c210921
        if((ipden.eq.0 .and. itden.eq.0) .or.
     c   (ipden.eq.2 .and. itden.eq.2))then   ! if pp e+e-

c       for pp,ppbar,pbarp,and e+e-
        n=0
        nbe=0   
        naf=0
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
        ifcom(i1)=0   
        ishp(i1)=0
        tau(i1)=0.
        enddo
        if(ipden.eq.0 .and. itden.eq.0)then
        mstp(111)=mstptj   ! 151021 230722
        if(ifram.eq.1)then   
        if(nzp.eq.1 .and. nzt.eq.1)call pyinit('cms','p','p',win)   
        if(nzp.eq.1 .and. nzt.eq.-1)call pyinit('cms','p','pbar',win)
        if(nzp.eq.-1 .and. nzt.eq.1)call pyinit('cms','pbar','p',win)   
        if(nzp.eq.0 .and. nzt.eq.0)call pyinit('cms','n','n',win)   
        if(nzp.eq.1 .and. nzt.eq.0)call pyinit('cms','p','n',win)
        if(nzp.eq.0 .and. nzt.eq.1)call pyinit('cms','n','p',win)   
        endif
        if(ifram.eq.0)then
        if(nzp.eq.1 .and. nzt.eq.1)call pyinit('fixt','p','p',win)   
        if(nzp.eq.1 .and. nzt.eq.-1)call pyinit('fixt','p','pbar',win)  
        if(nzp.eq.-1 .and. nzt.eq.1)call pyinit('fixt','pbar','p',win) 
        if(nzp.eq.0 .and. nzt.eq.0)call pyinit('fixt','n','n',win)   
        if(nzp.eq.1 .and. nzt.eq.0)call pyinit('fixt','p','n',win)
        if(nzp.eq.0 .and. nzt.eq.1)call pyinit('fixt','n','p',win)   
        endif        
c151021 endif
        if(itden.eq.0 .and. itorw.eq.1)call pyevnt  
        if(itden.eq.0 .and. itorw.eq.2)call pyevnw 
c151021 
        endif 
        if(ipden.eq.2 .and. itden.eq.2)then
        mstj(1)=mstptj   ! 230722
        call pyeevt(0,win)   ! for e+e- 
        endif
c151021
        call pyedit(2)
c       call pylist(1)


c230722
        if(mstptj.eq.1)then   !! 230722 
c       give four position to the particles generated in pythia ('pyjets')
        call ptcre(1,2,time)
c       write(22,*)'pp mstptj=',mstptj
c       call pylist(1)
        goto 998   ! toward hadron rescattering ('call hadcas')
        else   !! 230722
c230722

c       remove gamma from 'pyjets' to 'sgam'
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
c       give four position to the particles generated in pythia (in pyjets)
        call ptcre(1,2,time) ! arguments 1 and 2 make no sense indeed
c       if(ipden.eq.2 .and. itden.eq.2)print*,'af. ptcre'
c       remove hadrons from 'pyjets' to 'sbh' 
        call remo
        call pyedit(2)
c       if(ipden.eq.2 .and. itden.eq.2)print*,'af. remo'
c       preparing parton rescattering
c       'pyjets' to 'sbe'. etc.
        if(n.ge.1)then   ! 1
        do i1=1,n
        i3=i1   ! 050805
        kf=k(i1,2)
        kfab=iabs(kf)
c       identify diquarks        
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then   ! 2
c     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then   ! 2
        idi=idi+1
        ndiq(i1)=idi  
        endif   ! 2
c       'pyjets' to 'sbe'        
        do i2=1,5
        kbe(i3,i2)=k(i1,i2)
        pbe(i3,i2)=p(i1,i2)
        vbe(i3,i2)=v(i1,i2)
        enddo
        enddo
        nbeo=0
        nbe=i3
        endif   ! 1
c       if(ipden.eq.2 .and. itden.eq.2)print*,'af. pyjets to sbe'
c       break up diquark and give four momentum and four position
c        to broken quarks (working in 'pyjets')
        call break
        call pyedit(2)
c       if(ipden.eq.2 .and. itden.eq.2)print*,'af. break'
c       find number of strings and line number of first and last components
c        of each string
        nstr1=0
        jb=0
10000   do i1=jb+1,n
        if(k(i1,1).eq.2)then   ! i1 is 'A'
        do i2=i1+1,n
        if(k(i2,1).eq.1)then   ! i2 is 'V'
        nstr1=nstr1+1
        nstr1a(nstr1)=i1   ! line number of first component of nstr-th string
        nstr1v(nstr1)=i2   ! line number of first component of nstr-th string
        jb=i2
        if(jb.lt.n)goto 10000
        if(jb.eq.n)goto 20000
        endif
        enddo
        endif
        enddo
20000   continue
        nstr0=nstr1   ! 090620
c       if(ipden.eq.2 .and. itden.eq.2)print*,'af. count strings'
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
        idio=idi
c       if(ipden.eq.2 .and. itden.eq.2)print*,'be. 999,n=',n
        goto 999   ! toward parton rescattering ('call parcas')

        endif   !! 230722

        endif   ! if pp e+e- 
c210921 

	call parini(time_ini,parp21,parp22,win,psno,ijk) ! 081010 240513  
c120620  mstp(111)=1   ! 050620
        if(ijk.eq.1)goto 300   ! to avoide infinite loop in parcas 060813
        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813
csa     print*,'af. parini iii,n,adj140=',iii,n,adj140   ! sa
c260620 write(22,*)'af. parini iii,n,adj140=',iii,n,adj140
c       call pylist(1)
c260620 call prt_sbh(nbh,cc)
c140219
c	if(kjp22.eq.2 .or. kjp22.eq.3)then   ! 1
c	print*,'in parini nbh,win,kjp22=',nbh,win,kjp22
c	call prt_sbh(nbh,cc)

c       calculates wounded nucleons in colliding system
        unwoun=0.
        if(ifram.eq.1)then   ! 2 collider
        do i1=1,nbh
        if(dabs(pbh(i1,4)-0.5*win).le.1.)unwoun=unwoun+1.
        enddo
        woun=nap+nat-unwoun
c       woun: # of wounded nucleons
c       unwoun: # of unwounded nucleons
        endif   ! 2
        if(ifram.eq.0)then   ! 3 fixed target
        do i1=1,nbh
        ppp=pbh(i1,1)**2+pbh(i1,2)**2+pbh(i1,3)**2
        ppp=sqrt(ppp)
        if((ppp-win).le.0.1 .or. ppp.le.1.e-20)unwoun=unwoun+1.
        enddo
        woun=nap+nat-unwoun
        endif   ! 3
c	print*,'in parini woun,unwoun=',woun,unwoun
c	endif   ! 1
c230722
c       if(((ipden.eq.1.and.itden.eq.1).or.(ipden.eq.0.and.itden.eq.1)
c       c   .or.(ipden.eq.1.and.itden.eq.0)).and.mstptj.eq.1)then
        if(mstptj.eq.1)then
c       'sbh' to 'pyjets'
        n=nbh
        if(n.ge.1)then
        do i1=1,n
        do i2=1,5
        k(i1,i2)=kbh(i1,i2)
        p(i1,i2)=pbh(i1,i2)
        v(i1,i2)=vbh(i1,i2)
        enddo
        enddo
        endif
        do i1=n+1,kszj
        do i2=1,5
        k(i1,i2)=0
        p(i1,i2)=0.
        v(i1,i2)=0.
        enddo
        enddo
        goto 998   ! toward hadron rescattering ('call hadcas')
        endif
c230722                
999     continue   ! 230722
        if(mstptj.eq.0)then   ! 230722

c	no parton produced at all
	if(n.le.0)then
c	write(9,*)'neve,nncoll,n=',neve,nncoll,n   ! sa
	nncoll=nncoll+1
c060814	if(nncoll.gt.neve)then
c	stop 8888
c060814	endif
	iii=iii-1
	goto 300
	endif
c       print*,'af. parini n,adj140=',n,adj140   ! sa

c       throw away event with junction if iparres=1
        if(iparres.eq.1)then
        do i1=1,n   
        kf=k(i1,2)
        if(kf.eq.88)then
        iii=iii-1
        goto 300
        endif
        enddo
        endif
c180520

	if(adj140.eq.1 .or. adj140.eq.3)then   ! 290505 271205 020718 
c271205
        eevp=0.
        do i1=1,n
        eevp=eevp+p(i1,4)
        enddo
        eevh=0.
        do i1=1,nbh
        eevh=eevh+pbh(i1,4)
        enddo
c020718
	if(adj140.eq.3)then
        ithroq=0
        ithrob=0
        ithroc=0
        do i=1,4
        throe(i)=0.
        enddo
c	'pyjets' to 'parlist'
        iprl=n
        do i=1,n
        idp(i)=k(i,2)
        rp(4,i)=0.
c       parton cascade process is assumed to start at time 0.
        eee=p(i,4)
        pp(4,i)=eee
        if(eee.lt.1.e-12)eee=1.e-12   ! 180420
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
	call coales(iii,neve,nout,nap,nat,nzp,nzt)
	endif

        n44=0
        do j=1,nbh
        kf=kbh(j,2)
        if(kf.eq.22)then
        kbh(j,2)=44   ! '44': prompt direct photon
        n44=n44+1
        endif
        enddo
c       move "44" from 'sbh' to 'sgam'
        if(n44.gt.0)call remo_gam_sbh(44)
      
c       'sbh' to 'pyjets'
        if(nbh.eq.0)goto 5000   ! 261103
        do l=1,nbh
        l1=n+l
        do m=1,5
        k(l1,m)=kbh(l,m)
        p(l1,m)=pbh(l,m)
        v(l1,m)=vbh(l,m)
        enddo
        enddo
        n=n+nbh
        do i=n+1,kszj   ! 261103
        do j=1,5
        k(i,j)=0
        p(i,j)=0.
        v(i,j)=0.
        enddo
        enddo
5000    continue   ! 300407
c        if(iii.eq.3)call pylist(1)
        goto 888
        endif   ! 290505 271205 020718
c271205  
c	goto 889   ! temporal

c250209
        n44=0
        do j=1,nbh
        kf=kbh(j,2)
        if(kf.eq.22)then
        kbh(j,2)=44   ! '44': prompt direct photon
        n44=n44+1
        endif
        enddo
c       move "44" from 'sbh' to 'sgam'
        if(n44.gt.0)call remo_gam_sbh(44)
c250209
        call prt_sgam(n44,egam1,1)   ! 080419 160919 270220
c080419 egam1: energy of gamma after partonic initiation


c	partonic rescattering   
	if(n.lt.2)goto 889   ! 151302
c140718	if(itden.ne.1)goto 890   ! for e+e-,p+p,pbar_p, or p+pbar 080806 
c	write(9,*)'before parcas iii=',iii   ! sa

c	'pyjets' to 'parlist'
        iprl=n
        do i=1,n
        idp(i)=k(i,2)
        rp(4,i)=0.
c       parton cascade process is assumed to start at time 0.
        eee=p(i,4)
        pp(4,i)=eee
        if(eee.lt.1.e-12)eee=1.e-12   ! 180420
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
        time_par=0.d0   ! 081010
	iijk=0   ! 151203
csa     print*,'be. parcas iii,iijk,n=', iii,iijk,n
	call parcas(time_par,nnn,iijk,win,nap,rnt,rnp)   ! 120603 220110
c220110 nnn: nnn-th parton-parton interacion in a nucleus-nucleus collision
csa     print*,'af. parcas iii,iijk,n=',iii,iijk,n   ! sa
c2	write(9,*)'after parcas iii,iijk,time_par=',
c2     c	iii,iijk,time_par! sa
c120603
	if(iijk.eq.1)then
c	write(9,*)'iii,iijk=',iii,iijk
	goto 300   ! give up current event avoiding infinite collision loop
	endif
c120603
	if(iijk.eq.2)siijk=siijk+1   ! 201203
c       write(22,*)'af. parcas iii,n,iijk=',iii,n,iijk   ! bh
c       if(ipden.lt.11)call pyedit(2)
c       if(ipden.lt.11)call pyedit(1)
c       call pylist(1)
c       call prt_sbh(nbh,cc)
c       write(9,*)'af. parcas iii=',iii   ! 180520
c260620 write(9,100)(nreac(i1),i1=1,9)   ! 180520
c201203
	if(n.eq.0)goto 300   ! no parton at all, give up current event

        endif   ! 230722

	if(adj140.eq.2)then   ! 290505 271205  
c271205
        eevp=0.
        do i1=1,n
        eevp=eevp+p(i1,4)
        enddo
        eevh=0.
        do i1=1,nbh
        eevh=eevh+pbh(i1,4)
        enddo
c       'sbh' to 'pyjets'
        if(nbh.eq.0)goto 6000   ! 261103
        do l=1,nbh
        l1=n+l
        do m=1,5
        k(l1,m)=kbh(l,m)
        p(l1,m)=pbh(l,m)
        v(l1,m)=vbh(l,m)
        enddo
        enddo
        n=n+nbh
        do i=n+1,kszj   ! 261103
        do j=1,5
        k(i,j)=0
        p(i,j)=0.
        v(i,j)=0.
        enddo
        enddo
6000    continue   ! 300407
        goto 888
        endif   ! 290505 271205
c271205
c250209
        n55=0
        do i1=1,n
        kf=k(i1,2)
        if(kf.eq.22)then
        k(i1,2)=55
        n55=n55+1
        endif
        enddo
c       write(9,*)'iiii,iii,n55=',iiii,iii,n55   ! sa
c       move '55' from 'pyjets' to 'sgam'
        if(n55.gt.0)call remo_gam_par(55)
c250209
        call prt_sgam(n55,egam2,2)   ! 080419 160919 270220
c080419 egam2: energy of gammas after parton cascade 

890	continue   ! 020512
        if(adj12.ne.0)goto 889   ! coalescence 120520
c120520        if(adj12.ne.0.or.(adj12.eq.0.and.(nreac(4).gt.nreaco(4).or.
c120520     c   nreac(6).gt.nreaco(6).or.nreac(7).gt.nreaco(7))))goto 889 
! 020512 ->coalescence


c	recover parton configuration in 'sbe' (having diquark)
c       print*,'be. recover iii,adj12=',iii,adj12   ! sa 
c	loop over 'sbe'   ! 'pyjets' to 'sbe' in parini_23.f
	idii=0
c	write(9,*)'sbe nbe=',nbe
	do i=1,nbe
	kf=kbe(i,2)
        kfab=iabs(kf)
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
c060805     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then   ! 060805
	idii=idii+1
	do j=1,5
	kdiqq=kbe(i,j)
	kdiq(idii,j)=kdiqq
	enddo
	dgmas(idii)=pbe(i,5)
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
	j=npt(ndiqi) 
	do i1=1,5
        k(i,i1)=kdiq(idij,i1)
	enddo
c	write(9,*)'i,p=',i,(p(i,i1),i1=1,5)
c	write(9,*)'j,p=',j,(p(j,i1),i1=1,5)
	do i1=1,3   
	p(i,i1)=p(i,i1)+p(j,i1)
	enddo
	dimass=dgmas(idij)
	pi1=p(i,1)
	pi2=p(i,2)
	pi3=p(i,3)
	pi4=dsqrt(pi1*pi1+pi2*pi2+pi3*pi3+dimass*dimass)
	dele=dele+p(i,4)+p(j,4)-pi4
c	write(9,*)'dele,p(i,4),p(j,4),pi4=',dele,p(i,4),p(j,4),pi4
	p(i,4)=pi4
	p(i,5)=dimass
c	write(9,*)'ndiqi,pi1,pi2,pi3,pi4,dimass=',ndiqi,pi1,pi2,
c     c	pi3,pi4,dimass
c060805
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
c       move particle list,'pyjets' ('ndiq') one step downward from 
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
c       print*,'af. recover ,iii,n=',iii,n   ! sa

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
c	write(22,*)'af share jb,n=',jb,n
c	call prt_pyj(n,cc)
c180520
c260620 write(9,*)'be. sbe->pyjets nreaco iii=',iii
c       write(9,100)(nreaco(i1),i1=1,9)   
c       write(9,*)'be. sbe->pyjets nreac iii=',iii
c       write(9,100)(nreac(i1),i1=1,9)   
c180520

889	continue
c010518
        if(ipden.lt.11)call pyedit(2) 
        if(ipden.ge.11)call pyedit(1)
c       write(22,*)'af. recover iii,n,nbh=',iii,n,nbh   ! bh
c       call prt_pyj(n,cc)
c       call prt_sbh(nbh,cc) 
c010518


c	hadronization
c	write(9,*)'before hadronization iii=',iii   ! sa
csa     print*,'be. hadroni. iii,n,adj12=',iii,n,adj12   ! sa
c230618
        n77s=0   ! 270220
	nbe=0
        ithroq=0
        ithrob=0
        ithroc=0
        do i=1,4
        throe(i)=0.
        enddo
c	'pyjets' to 'parlist'
        iprl=n
        do i=1,n
        idp(i)=k(i,2)
        rp(4,i)=0.
c       parton cascade process is assumed to start at time 0.
        eee=p(i,4)
        pp(4,i)=eee
        if(eee.lt.1.e-12)eee=1.e-12   ! 170420
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
c230618
c220110
c       print*,'be. adj12=0'
	if(adj12.eq.0)then   !! 010518 
csa     print*,'iii,4,4o,6,6o,7,7o=',iii,nreac(4),nreaco(4),nreac(6)
csa     c   ,nreaco(6),nreac(7),nreaco(7)
c120520 if(nreac(4).gt.nreaco(4) .or. nreac(6).gt.nreaco(6)
c120520      c	 .or. nreac(7).gt.nreaco(7))then   ! 020512 010518
c020512	for inela. processes 4,6,and 7
c120520 call coales(iii,neve,nout,nap,nat,nzp,nzt)
csa     print*,'af. coales'
c120520 else   ! 020512 010518
c       otherwise
csa     print*,'be. sfm kjp22',kjp22

c	calculate the multiple effective string tension and parj(1) etc. 
	if(kjp22.eq.2 .or. kjp22.eq.3)then   ! 2 
	ampi=mint(31)
c250119 note mint(31)=0 for a low_pT event
c140219
c140219	pathn=(npinel(1)+npinel(592))/pirr
c180219	pathn=(npinel(1)+npinel(592))/(0.5*woun)  
c	numerator: NN collision # calculated
c	denominator: N_part of projectile nucleus calculated (pirr: 
c        from Glauber model)
c	print*,'bp,ampi=',bp,ampi
c	print*,'evbin,pirr,tirr=',evbin,pirr,tirr
c140219
        if(ampi.le.0)then   !
        ckapa=1.
c	print*,'ckapa=',ckapa
        elseif(ampi.gt.0)then   !
c120219	ampi=ampi*evbin
c180219
	if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.   ! 180921 yan
     c   itden.eq.2)) then   ! 180219 
	pathn=1.
	else   !!
	pathn=evbin/(pirr+tirr)   ! 180219
	endif   !!
c	evbin: N_bin of collision system (A+B)
c	pirr: N_part of projectile nucleus (Glauber)
c	tirr: N_part of target nucleus (Glauber)
c180219
	ampi=ampi*pathn   ! 120219
	ckapa=(1.+(ampi-1.)/(1.+1/(cp0**2)))**cr0
c       print*,'pathn,ckapa=',pathn,ckapa
c140219	ckapa=1.+(ampi-1.)
c	denoc=1.+1./(cp0**2)
c	if(denoc.lt.1.d-20)denoc=1.d-20
c	ckapa=ckapa/denoc
c	if(ckapa.lt.1.d-20)ckapa=1.d-20
c140219	ckapa=(ckapa)**cr0
        else   !
        endif   !
c250119
c       ckapa is multiple string tension
c       string tension of the pure qqbar string, kapa0, is assumed to be 1
c	write(9,*)'iii,jjj,bp,evbin,ampi=',iii,jjj,bp,evbin,ampi
c	write(9,*)'ampi,ckapa=',ampi,ckapa
c       print*,'iii,ampi,ckapa',iii,ampi,ckapa   ! sa
c020718 
c120718
	if(kjp22.eq.2)then
        parj_2=parj2**(1./ckapa)
        parj_21=parj21*((ckapa/1.)**(0.5))
        parj_1=parj1**(1./ckapa)
        parj_3=parj3**(1./ckapa)
        parj_4=parj4**(1./ckapa)
        parj(1)=parj_1
        parj(2)=parj_2
        parj(3)=parj_3
        parj(4)=parj_4
        parj(21)=parj_21
c270219 recalculate parj(1) with popcorn mechanism correction
        wxef=parj_3
        wyef=parj_4
        wrhoef=parj_2
        wnumer=1.+2.*wxef*wrhoef+9.*wyef+6.*wxef*wyef*wrhoef+
     c   3.*wxef*wxef*wyef*wrhoef*wrhoef
        wdenom=2.+wrhoef
        wweff=wnumer/wdenom
        parj(1)=seco*wweff*(parj10/seco/ww0)**(1./ckapa)
c270219
c290119
        akapa(1)=parj(1)
        akapa(2)=parj(2)
        akapa(3)=parj(3)
        akapa(4)=parj(4)
        akapa(5)=parj(21)
        akapa(6)=ckapa   
c290119
	endif
c120718 
	endif   ! 2

c       fragments strings all at once  
        if(kjp22.eq.2 .or. kjp22.eq.4)then
        kkii=0   ! 050920
        call sfm
        if(kkii.eq.2)then   ! 050920
        iii=iii-1
        goto 300   ! through away current event
        endif
        if(ipden.lt.11)call pyedit(2)
        if(ipden.ge.11)call pyedit(1)
c       write(22,*)'af. all at once n,nbh,mstp=',n,nbh,mstp(111)   ! bh
c       call pylist(1)
cc      call prt_pyj(n,cc)
c       call prt_sbh(nbh,cc)

c       removes gammas ('77') after hadronization from 'pyjets' to 'sgam'
        if(n.gt.0)then   
        n77=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=77   ! '77': photons after hadronization of current string
        n77=n77+1
        endif
        enddo
c       move "77" from 'pyjets' to 'sgam'
        if(n77.gt.0)call remo_gam_hadro(77)
        endif   
        goto 30001
        endif

c	calculate the single (single + multiple) effective string tension 
c	 and parj(1) etc. 
c010518
        if(kjp22.eq.1 .or. kjp22.eq.3)then   ! 070219
c070219 nstr=0
        itime=0
        gtime=0.
        adiv=0.
        gpmax=0.
        do i1=1,6
        akapa(i1)=0.
        enddo
        vfr24=3.5   ! parameter alpha
        vfr25=0.8   ! \sqrt(s_0) in GeV
	vfr252=vfr25*vfr25
        endif   ! 070219

c	find string and line number of its first and last components, 
c	 calculate parj(1) etc. if kjp22=1 or 3, and then hadronize string
csa     print*,'enter find string iii,n,kjp22',iii,n,kjp22

	nstr=0
c	'pyjets' to 'sin'
        nsin=n
        do i1=1,n
        do i2=1,5
        ksin(i1,i2)=k(i1,i2)
        psin(i1,i2)=p(i1,i2)
        vsin(i1,i2)=v(i1,i2)
        enddo
        enddo
        do i1=n+1,kszj
        do i2=1,5
        ksin(i1,i2)=0
        psin(i1,i2)=0.
        vsin(i1,i2)=0.
        enddo
        enddo
	do i1=1,kszj
        do i2=1,5
        k(i1,i2)=0
        p(i1,i2)=0.
        v(i1,i2)=0.
        enddo
        enddo
	naff=0
	do i1=1,kszj
	do i2=1,5
	kaff(i1,i2)=0
	paff(i1,i2)=0.
	vaff(i1,i2)=0.
	enddo
	enddo  

csa     print*,'begin loop over strings iii,n,kjp22',iii,n,kjp22
c       loop over string (begin)
c090219  jb=0
c090219 10000   do i1=jb+1,nsin
10001   do i1=1,nsin   ! 090219
c	find a string
        if(ksin(i1,1).eq.2)then   ! i1 is 'A'
        do i2=i1+1,nsin
        if(ksin(i2,1).eq.1)then   ! i2 is 'V'
        nstr=nstr+1
        nstra(nstr)=i1   ! line number of first component of nstr-th string
        nstrv(nstr)=i2   ! line number of last component of nstr-th string
csa     print*,'line # of first & last string, nsin=',i1,i2,nsin

        if(kjp22.eq.1 .or. kjp22.eq.3)then   ! 3
        toteng=0.0
        toten=0.0
        totglukt=0.0
        pmax=0.d0
        ggg=0.
	do i3=i1,i2
        toten=toten+psin(i3,4)   ! root_s, string total energy 
        pp2=psin(i3,1)**2+psin(i3,2)**2   ! pt*pt 
        ppp=dsqrt(pp2)
	if(ksin(i3,2).eq.21)then  
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
c       div: factor related to number of gluons and hardest gluon in 
c	 current string
c       pmax: transverse momentum of hardest gluon in current string
        adiv=adiv+div
        gpmax=gpmax+pmax
c       string tension of the pure qqbar string, kapa0, is assumed to be 1
c       calculate kapa and parj(1) etc. of current string
        effk2=(1.-div)**(-vfr24)
        itime=itime+1
        gtime=gtime+ggg
c       single
        if(kjp22.eq.1)then   ! 4
        parj_2=parj2**(1./effk2)   ! 210218
        parj_21=parj21*((effk2/1.)**(0.5))   ! 210218
        parj_1=parj1**(1./effk2)   ! 210218
        parj_3=parj3**(1./effk2)   ! 210218     
        parj_4=parj4**(1./effk2)   ! 210218  
        parj(1)=parj_1
        parj(2)=parj_2
        parj(3)=parj_3
        parj(4)=parj_4
        parj(21)=parj_21
c270219 recalculate parj(1) with popcorn mechanism correction
        wxef=parj_3
        wyef=parj_4
        wrhoef=parj_2
        wnumer=1.+2.*wxef*wrhoef+9.*wyef+6.*wxef*wyef*wrhoef+
     c   3.*wxef*wxef*wyef*wrhoef*wrhoef
        wdenom=2.+wrhoef
        wweff=wnumer/wdenom
        parj(1)=seco*wweff*(parj10/seco/ww0)**(1./effk2)
c270219
c       sum over strings
        akapa(1)=akapa(1)+parj(1)   ! 270219
        akapa(2)=akapa(2)+parj_2
        akapa(3)=akapa(3)+parj_3
        akapa(4)=akapa(4)+parj_4
        akapa(5)=akapa(5)+parj_21
        akapa(6)=akapa(6)+effk2
        endif   ! 4
c       single + multiple
        if(kjp22.eq.3)then   ! 5
        parj_2=parj2**(1./(ckapa*effk2))   ! 210218
        parj_21=parj21*(((effk2*ckapa)/1.)**(0.5))   ! 210218
        parj_1=parj1**(1./(ckapa*effk2))   ! 210218
        parj_3=parj3**(1./(ckapa*effk2))   ! 210218     
        parj_4=parj4**(1./(ckapa*effk2))   ! 210218 
        parj(1)=parj_1
        parj(2)=parj_2
        parj(3)=parj_3
        parj(4)=parj_4
        parj(21)=parj_21
c270219 recalculate parj(1) with popcorn mechanism correction
        wxef=parj_3
        wyef=parj_4
        wrhoef=parj_2
        wnumer=1.+2.*wxef*wrhoef+9.*wyef+6.*wxef*wyef*wrhoef+
     c   3.*wxef*wxef*wyef*wrhoef*wrhoef
        wdenom=2.+wrhoef
        wweff=wnumer/wdenom
        parj(1)=seco*wweff*(parj10/seco/ww0)**(1./(ckapa+effk2))
c270219
c       sum over strings
        akapa(1)=akapa(1)+parj(1)   ! 270219
        akapa(2)=akapa(2)+parj_2
        akapa(3)=akapa(3)+parj_3
        akapa(4)=akapa(4)+parj_4
        akapa(5)=akapa(5)+parj_21
        akapa(6)=akapa(6)+ckapa*effk2
cc      write(9,*)'af. parj(1) etc.=',parj(1),parj(2),parj(3),parj(21)
        endif   ! 5
        endif   ! 3

c       'sin' to 'pyjets' (the part of current string in 'sin' to 'pyjets')
        n=i2-i1+1
        do ii1=i1,i2
        ii3=ii1-i1+1
        do ii2=1,5
        k(ii3,ii2)=ksin(ii1,ii2)
        p(ii3,ii2)=psin(ii1,ii2)
        v(ii3,ii2)=vsin(ii1,ii2)
        enddo
        enddo
        do ii1=n+1,kszj
        do ii2=1,5
        k(ii1,ii2)=0
        p(ii1,ii2)=0.
        v(ii1,ii2)=0.
        enddo
        enddo
c090219	jb=i2   ! 160317
c       print*,'be. hadroni. of current string, nstr=',nstr
c       if(nstr.le.5)then
c       write(22,*)'be. current string har. nstr10,nstr=',nstr10,nstr ! bh
c       call pylist(1)
c       endif

c       hadronization of current string
	iikk=0
        kkii=0   ! 050920
	call sfm
        if(ipden.lt.11)call pyedit(2) 
        if(ipden.ge.11)call pyedit(1)
c280618
	if(iikk.eq.2 .and. ((ipden.eq.0 .and. itden.eq.0) .or.
     c  (ipden.eq.2 .and. itden.eq.2) .or. ipden.ge.11))then
	iii=iii-1
	goto 300
	elseif(iikk.eq.2)then
c	moves the part of current string in 'sin' to 'sbe
	do ii1=i1,i2
	nbe=nbe+1
        do ii2=1,5
        kbe(nbe,ii2)=ksin(ii1,ii2)
        pbe(nbe,ii2)=psin(ii1,ii2)
        vbe(nbe,ii2)=vsin(ii1,ii2)
        enddo
        enddo
	n=0   ! 230618
	else
	endif

        if(kkii.eq.2)then   ! 050920
c       moves the part of current string in 'sin' to 'sbe
        do ii1=i1,i2
        nbe=nbe+1
        do ii2=1,5
        kbe(nbe,ii2)=ksin(ii1,ii2)
        pbe(nbe,ii2)=psin(ii1,ii2)
        vbe(nbe,ii2)=vsin(ii1,ii2)
        enddo
        enddo
        n=0   
        endif
c       if(nstr.le.5)then
c       write(22,*)'af. hadroni. of current string, n,iikk=',n,iikk   ! bh
c       call pylist(1)
c       endif
csa     print*,'af. moves the part of current string n=',n
c280618
c240219
        if(n.gt.0)then   ! 240219
c       removes gammas ('77') after hadronization from 'pyjets' to 'sgam'
        n77=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=77   ! '77': photons after hadronization of current string 
        n77=n77+1
        n77s=n77s+1   ! 270220
        endif
        enddo
c       move "77" from 'pyjets' to 'sgam'
        if(n77.gt.0)call remo_gam_hadro(77)
        endif   ! 240219
csa     print*,'af. moves 77,n=',n

c       'pyjets' to 'aff'
	if(n.gt.0)then   ! 230618
c260620 if(nstr.le.5)then
c       if(ipden.lt.11)call pyedit(2)
c       if(ipden.ge.11)call pyedit(1)                
c260620 write(22,*)'af. current string had. nstr1,nstr=',nstr1,nstr
c       call pylist(1)
c       write(22,*)'summary and list of gammas'
c       call prt_sgam(ngam,egam,8)
c260620 endif
        do ii1=1,n
        ii3=naff+ii1
c190420
        if(ii3.ge.kszj)then
        print*,'iii,nstr,i1,i2,ii3=',iii,nstr,i1,i2,ii3
        stop 
        endif
c190420
        do ii2=1,5
        kaff(ii3,ii2)=k(ii1,ii2)
        paff(ii3,ii2)=p(ii1,ii2)
        vaff(ii3,ii2)=v(ii1,ii2)
        enddo
        enddo
        naff=ii3
	endif   ! 230618
c290518
	do ii1=naff+1,kszj
	do ii2=1,5
	kaff(ii1,ii2)=0
	paff(ii1,ii2)=0.
	vaff(ii1,ii2)=0.
	enddo
	enddo
cc	write(22,*)'af. sfm nstr,itime,naff,n,iikk=',nstr,itime,naff
cc     c	 ,n,iikk
cc	call pylist(1)  
c	call prt_aaff(naff)
c290518
		
        if(i2.lt.nsin)then   
c230618
c	revamps 'sin', i.e. moves parton list 'sin' ii (=i2-i1+1) steps 
c	 downward from i2+1 to nsin
	ii=i2-i1+1
	do j=i2+1, nsin
	do jj=1,5
        ksin(j-ii,jj)=ksin(j,jj)
        psin(j-ii,jj)=psin(j,jj)
        vsin(j-ii,jj)=vsin(j,jj)
        enddo
        enddo
	nsin=nsin-ii
c230618
	goto 10001
	endif

        if(i2.eq.nsin)goto 20001   ! without rest partons
        endif   ! i2
        enddo   ! i2
        endif   ! i1
        enddo   ! i1

c       rest partons which can not compose a string 
csa     print*,'be. rest partons iii,kjp22,nbe,nsin=',iii,kjp22,nbe,nsin
	if(nsin.ge.1)then   !!
c       'sin' to 'sbe'
        do ii1=1,nsin
        nbe=nbe+1
        do ii2=1,5
        kbe(nbe,ii2)=ksin(ii1,ii2)     
        pbe(nbe,ii2)=psin(ii1,ii2)
        vbe(nbe,ii2)=vsin(ii1,ii2)
        enddo
        enddo
c290518
        do ii1=nbe+1,kszj
        do ii2=1,5
        kbe(ii1,ii2)=0
        pbe(ii1,ii2)=0.
        vbe(ii1,ii2)=0.
        enddo
        enddo
c	write(22,*)'rest partons iii=',iii
c	call prt_sbe(nbe,cc)
c290518
	endif  !!
csa     print*,'af. rest partons iii,nbe,nsin=',iii,nbe,nsin
c       loop over string endded
c	fragments string by string endded

20001   continue   ! from 'nstr00=nstr' to 'continue' on 030620
csa     print*,'loop over string end  naff,nbe,kjp22=',naff,nbe,kjp22

        if(kjp22.eq.1 .or. kjp22.eq.3)then   ! 070219
c	average over strings in current hh collision 010518
        atime=dfloat(itime)
c	itime, # of strings in current hh collision 010518
        if(atime.gt.0.)then   ! 7
	do i1=1,6   ! 010518
        akapa(i1)=akapa(i1)/atime
	enddo
        gtime=gtime/atime
        adiv=adiv/atime
        gpmax=gpmax/atime
c gtime: averaged # of gluons in a string in current hh collision 
        endif   ! 7
        endif   ! 070219
c200420 iiire=iii   ! 230219
c200420 kjp211=kjp21   ! 230219

c       'aff' to 'pyjets'
csa     print*,'enter aff to pyjets'
        n=naff
        do ii1=1,n
        do ii2=1,5
        k(ii1,ii2)=kaff(ii1,ii2)
        p(ii1,ii2)=paff(ii1,ii2)
        v(ii1,ii2)=vaff(ii1,ii2)
        enddo
        enddo
c200420 iii=iiire   !   230219
c200420 kjp21=kjp211   ! 230219
        do ii1=n+1,kszj
        do ii2=1,5
        k(ii1,ii2)=0
        p(ii1,ii2)=0.
        v(ii1,ii2)=0.
        enddo
        enddo
c	write(22,*)'hadron list event=',iii
c	call pylist(1)

c       write(22,*)'af. frag. s-by-s naff,n,nbh,nbe=',naff,n,nbh,nbe   ! bh
c       call prt_pyj(n,cc)
c       call prt_sbh(nbh,cc)
c       call prt_sbe(nbe,cc)
c120520        endif   ! 020512 010518 
30001	continue
	endif   !! 010518 (adj12.eq.0)
c220110
c       print*,'be. adj12=1'
	if(adj12.ne.0)then
c	write(9,*)'af parcas iprl=',iprl   ! sa
c	write(9,887)(idp(i2),i2=1,iprl)   ! sa
	call coales(iii,neve,nout,nap,nat,nzp,nzt)   ! coalescence
        goto 333
	endif
887	format(20(1x,i3))
        if(ipden.lt.11)call pyedit(2)   
        if(ipden.ge.11)call pyedit(1)
csa     print*,'af. call pyedit(2)  iii,naff,n,nbe=',iii,naff,n,nbe
c141208

c030820 removes g,q,qbar,qq & (qq)bar from 'pyjets' to 'sbe'
        call remop   ! 030820
csa     print*,'af. remove iii,n,nbe=',iii,n,nbe
c       write(22,*)'af. remop n,nbh,nbe=',n,nbh,nbe   ! bh
c       call prt_pyj(n,cc)
c       call prt_sbh(nbh,cc)
c       call prt_sbe(nbe,cc)

333     continue
c       print*,'af. call coales'
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
c010518
csa     print*,'af. sbh->pyjets iii,n,nbh,nbe=',iii,n,nbh,nbe
        call prt_sgam(n77s,egam3,3)   ! 080419 100919 270220
c080419 egam3: gamma energy after hadronization
c       write(22,*)'be. hadr. rescattering iii,n,nbe=',iii,n,nbe
c       call pylist(1) 
c       call prt_sbe(nbe,cc)
c       call prt_sbh(nbh,cc)

c230722
998     continue
        if(mstptj.eq.1)then
c       removes gammas ('77') after hadronization from 'pyjets' to 'sgam'
        n77=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=77   ! '77': photons after hadronization of current string
        n77=n77+1
        n77s=n77s+1   ! 270220
        endif
        enddo
c       move "77" from 'pyjets' to 'sgam'
        if(n77.gt.0)call remo_gam_hadro(77)
        endif
c230722

c	hadronic cascade (rescattering, HRS)
csa     print*,'prepare for call hadcas'
	if(kjp21.eq.1)then   ! 1 241103
        call filt
        do i=1,kfmax
        nup=numbs(i)
        enddo
        nbh1=n-nup
c       nup is the number of particles kept in 'pyjets' (joints HRS)
c       nbh1 is the number of particles storing in 'sbh' (not joints HRS)
c060813	lepton is not rescattering with hadrons
        if(nbh1.eq.0)goto 7000
        do i=nup+1,n
        nbh=i-nup
        do j=1,5
        kbh(nbh,j)=k(i,j)
        pbh(nbh,j)=p(i,j)
        vbh(nbh,j)=v(i,j)
        enddo
        enddo
7000    continue
	n=nup
	nbh=nbh1   ! 261103
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
c241103
c	if(nout.eq.1 .or. iii.eq.1 .or. iii.eq.neve)then
c	write(9,*)'be hadcas iii,n,nbh=',iii,n,nbh   ! sa
c	write(22,*)'be hadcas iii=',iii   !sa
c	call prt_sbh(nbh,cc)
c	call prt_sa1_h(nn)
c	endif  
        time_had=0.d0   ! 081010
csa     print*,'be. hadcas iii,n,nbh,nbe=',iii,n,nbh,nbe   !sa
	call hadcas(iii,neve,nout,time_had,ijkk)   ! 241103
csa     print*,'af. hadcas iii,n,nbh,nbe,ijkk=',iii,n,nbh,nbe,ijkk
	if(ijkk.eq.1)then   ! 161203   
c110603	iii=iii-1   ! it has been executed in 'scat' in 'hadcas' 
	goto 300   ! give up current event avoiding infinite collision loop
	endif
c241103
c	if(nout.eq.1 .or. iii.eq.1 .or. iii.eq.neve)then
c	write(22,*)'af hadcas iii=',iii   !sa
c	write(9,*)'af hadcas iii=',iii   ! sa
c	call prt_sbh(nbh,cc)
c	call prt_sa1_h(nn)
c	endif  
c201203
c	'sa1_h' to 'pyjets'
c	call trans_h
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
c201203
c241103
c250209
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
c250209

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
	do i=n+1,kszj   ! 261103
	do j=1,5
	k(i,j)=0
	p(i,j)=0.
	v(i,j)=0.
	enddo
	enddo
9000	continue
csa     print*,'af. sbh->pyjets iii,ijkk,n,nbe=',iii,ijkk,n,nbe
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

c       perform particle, declared unstable in the 'mdcy' array, decay
c        and remove hadronic decay gammas ('22') from 'pyjets' to 'sgam'
c130205	call pyexec
	rrp=1.16   ! 130205
csa     print*,'be. decayh iii,n,nbe=',iii,n,nbe   
	call decayh(rrp)   ! 130205
csa     print*,'af. decayh iii,n,nbe=',iii,n,nbe   ! sa
c181003
888	continue
c       write(22,*)'af. hadr. rescattering iii,n,nbe=',iii,n,nbe   ! bh
c       call pylist(1) 
c       call prt_sbe(nbe,cc)
c       call prt_sbh(nbh,cc)
c	write(9,*)'paciae, nnstop,zstop=',nnstop,zstop
        if(nnstop.ne.0)then
        sstop1=zstop/dfloat(nnstop)
        sstop=sstop+sstop1
        else
        nzstop=nzstop+1
        endif
c181003
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

c       write(22,*)'be. 80004 iii,nbe=',iii,nbe   ! bh

	if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.   ! 180921 yan
     c   itden.eq.2) .or. ipden.ge.11)goto 80004   ! 280718, 180921 yan
c291218        
	if(adj140.ne.4)goto 80004   
        if(adj140.eq.4)then 
        if(((ipden.eq.0 .and. itden.eq.1) .or. 
     c   (ipden.eq.1 .and. itden.eq.0)) .and. nbe.le.3)goto 80004   ! pA or Ap
c       if(ipden.eq.1 .and. itden.eq.1)goto 80004 ! AB
        if((ipden.eq.1 .and. itden.eq.1) .and. nbe.le.18)goto 80004 ! AB
        endif
c       'goto 80004': not fragments rest partons
c291218
c       write(22,*)'be. hadro. rest partons,iii,nbe,naff=',iii,nbe,naff   ! bh
c       call prt_aaff(naff,cc)
c       call prt_sbe(nbe,cc)
c	'pyjets' to 'aff' (note, 'pyjets' here is result of simulation above)
c180720
        naff=0
        call tran_pyjets

c       'trs' to 'sbe'
        if(iparres.eq.1 .and. ntrs.gt.0)then
        do ii1=1,ntrs
        ii3=nbe+ii1
        do ii2=1,5
        kbe(ii3,ii2)=ktrs(ii1,ii2)
        pbe(ii3,ii2)=ptrs(ii1,ii2)
        vbe(ii3,ii2)=vtrs(ii1,ii2)
        enddo
        enddo
        nbe=nbe+ntrs        
        endif

c050920 breaks up diquarls in 'sbe'
        if(adj12.eq.1.or.(adj12.eq.0.and.(nreac(4).ne.0.or.nreac(6)
     c   .ne.0.or.(nreac(7).ne.0.and.mstj1_2.eq.1))))call break_sbe        

c       'sbe' to 'sbh'   ! 260720
        do i1=1,nbe
        do i2=1,5
        kbh(i1,i2)=kbe(i1,i2)
        pbh(i1,i2)=pbe(i1,i2)
        vbh(i1,i2)=vbe(i1,i2)
        enddo
        enddo
        nbh=nbe   ! 260720


c       hadronizes rest partons in 'sbh' by 'py2ent'
        if(adj18.eq.0.and.nbh.gt.0)call pa2ent   ! 050920
        goto 80006
       
c       hadronizes rest partons in 'sbh' by coalescence        
        if(adj18.eq.1.and.nbh.gt.0)then
        ithroq=0
        ithrob=0
        ithroc=0
        do i=1,4
        throe(i)=0.
        enddo
c	'sbh' to 'parlist' for hadronizes rest partons by coalescence
        iprl=nbh
        do i=1,nbh
        idp(i)=kbh(i,2)
        rp(4,i)=0.
        eee=pbh(i,4)
        pp(4,i)=eee
        if(eee.lt.1.e-12)eee=1.e-12   ! 180420
        do j=1,3
        rp(j,i)=vbh(i,j)
        ppp=pbh(i,j)
        pp(j,i)=ppp
        vp(j,i)=ppp/eee
        enddo
        rmp(i)=pbh(i,5)
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
csa     print*,'be. coales (rest), iii,nbe,n,naff',iii,nbe,n,naff   
        call coales(1,neve,nout,0,0,0,0)   ! 060119 080119
c100119 seting first component equal to 1 in order to have 'call tabhb'
c100119 seting last four components equal to 0 in order to have netba=0
csa     print*,'af. coales (rest), iii,nbe,n,naff',iii,nbe,n,naff  
        endif
80006   continue   ! 120720
c       write(9,*)'af. hadro. rest partons iii,naff,n',iii,naff,n   ! sa

c	'aff' to 'pyjets'
        do ii1=1,naff
	ii3=n+ii1
        do ii2=1,5
        k(ii3,ii2)=kaff(ii1,ii2)
        p(ii3,ii2)=paff(ii1,ii2)
        v(ii3,ii2)=vaff(ii1,ii2)
        enddo
        enddo
	n=n+naff
        do ii1=n+1,kszj
        do ii2=1,5
        k(ii1,ii2)=0
        p(ii1,ii2)=0.
        v(ii1,ii2)=0.
        enddo
        enddo
c       moves partons from 'pyjets' to 'sbe'
        call remop
c       write(22,*)'af. hadro. rest partons iii,n,nbe',iii,n,nbe   ! bh
c       call prt_pyj(n,cc)
c       call prt_sbe(nbe,cc)


c	pythia output 
80004	continue
        if(nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
     c   .eq.neve)then
        write(mstu(11),*)'event=',iii
        call pylist(1)
        write(mstu(11),*)'ppsa=',(ppsa(i1),i1=1,5)
        write(22,*)'throe_p=',throe_p
        write(22,*)'ithroq_p,ithrob_p,ich_p=',ithroq_p,ithrob_p,
     c   dfloat(ich_p)/3.
        write(22,*)'throe=',throe
        write(22,*)'ithroq,ithrob,ithroc=',ithroq,ithrob,
     c   dfloat(ithroc)/3.
c       call prt_sbh(nbh,cc)   
c       call prt_aaff(naff)   ! 010518
c240720 if((ipden.eq.0 .and. itden.eq.1) .or.
c240720 c   (ipden.eq.1 .and. itden.eq.0))then
        write(22,*)'summary of rest partons nbe=',nbe   ! 260720
        call prt_sbe(nbe,cc)   ! 291218
c240720 endif
c170720 if(ipden.eq.1 .and. itden.eq.1)then
c       write(22,*)'summary and list of rest partons'
c       call prt_sbe(nbe,cc)
c       endif
        write(22,*)'summary and list of gammas'
        call prt_sgam(ngam,egam,8)   ! 250209 080419 1609019
        write(9,*)'# of parton scaterring process undergone in an event'
        write(9,108)iii,(nreac(i1),i1=1,9)   
        endif
108     format(10(1x,i4))
c       ??????????????? oscar stander output ??????????????????????
        if((nosc.eq.1.and.(mod(iii,nout).eq.0)).or.nosc.eq.2)
     c	 call oscar(iii)   ! 160705   
c       ??????????????? oscar stander output ??????????????????????


c	analyse an event
600	format(25(1x,i2)/)
	if(mod(iii,nout).eq.0)then
	print*,'event=',iii
	endif

c       records the event for off-line analyse
c       records the hadrons
c       write(99,190)bp,n
c       do 400 j=1,n
c       p1=p(j,1)
c       p2=p(j,2)
c       p3=p(j,3)
c       p4=p(j,4)
c       write(99,191) ik, p1, p2, p3, p4
c400    continue
c       records the gammas
c	write(98,190)bp,ngam   
c	do j=1,ngam
c	ik=kgam(j,2)
c	p1=pgam(j,1)
c	p2=pgam(j,2)
c	p3=pgam(j,3)
c	p4=pgam(j,4)
c	write(98,191) ik, p1, p2, p3, p4
c	enddo
c190    format(f10.4,I7)   ! sa
c191    format(I6,4(1x,e15.7))   ! sa

c       analyses the event on-line
        call analy(nmin,nminf,ncha,nchaf)   ! 281219
c       sum over events        
        dnmin=dnmin+nmin
        dnminf=dnminf+nminf
        dncha=dncha+ncha
        dnchaf=dnchaf+nchaf
	do kk=1,ispmax
	sbn(kk)=sbn(kk)+bn(kk)
	sbnf(kk)=sbnf(kk)+bnf(kk)
	do i1=1,40   ! 070419
	do i2=1,isdmax
	san(i1,i2,kk)=san(i1,i2,kk)+an(i1,i2,kk)
	sanf(i1,i2,kk)=sanf(i1,i2,kk)+anf(i1,i2,kk)
	enddo
	enddo
	enddo
c120119        
        stime_ini=stime_ini+time_ini
        stime_par=stime_par+time_par
        stime_had=stime_had+time_had
c120119   
        segam=segam+egam   ! 080419        
        segam1=segam1+egam1   ! 080419        
        segam2=segam2+egam2   ! 080419        
        segam3=segam3+egam3   ! 080419        
c070417
        snnc=snnc+nnc
        sadiv=sadiv+adiv
        sgpmax=sgpmax+gpmax
	do i1=1,6
        skapa(i1)=skapa(i1)+akapa(i1)
	enddo
        sgtime=sgtime+gtime
        sitime=sitime+itime
c070417
	spathn=spathn+pathn   ! 140219
	sevbin=sevbin+evbin   ! 260219
	sel=sel+noel
	do i1=1,600
	sinel=sinel+noinel(i1)
	dinel(i1)=dinel(i1)+noinel(i1)
        einel(i1)=einel(i1)+npinel(i1)   ! 140820
	enddo
        swoun=swoun+woun   ! 140820
        snpctl0=snpctl0+npctl0   ! 170121
        snpctlm=snpctlm+npctlm   ! 180121
        snpar=snpar+ncpart   ! 280722
c       print*,'ncpart,snpar=',iii,ncpart,snpar   ! 280722
c140820 swoun: # of wounded nucleons sumed up over enents
c170121 snpctl0: # of nn collision pairs (in parini.f) sumed up over enents 
c280722 snpar: # of collided nucleons (in parini.f) sumed up over enents
c071103
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
	sthroq=sthroq+ithroq+ithroq_p
	sthrob=sthrob+ithrob+ithrob_p
        sthroc=sthroc+ithroc+ich_p
	do i1=1,4
	sthroe(i1)=sthroe(i1)+throe(i1)+throe_p(i1)
	enddo

        open(8,file='nout.out',status='unknown')
        write(8,*)'iii=',iii
        close(8)
c       print*,'be. internal and final printing iii,n=',iii,n
c       internal and final printing and controled return
        if(mod(iii,nout).eq.0 .or. iii.eq.neve)then
        open(10,file='rms.out',status='unknown')
	flaa=dfloat(iii-ich)   ! July/20/98
	if(flaa.le.1.e-20)goto 1200
c       averaged ovr events     
        if(psno.eq.2)then
        averbo=averb/flaa
        psnono=psnon/flaa
        psnopo=psnop/flaa
        psnoto=psnot/flaa
        endif
	dnmino=dnmin/flaa
	dnminfo=dnminf/flaa
        dnchao=dncha/flaa
        dnchafo=dnchaf/flaa
c010220
	do kk=1,ispmax
	sbo(kk)=sbn(kk)/flaa
	sbof(kk)=sbnf(kk)/flaa
	do i1=1,40   ! 070419
	do i2=1,isdmax
	sao(i1,i2,kk)=san(i1,i2,kk)/flaa
	saof(i1,i2,kk)=sanf(i1,i2,kk)/flaa
	enddo
	enddo
	enddo
        segamo=segam/flaa   ! 080419
        segam1o=segam1/flaa   ! 080419
        segam2o=segam2/flaa   ! 080419
        segam3o=segam3/flaa   ! 080419
c120119
        stime_ini=stime_ini/flaa
        stime_par=stime_par/flaa
        stime_had=stime_had/flaa
c120119
c070417
        snnco=snnc/flaa
        sgtimeo=sgtime/flaa
c       sgtimeo: average number of gluons in a string
        do i1=1,6   ! 010518
        skapao(i1)=skapa(i1)/flaa
        enddo
        sadivo=sadiv/flaa
        sgpmaxo=sgpmax/flaa
        sitimeo=sitime/flaa
c070417
c061103
	spathni=spathn/flaa   ! 140219
	sevbini=sevbin/flaa   ! 260219
c061103
c071103
	reli=nrel/flaa   ! original =rel/flaa 220110
	rineli=rinel/flaa
c071103
	seli=sel/flaa
	sineli=sinel/flaa
	do i1=1,600
	dineli(i1)=dinel(i1)/flaa
        eineli(i1)=einel(i1)/flaa   ! 140820
	enddo
        swouni=swoun/flaa   ! 140820
        snpctl0i=snpctl0/flaa   ! 170121
        snpctlmi=snpctlm/flaa   ! 180121
        snpari=snpar/flaa   ! 280722
c140820 swouni: # of wounded nucleons averaged over enents after parini.f
c200601
c170121 snpctl0i: # of nn collision pairs averaged over enents 
c180121 snpctlmi: maximum # of nn collision pairs averaged over enents 
c280722 snpari: MC Glauber-like <N_part>
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
c220110
        do i1=1,9
        nrea=nrea+nreac(i1)
        snreac(i1)=nreac(i1)/flaa
        enddo
        srea=float(nrea)/flaa
c220110
1200	continue
c       user out put        
	write(10,*)'parp81,bp,mstp82=',parp(81),bp,mstp(82)    ! 291207  
        write(10,*)'MC Glauber-like <N_coll>,<N_part> =',
     c   snpctl0i,snpari   ! 180121 280722
        write(10,*)'largest ave. # of nn collision pairs =',snpctlmi ! 280722
c140820
        write(10,*)'ave. # of nn collision pairs calling pythia, not 
     c   calling pythia=',eineli(592),eineli(593)
        write(10,*)'ave. # of wounded nucleons in parini.f =',swouni  
c140820 

	write(10,*)'colli. # suffered by projectile nucleon 
     c	 in target nucleus',spathni   ! 140219
	write(10,*)'event averaged N_bin',sevbini   ! 260219
        write(10,*)'event averaged energy of gamma after 
     c   partonic initiation, partonic cascade, hadronization 
     c   and end of event=',segam1o,segam2o,segam3o,segamo ! 080419

c071103
	write(10,*)'# of successful and blocked collision in parton 
     c	 cascade=',rineli,reli,reli+rineli
c071103
c220110
        write(10,*)'average collision # in parton cascade=',srea
        write(10,*)'total # of scaterring processes in parton cascade'
        write(10,*)(snreac(i1),i1=1,9)
c220110

	write(10,*)'average frequency of the occurring of each inela. 
     c   in hadron cascade='
	write(10,*)dineli
	write(10,*)'el. and inela. colli. # and sum in hadron cascade=',
     c	 seli,sineli,seli+sineli

c200601
	write(10,*)'(Npart)mini-jet,Nnn,Npp=',skparo,sknno,skppo
	write(10,*)'Nnp,Ntot,Nep=',sknpo,sknno+skppo+sknpo,skepo   ! 060813
c200601
        if(psno.eq.2)write(10,*)'psno, ave. b,N_part and N_bin=',
     c   psno,averbo,psnopo,psnoto,psnono*csnn   ! 280113
	if(ipden.ge.11.and.ipden.le.16)
     c   write(10,*)'event average number of lepton studied=',vnlep/flaa !260314
c070417
        write(10,*)'default par1,par2,par3,parj4,par21=',
     c   parj1,parj2,parj3,parj4,parj21   ! 010518
        write(10,*)'parj1,parj2,parj3,parj4,parj21,keff=',
     c   (skapao(i1),i1=1,6)   ! 010518
        write(10,*)'averaged # of gluon in a string',sgtimeo
        write(10,*)'event averaged value of the factor related to # of'
        write(10,*)'gluons and hardest gluon in a string,event averaged'
        write(10,*)'transverse momentum of hardest gluon,event averaged'
        write(10,*)'# strings=',sadivo,sgpmaxo,sitimeo
c070417
c120119
        write(10,*)'times & sum=',stime_ini,stime_par,stime_had,
     c   stime_ini+stime_par+stime_had        
c120119        
        write(10,*)'multiplicity of negative particles=',dnmino
	write(10,*)'multiplicity of negative particles=',dnminfo
	write(10,*)'multiplicity of positive particles,partial=',dnchao
	write(10,*)'multiplicity of positive particles,full=',dnchafo
	write(10,*)'throw away ithroq,ithrob,ithroc=',
     c	 wthroq,wthrob,wthroc/3.
	write(10,*)'throe=',wthroe
	write(10,*)'avb,avneu,astbp,astbt,aanbin=',
     c   avb,avneu,astbp,astbt,aanbin   ! 280113
        write(10,*)'particle multiplicity=',(sbof(ll),ll=1,ispmax)
	write(10,*)'particle multiplicity=',(sbo(ll),ll=1,ispmax)
601	format(8(1x,f7.4))
csa****************************************************************

	do m2=1,isdmax
	write(10,*)'ID of distribution m2=',m2
	do m3=1,ispmax
	write(10,*)'distribution belong to m3=',m3
	write(10,*)(sao(m1,m2,m3),m1=1,40)   ! 070419
	write(10,*)(saof(m1,m2,m3),m1=1,40)   ! 070419
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
        do i1=1,40   ! 070419
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
        write(10,*)(sao(m1,m2,m3),m1=1,40)   ! 070419  
	write(10,*)(saof(m1,m2,m3),m1=1,40)   ! 070419
        enddo
        enddo
	endif
c260314

	close(10)
	endif

1000	if(iii.lt.neve)then
c260718	if(dabs(bmin-bmax).lt.10d-4)goto 300
	if(psno.eq.0 .or. psno.eq.2)goto 300   ! 280113 260718
	if(psno.eq.1 .and. jjj.eq.10)then   ! 260718
c	10: total number of impact paremeters in systematic sampling for impact
c           parameter
	jjj=1
	goto 300
	elseif(psno.eq.1 .and. jjj.lt.10)then   ! 260718
	jjj=jjj+1
	goto 300
	else   ! 260718
	endif   ! 260718
	endif
	
	write(9,*)'nncoll=',nncoll   ! sa 060814
c060813 statistics of processes generated
        call pystat(1)   ! 060813
	close(2)
	close(3)
	close(5)
	close(9)
	close(22)
	close(34)
c	close(98)   ! 260219
c       close(99)
900 	format(i5,8(1x,f9.3),1x,f5.3)
c	timeb=dtime(ty)
c	write(9,*)'time consuming =',timeb
	stop
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remop   
c	moves q,qbar,g,diquark, and anti-diquark from 'pyjets' to 'sbe'  
      PARAMETER (KSZJ=80000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/pyjets/n,nonj,k(kszj,5),p(kszj,5),v(kszj,5)   
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        jb=0
201     do i1=jb+1,n   ! i1 loop
        kf=k(i1,2)
        kfab=iabs(kf)
        if(kfab.gt.6 .and. kfab.ne.2101 .and. kfab.ne.3101
     c   .and. kfab.ne.3201 .and. kfab.ne.1103 .and. kfab.ne.2103
     c   .and. kfab.ne.2203 .and. kfab.ne.3103 .and. kfab.ne.3203
     c   .and. kfab.ne.3303 .and. kfab.ne.21)then
        jb=jb+1
        goto 202
        endif
        nbe=nbe+1
        do i2=1,5
        kbe(nbe,i2)=k(i1,i2)
        pbe(nbe,i2)=p(i1,i2)
        vbe(nbe,i2)=v(i1,i2)
        enddo
        if(i1.eq.n)then
        n=n-1
        goto 203
        endif
c	moves particle list 'pyjets' one step downward from i1+1 to n
        do j=i1+1,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 201
202     enddo   ! i1 loop

203     continue
	return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pa2ent
c       hadronizes with 'py2ent'
        parameter(kszj=80000)
        implicit double precision (a-h,o-z)
        implicit integer (i-n)
        integer pyk,pychge,pycomp
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/aaff/naff,nonaf,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
c       not write header 
        mstp(127)=0
        mstp(122)=0
        mstu(11)=22

        nbe=0
        n=0  
c       making 'sbh' in order of g,qbar, & q        
        jj=0
        do i1=1,nbh
        kff=kbh(i1,2)
        if(kff.eq.21)then
        ii=i1
        jj=jj+1
c       moves 'ii' to position of 'jj' in 'sbh'
        kip1=kbh(jj,1)
        kip2=kbh(jj,2)
        kip3=kbh(jj,3)
        kip4=kbh(jj,4)
        kip5=kbh(jj,5)
        pip1=pbh(jj,1)
        pip2=pbh(jj,2)
        pip3=pbh(jj,3)
        pip4=pbh(jj,4)
        pip5=pbh(jj,5)
        vip1=vbh(jj,1)
        vip2=vbh(jj,2)
        vip3=vbh(jj,3)
        vip4=vbh(jj,4)
        vip5=vbh(jj,5)

        do i2=1,5
        kbh(jj,i2)=kbh(ii,i2)
        pbh(jj,i2)=pbh(ii,i2)
        vbh(jj,i2)=vbh(ii,i2)
        enddo

        kbh(ii,1)=kip1
        kbh(ii,2)=kip2
        kbh(ii,3)=kip3
        kbh(ii,4)=kip4
        kbh(ii,5)=kip5
        pbh(ii,1)=pip1
        pbh(ii,2)=pip2
        pbh(ii,3)=pip3
        pbh(ii,4)=pip4
        pbh(ii,5)=pip5
        vbh(ii,1)=vip1
        vbh(ii,2)=vip2
        vbh(ii,3)=vip3
        vbh(ii,4)=vip4
        vbh(ii,5)=vip5
        endif
        enddo
        ngg=i1   ! order # of last gluon in 'sbh'

        if(nbh.gt.ngg)then   ! upto here jj=ngg

        do i2=ngg+1,nbh  
        kff=kbh(i2,2)
        if(kff.lt.0.)then   ! diquark
        ii=i2
        jj=jj+1
c       moves 'ii' to position of 'jj' in 'sbh'
        kip1=kbh(jj,1)
        kip2=kbh(jj,2)
        kip3=kbh(jj,3)
        kip4=kbh(jj,4)
        kip5=kbh(jj,5)
        pip1=pbh(jj,1)
        pip2=pbh(jj,2)
        pip3=pbh(jj,3)
        pip4=pbh(jj,4)
        pip5=pbh(jj,5)
        vip1=vbh(jj,1)
        vip2=vbh(jj,2)
        vip3=vbh(jj,3)
        vip4=vbh(jj,4)
        vip5=vbh(jj,5)
        do i3=1,5
        kbh(jj,i3)=kbh(ii,i3)
        pbh(jj,i3)=pbh(ii,i3)
        vbh(jj,i3)=vbh(ii,i3)
        enddo
        kbh(ii,1)=kip1
        kbh(ii,2)=kip2
        kbh(ii,3)=kip3
        kbh(ii,4)=kip4
        kbh(ii,5)=kip5
        pbh(ii,1)=pip1
        pbh(ii,2)=pip2
        pbh(ii,3)=pip3
        pbh(ii,4)=pip4
        pbh(ii,5)=pip5
        vbh(ii,1)=vip1
        vbh(ii,2)=vip2
        vbh(ii,3)=vip3
        vbh(ii,4)=vip4
        vbh(ii,5)=vip5
        endif
        enddo
        nqb=i2   ! order # of last diquark

        endif
c       making 'sbh' in order of g, qbar, & q finished

c       fragmenting 'sbh'

c       fragmenting gg in 'sbh'
        nggg=0
        if(ngg.gt.0)then   ! 1
        do i1=1,ngg,2   ! 2
        kf1=kbh(i1,2)
        er1=pbh(i1,4)
        kf2=kbh(i1+1,2)
        er2=pbh(i1+1,4)
        ss=er1+er2
        call py2ent(0,kf1,kf2,ss)
        if(ipden.lt.11)call pyedit(2)
        if(ipden.ge.11)call pyedit(1)
c       moves "77" from 'pyjets' to 'sgam'
        if(n.gt.0)then
        n77=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=77   ! '77': photons after hadronization of current string
        n77=n77+1
        endif
        enddo
        if(n77.gt.0)call remo_gam_hadro(77)
        endif
c       'pyjets' to 'aff'
        call tran_pyjets
        nggg=nggg+2        
        enddo   ! 2
        if(nggg.lt.ngg)then
c       moves last g to 'sbe'
        nbe=nbe+1
        do jj=1,5
        kbe(nbe,jj)=kbh(ngg,jj)
        pbe(nbe,jj)=pbh(ngg,jj)
        vbe(nbe,jj)=vbh(ngg,jj)
        enddo
        endif

        endif   ! 1   

c       fragmenting qbarq in 'sbh'        
100     continue
        if(nqb.gt.ngg)then   ! 3 having antiquark
        kf1=kbh(ngg+1,2)
        er1=pbh(ngg+1,4)
        if(nbh.gt.nqb)then   ! 4 having quark
        kf2=kbh(nqb+1,2)
        er2=pbh(nqb+1,4)
        ss=er1+er2
        call py2ent(0,kf1,kf2,ss)
        if(ipden.lt.11)call pyedit(2)
        if(ipden.ge.11)call pyedit(1)
c       moves "77" from 'pyjets' to 'sgam'
        if(n.gt.0)then
        n77=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=77   ! '77': photons after hadronization of current string
        n77=n77+1
        endif
        enddo
        if(n77.gt.0)call remo_gam_hadro(77)
        endif
c       'pyjets' to 'aff'
        call tran_pyjets
        n=0
c       moves 'sbh' list one step downward since nqb+1 to nbh 
        do i1=nqb+1,nbh
        i2=i1-1
        do i3=1,5
        kbh(i2,i3)=kbh(i1,i3)
        kbh(i2,i3)=kbh(i1,i3)
        kbh(i2,i3)=kbh(i1,i3)
        enddo
        enddo
        nbh=nbh-1
c       moves 'sbh' list one step downward since ngg+1 to nbh
        do i1=ngg+1,nbh
        i2=i1-1
        do i3=1,5
        kbh(i2,i3)=kbh(i1,i3)
        kbh(i2,i3)=kbh(i1,i3)
        kbh(i2,i3)=kbh(i1,i3)
        enddo
        enddo
        nqb=nqb-1
        nbh=nbh-1
        n=0
        endif   ! 4
        endif   ! 3        
        if(nqb.gt.ngg.and.nbh.gt.nqb)goto 100
c       moves rest qbar and/or q to 'sbe'
        if(nbh.gt.0)then   
        do i1=ngg+1,nbh
        nbe=nbe+1
        do jj=1,5
        kbe(nbe,jj)=kbh(i1,jj)
        pbe(nbe,jj)=pbh(i1,jj)
        vbe(nbe,jj)=vbh(i1,jj)
        enddo
        enddo
        endif
c       fragmenting 'sbh' finished

        return
        end



c********************************************************************
        subroutine tran_pyjets
c       'pyjets' to 'aff'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (kszj=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/aaff/naff,nonaf,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)   
        do l=1,n
        l1=naff+l
        do m=1,5
        kaff(l1,m)=k(l,m)
        paff(l1,m)=p(l,m)
        vaff(l1,m)=v(l,m)
        enddo
        enddo
        naff=naff+n
	do l=naff+1,kszj
	do m=1,5
	kaff(l,m)=0
	paff(l,m)=0.
	vaff(l,m)=0.
	enddo
	enddo
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine oscar(iii)
c       record history of spatial and momentum coordinates due to
c        OSC1999A
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=80000)
        common/pyjets/nsa,npad,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
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


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine remo_gam_sbh(ii)   ! 250209
c	move particles with flavor code ii ('44') from  'sbh' to 'sgam'
	parameter (kszj=80000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	COMMON/sbh/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine remo_gam_par(ii)   ! 250209
c	move particles with flavor code ii ('55') from 'pyjets' to 'sgam'
	parameter (kszj=80000)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine remo_gam(ii)   ! 250209
c	move particles with flavor code ii ('66') from  'pyjets' to 'sgam'
	parameter (kszj=80000)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine remo_gam_hadro(ii)   ! 240219
c	move particles with flavor code ii ('77') from  'pyjets' to 'sgam'
	parameter (kszj=80000)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sgam(nn,egam,iprt)   ! 250209 080419 160919
c       print particle list 'sgam' and sum of momentum and energy
        parameter (kszj=80000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc   ! 260419
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        dimension peo(4)
        call psum(pgam,1,ngam,peo)
        egam=peo(4)   ! 080419
        ch1=0.
        do i1=1,nn
        kf=kgam(i1,2)
        if(kf.eq.22 .or. kf.eq.44 .or. kf.eq.55 .or. kf.eq.66 
     c   .or. kf.eq.77)goto 100
        ch1=ch1+pychge(kf)
100     enddo
        if((nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
     c   .eq.neve) .and. iprt.eq.8)then   ! 260419 160919
        write(22,*)'c & p sum=',ch1/3.,peo   !
        do i=1,nn
        write(22,*)i,kgam(i,1),kgam(i,2),(pgam(i,j),j=1,4)
        enddo
        endif   ! 260419
        return
        end

                                                                        

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_aaff(nn)
c       print particle list and sum of momentum and energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/aaff/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
c       do i=1,nn
c       write(22,*)i,kbh(i,2),(pbh(i,j),j=1,4)
c       enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3.
        write(22,*)'aaff nn=',nn
        write(22,*)'c & p sum=',cc,peo   !
        return
        end


 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break_sbe   ! 030920
c       breaks up diquark (anti-diquark), gives four momenta 
c	 and four positions to the broken partons
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/sbe/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa24/adj1(40),nnstop,non24,zstop   
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
	jb=0
	ii=idio   

100     do i1=jb+1,n
        kf=k(i1,2)
	kfab=iabs(kf)
        if(kfab.ne.2101 .and. kfab.ne.3101
     c   .and. kfab.ne.3201 .and. kfab.ne.1103 .and. kfab.ne.2103
     c   .and. kfab.ne.2203 .and. kfab.ne.3103 .and. kfab.ne.3203
     c   .and. kfab.ne.3303)then
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
c       give four momentum to the breaked quarks
	call bream_sbe(i1,kf1,kf2)
c       give four coordinate to the breaked quarks
        call coord_sbe(i1)
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
        subroutine bream_sbe(ii,kf1,kf2)
c       give four momentum to the broken quarks
c       ii: line number of diquark in 'sbe'
c       kf1,kf2: flavor codes of broken quarks
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/sbe/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)   ! 260503
        am1=pymass(kf1)
        am2=pymass(kf2)
        pp(1,5)=am1
        pp(2,5)=am2
c       pp : four momentum & mass of broken quarks, local variable 
        do i1=1,4
        ps(i1)=p(ii,i1)
        enddo
c       ps : four momentum of diquark, local variable 
	goto 400   ! activate it for 'decay method'
c       broken quarks share diquark four momentum randomly,
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
	call decmom_sbe(ps,pp,am1,am2,decsuc)
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
	subroutine decmom_sbe(ps,pp,am1,am2,decsuc)
c	calculate four momentum of decayed particles
c	ps: four momentum of decaying particle
c	am1 and am2: mass of decayed pair
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000)
      COMMON/sbe/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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
        subroutine coord_sbe(ii)
c       give four position to broken quarks
c       first broken quark takes the four position of diquark
c       second broken quark is arranged around first ones within
c        0.5 fm randumly in each of three position coordinates and has same
c        fourth position coordinate as diquark
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/sbe/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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



*******************************************************************************
c This code calculates the nuclear overlap functions which are needed 
c to scale pp to pA and AB.
c Units are fm
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
     c  part1(200),part2(200),binn(200)   ! 020511 020718
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
c Calculate npart and Nbin using thickness functions
        do i=1,200
        bbb=bb(i)
        if(kjp23.eq.2)then
        call PART(A1,A2,bbb,sigma_NN,denflag,density0,nshot,
     c  part1(i),part2(i))
	binn(i)=(sigma_NN/10.)*TA1A2(i)   ! 020718 180219
c        if(i.eq.36)write(9,*)'i,bbb,pir,tir=',i,bbb,part1(i),part2(i),
c     c	 binn(i)   ! 020718
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
      if(kjp23.eq.2)write(9,901)'i','b','TA','TB','TAB','Apart','Bpart',
     c	 'Nbin'
        if(kjp23.eq.1)write(9,902)'i','b','Apart','Bpart'
      do i=1,200
        if(kjp23.eq.2)write(9,905)i,bb(i),TA1(i),TA2(i),TA1A2(i),
     c  part1(i),part2(i),binn(i)   ! 020718
        if(kjp23.eq.1)write(9,906)i,bb(i),part1(i),part2(i)
      enddo
 901  format(2a6,6a10)   ! 020718
 902    format(2a6,2a10)
 905  format(i5,2x,f6.2,6(f10.3))   ! 020718
 906  format(i5,2x,f6.2,2(f10.3))
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
     c  part1(200),part2(200),binn(200)   ! 020511 020718
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
     c  part1(200),part2(200),binn(200)   ! 020718
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
     c  part1(200),part2(200),binn(200)   ! 020511 020718
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
      REAL ABI, A2B2, GAMMA, PI, DSQRT, AB   ! 081010                                     
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
     X BETA + 1.0E0) / GAMMA(ABI)   ! 081010                                              
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
   60 MUZERO = GAMMA(ALPHA + 1.0E0)   ! 081010                                            
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
      FUNCTION GAMMA(X)   ! 081010                                                      
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
    7 GAMMA=VALUE   ! 081010                                                           
      RETURN                                                                    
      END 
