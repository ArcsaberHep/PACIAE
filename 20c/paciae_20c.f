        program paciae   ! 020708 241110
c	parton and hadron cascade model for relativistic nucleus-nucleus
c	 collision 
c	it is composed of paciae_20c.f, parini_20c.f, parcas_20c.f, sfm_20c.f,  
c	 coales_20c.f, hadcas_20c.f, and p20c.f
c       paciae_20c.f: user program
c       parini_20c.f: generate a partonic initial state for a
c        nucleus-nucleus collision
c       parcas_20c.f: perform parton rescattering, where only 2->2 processes
c        are considered and LO pQCD cross section or its regularized
c        approximation is used
c       sfm_20c.f: hadronization according to LUND string fragmentation model
c       coales_20c.f: hadronization according to Monte Carlo coalescence model
c       select sfm_20c.f or coales_20c.f by parameter adj1(12) 
c       hadcas_20c.f: perform hadronic rescattering
c       read ``paciae_guide" for the details
c	note: the statistics made here is either for parton or   
c	 hadron according to purpose
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000,mplis=40000,KSZ1=30)
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
     &  iabsb,iabsm,non10,ajpsi,csspn,csspm
	common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
        common/sa16/dtt,dni(10),dpi(10),edi(10),bmin,bmax
     &   ,bar(10),abar(10),barf(10),abarf(10)   ! 033101
     &	 ,emin(10),eminf(10),eplu(10),epluf(10)   ! 033101
	common/sa18/tdh,itnum,non18,cptl,cptu,cptl2,cptu2,snum(4,20),
     &	 v1(4,20),v2(4,20),v12(4,20),v22(4,20)
	common/sa18_pt/snum_pt(21,20),v1_pt(21,20),v2_pt(21,20),
     c   v12_pt(21,20),v22_pt(21,20)   ! 280607
        common/sa18_eta/snum_eta(21,20),v1_eta(21,20),v2_eta(21,20),
     c   v12_eta(21,20),v22_eta(21,20)   ! 280607
	common/sa21/rots,degd   ! 120299
	common/sa23/kpar,knn,kpp,knp   ! 200601
	common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003
	common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 220110
        common/sa27/itime,kjp22,gtime,astr,akapa(5),parj1,parj2,parj3,
     c   parj21   ! 020708
	common/sa30/vneump,vneumt   ! 241110
        common/sa31/rmax,bbb(200),TA1(200),TA2(200),TA1A2(200),
     c  part1(200),part2(200)   ! 020511
        common/sa6_c/ithroq,ithrob,ithroc,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa6_t/ithroq_t,ithrob_t,ich_t,non6_t,throe_t(4)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)   ! 050603
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)   ! 050603
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	common/show/vip(mplis),xap(mplis)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	common/count/isinel(600)
	common/ctllist/npctl,npinel(600),npctl0,npel   ! 061103
        common/ctllist_p/nreac(9),nrel   ! 071103
	common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 250209
        common/ssin/nsin,nonsin,ksin(kszj,5),psin(kszj,5),vsin(kszj,5),
     c   bsin(kszj)   ! 250209
        dimension an(20,5,20),bn(20),san(20,5,20),sbn(20),
     c	 anf(20,5,20),bnf(20),sanf(20,5,20),sbnf(20)
	dimension sao(20,5,20),sbo(20),saof(20,5,20),sbof(20)
        dimension skapa(5),skapao(5),snreac(9)   ! 020708 220110
	dimension c(5),dinel(600),dineli(600),sthroe(4),wthroe(4)
	dimension bpp(20),kdiq(kszj,5),dgmas(kszj)
	dimension acoll(20),acollp(20),acollt(20)
	dimension pj(20),ej(20),pjp(20),ejp(20),ett(20),sbp(20) 
	dimension pji(20),eji(20),pjpi(20),ejpi(20),etti(20)
	dimension fpj(20),fej(20),fpjp(20),fejp(20),fett(20) 
	dimension fpji(20),feji(20),fpjpi(20),fejpi(20),fetti(20)
        dimension fettm(20),ettm(20),fettmi(20),ettmi(20)
	dimension sndi(10),spdi(10),sedi(10),sndio(10),spdio(10),
     c	 sedio(10)
c033101
	dimension sbar(10),sabar(10),sbarf(10),sabarf(10),
     c	 sbar2(10),sabar2(10),sbarf2(10),sabaf2(10) 
	dimension obar(10),oabar(10),obarf(10),oabarf(10),
     c	 obar2(10),oabar2(10),obarf2(10),oabaf2(10)
	dimension erat(10),eratf(10),erat2(10),eratf2(10),
     c	 orat(10),oratf(10),orat2(10),oratf2(10)
c033101 
     	dimension ssnum(4,20),sv1(4,20),sv2(4,20),
     c	 sv12(4,20),sv22(4,20),
     c	 sv1o(4,20),sv2o(4,20),sv12o(4,20),sv22o(4,20)
     	dimension ssnum_pt(21,20),sv1_pt(21,20),sv2_pt(21,20),
     c	 sv12_pt(21,20),sv22_pt(21,20),sv1o_pt(21,20),sv2o_pt(21,20),
     c   sv12o_pt(21,20),sv22o_pt(21,20)   ! 280607
     	dimension ssnum_eta(21,20),sv1_eta(21,20),sv2_eta(21,20),
     c	 sv12_eta(21,20),sv22_eta(21,20),sv1o_eta(21,20),sv2o_eta(21,20)
     c   ,sv12o_eta(21,20),sv22o_eta(21,20)   ! 280607
c021207
        dimension ssv1_pt(21,20),ssv2_pt(21,20),ssv12_pt(21,20),
     c   ssv22_pt(21,20),ssv1o_pt(21,20),ssv2o_pt(21,20),
     c   ssv12o_pt(21,20),ssv22o_pt(21,20)
        dimension ssv1_eta(21,20),ssv2_eta(21,20),ssv12_eta(21,20),
     c   ssv22_eta(21,20),ssv1o_eta(21,20),ssv2o_eta(21,20),
     c   ssv12o_eta(21,20),ssv22o_eta(21,20)
c021207
c	dimension uds(3),udsb(3),duds(3),dudsb(3),
c     c   uuds(3),uudsb(3)   ! 261002
        dimension fuds(3),fudsb(3),fduds(3),fdudsb(3)   ! 261002
        dimension avpt(20),savpt(20),avpto(20),navpt(20),nsavpt(20)  ! 170705
        dimension fvpt(20),sfvpt(20),fvpto(20),nfvpt(20),nsfvpt(20)  ! 170705
	dimension asdd(20),ptrus(3)   ! 300404 131204 130605
        dimension nreaco(9)   ! 220110
        real nmin,nminf,ncha,nchaf,nmine,nminef,nplue,npluef!020203
c210803 avpt: average transverse momentum of produced particle
c210803 navpt: number of given piece of particle in a event
c210803 
c	common block sa24: adjustable variables
c	adj1(i), i=
c	1: k factor used in parton cascade
c	2: parameter \alpha_s (as in program) in parton cascade
c	3: parameter tcut in program, to avoid divergence in calculating 
c          parton-parton differential cross section in parton cascade 
c	4: parameter idw, number of intervals in numerical integration in parton cascade
c	5: =1 and 0, with and without nuclear shadowing, respectively
c	6: parameter a, i.e. parj(41) in Lund string fragmentation function
c	   if adj1(12)=0
c	   parameter a in FF string fragmentation function if adj1(12)=1
c       7: parameter b, i.e. parj(42) in Lund string fragmentation function
c       8: mstp(82) in PYTHIA 6.4
c       9: parp(81), (D=1.4 GeV/c), effective minimum transverse momentum for 
c	   multiple interactions with mstp(82)=1                      
c	10: parp(31),k factor in pythia (parp(31))
c	11: time accuracy used in hadron cascade (time accuracy used in 
c           parton cascade is dddt)
c	12: model of hadronization: =0 string fragmentation; =1: coalescence
c	13: dimension of meson table considered
c       14: dimension of baryon table considered
c	15: string tension  
c	16: number of loops in deexcitation of energetic quark in coalescence
c	17: the threshold energy in deexcitation of energetic quark in coalescence
c	18: =0 and 1 without and with partonic Pauli in parton cascade, respectively
c	19: time accuracy used in parton cascade (dddt in program)
c	20: =0 exact pQCD parton-parton cross section
c	    =1 limited and regularized parton-parton cross section (B. Zhang)
c	    =2 the same as 0 but flat scattering angle distribution is assumed
c           =3 the same as 1 but flat scattering angle distribution is assumed
c	21: =0 and 1 without and with phase space adjudgment, respectively, 
c	 in coalescence
c	22: critical value of the product of radii both  in coordinate and 
c	 momentum phase space (4 is assumed) 
c	23: =0 LUND fragmentation function used in subroutine 'ffm' in coalescence
c	    =1 IF fragmentation function used 
c	24: the virtuality cut ('tl0') in time-like radiation in parton cascade
c	25: \Lambda_QCD in parton cascade
c	26: number of random number thrown away
c       27: largest momentum allowed for particle('dpmax')
c       28: largest position allowed for particle (drmax=para10*max(rnt,rnp), giving 1 
c           to it in usu.dat and calculating it in the running)  
c	29: width of two dimension Gaussian distribution sampling px and py of
c           produced quark pair in deexcitation of the energetic quark in coalescence
c	30: maximum $p_T^2$ in above two dimension Gaussian distribution 
c	31: parj(1) in pythia
c	32: parj(2) in pythia
c	33: parj(3) in pythia
c	34: parj(21) in pythia
c	35: mstp(91) in pythia,parton transverse momentum (k_{\perp}) distribution 
c           inside hadron; 
c	    =1, Gaussian; 
c           =2, exponential
c	36: =0 without phenomenological parton energy loss in parton cascade
c	    =1 with phenomenological parton energy loss
c	37: the coefficient ('c') in phenomenological parton energy loss
c	38: pt cut in phenomenological parton energy loss 
c	39: width of Gaussian k_{\perp} distribution in hadron if mstp(91)=1
c           width of exponential k_{\perp} distribution in hadron if mstp(91)=2
c	40: =1 transport processes endded after parton initiation
c	    =2 after parton rescattering
c	    =4 after hadron rescattering
c210803
c       para1_1: total cross section of nn, used in parton initiation
c        for nuclus-nucleus collision, this is irrelevant
c	 to the inelastic cross section of nn used in PYTHIA
c       para1_2: total cross section of nn, used in hadron cascade
c	dni: nucleon number density
c	dpi: particle (nucleon, pion, kaon ...) number density
c	edi: energy density of particle (nucleon, pion, kaon ...)
c	all of above three densities are calculated in r less or equal
c	 2 fm and t less or equal 10 fm/c
c	pj and ej: (pt)**2 and Et of J/psi 
c	pjp and ejp: (pt)**2 and Et of (J/psi)
c	ett: total Et 
c	acoll: array, the demension of which should be larger
c	 than or equal to 'nmax' 
c	note: the dimension of 'bpp' must be < or = nmax
c       ipden: =0,if projectile is proton
c              =1, projectile is nucleus
c	       =2, projectile is e+   ! 060605
c       itden: =0, if target is proton
c              =1, target is nucleus
c              =2, target is e-   ! 060605
c       suppm: the upper bound in sampling the radius of projectile nucleon
c       suptm: the upper bound in sampling the radius of target nucleon
c       suppc: the maximum radius in sampling for projectile
c       suptc: the maximum radius in sampling for target
c       r0p: projectile radius 
c       r0t: target radius 
c       pio: 3.1416
c	bp: impact parameter
c	iii: current run number
c	coor: position of CM
c	ispmax: maximum # of kinds of particles wanted to statistics
c	ispkf(i): flavor of i-th kind of particle wanted to statistics
c	kfmax: the maximum # of particles with different flavor considered 
c	kfaco(i): i-th particle flavor
c	numb(i): sum of particles up to the last one with flavor code of kfaco(i) 
c                in particle list
c	an(l,i,j) (san(l,i,j)):
c	 l: value of distribution argument (e. g. value of y or pt) 
c	 i: identify the distribution (e. g. i=2 is y distribution)
c	 j: order # of kf code of particle wanded to statistics (1 to ispmax)
c	bn,sbn,sbo: record the multiplicity of particles
c       anf(l,i,j) (sanf(l,i,j)), for instance, is corresponding to an(l,i,j) 
c        (san(l,i,j)) but is statistics of full phase space instead of partical
c        phase space
c	isdmax: maximum # of distributions wanted to calculate
c	asd(i): interval segmented for i-th distribution  
c	 i=1: for y
c	 i=2: for pt
c         .      . 
c         .      .
c         .      .
c	iflmax: maximum # of filters,=0 means no filter at all
c	afl(j,i,1): lower limit of i-th filter for the j-th particle
c	afl(j,i,2): upper limit of i-th filter for the j-th particle
c 	 i=1: y filter	
c	 i=2: pt filter
c         .        .
c         .        . 
c         .        .
c	ifram: = 0 for fixed target, = 1 for collider 
c	cspipi (fm^2): total cross section of pion + pion
c	sig (fm^2): cross section of pion + pion to kaon + kaon
c	cspin (fm^2): total cross section of pion + nucleon interaction
c	cskn (fm^2): total cross section of kaon + nucleon interaction
c	csnn (fm^2): total cross section of n + n interaction
c	rcsit: ratio of inelastic to total cross section
c	disbe(i,j): allowable minimum distance between two particles of
c        kfaco(i) & kfaco(j).
c	c17(i,1-3): position of particle i
c	tp(i): time of particle i
c	ishp(i): =1 if i-th particle inside the simulated volume
c	         =0 if i-th particle outside 
c	tau(i): formation time of particle i.
c	isinel(i): = 0 without i-th inelastic process
c	           = 1 with i-th inelastic process
c       nreac(i): statistics of successful i-th collision in parton cascade
c       nrel: statistics of blocked collisions in parton cascade
c	npel: statistics of bolcked nn collisions in parton initiation
c       npinel(600): statistics of successful nn collisions in parton initiation
c	noel : statistics of elastic collisions in hadron cascade
c       noinel(i): statistics the i-th inelastic channel in hadron cascade 
c	nosc = 1 : pythia type output only
c	       2 : also OSC1999A standard output event-by-event
c	       3 : no used 
c020708
c       itime: number of strings in current event
c       astr: number of strings in current event
c       gtime: number of gluon in current event
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
	read(11,*)ddt,dtt,bmin,bmax,nmax   ! 241108
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
	read(11,*)ttaup,taujp,iabsb,iabsm,nchan   ! 241108
	read(11,*)para13,para14,psno,para15,para16,ajpsi,vneum
	read(11,*)para1_1,para1_2,para2,para4   
	read(11,*)tdh,cptl,cptu,cptl2,cptu2,itnum   ! 241108
	read(11,*)mstu21,mstj1_1,mstj1_2,mstj2,mstj3
c210803
	read(11,*)(adj1(i),i=1,10)
	read(11,*)(adj1(i),i=11,20)
        read(11,*)(adj1(i),i=21,30)
        read(11,*)(adj1(i),i=31,40)
c210803
        read(11,*)kjp22,kjp23,kjp24   ! 020708 020511
	close(11)
c	tdh and itnum: time step and number of time steps used in subroutine 'flow_t'
c	cptl,cptu;cptl2,cptu2 : pt cut in 'flow_t' for particle 1;particle 2
c       nchan=0: inelastic (INEL)
c	nchan=1: Non Single Difractive (NSD) 
c	nchan=2: qqb --> gamma^*/Z^0, used to generate Drell-Yen
c	nchan=3: J/psi production
c	nchan=4: heavy-flavor production
c	nchan=5: direct photon
c	nchan=6: soft only
c	nchan=7: pythia 

c	print*,'after reading'
c	neve : # of events to be generate
c	nap (nzp) : # of nucleons (protons) in projectile nucleus
c	nat (nzt) : # of nucleons (protons) in target nucleus
c	t0 : average proper formation time at rest
c	ddt : time accuracy
c	dep : the accuracy in four momentum conservation
c	rou0 : normal nuclear density
c	rao : enlarge factor in the radius of simulated volume
c	bmin,bmax : minimum and maximum impact parameters, bmin=bmax means
c	 definite impact parameter, 2*nmax: the number of 
c	 intervals segmented in [bmin,bmax]	
c	kjp20: =1 constant cross sections 
c	       =0 energy dependent cross sections
c	kjp21: = 0 without hadron rescattering, 
c	       = 1 with hadron rescattering
c020708
c       kjp22: = 0 constant string tension with calculation for string effective tension
c       kjp22: = 1 variable effective string tension
c       kjp22: = 2 original case (constant string tension)
c020511 kjp23: = 1 npart calculated by geometric model
c020511	kjp23: = 2 npart calculated by Glauber model
c020511 kjp24: = 1 sharp sphere in Glauber model
c020511	kjp24: = 2 Woods-Saxon in Glauber model
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
	mstp(82)=adj1(8)
c	 =0: soft (two-string) only; =1: both of soft and hard
	parp(81)=adj1(9)
c	effective minimum transverse momentum  
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
c	mstj1_1: =0, no jet fragmentation at all used in parini.f
c	mstj1_2: =1, Lund string fragmentation used in sfm.f
c	no writing of header
	mstp(127)=0
	mstp(122)=0
c
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

 	ept=dsqrt(win*win+0.938*0.938)
        rots=dsqrt((ept+0.938)*(ept+0.938)-win*win)
        if(ifram.eq.1)rots=win
        degd=(1.-3.097/rots)**12
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
c020708
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
c020708
c061103
	pinel=0.
	pel=0.
c061103
c071103
	rinel=0.
	rel=0.
c071103
	sinel=0.
	sel=0.
	ssjp=0.
	nssjp=0
	do i1=1,600
	dinel(i1)=0.
	enddo
	do i1=1,20
	acoll(i1)=0.
	acollp(i1)=0.
	acollt(i1)=0.
	ett(i1)=0.
	fett(i1)=0.
        ettm(i1)=0.
        fettm(i1)=0.
	sbp(i1)=0.
	enddo
	do i1=1,10
	sndi(i1)=0.
	spdi(i1)=0.
	sedi(i1)=0.
c033101
	erat(i1)=0.
	eratf(i1)=0.
	erat2(i1)=0.
	eratf2(i1)=0.
	sbar(i1)=0.
	sabar(i1)=0.
	sbarf(i1)=0.
	sabarf(i1)=0.
	sbar2(i1)=0.
	sabar2(i1)=0.
	sbarf2(i1)=0.
	sabaf2(i1)=0.
c033101
	enddo
c081010	stime=0.
        stime_ini=0.d0   ! 081010
	stime_par=0.
	stime_had=0.
        stimei=0.d0   ! 081010
        stimep=0.d0   ! 081010
        stimeh=0.d0   ! 081010
	snspe=0.   ! 111899

csa*********************************************************************
	do i1=1,20
	do i2=1,4
	ssnum(i2,i1)=0.
	sv1(i2,i1)=0.
	sv2(i2,i1)=0.
	sv1o(i2,i1)=0.
	sv2o(i2,i1)=0.
	sv12(i2,i1)=0.
	sv22(i2,i1)=0.
	sv12o(i2,i1)=0.
	sv22o(i2,i1)=0.
	enddo
	enddo
	do i2=1,21
	do i1=1,20
        ssnum_pt(i2,i1)=0.   ! 280607
        ssnum_eta(i2,i1)=0.   ! 280607
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
c021207
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
c021207
	enddo
	enddo
csa**********************************************************************
c210803
        do i1=1,20
        savpt(i1)=0.
        sfvpt(i1)=0.
        nsavpt(i1)=0   ! 170705
        nsfvpt(i1)=0   ! 170705
        enddo
        savptn=0.
        sfvptn=0.
        savptp=0.
        sfvptp=0.
c210803

	volum=4.*3.1416/3.*2.**3
c261002	volume of sphere with radius of 2 fm in position space 
	qvolum=0.125   ! 261002   
c010220
        ncha=0
        nchaf=0
c010220
	if(nchan.eq.3)then
	nspsf=0
	nspsif=0
	nsps=0
	nspsi=0
	do i1=1,20
	pj(i1)=0.
	ej(i1)=0.
	pjp(i1)=0.
	ejp(i1)=0.
	fpj(i1)=0.
	fej(i1)=0.
	fpjp(i1)=0.
	fejp(i1)=0.
	enddo
	endif
c	nspsf: the number of primary J/psi in full phase space
c	nspsif: the number primary psi' in full phase space
c	nsps: the number primary J/psi in partial phase space
c	nspsi: the number primary psi' in partial phase space

c	ibp=0   ! July/20/98
c200601
	skpar=0.
	sknn=0.
	skpp=0.
	sknp=0.
c200601
	sthroq=0.
	sthrob=0.
        sthroc=0.
	do i1=1,4
	sthroe(i1)=0.
	enddo
c261002
        do i1=1,3
c       uds(i1)=0.
c       udsb(i1)=0.
        fuds(i1)=0.
        fudsb(i1)=0.
        enddo
	fg=0.   
	sxmax=0.
	symax=0.
	szmax=0.
	stmax=0.
	spxmax=0.
	spymax=0.
	spzmax=0.
	semax=0.
c261002
c181003
	sstop=0.
	nzstop=0   ! statistics of event w/o string 
c181003
	adj12=adj1(12)
	dpmax=adj1(27)   ! 281194
	adj140=adj1(40)   ! 290505
	open(5,file='sxp.out',status='unknown')   ! sa 26/05/99
	open(9,file='rms0.out',status='unknown')
	open(34,file='oscar.out',status='unknown')
c	give values to some important variables
	call sysini
c	write(9,*)'af. sysini csnn,kjp23,rou0,idw=',csnn,kjp23,rou0,idw !310805
c020511
        csnn1=csnn*10   ! csnn in fm^2 csnn1 in mb^2
        idw1=idw/50   ! *100
c        write(9,*)'csnn1,kjp24,idw1=',csnn1,kjp24,idw1
        call overlap(nap,nat,rnp,rnt,csnn1,kjp23,kjp24,rou0,idw1)
c020511
	adj1(28)=para10*dmax1(rnt,rnp)
	iii=0
c	write(9,100)frame
c	write(9,100)beam
c	write(9,100)target
	write(9,*)'win,nap,nzp,nat,nzt=',win,nap,nzp,nat,nzt
	write(9,*)'neve,nout,nosc=',neve,nout,nosc
	write(9,*)'dtt,bmin,bmax,nmax=',dtt,bmin,bmax,nmax
	write(9,*)'kjp21,ifram,para7,para10,kjp20,kjp22,kjp23,kjp24=',
     c   kjp21,ifram,para7,para10,kjp20,kjp22,kjp23,kjp24   ! 020511
	write(9,*)ispmax,isdmax,iflmax
	write(9,*)(ispkf(i),i=1,ispmax)
	write(9,*)(asd(i),i=1,isdmax)
c020708
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        write(9,*)'par1,par2,par3,par21,kjp22=',
     c   parj(1),parj(2),parj(3),parj(21),kjp22
        endif
c020708
	write(9,*)'parp21,parp22,ttaup,taujp=',parp21,parp22,ttaup,taujp
	write(9,*)'iabsb,iabsm,nchan=',iabsb,iabsm,nchan
        write(9,*)'para13 - 16 =',para13,para14,para15,para16
	write(9,*)'psno,rots,ajpsi,vneum=',psno,rots,ajpsi,vneum
        write(9,*)'para1_1,para1_2,para2,para4=',para1_1,para1_2,para2,
     c	 para4   
	write(9,*)'tdh,cptl,itnum=',tdh,cptl,itnum
	write(9,*)'cptu,cptl2,cptu2=',cptu,cptl2,cptu2
	write(9,*)'mstu21,mstj1_1,mstj1_2,mstj2,mstj3=',
     c	 mstu21,mstj1_1,mstj1_2,mstj2,mstj3   
c210803
	write(9,*)'adj1=',(adj1(i),i=1,10)
	write(9,*)'adj1=',(adj1(i),i=11,20)
        write(9,*)'adj1=',(adj1(i),i=21,30)
        write(9,*)'adj1=',(adj1(i),i=31,40)
c210803
c	write(9,*)'41-45=',parj(41),parj(42),parj(43),parj(44),
c     c	 parj(45)
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
c033101
	ratt=0.
	ratt2=0.
	rattf=0.
	ratt2f=0.
c033101
c041102
	qch=0.
	qch2=0.
	qchf=0.
	qchf2=0.
c041102
c010607
	tbar3=0
	tbar3f=0
	tabar3=0
	taba3f=0
	tbar1=0
	tbar1f=0
	tabar1=0
	taba1f=0
	tbare1=0
	tbre1f=0
	tbare3=0
	tbre3f=0
c010607
	ich=0	
	time=0.
c       write(9,*)'rnt,para10,adj1(28)=',rnt,para10,adj1(28)   ! sa
	jjj=1
c	for given b (impact parameter)
	if(dabs(bmin-bmax).lt.10d-4)then
	bp=bmin
        r4=rnp
        if(rnt.gt.rnp)r4=rnt
        rr4=bp/r4
        vneu=dexp(-rr4*rr4)
c	calculate the overlap region of two nuclei at given b by geometric 
c020511	 or Glauber method
c060605
	if(nap.ne.1 .and. nat.ne.1)then   ! 020511   
c020511
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
        pir=part1(ibpp)
        tir=part2(ibpp)
c       write(9,*)'bp,ibpp,part1,part2=',bp,ibpp,pir,tir
        endif
c020511
	if(nap.eq.1)pir=1.   
	if(nat.eq.1)tir=1.   
c060605
        vneump=pir   ! 111399
        vneumt=tir   ! 111399
	write(9,*)'bp,vneump,vneumt=',bp,vneump,vneumt   ! 190309
	goto 300
	endif
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
c060605
	if(nap.ne.1 .and. nat.ne.1)then   ! 020511
c020511
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
        pir=part1(ibpp)
        tir=part2(ibpp)
        endif
c020511
	if(nap.eq.1)pir=1.
	if(nat.eq.1)tir=1.
c060605
	acollp(i1)=pir
	acollt(i1)=tir
c	write(*,*)'bp,din,vneu,tb=',bp,din,acoll(i1),tb
	stbp=stbp+acollp(i1)
	stbt=stbt+acollt(i1)
	stb=stb+acoll(i1)
	enddo
	aneump=stbp/dfloat(i2)   ! 241110
	aneumt=stbt/dfloat(i2)   ! 241110
	vneum=stb/dfloat(i2)
	write(9,*)'neu=',(acoll(i1),i1=1,i2)   !!
	write(9,*)'(N_part)_p=',(acollp(i1),i1=1,i2)   ! 241110
        write(9,*)'(N_part)_t=',(acollt(i1),i1=1,i2)   ! 241110
	write(9,*)'(ave. participant)geo.=',aneump,aneumt,vneum   ! 241110
c	average b in [bmin,bmax]
        avb=2./3.*(bmin+bmax)
c	above equation is correct when bmin=0 only
	r4=rnp
	if(rnt.gt.rnp)r4=rnt
	rr4=avb/r4
	avneu=dexp(-rr4*rr4)
c	calculate the overlap region of two nuclei at given b 
c060605
	if(nap.ne.1 .and. nat.ne.1)then   ! 020511
c020511
        ibpp=int(avb/0.1+1.0)
        ibpp=min(ibpp,200)
        pir=part1(ibpp)
        tir=part2(ibpp)
        endif
c020511
	if(nap.eq.1)pir=1.
	if(nat.eq.1)tir=1.
c060605
	astbp=pir
	astbt=tir
c220110
        nrel=0
        nrea=0
        do i=1,9
        nreac(i)=0
        enddo
c220110
c	generate a event for nucleus-nucleou collision
	iran=adj1(26)
	if(iran.eq.0)goto 300
	do i1=1,iran
	thrr=pyr(1)
	enddo
	nncoll=0
300	iii=iii+1
c220110
        do i1=1,9
        nreaco(i1)=nreac(i1)
        enddo
c220110
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
        enddo
        bsin(i1)=0.   
        enddo
c250209
        do i1=1,21
        do i2=1,20
        snum_pt(i1,i2)=0.   ! 280607
        snum_eta(i1,i2)=0.   ! 280607
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
	siijk=0   ! 201203
c061103
	npel=0
        do i=1,600
        npinel(i)=0
        enddo
	noel=0
        do i=1,600
        noinel(i)=0
        enddo
c061103
	nspe=0   ! 111899
	nmine=0
	nplue=0
c033101
	nminef=0
	npluef=0
c033101

c010607
	tbar=0
	tbarf=0
	tabar=0
	tabarf=0 
c010507
        non6_c=123456   ! 141208
	sbp(jjj)=sbp(jjj)+1
c        open(22,file='main.out',status='unknown')
	if(dabs(bmin-bmax).lt.10d-4)goto 800
	bp=bpp(jjj) 
	vneump=acollp(jjj)   ! 241110
	vneumt=acollt(jjj)   ! 241110  
c241110 if(psno.eq.2)bp=dsqrt(pyr(1)*(bmax*bmax-bmin2)+bmin2)   ! 291207
	if(psno.eq.2)then
        bp=sqrt(pyr(1)*(bmax*bmax-bmin2)+bmin2)
c       calculate the overlap region of two nuclei at given bp 
        if(nap.ne.1 .and. nat.ne.1)then   ! 020511
c020511
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
        pir=part1(ibpp)
        tir=part2(ibpp)
        endif
c020511
        if(nap.eq.1)pir=1.
        if(nat.eq.1)tir=1.
        vneump=pir
        vneumt=tir
        write(9,*)'bp,vneump,vneumt=',bp,vneump,vneumt
        endif
c241110
c291207 psno: =1 systematic sampling method
c291207       =2 random sampling method
c       in case of psno=2 the outputs concerning jjj are nonsense
	vneu=acoll(jjj)
c	write(9,*)'in main vneu=',vneu   !!!
c	if(neve.le.40)write(9,*)'iii,jjj,bp,vneu=',iii,jjj,bp,vneu   ! sa
800     continue
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
c210803
        do i1=1,20
        avpt(i1)=0.
        navpt(i1)=0
        fvpt(i1)=0.
        nfvpt(i1)=0
        enddo
        navptn=0
        navptp=0
        avptn=0.
        avptp=0.
        nfvptn=0
        nfvptp=0
        fvptn=0.
        fvptp=0.
c210803
c	forbiden decay of particle, if mdcy(...)=0
c	mdcy(pycomp(111),1)=0
        mdcy(pycomp(310),1)=0   ! k0_S
        mdcy(pycomp(333),1)=0   ! phi
        mdcy(pycomp(3122),1)=0   ! Lambda
        mdcy(pycomp(-3122),1)=0
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
c	creat an event
c	parton initiation   
	ijk=0
        time_ini=0.d0   ! 081010 
c	write(9,*)'before call parini iii=',iii   ! sa
	if(itden.eq.1)call parini(time_ini,parp21,parp22,win,psno,ijk) ! 081010   
c060605	 for p+A and A+A
        if(ijk.eq.1)goto 300
	call pyedit(2)
c	no parton produced at all
	if(n.le.0)then
c	write(9,*)'neve,nncoll,n=',neve,nncoll,n   ! sa
	nncoll=nncoll+1
	if(nncoll.gt.neve)then
	write(9,*)'nncoll=',nncoll   ! sa
	stop 8888
	endif
	iii=iii-1
	goto 300
	endif
c        if(iii.eq.3)then 
c        write(mstu(11),*)'af parini'
c        call pylist(1)
c        write(mstu(11),*)'ppsa=',(ppsa(i1),i1=1,5)
c        call prt_sbh(nbh,cc)   
c        endif
	if(adj140.eq.1)then   ! 290505 271205 
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
        k(i,j)=0
        p(i,j)=0.
        v(i,j)=0.
        enddo
        enddo
5000    continue   ! 300407
c        if(iii.eq.3)call pylist(1)
        goto 888
        endif
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
	if(n.lt.2)goto 889   ! 151302
	if(itden.ne.1)goto 890   ! for e+e-,p+p,pbar_p, or p+pbar 080806 
c	partonic cascade
c	write(9,*)'before parcas iii=',iii   ! sa
c201203
c	'pyjets' to 'parlist'
        iprl=n
        do i=1,n
        idp(i)=k(i,2)
        rp(4,i)=0.
c       parton cascade process is assumed to start at time 0.
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
        time_par=0.d0   ! 081010
	iijk=0   ! 151203
	call parcas(time_par,nnn,iijk,win,nap,rnt,rnp,n00)   ! 120603 220110
c220110 nnn: nnn-th parton-parton interacion in a nucleus-nucleus collision
c2	write(9,*)'after parcas iii,iijk,time_par=',
c2     c	iii,iijk,time_par! sa
c120603
	if(iijk.eq.1)then
c	write(9,*)'iii,iijk=',iii,iijk
	goto 300   ! give up current event avoiding infinite collision loop
	endif
c120603
	if(iijk.eq.2)siijk=siijk+1   ! 201203
c201203
c       'parlist' to 'pyjets'
        n=iprl
	if(n.eq.0)goto 300   ! no parton at all, give up current event
        do i=1,n
c220110 k(i,1)=1
        k(i,2)=idp(i)
c220110 do j=3,5
c       k(i,j)=0
c220110 enddo
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
c       if(iii.eq.5)then 
c       write(22,*)'af parcas'
c       call pyedit(2)
c       call pylist(1)
c       call prt_sbh(nbh,cc)   
c       endif
c201203 
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
        endif
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
c       move '55' from 'lujets' to 'sgam'
        if(n55.gt.0)call remo_gam_par(55)
c250209
890	if(adj12.ne.0)goto 889   ! hadronization  
c	recover parton configuration in 'sbe' 
c220110
c       if(iii.eq.12.or.iii.eq.29.or.iii.eq.151.or.iii.eq.213)then
c       write(22,*)'be. recover n,nbe,n00,idi=',n,nbe,n00,idi
c       call prt_sbe(nbe,cc)
c       endif
c220110
c	loop over 'sbe'
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
c	write(22,*)'af recover jb,n=',jb,n
c	call prt_pyj(n,cc)

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

c220110do i1=1,n
c	do j1=1,5
c	pbe(i1,j1)=p(i1,j1)
c	enddo
c	enddo
c       'sbe' to 'pyjets'
c220110 call tran_sbe 
c220110
c       if(iii.eq.12.or.iii.eq.29.or.iii.eq.151.or.iii.eq.213)then
c       write(22,*)'af. recover n,nbh,nbe,noo,idi=',n,nbh,nbe,n00,idi   ! sa
cs      write(9,*)'af. recover n,nbh=',n,nbh   ! sa
c       call luedit(2)
c       call lulist(1)
cs      call prt_sbh(nbh,cc)
c       endif
c220110
889	continue
c	hadronization
c	write(9,*)'before hadronization iii=',iii   ! sa
c220110
	if(adj12.eq.0)then   ! 
c       write(9,*)'iii,4,4o,6,6o=',
c     c   iii,nreac(4),nreaco(4),nreac(6),nreaco(6)
        if(nreac(4).gt.nreaco(4) .or. nreac(6).gt.nreaco(6))then
c       for inela. processes 4 and 6
        call coales(iii,neve,nout,nap,nat,nzp,nzt)
        else
c       otherwise
        call sfm   ! string fragmentation 110604
        endif
        endif   !
c220110
	if(adj12.ne.0)then
c	write(9,*)'af parcas iprl=',iprl   ! sa
c	write(9,887)(idp(i2),i2=1,iprl)   ! sa
	call coales(iii,neve,nout,nap,nat,nzp,nzt)   ! coalescence
	endif
887	format(20(1x,i3))
c2	write(9,*)'after hadniz iii,nn,nbh=',iii,nn,nbh   ! sa
c161203
c010305	if(nout.eq.1 .or. iii.eq.1 .or. iii.eq.neve)then
c300407	if(nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
c     c	 .eq.neve)then 
c       if(iii.eq.5)then  
c	write(mstu(11),*)'af hadniz.'
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
c300407	endif
c       endif
c141208
        jb=1
111     do i1=jb,n   !
        kf=k(i1,2)
        kfab=iabs(kf)
        if(kfab.le.10.or.kfab.eq.21.or.kf.eq.2101.or.kf.eq.3101.or.
     c   kf.eq.3201.or.kf.eq.1103.or.kf.eq.2103.or.kf.eq.2203.or.kf
     c   .eq.3103.or.kf.eq.3203.or.kf.eq.3303)then   !!
        if(i1.eq.n)then
        n=n-1
        goto 222
        endif
        goto 333
        endif   !!
        enddo   !
        goto 222
333     continue
c       move particle list 'pyjets' one step downward since i1+1
        call updad_pyj(n,i1+1,1)
        n=n-1
        jb=i1
        goto 111
222     continue
c141208
c161203
	if(kjp21.eq.1)then   ! 1
        call filt
        do i=1,kfmax
        nup=numbs(i)
        enddo
        nbh1=n-nup
c       nup is the number of particles kept in 'pyjets'
c       nbh1 is the number of particles storing in 'sbh'
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
c	hadron cascade (rescattering)
c	if(nout.eq.1 .or. iii.eq.1 .or. iii.eq.neve)then
c	write(9,*)'be hadcas iii,n,nbh=',iii,n,nbh   ! sa
c	write(22,*)'be hadcas iii=',iii   !sa
c	call prt_sbh(nbh,cc)
c	call prt_sa1_h(nn)
c	endif  
        time_had=0.d0   ! 081010
	call hadcas(iii,neve,nout,time_had,ijkk)   ! 241103
c	write(9,*)'after hadcas iii,ijkk,time_had=',
c     c	 iii,ijkk,time_had   ! sa
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
c	if(iii.eq.7)then
c	write(mstu(11),*)'be luexec'
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
c	call prt_sbh(nbh,cc)
c	endif  
c250209
        n66=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=66
        n66=n66+1
        endif
        enddo
c       move "66" from 'lujets' to 'sgam'
        if(n66.gt.0)call remo_gam(66)
c250209
c       perform particle, declared unstable in the 'mdcy' array, decay
c130205	call pyexec
	rrp=1.16   ! 130205
	call decayh(rrp)   ! 130205
c181003
888	continue
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
	write(9,*)'ifram,rcsit,kfmax=',ifram,rcsit,kfmax
	write(9,*)(kfaco(i),i=1,kfmax)
	write(9,*)(disbe(i,i),i=1,kfmax)
	write(9,*)(disbe(1,i),i=1,8)
	write(9,*)'isinel='
	write(9,600)isinel
	endif
c050603
	call pyedit(2) 
	if(nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
     c	 .eq.neve)then   
	write(mstu(11),*)'event=',iii
	call pylist(1)
        write(mstu(11),*)'ppsa=',(ppsa(i1),i1=1,5)
	write(22,*)'throe_t=',throe_t
	write(22,*)'ithroq_t,ithrob_t,ich_t=',ithroq_t,ithrob_t,
     c	 dfloat(ich_t)/3.
	write(22,*)'throe_p=',throe_p
	write(22,*)'ithroq_p,ithrob_p,ich_p=',ithroq_p,ithrob_p,
     c	 dfloat(ich_p)/3.
	write(22,*)'throe=',throe
	write(22,*)'ithroq,ithrob,ithroc=',ithroq,ithrob,
     c	 dfloat(ithroc)/3.
c	call prt_sbh(nbh,cc)
        call prt_sgam(ngam)   ! 250209
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

	ibjf=0
	ibj=0
	ibpf=0
	ibp=0
	pbjf=0.
	pbj=0.
	pbpf=0.
	pbp=0.
	ssjp=sjp+ssjp
	nssjp=nsjp+nssjp

c????? following four statements are increaded on 26/05/99 by sa ******
csa	coox=coor(1)
csa	cooy=coor(2)
csa	cooz=coor(3)
c	write(5,*)coox,cooy,cooz
c**********************************************************************
c261002
	xmax=0.
	ymax=0.
	zmax=0.
	tmax=0.
	pxmax=0.
	pymax=0.
	pzmax=0.
	emax=0.
c261002

c	analyse an event
c130605	calculate trust axis for e+e-
        if(ipden.eq.2.and.itden.eq.2)then
	do j=1,3
	ptrus(j)=0.
	enddo
	do j=1,n
	ik=k(j,2)
	ikabs=iabs(ik)
	p1=p(j,1)
	if(p1.gt.0..and.(ikabs.eq.211.or.ikabs.eq.321.or.ikabs.eq.2212
     c	 .or.ikabs.eq.1.or.ikabs.eq.2.or.ikabs.eq.3))then
c	upper hemisphere
	p2=p(j,2)
        p3=p(j,3)
	ptrus(1)=ptrus(1)+p1
	ptrus(2)=ptrus(2)+p2
	ptrus(3)=ptrus(3)+p3
	endif
	enddo
	endif
c130605	
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
c021207 print p1,p2,p3 of charged particles
c        abplu6=dabs(plu6)
c        if(adj140.eq.4. .and. abplu6.gt.0.9)then
c        write(5,*)p1,p2,p3
c        endif
c        if(adj140.ne.4. .and. abplu6.gt.0.2)then
c        write(5,*)p1,p2,p3
c        endif
c021207
	ppt=pyp(j,10)
c060605
c	calculate the y_T for e+e-
	if(ipden.eq.2 .and. itden.eq.2)then
	ppm=p1*ptrus(1)+p2*ptrus(2)+p3*ptrus(3)
	pptru=ptrus(1)*ptrus(1)+ptrus(2)*ptrus(2)+ptrus(3)*ptrus(3)
	if(pptru.gt.1.e20)pptru=1.e20
	if(pptru.lt.1.e-20)pptru=1.e-20
	pptru=dsqrt(pptru)
	ppm=ppm/pptru
	epm=ppm*ppm+0.0196
c	0.0196=0.140*0.140
	if(epm.gt.1.e20)epm=1.e20
        if(epm.lt.1.e-20)epm=1.e-20
	epm=dsqrt(epm)
	denom=epm-ppm
	epnum=epm+ppm
	if(denom.lt.1.e-20)denom=1.e-20
	yy=0.5*dlog(epnum/denom)
	endif
c060605
c281104
	ppm=dsqrt(p1*p1+p2*p2+p3*p3)
	if(ppm.le.dpmax.and.p4.le.dpmax)then
	goto 3000
	else
	goto 400
	endif
c281104
c	write(9,*)'iii,n,j,ik,plu6=',iii,n,j,ik,plu6   ! 070802
c	analyse parton
c261002
3000	ikk=iabs(ik)   ! 281104
	if(ikk.eq.21)fg=fg+1   
        if(ikk.eq.1 .or. ikk.eq.2 .or. ikk.eq.3)then
        if(ik.gt.0)fuds(ik)=fuds(ik)+1   ! statistic u,d,s
        if(ik.lt.0)fudsb(ikk)=fudsb(ikk)+1   ! statistic anti-quark
        r1=v(j,1)
        r2=v(j,2)
        r3=v(j,3)
	r4=v(j,4)
	ar1=dabs(r1)
	ar2=dabs(r2)
	ar3=dabs(r3)
	ap1=dabs(p1)
	ap2=dabs(p2)
	ap3=dabs(p3)
	if(ar1.gt.xmax)xmax=ar1
	if(ar2.gt.ymax)ymax=ar2
	if(ar3.gt.zmax)zmax=ar3
	if(r4.gt.tmax)tmax=r4
	if(ap1.gt.pxmax)pxmax=ap1
	if(ap2.gt.pymax)pymax=ap2
	if(ap3.gt.pzmax)pzmax=ap3
	if(p4.gt.emax)emax=p4
        iadj=0
c	does the q is inside the cell with size 1 fm in position 
c	 space and 0.5 GeV/c in momentum sapce 
c	call adjins(r1,r2,r3,p1,p2,p3,iadj)
c        if(iadj.eq.1)then   ! yes
c        if(ik.gt.0)uds(ik)=uds(ik)+1
c        if(ik.lt.0)udsb(ikk)=udsb(ikk)+1
c	endif
	endif
c261002
c????? following four statements is increaded on 26/05/99 by sa ******
csa        x1=v(j,1)-coox
csa        y1=v(j,2)-cooy
csa        z1=v(j,3)-cooz
csa        t1=v(j,4)
c*********************************************************************

	c(1)=yy   
	if(ifram.eq.1)c(1)=eta   
	c(2)=ppt
c	.
c	.
c	.
	kkk=1
c033101
c	statistics of negative multiplicity
	if(adj140.eq.4. .and. plu6.lt.-0.9)then   ! for hadron
	nminf=nminf-plu6
	nminef=nminef-plu6
c210803
        nfvptn=nfvptn+1
        fvptn=fvptn+ppt
c210803
	do i=1,iflmax
	if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 700
	enddo
	nmin=nmin-plu6	
	nmine=nmine-plu6
c210803
        navptn=navptn+1
        avptn=avptn+ppt
c210803
700	endif
	if(adj140.ne.4. .and. plu6.lt.-0.2)then   ! for parton
        nminf=nminf+1.   ! -plu6 230206
        nminef=nminef-plu6
c210803
        nfvptn=nfvptn+1
        fvptn=fvptn+ppt
c210803
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 702
        enddo
        nmin=nmin+1.   ! -plu6 230206
        nmine=nmine-plu6
c210803
        navptn=navptn+1
        avptn=avptn+ppt
c210803
702     endif
c010220
c       statistics of charged particles multiplicity
c033101	statistics positive multiplicity
c070802
c       if(dabs(plu6).gt.0.9 .and. iabs(ik).ne.11)then
        if(adj140.eq.4. .and. plu6.gt.0.9)then   ! for hadron
c070802
        nchaf=nchaf+plu6 
	npluef=npluef+plu6
c210803
        nfvptp=nfvptp+1
        fvptp=fvptp+ppt
c210803
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 701
        enddo
        ncha=ncha+plu6 
	nplue=nplue+plu6
c210803
        navptp=navptp+1
        avptp=avptp+ppt
c210803
701     endif
	if(adj140.ne.4. .and. plu6.gt.0.2)then   ! for parton
c070802
        nchaf=nchaf+1.   ! +plu6 230206
        npluef=npluef+plu6
c210803
        nfvptp=nfvptp+1
        fvptp=fvptp+ppt
c210803
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 703
        enddo
        ncha=ncha+1.   ! +plu6 230206
        nplue=nplue+plu6
c210803
        navptp=navptp+1
        avptp=avptp+ppt
c210803
703     endif
c010220

c033101
c	calculate final state baryon (antibaryon) number fluctuation
c	 per baryon (antibaryon)
c	if(ik.eq.2212.or.ik.eq.2112.or.ik.eq.3122.or.
c     &   ik.eq.3212.or.ik.eq.3222)then
	if(ik.eq.2212)then
        tbarf=tbarf+1.
        do i1=1,iflmax
        if(c(i1).lt.afl(kkk,i1,1) .or. c(i1).gt.afl(kkk,i1,2))
     &   goto 402
        enddo
        tbar=tbar+1.
402     continue
        endif
c        if(ik.eq.-2212.or.ik.eq.-2112.or.ik.eq.-3122.or.
c     &   ik.eq.-3212.or.ik.eq.3112)then
	if(ik.eq.-2212)then
        tabarf=tabarf+1.
        do i1=1,iflmax
        if(c(i1).lt.afl(kkk,i1,1) .or. c(i1).gt.afl(kkk,i1,2))
     &   goto 403
        enddo
        tabar=tabar+1.
403     continue
        endif
c033101
c?????????????????????????????????????????????????????????????????????
c       statistics of negative and positive multiplicity for each event
c       if(plu6.lt.-0.2 .and. ik.ne.11)nmine=nmine+1   !070802
c       if(plu6.gt.0.2 .and. ik.ne.-11)nplue=nplue+1   !070802
c        if(nmine+nplue.gt.100)then
c        ich=ich+1
c        goto 1000
c        endif
c?????????????????????????????????????????????????????????????????????


c	statistics of transverse energy
c	if(iabs(ik).eq.2212 .or. iabs(ik).eq.211 .or. iabs(ik).eq.321)
c     c	 then
c	if(ik.eq.2112 .or. ik.eq.111 .or. ik.eq.310 .or. ik.eq.130 .or. 
c     c	 iabs(ik).eq.3122)then   ! neutral Et
	if(plu6.gt.-1.e-5 .and. plu6.lt.1.e-5)then   ! neutral Et
cc	if(ik.eq.22 .or. ik.eq.111)then   ! electro-magnetic Et (NA50)
	eii=pyp(j,4)
	set=pyp(j,13)
	sset=dsin(set)
	eiis=eii*sset
c	if(ik.eq.310)then
c	fettm(jjj)=fettm(jjj)+2.   ! counted due to neutral particles
c	fett(jjj)=fett(jjj)+2*eiis   ! counted due to neutral particles
c	fettm(jjj)=fettm(jjj)+2.*0.333   ! counted due to gamma
c	fett(jjj)=fett(jjj)+2*eiis*0.333   ! counted due to gamma
c	endif
c	assume K_L0 has the same transverse energy as K_S0
c	0.333: branch ratio of K_S0 -> gamma
	fettm(jjj)=fettm(jjj)+1.
	fett(jjj)=fett(jjj)+eiis
c	if(eta.ge.1.7 .and. eta.le.4.1)then   ! NA38 (PL,B262,362,1991)
	if(eta.ge.1.1 .and. eta.le.2.3)then   ! NA50 (L. Ramello)
c	if(ik.eq.310)then
c	ettm(jjj)=ettm(jjj)+2.   ! counted according to neutral particles
c	ett(jjj)=ett(jjj)+2*eiis   ! counted according to neutral particles
c       ettm(jjj)=ettm(jjj)+2.*0.333   ! counted according to gamma
c	ett(jjj)=ett(jjj)+2*eiis*0.333   ! counted according to gamma
c	endif
        ettm(jjj)=ettm(jjj)+1.
	ett(jjj)=ett(jjj)+eiis
	endif
810	endif

c	statistics of transverse momentum of J/psi
	if(nchan.eq.3 .and. ik.eq.443)then
	ibjf=ibjf+1
	eii=pyp(j,4)
	set=pyp(j,13)
	sset=dsin(set)
	eiis=eii*sset
	pluj9=pyp(j,9)
c        write(9,*)'iii,J/psi Pt=',iii,pluj9   ! sa
	fej(jjj)=fej(jjj)+eiis
	pbjf=pbjf+pluj9
	do i=1,iflmax
	if(c(i).lt.afl(6,i,1) .or. c(i).gt.afl(6,i,2))goto 820
	enddo
	ibj=ibj+1
	ej(jjj)=ej(jjj)+eiis
	pbj=pbj+pluj9
820	endif

c	statistics of transverse momentum of psi prime
	if(nchan.eq.3 .and. ik.eq.30443)then
	ibpf=ibpf+1
	eii=pyp(j,4)
	set=pyp(j,13)
	sset=dsin(set)
	eiis=eii*sset
	pluj9=pyp(j,9)
	fejp(jjj)=fejp(jjj)+eiis
	pbpf=pbpf+pluj9
	do i=1,iflmax
	if(c(i).lt.afl(7,i,1) .or. c(i).gt.afl(7,i,2))goto 830
	enddo
	ibp=ibp+1
	ejp(jjj)=ejp(jjj)+eiis
	pbp=pbp+pluj9
830	endif

c	statistics of y, pt, ect. distributions
	do 500 kk=1,ispmax
	kf=ispkf(kk)
	if(ik.ne.kf)goto 500

c????????????????????????????????????????????????????????????????????
        if(ik.eq.2212 .or. ik.eq.2112)then
c	eep=dsqrt(win*win+0.938*0.938)
c	if(ifram.eq.1)eep=0.5*win
c	eepd=eep-0.001
c	eepu=eep+0.001
c	if(ifram.eq.0 .and. ((p4.gt.eepd .and. p4.lt.eepu) .or.
c     c   p4.le.0.940))then   ! 111899
c	nspe=nspe+1   ! 111899
c	goto 500
c	endif   ! 111899
c	if(ifram.eq.1 .and. (p4.gt.eepd .and. p4.lt.eepu))then
	if(ppt.le.1.e-4)then   
	nspe=nspe+1   ! 111899
	goto 500
	endif   ! 111899 
        endif
c       exclude the projectile and the target spectator nucleons
c????????????????????????????????????????????????????????????????????
c210803
        nfvpt(kk)=nfvpt(kk)+1
        fvpt(kk)=fvpt(kk)+ppt
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 840
        enddo
        navpt(kk)=navpt(kk)+1
        avpt(kk)=avpt(kk)+ppt
c210803

840	call stati(yy,ppt,eta,p5,ik,kk,w,bn,an,bnf,anf)   ! 010600
	call flow_f(p1,p2,ppt,eta,p5,ik,kk)
c	write(9,*)'p1,p2,pt,p5,ik,kk=',p1,p2,ppt,p5,ik,kk   ! sa
c	write(9,*)'final state v1(pt)'
c	do i1=1,20
c	write(9,*)'i1=',i1
c	write(9,602)(v1_pt(i1,i2),i2=1,20)
c	enddo
	goto 901   ! sa 26/05/99 (originally  'goto 400')
500	continue
c?????? follow three statements are increaded on 26/05/99 by sa ******
901     if(nout.ne.1)then
csa        write(5,900)ik,x1,y1,z1,t1,p1,p2,p3,p4,p5
	endif
c********************************************************************
400	continue
c261002
	sxmax=sxmax+xmax
	symax=symax+ymax
	szmax=szmax+zmax
	stmax=stmax+tmax
	spxmax=spxmax+pxmax
	spymax=spymax+pymax
	spzmax=spzmax+pzmax
	semax=semax+emax
c261002

c010607
	tbar3=tbar3+tbar
	tbar3f=tbar3f+tbarf
	tabar3=tabar3+tabar
	taba3f=taba3f+tabarf
	tbar1=tbar1+tbar*tbar
	tbar1f=tbar1f+tbarf*tbarf
	tabar1=tabar1+tabar*tabar
	taba1f=taba1f+tabarf*tabarf
	tbare=tbar-tabar
	tbaref=tbarf-tabarf
	tbare3=tbare3+tbare
	tbre3f=tbre3f+tbaref
	tbare1=tbare1+tbare*tbare
	tbre1f=tbre1f+tbaref*tbaref
c010607

c033101
c	statistics of multiplicity distributions,
c	 spectator nucleons are excluded
	do kkk=1,ispmax
c	ik=ispkf(kkk)
c        if(iabs(ik).eq.2212.or.iabs(ik).eq.2112.or.iabs(ik).eq.3122.or.
c     &   iabs(ik).eq.3212.or.ik.eq.3222.or.ik.eq.3112)then
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
c033101
csa*******************************************************************
	do i1=1,2
	do i2=1,itnum
	ssnum(i1,i2)=ssnum(i1,i2)+snum(i1,i2)
	sv1(i1,i2)=sv1(i1,i2)+v1(i1,i2)
	sv2(i1,i2)=sv2(i1,i2)+v2(i1,i2)
	sv12(i1,i2)=sv12(i1,i2)+v12(i1,i2)
	sv22(i1,i2)=sv22(i1,i2)+v22(i1,i2)
	enddo
	enddo
	do i1=1,21
	do i2=1,20
        ssnum_pt(i1,i2)=ssnum_pt(i1,i2)+snum_pt(i1,i2)   ! 280607
        ssnum_eta(i1,i2)=ssnum_eta(i1,i2)+snum_eta(i1,i2)   ! 280607
	sv1_pt(i1,i2)=sv1_pt(i1,i2)+v1_pt(i1,i2)
	sv2_pt(i1,i2)=sv2_pt(i1,i2)+v2_pt(i1,i2)
	sv12_pt(i1,i2)=sv12_pt(i1,i2)+v12_pt(i1,i2)
	sv22_pt(i1,i2)=sv22_pt(i1,i2)+v22_pt(i1,i2)
	sv1_eta(i1,i2)=sv1_eta(i1,i2)+v1_eta(i1,i2)
	sv2_eta(i1,i2)=sv2_eta(i1,i2)+v2_eta(i1,i2)
	sv12_eta(i1,i2)=sv12_eta(i1,i2)+v12_eta(i1,i2)
	sv22_eta(i1,i2)=sv22_eta(i1,i2)+v22_eta(i1,i2)
c021207
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
c021207
	enddo
	enddo
csa*******************************************************************
	fibjf=dfloat(ibjf)
	if(fibjf.gt.1.e-15)fpj(jjj)=fpj(jjj)+pbjf/fibjf
	fibj=dfloat(ibj)
	if(fibj.gt.1.e-15)pj(jjj)=pj(jjj)+pbj/fibj
	fibpf=dfloat(ibpf)
	if(fibpf.gt.1.e-15)fpjp(jjj)=fpjp(jjj)+pbpf/fibpf
	fibp=dfloat(ibp)
	if(fibp.gt.1.e-15)pjp(jjj)=pjp(jjj)+pbp/fibp

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
        eevpav=eevpav+eevp   ! 271205
        eevhav=eevhav+eevh   ! 271205
c	write(9,*)'npel,npinel(592)=',npel,npinel(592)   ! 061103
	pel=pel+npel   ! 061103
	pinel=pinel+npinel(592)   
c220110 rel=rel+nrel   ! 071103
	sel=sel+noel
	do i1=1,600
c	pinel=pinel+npinel(i1)   ! 061103
	sinel=sinel+noinel(i1)
	dinel(i1)=dinel(i1)+noinel(i1)
	enddo
c071103
	do i1=1,9
	rinel=rinel+nreac(i1)
	enddo
c071103
	do i1=1,10
	sndi(i1)=sndi(i1)+dni(i1)
	spdi(i1)=spdi(i1)+dpi(i1)
	sedi(i1)=sedi(i1)+edi(i1)
c033101
	emini=emin(i1)
	eminfi=eminf(i1)
	eplui=eplu(i1)
	eplufi=epluf(i1)
	if(emini.ne.0.)then
	erat(i1)=erat(i1)+eplui/emini
	erat2(i1)=erat2(i1)+eplui*eplui/emini/emini
	endif
	if(eminfi.ne.0.)then
	eratf(i1)=eratf(i1)+eplufi/eminfi
	eratf2(i1)=eratf2(i1)+eplufi*eplufi/eminfi/eminfi
	endif
	bari=bar(i1)
	abari=abar(i1)
	barfi=barf(i1)
	abarfi=abarf(i1)
	sbar(i1)=sbar(i1)+bari
	sabar(i1)=sabar(i1)+abari
	sbarf(i1)=sbarf(i1)+barfi
	sabarf(i1)=sabarf(i1)+abarfi
	sbar2(i1)=sbar2(i1)+bari*bari
	sabar2(i1)=sabar2(i1)+abari*abari
	sbarf2(i1)=sbarf2(i1)+barfi*barfi
	sabaf2(i1)=sabaf2(i1)+abarfi*abarfi
c033101
	enddo
c033101
	if(nmine.ne.0)then
	ratte=nplue/nmine
	ratt=ratt+ratte
	ratt2=ratt2+ratte*ratte
	endif
	if(nminef.ne.0)then
	ratte=npluef/nminef
	rattf=rattf+ratte
	ratt2f=ratt2f+ratte*ratte
	endif
c033101
c041102
	qche=nplue-nmine
	qch=qch+qche
	qch2=qch2+qche*qche
	qche=npluef-nminef
	qchf=qchf+qche
	qchf2=qchf2+qche*qche
c041102
c200601
	skpar=skpar+kpar
	sknn=sknn+knn
	skpp=skpp+kpp
	sknp=sknp+knp
c200601
	sthroq=sthroq+ithroq+ithroq_p+ithroq_t
	sthrob=sthrob+ithrob+ithrob_p+ithrob_t
        sthroc=sthroc+ithroc+ich_p+ich_t
	do i1=1,4
	sthroe(i1)=sthroe(i1)+throe(i1)+throe_p(i1)+throe_t(i1)
	enddo
c081010	stime=stime+time_par+time_had
        stime_ini=stime_ini+time_ini   ! 081010
	stime_par=stime_par+time_par
	stime_had=stime_had+time_had
        if(time_ini.gt.0.d0)stimei=stimei+1.d0   ! 081010
        if(time_par.gt.0.d0)stimep=stimep+1.d0   ! 081010
        if(time_had.gt.0.d0)stimeh=stimeh+1.d0   ! 081010
	snspe=snspe+nspe   ! 111899
c210803
        do i1=1,20
        davpt=dfloat(navpt(i1))
        if(davpt.gt.0.)then   ! 170705
	savpt(i1)=savpt(i1)+avpt(i1)/davpt
	nsavpt(i1)=nsavpt(i1)+1 ! statistics of event with particle i1 (partial)
	endif   ! 170705
        davpt=dfloat(nfvpt(i1))
        if(davpt.gt.0.)then   ! 170705
	sfvpt(i1)=sfvpt(i1)+fvpt(i1)/davpt
        nsfvpt(i1)=nsfvpt(i1)+1 ! statistics of event with particle i1 (full)
	endif   ! 170705
        enddo
        davpt=dfloat(navptn)
        if(davpt.gt.0.)savptn=savptn+avptn/davpt
        davpt=dfloat(navptp)
        if(davpt.gt.0.)savptp=savptp+avptp/davpt
        davpt=dfloat(nfvptn)
        if(davpt.gt.0.)sfvptn=sfvptn+fvptn/davpt
        davpt=dfloat(nfvptp)
        if(davpt.gt.0.)sfvptp=sfvptp+fvptp/davpt
c210803

	if(nchan.eq.3)then
	call stapsi(nspsf,nspsif,nsps,nspsi)
        if(nssjp.ge.1 .and. ssjp.ge.1.e-15)then
        edeg=ssjp/nssjp
        degn=(1.-3.097/edeg)**12
        degr=degn/degd
	endif
c       degr : suppression factor caused by degradation
        endif

c020708
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
c020708
	open(8,file='nout.out',status='unknown')
	write(8,*)'iii=',iii
	close(8)
c	internal and final printing and controled return
	if(mod(iii,nout).eq.0 .or. iii.eq.neve)then
	open(10,file='rms.out',status='unknown')
c300404
c       pt1.dat, for example, stores pt dis. (two column) of 1-th particle
	open(41,file='pt1.dat',status='unknown')
	open(42,file='pt2.dat',status='unknown')
	open(43,file='pt3.dat',status='unknown')
	open(44,file='pt4.dat',status='unknown')
	open(45,file='pt5.dat',status='unknown')
	open(46,file='pt6.dat',status='unknown')
	open(47,file='pt7.dat',status='unknown')
	open(48,file='pt8.dat',status='unknown')
	open(49,file='pt9.dat',status='unknown')
	open(50,file='pt10.dat',status='unknown')
	open(51,file='pt11.dat',status='unknown')
	open(52,file='pt12.dat',status='unknown')
	open(53,file='pt13.dat',status='unknown')
	open(54,file='pt14.dat',status='unknown')
	open(55,file='pt15.dat',status='unknown')
	open(56,file='pt16.dat',status='unknown')
	open(57,file='pt17.dat',status='unknown')
	open(58,file='pt18.dat',status='unknown')
	open(59,file='pt19.dat',status='unknown')
	open(60,file='pt20.dat',status='unknown')
	open(61,file='pt1f.dat',status='unknown')
	open(62,file='pt2f.dat',status='unknown')
	open(63,file='pt3f.dat',status='unknown')
	open(64,file='pt4f.dat',status='unknown')
	open(65,file='pt5f.dat',status='unknown')
	open(66,file='pt6f.dat',status='unknown')
	open(67,file='pt7f.dat',status='unknown')
	open(68,file='pt8f.dat',status='unknown')
	open(69,file='pt9f.dat',status='unknown')
	open(70,file='pt10f.dat',status='unknown')
	open(71,file='pt11f.dat',status='unknown')
	open(72,file='pt12f.dat',status='unknown')
	open(73,file='pt13f.dat',status='unknown')
	open(74,file='pt14f.dat',status='unknown')
	open(75,file='pt15f.dat',status='unknown')
	open(76,file='pt16f.dat',status='unknown')
	open(77,file='pt17f.dat',status='unknown')
	open(78,file='pt18f.dat',status='unknown')
	open(79,file='pt19f.dat',status='unknown')
	open(80,file='pt20f.dat',status='unknown')
c300404
	flaa=dfloat(iii-ich)   ! July/20/98
c181003
	flaaa=dfloat(iii-ich-nzstop)
	if(flaaa.gt.1.d-20)dsstop=sstop/flaaa
c181003
	if(flaa.le.1.e-20)goto 1200
c261002
	dxmax=sxmax/flaa
	dymax=symax/flaa
	dzmax=szmax/flaa
	dtmax=stmax/flaa
	dpxmax=spxmax/flaa
	dpymax=spymax/flaa
	dpzmax=spzmax/flaa
	demax=semax/flaa
        do i1=1,3
        fduds(i1)=fuds(i1)/flaa
        fdudsb(i1)=fudsb(i1)/flaa
c	uuds(i1)=uds(i1)/flaa
c	uudsb(i1)=udsb(i1)/flaa
c       duds(i1)=uds(i1)/flaa/qvolum
c       dudsb(i1)=udsb(i1)/flaa/qvolum
        enddo
	fdg=fg/flaa   
c261002
	dnmino=nmin/flaa
	dnminfo=nminf/flaa
c010200
        dncha=ncha/flaa
        dnchaf=nchaf/flaa
c010220

c033101
	ratpm=ratt/flaa
	ratpmf=rattf/flaa
	ratpm2=ratt2/flaa
	rapmf2=ratt2f/flaa
c041102
	aqch=qch/flaa
	aqch2=qch2/flaa
	aqchf=qchf/flaa
	aqchf2=qchf2/flaa
c041102
c010607	tbaro=tbar/flaa
c010607	tbarfo=tbarf/flaa
c010607	tabaro=tabar/flaa
c010607	tabafo=tabarf/flaa
c010607	tbar2=tbar*tbar/flaa
c010607	tbarf2=tbarf*tbarf/flaa
c010607	tabar2=tabar*tabar/flaa
c010607	tabaf2=tabarf*tabarf/flaa
c010607
	tbaro=tbar3/flaa
	tbarfo=tbar3f/flaa
	tabaro=tabar3/flaa
	tabafo=taba3f/flaa
	tbar2=tbar1/flaa
	tbarf2=tbar1f/flaa
	tabar2=tabar1/flaa
	tabaf2=taba1f/flaa
c010607
c010607	tbare=tbar-tabar
c010607	tbaref=tbarf-tabarf
c010607	tbard=tbare/flaa
c010607	tbardf=tbaref/flaa
c010607 tbard2=tbare*tbare/flaa
c010607 tbadf2=tbaref*tbaref/flaa

c010607
	tbard=tbare3/flaa
	tbardf=tbre3f/flaa
	tbard2=tbare1/flaa
	tbadf2=tbre1f/flaa
c010607
c033101
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
c020708
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
c020708
        eevpao=eevpav/flaa   ! 271205
        eevhao=eevhav/flaa   ! 271205

csa*******************************************************************
	do i1=1,2
	do i2=1,itnum
	sss=ssnum(i1,i2)
	if(sss.gt.1.e-5)then
	sss1=sv1(i1,i2)/sss
	sss2=sv2(i1,i2)/sss
	sv1o(i1,i2)=sss1
        sv2o(i1,i2)=sss2
c121207 sv12o(i1,i2)=dsqrt((sv12(i1,i2)/sss-sss1*sss1)/(sss-1))
c121207 sv22o(i1,i2)=dsqrt((sv22(i1,i2)/sss-sss2*sss2)/(sss-1))
        sv12o(i1,i2)=dsqrt(sv12(i1,i2)/sss-sss1*sss1)   ! 121207
        sv22o(i1,i2)=dsqrt(sv22(i1,i2)/sss-sss2*sss2)   ! 121207
	endif
	enddo
	enddo   
c220607
        do i1=1,21
        do i2=1,20
        sss=ssnum_pt(i1,i2)   ! 280607
        if(sss.gt.1.e-5)then   ! 280607
        sv1op=sv1_pt(i1,i2)/sss   ! direct, average
        sv1o_pt(i1,i2)=sv1op
        sv12op=sv12_pt(i1,i2)/sss   ! direct, squared average
c121207 if(sss-1.gt.1.e-5)then   ! 280607
c121207 sv12op=(sv12op-sv1op*sv1op)/(sss-1.)
        sv12op=sv12op-sv1op*sv1op   ! 121207
        if(sv12op.le.1.e-20)sv12op=1.e-20
        sv12o_pt(i1,i2)=dsqrt(sv12op)   !  direct, error
c121207 endif   ! 280607
        sv2op=sv2_pt(i1,i2)/sss   ! elliptic, average
        sv2o_pt(i1,i2)=sv2op
        sv22op=sv22_pt(i1,i2)/sss   ! elliptic, squared average
c121207 if(sss-1.gt.1.e-5)then   ! 280607
c121207 sv22op=(sv22op-sv2op*sv2op)/(sss-1.)
        sv22op=sv22op-sv2op*sv2op   ! 121207
        if(sv22op.le.1.e-20)sv22op=1.e-20
        sv22o_pt(i1,i2)=dsqrt(sv22op)   ! elliptic, error
c121207 endif   ! 280607
        endif   ! 280607
        sss=ssnum_eta(i1,i2)   ! 280607
        if(sss.gt.1.e-5)then   ! 280607
        sv1oe=sv1_eta(i1,i2)/sss   ! direct, average
        sv1o_eta(i1,i2)=sv1oe   
        sv12oe=sv12_eta(i1,i2)/sss   ! direct, squared average
c121207 if(sss-1.gt.1.e-5)then   ! 280607
c121207 sv12oe=(sv12oe-sv1oe*sv1oe)/(sss-1.)
        sv12oe=sv12oe-sv1oe*sv1oe   ! 121207
        if(sv12oe.le.1.e-20)sv12oe=1.e-20
        sv12o_eta(i1,i2)=dsqrt(sv12oe)   !  direct, error
c121207 endif   ! 280607   
        sv2oe=sv2_eta(i1,i2)/sss   ! elliptic, average
        sv2o_eta(i1,i2)=sv2oe   
        sv22oe=sv22_eta(i1,i2)/sss   ! elliptic, squared average
c121207 if(sss-1.gt.1.e-5)then   ! 280607
c121207 sv22oe=(sv22oe-sv2oe*sv2oe)/(sss-1.)
        sv22oe=sv22oe-sv2oe*sv2oe   ! 121207
        if(sv22oe.le.1.e-20)sv22oe=1.e-20
        sv22o_eta(i1,i2)=dsqrt(sv22oe)   ! elliptic, error
c121207 endif   ! 280607
        endif   ! 280607
        enddo
        enddo
c220607
c061103
	peli=pel/flaa
	pineli=pinel/flaa
c061103
c071103
	reli=nrel/flaa   ! original =rel/flaa 220110
	rineli=rinel/flaa
c071103
	seli=sel/flaa
	sineli=sinel/flaa
	do i1=1,600
	dineli(i1)=dinel(i1)/flaa
	enddo
	do i1=1,10
	sndio(i1)=sndi(i1)/flaa/volum
	spdio(i1)=spdi(i1)/flaa/volum
	sedio(i1)=sedi(i1)/flaa/volum
c033101
	orat(i1)=erat(i1)/flaa
	oratf(i1)=eratf(i1)/flaa
	orat2(i1)=erat2(i1)/flaa
	oratf2(i1)=eratf2(i1)/flaa
	obar(i1)=sbar(i1)/flaa
	oabar(i1)=sabar(i1)/flaa
	obarf(i1)=sbarf(i1)/flaa
	oabarf(i1)=sabarf(i1)/flaa
	obar2(i1)=sbar2(i1)/flaa
	oabar2(i1)=sabar2(i1)/flaa
	obarf2(i1)=sbarf2(i1)/flaa
	oabaf2(i1)=sabaf2(i1)/flaa
c033101
	enddo
c200601
	skparo=skpar/flaa
	sknno=sknn/flaa
	skppo=skpp/flaa
	sknpo=sknp/flaa
c200601
c021207
        do i1=1,21
        do i2=1,20
        sv1op=ssv1_pt(i1,i2)/flaa
        ssv1o_pt(i1,i2)=sv1op   ! direct, average
        sv12op=ssv12_pt(i1,i2)/flaa   ! direct, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv12op=(sv12op-sv1op*sv1op)/(flaa-1.)
        sv12op=sv12op-sv1op*sv1op   ! 121207
        if(sv12op.le.1.e-20)sv12op=1.e-20
        ssv12o_pt(i1,i2)=dsqrt(sv12op)   !  direct, error
c121207 endif
        sv2op=ssv2_pt(i1,i2)/flaa
        ssv2o_pt(i1,i2)=sv2op   ! elliptic, average
        sv22op=ssv22_pt(i1,i2)/flaa   ! elliptic, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv22op=(sv22op-sv2op*sv2op)/(flaa-1.)
        sv22op=sv22op-sv2op*sv2op   ! 121207
        if(sv22op.le.1.e-20)sv22op=1.e-20
        ssv22o_pt(i1,i2)=dsqrt(sv22op)   ! elliptic, error
c121207 endif
        sv1op=ssv1_eta(i1,i2)/flaa
        ssv1o_eta(i1,i2)=sv1op   ! direct, average
        sv12op=ssv12_eta(i1,i2)/flaa   ! direct, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv12op=(sv12op-sv1op*sv1op)/(flaa-1.)
        sv12op=sv12op-sv1op*sv1op   ! 121207
        if(sv12op.le.1.e-20)sv12op=1.e-20
        ssv12o_eta(i1,i2)=dsqrt(sv12op)   !  direct, error
c121207 endif
        sv2op=ssv2_eta(i1,i2)/flaa
        ssv2o_eta(i1,i2)=sv2op   ! elliptic, average
        sv22op=ssv22_eta(i1,i2)/flaa   ! elliptic, squared average
c121207 if(flaa-1.gt.0.)then
c121207 sv22op=(sv22op-sv2op*sv2op)/(flaa-1.)
        sv22op=sv22op-sv2op*sv2op   ! 121207
        if(sv22op.le.1.e-20)sv22op=1.e-20
        ssv22o_eta(i1,i2)=dsqrt(sv22op)   ! elliptic, error
c121207 endif
        enddo
        enddo
c021207
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
c081010	timeo=stime/flaa
c081010	timeo_par=stime_par/flaa
        if(stimei.gt.0.d0)then
        timeo_ini=stime_ini/stimei
        else
        timeo_ini=0.d0
        endif
        if(stimep.gt.0.d0)then
        timeo_par=stime_par/stimep
        else
        timeo_par=0.d0
        endif
c081010	timeo_had=stime_had/flaa
        if(stimeh.gt.0.d0)then
        timeo_had=stime_had/stimeh
        else
        timeo_had=0.d0
        endif
        timeo=timeo_ini+timeo_par+timeo_had
c081010
	snspeo=snspe/flaa   ! 111899
c210803
        do i1=1,20
	avpti1=dfloat(nsavpt(i1)) !   170705
        if(avpti1.gt.0.)avpto(i1)=savpt(i1)/avpti1 !   flaa 170705
	avpti1=dfloat(nsfvpt(i1)) !   170705
        if(avpti1.gt.0.)fvpto(i1)=sfvpt(i1)/avpti1 !   flaa 170705
        enddo
        avptno=savptn/flaa
        avptpo=savptp/flaa
        fvptno=sfvptn/flaa
        fvptpo=sfvptp/flaa
c210803
	do i1=1,10
	sbpi=sbp(i1)
	if(sbpi.gt.1.e-15)then
	etti1=ett(i1)/sbpi
	etti(i1)=etti1
	fetti1=fett(i1)/sbpi
	fetti(i1)=fetti1
        if(etti1.ge.1.e-15)ettmi(i1)=ettm(i1)/sbpi/etti1
        if(fetti1.ge.1.e-15)fettmi(i1)=fettm(i1)/sbpi/fetti1
	endif
	enddo

	if(nchan.eq.3)then
	spsfo=nspsf/flaa
	spsifo=nspsif/flaa
	spso=nsps/flaa
	spsio=nspsi/flaa
	do i1=1,10
	sbpi=sbp(i1)
	if(sbpi.gt.1.e-15)then
	eji(i1)=ej(i1)/sbpi
	pji(i1)=pj(i1)/sbpi
	ejpi(i1)=ejp(i1)/sbpi
	pjpi(i1)=pjp(i1)/sbpi
	feji(i1)=fej(i1)/sbpi
	fpji(i1)=fpj(i1)/sbpi
	fejpi(i1)=fejp(i1)/sbpi
	fpjpi(i1)=fpjp(i1)/sbpi
	endif
	enddo
	endif
1200	continue
c061103
	write(10,*)'msel,mint(48),msub =',msel,mint(48)
        write(10,*)(msub(i),i=1,25) 
	write(10,*)(msub(i),i=26,50)
        write(10,*)(msub(i),i=51,75) 
	write(10,*)(msub(i),i=76,100)
        write(10,*)(msub(i),i=101,125) 
	write(10,*)(msub(i),i=126,150)
        write(10,*)(msub(i),i=151,175) 
	write(10,*)(msub(i),i=176,200)
	write(10,*)'parp81,bp,mstp82=',parp(81),bp,mstp(82)    ! 291207  
	write(10,*)'nn colli. # and blocked # in parton initialization=',
     c	 pineli,peli,peli+pineli
c071103
	write(10,*)'# of successful and blocked collision in parton 
     c	 cascade=',rineli,reli,reli+rineli
c071103
c220110
        write(10,*)'average collision # in parton cascade=',srea
        write(10,*)'total # of scaterring processes in parton cascade'
        write(10,*)(snreac(i1),i1=1,9)
c220110
	write(10,*)'el. and inela. colli. # and sum in hadron cascade=',
     c	 seli,sineli,seli+sineli
	write(10,*)'time_ini,par,had, sum=',timeo_ini,timeo_par,
     c   timeo_had,timeo   ! 081010
c061103
	write(10,*)'statistics case of no colli. in parcas =',siijk!201203  
c200601
	write(10,*)'(Npart)mini-jet,Nnn,Npp=',skparo,sknno,skppo
	write(10,*)'Nnp,Ntot=',sknpo,sknno+skppo+sknpo
c200601
	write(10,*)'(Npart)spec. and z=',nap+nat-snspeo,dsstop 
c020708
        if(kjp22.eq.0 .or. kjp22.eq.1)then
        write(10,*)'kjp22=0, par1,par2,par3,par21=',
     c   parj1,parj2,parj3,parj21
        write(10,*)'keff2,par2,par21,par1,par3=',
     c   (skapao(i1),i1=1,5)
        write(10,*)'average # of gluons in a string,average # of'
        write(10,*)'strings in an event and average # of events' 
        write(10,*)'with string',sgtimeo,sastro,mkapa
        endif
c020708
        write(10,*)'multiplicity of negative particles=',dnmino
	write(10,*)'multiplicity of negative particles=',dnminfo
c010220
c033101	write(10,*)'multiplicity of charged particles,partial=',dncha
c033101	write(10,*)'multiplicity of charged particles,full=',dnchaf
c033101
	write(10,*)'multiplicity of positive particles,partial=',dncha
	write(10,*)'multiplicity of positive particles,full=',dnchaf
c033101
	write(10,*)'throw away ithroq,ithrob,ithroc=',
     c	 wthroq,wthrob,wthroc/3.
	write(10,*)'throe=',wthroe
c210803
c        write(10,*)'neg. & posi. cha. par. <pt> partial',avptno,avptpo
c        write(10,*)'neg. & posi. cha. par. <pt> full',fvptno,fvptpo
c        write(10,*)'<pt> of produced particle, partial'
c        write(10,*)(avpto(i1),i1=1,5)
c        write(10,*)(avpto(i1),i1=6,10)
c        write(10,*)(avpto(i1),i1=11,15)
c        write(10,*)(avpto(i1),i1=16,20)
c        write(10,*)'<pt> of produced particle, full'
c        write(10,*)(fvpto(i1),i1=1,5)
c        write(10,*)(fvpto(i1),i1=6,10)
c        write(10,*)(fvpto(i1),i1=11,15)
c        write(10,*)(fvpto(i1),i1=16,20)
c210803
c261002
c	write(10,*)'average largest x,y,z,t=',dxmax,dymax,dzmax,dtmax
c	write(10,*)'average largest px,py,pz,e=',dpxmax,dpymax,dpzmax,
c     c	 demax
c        write(10,*)'average total d,u,s=',(fduds(i1),i1=1,3)
c        write(10,*)'average toatl d,u,s bar=',(fdudsb(i1),i1=1,3)
c	write(10,*)'average toatl g=',fdg   
c	write(10,*)'volum=',volum
c        write(10,*)'average number of d,u,s=',(uuds(i1),i1=1,3)
c        write(10,*)'average number of d,u,s bar=',(uudsb(i1),i1=1,3)
c        write(10,*)'average d,u,s density=',(duds(i1),i1=1,3)
c        write(10,*)'average d,u,s bar density=',(dudsb(i1),i1=1,3)
c261002
c033101
c	write(10,*)'final (+) to (-) charge ratio,partial',ratpm,ratpm2
c	write(10,*)'final (+) to (-) charge ratio,full',ratpmf,rapmf2
c	write(10,*)'final (+) to (-) charge ratio fluctuation,partial'
c	write(10,*)(ratpm2-ratpm*ratpm)*(dncha+dnmino)
c	write(10,*)'final (+) to (-) charge ratio fluctuation,full'
c	write(10,*)(rapmf2-ratpmf*ratpmf)*(dnchaf+dnminfo)
c041102
c	write(10,*)'net charge,partial',aqch,aqch2,aqch2-aqch*aqch
c	write(10,*)'net charge,full',aqchf,aqchf2,aqchf2-aqchf*aqchf
c041102
c	write(10,*)'final baryon number fluctuation per baryon,partial'
c	write(10,*)'tbaro,tbar2=',tbaro,tbar2
c	if(tbaro.ne.0.)write(10,*)(tbar2-tbaro*tbaro)/tbaro
c	write(10,*)'final baryon number fluctuation per baryon,full'
c	write(10,*)'tbarfo,tbarf2=',tbarfo,tbarf2
c	if(tbarfo.ne.0.)write(10,*)(tbarf2-tbarfo*tbarfo)/tbarfo
c	write(10,*)'final antibaryon fluctuation per antibaryon,partial'
c	write(10,*)'tabaro,tabar2=',tabaro,tabar2
c	if(tabaro.ne.0.)write(10,*)(tabar2-tabaro*tabaro)/tabaro
c	write(10,*)'final antibaryon fluctuation per antibaryon,full'
c	write(10,*)'tabafo,tabaf2=',tabafo,tabaf2
c	if(tabafo.ne.0.)write(10,*)(tabaf2-tabafo*tabafo)/tabafo
c	write(10,*)'final net baryon fluctuation per ...,partial'
c	write(10,*)'tbard,tbard2=',tbard,tbard2
c	if(tbaro+tabaro.ne.0.)write(10,*)(tbard2-tbard*tbard)/
c    &	 (tbaro+tabaro)
c	write(10,*)'final net baryon fluctuation per ...,full'
c	write(10,*)'tbardf,tbadf2=',tbardf,tbadf2	
c	if(tbarfo+tabafo.ne.0.)write(10,*)(tbadf2-tbardf*tbardf)/
c     &	 (tbarfo+tabafo)

c033101
c        write(10,*)'avb,avneu,astbp,astbt=',avb,avneu,astbp,astbt
c	write(10,*)'ave. participant=',aneump,aneumt,vneum   ! 241110
c        write(10,*)'ich=',ich   ! July/20/98
        write(10,*)'particle multiplicity=',(sbof(ll),ll=1,ispmax)
	write(10,*)'particle multiplicity=',(sbo(ll),ll=1,ispmax)
c	write(10,*)'nucleon number density in center'
c	write(10,*)(sndio(i1),i1=1,10)
c	write(10,*)'particle number density in center'
c	write(10,*)(spdio(i1),i1=1,10)
c	write(10,*)'energy density in center'
c	write(10,*)(sedio(i1),i1=1,10)
c033101
c	write(10,*)'baryon number fluctuation per baryon,partial'
c	do i1=1,10
c	bb=obar(i1)
c	bb2=obar2(i1)
c	if(bb.ne.0.)then
c	write(10,*)'<o2>,<o>2,vari.=',bb,bb2,(bb2-bb*bb)/bb
c	else
c	write(10,*)'<o2>,<o>2=',bb,bb2
c	endif
c	enddo
c	write(10,*)'baryon number fluctuation per baryon,full'
c	do i1=1,10
c	bb=obarf(i1)
c	bb2=obarf2(i1)
c	if(bb.ne.0.)then
c	write(10,*)'<o2>,<o>2,vari.=',bb,bb2,(bb2-bb*bb)/bb
c	else
c	write(10,*)'<o2>,<o>2=',bb,bb2
c	endif
c	enddo
c      write(10,*)'antibaryon number fluctuation per antibaryon,partial'
c	do i1=1,10
c	bb=oabar(i1)
c	bb2=oabar2(i1)
c	if(bb.ne.0.)then
c	write(10,*)'<o2>,<o>2,vari.=',bb,bb2,(bb2-bb*bb)/bb
c	else
c	write(10,*)'<o2>,<o>2=',bb,bb2
c	endif
c	enddo
c	write(10,*)'anaibaryon number fluctuation per antibaryon,full'
c	do i1=1,10
c	bb=oabarf(i1)
c	bb2=oabaf2(i1)
c	if(bb.ne.0.)then
c	write(10,*)'<o2>,<o>2,vari.=',bb,bb2,(bb2-bb*bb)/bb
c	else
c	write(10,*)'<o2>,<o>2=',bb,bb2
c	endif
c	enddo
c	write(10,*)'(+) to (-) charge ratio fluctuation,partial'
c	do i1=1,10
c	bb=orat(i1)
c	bb2=orat2(i1)
c	if(bb.ne.0.)then
c	write(10,*)'<o2>,<o>2,vari.=',bb,bb2,(bb2-bb*bb)/bb
c	else
c	write(10,*)'<o2>,<o>2=',bb,bb2
c	endif
c	enddo
c	write(10,*)'(+) to (-) charge ratio fluctuation,full'
c	do i1=1,10
c	bb=oratf(i1)
c	bb2=oratf2(i1)
c	if(bb.ne.0.)then
c	write(10,*)'<o2>,<o>2,vari.=',bb,bb2,(bb2-bb*bb)/bb
c	else
c	write(10,*)'<o2>,<o>2=',bb,bb2
c	endif
c	enddo
c033101
c	write(10,*)'total transverse energy='
c	write(10,*)(fetti(i1),i1=1,10)
c	write(10,*)(etti(i1),i1=1,10)
c       write(10,*)'total transverse energy distribution='
c        write(10,*)(fettmi(i1),i1=1,10)
c        write(10,*)(ettmi(i1),i1=1,10)
c        if(nchan.eq.3)then
c        write(10,*)'iii,njp,edeg,degr=',iii,nssjp,edeg,degr
c        write(10,*)'J/psi and psi=',spsfo,spsifo,spso,spsio
c	if(spsfo.gt.1.e-10)sbof6=sbof(6)/spsfo
c	write(10,*)'J/psi S fac. without and with D=',sbof6,sbof6*degr
c	if(spsifo.gt.1.e-10)sbof7=sbof(7)/spsifo
c	write(10,*)'psi S fac. without and with D=',sbof7,sbof7*degr
c        sbof67=sbof(6)+sbof(7)
c        spsfot=spsfo+spsifo
c        if(spsfot.gt.1.e-10)sboft=sbof67/spsfot
c        write(10,*)'S fac. without and with D=',sboft,sboft*degr
c	if(spso.gt.1.e-10)sbo6=sbo(6)/spso
c	write(10,*)'J/psi S fac. without and with D=',sbo6,sbo6*degr
c	if(spsio.gt.1.e-10)sbo7=sbo(7)/spsio
c	write(10,*)'psi S fac. without and with D=',sbo7,sbo7*degr
c       sbo67=sbo(6)+sbo(7)
c        spsot=spso+spsio
c        if(spsot.gt.1.e-10)sbot=sbo67/spsot
c        write(10,*)'S fac. without and with D=',sbot,sbot*degr
c	write(10,*)'J/psi transverse energy='
c	write(10,*)(feji(i1),i1=1,10)
c	write(10,*)(eji(i1),i1=1,10)
c	write(10,*)'J/psi transverse momentum squared='
c	write(10,*)(fpji(i1),i1=1,10)
c	write(10,*)(pji(i1),i1=1,10)
c	write(10,*)'psi transverse energy='
c	write(10,*)(fejpi(i1),i1=1,10)
c	write(10,*)(ejpi(i1),i1=1,10)
c	write(10,*)'psi transverse momentum squared='
c	write(10,*)(fpjpi(i1),i1=1,10)
c	write(10,*)(pjpi(i1),i1=1,10)
c	endif

csa****************************************************************
c	write(10,*)'time dependent gluon v1'
c	write(10,601)(sv1o(1,i2),i2=1,itnum)
c	write(10,*)'error of gluon v1'
c	write(10,601)(sv12o(1,i2),i2=1,itnum)
c	write(10,*)'time dependent gluon v2'
c	write(10,601)(sv2o(1,i2),i2=1,itnum)
c	write(10,*)'error of gluon v2'
c	write(10,601)(sv22o(1,i2),i2=1,itnum)
c	write(10,*)'time dependent u pair v1'
c	write(10,601)(sv1o(2,i2),i2=1,itnum)
c	write(10,*)'error of u pair v1'
c	write(10,601)(sv12o(2,i2),i2=1,itnum)
c	write(10,*)'time dependent u pair v2'
c	write(10,601)(sv2o(2,i2),i2=1,itnum)
c	write(10,*)'error of u pair v2'
c	write(10,601)(sv22o(2,i2),i2=1,itnum)
	write(10,*)'final state v1(pt)'
	do i1=1,21
	write(10,*)'i1=',i1
	write(10,602)(sv1o_pt(i1,i2),i2=1,20)
	enddo
	write(10,*)'error of v1(pt)'
	do i1=1,21
	write(10,*)'i1=',i1
	write(10,602)(sv12o_pt(i1,i2),i2=1,20)
	enddo
	write(10,*)'final state v2(pt)'
	do i1=1,21
	write(10,*)'i1=',i1
	write(10,602)(sv2o_pt(i1,i2),i2=1,20)
	enddo
	write(10,*)'error of v2(pt)'
	do i1=1,21
	write(10,*)'i1=',i1
	write(10,602)(sv22o_pt(i1,i2),i2=1,20)
	enddo
	write(10,*)'final state v1(eta)'
	do i1=1,21
	write(10,*)'i1=',i1
	write(10,602)(sv1o_eta(i1,i2),i2=1,20)
	enddo
	write(10,*)'error of v1(eta)'
	do i1=1,21
	write(10,*)'i1=',i1
	write(10,602)(sv12o_eta(i1,i2),i2=1,20)
	enddo
	write(10,*)'final state v2(eta)'
	do i1=1,21
	write(10,*)'i1=',i1
	write(10,602)(sv2o_eta(i1,i2),i2=1,20)
	enddo
	write(10,*)'error of v2(eta)'
	do i1=1,21
	write(10,*)'i1=',i1
	write(10,602)(sv22o_eta(i1,i2),i2=1,20)
	enddo
601	format(8(1x,f7.4))
602     format(6(1x,e10.3))
c021207
        write(10,*)'final state v1(pt) s'
        do i1=1,21
        write(10,*)'i1=',i1
        write(10,602)(ssv1o_pt(i1,i2),i2=1,20)
        enddo
        write(10,*)'error of v1(pt) s'
        do i1=1,21
        write(10,*)'i1=',i1
        write(10,602)(ssv12o_pt(i1,i2),i2=1,20)
        enddo
        write(10,*)'final state v2(pt) s'
        do i1=1,21
        write(10,*)'i1=',i1
        write(10,602)(ssv2o_pt(i1,i2),i2=1,20)
        enddo
        write(10,*)'error of v2(pt) s'
        do i1=1,21
        write(10,*)'i1=',i1
        write(10,602)(ssv22o_pt(i1,i2),i2=1,20)
        enddo
        write(10,*)'final state v1(eta) s'
        do i1=1,21
        write(10,*)'i1=',i1
        write(10,602)(ssv1o_eta(i1,i2),i2=1,20)
        enddo
        write(10,*)'error of v1(eta) s'
        do i1=1,21
        write(10,*)'i1=',i1
        write(10,602)(ssv12o_eta(i1,i2),i2=1,20)
        enddo
        write(10,*)'final state v2(eta) s'
        do i1=1,21
        write(10,*)'i1=',i1
        write(10,602)(ssv2o_eta(i1,i2),i2=1,20)
        enddo
        write(10,*)'error of v2(eta) s'
        do i1=1,21
        write(10,*)'i1=',i1
        write(10,602)(ssv22o_eta(i1,i2),i2=1,20)
        enddo
c021207
csa****************************************************************

	do m2=1,isdmax
	write(10,*)'ID of distribution m2=',m2
	do m3=1,ispmax
	write(10,*)'distribution belong to m3=',m3
	write(10,*)(sao(m1,m2,m3),m1=1,20)
	write(10,*)(saof(m1,m2,m3),m1=1,20)
	enddo
	enddo
	write(10,*)'average frequency of the occurring of each inela.'
	write(10,*)dineli

c	sumine=0.
c	do i1=1,600
c	sumine=sumine+dineli(i1)
c	enddo
c	write(10,*)'average total number of inela. colli.=',sumine
c300404
        write(10,*)'eevpav,eevhav=',eevpao,eevhao   ! 271205
	asd2=asd(2)
	asd22=asd(2)/2.
	do m1=1,20
	asdd(m1)=asd22+(m1-1)*asd2
	enddo
	do m1=1,20
	write(41,*)asdd(m1),sao(m1,2,1)
	enddo
	do m1=1,20
	write(42,*)asdd(m1),sao(m1,2,2)
	enddo
	do m1=1,20
	write(43,*)asdd(m1),sao(m1,2,3)
	enddo
	do m1=1,20
	write(44,*)asdd(m1),sao(m1,2,4)
	enddo
	do m1=1,20
	write(45,*)asdd(m1),sao(m1,2,5)
	enddo
	do m1=1,20
	write(46,*)asdd(m1),sao(m1,2,6)
	enddo
	do m1=1,20
	write(47,*)asdd(m1),sao(m1,2,7)
	enddo
	do m1=1,20
	write(48,*)asdd(m1),sao(m1,2,8)
	enddo
	do m1=1,20
	write(49,*)asdd(m1),sao(m1,2,9)
	enddo
	do m1=1,20
	write(50,*)asdd(m1),sao(m1,2,10)
	enddo
	do m1=1,20
	write(51,*)asdd(m1),sao(m1,2,11)
	enddo
	do m1=1,20
	write(52,*)asdd(m1),sao(m1,2,12)
	enddo
	do m1=1,20
	write(53,*)asdd(m1),sao(m1,2,13)
	enddo
	do m1=1,20
	write(54,*)asdd(m1),sao(m1,2,14)
	enddo
	do m1=1,20
	write(55,*)asdd(m1),sao(m1,2,15)
	enddo
	do m1=1,20
	write(56,*)asdd(m1),sao(m1,2,16)
	enddo
	do m1=1,20
	write(57,*)asdd(m1),sao(m1,2,17)
	enddo
	do m1=1,20
	write(58,*)asdd(m1),sao(m1,2,18)
	enddo
	do m1=1,20
	write(59,*)asdd(m1),sao(m1,2,19)
	enddo
	do m1=1,20
	write(60,*)asdd(m1),sao(m1,2,20)
	enddo
	do m1=1,20
	write(61,*)asdd(m1),saof(m1,2,1)
	enddo
	do m1=1,20
	write(62,*)asdd(m1),saof(m1,2,2)
	enddo
	do m1=1,20
	write(63,*)asdd(m1),saof(m1,2,3)
	enddo
	do m1=1,20
	write(64,*)asdd(m1),saof(m1,2,4)
	enddo
	do m1=1,20
	write(65,*)asdd(m1),saof(m1,2,5)
	enddo
	do m1=1,20
	write(66,*)asdd(m1),saof(m1,2,6)
	enddo
	do m1=1,20
	write(67,*)asdd(m1),saof(m1,2,7)
	enddo
	do m1=1,20
	write(68,*)asdd(m1),saof(m1,2,8)
	enddo
	do m1=1,20
	write(69,*)asdd(m1),saof(m1,2,9)
	enddo
	do m1=1,20
	write(70,*)asdd(m1),saof(m1,2,10)
	enddo
	do m1=1,20
	write(71,*)asdd(m1),saof(m1,2,11)
	enddo
	do m1=1,20
	write(72,*)asdd(m1),saof(m1,2,12)
	enddo
	do m1=1,20
	write(73,*)asdd(m1),saof(m1,2,13)
	enddo
	do m1=1,20
	write(74,*)asdd(m1),saof(m1,2,14)
	enddo
	do m1=1,20
	write(75,*)asdd(m1),saof(m1,2,15)
	enddo
	do m1=1,20
	write(76,*)asdd(m1),saof(m1,2,16)
	enddo
	do m1=1,20
	write(77,*)asdd(m1),saof(m1,2,17)
	enddo
	do m1=1,20
	write(78,*)asdd(m1),saof(m1,2,18)
	enddo
	do m1=1,20
	write(79,*)asdd(m1),saof(m1,2,19)
	enddo
	do m1=1,20
	write(80,*)asdd(m1),saof(m1,2,20)
	enddo
	close(41)
	close(42)
	close(43)
	close(44)
	close(45)
	close(46)
	close(47)
	close(48)
	close(49)
	close(50)
	close(51)
	close(52)
	close(53)
	close(54)
	close(55)
	close(56)
	close(57)
	close(58)
	close(59)
	close(60)
	close(61)
	close(62)
	close(63)
	close(64)
	close(65)
	close(66)
	close(67)
	close(68)
	close(69)
	close(70)
	close(71)
	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
	close(77)
	close(78)
	close(79)
	close(80)
c300404
	close(10)
	endif

1000	if(iii.lt.neve)then
	if(dabs(bmin-bmax).lt.10d-4)goto 300
	if(jjj.ge.10)then
c	10: total number of impact paremeters in systematic sampling for impact
c           parameter
	jjj=1
	goto 300
	endif
	jjj=jjj+1
	goto 300
	endif
	
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
	subroutine stati(y,pt,eta,p5,ik,kk,ww,a,b,af,bf)
c	on line statistics
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &  iabsb,iabsm,non10,ajpsi,csspn,csspm
	dimension a(20),b(20,5,20),c(5),af(20),bf(20,5,20),id(5)
	amass=p5   ! 010600
        amass2=amass*amass
        pt2=pt*pt
        et=dsqrt(pt2+amass2)
        do 10000 i=1,iflmax
        goto (10,20,30,40,50) i
10      c(i)=y   
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
100	ii=dabs(y)/asd(i)+1
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
300	ii=dabs(eta)/asd(i)+1
	if(ifram.eq.1 .and. eta.gt.0.)ii=ii+10
        if(ifram.eq.1 .and. eta.lt.0.)ii=10-ii+1
c       note: 10 here should be change together with the dimension
c        20
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine flow_f(px,py,pt,eta,p5,ik,kk)
c	calculate pt (eta) dependent direct and elliptic flow at final
c        state (partonic or hadronic)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=40000,KSZ1=30)
c210607        common/pyjets/nsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     c  iabsb,iabsm,non10,ajpsi,csspn,csspm
        common/sa18_pt/snum_pt(21,20),v1_pt(21,20),v2_pt(21,20),
     c   v12_pt(21,20),v22_pt(21,20)   ! 280607
	common/sa18_eta/snum_eta(21,20),v1_eta(21,20),v2_eta(21,20),
     c   v12_eta(21,20),v22_eta(21,20)   ! 280607
        common/sa24/adj1(40),nnstop,non24,zstop   ! 120707
c	snum(1,i1)=snum(1,i1)+1.
c	write(9,*)'px,py,pt,p5,ik,kk=',px,py,pt,p5,ik,kk   ! sa
	amass=p5   
        amass2=amass*amass
        pt2=pt*pt
	iik=iabs(ik)
        adj140=adj1(40)

c       pt is located at which interval ?
	idpt=pt/asd(2)+1
c       eta is located at which interval ?
        ii=dabs(eta)/asd(3)+1
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
        snum_pt(kk,idpt)=snum_pt(kk,idpt)+1   ! 280607
	v1_pt(kk,idpt)=v1_pt(kk,idpt)+pxt
	v2_pt(kk,idpt)=v2_pt(kk,idpt)+pxt2
	v12_pt(kk,idpt)=v12_pt(kk,idpt)+pxt*pxt
	v22_pt(kk,idpt)=v22_pt(kk,idpt)+pxt2*pxt2
c120607 direct and elliptic flows of charged particle in hadronic final state
        if(adj140.eq.4)then   ! 120707
	if(iik.eq.211 .or. iik.eq.321 .or. iik.eq.2212)then
        snum_pt(21,idpt)=snum_pt(21,idpt)+1   ! 280607
	v1_pt(21,idpt)=v1_pt(21,idpt)+pxt
	v2_pt(21,idpt)=v2_pt(21,idpt)+pxt2
	v12_pt(21,idpt)=v12_pt(21,idpt)+pxt*pxt
	v22_pt(21,idpt)=v22_pt(21,idpt)+pxt2*pxt2
        endif
        endif   ! 120707
c120707
        if(adj140.ne.4)then
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
        snum_eta(kk,idmt)=snum_eta(kk,idmt)+1   ! 280607
	v1_eta(kk,idmt)=v1_eta(kk,idmt)+pxt
        v2_eta(kk,idmt)=v2_eta(kk,idmt)+pxt2
        v12_eta(kk,idmt)=v12_eta(kk,idmt)+pxt*pxt
        v22_eta(kk,idmt)=v22_eta(kk,idmt)+pxt2*pxt2
        if(adj140.eq.4)then   ! 120707
        if(iik.eq.211 .or. iik.eq.321 .or. iik.eq.2212)then
        snum_eta(21,idmt)=snum_eta(21,idmt)+1   ! 280607
	v1_eta(21,idmt)=v1_eta(21,idmt)+pxt
        v2_eta(21,idmt)=v2_eta(21,idmt)+pxt2
        v12_eta(21,idmt)=v12_eta(21,idmt)+pxt*pxt
        v22_eta(21,idmt)=v22_eta(21,idmt)+pxt2*pxt2
	endif
        endif   ! 120707
c120707
        if(adj140.ne.4)then   
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine stapsi(nspsf,nspsif,nsps,nspsi)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &  iabsb,iabsm,non10,ajpsi,csspn,csspm
        common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
	dimension c(5)
	nspsf=nspsf+nps
	nspsif=nspsif+npsi
cc	goto 100   ! 28/04/98
c        if(iflmax.eq.0)goto 100
c       iflmax = 0 : no filter at all
	amass=pmas(pycomp(443),1)
	do i=1,nps
	px=pps(i,1)
	py=pps(i,2)
	pz=pps(i,3)
	pp=px*px+py*py+pz*pz
        e=dsqrt(pp)
        ee=dsqrt(amass*amass+pp)
        ppp=e-pz
        if(ppp.le.1.e-20)ppp=1.e-20
	ppy=ee-pz   ! 081102
        if(ppy.le.1.e-20)ppy=1.e-20   !081102
        eee=e+pz
	eey=ee+pz   !081102
	if(eee.lt.1.e-20)then
	write(9,*)'in stapsi eee,e=',eee,e
	write(9,*)'some thing wrong p=',px,py,pz
	eee=1.e-20
	endif
        y=0.5*dlog(eey/ppy)   ! 081102
        eta=0.5*dlog(eee/ppp)   ! 081102
	pt=dsqrt(px**2+py**2)
	c(1)=y
	if(ifram.eq.1)c(1)=eta   !081102
	c(2)=pt
c	.
c	.
c	.
        do i1=1,iflmax
        if(c(i1).lt.afl(6,i1,1) .or. c(i1).gt.afl(6,i1,2))goto 200
        enddo
	nsps=nsps+1
200	enddo
	amass=pmas(pycomp(30443),1)
	do i=1,npsi
	px=ppsi(i,1)
	py=ppsi(i,2)
	pz=ppsi(i,3)
	pp=px*px+py*py+pz*pz
        e=dsqrt(pp)
        ee=dsqrt(amass*amass+pp)
        ppp=e-pz
        if(ppp.le.1.e-20)ppp=1.e-20
	ppy=ee-pz
        if(ppy.le.1.e-20)ppy=1.e-20   ! 081102
        eee=e+pz
        eey=ee+pz   !081102
	if(eee.lt.1.e-20)then
	write(9,*)'in stapsi eee,e=',eee,e
	write(9,*)'some thing wrong p=',px,py,pz
	eee=1.e-20
	endif
        y=0.5*dlog(eey/ppy)
        eta=0.5*dlog(eee/ppp)   ! 081102
	pt=dsqrt(px**2+py**2)
	c(1)=y
	if(ifram.eq.1)c(1)=eta   ! 081102
	c(2)=pt
c	.
c	.
c	.
        do i1=1,iflmax
        if(c(i1).lt.afl(7,i1,1) .or. c(i1).gt.afl(7,i1,2))goto 300
        enddo
	nspsi=nspsi+1
300	enddo
100	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine adjins(r1,r2,r3,p1,p2,p3,iadj)
c       does the point of (r1,r2,r3) inside the sphere with radius 2
c        fm and located at coor(3)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter(kszj=40000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
c        x1=r1-coor(1)
c        x2=r2-coor(2)
c        x3=r3-coor(3)
        x1=r1-coor(1)
        x2=r2-coor(2)
        x3=r3-coor(3)
	ax1=dabs(x1)
	ax2=dabs(x2)
	ax3=dabs(x3)
	ap1=dabs(p1)
	ap2=dabs(p2)
	ap3=dabs(p3)
        if(ax1.le.0.5 .and. ax2.le.0.5 .and. ax3.le.0.5 .and. 
     c	 ap1.le.0.25 .and. ap2.le.0.25 .and. ap3.le.0.25)iadj=1
        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shanul(kf,x,q2,xpq)
c       calculate nuclear ratio R^A_i=f_{i/A}(x,Q2)/f_i(x,Q2) according
c        to Xin-Nian Wang's paper (PL, B527(2002)85), multiply it to
c        the parton distribution function in pythia, resulted parton distribution 
c        function is including nuclear shadowing effect 
c       it was proved in Eur. Phys. J. C9(1999)61 that nuclear ratio does
c        not depend strongly on the choice for the parton distribution
c        function in nucleon f_i(x,Q2)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        dimension xpq(-25:25)
        sq=0.1
        sg=0.26
        x35=x**0.35
        xsqr=dsqrt(x)
        x2=x*x
        x3=x2*x
        xx=x3-1.2*x2+0.21*x
        a=nap
c       what is the definition of "a" for asymmetry reaction system ?
        a13=a**0.3333
        aa=(a13-1)**0.6
        coa=dlog(a)
        coa16=coa**0.16666
        bbq=1.-3.5*xsqr
        bbg=1.-1.5*x35
        eq=dexp(-x2/0.01)
        eg=dexp(-x2/0.004)
c       raq=a*(1.+1.19*coa16*xx-sq*aa*bbq*eq)
c       rag=a*(1.+1.19*coa16*xx-sg*aa*bbg*eg)
        raq=1.+1.19*coa16*xx-sq*aa*bbq*eq
        rag=1.+1.19*coa16*xx-sg*aa*bbg*eg
c	write(9,*)'nap,x,q2,coa=',a,x,q2,coa   ! sa
c	write(9,*)'coa16,xx,aa=',coa16,xx,aa   ! sa
c	write(9,*)'bbq,eq,bbg=',bbq,eq,bbg   ! sa
c	write(9,*)'eg,raq,rag=',eg,raq,rag   ! sa
c       write(9,*)'before corre., R=',(xpq(i),i=-2,2)   ! sa
        do i=-6,6
	if(i.eq.0)then
	xpq(i)=rag*xpq(i)
	goto 100
	endif
	xpq(i)=raq*xpq(i)
100	enddo
c       write(9,*)'after corre., R=',(xpq(i),i=-2,2)   ! sa
        return
        end



c********************************************************************
        subroutine trans_h
c       'sa1_h' to 'pyjets' after finish calculation
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (kszj=40000,KSZ1=30)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nsa,non1_h,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        do l=1,nsa
        do m=1,5
        p(l,m)=psa(l,m)
        v(l,m)=vsa(l,m)
        enddo
        do m=1,2
        k(l,m)=ksa(l,m)
        enddo
        do m=3,5
        k(l,m)=0
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine oscar(iii)
c       record history of spatial and momentum coordinates due to
c        OSC1999A
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=40000,KSZ1=30)
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
	parameter (kszj=40000)
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
c	move particles with flavor code ii ('55') from 'lujets' to 'sgam'
	parameter (kszj=40000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
	common/ssin/nsin,nonsin,ksin(kszj,5),psin(kszj,5),vsin(kszj,5),
     c   bsin(kszj)   ! 140506
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
c	bsin(n)=0.
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
	bsin(j1)=bsin(j)
        enddo
        n=n-1
        goto 201
202     enddo
203     continue
        return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine remo_gam(ii)   ! 250209
c	move particles with flavor code ii ('66') from  'lujets' to 'sgam'
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sgam(nn)   ! 250209
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
        write(22,*)'c & p sum=',ch1/3.,peo   !
        do i=1,nn
        write(22,*)i,kgam(i,1),kgam(i,2),(pgam(i,j),j=1,4)
        enddo
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
c D.Miskowiec 1997, updated in 2001, http://www-linux.gsi.de/~misko/overlap/ . 
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
        RA=1.12*A**0.333333   ! -0.86/A**0.333333   ! 020511
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
