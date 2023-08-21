        program main_23
c       administrates MC simulation for relativistic pp(lp), pA (Ap), 
c        AB (lA), & e+e- collisions 
c100223 the program composes of main_23.f,parini_23.f,parcas_23.f, 
c        p_23.f,sfm_23.f,coales_23.f,hadcas_23.f,and analy.f
c       main_23.f: administrates the MC simulation
c       parini_23.f: generates partonic initial state of colliding system
c       parcas_23.f: performs parton rescattering, where 2->2 processes
c        are considered only and LO pQCD cross section or its regularized
c        approximation is used
c       p_23.f (pythia 6.4.28): for generation of partonic initial state 
c        and/or string hadronization
c       sfm_23.f: hadronization with LUND string fragmentation model
c       coales_23.f: hadronization with Monte Carlo coalescence model
c       hadcas_23.f: performs hadronic rescattering
c       analy.f: an example of event analysis subroutine, user is free
c        to replace it with your own one

c250420 note: because of history reason, all of the Fortran programs are 
c        written with mode of regardless either capital lette or small 
c        letter. typing ": set ic + enter key" before searching in Vi/Vim.

        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MPLIS=80000)
        PARAMETER (NSIZE=280000)   !Lei2023060
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
        COMMON/PYDATR/MRPY(6),RRPY(100)   ! 171022 Lei
        common/pyint1/mint(400),vint(400)
        common/pycidat1/kfacot(100),disdet(100),isinelt(600)
        common/pycidat2/kfmaxt,nont2,param(20),weigh(600)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(6),
     c   afl(20,6,2)   !Lei2023060 5 -> 6
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
        common/sa16/x_ratio,dni(10),dpi(10),edi(10),bmin,bmax
     &   ,bar(10),abar(10),barf(10),abarf(10)   ! 033101 Lei202307
     &   ,emin(10),eminf(10),eplu(10),epluf(10)   ! 033101
c00623 common/sa18/tdh,itnum,non18,cptl,cptu,cptl2,cptu2,snum(4,20),
c00623 &   v1(4,20),v2(4,20),v12(4,20),v22(4,20)   !Lei2023060
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph
     c   ,vnlep   ! 260314
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
     c   part1(200),part2(200),binn(200)   ! 020511 020718
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 220312 240412 131212
        common/sa34/itorw,iikk,cp0,cr0,kkii   ! 010418 010518 050920
        common/sa35/ncpart,ncpar(kszj)   ! 280722
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)   ! 150922
        common/sa38/csp_31,csp_32,csp_41,csp_42,csp_43,csp_51,csp_52,
     c   csp_53,csp_54,csp_61,csp_62,csp_63,csp_64,csp_65   ! 161022
        common/sa6_c/ithroq,ithrob,ithroc,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)   ! 210921
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5) ! 010518
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)   ! 050603
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)   ! 050603
        common/delt/ndel,nodel,kdel(kszj,5),pdel(kszj,5),vdel(kszj,5) ! 150323
        common/show/vip(mplis),xap(mplis)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
        common/count/isinel(600)
        common/ctllist/npctl,npinel(600),npctl0,npctlm ! 061103 180121 230121
        common/ctllist_p/nreac(9),nrel   ! 071103
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5) ! 250209
        common/anly1/an(40,6,20),bn(20),anf(40,6,20),bnf(20)   ! 281219 Lei2023060 5 -> 6
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5) ! 280620
        common/work7/reac(9),crose(9)   !Lei2023060 in parcas
        common/ancoal/icoal1,icoal2,xkappa1,xkappa2 ! 200820 yan
        common/schuds/schun,schudn,schudsn,sfra   ! 211022 She and Lei
        dimension san(40,6,20),sbn(20),sanf(40,6,20),sbnf(20)   ! 070419 Lei2023060 5 -> 6
        dimension sao(40,6,20),sbo(20),saof(40,6,20),sbof(20)   ! 070419 Lei2023060 5 -> 6
        dimension skapa(6),skapao(6),snreac(9)   ! 020708 220110 010518
        dimension c(5),dinel(600),dineli(600),sthroe(4),wthroe(4)
        dimension einel(600),eineli(600)   ! 140820
        dimension bpp(20),kdiq(kszj,5),dgmas(kszj)
        dimension acoll(20),acollp(20),acollt(20)
        dimension sbp(20),numbth(3)   ! 010718
        dimension nreaco(9),pl(100,5)   ! 220110 260314
        dimension ksin(kszj,5),psin(kszj,5),vsin(kszj,5)   ! 010518 230618
        dimension pi(4),pj(4),b(3)
        real*8 nmin,nminf,ncha,nchaf   ! 020203
        dimension lc(nsize,5),tc(nsize),tw(nsize)   !Lei202305

c171022 MRPY(1) is the seed of PYTHIA random number generator. (D=19780503)
c       For the intrinsic subroutine DATE_AND_TIME calling.
        dimension n_current_date_and_time(8)   ! 171022 Lei
        character*4 c_date_and_time(8)   !Lei2023060
        dimension ps0(6),ps1(6)   !Lei2023060

!Lei20230303B
        dimension sum_pT2(40,2),sum_pT(40,2),sum_y(40,2),sum_eta(40,2),
     &            i_h(6)
        dimension i_hist(100)
!Lei20230303D

c010418 itorw: =1 executing pyevnt, =2 executing pyevnw
c260314 pl(ii,5): four momentum and mass of ii-th lepton
c       adj1(i), i=
c       1: k factor used in parton rescattering,k=0: no parton rescattering
c       2: parameter \alpha_s ('as' in program) in parton rescattering
c       3: parameter (\mu)^2 ('tcut' in program): to avoid divergence in 
c          calculating parton-parton differential cross section in parcas_23.f
c       4: parameter idw, # of integration intervals in parcas_23.
c       5: =0, w/o nuclear shadowing,
c          =1, Wang's nuclear shadowing (PLB 527(2002)85).
c       6: parameter 'a', i.e. parj(41), in Lund string fragmentation function
c       7: parameter 'b', i.e. parj(42), in Lund string fragmentation function
c       8: i.e. MSTP(82) in PYTHIA 6.4
c       9: = PARP(81) (D=1.9 GeV/c), effective minimum transverse momentum 
c            p_erp_min of multiple interactions if MSTP(82)=1, for old MPI;
c          = PARP(82) (D=2.0 GeV/c), regularization scale p_erp_0 of the 
c            transverse-momentum spectrum for multiple interactions with 
c            MSTP(82) >= 2, for new MPI.   !Lei2023060
c       10: parp(31),k factor in pythia64
c       11: time accuracy used in hadas_23.f
c       12: model of hadronization: 
c           =0 string fragmentation; 
c           =1 coalescence mode 1; Lei2023060
c           =2 coalescence mode 2; with gluon splitting & deexcitation 
c                                  before parcas. Lei2023060
c       13: dimension of meson table considered if adj1(12) > 0
c       14: dimension of baryon table considered if adj1(12) > 0
c       15: default string tension
c       16: # of allowable generations in q (qbar) deexcitation in coales_23.f
c       17: threshold energy in deexcitation of energetic quark in coales_23.f
c       18: =0, rest partons hadronize by string fragmentation
c           =1, rest partons hadronize by coalescence
c       19: time accuracy used in parcas_23.f ('dddt' in program)
c       20: =0 exact pQCD parton-parton cross section
c           =1 limited and regularized parton-parton cross section (B. Zhang,
c              Comput. Phys. Commun. 109(1998)193)
c           =2 the same as 0 but flat scattering angle distribution is assumed
c           =3 the same as 1 but flat scattering angle distribution is assumed
c       21: =0 without phase space requirement in coales_23.f
c070223     =1 with complete phase space requirement
c070223     =2 with phase space requirement in spatial phase space only 
c110223     =3 with phase space requirement in momentum phase space only   ! 110223 Lei
c       22: critical value of the product of radii both in coordinate and 
c        momentum phase space (4 is assumed) used in coales_23.f
c       23: switch for chiral magnetic effect (CME)
c           =0: CME off
c           =1: CME on
c       24: =tl0,the virtuality cut in time-like radiation in parcas_23.f, 
c            4*tl0 is assumed
c       25: \Lambda_QCD in parcas_23.f
c171022 26: selection of random number seed   ! 171022 Lei
c           =0, default PYTHIA seed (19780503), can be used for debug
c           =other, seed from the real-time colock
c       27: largest momentum allowed for particle ('dpmax')
c       28: largest position allowed for particle in hadcas ( which is 1. in 
c            usu.dat, but is recalculated in the running, 
c            drmax=para10*max(rnt,rnp)*adj1(28) )   Lei2023060
c140223 29: For sfm in PYTHIA, it is MSTJ(11). Choice of longitudinal   ! 140223 Lei
c            fragmentation function, i.e. how large a fraction of the energy
c            available a newly-created hadron takes.
c           For coal, sampling deexcited daughter qqbar-pair energy fraction z 
c            taking from mother by 'PYZDIS' in PYTHIA or 'funcz' or random z.
c           =1: Lund symmetric fragmentation function, see PARJ(41) - PARJ(45)
c           =2: Field–Feynman + Peterson/SLAC, see PARJ(51) PARJ(59)
c           =3: Lund + Peterson/SLAC (light flavor + heavier)
c           =4: default PYTHIA. Lund + Bowler
c           =5: as = 4, but interpolate for c and b; see PARJ(46) and PARJ(47)
c           =0: as = 4 in sfm / random z in 'coal'   !Lei2023060
c           Note: for coal, when using 11, 12 or 13, it is LUND/FF/PS
c140223      fragmentation function in 'funcz' which comes up with PACIAE.
c111222 30: =1: distributes participant nucleons in overlapping areas forcely
c           =0: no more requirements   ! 111222 Lei
c       31: parj(1) in pythia64
c       32: parj(2) in pythia64
c       33: parj(3) in pythia64
c       34: parj(21) in pythia64
c       35: mstp(91) in pythia64, selects parton transverse momentum 
c           (k_{\perp}) distribution inside hadron; 
c           =1, Gaussian; 
c           =2, exponential
c       36: =0 without phenomenological parton energy loss in parcas_23.f
c           =1 otherwise
c       37: the coefficient ('c') in phenomenological parton energy loss
c       38: pt cut in phenomenological parton energy loss 
c       39: width of Gaussian k_{\perp} distribution in hadron if mstp(91)=1
c           width of exponential k_{\perp} distribution in hadron if mstp(91)=2
c       40: =1 simulation ends after partonic initiation
c           =2 simulation ends after partonic rescattering
c           =3 hadronization by coalesce model follows partonic initiation
c           =4 simulation ends after hadron rescattering

c       smadel: small perpurbation of ellipse relative to circle
c       parecc: a parameter converting initial spatial space eccentricity 
c        to final momentum space ellipticity
c       iparres: =0 considers elastic parton-parton collisions only in 
c                   parcas_23.f
c                =1 with inelastic parton-parton collisions as well

c       pincl (pscal): four momentum and mass of incident (scatterd) lepon
c       pinch: four momentum and mass of incident hadron
c        vnu: \nu; fq2: Q^2=-q^2; w2l: W^2; yyl: y; zl: z; xb: x_B; pph: P_h

c       para1_1: nn total cross section, used in parini_23.f 
c       para1_2: nn total cross section  used in hadcas_23.f
c       dni: nucleon number density

c       pj and ej: (pt)**2 and transverse energy of J/psi 
c       pjp and ejp: (pt)**2 and transverse energy of (J/psi) prime

c       acoll: a array, demension of which should be larger
c        than or equal to 'nmax' (interval # of impact parameter)
c       the dimension of 'bpp' must be small than or equal to 'nmax'

c       ipden: =0, if projectile is proton
c              =2, for e+e- collisions  ! 180921 yan
c              =1, if projectile is nucleus
c              =11, projectile is e- (e+)  
c              =12, projectile is nu_e (nu_ebar)  
c              =13, projectile is mu- (mu+)  
c              =14, projectile is nu_mu (nu_mubar)
c              =15, projectile is tau- (tau+)  
c              =16, projectile is nu_tau (nu_taubar)  
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
c       bp: impact parameter (=0 for pp,lp and lA)
c       iii: current run number
c       coor: CM position of collision system
c       ispmax: maximum # of particle KF code wanted to statistics
c       ispkf(i): KF code of i-th particle wanted to statistics
c       kfmax: the maximum # of particle KF code considered 
c       kfaco(i): KF code of i-th particle
c       numb(i): order number of last particle among particles with same flavor 
c        code of kfaco(i) in particle list

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
c       isinel(i): = 0 without i-th inelastic process in hadcas_23.f
c                  = 1 with i-th inelastic process in hadcas_23.f

c       nrel: statistics of blocked parton-parton scattering process in
c        parcas_23.f
c       nreac(i): statistics of successful i-th parton-parton scattering 
c        process in parcas_23.f
c       npinel(592): # of nn collision calling pythia' in parini_23.f 
c       npinel(593): # of nn collision not calling pythia' in parini_23.f
c       noel : statistics of elastic collisions in hadcas_23.f
c       noinel(i): statistics the i-th inelastic collisions in hadcas_23.f
c140223 nosc = 0 : no OSCAR output, see subroutine oscar   ! 140223 Lei
c            = 1 : OSCAR1997A (final_id_p_x, just PACIAE final output)
c            = 2 : OSCAR1999A (full_event_history)
c140223      = 3 : OSCAR2013A (full_event_history), dummy now
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

c161022 cumulent sum of flavor generation probability, for an example: 
c        p_31=(1-amd/sms)/2., p_32=(1-amu/sms)/2., p_33=(1-ams/sms)/2.; 
c        amd (amu, ams): d (u,s) quark constituent mass, sms=amd+amu+ams; 
c        p_31+p_32+p33=1;
c        cumulent sum: csp_31=p_31, csp_32=p_31+p_32
c       u and d have same mass and probability

c       ?????????????????????????????????????
c       0. is the hard distance between two pions
c       0.5 is the hard distance between two nucleons
c       0. is the hard distance between pion and nucleon
c       ??????????????????????????????????????

c       win: =incident momentum if ifram=0 
c            =sqrt(s) if ifram=1

c250209 flavor code 22: hardonic decay photon
c                   44: prompt direct photon (<- pythia)
c                   55: photon from parton-parton scattering
c                       qg->q(gamma) and q(-q)->g(gamma)
c                   66: hardonic direct photon
c250209                 pi+pi->rho+(gamma) and pi+rho->pi+(gamma)
c240219             77: photons from hadronization


2000    open(2,file='prx.out',status='unknown')
        open(3,file='pix.out',status='unknown')
        open(11,file='usu.dat',status='unknown')
        mstu(11)=22
        open(22,file='main.out',status='unknown')
c       open(98,file='databs_g.dat',status='unknown')   ! 260219
c260219 databs_g.dat: record gammas E-by-E
c       open(99,file='databs_h.dat',status='unknown')   ! 240119 260219
c260219 databs_h.dat: record hadrons E-by-E
c120214 if projectile is lepton, set nap =1


c-------------------------------------------------------------------------------
c------------------------------   Input Reading   ------------------------------
c       reads input variables for event generation
        read(11,*)neve,nout,nosc
        read(11,*)nap,nzp,nat,nzt
        read(11,*)ddt,x_ratio,bmin,bmax,nmax   ! 241108 Lei202307
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
c00623 read(11,*)tdh,cptl,cptu,cptl2,cptu2,itnum   ! 241108 Lei2023060
        read(11,*)i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,a_FF,aPS_c,aPS_b,parj23,parj24,
     &   parp82,i_coord_recover,i_tune   !Lei2023060
        read(11,*)mstu21,mstj1_1,mstj1_2,mstj2,decpro,itorw   ! 160617 010418 Lei202307
        read(11,*)(adj1(i),i=1,10)
        read(11,*)(adj1(i),i=11,20)
        read(11,*)(adj1(i),i=21,30)
        read(11,*)(adj1(i),i=31,40)
        read(11,*)kjp22,kjp23,kjp24,parp78,mstptj   !  100821 230722
        read(11,*)parecc,iparres,smadel,dparj4,cp0,cr0,seco   ! 120219 260219
        read(11,*)csp_31,csp_32   ! 161022
        read(11,*)csp_41,csp_42,csp_43   ! 161022
        read(11,*)csp_51,csp_52,csp_53,csp_54   ! 161022
        read(11,*)csp_61,csp_62,csp_63,csp_64,csp_65   ! 161022
!       h+- histogram
        read(11,*) pT_low, pT_upp, y_low, y_upp, eta_low, eta_upp   !Lei20230313
        close(11)
c------------------------------   Input Reading   ------------------------------
c-------------------------------------------------------------------------------

!Lei20230303B
        do i=1,100,1
            i_hist(i) = i
        end do
        i_h = 0
        i_h(1)= 211
        i_h(2)=-211
        i_h(3)= 321
        i_h(4)=-321
        i_h(5)= 2212
        i_h(6)=-2212
    !     pt_max = 5.
    !     CALL PYBOOK(1, 'px distribution', 100,  -pt_max, pt_max)
    !     CALL PYBOOK(2, 'py distribution', 100,  -pt_max, pt_max)
    !     CALL PYBOOK(3, 'pT distribution', 100,  0D0, pt_max)
    !     CALL PYBOOK(4, 'z distribution of uds' , 100,  0D0, 1D0 )
    !     CALL PYBOOK(5, 'z distribution of c' , 100,  0D0, 1D0 )
    !     CALL PYBOOK(6, 'z distribution of b' , 100,  0D0, 1D0 )
    !     CALL PYBOOK(66, 'z distribution, total' , 100,  0D0, 1D0 )
    !     CALL PYBOOK(7, 'KF distribution af. gluon splitting' , 
    !  &                                      6,  0D0, 6.5D0 )
    !     CALL PYBOOK(8, 'KF distribution af. deexcitation' , 
    !  &                                      6,  0D0, 6.5D0 )
        CALL PYBOOK(9,  'h+-/2 inv. dN/dpT, partial' , 100,  0D0, 20D0 )
        CALL PYBOOK(10, 'h+-/2 inv. dN/dpT, full' ,    100,  0D0, 20D0 )
        CALL PYBOOK(11, 'h+- dN/y, partial' ,          100, -6D0,  6D0 )
        CALL PYBOOK(12, 'h+- dN/y, full' ,             100, -6D0,  6D0 )
        CALL PYBOOK(13, 'h+- dN/eta, partial' ,        100, -6D0,  6D0 )
        CALL PYBOOK(14, 'h+- dN/eta, full' ,           100, -6D0,  6D0 )
        CALL PYBOOK(15, 'h+-/2 dN/dpT, partial' , 100, 0D0, 20D0 )
        CALL PYBOOK(16, 'h+-/2 dN/dpT, full' ,    100, 0D0, 20D0 )
        CALL PYBOOK(17, 'time parini' ,           100, 0D0, 10D0 )
        CALL PYBOOK(18, 'time parcas' ,           100, 0D0, 100D0 )
        CALL PYBOOK(19, 'time hadcas' ,           100, 0D0, 100D0 )
        CALL PYBOOK(20, 'Ncoll parcas' ,          100, 0D0, 1000D0 )
        CALL PYBOOK(21, 'Ncoll hadcas' ,          100, 0D0, 1000D0 )
!Lei20230303D

c       bmin,bmax : minimum and maximum impact parameters, bmin=bmax means
c        definite impact parameter, 2*nmax: the number of 
c        intervals segmented in [bmin,bmax]
c       nchan=0: inelastic (INEL)
c       nchan=1: Non Single Difractive (NSD)
c       nchan=2:
c       nchan=3: J/psi production
c       nchan=4: heavy-flavor production
c       nchan=5: direct photon
c       nchan=6: soft only
c       nchan=7: W+/- production
c       nchan=8: PYTHIA default (msel=1)
c       nchan=9: Z0 production
c       setting of nchan=0,1,3,6,7,8 and 9 is ready

c       neve : # of events to be generate
c       nap (nzp) : # of nucleons (protons) in projectile nucleus
c       nat (nzt) : # of nucleons (protons) in target nucleus
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

c       t0 : average proper formation time at rest
c       ddt : time accuracy
c       dep : the accuracy in four momentum conservation
c       rou0 : normal nuclear density
c       rao : enlarge factor in the radius of simulated volume
c       kjp20: =1 constant cross sections in hadcas
c              =0 energy dependent cross sections in hadcas
c       kjp21: = 0 without hadron rescattering, 
c              = 1 with hadron rescattering
c
c       kjp22 = 1 : variable single string tension and parj(1) etc. 
c       kjp22 = 2 : variable multiple string tension and parj(1) etc.
c       kjp22 = 3 : variable (single+multiple) string tension and parj(1) etc.
c       kjp22 = 4 : default string tension and parj(1) etc.
c       kjp23: = 1 npart is calculated by geometric model
c       kjp23: = 2 npart is calculated by Glauber model
c       kjp24: = 1 sharp sphere in Glauber model
c       kjp24: = 2 Woods-Saxon in Glauber model
c260223 mstj2: =1, low energy simulation A-loop
c              =2, PYTHIA-like simulation B-loop
c              =3, PACIAE simulation C-loop
c260223 decpro is Delta decay probability in low energy A-loop   !Lei202307

c       pathn: collision numer suffered by projectile nucleon in target nucleus

c00623 i_deex: (D=1) the deexcitation mode used in coal   !Lei2023060
c               = 1, light-cone variable mode
c               = 2, energy mode
c       i_deex_gen: (D=0) the deexcitation generation of newly produced qqbar
c                   = 0, means no deexcitation for any newly produced qqbar
c                   = 1, means just do deexcitation for the directly proudced 
c                        qqbar pairs (1-st daughters) from original mother 
c                        quarks (Orig mothers)
c                   = 2, means do deexcitation for "1-st daughters" from 
c                        "Orig mothers" and the subsequent qqbar pairs 
c                        produced from "1-st daughters". (2-nd daughters)
c                     ...
c                   = 999, always do deexcitation for newly produced qqbar
c       i_pT: (D=3) the pT sampling method of the daughter qqbar pair in coal
c             = 1, Gaussian px and py with width PARJ(21)
c             = 2, Exponential px and py with width PARJ(21)
c             = 3, Exponential pT with width PARJ(21)
c             = 4, random pT from mother
c             = 5, random px and random py from mother, different random factors
c             = 6, random (px and py) from mother, the same random factor
c             = 7, random (px and py) from mother, the same random factor as 
c                   z which related to adj1(29).
c       i_pt_max: (D=0) whether the sampled pT in coal deexitation is greater 
c                       than the mother quark or not.
c       i_split_diq: (D=1) momentum splitting/allocating mode of diquark breakup
c       i_split_qqb: (D=2) momentum splitting mode of qqbar -> q + qbar in coal
c       i_split_g: (D=1) momentum splitting mode of g -> q + qbar in coal
c             = 1, decay mode, i.e. decmom + random 3-momentum method
c             = 2, random 3-momentum method with the different factors
c             = 3, random 3-momentum method  with the same factor
c             = 4, divided equally
c       a_FF: (D=0.77) parameter for light hadron in Field-Feynman function, 
c             i.e. u, d, and s --PARJ(51), (52), and (53)--, set them equal
c       aPS_c: (D=0.05) -PARJ(54), parameter for charm-hadron in Petersono/SLAC
c       aPS_b: (D=0.005) -PARJ(55), parameter for bottom-hadron in P/S function
c       parj23: (D=0.01) PARJ(23), non-uniform tail fraction
c00623 parj24: (D=2.) PARJ(24), width of tail = PARJ(21)*PARJ(24)   !Lei2023060


c-------------------------------------------------------------------------------
c-------------------------   Parameter Value Giving   --------------------------
c       param(1)=para1   ! given in 'sysini' in parini_23.f
        param(2)=para2   ! read from usu.dat
        param(4)=para4
        param(6)=x_ratio ! 240323 Lei202307
c240323 x_ratio: ratio of inela. cross section to total cross section, 
c        i.e. param(6) in hadcas_23.f, with default value of 
c        D=0.85 for high energy and D=0.1 for low energy
        if(mstj2.eq.1)then   !Lei2023060
        if(win.gt.2.015.and.win.lt.3.0)then
        param(6)=1.35*(win-2.015)**2/(0.015+(win-2.015))*x_ratio*10.0 ! 010423 yan Lei202307
        endif
c       x_ratio=1.35*(win-2.015)**2/(0.015+(win-2.015)) revised by yan
c       default x_ratio=0.1 in usu.dat, so there are x_ratio*10.0, 010423 yan
        endif

        param(7)=para7
        param(8)=ddt
        param(10)=para10
c       totle cross-section of J/Psi + n
        param(13)=para13
c       totle cross-section of J/Psi + meson
        param(14)=para14
c       totle cross-section of Psi' + n
        param(15)=para15
c       totle cross-section of Psi' + meson
        param(16)=para16
c020511 # of segments used in integration
        idw=adj1(4)
        mstp(82)=adj1(8)
        parp(81)=adj1(9)
        parp(82)=parp82   !Lei2023060 For new MPI model.
        parj(2)=adj1(32)
        parj2=parj(2)
        parj(1)=adj1(31)
        parj(3)=adj1(33)
        parj(4)=dparj4
        parj(21)=adj1(34)
        PARJ(23)=parj23   !Lei2023060
        PARJ(24)=parj24   !Lei2023060
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
c       parp(93)=adj1(30) ! upper cut-off for k_perp distribution in hadron
c       parj(21)=adj1(29)
        parp(2)=parp21
c       parp21: lowest CM energy for calling 'pythia' (D=10.), for 
c        the case of nchan not equal to 3
c140223 parp22: select y or eta in partial phase-space statistics.   ! 140223 Lei
c               = 0 , y
c140223         = 1 , eta

c       inclusion of k factor multiplying the differential cross sections for 
c       hard parton-parton process
        mstp(33)=1
        parp(31)=adj1(10)   ! D=1.5

c070417 control the strength of colour reconnection
        parp(78)=parp78

c       mstj1_1: =6, with inelastic processes 4, 6, and 7 (parcas_23.f)
c                =7, with inelastic process 7 only (parcas_23.f)
c       mstj1_2: =0, w/o final state time-like parton shower if iparres=1
c                =1, w/ final state time-like parton shower if iparres=1
c230722 mstptj: =0, input mstp(111) (mstj(1)) for pp,pA(Ap),and AA (for e+e-) 
c121222          in PACIAE simulation developed from partonic initial stage, 
c                to partonic rescattering, hadronization, and to hadronic 
c                rescttering stage
c       mstptj: =1, PYTHIA-like simulation without partonic & hadronic
c121222           rescatterings and low energy simulation   ! 260223

c       check on possible errors during execution
        mstu(21)=mstu21   ! 120603
c
c140223 Choses fragmentation function in PYTHIA for sfm and 'PYZDIS' for coal.
        MSTJ(11)=INT( adj1(29) )   ! 140223 Lei D=4 in PYTHIA
c00623 When use 'funcz' or random z in coal.   !Lei2023060
        if(MSTJ(11).lt.1 .OR. MSTJ(11).gt.5) MSTJ(11)=4   !Lei2023060
c       parameters in Lund string fragmentation function
        parj(41)=adj1(6)   ! D=0.3
        parj(42)=adj1(7)   ! D=0.58
c00623 parameters in Field-Feynman fragmentation function   !Lei2023060
        PARJ(51)=a_FF   ! u-, D=0.77
        PARJ(52)=a_FF   ! d-
        PARJ(53)=a_FF   ! s-, set them equal
c00623 parameters in Peterson/SLAC fragmentation function   !Lei2023060
        PARJ(54)=-aPS_c   ! c-, D=-0.05, note the minus sign
        PARJ(55)=-aPS_b   ! b-, D=-0.005

        adj12=adj1(12)
        adj18=adj1(18)
        adj140=adj1(40)
c-------------------------   Parameter Value Giving   --------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c-----------------------   Global Variable Initializing   ----------------------
c211022 fraction for charge separation    ! 211022 She and Lei
        schun=0.
        schudn=0.
        schudsn=0.
        sfra=0.

c       initiation of event averaged variales
        dnmin=0.
        dnminf=0.
        dncha=0.
        dnchaf=0.
        sbn=0.
        sbnf=0.
        san=0.
        sanf=0.

        stime_ini=0.
        stime_par=0.
        stime_had=0.
        snnc=0.
        sgtime=0.
        sgtimeo=0.
        sitime=0.
        sadiv=0.
        sgpmax=0.
        skapa=0.
        skapao=0.
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
        dinel=0.
        einel=0.
        acoll=0.
        acollp=0.
        acollt=0.
        sbp=0.

        volum=4.*3.1416/3.*2.**3
c261002 volume of sphere with radius of 2 fm in position phase space 
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
        sthroe=0.
c220110
        nrel=0
        nrea=0
        nreac=0
c220110
c280113
        if(INT(psno).eq.2)then
        averb=0.
        psnon=0.   ! N_bin in case of psno=2
        psnop=0.   ! parojectile N_part in case of psno=2
        psnot=0.   ! target N_part in case of psno=2
        endif
        nncoll=0
        vnlep=0.d0 ! statistics of the number of studied leptons 260314
c280113
        iii=0
        ich=0
        time=0.
        jjj=1
c----------------------   Global Variable Initializing   -----------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c----------------------------   Mistake-proofing   -----------------------------
c140223 Lei
c       Mistake-proofing.
c       Values given automatically avoiding some manual input mistakes in usu.dat.
        if( ipden.eq.2 .AND. itden.eq.2 )then
c       For e+e- , one just needs to specify ipden=itden=2 in usu.dat.
            nap =  1   ! nominal value
            nzp =  1   ! nominal value
            nat =  1   ! nominal value
            nzt = -1   ! nominal value
            bmin= 0D0
            bmax= 0D0
            psno= 0D0
            adj1(5) = 0D0   !Lei2023061
        else if( ipden.gt.1 )then
c       For lN and lA, one just needs to specify ipden >= 11 and nzp/nat/nzt.
            nap = 1   ! nominal value
            ! nzp = -1   ! nominal value for e-/mu-/tau-/nu_e...
            ! nzp = +1   ! nominal value for e+/mu+/tau+/nu_ebar...
            itden=0   ! lN
            if(nat.gt.1) itden=1 !lA
            if( itden.eq.0 ) then   ! lN
                bmin = 0D0
                bmax = 0D0
                psno = 0D0
            end if
        else
c       For NN , NA(AN) and AA, one just needs to specify nap/nzp/nat/nzt.
            ipden=0
            itden=0
            if(nap.gt.1) ipden=1
            if(nat.gt.1) itden=1
            if( ipden.eq.0 .AND. itden.eq.0 ) then   ! NN
                bmin = 0D0
                bmax = 0D0
                psno = 0D0
                adj1(5) = 0D0   !Lei2023061
            end if
        end if
        mstptj=0   !Lei2023060
        if(mstj2.eq.1 .OR. mstj2.eq.2) mstptj=1   !Lei2023060
        ! if(mstj2.eq.2 .OR. mstj2.eq.3) x_ratio=0.85   !Lei2023060 hardcode
        ! if(mstj2.eq.2 .OR. mstj2.eq.3) param(6)=0.85   !Lei2023060
c       if win.lt.parp21 & mstptj=0, dead loop in low-energy AA collision.
c260223 changed from (win.lt.parp21)
c140223 Lei
c----------------------------   Mistake-proofing   -----------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c--------------------------   System Initializing   ----------------------------
c       gives values to some important variables
        call sysini(win)   ! in parini.f
        adj1(28)=para10*dmax1(rnt,rnp)*adj1(28)   !Lei2023060 "*adj1(28)"

        open(5,file='sxp.out',status='unknown')   ! sa 26/05/99
        open(9,file='rms0.out',status='unknown')

c140223 nosc=0 : no oscar output
        if(nosc.gt.0) open(34,file='oscar.out',status='unknown')   ! 140223 Lei
        call oscar(win,-1)   ! 140223 Lei Prints oscar file header if nosc > 0.
        call PASTAT(-2,0)  !Lei202306 Initializes the variables couting cross sections.
c--------------------------   System Initializing   ----------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c-------------   Independent Optical Initial Geometry Calculation  -------------
c020511 calculates the nuclear overlap function etc.
c       if((ipden.eq.0 .and. itden.eq.0) .or. ipden.ge.11)goto 80001!180921 yan
        if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.
     c   itden.eq.2) .or. ipden.ge.11)goto 80001   ! 180921 yan
        csnn1=csnn*10   ! csnn in fm^2 csnn1 in mb
        idw1=idw/50   ! *100
        if(ipden.lt.2)call overlap(nap,nat,rnp,rnt,csnn1,kjp23,kjp24,
     c   rou0,idw1)   ! 060813 120214 changed from .ne. to .lt.
c020511

80001   continue   ! 310518
c       psno: =0 fixed impact parameter
c       psno: =1 systematic sampling method
c       psno: =2 random sampling method

c       calculate the impact parameter etc.
c       for given b (impact parameter)

c171022 Initialies them avoiding bad display in rms.out.   ! 171022 Lei
        avb    = 0D0
        avneu  = 0D0
        astbp  = 0D0
        astbt  = 0D0
        aanbin = 0D0
        swouni = 2D0

c260718 if(dabs(bmin-bmax).lt.10d-4)then   ! i. e. case of psno=0. 280113
        if(psno.eq.0)then   ! 260718
        bp=bmin
        r4=rnp
        if(rnt.gt.rnp)r4=rnt
        rr4=bp/r4
        vneu=dexp(-rr4*rr4)
c       calculates the overlap region of two nuclei at given b by 
c020511  interpolation method
        if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.   ! 180921 yan
     c   itden.eq.2) .or. ipden.ge.11)goto 80002   ! 020718, 180921 yan
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
c291118
c180219 if(ipden.eq.1 .and. itden.eq.1)then   ! A+B
        anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
c180219 elseif(ipden.eq.0 .and. itden.eq.1)then   ! p+A
c       anbin=ta2(ibpp)   ! overlap function of B nucleus (1/fm^2)
c       elseif(ipden.eq.1 .and. itden.eq.0)then   ! A+p
c       anbin=ta1(ibpp)   ! overlap function of A nucleus (1/fm^2) 020718
c       else
c180219 endif
c291118
        pir=part1(ibpp)
        tir=part2(ibpp)
        evbin=anbin*csnn   !  Nbin in current event 020718
        pirr=pir   ! 140219
        tirr=tir   ! 180219
c       write(9,*)'bp,ibpp,part1,part2,evbin=',bp,ibpp,pir,tir,evbin   ! 020718
c280113 endif
c020511

80002   if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.   ! 180921 yan
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
        write(9,*)
        write(9,*)'psno,b,N_part_p,N_part_t,N_bin=',
     c   psno,bp,vneump,vneumt,evbin   ! 190309 280113
        goto 80003   ! 010518 changed from 300 to 80003
        endif

        if(psno.eq.1.)then   ! 280113
c       systematic sampling method for given interval of b according to b**2 law
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
c       calculate the overlap region of two nuclei at given b 
        ibpp=int(bp/0.1+1.0)
        ibpp=min(ibpp,200)
c291118
c180219 if(ipden.eq.1 .and. itden.eq.1)then   ! A+B
        anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
c180219 elseif(ipden.eq.0 .and. itden.eq.1)then   ! p+A
c       anbin=ta2(ibpp)   ! overlap function of B nucleus (1/fm^2)
c       elseif(ipden.eq.1 .and. itden.eq.0)then   ! A+p
c       anbin=ta1(ibpp)   ! overlap function of A nucleus (1/fm^2) 020718
c       else
c180219 endif
c291118
        pir=part1(ibpp)
        tir=part2(ibpp)
c060605
        stab=stab+bp   ! 280113
        acoll(i1)=anbin   ! 280113
        acollp(i1)=pir
        acollt(i1)=tir
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
     c   aneump,aneumt,vneum*csnn ! 191110      280113

c       average b in [bmin,bmax]
        avb=2./3.*(bmin+bmax)
c       above equation is correct when bmin=0 only
        r4=rnp
        if(rnt.gt.rnp)r4=rnt
        rr4=avb/r4
        avneu=dexp(-rr4*rr4)
c       calculate the overlap region of two nuclei at given b
        ibpp=int(avb/0.1+1.0)
        ibpp=min(ibpp,200)
c291118
c180219 if(ipden.eq.1 .and. itden.eq.1)then   ! A+B
        anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
c180219 elseif(ipden.eq.0 .and. itden.eq.1)then   ! p+A
c       anbin=ta2(ibpp)   ! overlap function of B nucleus (1/fm^2)
c       elseif(ipden.eq.1 .and. itden.eq.0)then   ! A+p
c       anbin=ta1(ibpp)   ! overlap function of A nucleus (1/fm^2) 020718
c       else
c180219 endif
c291118
        pir=part1(ibpp)
        tir=part2(ibpp)
        aanbin=anbin   ! 280113
        astbp=pir
        astbt=tir
c220110
        endif   ! 280113
c-------------   Independent Optical Initial Geometry Calculation  -------------
c-------------------------------------------------------------------------------


80003   continue   ! 310518


c-------------------------------------------------------------------------------
c----------------------------   Input Recording   ------------------------------
        write(9,*)'win,nap,nzp,nat,nzt=',win,nap,nzp,nat,nzt
        write(9,*)'pio,ipden,itden=',pio,ipden,itden   !Lei2023060
        write(9,*)'neve,nout,nosc=',neve,nout,nosc
        write(9,*)'x_ratio,bmin,bmax,nmax,parp78,mstptj=',x_ratio,bmin,
     &   bmax,nmax,parp78,mstptj   ! 150612 yan 070417 100821 230722 Lei202307
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
     c   para4
        write(9,*)'i_deex,i_deex_gen,i_pT,i_pt_max,i_split_diq,'//
     c   'i_split_qqb,i_split_g,a_FF,aPS_c,aPS_b,parj23,parj24='   !Lei2023060
        write(9,*) i_deex,i_deex_gen,i_pT,i_pt_max,i_split_diq,
     c   i_split_qqb,i_split_g,a_FF,aPS_c,aPS_b,parj23,parj24    !Lei2023060
        write(9,*)'mstu21,mstj1_1,mstj1_2,mstj2,decpro,itorw=',
     c   mstu21,mstj1_1,mstj1_2,mstj2,decpro,itorw   ! 160617 010418 Lei202307
c210803
        write(9,*)'adj1=',(adj1(i),i=1,10)
        write(9,*)'adj1=',(adj1(i),i=11,20)
        write(9,*)'adj1=',(adj1(i),i=21,30)
        write(9,*)'adj1=',(adj1(i),i=31,40)
c210803
        write(9,*)'parecc,iparres,smadel,dparj4,cp0,cr0,seco=',
     c   parecc,iparres,smadel,dparj4,cp0,cr0,seco   ! 120219 260219
        if(iflmax.ne.0)then
        do kk=1,ispmax
        do i=1,iflmax
        write(9,*)(afl(kk,i,j),j=1,2)
        enddo
        enddo
        endif
c----------------------------   Input Recording   ------------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c------------------------   Random Number Seed Giving   ------------------------
c171022 Note it is an array. One doesn't need to initialize it using do-loop.
        n_current_date_and_time = 0   ! 171022 Lei
        i_seed = MRPY(1)

        iran = INT(adj1(26))
        if(iran.eq.0)goto 300
c171022 do i1=1,iran
c171022 thrr=pyr(1)
c171022 enddo

c171022 A new method is introduced now.   ! 171022 Lei
c       The initial seed of the random number generator will be given by 
c        real-time clock from computer.
c       Get the current real-date_and_time from the computer.
        call DATE_AND_TIME(VALUES=n_current_date_and_time)   ! Fortran intrinsic
C       Only use hour-to-milliseconds because of the limitation of INTEGER type.
c                n_current_date_and_time(1)    ! year   (CCYY)
c                n_current_date_and_time(2)    ! month  (1-12)
c                n_current_date_and_time(3)    ! day    (1-31)
c                n_current_date_and_time(4)    ! The time difference, in minutes, 
c                                              !  with respect to UTC.
        i_seed = n_current_date_and_time(5) * 10000000 +   ! hour   (1-23)
     &           n_current_date_and_time(6) * 100000   +   ! minute (1-59)
     &           n_current_date_and_time(7) * 1000     +   ! second (0-60)
     &           n_current_date_and_time(8) * 1            ! milliseconds (0-999)
c       The seed will be 8-digit HH-MM-SS-ms-ms.
        MRPY(1) = i_seed
c------------------------   Random Number Seed Giving   ------------------------
c-------------------------------------------------------------------------------




c*******************************************************************************
c-------------------------------------------------------------------------------
c----------------------------   Event Generating   -----------------------------
c       loop over event
300     iii=iii+1
c       print*,'loop over event, iii=',iii   ! sa
c220110
        do i1=1,9
        nreaco(i1)=nreac(i1)
        enddo
c220110


c-------------------------------------------------------------------------------
c-------------------   Single Event Variable Initializing   --------------------
c250209
        ngam=0
        nsin=0
        ntrs=0   !Lei2023060
        nth=0    !Lei2023060
        naff=0   !Lei2023060
        if(iii.eq.1)then   !Lei2023060
        kgam=0
        pgam=0.
        vgam=0.
        ksin=0
        psin=0.
        vsin=0.
c280620
        ktrs=0
        ptrs=0.
        vtrs=0.
c280620
        kaff=0   !Lei2023060
        paff=0.
        vaff=0.
        endif   !Lei2023060
        siijk=0   ! 201203
c061103
        noel=0
        noinel=0
        npinel=0   ! 140820
c061103
        ppsa=0.
        throe_p=0.   !Lei2023060 For global 4-momentum adjustment.
c250209

        ijk=0   ! 10/08/98
c181003
        nnstop=0
        zstop=0.
c181003
c-------------------   Single Event Variable Initializing   --------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c------------------------   Impact Parameter Sampling   ------------------------
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
c180219 if(ipden.eq.1 .and. itden.eq.1)then   ! A+B
        anbin=ta1a2(ibpp)   ! overlap function of A+B (1/fm^2) 280113
c180219 elseif(ipden.eq.0 .and. itden.eq.1)then   ! p+A
c       anbin=ta2(ibpp)   ! overlap function of B nucleus (1/fm^2)
c       elseif(ipden.eq.1 .and. itden.eq.0)then   ! A+p
c       anbin=ta1(ibpp)   ! overlap function of A nucleus (1/fm^2) 020718
c       else
c       endif
c291118
        pir=part1(ibpp)
        tir=part2(ibpp)
        evbin=anbin*csnn   ! 020718
        pirr=pir   ! 140219
        tirr=tir   ! 180219
c280113 endif
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
c------------------------   Impact Parameter Sampling   ------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c------------------------   Particle Decay Specifying   ------------------------
c       forbiden decay of particle (sets mdcy(...)=0)
c       mdcy(pycomp(111),1)=0    ! pi0
        mdcy(pycomp(310),1)=0    ! K0_S
        mdcy(pycomp(333),1)=0    ! phi
        mdcy(pycomp(3122),1)=0   ! Lambda0
        mdcy(pycomp(-3122),1)=0  ! Lambdabar0
        mdcy(pycomp(443),1)=0    ! J/psi
c       mdcy(pycomp(10441),1)=0  ! chi_0c
c       mdcy(pycomp(20443),1)=0  ! chi_1c
c       mdcy(pycomp(445),1)=0    ! chi_2c
c       mdcy(pycomp(411),1)=0    ! Sigma_c0
c       mdcy(pycomp(-411),1)=0   ! Sigma_cbar0
c       mdcy(pycomp(421),1)=0    ! Sigma_c+
c       mdcy(pycomp(-421),1)=0   ! Sigma_cbar-
c       mdcy(pycomp(4122),1)=0   ! Lambda_c+
c       mdcy(pycomp(4112),1)=0   ! Sigma_c0Sigma_c0
c       mdcy(pycomp(4212),1)=0   ! Sigma_c+
c       mdcy(pycomp(4222),1)=0   ! Sigma_cbar-
c       mdcy(pycomp(3212),1)=0   ! Sigma0
c       mdcy(pycomp(-3212),1)=0  ! Sigmabar0
        mdcy(pycomp(3112),1)=0   ! Sigma-
c       mdcy(pycomp(-3112),1)=0  ! Sigmabar+
        mdcy(pycomp(3222),1)=0   ! Sigma+
c       mdcy(pycomp(-3222),1)=0  ! Sigmabar-
        mdcy(pycomp(3312),1)=0   ! Xi-
        mdcy(pycomp(-3312),1)=0  ! Xibar+
c       mdcy(pycomp(3322),1)=0   ! Xi0
c       mdcy(pycomp(-3322),1)=0  ! Xibar0
        mdcy(pycomp(3334),1)=0   ! Omega-
        mdcy(pycomp(-3334),1)=0  ! Omegabar+
c       mdcy(pycomp(1114),1)=0   ! Delta-
c       mdcy(pycomp(2114),1)=0   ! Delta0
c       mdcy(pycomp(2214),1)=0   ! Delta+
c       mdcy(pycomp(2224),1)=0   ! Delta++
c       mdcy(pycomp(213),1)=0    ! rho+
c       mdcy(pycomp(-213),1)=0   ! rho-
c       mdcy(pycomp(113),1)=0    ! rho0
c       mdcy(pycomp(223),1)=0    ! omega
c       mdcy(pycomp(413),1)=0    ! D*+
c       mdcy(pycomp(-413),1)=0   ! D*-
c       mdcy(pycomp(423),1)=0    ! D*0
c       mdcy(pycomp(-423),1)=0   ! D*bar0
c       mdcy(pycomp(13),1)=0     ! mu-
c       mdcy(pycomp(-13),1)=0    ! mu+
c00623 Lei2023060B
c       MDCY( KC , 1 ) = 1 means the particle will decay (unstable)
        if(mstj2.eq.1)then   ! For A-loop, hardcode
            MDCY( PYCOMP(111),  1 ) = 1   ! pi0
            MDCY( PYCOMP(1114), 1 ) = 1   ! Delta-
            MDCY( PYCOMP(2114), 1 ) = 1   ! Delta0
            MDCY( PYCOMP(2214), 1 ) = 1   ! Delta+
            MDCY( PYCOMP(2224), 1 ) = 1   ! Delta++
        end if
c00623 Lei2023060E
c------------------------   Particle Decay Specifying   ------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c-------------------------   Subprocesses Selecting   --------------------------
c060620
c       default pythia or user own selection for subprocesses
        if(nchan.eq.0)then
c       Inelastic (INEL)
        msel=0
        msub(11)=1   ! Hard QCD, f_i + f_j -> f_i + f_j
        msub(12)=1   ! Hard QCD, f_i + f_i^bar -> f_k + f_k^bar
        msub(13)=1   ! Hard QCD, f_i + f_i^bar -> g + g
        msub(28)=1   ! Hard QCD, f_i + g -> f_i + g
        msub(53)=1   ! Hard QCD, g + g -> f_k + f_k^bar
        msub(68)=1   ! Hard QCD, g + g -> g + g
c       msub(91)=1   ! Soft QCD, elastic scattering
        msub(92)=1   ! Soft QCD, single diffraction (AB -> XB)
        msub(93)=1   ! Soft QCD, single diffraction (AB -> AX)
        msub(94)=1   ! Soft QCD, double diffraction
        msub(95)=1   ! Soft QCD, low_pT production
        endif
        if(nchan.eq.1)then
c       Non Single Difractive (NSD)
        msel=0
        msub(11)=1   ! Hard QCD, f_i + f_j -> f_i + f_j
        msub(12)=1   ! Hard QCD, f_i + f_i^bar -> f_k + f_k^bar
        msub(13)=1   ! Hard QCD, f_i + f_i^bar -> g + g
        msub(28)=1   ! Hard QCD, f_i + g -> f_i + g
        msub(53)=1   ! Hard QCD, g + g -> f_k + f_k^bar
        msub(68)=1   ! Hard QCD, g + g -> g + g
c       msub(91)=1   ! Soft QCD, elastic scattering
c       msub(92)=1   ! Soft QCD, single diffraction (AB -> XB)
c       msub(93)=1   ! Soft QCD, single diffraction (AB -> AX)
        msub(94)=1   ! Soft QCD, double diffraction
        msub(95)=1   ! Soft QCD, low_pT production
        endif
c090921
        if(nchan.eq.3)then
c       J/psi production
        MSEL=61
c       Octec radiation.
        ! MSTP(148)=1
        ! MSTP(149)=1
        endif
c00623   Lei2023060B---
        if(nchan.eq.4)then
c       Heavy-flavor production
        msel=0
        msub(81)=1   ! Open HF, f_i + f_i^bar -> Q_k + Q_k^bar
        msub(82)=1   ! Open HF, g + g -> Q_k + Q_k^bar
        ! msub(83)=1 ! Open HF, q_i + f_i -> Q_k + f_l
        endif
        if(nchan.eq.5)then
c       Direct photon
        msel=0
        msub(14)=1    ! Prompt photon, f_i + f_i^bar -> g + gamma
        msub(18)=1    ! Prompt photon, f_i + f_i^bar -> gamma + gamma
        msub(29)=1    ! Prompt photon, f_i + g -> f_i + gamma
        msub(114)=1   ! Prompt photon, g + g -> gamma + gamma
        msub(115)=1   ! Prompt photon, g + g -> g + gamma
        endif
        if(nchan.eq.6)then
c       Soft QCD (Minimum Bias)
        msel=0
        msub(91)=1   ! Soft QCD, elastic scattering
        msub(92)=1   ! Soft QCD, single diffraction (AB -> XB)
        msub(93)=1   ! Soft QCD, single diffraction (AB -> AX)
        msub(94)=1   ! Soft QCD, double diffraction
        msub(95)=1   ! Soft QCD, low_pT production
        endif
c00623   Lei2023060E---
c090921
        if(nchan.eq.7)then
c       W+/- production (nchan=7,isub=2,16,20,23,25,31)
        msel=0
        msub(2)=1      ! Single W, f_i + f_j^bar -> W
        ! msub(16)=1   ! Single W, f_i + f_j^bar -> W + g
        ! msub(20)=1   ! Single W, f_i + f_j^bar -> W + gamma
        ! msub(23)=1   ! W/Z pair, f_i + f_j^bar -> W + Z0
        ! msub(25)=1   ! W pair,   f_i + f_i^bar -> W+ + W-
        ! msub(31)=1   ! Single W, f_i + g -> W + f_k
        mdcy(pycomp(24),1)=0    ! W+
        mdcy(pycomp(-24),1)=0   ! W-
        endif
        if(nchan.eq.8) MSEL=1   !Lei2023060 PYTHIA default
        if(nchan.eq.9)then
c       Z0 production (nchan=9,isub=1,15,19,22,23,30)
        msel=0
        msub(1)=1      ! Single Z, f_i + f_i^bar -> gamma*/Z0
        ! msub(15)=1   ! Single Z, f_i + f_i^bar -> gamma*/Z0 + g
        ! msub(19)=1   ! Single Z, f_i + f_i^bar -> gamma*/Z0 + gamma
        ! msub(22)=1   ! Single Z, f_i + f_j^bar -> W
        ! msub(23)=1   ! W/Z pair, f_i + f_j^bar -> W + Z0
        ! msub(30)=1   ! Single Z, f_i + g -> gamma*/Z0 + f_k
        mdcy(pycomp(23),1)=0   ! Z0
        endif
c060620
c-------------------------   Subprocesses Selecting   --------------------------
c-------------------------------------------------------------------------------


c       creats pp (pA,Ap,AB,lp, lA, & e+e-) collision events


c-------------------------------------------------------------------------------
c--------------------   Independent NN & e+e- Generating   ---------------------
c210921
        if((ipden.eq.0 .and. itden.eq.0) .or.
     c   (ipden.eq.2 .and. itden.eq.2))then   ! if
c131019 writing initialization and differential cross section maximum
        if(iii.eq.1)then
        mstp(122)=1
        else
        mstp(122)=0
        endif
c131019
c       for nucleon (antinucleon)-nucleon (antinucleon) or e+e-
        n=0
        nbe=0
        naf=0
        nsa=0
        idi=0
        idio=0
        if(iii.eq.1)then   !Lei2023060
        k=0
        p=0.
        v=0.
        kbe=0
        pbe=0.
        vbe=0.
        kaf=0
        paf=0.
        vaf=0.
        ksa=0
        psa=0.
        vsa=0.
        endif   !Lei2023060
        ndiq=0
        npt=0
        ifcom=0
        ishp=0
        tau=0.

c--------------------------   NN Low Energy A-loop   ---------------------------
c161222
        if((ipden.eq.0.and.itden.eq.0).and.mstj2.eq.1)then  !!! 161222
c260223 changed from 'win.lt.parp21' to mstj2.eq.1, for low energy A-loop

c       initiation for a nucleon-nucleon collision in 'pyjets'
        n=2
c       initiation array 'k'
        k(1,1)=2   ! A
        k(2,1)=1   ! V

        if(nzp.eq.1)then
        k(1,2)=2212
        elseif(nzp.eq.-1)then
        k(1,2)=-2212
        elseif(nzp.eq.0)then
        k(1,2)=2112
        endif

        if(nzt.eq.1)then
        k(2,2)=2212
        elseif(nzt.eq.-1)then
        k(2,2)=-2212
        elseif(nzt.eq.0)then
        k(2,2)=2112
        endif

c       initiation in spatial phase space
c       v=0.

c       initiation in momentum phase space

        if(ifram.eq.1)then
        ep1=0.5*win   ! energy of projetile particle (if it is proton)
        et1=ep1   ! energy of target particle (if proton)
        ep2=0.5*win   ! energy of projetile particle (if neutron)
        et2=ep2   ! energy of target particle (if neutron)
        pm2=pymass(2212)**2   ! square mass of proton   
        pp1=dsqrt(ep1*ep1-pm2)   ! momentum of projetile particle (if proton)
        pt1=-dsqrt(et1*et1-pm2)  ! momentum of target particle (if proton)
        pm2=pymass(2112)**2   ! square mass of nucleon  
        pp2=dsqrt(ep2*ep2-pm2)   ! momentum of projetile particle (if neutron)
        pt2=-dsqrt(et2*et2-pm2)  ! momentum of target particle (if neutron)
        endif

        if(ifram.eq.0)then
        pp1=win   ! momentum of projetile particle (if proton)
        pt1=1.e-20   ! momentum of target particle (if proton)
        pp2=win   ! momentum of projetile particle (if neutron)
        pt2=1.e-20   ! momentum of target particle (if neutron)
        pm2=pymass(2212)**2   ! square mass of proton
        ep1=dsqrt(pp1*pp1+pm2)   ! energy of projetile particle (if proton)
        et1=dsqrt(pt1*pt1+pm2)   ! energy of target particle (if proton)
        pm2=pymass(2112)**2   ! square mass of neutron
        ep2=dsqrt(pp2*pp2+pm2)   ! energy of projetile particle (if neutron)
        et2=dsqrt(pt2*pt2+pm2)   ! energy of target particle (if neutron)
        endif

        if(ifram.eq.1)then   !!!

c       four momenta of projectile particle
        p(1,1)=0.
        p(1,2)=0.
        if(iabs(nzp).eq.1)then
        p(1,3)=pp1   ! projectile particle is p (pba)
        p(1,4)=ep1
        p(1,5)=pymass(2212)
        elseif(nzp.eq.0)then
        p(1,3)=pp2   ! projectile particle is neutron
        p(1,4)=ep2
        p(1,5)=pymass(2112)
        endif

c       four momenta of target particle
        p(2,1)=0.
        p(2,2)=0.
        if(iabs(nzt).eq.1)then
        p(2,3)=pt1   ! target particle is p (pba)
        p(2,4)=et1
        p(2,5)=pymass(2212)
        elseif(nzt.eq.0)then
        p(2,3)=pt2   ! target particle is neutron
        p(2,4)=et2
        p(2,5)=pymass(2112)
        endif

        endif   !!!


        if(ifram.eq.0)then   !!!!

c       four momenta of projectile particle
        p(1,1)=0.
        p(1,2)=0.
        if(iabs(nzp).eq.1)then
        p(1,3)=pp1   ! projectile particle is p (pba)
        p(1,4)=ep1
        p(1,5)=pymass(2212)
        elseif(nzp.eq.0)then
        p(1,3)=pp2   ! projectile particle is neutron
        p(1,4)=ep2
        p(1,5)=pymass(2112)
        endif

c       four momenta of target particle
        p(2,1)=0.
        p(2,2)=0.
        if(iabs(nzt).eq.1)then
        p(2,3)=pt1   ! target particle is p (pba) 
        p(2,4)=et1
        p(2,5)=pymass(2212)
        elseif(nzt.eq.0)then
        p(2,3)=pt2   ! target particle is neutron
        p(2,4)=et2
        p(2,5)=pymass(2112)
        endif

        endif   !!!!

c       initiation finished

c       'pyjets' to 'sa2'
        nsa=n
        do i2=1,5
        do i1=1,n
        ksa(i1,i2)=k(i1,i2)
        psa(i1,i2)=p(i1,i2)
        vsa(i1,i2)=v(i1,i2)
        enddo
        enddo

        do i2=1,4
        pi(i2)=psa(1,i2)   ! four momentum of one colliding particle
        pj(i2)=psa(2,i2)   ! four momentum of another one colliding particle
        enddo

        if(ifram.eq.0)then
c       boost to CMS frame of colliding pair
        do i2=1,3
        b(i2)=(pi(i2)+pj(i2))/(pi(4)+pj(4))
        enddo
        ilo=0
        call lorntz(ilo,b,pi,pj)   ! in parini.f
        endif
c260223
        if(ifram.eq.1)b=0.    ! ifram=1, b(i)=0 & \gamma=1
c260223

c150323 following is for case of NN collision only
        ww=rcsit
c       the cross section ratio of (ela.)/tot =1- rcsit (rcsit=inela./tot)
        rrlu=pyr(1)

        if(rrlu.gt.ww)then   ! ela.  !! 161222
c       perform nucleon-nucleon (NN) elastic scattering
        call coelas_nn(1,2,pi,pj,lc,tc,tw,time,b,2)   ! in parini_23.f
        inorex=1
        goto 998   ! toward hadron rescattering (call 'hadcas')
        elseif(rrlu.le.ww)then   ! inela.   !! 161222
c       perform NN inelastic scattering
c     p: 2212; n:2112; delta0: 2114; delta+: 2214; delta-: 1114; delta++: 2224
c       consider following 2->2 NN ielastic channels
c     1       p + p to delta+ + p
c     2       p + p to delta++ + n
c     3       p + n to delta+ + n
c     4       p + n to delta0 + p
c     5       n + n to delta0 + n
c     6       n + n to delta- + p
c       reverse scattering is not considered

c260223 inorex: endothermic or exothermic
c       endothermic reaction, inorex=1 & ela. (cf. 'coinel_nn')
c       exothermic reaction, inorex=2 & inela. (cf. 'coinel_nn')

c       give flavor code to inelastically scattered particles
c150323 inorex=1   ! 260223
        rpy=pyr(1)
        if(nzp.eq.1 .and. nzt .eq.1)then   ! pp   !!
        if(rpy .gt. 0.5)then
        ksa(1,2)=2214
        ksa(2,2)=2212
        call adjudg_nn(1,2,2214,2212,pi,pj,inorex)   ! in parini.f
        if(inorex.eq.1)then
        call coelas_nn(1,2,pi,pj,lc,tc,tw,time,b,2)   ! 120323
        goto 998
        endif
        if(inorex.eq.2)then
c       give four momentum to inelastically scattered particles 
        call coinel_nn(1,2,2214,2212,pi,pj)   ! in parini.f
c150323 call padecy(1,2)   ! in parini.f
        goto 998
        endif
        else
        ksa(1,2)=2224
        ksa(2,2)=2112
        call adjudg_nn(1,2,2224,2112,pi,pj,inorex)
        if(inorex.eq.1)then
        call coelas_nn(1,2,pi,pj,lc,tc,tw,time,b,2)   ! 120323
        goto 998
        endif
        if(inorex.eq.2)then
        call coinel_nn(1,2,2224,2112,pi,pj)
c150323 call padecy(1,2)
        goto 998
        endif
        endif

        elseif((nzp.eq.1 .and. nzt .eq.0).or.   ! pn   !!
     c   (nzt.eq.1 .and. nzp .eq.0))then
        if(rpy .gt. 0.5)then
        ksa(1,2)=2214
        ksa(2,2)=2112
        call adjudg_nn(1,2,2214,2112,pi,pj,inorex)
        if(inorex.eq.1)then
        call coelas_nn(1,2,pi,pj,lc,tc,tw,time,b,2)   ! 120323
        goto 998
        endif
        if(inorex.eq.2)then
        call coinel_nn(1,2,2214,2112,pi,pj)
c150323 call padecy(1,2)
        goto 998
        endif
        else
        ksa(1,2)=2114
        ksa(2,2)=2212
        call adjudg_nn(1,2,2114,2212,pi,pj,inorex)
        if(inorex.eq.1)then
        call coelas_nn(1,2,pi,pj,lc,tc,tw,time,b,2)   ! 120323
        goto 998
        endif
        if(inorex.eq.2)then
        call coinel_nn(1,2,2114,2212,pi,pj)
c150323 call padecy(1,2)
        goto 998
        endif
        endif

        elseif(nzp.eq.0 .and. nzt .eq.0)then   ! nn   !!
        if(rpy .gt. 0.5)then
        ksa(1,2)=2114
        ksa(2,2)=2112
        call adjudg_nn(1,2,2114,2112,pi,pj,inorex)
        if(inorex.eq.1)then
        call coelas_nn(1,2,pi,pj,lc,tc,tw,time,b,2)   ! 120323
        goto 998
        endif
        if(inorex.eq.2)then
        call coinel_nn(1,2,2114,2112,pi,pj)
c150323 call padecy(1,2)
        goto 998
        endif
        else
        ksa(1,2)=1114
        ksa(2,2)=2212
        call adjudg_nn(1,2,1114,2212,pi,pj,inorex)
        if(inorex.eq.1)then
        call coelas_nn(1,2,pi,pj,lc,tc,tw,time,b,2)   ! 120323
        goto 998
        endif
        if(inorex.eq.2)then
        call coinel_nn(1,2,1114,2212,pi,pj)
c150323 call padecy(1,2)
        goto 998
        endif
        endif
        endif   !!

        endif   !! 161222

c       low energy NN collision A-loop finished

        endif   !!! 161222
c161222
c--------------------------   NN Low Energy A-loop   ---------------------------

c------------------------   NN B- & C-loop Generating  -------------------------
        if(ipden.eq.0 .and. itden.eq.0)then

99999   continue   !Lei2023060
        mstp(111)=mstptj   ! 151021 230722
        MSTP(5)=i_tune   !Lei202306

        if(ifram.eq.1)then
        if(nzp.eq.1 .and. nzt.eq.1)call pyinit('cms','p','p',win)   ! in p_23.f
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
c00623 Lei2023060B
c       Sums of incident px, py, pz, E, inv. m, and charge.
        ps0=0.
        do i=1,6,1
            ps0(i)=PYP(0,i)
        end do
!20523 Lei2023060E
c151021 endif
        if(itden.eq.0 .and. itorw.eq.1)call pyevnt   ! in p_23.f
        if(itden.eq.0 .and. itorw.eq.2)call pyevnw   ! in p_23.f
c151021
c00623 Lei2023060B
c       Sums of px, py, pz, E, inv. m, and charge after the excution.
        ps1=0.
        do i=1,6,1
            ps1(i)=PYP(0,i)
        end do
c       Charge is not conserved. Re-generate the event.
        if( ABS(ps0(6)-ps1(6)).gt.1D-10 ) goto 99999   ! Charge.
c       4-momentum is not conserved. Re-generate the event.
        do i=1,4,1   ! px, py, pz, E
            if( ABS(ps0(i)-ps1(i)).gt.1D-10 ) goto 99999
        end do
c       Error count in PYTHIA. Re-generate the event if any errors exit.
        if( MSTU(23).gt.0 .OR. MSTU(30).gt.0 ) goto 99999
        call PASTAT(-1,1)  !Lei202306 Counts single-event cross sections.
c00623 Lei2023060E

        endif
c------------------------   NN B- & C-loop Generating  -------------------------

c-----------------------------   e+e- Generating   -----------------------------
        if(ipden.eq.2 .and. itden.eq.2)then
        mstj(1)=mstptj   ! 230722
        call pyeevt(0,win)   ! for e+e-,   ! in p_23.f
        endif
c-----------------------------   e+e- Generating   -----------------------------

c151021
        call pyedit(2)   ! in p_23.f
c       call pylist(1)   ! in p_23.f

c-----------------------------   B-loop Treating  ------------------------------
c230722
        if(mstptj.eq.1)then   !! 230722 PYTHIA-like simulation for pp & e+e-
c       give four position to the particles generated in pythia ('pyjets')
        call ptcre(1,2,time)   ! in parini.f
        goto 998   ! toward hadron rescattering ('call hadcas') for pp & e+e-
        else   !! 230722
c230722
c-----------------------------   B-loop Treating  ------------------------------

c*****************************   C-loop Treating  ******************************
c----------------------------   Gamma 66 Removing   ----------------------------
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
c       remove hadrons from 'pyjets' to 'sbh'
        call remo   ! in parini.f
        call pyedit(2)
c----------------------------   Gamma 66 Removing   ----------------------------

c----------------------------   Diquark Locating   -----------------------------
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
c----------------------------   Diquark Locating   -----------------------------

c---------------------------   Diquark Breaking-up   ---------------------------
c       break up diquark and give four momentum and four position
c        to broken quarks (working in 'pyjets')
        call break   ! in parini.f
        call pyedit(2)
c---------------------------   Diquark Breaking-up   ---------------------------

c140223 Lei full_events_history of OSC1999A
        call oscar(win,1)

c-----------------------------   String Locating   -----------------------------
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
c-----------------------------   String Locating   -----------------------------

        goto 999   ! toward parton rescattering ('call parcas')

        endif   !! 230722
c*****************************   C-loop Treating  ******************************

        endif   ! if
c210921
c--------------------   Independent NN & e+e- Generating   ---------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c---------------------------   Partonic Initiation   ---------------------------
        time_ini=0.d0   ! 081010
c120620 mstp(111)=0   ! 050620 
        npctlm=0   ! 180121
c       partonic initiation for a nuclear-nuclear collision   ! 161222
        call parini(time_ini,parp21,parp22,win,psno,ijk,mstj2,decpro,
     &   i_tune)   ! in parini.f Lei2023061 added i_tune
c       081010 240513 260223
c260223 added the mstj2 & decpro
c120620  mstp(111)=1   ! 050620
        if(ijk.eq.1)goto 300   ! to avoide infinite loop in parcas 060813
        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813
c---------------------------   Partonic Initiation   ---------------------------
c-------------------------------------------------------------------------------


c140223 Lei full_events_history of OSC1999A
        call oscar(win,1)


c-------------------------------------------------------------------------------
c------------------------   Wounded Nucleons Counting   ------------------------
        if(mstj2.eq.2.or.mstj2.eq.3)then   ! 150323
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
c       endif   ! 1
c230722
c       if(((ipden.eq.1.and.itden.eq.1).or.(ipden.eq.0.and.itden.eq.1)
c       c   .or.(ipden.eq.1.and.itden.eq.0)).and.mstptj.eq.1)then
        endif   ! 150323
c------------------------   Wounded Nucleons Counting   ------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c-----------------------------   B-loop Treating  ------------------------------
        if(mstptj.eq.1)then   !! 100223
c121222 PYTHIA like simulation for pA (Ap) & AB, and low energy simulation
c       'sbh' to 'pyjets'
        n=nbh
        if(n.ge.1)then
        do i2=1,5
        do i1=1,n
        k(i1,i2)=kbh(i1,i2)
        p(i1,i2)=pbh(i1,i2)
        v(i1,i2)=vbh(i1,i2)
        enddo
        enddo
        endif
        goto 998  !toward hadron rescattering ('call hadcas'), pA (Ap) & AB
        endif   !!

c240123 if(win.lt.parp21)goto 998!toward hadron rescattering,pA(Ap),&,AB,161222
c-----------------------------   B-loop Treating  ------------------------------
c-------------------------------------------------------------------------------



c*******************************************************************************
c-----------------------------   C-loop Treating  ------------------------------
c230722
999     continue   ! 230722
        if(mstptj.eq.0)then   ! 230722


c-------------------------------------------------------------------------------
c-----------------------   Diffractive Event Treating  -------------------------
c       no parton produced at all
        if(n.le.0)then
        nncoll=nncoll+1
        if(nbh.gt.0) goto 333   !Lei2023060 Diffractive event in PYTHIA. Do not throw it away.
c060814 if(nncoll.gt.neve)then
c       stop 8888
c060814 endif
        iii=iii-1
        goto 300
        endif
c-----------------------   Diffractive Event Treating  -------------------------
c-------------------------------------------------------------------------------


c       throw away event with junction if iparres=1
c040223 if(iparres.eq.1)then
c       do i1=1,n
c       kf=k(i1,2)
c       if(kf.eq.88)then
c       iii=iii-1
c       goto 300
c       endif
c       enddo
c       endif
c180520


c-------------------------------------------------------------------------------
c------------------------   Simulation Stopping 1 & 3  -------------------------
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
        if(adj140.eq.3)call coales(iii,neve,nout,nap,nat,nzp,nzt)   ! in coales.f

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
5000    continue   ! 300407
c       if(iii.eq.3)call pylist(1)
        goto 888
        endif   ! 290505 271205 020718
c271205
c       goto 889   ! temporal
c------------------------   Simulation Stopping 1 & 3  -------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c----------------------------   Gamma 44 Removing   ----------------------------
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
c----------------------------   Gamma 44 Removing   ----------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c---------------------------   Coalescence Mode 2   ----------------------------
!Lei2023060B---
c Special version of coalescence, gluon splitting and energetic quark 
c  deexcitation before 'parcas'. It means that there will be no gluons 
c  into 'parcas', similar to the ZPC of AMPT string-melting in some sence.
        if( INT(adj1(12)).eq.2 )then
c           Moves gluons from "PYJEST" to "sa36".
            call remo_glu   ! in coales.f
c           Breaks gluons up (with E_g > 2E_u in "sa36") -> qqbar string 
c            (filling in "PYJETS")
            call break_glu   ! in coales.f
c           Shares 4-momentum in "throe_p" among partons.
            call share_p_PYJETS
c           Energetic q (qbar) de-excitation.
            n00   = n   ! Original total entries in PYJETS
            igens = 0
            i_daught_gen = 0   ! the #-th newly produced daughter qqbar
            n_deex = 0   ! number of successful deexcitation
            jb = 0
            n0 = n   ! Current total entries in PYJETS
700         continue
            do i1=jb+1,n0,1
                kf0   = k(i1,2)
                ee    = p(i1,4)
                iflav =1
                if(kf0.lt.0) iflav = -1
c               iflav = 1 : if source parton is quark
c                     =-1 : if source parton is antiquark
                if(ee.gt.adj1(17))then
                if(i_deex.eq.1) call deexcitation_EP(i1,kf0,igen,iflav)   ! in coales.f
                if(i_deex.eq.2) call deexcitation_E(i1,kf0,igen,iflav)    ! in coales.f
      if(i_deex.eq.3) call deexcitation_EP_comp_pT_1(i1,kf0,igen,iflav)    !Lei2023071
      if(i_deex.eq.4) call deexcitation_EP_comp_pT_2(i1,kf0,igen,iflav)    !Lei2023071
                if(igen.gt.0) n_deex = n_deex + 1
                igens = igens + 1   ! Number of "call deexcitation"
                endif
c               igen : number of generations per source q (qbar)
c               Updates n0 and does deexcitation for newly produced qqbar pair
        if(i1.eq.n0 .AND. n.gt.n0 .AND. i_daught_gen.lt.i_deex_gen)then
c         i_deex_gen=0 means no deexcitation for any newly produced qqbar pairs.
c         i_deex_gen=1 means just do deexcitation for the directly proudced qqbar 
c                       pairs (1-st daughters) from original PYJETS (Orig mothers).
c         i_deex_gen=2 means do deexcitation for "1-st daughters" from "Orig mothers" 
c                       and the subsequent qqbar pairs produced from "1-st daughters".
c                       (2-nd daughters).
c         i_deex_gen=3,4,...
c         ...
c         i_deex_gen=999 means always do deexcitation for newly produced qqbar pair
            jb = i1
            i_daught_gen = i_daught_gen + 1
            n0 = n
            goto 700
        end if
            end do
c           Shares the 4-momentum in 'throe_p' among partons.
            call share_p_PYJETS
        end if
!Lei2023060E---
c---------------------------   Coalescence Mode 2   ----------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c--------------------------   Partonic Rescattering   --------------------------
c       partonic rescattering
        if(n.lt.2)goto 889   ! 151302
c140718 if(itden.ne.1)goto 890   ! for e+e-,p+p,pbar_p, or p+pbar 080806
c201203
c       goto 890
        time_par=0.d0   ! 081010
        iijk=0   ! 151203
        call parcas(time_par,nnn,iijk,win,nap,rnt,rnp)   ! 120603 220110   ! in parcas.f
c220110 nnn: nnn-th parton-parton interacion in a nucleus-nucleus collision
c120603
        if(iijk.eq.1)then
        goto 300   ! give up current event avoiding infinite collision loop
        endif
c120603
        if(iijk.eq.2)siijk=siijk+1   ! 201203
c       if(ipden.lt.11)call pyedit(2)
c       if(ipden.ge.11)call pyedit(1)
c201203
        if(n.eq.0)goto 300   ! no parton at all, give up current event
c--------------------------   Partonic Rescattering   --------------------------
c-------------------------------------------------------------------------------


        endif   ! 230722
c-----------------------------   C-loop Treating  ------------------------------
c*******************************************************************************



c140223 Lei full_events_history of OSC1999A
        call oscar(win,2)


c-------------------------------------------------------------------------------
c--------------------------   Simulation Stopping 2  ---------------------------
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
6000    continue   ! 300407
        goto 888
        endif   ! 290505 271205
c271205
c--------------------------   Simulation Stopping 2  ---------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c----------------------------   Gamma 55 Removing   ----------------------------
c250209
        n55=0
        do i1=1,n
        kf=k(i1,2)
        if(kf.eq.22)then
        k(i1,2)=55
        n55=n55+1
        endif
        enddo
c       move '55' from 'pyjets' to 'sgam'
        if(n55.gt.0)call remo_gam_par(55)
c250209
        call prt_sgam(n55,egam2,2)   ! 080419 160919 270220
c080419 egam2: energy of gammas after parton cascade 
c----------------------------   Gamma 55 Removing   ----------------------------
c-------------------------------------------------------------------------------


890     continue   ! 020512



c*******************************************************************************
c------------------------------   Hadronization  -------------------------------
        if(adj12.ne.0)goto 889   ! coalescence 120520
c120520        if(adj12.ne.0.or.(adj12.eq.0.and.(nreac(4).gt.nreaco(4).or.
c120520     c   nreac(6).gt.nreaco(6).or.nreac(7).gt.nreaco(7))))goto 889 
! 020512 ->coalescence


c-------------------------------------------------------------------------------
c---------------------------   Diquark Recovering   ----------------------------
c       recover parton configuration in 'sbe' (having diquark)
c       loop over 'sbe'   ! 'pyjets' to 'sbe' in parini_23.f
        idii=0
        do i=1,nbe
        kf=kbe(i,2)
        kfab=iabs(kf)
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then   ! 060805
c060805     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then
        idii=idii+1
        do j=1,5
        kdiqq=kbe(i,j)
        kdiq(idii,j)=kdiqq
        enddo
        dgmas(idii)=pbe(i,5)
        endif
        enddo

c       loop over 'pyjets'
        idij=0
        jb=0
c00623 dele=0.   !Lei2023060
880     do 980 i=jb+1,n
        jb=jb+1
        ndiqi=ndiq(i)
        if(ndiqi.ne.0)then   ! diquark (anti-diquark)
        idij=idij+1
        j=npt(ndiqi)
        do i1=1,5
        k(i,i1)=kdiq(idij,i1)
        enddo
        do i1=1,3
        p(i,i1)=p(i,i1)+p(j,i1)
        enddo
        dimass=dgmas(idij)
        pi1=p(i,1)
        pi2=p(i,2)
        pi3=p(i,3)
        pi4=dsqrt(pi1*pi1+pi2*pi2+pi3*pi3+dimass*dimass)
c00623 dele=dele+p(i,4)+p(j,4)-pi4   !Lei2023060
        throe_p(4)=throe_p(4)+p(i,4)+p(j,4)-pi4   !Lei2023060
        p(i,4)=pi4
        p(i,5)=dimass
c00623 Lei2023060
c       Sets 3-coordinate of the diquark as one of the cooresponding quarks.
        ! do i1=1,3
        ! v(i,i1)=v(i,i1)
        ! enddo
        if(PYR(1).gt.0.5)then
        do i1=1,3
        v(i,i1)=v(j,i1)
        enddo
        endif
c00623 Lei2023060
c060805
        if(j.eq.n)then
        n=n-1
        goto 1800
        endif
c060805
        goto 1100
        endif
980     continue
        goto 1800
1100    continue
c       move particle list,'pyjets' ('ndiq') one step downward from 
c       j+1 to n
        do j1=j+1,n
        ndiq(j1-1)=ndiq(j1)
        do jj=1,5
        k(j1-1,jj)=k(j1,jj)
        p(j1-1,jj)=p(j1,jj)
        v(j1-1,jj)=v(j1,jj)
        enddo
        enddo
        n=n-1
c       subtract 'npt' by one from idij+1 to idi
        if(idij.lt.idi)then
        do j1=idij+1,idi
        npt(j1)=npt(j1)-1
        enddo
        endif
        goto 880
1800    continue
c060805 n=jb
c00623 Assigns 4-coordinates of one of the rest partons (or last one?) to the 
c       first parton in the string randomly. This treatment would give random 
c       3-coordinates to produced hadrons that surround the first parton after 
c       PYTHIA sfm, i.e. more random position distribution for produced hadrons. !Lei2023060
        do i1=1,nstr1,1
            i_string_A=nstr1a(i1)
            i_string_V=nstr1v(i1)
            do while(.true.)
                i_A = INT( PYR(1)*(i_string_V-i_string_A+1)+i_string_A )
                if(K(i_A,2).ne.88) exit   ! Excludes junction
            end do
            ! if( PYR(1).ge.0.5 ) i_A = i_string_V
            do j1=1,5,1
                V(i_string_A,j1) = V(i_A,j1)
            end do
        end do
c00623 Recovers the 4-coordinate from 'parini'. (not ?) This treatment is 
c       equivalent to giving medium correction in momentum space.   !Lei2023060
        if(INT(adj12).eq.0 .AND. i_coord_recover.eq.1)then
            do i2=1,5,1
                do i1=1,n,1
                    v(i1,i2)=vbe(i1,i2)
                end do
            end do
        end if
c00623
c---------------------------   Diquark Recovering   ----------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c---------------------------   4-momentum Sharing   ----------------------------
c       share energy 'dele' into particles
c00623 del=dele/dfloat(n)   !Lei2023060
c       do j3=1,n
c       p(j3,4)=p(j3,4)+del
c       if(del.lt.0.)then
c       if(p(j3,4).lt.0.)p(j3,4)=p(j3,4)-del
c       pabs=dabs(p(j3,3))
c       if(pabs.ge.p(j3,4))p(j3,4)=p(j3,4)-del
c       endif
c00623 enddo
c00623 Shares 4-momentum 'throe_p' into partons.   !Lei2023060
        call share_p_PYJETS   !Lei2023060
c---------------------------   4-momentum Sharing   ----------------------------
c-------------------------------------------------------------------------------


889     continue
c010518
        if(ipden.lt.11)call pyedit(2)
        if(ipden.ge.11)call pyedit(1)
c010518


c-------------------------------------------------------------------------------
c-----------------------   String Fragmentation Model  -------------------------
c       hadronization
c230618
        n77s=0   ! 270220
        nbe=0
c230618
        if(adj12.eq.0)then   !! 010518
c120520 if(nreac(4).gt.nreaco(4) .or. nreac(6).gt.nreaco(6)
c120520      c   .or. nreac(7).gt.nreaco(7))then   ! 020512 010518
c020512 for inela. processes 4,6,and 7
c120520 call coales(iii,neve,nout,nap,nat,nzp,nzt)
c120520 else   ! 020512 010518
c       otherwise

c--------------------   String Tension 2 & 3 Calculating  ----------------------
c       calculate the multiple effective string tension and parj(1) etc.
        if(kjp22.eq.2 .or. kjp22.eq.3)then   ! 2
        ampi=mint(31)
c250119 note mint(31)=0 for a low_pT event
c140219
c140219 pathn=(npinel(1)+npinel(592))/pirr
c180219 pathn=(npinel(1)+npinel(592))/(0.5*woun)
c       numerator: NN collision # calculated
c       denominator: N_part of projectile nucleus calculated (pirr: 
c        from Glauber model)
c140219
        if(ampi.le.0)then   !
        ckapa=1.
        elseif(ampi.gt.0)then   !
c120219 ampi=ampi*evbin
c180219
        if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and.   ! 180921 yan
     c   itden.eq.2)) then   ! 180219 
        pathn=1.
        else   !!
        pathn=evbin/(pirr+tirr)   ! 180219
        endif   !!
c       evbin: N_bin of collision system (A+B)
c       pirr: N_part of projectile nucleus (Glauber)
c       tirr: N_part of target nucleus (Glauber)
c180219
        ampi=ampi*pathn   ! 120219
        ckapa=(1.+(ampi-1.)/(1.+1/(cp0**2)))**cr0
c140219 ckapa=1.+(ampi-1.)
c       denoc=1.+1./(cp0**2)
c       if(denoc.lt.1.d-20)denoc=1.d-20
c       ckapa=ckapa/denoc
c       if(ckapa.lt.1.d-20)ckapa=1.d-20
c140219 ckapa=(ckapa)**cr0
        endif   !
c250119
c       ckapa is multiple string tension
c       string tension of the pure qqbar string, kapa0, is assumed to be 1
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
c--------------------   String Tension 2 & 3 Calculating  ----------------------

c***********************   String Fragmentation 2 & 4  *************************
c       fragments strings all at once
        if(kjp22.eq.2 .or. kjp22.eq.4)then
        ps0=0.   !Lei2023060
        ps1=0.   !Lei2023060
        do ii=1,4,1   !Lei2023060
            ps0(ii) = PYP(0,ii)   ! 4-momentum before 'sfm'
        end do
        kkii=0   ! 050920
        call sfm   ! in sfm.f
c00623 Lei2023060 Throw away the current event if any errors exit.
        ! if( MSTU(23).gt.0 .OR. MSTU(30).gt.0 )then
        !     write(22,*) "iii, MSTU23, MSTU30", iii, MSTU(23), MSTU(30)   !Lei_debug
        !     iii = iii - 1
        !     goto 300
        ! end if
c00623 Lei2023060
        if(kkii.eq.2)then   ! 050920
        iii=iii-1
        goto 300   ! throw away current event
        endif
        do i1=1,4,1   !Lei2023060
            ps1(i1) = PYP(0,i1)   ! 4-momentum after 'sfm'
            throe_p(i1) = throe_p(i1) + ps0(i1) - ps1(i1)
        end do
        if(ipden.lt.11)call pyedit(2)
        if(ipden.ge.11)call pyedit(1)

c----------------------------   Gamma 77 Removing   ----------------------------
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
c----------------------------   Gamma 77 Removing   ----------------------------

        goto 30001
        endif
c***********************   String Fragmentation 2 & 4  *************************

c--------------------   String Tension 1 & 3 Calculating  ----------------------
c       calculate the single (single + multiple) effective string tension 
c        and parj(1) etc. 
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
c--------------------   String Tension 1 & 3 Calculating  ----------------------

c***********************   String Fragmentation 1 & 3  *************************
c       find string and line number of its first and last components, 
c        calculate parj(1) etc. if kjp22=1 or 3, and then hadronize string

        nstr=0
c       'pyjets' to 'sin'
        nsin=n
        do i2=1,5
        do i1=1,n
        ksin(i1,i2)=k(i1,i2)
        psin(i1,i2)=p(i1,i2)
        vsin(i1,i2)=v(i1,i2)
        enddo
        enddo
        n=0
        naff=0
        if(iii.eq.1)then   !Lei2023060
        kaff=0
        paff=0.
        vaff=0.
        endif   !Lei2023060

c       loop over string (begin)
c090219  jb=0
c090219 10000   do i1=jb+1,nsin
10001   do i1=1,nsin   ! 090219
c       find a string
        if(ksin(i1,1).eq.2)then   ! i1 is 'A'
        do i2=i1+1,nsin
        if(ksin(i2,1).eq.1)then   ! i2 is 'V'
        nstr=nstr+1
        nstra(nstr)=i1   ! line number of first component of nstr-th string
        nstrv(nstr)=i2   ! line number of last component of nstr-th string

c------------------------   Gluon Wrinkle Correction  --------------------------
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
c        current string
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
        endif   ! 5
        endif   ! 3
c------------------------   Gluon Wrinkle Correction  --------------------------

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
c090219 jb=i2   ! 160317

c--------------------------   String Fragmentation - ---------------------------
c       hadronization of current string
        iikk=0
        kkii=0   ! 050920
        ps0=0.   !Lei2023060
        ps1=0.   !Lei2023060
        do ii1=1,4,1   !Lei2023060
            ps0(ii1) = PYP(0,ii1)   ! 4-momentum before 'sfm'
        end do
        call sfm
c00623 Lei2023060 Throw away the current event if any errors exit.
        ! if( MSTU(23).gt.0 .OR. MSTU(30).gt.0 )then
        !     write(22,*) "iii, MSTU23, MSTU30", iii, MSTU(23), MSTU(30)   !Lei_debug
        !     iii = iii - 1
        !     goto 300
        ! end if
c00623 Lei2023060
c--------------------------   String Fragmentation   ---------------------------

        if(ipden.lt.11)call pyedit(2)
        if(ipden.ge.11)call pyedit(1)

c-------------------------   Failed String Treating   --------------------------
c280618
        if(iikk.eq.2 .and. ((ipden.eq.0 .and. itden.eq.0) .or.
     c   (ipden.eq.2 .and. itden.eq.2) .or. ipden.ge.11))then
        iii=iii-1
        goto 300
        elseif(iikk.eq.2)then
c       moves the part of current string in 'sin' to 'sbe
        do ii1=i1,i2
        nbe=nbe+1
        do ii2=1,5
        kbe(nbe,ii2)=ksin(ii1,ii2)
        pbe(nbe,ii2)=psin(ii1,ii2)
        vbe(nbe,ii2)=vsin(ii1,ii2)
        enddo
        enddo
        n=0   ! 230618
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
c280618
c-------------------------   Failed String Treating   --------------------------

c240219
        if(n.gt.0)then   ! 240219
        do ii1=1,4,1   !Lei2023060
            ps1(ii1) = PYP(0,ii1)   ! 4-momentum after 'sfm'
            throe_p(ii1) = throe_p(ii1) + ps0(ii1) - ps1(ii1)
        end do

c---------------------------   Gamma 77 Removing   -----------------------------
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
c---------------------------   Gamma 77 Removing   -----------------------------

c       'pyjets' to 'aaff'
        if(n.gt.0)then   ! 230618
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

        if(i2.lt.nsin)then
c230618
c       revamps 'sin', i.e. moves parton list 'sin' ii (=i2-i1+1) steps 
c        downward from i2+1 to nsin
        ii=i2-i1+1
        do jj=1,5
        do j=i2+1, nsin
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

c--------------------------   Rest Parton Dumping   ----------------------------
c       rest partons which can not compose a string 
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
        endif  !!
c--------------------------   Rest Parton Dumping   ----------------------------

c       loop over string endded
c       fragments string by string endded

20001   continue   ! from 'nstr00=nstr' to 'continue' on 030620

c----------------------------   Tension Counting   -----------------------------
        if(kjp22.eq.1 .or. kjp22.eq.3)then   ! 070219
c       average over strings in current hh collision 010518
        atime=dfloat(itime)
c       itime, # of strings in current hh collision 010518
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
c----------------------------   Tension Counting   -----------------------------

c       'aff' to 'pyjets'
        n=naff
        do ii2=1,5
        do ii1=1,n
        k(ii1,ii2)=kaff(ii1,ii2)
        p(ii1,ii2)=paff(ii1,ii2)
        v(ii1,ii2)=vaff(ii1,ii2)
        enddo
        enddo
c200420 iii=iiire   !   230219
c200420 kjp21=kjp211   ! 230219
c***********************   String Fragmentation 1 & 3  *************************

c120520 endif   ! 020512 010518
30001   continue
        endif   !! 010518 (adj12.eq.0)
c-----------------------   String Fragmentation Model  -------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c---------------------------   Coalescence  Model  -----------------------------
        if(adj12.ne.0)then
        call coales(iii,neve,nout,nap,nat,nzp,nzt)   ! coalescence
c00623 goto 333   !Lei2023060
        endif
c---------------------------   Coalescence  Model  -----------------------------
c-------------------------------------------------------------------------------


        if(ipden.lt.11)call pyedit(2)
        if(ipden.ge.11)call pyedit(1)
c141208


c-------------------------------------------------------------------------------
c----------------------   Rest Parton Re-hadronization  ------------------------
c030820 removes g,q,qbar,qq & (qq)bar from 'pyjets' to 'sbe'
        call remop   ! 030820
        call rest_hadronization   !Lei2023060
c----------------------   Rest Parton Re-hadronization  ------------------------
c-------------------------------------------------------------------------------


c------------------------------   Hadronization  -------------------------------
c*******************************************************************************



333     continue


c-------------------------------------------------------------------------------
c-------------------------------   sbh Moving  ---------------------------------
c       'sbh' to 'pyjets'
        if(nbh.ge.1)then
c00623 Lei2023060
c       Spectators in 'sbh' moved during time of 'parcas'.
        ! if(INT(adj1(12)).ne.0)then
        ! do l=1,nbh,1
        ! kf=kbh(l,2)
        ! pT2=pbh(l,1)**2+pbh(l,2)**2
        ! if( pT2.le.1D-15 .AND. (kf.eq.2212 .OR. kf.eq.2112) )then
        ! vbh(l,1) = vbh(l,1) + time_par * pbh(l,1)/pbh(l,4)
        ! vbh(l,2) = vbh(l,2) + time_par * pbh(l,2)/pbh(l,4)
        ! vbh(l,3) = vbh(l,3) + time_par * pbh(l,3)/pbh(l,4)
        ! vbh(l,4) = time_par
        ! end if
        ! end do
        ! end if
c00623 Lei2023060
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
c-------------------------------   sbh Moving  ---------------------------------
c-------------------------------------------------------------------------------


        call prt_sgam(n77s,egam3,3)   ! 080419 100919 270220
c080419 egam3: gamma energy after hadronization


c230722
998     continue


c-------------------------------------------------------------------------------
c---------------------------   NN A-loop Treating  -----------------------------
c150323
        if((ipden.eq.0.and.itden.eq.0).and.mstj2.eq.1)then
        if(inorex.eq.1)then
c       'sa2' to 'pyjets'
        n=nsa
        do i2=1,5,1
        do i1=1,nsa,1
        k(i1,i2)=ksa(i1,i2)
        p(i1,i2)=psa(i1,i2)
        v(i1,i2)=vsa(i1,i2)
        enddo
        enddo
        endif
        if(inorex.eq.2)call padecy(1,2)   ! in parini.f
        endif
c150323
c---------------------------   NN A-loop Treating  -----------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c------------------------   B-loop Gamma 77 Removing  --------------------------
        if(mstptj.eq.1.and.mstj2.eq.2)then   !!! 161222 260223
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
        endif  !!! 161222
c230722
c------------------------   B-loop Gamma 77 Removing  --------------------------
c-------------------------------------------------------------------------------

!Lei2023060B--
!       Share or not? Hadron should be on-shell?
        ! IF(.TRUE.)THEN
        IF(.FALSE.)THEN
            call share_p_PYJETS
        END IF
!Lei2023060E--

c140223 Lei full_events_history of OSC1999A
        call oscar(win,3)


c-------------------------------------------------------------------------------
c--------------------------   Hadronic Rescattering  ---------------------------
c       hadronic cascade (rescattering, HRS)
        if(kjp21.eq.1)then   ! 1 241103

c----------------------------   Hadron Specifying  -----------------------------
        call filt   ! in parini.f
        do i=1,kfmax
        nup=numbs(i)
        enddo
        nbh1=n-nup
c       nup is the number of particles kept in 'pyjets' (joints HRS)
c       nbh1 is the number of particles storing in 'sbh' (not joints HRS)
c060813 lepton is not rescattering with hadrons
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
c      'pyjets' to 'sa1_h'
        nn=n
        do i2=1,5
        do i1=1,n
        kn(i1,i2)=k(i1,i2)
        pn(i1,i2)=p(i1,i2)
        rn(i1,i2)=v(i1,i2)
        enddo
        enddo
c----------------------------   Hadron Specifying  -----------------------------

c241103
        time_had=0.d0   ! 081010
        call hadcas(iii,neve,nout,time_had,ijkk)   ! 241103   ! in hadcas.f
        if(ijkk.eq.1)then   ! 161203
c110603 iii=iii-1   ! it has been executed in 'scat' in 'hadcas'
        goto 300   ! give up current event avoiding infinite collision loop
        endif
c241103

c201203
c       'sa1_h' to 'pyjets'
        n=nn
        do i2=1,5
        do i1=1,n
        k(i1,i2)=kn(i1,i2)
        p(i1,i2)=pn(i1,i2)
        v(i1,i2)=rn(i1,i2)
        enddo
        enddo
c201203
c241103

c----------------------------   Gamma 66 Removing  -----------------------------
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
c----------------------------   Gamma 66 Removing  -----------------------------

c-------------------------------   sbh Moving  ---------------------------------
c       'sbh' to 'pyjets'
        if(nbh.eq.0)goto 9000   ! 261103
c00623 Lei2023060
c       Spectators in 'sbh' moved during time of 'hadcas'.
        ! do l=1,nbh,1
        ! kf=kbh(l,2)
        ! pT2=pbh(l,1)**2+pbh(l,2)**2
        ! if( pT2.le.1D-15 .AND. (kf.eq.2212 .OR. kf.eq.2112) )then
        ! vbh(l,1) = vbh(l,1) + time_had * pbh(l,1)/pbh(l,4)
        ! vbh(l,2) = vbh(l,2) + time_had * pbh(l,2)/pbh(l,4)
        ! vbh(l,3) = vbh(l,3) + time_had * pbh(l,3)/pbh(l,4)
        ! vbh(l,4) = time_had
        ! end if
        ! end do
c00623 Lei2023060
        do l=1,nbh
        l1=n+l
        do m=1,5
        k(l1,m)=kbh(l,m)
        p(l1,m)=pbh(l,m)
        v(l1,m)=vbh(l,m)
        enddo
        enddo
        n=n+nbh
9000    continue
c-------------------------------   sbh Moving  ---------------------------------

        endif   ! 1 241103
c--------------------------   Hadronic Rescattering  ---------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c-----------------------------   Kaon Oscillation  -----------------------------
c       change K0,K0ba to K0L and K0S
        do j=1,n
        kf=k(j,2)
        if(kf.eq.311 .or. kf.eq.-311)then
        rrlu=pyr(1)
        k(j,2)=130
        if(rrlu.gt.0.5)k(j,2)=310
        endif
        enddo
c-----------------------------   Kaon Oscillation  -----------------------------
c-------------------------------------------------------------------------------


c       perform particle, declared unstable in the 'mdcy' array, decay
c        and remove hadronic decay gammas ('22') from 'pyjets' to 'sgam'
c130205 call pyexec   ! in p_23.f


c-------------------------------------------------------------------------------
c-----------------------------   A-loop Treating   -----------------------------
c150323 
        if(mstj2.eq.1)then

c       decay the delta
c       moves delta from 'pyjets' to 'delt'
        call remo_delt   ! in 'parini_23.f'
c       'pyjets' to 'saf'
        naf=n
        do i2=1,5,1
        do i1=1,n,1
        kaf(i1,i2)=k(i1,i2)
        paf(i1,i2)=p(i1,i2)
        vaf(i1,i2)=v(i1,i2)
        enddo
        enddo
c       delta decay one by one
        do i1=1,ndel
        call padecy_delt(i1)   ! in parini.f
c       move decayed particle from 'pyjets' to 'saf'
        do i3=1,n
        naf=naf+1
        do i4=1,5
        kaf(naf,i4)=k(i3,i4)
        paf(naf,i4)=p(i3,i4)
        vaf(naf,i4)=v(i3,i4)
        enddo
        enddo
        enddo
c       'saf' to 'pyjets'
        call tran_saf   ! in 'parini_23.f

c       decay the pi^0
        i11=0
        n0=n
907     continue
        do i1=i11+1,n0
        kf=k(i1,2)
        if(kf.eq.111)then
        call pydecy(i1)
        call pyedit(2)
        i11=i1
        goto 907
        endif
        enddo

        endif
c150323
c-----------------------------   A-loop Treating   -----------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c------------------------------   Hadron Decay   -------------------------------
        rrp=1.16   ! 130205
c       if(win.lt.parp21)then   ! 260123
        if(mstj2.eq.2.or.mstj2.eq.3)call decayh(rrp)   ! in sfm.f
c       in sfm_23.f, 130205 260123 150323
c150323 call decayh(rrp)
c       endif   ! 260123
c------------------------------   Hadron Decay   -------------------------------
c-------------------------------------------------------------------------------

!Lei2023060B-- Share or not? Hadron should be on-shell?
        ! IF(.TRUE.)THEN
        IF(.FALSE.)THEN
            call share_p_PYJETS
        END IF
!Lei2023060B--

c181003
888     continue
        if(nnstop.ne.0)then
        sstop1=zstop/dfloat(nnstop)
        sstop=sstop+sstop1
        else
        nzstop=nzstop+1
        endif
c181003


c-------------------------------------------------------------------------------
c--------------------------   Information Recording   --------------------------
        if(iii.eq.1)then
        write(9,*)
        write(9,*)'iii,neve=',iii,neve
        write(9,*)'nap,nzp,nat,nzp,bp=',nap,nzp,nat,nzt,bp
        write(9,*)'sig,t0,ddt,dep=',sig,t0,ddt,dep
        write(9,*)'rou0,rao,rnp,rnt=',rou0,rao,rnp,rnt
        write(9,*)'csnn,cspin,cskn=',csnn,cspin,cskn
        write(9,*)'cspipi,cspsn,cspsm=',cspipi,cspsn,cspsm
        write(9,*)'ifram,rcsit,kfmax,ipden,itden=',
     c   ifram,rcsit,kfmax,ipden,itden   ! 060813
c00623 !Lei2023060
        n_kfmax=kfmax/10
        do i_kfmax=1,n_kfmax+1,1
            i_low=i_kfmax*10-9
            i_upp=i_kfmax*10
            if(i_kfmax.eq.(n_kfmax+1)) i_upp=kfmax
        write(9,*) (kfaco(i),i=i_low,i_upp,1)
        enddo
        do i_kfmax=1,n_kfmax+1,1
            i_low=i_kfmax*10-9
            i_upp=i_kfmax*10
            if(i_kfmax.eq.(n_kfmax+1)) i_upp=kfmax
        write(9,*) (disbe(i,i),i=i_low,i_upp,1)
        enddo
c00623 !Lei2023060
c00623 write(9,*)(kfaco(i),i=1,kfmax)
c00623 write(9,*)(disbe(i,i),i=1,kfmax)
        write(9,*)(disbe(1,i),i=1,10)
        write(9,*)'isinel='
        write(9,600)isinel
        endif
600     format(25(1x,i2)/)
c050603
c--------------------------   Information Recording   --------------------------
c-------------------------------------------------------------------------------


        if(ipden.lt.11)call pyedit(2)   ! 060813
        if(ipden.ge.11)call pyedit(1)   ! 060813


c291218
        if(adj140.lt.3)goto 80004   !Lei2023060 .ne.4 -> .lt.3
        if(adj140.eq.4.and.mstj2.eq.1)goto 80004   ! 150323
80006   continue   ! 120720
c       moves partons from 'pyjets' to 'sbe'
        call remop   !Lei2023060 Keep it here for error checking.


c-------------------------------------------------------------------------------
c------------------------   Final Information Output   -------------------------
c       pythia output
80004   continue
        if(nout.eq.1 .or. iii.eq.1 .or. mod(iii,nout).eq.0 .or. iii
     c   .eq.neve)then
        write(mstu(11),*)
        write(mstu(11),*)
        write(mstu(11),*)"#!-------------------------------------"//
     &                   "----------------------------------------"
        write(mstu(11),*)'event=',iii
c       call prt_pyj(n,cc)   ! 190922   ! in parini.f
        call pylist(1)
        call prt_final_info(win)   !Lei2023060
        if(iii.eq.1) write(9,*)'# of parton scaterring process '//
     &                          'undergone in an event'
        write(9,108)iii,(INT(reac(i1)),i1=1,9)
        endif
108     format(10(1x,i4))
c------------------------   Final Information Output   -------------------------
c-------------------------------------------------------------------------------


c       ??????????????? OSCAR standard output ??????????????????????
        call oscar(win,4)   ! 140223 Lei Outputs if nosc > 0
c       ??????????????? OSCAR standard output ??????????????????????
        call PASTAT(0,0)  !Lei202306 Counts total-event cross sections.


c-------------------------------------------------------------------------------
c-----------------------------   Event Analysis   ------------------------------

!Lei20230311B
        IF(.TRUE.)THEN

        if(iii.eq.1)then
            mult_h_full = 0
            mult_h_partial_pT_rap = 0
            mult_h_partial_pT     = 0
            mult_h_partial_rap    = 0
        end if
        d_pT = 20./100.
        d_y  = 12./100.
        d_eta= 12./100.
        do i=1,N,1
            kf_a = ABS(K(i,2))
            pT   = PYP(i,10)
            y    = PYP(i,17)
            eta  = PYP(i,19)
            rap  = y
            rap_low = y_low
            rap_upp = y_upp
            if( INT(parp22).eq.1 )then
                rap = eta
                rap_low = eta_low
                rap_upp = eta_upp
            end if
            d_rap = rap_upp - rap_low
            ! h+- only, i.e. pi+-, K+- and p+-
            if(kf_a.ne.211 .AND. kf_a.ne.321 .AND. kf_a.ne.2212) cycle
            ! Excludes spectators.
            ! if(kf_a.eq.2212 .AND. pT.le.1D-15) cycle
            ! Does not counts h+- with too small pT.
            if(pT.le.1D-10) cycle
            ! Multiplicity
            mult_h_full = mult_h_full + 1
            if( rap.gt.rap_low .AND. rap.lt.rap_upp .AND. 
     &          pT.gt.pT_low   .AND. pT.lt.pT_upp )then
                mult_h_partial_pT_rap = mult_h_partial_pT_rap + 1
            end if
            ! Inv. dN/dpT
            if(rap.gt.rap_low .AND. rap.lt.rap_upp)then
                mult_h_partial_rap = mult_h_partial_rap + 1
                CALL PYFILL(9, pT, 1./pT/PARU(2)/d_pT/d_rap/2.)
                ! dN/dpT
                CALL PYFILL(15, pT, 1./d_pT/d_rap/2.)
            end if
            CALL PYFILL(10, pT, 1./pT/PARU(2)/d_pT/d_rap/2.)
            ! dN/dpT
            CALL PYFILL(16, pT, 1./d_pT/d_rap/2.)
            ! dN/drap
            if(pT.gt.pT_low .AND. pT.lt.pT_upp)then
                mult_h_partial_pT = mult_h_partial_pT + 1
                CALL PYFILL(11,    y, 1./d_y)
                CALL PYFILL(13,  eta, 1./d_eta)
            end if
            CALL PYFILL(12,    y, 1./d_y)
            CALL PYFILL(14,  eta, 1./d_eta)
        end do
        CALL PYFILL(17, time_ini, 1D0)
        CALL PYFILL(18, time_par, 1D0)
        CALL PYFILL(19, time_had, 1D0)
        sum_reac=0.
        do i=1,9,1
        sum_reac = sum_reac + reac(i)
        end do
        CALL PYFILL(20, sum_reac, 1D0)
        sum_reac=0.
        do i=1,600,1
        sum_reac = sum_reac + noinel(i)*1.
        end do
        sum_reac = sum_reac + noel*1.
        CALL PYFILL(21, sum_reac, 1D0)

        END IF
!Lei20230311E

c       analyse an event
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
c       write(98,190)bp,ngam
c       do j=1,ngam
c       ik=kgam(j,2)
c       p1=pgam(j,1)
c       p2=pgam(j,2)
c       p3=pgam(j,3)
c       p4=pgam(j,4)
c       write(98,191) ik, p1, p2, p3, p4
c       enddo
c190    format(f10.4,I7)   ! sa
c191    format(I6,4(1x,e15.7))   ! sa

c       analyses the event on-line
        call analy(nmin,nminf,ncha,nchaf,parp22)   ! 281219 140223 Lei ! in analy.f
        call analy_parton(parp22)   !Lei2023060   ! in analy.f
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
        call prt_sgam(0,egam,4)   !Lei20230819
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
c140820 swoun: # of wounded nucleons sumed up over enents
c170121 snpctl0: # of nn collision pairs (in parini.f) sumed up over enents
c280722 snpar: # of collided nucleons (in parini.f) sumed up over enents
c071103
        rinel=0   !Lei2023060
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
c-----------------------------   Event Analysis   ------------------------------
c-------------------------------------------------------------------------------


        open(8,file='nout.out',status='unknown')
        write(8,*)'iii=',iii
        close(8)


c-------------------------------------------------------------------------------
c--------------------------   Event Averaged Output   --------------------------
c       internal and final printing and controled return
        if(mod(iii,nout).eq.0 .or. iii.eq.neve)then

        open(10,file='rms.out',status='unknown')
        flaa=dfloat(iii-ich)   ! July/20/98 ich=0
        if(flaa.le.1.e-20)goto 1200

c-----------------------------   Event Averaging   -----------------------------
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
        nrea=0   !Lei2023060
        do i1=1,9
        nrea=nrea+nreac(i1)
        snreac(i1)=nreac(i1)/flaa
        enddo
        srea=float(nrea)/flaa
c220110
c-----------------------------   Event Averaging   -----------------------------

1200    continue

c-------------------------------   User Output   -------------------------------
c00623 Lei2023060
c       User output
c       Header.   Lei2023060
        write(10,*)'#!**************************************'//
     &             '*************************************!#'
        write(10,*)'#!*********************|'//
     &             '    PACIAE  Analysis  Output    '//
     &             '|********************!#'
        write(10,*)'#!**************************************'//
     &             '*************************************!#'

c       Records the current real date and time.   !Lei2023060
        call DATE_AND_TIME(VALUES=n_current_date_and_time)
        do i1=1,8
        write(c_date_and_time(i1),"(I4)") n_current_date_and_time(i1)
        enddo
        write(10,*)
        write(10,*)'#! Now is   ' //TRIM(ADJUSTL(c_date_and_time(5))) //
     &                        ':' //TRIM(ADJUSTL(c_date_and_time(6))) //
     &                        ':' //TRIM(ADJUSTL(c_date_and_time(7))) //
     &                      '   ' //TRIM(ADJUSTL(c_date_and_time(3))) //
     &                        '/' //TRIM(ADJUSTL(c_date_and_time(2))) //
     &                        '/' //TRIM(ADJUSTL(c_date_and_time(1)))
c171022 Records the seed of random number generator.
        write(10,*)'#! Seed (PYTHIA default=19780503) =', MRPY(1)   ! 171022 Lei
        write(10,*)

        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! parp81, parp82, bp, mstp82 ='    ! 291207 Lei2023060 parp82
        write(10,*) parp81,parp82,bp,mstp(82)    ! 291207 Lei2023060 parp82

        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! MC Glauber-like <N_coll>, <N_part> ='
        write(10,*) snpctl0i,snpari   ! 180121 280722
        write(10,*)'#! largest ave. # of NN collision pairs =' ! 280722
        write(10,*) snpctlmi   ! 280722
c140820
        write(10,*)'#! ave. # of NN collision pairs calling pythia, '//
     c             'not calling pythia ='
        write(10,*) eineli(592),eineli(593)
        write(10,*)'#! ave. # of wounded nucleons in parini ='
        write(10,*) swouni
c140820 
        write(10,*)'#! colli. # suffered by projectile nucleon '//
     c             'in target nucleus'   ! 140219
        write(10,*) spathni   ! 140219
        write(10,*)'#! event averaged N_bin'   ! 260219
        write(10,*) sevbini   ! 260219
c200601
        write(10,*)'#! (Npart)mini-jet, Nnn, Npp='
        write(10,*) skparo,sknno,skppo
        write(10,*)'#! Nnp, Ntot, Nep='   ! 060813
        write(10,*) sknpo,sknno+skppo+sknpo,skepo   ! 060813
c200601
!Lei20230820B-
        if( INT(psno).eq.1 )then
            write(10,*)'#! event averaged b, avneu, Npart_p, '
     &               //'Npart_t, T_pt='   !Lei2023060 Lei20230820
            write(10,*) avb,avneu,astbp,astbt,aanbin   ! 280113
        else if( INT(psno).eq.2 )then
            write(10,*)'#! psno, ave. b, N_part and N_bin ='   ! 280113
            write(10,*) psno,averbo,psnopo,psnoto,psnono*csnn   ! 280113
        else
            write(10,*)
            write(10,*)
        end if
!Lei20230820E-

        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! event averaged energy of gamma after '//
     c             'partonic initiation, partonic cascade,'
        write(10,*)'#!  hadronization and end of event =' ! 080419
        write(10,*) segam1o,segam2o,segam3o,segamo ! 080419
        if(ipden.ge.11.and.ipden.le.16)
     c   write(10,*)'#! event average number of lepton studied ='   !260314
        if(ipden.ge.11.and.ipden.le.16) write(10,*) vnlep/flaa   !260314

c071103
        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! # of successful, blocked and all collision '//
     c   'in parton cascade ='
        write(10,*) rineli,reli,reli+rineli
c071103
c220110
        write(10,*)'#! average collision # in parton cascade ='
        write(10,*) srea
        write(10,*)'#! # of scaterring processes in parton cascade'
        write(10,*) (snreac(i1),i1=1,3)
        write(10,*) (snreac(i1),i1=4,6)
        write(10,*) (snreac(i1),i1=7,9)
c220110

        write(10,*)'#! average frequency of the occurring of each '//
     c   'inela. in hadron cascade (at the end of the file)'   ! 140223 Lei
c140223 write(10,*)dineli   ! 140223 Lei
        write(10,*)'#! el. and inel. coll. # and sum in hadron cascade='
        write(10,*) seli,sineli,seli+sineli

c070417
        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! default parj1, parj2, parj3, parj4, par21 ='   ! 010518
        write(10,*) parj1,parj2,parj3,parj4,parj21   ! 010518
        write(10,*)'#! Eff-parj1, parj2, parj3, parj4, parj21, keff ='   ! 010518
        write(10,*) (skapao(i1),i1=1,6)   ! 010518
        write(10,*)'#! averaged # of gluon in a string when kjp22=1,3'
        write(10,*) sgtimeo
        write(10,*)'#! event averaged value of the factor related to # '
        write(10,*)'#!  of gluons and hardest gluon in a string, event '
        write(10,*)'#!  averaged transverse momentum of hardest gluon,'
        write(10,*)'#!  event averaged # strings when kjp22=1,3 ='
        write(10,*) sadivo,sgpmaxo,sitimeo
c070417
c120119
        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! times & sum='
        write(10,*) stime_ini,stime_par,stime_had,
     c              stime_ini+stime_par+stime_had
c120119

        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! q, qbar, charge thrown away ='
        write(10,*) wthroq,wthrob,wthroc/3.
        write(10,*)'#! 3-momentum and energy thrown away ='
        write(10,*) wthroe

        write(10,*)
        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! multiplicity of negative, positive particles '//
     &             'and sums, partial & full ='
        write(10,*) dnmino, dnchao, dnmino+dnchao
        write(10,*) dnminfo, dnchafo, dnminfo+dnchafo
c00623 write(10,*)"#!-------------------------------------"//
c    &             "----------------------------------------"
c       write(10,*)'#! particle multiplicity, partial ='
c       write(10,*) (sbo(ll),ll=1,ispmax)
c       write(10,*)'#! particle multiplicity, full    ='
c00623 write(10,*) (sbof(ll),ll=1,ispmax)   !Lei2023060 in output_hadron_distribution now
csa****************************************************************

c00623 do m2=1,isdmax   !Lei2023060 in output_hadron_distribution now
c       write(10,*)'ID of distribution m2=',m2
c       do m3=1,ispmax
c       write(10,*)'distribution belong to m3=',m3
c       write(10,*)(sao(m1,m2,m3),m1=1,40)   ! 070419
c       write(10,*)(saof(m1,m2,m3),m1=1,40)   ! 070419
c       enddo
c00623 enddo

c00623 Lei2023060B-
c       Outputs multiplicities, abscissa, 5-distributions of pi+K+p and 20 
c        particles specified in usu.dat.
        call output_hadron_distribution(sao,sbo,saof,sbof)   ! in analy.f
c       Outputs multiplicities, abscissa, 5-distributions of g, u+d+s +anti- 
c        and q/qbar.
        call output_parton_distribution   ! in analy.f
c00623 Lei2023060E-

c140223 Lei
        write(10,"(2/)")   ! 140223 Lei Empty lines.
        write(10,*)'#! average frequency of the occurring of each '//
     c   'inela. in hadron cascade ='
        do i=1,60,1
            j=(i-1)*6
            write(10,*) (dineli(j+kt),kt=1,6,1)
        end do
c140223 Lei

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
        write(10,*)'#! relative multiplicity, p =',(sbo(ll),ll=1,10)
        write(10,*)'#! relative multiplicity, f =',(sbof(ll),ll=1,10)
        do m2=1,isdmax
        write(10,*)'#! ID of relative distribution m2 =',m2
        do m3=1,10
        write(10,*)'#! distribution belong to m3 =',m3
        write(10,*)(sao(m1,m2,m3),m1=1,40)   ! 070419
        write(10,*)(saof(m1,m2,m3),m1=1,40)   ! 070419
        enddo
        enddo
        endif
c260314

!Lei20230303B--
        sum_y = 0.
        sum_pT = 0.
        sum_eta = 0.
        sum_pT2 = 0.

        do i_kf=1,20,1
            if(ispkf(i_kf).eq. 211)  i_h(1) = i_kf
            if(ispkf(i_kf).eq.-211)  i_h(2) = i_kf
            if(ispkf(i_kf).eq. 321)  i_h(3) = i_kf
            if(ispkf(i_kf).eq.-321)  i_h(4) = i_kf
            if(ispkf(i_kf).eq. 2212) i_h(5) = i_kf
            if(ispkf(i_kf).eq.-2212) i_h(6) = i_kf
        end do

        do m1=1,40,1
            sum_y(m1,1) = sao(m1,1,i_h(1)) + sao(m1,1,i_h(2))
     &                  + sao(m1,1,i_h(3)) + sao(m1,1,i_h(4))
     &                  + sao(m1,1,i_h(5)) + sao(m1,1,i_h(6))
            sum_pT(m1,1) = sao(m1,2,i_h(1)) + sao(m1,2,i_h(2))
     &                   + sao(m1,2,i_h(3)) + sao(m1,2,i_h(4))
     &                   + sao(m1,2,i_h(5)) + sao(m1,2,i_h(6))
            sum_eta(m1,1) = sao(m1,3,i_h(1)) + sao(m1,3,i_h(2))
     &                    + sao(m1,3,i_h(3)) + sao(m1,3,i_h(4))
     &                    + sao(m1,3,i_h(5)) + sao(m1,3,i_h(6))
            sum_pT2(m1,1) = sao(m1,6,i_h(1)) + sao(m1,6,i_h(2))
     &                    + sao(m1,6,i_h(3)) + sao(m1,6,i_h(4))
     &                    + sao(m1,6,i_h(5)) + sao(m1,6,i_h(6))

            sum_y(m1,2) = saof(m1,1,i_h(1)) + saof(m1,1,i_h(2))
     &                  + saof(m1,1,i_h(3)) + saof(m1,1,i_h(4))
     &                  + saof(m1,1,i_h(5)) + saof(m1,1,i_h(6))
            sum_pT(m1,2) = saof(m1,2,i_h(1)) + saof(m1,2,i_h(2))
     &                   + saof(m1,2,i_h(3)) + saof(m1,2,i_h(4))
     &                   + saof(m1,2,i_h(5)) + saof(m1,2,i_h(6))
            sum_eta(m1,2) = saof(m1,3,i_h(1)) + saof(m1,3,i_h(2))
     &                    + saof(m1,3,i_h(3)) + saof(m1,3,i_h(4))
     &                    + saof(m1,3,i_h(5)) + saof(m1,3,i_h(6))
            sum_pT2(m1,2) = saof(m1,6,i_h(1)) + saof(m1,6,i_h(2))
     &                    + saof(m1,6,i_h(3)) + saof(m1,6,i_h(4))
     &                    + saof(m1,6,i_h(5)) + saof(m1,6,i_h(6))
        end do

        sum_pT  = sum_pT  / 2. / PARU(2)
        sum_pT2 = sum_pT2 / 2.
        do m1=1,40,1
            sum_pT(m1,1)  = sum_pT(m1,1)  / (afl(1,1,2)-afl(1,1,1))
            sum_pT2(m1,1) = sum_pT2(m1,1) / (afl(1,1,2)-afl(1,1,1))
        end do
        sum_y = sum_y
        sum_eta = sum_eta
        write(10,"(2/)")
        write(10,*) "#! Inv. dN/dpT, dN/dpT, dN/dy, dN/deta of h+-,"//
     &              " p and f"
        do m1=1,40,1
            write(10,*) sum_pT(m1,1),sum_pT(m1,2),
     &                  sum_pT2(m1,1),sum_pT2(m1,2),
     &                  sum_y(m1,1), sum_y(m1,2),
     &                  sum_eta(m1,1),sum_eta(m1,2)
        end do
        sum_mul_partial = 0.
        sum_mul_full = 0.
        do i_kf=1,6,1
            sum_mul_partial = sum_mul_partial + sbo(i_h(i_kf))
            sum_mul_full = sum_mul_full + sbof(i_h(i_kf))
        end do
        write(10,*)
        write(10,*) "#! partial/full multiplicity of h+-:"
        write(10,*) sum_mul_partial, sum_mul_full
!Lei20230303E--

        close(10)
c00623 Lei2023060
c-------------------------------   User Output   -------------------------------

!Lei20230303B
        mstu(11)=99
        open(99,file='h+-.out',status='unknown')
        if(iii.eq.neve)then
            CALL PYFACT(9, 1D0/neve)
            CALL PYFACT(10,1D0/neve)
            CALL PYFACT(11,1D0/neve)
            CALL PYFACT(12,1D0/neve)
            CALL PYFACT(13,1D0/neve)
            CALL PYFACT(14,1D0/neve)
            CALL PYFACT(15,1D0/neve)
            CALL PYFACT(16,1D0/neve)
        end if
        write(99,*) "#!------------------------------------------------"
        write(99,*) "#! Multiplicity, full, (-pT and -rap), -pT, -rap"
        write(99,*) mult_h_full*1D0/iii, mult_h_partial_pT_rap*1D0/iii,
     &           mult_h_partial_pT*1D0/iii, mult_h_partial_rap*1D0/iii
        write(99,*) "#!------------------------------------------------"
        write(99,*) "#! (1/2*pi*pT)*d2N/dpT/drap, partial rap, full",
     &              rap_low, "< rap <", rap_upp
        CALL PYDUMP(3,99,2,i_hist(9))
        write(99,*) "#!------------------------------------------------"
        write(99,*) "#! d2N/dpT/drap, partial rap, full",
     &              rap_low, "< rap <", rap_upp
        CALL PYDUMP(3,99,2,i_hist(15))
        write(99,*) "#!------------------------------------------------"
        write(99,*) "#! dN/dy, dN/deta, partial pT, full"
        CALL PYDUMP(3,99,4,i_hist(11))
        write(99,*) "#!------------------------------------------------"
        write(99,*) "#! time parini"
        CALL PYDUMP(3,99,1,i_hist(17))
        write(99,*) "#! time parcas"
        CALL PYDUMP(3,99,1,i_hist(18))
        write(99,*) "#! time hadcas"
        CALL PYDUMP(3,99,1,i_hist(19))
        write(99,*) "#! Ncoll parcas, hadcas"
        CALL PYDUMP(3,99,2,i_hist(20))

        CALL PYPLOT(9)
        CALL PYPLOT(10)
        CALL PYPLOT(11)
        CALL PYPLOT(12)
        CALL PYPLOT(13)
        CALL PYPLOT(14)
        CALL PYPLOT(15)
        CALL PYPLOT(16)
        CALL PYPLOT(17)
        CALL PYPLOT(18)
        CALL PYPLOT(19)
        CALL PYPLOT(20)
        CALL PYPLOT(21)
        close(99)
        mstu(11)=22
!Lei20230303E

        endif
c--------------------------   Event Averaged Output   --------------------------
c-------------------------------------------------------------------------------


1000    if(iii.lt.neve)then
c260718 if(dabs(bmin-bmax).lt.10d-4)goto 300
        if(psno.eq.0 .or. psno.eq.2)goto 300   ! 280113 260718
        if(psno.eq.1 .and. jjj.eq.10)then   ! 260718
c       10: total number of impact paremeters in systematic sampling for impact
c           parameter
        jjj=1
        goto 300
        elseif(psno.eq.1 .and. jjj.lt.10)then   ! 260718
        jjj=jjj+1
        goto 300
        endif   ! 260718
        endif
c----------------------------   Event Generating   -----------------------------
c-------------------------------------------------------------------------------
c*******************************************************************************




        write(9,*)'nncoll (total diffractive NN)=',nncoll   ! sa 060814 Lei2023060
c060813 statistics of processes generated
        write(MSTU(11),*)
        write(MSTU(11),*)
c0623 call pystat(1)   ! 060813   in p_23.f !Lei202306
        call PASTAT(1,0)   !Lei202306 in analy.f

        close(2)
        close(3)
        close(5)
        close(9)
        close(22)
        if(nosc.gt.0) close(34)   ! 140223 Lei
c       close(98)   ! 260219
c       close(99)
c       timeb=dtime(ty)
c       write(9,*)'time consuming =',timeb


        stop
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remop
c       moves q,qbar,g,diquark, and anti-diquark from 'pyjets' to 'sbe'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
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
c       moves particle list 'pyjets' one step downward from i1+1 to n
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
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
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
     c   nap,nat,nzp,nzt,pio
        common/sa34/itorw,iikk,cp0,cr0,kkii   !Lei2023060
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
        ngg=jj   ! order # of last gluon in 'sbh' !Lei2023060 i1 -> jj

        if(nbh.gt.ngg)then   ! upto here jj=ngg

        do i2=ngg+1,nbh
        kff=kbh(i2,2)
        if(kff.lt.0)then   ! anti-quark
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
        nqb=jj   ! order # of last anti-quark !Lei2023060 i2 -> jj

        endif
c       making 'sbh' in order of g, qbar, & q finished

c       fragmenting 'sbh'

c       fragmenting gg in 'sbh'
        nggg=0
        if(ngg.gt.0)then   ! 1
        n_end = ngg   !Lei2023060
        if( MOD(ngg,2).eq.1 ) n_end = ngg - 1   !Lei2023060 Odd ngg
        do i1=1,n_end,2   ! 2 !Lei2023060 ngg -> n_end
        kf1=kbh(i1,2)
        er1=pbh(i1,4)
        kf2=kbh(i1+1,2)
        er2=pbh(i1+1,4)
        ss=er1+er2
        call py2ent(0,kf1,kf2,ss)   ! in p_23.f
c00623 Lei2023060
        if(iikk.eq.2 .OR. kkii.eq.2)then   ! Failed, recover PYJETS
            n = 0
            ! Throw away them to "sbe"
            nbe = nbe+1
            do jj=1,5
                kbe(nbe,jj) = kbh(i1,jj)
                pbe(nbe,jj) = pbh(i1,jj)
                vbe(nbe,jj) = vbh(i1,jj)
            enddo
            nbe = nbe+1
            do jj=1,5
                kbe(nbe,jj) = kbh(i1+1,jj)
                pbe(nbe,jj) = pbh(i1+1,jj)
                vbe(nbe,jj) = vbh(i1+1,jj)
            enddo
            cycle
        end if
c00623 Lei2023060
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
c00623 if(nggg.lt.ngg)then   !Lei2023060
        if( MOD(ngg,2).eq.1 )then   !Lei2023060
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
c00623 Lei2023060
        if(iikk.eq.2 .OR. kkii.eq.2)then   !Failed, recover PYJETS
            n = 0
            ! Throw away them to "sbe"
            nbe = nbe+1
            do jj=1,5
                kbe(nbe,jj) = kbh(ngg+1,jj)
                pbe(nbe,jj) = pbh(ngg+1,jj)
                vbe(nbe,jj) = vbh(ngg+1,jj)
            enddo
            nbe = nbe+1
            do jj=1,5
                kbe(nbe,jj) = kbh(nqb+1,jj)
                pbe(nbe,jj) = pbh(nqb+1,jj)
                vbe(nbe,jj) = vbh(nqb+1,jj)
            enddo
            ! cycle
            goto 200
        end if
c00623 Lei2023060
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
200     continue   !Lei2023060
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
        n_begin = ngg + 1   !Lei2023060
        if(ngg.eq.0) n_begin = 2  !Lei2023060
        do i1=n_begin,nbh,1   !Lei2023060 ngg+1 -> n_begin, explicit 1 if N: qbar = q = 1
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
        nbh=0   !Lei2023060
c       fragmenting 'sbh' finished

        return
        end



        subroutine rest_hadronization
c00623 Lei2023060
c       Hadronize the rest partons ("sbe") those failed in sfm/coal by 
c        calling two subroutines "rest_sfm" and/or "rest_coal".
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
c       PYDAT1,PYDAT2,PYDAT3 and PYJETS are the subroutines in PYTHIA.
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_c/ithroq,ithrob,ithroc,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5)
        common/ctllist_p/nreac(9),nrel

c       Dumps "PYJETS" into "aaff".
        naff = 0
        call tran_pyjets   ! in main.f
        n = 0

c       Appends "trs" to "sbe". (partons of inel. coll. in parcas with sfm)
        if( INT(adj1(12)).eq.0. .AND. iparres.eq.1 .AND. ntrs.gt.0 )then
            do ii1=1,ntrs,1
                ii3 = nbe + ii1
                do ii2=1,5,1
                    kbe(ii3,ii2) = ktrs(ii1,ii2)
                    pbe(ii3,ii2) = ptrs(ii1,ii2)
                    vbe(ii3,ii2) = vtrs(ii1,ii2)
                enddo
            enddo
            nbe  = nbe + ntrs
            ntrs = 0
        endif

c       Appends "sa37" to "sbe". (partons failed in coal)
        if( INT(adj1(12)).ne.0 .AND. nth.gt.0 )then   !Lei20230810 .eq.1 -> .ne.0
            do ii1=1,nth,1
                ii3 = nbe + ii1
                do ii2=1,5,1
                    kbe(ii3,ii2) = kth(ii1,ii2)
                    pbe(ii3,ii2) = pth(ii1,ii2)
                    vbe(ii3,ii2) = vth(ii1,ii2)
                enddo
            enddo
            nbe = nbe + nth
            nth = 0
        end if

c       Breaks up potential diquarks in "sbe" if sfm.
        if( INT(adj1(12)).eq.0 .AND. (nreac(4).ne.0.or.nreac(6).ne.0 
     &      .OR. (nreac(7).ne.0 .AND. mstj1_2.eq.1)) ) call break_sbe   ! in main.f

        if( INT(adj1(18)).eq.0 .AND. nbe.gt.1 ) call rest_sfm   ! in main.f

        if( INT(adj1(18)).ge.0 .AND. ( (INT(adj1(12)).eq.0 .AND. 
     &      nbe.gt.0) .OR. INT(adj1(12)).ne.0) )
     &  call rest_coal   ! Do coalescence again even though adj18=0.   ! in main.f

c       Appends "aaff" to "PYJETS". (dumps "aaff" into "PYJETS" in fact)
        do ii1=1,naff,1
            ii3 = N + ii1
            do ii2=1,5,1
                K(ii3,ii2) = kaff(ii1,ii2)
                P(ii3,ii2) = paff(ii1,ii2)
                V(ii3,ii2) = vaff(ii1,ii2)
            enddo
        enddo
        N = N + naff
        naff = 0


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rest_sfm
c00623 Lei2023060
c       Hadronizes "sbe" with "PYEXEC".
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
c       PYDAT1,PYDAT2,PYDAT3 and PYJETS are the subroutines in PYTHIA.
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/aaff/naff,nonaf,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
        common/sa34/itorw,iikk,cp0,cr0,kkii
c       Local variables.
        integer n_g, n_q, n_b   ! g, q, qbar
        dimension p0(4), dp(4), ps1(6), ps0(6), rc(3)
        dimension k_g(100,5), p_g(100,5), v_g(100,5)   ! n_g for gluon
        dimension k_q(100,5), p_q(100,5), v_q(100,5)   ! n_q for quark/di-
        dimension k_b(100,5), p_b(100,5), v_b(100,5)   ! n_b for anti-quark/di-

        if(nbe.gt.100) write(*,*) "Warning! nbe overflows in rest_sfm!"
        if(nbe.gt.100) stop

c       Shares the 4-momentum in 'throe_p' among partons.
        call share_p_sbe   ! in main.f

c       Dumps ("sbe") g to "n_g", q to "n_q", and qbar to "n_b", respectively.
        n_g = 0
        k_g = 0
        p_g = 0.
        v_g = 0.
        do i=1,nbe,1
            if( kbe(i,2).eq.21 )then   ! g
                n_g = n_g + 1
                do j=1,5,1
                    k_g(n_g,j) = kbe(i,j)
                    p_g(n_g,j) = pbe(i,j)
                    v_g(n_g,j) = vbe(i,j)
                end do
            end if
        end do

        n_q = 0
        k_q = 0
        p_q = 0.
        v_q = 0.
        do i=1,nbe,1
            if( kbe(i,2).gt.0 .AND. kbe(i,2).ne.21 )then   ! q/diquark
                n_q = n_q + 1
                do j=1,5,1
                    k_q(n_q,j) = kbe(i,j)
                    p_q(n_q,j) = pbe(i,j)
                    v_q(n_q,j) = vbe(i,j)
                end do
            end if
        end do

        n_b = 0
        k_b = 0
        p_b = 0.
        v_b = 0.
        do i=1,nbe,1
            if( kbe(i,2).lt.0 )then   ! qbar/anti-diquark
                n_b = n_b + 1
                do j=1,5,1
                    k_b(n_b,j) = kbe(i,j)
                    p_b(n_b,j) = pbe(i,j)
                    v_b(n_b,j) = vbe(i,j)
                end do
            end if
        end do

        nbe = 0

c       Fragments pure single g-...-g string in "n_g".
        iteration = 0   ! No more than 50 times.
100     continue
        if( n_g.gt.1 )then
            do i=1,n_g,1
                do kk=1,5,1
                    P(i,kk) = p_g(i,kk)
                    V(i,kk) = v_g(i,kk)
                end do
                K(i,1) = 2   ! 'A'-'I'-
                K(i,2) = k_g(i,2)
                K(i,3) = 0
                K(i,4) = 0
                K(i,5) = 0
            end do
            K(n_g,1) = 1   ! 'V'
            N = n_g
            iikk = 0
            kkii = 0
c           Sums of px, py, pz, E, inv. m, and charge before the "PYEXEC".
            ps0  = 0.
            do i=1,6,1
                ps0(i) = PYP(0,i)
            end do
c           Sets fragmentation flag in PYEXEC.
            if(ipden.eq.2 .AND. itden.eq.2)then
                mstj(1)   = 1
            else
                mstp(111) = 1
            endif

            call PYEXEC   ! Consider px, py, pz, E.   ! in p_23.f

c           Sums of px, py, pz, E, inv. m, and charge after the "PYEXEC".
            ps1 = 0.
            do i=1,6,1
                ps1(i) = PYP(0,i)
            end do
c           Charge is not conserved or errors occur. Re-generate the event.
            if( ABS(ps0(6)-ps1(6)).gt.1D-10 .OR.
     &          iikk.eq.2 .OR. kkii.eq.2 )then
                N = 0
                iteration = iteration + 1
                if(iteration.gt.50) then   ! Throws away these gluons.
                    do kk=1,5,1
                        do i=1,n_g,1
                            kbe(i,kk) = k_g(i,kk)
                            pbe(i,kk) = p_g(i,kk)
                            vbe(i,kk) = v_g(i,kk)
                        end do
                    end do
                    nbe = n_g
                    n_g = 0
                    goto 200
                end if
                goto 100
            end if

            ! Success.
            do i=1,4,1
            throe_p(i) = throe_p(i) + ps0(i) - ps1(i)
            end do

            if(ipden.lt.11) call pyedit(1)   ! in p_23.f
            if(ipden.ge.11) call pyedit(1)   ! in p_23.f

c           Moves "77" from "PYJETS" to "sgam".
            if(N.gt.0)then
                n77 = 0
                do jj=1,N,1
                    kf = K(jj,2)
                    if(kf.eq.22)then
                        ! '77': photons after hadronization of current string
                        K(jj,2) = 77
                        n77     = n77 + 1
                    endif
                enddo
                if(n77.gt.0) call remo_gam_hadro(77)   ! in main.f
            endif

c           Dumps produced hadrons (from PYEXEC, no gamma) into "sa1_h".
            do jj=1,5,1
                do ii=1,N,1
                    kn(ii,jj) = K(ii,jj)
                    pn(ii,jj) = P(ii,jj)
                    rn(ii,jj) = V(ii,jj)
                enddo
            enddo
            nn  = N
c           The first position in "sa1_h" is selected from original one gluon randomly.
            i_g = INT( PYR(1)*n_g + 1 )
            do jj=1,5,1
                rn(1,jj) = v_g(i_g,jj)
            end do

c           Arranges produced particles on the surface of sphere (with 
c            radius rrp and centred on parent position), first produced 
c            particle is put on its parent position originally.
            rrp = 1.16
            rc(1)    = rn(1,1)
            rc(2)    = rn(1,2)
            rc(3)    = rn(1,3)
            rn(1,4)  = 0.   !TODO(Lei2023060): may need to be extended.
c           In corresponding with the time set in >posi".
            call posi(1,nn,rc,rrp)   ! in sfm.f

c           Transfers four position messages from "sa1_h" to "PYJETS".
            do jj=1,4,1
                do ii=1,nn,1
                    v(ii,jj) = rn(ii,jj)
                enddo
            enddo

            n_g = 0
c           Appends "PYJETS" to "aaff".
            call tran_pyjets   ! in main.f
            N = 0
        end if

200     continue

        sn = (n_q + n_b)*1.
        if( sn.le.0.9 ) goto 400   ! No q and qbar.

c       Share the 4-momentum in 'throe_p' among q and qbar.
501     continue
        dp = throe_p/sn
        p0 = throe_p   !  Old
        do i1=1,n_q,1   ! Share with q.
            p1  = p_q(i1,1)
            p2  = p_q(i1,2)
            p3  = p_q(i1,3)
            p4  = p_q(i1,4)
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2 = (p4+dp(4))**2 - (p1+dp(1))**2 - (p2+dp(2))**2
     &                                          - (p3+dp(3))**2
c           inv. m^2 >= 0 or < 0 but closer to 0, E > 0, then share.
            if( (sm2.ge.0. .OR. (sm2.lt.0. .AND. sm2.ge.sm2_org)) .AND. 
     &          (p4+dp(4)).gt.0. )then
                do i2=1,4,1
                    p_q(i1,i2) = p_q(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do   ! Compare the old and new one
        do i1=1,n_b,1  ! Share with qbar.
            p1  = p_b(i1,1)
            p2  = p_b(i1,2)
            p3  = p_b(i1,3)
            p4  = p_b(i1,4)
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2 = (p4+dp(4))**2 - (p1+dp(1))**2 - (p2+dp(2))**2
     &                                          - (p3+dp(3))**2
c           inv. m^2 >= 0 or < 0 but closer to 0, E > 0, then share.
            if( (sm2.ge.0. .OR. (sm2.lt.0. .AND. sm2.ge.sm2_org)) .AND. 
     &          (p4+dp(4)).gt.0. )then
                do i2=1,4,1
                    p_b(i1,i2) = p_b(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do   ! Compare the old and new one
        do i1=1,4,1
            if( ABS(throe_p(i1)).gt.1D-15 .AND. 
     &          ABS(throe_p(i1)-p0(i1)) .gt. 1D-15 ) goto 501
        end do

c       Fragments simple (di-)qbar-(di-)q string.
300     continue
        if( n_b.gt.0 .AND. n_q.gt.0 )then   ! 1
            do i=1,n_b,1   ! 2
                do kk=1,5,1
                    P(1,kk) = p_b(i,kk)
                    V(1,kk) = v_b(i,kk)
                end do
                K(1,1) = 2   ! 'A'
                K(1,2) = k_b(i,2)
                K(1,3) = 0
                K(1,4) = 0
                K(1,5) = 0
                do j=1,n_q,1   ! 3
                
                    do kk=1,5,1
                        P(2,kk) = p_q(j,kk)
                        V(2,kk) = v_q(j,kk)
                    end do
                    K(2,1) = 1   ! 'V'
                    K(2,2) = k_q(j,2)
                    K(2,3) = 0
                    K(2,4) = 0
                    K(2,5) = 0
                    N = 2
                    iikk = 0
                    kkii = 0
c                   Sums of px, py, pz, E, inv. m, and charge before the "PYEXEC".
                    ps0=0.
                    do ii=1,6,1
                        ps0(ii)=PYP(0,ii)
                    end do
c                   Sets fragmentation flag in PYEXEC.
                    if(ipden.eq.2 .AND. itden.eq.2)then
                        mstj(1)   = 1
                    else
                        mstp(111) = 1
                    endif
                    
                    call PYEXEC   ! Consider px, py, pz, E   ! in p_23.f
                    
c                   Sums of px, py, pz, E, inv. m, and charge after the "PYEXEC".
                    ps1=0.
                    do ii=1,6,1
                        ps1(ii)=PYP(0,ii)
                    end do
c                   Charge is not conserved or errors occur. Re-generate the string.
                    if( ABS(ps0(6)-ps1(6)).gt.1D-10 .OR.
     &                  iikk.eq.2 .OR. kkii.eq.2 )then
                        N = 0
                        if( j.eq.n_q )then
                            nbe = nbe + 1
                            do kk=1,5,1
                                ! Throws away i-qbar/anti-diq to "sbe".
                                kbe(nbe,kk) = k_b(i,kk)
                                pbe(nbe,kk) = p_b(i,kk)
                                vbe(nbe,kk) = v_b(i,kk)
                                ! Moves one step down.
                                do ii=i+1,n_b,1
                                    k_b(i,kk) = k_b(ii,kk)
                                    p_b(i,kk) = p_b(ii,kk)
                                    v_b(i,kk) = v_b(ii,kk)
                                end do
                            enddo
                            n_b = n_b - 1
                            if(n_b.eq.0) goto 400
                            goto 300
                        end if
                        cycle
                    end if
                    
                    ! Success.
                    do ii=1,4,1
                    throe_p(ii) = throe_p(ii) + ps0(ii) - ps1(ii)
                    end do
                    
                    if(ipden.lt.11) call pyedit(1)   ! in p_23.f
                    if(ipden.ge.11) call pyedit(1)   ! in p_23.f

c                   moves "77" from "PYJETS" to "sgam"
                    if(n.gt.0)then
                        n77 = 0
                        do jj=1,N,1
                            kf = K(jj,2)
                            if(kf.eq.22)then
                                ! '77': photons after hadronization of current string
                                K(jj,2) = 77
                                n77    = n77 + 1
                            endif
                        enddo
                        if(n77.gt.0) call remo_gam_hadro(77)   ! in main.f
                    endif

c                   Dumps produced hadrons (from PYEXEC, no gamma) into "sa1_h".
                    do jj=1,5
                        do ii=1,N,1
                            kn(ii,jj) = k(ii,jj)
                            pn(ii,jj) = p(ii,j)
                            rn(ii,jj) = v(ii,jj)
                        enddo
                    enddo
                    nn  = N
c                   The first position in "sa1_h" is selected from one of partons.
                    i_p = 0
                    if( PYR(1).gt.0.5 ) i_p = 1
                    do jj=1,5
                        rn(1,jj) = v_b(i,jj)
                        if(i_p.eq.1) rn(1,jj) = v_q(j,jj)
                    end do

c                   Arranges produced particles on the surface of sphere (with 
c                    radius rrp and centred on parent position), first produced 
c                    particle is put on its parent position originally.
                    rrp = 1.16
                    rc(1)    = rn(1,1)
                    rc(2)    = rn(1,2)
                    rc(3)    = rn(1,3)
                    rn(1,4)=0.   !TODO(Lei2023060): may be extended.
c                   In corresponding with the time set in >posi".
                    call posi(1,nn,rc,rrp)   ! in sfm.f

c                   Transfers four position messages from "sa1_h" to "PYJETS".
                    do jj=1,4,1
                        do ii=1,nn,1
                            v(ii,jj) = rn(ii,jj)
                        enddo
                    enddo
c                   "PYJETS" to "aff"
                    call tran_pyjets   ! in main.f
                    N = 0
                    ! Moves one step down. (j-q/diq)
                    do kk=1,5
                        do ii=j+1,n_q,1
                            k_q(j,kk) = k_q(ii,kk)
                            p_q(j,kk) = p_q(ii,kk)
                            v_q(j,kk) = v_q(ii,kk)
                        end do
                    end do
                    n_q = n_q - 1
                    ! Moves one step down. (i-qbar/anti-diq)
                    do kk=1,5
                        do ii=i+1,n_b,1
                            k_b(i,kk) = k_b(ii,kk)
                            p_b(i,kk) = p_b(ii,kk)
                            v_b(i,kk) = v_b(ii,kk)
                        end do
                    end do
                    n_b = n_b - 1
                    goto 300
                end do   ! 3
            enddo   ! 2
        endif   ! 1

400     continue

c       Moves rest qbar/anti-diqaruk to 'sbe'. (for the case only one gluon)
        if(n_g.gt.0)then
            do i=1,n_g,1
                nbe = nbe + 1
                do j=1,5
                    kbe(nbe,j) = k_g(i,j)
                    pbe(nbe,j) = p_g(i,j)
                    vbe(nbe,j) = v_g(i,j)
                enddo
            enddo
        endif
        n_g = 0

c       Moves rest qbar/anti-diqaruk to 'sbe'.
        if(nbh.gt.0)then
            do i=1,n_b,1
                nbe = nbe + 1
                do j=1,5
                    kbe(nbe,j) = k_b(i,j)
                    pbe(nbe,j) = p_b(i,j)
                    vbe(nbe,j) = v_b(i,j)
                enddo
            enddo
        endif
        nbh = 0

c       Moves rest q/diqaruk to 'sbe'.
        if(n_q.gt.0)then
            do i=1,n_q,1
                nbe = nbe + 1
                do j=1,5
                    kbe(nbe,j) = k_q(i,j)
                    pbe(nbe,j) = p_q(i,j)
                    vbe(nbe,j) = v_q(i,j)
                enddo
            enddo
        endif
        n_q = 0

c       Gives proper status code, etc.
        if(nbe.gt.0)then
            do i=1,nbe,1
                kbe(i,1) = 2
                kbe(i,3) = 0
                kbe(i,4) = 0
                kbe(i,5) = 0
            end do
        end if
        N = 0


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rest_coal
c00623 Lei2023060
c       Hadronizes/Coalesces "sbe" with "coal".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
c       PYDAT1,PYDAT2,PYDAT3 and PYJETS are the subroutines in PYTHIA.
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/aaff/naff,nonaf,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     &   napp,natt,nzpp,nztt,pio
c       Local variables.
        dimension numb(3)   ! g, qbar, q

        iteration = 0
12345   continue   ! Process the rest partons until empty.
        iteration = iteration + 1  ! No more than 50 times.
        if(nbe.eq.0) goto 54231   ! No rest partons.

c       Shares the 4-momentum in 'throe_p' among partons.
        call share_p_sbe   ! in main.f

c       Dumps "sbe" to "PYJETS".
        N = nbe
        do i2=1,5,1
            do i1=1,nbe,1
                K(i1,i2) = kbe(i1,i2)
                P(i1,i2) = pbe(i1,i2)
                V(i1,i2) = vbe(i1,i2)
            enddo
        enddo
        nbe = 0

        rrp = 1.16
        nn  = 0

        imc   = INT(adj1(13))
        ibc   = INT(adj1(14))
        iphas = INT(adj1(21))

c       Removes potential junctions if sfm in main.
        if(INT(adj1(12)).eq.0 .AND. iteration.eq.1)then
            jb = 0
100         continue
            do i1=jb+1,N,1
                kf   = K(i1,2)
                kfab = iabs(kf)
                if(kfab.ne.88)then
                    jb = jb + 1
                    cycle
                endif
                call updad_pyj(n,i1+1,1)   ! in sfm.f
                N = N - 1
                goto 100
            enddo
        endif

c       Moves gluons from "PYJETS" to "sa36".
        call remo_glu   ! in coales.f

c       Breaks up gluons (with E_g > 2*E_u in "sa36") -> qqbar string 
c        (filling in "PYJETS").
        call break_glu   ! in coales.f

c       Shares 4-momentum in 'throe_p' among partons.
        call share_p_PYJETS   ! in main.f

c       Energetic q (qbar) de-excitation.
        n00   = N   ! Original total entries in PYJETS
        igens = 0
        i_daught_gen = 0   ! the #-th newly produced daughter qqbar
        n_deex = 0   ! number of successful deexcitation
        jb = 0
        n0 = N   ! Current total entries in PYJETS
700     continue
        do i1=jb+1,n0,1
            kf0   = K(i1,2)
            ee    = P(i1,4)
            iflav = 1
            if(kf0.lt.0) iflav = -1
c           iflav = 1 : if source parton is quark
c                 =-1 : if source parton is antiquark
            if(ee.gt.adj1(17))then
                if(i_deex.eq.1) call deexcitation_EP(i1,kf0,igen,iflav)   ! in coales.f
                if(i_deex.eq.2) call deexcitation_E(i1,kf0,igen,iflav)    ! in coales.f
      if(i_deex.eq.3) call deexcitation_EP_comp_pT_1(i1,kf0,igen,iflav)    !Lei2023071
      if(i_deex.eq.4) call deexcitation_EP_comp_pT_2(i1,kf0,igen,iflav)    !Lei2023071
                if(igen.gt.0) n_deex = n_deex + 1
                igens = igens + 1   ! Number of "call deexcitation"
            endif
c           igen : number of generations per source q (qba)
c               Updates n0 and does deexcitation for newly produced qqbar pair
        if(i1.eq.n0 .AND. n.gt.n0 .AND. i_daught_gen.lt.i_deex_gen)then
c         i_deex_gen=0 means no deexcitation for any newly produced qqbar pairs.
c         i_deex_gen=1 means just do deexcitation for the directly proudced qqbar 
c                       pairs (1-st daughters) from original PYJETS (Orig mothers).
c         i_deex_gen=2 means do deexcitation for "1-st daughters" from "Orig mothers" 
c                       and the subsequent qqbar pairs produced from "1-st daughters".
c                       (2-nd daughters).
c         i_deex_gen=3,4,...
c         ...
c         i_deex_gen=999 means always do deexcitation for newly produced qqbar pair
            jb = i1
            i_daught_gen = i_daught_gen + 1
            n0 = n
            goto 700
        end if
        enddo
c       Share the 4-momentum in 'throe_p' among partons.
        call share_p_PYJETS   ! in main.f
c       Energetic q (qbar) de-excitation, finished.

        numb = 0
c       numb(1),(2), and (3): the order # of last g,qba & q in "PYJETS" 
c       1: refers to g (no gluon at all now), 2: qba, 3: q

c       Make the partons in order of qba and q, i.e. move q to the end.
        jh = N
        jl = 0
200     continue
        do j=jl+1,jh,1
            kf   = K(j,2)
            kfab = abs(kf)
            if(kfab.lt.7 .and. kf.gt.0)then   ! q, consider d,u,s,c,b,t
                N = N + 1
                numb(3) = numb(3)+1
                do i4=1,5,1
                    K(n,i4) = K(j,i4)
                    P(n,i4) = P(j,i4)
                    V(n,i4) = V(j,i4)
                enddo
                call updad_pyj(n,j+1,1)   ! in sfm.f
                N  = N  - 1
                jh = jh - 1
                jl = j  - 1
                goto 200
            endif
        enddo
        numb(1) = 0
        numb(2) = N - numb(3)
        numb(3) = n
        n1 = numb(1)
        n2 = numb(2)
        n3 = N

c       Read the table of hadron (meson: pseudoscalar-spin 0 & vector-spin 1 
c        only, baryon: octet-spin 1/2 & decuplet-spin 3/2 only).
        if(iteration.eq.1) call tabhb   ! in coales.f
c       iteration is the calling number

c       Parton coalescence
        netba = 0
        iqba = n2
        call coal(n3,iqba,iii,rrp,iphas,netba)   ! in coales.f
c       n3: total number of partons (qba and q)
c       iqba: total # of qba (qba is ordered before q)
c       ithroq : the total # of quarks thrown away
c       ithrob : the total # of antiquarks thrown away
c       throe : total 4-momentum of the partons thrown away
c       ich : total charge of the partons thrown away

c       Re-coalesces if iphase /= 0.
        iqba = ithrob
        n3   = ithroq + ithrob
        if(iphas.ne.0 .and. n3.ge.2)then
            call coal(n3,iqba,iii,rrp,0,0)
        endif

c       Trasfers "sa1_h" to "PYJETS".
        n = nn
        do j2=1,5,1
            do j1=1,nn,1
                K(j1,j2) = kn(j1,j2)
                P(j1,j2) = pn(j1,j2)
                V(j1,j2) = rn(j1,j2)
            enddo
        enddo

c       Decay of unstable hadrons
        call decayh(rrp)   ! in sfm.f
        nn = 0

c       Moves partons from "PYJETS" to "sbe" after decayh.
        call remop   ! in main.f

c       Appends "PYJETS" to "aaff".
        call tran_pyjets   ! in main.f
        N = 0

c       Appends "sa37" to "sbe".
        if(nth.gt.0)then
            do ii1=1,nth,1
                ii3 = nbe + ii1
                do ii2=1,5,1
                    kbe(ii3,ii2) = kth(ii1,ii2)
                    pbe(ii3,ii2) = pth(ii1,ii2)
                    vbe(ii3,ii2) = vth(ii1,ii2)
                enddo
            enddo
            nbe = nbe + nth
            nth = 0
            ithroq = 0
            ithrob = 0
            ich    = 0
            throe  = 0.
        endif

54231   if(iteration.gt.50) return

        if(nbe.gt.0) goto 12345

        RETURN   !Lei2023060 Do not use the following method.

c       Assume the lost 4-momentum (positive energy) excites one gluon 
c        (off-shell), then try re-coalescence again.
c       The inv. mass should greater than mass of pion (~140 MeV).
        p1  = throe_p(1)
        p2  = throe_p(2)
        p3  = throe_p(3)
        p4  = throe_p(4)
        sm2 = p4**2 - p1**2 - p2**2 - p3**2
        if( p4.gt.0. .AND. nbe.eq.0 .AND. sm2.gt.PYMASS(211)**2 )then
            ! Status.
            kbe(1,1) = 2
            kbe(1,2) = 21
            kbe(1,3) = 0
            kbe(1,4) = 0
            kbe(1,5) = 0
            ! Momentum.
            pbe(1,1) = p1
            pbe(1,2) = p2
            pbe(1,3) = p3
            pbe(1,4) = p4
            pbe(1,5) = 0.
            ! Position. Assume 10*10*10 fm^3 volume. 
            !           (1~10 fm/c of QGP evolution)
            vbe(1,1) = PYR(1)*10.
            vbe(1,2) = PYR(1)*10.
            vbe(1,3) = PYR(1)*10.
            vbe(1,4) = PYR(1)*10.   !TODO(Lei2023060): time of parcas or 0. ?
            vbe(1,5) = 0.
            nbe = 1
            throe_p = 0.
            goto 12345
        end if

        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine share_p_PYJETS
c00623 Lei2023060
c       Share the lost 4-momentum in "throe_p" with particles in "PYJETS".
c       The critirion is inv. m^2 > 0 and E > 0, or original inv. m^2 < 0 
c         but the new one is closer to 0 and E > 0, after sharing.
c       Sometimes there are junctions, which should be excluded.
c       "throe_p" is 4-momentum accumulated before this calling.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        dimension p0(4),dp(4)

        if(N.eq.0) return
c       Potential number of junctions and spectators.
        n_junc = 0
        do i=1,N,1
            if( K(i,2).eq.88 .OR. (P(i,1)**2+P(i,2)**2).le.1D-15 ) 
     &          n_junc = n_junc + 1
        end do
        if(N.eq.n_junc) return

        sn = (N - n_junc)*1.
100     continue
        dp = throe_p/sn
        p0 = throe_p   ! Old
        do i1=1,N,1
            kf = K(i1,2)
            p1 = P(i1,1)
            p2 = P(i1,2)
            p3 = P(i1,3)
            p4 = P(i1,4)
            p5 = P(i1,5)
            pT2= p1**2 + p2**2
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2=(p4+dp(4))**2-(p1+dp(1))**2-(p2+dp(2))**2-(p3+dp(3))**2
c           Inv. m^2 >= 0 or < 0 but closer to 0, E >= m, and not junction or spectator.
            if( (sm2.ge.0. .OR. (sm2.lt.0. .AND. sm2.ge.sm2_org)) .AND. 
     &          (p4+dp(4)).ge.p5 .AND. kf.ne.88 .AND. pT2.gt.1D-15 )then
    !  &          (p4+dp(4)).gt.0. .AND. kf.ne.88 .AND. pT2.gt.1D-15 )then
                do i2=1,4,1
                    P(i1,i2) = P(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do
c       If remaining 4-momentum > sigma (accuracy, 1D-15 here) and sharing 
c        suceeded, try another sharing (iteration).
        do i1=1,4,1   ! Compares the old with new one.
            if( ABS(throe_p(i1)).gt.1D-15 .AND. 
     &          ABS(throe_p(i1)-p0(i1)) .gt. 1D-15 ) goto 100
        end do

        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine share_p_sbe
c00623 Lei2023060
c       Share the lost 4-momentum in "throe_p" with particles in "PYJETS".
c       The critirion is inv. m^2 > 0 and E > 0, or original inv. m^2 < 0 
c         but the new one is closer to 0 and E > 0, after sharing.
c       Since sometimes there are junctions, which should be excluded.
c       "throe_p" is 4-momentum accumulated before this calling.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        common/sbe/nbe,npad,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        dimension p0(4),dp(4)

        if(nbe.eq.0) return
c       Potential number of junctions and spectators.
        n_junc = 0
        do i=1,nbe,1
            if( kbe(i,2).eq.88 .OR. (pbe(i,1)**2+pbe(i,2)**2).le.1D-15 )
     &          n_junc = n_junc + 1
        end do
        if(nbe.eq.n_junc) return

        sn = (nbe - n_junc)*1.
100     continue
        dp = throe_p/sn
        p0 = throe_p   ! Old
        do i1=1,nbe,1
            kf = kbe(i1,2)
            p1 = pbe(i1,1)
            p2 = pbe(i1,2)
            p3 = pbe(i1,3)
            p4 = pbe(i1,4)
            p5 = pbe(i1,5)
            pT2= p1**2 + p2**2
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2=(p4+dp(4))**2-(p1+dp(1))**2-(p2+dp(2))**2-(p3+dp(3))**2
c           Inv. m^2 >= 0 or < 0 but closer to 0, E >= m, and not junction or spectator.
            if( (sm2.ge.0. .OR. (sm2.lt.0. .AND. sm2.ge.sm2_org)) .AND. 
     &          (p4+dp(4)).ge.p5 .AND. kf.ne.88 .AND. pT2.gt.1D-15 )then
    !  &          (p4+dp(4)).gt.0. .AND. kf.ne.88 .AND. pT2.gt.1D-15 )then
                do i2=1,4,1
                    pbe(i1,i2) = pbe(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do
c       If remaining 4-momentum > sigma (accuracy, 1D-15 here) and sharing 
c        suceeded, try another sharing (iteration).
        do i1=1,4,1   ! Compares the old with new one.
            if( ABS(throe_p(i1)).gt.1D-15 .AND. 
     &          ABS(throe_p(i1)-p0(i1)) .gt. 1D-15 ) goto 100
        end do

        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine share_p_PYJETS_sa1h
c00623 Lei2023060
c       Share the lost 4-momentum in "throe_p" with particles in "PYJETS".
c       The critirion is inv. m^2 > 0 and E > 0, or original inv. m^2 < 0 
c         but the new one is closer to 0 and E > 0, after sharing.
c       Sometimes there are junctions, which should be excluded.
c       "throe_p" is 4-momentum accumulated before this calling.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        dimension p0(4),dp(4)

        if( (N+nn).eq.0 ) return
c       Potential number of junctions and spectators.
        n_junc = 0
    !     do i=1,N,1
    !         if( K(i,2).eq.88 .OR. (P(i,1)**2+P(i,2)**2).le.1D-15 ) 
    !  &          n_junc = n_junc + 1
    !     end do
    !     do i=1,nn,1
    !         if( kn(i,2).eq.88 .OR. (pn(i,1)**2+pn(i,2)**2).le.1D-15 ) 
    !  &          n_junc = n_junc + 1
    !     end do
    !     if( (N+nn).eq.n_junc ) return

        sn = (N + nn - n_junc)*1.
100     continue
        dp = throe_p/sn
        p0 = throe_p   ! Old
c       Share with PYJETS.
        do i1=1,N,1
            kf = K(i1,2)
            p1 = P(i1,1)
            p2 = P(i1,2)
            p3 = P(i1,3)
            p4 = P(i1,4)
            p5 = P(i1,5)
            pT2= p1**2 + p2**2
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2=(p4+dp(4))**2-(p1+dp(1))**2-(p2+dp(2))**2-(p3+dp(3))**2
c           Inv. m^2 >= 0 or < 0 but closer to 0, E >= m, and not junction or spectator.
            if( (sm2.ge.0. .OR. (sm2.lt.0. .AND. sm2.ge.sm2_org)) .AND. 
     &          (p4+dp(4)).ge.p5 .AND. kf.ne.88 .AND. pT2.gt.1D-15 )then
    !  &          (p4+dp(4)).gt.0. .AND. kf.ne.88 .AND. pT2.gt.1D-15 )then
                do i2=1,4,1
                    P(i1,i2) = P(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do
c       Share with sa1_h.
        do i1=1,nn,1
            kf = kn(i1,2)
            p1 = pn(i1,1)
            p2 = pn(i1,2)
            p3 = pn(i1,3)
            p4 = pn(i1,4)
            p5 = pn(i1,5)
            pT2= p1**2 + p2**2
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2=(p4+dp(4))**2-(p1+dp(1))**2-(p2+dp(2))**2-(p3+dp(3))**2
c           Inv. m^2 >= 0 or < 0 but closer to 0, E >= m, and not junction or spectator.
            if( (sm2.ge.0. .OR. (sm2.lt.0. .AND. sm2.ge.sm2_org)) .AND. 
     &          (p4+dp(4)).ge.p5 .AND. kf.ne.88 .AND. pT2.gt.1D-15 )then
    !  &          (p4+dp(4)).gt.0. .AND. kf.ne.88 .AND. pT2.gt.1D-15 )then
                do i2=1,4,1
                    pn(i1,i2) = pn(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do
c       If remaining 4-momentum > sigma (accuracy, 1D-15 here) and sharing 
c        suceeded, try another sharing (iteration).
        do i1=1,4,1   ! Compares the old with new one.
            if( ABS(throe_p(i1)).gt.1D-15 .AND. 
     &          ABS(throe_p(i1)-p0(i1)) .gt. 1D-15 ) goto 100
        end do

        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_final_info(win)
c00623 Lei2023060
c       Print the sum of momentum and energy for NN, NA(AN) and AA.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc   ! 260419
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
        dimension p_init(4),p_h(4),p_l(4),p_rest(4),p_lost(4), p_tot(4)

        p_init = 0.
        p_h    = 0.
        p_l    = 0.
        p_rest = 0.
        p_lost = throe_p

c       Initial c & p.
        ic_init = (nzp+nzt)*PYCHGE(2212)
        p_init(4) = win/2.*(nap+nat)

c       Final hadron ( + part of lepton e/mu/tau ), "PYJETS"
        ic_h = 0
        do i=1,N,1
            kf   = K(i,2)
            ic_h = ic_h + PYCHGE(kf)
            do j=1,4,1
                p_h(j) = p_h(j) + P(i,j)
            end do
        enddo

c       Final lepton ( gamma ), "sgam"
        ic_l = 0
        do i=1,ngam,1
            kf   = kgam(i,2)
            ic_l = ic_l + PYCHGE(kf)
            do j=1,4,1
                p_l(j) = p_l(j) + pgam(i,j)
            end do
        enddo

c       Rest parton, "sbe"
        ic_rest = 0
        do i=1,nbe,1
            kf   = kbe(i,2)
            ic_rest = ic_rest + PYCHGE(kf)
            do j=1,4,1
                p_rest(j) = p_rest(j) + pbe(i,j)
            end do
        enddo

c       Total.
        ic_tot = ic_h + ic_l + ic_rest   ! + ic_lost=0
        p_tot  = p_h + p_l + p_rest + p_lost
    !     if( ic_init.ne.ic_tot .OR. ABS(p_tot(1)).gt.1D-10 .OR. 
    !  &      ABS(p_tot(2)).gt.1D-10 .OR. ABS(p_tot(3)).gt.1D-10 .OR.
    !  &      ABS(p_tot(4)-p_init(4)).gt.1D-10 )then
    !         write(777,*) "iii, c & p sum tot=", iii, ic_tot/3., p_tot
    !     end if
    !     if(iii.eq.neve) write(777,*) "Check! iii=", iii
    !     write(888,*) "iii, c & p sum tot=", iii, ic_tot/3., p_tot

c       h + l.
        if( ic_init.ne.(ic_h + ic_l) .OR. ABS(p_h(1) + p_l(1)).gt.1D-10 
     &      .OR. ABS(p_h(2) + p_l(2)).gt.1D-10 
     &      .OR. ABS(p_h(3) + p_l(3)).gt.1D-10 
     &      .OR. ABS(p_h(4) + p_l(4) - p_init(4)).gt.1D-10 )then
        write(2,*) "iii, c & p sum h+l=",iii,(ic_h + ic_l)/3., p_h + p_l
        end if
        if(iii.eq.neve) write(2,*) "Check! iii=", iii
    !    write(77,*) "iii, c & p sum h+l=",iii,(ic_h + ic_l)/3., p_h + p_l

c       Prints information of 3-momentum and energy.
        write(22,*) "------------------------------------------------"//
     &              "------------------------------------------------"//
     &              "------------------------------------------------"
        write(22,*) "Initial              c & p =",
     &                            ic_init / 3., p_init
        write(22,*) "Final   ( h + l )    c & p =",
     &                      (ic_h + ic_l) / 3., (p_h + p_l)
        write(22,*) "Rest parton          c & p =",
     &                            ic_rest / 3., p_rest
        write(22,*) "Final + rest         c & p =",
     &            (ic_h + ic_l + ic_rest) / 3., (p_h + p_l + p_rest)
        write(22,*) "Lost                 c & p =",
     &                                 0. /3., p_lost
        write(22,*) "Final + rest + lost  c & p =", 
     &       (ic_h + ic_l + ic_rest) / 3., (p_h + p_l + p_rest + p_lost)
        write(22,*) "------------------------------------------------"//
     &              "------------------------------------------------"//
     &              "------------------------------------------------"

c       Failure and gamma information.
        if( ABS(ppsa(5)).gt.1D-10 ) write(mstu(11),*)"ppsa=", ppsa
        if( INT(adj1(40)).gt.2 )then
c           Prints the partons those failed in coales.
            if(nth.gt.0)then
        write(22,*) "Summary and list of partons thrown away in coales"
                call prt_sa37(nth,cc)   ! in coales.f
            endif
c           Prints the rest partons those cannot be processed.
            if(nbe.gt.0)then
                write(22,*) "Summary and list of rest partons"
                call prt_sbe(nbe,cc)   ! in parini.f
            endif
c           Prints gamma.
            if(ngam.gt.0)then
                write(22,*)'summary and list of gammas'
                call prt_sgam(ngam,egam,8)
            endif
        endif


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_pyj_sbe_sbh_sgam(win)   !Lei2023060
c       print particle list and sum of momentum and energy
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc   ! 260419
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
        dimension peo(4), p_org(4)

        ich1=0
        peo=0.
        p_org=0.

        ! PYJETS
        do i1=1,N,1
            kf   = K(i1,2)
            ich1 = ich1 + PYCHGE(kf)
            do i2=1,4,1
                peo(i2) = peo(i2) + P(i1,i2)
            end do
        enddo
        ! sbe
        do i1=1,nbe,1
            kf   = kbe(i1,2)
            ich1 = ich1 + PYCHGE(kf)
            do i2=1,4,1
                peo(i2) = peo(i2) + pbe(i1,i2)
            end do
        enddo
        ! sbh
        do i1=1,nbe,1
            kf   = kbh(i1,2)
            ich1 = ich1 + PYCHGE(kf)
            do i2=1,4,1
                peo(i2) = peo(i2) + pbh(i1,i2)
            end do
        enddo
        ! sgam
        do i1=1,ngam,1
            kf   = kgam(i1,2)
            ich1 = ich1 + PYCHGE(kf)
            do i2=1,4,1
                peo(i2) = peo(i2) + pgam(i1,i2)
            end do
        enddo

        cc = ich1/3.

        i_org=(nzp+nzt)*PYCHGE(2212)
        c_org = i_org/3.
        e_org = win/2.*(nap+nat)
    !     if( c_org.ne.cc .OR. ABS(peo(1)).gt.1D-10 .OR. 
    !  &      ABS(peo(2)).gt.1D-10 .OR. ABS(peo(3)).gt.1D-10 .OR.
    !  &      ABS(peo(4)-e_org).gt.1D-10 )then
    !     write(2,*) "iii, c & p sum=", iii, cc, peo
    !     end if
    !     if(iii.eq.neve) write(2,*) "Check! iii=", iii
    !     write(66,*) "iii, c & p sum=", iii, cc, peo


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_pyj_sa1h_sa37!(win)   !Lei2023060
c       print particle list and sum of momentum and energy
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc   ! 260419
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)   ! 150922
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
        dimension peo(4,3), p_org(4), pe(4)

        ich1=0
        ich2=0
        ich3=0
        peo=0.
        p_org=0.
        pe=0.

        ! PYJETS
        write(88,*)
        write(88,*) "iii, N", N
        do i1=1,N,1
            kf   = K(i1,2)
            write(88,*) i1, kf, (p(i1,i3),i3=1,5,1)
            ich1 = ich1 + PYCHGE(kf)
            do i2=1,4,1
                peo(i2,1) = peo(i2,1) + P(i1,i2)
                pe(i2) = pe(i2) + P(i1,i2)
            end do
        enddo
        cc1 = ich1/3.
        write(88,*) "iii, p sum and c", iii, (peo(i3,1),i3=1,4,1), cc1
        ! sa1_h
        write(88,*)
        write(88,*) "iii, nn", nn
        do i1=1,nn,1
            kf   = kn(i1,2)
            write(88,*) i1, kf, (pn(i1,i2),i2=1,5,1)
            ich2 = ich2 + PYCHGE(kf)
            do i2=1,4,1
                peo(i2,2) = peo(i2,2) + pn(i1,i2)
                pe(i2) = pe(i2) + pn(i1,i2)
            end do
        enddo
        cc2 = ich2/3.
        write(88,*) "iii, p sum & c", iii, (peo(i3,2),i3=1,4,1), cc2
        ! sa37
        write(88,*)
        write(88,*) "iii, nth", nth
        do i1=1,nth,1
            kf   = kth(i1,2)
            write(88,*) i1, kf, (pth(i1,i3),i3=1,5,1)
            ich3 = ich3 + PYCHGE(kf)
            do i2=1,4,1
                peo(i2,3) = peo(i2,3) + pth(i1,i2)
                pe(i2) = pe(i2) + pth(i1,i2)
            end do
        enddo
        cc3 = ich3/3.
        !     i_org=(nzp+nzt)*PYCHGE(2212) + (nap+nat-nzp-nzt)*PYCHGE(-2212)
    !     c_org = i_org/3.
    !     e_org = win/2.*(nap+nat)
    !     if( c_org.ne.cc .OR. ABS(peo(1)).gt.1D-10 .OR. 
    !  &      ABS(peo(2)).gt.1D-10 .OR. ABS(peo(3)).gt.1D-10 .OR.
    !  &      ABS(peo(4)-e_org).gt.1D-10 )then
        write(88,*)
        write(88,*) "iii, p sum & c=", iii,  pe, cc1+cc2+cc3
        write(88,*)
    !         end if


        return
        end



c********************************************************************
        subroutine tran_pyjets
c       'pyjets' to 'aff'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
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
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       ??????????????? OSCAR standard output ??????????????????????
        subroutine oscar(win,i_stage)
c140223 The subroutine has been re-wrote.   ! 140223 Lei
c       Records final particle information or event history
c       Obey the OSC1997A/1999A/2013A standard format.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c       common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     &   nap,nat,nzp,nzt,pio
        common/oscar0/ n0, npad0, k0(1000,5), p0(1000,5), v0(1000,5)   !Lei2023060
        common/oscar1/ n1, npad1, k1(kszj,5), p1(kszj,5), v1(kszj,5)   !Lei2023060
        common/oscar2/ n2, npad2, k2(kszj,5), p2(kszj,5), v2(kszj,5)   !Lei2023060
        common/oscar3/ n3, npad3, k3(kszj,5), p3(kszj,5), v3(kszj,5)   !Lei2023060
        character(8) PACIAE, version, frame
c       win : energy
c       i_stage =-1 : prints file header.
c               = 0 : initial nucleons configuration.
c               = 1 : initial parton state (IPS), after parini calling.
c               = 2 : final parton state (FPS), after parcas calling.
c               = 3 : initial hadron state (FPS), after sfm/coal calling.
c               = 4 : final program output, after PACIAE running. (Final state FS)
c                     Normally final hadron state (FHS), after hadcas calling.
c       nosc = 0 : no OSCAR output.
c            = 1 : OSCAR1997A (final_id_p_x, just PACIAE final output)
c            = 2 : OSCAR1999A (full_event_history)
c            = 3 : OSCAR2013A (full_event_history, dummy now)


c*************************************************************************
c****************************** File header ******************************
        if(i_stage.eq.-1)then
            ntest = 1
            PACIAE = "PACIAE"
            version= "2.3"
            if(ifram.eq.0)then
                frame = "lab"
                ebeam = SQRT(win*win+0.938*0.938)
            else if(ifram.eq.1)then
                frame = "nncm"
                ebeam = win*win/2./0.938-1.
            else
                frame = "unknown"
                ebeam = 0D0
                ! stop
            end if
ccccccccccccccccccccccccc   No OSCAR output.   ccccccccccccccccccccccccc
            if(nosc.eq.0)then
cccccccccccccccccccccccccccccc   OSC1997A   cccccccccccccccccccccccccccccc
            else if(nosc.eq.1)then
                write(34,"(A8)")  "OSC1997A"
                write(34,"(A12)") "final_id_p_x"
                write(34,100) PACIAE, version, nap, nzp, nat, nzt, 
     &                        frame, ebeam, ntest
100             format(2(A8,2X),"(", I3, ",", I6, ")+(", I3, ",", 
     &                 I6, ")", 2X, A4, 2X, E10.4, 2X, I8)
cccccccccccccccccccccccccccccc   OSC1999A   cccccccccccccccccccccccccccccc
            else if(nosc.eq.2)then
                write(34,"(A10)") "# OSC1999A"
                write(34,"(A20)") "# full_event_history"
                write(34,"(A12)") "# PACIAE 2.3"
                write(34,200) nap, nzp, nat, nzt, frame, ebeam, ntest
200             format("# (", I3, ",", I6, ")+(", I3, ",", 
     &                 I6, ")", 2X, A4, 2X, E10.4, 2X, I8)
cccccccccccccccccccccccccccccc   OSC2013A   cccccccccccccccccccccccccccccc
            else if(nosc.eq.3)then   !TODO(Lei20230214): need to extend.
                write(34,*)"#!OSCAR2013 full_event_history ID "//
     &                     "px py pz E m x y z t"
                write(34,*)"# Units: none "//
     &                     "GeV GeV GeV GeV GeV fm fm fm fm "
                write(34,"(A12)") "# PACIAE 2.3"
300             continue
            end if
c00623 Initialization. Lei2023060
            n0 = 0
            n1 = 0
            n2 = 0
            n3 = 0
            k0 = 0
            k1 = 0
            k2 = 0
            k3 = 0
            p0 = 0.
            p1 = 0.
            p2 = 0.
            p3 = 0.
            v0 = 0.
            v1 = 0.
            v2 = 0.
            v3 = 0.
c00623 Initialization. Lei2023060
            return
        endif
c****************************** File header ******************************
c*************************************************************************


c*************************************************************************
c****************************** Data dumping *****************************
        if(i_stage.eq.0)then   ! Initial nucleon configuration.
            n0 = N
            do j=1,5,1
                do i=1,N,1
                    k0(i,j) = K(i,j)
                    p0(i,j) = P(i,j)
                    v0(i,j) = V(i,j)
                end do
            end do
        else if(i_stage.eq.1)then   ! Initial parton state (IPS).
            n1 = N
            do j=1,5,1
                do i=1,N,1
                    k1(i,j) = K(i,j)
                    p1(i,j) = P(i,j)
                    v1(i,j) = V(i,j)
                end do
            end do
            if(nbh.gt.0)then
                do j=1,5,1
                    do i=1,nbh,1
                        k1(n1+i,j) = kbh(i,j)
                        p1(n1+i,j) = pbh(i,j)
                        v1(n1+i,j) = vbh(i,j)
                    end do
                end do
                n1 = n1 + nbh
            end if
        else if(i_stage.eq.2)then   ! Final parton state (FPS).
            n2 = N
            do j=1,5,1
                do i=1,N,1
                    k2(i,j) = K(i,j)
                    p2(i,j) = P(i,j)
                    v2(i,j) = V(i,j)
                end do
            end do
            if(nbh.gt.0)then
                do j=1,5,1
                    do i=1,nbh,1
                        k2(n2+i,j) = kbh(i,j)
                        p2(n2+i,j) = pbh(i,j)
                        v2(n2+i,j) = vbh(i,j)
                    end do
                end do
                n2 = n2 + nbh
            end if
        else if(i_stage.eq.3)then   ! Hadronization.
            n3 = N
            do j=1,5,1
                do i=1,N,1
                    k3(i,j) = K(i,j)
                    p3(i,j) = P(i,j)
                    v3(i,j) = V(i,j)
                end do
            end do
        end if
        if(i_stage.lt.4) return
c****************************** Data dumping *****************************
c*************************************************************************


c*************************************************************************
c****************************** Event block ******************************
cccccccccccccccccccccccccccccc   OSC1997A   cccccccccccccccccccccccccccccc
        phi = 0D0
        if(nosc.eq.1)then
            write(34,101) iii, neve, N, bp, phi
            do i=1,N,1
                write(34,102) i, K(i,2), (P(i,j),j=1,5), (V(i,j),j=1,4)
            end do
101         format(I10, 2(2X, I10), 2X, F8.3, 2X, F8.3)
102         format(I10, 2X, I10, 2X, 9(E12.6, 2X))
cccccccccccccccccccccccccccccc   OSC1999A   cccccccccccccccccccccccccccccc
        else if(nosc.eq.2)then
c00623      Prints stage 0, 1, 2, 3 and 4 at once. Lei2023060
            write(34,201) iii, neve, n0, bp, phi, 0   ! Initial nucleons configuration.
            do i=1,n0,1
              write(34,202) i, k0(i,2), (p0(i,j),j=1,5), (v0(i,j),j=1,4)
            end do
            write(34,201) iii, neve, n1, bp, phi, 1   ! IPS
            do i=1,n1,1
              write(34,202) i, k1(i,2), (p1(i,j),j=1,5), (v1(i,j),j=1,4)
            end do
            write(34,201) iii, neve, n2, bp, phi, 2   ! FPS
            do i=1,n2,1
              write(34,202) i, k2(i,2), (p2(i,j),j=1,5), (v2(i,j),j=1,4)
            end do
            write(34,201) iii, neve, n3, bp, phi, 3   ! IHS
            do i=1,n3,1
              write(34,202) i, k3(i,2), (p3(i,j),j=1,5), (v3(i,j),j=1,4)
            end do
c00623 Lei2023060
            write(34,201) iii, neve, N, bp, phi, i_stage   ! FS (FHS)
            do i=1,N,1
                write(34,202) i, K(i,2), (P(i,j),j=1,5), (V(i,j),j=1,4)
            end do
201         format(I10, 2(2X, I10), 2X, F8.3, 2X, F8.3, 2X, I10)
202         format(I10, 2X, I10, 2X, 9(E12.6, 2X))
            if(i_stage.eq.4) write(34,"(11(I10,2X))") 
     &                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
cccccccccccccccccccccccccccccc   OSC2013A   cccccccccccccccccccccccccccccc
        else if(nosc.eq.3)then  !TODO(Lei20230214): need to extend.

        end if


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_gam_sbh(ii)   ! 250209
c       move particles with flavor code ii ('44') from  'sbh' to 'sgam'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        common/sbh/n,nonj,k(kszj,5),p(kszj,5),v(kszj,5)
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
c       move particles with flavor code ii ('55') from 'pyjets' to 'sgam'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
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
c       move particles with flavor code ii ('66') from  'pyjets' to 'sgam'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
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
c       move particles with flavor code ii ('77') from  'pyjets' to 'sgam'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
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
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc   ! 260419
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        dimension peo(4)
        call psum(pgam,1,ngam,peo)   ! in parini.f
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
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/aaff/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
c       do i=1,nn
c       write(22,*)i,kbh(i,2),(pbh(i,j),j=1,4)
c       enddo
        call psum(pbh,1,nbh,peo)   ! in parini.f
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
c        and four positions to the broken partons
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        common/sbe/n,npad,k(kszj,5),p(kszj,5),v(kszj,5)
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
200     k(i1,2)=kf1
        k(n+1,2)=kf2
c221203 k(i1,1)=1
        k(n+1,1)=1
c221203 k(i1,3)=0
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
400     return
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
        common/sbe/n,npad,k(kszj,5),p(kszj,5),v(kszj,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 Lei2023060
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
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

        i_p_split = i_split_diq   !Lei2023060
c       i_p_split = 1: decay method   !Lei2023060
c                 = 2, random 3-momentum method with the different factors
c                 = 3, random 3-momentum method with the same factor for 3-mom
c                 = 4: equal division method

        if( i_p_split.eq.1 ) goto 400   ! activate it for 'decay method'   Lei2023060

c       broken quarks share diquark four momentum randomly,
c        denoted as 'random four momentum method'
c       do i1=1,4   ! activate it for 'random four momentum method'
c       broken quarks share out diquark three momentum randomly,
c        denoted as 'random three momentum method'

c       Random three momentum method with the different factors. Lei2023060
401     do i1=1,3   ! activate it for 'random three momentum method'
        pi(i1)=pyr(1)*p(ii,i1)
        pp(1,i1)=pi(i1)
        pp(2,i1)=ps(i1)-pi(i1)
        enddo

c00623 Lei2023060B--
c       Random three momentum method with the same factor.
        if( i_p_split.eq.3 )then
            factor=PYR(1)
            do i1=1,3,1
                pi(i1)=factor*p(ii,i1)
                pp(1,i1)=pi(i1)
                pp(2,i1)=ps(i1)-pi(i1)
            enddo
        end if

c       Equal division method.
        if( i_p_split.eq.4 )then
            do i1=1,3,1
                pi(i1)=0.5*p(ii,i1)
                pp(1,i1)=pi(i1)
                pp(2,i1)=ps(i1)-pi(i1)
            enddo
        end if
c00623 Lei2023060E--

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
400     continue

c       Decay method.
        if( i_p_split.eq.1 )then   !Lei2023060
        decsuc=1
        call decmom_sbe(ps,pp,am1,am2,decsuc)
        if(decsuc.eq.0)goto 401   ! return to random three momentum method
        endif   !Lei2023060

300     continue
c       adjust four momentum conservation by iteration,no more than
c        4000 iterations
c       call conser(2,pp,ps)   ! in parini.f
        do i1=1,4
        p(ii,i1)=pp(1,i1)
        enddo
        p(ii,5)=am1
        do i1=1,4
        p(n+1,i1)=pp(2,i1)
        enddo
        p(n+1,5)=am2

        do i2=1,4   !Lei2023060 Collects lost 4-momentum.
        throe_p(i2) = throe_p(i2) + ( ps(i2) - p(ii,i2) - p(n+1,i2) )
        enddo


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine decmom_sbe(ps,pp,am1,am2,decsuc)
c       calculate four momentum of decayed particles
c       ps: four momentum of decaying particle
c       am1 and am2: mass of decayed pair
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
        common/sbe/n,npad,k(kszj,5),p(kszj,5),v(kszj,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)   
c       calculate the E and |p| of broken quark in rest frame of diquark
        sm2=ps(4)*ps(4)-ps(1)*ps(1)-ps(2)*ps(2)-ps(3)*ps(3)
c       one problem here is that 'sm2' may not equal to square of diquark 
c        (gluon) rest mass,'bream' is called for spliting g especially
c030603
c1      if(sm2.lt.1.d-10)then
c1      sm2=1.d-10
c1      endif
c030603
        if(sm2.lt.0.005)then   ! 110211
        decsuc=0   ! go back to random three momentum method
        return
        endif
        sm=dsqrt(sm2)   ! M (should be diquark mass)
c       pp(1,4)=(sm2-am2*am2+am1*am1)/2./sm
c       pp(2,4)=(sm2-am1*am1+am2*am2)/2./sm
        ppp=(sm2-(am1+am2)*(am1+am2))*(sm2-(am1-am2)*(am1-am2))
        if(ppp.lt.0.)then   !Lei2023060
        decsuc=0   ! go back to random three momentum method
        return
        endif   !Lei2023060
c161204 ppp=dabs(ppp)   ! 030603 ?
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
c       boost to moving frame of diquark
        ee=ps(4)
        if(ee.lt.1.d-14)ee=1.d-14   ! 021005
        do i1=1,3
        bb(i1)=ps(i1)/ee
        enddo
        call lorntz(1,bb,pi,pj)   ! in parini.f
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
        common/sbe/n,npad,k(kszj,5),p(kszj,5),v(kszj,5)
        dimension rr(3)
        do i1=1,3
c261002 rr(i1)=pyr(1)*v(ii,i1)
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
        subroutine overlap(A1,A2,rnp,rnt,sigma_NN,kjp23,denflag,density0
     c   ,nshot)   ! 020511
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        integer denflag,A1,A2
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200),
     c   part1(200),part2(200),binn(200)   ! 020511 020718
        rmax=10.0
c        write(9,*)'sigma_NN,kjp23,denflag,density0,nshot=',
c     c   sigma_NN,kjp23,denflag,density0,nshot
c       calculate thickness functions for A1 and A2
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
c       Calculate overlap function for A1+A2 collision using densities
        if(kjp23.eq.2)then   ! 190511
        do i=1,200
        bbb=bb(i)
        TA1A2(i)=TAB(A1,A2,bbb,denflag,density0,nshot)
        enddo
        endif
c       Calculate Npart and Nbin using thickness functions
        do i=1,200
        bbb=bb(i)
        if(kjp23.eq.2)then
        call PART(A1,A2,bbb,sigma_NN,denflag,density0,nshot,
     c   part1(i),part2(i))
        binn(i)=(sigma_NN/10.)*TA1A2(i)   ! 020718 180219
c       if(i.eq.36)write(9,*)'i,bbb,pir,tir=',i,bbb,part1(i),part2(i),
c     c   binn(i)   ! 020718
        endif
        if(kjp23.eq.1)then
        call irpt1(bbb,1,a1,a2,rou,rnp,rnt,1,0.d0,pir,tir,1.d0,
     c   1.d0,1.d0,1.d0)
        part1(i)=pir
        part2(i)=tir
c       if(i.eq.36)write(9,*)'i,bbb,pir,tir=',i,bbb,pir,tir
        endif
        enddo
c       Print data
c-----------
        if(kjp23.eq.2)write(9,901)'i','b','TA','TB','TAB','Apart',
     c   'Bpart','Nbin'
        if(kjp23.eq.1)write(9,902)'i','b','Apart','Bpart'
        do i=1,200
        if(kjp23.eq.2)write(9,905)i,bb(i),TA1(i),TA2(i),TA1A2(i),
     c   part1(i),part2(i),binn(i)   ! 020718
        if(kjp23.eq.1)write(9,906)i,bb(i),part1(i),part2(i)
        enddo
 901    format(2a6,6a10)   ! 020718
 902    format(2a6,2a10)
 905    format(i5,2x,f6.2,6(f10.3))   ! 020718
 906    format(i5,2x,f6.2,2(f10.3))
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
     c   part1(200),part2(200),binn(200)   ! 020511 020718
c       write(9,*)'ta a,b,denflag,density0,nshot=',
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
c       write(9,*)'out of ta ta=',ta
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
     c   part1(200),part2(200),binn(200)   ! 020718
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
     c   density(A2,x,y,z1,denflag,density0)
cc      TAB=TAB+o
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
        subroutine PART(A1,A2,b,sigma_NN,denflag,density0,nshot,art1,
     &   art2)! 020511
c       two dimensions trapezoid integral
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        integer denflag,A1,A2
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200),
     c   part1(200),part2(200),binn(200)   ! 020511 020718
c      write(9,*)'part A1,A2,b,sigma_NN=',A1,A2,b,sigma_NN
        art1=0.0
        art2=0.0
cc      do j=1,nshot
cc      x=(pyr(1)-0.5)*2*rmax
cc      y=(pyr(1)-0.5)*2*rmax
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
c       art1=art1+TA1*(1-exp(-sigma_NN*TA2/10))
c       art2=art2+TA2*(1-exp(-sigma_NN*TA1/10))
c       art1=art1+TA1(ib1)*(1-(1-sigma_NN*TA2(ib2)/10/A2)**A2)
c       art2=art2+TA2(ib2)*(1-(1-sigma_NN*TA1(ib1)/10/A1)**A1)
cc      enddo
c       write(9,*)'part1,part2=',part1,part2
cc      art1=art1/nshot*(2*rmax)**2
cc      art2=art2/nshot*(2*rmax)**2
        art1=art1*delta**2
        art2=art2*delta**2
c       write(9,*)'out of part part1,part2=',part1,part2
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shanul(x,q2,rag,raq)   ! 181213
c       calculate nuclear ratio R^A_i=f_{i/A}(x,Q2)/f_i(x,Q2) according
c        to Xin-Nian Wang's paper (PLB 527 (2002) 85), multiply it to
c        the parton distribution function in pythia, resulted parton
c        distribution function is including nuclear shadowing effect
c       it was proved in Eur. Phys. J. C9(1999)61 that nuclear ratio does
c        not depend strongly on the choice for the parton distribution
c        function in nucleon f_i(x,Q2)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
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
c171022      DO 101 I=1,N1   ! 171022 Lei
      DO I=1,N1 
 101  X1(I)=0.5*(U1+D1+(U1-D1)*X11(I))
      END DO
c171022      DO 102 I=1,N2   ! 171022 Lei
      DO I=1,N2 
 102  X2(I)=0.5*(U2+D2+(U2-D2)*X22(I))
      END DO
c171022      DO 103 I=1,N3   ! 171022 Lei
      DO I=1,N3 
 103  X3(I)=0.5*(U3+D3+(U3-D3)*X33(I))
      END DO
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
c171022      DO 105 I1=1,N1   ! 171022 Lei
      DO I1=1,N1
      XX1=X1(I1)
      SUM2=0.
c171022      DO 106 I2=1,N2   ! 171022 Lei
      DO I2=1,N2
      XX2=X2(I2)
      SUM3=0.
c171022      DO 107 I3=1,N3   ! 171022 Lei
      DO I3=1,N3
      XX3=X3(I3)
 107  SUM3=SUM3+W3(I3)*
     $ FUNC(XX1,XX2,XX3,II,IAP,IAT,RP,RT,ROU,BB,AP,BP,AT,BT,CSETA)
      END DO
C
C
C
 106  SUM2=SUM2+W2(I2)*SUM3
      END DO
 105  SUM(II)=SUM(II)+W1(I1)*SUM2 
      END DO
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
C     CALCULATE THE VOLUME OF INTERSECTION FOR PROJECTILE
      ARG1=SQRT(X12+BY2+X32)
      F1=0.
      IF(RP.GE.ARG1)F1=1.
      ARG2=SQRT(X12+X22)
      F2=0.
      IF(RT.GE.ARG2)F2=1.
      GOTO 300
C     CALCULATE THE VOLUME OF INTERSECTION FOR TARGET
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
c171022      DO 105 I = 2, N   ! 171022 Lei                                     
      DO I = 2, N                                                               
  105    W(I) = 0.0E0                                                           
      END DO                                                                    
C                                                                               
      CALL IMTQL2 (N, T, B, W, IERR)                                            
c171022      DO 110 I = 1, N   ! 171022 Lei                                     
      DO I = 1, N                                                               
  110    W(I) = MUZERO * W(I) * W(I)                                            
      END DO                                                                    
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
c171022      DO 10 I = 2, NM1   ! 171022 Lei                                    
      DO I = 2, NM1                                                             
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA                                 
      END DO                                                                    
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
c171022      DO 11 I = 1, NM1   ! 171022 Lei                                    
      DO I = 1, NM1                                                             
         A(I) = 0.0E0                                                           
         ABI = I                                                                
   11    B(I) = ABI/SQRT(4*ABI*ABI - 1.0E0)                                     
      END DO                                                                    
      A(N) = 0.0E0                                                              
      RETURN                                                                    
C                                                                               
C              KIND = 2:  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)          
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)                            
C                                                                               
   20 MUZERO = PI                                                               
c171022      DO 21 I = 1, NM1   ! 171022 Lei                                    
      DO I = 1, NM1                                                             
         A(I) = 0.0E0                                                           
   21    B(I) = 0.5E0                                                           
      END DO                                                                    
      B(1) = SQRT(0.5E0)                                                        
      A(N) = 0.0E0                                                              
      RETURN                                                                    
C                                                                               
C              KIND = 3:  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)         
C              ON (-1, +1), W(X) = SQRT(1 - X*X)                                
C                                                                               
   30 MUZERO = PI/2.0E0                                                         
c171022      DO 31 I = 1, NM1   ! 171022 Lei                                    
      DO I = 1, NM1                                                             
         A(I) = 0.0E0                                                           
   31    B(I) = 0.5E0                                                           
      END DO                                                                    
      A(N) = 0.0E0                                                              
      RETURN                                                                    
C                                                                               
C              KIND = 4:  HERMITE POLYNOMIALS H(X) ON (-INFINITY,               
C              +INFINITY), W(X) = EXP(-X**2)                                    
C                                                                               
   40 MUZERO = SQRT(PI)                                                         
c171022      DO 41 I = 1, NM1   ! 171022 Lei                                    
      DO I = 1, NM1                                                             
         A(I) = 0.0E0                                                           
   41    B(I) = SQRT(I/2.0E0)                                                   
      END DO                                                                    
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
c171022      DO 51 I = 2, NM1   ! 171022 Lei                                    
      DO I = 2, NM1                                                             
         ABI = 2.0E0*I + AB                                                     
         A(I) = A2B2/((ABI - 2.0E0)*ABI)                                        
   51    B(I) = SQRT (4.0E0*I*(I + ALPHA)*(I + BETA)*(I + AB)/                  
     1   ((ABI*ABI - 1)*ABI*ABI))                                               
      END DO                                                                    
      ABI = 2.0E0*N + AB                                                        
      A(N) = A2B2/((ABI - 2.0E0)*ABI)                                           
      RETURN                                                                    
C                                                                               
C              KIND = 6:  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON                   
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER         
C              THAN -1.                                                         
C                                                                               
   60 MUZERO = GAMMA(ALPHA + 1.0E0)   ! 081010                                  
c171022      DO 61 I = 1, NM1   ! 171022 Lei                                    
      DO I = 1, NM1                                                             
         A(I) = 2.0E0*I - 1.0E0 + ALPHA                                         
   61    B(I) = SQRT(I*(I + ALPHA))                                             
      END DO                                                                    
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
c171022      DO 200 II=1,MML   ! 171022 Lei                                     
      DO II=1,MML                                                               
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
      END DO                                                                    
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
c171022      IF(Y-2.0E0)1,2,2   ! 171022 Lei                                    
      IF((Y-2.0E0).LT.0E0) GOTO 1                                               
c171022 2    IF(Y-3.0E0)3,3,4   ! 171022 Lei                                    
 2    IF((Y-3.0E0).LE.0E0) GOTO 3                                               
c171022 4    IF(Y-10.0E0)5,5,6   ! 171022 Lei                                   
 4    IF((Y-10.0E0).GT.0E0) GOTO 6                                              
 5    Y=Y-1.0E0                                                                 
      C=C*Y                                                                     
c171022      IF(Y-3.0E0)3,3,5   ! 171022 Lei                                    
      IF((Y-3.0E0).LE.0E0)THEN                                                  
         GOTO 3                                                                 
      ELSE                                                                      
         GOTO 5                                                                 
      END IF                                                                    
 1    C=Y*C                                                                     
      Y=Y+1.0E0                                                                 
c171022      IF(Y-2.0E0)1,9,9   ! 171022 Lei                                    
      IF((Y-2.0E0).LT.0E0) GOTO 1                                               
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
