        subroutine coales(ijk,neve,nnout,nap,nat,nzp,nzt)
c       A phenomenological coalescence model writen by Sa Ben-Hao on 04/06/2004
c       Its input messages are in 'pyjets'   ! 220822
c       Its storing array is 'pyjets'
c       Its output message is in 'sa1_h' (in 'pyjets' either)
c       ijk: the run number
c       neve: total number of runs
c       nnout: a internal printing per nnout runs
c       nap and nzp: atomic and charge number of projectile
c       nat and nzt: atomic and charge number of target
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        PARAMETER (MPLIS=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c       Those variables in above four statements are only used here and
c        in subroutine 'decayh','findb' and 'thephi'.
c       PYDAT1,PYDAT2,PYDAT3 and PYJETS are the subroutines in PYTHIA
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 Lei2023060
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) ! 220822
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)   ! 150922
        common/sbh/nbh,nonh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   napp,natt,nzpp,nztt,pio
        dimension rc(3),b(3),pstr(3),rstr(3),numb(3),jk(20)
        dimension p0(4),pf(20,4),ppt(50,2),peo(5)
        dimension ppr(kszj,4),peoh(5),pnn(kszj,5),rr(3)
        dimension spglu(4) ! 080512, e (p) sum of unsplitted gluons

        rrp=1.16
        nn=0
        if(ijk.eq.1)then   !Lei2023060
            kn=0
            pn=0.
            rn=0.
        endif   !Lei2023060

        nout=nnout
        imc=INT(adj1(13))
        ibc=INT(adj1(14))
        iphas=INT(adj1(21))

        if(INT(adj1(40)).eq.3)then   ! 070223
c220822     remove junctions
            jb=0
2010        do i1=jb+1,n  ! i1 loop
                kf=k(i1,2)
                kfab=iabs(kf)
                if(kfab.ne.88)then
                    jb=jb+1
                    goto 2020
                endif
                call updad_pyj(n,i1+1,1)   ! 090922 'updad_pyj' in sfm_23.f
                n=n-1
                goto 2010
2020        enddo   ! i1 loop
        endif   ! 070223

c       The baryon number conservation processing is unnecessary in fact.   !Lei2023060
c       It will be satisfied automatically according to the number of q and qbar.
c       Conservation of net baryon.
        netba=0   !Lei2023060
        do i1=1,nbh,1 ! note: 'sbh' is hadron list before hadronization 060119
            kf=kbh(i1,2)
            kfab=iabs(kf)
            if(kf.gt.0.and.(kf.gt.1000 .and. kf.lt.10000)) netba=netba+1
        if(kf.lt.0.and.(kfab.gt.1000 .and. kfab.lt.10000)) netba=netba-1
        enddo
c060813 120214
        if(nap.gt.1 .and. nat.gt.1)then
            netba=nap+nat-netba   ! AA
        elseif(nap.eq.1.and.nat.gt.1.and.ipden.eq.0)then
            netba=nzp+nat-netba   ! p-A (pbar-A)
        elseif(nap.gt.1.and.nat.eq.1)then
            netba=nap+nzt-netba   ! A-p (A-pbar)
        elseif(nap.eq.1.and.nat.eq.1)then
            netba=nzp+nzt-netba   ! pp (p-pbar,pbar-p)
        else
            netba=nat-netba   ! lepton-A
        endif
c060813 120214
c060119 if(ijk.eq.1)then
c       napt=nap+nat
c       nzpt=nzp+nzt
c       sbaryi=float(napt)
c       schgei=3.*float(nzpt)
c060119 endif

888     continue

        if( INT(adj1(12)).eq.2 ) goto 1000   !Lei2023060 g splitting and deexc. bf. "parcas".

c220122
c       move gluons from 'pyjest' to 'sa36'
        call remo_glu
c       break-up gluon (with E_g>2E_u in 'sa36') -> qqbar string 
c        (filling in 'pyjets')
        call break_glu
c       so far, the parton list ('pyjets') is composed of q and qba only

        adj12=adj1(12)
        adj16=adj1(16)
        adj17=adj1(17)
c200222 adj17=max(4.0,adj17) ! 070612, yan

c00623 Shares 4-momentum in "throe_p" among partons.   !Lei2023060
        call share_p_PYJETS   !Lei2023060

        if( INT(adj1(12)).eq.3 ) goto 1000   !Lei2023070 w/ g splitting and w/o deexc.

c280822 energetic q (qbar) de-excitation
        n00=n   ! Original total entries in PYJETS
        igens=0
        i_daught_gen=0   ! the #-th newly produced daughter qqbar
        n_deex=0   ! the # of successful deexcitation
        jb=0
        n0=n   ! Current total entries in PYJETS
700     continue
        do i1=jb+1,n0,1
            kf0=k(i1,2)
            ee=p(i1,4)
            iflav=1
            if(kf0.lt.0)iflav=-1
c           iflav = 1 : if source parton is quark
c                 =-1 : if source parton is antiquark
            if(ee.gt.adj17)then
                if(i_deex.eq.1) call deexcitation_EP(i1,kf0,igen,iflav)   !Lei2023060
                if(i_deex.eq.2) call deexcitation_E(i1,kf0,igen,iflav)    !Lei2023060
      if(i_deex.eq.3) call deexcitation_EP_comp_pT_1(i1,kf0,igen,iflav)    !Lei2023071
      if(i_deex.eq.4) call deexcitation_EP_comp_pT_2(i1,kf0,igen,iflav)    !Lei2023071
                if(igen.gt.0) n_deex = n_deex + 1   !Lei2023060
                igens=igens+1   ! # of 'call deexcitation'
            endif
c           igen : # of generations per source q (qba)
c00623 Lei2023060B---
c           Updates n0 and does deexcitation for newly produced qqbar pair
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
c00623 Lei2023060E---
800     enddo   ! 280822 continue->enddo
900     continue
c       energetic q (qbar) de-excitation, finished.

c00623 Shares the 4-momentum in 'throe_p' among partons.   !Lei2023060
        call share_p_PYJETS   !Lei2023060

c220122

1000    continue   !Lei2023060 For adj12 > 1

        numb(1)=0   ! 220822
        numb(2)=0
        numb(3)=0
c       numb(1),(2), and (3): the order # of last g,qba & q in 'pyjets'
c220822 1: refers to g (no gluon at all now), 2: qba, 3: q

c       make the partons in order of qba and q   ! 220822
c        i.e. move q to the end   ! 220822
        jh=n
        jl=0
2030    continue
        do j=jl+1,jh
        kf=k(j,2)
        kfab=iabs(kf)
        if(kfab.lt.7 .and. kf.gt.0)then   ! q, consider d, u, s, c, b, t
        n=n+1
        numb(3)=numb(3)+1
        do i4=1,5   !Lei2023060 Lei 4 -> 5
        k(n,i4)=k(j,i4)
        p(n,i4)=p(j,i4)
        v(n,i4)=v(j,i4)
        enddo
        call updad_pyj(n,j+1,1)
        n=n-1
        jh=jh-1
        jl=j-1
        goto 2030
        endif
        enddo
        numb(1)=0
        numb(2)=n-numb(3)
        numb(3)=n
        n1=numb(1)
        n2=numb(2)
        n3=n

c       Order the qba according to energy from the maximal to minimal.
c00623 call eord(n1+1,n2)   !Lei2023060
c       Order the q according to energy from the maximal to minimal.
c00623 call eord(n2+1,n3)   !Lei2023060

c       Parton coalescence
        if(ijk.eq.1)call tabhb
c       ijk is the event number
c       Read the table of hadron (meson: pseudoscalar-spin 0 & vector-spin 1 
c        only, baryon: octet-spin 1/2 & decuplet-spin 3/2 only)

        iqba=n2
        call coal(n3,iqba,ijk,rrp,iphas,netba)
c       n3: total number of partons (qba and q)
c       iqba: total # of qba (qba is ordered before q)
c       ithroq : the total # of quarks thrown away
c       ithrob : the total # of antiquarks thrown away
c       throe : total 4-momentum of the partons thrown away
c       ich : total charge of the partons thrown away

        iqba=ithrob   ! not active originally
        n3=ithroq+ithrob   ! not active originally
        if(iphas.ne.0 .and. n3.ge.2)then   !070223
            call coal(n3,iqba,ijk,rrp,0,0)   ! 110905
c070223 re-coalesce rest parton after last 'call coal'
        endif

c150922 ichth=ich   ! 092600

c       'sa1_h' to 'pyjets'.
        n=nn
        do j2=1,5
            do j1=1,nn
                k(j1,j2)=kn(j1,j2)
                p(j1,j2)=pn(j1,j2)
                v(j1,j2)=rn(j1,j2)
            enddo
        enddo

c       Decay of unstable hadrons
        call decayh(rrp)    ! 060119

        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine pes(pei,il,ih,peo)
c       sum up momentum and energy.
c       pei : two dimension array of input momentum and energy
c       il and ih : lower and higher limits of summation
c       peo : one dimension array of output momentum,energy & sqrt(s)
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        PARAMETER (MPLIS=80000)
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
          peo5=abs(peo5)
        endif

100     format('            px           py          pz         e    '//
     c   '  sqrt(s)')
200     format(4x,5(1x,f9.3))
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tabhb
c       Table of primary meson (pseudoscalar (spin 0) and vector (spin 1) 
c        only) and baryon (octet (spin 1/2) and decuplet (spin 3/2) only)
c00623 Added meson: B_s0, B*_s0, B_c+, B*_c+, B*+, and anti-ones, along 
c00623  with eta_b.
c00623 More baryons, especially the heavy-flavor.
c00623 imc: 26 -> 30; ibc: 18 -> 45; total 200 hadrons.   !Lei2023060
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
c       imc (ibc): dimension of meson (baryon) table considered
c       kqh(i,j): flavor code of j-th constituent q (or qba) in i-th quark 
c        configuration of meson
c       kfh(i,1): flavor code of pseudoscalar meson with i-th quark
c        configuration of meson
c       kfh(i,2): flavor code of vector meson with i-th quark
c        configuration of meson
c       proh(i,1): probability of pseudoscalar meson with i-th quark
c        configuration of meson
c       proh(i,2): probability of vector meson with i-th quark
c        configuration of meson
c       amash(i,1): mass of pseudoscalar meson with i-th quark
c        configuration of meson
c       amash(i,2): mass of vector meson with i-th quark
c        configuration of meson
c       kqb(i,j): flavor code of j-th constituent q (or qba) in i-th quark
c        configuration of baryon
c       kfb(i,1): flavor code of octet baryon with i-th quark
c        configuration of baryon
c       kfb(i,2): flavor code of decuplet baryon with i-th quark
c        configuration of baryon
c       prob(i,1): probability of octet baryon with i-th quark
c        configuration of baryon
c       prob(i,2): probability of decuplet baryon with i-th quark
c        configuration of baryon
c       amasb(i,1): mass of octet baryon with i-th quark
c        configuration of baryon
c       amasb(i,2): mass of decuplet baryon with i-th quark
c        configuration of baryon
c00623 Old table.   !Lei2023060
c        data (kqh(i,1),i=1,80)/2*2,2*1,2*3,3*2,3*1,2*3,4,1,4,2,4,3,4,
c     c       2,5,1,5,5,54*0/
c        data (kqh(i,2),i=1,80)/-1,-3,-2,-3,-2,-1,-2,-2,-2,-1,-1,-1,-3,
c     c       -3,-1,-4,-2,-4,-3,-4,-4,-5,-2,-5,-1,-5,54*0/
c        data (kfh(i,1),i=1,80)/211,321,-211,311,-321,-311,111,221,331,
c     c       111,221,331,221,331,411,-411,421,-421,431,-431,441,521,
c     c       -521,511,-511,553,54*0/
c        data (kfh(i,2),i=1,80)/213,323,-213,313,-323,-313,113,223,0,
c     c       113,223,0,333,0,413,-413,423,-423,433,-433,443,2*0,513,
c     c       -513,0,54*0/
c        data (proh(i,1),i=1,80)/6*1.,0.5,2*0.25,0.5,2*0.25,2*0.5,12*1,
c     c       54*0./
c        data (proh(i,2),i=1,80)/6*1.,2*0.5,0.,2*0.5,0.,1.,0.,7*1,
c     c       2*0,2*1,0,54*0./
c
c        data (kqb(i,1),i=1,80)/5*2,3*1,3,4*2,1,2,1,3,2,62*0/
c        data (kqb(i,2),i=1,80)/3*2,1,3,2*1,2*3,3*1,2,1,3*3,1,62*0/
c        data (kqb(i,3),i=1,80)/2,1,3,1,3,1,5*3,6*4,5,62*0/
c        data (kfb(i,1),i=1,80)/0,2212,3222,2112,3322,0,3112,3312,0,
c     c       3122,3212,4122,4222,4112,4232,4132,4332,5122,62*0/
c        data (kfb(i,2),i=1,80)/2224,2214,3224,2114,3324,1114,3114,3314,
c     c       3334,0,3214,4212,4222,0,4232,4132,2*0,62*0/
c        data (prob(i,1),i=1,80)/0.,4*1,0.,2*1.,0.,2*0.5,7*1.,62*0./
c        data (prob(i,2),i=1,80)/9*1.,0.,1.,2*1.,0.,2*1.,2*0.,62*0./
c00623 Old table.   !Lei2023060

c00623 New table.   !Lei2023060
c***********************************************************************
c       imc: 26 -> 30; ibc: 18 -> 45
c       Meson  : 20 + 20 + 6 + 4 = 50  (meson  + anti- + onium)
c       Baryon : 75 + 75         = 150 (baryon + anti-)
c       Total  : 50 + 150        = 200
c-----------------------------------------------------------------------
c       Onium (no anti-): 6 + 4
c       Light mixing onium: pi0, rho0, eta, omega, eta', phi
c       Heavy onium : eta_c, J/psi, eta_b, Upsilon
c***********************************************************************
c       Meson table.
c-----------------------------------------------------------------------
c       KF code of quark and anti-quark.
        data (kqh(i,1),i=1,80)
     &       /  1,  1,  1,  1,  2,  1,  3,  1,  4,  1,   ! 10
     1          5,  2,  2,  2,  2,  3,  2,  4,  2,  5,   ! 20
     2          3,  3,  3,  4,  3,  5,  4,  4,  5,  5,   ! 30
     3          50*0 /                                   ! 80
        data (kqh(i,2),i=1,80)
     &       / -1, -1, -1, -2, -1, -3, -1, -4, -1, -5,   ! 10
     1         -1, -2, -2, -2, -3, -2, -4, -2, -5, -2,   ! 20
     2         -3, -3, -4, -3, -5, -3, -4, -5, -4, -5,   ! 30
     3          50*0 /                                   ! 80
c-----------------------------------------------------------------------
c       KF code of hadron (meson).
c       Pseudoscalar meson, s = 0.
        data (kfh(i,1),i=1,80)
     &   /  111,  221,  331, -211,  211,  311, -311, -411,  411,  511,   ! 10
     1     -511,  111,  221,  331,  321, -321, -421,  421,  521, -521,   ! 20
     2      221,  331, -431,  431,  531, -531,  441,  541, -541,  551,   ! 30
     3      50*0 /                                                       ! 80
c       Vector meson, s = 1.
        data (kfh(i,2),i=1,80)
     &   /  113,  223,    0, -213,  213,  313, -313, -413,  413,  513,   ! 10
     1     -513,  113,  223,    0,  323, -323, -423,  423,  523, -523,   ! 20
     2      333,    0, -433,  433,  533, -533,  443,  543, -543,  553,   ! 30
     3      50*0 /                                                       ! 80
c-----------------------------------------------------------------------
c       Probability of hadron (meson).
        data (proh(i,1),i=1,80)
     &   /  0.5,0.167,0.333,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   ! 10
     1       1.,  0.5,0.167,0.333,   1.,   1.,   1.,   1.,   1.,   1.,   ! 20
     2    0.667,0.333,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   ! 30
     3      50*0. /                                                      ! 80
        data (proh(i,2),i=1,80)
     &   /  0.5,  0.5,   0.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   ! 10
     1       1.,  0.5,  0.5,   0.,   1.,   1.,   1.,   1.,   1.,   1.,   ! 20
     2       1.,   0.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   ! 30
     3      50*0. /                                                      ! 80
c***********************************************************************
c       Baryon table.
c-----------------------------------------------------------------------
c       KF code of 3-quarks.
        data (kqb(i,1),i=1,80)
     &       /  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   ! 10
     1          1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   ! 20
     2          1,  2,  2,  2,  2,  2,  2,  2,  2,  2,   ! 30
     3          2,  2,  2,  2,  3,  3,  3,  3,  3,  3,   ! 40
     4          3,  4,  4,  4,  5,                       ! 45
     4          0,  0,  0,  0,  0,                       ! 50
     5          30*0 /                                   ! 80
        data (kqb(i,2),i=1,80)
     &       /  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,   ! 10
     1          2,  2,  3,  3,  3,  3,  3,  4,  4,  4,   ! 20
     2          5,  2,  2,  2,  2,  3,  3,  3,  3,  3,   ! 30
     3          4,  4,  4,  5,  3,  3,  3,  4,  4,  4,   ! 40
     4          5,  4,  4,  5,  5,                       ! 45
     4          0,  0,  0,  0,  0,                       ! 50
     5          30*0 /                                   ! 80
        data (kqb(i,3),i=1,80)
     &       /  1,  2,  3,  4,  5,  2,  3,  3,  4,  4,   ! 10
     1          5,  5,  3,  4,  4,  5,  5,  4,  5,  5,   ! 20
     2          5,  2,  3,  4,  5,  3,  4,  4,  5,  5,   ! 30
     3          4,  5,  5,  5,  3,  4,  5,  4,  5,  5,   ! 40
     4          5,  4,  5,  5,  5,                       ! 45
     4          0,  0,  0,  0,  0,                       ! 50
     5          30*0 /                                   ! 80
c-----------------------------------------------------------------------
c       KF code of baryon.
        data (kfb(i,1),i=1,80)
     &   /    0, 2112, 3112, 4112, 5112, 2212, 3122, 3212, 4122, 4212,   ! 10
     1     5122, 5212, 3312, 4132, 4312, 5132, 5312, 4412, 5142, 5412,   ! 20
     2     5512,    0, 3222, 4222, 5222, 3322, 4232, 4322, 5232, 5322,   ! 30
     3     4422, 5242, 5422, 5522,    0, 4332, 5332, 4432, 5342, 5432,   ! 40
     4     5532,    0, 5442, 5542,    0,                                 ! 45
     4         0,   0,    0,    0,    0,                                 ! 50
     5     30*0 /                                                        ! 80
        data (kfb(i,2),i=1,80)
     &   / 1114, 2114, 3114, 4114, 5114, 2214,    0, 3214,    0, 4214,   ! 10
     1        0, 5214, 3314,    0, 4314,    0, 5314, 4414,    0, 5414,   ! 20
     2     5514, 2224, 3224, 4224, 5224, 3324,    0, 4324,    0, 5324,   ! 30
     3     4424,    0, 5424, 5524, 3334, 4334, 5334, 4434,    0, 5434,   ! 40
     4     5534, 4444, 5444, 5544, 5554,                                 ! 45
     4         0,   0,    0,    0,    0,                                 ! 50
     5     30*0 /                                                        ! 80
c-----------------------------------------------------------------------
c       Probability of baryon.
        data (prob(i,1),i=1,80)
     &   /   0.,   1.,   1.,   1.,   1.,   1.,  0.5,  0.5,  0.5,  0.5,   ! 10
     1      0.5,  0.5,   1.,  0.5,  0.5,  0.5,  0.5,   1.,  0.5,  0.5,   ! 20
     2       1.,   0.,   1.,   1.,   1.,   1.,  0.5,  0.5,  0.5,  0.5,   ! 30
     3       1.,  0.5,  0.5,   1.,   0.,   1.,   1.,   1.,  0.5,  0.5,   ! 40
     4       1.,   0.,   1.,   1.,   0.,                                 ! 45
     4       0.,   0.,   0.,   0.,   0.,                                 ! 50
     5      30*0. /                                                      ! 80
        data (prob(i,2),i=1,80)
     &    /  1.,   1.,   1.,   1.,   1.,   1.,   0.,   1.,   0.,   1.,   ! 10
     1       0.,   1.,   1.,   0.,   1.,   0.,   1.,   1.,   0.,   1.,   ! 20
     2       1.,   1.,   1.,   1.,   1.,   1.,   0.,   1.,   0.,   1.,   ! 30
     3       1.,   0.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,   1.,   ! 40
     4       1.,   1.,   1.,   1.,   1.,                                 ! 45
     4       0.,   0.,   0.,   0.,   0.,                                 ! 50
     5      30*0. /                                                      ! 80
c***********************************************************************
c00623 New table.   !Lei2023060

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
        subroutine eord(ni,nc)
c   Order particle set (ni to nc) according to energy
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MPLIS=80000)   ! 280822
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 280822
        dimension rr(4),p1(4)
        do 100 i1=ni,nc
          ii=i1
          if(i1.eq.nc)goto 100
c280822   j=ii
          alar=p(ii,4)
          do 200 i2=i1+1,nc
c280822 communication between i1 and i2 of which the energy is lagest
c        among i1+1, i1+2, ..., nc
            ee=p(i2,4)
            if(ee.gt.alar)then   ! 280822 .ge. -> .gt.
                j=i2
                alar=ee
            endif
200       enddo   ! continue-> enddo 280822
c280822 now, j: order # of particle with largest energy
          kii2=k(ii,2)
          do jj=1,4
          p1(jj)=p(ii,jj)
          rr(jj)=v(ii,jj)
          enddo
          pii5=p(ii,5)

          k(ii,2)=k(j,2)
          do jj=1,4
          p(ii,jj)=p(j,jj)
          v(ii,jj)=v(j,jj)
          enddo
          p(ii,5)=p(j,5)

          k(j,2)=kii2
          do jj=1,4
            p(j,jj)=p1(jj)
            v(j,jj)=rr(jj)
          enddo
          p(j,5)=pii5
100     enddo

        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coal(n1,nqb,ijk,rrp,iphas,netba)
c       Parton coalescence (hadronization)
c       n1 : total # of partons (q & qba only)
c       nqb : total # of qba (qba is ordered before q)
c       ijk : the run number
c       iphas: = 1, complete phase space constraint  !Lei2023060
c              = 2, position constraint only
c              = 3, momentum constraint only
c       netba: number of baryons to be first generated keeping
c              baryon conservation 
c       ithroq : total # of quarks thrown away
c       ithrob : total # of antiquarks thrown away
c       throe : total four momentum of the partons thrown away   ! 090922
c       ich : total charge of the partons thrown away
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MPLIS=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104   Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)   ! 150922
        dimension pc(4),rc(3),iar(3),rcp(3)
        dimension psu(3),peo(5),pnn(kszj,5)
        dimension numb(3)   ! 110905

c00623 Moved from 'coales' to here.   !Lei2023060
        ithroq=0
        ithrob=0
        ich=0
        throe=0.

        nth=0
        if(ijk.eq.1)then   !Lei2023060
            kth=0
            pth=0.
            vth=0.
        endif   !Lei2023060

c       Coalescence
        ibarp=0
        ibarm=0
        imes=0

c       The baryon number conservation processing is unnecessary in fact.   !Lei2023060
c       It will be satisfied automatically according to the number of q and qbar.
c       Generate first 'netba' baryons (if 'netba'>0) or -'netba' antibaryons 
c        (if 'netba'<0) keeping baryon conservation 
        if(netba.gt.0)call barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,1)
        if(netba.lt.0)then
302       continue
          do ii1=1,nqb
            call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,1)
            if(isu.eq.1 .and. ibarm.lt.-netba)goto 303
              if(isu.eq.1 .and. ibarm.eq.-netba)goto 304
              if(isu.eq.0) stop 6666
          enddo
303       goto 302
304       continue
        endif

c       Compose first the antiquark, one by one, into hadron 150922
c       A antiquark can be a component of antibaryon or meson.
        jj=0
300     continue
        jj=jj+1
        if(nqb.eq.0)goto 100   ! composing baryon
        if(nqb.lt.3 .and. n1.eq.nqb)goto 301
c       No q, cannot compose meson, throw away nqb.

        if(nqb.eq.1)ii1=1
        if(nqb.gt.1)ii1=nqb*pyr(1)+1

        kf=k(ii1,2)
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

          if(isu.eq.1 .and. nba.eq.1)then   ! one antibaryon produced
            goto 300
          endif

          if(isu.eq.0)then   ! ii1 can not produce antibryon
            call mespro(n1,nqb,ii1,iphas,nme,imes,rrp,isu)
            if(isu.eq.1 .and. nme.eq.1)then   ! one hadron produced
              goto 300
            endif

            if(isu.eq.0)then   !!! ii1 can not produce meson either, throw away
              nth=nth+1
              ithrob=ithrob+1
              kk=k(ii1,2)
              ich=ich+pychge(kk)
              kth(nth,2)=kk
              pth(nth,5)=p(ii1,5)
              do i2=1,4
                throe(i2)=throe(i2)+p(ii1,i2)
                pth(nth,i2)=p(ii1,i2)
                vth(nth,i2)=v(ii1,i2)
              enddo

c       Move parton list one steps downward since ii1+1
              call updad_pyj(n1,ii1+1,1)
              n=n-1   !Lei2023060
              n1=n1-1
              nqb=nqb-1
              goto 300
            endif
          endif   ! ii1 can not produce antibryon
        elseif(rand.gt.relpr .or. nqb.lt.3)then   !
          call mespro(n1,nqb,ii1,iphas,nme,imes,rrp,isu)
          if(isu.eq.1 .and. nme.eq.1)then   ! one hadron produced
            goto 300
          endif
          if(isu.eq.0)then   !! ii1 can not produce meson
            call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,2)
            if(isu.eq.1 .and. nba.eq.1)then   ! one antibryon produced
              goto 300
            endif
            if(isu.eq.0)then !!! can not produce antibryon either, throw away
              nth=nth+1
              ithrob=ithrob+1
              kk=k(ii1,2)
              ich=ich+pychge(kk)
              kth(nth,2)=kk
              pth(nth,5)=p(ii1,5)
              do i2=1,4
                throe(i2)=throe(i2)+p(ii1,i2)
                pth(nth,i2)=p(ii1,i2)
                vth(nth,i2)=v(ii1,i2)
              enddo

c       Move parton list one steps downward since ii1+1
              call updad_pyj(n1,ii1+1,1)
              n=n-1   !Lei2023060
              n1=n1-1
              nqb=nqb-1
              goto 300
            endif   !!!
          endif   !!
        else   !
        endif  !

c       Throw away those qba remained.
301     continue
        if(nqb.eq.0)goto 100
        do i1=1,nqb
          nth=nth+1   ! 110905
          ithrob=ithrob+1
          kk=k(i1,2)
          ich=ich+pychge(kk)
          kth(nth,2)=kk
          v(nth,5)=v(i1,5)
          do i2=1,4
            throe(i2)=throe(i2)+p(i1,i2)
            pth(nth,i2)=p(i1,i2)
            vth(nth,i2)=v(i1,i2)
          enddo
        enddo

c       Move parton list nqb steps downward since nqb+1
        call updad_pyj(n1,nqb+1,nqb)
        n=n-nqb   !Lei2023060
        n1=n1-nqb
        nqb=0
100     continue

600     format(20(1x,i3))
400     if(n1-nqb.eq.0)goto 505   ! get out of 'coal'
        if(n1-nqb.lt.3)goto 500   ! throw away those quarks

c       Proceed for baryon production
        call barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,0)
        if(n1-nqb.eq.0)goto 505
500     continue

c     Throw away remained q   150922
503     continue
        if(n1-nqb.eq.0)goto 505
        ithroq=ithroq+n1-nqb
        do i1=nqb+1,n1
          nth=nth+1   ! 110905 171122 Lei corrects it from nth=nth+i1
          j=nth
          kk=k(i1,2)
          ich=ich+pychge(kk)
          kth(nth,2)=k(i1,2)
          vth(nth,5)=v(i1,5)
          do i2=1,4
            throe(i2)=throe(i2)+p(i1,i2)
            pth(nth,i2)=p(i1,i2)
            vth(nth,i2)=v(i1,i2)
          enddo
        enddo
        n1=0   ! 110905
505     continue

c150922 Reconstruct parton list for calling 'coal' again
c00623 if(iphas.ne.0)then   !Lei2023060 .eq.1 -> .ne.0
c110905
c     Make the parton list ('sa37') in order of qba and q
c150922
        jh = nth
        jl = 0
        n_q = 0
2030    continue
        do j=jl+1,jh,1
            kf = kth(j,2)
            if( kf.gt.0 .AND. kf.lt.7 )then   ! q
                nth = nth + 1
                n_q = n_q + 1
                do i4=1,5,1
                    kth(nth,i4) = kth(j,i4)
                    pth(nth,i4) = pth(j,i4)
                    vth(nth,i4) = vth(j,i4)
                enddo
c               move particle list 'sa37' one step downward since j+1 to nth
                do jjj=1,5,1
                    do ii2=j+1,nth,1
                        kth(ii2-1,jjj) = kth(ii2,jjj)
                        pth(ii2-1,jjj) = pth(ii2,jjj)
                        vth(ii2-1,jjj) = vth(ii2,jjj)
                    enddo
                enddo
                nth = nth - 1
                jh  = jh - 1
                jl  = j  - 1
                goto 2030
            endif
        enddo
        n_tot  = nth
        i_qbar = n_tot - n_q
c       'sa37' -> 'PYJETS'
        n = nth
        do j=1,5,1
            do i=1,nth,1
                k(i,j) = kth(i,j)
                p(i,j) = pth(i,j)
                v(i,j) = vth(i,j)
            end do
        end do
        ithroq = 0
        ithrob = 0
        ich    = 0
        throe  = 0.
c00623 endif   ! 150922 Lei2023060

c00623 Final try for meson production. No more than 50 times.   !Lei2023060
        if(i_qbar.gt.0 .AND. n_q.gt.0)then
            do i=1,50,1
                ii1 = INT( i_qbar*pyr(1) + 1 )
                call mespro(n_tot,i_qbar,ii1,iphas,nme,imes,rrp,isu)
                    if(isu.eq.1 .and. nme.eq.1)then   ! one hadron produced
                        n_q = n_q - 1
                        if( i_qbar.eq.0 .OR. n_q.eq.0) exit
                    end if
            end do
            if(i.eq.50) write(*,*) "Final mespro failed! iii, n" ,ijk,n
        end if
c       Final try for baryon production. No more than 50 times.
        if(n_q.gt.2)then
            do i=1,50,1
                call barpro(n_tot,i_qbar,iphas,nba,ibarp,netba,rrp,0)
                if(n_tot.eq.i_qbar) exit
            end do
            if(i.eq.50) write(*,*) "Final barpro failed! iii, n" ,ijk,n
        end if
        n_q = n_tot - i_qbar
c       Final try for anti-baryon production. No more than 50 times.
        if(i_qbar.gt.2)then
            do i=1,50,1
                ii1 = INT( i_qbar*pyr(1) + 1 )
        call an_barpro(n_tot,i_qbar,ii1,iphas,nba,ibarm,netba,rrp,isu,2)
                if(isu.eq.1 .and. nba.eq.1)then   ! one hadron produced
                    if( i_qbar.eq.0) exit
                end if
            end do
        if(i.eq.50) write(*,*) "Final anbarpro failed! ijk, n" ,ijk,n
        end if
c       'PYJETS' -> 'sa37'
        if(n.gt.0)then
            do i=1,n,1
                do j=1,5,1
                    kth(i,j) = k(i,j)
                    pth(i,j) = p(i,j)
                    vth(i,j) = v(i,j)
                end do
                if( k(i,2).gt.0 )then
                    ithroq = ithroq + 1
                    ich    = ich    + PYCHGE(k(i,2))
                    throe  = throe  + p(i,4)
                elseif( k(i,2).lt.0 )then
                    ithrob = ithrob + 1
                    ich    = ich    + PYCHGE(k(i,2))
                    throe  = throe  + p(i,4)
                end if
            end do
            nth = n
        else
            nth = 0
        end if
c00623 Lei2023060

        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine findm(kf1,kf2,cm,kfii,amasi,isucc,iflav)
c       Find out the primary meson from mesonic table according to kf1 & kf2
c       cm : invariant mass of kf1 & kf2
c       kfii : flavor code of the primary meson
c       amasi : mass of the primary meson
c       isucc = 1 : success
c             = 0 : fail
c       iflav = 1 : kf1>0,do not need to permute kf1 & kf2
c             = -1 : kf1<0,need to permute kf1 & kf2 (never used now, 241022)
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
     c      goto 506   ! vector

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
500       continue
          isucc=0   ! fail
          return
        endif   ! 1

        ikf(1)=kf1
        ikf(2)=kf2
c       Two body permutation = arrangement (2,2)
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
     c      goto 505   ! vector

c           Proceed for pseudoscalar
            kfii=kfh(i4,1)
            amasi=amash(i4,1)
            proi=proh(i4,1)
            ran1=pyr(1)
            if(ran1.gt.proi)goto 503
            goto 504   ! success

505         kfii=kfh(i4,2)   ! vector
            amasi=amash(i4,2)
            proi=proh(i4,2)
            ran1=pyr(1)
            if(ran1.gt.proi)goto 503
            goto 504

            endif
503         continue
502       continue
501     continue
        isucc=0   ! fail
        return
504     isucc=1   ! success
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine findb(kf0,kf1,kf2,cm,kfii,amasi,isucc,iflav)
c       Find out the primary baryon (antibaryon) from baryonic table 
c        according to kf0,kf1,and kf2,these flavor codes are all > 0
c       cm: invariant mass of kf0,kf1 & kf2
c       kfii : flavor code of the primary baryon
c       amasi : mass of the primary baryon
c       isucc = 1 : success
c       isucc = 0 : fail
c       iflav = 1 : if composing parton is quark
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

c       Three body permutation = arrangement (3,3)
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
                amas1=abs(cm-amas1)
                amas2=amasb(i4,2)
                amas2=abs(cm-amas2)
                if(kfi.eq.if1 .and. kfj.eq.if2 .and. kfk.eq.if3)then ! success
                  if(prob(i4,1).eq.0.and.prob(i4,2).eq.0)goto 107 ! fail and try again
                  if(prob(i4,1).eq.0.and.prob(i4,2).ne.0.)goto 108   ! 3/2
        if((prob(i4,1).ne.0.and.prob(i4,2).ne.0.).and.amas2.le.amas1)
     c    goto 108   ! 3/2
c       Goto 108, for spin 3/2 decuplet.
c       Proceed for spin 1/2 octet.
                  kfii=kfb(i4,1)
                  amasi=amasb(i4,1)
                  if(iflav.eq.-1)then
                    kfii=-kfb(i4,1)
                  endif
                  proi=prob(i4,1)
                  ran1=pyr(1)
                  if(ran1.gt.proi)goto 107 ! fail and try again
                  goto 109   ! success

108               kfii=kfb(i4,2)   ! spin 3/2 decuplet
                  amasi=amasb(i4,2)
                  if(iflav.eq.-1)then
                    kfii=-kfb(i4,2)
                  endif
                  proi=prob(i4,2)
                  ran1=pyr(1)
                  if(ran1.gt.proi)goto 107 ! fail and try again
                  goto 109   ! success
                endif

107           continue
106         continue
105       continue
104     continue
        isucc=0   ! fail
        return
109     isucc=1   ! success

        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,iway)
c       To compose baryon
c       n1 : total # of partons (q & qba)
c       nqb : total # of qba (qba is ordered before q)
c       iphas: = 1, complete phase space constraint  !Lei2023060
c              = 2, position constraint only
c              = 3, momentum constraint only
c       ibarp: statistic number of baryon
c       netba: number of baryons keeping baryon conservation
c       iway: a switch
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MPLIS=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104   Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension rcp(3)
        dpmax=adj1(27)
        drmax=adj1(28)
        nba=0
400     continue

        if(n1-nqb.lt.3)return   ! number of quarks less than three
        do 404 i1=nqb+1,n1-2
          kf1=k(i1,2)   ! 150922
          do 405 i2=i1+1,n1-1
            kf2=k(i2,2)
            do 406 i3=i2+1,n1
              kf3=k(i3,2)

              sume=p(i1,4)+p(i2,4)+p(i3,4)
              sump1=p(i1,1)+p(i2,1)+p(i3,1)
              sump2=p(i1,2)+p(i2,2)+p(i3,2)
              sump3=p(i1,3)+p(i2,3)+p(i3,3)
              cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
              if(cm.gt.1.e20)cm=1.e20
              if(cm.le.0.)goto 406   ! fail 071204
              cm=sqrt(cm)

c       Find out the primary baryon from hadron table due to kf1,kf2 & kf3
              call findb(kf1,kf2,kf3,cm,kfii,amasi,isucc,1)

              if(isucc.eq.0)goto 406   ! fail, and keep on cycle, try again.

c       Proceed for success
c070223
c       Phase space adjudgment
              if(iphas.ne.0)then   ! change from '.eq.1' to '.ne.0'
                call phas(i1,i2,i3,isucc,3,iphas)   ! phase spase adjudgment 140223 Lei added iphas
                if(isucc.eq.0)goto 406   ! fail
              endif   !!1
c070223

c       Proceed for success
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
              if(pnnmm.gt.1.e20)pnnmm=1.e20
              if(pnnmm.le.0.)pnnmm=1.e-20
              pnnn=sqrt(pnnmm)
              pn(nn,4)=pnnn
              ! dele=sume-pnnn
              throe_p(4)=throe_p(4)+sume-pnnn   !Lei2023060

c       Produced hadron is arranged among constituent partons randomly.
              pyrx=pyr(1)
              pyry=pyr(1)
              pyrz=pyr(1)
              rn(nn,1)=pyrx*v(i1,1)+pyry*v(i2,1)+pyrz*v(i3,1)
              rn(nn,2)=pyrx*v(i1,2)+pyry*v(i2,2)+pyrz*v(i3,2)
              rn(nn,3)=pyrx*v(i1,3)+pyry*v(i2,3)+pyrz*v(i3,3)

c       Move parton list one step downward from i3+1 to n1
411           call updad_pyj(n1,i3+1,1)
              n=n-1   !Lei2023060
              n1=n1-1

c       Move parton list one step downward from i2+1 to n1
              call updad_pyj(n1,i2+1,1)
              n=n-1   !Lei2023060
              n1=n1-1

c       Move parton list one step downward from i1+1 to n1
              call updad_pyj(n1,i1+1,1)
              n=n-1   !Lei2023060
              n1=n1-1

c00623 Share the surplus 4-momentum in throe_p.   !Lei2023060
              ! call share_p_PYJETS   !Lei2023060
              call share_p_PYJETS_sa1h   !Lei2023060
              if(iway.eq.1)then
                if(nba.lt.netba)goto 400 ! recycle all the partons remained, do again.
                if(nba.eq.netba)return
              endif
              if(iway.eq.2 .and. nba.eq.1)return
              if(iway.eq.0)goto 400   ! 121204
c   iway=1: when the baryon number equal net baryon, return. Used in creat net baryon
c   iway=2: it can return when there generate one baryon. 
c   iway=0: check all the probability of parton constituent baryon, then return.
406         continue   ! fail
405       continue
404     continue
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine an_barpro(n1,nqb,i1,iphas,nba,ibarm,netba,rrp,isu,
     c   iway)
c       To compose an anti-baryon
c       n1 : total # of partons (q & qba)
c       nqb : total # of qba (qba is ordered before q)
c       i1: antiquark wanted to compose antibaryon
c       iphas: = 1, complete phase space constraint  !Lei2023060
c              = 2, position constraint only
c              = 3, momentum constraint only
c       ibarm: statistic number of anti-baryon
c       -netba: number of anti-baryons keeping baryon conservation
c       isu: =1 success
c            =0 fail
c       iway: a switch
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MPLIS=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 150922
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104   Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension rcp(3)
        dpmax=adj1(27)
        drmax=adj1(28)
        isu=1
        nba=0

        if(nqb.lt.3)goto 100   ! number of qbar less than three,return
        kf1=k(i1,2)
        do 405 i2=1,nqb
          if(i2.eq.i1)goto 405
          kf2=k(i2,2)
          do 406 i3=1,nqb
            if(i3.eq.i1)goto 406
            if(i3.eq.i2)goto 406
            kf3=k(i3,2)

            sume=p(i1,4)+p(i2,4)+p(i3,4)
            sump1=p(i1,1)+p(i2,1)+p(i3,1)
            sump2=p(i1,2)+p(i2,2)+p(i3,2)
            sump3=p(i1,3)+p(i2,3)+p(i3,3)
            cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
            if(cm.gt.1.e20)cm=1.e20
            if(cm.le.0.)goto 406   ! fail 071204
            cm=sqrt(cm)

c       Find out the primary antibaryon from hadron table according to kf1,kf2 
c        & kf3
            call findb(-kf1,-kf2,-kf3,cm,kfii,amasi,isucc,-1)

            if(isucc.eq.0)goto 406   ! fail

c       Proceed for success
c070223
c       Phase space adjudgment
            if(iphas.ne.0)then   ! change from '.eq.1' to '.ne.0'
              call phas(i1,i2,i3,isucc,3,iphas)   ! phase spase adjudgment 140223 Lei added iphas
              if(isucc.eq.0)goto 406   ! fail
            endif   !!1
c070223

c       Proceed for success
            goto 400
406       continue   ! fail

405     continue
        goto 100
400     continue
        ibarm=ibarm+1
        nba=nba+1

c       Give proper variables to the primary antibaryon
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
        if(pnnmm.gt.1.e20)pnnmm=1.e20
        if(pnnmm.le.0.)pnnmm=1.e-20
        pnnn=sqrt(pnnmm)
        pn(nn,4)=pnnn
        ! dele=sume-pnnn
        throe_p(4)=throe_p(4)+sume-pnnn   !Lei2023060

c       Arrange produced particle on the surface of sphere with radius
c        rrp and centered at the center of mass

c       Produced hadron is arranged among contituent partons randomly
        pyrx=pyr(1)
        pyry=pyr(1)
        pyrz=pyr(1)
        rn(nn,1)=pyrx*v(i1,1)+pyry*v(i2,1)+pyrz*v(i3,1)
        rn(nn,2)=pyrx*v(i1,2)+pyry*v(i2,2)+pyrz*v(i3,2)
        rn(nn,3)=pyrx*v(i1,3)+pyry*v(i2,3)+pyrz*v(i3,3)

c       Move parton list one step downward from i3+1 to n1
411     call updad_pyj(n1,i3+1,1)
        if(i1.gt.i3)i1=i1-1
        if(i2.gt.i3)i2=i2-1
        nqb=nqb-1
        n1=n1-1
        n=n-1   !Lei2023060

c       Move parton list one step downward from i2+1 to n1
        call updad_pyj(n1,i2+1,1)
        if(i1.gt.i2)i1=i1-1
        nqb=nqb-1
        n1=n1-1
        n=n-1   !Lei2023060

c       Move parton list one step downward from i1+1 to n1
        call updad_pyj(n1,i1+1,1)
        nqb=nqb-1
        n1=n1-1
        n=n-1   !Lei2023060

c00623 Share the surplus 4-momentum in throe_p.   !Lei2023060
        ! call share_p_PYJETS
        call share_p_PYJETS_sa1h
c   iway=1: creat an antibaryon and return
c   iway=2: creat an antibaryon and a baryon, then return
        if(iway.eq.1 .and. nba.eq.1)return
        if(iway.eq.2 .and. nba.eq.1)then
c       An antibaryon generation must be followed immediately a baryon
c        generation keeping baryon conservation
          call barpro(n1,nqb,iphas,nbaa,ibarp,netba,rrp,2)
          return
        endif

100     isu=0
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mespro(n1,nqb,i1,iphas,nme,imes,rrp,isu)
c       Compose a meson
c       n1 : total # of partons (q & qba)
c       nqb : total # of qba (qba is ordered before q)
c       i1: antiquark wanted to compose a meson
c       iphas: = 1, complete phase space constraint  !Lei2023060
c              = 2, position constraint only
c              = 3, momentum constraint only
c       imes: statistic number of meson
c       isu: =1 success
c            =0 fail
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MPLIS=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 150922
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104   Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension rcp(3)
        dpmax=adj1(27)
        drmax=adj1(28)
        isu=1
        nme=0

        if(n1.eq.nqb)goto 100   ! 300105, no quark, return
        kf1=k(i1,2)
        do 102 i2=nqb+1,n1
        kf2=k(i2,2)
        sume=p(i1,4)+p(i2,4)
        sump1=p(i1,1)+p(i2,1)
        sump2=p(i1,2)+p(i2,2)
        sump3=p(i1,3)+p(i2,3)
        cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
        if(cm.gt.1.e20)cm=1.e20
        if(cm.le.0.)goto 102   ! fail 071204
        cm=sqrt(cm)

c       Find out primary meson from hadronic table according to kf2 & kf1
        call findm(kf2,kf1,cm,kfii,amasi,isucc,1)

        if(isucc.eq.0)goto 102   ! fail
c       Proceed for success
c070223
c       Phase space adjudgment
        if(iphas.ne.0)then   ! change from '.eq.1' to '.ne.0'
        call phas(i1,i2,i3,isucc,2,iphas)   ! phase spase adjudgment 140223 Lei added iphas
        if(isucc.eq.0)goto 102   ! fail
        endif   !
c070223 

c       Proceed for success
        imes=imes+1
        nme=nme+1

c       Give proper variables to the primary meson
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
        if(pnnmm.gt.1.e20)pnnmm=1.e20
        if(pnnmm.le.0.)pnnmm=1.e-20
        pnnn=sqrt(pnnmm)
        pn(nn,4)=pnnn
c090922
c       do j1=1,5
c       pnnnj1=pn(nn,j1)
c       aabbss=abs(pnnnj1)
c       pnnnj1=pnnnj1/aabbss   ! symbol, i. e. (+ or -)1
c       dpmax=dpmax/1000.0
c       if(aabbss.ge.dpmax)then
c       aabbss=dpmax
c       pn(nn,j1)=aabbss*pnnnj1
c       endif
c       enddo
c090922
        ! dele=sume-pnnn
        throe_p(4)=throe_p(4)+sume-pnnn   !Lei2023060

c       Produced hadron is seded in between contituent partons randomly
        pyrx=pyr(1)
        pyry=pyr(1)
        rn(nn,1)=pyrx*v(i1,1)+pyry*v(i2,1)
        rn(nn,2)=pyrx*v(i1,2)+pyry*v(i2,2)
        rn(nn,3)=pyrx*v(i1,3)+pyry*v(i2,3)

c       Move parton list one step downward since i2+1
111     call updad_pyj(n1,i2+1,1)
        n1=n1-1
        n=n-1   !Lei2023060

c       Move parton list one step downward since i1+1
        call updad_pyj(n1,i1+1,1)
        nqb=nqb-1
        n1=n1-1
        n=n-1   !Lei2023060

c00623 Share the surplus 4-momentum in throe_p.   !Lei2023060
        ! call share_p_PYJETS
        call share_p_PYJETS_sa1h
        if(nme.eq.1)return
102     enddo   ! 090922 continue-> enddo
100     isu=0   ! 300105
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conse_c(pp,ps,npl,np,nstep)
c       Adjust four momentum conservation by iteration,no more than
c        5000 iterations
c       pp : four momentum of particle
c       ps : the four momentum should be conserved to
c       npl : order # of the first particle
c       np : order # of last particle
c       nstep : interval of the step
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 140604 060813
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
        subroutine deexcitation_E(ii,kf0,igen,iflav)
c00623 Renamed from orginal 'ffm' and some modifications were made.   !Lei2023060
c       qqbar pair generation according to energy conservation ! 200223
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order # of source quark (or antiquark)
c       kf0 : flavor code of source quark (or antiquark)
c       igen : # of generations per source q (qbar)  ! 280822
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is antiquark (kf0<0)
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj16 = adj1(16)   ! # of allowed generations 280822
        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = n

        do i1=1,3
            rc(i1)  = v(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = p(ii,i2)  ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qbar) 260223
        pm0    = p(ii,5) ! m, mass of source q (qba), maybe 0 from PYTHIA.
        ! pm0    = PYMASS(kf0) ! m, mass of source q (qba)

        igen = 0   ! Count # of generation

        if(e0.lt.0.) return   ! Stop generation 280822

c       qqbar creation from vacuum
100     continue

        ie1 = 0
c       Sample the flavor (mass) of generated q (qbar)   260223
        eg = e0
        call break_f(eg,kf,amq)
        amasi =  2*amq    ! 220222 Mass of created qqbar pair (an object) 280822
        kf1   =  kf
        kf2   = -kf
        if(iflav.eq.1)then
            kf1 = -kf
            kf2 =  kf
        endif

104     continue

        pT2_0 = p0(1)**2 + p0(2)**2
        pT_0  = SQRT(pT2_0)
        px_0  = p0(1)
        py_0  = p0(2)
c       Sample transverse momentum of created qqbar pair by "PAPTDI" in coalles.f .
        call PAPTDI(kf0,kf1,kf2,px_0,py_0,pT_0,px,py,pT,i_pT,i_pT_max,x)
        p1(1) = px   ! of created qqbar pair
        p1(2) = py   ! of created qqbar pair
        ppt   = px*px + py*py   ! pT square of created qqbar pair

c       Sample z (energy fraction of created qqbar pair taking from source 
c        q (qbar) ) by Field-Feymman fragmentation function, etc. 081022 240223
        p_mT2 = amasi*amasi + ppt
        if(i_z.gt.0  .AND. i_z.lt.6)  call PYZDIS(kf1,kf2,p_mT2,z1)
        if(i_z.gt.10 .AND. i_z.lt.16) call funcz(z1,ie1)
        if(i_z.eq.0) z1 = PYR(1)
        if(i_pT.eq.7 .OR. i_pT.eq.8) z1 = x

        e1  = z1 * e0   ! Energy of first generation qqbar pair
        p13 = -(amasi*amasi+ppt) + e1*e1   ! Sqruare p1(3) of created qqbar pair

        if(e1.gt.e0 .OR. p13.lt.0.)then ! Something goes wrong, re-samlpling
        ! if((e1-e0).lt.pm0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
            ie1 = ie1 + 1   ! No more than 100 iterations
            if(ie1.gt.100) goto 106   ! Stop re-samlpling
            goto 104
        endif

        p1(3) = sqrt(p13)
        p1(4) = e1   ! E

c       Fill the generated q & qbar into parton list ('PYJETS') after N.

c       Give four position to generated q & qbar 
c        generated q and qbar are arranged around source parton within 0.5 fm 
c        randomly in each one of the three coordinates and has same fourth 
c        coordinate as source parton.
        do i=1,3
            rr(i)    = pyr(1)*0.5
            v(n+1,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+1,i) = rc(i) - rr(i)
            rr(i)    = pyr(1)*0.5
            v(n+2,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+2,i) = rc(i) - rr(i)
        enddo
        v(n+1,4) = v(ii,4) ! ii: order # of sourse parton
        v(n+2,4) = v(ii,4) ! n+1 (n+2): order # of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

        i_p_split = i_split_qqb
c       i_p_split = 1: decay method
c                 = 2, random 3-momentum method with the different factors
c                 = 3, random 3-momentum method with the same factor for 3-mom
c                 = 4: equal division method

c       Decay method.
        if( i_p_split.eq.1 )then
            decsuc = 1
            call decmom(p1,pdec,amq,amq,decsuc)
c090922 p1: four momentum of decaying particle
c       pdec: four momentum of decayed particles   ! 090922
c       amq: mass of decayed particle
c090922 'decmom': in parini_23.f
            if(decsuc.eq.1)then
                do i1=1,4
                    p(n+1,i1) = pdec(1,i1)
                enddo
                do i1=1,4
                    p(n+2,i1) = pdec(2,i1)
                enddo
                if(p00(3).lt.0.)then  ! For the negative direction p_z
                    p(n+1,3) = -p(n+1,3)
                    p(n+2,3) = -p(n+2,3)
                endif
                goto 300
            endif
        end if

c       Random three momentum method with the different factors.
        do i1=1,3,1
            pii = pyr(1)*p1(i1)
            p(n+1,i1) = pii
            p(n+2,i1) = p1(i1) - pii
        enddo

c       Random three momentum method with the same factor.
        if( i_p_split.eq.3 )then
            factor=PYR(1)
            do i1=1,3,1
                pii = factor*p1(i1)
                p(n+1,i1) = pii
                p(n+2,i1) = p1(i1) - pii
            enddo
        end if

c       Equal division method.
        if( i_p_split.eq.4 )then
            do i1=1,3,1
                p(n+1,i1) = p1(i1) * 0.5
                p(n+2,i1) = p(n+1,i1)
            enddo
        end if

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            p(n+1,3) = -p(n+1,3)
            p(n+2,3) = -p(n+2,3)
        endif

c       Recalculates E.
        pn11 = p(n+1,1)   ! pnn(1,1) 280822
        pn12 = p(n+1,2)   ! pnn(1,2) 280822
        pn13 = p(n+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+1,4) = sqrt(agsq)   ! 280822
        pn21 = p(n+2,1)   ! pnn(2,1) 280822
        pn22 = p(n+2,2)   ! pnn(2,2) 280822
        pn23 = p(n+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - p(n+1,1) - p(n+2,1)
        d_py = p1(2) - p(n+1,2) - p(n+2,2)
        d_pz = p1(3) - p(n+1,3) - p(n+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - p(n+1,3) - p(n+2,3)
        d_pE = p1(4) - p(n+1,4) - p(n+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        p(n+1,5) = amq   ! Mass
        p(n+2,5) = amq   ! Mass
        k(n+1,1) = 2     ! 'A'
        k(n+2,1) = 1     ! 'V'
        k(n+1,2) = kf1   ! KF code
        k(n+2,2) = kf2   ! KF code
        k(n+1,3) = ii    ! Mother
        k(n+2,3) = ii    ! Mother
        k(n+1,4) = 0     ! First daughter
        k(n+2,4) = 0     ! First aughter
        k(n+1,5) = 0     ! Last daugter
        k(n+2,5) = 0     ! Last daugter

c       Give proper variables to the remnant parton.
        e1c = e0 - e1   ! Conservation
        do i3=1,3   ! Three-momentum conservation
            p1c(i3) = p0(i3) - p1(i3)
        enddo
        if(p00(3).lt.0.) p1c(3) = p0(3) + p1(3)
c       p0 refers to original q (qbar), p1 refers to generated qqbar pair, 
c       p1c refers to remnant (original parion after generating qqbr pair)

        p1c(4) = e1c

        kf0 = kf00
        if(kf0.gt.0) iflav =  1
        if(kf0.lt.0) iflav = -1
        e0 = e1c
        do i3=1,4
            p0(i3) = p1c(i3)
        enddo
        igen = igen + 1

        n = n + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if(igen.ge.INT(adj16)) goto 106   ! Stop generation 280822
        if( e0.le.adj17) goto 106   ! Stop generation

        goto 100

106     continue
        if(igen.eq.0)return
c       Update four momentum of the remnant of ii-th source q (qbar) 280822
        do i3=1,4   ! 280822
            p(ii,i3) = p0(i3)
        enddo

c       Re-calculate E for the remenant, let its inv. mass >= 0.
        sm2 = p(ii,4)**2 - p(ii,1)**2 - p(ii,2)**2 - p(ii,3)**2
        if( sm2.lt.0. )then
        p(ii,4) = SQRT(p(ii,5)**2 + p(ii,1)**2 + p(ii,2)**2 +p(ii,3)**2)
        p(ii,4) = p(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - p(ii,4)
        end if


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
        subroutine deexcitation_EP(ii,kf0,igen,iflav)
c00623 Rename from orginal 'ffm' and some modifications were made.   !Lei2023060
c       qqbar pair generation according to light-cone variable !
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order # of source quark (or antiquark)
c       kf0 : flavor code of source quark (or antiquark)
c       igen : # of generations per source q (qbar)  ! 280822
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is antiquark (kf0<0)
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj16 = adj1(16)   ! # of allowed generations 280822
        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = n

        do i1=1,3
            rc(i1)  = v(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = p(ii,i2)   ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qba)   260223
        pm0    = p(ii,5) ! m, mass of source q (qba), maybe 0 from PYTHIA.
        ! pm0    = PYMASS(kf0) ! m, mass of source q (qba)

        if(p00(3).lt.0.) p0(3) = -p0(3)  ! Converts E-p_z as E+p_z
        w0 = e0 + p0(3) ! E+p_z

        igen = 0   ! Count # of generation

        if(w0.lt.0.) return   ! Stop generation

c       qqbar creation from vacuum
100     continue

        ie1 = 0
c       Sample the flavor (mass) of generated  q (qbar).   260223
        eg  = e0
        call break_f(eg,kf,amq)
        amasi = 2*amq    ! 220222 Mass of created qqbar pair (a object) 280822
        kf1   =  kf
        kf2   = -kf
        if(iflav.eq.1)then
            kf1 = -kf
            kf2 =  kf
        endif

104     continue

        pT2_0 = p0(1)**2 + p0(2)**2
        pT_0  = SQRT(pT2_0)
        px_0  = p0(1)
        py_0  = p0(2)
c       Sample transverse momentum of created qqbar pair by "PAPTDI" in coales.f .
        call PAPTDI(kf0,kf1,kf2,px_0,py_0,pT_0,px,py,pT,i_pT,i_pT_max,x)
        p1(1) = px   ! of created qqbar pair
        p1(2) = py   ! of created qqbar pair
        ppt   = px*px + py*py   ! pT square of created qqbar pair

c       Sample z (energy fraction of created qqbar pair taking from source 
c        q (qbar) ) by Field-Feymman fragmentation function, etc. 081022 240223
        p_mT2 = amasi*amasi + ppt
        if(i_z.gt.0  .AND. i_z.lt.6)  call PYZDIS(kf1,kf2,p_mT2,z1)
        if(i_z.gt.10 .AND. i_z.lt.16) call funcz(z1,ie1)
        if(i_z.eq.0) z1 = PYR(1)
        if(i_pT.eq.7 .OR. i_pT.eq.8) z1 = x

        w1    = z1 * w0             ! E1 + p_z1 = z * ( E0 + p_z0 )
        amt12 = amasi*amasi + ppt   ! m_qqbar^2 + p_T^2
        e1    = 0.5*( w1 + amt12/w1 )
        p13   = 0.5*( w1 - amt12/w1 )

        if(e1.gt.e0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
        ! if((e0-e1).lt.pm0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
            ie1 = ie1 + 1   ! No more than 100 iterations
            if(ie1.gt.100) goto 106   ! Stop re-samlpling
            goto 104
        endif

        p1(3) = p13
        p1(4) = e1   ! E

c       Fill the generated q & qbar into parton list ('PYJETS') after N.

c       Give four position to generated q & qbar 
c        generated q and qbar are arranged around source parton within 0.5 fm 
c        randomly in each one of the three coordinates and has same fourth 
c        coordinate as source parton
        do i=1,3
            rr(i)    = pyr(1)*0.5
            v(n+1,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+1,i) = rc(i) - rr(i)
            rr(i)    = pyr(1)*0.5
            v(n+2,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+2,i) = rc(i) - rr(i)
        enddo
        v(n+1,4) = v(ii,4) ! ii: order # of sourse parton
        v(n+2,4) = v(ii,4) ! n+1 (n+2): order # of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

        i_p_split = i_split_qqb
c       i_p_split = 1: decay method
c                 = 2, random 3-momentum method with the different factors
c                 = 3, random 3-momentum method with the same factor for 3-mom
c                 = 4: equal division method

c       Decay method.
        if( i_p_split.eq.1 )then
            decsuc = 1
            call decmom(p1,pdec,amq,amq,decsuc)
c090922 p1: four momentum of decaying particle
c       pdec: four momentum of decayed particles   ! 090922
c       amq: mass of decayed particle
c090922 'decmom': in parini_23.f
            if(decsuc.eq.1)then
                do i1=1,4
                    p(n+1,i1) = pdec(1,i1)
                enddo
                do i1=1,4
                    p(n+2,i1) = pdec(2,i1)
                enddo
                if(p00(3).lt.0.)then  ! for the negative direction p_z
                    p(n+1,3) = -p(n+1,3)
                    p(n+2,3) = -p(n+2,3)
                endif
                goto 300
            endif
        end if

c       Random three momentum method with the different factors.
        do i1=1,3,1
            pii = pyr(1)*p1(i1)
            p(n+1,i1) = pii
            p(n+2,i1) = p1(i1) - pii
        enddo

c       Random three momentum method with the same factor.
        if( i_p_split.eq.3 )then
            factor=PYR(1)
            do i1=1,3,1
                pii = factor*p1(i1)
                p(n+1,i1) = pii
                p(n+2,i1) = p1(i1) - pii
            enddo
        end if

c       Equal division method.
        if( i_p_split.eq.4 )then
            do i1=1,3,1
                p(n+1,i1) = p1(i1) * 0.5
                p(n+2,i1) = p(n+1,i1)
            enddo
        end if

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            p(n+1,3) = -p(n+1,3)
            p(n+2,3) = -p(n+2,3)
        endif

c       Recalculates E.
        pn11 = p(n+1,1)   ! pnn(1,1) 280822
        pn12 = p(n+1,2)   ! pnn(1,2) 280822
        pn13 = p(n+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+1,4) = sqrt(agsq)   ! 280822
        pn21 = p(n+2,1)   ! pnn(2,1) 280822
        pn22 = p(n+2,2)   ! pnn(2,2) 280822
        pn23 = p(n+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - p(n+1,1) - p(n+2,1)
        d_py = p1(2) - p(n+1,2) - p(n+2,2)
        d_pz = p1(3) - p(n+1,3) - p(n+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - p(n+1,3) - p(n+2,3)
        d_pE = p1(4) - p(n+1,4) - p(n+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        p(n+1,5) = amq   ! Mass
        p(n+2,5) = amq   ! Mass
        k(n+1,1) = 2     ! 'A'
        k(n+2,1) = 1     ! 'V'
        k(n+1,2) = kf1   ! KF code
        k(n+2,2) = kf2   ! KF code
        k(n+1,3) = ii    ! Mother
        k(n+2,3) = ii    ! Mother
        k(n+1,4) = 0     ! First daughter
        k(n+2,4) = 0     ! First aughter
        k(n+1,5) = 0     ! Last daugter
        k(n+2,5) = 0     ! Last daugter

c       Give proper variables to the remnant parton.
        w1c = w0 - w1   ! Conservation
        do i3=1,3   ! Three-momentum conservation
            p1c(i3) = p0(i3) - p1(i3)
        enddo
c       p0 refers to original q (qbar), p1 refers to generated qqbar pair, 
c       p1c refers to remnant (original parion after generating qqbr pair)

        e1c = w1c - p1c(3)
        p1c(4) = e1c

        kf0 = kf00
        if(kf0.gt.0) iflav =  1
        if(kf0.lt.0) iflav = -1
        w0 = w1c
        e0 = e1c
        do i3=1,4
            p0(i3) = p1c(i3)
        enddo
        igen = igen + 1

        n = n + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if( igen.ge.INT(adj16) ) goto 106   ! Stop generation 280822
        if( e0.le.adj17 ) goto 106   ! Stop generation
        if( w0.le.0. )    goto 106   ! Stop generation

        goto 100

106     continue
        if(igen.eq.0)return
c       Update four momentum of the remnant of ii-th source q (qbar). 280822
        do i3=1,4   ! 280822
            p(ii,i3) = p0(i3)
        enddo
        if(p00(3).lt.0.) p(ii,3) = -p0(3)   ! For the negative direction p_z.

c       Re-calculates E for the remenant, let its inv. mass >= 0.
        sm2 = p(ii,4)**2 - p(ii,1)**2 - p(ii,2)**2 - p(ii,3)**2
        if( sm2.lt.0. )then
        p(ii,4) = SQRT(p(ii,5)**2 + p(ii,1)**2 + p(ii,2)**2 +p(ii,3)**2)
        p(ii,4) = p(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - p(ii,4)
        end if


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
        subroutine deexcitation_EP_comp_pT_1(ii,kf0,igen,iflav)
c00623 Rename from orginal 'ffm' and some modifications were made.   !Lei2023060
c       qqbar pair generation according to light-cone variable !
c      Assuming local pT compensation, i.e. px(-px) and py(-py) for q(qbar), 
c       vice versa. Sample z for qqbar.
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order # of source quark (or antiquark)
c       kf0 : flavor code of source quark (or antiquark)
c       igen : # of generations per source q (qbar)  ! 280822
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is antiquark (kf0<0)
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj16 = adj1(16)   ! # of allowed generations 280822
        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = n

        do i1=1,3
            rc(i1)  = v(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = p(ii,i2)   ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qba)   260223
        pm0    = p(ii,5) ! m, mass of source q (qba), maybe 0 from PYTHIA.
        ! pm0    = PYMASS(kf0) ! m, mass of source q (qba)

        if(p00(3).lt.0.) p0(3) = -p0(3)  ! Converts E-p_z as E+p_z
        w0 = e0 + p0(3) ! E+p_z

        igen = 0   ! Count # of generation

        if(w0.lt.0.) return   ! Stop generation

c       qqbar creation from vacuum
100     continue

        ie1 = 0
c       Sample the flavor (mass) of generated  q (qbar).   260223
        eg  = e0
        call break_f(eg,kf,amq)
        amasi = 2*amq    ! 220222 Mass of created qqbar pair (a object) 280822
        kf1   =  kf
        kf2   = -kf
        if(iflav.eq.1)then
            kf1 = -kf
            kf2 =  kf
        endif

104     continue

        pT2_0 = p0(1)**2 + p0(2)**2
        pT_0  = SQRT(pT2_0)
        px_0  = p0(1)
        py_0  = p0(2)
c       Sample transverse momentum of created qqbar pair by "PAPTDI" in coales.f .
        call PAPTDI(kf0,kf1,kf2,px_0,py_0,pT_0,px,py,pT,i_pT,i_pT_max,x)
        px_1 =  px   ! Local pT compensation.
        py_1 =  py
        px_2 = -px
        py_2 = -py
        p1(1) = px_1 + px_2   ! of created qqbar pair
        p1(2) = py_1 + py_2   ! of created qqbar pair
        ppt   = p1(1)*p1(1) + p1(2)*p1(2)   ! pT square of created qqbar pair

c       Sample z (energy fraction of created qqbar pair taking from source 
c        q (qbar) ) by Field-Feymman fragmentation function, etc. 081022 240223
        p_mT2 = amasi*amasi + ppt
        if(i_z.gt.0  .AND. i_z.lt.6)  call PYZDIS(kf1,kf2,p_mT2,z1)
        if(i_z.gt.10 .AND. i_z.lt.16) call funcz(z1,ie1)
        if(i_z.eq.0) z1 = PYR(1)
        if(i_pT.eq.7 .OR. i_pT.eq.8) z1 = x

        w1    = z1 * w0             ! E1 + p_z1 = z * ( E0 + p_z0 )
        amt12 = amasi*amasi + ppt   ! m_qqbar^2 + p_T^2
        e1    = 0.5*( w1 + amt12/w1 )
        p13   = 0.5*( w1 - amt12/w1 )

        if(e1.gt.e0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
        ! if((e0-e1).lt.pm0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
            ie1 = ie1 + 1   ! No more than 100 iterations
            if(ie1.gt.100) goto 106   ! Stop re-samlpling
            goto 104
        endif

        p1(3) = p13
        p1(4) = e1   ! E

c       Fill the generated q & qbar into parton list ('PYJETS') after N.

c       Give four position to generated q & qbar 
c        generated q and qbar are arranged around source parton within 0.5 fm 
c        randomly in each one of the three coordinates and has same fourth 
c        coordinate as source parton
        do i=1,3
            rr(i)    = pyr(1)*0.5
            v(n+1,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+1,i) = rc(i) - rr(i)
            rr(i)    = pyr(1)*0.5
            v(n+2,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+2,i) = rc(i) - rr(i)
        enddo
        v(n+1,4) = v(ii,4) ! ii: order # of sourse parton
        v(n+2,4) = v(ii,4) ! n+1 (n+2): order # of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

        i_p_split = i_split_qqb
c       i_p_split = 1: decay method
c                 = 2, random 3-momentum method with the different factors
c                 = 3, random 3-momentum method with the same factor for 3-mom
c                 = 4: equal division method

c       Decay method.
        if( i_p_split.eq.1 )then
            decsuc = 1
            call decmom(p1,pdec,amq,amq,decsuc)
c090922 p1: four momentum of decaying particle
c       pdec: four momentum of decayed particles   ! 090922
c       amq: mass of decayed particle
c090922 'decmom': in parini_23.f
            if(decsuc.eq.1)then
                do i1=1,4
                    p(n+1,i1) = pdec(1,i1)
                enddo
                do i1=1,4
                    p(n+2,i1) = pdec(2,i1)
                enddo
                if(p00(3).lt.0.)then  ! for the negative direction p_z
                    p(n+1,3) = -p(n+1,3)
                    p(n+2,3) = -p(n+2,3)
                endif
                goto 300
            endif
        end if

c       Random three momentum method with the different factors.
        do i1=1,3,1
            pii = pyr(1)*p1(i1)
            p(n+1,i1) = pii
            p(n+2,i1) = p1(i1) - pii
        enddo

c       Random three momentum method with the same factor.
        if( i_p_split.eq.3 )then
            factor=PYR(1)
            do i1=1,3,1
                pii = factor*p1(i1)
                p(n+1,i1) = pii
                p(n+2,i1) = p1(i1) - pii
            enddo
        end if

c       Equal division method.
        if( i_p_split.eq.4 )then
            do i1=1,3,1
                p(n+1,i1) = p1(i1) * 0.5
                p(n+2,i1) = p(n+1,i1)
            enddo
        end if

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            p(n+1,3) = -p(n+1,3)
            p(n+2,3) = -p(n+2,3)
        endif

c       Local pT compensation.
        p(n+1,1) = px_1
        p(n+1,2) = py_1
        p(n+2,1) = px_2
        p(n+2,2) = py_2

c       Recalculates E.
        pn11 = p(n+1,1)   ! pnn(1,1) 280822
        pn12 = p(n+1,2)   ! pnn(1,2) 280822
        pn13 = p(n+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+1,4) = sqrt(agsq)   ! 280822
        pn21 = p(n+2,1)   ! pnn(2,1) 280822
        pn22 = p(n+2,2)   ! pnn(2,2) 280822
        pn23 = p(n+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - p(n+1,1) - p(n+2,1)
        d_py = p1(2) - p(n+1,2) - p(n+2,2)
        d_pz = p1(3) - p(n+1,3) - p(n+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - p(n+1,3) - p(n+2,3)
        d_pE = p1(4) - p(n+1,4) - p(n+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        p(n+1,5) = amq   ! Mass
        p(n+2,5) = amq   ! Mass
        k(n+1,1) = 2     ! 'A'
        k(n+2,1) = 1     ! 'V'
        k(n+1,2) = kf1   ! KF code
        k(n+2,2) = kf2   ! KF code
        k(n+1,3) = ii    ! Mother
        k(n+2,3) = ii    ! Mother
        k(n+1,4) = 0     ! First daughter
        k(n+2,4) = 0     ! First aughter
        k(n+1,5) = 0     ! Last daugter
        k(n+2,5) = 0     ! Last daugter

c       Give proper variables to the remnant parton.
        w1c = w0 - w1   ! Conservation
        do i3=1,3   ! Three-momentum conservation
            p1c(i3) = p0(i3) - p1(i3)
        enddo
c       p0 refers to original q (qbar), p1 refers to generated qqbar pair, 
c       p1c refers to remnant (original parion after generating qqbr pair)

        e1c = w1c - p1c(3)
        p1c(4) = e1c

        kf0 = kf00
        if(kf0.gt.0) iflav =  1
        if(kf0.lt.0) iflav = -1
        w0 = w1c
        e0 = e1c
        do i3=1,4
            p0(i3) = p1c(i3)
        enddo
        igen = igen + 1

        n = n + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if( igen.ge.INT(adj16) ) goto 106   ! Stop generation 280822
        if( e0.le.adj17 ) goto 106   ! Stop generation
        if( w0.le.0. )    goto 106   ! Stop generation

        goto 100

106     continue
        if(igen.eq.0)return
c       Update four momentum of the remnant of ii-th source q (qbar). 280822
        do i3=1,4   ! 280822
            p(ii,i3) = p0(i3)
        enddo
        if(p00(3).lt.0.) p(ii,3) = -p0(3)   ! For the negative direction p_z.

c       Re-calculates E for the remenant, let its inv. mass >= 0.
        sm2 = p(ii,4)**2 - p(ii,1)**2 - p(ii,2)**2 - p(ii,3)**2
        if( sm2.lt.0. )then
        p(ii,4) = SQRT(p(ii,5)**2 + p(ii,1)**2 + p(ii,2)**2 +p(ii,3)**2)
        p(ii,4) = p(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - p(ii,4)
        end if


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
        subroutine deexcitation_EP_comp_pT_2(ii,kf0,igen,iflav)
c00623 Rename from orginal 'ffm' and some modifications were made.   !Lei2023060
c       qqbar pair generation according to light-cone variable !
c      Assuming local pT compensation, i.e. px(-px) and py(-py) for q(qbar), 
c       vice versa. Sample z for q0-q (q0-qbar).
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order # of source quark (or antiquark)
c       kf0 : flavor code of source quark (or antiquark)
c       igen : # of generations per source q (qbar)  ! 280822
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is antiquark (kf0<0)
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj16 = adj1(16)   ! # of allowed generations 280822
        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = n

        do i1=1,3
            rc(i1)  = v(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = p(ii,i2)   ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qba)   260223
        ! pm0    = p(ii,5) ! m, mass of source q (qba), maybe 0 from PYTHIA.
        pm0    = PYMASS(kf0) ! m, mass of source q (qba)

        if(p00(3).lt.0.) p0(3) = -p0(3)  ! Converts E-p_z as E+p_z
        w0 = e0 + p0(3) ! E+p_z

        igen = 0   ! Count # of generation

        if(w0.lt.0.) return   ! Stop generation

c       qqbar creation from vacuum
100     continue

        ie1 = 0
c       Sample the flavor (mass) of generated  q (qbar).   260223
        eg  = e0
        call break_f(eg,kf,amq)
        amasi = 2*amq    ! 220222 Mass of created qqbar pair (a object) 280822
        kf1   =  kf
        kf2   = -kf
        if(iflav.eq.1)then
            kf1 = -kf
            kf2 =  kf
        endif

104     continue

        pT2_0 = p0(1)**2 + p0(2)**2
        pT_0  = SQRT(pT2_0)
        px_0  = p0(1)
        py_0  = p0(2)
c       Sample transverse momentum of created qqbar pair by "PAPTDI" in coales.f .
        call PAPTDI(kf0,kf1,kf2,px_0,py_0,pT_0,px,py,pT,i_pT,i_pT_max,x)
        px_1 =  px   ! Local pT compensation.
        py_1 =  py
        px_2 = -px
        py_2 = -py
        p1(1) = px_1 + px_0   ! of q0-q or q0-qbar pair
        p1(2) = py_1 + py_0   ! of q0-q or q0-qbar pair
        ppt   = p1(1)*p1(1) + p1(2)*p1(2)   ! pT square of q0-q or q0-qbar pair

c       Sample z (energy fraction of q0-q or q0-qbar pair taking from source 
c        q (qbar) ) by Field-Feymman fragmentation function, etc. 081022 240223
        amasi = amq + pm0    ! 220222 Mass of created qqbar pair (a object) 280822
        p_mT2 = amasi*amasi + ppt
        if(i_z.gt.0  .AND. i_z.lt.6)  call PYZDIS(kf1,kf2,p_mT2,z1)
        if(i_z.gt.10 .AND. i_z.lt.16) call funcz(z1,ie1)
        if(i_z.eq.0) z1 = PYR(1)
        if(i_pT.eq.7 .OR. i_pT.eq.8) z1 = x

        w1    = z1 * w0             ! E1 + p_z1 = z * ( E0 + p_z0 )
        amt12 = amasi*amasi + ppt   ! m_qqbar^2 + p_T^2
        e1    = 0.5*( w1 + amt12/w1 )
        p13   = 0.5*( w1 - amt12/w1 )

        if(e1.gt.e0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
        ! if((e0-e1).lt.pm0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
            ie1 = ie1 + 1   ! No more than 100 iterations
            if(ie1.gt.100) goto 106   ! Stop re-samlpling
            goto 104
        endif

        p1(3) = p13
        p1(4) = e1   ! E

c       Fill the generated q & qbar into parton list ('PYJETS') after N.

c       Give four position to generated q & qbar 
c        generated q and qbar are arranged around source parton within 0.5 fm 
c        randomly in each one of the three coordinates and has same fourth 
c        coordinate as source parton
        do i=1,3
            rr(i)    = pyr(1)*0.5
            v(n+1,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+1,i) = rc(i) - rr(i)
            rr(i)    = pyr(1)*0.5
            v(n+2,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+2,i) = rc(i) - rr(i)
        enddo
        v(n+1,4) = v(ii,4) ! ii: order # of sourse parton
        v(n+2,4) = v(ii,4) ! n+1 (n+2): order # of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

        i_p_split = i_split_qqb
c       i_p_split = 1: decay method
c                 = 2, random 3-momentum method with the different factors
c                 = 3, random 3-momentum method with the same factor for 3-mom
c                 = 4: equal division method

c       Decay method.
        if( i_p_split.eq.1 )then
            decsuc = 1
            call decmom(p1,pdec,amq,pm0,decsuc)
c090922 p1: four momentum of decaying particle
c       pdec: four momentum of decayed particles   ! 090922
c       amq: mass of decayed particle
c090922 'decmom': in parini_23.f
            if(decsuc.eq.1)then
                do i1=1,4
                    p(n+1,i1) = pdec(1,i1)
                enddo
                do i1=1,4
                    p(n+2,i1) = pdec(2,i1)
                enddo
                if(p00(3).lt.0.)then  ! for the negative direction p_z
                    p(n+1,3) = -p(n+1,3)
                    p(n+2,3) = -p(n+2,3)
                endif
                goto 300
            endif
        end if

c       Random three momentum method with the different factors.
        do i1=1,3,1
            pii = pyr(1)*p1(i1)
            p(n+1,i1) = pii
            p(n+2,i1) = p1(i1) - pii
        enddo

c       Random three momentum method with the same factor.
        if( i_p_split.eq.3 )then
            factor=PYR(1)
            do i1=1,3,1
                pii = factor*p1(i1)
                p(n+1,i1) = pii
                p(n+2,i1) = p1(i1) - pii
            enddo
        end if

c       Equal division method.
        if( i_p_split.eq.4 )then
            do i1=1,3,1
                p(n+1,i1) = p1(i1) * 0.5
                p(n+2,i1) = p(n+1,i1)
            enddo
        end if

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            p(n+1,3) = -p(n+1,3)
            p(n+2,3) = -p(n+2,3)
        endif

c       Local pT compensation.
        p(n+1,1) = px_1
        p(n+1,2) = py_1
        p(n+2,1) = px_2
        p(n+2,2) = py_2

c       Recalculates E.
        pn11 = p(n+1,1)   ! pnn(1,1) 280822
        pn12 = p(n+1,2)   ! pnn(1,2) 280822
        pn13 = p(n+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+1,4) = sqrt(agsq)   ! 280822
        pn21 = p(n+2,1)   ! pnn(2,1) 280822
        pn22 = p(n+2,2)   ! pnn(2,2) 280822
        pn23 = p(n+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - p(n+1,1) - p(n+2,1)
        d_py = p1(2) - p(n+1,2) - p(n+2,2)
        d_pz = p1(3) - p(n+1,3) - p(n+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - p(n+1,3) - p(n+2,3)
        d_pE = p1(4) - p(n+1,4) - p(n+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        p(n+1,5) = amq   ! Mass
        p(n+2,5) = amq   ! Mass
        k(n+1,1) = 2     ! 'A'
        k(n+2,1) = 1     ! 'V'
        k(n+1,2) = kf1   ! KF code
        k(n+2,2) = kf2   ! KF code
        k(n+1,3) = ii    ! Mother
        k(n+2,3) = ii    ! Mother
        k(n+1,4) = 0     ! First daughter
        k(n+2,4) = 0     ! First aughter
        k(n+1,5) = 0     ! Last daugter
        k(n+2,5) = 0     ! Last daugter

c       Give proper variables to the remnant parton.
        w1c = w0 - w1   ! Conservation
        do i3=1,3   ! Three-momentum conservation
            p1c(i3) = p0(i3) - p1(i3)
        enddo
c       p0 refers to original q (qbar), p1 refers to generated qqbar pair, 
c       p1c refers to remnant (original parion after generating qqbr pair)

        e1c = w1c - p1c(3)
        p1c(4) = e1c

        kf0 = kf00
        if(kf0.gt.0) iflav =  1
        if(kf0.lt.0) iflav = -1
        w0 = w1c
        e0 = e1c
        do i3=1,4
            p0(i3) = p1c(i3)
        enddo
        igen = igen + 1

        n = n + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if( igen.ge.INT(adj16) ) goto 106   ! Stop generation 280822
        if( e0.le.adj17 ) goto 106   ! Stop generation
        if( w0.le.0. )    goto 106   ! Stop generation

        goto 100

106     continue
        if(igen.eq.0)return
c       Update four momentum of the remnant of ii-th source q (qbar). 280822
        do i3=1,4   ! 280822
            p(ii,i3) = p0(i3)
        enddo
        if(p00(3).lt.0.) p(ii,3) = -p0(3)   ! For the negative direction p_z.

c       Re-calculates E for the remenant, let its inv. mass >= 0.
        sm2 = p(ii,4)**2 - p(ii,1)**2 - p(ii,2)**2 - p(ii,3)**2
        if( sm2.lt.0. )then
        p(ii,4) = SQRT(p(ii,5)**2 + p(ii,1)**2 + p(ii,2)**2 +p(ii,3)**2)
        p(ii,4) = p(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - p(ii,4)
        end if


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
        subroutine deexcitation_EP_comp_pT_3(ii,kf0,igen,iflav)
c00623 Rename from orginal 'ffm' and some modifications were made.   !Lei2023060
c       qqbar pair generation according to light-cone variable !
c      Assuming local pT compensation, i.e. px(-px) and py(-py) for q(qbar), 
c       vice versa. Sample z for q0-q (q0-qbar). Using PYKFDI (PYINDF-like).
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order # of source quark (or antiquark)
c       kf0 : flavor code of source quark (or antiquark)
c       igen : # of generations per source q (qbar)  ! 280822
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is antiquark (kf0<0)
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj16 = adj1(16)   ! # of allowed generations 280822
        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = n

        do i1=1,3
            rc(i1)  = v(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = p(ii,i2)   ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qba)   260223
        ! pm0    = p(ii,5) ! m, mass of source q (qba), maybe 0 from PYTHIA.
        pm0    = PYMASS(kf0) ! m, mass of source q (qba)

        if(p00(3).lt.0.) p0(3) = -p0(3)  ! Converts E-p_z as E+p_z
        w0 = e0 + p0(3) ! E+p_z

        igen = 0   ! Count # of generation

        if(w0.lt.0.) return   ! Stop generation

c       qqbar creation from vacuum
100     continue

        ie1 = 0
c       Sample the flavor (mass) of generated  q (qbar).   260223
        eg  = e0
        call break_f(eg,kf,amq)
        amasi = 2*amq    ! 220222 Mass of created qqbar pair (a object) 280822
        kf1   =  kf
        kf2   = -kf
        if(iflav.eq.1)then
            kf1 = -kf
            kf2 =  kf
        endif

104     continue

        pT2_0 = p0(1)**2 + p0(2)**2
        pT_0  = SQRT(pT2_0)
        px_0  = p0(1)
        py_0  = p0(2)
c       Sample transverse momentum of created qqbar pair by "PAPTDI" in coales.f .
        call PAPTDI(kf0,kf1,kf2,px_0,py_0,pT_0,px,py,pT,i_pT,i_pT_max,x)
        px_1 =  px   ! Local pT compensation.
        py_1 =  py
        px_2 = -px
        py_2 = -py
        p1(1) = px_1 + px_2   ! of created qqbar pair
        p1(2) = py_1 + py_2   ! of created qqbar pair
        ppt   = p1(1)*p1(1) + p1(2)*p1(2)   ! pT square of created qqbar pair

c       Sample z (energy fraction of created qqbar pair taking from source 
c        q (qbar) ) by Field-Feymman fragmentation function, etc. 081022 240223
        p_mT2 = amasi*amasi + ppt
        if(i_z.gt.0  .AND. i_z.lt.6)  call PYZDIS(kf1,kf2,p_mT2,z1)
        if(i_z.gt.10 .AND. i_z.lt.16) call funcz(z1,ie1)
        if(i_z.eq.0) z1 = PYR(1)
        if(i_pT.eq.7 .OR. i_pT.eq.8) z1 = x

        w1    = z1 * w0             ! E1 + p_z1 = z * ( E0 + p_z0 )
        amt12 = amasi*amasi + ppt   ! m_qqbar^2 + p_T^2
        e1    = 0.5*( w1 + amt12/w1 )
        p13   = 0.5*( w1 - amt12/w1 )

        if(e1.gt.e0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
        ! if((e0-e1).lt.pm0 .OR. p13.lt.0.)then ! Something goes wrong, re-sampling.
            ie1 = ie1 + 1   ! No more than 100 iterations
            if(ie1.gt.100) goto 106   ! Stop re-samlpling
            goto 104
        endif

        p1(3) = p13
        p1(4) = e1   ! E

c       Fill the generated q & qbar into parton list ('PYJETS') after N.

c       Give four position to generated q & qbar 
c        generated q and qbar are arranged around source parton within 0.5 fm 
c        randomly in each one of the three coordinates and has same fourth 
c        coordinate as source parton
        do i=1,3
            rr(i)    = pyr(1)*0.5
            v(n+1,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+1,i) = rc(i) - rr(i)
            rr(i)    = pyr(1)*0.5
            v(n+2,i) = rc(i) + rr(i)
            if(pyr(1).gt.0.5) v(n+2,i) = rc(i) - rr(i)
        enddo
        v(n+1,4) = v(ii,4) ! ii: order # of sourse parton
        v(n+2,4) = v(ii,4) ! n+1 (n+2): order # of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

        i_p_split = i_split_qqb
c       i_p_split = 1: decay method
c                 = 2, random 3-momentum method with the different factors
c                 = 3, random 3-momentum method with the same factor for 3-mom
c                 = 4: equal division method

c       Decay method.
        if( i_p_split.eq.1 )then
            decsuc = 1
            call decmom(p1,pdec,amq,amq,decsuc)
c090922 p1: four momentum of decaying particle
c       pdec: four momentum of decayed particles   ! 090922
c       amq: mass of decayed particle
c090922 'decmom': in parini_23.f
            if(decsuc.eq.1)then
                do i1=1,4
                    p(n+1,i1) = pdec(1,i1)
                enddo
                do i1=1,4
                    p(n+2,i1) = pdec(2,i1)
                enddo
                if(p00(3).lt.0.)then  ! for the negative direction p_z
                    p(n+1,3) = -p(n+1,3)
                    p(n+2,3) = -p(n+2,3)
                endif
                goto 300
            endif
        end if

c       Random three momentum method with the different factors.
        do i1=1,3,1
            pii = pyr(1)*p1(i1)
            p(n+1,i1) = pii
            p(n+2,i1) = p1(i1) - pii
        enddo

c       Random three momentum method with the same factor.
        if( i_p_split.eq.3 )then
            factor=PYR(1)
            do i1=1,3,1
                pii = factor*p1(i1)
                p(n+1,i1) = pii
                p(n+2,i1) = p1(i1) - pii
            enddo
        end if

c       Equal division method.
        if( i_p_split.eq.4 )then
            do i1=1,3,1
                p(n+1,i1) = p1(i1) * 0.5
                p(n+2,i1) = p(n+1,i1)
            enddo
        end if

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            p(n+1,3) = -p(n+1,3)
            p(n+2,3) = -p(n+2,3)
        endif

c       Local pT compensation.
        p(n+1,1) = px_1
        p(n+1,2) = py_1
        p(n+2,1) = px_2
        p(n+2,2) = py_2

c       Recalculates E.
        pn11 = p(n+1,1)   ! pnn(1,1) 280822
        pn12 = p(n+1,2)   ! pnn(1,2) 280822
        pn13 = p(n+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+1,4) = sqrt(agsq)   ! 280822
        pn21 = p(n+2,1)   ! pnn(2,1) 280822
        pn22 = p(n+2,2)   ! pnn(2,2) 280822
        pn23 = p(n+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        p(n+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - p(n+1,1) - p(n+2,1)
        d_py = p1(2) - p(n+1,2) - p(n+2,2)
        d_pz = p1(3) - p(n+1,3) - p(n+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - p(n+1,3) - p(n+2,3)
        d_pE = p1(4) - p(n+1,4) - p(n+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        p(n+1,5) = amq   ! Mass
        p(n+2,5) = amq   ! Mass
        k(n+1,1) = 2     ! 'A'
        k(n+2,1) = 1     ! 'V'
        k(n+1,2) = kf1   ! KF code
        k(n+2,2) = kf2   ! KF code
        k(n+1,3) = ii    ! Mother
        k(n+2,3) = ii    ! Mother
        k(n+1,4) = 0     ! First daughter
        k(n+2,4) = 0     ! First aughter
        k(n+1,5) = 0     ! Last daugter
        k(n+2,5) = 0     ! Last daugter

c       Give proper variables to the remnant parton.
        w1c = w0 - w1   ! Conservation
        do i3=1,3   ! Three-momentum conservation
            p1c(i3) = p0(i3) - p1(i3)
        enddo
c       p0 refers to original q (qbar), p1 refers to generated qqbar pair, 
c       p1c refers to remnant (original parion after generating qqbr pair)

        e1c = w1c - p1c(3)
        p1c(4) = e1c

        kf0 = kf00
        if(kf0.gt.0) iflav =  1
        if(kf0.lt.0) iflav = -1
        w0 = w1c
        e0 = e1c
        do i3=1,4
            p0(i3) = p1c(i3)
        enddo
        igen = igen + 1

        n = n + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if( igen.ge.INT(adj16) ) goto 106   ! Stop generation 280822
        if( e0.le.adj17 ) goto 106   ! Stop generation
        if( w0.le.0. )    goto 106   ! Stop generation

        goto 100

106     continue
        if(igen.eq.0)return
c       Update four momentum of the remnant of ii-th source q (qbar). 280822
        do i3=1,4   ! 280822
            p(ii,i3) = p0(i3)
        enddo
        if(p00(3).lt.0.) p(ii,3) = -p0(3)   ! For the negative direction p_z.

c       Re-calculates E for the remenant, let its inv. mass >= 0.
        sm2 = p(ii,4)**2 - p(ii,1)**2 - p(ii,2)**2 - p(ii,3)**2
        if( sm2.lt.0. )then
        p(ii,4) = SQRT(p(ii,5)**2 + p(ii,1)**2 + p(ii,2)**2 +p(ii,3)**2)
        p(ii,4) = p(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - p(ii,4)
        end if


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine funcz(zz,ij)   ! 081222
c00623 This subroutine has been rewrote and fixed.   !Lei2023060
c       Sample daughter energy fraction, zz, from mother according to 
c        fragmentation function by selecting sample method.
c       adj1(29)=1 : Lund string fragmentation function   ! 140223 Lei
c               =2 : Field-Feymman fragmentation function
c200223         =3 : Perterson/SLAC fragmentation function
c               =4 : dummy, but set FF as default
c               =5 : dummy, but set FF as default
c       Lund string fragmentation function: f(z)=(1/z)*(1-z)^a
c        *exp(-b*m_{\per}^2/z) \sim (1/z)*(1-z)^a*exp(-2.36*b/z)
c        assume: m_{\per}^2 \sim <p_T>^2+m_u^2 \sim 1.5^2+0.333^2=2.36
c       FF fragmentation function: f(z)dz=[1-a+3*a*(1-z)**2]dz, 0<=z<=1,
c        and its largest value: fmax=f(0)=1-a+3*a=1+2*a.
c200223 P/S fragmentation function: f(z)dz=1/(z*(1-1/z-a/(1-z))^2)
c260223 P/S refers to Peterson/SLAC
c       ij: number of calling 'funcz'
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc   ! 081222
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/local_fmax/ fmax_value(5)   !Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop

c       Sets parameters.
        a   = adj1(6)   ! Lund a parameter
        b   = adj1(7)   ! Lund b parameter
        aFF = a_FF      ! Field-Feynman parameter
        c   = a_PS_c    ! Peterson/SLAC parameter
        i_z_samp_func = INT( adj1(29)-10. )   ! 11, 12, 13
c       Finds max values.
        if(iii.eq.1 .AND. ij.eq.0)then
            fmax_value = 0.
c           1 = LUND
            fmax = 0.
            z1   = 0.025
            do i1=1,20
                z   = z1 + (i1-1)*0.05
                fzz = (1./z) * (1.-z)**a * dexp(-2.36*b/z)
                if(fzz.gt.fmax) fmax = fzz
            enddo
            fmax_value(1) = fmax
c           2 = FF
            fmax = 1. + 2.*aFF
            if(aFF.lt.0) fmax = 1. - aFF
            fmax_value(2) = fmax
c           3 = P/S
            fmax = 0.
            z1   = 0.025
            do i1=1,20
                z   = z1 + (i1-1)*0.05
                fzz = 1./( z * (1.-1./z-c/(1.-z)) * (1.-1./z-c/(1.-z)) )
                if(fzz.gt.fmax) fmax = fzz
            enddo
            fmax_value(3) = fmax
c           4 = dummy, but set FF as default
            ! fmax = 
            fmax = fmax_value(2)
            fmax_value(4) = fmax
c           5 = dummy, but set FF as default
            ! fmax = 
            fmax = fmax_value(2)
            fmax_value(5) = fmax
        end if

!       Samples z.
        fmax = fmax_value( i_z_samp_func )

100     ran1 = pyr(1)
        ran2 = pyr(1)
        fm = ran1 * fmax

        select case( i_z_samp_func )
        case(1)   ! LUND
            fran2 = (1./ran2) * (1.-ran2)**a * dexp(-2.36*b/ran2)
        case(2)   ! FF
            fran2 = 1. - aFF + 3.*aFF*(1.-ran2)**2
            if(aFF.lt.0) fran2 = aFF + 3.*aFF*(1.-ran2)**2
        case(3)   ! P/S
            fran2 = 1./( ran2 * (1.-1./ran2-c/(1.-ran2))**2 )
        ! case(4)
        !     fran2 = 
        ! case(5)
        !     fran2 = 
        case default   ! Sets FF as default.
            fran2 = 1. - aFF + 3.*aFF*(1.-ran2)**2.
            if(aFF.lt.0) fran2 = aFF + 3.*aFF*(1.-ran2)**2
        end select

        if(fm.gt.fran2) goto 100
        zz = ran2


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine phas(i1,i2,i3,isucc,j,iphas)   ! 140223 Lei added iphas
c       The phase space judgement.
c       j=2 for meson
c       j=3 for baryon
c       iphas: = 1, complete phase space constraint  ! 140223 Lei
c              = 2, position constraint only
c              = 3, momentum constraint only
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,MPLIS=80000)   ! 150922
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 150922
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension ri1(3),ri2(3),ri3(3),pi1(3),pi2(3),pi3(3)
        delc=adj1(22)
        if( ABS(adj1(21) - 2.).lt.1D-15 .OR. 
     &      ABS(adj1(21) - 3.).lt.1D-15 ) delc=0.5*delc   ! 140223 Lei adj21=2 or 3

        if(j.eq.2)goto 100   ! for meson
c       proceed for baryon 
        ri1(1)=v(i1,1)
        ri1(2)=v(i1,2)
        ri1(3)=v(i1,3)
        ri2(1)=v(i2,1)
        ri2(2)=v(i2,2)
        ri2(3)=v(i2,3)
        ri3(1)=v(i3,1)
        ri3(2)=v(i3,2)
        ri3(3)=v(i3,3)
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
        pi1(1)=p(i1,1)
        pi1(2)=p(i1,2)
        pi1(3)=p(i1,3)
        pi2(1)=p(i2,1)
        pi2(2)=p(i2,2)
        pi2(3)=p(i2,3)
        pi3(1)=p(i3,1)
        pi3(2)=p(i3,2)
        pi3(3)=p(i3,3)
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

c140223 Lei
        if(iphas.eq.1)then   ! complete phase space constain
            del12=delr12*delp12
            del13=delr13*delp13
            del23=delr23*delp23
        elseif(iphas.eq.2)then   ! position constraint
            del12 = delr12
            del13 = delr13
            del23 = delr23
        elseif(iphas.eq.3)then   ! momentum constraint
            del12 = delp12
            del13 = delp13
            del23 = delp23
        else   ! no constraint
            del12 = 0D0
            del13 = 0D0
            del23 = 0D0
        end if
c140223 Lei

        if(del12.le.delc.and.del13.le.delc.and.del23.le.delc)then
        isucc=1
        else
        isucc=0
        endif
        return

100     continue   ! for meson
        ri1(1)=v(i1,1)
        ri1(2)=v(i1,2)
        ri1(3)=v(i1,3)
        ri2(1)=v(i2,1)
        ri2(2)=v(i2,2)
        ri2(3)=v(i2,3)
        ri121=ri1(1)-ri2(1)
        ri122=ri1(2)-ri2(2)
        ri123=ri1(3)-ri2(3)
        delr=sqrt(ri121*ri121+ri122*ri122+ri123*ri123)
        pi1(1)=p(i1,1)
        pi1(2)=p(i1,2)
        pi1(3)=p(i1,3)
        pi2(1)=p(i2,1)
        pi2(2)=p(i2,2)
        pi2(3)=p(i2,3)
        pi121=pi1(1)-pi2(1)
        pi122=pi1(2)-pi2(2)
        pi123=pi1(3)-pi2(3)
        delp=sqrt(pi121*pi121+pi122*pi122+pi123*pi123)

c140223 Lei
        if(iphas.eq.1)then   ! complete phase space constraint
            delrp = delr*delp
        elseif(iphas.eq.2)then   ! position consttaint
            delrp = delr
        elseif(iphas.eq.3)then   ! momentum constraint
            delrp = delp
        else   ! no constraint
            delrp = 0D0
        end if
c140223 Lei

        if(delrp.le.delc)then
        isucc=1
        else
        isucc=0
        endif
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break_f(eg,kf,amq)
c       sample flavor (mass) of generated qqbar pair
c       eg: energy of original q or qbar
c       kf (amq): flavor code (mass) of generated quark
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa38/csp_31,csp_32,csp_41,csp_42,csp_43,csp_51,csp_52,
     c   csp_53,csp_54,csp_61,csp_62,csp_63,csp_64,csp_65   ! 161022
        amd=pymass(1)   ! kinematical mass in GeV, 0.33
        amu=pymass(2)   ! amu=amd, 0.33
        ams=pymass(3)   ! 0.5
        amc=pymass(4)   ! 1.5
        amb=pymass(5)   ! 4.8
        amt=pymass(6)   ! 175
        ! amd=0.0099    ! current algebra mass
        ! amu=0.0056
        ! ams=0.199
        ! amc=1.23
        ! amb=4.17
        ! amt=165
        ! amd=0.325    ! constituent mass
        ! amu=0.325    ! amu=amd
        ! ams=0.5
        ! amc=1.6
        ! amb=5
        ! amt=xxx
        amuu=2*amu
        amdd=2*amd
        amss=2*ams
        amcc=2*amc
        ambb=2*amb
        amtt=2*amt
        aa=pyr(1)
c       if(eg.lt.amuu)goto 200   ! throw away amuu
c161022 if(eg.ge.amdd .and. eg.lt.amss)then   ! d,u
        if(eg.lt.amss)then   ! d,u (with same flavor generation probability)
          if(aa.le.0.5)then
          kf=1   ! d
          amq=amd
          else
          kf=2   ! u
          amq=amu
          endif
          goto 200
        endif

        ! if(eg.ge.amss .and. eg.lt.amcc)then   ! d,u,s
        if(eg.ge.amss)then   ! d,u,s
!Lei20230817B-
        !   if(aa.le.csp_31)then
        !   kf=1   ! d
        !   amq=amd
        !   elseif(aa.gt.csp_31 .and. aa.le.csp_32)then
        !   kf=2   ! u
        !   amq=amu
        !   else
        !   kf=3   ! s
        !   amq=ams
        !   endif
        !   goto 200
          kf  = 1 + INT( ( 2D0 + adj1(32) )*PYR(1) )   !Lei20230817
          amq = PYMASS(kf)    !Lei20230817
          return
!Lei20230817E-
        endif

        ! IF(.TRUE.)THEN   !Lei2023060
        IF(.FALSE.)THEN   !Lei2023060

        if(eg.ge.amcc .and. eg.lt.ambb)then ! d,u,s,c
          if(aa.le.csp_41)then
          kf=1   ! d
          amq=amd
          elseif(aa.gt.csp_41 .and. aa.le.csp_42)then
          kf=2   ! u
          amq=amu
          elseif(aa.gt.csp_42 .and. aa.le.csp_43)then
          kf=3
          amq=ams
          else
          kf=4
          amq=amc
          endif
          goto 200
        endif

c00623 if(eg.ge.ambb .and. eg.lt.amtt)then ! d,u,s,c,b
        if(eg.ge.ambb)then ! d,u,s,c,b !Lei2023060
          if(aa.le.csp_51)then
          kf=1
          amq=amd
          elseif(aa.gt.csp_51 .and. aa.le.csp_52)then
          kf=2
          amq=amu
          elseif(aa.gt.csp_52 .and. aa.le.csp_53)then
          kf=3
          amq=ams
          elseif(aa.gt.csp_53 .and. aa.le.csp_54)then
          kf=4
          amq=amc
          else
          kf=5
          amq=amb
          endif
          goto 200
        endif

c00623 Do not consider top.   !Lei2023060
c       if(eg.ge.amtt)then ! d,u,s,c,b,t
c         if(aa.le.csp_61)then
c         kf=1
c         amq=amd
c         elseif(aa.gt.csp_61 .and. aa.le.csp_62)then
c         kf=2
c         amq=amu
c         elseif(aa.gt.csp_62 .and. aa.le.csp_63)then
c         kf=3
c         amq=ams
c         elseif(aa.gt.csp_63 .and. aa.le.csp_64)then
c         kf=4
c         amq=amc
c         elseif(aa.gt.csp_64 .and. aa.le.csp_65)then
c         kf=5
c         amq=amb
c         else
c         kf=6
c         amq=amt
c         endif
c00623 endif
        ENDIF  !Lei2023060

200     continue
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tdgaus(v,pmax,np,pp)
c...... 2-d Gaussian distribution with width v, i.e., e^(-p2/v)dp2, 0<p2<pmax
c...... set pmax < 0 if pmax should be infinity.
c...... np : the total # of particles wanted to sample their transverse momenta
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 220312 240412 131212
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
10        p2 = -v*log(max(1.e-20,1. - a*pyr(1)))
          if(p2.LT.0.)goto 10
          ps=sqrt(p2)
          fi=2.*3.1415926*pyr(1)
c220312 randomly sample [px,py] on circle of sphere with radius ps
c220312   pp(i,1)=ps*cos(fi)
c220312   pp(i,2)=ps*sin(fi)
c220312 randomly sample [px,py] on circle of ellipsoid with half major axis
c220312 of ps*(1+smadel) and half minor axis of ps*(1-smadel)
          pp(i,1)=ps*cos(fi)*(1+smadel)   ! 220312
          pp(i,2)=ps*sin(fi)*(1-smadel)   ! 220312
c220312 note: ps is not in the dimension list
30      continue
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_glu   ! 160822
c       moves gluons from 'pyjets' to 'sa36'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        nglu=0
        if(iii.eq.1)then   !Lei2023060
        kglu=0
        pglu=0.
        vglu=0.
        endif   !Lei2023060
        jb=0
201     do i1=jb+1,n
        kf=k(i1,2)
        kfab=iabs(kf)
        eng=p(i1,4)
        if(kfab.ne.21)then   ! stay
        jb=jb+1
        goto 202
        endif

        nglu=nglu+1
        do i2=1,5
        kglu(nglu,i2)=k(i1,i2)
        pglu(nglu,i2)=p(i1,i2)
        vglu(nglu,i2)=v(i1,i2)
        enddo
        if(i1.eq.n)then
        n=n-1
        goto 203
        endif
c       move particle list one step downward from i1+1 to n
        do jj=1,5
        do j=i1+1,n
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 201
202     enddo   ! do loop
203     continue
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break_glu   ! 160822
c       break-up gluon (in 'sa36') and give flavor (mass) and four momentum 
c        (position) to broken objects, which are assumed to be a string and 
c        filling in 'pyjets'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104   Lei2023060
        common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 080104 220110
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) 
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        amu=pymass(2)   ! 271022
        amuu=2*amu   ! 271022
c       throw away g with energy<amuu
        jb=1       !Lei2023060
2010    continue   !Lei2023060
        do i1=jb,nglu,1   !Lei2023060 1 -> jb
        eg=pglu(i1,4)
        if(eg.lt.amuu)then
        do i2=1,4
        throe_p(i2)=throe_p(i2)+pglu(i1,i2)   !Lei2023060 Replace ppsa by throe_p.
        enddo
c00623 Lei2023060
        if(i1.eq.nglu)then
            nglu = nglu - 1
            goto 2020
        endif
c00623 Lei2023060
c       move particle list ('sa36') one step downward from i1+1 to nglu
        do jj=1,5
        do j=i1+1,nglu
        kglu(j-1,jj)=kglu(j,jj)
        pglu(j-1,jj)=pglu(j,jj)
        vglu(j-1,jj)=vglu(j,jj)
        enddo
        enddo
        nglu=nglu-1
        jb=i1       !Lei2023060
        goto 2010   !Lei2023060
        endif
        enddo
2020    continue    !Lei2023060

c       g (in 'sa36') -> qq_{bar} (as a string filling in 'pyjets')
100     do i1=1,nglu   ! do loop over gluons
        eg=pglu(i1,4)
        call break_f(eg,kf,amq)
        kf1=kf
        kf2=-kf
        am1=amq
        am2=amq

        k(n+1,1)=2   ! 'A'
        k(n+2,1)=1   ! 'V'
        k(n+1,2)=kf1
        k(n+2,2)=kf2
        k(n+1,3)=0
        k(n+2,3)=0
        k(n+1,4)=0
        k(n+2,4)=0
        k(n+1,5)=0
        k(n+2,5)=0
c       p(n+1,5)=am1
c       p(n+2,5)=am2

c       give four momentum to the breaked quarks
        call bream_glu(i1,kf1,kf2)
c       give four coordinate to the breaked quarks
        call coord_glu(i1)
        if(i1.eq.nglu)then
        n=n+2
        goto 200
        endif
        n=n+2
c       goto 100
        enddo   ! do loop over gluons ended
200     return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine bream_glu(ii,kf1,kf2)   ! 160822
c       give four momentum to the broken quarks
c       ii: line number of initial gluon in 'sa36'
c       kf1,kf2: flavor codes of broken quarks
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104   Lei2023060
        common/sa18/i_deex,i_deex_gen,i_pT,i_pT_max,i_split_diq,
     &   i_split_qqb,i_split_g,i_pad,a_FF,aPS_c,aPS_b,parj23,parj24   !Lei2023060
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)
        am1=pymass(kf1)
        am2=pymass(kf2)
        pp(1,5)=am1   ! mass of first broken quark
        pp(2,5)=am2   ! mass of second broken quark
c       pp : four momenta & mass of broken quarks, local variable
        do i1=1,4
        ps(i1)=pglu(ii,i1)
        enddo
c       ps : four momentum of initial gluon, local variable 

        i_p_split = i_split_g   !Lei2023060
c       i_p_split = 1: decay method   !Lei2023060
c                 = 2, random 3-momentum method with the different factors
c                 = 3, random 3-momentum method with the same factor for 3-mom
c                 = 4: equal division method

        if( i_p_split.eq.1 ) goto 400   ! for 'decay method' Lei2023060

c       Random three momentum method with the different factors. Lei2023060
401     do i1=1,3   ! for 'random three momentum method'
        pi(i1)=pyr(1)*pglu(ii,i1)
        pp(1,i1)=pi(i1)
        pp(2,i1)=ps(i1)-pi(i1)
        enddo

c00623 Lei2023060B--
c       Random three momentum method with the same factor.
        if( i_p_split.eq.3 )then
            factor=PYR(1)
            do i1=1,3,1
                pi(i1)=factor*pglu(ii,i1)
                pp(1,i1)=pi(i1)
                pp(2,i1)=ps(i1)-pi(i1)
            enddo
        end if

c       Equal division method.
        if( i_p_split.eq.4 )then
            do i1=1,3,1
                pi(i1)=0.5*pglu(ii,i1)
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
400     continue   ! for 'decay method'

c       Decay method.
        if( i_p_split.eq.1 )then   !Lei2023060
        decsuc=1
        call decmom(ps,pp,am1,am2,decsuc)
        if(decsuc.eq.0)goto 401   ! return to random three momentum method
        endif   !Lei2023060

300     continue
c       adjust four momentum conservation by iteration,no more than
c        4000 iterations
c       call conser(2,pp,ps)
c160822
c       do i1=1,2
c       do i2=1,4
c       ppi1i2=pp(i1,i2)
c       pp12=abs(ppi1i2)
c       if(pp12.gt.1.d6)pp12=1.d6
c       pp(i1,i2)=sign(pp12,ppi1i2)
c       enddo
c       enddo
c160822
        do i1=1,4
        p(n+1,i1)=pp(1,i1)
        enddo
        p(n+1,5)=am1
        do i1=1,4
        p(n+2,i1)=pp(2,i1)
        enddo
        p(n+2,5)=am2

        do i2=1,4   !Lei2023060 Collects lost 4-momentum.
        throe_p(i2) = throe_p(i2) + ( ps(i2) - p(n+1,i2) - p(n+2,i2) )
        enddo


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coord_glu(ii)   ! 160822
c       give four position to broken quarks
c       first broken quark takes the four position of gluon
c       second broken quark is arranged around first ones within
c        0.5 fm randumly in each of three positions and has same
c        fourth position as gluon
c       ii: order # of gluon in 'sa36'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) 
        dimension rr(3)

        do i1=1,3
        v(n+1,i1)=vglu(ii,i1)
        rr(i1)=pyr(1)*0.5   ! 261002
        v(n+2,i1)=vglu(ii,i1)+rr(i1)
        if(pyr(1).gt.0.5d0)v(n+2,i1)=vglu(ii,i1)-rr(i1)
        enddo
        v(n+1,4)=vglu(ii,4)
        v(n+2,4)=vglu(ii,4)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa36(nn,cc)
c       print particle list and sum of momentum and energy
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        dimension peo(4)
c       do i=1,nn
c       write(mstu(11),*)i,kglu(i,2),(pglu(i,j),j=1,4)
c       write(9,*)i,kglu(i,2),(pglu(i,j),j=1,4)
c       enddo
        call psum(pglu,1,nglu,peo)
        ich1=0
        do i1=1,nn
        kf=kglu(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3.
        write(22,*)'sa36 nn=',nn
        write(mstu(11),*)'c & p sum=',cc,peo   !
c       write(9,*)peo,ich1/3   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa37(nn,cc)   ! Lei2023060
c       print particle list and sum of momentum and energy
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa37/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
        call psum(pbh,1,nbh,peo)
        ich1=0
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3.
        ! write(22,*)'sa37 nn=',nn
        write(22,*)'c & p sum=',cc,peo   !
        do i=1,nn
        write(22,*)i,kbh(i,1),kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        write(22,*) "------------------------------------------------"//
     &              "------------------------------------------------"//
     &              "------------------------------------------------"
        return
        end



c*********************************************************************
        subroutine thptdi(pt,px,py)   ! 150822 
c       Generates transverse momentum according to thermodynamic distribution
c        (exponential distribution)
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        SAVE /PYDAT1/
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 220312 240412 131212
c       Generate p_T and azimuthal angle, gives p_x and p_y.
        pt=-parj(21)*log(max(1d-10,pyr(1)))
        phi=paru(2)*pyr(1)
c       randomly sample [px,py] on circle of sphere with radius pt
        px=pt*cos(phi)*(1+smadel)   !Lei2023060
        py=pt*sin(phi)*(1-smadel)   !Lei2023060
        pt=SQRT(px*px+py*py)        !Lei2023060
        return
        end



C*********************************************************************
C...PAPTDI
c00623 Lei2023060
C...Generates transverse momentum according to a Gaussian/Exponential/Random.
 
        SUBROUTINE PAPTDI(KF0,KF1,KF2,PX0,PY0,PT0,PX,PY,PT,i_pT,i_max,x)
 
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 112819
        SAVE /PYDAT1/
 
        PT=0D0
        x=0D0
100     continue
C...Generate p_T and azimuthal angle, gives Gaussian p_x and Gaussian p_y.
        if(i_pT.eq.1) PT=PARJ(21)*SQRT(-LOG(MAX(1D-10,PYR(1))))
 
C...Generate p_T and azimuthal angle, gives Exponential p_x and Exponential p_y.
        if(i_pT.eq.2) PT=-PARJ(21)*LOG(PYR(1)*PYR(1))
 
C...Generate p_T according to the Exponential distribution. (thermodynamic)
        if(i_pT.eq.3) PT=-PARJ(21)*LOG(MAX(1D-10,PYR(1)))
 
C...Non-uniform tail.
        if(PARJ(23).GT.PYR(1)) PT=PARJ(24)*PT
 
C...PT constrainting and re-sampling
        if(i_max.eq.1)then
            if(PT.gt.PT0) goto 100
        endif
 
C...Sampling p_T form given PT0 randomly.
        if(i_pT.eq.4 .OR. i_pT.eq.8)then
            factor=PYR(1)
            PT=factor*PT0
            x=factor
        endif
 
200     continue
C...Generate azimuthal angle.
        PHI=PARU(2)*PYR(1)
 
C...Randomly sample [PX,PY] on circle of sphere with radius PT
C       PX=PT*COS(PHI)
C       PY=PT*SIN(PHI)
C...Randomly sample [PX,PY] on circle of ellipsoid with half major axis 
C       of PT*(1+smadel) and half minor axis of PT*(1-smadel)
        PX=PT*COS(PHI)*(1.+smadel)
        PY=PT*SIN(PHI)*(1.-smadel)
 
C...PX, PY constrainting and re-sampling
        if(i_max.eq.1)then
            if(PX0.gt.0. .AND. PX.gt.PX0) goto 200
            if(PX0.lt.0. .AND. PX.lt.PX0) goto 200
            if(PY0.gt.0. .AND. PY.gt.PY0) goto 200
            if(PY0.gt.0. .AND. PY.gt.PY0) goto 200
        endif
 
C...Samples p_x and p_y from given PX0 and PY0 with different random factors.
        if(i_pT.eq.5) PX=PYR(1)*PX0
        if(i_pT.eq.5) PY=PYR(1)*PY0
 
C...Samples p_x and p_y from given PX0 and PY0 with a same random factor.
        if(i_pT.eq.6 .OR. i_pT.eq.7)then
            factor=PYR(1)
            PX=factor*PX0
            PY=factor*PY0
            x=factor
        endif
 
C...PT might deviate little bit from original one.
        PT=SQRT(PX*PX+PY*PY)
 
        RETURN
        END
