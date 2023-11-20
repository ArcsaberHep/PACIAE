        subroutine coales(ijk,neve,nnout,nap,nat,nzp,nzt,i_coal)   ! 250823 Lei Added i_coal.
c       A phenomenological coalescence model writen by Sa Ben-Hao on 04/06/2004
c       Its input messages are in 'pyjets'   ! 220822
c       Its storing array is 'pyjets'
c       Its output message is in 'sa1_h' (in 'pyjets' either)
c       ijk: the run number
c       neve: total number of runs
c       nnout: a internal printing per nnout runs
c       nap and nzp: atomic and charge number of projectile
c       nat and nzt: atomic and charge number of target
c       i_coal: with or without gluon splitting and deexcitation.    ! 250823 Lei
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
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 300623 Lei
        common/sa18/i_deex,n_deex_step,i_pT,i_pT_max,a_FF,aPS_c,aPS_b   ! 280823 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) ! 220822
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)   ! 150922
        common/sbh/nbh,nonh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   napp,natt,nzpp,nztt,pio
        dimension numb(3)


        if( i_coal.eq.0 ) goto 888  ! 250823 Lei


c-------------------------------------------------------------------------------
c-------------------------   Variables initialization   ------------------------
        rrp = 1.16
        nn  = 0
        if( ijk.eq.1 )then   ! 300623 Lei
            kn = 0
            pn = 0D0
            rn = 0D0
        endif   ! 300623 Lei
        nout  = nnout
        imc   = INT(adj1(13))
        ibc   = INT(adj1(14))
        iphas = INT(adj1(21))
c-------------------------   Variables initialization   ------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c----------------------------   Junctions removing   ---------------------------
        if( INT(adj1(40)).eq.3 )then   ! 070223
c220822     Remove junctions.
            jb = 0
2010        do i1=jb+1,N,1  ! i1 loop
                kf   = K(i1,2)
                kfab = ABS(kf)
                if(kfab.ne.88)then
                    jb = jb + 1
                    goto 2020
                endif
                call updad_pyj(N,i1+1,1)   ! 090922 'updad_pyj' in sfm_30.f
                N = N - 1
                goto 2010
2020        enddo   ! i1 loop
        endif   ! 070223
c----------------------------   Junctions removing   ---------------------------
c-------------------------------------------------------------------------------


!TODO(Lei20230825): The baryon number conservation process is unnecessary in fact.
c       It will be satisfied automatically according to the number of q and qbar.
c-------------------------------------------------------------------------------
c------------------------   Net baryon number counting   -----------------------
c       Conservation of net baryon.
        netba = 0   ! 300623 Lei
        do i1=1,nbh,1 ! note: 'sbh' is hadron list before hadronization 060119
            kf   = kbh(i1,2)
            kfab = ABS(kf)
            if(kf.gt.0.and.(kf.gt.1000 .and. kf.lt.10000)) netba=netba+1
        if(kf.lt.0.and.(kfab.gt.1000 .and. kfab.lt.10000)) netba=netba-1
        enddo
c060813 120214
        if( nap.gt.1 .and. nat.gt.1 )then
            netba = nap + nat - netba   ! AA
        elseif( nap.eq.1.and.nat.gt.1.and.ipden.eq.0 )then
            netba = nzp + nat - netba   ! p-A (pbar-A)
        elseif( nap.gt.1.and.nat.eq.1 )then
            netba = nap + nzt - netba   ! A-p (A-pbar)
        elseif( nap.eq.1.and.nat.eq.1 )then
            netba = nzp + nzt - netba   ! pp (p-pbar,pbar-p)
        else
            netba = nat - netba   ! lepton-A
        endif
c060813 120214
c060119 if(ijk.eq.1)then
c       napt=nap+nat
c       nzpt=nzp+nzt
c       sbaryi=float(napt)
c       schgei=3.*float(nzpt)
c060119 endif
c------------------------   Net baryon number counting   -----------------------
c-------------------------------------------------------------------------------


c300623 g-splitting and q-deexcitation before "parcas".   ! 300623
888     if( INT(adj1(12)).eq.2 .AND. i_coal.eq.1 ) goto 1000   ! 250823 Lei


c-------------------------------------------------------------------------------
c-----------------------------   Gluon splitting   -----------------------------
c220122
c       Move gluons from 'pyjest' to 'sa36'.
        call remo_glu
c       Break-up gluon (with E_g>2E_u in 'sa36') -> qqbar string 
c        (filling in 'pyjets').
        call break_glu
c       So far, the parton list ('pyjets') is composed of q and qbar only.
        adj12 = adj1(12)
        i_deex_gen = INT( adj1(16) )   ! 180923 Lei
        adj17 = adj1(17)
c200222 adj17=max(4.0,adj17) ! 070612, yan

c300623 Shares 4-momentum in "throe_p" among partons.   ! 300623 Lei
        call share_p_PYJETS   ! 300623 Lei
c-----------------------------   Gluon splitting   -----------------------------
c-------------------------------------------------------------------------------


c250823 Debug mode.   ! 250823 Lei
        if( INT(adj1(12)).eq.3 ) goto 1000   ! 310723 Lei w/ g-splitting and w/o deexc.


c-------------------------------------------------------------------------------
c---------------------------   Quark deexcitation   ----------------------------
c280822 energetic q (qbar) de-excitation
        n00   = N   ! Original total entries in PYJETS
        i_call_deex  = 0
        i_daught_gen = 1   ! the #-th newly produced daughter qqbar
        n_deex = 0   ! the number of successful deexcitation
        jb = 0
        n0 = N   ! Current total entries in PYJETS
700     continue
        do i1=jb+1,n0,1
            kf0   = K(i1,2)
            ee    = P(i1,4)
            iflav = 1
            if( kf0.lt.0 ) iflav = -1
c           iflav = 1 : if source parton is quark
c                 =-1 : if source parton is anti-quark
            if( ee.gt.adj17 )then
              if( i_deex.eq.1 ) call deexcitation_EP(i1,kf0,nstep,iflav)   ! 300623 Lei
              if( i_deex.eq.2 ) call deexcitation_E(i1,kf0,nstep,iflav)    ! 300623 Lei
        if(i_deex.eq.3) call deexcitation_EP_comp_pT(i1,kf0,nstep,iflav)   ! 310723 Lei
        if(i_deex.eq.4) call deexcitation_E_comp_pT(i1,kf0,nstep,iflav)    ! 310723 Lei
                if( nstep.gt.0 ) n_deex = n_deex + 1   ! 300623 Lei
                i_call_deex = i_call_deex + 1   ! times of 'call deexcitation'
            endif
c           nstep : number of deexcitation steps per source q (qbar)
c300623 Lei
c           Updates n0 and does deexcitation for newly produced qqbar pairs.
        if( i1.eq.n0 .AND. N.gt.n0 .AND. i_daught_gen.lt.i_deex_gen)then
            jb = i1
            i_daught_gen = i_daught_gen + 1
            n0 = N
            goto 700
        end if
c300623 Lei
800     enddo   ! 280822 continue->enddo
900     continue
c       energetic q (qbar) de-excitation, finished.

c300623 Shares the 4-momentum in 'throe_p' among partons.   ! 300623 Lei
        call share_p_PYJETS   ! 300623 Lei
c220122
c---------------------------   Quark deexcitation   ----------------------------
c-------------------------------------------------------------------------------


1000    if( i_coal.eq.0 ) return  ! 300623 Lei For adj12 = 2


c-------------------------------------------------------------------------------
c-----------------------------   Parton sorting   ------------------------------
        numb = 0   ! 220822
c       numb(1),(2), and (3): the order number of last g, qbar & q in 'pyjets'
c220822 1: refers to g (no gluon at all now), 2: qbar, 3: q

c       Sort the partons in order of qbar and q,   ! 220822
c        i.e. move q to the end.   ! 220822
        jh = N
        jl = 0
2030    continue
        do j=jl+1,jh,1
            kf   = K(j,2)
            kfab = ABS(kf)
            if( kfab.lt.7 .and. kf.gt.0 )then   ! q, consider d, u, s, c, b, t
                N = N + 1
                numb(3) = numb(3) + 1
                do i4=1,5,1   ! 300623 Lei Lei 4 -> 5
                    K(N,i4) = K(j,i4)
                    P(N,i4) = P(j,i4)
                    V(N,i4) = V(j,i4)
                enddo
                call updad_pyj(N,j+1,1)
                N  = N  - 1
                jh = jh - 1
                jl = j  - 1
                goto 2030
            endif
        enddo
        numb(1) = 0
        numb(2) = N - numb(3)
        numb(3) = N
        n1 = numb(1)
        n2 = numb(2)
        n3 = N
c       Order the qbar according to energy from the maximal to minimal.
c300623 call eord(n1+1,n2)   ! 300623 Lei
c       Order the q according to energy from the maximal to minimal.
c300623 call eord(n2+1,n3)   ! 300623 Lei
c-----------------------------   Parton sorting   ------------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c----------------------------   Parton coalescence   ---------------------------
c       Load the table of hadron (meson: pseudoscalar-spin 0 & vector-spin 1 
c        only, baryon: octet-spin 1/2 & decuplet-spin 3/2 only).
        if( ijk.eq.1 ) call tabhb
c       ijk is the event number.

c       Normal cpalescence process.
        iqba = n2
        call coal(n3,iqba,ijk,rrp,iphas,netba)
c       n3: total number of partons (qbar and q)
c       iqba: total number of qbar (qbar is ordered before q)
c       ithroq : the total number of quarks thrown away
c       ithrob : the total number of anti-quarks thrown away
c       throe : total 4-momentum of the partons thrown away
c       ich : total charge of the partons thrown away

c070223 Re-coalesce failed parton after last 'call coal'.
c       Try coalescence without phase-space constraint agian for those partons 
c        failed in normal coalescence process with constraint originally.
        iqba = ithrob   ! Not active originally.
        n3   = ithroq + ithrob   ! Not active originally.
        if( iphas.ne.0 .and. n3.ge.2 )then   ! 070223
            call coal(n3,iqba,ijk,rrp,0,0)   ! 110905
        endif

c150922 ichth=ich   ! 092600
c----------------------------   Parton coalescence   ---------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c------------------------------   Data dumping   -------------------------------
c       'sa1_h' to 'pyjets'.
        N = nn
        do j2=1,5,1
            do j1=1,nn,1
                K(j1,j2) = kn(j1,j2)
                P(j1,j2) = pn(j1,j2)
                V(j1,j2) = rn(j1,j2)
            enddo
        enddo
c------------------------------   Data dumping   -------------------------------
c-------------------------------------------------------------------------------


c       Decay of unstable hadrons.
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
C...Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
c       imc (ibc): dimension of meson (baryon) table considered
c       kqh(i,j): flavor code of j-th constituent q (or qbar) in i-th quark 
c        configuration of meson
c       kfh(i,1), kfh(i,2): flavor code of pseudoscalar meson and vector meson 
c        with i-th quark configuration of meson
c       proh(i,1), proh(i,2): proper probability of meson
c       amash(i,1), amash(i,2): mass of meson
c       kqb(i,j): flavor code of j-th constituent q (or qbar) in i-th quark
c        configuration of baryon
c       kfb(i,1), kfb(i,2): flavor code of octet baryon and decuplet baryon 
c        with i-th quark configuration of baryon
c       prob(i,1), prob(i,2): probability of baryon
c       amasb(i,1), amasb(i,2): mass of baryon

c300623 New table.   ! 300623 Lei
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
c       Note that kqh1 > 0 and kqh2 < 0.
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
! Note that kqb1 <= kqb2 <= kqb3.
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
c       Octet baryon, s = 1/2.
        data (kfb(i,1),i=1,80)
     &   /    0, 2112, 3112, 4112, 5112, 2212, 3122, 3212, 4122, 4212,   ! 10
     1     5122, 5212, 3312, 4132, 4312, 5132, 5312, 4412, 5142, 5412,   ! 20
     2     5512,    0, 3222, 4222, 5222, 3322, 4232, 4322, 5232, 5322,   ! 30
     3     4422, 5242, 5422, 5522,    0, 4332, 5332, 4432, 5342, 5432,   ! 40
     4     5532,    0, 5442, 5542,    0,                                 ! 45
     4         0,   0,    0,    0,    0,                                 ! 50
     5     30*0 /                                                        ! 80
c       Decuplet baryon, s = 3/2.
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
c300623 New table.   ! 300623 Lei

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
          alar=P(ii,4)
          do 200 i2=i1+1,nc
c280822 communication between i1 and i2 of which the energy is lagest
c        among i1+1, i1+2, ..., nc
            ee=P(i2,4)
            if(ee.gt.alar)then   ! 280822 .ge. -> .gt.
                j=i2
                alar=ee
            endif
200       enddo   ! continue-> enddo 280822
c280822 now, j: order number of particle with largest energy
          kii2=K(ii,2)
          do jj=1,4
          p1(jj)=P(ii,jj)
          rr(jj)=V(ii,jj)
          enddo
          pii5=P(ii,5)

          K(ii,2)=K(j,2)
          do jj=1,4
          P(ii,jj)=P(j,jj)
          V(ii,jj)=V(j,jj)
          enddo
          P(ii,5)=P(j,5)

          K(j,2)=kii2
          do jj=1,4
            P(j,jj)=p1(jj)
            V(j,jj)=rr(jj)
          enddo
          P(j,5)=pii5
100     enddo

        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coal(n1,nqb,ijk,rrp,iphas,netba)
c       Parton coalescence (hadronization)
c       n1 : total number of partons (q & qbar only)
c       nqb : total number of qbar (qbar is ordered before q)
c       ijk : the run number
c       iphas: = 1, complete phase space constraint  ! 300623 Lei
c              = 2, position constraint only
c              = 3, momentum constraint only
c       netba: number of baryons to be first generated keeping
c              baryon conservation 
c       ithroq : total number of quarks thrown away
c       ithrob : total number of anti-quarks thrown away
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
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 300623 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)   ! 150922
        dimension pc(4),rc(3),iar(3),rcp(3)
        dimension psu(3),peo(5),pnn(kszj,5)
        dimension numb(3)   ! 110905


c-------------------------------------------------------------------------------
c-------------------------   Variables initialization   ------------------------
c300623 Moved from 'coales' to here.   ! 300623 Lei
        ithroq = 0
        ithrob = 0
        ich    = 0
        throe  = 0D0
c       Arrays for partons thrown away.
        nth = 0
        if( ijk.eq.1 )then   ! 300623 Lei
            kth = 0
            pth = 0D0
            vth = 0D0
        endif   ! 300623 Lei
c-------------------------   Variables initialization   ------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c-------------------------------   Coalescence   -------------------------------
        ibarp = 0   ! Numer of baryon generated. (baryon plus)
        ibarm = 0   ! Number of anti-baryon generated. (baryon minus)
        imes  = 0   ! Number of meson generated.

!TODO(Lei20230825): The baryon number conservation process is unnecessary in fact.
c       It will be satisfied automatically according to the number of q and qbar.
c----------------------   Baryon conservation treatment   ----------------------
c       Generate first 'netba' baryons (if 'netba'>0) or -'netba' anti-baryons 
c        (if 'netba'<0) to keep baryon conservation.
        if( netba.gt.0 ) call barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,1)
        if( netba.lt.0 )then
302       continue
          do ii1=1,nqb
            call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,1)
            if( isu.eq.1 .and. ibarm.lt.-netba ) goto 303
            if( isu.eq.1 .and. ibarm.eq.-netba ) goto 304
            if( isu.eq.0 ) stop 6666
          enddo
303       goto 302
304       continue
        endif
c----------------------   Baryon conservation treatment   ----------------------

c----------------------------   Normal coalescence   ---------------------------
c       Compose first the anti-quark, one by one, into hadron. 150922
c       An anti-quark can be a component of anti-baryon or meson.
        jj = 0
300     continue   ! Start the first and next coalescence. Cycle beginning.
        jj = jj + 1

c       No qbar. Try to compose baryons.
        if( nqb.eq.0 ) goto 100
c       No q and number of qbar is lees than 3. One cannot compose 
c        mesons or anti-baryons. Throw away qbar.
        if( nqb.lt.3 .and. n1.eq.nqb ) goto 301

c       Select an anti-quark from list randomly.
        if( nqb.eq.1 ) ii1 = 1
        if( nqb.gt.1 ) ii1 = INT( nqb*PYR(1) + 1 )

        kf = K(ii1,2)
c       Probability of anti-baryon/meson.
        relpr = adj1(31)
        relpr = relpr/(1.+relpr)
c       Extra suppression of strange anti-baryon/meson.
        if(kf.eq.-3)then
            relpr = adj1(31)*adj1(33)
            relpr = relpr/(1.+relpr)
        endif

        rand = PYR(1)
c       Compose an anti-baryon.
        if( rand.le.relpr .and. nqb.ge.3 )then   !
          call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,2)
          if( isu.eq.1 .and. nba.eq.1 )then   ! One anti-baryon produced.
            goto 300
          endif

c       Fail to compose ii1 into an anti-baryon. Try to compose a meson.
          if( isu.eq.0 )then   !!
            call mespro(n1,nqb,ii1,iphas,nme,imes,rrp,isu)
            if( isu.eq.1 .and. nme.eq.1 )then   ! One meson produced.
              goto 300
            endif
c       Fail to compose ii1 into a meson, either. Throw away the qbar.
            if( isu.eq.0 )then   !!!
              nth = nth + 1
              ithrob = ithrob + 1
              kk  = K(ii1,2)
              ich = ich + PYCHGE(kk)
              kth(nth,2) = kk
              vth(nth,5) = V(ii1,5)   ! 250823 Lei
              pth(nth,5) = P(ii1,5)
              do i2=1,4
                throe(i2)   = throe(i2) + P(ii1,i2)
                pth(nth,i2) = P(ii1,i2)
                vth(nth,i2) = V(ii1,i2)
              enddo
c       Move parton list one step downward from ii1+1 to n1.
              call updad_pyj(n1,ii1+1,1)
              N   = N   - 1   ! 300623 Lei
              n1  = n1  - 1
              nqb = nqb - 1
              goto 300
            endif   !!!
          endif   !! Fail to compose ii1 into an anti-baryon.

c       Compose a meson.
        elseif( rand.gt.relpr .or. nqb.lt.3 )then   !
          call mespro(n1,nqb,ii1,iphas,nme,imes,rrp,isu)
          if( isu.eq.1 .and. nme.eq.1 )then   ! One meson produced.
            goto 300
          endif

c       Fail to compose ii1 into a meson. Try to compose an anti-baryon.
          if( isu.eq.0 )then   !!
            call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,2)
            if( isu.eq.1 .and. nba.eq.1 )then   ! One anti-baryon produced.
              goto 300
            endif
c       Fail to compose ii1 into an anti-baryon, either. Throw away the qbar.
            if( isu.eq.0 )then !!!
              nth = nth+1
              ithrob = ithrob + 1
              kk  = K(ii1,2)
              ich = ich + PYCHGE(kk)
              kth(nth,2) = kk
              vth(nth,5) = V(ii1,5)   ! 250823 Lei
              pth(nth,5) = P(ii1,5)
              do i2=1,4
                throe(i2)   = throe(i2) + P(ii1,i2)
                pth(nth,i2) = P(ii1,i2)
                vth(nth,i2) = V(ii1,i2)
              enddo
c       Move parton list one step downward from ii1+1 to n1.
              call updad_pyj(n1,ii1+1,1)
              N   = N   - 1   ! 300623 Lei
              n1  = n1  - 1
              nqb = nqb - 1
              goto 300
            endif   !!!
          endif   !!

        endif  !
c----------------------------   Normal coalescence   ---------------------------

c-------------------------------   Coalescence   -------------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c-------------------------   Failed quarks treatment   -------------------------
c       Throw away those qbar remained.
301     continue
        if( nqb.eq.0 ) goto 100   ! No qbar. Try to compose baryons.
        do i1=1,nqb,1
          nth = nth + 1   ! 110905
          ithrob = ithrob + 1
          kk  = K(i1,2)
          ich = ich + PYCHGE(kk)
          kth(nth,2) = kk
          vth(nth,5) = V(i1,5)   ! 250823 Lei Fixed from V(nth, to vth(nth
          pth(nth,5) = P(i1,5)   ! 250823 Lei
          do i2=1,4
            throe(i2)   = throe(i2) + P(i1,i2)
            pth(nth,i2) = P(i1,i2)
            vth(nth,i2) = V(i1,i2)
          enddo
        enddo
c       Move parton list nqb steps downward from nqb+1 to n1.
        call updad_pyj(n1,nqb+1,nqb)
        N   = N  - nqb   ! 300623 Lei
        n1  = n1 - nqb
        nqb = 0

100     continue
        if( n1-nqb.eq.0 ) goto 505   ! No quarks. Get out of 'coal'.
        if( n1-nqb.lt.3 ) goto 500   ! Throw away those quarks.
c       No qbar. Try to compose baryons.
        call barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,0)
        if( n1-nqb.eq.0 ) goto 505  ! No quarks. Get out of 'coal'.

500     continue
        if( n1-nqb.eq.0 ) goto 505   ! No quarks. Get out of 'coal'.
c       Throw away remained q.   150922
        ithroq = ithroq + n1-nqb
        do i1=nqb+1,n1,1
          nth = nth + 1   ! 110905 171122 Lei Corrected it from nth=nth+i1.
          kk  = K(i1,2)
          ich = ich + PYCHGE(kk)
          kth(nth,2) = K(i1,2)
          vth(nth,5) = V(i1,5)
          pth(nth,5) = P(i1,5)   ! 250823 Lei
          do i2=1,4
            throe(i2)   = throe(i2) + P(i1,i2)
            pth(nth,i2) = P(i1,i2)
            vth(nth,i2) = V(i1,i2)
          enddo
        enddo
        N   = nqb   ! 250823 Lei
        n1  = nqb
c-------------------------   Failed quarks treatment   -------------------------
c-------------------------------------------------------------------------------


505     continue


c-------------------------------------------------------------------------------
c---------------------------   List reconstructing   ---------------------------
c300623 Lei
c       Reconstruct parton list for calling 'coal' again.
c       Sort the parton list ('sa37') in order of qbar and q.
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
c       Move particle list 'sa37' one step downward from j+1 to nth.
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
        N = nth
        do j=1,5,1
            do i=1,nth,1
                K(i,j) = kth(i,j)
                P(i,j) = pth(i,j)
                V(i,j) = vth(i,j)
            end do
        end do
        ithroq = 0
        ithrob = 0
        ich    = 0
        throe  = 0.
c300623 Lei
c---------------------------   List reconstructing   ---------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c--------------------------   Final coalescence try   --------------------------
c300623 Lei
c       Final try for meson production. No more than 50 times.
        if( i_qbar.gt.0 .AND. n_q.gt.0 )then
            do i=1,50,1
                ii1 = INT( i_qbar*PYR(1) + 1 )
                call mespro(n_tot,i_qbar,ii1,iphas,nme,imes,rrp,isu)
                    if( isu.eq.1 .and. nme.eq.1 )then   ! One meson produced.
                        n_q = n_q - 1
                        if( i_qbar.eq.0 .OR. n_q.eq.0) exit
                    end if
            end do
            if(i.eq.50) write(*,*) "Final mespro failed! iii, N" ,ijk,N
        end if
c       Final try for baryon production. No more than 50 times.
        if( n_q.gt.2 )then
            do i=1,50,1
                call barpro(n_tot,i_qbar,iphas,nba,ibarp,netba,rrp,0)
                if( n_tot.eq.i_qbar ) exit
            end do
            if(i.eq.50) write(*,*) "Final barpro failed! iii, N" ,ijk,N
        end if
        n_q = n_tot - i_qbar
c       Final try for anti-baryon production. No more than 50 times.
        if( i_qbar.gt.2 )then
            do i=1,50,1
                ii1 = INT( i_qbar*PYR(1) + 1 )
        call an_barpro(n_tot,i_qbar,ii1,iphas,nba,ibarm,netba,rrp,isu,2)
                if( isu.eq.1 .and. nba.eq.1 )then   ! One anti-baryon produced.
                    if( i_qbar.eq.0 ) exit
                end if
            end do
        if(i.eq.50) write(*,*) "Final an_barpro failed! ijk, N" ,ijk,N
        end if
c300623 Lei
c--------------------------   Final coalescence try   --------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c------------------------------   Data dumping   -------------------------------
c300623 Lei
c       'PYJETS' -> 'sa37'
        if(N.gt.0)then
            do i=1,N,1
                do j=1,5,1
                    kth(i,j) = K(i,j)
                    pth(i,j) = P(i,j)
                    vth(i,j) = V(i,j)
                end do
                if( K(i,2).gt.0 )then
                    ithroq = ithroq + 1
                    ich    = ich    + PYCHGE(K(i,2))
                    throe  = throe  + P(i,4)
                elseif( K(i,2).lt.0 )then
                    ithrob = ithrob + 1
                    ich    = ich    + PYCHGE(K(i,2))
                    throe  = throe  + P(i,4)
                end if
            end do
            nth = N
        else
            nth = 0
        end if
c300623 Lei
c------------------------------   Data dumping   -------------------------------
c-------------------------------------------------------------------------------


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
c       To compose baryon.
c       n1 : total number of partons (q & qbar)
c       nqb : total number of qbar (qbar is ordered before q)
c       iphas: = 1, complete phase space constraint  ! 300623 Lei
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
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 300623 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension rcp(3)

        ! dpmax = adj1(27)
        ! drmax = adj1(28)
        nba = 0
400     continue

        if( n1-nqb.lt.3 ) return   ! Number of quarks is less than three.

        do 404 i1=nqb+1,n1-2,1
            kf1 = K(i1,2)   ! 150922
            do 405 i2=i1+1,n1-1,1
                kf2 = K(i2,2)
                do 406 i3=i2+1,n1,1
                kf3 = K(i3,2)

                sume  = P(i1,4) + P(i2,4) + P(i3,4)
                sump1 = P(i1,1) + P(i2,1) + P(i3,1)
                sump2 = P(i1,2) + P(i2,2) + P(i3,2)
                sump3 = P(i1,3) + P(i2,3) + P(i3,3)
                cm = sume*sume - sump1*sump1 - sump2*sump2 - sump3*sump3
                if( cm.gt.1D20 ) cm = 1D20
                if( cm.le.0D0  ) goto 406   ! fail 071204
                cm = SQRT(cm)

c       Find out the primary baryon from hadron table.
                call findb(kf1,kf2,kf3,cm,kfii,amasi,isucc,1)

                if( isucc.eq.0 ) goto 406   ! fail, and keep on cycle, try again.

c070223
c       Phase space adjudgment.
                if( iphas.ne.0 )then   ! Change from '.eq.1' to '.ne.0'.
                    call phas(i1,i2,i3,isucc,3,iphas)   ! 140223 Lei Added iphas.
                    if( isucc.eq.0 ) goto 406   ! fail
                endif   !!1
c070223

c       Proceed for success.
                ibarp = ibarp + 1
                nba   = nba   + 1

c       Give proper variables to the primary baryon.
                nnol  = nn
                nn    = nn+1
                kn(nn,1) = 1
                kn(nn,2) = kfii
                pn(nn,5) = amasi
                pn(nn,1) = sump1
                pn(nn,2) = sump2
                pn(nn,3) = sump3
                pnnm  = sump1*sump1 + sump2*sump2 + sump3*sump3
                pnnmm = amasi*amasi + pnnm
                if( pnnmm.gt.1D20 ) pnnmm = 1D20
                if( pnnmm.le.0D0  ) pnnmm = 1D-20
                pnnn = SQRT(pnnmm)
                pn(nn,4)   = pnnn
                ! dele     = sume - pnnn
                throe_p(4) = throe_p(4) + sume-pnnn   ! 300623 Lei

c       Produced hadron is arranged among constituent partons randomly.
                pyrx = PYR(1)
                pyry = PYR(1)
                pyrz = PYR(1)
                rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1) + pyrz*V(i3,1)
                rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2) + pyrz*V(i3,2)
                rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3) + pyrz*V(i3,3)

c       Move parton list one step downward from i3+1 to n1.
411             call updad_pyj(n1,i3+1,1)
                n1 = n1 - 1
                N  = N  - 1   ! 300623 Lei

c       Move parton list one step downward from i2+1 to n1.
                call updad_pyj(n1,i2+1,1)
                n1 = n1 - 1
                N  = N  - 1   ! 300623 Lei

c       Move parton list one step downward from i1+1 to n1.
                call updad_pyj(n1,i1+1,1)
                n1 = n1 - 1
                N  = N  - 1   ! 300623 Lei

c300623 Share the surplus 4-momentum in throe_p.   ! 300623 Lei
                ! call share_p_PYJETS   ! 300623 Lei
                call share_p_PYJETS_sa1h   ! 300623 Lei
                if( iway.eq.1 )then
                    if( nba.lt.netba ) goto 400 ! Recycle all the partons remained, do again.
                    if( nba.eq.netba ) return
                endif
                if( iway.eq.2 .and. nba.eq.1 ) return
                if( iway.eq.0 ) goto 400   ! 121204
c       iway=1: when the baryon number equal net baryon, return. Used in conserving net baryon number.
c       iway=2: generate one baryon. 
c       iway=0: check all the probability of parton constituent baryon, then return.
406             continue   ! fail
405         continue
404     continue


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine an_barpro(n1,nqb,i1,iphas,nba,ibarm,netba,rrp,isu,
     c   iway)
c       To compose an anti-baryon.
c       n1 : total number of partons (q & qbar)
c       nqb : total number of qbar (qbar is ordered before q)
c       i1: anti-quark wanted to compose anti-baryon
c       iphas: = 1, complete phase space constraint  ! 300623 Lei
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
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 300623 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension rcp(3)

        ! dpmax = adj1(27)
        ! drmax = adj1(28)
        isu = 1
        nba = 0

        if( nqb.lt.3 ) goto 100   ! Number of qbar is less than three, return.

        kf1 = K(i1,2)
        do 405 i2=1,nqb,1
            if( i2.eq.i1 ) goto 405
            kf2 = K(i2,2)
            do 406 i3=1,nqb,1
                if( i3.eq.i1 ) goto 406
                if( i3.eq.i2 ) goto 406
                kf3 = K(i3,2)

                sume  = P(i1,4) + P(i2,4) + P(i3,4)
                sump1 = P(i1,1) + P(i2,1) + P(i3,1)
                sump2 = P(i1,2) + P(i2,2) + P(i3,2)
                sump3 = P(i1,3) + P(i2,3) + P(i3,3)
                cm = sume*sume - sump1*sump1 - sump2*sump2 - sump3*sump3
                if(cm.gt.1D20) cm = 1D20
                if( cm.le.0D0 ) goto 406   ! fail 071204
                cm = SQRT(cm)

c       Find out the primary anti-baryon from hadron table.
                call findb(-kf1,-kf2,-kf3,cm,kfii,amasi,isucc,-1)

                if( isucc.eq.0 ) goto 406   ! fail

c070223
c       Phase space adjudgment.
                if( iphas.ne.0 )then   ! Changed from '.eq.1' to '.ne.0'.
                    call phas(i1,i2,i3,isucc,3,iphas)   ! 140223 Lei Added iphas.
                    if( isucc.eq.0 ) goto 406   ! fail
                endif   !!1
c070223

c       Proceed for success.
                goto 400

406         continue   ! fail
405     continue
        goto 100

400     continue
        ibarm = ibarm + 1
        nba   = nba + 1

c       Give proper variables to the primary anti-baryon.
        nnol = nn
        nn   = nn + 1
        kn(nn,1) = 1
        kn(nn,2) = kfii
        pn(nn,5) = amasi
        pn(nn,1) = sump1
        pn(nn,2) = sump2
        pn(nn,3) = sump3
        pnnm  = sump1*sump1 + sump2*sump2 + sump3*sump3
        pnnmm = amasi*amasi + pnnm
        if( pnnmm.gt.1D20 ) pnnmm = 1D20
        if( pnnmm.le.0D0 )  pnnmm = 1D-20
        pnnn = SQRT(pnnmm)
        pn(nn,4)   = pnnn
        ! dele     = sume - pnnn
        throe_p(4) = throe_p(4) + sume-pnnn   ! 300623 Lei

c       Produced hadron is arranged among contituent partons randomly.
        pyrx = PYR(1)
        pyry = PYR(1)
        pyrz = PYR(1)
        rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1) + pyrz*V(i3,1)
        rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2) + pyrz*V(i3,2)
        rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3) + pyrz*V(i3,3)

c       Move parton list one step downward from i3+1 to n1.
411     call updad_pyj(n1,i3+1,1)
        if( i1.gt.i3 ) i1 = i1 - 1
        if( i2.gt.i3 ) i2 = i2 - 1
        nqb = nqb - 1
        n1  = n1  - 1
        N   = N   - 1   ! 300623 Lei

c       Move parton list one step downward from i2+1 to n1.
        call updad_pyj(n1,i2+1,1)
        if( i1.gt.i2 ) i1 = i1 - 1
        nqb = nqb - 1
        n1  = n1  - 1
        N   = N   - 1   ! 300623 Lei

c       Move parton list one step downward from i1+1 to n1.
        call updad_pyj(n1,i1+1,1)
        nqb = nqb - 1
        n1  = n1  - 1
        N   = N   - 1   ! 300623 Lei

c300623 Share the surplus 4-momentum in throe_p.   ! 300623 Lei
        ! call share_p_PYJETS
        call share_p_PYJETS_sa1h
c       iway=1: creat an anti-baryon and return
c       iway=2: creat an anti-baryon and a baryon, then return
        if( iway.eq.1 .and. nba.eq.1 ) return
        if( iway.eq.2 .and. nba.eq.1 ) then
c       An anti-baryon generation must be followed immediately a baryon
c        generation to keep baryon conservation.
            call barpro(n1,nqb,iphas,nbaa,ibarp,netba,rrp,2)
            return
        endif

100     isu=0


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mespro(n1,nqb,i1,iphas,nme,imes,rrp,isu)
c       Compose a meson.
c       n1 : total number of partons (q & qbar)
c       nqb : total number of qbar (qbar is ordered before q)
c       i1: anti-quark wanted to compose a meson
c       iphas: = 1, complete phase space constraint  ! 300623 Lei
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
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 300623 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension rcp(3)

        ! dpmax = adj1(27)
        ! drmax = adj1(28)
        isu = 1
        nme = 0

        if( n1.eq.nqb ) goto 100   ! 300105 No quark, return.

        kf1 = K(i1,2)
        do 102 i2=nqb+1,n1,1
            kf2 = K(i2,2)

            sume  = P(i1,4) + P(i2,4)
            sump1 = P(i1,1) + P(i2,1)
            sump2 = P(i1,2) + P(i2,2)
            sump3 = P(i1,3) + P(i2,3)
            cm = sume*sume - sump1*sump1 - sump2*sump2 - sump3*sump3
            if( cm.gt.1.e20 ) cm = 1D20
            if( cm.le.0D0   ) goto 102   ! fail 071204
            cm = SQRT(cm)

c       Find out primary meson from hadronic table according to kf2 & kf1
            call findm(kf2,kf1,cm,kfii,amasi,isucc,1)

            if( isucc.eq.0 ) goto 102   ! fail

c070223
c       Phase space adjudgment.
        if( iphas.ne.0 )then   ! Changed from '.eq.1' to '.ne.0'.
            call phas(i1,i2,i3,isucc,2,iphas)   ! 140223 Lei Added iphas.
            if( isucc.eq.0 ) goto 102   ! fail
        endif   !
c070223 

c       Proceed for success
        imes = imes + 1
        nme  = nme  + 1

c       Give proper variables to the primary meson.
        nnol = nn
        nn   = nn+1
        kn(nn,1) = 1
        kn(nn,2) = kfii
        kn(nn,3) = 0
        pn(nn,5) = amasi
        pn(nn,1) = sump1
        pn(nn,2) = sump2
        pn(nn,3) = sump3
        pnnm  = sump1*sump1 + sump2*sump2 + sump3*sump3
        pnnmm = amasi*amasi + pnnm
        if( pnnmm.gt.1D20 ) pnnmm = 1D20
        if( pnnmm.le.0D0  ) pnnmm = 1D-20
        pnnn = SQRT(pnnmm)
        pn(nn,4)   = pnnn
        ! dele     = sume - pnnn
        throe_p(4) = throe_p(4) + sume-pnnn   ! 300623 Lei

c       Produced hadron is set in between contituent partons randomly.
        pyrx = PYR(1)
        pyry = PYR(1)
        rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1)
        rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2)
        rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3)

c       Move parton list one step downward since i2+1.
111     call updad_pyj(n1,i2+1,1)
        n1 = n1 - 1
        N  = N  - 1   ! 300623 Lei

c       Move parton list one step downward since i1+1.
        call updad_pyj(n1,i1+1,1)
        nqb = nqb - 1
        n1  = n1  - 1
        N   = N   - 1   ! 300623 Lei

c300623 Share the surplus 4-momentum in throe_p.   ! 300623 Lei
        ! call share_p_PYJETS
        call share_p_PYJETS_sa1h
        if( nme.eq.1 ) return
102     continue

100     isu=0   ! 300105


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conse_c(pp,ps,npl,np,nstep)
c       Adjust four momentum conservation by iteration,no more than
c        5000 iterations
c       pp : four momentum of particle
c       ps : the four momentum should be conserved to
c       npl : order number of the first particle
c       np : order number of last particle
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
        subroutine deexcitation_E(ii,kf0,nstep,iflav)
c300623 Renamed from orginal 'ffm' and some modifications were made.   ! 300623 Lei
c       qqbar pair generation according to energy conservation ! 200223
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order number of source quark (or anti-quark)
c       kf0 : flavor code of source quark (or anti-quark)
c       nstep : number of deexcitation steps per source q (qbar)
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is anti-quark (kf0<0)
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
        common/sa18/i_deex,n_deex_step,i_pT,i_pT_max,a_FF,aPS_c,aPS_b   ! 280823 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = N

        do i1=1,3
            rc(i1)  = V(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = P(ii,i2)  ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qbar) 260223
        pm0    = P(ii,5) ! m, mass of source q (qbar), maybe 0 from PYTHIA.
        ! pm0    = PYMASS(kf0) ! m, mass of source q (qbar)

        nstep = 0   ! Counts number of deexcitation step.

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
        ! i_pT_max = 0   ! 280823 Lei
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
            rr(i)    = PYR(1)*0.5
            V(N+1,i) = rc(i) + rr(i)
            if(PYR(1).gt.0.5) V(N+1,i) = rc(i) - rr(i)
            rr(i)    = PYR(1)*0.5
            V(N+2,i) = rc(i) + rr(i)
            if(PYR(1).gt.0.5) V(N+2,i) = rc(i) - rr(i)
        enddo
        V(N+1,4) = V(ii,4) ! ii: order number of sourse parton
        V(N+2,4) = V(ii,4) ! N+1 (N+2): order number of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

c       Random three momentum method.
        do i1=1,3,1
            pii = PYR(1)*p1(i1)
            P(N+1,i1) = pii
            P(N+2,i1) = p1(i1) - pii
        enddo

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            P(N+1,3) = -P(N+1,3)
            P(N+2,3) = -P(N+2,3)
        endif

c       Recalculates E.
        pn11 = P(N+1,1)   ! pnn(1,1) 280822
        pn12 = P(N+1,2)   ! pnn(1,2) 280822
        pn13 = P(N+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        P(N+1,4) = sqrt(agsq)   ! 280822
        pn21 = P(N+2,1)   ! pnn(2,1) 280822
        pn22 = P(N+2,2)   ! pnn(2,2) 280822
        pn23 = P(N+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        P(N+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - P(N+1,1) - P(N+2,1)
        d_py = p1(2) - P(N+1,2) - P(N+2,2)
        d_pz = p1(3) - P(N+1,3) - P(N+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - P(N+1,3) - P(N+2,3)
        d_pE = p1(4) - P(N+1,4) - P(N+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        P(N+1,5) = amq   ! Mass
        P(N+2,5) = amq   ! Mass
        K(N+1,1) = 2     ! 'A'
        K(N+2,1) = 1     ! 'V'
        K(N+1,2) = kf1   ! KF code
        K(N+2,2) = kf2   ! KF code
        K(N+1,3) = ii    ! Mother
        K(N+2,3) = ii    ! Mother
        K(N+1,4) = 0     ! First daughter
        K(N+2,4) = 0     ! First aughter
        K(N+1,5) = 0     ! Last daugter
        K(N+2,5) = 0     ! Last daugter

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
        nstep = nstep + 1

        N = N + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if( nstep.ge.n_deex_step ) goto 106   ! Stop generation 280822
        if( e0.le.adj17) goto 106   ! Stop generation

        goto 100

106     continue
        if( nstep.eq.0 ) return
c       Update four momentum of the remnant of ii-th source q (qbar) 280822
        do i3=1,4   ! 280822
            P(ii,i3) = p0(i3)
        enddo

c       Re-calculate E for the remenant, let its inv. mass >= 0.
        sm2 = P(ii,4)**2 - P(ii,1)**2 - P(ii,2)**2 - P(ii,3)**2
        if( sm2.lt.0. )then
        P(ii,4) = SQRT(P(ii,5)**2 + P(ii,1)**2 + P(ii,2)**2 +P(ii,3)**2)
        P(ii,4) = P(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - P(ii,4)
        end if


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
        subroutine deexcitation_EP(ii,kf0,nstep,iflav)
c300623 Rename from orginal 'ffm' and some modifications were made.   ! 300623 Lei
c       qqbar pair generation according to light-cone variable !
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order number of source quark (or anti-quark)
c       kf0 : flavor code of source quark (or anti-quark)
c       nstep : number of deexcitation steps per source q (qbar)
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is anti-quark (kf0<0)
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
        common/sa18/i_deex,n_deex_step,i_pT,i_pT_max,a_FF,aPS_c,aPS_b   ! 280823 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = N

        do i1=1,3
            rc(i1)  = V(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = P(ii,i2)   ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qbar)   260223
        pm0    = P(ii,5) ! m, mass of source q (qbar), maybe 0 from PYTHIA.
        ! pm0    = PYMASS(kf0) ! m, mass of source q (qbar)

        if(p00(3).lt.0.) p0(3) = -p0(3)  ! Converts E-p_z as E+p_z
        w0 = e0 + p0(3) ! E+p_z

        nstep = 0   ! Counts number of deexcitation step.

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
        ! i_pT_max = 0   ! 280823 Lei
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
            rr(i)    = PYR(1)*0.5
            V(N+1,i) = rc(i) + rr(i)
            if(PYR(1).gt.0.5) V(N+1,i) = rc(i) - rr(i)
            rr(i)    = PYR(1)*0.5
            V(N+2,i) = rc(i) + rr(i)
            if(PYR(1).gt.0.5) V(N+2,i) = rc(i) - rr(i)
        enddo
        V(N+1,4) = V(ii,4) ! ii: order number of sourse parton
        V(N+2,4) = V(ii,4) ! N+1 (N+2): order number of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

c       Random three momentum method.
        do i1=1,3,1
            pii = PYR(1)*p1(i1)
            P(N+1,i1) = pii
            P(N+2,i1) = p1(i1) - pii
        enddo

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            P(N+1,3) = -P(N+1,3)
            P(N+2,3) = -P(N+2,3)
        endif

c       Recalculates E.
        pn11 = P(N+1,1)   ! pnn(1,1) 280822
        pn12 = P(N+1,2)   ! pnn(1,2) 280822
        pn13 = P(N+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        P(N+1,4) = sqrt(agsq)   ! 280822
        pn21 = P(N+2,1)   ! pnn(2,1) 280822
        pn22 = P(N+2,2)   ! pnn(2,2) 280822
        pn23 = P(N+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        P(N+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - P(N+1,1) - P(N+2,1)
        d_py = p1(2) - P(N+1,2) - P(N+2,2)
        d_pz = p1(3) - P(N+1,3) - P(N+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - P(N+1,3) - P(N+2,3)
        d_pE = p1(4) - P(N+1,4) - P(N+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        P(N+1,5) = amq   ! Mass
        P(N+2,5) = amq   ! Mass
        K(N+1,1) = 2     ! 'A'
        K(N+2,1) = 1     ! 'V'
        K(N+1,2) = kf1   ! KF code
        K(N+2,2) = kf2   ! KF code
        K(N+1,3) = ii    ! Mother
        K(N+2,3) = ii    ! Mother
        K(N+1,4) = 0     ! First daughter
        K(N+2,4) = 0     ! First aughter
        K(N+1,5) = 0     ! Last daugter
        K(N+2,5) = 0     ! Last daugter

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
        nstep = nstep + 1

        N = N + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if( nstep.ge.n_deex_step ) goto 106   ! Stop generation 280822
        if( e0.le.adj17 ) goto 106   ! Stop generation
        if( w0.le.0. )    goto 106   ! Stop generation

        goto 100

106     continue
        if( nstep.eq.0 ) return
c       Update four momentum of the remnant of ii-th source q (qbar). 280822
        do i3=1,4   ! 280822
            P(ii,i3) = p0(i3)
        enddo
        if(p00(3).lt.0.) P(ii,3) = -p0(3)   ! For the negative direction p_z.

c       Re-calculates E for the remenant, let its inv. mass >= 0.
        sm2 = P(ii,4)**2 - P(ii,1)**2 - P(ii,2)**2 - P(ii,3)**2
        if( sm2.lt.0. )then
        P(ii,4) = SQRT(P(ii,5)**2 + P(ii,1)**2 + P(ii,2)**2 +P(ii,3)**2)
        P(ii,4) = P(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - P(ii,4)
        end if


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
        subroutine deexcitation_EP_comp_pT(ii,kf0,nstep,iflav)
c300623 Rename from orginal 'ffm' and some modifications were made.   ! 300623 Lei
c       qqbar pair generation according to light-cone variable !
c      Assuming local pT compensation, i.e. px(-px) and py(-py) for q(qbar), 
c       vice versa. Sample z for qqbar.
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order number of source quark (or anti-quark)
c       kf0 : flavor code of source quark (or anti-quark)
c       nstep : number of deexcitation steps per source q (qbar)
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is anti-quark (kf0<0)
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
        common/sa18/i_deex,n_deex_step,i_pT,i_pT_max,a_FF,aPS_c,aPS_b   ! 280823 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = N

        do i1=1,3
            rc(i1)  = V(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = P(ii,i2)   ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qbar)   260223
        pm0    = P(ii,5) ! m, mass of source q (qbar), maybe 0 from PYTHIA.
        ! pm0    = PYMASS(kf0) ! m, mass of source q (qbar)

        if(p00(3).lt.0.) p0(3) = -p0(3)  ! Converts E-p_z as E+p_z
        w0 = e0 + p0(3) ! E+p_z

        nstep = 0   ! Counts number of deexcitation step.

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
        ! i_pT_max = 0   ! 280823 Lei
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
            rr(i)    = PYR(1)*0.5
            V(N+1,i) = rc(i) + rr(i)
            if(PYR(1).gt.0.5) V(N+1,i) = rc(i) - rr(i)
            rr(i)    = PYR(1)*0.5
            V(N+2,i) = rc(i) + rr(i)
            if(PYR(1).gt.0.5) V(N+2,i) = rc(i) - rr(i)
        enddo
        V(N+1,4) = V(ii,4) ! ii: order number of sourse parton
        V(N+2,4) = V(ii,4) ! N+1 (N+2): order number of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

c       Random three momentum method.
        do i1=1,3,1
            pii = PYR(1)*p1(i1)
            P(N+1,i1) = pii
            P(N+2,i1) = p1(i1) - pii
        enddo

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            P(N+1,3) = -P(N+1,3)
            P(N+2,3) = -P(N+2,3)
        endif

c       Local pT compensation.
        P(N+1,1) = px_1
        P(N+1,2) = py_1
        P(N+2,1) = px_2
        P(N+2,2) = py_2

c       Recalculates E.
        pn11 = P(N+1,1)   ! pnn(1,1) 280822
        pn12 = P(N+1,2)   ! pnn(1,2) 280822
        pn13 = P(N+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        P(N+1,4) = sqrt(agsq)   ! 280822
        pn21 = P(N+2,1)   ! pnn(2,1) 280822
        pn22 = P(N+2,2)   ! pnn(2,2) 280822
        pn23 = P(N+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        P(N+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - P(N+1,1) - P(N+2,1)
        d_py = p1(2) - P(N+1,2) - P(N+2,2)
        d_pz = p1(3) - P(N+1,3) - P(N+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - P(N+1,3) - P(N+2,3)
        d_pE = p1(4) - P(N+1,4) - P(N+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        P(N+1,5) = amq   ! Mass
        P(N+2,5) = amq   ! Mass
        K(N+1,1) = 2     ! 'A'
        K(N+2,1) = 1     ! 'V'
        K(N+1,2) = kf1   ! KF code
        K(N+2,2) = kf2   ! KF code
        K(N+1,3) = ii    ! Mother
        K(N+2,3) = ii    ! Mother
        K(N+1,4) = 0     ! First daughter
        K(N+2,4) = 0     ! First aughter
        K(N+1,5) = 0     ! Last daugter
        K(N+2,5) = 0     ! Last daugter

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
        nstep = nstep + 1

        N = N + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if( nstep.ge.n_deex_step ) goto 106   ! Stop generation 280822
        if( e0.le.adj17 ) goto 106   ! Stop generation
        if( w0.le.0. )    goto 106   ! Stop generation

        goto 100

106     continue
        if( nstep.eq.0 ) return
c       Update four momentum of the remnant of ii-th source q (qbar). 280822
        do i3=1,4   ! 280822
            P(ii,i3) = p0(i3)
        enddo
        if(p00(3).lt.0.) P(ii,3) = -p0(3)   ! For the negative direction p_z.

c       Re-calculates E for the remenant, let its inv. mass >= 0.
        sm2 = P(ii,4)**2 - P(ii,1)**2 - P(ii,2)**2 - P(ii,3)**2
        if( sm2.lt.0. )then
        P(ii,4) = SQRT(P(ii,5)**2 + P(ii,1)**2 + P(ii,2)**2 +P(ii,3)**2)
        P(ii,4) = P(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - P(ii,4)
        end if


        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
        subroutine deexcitation_E_comp_pT(ii,kf0,nstep,iflav)
c300623 Renamed from orginal 'ffm' and some modifications were made.   ! 300623 Lei
c       qqbar pair generation according to energy conservation ! 200223
c      Assuming local pT compensation, i.e. px(-px) and py(-py) for q(qbar), 
c       vice versa. Sample z for qqbar.
c280822 i.e. energetic q (qbar) de-exciatation
c       ii : the order number of source quark (or anti-quark)
c       kf0 : flavor code of source quark (or anti-quark)
c       nstep : number of deexcitation steps per source q (qbar)
c       iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is anti-quark (kf0<0)
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
        common/sa18/i_deex,n_deex_step,i_pT,i_pT_max,a_FF,aPS_c,aPS_b   ! 280823 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c   peo(5),pdec(20,5)   ! 090922

        adj17 = adj1(17)   ! Threshold energy
        i_z   = INT( adj1(29) )   ! Function for selecting z
        n0    = N

        do i1=1,3
            rc(i1)  = V(ii,i1)   ! Three-coordinate of source q (qbar)
        enddo
        do i2=1,4
            p0(i2)  = P(ii,i2)  ! Four-momentum of source q (qbar)
            p00(i2) = p0(i2)
        enddo
        kf00   = kf0     ! KF code of source q (qbar)
        e0     = p0(4)   ! E, energy of source q (qbar) 260223
        pm0    = P(ii,5) ! m, mass of source q (qbar), maybe 0 from PYTHIA.
        ! pm0    = PYMASS(kf0) ! m, mass of source q (qbar)

        nstep = 0   ! Counts number of deexcitation step.

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
        ! i_pT_max = 0   ! 280823 Lei
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
            rr(i)    = PYR(1)*0.5
            V(N+1,i) = rc(i) + rr(i)
            if(PYR(1).gt.0.5) V(N+1,i) = rc(i) - rr(i)
            rr(i)    = PYR(1)*0.5
            V(N+2,i) = rc(i) + rr(i)
            if(PYR(1).gt.0.5) V(N+2,i) = rc(i) - rr(i)
        enddo
        V(N+1,4) = V(ii,4) ! ii: order number of sourse parton
        V(N+2,4) = V(ii,4) ! N+1 (N+2): order number of new generated parton

c       Give four momentum to generated q (qbar).   ! 090922
        amq = 0.5*amasi   ! Mass

c       Random three momentum method.
        do i1=1,3,1
            pii = PYR(1)*p1(i1)
            P(N+1,i1) = pii
            P(N+2,i1) = p1(i1) - pii
        enddo

c       For the negative direction p_z.
        if(p00(3).lt.0.)then
            P(N+1,3) = -P(N+1,3)
            P(N+2,3) = -P(N+2,3)
        endif

c       Local pT compensation.
        P(N+1,1) = px_1
        P(N+1,2) = py_1
        P(N+2,1) = px_2
        P(N+2,2) = py_2

c       Recalculates E.
        pn11 = P(N+1,1)   ! pnn(1,1) 280822
        pn12 = P(N+1,2)   ! pnn(1,2) 280822
        pn13 = P(N+1,3)   ! pnn(1,3) 280822
        agsq = amq*amq + pn11*pn11 + pn12*pn12 + pn13*pn13   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        P(N+1,4) = sqrt(agsq)   ! 280822
        pn21 = P(N+2,1)   ! pnn(2,1) 280822
        pn22 = P(N+2,2)   ! pnn(2,2) 280822
        pn23 = P(N+2,3)   ! pnn(2,3) 280822
        agsq = amq*amq + pn21*pn21 + pn22*pn22 + pn23*pn23   ! 280822
        ! if(agsq.lt.1.d-20) agsq = 1.d-20   ! 280822
        P(N+2,4) = sqrt(agsq)   ! 280822

300     continue
        d_px = p1(1) - P(N+1,1) - P(N+2,1)
        d_py = p1(2) - P(N+1,2) - P(N+2,2)
        d_pz = p1(3) - P(N+1,3) - P(N+2,3)
        if(p00(3).lt.0.) d_pz = -p1(3) - P(N+1,3) - P(N+2,3)
        d_pE = p1(4) - P(N+1,4) - P(N+2,4)
        throe_p(1) = throe_p(1) + d_px
        throe_p(2) = throe_p(2) + d_py
        throe_p(3) = throe_p(3) + d_pz
        throe_p(4) = throe_p(4) + d_pE

c       Give other properties to generated q and qbar.
        P(N+1,5) = amq   ! Mass
        P(N+2,5) = amq   ! Mass
        K(N+1,1) = 2     ! 'A'
        K(N+2,1) = 1     ! 'V'
        K(N+1,2) = kf1   ! KF code
        K(N+2,2) = kf2   ! KF code
        K(N+1,3) = ii    ! Mother
        K(N+2,3) = ii    ! Mother
        K(N+1,4) = 0     ! First daughter
        K(N+2,4) = 0     ! First aughter
        K(N+1,5) = 0     ! Last daugter
        K(N+2,5) = 0     ! Last daugter

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
        nstep = nstep + 1

        N = N + 2
        ! sm2 = e0**2 - p0(1)**2 - p0(2)**2 - p0(3)**2
        ! if( SQRT(sm2).gt.pm0 ) goto 100 ! Continue to another generation.
        if( nstep.ge.n_deex_step ) goto 106   ! Stop generation 280822
        if( e0.le.adj17) goto 106   ! Stop generation

        goto 100

106     continue
        if( nstep.eq.0 ) return
c       Update four momentum of the remnant of ii-th source q (qbar) 280822
        do i3=1,4   ! 280822
            P(ii,i3) = p0(i3)
        enddo

c       Re-calculate E for the remenant, let its inv. mass >= 0.
        sm2 = P(ii,4)**2 - P(ii,1)**2 - P(ii,2)**2 - P(ii,3)**2
        if( sm2.lt.0. )then
        P(ii,4) = SQRT(P(ii,5)**2 + P(ii,1)**2 + P(ii,2)**2 +P(ii,3)**2)
        P(ii,4) = P(ii,4) + 1D-10   ! Give small sigma for machine precision.
        throe_p(4) = throe_p(4) + p0(4) - P(ii,4)
        end if


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine funcz(zz,ij)   ! 081222
c300623 This subroutine has been rewrote and fixed.   ! 300623 Lei
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
        common/sa18/i_deex,n_deex_step,i_pT,i_pT_max,a_FF,aPS_c,aPS_b   ! 280823 Lei
        common/local_fmax/ fmax_value(5)   ! 300623 Lei
        common/sa24/adj1(40),nnstop,non24,zstop

c       Sets parameters.
        a   = adj1(6)   ! Lund a parameter
        b   = adj1(7)   ! Lund b parameter
        aFF = a_FF      ! Field-Feynman parameter
        c   = aPS_c     ! Peterson/SLAC parameter
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

100     ran1 = PYR(1)
        ran2 = PYR(1)
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
        ri1(1)=V(i1,1)
        ri1(2)=V(i1,2)
        ri1(3)=V(i1,3)
        ri2(1)=V(i2,1)
        ri2(2)=V(i2,2)
        ri2(3)=V(i2,3)
        ri3(1)=V(i3,1)
        ri3(2)=V(i3,2)
        ri3(3)=V(i3,3)
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
        pi1(1)=P(i1,1)
        pi1(2)=P(i1,2)
        pi1(3)=P(i1,3)
        pi2(1)=P(i2,1)
        pi2(2)=P(i2,2)
        pi2(3)=P(i2,3)
        pi3(1)=P(i3,1)
        pi3(2)=P(i3,2)
        pi3(3)=P(i3,3)
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
        ri1(1)=V(i1,1)
        ri1(2)=V(i1,2)
        ri1(3)=V(i1,3)
        ri2(1)=V(i2,1)
        ri2(2)=V(i2,2)
        ri2(3)=V(i2,3)
        ri121=ri1(1)-ri2(1)
        ri122=ri1(2)-ri2(2)
        ri123=ri1(3)-ri2(3)
        delr=sqrt(ri121*ri121+ri122*ri122+ri123*ri123)
        pi1(1)=P(i1,1)
        pi1(2)=P(i1,2)
        pi1(3)=P(i1,3)
        pi2(1)=P(i2,1)
        pi2(2)=P(i2,2)
        pi2(3)=P(i2,3)
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
c220823 Rewritten to use probability ratios directly.   ! 220823 Lei
c       sample flavor (mass) of generated qqbar pair
c       eg: energy of original q (qbar) or g   ! sa 111123
c       kf (amq): flavor code (mass) of generated quark
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 250823 Lei
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc   ! 081222
        common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa38/ prob_ratio_q(6), am(6), amqq(6)   ! 290823 Lei sa 111123
csa     prob_ratio_q(i): ratio probability of quark with kf=i  
csa     In PYTHIA the ratio probability of d,u,s,c is equal to 
csa      1,1,0.3,10^{11}, respectively.

        dimension prob_ratio_sum(6)   ! sa 111123
c       Constituent quark mass.
        do i=1,6,1
           am(i) = pymass(i) 
        end do
csa     Note: constituent quark mass is assumed in other parts of PACIAE3.0 . 
c       Mass of qqbar pair.   ! sa
            amqq = 2*am

c       Samples flavor of generated quark pair according to prob_ratio_q.
csa     General generated flavor candidates: kf=d,u,s,c,b,t.
        rand_num = PYR(1)
        amdd = amqq(1)
        amss = amqq(3)
csa     If energy ('eg') is not enough to excite dd_bar pair        
        if( eg .lt. amdd ) return

        if( eg .lt. amss)then   ! d and u only
          if(rand_num.le.0.5)then
          kf=1   ! d
          amq = am(1)
          else
          kf=2   ! u
          amq = am(2)
          endif
        endif  

        prob_ratio_tot = prob_ratio_q(1)+prob_ratio_q(2)
        do i=3,6,1
csa     Sum of ratio probability upto i ($\sum_{j=1}^{j=i}$).
        prob_ratio_tot = prob_ratio_tot + prob_ratio_q(i)
            if( eg .ge. amqq(i) )then
                n_flavor = i
csa     Possible candidate with kf<=n_flavor at a given 'eg'.
csa     Sum of those possible candidate ratio probability.
                prob_ratio_sum(i)=prob_ratio_tot   ! sa
            end if
        end do
        prob_interval_low = 0D0
        prob_interval_upp = 0D0
csa     Calculate probability value at upper end point of each probability 
c        interval and determine simulated flavor code.
        do j=1,n_flavor,1   
csa     Note: possible candidate with kf<=n_flavor at a given energy 'eg'.  
            ratio_interval = prob_ratio_q(j) / prob_ratio_sum(n_flavor) ! sa
            prob_interval_upp = prob_interval_upp + ratio_interval
            if( rand_num.gt.prob_interval_low .AND.
     &          rand_num.le.prob_interval_upp )then
                kf = j
                amq = am(j)
                return
            end if
            prob_interval_low = prob_interval_upp
        end do

csa     For example:         
c        in the case of ratio probability is u:d:s:c:b:t=1:1:0.5:*:*:* 
c        and has 'eg.ge.amqq(3)', the possible candidates are d,u,and s,
c        the corresponding probability intervals are [0,0.4), (0.4,0.80), 
c        and (0.8,1].
c        Ratio: u:d:s=          1          :         1          :   0.5     
c                     |____________________|____________________|__________|
c        Interval:    0                   0.4                  0.8         1
c        0.4=1/sum (sum=1+1+0.5); 0.8=(1+1)/sum.
csa

c       Mass defined in PYTHIA, in GeV.
c       Kinematical mass:
c           amd = 0.33D0
c           amu = 0.33D0   ! amu=amd
c           ams = 0.5D0
c           amc = 1.5D0
c           amb = 4.8D0
c           amt = 175D0
c       Current algebra mass
c           amd = 0.0099D0
c           amu = 0.0056D0   ! amu < amd
c           ams = 0.199D0
c           amc = 1.23D0
c           amb = 4.17D0
c           amt = 165D0
c       Constituent mass
c           amd = 0.333   ! sa
c           amu = 0.333   ! amu=amd
c           ams = 0.5D0
c           amc = 1.5D0   ! sa
c           amb = 4.8D0   ! sa
c           amt = 174.4   ! sa

        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tdgaus(v,pmax,np,pp)
c...... 2-d Gaussian distribution with width v, i.e., e^(-p2/v)dp2, 0<p2<pmax
c...... set pmax < 0 if pmax should be infinity.
c...... np : the total number of particles wanted to sample their transverse momenta
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
10        p2 = -v*log(max(1.e-20,1. - a*PYR(1)))
          if(p2.LT.0.)goto 10
          ps=sqrt(p2)
          fi=2.*3.1415926*PYR(1)
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
        if(iii.eq.1)then   ! 300623 Lei
        kglu=0
        pglu=0.
        vglu=0.
        endif   ! 300623 Lei
        jb=0
201     do i1=jb+1,N
        kf=K(i1,2)
        kfab=iabs(kf)
        eng=P(i1,4)
        if(kfab.ne.21)then   ! stay
        jb=jb+1
        goto 202
        endif

        nglu=nglu+1
        do i2=1,5
        kglu(nglu,i2)=K(i1,i2)
        pglu(nglu,i2)=P(i1,i2)
        vglu(nglu,i2)=V(i1,i2)
        enddo
        if(i1.eq.N)then
        N=N-1
        goto 203
        endif
c       move particle list one step downward from i1+1 to N
        do jj=1,5
        do j=i1+1,N
        K(j-1,jj)=K(j,jj)
        P(j-1,jj)=P(j,jj)
        V(j-1,jj)=V(j,jj)
        enddo
        enddo
        N=N-1
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
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 250823 Lei
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 300623 Lei
        common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 080104 220110
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) 
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)

       amu = pymass(2)   ! 271022

        amuu=2*amu   ! 271022
c       throw away g with energy<amuu
        jb=1       ! 300623 Lei
2010    continue   ! 300623 Lei
        do i1=jb,nglu,1   ! 300623 Lei 1 -> jb
        eg=pglu(i1,4)
        if(eg.lt.amuu)then
        do i2=1,4
        throe_p(i2)=throe_p(i2)+pglu(i1,i2)   ! 300623 Lei Replace ppsa by throe_p.
        enddo
c300623 Lei
        if(i1.eq.nglu)then
            nglu = nglu - 1
            goto 2020
        endif
c300623 Lei
c       move particle list ('sa36') one step downward from i1+1 to nglu
        do jj=1,5
        do j=i1+1,nglu
        kglu(j-1,jj)=kglu(j,jj)
        pglu(j-1,jj)=pglu(j,jj)
        vglu(j-1,jj)=vglu(j,jj)
        enddo
        enddo
        nglu=nglu-1
        jb=i1       ! 300623 Lei
        goto 2010   ! 300623 Lei
        endif
        enddo
2020    continue    ! 300623 Lei

c       g (in 'sa36') -> qq_{bar} (as a string filling in 'pyjets')
100     do i1=1,nglu   ! do loop over gluons
        eg=pglu(i1,4)
        call break_f(eg,kf,amq)
        kf1=kf
        kf2=-kf
        am1=amq
        am2=amq

        K(N+1,1)=2   ! 'A'
        K(N+2,1)=1   ! 'V'
        K(N+1,2)=kf1
        K(N+2,2)=kf2
        K(N+1,3)=0
        K(N+2,3)=0
        K(N+1,4)=0
        K(N+2,4)=0
        K(N+1,5)=0
        K(N+2,5)=0
c       P(N+1,5)=am1
c       P(N+2,5)=am2

c       give four momentum to the breaked quarks
        call bream_glu(i1,kf1,kf2)
c       give four coordinate to the breaked quarks
        call coord_glu(i1)
        if(i1.eq.nglu)then
        N=N+2
        goto 200
        endif
        N=N+2
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
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 250823 Lei
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104 300623 Lei
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
        goto 400   ! for 'decay method'
401     do i1=1,3   ! for 'random three momentum method'
        pi(i1)=PYR(1)*pglu(ii,i1)
        pp(1,i1)=pi(i1)
        pp(2,i1)=ps(i1)-pi(i1)
        enddo
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
        decsuc=1
        call decmom(ps,pp,am1,am2,decsuc)
        if(decsuc.eq.0)goto 401   ! return to random three momentum method
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
        P(N+1,i1)=pp(1,i1)
        enddo
        P(N+1,5)=am1
        do i1=1,4
        P(N+2,i1)=pp(2,i1)
        enddo
        P(N+2,5)=am2

        do i2=1,4   ! 300623 Lei Collects lost 4-momentum.
        throe_p(i2) = throe_p(i2) + ( ps(i2) - P(N+1,i2) - P(N+2,i2) )
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
c       ii: order number of gluon in 'sa36'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) 
        dimension rr(3)

        do i1=1,3
        V(N+1,i1)=vglu(ii,i1)
        rr(i1)=PYR(1)*0.5   ! 261002
        V(N+2,i1)=vglu(ii,i1)+rr(i1)
        if(PYR(1).gt.0.5d0)V(N+2,i1)=vglu(ii,i1)-rr(i1)
        enddo
        V(N+1,4)=vglu(ii,4)
        V(N+2,4)=vglu(ii,4)
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
        ich1=ich1+PYCHGE(kf)
        enddo
        cc=ich1/3.
        write(22,*)'sa36 nn=',nn
        write(mstu(11),*)'c & p sum=',cc,peo   !
c       write(9,*)peo,ich1/3   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa37(nn,cc)   ! 300623 Lei
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
        ich1=ich1+PYCHGE(kf)
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
        pt=-parj(21)*log(max(1d-10,PYR(1)))
        phi=paru(2)*PYR(1)
c       randomly sample [px,py] on circle of sphere with radius pt
        px=pt*cos(phi)*(1+smadel)   ! 300623 Lei
        py=pt*sin(phi)*(1-smadel)   ! 300623 Lei
        pt=SQRT(px*px+py*py)        ! 300623 Lei
        return
        end



C*********************************************************************
C...PAPTDI
C...Generates transverse momentum according to a Gaussian/Exponential/Random.
 
        SUBROUTINE PAPTDI(KF0,KF1,KF2,PX0,PY0,PT0,PX,PY,PT,i_pT,i_max,x)
c300623 Lei
c       i_pT: (D=1) the pT sampling method of the daughter qqbar pair in coal
c             = 1, Gaussian px and py with width PARJ(21)
c             = 2, Exponential px and py with width PARJ(21)
c             = 3, Exponential pT with width PARJ(21)
c             = 4, random pT from mother
c             = 5, random px and random py from mother, different random factors
c             = 6, random (px and py) from mother, the same random factor
c             = 7, random (px and py) from mother, the same random factor as 
c                   z which related to adj1(29).
c       i_max: (D=0) whether the sampled pT in coal deexitation is greater 
c                       than the mother quark or not.
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
            if(PX0.gt.0. .AND. PX.gt.PX0) goto 100
            if(PX0.lt.0. .AND. PX.lt.PX0) goto 100
            if(PY0.gt.0. .AND. PY.gt.PY0) goto 100
            if(PY0.lt.0. .AND. PY.lt.PY0) goto 100   ! 060923 Lei fixed .gt. -> .lt., 200 -> 100
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
