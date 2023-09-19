        subroutine analy(nmin,nminf,ncha,nchaf,yOrEta)   ! 140223 Lei
c       analyses an event based on the messages in 'pyjets'
c       it is composed by Ben-Hao Sa on 28/12/19
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(6),
     c   afl(20,6,2)   ! 300623 Lei 5 -> 6
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     c   iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph
     c   ,vnlep   ! 260314
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003
        common/sa30/vneump,vneumt,mstptj
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
        common/anly1/an(40,6,20),bn(20),anf(40,6,20),bnf(20)   ! 300623 Lei 5 -> 6
        real*8 nmin,nminf,ncha,nchaf
        dimension c(5),dinel(600),dineli(600),sthroe(4),wthroe(4)

c       ispmax: maximum # of particle KF code wanted to statistics
c       ispkf: array of particle KF code 
c       ispkf(i): KF code of i-th particle
c       kfmax: dimention of ispkf array, i.e. maximum # of KF code considered
c       bn(i) (bnf(i)): multiplicity of i-th particle in partial (full)
c        phase space
c       an(l,i,j) (anf(l,i,j)): three dimension array in partial (full)
c         phase space
c       l: order # of y (pt, eta, ...) interval
c       i: kind of distribution 
c       for pp, NA,AN and AB collisions
c        i=1: for y
c        i=2: for pt (inv. pT)
c        i=3: for eta
c        i=4: for mT
c        i=5: for event-wise multiplicity
c        i=6: for pT (dN/dpT)   ! 300623 Lei
c       for lp and lA collisions
c        i=1 : z
c        i=2 : \nu
c        i=3 : Q^2
c         .        .
c         .        .
c         .        .
c       j: order # of particle in array of 'ispkf' (1 to ispmax)
c       isdmax: maximum # of distributions considered
c       asd(i): interval segmented for i-th distribution
c       iflmax: maximum # of windows, =0 means no window at all
c       afl(j,i,1): lower bound of i-th window for j-th particle
c       afl(j,i,2): upper bound of i-th window for j-th particle
c260314 for pp,NA,AN and AB collisions
c        i=1: y/eta window (depends on yOrEta)
c        i=2: pt window
c         .        .
c         .        .
c         .        .
c260314 for lp and lA collisions
c        i=1 : Q^2=-q^2 (fq2 named in program) window
c        i=2 : W^2 (w2l) window
c        i=3 : y (yyl) window
c        i=4 : P_h (pph) window
c        i=5 : z (zl) window
c         .        .
c         .        .
c         .        .

c       for lp and lA collisions
c       pincl (pscal): four momentum and mass of incident (scatterd) lepon
c       pinch: four momentum and mass of incident hadron
c       vnu: \nu; fq2: Q^2=-q^2; w2l: W^2; yyl: y; zl: z; xb: x_B; pph: P_h

c140223 Lei
c       yOrEta: select y or eta in partial phase-space statistics.
c               = 0 , y
c               = 1 , eta
c140223 Lei

        dpmax=adj1(27)   ! largest momentum allowed for particle
        adj140=adj1(40)

c       initializes the variales
        nmin=0.
        nminf=0.
        ncha=0.
        nchaf=0.
        bn=0.
        bnf=0.
        an=0.
        anf=0.

c       analyses an event
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
        ppm=dsqrt(p1*p1+p2*p2+p3*p3)
        if(ppm.le.dpmax.and.p4.le.dpmax)then
        goto 3000
        else
        goto 400   ! throw away that particle
        endif
c281104
3000    continue
c060718 if((itden.eq.0.and.ipden.eq.1).or.(itden.eq.1.and.ipden.eq.0)
c060718     c   .or.(itden.eq.1.and.ipden.eq.1))then   ! 260314
        if(ipden.lt.11)then   !  260314 060718 for pp,pA (Ap) & AB
        c(1)=yy
        if( INT(yOrEta).eq.1 ) c(1)=eta   ! 140223 Lei
        c(2)=ppt
c       .
c       .
c       .
        kkk=1
c033101
c       statistics of negative multiplicity
        if(adj140.ge.3 .and. plu6.lt.-0.9)then   ! for hadron 140414
        nminf=nminf+1.
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 700
        enddo
        nmin=nmin+1.
700     endif
        if(adj140.lt.3 .and. plu6.lt.-0.2)then   ! for parton 140414
        nminf=nminf+1.   ! -plu6 230206
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 702
        enddo
        nmin=nmin+1.   ! -plu6 230206
702     endif
c010220
c033101 statistics of positive multiplicity
        if(adj140.ge.3 .and. plu6.gt.0.9)then   ! for hadron 140414
        nchaf=nchaf+1.   ! +plu6
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 701
        enddo
        ncha=ncha+1.   ! +plu6
701     endif
        if(adj140.lt.3 .and. plu6.gt.0.2)then   ! for parton 140414
c070802
        nchaf=nchaf+1.   ! +plu6 230206
        do i=1,iflmax
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 703
        enddo
        ncha=ncha+1.   ! +plu6 230206
703     endif
        endif   ! 260314

c       statistics of y, pt, ... distributions (for NA,AN,and BA); z, \nu,
c        ... distributions (for lp and -lA)
        do 500 kk=1,ispmax
        kf=ispkf(kk)
        if(ik.ne.kf)goto 500

c????????????????????????????????????????????????????????????????????
c       excludes the projectile and the target spectator nucleons
        if(ik.eq.2212 .or. ik.eq.2112)then
c       eep=dsqrt(win*win+0.938*0.938)
c       if(ifram.eq.1)eep=0.5*win
c       eepd=eep-0.001
c       eepu=eep+0.001
c       if(ifram.eq.0 .and. ((p4.gt.eepd .and. p4.lt.eepu) .or.
c     c   p4.le.0.940))then   ! 111899
c       goto 500
c       endif   ! 111899
c       if(ifram.eq.1 .and. (p4.gt.eepd .and. p4.lt.eepu))then
        if(ppt.le.1.e-5)goto 500 ! 310521 e-4->e-5
        endif
c????????????????????????????????????????????????????????????????????

c       analyses for pp,pA and AB
        if(ipden.lt.11)
     &   call stati_h(yy,ppt,eta,p5,ik,kk,w,bn,an,bnf,anf,yOrEta)   ! 140223 Lei
c       analyses for lp and lA
        if(ipden.ge.11.and.ipden.le.16)
     c   call stati_l(p1,p2,p3,p4,p5,ik,kk,w,bn,an,bnf,anf)
c260314
        goto 400
500     continue
400     continue

c       statistics of multiplicity distributions,
c        spectator nucleons are excluded
        do kkk=1,ispmax
c       ik=ispkf(kkk)
c       if(iabs(ik).eq.2212.or.iabs(ik).eq.2112.or.iabs(ik).eq.3122.or.
c     &   iabs(ik).eq.3212.or.ik.eq.3222.or.ik.eq.3112)then
c       multiplicity is located at which interval
c       if(ik.eq.2212)then
        idf=bnf(kkk)/asd(5)+1
c       the 5-th distribution is particle multiplicity distribution
        if(idf.lt.1 .or. idf.gt.40)goto 405   ! 131204 300623 Lei 20 -> 40
        anf(idf,5,kkk)=anf(idf,5,kkk)+1. !/asd(5)
c       active the window 
405     do i=1,iflmax   ! 131204
        if(c(i).lt.afl(kkk,i,1) .or. c(i).gt.afl(kkk,i,2))goto 404
        enddo
        idd=bn(kkk)/asd(5)+1
        if(idd.lt.1 .or. idd.gt.40)goto 404   ! 131204 300623 Lei 20 -> 40
        an(idd,5,kkk)=an(idd,5,kkk)+1. !/asd(5)
404     continue
c       endif
        enddo


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine stati_h(y,pt,eta,p5,ik,kk,ww,a,b,af,bf,yOrEta)   ! 260314 140223 Lei
c       on line statistics for NA,AN,and AA collisions   ! 260314
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(6),
     c   afl(20,6,2)   ! 300623 Lei 5 -> 6
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        dimension a(20),b(40,6,20),c(6),af(20),bf(40,6,20),id(6) ! 070419 300623 Lei 5 -> 6
        if(pt.le.1D-15) return   ! 300623 Lei Excludes particles with too small/zero pT.
        amass=p5   ! 010600
        amass2=amass*amass
        pt2=pt*pt
        pmt=dsqrt(pt2+amass2)   ! 300623 Lei m_T
        do 10000 i=1,iflmax
        goto (10,20,30,40,50,60) i   ! 300623 Lei Added 60
10      c(i)=y
        if( INT(yOrEta).eq.1 ) c(1)=eta   ! 140223 Lei
        goto 10000
20      c(i)=pt
        goto 10000
30      continue
40      continue
50      continue
60      continue   ! 300623 Lei added 60
10000   continue
c       calculate the abscissa one by one
40000   do 20000 i=1,isdmax
        goto (100,200,300,400,500,600) i   ! 300623 Lei Added 600
c       y is located in which interval?
100     ii=dabs(y)/asd(i)+1
        if(ifram.eq.1 .and. y.gt.0.)ii=ii+20   ! 311019
        if(ifram.eq.1 .and. y.lt.0.)ii=20-ii+1   ! 311019
c       note: 20 here should be change together with the dimension
c        40   ! 311019
        id(i)=ii
c100    id(i)=y/asd(i)+1
c       if(id(i).le.0)id(i)=1
        goto 20000
c       pt is located in which interval? (inv. pT)
200     id(i)=pt/asd(i)+1
c       if(id(i).le.0)id(i)=1
        goto 20000
c       eta is located in which interval?
300     ii=dabs(eta)/asd(i)+1
        if(ifram.eq.1 .and. eta.gt.0.)ii=ii+20   ! 311019
        if(ifram.eq.1 .and. eta.lt.0.)ii=20-ii+1   ! 311019
c       note: 20 here should be change together with the dimension
c        40   ! 311019
        id(i)=ii
        goto 20000
c       mt is counted in which mt interval ?
400     id(i)=pmt/asd(i)+1   ! 300623 Lei
        goto 20000
500     continue
        goto 20000   ! 300623 Lei
c300623 pt is located in which interval? (dN/dpT)   ! 300623 Lei
600     id(i)=pt/asd(i)+1
20000   continue

c       make statistics of particle yield and desired distributions in full space
        af(kk)=af(kk)+ww
        do i=1,isdmax
        ii=id(i)
        if(ii.lt.1 .or. ii.gt.40)goto 30000   ! 010218 070419
        if(i.eq.1)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)       ! dN/dy
        if(i.eq.2)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)/pt    ! (1/pT)dN/dpT
        if(i.eq.3)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)       ! dN/deta
        if(i.eq.4)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)/pmt   ! (1/mT)dN/dmT
        if(i.eq.6)bf(ii,i,kk)=bf(ii,i,kk)+ww/asd(i)       ! 300623 Lei dN/dpT
30000   enddo

c       if(iflmax.eq.0)return

c300623 make statistics for desired distributions in partial space   ! 300623 Lei
        do i=1,isdmax
            ii=id(i)
            if(ii.lt.1 .or. ii.gt.40)goto 50000   ! 010218  070419
            if( c(2).ge.afl(kk,2,1) .AND. c(2).le.afl(kk,2,2) )then   ! pT cut
                if(i.eq.1) b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)      ! dN/dy
                if(i.eq.3) b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)      ! dN/deta
            endif
            if( c(1).ge.afl(kk,1,1) .AND. c(1).le.afl(kk,1,2) )then   ! y/eta cut
                if(i.eq.2) b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)/pt   ! (1/pT)dN/dpT
                if(i.eq.4) b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)/pmt  ! (1/mT)dN/dmT
                if(i.eq.6) b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)      ! 300623 Lei dN/dpT
            endif
50000   enddo

c       put the filter to be effective
        do i=1,iflmax
        if(c(i).lt.afl(kk,i,1) .or. c(i).gt.afl(kk,i,2))return
        enddo
c       make statistics for particle yield in partial space
        a(kk)=a(kk)+ww

        return
        end



c260314cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine stati_l(p1,p2,p3,p4,p5,ik,kk,ww,a,b,af,bf)
c       on line statistics for case of lepton incidence
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(6),
     c   afl(20,6,2)   ! 300623 Lei 5 -> 6
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph
     c   ,vnlep
        dimension a(20),b(40,6,20),af(20),bf(40,6,20),c(6),id(6) ! 070419  300623 Lei 5 -> 6
c       vnu: \nu; fq2: Q^2=-q^2; w2l: W^2; yyl: y; zl: z; xb: x_B; pph: P_h
c       calculate kinematic variable relevant to the produced hadron
        pph=p1*p1+p2*p2+p3*p3
        pph=dmax1(pph,1.d-20)
        pph=dsqrt(pph)
        zln=pinch(4)*p4-pinch(1)*p1-pinch(2)*p2-pinch(3)*p3   ! numerator of z
        zld=pinch(5)*vnu   ! denominator of z
        zld=dmax1(zld,1.d-20)
        zl=zln/zld
c       write(9,*)'ik,kk,pph,zl=',ik,kk,pph,zl
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
c       write(9,*)'c(i)=',(c(i),i=1,5)
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
c       write(9,*)'asd(i)=',(asd(i),i=1,5)
c       write(9,*)'id(i)=',(id(i),i=1,3)
c       make statistics of particle yield and desired distributions
        af(kk)=af(kk)+ww
        do i=1,isdmax
        ii=id(i)
        if(ii.lt.1 .or. ii.gt.40)goto 30000   ! 070419
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
        if(ii.lt.1 .or. ii.gt.40)goto 50000   ! 070419
        if(i.eq.1)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
        if(i.eq.2)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
        if(i.eq.3)b(ii,i,kk)=b(ii,i,kk)+ww/asd(i)
50000   enddo
        return
        end



        subroutine analy_parton(yOrEta)
c300623 Lei
c       Analyses an event based on the messages in "pyjets"
c       Use the nominal cuts of 1-st setting in usu.dat for partons.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        ! COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        ! common/oscar0/ n0, npad0, k0(1000,5), p0(1000,5), v0(1000,5)
        ! common/oscar1/ n1, npad1, k1(kszj,5), p1(kszj,5), v1(kszj,5)
        ! common/oscar2/ n2, npad2, k2(kszj,5), p2(kszj,5), v2(kszj,5)
        COMMON/oscar2/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! Note name here.
        ! common/oscar3/ n3, npad3, k3(kszj,5), p3(kszj,5), v3(kszj,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(6),
     c   afl(20,6,2)   ! 300623 Lei 5 -> 6
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
c300623 For parton.   ! 300623 Lei
c       u, ubar, d, dbar, s, sbar, c, cbar, b, bbar, t, tbar, g, (u+d+s + anti-)
        common/anly2/ sn_min_p, sn_min_p_f, sn_cha_p, sn_cha_p_f, sn_g,
     &  sn_g_f, san_p(40,6,14), sbn_p(14), san_p_f(40,6,14), sbn_p_f(14)   ! 300623 Lei 5 -> 6
c300623 Lei
        dimension c(6),bn_p(14),bn_p_f(14),an_p(40,6,14),an_p_f(40,6,14)   ! 300623 Lei 5 -> 6

c       bn(i) (bnf(i)): multiplicity of i-th particle in partial (full)
c        phase space
c       an(l,i,j) (anf(l,i,j)): three dimension array in partial (full)
c         phase space
c       l: order # of y (pt, eta, ...) interval
c       i: kind of distribution 
c       for pp, NA,AN and AB collisions
c        i=1: for y
c        i=2: for pt (inv. pT)
c        i=3: for eta
c        i=4: for mT
c        i=5: for event-wise multiplicity
c        i=6: for pT (dN/dpT)   ! 300623 Lei
c       j: order # of particle in array of 'ispkf' (1 to ispmax)
c       isdmax: maximum # of distributions considered
c       asd(i): interval segmented for i-th distribution
c       iflmax: maximum # of windows, =0 means no window at all
c       afl(j,i,1): lower bound of i-th window for j-th particle
c       afl(j,i,2): upper bound of i-th window for j-th particle
c260314 for pp,NA,AN and AB collisions
c        i=1: y/eta window (depends on yOrEta)
c        i=2: pt window
c         .        .
c         .        .
c         .        .

c140223 Lei
c       yOrEta: select y or eta in partial phase-space statistics.
c               = 0 , y
c               = 1 , eta
c140223 Lei

        dpmax = adj1(27)   ! Largest momentum allowed for particle

c       Initializes the variales
        n_min_p   = 0
        n_min_p_f = 0
        n_cha_p   = 0
        n_cha_p_f = 0
        n_g       = 0
        n_g_f     = 0
        bn_p   = 0.
        bn_p_f = 0.
        an_p   = 0.
        an_p_f = 0.
        if(iii.eq.1)then
            sn_min_p   = 0.
            sn_min_p_f = 0.
            sn_cha_p   = 0.
            sn_cha_p_f = 0.
            sn_g       = 0.
            sn_g_f     = 0.
            sbn_p   = 0.
            sbn_p_f = 0.
            san_p   = 0.
            san_p_f = 0.
        end if

c       Analyses an event
        do 400 j=1,N,1
            kf   = K(j,2)
            plu6 = PYP(j,6)
            eta  = PYP(j,19)
            w    = 1.
            p1   = P(j,1)
            p2   = P(j,2)
            p3   = P(j,3)
            p4   = P(j,4)
            p5   = P(j,5)
            yy   = PYP(j,17)
            ppt  = PYP(j,10)
            ppm  = SQRT( p1*p1 + p2*p2 + p3*p3 )
            if( ppm.le.dpmax .and. p4.le.dpmax )then
                goto 3000
            else
                goto 400   ! Throw away that particle
            endif
3000        continue
            if( ipden.lt.11 )then   ! For pp,pA (Ap) & AB
                c(1) = yy
                if( INT(yOrEta).eq.1 ) c(1) = eta
                c(2) = ppt
c               .
c               .
c               .
c               Statistics of negative multiplicity
                if( plu6.lt.-0.2 )then   ! For parton
                    n_min_p_f = n_min_p_f + 1   ! -plu6
                    do i=1,iflmax,1
            if( c(i).lt.afl(1,i,1) .or. c(i).gt.afl(1,i,2) ) goto 702
                    enddo
                    n_min_p = n_min_p + 1   ! -plu6
702             endif
c               Statistics of positive multiplicity
                if( plu6.gt.0.2 )then   ! For parton
                    n_cha_p_f = n_cha_p_f + 1   ! +plu6
                    do i=1,iflmax,1
            if( c(i).lt.afl(1,i,1) .or. c(i).gt.afl(1,i,2) ) goto 703
                    enddo
                    n_cha_p = n_cha_p + 1   ! +plu6
703             endif
c               Statistics of gluon multiplicity
                if( kf.eq.21 )then
                    n_g_f = n_g_f + 1
                    do i=1,iflmax,1
            if( c(i).lt.afl(1,i,1) .or. c(i).gt.afl(1,i,2) ) goto 704
                    enddo
                    n_g = n_g + 1
704             endif
            endif

c           g, (u+d+s + anti-), u, ubar, d, dbar, s, sbar, c, cbar, b, bbar, t, tbar
            kk = 0
            select case(kf)
            case(21)
                kk = 1
            ! case ABS( 1~3 )
            !     kk = 2
            case(2)
                kk = 3
            case(-2)
                kk = 4
            case(1)
                kk = 5
            case(-1)
                kk = 6
            case(3)
                kk = 7
            case(-3)
                kk = 8
            case(4)
                kk = 9
            case(-4)
                kk = 10
            case(5)
                kk = 11
            case(-5)
                kk = 12
            case(6)
                kk = 13
            case(-6)
                kk = 14
            end select

c           Statistics of y, pt, ... distributions (for NA,AN,and AB)

c           Analyses for pp,pA and AB
            if(ipden.lt.11 .AND. kk.gt.0) call stati_parton(yy, ppt,
     &          eta, p5, kf, kk, w, bn_p, an_p, bn_p_f, an_p_f, yOrEta)

400     continue

c       Statistics of multiplicity distributions
        do kkk=1,14,1
c           Multiplicity is located at which interval
            idf = INT( bn_p_f(kkk)/asd(5) + 1 )
c           The 5-th distribution is particle multiplicity distribution
            if(idf.lt.1 .or. idf.gt.40) goto 405
            an_p_f(idf,5,kkk) = an_p_f(idf,5,kkk) + 1. !/asd(5)
c           Active the window
c           Use the nominal cuts of 1-st setting in usu.dat for partons.
405         do i=1,iflmax,1
                if(c(i).lt.afl(1,i,1) .or. c(i).gt.afl(1,i,2)) goto 404
            enddo
            idd = INT( bn_p(kkk)/asd(5) + 1 )
            if(idd.lt.1 .or. idd.gt.40) goto 404
            an_p(idd,5,kkk) = an_p(idd,5,kkk) + 1. !/asd(5)
404         continue
        enddo

c       Event accumulation
        sn_min_p   = sn_min_p   + n_min_p * 1.
        sn_min_p_f = sn_min_p_f + n_min_p_f * 1.
        sn_cha_p   = sn_cha_p   + n_cha_p * 1.
        sn_cha_p_f = sn_cha_p_f + n_cha_p_f * 1.
        sn_g       = sn_g       + n_g * 1.
        sn_g_f     = sn_g_f     + n_g_f * 1.

        sbn_p   = sbn_p   + bn_p
        sbn_p_f = sbn_p_f + bn_p_f
        san_p   = san_p   + an_p
        san_p_f = san_p_f + an_p_f


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine stati_parton(y,pt,eta,p5,kf,kk,ww,a,b,af,bf,yOrEta)
c300623 Lei
c       Online statistics for NA,AN,and AA collisions
c       Use the nominal cuts of 1-st setting in usu.dat for partons.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(6),
     c   afl(20,6,2)   ! 300623 Lei 5 -> 6
        dimension c(6), id(6), a(14), af(14), b(40,6,14), bf(40,6,14)   ! 300623 Lei 5 -> 6

        if(pt.le.1D-15) return   ! Excludes particles with too small/zero pT.
        amass  = p5
        amass2 = amass*amass
        pt2    = pt*pt
        pmt    = sqrt( pt2 + amass2 )
        c(1) = y
        if( INT(yOrEta).eq.1 ) c(1) = eta
        c(2)=pt
c       .
c       .
c       .
c       Calculate the abscissa one by one
c       y is located in which interval?
        id = 0
        ii = INT( ABS(y)/asd(1) + 1 )
        if( y.ge.0. ) ii = ii + 20
        if( y.lt.0. ) ii = 20 - ii + 1
        id(1) = ii
c       pt is located in which interval? (inv. pT)
        id(2) = INT( pt/asd(2) + 1 )
c       eta is located in which interval?
        ii = INT( ABS(eta)/asd(3) + 1 )
        if( eta.ge.0. ) ii = ii + 20
        if( eta.lt.0. ) ii = 20 - ii + 1
        id(3) = ii
c       mt is counted in which mt interval ?
        id(4) = INT( pmt/asd(4) + 1 )
c300623 pt is located in which interval? (dN/dpT)   ! 300623 Lei
        id(6) = INT( pt/asd(6) + 1 )

c       Make statistics of particle yield and desired distributions in full space
        af(kk) = af(kk) + ww
        if( ABS(kf).ge.1 .AND. ABS(kf).le.3 ) af(2) = af(2) + ww
        do i=1,isdmax,1
            ii = id(i)
            if(ii.lt.1 .or. ii.gt.40) goto 30000
            if(i.eq.1) bf(ii,i,kk) = bf(ii,i,kk) + ww/asd(i)
            if(i.eq.2) bf(ii,i,kk) = bf(ii,i,kk) + ww/asd(i)/pt
            if(i.eq.3) bf(ii,i,kk) = bf(ii,i,kk) + ww/asd(i)
            if(i.eq.4) bf(ii,i,kk) = bf(ii,i,kk) + ww/asd(i)/pmt
            if(i.eq.6) bf(ii,i,kk) = bf(ii,i,kk) + ww/asd(i)   ! 300623 Lei dN/dpT
            if( ABS(kf).ge.1 .AND. ABS(kf).le.3 )then   ! u+d+s + anti-
                if(i.eq.1) bf(ii,i,2) = bf(ii,i,2) + ww/asd(i)
                if(i.eq.2) bf(ii,i,2) = bf(ii,i,2) + ww/asd(i)/pt
                if(i.eq.3) bf(ii,i,2) = bf(ii,i,2) + ww/asd(i)
                if(i.eq.4) bf(ii,i,2) = bf(ii,i,2) + ww/asd(i)/pmt
                if(i.eq.6) bf(ii,i,2) = bf(ii,i,2) + ww/asd(i)   ! 300623 Lei
            endif
30000   enddo

c       Make statistics for desired distributions in partial space
c       Use the nominal cuts of 1-st setting in usu.dat for partons.
        do i=1,isdmax,1
            ii=id(i)
            if(ii.lt.1 .or. ii.gt.40) goto 50000
            if( c(2).ge.afl(1,2,1) .AND. c(2).le.afl(1,2,2) )then   ! pT cut
                if(i.eq.1) b(ii,i,kk) = b(ii,i,kk) + ww/asd(i)      ! dN/dy
                if(i.eq.3) b(ii,i,kk) = b(ii,i,kk) + ww/asd(i)      ! dN/deta
                if( ABS(kf).ge.1 .AND. ABS(kf).le.3 )then
                    if(i.eq.1) b(ii,i,2) = b(ii,i,2) + ww/asd(i)
                    if(i.eq.3) b(ii,i,2) = b(ii,i,2) + ww/asd(i)
                endif
            endif
            if( c(1).ge.afl(1,1,1) .AND. c(1).le.afl(1,1,2) )then   ! y/eta cut
                if(i.eq.2) b(ii,i,kk) = b(ii,i,kk) + ww/asd(i)/pt   ! (1/pT)dN/dpT
                if(i.eq.4) b(ii,i,kk) = b(ii,i,kk) + ww/asd(i)/pmt  ! (1/mT)dN/dmT
                if(i.eq.6) b(ii,i,kk) = b(ii,i,kk) + ww/asd(i)   ! 300623 Lei dN/dpT
                if( ABS(kf).ge.1 .AND. ABS(kf).le.3 )then
                    if(i.eq.2) b(ii,i,2) = b(ii,i,2) + ww/asd(i)/pt
                    if(i.eq.4) b(ii,i,2) = b(ii,i,2) + ww/asd(i)/pmt
                    if(i.eq.6) b(ii,i,2) = b(ii,i,2) + ww/asd(i)   ! 300623 Lei
                endif
            endif
50000   enddo

c       Put the filter to be effective
        do i=1,iflmax,1
            if( c(i).lt.afl(1,i,1) .or. c(i).gt.afl(1,i,2) ) return
        enddo
c       Make statistics for particle yield in partial space
        a(kk) = a(kk) + ww
        if( ABS(kf).ge.1 .AND. ABS(kf).le.3 ) a(2) = a(2) + ww

        return
        end



        subroutine output_hadron_distribution(sao,sbo,saof,sbof)
c300623 Lei
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(6),
     c   afl(20,6,2)
        dimension sao(40,6,20),sbo(20),saof(40,6,20),sbof(20)
c       The following statements are for user-output after event analyzing.
        dimension xcoor(40,6),xcoor_low(6)
        character*16 name_KF(20)
        character*5 cut_off(2,2)
        character*4 id_abscissa(6)
        data id_abscissa / "y", "pT", "eta", "mT", "mult", "pT" /
        character*25 id_distr(6)
        data id_distr / "dN/dy VS. y", "(1/pT)dN/dpT VS. pT", 
     &                  "dN/deta VS. eta", "(1/mT)dN/dmT VS. mT", 
     &                  "event-wise multiplicity", "dN/dpT VS. pT" /
        ! common/anly2/   ! For parton, in analy_parton.
        common/anly3/ xcoor, xcoor_low, cut_off, id_abscissa, id_distr
        dimension sum_h(40,6,2), i_h(6)

        if(iii.eq.nout .OR. iii.eq.neve)then
c           Calculates the abscissa of the distributions.
            do m2=1,isdmax
                xcoor_low(m2) = -0.5*asd(m2)
        if(m2.eq.1 .OR. m2.eq.3) xcoor_low(m2) = -21*asd(m2)+0.5*asd(m2)
                do m1=1,40
                    xcoor(m1,m2) = xcoor_low(m2) + m1*asd(m2)
                enddo
            enddo
c           Gets characters of the cut.
            write( cut_off(1,1),"(F5.2)" ) afl(1,2,1)
            write( cut_off(2,1),"(F5.2)" ) afl(1,2,2)
            write( cut_off(1,2),"(F5.2)" ) afl(1,1,1)
            write( cut_off(2,2),"(F5.2)" ) afl(1,1,2)
        endif
c       Gets names of the specified particles in usu.dat.
        do m3=1,ispmax,1
            call PYNAME(ispkf(m3),name_KF(m3))
        enddo
c       Finds the storing locations of pi, K and p.
        i_h = 0
        do i_kf=1,20,1
            if(ispkf(i_kf).eq. 211)  i_h(1) = i_kf
            if(ispkf(i_kf).eq.-211)  i_h(2) = i_kf
            if(ispkf(i_kf).eq. 321)  i_h(3) = i_kf
            if(ispkf(i_kf).eq.-321)  i_h(4) = i_kf
            if(ispkf(i_kf).eq. 2212) i_h(5) = i_kf
            if(ispkf(i_kf).eq.-2212) i_h(6) = i_kf
        end do
c       Calculates the multiplicities of pi, K and p, partial and full.
        sum_mul_h_partial = 0.
        sum_mul_h_full = 0.
        do ll=1,6,1
            sum_mul_h_partial = sum_mul_h_partial + sbo(i_h(ll))
            sum_mul_h_full    = sum_mul_h_full + sbof(i_h(ll))
        end do
c       Outputs particle multiplicities of 20 particles specified in usu.dat.
        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! particle multiplicity, partial ='
        write(10,*)"#! pi+K+p          ",
     &   ("          "//name_KF(m3),m3=1,20)
        write(10,*) sum_mul_h_partial, (sbo(ll),ll=1,ispmax)
        write(10,*)'#! particle multiplicity, full    ='
        write(10,*) sum_mul_h_full, (sbof(ll),ll=1,ispmax)
c       Calculates the distributions of pi, K and p, partial and full.
        sum_h = 0.
        do m2=1,isdmax,1
            do m1=1,40,1
                sum_h(m1,m2,1) = sao(m1,m2,i_h(1)) + sao(m1,m2,i_h(2))
     &                         + sao(m1,m2,i_h(3)) + sao(m1,m2,i_h(4))
     &                         + sao(m1,m2,i_h(5)) + sao(m1,m2,i_h(6))
                sum_h(m1,m2,2) = saof(m1,m2,i_h(1)) + saof(m1,m2,i_h(2))
     &                         + saof(m1,m2,i_h(3)) + saof(m1,m2,i_h(4))
     &                         + saof(m1,m2,i_h(5)) + saof(m1,m2,i_h(6))
            end do
        end do
c       Outputs abscissa, 6-distributions of pi+K+p and 20 particles specified in usu.dat.
        write(10,*)
        write(10,*)
        write(10,*)
        write(10,*) "#!*******************|"//
     &              "    Hadron  Distribution  Output    "//
     &              "|******************!#"

        do m2=1,isdmax,1

        write(10,*) "#!-------------------------------------"//
     &              "----------------------------------------"
        write(10,*) "#! ID of distribution m2=",m2,"   "//id_distr(m2)
        if(m2.eq.1 .OR. m2.eq.3)then
            write(10,*) "#! partial phase-space, "//
     &                  cut_off(1,1)//" < pT < "//cut_off(2,1)//
     &                  "  (nominal cuts of 1-st setting in usu.dat)"
        end if
        if(m2.eq.2 .OR. m2.eq.4 .OR. m2.eq.6)then
            write(10,*) "#! partial phase-space, "//
     &                  cut_off(1,2)//" < y/eta < "//cut_off(2,2)//
     &                  "  (nominal cuts of 1-st setting in usu.dat)"
        end if
        if(m2.eq.5)then
            write(10,*) "#! partial phase-space, "//
     &                  cut_off(1,1)//" < pT < "//cut_off(2,1)//";"//
     &                  cut_off(1,2)//" < y/eta < "//cut_off(2,2)//
     &                  "  (nominal cuts of 1-st setting in usu.dat)"
        end if
        write(10,*) "#!       "//id_abscissa(m2)//
     &              "                   pi+K+p        ",
     &              ("          "//name_KF(m3),m3=1,20,1)
        do m1=1,40,1
        write(10,*) xcoor(m1,m2),sum_h(m1,m2,1), (sao(m1,m2,m3),m3=1,20)
        enddo
        write(10,*) "#! ID of distribution m2=",m2,"   "//id_distr(m2)
        write(10,*) "#! full phase-space"
        write(10,*) "#!       "//id_abscissa(m2)//
     &              "                   pi+K+p        ",
     &              ("          "//name_KF(m3),m3=1,20,1)
        do m1=1,40,1
        write(10,*) xcoor(m1,m2),sum_h(m1,m2,2),(saof(m1,m2,m3),m3=1,20)
        enddo

        enddo


        return
        end



        subroutine output_parton_distribution
c300623 Lei
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(6),
     c   afl(20,6,2)
c       u, ubar, d, dbar, s, sbar, c, cbar, b, bbar, t, tbar, g, (u+d+s + anti-)
        common/anly2/ sn_min_p, sn_min_p_f, sn_cha_p, sn_cha_p_f, sn_g,
     &  sn_g_f, san_p(40,6,14), sbn_p(14), san_p_f(40,6,14), sbn_p_f(14)
        dimension xcoor(40,6),xcoor_low(6)
        character*16 name_KF(20)
        character*5 cut_off(2,2)
        character*4 id_abscissa(6)
        ! data id_abscissa / "y", "pT", "eta", "mT", "mult", "pT" /
        character*25 id_distr(6)
    !     data id_distr / "dN/dy VS. y", "(1/pT)dN/dpT VS. pT", 
    !  &                  "dN/deta VS. eta", "(1/mT)dN/dmT VS. mT", 
    !  &                  "event-wise multiplicity", "dN/dpT VS. pT" /
        common/anly3/ xcoor, xcoor_low, cut_off, id_abscissa, id_distr


c       Outputs abscissa, 6-distributions of pi+K+p and 20 particles specified in usu.dat.
        write(10,*)
        write(10,*)
        write(10,*)

        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)"#! multiplicity of negative, positive quark, "//
     &             "sums and gluon, partial & full ="
        write(10,*) sn_min_p/iii, sn_cha_p/iii, (sn_min_p+sn_cha_p)/iii, 
     &              sn_g/iii
        write(10,*) sn_min_p_f/iii, sn_cha_p_f/iii, 
     &              (sn_min_p_f+sn_cha_p_f)/iii, sn_g_f/iii

        write(10,*)"#!-------------------------------------"//
     &             "----------------------------------------"
        write(10,*)'#! particle multiplicity, partial ='
        write(10,*)   ! Hardcode
     &        "#! g                        u+d+s + anti-             "//
     &        "u                         ubar                      "//
     &        "d                         dbar                      "//
     &        "s                         sbar                      "//
     &        "c                         cbar                      "//
     &        "b                         bbar                      "//
     &        "t                         tbar"
        write(10,*) (sbn_p(ll)/iii,ll=1,14,1)
        write(10,*)'#! particle multiplicity, full    ='
        write(10,*) (sbn_p_f(ll)/iii,ll=1,14,1)

c       Outputs abscissa, 6-distributions of g, u+d+s + anti- and q.
        write(10,*)
        write(10,*)
        write(10,*)
        write(10,*) "#!*******************|"//
     &              "    Parton  Distribution  Output    "//
     &              "|******************!#"

        iii_org = iii   ! 220823 Lei
        do m2=1,isdmax,1

        write(10,*) "#!-------------------------------------"//
     &              "----------------------------------------"
        write(10,*)"#! ID of distribution m2=",m2,"   "//id_distr(m2)
        if(m2.eq.1 .OR. m2.eq.3)then
            write(10,*) "#! partial phase-space, "//
     &                  cut_off(1,1)//" < pT < "//cut_off(2,1)//
     &                  "  (nominal cuts of 1-st setting in usu.dat)"
        end if
        if(m2.eq.2 .OR. m2.eq.4 .OR. m2.eq.6)then
            write(10,*) "#! partial phase-space, "//
     &                  cut_off(1,2)//" < y/eta < "//cut_off(2,2)//
     &                  "  (nominal cuts of 1-st setting in usu.dat)"
        end if
        if(m2.eq.5)then
            write(10,*) "#! partial phase-space, "//
     &                  cut_off(1,1)//" < pT < "//cut_off(2,1)//";"//
     &                  cut_off(1,2)//" < y/eta < "//cut_off(2,2)//
     &                  "  (nominal cuts of 1-st setting in usu.dat)"
        end if
        write(10,*)"#!       "//id_abscissa(m2)//"               "//   ! Hardcode
     &        "g                         u+d+s + anti-             "//
     &        "u                         ubar                      "//
     &        "d                         dbar                      "//
     &        "s                         sbar                      "//
     &        "c                         cbar                      "//
     &        "b                         bbar                      "//
     &        "t                         tbar"
        do m1=1,40,1
            if(m2.eq.5) iii = 1   ! 220823 Lei
            write(10,*) xcoor(m1,m2), (san_p(m1,m2,m3)/iii,m3=1,14,1)
            if(m2.eq.5) iii = iii_org   ! 220823 Lei
        enddo
        write(10,*)"#! ID of distribution m2=",m2,"   "//id_distr(m2)
        write(10,*)"#! full phase-space"
        write(10,*)"#!       "//id_abscissa(m2)//"               "//   ! Hardcode
     &        "g                         u+d+s + anti-             "//
     &        "u                         ubar                      "//
     &        "d                         dbar                      "//
     &        "s                         sbar                      "//
     &        "c                         cbar                      "//
     &        "b                         bbar                      "//
     &        "t                         tbar"
        do m1=1,40,1
            write(10,*) xcoor(m1,m2),(san_p_f(m1,m2,m3),m3=1,14,1)
        enddo

        enddo


        return
        end
        


C***********************************************************************
 
C...PASTAT
C...Prints out information about cross-sections, decay widths, branching
C...ratios, kinematical limits, status codes and parameter values.
      SUBROUTINE PASTAT(I_STAT,I_CALL)
C300623 Lei
C...Modified from PYTHIA/PYSTAT for PACIAE. It just prints 
C    cross-sections now, i.e. hardcode of MSTAT=1.
C...The blocks used for cross section are /PYINT5/ and /PYINT5_S/
C    MSUB(I) stored switch of subprocesses.
C    NGEN(I,1) is the number of times the differential cross section has been 
C              evaluated for subprocess ISUB. NGEN(0,1) is the sum of these.
C    NGEN(I,3) is the number of times an event of subprocess 
C              type ISUB is generated. NGEN(0,3) is the sum of these.
C    XSEC(I,3) is the estimated integrated cross section for subproces ISUB.
C    XSEC(0,3) is the estimated the estimated total cross for section for 
C              all subprocesses included. (mb)
C    PROC(I) is the character strings for the different possible subprocesses.
C    PROC(0) is denotes all processes.

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
      PARAMETER (EPS=1D-3)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28, CHTMP*16
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYMSRV/RVLAM(3,3,3), RVLAMP(3,3,3), RVLAMB(3,3,3)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT2/,/PYINT4/,/PYINT5/,/PYINT6/,/PYMSSM/,/PYMSRV/
C...Local arrays, character variables and data.
      DIMENSION WDTP(0:400),WDTE(0:400,0:5),NMODES(0:20),PBRAT(10)
      CHARACTER PROGA(6)*28,CHAU*16,CHKF*16,CHD1*16,CHD2*16,CHD3*16,
     &CHIN(2)*12,STATE(-1:5)*4,CHKIN(21)*18,DISGA(2)*28,
     &PROGG9(13)*28,PROGG4(4)*28,PROGG2(2)*28,PROGP4(4)*28
      CHARACTER*24 CHD0, CHDC(10)
      CHARACTER*6 DNAME(3)
      DATA PROGA/
     &'VMD/hadron * VMD            ','VMD/hadron * direct         ',
     &'VMD/hadron * anomalous      ','direct * direct             ',
     &'direct * anomalous          ','anomalous * anomalous       '/
      DATA DISGA/'e * VMD','e * anomalous'/
      DATA PROGG9/
     &'direct * direct             ','direct * VMD                ',
     &'direct * anomalous          ','VMD * direct                ',
     &'VMD * VMD                   ','VMD * anomalous             ',
     &'anomalous * direct          ','anomalous * VMD             ',
     &'anomalous * anomalous       ','DIS * VMD                   ',
     &'DIS * anomalous             ','VMD * DIS                   ',
     &'anomalous * DIS             '/
      DATA PROGG4/
     &'direct * direct             ','direct * resolved           ',
     &'resolved * direct           ','resolved * resolved         '/
      DATA PROGG2/
     &'direct * hadron             ','resolved * hadron           '/
      DATA PROGP4/
     &'VMD * hadron                ','direct * hadron             ',
     &'anomalous * hadron          ','DIS * hadron                '/
      DATA STATE/'----','off ','on  ','on/+','on/-','on/1','on/2'/,
     &CHKIN/' m_hard (GeV/c^2) ',' p_T_hard (GeV/c) ',
     &'m_finite (GeV/c^2)','   y*_subsystem   ','     y*_large     ',
     &'     y*_small     ','    eta*_large    ','    eta*_small    ',
     &'cos(theta*)_large ','cos(theta*)_small ','       x_1        ',
     &'       x_2        ','       x_F        ',' cos(theta_hard)  ',
     &'m''_hard (GeV/c^2) ','       tau        ','        y*        ',
     &'cos(theta_hard^-) ','cos(theta_hard^+) ','      x_T^2       ',
     &'       tau''       '/
      DATA DNAME /'q     ','lepton','nu    '/
 
C300623 Lei
      COMMON/SA1/KJP21,NON1,BP,III,NEVE,NOUT,NOSC
      COMMON/SYSPAR/IPDEN,ITDEN,SUPPM,SUPTM,SUPPC,SUPTC,R0P,R0T,
     &NAP,NAT,NZP,NZT,PIO
      COMMON/PYINT5_1/NGENPD_1,NGEN_1(0:500,3),XSEC_1(0:500,3)
      COMMON/PYDAT1_1/MSTU_1(200),PARU_1(200),MSTJ_1(200),PARJ_1(200)
      COMMON/PYINT5_S/NGENPD_S,NGEN_S(0:500,3),XSEC_S(0:500,3)
      COMMON/PYDAT1_S/MSTU_S(200),PARU_S(200),MSTJ_S(200),PARJ_S(200)
C300623 Lei


C300623 Lei----------
      MSTAT = 1   ! 300623 Lei Hardcode!
C...Initializes global variables.
      IF(I_STAT.EQ.-2) THEN
        NGEN_S = 0
        XSEC_S = 0D0
        MSTU_S = 0
        PARU_S = 0D0
        MSTJ_S = 0
        PARJ_S = 0D0
C...Counts single-event cross sections.
      ELSEIF(I_STAT.EQ.-1) THEN
        IF(I_CALL.EQ.1) THEN
          NGEN_1 = 0
          XSEC_1 = 0D0
          MSTU_1 = 0
          PARU_1 = 0D0
          MSTJ_1 = 0
          PARJ_1 = 0D0
        ENDIF
        NGEN_1 = NGEN_1 + NGEN
        XSEC_1 = XSEC_1 + XSEC
C...Counts total-event cross sections and numbers of errors and warnings.
      ELSEIF(I_STAT.EQ.0) THEN
        NGEN_S = NGEN_S + NGEN_1
        XSEC_S = XSEC_S + XSEC_1
        MSTU_S(23)  = MSTU_S(23)  + MSTU(23)
        MSTU_S(27)  = MSTU_S(27)  + MSTU(27)
        MSTU_S(30)  = MSTU_S(30)  + MSTU(30)
      ELSEIF(I_STAT.EQ.1) THEN
        write(MSTU(11),*) "All of events finished, iii=", iii
        write(MSTU(11),"(/)")
!        write(MSTU(11),*) "****************** PACIAE-wise Event " //
!     &                    "Averaged Cross Sections *****************"
      ENDIF
      IF(MSTAT.NE.1 .OR. I_STAT.NE.1 .OR. III.NE.NEVE) RETURN
      NGEN = NGEN_S !/ III
      XSEC = XSEC_S !/ DBLE(III)
      MSTU(23)  = MSTU_S(23)
      MSTU(27)  = MSTU_S(27)
      MSTU(30)  = MSTU_S(30)
C300623 Lei----------

C...Cross-sections.
      IF(MSTAT.LE.1) THEN
        IF(MINT(121).GT.1) CALL PYSAVE(5,0)
        WRITE(MSTU(11),5000)
        WRITE(MSTU(11),5100)
        WRITE(MSTU(11),5200) 0,PROC(0),NGEN(0,3),NGEN(0,1),XSEC(0,3)
        DO 100 I=1,500
          IF(MSUB(I).NE.1) GOTO 100
          WRITE(MSTU(11),5200) I,PROC(I),NGEN(I,3),NGEN(I,1),XSEC(I,3)
  100   CONTINUE
        IF(MINT(121).GT.1) THEN
          WRITE(MSTU(11),5300)
          DO 110 IGA=1,MINT(121)
            CALL PYSAVE(3,IGA)
            IF(MINT(121).EQ.2.AND.MSTP(14).EQ.10) THEN
              WRITE(MSTU(11),5200) IGA,DISGA(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSEIF(MINT(121).EQ.9.OR.MINT(121).EQ.13) THEN
              WRITE(MSTU(11),5200) IGA,PROGG9(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSEIF(MINT(121).EQ.4.AND.MSTP(14).EQ.30) THEN
              WRITE(MSTU(11),5200) IGA,PROGP4(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSEIF(MINT(121).EQ.4) THEN
              WRITE(MSTU(11),5200) IGA,PROGG4(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSEIF(MINT(121).EQ.2) THEN
              WRITE(MSTU(11),5200) IGA,PROGG2(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSE
              WRITE(MSTU(11),5200) IGA,PROGA(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ENDIF
  110     CONTINUE
          CALL PYSAVE(5,0)
        ENDIF
        WRITE(MSTU(11),5400) MSTU(23),MSTU(30),MSTU(27),
     &  1D0-DBLE(NGEN(0,3))/MAX(1D0,DBLE(NGEN(0,2)))
 
C...Decay widths and branching ratios.
      ELSEIF(MSTAT.EQ.2) THEN
        WRITE(MSTU(11),5500)
        WRITE(MSTU(11),5600)
        DO 140 KC=1,500
          KF=KCHG(KC,4)
          CALL PYNAME(KF,CHKF)
          IOFF=0
          IF(KC.LE.22) THEN
            IF(KC.GT.2*MSTP(1).AND.KC.LE.10) GOTO 140
            IF(KC.GT.10+2*MSTP(1).AND.KC.LE.20) GOTO 140
            IF(KC.LE.5.OR.(KC.GE.11.AND.KC.LE.16)) IOFF=1
            IF(KC.EQ.18.AND.PMAS(18,1).LT.1D0) IOFF=1
            IF(KC.EQ.21.OR.KC.EQ.22) IOFF=1
          ELSE
            IF(MWID(KC).LE.0) GOTO 140
            IF(IMSS(1).LE.0.AND.(KF/KSUSY1.EQ.1.OR.
     &      KF/KSUSY1.EQ.2)) GOTO 140
          ENDIF
C...Off-shell branchings.
          IF(IOFF.EQ.1) THEN
            NGP=0
            IF(KC.LE.20) NGP=(MOD(KC,10)+1)/2
            IF(NGP.LE.MSTP(1)) WRITE(MSTU(11),5700) KF,CHKF(1:10),
     &      PMAS(KC,1),0D0,0D0,STATE(MDCY(KC,1)),0D0
            DO 120 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              NGP1=0
              IF(IABS(KFDP(IDC,1)).LE.20) NGP1=
     &        (MOD(IABS(KFDP(IDC,1)),10)+1)/2
              NGP2=0
              IF(IABS(KFDP(IDC,2)).LE.20) NGP2=
     &        (MOD(IABS(KFDP(IDC,2)),10)+1)/2
              CALL PYNAME(KFDP(IDC,1),CHD1)
              CALL PYNAME(KFDP(IDC,2),CHD2)
              IF(KFDP(IDC,3).EQ.0) THEN
                IF(MDME(IDC,2).EQ.102.AND.NGP1.LE.MSTP(1).AND.
     &          NGP2.LE.MSTP(1)) WRITE(MSTU(11),5800) IDC,CHD1(1:10),
     &          CHD2(1:10),0D0,0D0,STATE(MDME(IDC,1)),0D0
              ELSE
                CALL PYNAME(KFDP(IDC,3),CHD3)
                IF(MDME(IDC,2).EQ.102.AND.NGP1.LE.MSTP(1).AND.
     &          NGP2.LE.MSTP(1)) WRITE(MSTU(11),5900) IDC,CHD1(1:10),
     &          CHD2(1:10),CHD3(1:10),0D0,0D0,STATE(MDME(IDC,1)),0D0
              ENDIF
  120       CONTINUE
C...On-shell decays.
          ELSE
            CALL PYWIDT(KF,PMAS(KC,1)**2,WDTP,WDTE)
            BRFIN=1D0
            IF(WDTE(0,0).LE.0D0) BRFIN=0D0
            WRITE(MSTU(11),5700) KF,CHKF(1:10),PMAS(KC,1),WDTP(0),1D0,
     &      STATE(MDCY(KC,1)),BRFIN
            DO 130 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              NGP1=0
              IF(IABS(KFDP(IDC,1)).LE.20) NGP1=
     &        (MOD(IABS(KFDP(IDC,1)),10)+1)/2
              NGP2=0
              IF(IABS(KFDP(IDC,2)).LE.20) NGP2=
     &        (MOD(IABS(KFDP(IDC,2)),10)+1)/2
              BRPRI=0D0
              IF(WDTP(0).GT.0D0) BRPRI=WDTP(J)/WDTP(0)
              BRFIN=0D0
              IF(WDTE(0,0).GT.0D0) BRFIN=WDTE(J,0)/WDTE(0,0)
              CALL PYNAME(KFDP(IDC,1),CHD1)
              CALL PYNAME(KFDP(IDC,2),CHD2)
              IF(KFDP(IDC,3).EQ.0) THEN
                IF(NGP1.LE.MSTP(1).AND.NGP2.LE.MSTP(1))
     &          WRITE(MSTU(11),5800) IDC,CHD1(1:10),
     &          CHD2(1:10),WDTP(J),BRPRI,
     &          STATE(MDME(IDC,1)),BRFIN
              ELSE
                CALL PYNAME(KFDP(IDC,3),CHD3)
                IF(NGP1.LE.MSTP(1).AND.NGP2.LE.MSTP(1))
     &          WRITE(MSTU(11),5900) IDC,CHD1(1:10),
     &          CHD2(1:10),CHD3(1:10),WDTP(J),BRPRI,
     &          STATE(MDME(IDC,1)),BRFIN
              ENDIF
  130       CONTINUE
          ENDIF
  140   CONTINUE
        WRITE(MSTU(11),6000)
 
C...Allowed incoming partons/particles at hard interaction.
      ELSEIF(MSTAT.EQ.3) THEN
        WRITE(MSTU(11),6100)
        CALL PYNAME(MINT(11),CHAU)
        CHIN(1)=CHAU(1:12)
        CALL PYNAME(MINT(12),CHAU)
        CHIN(2)=CHAU(1:12)
        WRITE(MSTU(11),6200) CHIN(1),CHIN(2)
        DO 150 I=-20,22
          IF(I.EQ.0) GOTO 150
          IA=IABS(I)
          IF(IA.GT.MSTP(58).AND.IA.LE.10) GOTO 150
          IF(IA.GT.10+2*MSTP(1).AND.IA.LE.20) GOTO 150
          CALL PYNAME(I,CHAU)
          WRITE(MSTU(11),6300) CHAU,STATE(KFIN(1,I)),CHAU,
     &    STATE(KFIN(2,I))
  150   CONTINUE
        WRITE(MSTU(11),6400)
 
C...User-defined limits on kinematical variables.
      ELSEIF(MSTAT.EQ.4) THEN
        WRITE(MSTU(11),6500)
        WRITE(MSTU(11),6600)
        SHRMAX=CKIN(2)
        IF(SHRMAX.LT.0D0) SHRMAX=VINT(1)
        WRITE(MSTU(11),6700) CKIN(1),CHKIN(1),SHRMAX
        PTHMIN=MAX(CKIN(3),CKIN(5))
        PTHMAX=CKIN(4)
        IF(PTHMAX.LT.0D0) PTHMAX=0.5D0*SHRMAX
        WRITE(MSTU(11),6800) CKIN(3),PTHMIN,CHKIN(2),PTHMAX
        WRITE(MSTU(11),6900) CHKIN(3),CKIN(6)
        DO 160 I=4,14
          WRITE(MSTU(11),6700) CKIN(2*I-1),CHKIN(I),CKIN(2*I)
  160   CONTINUE
        SPRMAX=CKIN(32)
        IF(SPRMAX.LT.0D0) SPRMAX=VINT(1)
        WRITE(MSTU(11),6700) CKIN(31),CHKIN(15),SPRMAX
        WRITE(MSTU(11),7000)
 
C...Status codes and parameter values.
      ELSEIF(MSTAT.EQ.5) THEN
        WRITE(MSTU(11),7100)
        WRITE(MSTU(11),7200)
        DO 170 I=1,100
          WRITE(MSTU(11),7300) I,MSTP(I),PARP(I),100+I,MSTP(100+I),
     &    PARP(100+I)
  170   CONTINUE
 
C...List of all processes implemented in the program.
      ELSEIF(MSTAT.EQ.6) THEN
        WRITE(MSTU(11),7400)
        WRITE(MSTU(11),7500)
        DO 180 I=1,500
          IF(ISET(I).LT.0) GOTO 180
          WRITE(MSTU(11),7600) I,PROC(I),ISET(I),KFPR(I,1),KFPR(I,2)
  180   CONTINUE
        WRITE(MSTU(11),7700)
 
      ELSEIF(MSTAT.EQ.7) THEN
      WRITE (MSTU(11),8000)
      NMODES(0)=0
      NMODES(10)=0
      NMODES(9)=0
      DO 290 ILR=1,2
        DO 280 KFSM=1,16
          KFSUSY=ILR*KSUSY1+KFSM
          NRVDC=0
C...SDOWN DECAYS
          IF (KFSM.EQ.1.OR.KFSM.EQ.3.OR.KFSM.EQ.5) THEN
            NRVDC=3
            DO 190 I=1,NRVDC
              PBRAT(I)=0D0
              NMODES(I)=0
  190       CONTINUE
            CALL PYNAME(KFSUSY,CHTMP)
            CHD0=CHTMP//' '
            CHDC(1)=DNAME(3) // ' + ' // DNAME(1)
            CHDC(2)=DNAME(2) // ' + ' // DNAME(1)
            CHDC(3)=DNAME(1) // ' + ' // DNAME(1)
            KC=PYCOMP(KFSUSY)
            DO 200 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              ID1=IABS(KFDP(IDC,1))
              ID2=IABS(KFDP(IDC,2))
              IF (KFDP(IDC,3).EQ.0) THEN
                IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(1)=PBRAT(1)+BRAT(IDC)
                  NMODES(1)=NMODES(1)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &                 .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6)) THEN
                  PBRAT(2)=PBRAT(2)+BRAT(IDC)
                  NMODES(2)=NMODES(2)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ELSE IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND
     &                 .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(3)=PBRAT(3)+BRAT(IDC)
                  NMODES(3)=NMODES(3)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
              ENDIF
  200       CONTINUE
          ENDIF
C...SUP DECAYS
          IF (KFSM.EQ.2.OR.KFSM.EQ.4.OR.KFSM.EQ.6) THEN
            NRVDC=2
            DO 210 I=1,NRVDC
              NMODES(I)=0
              PBRAT(I)=0D0
  210       CONTINUE
            CALL PYNAME(KFSUSY,CHTMP)
            CHD0=CHTMP//' '
            CHDC(1)=DNAME(2) // ' + ' // DNAME(1)
            CHDC(2)=DNAME(1) // ' + ' // DNAME(1)
            KC=PYCOMP(KFSUSY)
            DO 220 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              ID1=IABS(KFDP(IDC,1))
              ID2=IABS(KFDP(IDC,2))
              IF (KFDP(IDC,3).EQ.0) THEN
                IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(1)=PBRAT(1)+BRAT(IDC)
                  NMODES(1)=NMODES(1)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ELSE IF ((ID1.EQ.1.OR.ID1.EQ.3.OR.ID1.EQ.5).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(2)=PBRAT(2)+BRAT(IDC)
                  NMODES(2)=NMODES(2)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
              ENDIF
  220       CONTINUE
          ENDIF
C...SLEPTON DECAYS
          IF (KFSM.EQ.11.OR.KFSM.EQ.13.OR.KFSM.EQ.15) THEN
            NRVDC=2
            DO 230 I=1,NRVDC
              PBRAT(I)=0D0
              NMODES(I)=0
  230       CONTINUE
            CALL PYNAME(KFSUSY,CHTMP)
            CHD0=CHTMP//' '
            CHDC(1)=DNAME(3) // ' + ' // DNAME(2)
            CHDC(2)=DNAME(1) // ' + ' // DNAME(1)
            KC=PYCOMP(KFSUSY)
            DO 240 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              ID1=IABS(KFDP(IDC,1))
              ID2=IABS(KFDP(IDC,2))
              IF (KFDP(IDC,3).EQ.0) THEN
                IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &               .EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15)) THEN
                  PBRAT(1)=PBRAT(1)+BRAT(IDC)
                  NMODES(1)=NMODES(1)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
                IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(2)=PBRAT(2)+BRAT(IDC)
                  NMODES(2)=NMODES(2)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
              ENDIF
  240       CONTINUE
          ENDIF
C...SNEUTRINO DECAYS
          IF ((KFSM.EQ.12.OR.KFSM.EQ.14.OR.KFSM.EQ.16).AND.ILR.EQ.1)
     &         THEN
            NRVDC=2
            DO 250 I=1,NRVDC
              PBRAT(I)=0D0
              NMODES(I)=0
  250       CONTINUE
            CALL PYNAME(KFSUSY,CHTMP)
            CHD0=CHTMP//' '
            CHDC(1)=DNAME(2) // ' + ' // DNAME(2)
            CHDC(2)=DNAME(1) // ' + ' // DNAME(1)
            KC=PYCOMP(KFSUSY)
            DO 260 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              ID1=IABS(KFDP(IDC,1))
              ID2=IABS(KFDP(IDC,2))
              IF (KFDP(IDC,3).EQ.0) THEN
                IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND.(ID2
     &               .EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15)) THEN
                  PBRAT(1)=PBRAT(1)+BRAT(IDC)
                  NMODES(1)=NMODES(1)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
                IF ((ID1.EQ.1.OR.ID1.EQ.3.OR.ID1.EQ.5).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  NMODES(2)=NMODES(2)+1
                  PBRAT(2)=PBRAT(2)+BRAT(IDC)
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
              ENDIF
  260       CONTINUE
          ENDIF
          IF (NRVDC.NE.0) THEN
            DO 270 I=1,NRVDC
              WRITE (MSTU(11),8200) CHD0, CHDC(I), PBRAT(I), NMODES(I)
              NMODES(0)=NMODES(0)+NMODES(I)
  270       CONTINUE
          ENDIF
  280   CONTINUE
  290 CONTINUE
      DO 370 KFSM=21,37
        KFSUSY=KSUSY1+KFSM
        NRVDC=0
C...NEUTRALINO DECAYS
        IF (KFSM.EQ.22.OR.KFSM.EQ.23.OR.KFSM.EQ.25.OR.KFSM.EQ.35) THEN
          NRVDC=4
          DO 300 I=1,NRVDC
            PBRAT(I)=0D0
            NMODES(I)=0
  300     CONTINUE
          CALL PYNAME(KFSUSY,CHTMP)
          CHD0=CHTMP//' '
          CHDC(1)=DNAME(3) // ' + ' // DNAME(2) // ' + ' // DNAME(2)
          CHDC(2)=DNAME(3) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(3)=DNAME(2) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(4)=DNAME(1) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          KC=PYCOMP(KFSUSY)
          DO 310 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            ID1=IABS(KFDP(IDC,1))
            ID2=IABS(KFDP(IDC,2))
            ID3=IABS(KFDP(IDC,3))
            IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &           .EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15).AND.(ID3.EQ.11.OR
     &           .ID3.EQ.13.OR.ID3.EQ.15)) THEN
              PBRAT(1)=PBRAT(1)+BRAT(IDC)
              NMODES(1)=NMODES(1)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(2)=PBRAT(2)+BRAT(IDC)
              NMODES(2)=NMODES(2)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(3)=PBRAT(3)+BRAT(IDC)
              NMODES(3)=NMODES(3)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(4)=PBRAT(4)+BRAT(IDC)
              NMODES(4)=NMODES(4)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ENDIF
  310     CONTINUE
        ENDIF
C...CHARGINO DECAYS
        IF (KFSM.EQ.24.OR.KFSM.EQ.37) THEN
          NRVDC=5
          DO 320 I=1,NRVDC
            PBRAT(I)=0D0
            NMODES(I)=0
  320     CONTINUE
          CALL PYNAME(KFSUSY,CHTMP)
          CHD0=CHTMP//' '
          CHDC(1)=DNAME(3) // ' + ' // DNAME(3) // ' + ' // DNAME(2)
          CHDC(2)=DNAME(2) // ' + ' // DNAME(2) // ' + ' // DNAME(2)
          CHDC(3)=DNAME(3) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(4)=DNAME(2) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(5)=DNAME(1) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          KC=PYCOMP(KFSUSY)
          DO 330 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            ID1=IABS(KFDP(IDC,1))
            ID2=IABS(KFDP(IDC,2))
            ID3=IABS(KFDP(IDC,3))
            IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &           .EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15).AND.(ID3.EQ.12.OR
     &           .ID3.EQ.14.OR.ID3.EQ.16)) THEN
              PBRAT(1)=PBRAT(1)+BRAT(IDC)
              NMODES(1)=NMODES(1)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND
     &             .(ID2.EQ.12.OR.ID2.EQ.14.OR.ID2.EQ.16).AND.(ID3.EQ
     &             .11.OR.ID3.EQ.13.OR.ID3.EQ.15)) THEN
              PBRAT(1)=PBRAT(1)+BRAT(IDC)
              NMODES(1)=NMODES(1)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15).AND.(ID3.EQ
     &             .11.OR.ID3.EQ.13.OR.ID3.EQ.15)) THEN
              PBRAT(2)=PBRAT(2)+BRAT(IDC)
              NMODES(2)=NMODES(2)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ
     &             .2.OR.ID3.EQ.4.OR.ID3.EQ.6)) THEN
              PBRAT(3)=PBRAT(3)+BRAT(IDC)
              NMODES(3)=NMODES(3)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ
     &             .1.OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(3)=PBRAT(3)+BRAT(IDC)
              NMODES(3)=NMODES(3)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ
     &             .2.OR.ID3.EQ.4.OR.ID3.EQ.6)) THEN
              PBRAT(4)=PBRAT(4)+BRAT(IDC)
              NMODES(4)=NMODES(4)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ
     &             .1.OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(4)=PBRAT(4)+BRAT(IDC)
              NMODES(4)=NMODES(4)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ
     &             .1.OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(5)=PBRAT(5)+BRAT(IDC)
              NMODES(5)=NMODES(5)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.1.OR.ID1.EQ.3.OR.ID1.EQ.5).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ
     &             .1.OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(5)=PBRAT(5)+BRAT(IDC)
              NMODES(5)=NMODES(5)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ENDIF
  330     CONTINUE
        ENDIF
C...GLUINO DECAYS
        IF (KFSM.EQ.21) THEN
          NRVDC=3
          DO 340 I=1,NRVDC
            PBRAT(I)=0D0
            NMODES(I)=0
  340     CONTINUE
          CALL PYNAME(KFSUSY,CHTMP)
          CHD0=CHTMP//' '
          CHDC(1)=DNAME(3) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(2)=DNAME(2) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(3)=DNAME(1) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          KC=PYCOMP(KFSUSY)
          DO 350 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            ID1=IABS(KFDP(IDC,1))
            ID2=IABS(KFDP(IDC,2))
            ID3=IABS(KFDP(IDC,3))
            IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &           .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ.1.OR
     &           .ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(1)=PBRAT(1)+BRAT(IDC)
              NMODES(1)=NMODES(1)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(2)=PBRAT(2)+BRAT(IDC)
              NMODES(2)=NMODES(2)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(3)=PBRAT(3)+BRAT(IDC)
              NMODES(3)=NMODES(3)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ENDIF
  350     CONTINUE
        ENDIF
 
        IF (NRVDC.NE.0) THEN
          DO 360 I=1,NRVDC
            WRITE (MSTU(11),8200) CHD0, CHDC(I), PBRAT(I), NMODES(I)
            NMODES(0)=NMODES(0)+NMODES(I)
  360     CONTINUE
        ENDIF
  370 CONTINUE
      WRITE (MSTU(11),8100) NMODES(0), NMODES(10), NMODES(9)
 
      IF (IMSS(51).GE.1.OR.IMSS(52).GE.1.OR.IMSS(53).GE.1) THEN
        WRITE (MSTU(11),8500)
        DO 400 IRV=1,3
          DO 390 JRV=1,3
            DO 380 KRV=1,3
              WRITE (MSTU(11),8700) IRV,JRV,KRV,RVLAM(IRV,JRV,KRV)
     &             ,RVLAMP(IRV,JRV,KRV),RVLAMB(IRV,JRV,KRV)
  380       CONTINUE
  390     CONTINUE
  400   CONTINUE
        WRITE (MSTU(11),8600)
      ENDIF
      ENDIF
 
C...Formats for printouts.
!  5000 FORMAT('1',9('*'),1X,'PYSTAT:  Statistics on Number of ',
!     &'Events and Cross-sections',1X,9('*'))
 5000 FORMAT('1',5('*'),1X,'PASTAT: Statistics on Number of ',
     &'NN Events and Total Cross-sections',1X,6('*'))
 5100 FORMAT(/1X,78('=')/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',12X,
     &'Subprocess',12X,'I',6X,'Number of points',6X,'I',4X,'Sigma',3X,
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',34('-'),'I',28('-'),
     &'I',4X,'(mb)',4X,'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',1X,
     &'N:o',1X,'Type',25X,'I',4X,'Generated',9X,'Tried',1X,'I',12X,
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/1X,'I',34X,'I',28X,
     &'I',12X,'I')
 5200 FORMAT(1X,'I',1X,I3,1X,A28,1X,'I',1X,I12,1X,I13,1X,'I',1X,1P,
     &D10.3,1X,'I')
 5300 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/
     &1X,'I',34X,'I',28X,'I',12X,'I')
 5400 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')//
     &1X,'********* Total number of errors, excluding junctions =',
     &1X,I8,' *************'/
     &1X,'********* Total number of errors, including junctions =',
     &1X,I8,' *************'/
     &1X,'********* Total number of warnings =                   ',
     &1X,I8,' *************'/
     &1X,'********* Fraction of events that fail fragmentation ',
     &'cuts =',1X,F8.5,' *********'/)
!  5500 FORMAT('1',27('*'),1X,'PYSTAT:  Decay Widths and Branching ',
 5500 FORMAT('1',27('*'),1X,'PASTAT:  Decay Widths and Branching ',
     &'Ratios',1X,27('*'))
 5600 FORMAT(/1X,98('=')/1X,'I',49X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/
     &1X,'I',5X,'Mother  -->  Branching/Decay Channel',8X,'I',1X,
     &'Width (GeV)',1X,'I',7X,'B.R.',1X,'I',1X,'Stat',1X,'I',2X,
     &'Eff. B.R.',1X,'I'/1X,'I',49X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/
     &1X,98('='))
 5700 FORMAT(1X,'I',49X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,'I',1X,
     &I8,2X,A10,3X,'(m =',F10.3,')',2X,'-->',5X,'I',2X,1P,D10.3,0P,1X,
     &'I',1X,1P,D10.3,0P,1X,'I',1X,A4,1X,'I',1X,1P,D10.3,0P,1X,'I')
 5800 FORMAT(1X,'I',1X,I8,2X,A10,1X,'+',1X,A10,15X,'I',2X,
     &1P,D10.3,0P,1X,'I',1X,1P,D10.3,0P,1X,'I',1X,A4,1X,'I',1X,
     &1P,D10.3,0P,1X,'I')
 5900 FORMAT(1X,'I',1X,I8,2X,A10,1X,'+',1X,A10,1X,'+',1X,A10,2X,'I',2X,
     &1P,D10.3,0P,1X,'I',1X,1P,D10.3,0P,1X,'I',1X,A4,1X,'I',1X,
     &1P,D10.3,0P,1X,'I')
 6000 FORMAT(1X,'I',49X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,98('='))
!  6100 FORMAT('1',7('*'),1X,'PYSTAT: Allowed Incoming Partons/',
 6100 FORMAT('1',7('*'),1X,'PASTAT: Allowed Incoming Partons/',
     &'Particles at Hard Interaction',1X,7('*'))
 6200 FORMAT(/1X,78('=')/1X,'I',38X,'I',37X,'I'/1X,'I',1X,
     &'Beam particle:',1X,A12,10X,'I',1X,'Target particle:',1X,A12,7X,
     &'I'/1X,'I',38X,'I',37X,'I'/1X,'I',1X,'Content',6X,'State',19X,
     &'I',1X,'Content',6X,'State',18X,'I'/1X,'I',38X,'I',37X,'I'/1X,
     &78('=')/1X,'I',38X,'I',37X,'I')
 6300 FORMAT(1X,'I',1X,A9,5X,A4,19X,'I',1X,A9,5X,A4,18X,'I')
 6400 FORMAT(1X,'I',38X,'I',37X,'I'/1X,78('='))
!  6500 FORMAT('1',12('*'),1X,'PYSTAT: User-Defined Limits on ',
 6500 FORMAT('1',12('*'),1X,'PASTAT: User-Defined Limits on ',
     &'Kinematical Variables',1X,12('*'))
 6600 FORMAT(/1X,78('=')/1X,'I',76X,'I')
 6700 FORMAT(1X,'I',16X,1P,D10.3,0P,1X,'<',1X,A,1X,'<',1X,1P,D10.3,0P,
     &16X,'I')
 6800 FORMAT(1X,'I',3X,1P,D10.3,0P,1X,'(',1P,D10.3,0P,')',1X,'<',1X,A,
     &1X,'<',1X,1P,D10.3,0P,16X,'I')
 6900 FORMAT(1X,'I',29X,A,1X,'=',1X,1P,D10.3,0P,16X,'I')
 7000 FORMAT(1X,'I',76X,'I'/1X,78('='))
!  7100 FORMAT('1',12('*'),1X,'PYSTAT: Summary of Status Codes and ',
 7100 FORMAT('1',12('*'),1X,'PASTAT: Summary of Status Codes and ',
     &'Parameter Values',1X,12('*'))
 7200 FORMAT(/3X,'I',4X,'MSTP(I)',9X,'PARP(I)',20X,'I',4X,'MSTP(I)',9X,
     &'PARP(I)'/)
 7300 FORMAT(1X,I3,5X,I6,6X,1P,D10.3,0P,18X,I3,5X,I6,6X,1P,D10.3)
!  7400 FORMAT('1',13('*'),1X,'PYSTAT: List of implemented processes',
 7400 FORMAT('1',13('*'),1X,'PASTAT: List of implemented processes',
     &1X,13('*'))
 7500 FORMAT(/1X,65('=')/1X,'I',34X,'I',28X,'I'/1X,'I',12X,
     &'Subprocess',12X,'I',1X,'ISET',2X,'KFPR(I,1)',2X,'KFPR(I,2)',1X,
     &'I'/1X,'I',34X,'I',28X,'I'/1X,65('=')/1X,'I',34X,'I',28X,'I')
 7600 FORMAT(1X,'I',1X,I3,1X,A28,1X,'I',1X,I4,1X,I10,1X,I10,1X,'I')
 7700 FORMAT(1X,'I',34X,'I',28X,'I'/1X,65('='))
 8000 FORMAT(1X/ 1X/
     &     17X,'Sums over R-Violating branching ratios',1X/ 1X
     &     /1X,70('=')/1X,'I',50X,'I',11X,'I',5X,'I'/1X,'I',4X
     &     ,'Mother  -->  Sum over final state flavours',4X,'I',2X
     &     ,'BR(sum)',2X,'I',2X,'N',2X,'I'/1X,'I',50X,'I',11X,'I',5X,'I'
     &     /1X,70('=')/1X,'I',50X,'I',11X,'I',5X,'I')
 8100 FORMAT(1X,'I',50X,'I',11X,'I',5X,'I'/1X,70('=')/1X,'I',1X
     &     ,'Total number of R-Violating modes :',3X,I5,24X,'I'/
     &     1X,'I',1X,'Total number with non-vanishing BR :',2X,I5,24X
     &     ,'I'/1X,'I',1X,'Total number with BR > 0.001 :',8X,I5,24X,'I'
     &     /1X,70('='))
 8200 FORMAT(1X,'I',1X,A9,1X,'-->',1X,A24,11X,
     &     'I',2X,1P,D8.2,0P,1X,'I',2X,I2,1X,'I')
 8300 FORMAT(1X,'I',50X,'I',11X,'I',5X,'I')
 8500 FORMAT(1X/ 1X/
     &     1X,'R-Violating couplings',1X/ 1X /
     &     1X,55('=')/
     &     1X,'I',1X,'IJK',1X,'I',2X,'LAMBDA(IJK)',2X,'I',2X
     &     ,'LAMBDA''(IJK)',1X,'I',1X,"LAMBDA''(IJK)",1X,'I'/1X,'I',5X
     &     ,'I',15X,'I',15X,'I',15X,'I')
 8600 FORMAT(1X,55('='))
 8700 FORMAT(1X,'I',1X,I1,I1,I1,1X,'I',1X,1P,D13.3,0P,1X,'I',1X,1P
     &     ,D13.3,0P,1X,'I',1X,1P,D13.3,0P,1X,'I')
 
      RETURN
      END
 