        subroutine hadcas(ijk,neve,nout,time_had,ijkk)
c       deals with the hadronic rescattering, 
c        composed by Ben-Hao Sa, 20/09/2000
c       input message is in 'sa1_h', which plays working block as well
c       output message is in 'sa1_h'
c       ijk: the event number
c       neve: total number of events
c       nout: a internal output per nout events
c060112 if ijkk=1 give up current event avoiding infinite loop
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/sa19_h/coor(3)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,ecsspn,ecsspm
        common/count_h/isinel(600)
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/syspar_h/pio
        common/sa24/adj1(40),nnstop,non24,zstop   ! 231104
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2   ! 221203 250204
        common/pycidat2/kfmaxt,nont2,PARAM(20),weigh(600)   ! 250204
c       ifram = 0 for fixed target, = 1 for collider 
c       cspipi (fm^2): total cross section of pion + pion
c       cspin (fm^2): total cross section of pion + nucleon 
c       cskn (fm^2): total cross section of kaon + nucleon 
c       csnn (fm^2): total cross section of n + n 
c       cspsn: total cross section of J/Psi (Psi') + n
c       cspsm: total cross section of J/Psi (Psi') + meson
c       rcsit: ratio of inelastic to total cross section
c       kfmax: the maximum # of particles with given flavor code 
c       kfaco(i): flavor code of i-th particle among kfmax
c       numb(i): order # of last particle of particles with same flavor of 
c        kfaco(i) in particle list
c       disbe(i,j): allowable minimum approaching distance between particles
c                   kfaco(i) & kfaco(j)
c       sig (fm^2): cross section of pion + pion to kaon + kaon
c       edipi: largest interaction distance between two pions.
c       epin: largest interaction distance between pion and nucleon.
c       ekn: largest interaction distance between kaon and nucleon.
c       ecsnn: largest interaction distance between two nucleons.
c       t0: average proper formation time at rest.
c       ddt: time accuracy 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       0. is hard core distance between two pions
c       0.5 is hard core distance between two nucleons
c       0. is hard core distance between pion and nucleon
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       tau(i) : formation time of particle i.
c       ishp(i)=1 if i-th particle inside the simulated volume
c              =0 if i-th particle outside the simulated volume
c       isinel(i) = 0 without i-th inelastic process
c                 = 1 with i-th inelastic process
c       lc(i,1) and lc(i,2) are, respectively, line # in particle
c        list of the colliding particles of i-th collision pair
c       lc(i,3) and lc(i,4) are the flavor codes of scattered particles
c        of i-th collision pair
c       lc(i,5) identifies the different inelastic process,
c        lc(i,5)=592, not used
c       tc(i): collision time of i-th colliding pair
c       tw(i): cross section ratio of (i-th inelas.)/tot
c       pio : 3.1416
c       nctl: number of collision pairs in the current collision list
c       nctl0: number of collision pairs in last collision list
c       noinel(i): statistics of the occurring of i-th inelastic process
c       noel: statistics of the occurring of elastic process
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension peo(4)
        time=time_had   ! recovered on 280910
        PARAM(1)=para1_2   ! 250204
        pio=3.1416
c110923 iabsb=0   ! 110923 Lei
c110923 iabsm=0   ! 110923 Lei
c       iabsb = 0 : without J/Psi (Psi') + baryon
c             = 1 : with J/Psi (Psi') + baryon
c       iabsm = 0 : without J/Psi (Psi') + meson
c             = 1 : with J/Psi (Psi') + meson

c280910 if(ijk.eq.1)then
c       give initial value to quantities needed in hardon rescattering
        call sysini_h   ! it has been called in main.f
c280910 endif

c       initiation
        nctl=0
        lc=0
        tc=0.
        tw=0.
        numb=0
        do i=1,nsa
        vsa(i,4)=0.   ! 231104
        tau(i)=0.
        enddo
c231104
        dpmax=adj1(27)
        drmax=adj1(28)
        do i1=1,nsa
        pnn1=psa(i1,1)
        pnn2=psa(i1,2)
        pnn3=psa(i1,3)
        pnn4=psa(i1,4)
        rnn1=vsa(i1,1)
        rnn2=vsa(i1,2)
        rnn3=vsa(i1,3)
        pnnm=pnn1*pnn1+pnn2*pnn2+pnn3*pnn3
        if(pnnm.lt.1.e-28)pnnm=1.e-28
        if(pnnm.gt.1.e28)then
        ishp(i1)=0
        goto 200
        endif
        pnnm=sqrt(pnnm)
        rnnm=rnn1*rnn1+rnn2*rnn2+rnn3*rnn3
        if(rnnm.lt.1.e-28)rnnm=1.e-28
        if(rnnm.gt.1.e28)then
        ishp(i1)=0
        goto 200
        endif
        rnnm=sqrt(rnnm)
        if((pnnm.le.dpmax.and.pnn4.le.dpmax).and.rnnm.le.drmax)then
        ishp(i1)=1
        else
        ishp(i1)=0
        endif
200     enddo
c231104
        noel=0
        noinel=0

c       change K0S, K0L to K0, K0ba
        do j=1,nsa
        kf=ksa(j,2)
        if(kf.eq.130 .or. kf.eq.310)then
        rrlu=pyr(1)
        ksa(j,2)=311
        if(rrlu.gt.0.5)ksa(j,2)=-311
        endif
        enddo

c       filter out particles wanted to study and make in order of proton, 
c        neutron, ...
c       initial particle list is compsed of the arraies in common block 
c        'sa1_h', 'tau' and 'ishp' in 'sa8_h', and 'numb' in 'sa9_h'
c       call prt_sa1_h(nsa) ! sa 
        call filt_h
c       call prt_sa1_h(nsa) ! sa
c280910 time=0.

c       calculate position of center of mass of the system. distance of a 
c        particle from this cms is used to check whether it is freezes out 
c        or not
        call copl_h(time)

c       creat the initial collision list, note: be sure that the initial
c        collision list must not be empty
        call ctlcre_h(lc,tc,tw,time)

c       administrate hadron rescattering
        call scat_h(time,lc,tc,tw,ijkk,ijk,neve)
        if(ijkk.eq.1)return   ! 100603
c       if ijkk=1 (infinite loops) give up the event 
c       call prt_sa1_h(nsa) ! sa
        time_had=time

c       change K0,K0ba to K0L and K0S
        do j=1,nsa
        kf=ksa(j,2)
        if(kf.eq.311 .or. kf.eq.-311)then
        rrlu=pyr(1)
        ksa(j,2)=130
        if(rrlu.gt.0.5)ksa(j,2)=310
        endif
        enddo
c       call prt_sa1_h(nsa) ! sa
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sysini_h
c       give the initial values to the quantities needed in hardon rescattering
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,ecsspn,ecsspm
        common/sa24/adj1(40),nnstop,non24,zstop   ! 060404
        common/syspar_h/pio
        common/count_h/isinel(600)
c       cross sections are given in fm.
        csnn=PARAM(1)*0.1
        cspin=PARAM(2)*0.1
        cskn=PARAM(3)*0.1
        cspipi=PARAM(4)*0.1
        cspsn=PARAM(13)*0.1
        cspsm=PARAM(14)*0.1
        csspn=PARAM(15)*0.1
        csspm=PARAM(16)*0.1
c       largest interaction distance of two colliding particles.
        edipi=sqrt(cspipi/3.1416)
        epin=sqrt(cspin/3.1416)
        ekn=sqrt(cskn/3.1416)
        ecsnn=sqrt(csnn/3.1416)
        ecspsn=sqrt(cspsn/3.1416)
        ecspsm=sqrt(cspsm/3.1416)
        ecsspn=sqrt(csspn/3.1416)
        ecsspm=sqrt(csspm/3.1416)

        sig=PARAM(5)*0.1
        rcsit=PARAM(6)
        t0=PARAM(7)
c       dep=PARAM(9)
c060404 ddt=PARAM(8)
        ddt=adj1(11)   ! 060404
c       rao=PARAM(10)
        kfmax=KFMAXT
        kfaco=KFACOT
        isinel=ISINELT
        disbe=0.
        do j=1,kfmax
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
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine filt_h
c       filter out particles wanted to study and make
c       in order of proton,neutron,pba,nba,pi+,pi-,pi0,k-,k0-,sigma0,
c        sigma-,sigma+,sigma0ba,sigma-ba,sigma+ba,lamda,lamdaba,k0,k+,
c        cascade-,cascade-ba,cascade0,cascade0ba,Omega-,Omega+,Delta-,
c        Delta0,Delta+,Delta++,rho+,rho-,rho0,J/Psi,Psi',Upsilon,
c        Upsilon',x0c,x1c,x2c,D+,D+ba,D0,D0ba,lamdac+,sigmac0,sigmac+,
c        sigmac++,omega,k*+,K*0,D*+,D*+ba,D*0,D*0ba,Ds+,Ds+ba,B0,B0ba,
c        B+,B+ba,Bs0,Bs0ba,Bc+,Bc+ba,B*0,B*0ba,B*+,B*+ba,Bs*0,Bs*0ba,
c        Bc*+,Bc*+ba   ! 250420
c        (72 kinds of particle altogether)   ! 250420
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        iii=0
        jjj=0
        do i=1,kfmax
        kf=kfaco(i)
        do j=iii+1,n
        call ord_h(jjj,j,kf)
        enddo
        iii=jjj
        numb(i)=jjj
        enddo
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ord_h(ipi,j,kf)
c       make in order for particles with flavor code kf.
c       j : the particle wanted to order
c       kf: flavor code of j-th particle
c       ipi : j-th particle should order after ipi
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine copl_h(tt)
c       calculate position of center of mass of the non-freeze-out system 
c       distance of a particle from this cms is used to checke whether
c        it freezes out or not 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa19_h/coor(3)
        do ii=1,3
        coor(ii)=0.
        enddo
        samass=0.
        do 110 ii=1,n
        if(ishp(ii).eq.0)goto 110
        kf=k(ii,2)
        amass=pmas(pycomp(kf),1)
        samass=samass+amass
        do 100 jj=1,3
        coor(jj)=coor(jj)+amass*v(ii,jj)
100     continue
110     continue
        do ii=1,3
        coor(ii)=coor(ii)/max(0.14,samass)
        enddo
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ctlcre_h(lc,tc,tw,time)
c       create the initial collision list
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,ecsspn,ecsspm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/syspar_h/pio
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        m2=numb(2)   ! up to n (neutron)
        m4=numb(4)   ! up to nbar
        m7=numb(7)   ! up to pi0
        m9=numb(9)   ! up to k0-
        m17=numb(17)   ! up to Lambdabar
        m19=numb(19)   ! up to k+
        m25=numb(25)   ! up to omega+
        m29=numb(29)   ! up to Delta++
        m32=numb(32)   ! up to rho0
        m34=numb(34)   ! up to Psi'
c       m34=numb(kfmax-11)
c       subtract 11, since do not consider the rescattering of x0c, etc

        nctl=1
        do 500 l=1,nsa-1
        if(l.le.m19
     c   .or. (l.gt.m25 .and. l.le.m34))goto 300
c       consider only the reinteraction among nucleon, pion, kaon,
c        sigma, lambda, delta, rho and psi
        goto 500
300     do 600 l1=l+1,nsa

        if(nctl.gt.nsize)then
        write(9,*)'1 nsa,nctl,ddt=',nsa,nctl,ddt
        write(9,*)'size of array "nsize" needs to be extended'
        write(9,*)'error is serious,stop running'
        stop 30000
        endif

        if(l1.le.m19 
     c   .or. (l1.gt.m25 .and. l1.le.m34))goto 700
c       consider only the reinteraction among nucleon, pion, kaon,
c        sigma, lambda, delta,rho and psi
        goto 600
700     iflag=0
        call rsfilt_h(l,l1,iflag)
        if(iflag.eq.0)goto 600
        tc(nctl)=0.0
        call tcolij_h(l,l1,time,nctl,lc,tc,tw)
c170204 if(tc(nctl).gt.1.0e-7) nctl=nctl+1
c170204
        tci=tc(nctl)
c110504
c       if(tci.gt.1.0e-7)then
        if(tci.eq.0.0)goto 600
c       from 'tcolij' unsuccessfully, goto 600
        if(nctl.eq.1)then
        nctl=nctl+1
        goto 600
        endif
c110504
        do j1=1,nctl-1
        tcj=tc(j1)
        if(abs(tcj-tci).lt.ddt)goto 600
        enddo
        nctl=nctl+1
c       endif
c170204
600     continue
500     continue
c00623 if(tc(nctl).le.1.e-7) nctl=nctl-1   !Lei2023060
        nctl=nctl-1   !Lei2023060
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine scat_h(time,lc,tc,tw,ijkk,ijk,neve)
c       administrate hadron rescattering 
c       ijk: the run number
c       if ijkk=1 (infinite loop) give up the event 
c       neve: total number of runs
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/sa19_h/coor(3)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,ecsspn,ecsspm
        common/syspar_h/pio
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension pi(4),pj(4),pii(4),pjj(4),peo(4),pint(4),b(3)
        integer winel
        ijkk=0
        nctl0=nctl
        iii=1
10      if(iii.eq.1)goto 1000
        call copl_h(time)

c       find out the binary colli. with minimum collsion time
1000    call find_h(icp,tcp,lc,tc,tw,1)
        if(icp.eq.0)goto 100
c       icp=0 means the collision list is empty
c141104
        if(tcp.gt.1.e10)then
        do i1=icp+1,nctl
        do j1=1,5
        lc(i1-1,j1)=lc(i1,j1)
        enddo
        tc(i1-1)=tc(i1)
        tw(i1-1)=tw(i1)
        enddo
        nctl=nctl-1
        goto 10
        endif
c141104
        time0=time
        l=lc(icp,1)
        l1=lc(icp,2)
        kfa=ksa(l,2)    ! 060603
        kfb=ksa(l1,2)   ! 060603
        time=tcp
c       record this collision time
20      continue

c       perform classical Newton motion in Lab. system
        call his_h(time,lc,tc,tw,istop)
        if(istop.eq.1)goto 100
c       istop=1 means all particles have get out of considered volume

        pi(4)=psa(l,4)
        pj(4)=psa(l1,4)
        if(pi(4).lt.1.e-10)pi(4)=1.e-10   ! 031204
        if(pj(4).lt.1.e-10)pj(4)=1.e-10   ! 031204
        do i=1,3
        pi(i)=psa(l,i)
        pj(i)=psa(l1,i)
        b(i)=(pi(i)+pj(i))/(pi(4)+pj(4))
        enddo

c       boost to CMS frame of colliding pair
        ilo=0
        call lorntz(ilo,b,pi,pj)
        ss=pi(4)+pj(4)

        ww=rcsit
c       the cross section ratio of (ela.)/tot =1- rcsit
        rrlu=pyr(1)
        if(rrlu.gt.ww)then
        winel=0   ! ela.
        goto 640
        endif
700     winel=1   ! inela.
c       two particles with four-momentum pi and pj in CMS frame and flavor
c       ksa(l,2),ksa(l1,2) might go through inelastic reaction

        call coinel(l,l1,ss,b,pi,pj,icp,pii,pjj,lc,tc,tw,winel,ik3,ik4)
c       if winel=0 the inelastic reaction has not really happened 
c       if winel=1 inelastic collision happens, ik3 and ik4 are  
c        flavors of the scattered particles, pii and pjj are four-momentum of 
c        scattered particles in Lab frame, the two colliding particles are 
c        still with line numbers of l and l1 in the particle list

640     if(winel.ne.0)then   !!!

c       treat the inelastic collision 
        icp5=lc(icp,5)
        noinel(icp5)=noinel(icp5)+1
c       update particle list after inelastic collision
c        and truncates collision list correspondingly 
        call updpli(l,l1,icp,ss,pii,pjj,lc,tc,tw,winel,time,icp5)
        l=lc(icp,1)
        l1=lc(icp,2)
c       l and l1 are now the line numbers of scattered particles in 
c        particle list.

        else   !!!
c       calculate four-momentum of scattered particles after elastic 
c        collistion (pi and pj in CMS frame)
        call coelas_h(l,l1,ss,pi,pj)

c       update particle list for elastic scattering, pi and pj have been
c        boosted back to Lab fram 
        call updple_h(l,l1,b,pi,pj,time)
        noel=noel+1

        endif   !!!

c       update the collision list
        call updatl_h(l,l1,time,lc,tc,tw,winel,iii)
301     continue

c       if(nctl.le.1)goto 300
c        deltt=time-time0
c       deal with particle decay in transport processes
c        call decay(time,deltt,lc,tc,tw,iii,ijk)
c       ich2=0.   !!
c       do i1=1,nsa   !!
c       kf=ksa(i1,2)   !!
c       ich2=ich2+pychge(kf)   !!
c       enddo   !!

c       goto 100   ! it is actived temporally

300     iii=iii+1

        if(iii.gt.100*(nctl0))then
        write(9,*)'infinite loop may have happened in'
        write(9,*)'subroutine scat ijk=',ijk
c10/08/98       stop 'infinite loop occurs'
        ijk=ijk-1   ! 10/08/98
        ijkk=1  ! 10/08/98
c       if ijkk=1 (infinite loop) give up the event
        return   ! 10/08/98
        endif

        goto 10
100     continue
1200    continue
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine find_h(icp,tcp,lc,tc,tw,ico)
c       find out the binary collision with minimum collision time
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        icp=0
        tcp=20000.
        do i=1,nctl
        if(ico.eq.0)goto 100
        if(tc(i).le.1.0e-7)goto 241
100     if(tcp.lt.tc(i))goto 241
        icp=i
        tcp=tc(i)
241     continue
        enddo
        if(nctl.eq.0)icp=0
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine his_h(t1,lc,tc,tw,istop)   ! 231104
c       classical Newton motion in Lab. system
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa19_h/coor(3)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,ecsspn,ecsspm
        common/sa24/adj1(40),nnstop,non24,zstop
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        istop=1   ! 231104
        in=0   ! 231104
        r0=adj1(28)
        do 200 i=1,nsa
c       if(t1.le.tau(i))goto 100
c       do move particles which have not produced
        if(ishp(i).eq.1) goto 10
        in=in+1
        goto 100
10      aa=0.
        pp4=psa(i,4)
        do j=1,3
        vp=psa(i,j)/pp4
        vsa(i,j)=vsa(i,j)+vp*(t1-vsa(i,4))
        aa=aa+(vsa(i,j)-coor(j))**2
        enddo
c100505 vsa(i,4)=t1
        aa=sqrt(aa)
        if(aa.lt.r0) goto 300   ! 100 originally, 100505
c       if freeze-out already, deduct the distance between the last and 
c        current collisions
        do j=1,3
        vp=psa(i,j)/pp4
        vsa(i,j)=vsa(i,j)-vp*(t1-vsa(i,4))
        enddo
        ishp(i)=0
        do il=1,nctl
        if(lc(il,1).eq.i.or.lc(il,2).eq.i) tc(il)=0.
        enddo
        goto 200   !Lei2023060
300     vsa(i,4)=t1   ! 100505
100     continue
200     continue
        if(in.eq.nsa) return
        istop=0
        return
        end



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine coinel(l,l1,ss,b,pi,pj,icp,pii,pjj,lc,tc,tw,winel,
     &   ik1,ik2)

c       treat the inelastic collision
c       the inelastic processes considered are:

c        strangeness production reactions:
c        1. pion+ + pion- to k+ + k-
c        2. pion+ + pion- to k0 + k0-
c        3. pion+ + pion0 to k+ + k0-
c        4. pion- + pion0 to k- + k0
c        5. pion0 + pion0 to k+ + k-
c        6. pion0 + pion0 to k0 + k0-

c        7. pion+ + p to k+ + sigma+
c        8. pion+ + n to k+ + sigma0
c        9. pion+ + n to k+ + lambda
c        10. pion+ + n to k0 + sigma+
c        11. pion- + p to k+ + sigma-
c        12. pion- + p to k0 + lambda
c        13. pion- + p to k0 + sigma0
c        14. pion- + n to k0 + sigma-
c        15. pion0 + p to k+ + sigma0
c        16. pion0 + p to k+ + lambda
c        17. pion0 + p to k0 + sigma+
c        18. pion0 + n to k+ + sigma-
c        19. pion0 + n to k0 + lambda
c        20. pion0 + n to k0 + sigma0

c        21. pion+ + pba to k0- + lambdaba
c        22. pion+ + pba to k0- + sigma0ba
c        23. pion+ + pba to k- + sigma-ba
c        24. pion+ + nba to k0- + sigma-ba
c        25. pion- + pba to k- + sigma+ba
c        26. pion- + nba to k- + lambdaba
c        27. pion- + nba to k- + sigma0ba
c        28. pion- + nba to k0- + sigma+ba
c        29. pion0 + pba to k- + lambdaba
c        30. pion0 + pba to k- + sigma0ba
c        31. pion0 + pba to k0- + sigma+ba
c        32. pion0 + nba to k- + sigma-ba
c        33. pion0 + nba to k0- + lambdaba
c        34. pion0 + nba to k0- + sigma0ba

c        35. pion+ + sigma- to k+ + cascade-
c        36. pion- + lambda to k0 + cascade-
c        37. pion- + sigma+ to k+ + cascade-
c        38. pion- + sigma0 to k0 + cascade-
c        39. pion0 + lambda to k+ + cascade-
c        40. pion0 + sigma- to k0 + cascade-
c        41. pion0 + sigma0 to k+ + cascade-
c        42. pion+ + lambdaba to k0- + cascade-ba
c        43. pion+ + sigma+ba to k- + cascade-ba
c        44. pion+ + sigma0ba to k0- + cascade-ba
c        45. pion- + sigma-ba to k- + cascade-ba
c        46. pion0 + lambdaba to k- + cascade-ba
c        47. pion0 + sigma-ba to k0- + cascade-ba
c        48. pion0 + sigma0ba to k- + cascade-ba

c        strangeness exchange reactions:
c        49. k- + p to pion0 + lambda
c        50. k- + p to pion0 + sigma0
c        51. k- + p to pion- + sigma+
c        52. k- + p to pion+ + sigma-
c        53. k- + p to k+ + cascade- (strangeness production)
c        54. k- + n to pion- + sigma0
c        55. k- + n to pion- + lambda
c        56. k- + n to pion0 + sigma-
c        57. k- + n to k0 + cascade- (strangeness production)

c        58. k0- + p to pion+ + lambda
c        59. k0- + p to pion+ + sigma0
c        60. k0- + p to pion0 + sigma+
c        61. k0- + n to pion+ + sigma-
c        62. k0- + n to pion- + sigma+
c        63. k0- + n to pion0 + sigma0
c        64. k0- + n to pion0 + lambda
c        65. k0- + n to k+ + cascade- (strangeness production)

c        66. k+ + p- to pion0 + lambda-
c        67. k+ + p- to pion0 + sigma0-
c        68. k+ + p- to pion+ + sigma+ba
c        69. k+ + p- to pion- + sigma-ba
c        70. k+ + p- to k- + cascade-ba (strangeness production)

c        71. k+ + n- to pion+ + sigma0-
c        72. k+ + n- to pion+ + lambda-
c        73. k+ + n- to pion0 + sigma-ba
c        74. k+ + n- to k0- + cascade-ba(strangeness production)

c        75. k0 + p- to pion- + lambda-
c        76. k0 + p- to pion- + sigma0-
c        77. k0 + p- to pion0 + sigma+ba
c        78. k0 + n- to pion- + sigma-ba
c        79. k0 + n- to pion+ + sigma+ba
c        80. k0 + n- to pion0 + sigma0-
c        81. k0 + n- to pion0 + lambda-
c        82. k0 + n- to k- + cascade-ba (strangeness production)

c        83. k- + lambda to pion0 + cascade-
c        84. k- + sigma+ to pion+ + cascade-
c        85. k- + sigma- to pion- + cascade-
c        86. k- + sigma0 to pion0 + cascade-
c        87. k0- + lambda to pion+ + cascade-
c        88. k0- + sigma0 to pion+ + cascade-
c        89. k0- + sigma- to pion0 + cascade-
c        90. k+ + lambda- to pion0 + cascade-ba
c        91. k+ + sigma+ba to pion- + cascade-ba
c        92. k+ + sigma-ba to pion+ + cascade-ba
c        93. k+ + sigma0- to pion0 + cascade-ba
c        94. k0 + lambda- to pion- + cascade-ba
c        95. k0 + sigma0- to pion- + cascade-ba
c        96. k0 + sigma-ba to pion0 + cascade-ba
c        97. pion+ + sigma- to k0 + cascade0
c        98. pion+ + sigma0 to k+ + cascade0
c        99. pion+ + lambda0 to k+ + cascade0
c       100. pion- + sigma+ to k0 + cascade0
c       101. pion0 + sigma+ to k+ + cascade0
c       102. pion0 + sigma0 to k0 + cascade0
c       103. pion0 + lambda to k0 + cascade0
c       104. pion+ + sigma+ba to k0- + cascade0-
c       105. pion- + sigma-ba to k0- + cascade0-
c       106. pion- + sigma0- to k- + cascade0-
c       107. pion- + lambda- to k- + cascade0-
c       108. pion0 + sigma+- to k- + cascade0-
c       109. pion0 + sigma0- to k0- + cascade0-
c       110. pion0 + lambda- to k0- + cascade0-
c       111. k- + sigma+ to pion0 + cascade0
c       112. k- + sigma0 to pion0 + cascade-
c       113. k- + lambda to pion- + cascade0
c       114. k0- + sigma+ to pion+ + cascade0
c       115. k0- + sigma- to pion- + cascade0
c       116. k0- + sigma0 to pion0 + cascade0
c       117. k0- + lambda to pion0 + cascade0
c       118. k+ + sigma+ba to pion0 + cascade0-
c       119. k+ + sigma0- to pion+ + cascade0-
c       120. k+ + lambda- to pion+ + cascade0-
c       121. k+ + cascade-ba to pion+ + omiga-ba
c       122. k0 + sigma-ba to pion+ + cascade0ba
c       123. k0 + sigma0- to pion0 + cascade0-
c       124. k0 + lambda- to pion0 + cascade0ba
c       125. k- + p to k0 + cascade0
c       126. k0- + p to k+ + cascade0
c       127. k0- + n to k0 + cascade0
c       128. k+ + p- to k0- + cascade0ba
c       129. k0 + p- to k- + cascade0-
c       130. k0 + n- to k0- + cascade0-
c       131. pion+ + cascade- to k+ + omiga-
c       132. pion0 + cascade- to k0 + omiga-
c       133. pion- + cascade-ba to k- + omiga-ba
c       134. pion0 + cascade-ba to k0- + omiga-ba
c       135. pion- + cascade0 to k0 + omiga-
c       136. pion0 + cascade0 to k+ + omiga-
c       137. pion+ + cascade0- to k0- + omiga-ba
c       138. pion0 + cascade0- to k- + omiga-ba
c       139. k- + cascade- to pion- + omiga-
c       140. k0- + cascade- to pion0 + omiga-
c       141. k- + cascade0 to pion0 + omiga-
c       142. k0- + cascade0 to pion+ + omiga-
c       143. k+ + cascade-ba to pion+ + omiga-ba
c       144. k0 + cascade-ba to pion0 + omiga-ba
c       145. k+ + cascade0- to pion0 + omiga-ba
c       146. k0 + cascade0- to pion- + omiga-ba

c       147. pion- + p to delta- + pion+
c       148. pion- + p to rho0 + n
c       149. pion- + p to rho- + p
c       150. pion- + p to delta+ + pion-
c       151. pion- + p to delta0 + pion0
c       152. pion- + n to delta- + pion0
c       153. pion- + n to rho- + n
c       154. pion- + n to delta0 + pion-
c       155. pion+ + p to delta++ + pion0
c       156. pion+ + p to delta+ + pion+
c       157. pion+ + p to rho+ + p
c       158. pion+ + n to delta++ + pion-
c       159. pion+ + n to delta0 + pion+
c       160. pion+ + n to delta+ + pion0
c       161. pion+ + n to rho0 + p
c       162. pion+ + n to rho+ + n
c       163. pion0 + p to delta0 + pion+
c       164. pion0 + p to delta++ + pion-
c       165. pion0 + p to rho+ + n
c       166. pion0 + p to rho0 + p
c       167. pion0 + p to delta+ + pion0
c       168. pion0 + n to delta+ + pion-
c       169. pion0 + n to delta- + pion+
c       170. pion0 + n to delta0 + pion0
c       171. pion0 + n to rho0 + n
c       172. pion0 + n to rho- + p
c       173. p + p to delta+ + p
c       174. p + p to delta++ + n
c       175. p + n to delta+ + n
c       176. p + n to delta0 + p
c       177. n + n to delta0 + n
c       178. n + n to delta- + p

c       179. J/psi + n to lamdac + Dba
c       180. J/psi + n to sigmac + Dba
c       181. J/psi + n to sigmac0 + D0ba
c       182. J/psi + p to lamdac + D0ba
c       183. J/psi + p to sigmac + D0ba
c       184. J/psi + p to sigmac++ + Dba
c       185. J/psi + pion+ to D + D*0ba
c       186. J/psi + pion0 to D0 + D*0ba
c       187. J/psi + pion0 to D + D*ba
c       188. J/psi + pion- to D0 + D*ba
c       189. J/psi + rho+ to D + D0ba
c       190. J/psi + rho0 to D0 + D0ba
c       191. J/psi + rho0 to D + Dba
c       192. J/psi + rho- to D0 + Dba
c       193. psi' + n to lamdac + Dba
c       194. psi' + n to sigmac + Dba
c       195. psi' + n to sigmac0 + D0ba
c       196. psi' + p to lamdac + D0ba
c       197. psi' + p to sigmac + D0ba
c       198. psi' + p to sigmac++ + Dba
c       199. psi' + pion+ to D + D*0ba
c       200. psi' + pion0 to D0 + D*0ba

c       reverse reactions, after '201'
c       201. pion+ + pion- to k+ + k- (r)
c       202. pion0 + pion0 to k+ + k-(r)
c       203. pion+ + pion0 to k+ + k0-(r)
c       204. pion- + pion0 to k- + k0(r)
c       205. pion+ + pion- to k0 + k0-(r)
c       206. pion0 + pion0 to k0 + k0-(r)
c       207. k+ + sigma+ to pion+ + p
c       208. k+ + sigma- to pion- + p
c       209. k+ + sigma- to pion0 + n
c       210. k+ + sigma0 to pion+ + n
c       211. k+ + sigma0 to pion0 + p
c       212. k+ + lambda0 to pion+ + n
c       213. k+ + lambda0 to pion0 + p
c       214. k+ + sigma+ to pion+ + n
c       215. k+ + sigma+ to pion0 + p
c       216. k0 + sigma- to pion- +n
c       217. k0 + sigma0 to pion- + p
c       218. k0 + sigma0 to pion0 + n
c       219. k0 + lambda0 to pion- + p
c       220. k0 + lambda0 to pion0 + n
c       221. k- + sigma+ba to pion- + pba
c       222. k- + sigma-ba to pion+ + pba
c       223. k- + sigma-ba to pion0 + nba
c       224. k- + sigma0ba to pion- + nba
c       225. k- + sigma0ba to pion0 + pba
c       226. k- + lambda0ba to pion- + nba
c       227. k- + lambda0ba to pion0 + pba
c       228. k0ba + sigma+ba to pion- + nba
c       229. k0ba + sigma+ba to pion0 + pba
c       230. k0- + sigma-ba to pion+ + nba
c       231. k0- + sigma0- to pion+ + pba
c       232. k0- + sigma0- to pion0 + nba
c       233. k0- + lambda0- to pion+ + pba
c       234. k0- + lambda0- to pion0 + nba
c       235. k++ cascade- to pi+ sigma-
c       236. k+ + cascade- to pion- + sigma+
c       237. k+ + cascade- to pion0 + lambda0
c       238. k+ + cascade- to pion0 + sigma0
c       239. k- + cascade-ba to pion+ + sigma+ba
c       240. k- + cascade-ba to pion- + sigma-ba
c       241. k- + cascade-ba to pion0 + lambda0-
c       242. k- + cascade-ba to pion0 + sigma0-
c       243. k0 + cascade- to pion- + lambda0
c       244. k0 + cascade- to pion- + sigma0
c       245. k0 + cascade- to pion0 + sigma-
c       246. k0- + cascade-ba to pion+ + lambda0-
c       247. k0- + cascade-ba to pion+ + sigma0-
c       248. k0- + cascade-ba to pion0 + sigma-ba
c       249. pion+ + sigma- to k- + p
c       250. pion+ + sigma- to k0- + n
c       251. pion+ + sigma0 to k0- + p
c       252. pion+ + lambda0 to k0- + p
c       253. k+ + cascade- to k- + p
c       254. pion- + sigma+ to k- + p
c       255. pion- + sigma+ to k0- + n
c       256. pion- + sigma0 to k- + n
c       257. k0 + cascade- to k- + n
c       258. pion- + lambda0 to k- + n
c       259. pion0 + sigma+ to k0- + p
c       260. pion0 + sigma- to k- + n
c       261. pion0 + sigma0 to k- + p
c       262. pion0 + sigma0 to k0- + n
c       263. pion0 + lambda0 to k- + p
c       264. pion0 + lambda0 to k0- + n
c       265. k+ + cascade- to k0- + n
c       266. pion+ + sigma+ba to k+ + pba
c       267. pion+ + sigma+ba to k0 + nba
c       268. pion+ + sigma0- to k+ + nba
c       269. pion+ + lambda0- to k+ + nba
c       270. k- + cascade-ba to k+ + pba
c       271. pion- + sigma-ba to k+ + pba
c       272. pion- + sigma-ba to k0 + nba
c       273. pion- + sigma0- to k0 + pba
c       274. k0- + cascade-ba to k+ + nba
c       275. pion- + lambda0- to k0 + pba
c       276. pion0 + sigma+- to k0 + pba
c       277. pion0 + sigma-ba to k+ + nba
c       278. pion0 + sigma0- to k+ + pba
c       279. pion0 + sigma0- to k0 + nba
c       280. pion0 + slambda0- to k+ + pba
c       281. pion0 + lambda0- to k0 + nba
c       282. k- + cascade-ba to k0 + nba
c       283. pion+ + cascade- to k- + sigma+
c       284. pion+ + cascade- to k0- + lambda0
c       285. pion+ + cascade- to k0- + sigma0
c       286. pion- + cascade- to k- + sigma-
c       287. pion0 + cascade- to k- + lambda0
c       288. pion0 + cascade- to k- + gigma0
c       289. pion0 + cascade- to k0- + sigma-
c       290. pion+ + cascade-ba to k+ + sigma-ba
c       291. pion- + cascade-ba to k+ + sigma+-
c       292. pion- + cascade-ba to k0 + lambda0-
c       293. pion- + cascade-ba to k0 + sigma0-
c       294. pion0 + cascade-ba to k+ + lambda0-
c       295. pion0 + cascade-ba to k+ + gigma0-
c       296. pion0 + cascade-ba to k0 + sigma-ba
c       297. k0 + cascade0 to pion+ + sigma-
c       298. k+ + cascade0 to pion+ + sigma0
c       299. k+ + cascade0 to pion+ + lambda
c       300. k0 + cascade0 to pion- + sigma+
c       301. k+ + cascade0 to pion0 + sigma+
c       302. k0 + cascade0 to pion0 + sigma0
c       303. k0 + cascade0 to pion0 + lambda
c       304. k-0 + cascade0- to pion+ + sigma-ba
c       305. k0- + cascade0- to pion- + sigma-ba
c       306. k- + cascade0- to pion- + sigma0-
c       307. k- + cascade0- to pion- + lambda-
c       308. k- + cascade0- to pion0 + sigma+ba
c       309. k0- + cascade0- to pion0 + sigma0-
c       310. k0- + cascade0- to pion0 + lambda-
c       311. pion0 + cascade0 to k- + sigma+
c       312. pion- + cascade0 to k- + sigma0
c       313. pion- + cascade0 to k- + lambda
c       314. pion+ + cascade0 to k0- + sigma+
c       315. pion- + cascade0 to k0- + sigma-
c       316. pion0 + cascade0 to k0- + sigma0
c       317. pion0 + cascade0 to k0- + lambda
c       318. pion0 + cascade0- to k+ + sigma+ba
c       319. pion+ + cascade0- to k+ + sigma0-
c       320. pion+ + cascade0- to k+ + lambda-
c       321. pion- + cascade0- to k0 + sigma+ba
c       322. pion+ + cascade0- to k0 + sigma-ba
c       323. pion0 + cascade0- to k0 + sigma0-
c       324. pion0 + cascade0- to k0 + lambda-
c       325. k0 + cascade0 to k- + p
c       326. k+ + cascade0 to k0- + p
c       327. k0 + cascade0 to k0- + n
c       328. k0- + cascade0- to k+ + p-
c       329. k- + cascade0- to k0 + p-
c       330. k0- + cascade0- to k0 + n-
c       331. k+ + omiga- to pion+ + cascade-
c       332. k0 + omiga- to pion0 + cascade-
c       333. k- + omiga-ba to pion- + cascade-ba
c       334. k0- + omiga-ba to pion0 + cascade-ba
c       335. k0 + omiga- to pion- + cascade0
c       336. k+ + omiga- to pion0 + cascade0
c       337. k0- + omiga-ba to pion+ + cascade0-
c       338. k- + omiga-ba to pion0 + cascade0-
c       339. pion- + omiga- to k- + cascade-
c       340. pion0 + omiga- to k0- + cascade-
c       341. pion0 + omiga- to k- + cascade0
c       342. pion+ + omiga- to k0- + cascade0
c       343. pion+ + omiga-ba to k+ + cascade-ba
c       344. pion0 + omiga-ba to k0 + cascade-ba
c       345. pion0 + omiga-ba to k+ + cascade0-
c       346. pion- + omiga-ba to k0 + cascade0-
c       347. pion+ + delta- to pion- + p
c       348. pion+ + delta- to pion0 + n
c       349. pion+ + delta0 to pion+ + n
c       350. pion+ + delta0 to pion0 + p
c       351. pion+ + delta+ to pion+ + p
c       352. pion0 + delta++ to pion+ p
c       353. pion0 + delta+ to pion0 + p
c       354. pion0 + delta+ to pion+ + n
c       355. pion0 + delta0 to pion0 + n
c       356. pion0 + delta0 to pion- + p
c       357. pion0 + delta- to pion- + n
c       358. pion- + delta++ to pion0 + p
c       359. pion- + delta++ to pion+ + n
c       360. pion- + delta+ to pion- + p
c       361. pion- + delta+ to pion0 + n
c       362. pion- + delta0 to pion- + n
c       363. rho0 + n to pion- + p
c       364. rho0 + n to pion0 + n
c       365. rho- + n to pion- + n
c       366. rho+ + n to pion0 + p
c       367. rho+ + n to pion+ + n
c       368. rho0 + p to pion0 + p
c       369. rho0 + p to pion+ + n
c       370. rho- + p to pion0 + n
c       371. rho- + p to pion- + p
c       372. rho+ + p to pion+ + p
c       373. delta++ + n to p + p
c       374. delta+ + n to p + n
c       375. delta+ + p to p + p
c       376. delta0 + p to p + n
c       377. delta0 + n to n + n
c       378. delta- + p to n + n

c       follows are direct reactions again
c       379. psi' + pion0 to D + D*ba
c       380. psi' + pion- to D0 + D*ba
c       381. psi' + rho+ to D + D0ba
c       382. psi' + rho0 to D0 + D0ba
c       383. psi' + rho0 to D + Dba
c       384. psi' + rho- to D0 + Dba

c       593. lambda- + p to K*+ + omiga
c       594. lambda- + n to K*0 + omiga
c       595. sigma0- + p to K*+ + omiga
c       596. sigma0- + n to K*0 + omiga
c       597. p- + p to rho0 + omiga
c       598. p- + n to rho- + omiga
c       599. n- + p to rho+ + omiga
c       600. n- + n to rho0 + omiga
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
c       note the name of the arraies in 'sa1' in this subroutine
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension pi(4),pj(4),pii(4),pjj(4),b(3)
        integer winel
        jjj=0
        ww=tw(icp)
        kl=k(l,2)
        kl1=k(l1,2)
        am01=pmas(pycomp(kl),1)
        am02=pmas(pycomp(kl1),1)
c       the following statements guarantees the line number in particle 
c        list of pion, kion, rho, and J/psi should be first component 
c        of prod()
c       tw(i): the cross section ratio of (i-th inela.)/tot
        if(abs(kl).eq.211 .or. abs(kl).eq.111 .or. abs(kl).eq.321
     &   .or. abs(kl).eq.311.or. abs(kl).eq.213.or. abs(kl).eq.113
     &   .or. kl.eq.443 .or. kl.eq.100443)then
        idpl=1   ! meson
        else
        idpl=3   ! baryon
        endif
c       the meson in a colliding pair should be first component of prod(), 
c        but not for the case of collision occurs between two baryons (or 
c        two mesons)
        if(idpl.eq.1)call prod(l,l1,kl,kl1,ss,icp,lc,tc,tw)
        if(idpl.eq.3)call prod(l1,l,kl1,kl,ss,icp,lc,tc,tw)
        ik1=lc(icp,3)
        ik2=lc(icp,4)
        w1=tw(icp)/rcsit
c       1/rcsit : the cross section ratio of tot/(inela.)
        ww=1.
        if(pyr(1).gt.w1)then
        winel=0
c       treated as elastic then
        return
        endif
        icp5=lc(icp,5)
        if(icp5.ge.1593)goto 400   ! nonsense indeed
c100    fi1=atan2(pi(2),pi(1))
c       cta1=atan2(sqrt(pi(1)**2+pi(2)**2),pi(3))
100     fi1=pyangl(pi(1),pi(2))
        cta1=pyangl(pi(3),sqrt(pi(1)**2+pi(2)**2))
        cfi1=cos(fi1)
        sfi1=sin(fi1)
        ccta1=cos(cta1)
        scta1=sin(cta1)
        am1=pmas(pycomp(ik1),1)
        am2=pmas(pycomp(ik2),1)
        pp=(ss*ss-(am1+am2)**2)*(ss*ss-(am1-am2)**2)/(4.*ss*ss)
        if(pp.lt.0.)pp=1.e-10
        pp=sqrt(pp)
        pii(4)=(ss*ss+am1**2-am2**2)/(2.*ss)
c       energy of one particle (between two) after scattering
        fis=2.*3.1415926*pyr(1)
        cfis=cos(fis)
        sfis=sin(fis)
        call cosin(am01,am1,ss,pi,pp,cctas)
c200    cctas=2*pyr(1)-1.
        sctas=sqrt(1.-cctas*cctas)
111     continue
        pii(1)=cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis
        pii(2)=sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis
        pii(3)=ccta1*cctas-scta1*sctas*cfis
        pii(1)=pp*pii(1)
        pii(2)=pp*pii(2)
        pii(3)=pp*pii(3)
        do i=1,3
        pjj(i)=0.-pii(i)
        enddo
        pjj(4)=ss-pii(4)
        if(pii(4).lt.0. .or. pjj(4).lt.0.)then
        write(9,*)'error may happen here in subroutine'
        write(9,*)'coinel(),energy is negative'
        if(pii(4).lt.0.)pii(4)=1.e-18   ! 050105
        if(pjj(4).lt.0.)pjj(4)=1.e-18   ! 050105
        endif
        ilo=1
        call lorntz(ilo,b,pii,pjj)
        if(pii(4).lt.0. .or. pjj(4).lt.0.)then
        write(9,*)'error may happen here in subroutine'
        write(9,*)'coinel(),energy is negative'
        if(pii(4).lt.0.)pii(4)=1.e-18   ! 050105
        if(pjj(4).lt.0.)pjj(4)=1.e-18   ! 050105
        endif
c       if(jjj.eq.1)goto 300
400     return
        end



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine cosin(am01,am1,eij,pi,pp,cctas)
c       calculate cos(theta_s)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        dimension pi(4)
c       d=3.65*(eij-am1-am2)
c       if(d.lt.1.e-10)return
        pt=0.5   ! 100111
        a=min(10.3,1./(1.12*pt)/(1.12*pt))
c       d6=d**6
c       b=d6*a/(1.+d6)
        b=a   ! 100111 changed from 10.3 to 'a'
c       if(b.lt.1.e-20)then
c       b=1.e-20
c       endif
        pm2=pi(1)**2+pi(2)**2+pi(3)**2
        pm=dsqrt(pm2)
        em=dsqrt(pm*pm+dble(am01*am01))
        em1=dsqrt(dble(pp*pp)+dble(am1*am1))
        tmin=am01**2+am1**2-2*(em*em1+pm*pp)
        tmax=am01**2+am1**2-2*(em*em1-pm*pp)

        cc=pyr(1)

        abt=dexp(dmax1(-7.0D2,b*tmin))
        abu=dexp(dmax1(-7.0D2,b*tmax))

        tt1=dlog(cc*abu+(1.-cc)*abt)

        tt=tt1/b
        cctas=(0.5*(tt-am01**2-am1**2)+em*em1)/(pm*pp)

c050510 if(abs(cctas).gt.1)cctas=sign(1.,cctas)
        if(dabs(cctas).gt.1.d0)cctas=dsign(1.d0,cctas)
        return
        end



c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine prod(l,l1,kl,kl1,ss,icp,lc,tc,tw)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (given inela.)/tot
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,ecsspn,ecsspm
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/count_h/isinel(600)
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        integer fact1,fact2
        para13=PARAM(1)*0.1*PARAM(6)*0.8
c       p-p annihilation cross section is assumed to be equal to 0.8*(
c        total inelastic cross section)
        ioo=0
        ilo=1
        ilo1=1
        ilo2=1
        ilo3=1
        ilo4=1
        ilo5=1
        ilo6=1
        ilo7=1
        ilo8=1
        ilo9=1

        nchargei=PYCHGE(kl)+PYCHGE(kl1)
c p-p ------------------------>
        if(kl.eq. 2212.and. kl1.eq.2212)then
        call ppdelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
c p+n ----------------------->
        elseif((kl.eq. 2212.and. kl1.eq.2112).or.
     &         (kl.eq. 2112.and. kl1.eq.2212))then
        call pndelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
c n+n ------->
        elseif(kl.eq. 2112.and. kl1.eq.2112)then
        call nndelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion+ + pion-
        if((kl.eq.211 .and. kl1.eq.-211)
     c   .or.(kl.eq.-211 .and. kl1.eq.211))then
        if(isinel(1).eq.0)then
        fact1=0
        goto 101
        endif
        fact1=1
        ik3=-321
        ik4=321
        ic3=1
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1.eq.0)fact1=0.
101     if(isinel(2).eq.0)then
        fact2=0
        goto 102
        endif
        fact2=1
        ik5=-311
        ik6=311
        ic5=2
        call spipi(ik5,ik6,ss,ilo2)
        if(ilo2.eq.0)fact2=0.
102     fact=fact1+fact2
        if(fact1.eq.0. .and. fact2.eq.0.)goto 13
c       if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(pyr(1).gt.fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*sig/cspipi
        goto 10
        endif

c       pion+ + pion0
        if((kl.eq.211 .and. kl1.eq.111)
     c   .or.(kl.eq.111 .and. kl1.eq.211))then
        if(isinel(3).eq.0)goto 13
        ik3=-311
        ik4=321
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=3
        tw(icp)=sig/cspipi
        goto 10
        endif
        
c       pion- + pion0
        if((kl.eq.-211 .and. kl1.eq.111)
     c   .or.(kl.eq.111 .and. kl1.eq.-211))then
        if(isinel(4).eq.0)goto 13
        ik3=311
        ik4=-321
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=4
        tw(icp)=sig/cspipi
        goto 10
        endif

c       pion0 + pion0
        if(kl.eq.111 .and. kl1.eq.111)then
        if(isinel(5).eq.0)then
        fact1=0.
        goto 103
        endif
        fact1=1.
        ik3=-321
        ik4=321
        ic3=5
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1.eq.0)fact1=0.
103     if(isinel(6).eq.0)then
        fact2=0.
        goto 104
        endif
        fact2=1.
        ik5=-311
        ik6=311
        ic5=6
        call spipi(ik5,ik6,ss,ilo2)
        if(ilo2.eq.0)fact2=0.
104     fact=fact1+fact2
        if(fact1.eq.0. .and. fact2.eq.0.)goto 13
c       if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(pyr(1).gt.fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*sig/cspipi
        goto 10
        endif

c       pion+ + p

500     if(kl.eq.211 .and. kl1.eq.2212)then
        call pip2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion+ + n
        if(kl.eq.211 .and. kl1.eq.2112)then
        call pin1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif


c       pion- + p
        if(kl.eq.-211 .and. kl1.eq.2212)then
        call pip1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion- + n
        if(kl.eq.-211 .and. kl1.eq.2112)then
        call pin3(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion0 + p
        if(kl.eq.111 .and. kl1.eq.2212)then
        call pip3(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion0 + n
        if(kl.eq.111 .and. kl1.eq.2112)then
        call pin2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion+ + pba
        if(kl.eq.211 .and. kl1.eq.-2212)then
        if(isinel(21).eq.0)then
        si1=0.
        goto 319
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of pion+ + pba to k0- + lambdaba
c       cross section of pion+ + pba to k0- + lambdaba is assumed to be 
c        equal to pion- + p to k0 + lambda = s1724
319     if(isinel(22).eq.0)then
        si2=0.
        goto 320
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
c       cross section of pion+ + pba to k0- + sigma0-
320     if(isinel(23).eq.0)then
        si3=0.
        goto 321
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)
c       cross section of pion+ + pba to k- + sigma-ba
321     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-311
        ik2=-3122
        ic=21
c       pion+ + pba to k0- + lambda-
        goto 322
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=-3212
        ic=22
c       pion+ + pba to k0- + sigma0-
        goto 322
        endif
        ik1=-321
        ik2=-3112
        ic=23
c       pion+ + pba to k- + sigma-ba
322     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion+ + nba to k0- + sigma-ba
        if(kl.eq.211 .and. kl1.eq.-2112)then
        if(isinel(24).eq.0)goto 13
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3112),1)
        tw(icp)=s1724(ss,ilo,0,the)/cspin/10.
c       cross section of pion+ + nba to k0- + sigma-ba is assumed to be 
c        equal to pion- + n to k + y (isotropic averaged)
        if(ilo.eq.0)goto 13
        lc(icp,3)=-311
        lc(icp,4)=-3112
        lc(icp,5)=24
        goto 10
        endif

c       pion- + pba to k- + sigma+ba
        if(kl.eq.-211 .and. kl1.eq.-2212)then
        if(isinel(25).eq.0)goto 13
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3222),1)
        tw(icp)=s1724(ss,ilo,0,the)/cspin/10.
c       cross section of pion- + pba to k- + sigma+ba is assumed to be 
c        equal to pion+ + p to k + y (isotropic averaged)
        if(ilo.eq.0)goto 13
        lc(icp,3)=-321
        lc(icp,4)=-3222
        lc(icp,5)=25
        goto 10
        endif

c       pion- + nba
        if(kl.eq.-211 .and. kl1.eq.-2112)then
        if(isinel(26).eq.0)then
        si1=0.
        goto 323
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of pion- + nba to k- + y- is assumed to be 
c        equal to pion+ + n to k + y
323     if(isinel(27).eq.0)then
        si2=0.
        goto 324
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
324     if(isinel(28).eq.0)then
        si3=0.
        goto 325
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
325     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3122
        ic=26
c       pion- + nba to k- + lambdaba
        goto 326
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-321
        ik2=-3212
        ic=27
c       pion- + nba to k- + sigma0ba
        goto 326
        endif
        ik1=-311
        ik2=-3222
        ic=28
c       pion- + nba to k0- + sigma+ba
326     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion0 + pba
        if(kl.eq.111 .and. kl1.eq.-2212)then
        if(isinel(29).eq.0)then
        si1=0.
        goto 327
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
327     if(isinel(30).eq.0)then
        si2=0.
        goto 328
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
328     if(isinel(31).eq.0)then
        si3=0.
        goto 329
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
329     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3122
        ic=29
c       pion0 + pba to k- + lambda-
        goto 330
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-321
        ik2=-3212
        ic=30
c       pion0 + pba to k- + sigma0-
        goto 330
        endif
        ik1=-311
        ik2=-3222
        ic=31
c       pion0 + pba to k0- + sigma+ba
330     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion0 + nba
        if(kl.eq.111 .and. kl1.eq.-2112)then
        if(isinel(32).eq.0)then
        si1=0.
        goto 331
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3112),1)
        si1=s1724(ss,ilo1,0,the)
331     if(isinel(33).eq.0)then
        si2=0.
        goto 332
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3122),1)
        si2=s1724(ss,ilo2,0,the)
332     if(isinel(34).eq.0)then
        si3=0.
        goto 333
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)
333     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3112
        ic=32
c       pion0 + nba to k- + sigma-ba
        goto 334
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=-3122
        ic=33
c       pion0 + nba to k0- + lambdaba
        goto 334
        endif
        ik1=-311
        ik2=-3212
        ic=34
c       pion0 + nba to k0- + sigma0ba
334     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion+ + sigma-
        if(kl.eq.211 .and. kl1.eq.3112)then
        if(isinel(35).eq.0)then
        si1=0.
        goto 701
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of pion + y to kaon + cascade is assumed to be 
c        equal to pion + n to kaon + y,but take the different of
c        threshold energy into account

701     if(isinel(97).eq.0)then
        si2=0.
        goto 660
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

660     if(isinel(249).eq.0)then
        si3=0.
        goto 661
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
661     if(isinel(250).eq.0)then
        si4=0.
        goto 702
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
702     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3312
        ic=35
c       pion+ + sigma- to k+ + cascade-
        goto 703
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=311
        ik2=3322
        ic=97
c       pion+ + sigma- to k0 + cascade0
        goto 703
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=-321
        ik2=2212
        ic=249
c       pion+ + sigma- to k- + p
        goto 703
        endif
        ik1=-311
        ik2=2112
        ic=250
c       pion+ + sigma- to k0- + n
703     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + sigma-ba
        if(kl.eq.-211 .and. kl1.eq.-3112)then
        if(isinel(45).eq.0)then
        si1=0.
        goto 7011
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of pion + y to kaon + cascade is assumed to be 
c        equal to pion + n to kaon + y,but take the different of
c        threshold energy into account

7011    if(isinel(105).eq.0)then
        si2=0.
        goto 6601
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

6601    if(isinel(271).eq.0)then
        si3=0.
        goto 6611
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
6611    if(isinel(272).eq.0)then
        si4=0.
        goto 7021
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
7021    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0
     &   .and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3312
        ic=45
c       pion- + sigma-ba to k- + cascade-ba
        goto 7031
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=-3322
        ic=105
c       pion- + sigma-ba to k0- + cascade0-
        goto 7031
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=321
        ik2=-2212
        ic=271
c       pion- + sigma-ba to k+ + pba
        goto 7031
        endif
        ik1=311
        ik2=-2112
        ic=272
c       pion- + sigma-ba to k0 + nba
7031    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + sigma+ 
        if(kl.eq.-211 .and. kl1.eq.3222)then
        if(isinel(37).eq.0)then
        si1=0.
        goto 704
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)

704     if(isinel(100).eq.0)then
        si2=0.
        goto 663
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

663     if(isinel(254).eq.0)then
        si3=0.
        goto 664
        endif
        ik1=-321
        ik2=2212
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
664    if(isinel(255).eq.0)then
        si4=0.
        goto 705
        endif
        ik1=-311
        ik2=2112
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
705     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3312
        ic=37
c       pion- + sigma+ to k+ + cascade-
        goto 706
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=311
        ik2=3322
        ic=100
c       pion- + sigma+ to k0 + cascade0
        goto 706
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=-321
        ik2=2212
        ic=254
c       pion- + sigma+ to k- + p
        goto 706
        endif
        ik1=-311
        ik2=2112
        ic=255
c       pion- + sigma+ to k0- + n
706     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion+ + sigma+bar
        if(kl.eq.211 .and. kl1.eq.-3222)then
        if(isinel(43).eq.0)then
        si1=0.
        goto 7041
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)

7041    if(isinel(104).eq.0)then
        si2=0.
        goto 6631
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

6631    if(isinel(266).eq.0)then
        si3=0.
        goto 6641
        endif
        ik1=321
        ik2=-2212
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
6641    if(isinel(267).eq.0)then
        si4=0.
        goto 7051
        endif
        ik1=311
        ik2=-2112
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
7051    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3312
        ic=43
c       pion+ + sigma+bar to k- + cascade-bar
        goto 7061
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=-3322
        ic=104
c       pion+ + sigma+bar to k0- + cascade0bar
        goto 7061
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=321
        ik2=-2212
        ic=266
c       pion+ + sigma+bar to k+ + pbar
        goto 7061
        endif
        ik1=311
        ik2=-2112
        ic=267
c       pion+ + sigma+bar to k0 + nbar
7061    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif
c       pion- + sigma0 
        if(kl.eq.-211 .and. kl1.eq.3212)then
        if(isinel(38).eq.0)then
        si1=0.
        goto 1681
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
1681    if(isinel(256).eq.0)then
        si2=0.
        goto 1682
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(2112),1)
        ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
1682    if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=311
        ik2=3312
        ic=38
c       pion- + sigma0 to k0 + cascade-
        goto 683
        endif
        ik1=-321
        ik2=2112
        ic=256
c       pion- + sigma0 to k- + n
683     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif


c       pion+ + sigma0 bar
        if(kl.eq.211 .and. kl1.eq.-3212)then
        if(isinel(44).eq.0)then
        si1=0.
        goto 1781
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
1781    if(isinel(268).eq.0)then
        si2=0.
        goto 1782
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(-2112),1)
        ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
1782    if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-311
        ik2=-3312
        ic=44
c       pion+ + sigma0bar to k0- + cascade-bar
        goto 6831
        endif
        ik1=321
        ik2=-2112
        ic=268
c       pion+ + sigma0bar to k+ + nbar
6831    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif


c       pion0 + sigma0 
        if(kl.eq.111 .and. kl1.eq.3212)then
        if(isinel(41).eq.0)then
        si1=0.
        goto 710
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)

710     if(isinel(102).eq.0)then
        si2=0.
        goto 666
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

666     if(isinel(261).eq.0)then
        si3=0.
        goto 667
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
667     if(isinel(262).eq.0)then
        si4=0.
        goto 711
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
711     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3312
        ic=41
c       pion0 + sigma0 to k+ + cascade-
        goto 712
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=311
        ik2=3322
        ic=102
c       pion0 + sigma0 to k0 + cascade0
        goto 712
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=-321
        ik2=2212
        ic=261
c      pion0 + sigma0 to k- + p
        goto 712
        endif
       ik1=-311
        ik2=2112
        ic=262
c       pion0 + sigma0 to k0- + n
712     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif
c       pion0 + sigma0bar
        if(kl.eq.111 .and. kl1.eq.-3212)then
        if(isinel(48).eq.0)then
        si1=0.
        goto 1910
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)

1910    if(isinel(109).eq.0)then
        si2=0.
        goto 1666
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

1666    if(isinel(278).eq.0)then
        si3=0.
        goto 1667
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
1667    if(isinel(279).eq.0)then
        si4=0.
        goto 5711
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
5711    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3312
        ic=48
c       pion0 + sigma0bar to k- + cascade-bar
        goto 7125
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=-3322
        ic=109
c       pion0 + sigma0bar to k0- + cascade0bar
        goto 7125
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=321
        ik2=-2212
        ic=278
c      pion0 + sigma0bar to k+ + pbar
        goto 7125
        endif
        ik1=311
        ik2=-2112
        ic=279
c       pion0 + sigma0bar to k0 + nbar
7125    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion+ + sigma0
        if(kl.eq.211 .and. kl1.eq.3212)then
        if(isinel(98).eq.0)then
        si1=0.
        goto 2681
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3322),1)
        si1=s1724(ss,ilo1,0,the)
2681    if(isinel(251).eq.0)then
        si2=0.
        goto 2682
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(2212),1)
        ik3=-311
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
2682    if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3322
        ic=98
c       pion+ + sigma0 to k+ + cacade0
        goto 1683
        endif
        ik1=-311
        ik2=2212
        ic=251
c       pion+ + sigma0 to k0- + p
1683    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif
c       pion- + sigma0-
        if(kl.eq.-211 .and. kl1.eq.-3212)then
        if(isinel(106).eq.0)then
        si1=0.
        goto 2781
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3322),1)
        si1=s1724(ss,ilo1,0,the)
2781    if(isinel(273).eq.0)then
        si2=0.
        goto 2782
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(-2212),1)
        ik3=311
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
2782    if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3322
        ic=106
c       pion- + sigma0- to k- + cacade0-
        goto 1783
        endif
        ik1=311
        ik2=-2212
        ic=273
c       pion- + sigma0- to k0 + pba
1783    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion0 + sigma+
        if(kl.eq.111 .and. kl1.eq.3222)then
        if(isinel(101).eq.0)then
        si1=0.
        goto 3681
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3322),1)
        si1=s1724(ss,ilo1,0,the)
3681    if(isinel(259).eq.0)then
        si2=0.
        goto 3682
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(2212),1)
        ik3=-311
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
3682    if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3322
        ic=101
c       pion0 + sigma+ to k+ + cascade0
        goto 7683
        endif
        ik1=-311
        ik2=2212
        ic=259
c       pion0 + sigma+ to k0- + p
7683    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif
c       pion0 + sigma+bar
        if(kl.eq.111 .and. kl1.eq.-3222)then
        if(isinel(108).eq.0)then
        si1=0.
        goto 3781
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3322),1)
        si1=s1724(ss,ilo1,0,the)
3781    if(isinel(276).eq.0)then
        si2=0.
        goto 3782
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(-2212),1)
        ik3=311
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
3782    if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3322
        ic=108
c       pion0 + sigma+bar to k- + cascade0bar
        goto 6783
        endif
        ik1=311
        ik2=-2212
        ic=276
c       pion0 + sigma+bar to k0 + pbar
6783    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion0 + sigma-
        if(kl.eq.111 .and. kl1.eq.3112)then
        if(isinel(40).eq.0)then
        si1=0.
        goto 4681
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
4681    if(isinel(260).eq.0)then
        si2=0.
        goto 4682
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(2112),1)
        ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
4682    if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=311
        ik2=3312
        ic=40
c       pion0 + sigma- to k0+cascade-
        goto 2683
        endif
        ik1=-321
        ik2=2112
        ic=260
c       pion0 + sigma- to k- + n
2683    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif
c       pion0 + sigma-bar
        if(kl.eq.111 .and. kl1.eq.-3112)then
        if(isinel(47).eq.0)then
        si1=0.
        goto 7681
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
7681    if(isinel(277).eq.0)then
        si2=0.
        goto 7682
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(-2112),1)
        ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
7682    if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-311
        ik2=-3312
        ic=47
c       pion0 + sigma-bar to k0- +cascade-bar
        goto 7688
        endif
        ik1=321
        ik2=-2112
        ic=277
c       pion0 + sigma-bar to k+ + nbar
7688    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       k- + p
        if(kl.eq.-321 .and. kl1.eq.2212)then
        if(isinel(49).eq.0)then
        si1=0.
        goto 407
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3122),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of k + n to pion + y is assumed to be equal to 
c        ten times of pion + n to k + y
407     if(isinel(50).eq.0)then
        si2=0.
        goto 408
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3212),1)
        si2=s1724(ss,ilo2,0,the)
408     if(isinel(51).eq.0)then
        si3=0.
        goto 411
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3222),1)
        si3=s1724(ss,ilo3,0,the)
411     if(isinel(52).eq.0)then
        si4=0.
        goto 431
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3112),1)
        si4=s1724(ss,ilo4,0,the)
431     if(isinel(53).eq.0)then
        si5=0.
        goto 412
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si5=s1724(ss,ilo5,0,the)/10.
c       cross section of kaon + n to kaon + cascade is assumed to be 
c        equal to pion + n to kaon + y,but take the different of 
c        threshold energy into account

412     if(isinel(125).eq.0)then
        si6=0.
        goto 725
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si6=s1724(ss,ilo6,0,the)/10.

725     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0 
     c   .and. ilo5.eq.0 .and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6 .and. si5.eq.1.e-6 .and. si6.eq.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=3122
        ic=49
c       k- + p to pion0 + lambda
        goto 414
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=111
        ik2=3212
        ic=50
c       k- + p to pion0 + sigma0
        goto 414
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-211
        ik2=3222
        ic=51
c       k- + p to pion- + sigma+
        goto 414
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=211
        ik2=3112
        ic=52
c       k- + p to pion+ + sigma-
        goto 414
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=321
        ik2=3312
        ic=53
c       k- + p to k+ + cascade-
        goto 414
        endif
        ik1=311
        ik2=3322
        ic=125
c       k- + p to k0 + cascade0
414     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

c       k- + n
        if(kl.eq.-321 .and. kl1.eq.2112)then
        if(isinel(54).eq.0)then
        si1=0.
        goto 415
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3212),1)
        si1=s1724(ss,ilo1,0,the)
415     if(isinel(55).eq.0)then
        si2=0.
        goto 416
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3122),1)
        si2=s1724(ss,ilo2,0,the)
416     if(isinel(56).eq.0)then
        si3=0.
        goto 417
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3112),1)
        si3=s1724(ss,ilo3,0,the)
417     if(isinel(57).eq.0)then
        si4=0.
        goto 432
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3312),1)
        si4=s1724(ss,ilo4,0,the)/10.

432     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     $   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=3212
        ic=54
c       k- + n to pion- + sigma0
        goto 418
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-211
        ik2=3122
        ic=55
c       k- + n to pion- + lambda
        goto 418
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=111
        ik2=3112
        ic=56
c       k- + n to pion0 + sigma-
        goto 418
        endif
        ik1=311
        ik2=3312
        ic=57
c       k- + n to k0 + cascade-
418     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn
        goto 10
        endif

c       k0- + p
        if(kl.eq.-311 .and. kl1.eq.2212)then
        if(isinel(58).eq.0)then
        si1=0.
        goto 419
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3122),1)
        si1=s1724(ss,ilo1,0,the)
419     if(isinel(59).eq.0)then
        si2=0.
        goto 420
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3212),1)
        si2=s1724(ss,ilo2,0,the)
420     if(isinel(60).eq.0)then
        si3=0.
        goto 421
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3222),1)
        si3=s1724(ss,ilo3,0,the)
421     if(isinel(126).eq.0)then
        si4=0.
        goto 726
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3322),1)
        si4=s1724(ss,ilo4,0,the)/10.

726     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.si4.lt.
     &   1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        sit=si13+si4
        s1=si1/sit
        s2=si12/sit
        s3=si13/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3122
        ic=58
c       k0- + p to pion+ + lambda
        goto 422
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=211
        ik2=3212
        ic=59
c       k0- + p to pion+ + sigma0
        goto 422
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=111
        ik2=3222
        ic=60
c       k0- + p to pion0 + sigma+
        goto 422
        endif
        ik1=321
        ik2=3322
        ic=126
c       k0- + p to k+ + cascade0
422     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + n
        if(kl.eq.-311 .and. kl1.eq.2112)then
        if(isinel(61).eq.0)then
        si1=0.
        goto 423
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3112),1)
        si1=s1724(ss,ilo1,0,the)
423     if(isinel(62).eq.0)then
        si2=0.
        goto 424
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3222),1)
        si2=s1724(ss,ilo2,0,the)
424     if(isinel(63).eq.0)then
        si3=0.
        goto 425
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3212),1)
        si3=s1724(ss,ilo3,0,the)
425     if(isinel(64).eq.0)then
        si4=0.
        goto 426
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3122),1)
        si4=s1724(ss,ilo4,0,the)
426     if(isinel(65).eq.0)then
        si5=0.
        goto 433
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si5=s1724(ss,ilo5,0,the)/10.
433     if(isinel(127).eq.0)then
        si6=0.
        goto 729
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si6=s1724(ss,ilo6,0,the)/10.
729     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0 
     c   .and. ilo5.eq.0 .and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=3112
        ic=61
c       k0- + n to pion+ + sigma-
        goto 427
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=3222
        ic=62
c       k0- + n to pion- + sigma+
        goto 427
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=3212
        ic=63
c       k0- + n to pion0 + sigma0
        goto 427
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=111
        ik2=3122
        ic=64
c       k0- + n to pion0 + lambda
        goto 427
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=321
        ik2=3312
        ic=65
c       k0- + n to k+ + cascade-
        goto 427
        endif
        ik1=311
        ik2=3322
        ic=127
c       k0- + n to k0 + cascade0
427     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

c       k+ + p-
        if(kl.eq.321 .and. kl1.eq.-2212)then
        if(isinel(66).eq.0)then
        si1=0.
        goto 510
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of k + n- to pion + y- is assumed to be equal to 
c        ten times of pion + n to k + y
510     if(isinel(67).eq.0)then
        si2=0.
        goto 511
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
511     if(isinel(68).eq.0)then
        si3=0.
        goto 512
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
512     if(isinel(69).eq.0)then
        si4=0.
        goto 513
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3112),1)
        si4=s1724(ss,ilo4,0,the)
513     if(isinel(70).eq.0)then
        si5=0.
        goto 514
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si5=s1724(ss,ilo5,0,the)/10.

514     if(isinel(128).eq.0)then
        si6=0.
        goto 730
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)

        si6=s1724(ss,ilo6,0,the)/10.

730     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   .and. ilo5.eq.0 .and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6 .and. si5.eq.1.e-6.and.si6.eq.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=-3122
        ic=66
c       k+ + p- to pion0 + lambda-
        goto 515
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=111
        ik2=-3212
        ic=67
c       k+ + p- to pion0 + sigma0-
        goto 515
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=211
        ik2=-3222
        ic=68
c       k+ + p- to pion+ + sigma+ba
        goto 515
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=-211
        ik2=-3112
        ic=69
c       k+ + p- to pion- + sigma-ba
        goto 515
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=-321
        ik2=-3312
        ic=70
c       k+ + p- to k- + cascade-ba
        goto 515
        endif
        ik1=-311
        ik2=-3322
        ic=128
c       k+ + p- to k0- + cascade0ba
515     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

c       k+ + n-
        if(kl.eq.321 .and. kl1.eq.-2112)then
        if(isinel(71).eq.0)then
        si1=0.
        goto 516
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3212),1)
        si1=s1724(ss,ilo1,0,the)
516     if(isinel(72).eq.0)then
        si2=0.
        goto 517
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3122),1)
        si2=s1724(ss,ilo2,0,the)
517     if(isinel(73).eq.0)then
        si3=0.
        goto 518
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)
518     if(isinel(74).eq.0)then
        si4=0.
        goto 519
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3312),1)
        si4=s1724(ss,ilo4,0,the)/10.

519     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     $   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=-3212
        ic=71
c       k+ + n- to pion+ + sigma0-
        goto 520
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=211
        ik2=-3122
        ic=72
c       k+ + n- to pion+ + lambda-
        goto 520
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=111
        ik2=-3112
        ic=73
c       k+ + n- to pion0 + sigma-ba
        goto 520
        endif
        ik1=-311
        ik2=-3312
        ic=74
c       k+ + n- to k0- + cascade-ba
520     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn
        goto 10
        endif

c       k0 + p-
        if(kl.eq.311 .and. kl1.eq.-2212)then
        if(isinel(75).eq.0)then
        si1=0.
        goto 521
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
521     if(isinel(76).eq.0)then
        si2=0.
        goto 522
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
522     if(isinel(77).eq.0)then
        si3=0.
        goto 523
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
523     if(isinel(129).eq.0)then
        si4=0.
        goto 731
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3322),1)
        si4=s1724(ss,ilo4,0,the)/10.

731     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        sit=si13+si4
        s1=si1/sit
        s2=si12/sit
        s3=si13/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3122
        ic=75
c       k0 + p- to pion- + lambda-
        goto 524
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-211
        ik2=-3212
        ic=76
c       k0 + p- to pion- + sigma0-
        goto 524
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=111
        ik2=-3222
        ic=77
c       k0 + p- to pion0 + sigma+ba
        goto 524
        endif
        ik1=-321
        ik2=-3322
        ic=129
c       k0 + p- to k- + cascade0-
524     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0 + n-
        if(kl.eq.311 .and. kl1.eq.-2112)then
        if(isinel(78).eq.0)then
        si1=0.
        goto 525
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3112),1)
        si1=s1724(ss,ilo1,0,the)
525     if(isinel(79).eq.0)then
        si2=0.
        goto 526
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3222),1)
        si2=s1724(ss,ilo2,0,the)
526     if(isinel(80).eq.0)then
        si3=0.
        goto 527
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)
527     if(isinel(81).eq.0)then
        si4=0.
        goto 528
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3122),1)
        si4=s1724(ss,ilo4,0,the)
528     if(isinel(82).eq.0)then
        si5=0.
        goto 529
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si5=s1724(ss,ilo5,0,the)/10.

        if(isinel(130).eq.0)then
        si6=0.
        goto 529
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si6=s1724(ss,ilo6,0,the)/10.

529     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   .and. ilo5.eq.0.and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=-3112
        ic=78
c       k0 + n- to pion- + sigma-ba
        goto 530
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=211
        ik2=-3222
        ic=79
c       k0 + n- to pion+ + sigma+ba
        goto 530
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=-3212
        ic=80
c       k0 + n- to pion0 + sigma0-
        goto 530
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=111
        ik2=-3122
        ic=81
c       k0 + n- to pion0 + lambda-
        goto 530
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=-321
        ik2=-3312
        ic=82
c       k0 + n- to k- + cascade-ba
        goto 530
        endif
        ik1=-311
        ik2=-3322
        ic=130
c       k0 + n- to k0- + cascade0-
530     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

c       k- + lambda 
        if(kl.eq.-321 .and. kl1.eq.3122)then
        if(isinel(83).eq.0)then
        si1=0.
        goto 732
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of k + y to pion + cascade is assumed to be equal to 
c        k + n to pion + y
732     if(isinel(113).eq.0)then
        si2=0.
        goto 733
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
733     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=3312
        ic=83
c       k- + lambda to pion0 + cascade-
        goto 734
        endif
        ik1=-211
        ik2=3322
        ic=113
c       k- + lambda to pion- + cascade0
734     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k- + sigma+
        if(kl.eq.-321 .and. kl1.eq.3222)then
        if(isinel(84).eq.0)then
        si1=0.
        goto 735
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
735     if(isinel(111).eq.0)then
        si2=0.
        goto 736
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
736     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3312
        ic=84
c       k- + sigma+ to pion+ + cascade-
        goto 737
        endif
        ik1=111
        ik2=3322
        ic=111
c       k- + sigma+ to pion0 + cascade0
737     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k- + sigma- to pion- + cascade-
        if(kl.eq.-321 .and. kl1.eq.3112)then
        if(isinel(85).eq.0)goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(3312),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=3312
        lc(icp,5)=85
        goto 10
        endif

c       k- + sigma0 
        if(kl.eq.-321 .and. kl1.eq.3212)then
        if(isinel(86).eq.0)then
        si1=0.
        goto 738
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
738     if(isinel(112).eq.0)then
        si2=0.
        goto 739
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
739     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=3312
        ic=86
c       k- + sigma0 to pion0 + cascade-
        goto 740
        endif
        ik1=-211
        ik2=3322
        ic=112
740     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + lambda
        if(kl.eq.-311 .and. kl1.eq.3122)then
        if(isinel(87).eq.0)then
        si1=0.
        goto 741
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
741     if(isinel(117).eq.0)then
        si2=0.
        goto 742
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
742     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3312
        ic=87
c       k0- + lambda to pion+ + cascade-
        goto 743
        endif
        ik1=111
        ik2=3322
        ic=117
c       k0- + lambda to pion0 + cascade0
743     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + sigma0
        if(kl.eq.-311 .and. kl1.eq.3212)then
        if(isinel(88).eq.0)then
        si1=0.
        goto 744
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
744     if(isinel(116).eq.0)then
        si2=0.
        goto 745
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
745     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3312
        ic=88
c       k0- + sigma0 to pion+ + cascade-
        goto 746
        endif
        ik1=111
        ik2=3322
        ic=116
c       k0- + sigma0 to pion0 + cascade0
746     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + sigma- 
        if(kl.eq.-311 .and. kl1.eq.3112)then
        if(isinel(89).eq.0)then
        si1=0.
        goto 747
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
747     if(isinel(115).eq.0)then
        si2=0.
        goto 748
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
748     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=3312
        ic=89
c       k0- + sigma- to pion0 + cascade-
        goto 749
        endif
        ik1=-211
        ik2=3322
        ic=115
c       k0- + sigma- to pion- + cascade0
749     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k+ + lambda- 
        if(kl.eq.321 .and. kl1.eq.-3122)then
        if(isinel(90).eq.0)then
        si1=0.
        goto 750
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
750     if(isinel(120).eq.0)then
        si2=0.
        goto 751
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
751     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=-3312
        ic=90
c       k+ + lambda- to pion0 + cascade-ba
        goto 752
        endif
        ik1=211
        ik2=-3322
        ic=120
c       k+ + lambda- to pion+ + cascade0-
752     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k+ + sigma+ba 
        if(kl.eq.321 .and. kl1.eq.-3222)then
        if(isinel(91).eq.0)then
        si1=0.
        goto 753
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
753     if(isinel(118).eq.0)then
        si2=0.
        goto 754
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
754     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3312
        ic=91
c       k+ + sigma+ba to pion- + cascade-ba
        goto 755
        endif
        ik1=111
        ik2=-3322
        ic=118
c       k+ + sigma+ba to pion0 + cascade0-
755     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k+ + sigma-ba to pion+ + cascade-ba
        if(kl.eq.321 .and. kl1.eq.-3112)then
        if(isinel(92).eq.0)goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(-3312),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-3312
        lc(icp,5)=92
        goto 10
        endif

c       k+ + sigma0ba
        if(kl.eq.321 .and. kl1.eq.-3212)then
        if(isinel(93).eq.0)then
        si1=0.
        goto 756
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
756     if(isinel(119).eq.0)then
        si2=0.
        goto 757
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
757     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=-3312
        ic=93
c       k+ + sigma0- to pion0 + cascade-ba
        goto 758
        endif
        ik1=211
        ik2=-3322
        ic=119
c       k+ + sigma0- to pion+ + cascade0-
758     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0+ + lambda-
        if(kl.eq.311 .and. kl1.eq.-3122)then
        if(isinel(94).eq.0)then
        si1=0.
        goto 759
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
759     if(isinel(124).eq.0)then
        si2=0.
        goto 760
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
760     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3312
        ic=94
c       k0 + lambda- to pion- + cascade-ba
        goto 761
        endif
        ik1=111
        ik2=-3322
        ic=124
c       k0 + lambda- to pion0 + cascade0ba

761     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0 + sigma0-
        if(kl.eq.311 .and. kl1.eq.-3212)then
        if(isinel(95).eq.0)then
        si1=0.
        goto 762
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
762     if(isinel(123).eq.0)then
        si2=0.
        goto 763
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
763     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3312
        ic=95
c       k0 + sigma0- to pion- + cascade-ba
        goto 764
        endif
        ik1=111
        ik2=-3322
        ic=123
c       k0 + sigma0- to pion0 + cascade0-
764     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0 + sigma-ba
        if(kl.eq.311 .and. kl1.eq.-3112)then
        if(isinel(96).eq.0)then
        si1=0.
        goto 765
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
765     if(isinel(122).eq.0)then
        si2=0.
        goto 766
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
766     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=-3312
        ic=96
c       k0 + sigma-ba to pion0 + cascade-ba
        goto 767
        endif
        ik1=211
        ik2=-3322
        ic=122
c       k0 + sigma-ba to pion+ + cascade0ba
767     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + sigma+ to pion+ + cascade0
        if(kl.eq.-311 .and. kl1.eq.3222)then
        if(isinel(114).eq.0)goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(3322),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=3322
        lc(icp,5)=114
        goto 10
        endif

c       k0 + sigma+bar to pion- + cascade0ba
        if(kl.eq.311 .and. kl1.eq.-3222)then
        if(isinel(121).eq.0)goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3322),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-3322
        lc(icp,5)=121
        goto 10
        endif

c       k+ + cascade-ba to pion+ + omiga-ba
        if(kl.eq.321 .and. kl1.eq.-3312)then
        if(isinel(143).eq.0)goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-3334
        lc(icp,5)=143
        goto 10
        endif

c       k+ + cascade0- to pion0 + omiga-ba
        if(kl.eq.321 .and. kl1.eq.-3322)then
        if(isinel(145).eq.0)goto 13
        the=pmas(pycomp(111),1)+pmas(pycomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=111
        lc(icp,4)=-3334
        lc(icp,5)=145
        goto 10
        endif

c       k- + cascade- to pion- + omiga-
        if(kl.eq.-321 .and. kl1.eq.3312)then
        if(isinel(139).eq.0)goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=3334
        lc(icp,5)=139
        goto 10
        endif
       
c       k- + cascade0 to pion0 + omiga-
        if(kl.eq.-321 .and. kl1.eq.3322)then
        if(isinel(141).eq.0)goto 13
        the=pmas(pycomp(111),1)+pmas(pycomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=111
        lc(icp,4)=3334
        lc(icp,5)=141
        goto 10
        endif

c       k0 + cascade-ba to pion0 + omiga-ba
        if(kl.eq.311 .and. kl1.eq.-3312)then
        if(isinel(144).eq.0)goto 13
        the=pmas(pycomp(111),1)+pmas(pycomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=111
        lc(icp,4)=-3334
        lc(icp,5)=144
        goto 10
        endif

c       k0 + cascade0- to pion- + omiga-ba
        if(kl.eq.311 .and. kl1.eq.-3322)then
        if(isinel(146).eq.0)goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-3334
        lc(icp,5)=146
        goto 10
        endif

c       k0- + cascade- to pion0 + omiga-
        if(kl.eq.-311 .and. kl1.eq.3312)then
        if(isinel(140).eq.0)goto 13
        the=pmas(pycomp(111),1)+pmas(pycomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=111
        lc(icp,4)=3334
        lc(icp,5)=140
        goto 10
        endif

c       k0- + cascade0 to pion+ + omiga-
        if(kl.eq.-311 .and. kl1.eq.3322)then
        if(isinel(142).eq.0)goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=3334
        lc(icp,5)=142
        goto 10
        endif

c       follows are for reverse reactions
c       k+ + k-
        if((kl.eq.321 .and. kl1.eq.-321)
     c   .or.(kl.eq.-321 .and. kl1.eq.321))then
        if(isinel(201).eq.0)then
        fact1=0.
        goto 601
        endif
        ik3=211
        ik4=-211
        ic3=201
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
c       k+ + k- to pi+ + pi-
        fact1=fac
        if(ilo1.eq.0)fact1=0.
601     if(isinel(202).eq.0)then
        fact2=0.
        goto 602
        endif
        ik5=111
        ik6=111
        ic5=202
        call srev(kl,kl1,ik5,ik6,ss,ilo2,fac,
     &   0.5,0.5,0.,0.,1.,1.,0.,0.,0.5)
c       k+ + k- to pi0 + pi0
        fact2=fac
        if(ilo2.eq.0)fact2=0.
602     fact=fact1+fact2
        if(fact1.eq.0. .and. fact2.eq.0.)goto 13
c       if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(pyr(1).gt.fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*sig/cspipi
        goto 10
        endif

c       k+ + k0-
        if((kl.eq.321 .and. kl1.eq.-311)
     c   .or.(kl.eq.-311 .and. kl1.eq.321))then
        if(isinel(203).eq.0)goto 13
        ik3=211
        ik4=111
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
c       k+ + k0- to pi+ + pi0
        if(ilo1.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=203
        tw(icp)=fac*sig/cspipi
        goto 10
        endif

c       k- + k0
        if((kl.eq.-321 .and. kl1.eq.311)
     c   .or.(kl.eq.311 .and. kl1.eq.-321))then
        if(isinel(204).eq.0)goto 13
        ik3=-211
        ik4=111
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
c       k- + k0 to pi- + pi0
        if(ilo1.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=204
        tw(icp)=fac*sig/cspipi
        goto 10
        endif 

c       k0 + k0-
        if((kl.eq.311 .and. kl1.eq.-311)
     c   .or.(kl.eq.-311 .and. kl1.eq.311))then
        if(isinel(205).eq.0)then
        fact1=0.
        goto 603
        endif
        ik3=211
        ik4=-211
        ic3=205
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
c       k0 + k0- to pi- + pi+
        fact1=fac
        if(ilo1.eq.0)fact1=0.
603     if(isinel(206).eq.0)then
        fact2=0.
        goto 604
        endif
        ik5=111
        ik6=111
        ic5=206
        call srev(kl,kl1,ik5,ik6,ss,ilo2,fac,
     &   0.5,0.5,0.,0.,1.,1.,0.,0.,0.5)
c       k0 + k0- to pi0 + pi0
        fact2=fac
        if(ilo2.eq.0)fact2=0.
604     fact=fact1+fact2
        if(fact1.eq.0. .and. fact2.eq.0.)goto 13
c       if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(pyr(1).gt.fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*sig/cspipi
        goto 10
        endif


c       k+ + sigma+ to pion+ + p
        if(kl.eq.321 .and. kl1.eq.3222)then
c       3222 is the flavor code of sigma+
        if(isinel(207).eq.0)goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(2212),1)
        ww=s1713(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=207
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif


c       k+ + sigma-
        if(kl.eq.321 .and. kl1.eq.3112)then
        if(isinel(208).eq.0)then
        si1=0.
        goto 606
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(2212),1)
        si1=s0715(ss,ilo1,0,the)*fac
606     if(isinel(209).eq.0)then
        si2=0.
        goto 607
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2112),1)
        si2=s2325(ss,ilo2,0,the)*fac
607     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=2212
        ic=208
c       k+ + sigma- to pion- + p
        goto 609
        endif
        ik1=111
        ik2=2112
        ic=209
c       k+ + sigma- to pion0 + n
609     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif



c       k+ + sigma0
        if(kl.eq.321 .and. kl1.eq.3212)then
        if(isinel(210).eq.0)then
        si1=0.
        goto 610
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
610     if(isinel(211).eq.0)then
        si2=0.
        goto 611
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2212),1)
        si2=s2314(ss,ilo2,0,the)*fac
611     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=2112
        ic=210
c       k+ + sigma0 to pion+ + n
        goto 612
        endif
        ik1=111
        ik2=2212
        ic=211
c       k+ + sigma0 to pion0 + p
612     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif


c       k+ + lambda0
        if(kl.eq.321 .and. kl1.eq.3122)then
        if(isinel(212).eq.0)then
        si1=0.
        goto 613
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(2112),1)
        si1=s1727(ss,ilo1,0,the)*fac
613     if(isinel(213).eq.0)then
        si2=0.
        goto 614
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2212),1)
        si2=s2317(ss,ilo2,0,the)*fac
614     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=2112
        ic=212
c       k+ + lambda0 to pion+ + n
        goto 615
        endif
        ik1=111
        ik2=2212
        ic=213
c       k+ + lambda0 to pion0 + p
615     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif


c       k0 + sigma+
        if(kl.eq.311 .and. kl1.eq.3222)then
        if(isinel(214).eq.0)then
        si1=0.
        goto 616
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
616     if(isinel(215).eq.0)then
        si2=0.
        goto 617
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
617     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=2112
        ic=214
c       k+ + sigma+ to pion+ + n
        goto 618
        endif
        ik1=111
        ik2=2212
        ic=215
c       k+ + sigma+ to pion0 + p
618     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0 + sigma- to pion- +n 
        if(kl.eq.311 .and. kl1.eq.3112)then
c       3112 is the flavor code of sigma-
        if(isinel(216).eq.0)goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=216
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       k0 + sigma0
        if(kl.eq.311 .and. kl1.eq.3212)then
        if(isinel(217).eq.0)then
        si1=0.
        goto 619
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(2212),1)
        si1=s07123(ss,ilo1,0,the)*fac
619     if(isinel(218).eq.0)then
        si2=0.
        goto 620
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
620     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=2212
        ic=217
c       k0 + sigma0 to pion- + p
        goto 621
        endif
        ik1=111
        ik2=2112
        ic=218
c       k0 + sigma0 to pion0 + n
621     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0 + lambda0
        if(kl.eq.311 .and. kl1.eq.3122)then
        if(isinel(219).eq.0)then
        si1=0.
        goto 622
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(2212),1)
        si1=s07122(ss,ilo1,0,the)*fac
622     if(isinel(220).eq.0)then
        si2=0.
        goto 623
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
623     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=2212
        ic=219
c       k0 + lambda0 to pion- + p
        goto 624
        endif
        ik1=111
        ik2=2112
        ic=220
c       k0 + lambda0 to pion0 + n
624     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k- + sigma+ba to pion- + pba
        if(kl.eq.-321 .and. kl1.eq.-3222)then
        if(isinel(221).eq.0)goto 13
        the=pmas(pycomp(-211),1)+pmas(pycomp(-2212),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-2212
        lc(icp,5)=221
        ik3=-211
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       k- + sigma-ba
        if(kl.eq.-321 .and. kl1.eq.-3112)then
        if(isinel(222).eq.0)then
        si1=0.
        goto 625
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
625     if(isinel(223).eq.0)then
        si2=0.
        goto 626
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
626     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=-2212
        ic=222
c       k- + sigma-ba to pion+ + pba
        goto 627
        endif
        ik1=111
        ik2=-2112
        ic=223
c       k- + sigma-ba to pion0 + nba
627     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif


c       k- + sigma0ba
        if(kl.eq.-321 .and. kl1.eq.-3212)then
        if(isinel(224).eq.0)then
        si1=0.
        goto 628
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
628     if(isinel(225).eq.0)then
        si2=0.
        goto 629
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
629     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-2112
        ic=224
c       k- + sigma0ba to pion- + nba
        goto 630
        endif
        ik1=111
        ik2=-2212
        ic=225
c       k- + sigma0ba to pion0 + pba
630     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k- + lambda0-
        if(kl.eq.-321 .and. kl1.eq.-3122)then
        if(isinel(226).eq.0)then
        si1=0.
        goto 631
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
631     if(isinel(227).eq.0)then
        si2=0.
        goto 632
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
632     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-2112
        ic=226
c       k- + lambda0ba to pion- + nba
        goto 633
        endif
        ik1=111
        ik2=-2212
        ic=227
c       k- + lambda0ba to pion0 + pba
633     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0- + sigma+ba
        if(kl.eq.-311 .and. kl1.eq.-3222)then
        if(isinel(228).eq.0)then
        si1=0.
        goto 634
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
634     if(isinel(229).eq.0)then
        si2=0.
        goto 635
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
635     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-2112
        ic=228
c       k0ba + sigma+ba to pion- + nba
        goto 636
        endif
        ik1=111
        ik2=-2212
        ic=229
c       k0ba + sigma+ba to pion0 + pba
636     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0- + sigma-ba to pion+ + nba
        if(kl.eq.-311 .and. kl1.eq.-3112)then
        if(isinel(230).eq.0)goto 13
        the=pmas(pycomp(211),1)+pmas(pycomp(-2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-2112
        lc(icp,5)=230
        ik3=211
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       k0- + sigma0-
        if(kl.eq.-311 .and. kl1.eq.-3212)then
        if(isinel(231).eq.0)then
        si1=0.
        goto 637
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
637     if(isinel(232).eq.0)then
        si2=0.
        goto 638
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
638     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=-2212
        ic=231
c       k0- + sigma0- to pion+ + pba
        goto 639
        endif
        ik1=111
        ik2=-2112
        ic=232
c       k0- + sigma0- to pion0 + nba
639     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0- + lambda0-
        if(kl.eq.-311 .and. kl1.eq.-3122)then
        if(isinel(233).eq.0)then
        si1=0.
        goto 640
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
640     if(isinel(234).eq.0)then
        si2=0.
        goto 641
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
641     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=-2212
        ic=233
c       k0- + lambda0- to pion+ + pba
        goto 642
        endif
        ik1=111
        ik2=-2112
        ic=234
c       k0- + lambda0- to pion0 + nba
642     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k+ + cascade-
        if(kl.eq.321 .and. kl1.eq.3312)then
        if(isinel(235).eq.0)then
        si1=0.
        goto 699
        endif
        ik1=211
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3112),1)
        si1=s1724(ss,ilo1,0,the)*fac
699     if(isinel(236).eq.0)then
        si2=0.
        goto 643
        endif
        ik1=-211
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3222),1)
        si2=s1724(ss,ilo2,0,the)*fac
643     if(isinel(237).eq.0)then
        si3=0.
        goto 644
        endif
        ik1=111
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3122),1)
        si3=s1724(ss,ilo3,0,the)*fac
644     if(isinel(238).eq.0)then
        si4=0.
        goto 6430
        endif
        ik1=111
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3212),1)
        si4=s1724(ss,ilo4,0,the)*fac
6430    if(isinel(253).eq.0)then
        si5=0.
        goto 6440
        endif
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
6440    if(isinel(265).eq.0)then
        si6=0.
        goto 645
        endif
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
645     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   .and. ilo5.eq.0 .and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6 .and.
     c   si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=3112
        ic=235
c       k+ + cascade- to pion+ + sigma-
        goto 646
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=3222
        ic=236
c       k+ + cascade- to pion- + sigma+
        goto 646
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=3122
        ic=237
c       k+ + cascade- to pion0 + lambda0
        goto 646
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=111
        ik2=3212
        ic=238
c       k+ + cascade- to pion0 + sigma0
        goto 646
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=-321
        ik2=2212
        ic=253
c       k+ + cascade- to k- + p
        goto 646
        endif
        ik1=-311
        ik2=2112
        ic=265
c       k+ + cascade- to k0- + n
646     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10.
        goto 10
        endif




c       k- + cascade-ba
        if(kl.eq.-321 .and. kl1.eq.-3312)then
        if(isinel(239).eq.0)then
        si1=0.
        goto 647
        endif
        ik1=211
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3222),1)
        si1=s1724(ss,ilo1,0,the)*fac
647     if(isinel(240).eq.0)then
        si2=0.
        goto 648
        endif
        ik1=-211
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3112),1)
        si2=s1724(ss,ilo2,0,the)*fac
648     if(isinel(241).eq.0)then
        si3=0.
        goto 649
        endif
        ik1=111
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3122),1)
        si3=s1724(ss,ilo3,0,the)*fac
649     if(isinel(242).eq.0)then
        si4=0.
        goto 6500
        endif
        ik1=111
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3212),1)
        si4=s1724(ss,ilo4,0,the)*fac
6500    if(isinel(270).eq.0)then
        si5=0.
        goto 6501
        endif
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
6501    if(isinel(282).eq.0)then
        si6=0.
        goto 650
        endif
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
650     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   .and. ilo5.eq.0 .and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     c   si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=-3222
        ic=239
c       k- + cascade-ba to pion+ + sigma+ba
        goto 651
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=-3112
        ic=240
c       k- + cascade-ba to pion- + sigma-ba
        goto 651
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=-3122
        ic=241
c       k- + cascade-ba to pion0 + lambda0-
        goto 651
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=111
        ik2=-3212
        ic=242
c       k- + cascade-ba to pion0 + sigma0-
        goto 651
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=321
        ik2=-2212
        ic=270
c       k- + cascade-ba to k+ + pba
        goto 651
        endif
        ik1=311
        ik2=-2112
        ic=282
c       k- + cascade-ba to k0 +nba
651     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10.
        goto 10
        endif




c       k0 + cascade-
        if(kl.eq.311 .and. kl1.eq.3312)then
        if(isinel(243).eq.0)then
        si1=0.
        goto 652
        endif
        ik1=-211
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3122),1)
        si1=s1724(ss,ilo1,0,the)*fac
652     if(isinel(244).eq.0)then
        si2=0.
        goto 653
        endif
        ik1=-211
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3212),1)
        si2=s1724(ss,ilo2,0,the)*fac
653     if(isinel(245).eq.0)then
        si3=0.
        goto 6540
        endif
        ik1=111
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3112),1)
        si3=s1724(ss,ilo3,0,the)*fac
6540    if(isinel(257).eq.0)then
        si4=0.
        goto 654
        endif
        ik1=-321
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(2112),1)
        si4=s1724(ss,ilo4,0,the)*fac
654     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &   .and. ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=3122
        ic=243
c       k0 + cascade- to pion- + lambda0
        goto 655
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=3212
        ic=244
c       k0 + cascade- to pion- + sigma0
        goto 655
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=3112
        ic=245
c       k0 + cascade- to pion0 + sigma-
        goto 655
        endif
        ik1=-321
        ik2=2112
        ic=257
c       k0 + cascade- to k- + n
655     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10.
        goto 10
        endif



c       k0- + cascade-ba
        if(kl.eq.-311 .and. kl1.eq.-3312)then
        if(isinel(246).eq.0)then
        si1=0.
        goto 656
        endif
        ik1=211
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)*fac
656     if(isinel(247).eq.0)then
        si2=0.
        goto 657
        endif
        ik1=211
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)*fac
657     if(isinel(248).eq.0)then
        si3=0.
        goto 6580
        endif
        ik1=111
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)*fac
6580    if(isinel(274).eq.0)then
        si4=0.
        goto 658
        endif
        ik1=321
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(2112),1)
        si4=s1724(ss,ilo4,0,the)*fac
658     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &   .and. ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=-3122
        ic=246
c       k0- + cascade-ba to pion+ + lambda0-
        goto 659
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=211
        ik2=-3212
        ic=247
c       k0- + cascade-ba to pion+ + sigma0-
        goto 659
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=-3112
        ic=248
c       k0- + cascade-ba to pion0 + sigma-ba
        goto 659
        endif
        ik1=321
        ik2=-2112
        ic=274
c       k0- + cascade-ba to k+ + nba
659     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10.
        goto 10
        endif




c       pion+ + lambda0 
        if(kl.eq.211 .and. kl1.eq.3122)then
        if(isinel(99).eq.0)then
        sigma1=0.
        goto 2522
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3322),1)
        sigma1=s1724(ss,ilo1,0,the)
c       pion+ + lambda0 to k+ + cascade
2522    if(isinel(252).eq.0)then
        sigma2=0.
        goto 2523
        endif
        ik3=-311
         ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(2212),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10.

c       cross section of     pion+ + lambda0 to k0- + p
2523    if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6) goto 13
        ik1=321
        ik2=3322
        ic=99
c       pion+ + lambda0 to k+ + cascade0
        sigm12=sigma1+sigma2
        if(pyr(1).gt.sigma1/sigm12)then
        ik1=-311
        ik2=2212
        ic=252
c       cross section of     pion+ + lambda0 to k0- + p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10.
        goto 10
        endif 

c       pion- + lambda0-
        if(kl.eq.-211 .and. kl1.eq.-3122)then
        if(isinel(107).eq.0)then
        sigma1=0.
        goto 3522
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3322),1)
        sigma1=s1724(ss,ilo1,0,the)
c       pion- + lambda- to k- + cascade0-
3522        if(isinel(275).eq.0)then
        sigma2=0.
        goto 3523
        endif
        ik3=311
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-2212),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10.
c       pion- + lambda0- to k0 + pba

3523    if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6) goto 13
        ik1=-321
        ik2=-3322
        ic=107
c       pion- + lambda- to k- + cascade0-
        sigm12=sigma1+sigma2
        if(pyr(1).gt.sigma1/sigm12)then
        ik1=311
        ik2=-2212
        ic=275
c       pion- + lambda0- to k0 + pba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10.
        goto 10
        endif 



c       pion- + lambda0 
        if(kl.eq.-211 .and. kl1.eq.3122)then
        if(isinel(36).eq.0)then
        sigma1=0.
        goto 1522
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3312),1)
        sigma1=s1724(ss,ilo1,0,the)
c       pion- + lambda to k0 + cascade-
1522    if(isinel(258).eq.0)then
        sigma2=0.
        goto 1523
        endif
        ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(2112),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10.

c       cross section of   pion- + lambda0 to k- + n
1523    if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)  goto 13
        ik1=311
        ik2=3312
        ic=36
c       pion- + lambda to k0 + cascade-
        sigm12=sigma1+sigma2
        if(pyr(1).gt.sigma1/sigm12)then
        ik1=-321
        ik2=2112
        ic=258
c       cross section of pion- + lambda0 to k- + n
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10.
        goto 10
        endif 
c       pion+ + lambda-
        if(kl.eq.211 .and. kl1.eq.-3122)then
        if(isinel(42).eq.0)then
        sigma1=0.
        goto 6522
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3312),1)
        sigma1=s1724(ss,ilo1,0,the)
c       pion+ + lambdaba to k0- + cascade-ba
6522    if(isinel(269).eq.0)then
        sigma2=0.
        goto 6523
        endif
        ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &   1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-2112),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10.

c       pion+ + lambda0- to k+ + nba
6523    if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6) goto 13
        ik1=-311
        ik2=-3312
        ic=42
c       pion+ + lambdaba to k0- + cascade-ba
        sigm12=sigma1+sigma2
        if(pyr(1).gt.sigma1/sigm12)then
        ik1=321
        ik2=-2112
        ic=269
c       pion+ + lambda0- to k+ + nba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10.
        goto 10
        endif 


c       pion0 + lambda0
        if(kl.eq.111 .and. kl1.eq.3122)then
        if(isinel(263).eq.0)then
        si1=0.
        goto 669
        endif
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
669     if(isinel(264).eq.0)then
        si2=0.
        goto 670
        endif
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
670      if(isinel(39).eq.0)then
        si3=0.
        goto 1670
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3312),1)
        si3=s1724(ss,ilo3,0,the)
1670      if(isinel(103).eq.0)then
        si4=0.
        goto 1671
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3322),1)
        si4=s1724(ss,ilo4,0,the)
1671    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=2212
        ic=263
c       pion0 + lambda0 to k- + p
        goto 671
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=2112
        ic=264
c       pion0 + lambda0 to k0- + n
        goto 671
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
      ik1=321
        ik2=3312
        ic=39
c       pion0 + lambda0 to  k+ + cascade-
        goto 671
        endif
      ik1=311
        ik2=3322
        ic=103
c       pion0 + lambda to k0 + cascade0
671     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif
c       pion0 + lambda0-
        if(kl.eq.111 .and. kl1.eq.-3122)then
        if(isinel(280).eq.0)then
        si1=0.
        goto 2669
        endif
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
2669    if(isinel(281).eq.0)then
        si2=0.
        goto 2677
        endif
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
2677      if(isinel(46).eq.0)then
        si3=0.
        goto 2670
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3312),1)
        si3=s1724(ss,ilo3,0,the)
c       pion0 + lambda- to k- + cascade-ba
2670      if(isinel(110).eq.0)then
        si4=0.
        goto 2671
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3322),1)
        si4=s1724(ss,ilo4,0,the)
2671    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     &   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=-2212
        ic=280
c       pion0 + lambda0- to k+ + pbar
        goto 6711
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=311
        ik2=-2112
        ic=281
c       pion0 + lambda0- to k0 + nba
        goto 6711
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
      ik1=-321
        ik2=-3312
        ic=46
c       pion0 + lambda- to k- + cascade-ba
        goto 6711
        endif
      ik1=-311
        ik2=-3322
        ic=110
c          pion0 + lambda- to k0- + cascade0-
6711    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif



c       pion+ + cascade-
        if(kl.eq.211 .and. kl1.eq.3312)then
        if(isinel(283).eq.0)then
        si1=0.
        goto 684
        endif
        ik1=-321
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3222),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
684     if(isinel(284).eq.0)then
        si2=0.
        goto 685
        endif
        ik1=-311
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3122),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
685     if(isinel(285).eq.0)then
        si3=0.
        goto 686
        endif
        ik1=-311
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3212),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
686     if(isinel(131).eq.0)then
        si4=0.
        goto 1686
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3334),1)
        si4=s1724(ss,ilo4,0,the)

1686    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &   .and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-321
        ik2=3222
        ic=283
c       pion+ + cascade- to k- + sigma+
        goto 687
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-311
        ik2=3122
        ic=284
c       pion+ + cascade- to k0- + lambda0
        goto 687
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-311
        ik2=3212
        ic=285
c       pion+ + cascade- to k0- + sigma0
        goto 687
        endif
       ik1=321
        ik2=3334
        ic=131
c       pion+ + cascade- to k+ + omiga-
687     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + cascade- to k- + sigma-
        if(kl.eq.-211 .and. kl1.eq.3312)then
        if(isinel(286).eq.0)goto 13
        the=pmas(pycomp(-321),1)+pmas(pycomp(3112),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=3112
        lc(icp,5)=286
        ik3=-321
        ik4=3112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif 

c       pion0 + cascade-
        if(kl.eq.111 .and. kl1.eq.3312)then
        if(isinel(287).eq.0)then
        si1=0.
        goto 688
        endif
        ik1=-321
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3122),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
688     if(isinel(288).eq.0)then
        si2=0.
        goto 689
        endif
        ik1=-321
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3212),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
689     if(isinel(289).eq.0)then
        si3=0.
        goto 690
        endif
        ik1=-311
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3112),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
690     if(isinel(132).eq.0)then
        si4=0.
        goto 1690
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3334),1)
        si4=s1724(ss,ilo4,0,the)

c       pion0 + cascade- to k0 + omiga-
1690    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &   ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-321
        ik2=3122
        ic=287
c       pion0 + cascade- to k- + lambda0
        goto 691
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-321
        ik2=3212
        ic=288
c       pion0 + cascade- to k- + sigma0
        goto 691
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-311
        ik2=3112
        ic=289
c       pion0 + cascade- to k0- + sigma-
        goto 691
        endif
        ik1=311
        ik2=3334
        ic=132
c       pion0 + cascade- to k0 + omiga-
691     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion+ + cascade-ba to k+ + sigma-ba
        if(kl.eq.211 .and. kl1.eq.-3312)then
        if(isinel(290).eq.0)goto 13
        the=pmas(pycomp(321),1)+pmas(pycomp(-3112),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=321
        lc(icp,4)=-3112
        lc(icp,5)=290
        ik3=321
        ik4=-3112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)

        tw(icp)=ww*fac
        goto 10
        endif 

c       pion- + cascade-ba
        if(kl.eq.-211 .and. kl1.eq.-3312)then
        if(isinel(291).eq.0)then
        si1=0.
        goto 692
        endif
        ik1=321
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3222),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
692     if(isinel(292).eq.0)then
        si2=0.
        goto 693
        endif
        ik1=311
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3122),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
693     if(isinel(293).eq.0)then
        si3=0.
        goto 694
        endif
        ik1=311
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3212),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
694     if(isinel(133).eq.0)then
        si4=0.
        goto 1694
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3334),1)
        si4=s1724(ss,ilo4,0,the)
1694    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &   ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=321
        ik2=-3222
        ic=291
c       pion- + cascade-ba to k+ + sigma+-
        goto 695
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=311
        ik2=-3122
        ic=292
c       pion- + cascade-ba to k0 + lambda0-
        goto 695
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=311
        ik2=-3212
        ic=293
c       pion- + cascade-ba to k0 + sigma0-
        goto 695
        endif
        ik1=-321
        ik2=-3334
        ic=133
c       pion- + cascade-ba to k- + omiga-ba
695     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion0 + cascade-ba
        if(kl.eq.111 .and. kl1.eq.-3312)then
        if(isinel(294).eq.0)then
        si1=0.
        goto 696
        endif
        ik1=321
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3122),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
696     if(isinel(295).eq.0)then
        si2=0.
        goto 697
        endif
        ik1=321
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3212),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
697     if(isinel(296).eq.0)then
        si3=0.
        goto 698
        endif
        ik1=311
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3112),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
698     if(isinel(134).eq.0)then
        si4=0.
        goto 1698
        endif
         the=pmas(pycomp(-311),1)+pmas(pycomp(-3334),1)
        si4=s1724(ss,ilo4,0,the)
1698    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &   ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=321
        ik2=-3122
        ic=294
c       pion0 + cascade-ba to k+ + lambda0-
        goto 700
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=321
        ik2=-3212
        ic=295
c       pion0 + cascade-ba to k+ + sigma0-
        goto 700
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=311
        ik2=-3112
        ic=296
c       pion0 + cascade-ba to k0 + sigma-ba
        goto 700
        endif
       ik1=-311
        ik2=-3334
        ic=134
c       pion0 + cascade-ba to k0- + omiga-ba
700     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       k+ + cascade0
        if(kl.eq.321 .and. kl1.eq.3322)then
        if(isinel(298).eq.0)then
        si1=0.
        goto 768
        endif
        ik1=211
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3212),1)
        si1=s1724(ss,ilo1,0,the)*fac
768     if(isinel(299).eq.0)then
        si2=0.
        goto 769
        endif
        ik1=211
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3122),1)
        si2=s1724(ss,ilo2,0,the)*fac
769     if(isinel(301).eq.0)then
        si3=0.
        goto 770
        endif
        ik1=111
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3222),1)
        si3=s1724(ss,ilo3,0,the)*fac
770     if(isinel(326).eq.0)then
        si4=0.
        goto 771
        endif
        ik1=-311
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(2212),1)
        si4=s1724(ss,ilo4,0,the)*fac
771     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   )goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=3212
        ic=298
c       k+ + cascade0 to pion+ + sigma0
        goto 772
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=211
        ik2=3122
        ic=299
c       k+ + cascade0 to pion+ + lambda
        goto 772
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=3222
        ic=301
c       k+ + cascade0 to pion0 + sigma+
        goto 772
        endif
        ik1=-311
        ik2=2212
        ic=326
c       k+ + cascade0 to k0- + p
772     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10.
        goto 10
        endif

c       k- + cascad0-
        if(kl.eq.-321 .and. kl1.eq.-3322)then
        if(isinel(306).eq.0)then
        si1=0.
        goto 773
        endif
        ik1=-211
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3212),1)
        si1=s1724(ss,ilo1,0,the)*fac
773     if(isinel(307).eq.0)then
        si2=0.
        goto 774
        endif
        ik1=-211
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3122),1)
        si2=s1724(ss,ilo2,0,the)*fac
774     if(isinel(308).eq.0)then
        si3=0.
        goto 775
        endif
        ik1=111
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)*fac
775     if(isinel(329).eq.0)then
        si4=0.
        goto 776
        endif
        ik1=311
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-2212),1)
        si4=s1724(ss,ilo4,0,the)*fac
776     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   )goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=-3212
        ic=306
c       k- + cascade0- to pion- + sigma0-
        goto 777
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=-3122
        ic=307
c       k- + cascade0- to pion- + lambda-
        goto 777
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=-3222
        ic=308
c       k- + cascade0- to pion0 + sigma+ba
        goto 777
        endif
        ik1=311
        ik2=-2212
        ic=329
c       k- + cascade0- to k0 + p-
777     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10.
        goto 10
        endif

c       k0 + cascad0
        if(kl.eq.311 .and. kl1.eq.3322)then
        if(isinel(297).eq.0)then
        si1=0.
        goto 778
        endif
        ik1=211
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3112),1)
        si1=s1724(ss,ilo1,0,the)*fac
778     if(isinel(300).eq.0)then
        si2=0.
        goto 779
        endif
        ik1=-211
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3222),1)
        si2=s1724(ss,ilo2,0,the)*fac
779     if(isinel(302).eq.0)then
        si3=0.
        goto 780
        endif
        ik1=111
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3212),1)
        si3=s1724(ss,ilo3,0,the)*fac
780     if(isinel(303).eq.0)then
        si4=0.
        goto 781
        endif
        ik1=111
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3122),1)
        si4=s1724(ss,ilo4,0,the)*fac
781     if(isinel(325).eq.0)then
        si5=0.
        goto 782
        endif
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac,
     &   0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
782     if(isinel(327).eq.0)then
        si6=0.
        goto 783
        endif
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac,
     &   0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
783     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   .and. ilo5.eq.0.and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=3112
        ic=297
c       k0 + cascade0 to pion+ + sigma-
        goto 784
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=3222
        ic=300
c       k0 + cascade0 to pion- + sigma+
        goto 784
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=3212
        ic=302
c       k0 + cascade0 to pion0 + sigma0
        goto 784
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=111
        ik2=3122
        ic=303
c       k0 + cascade0 to pion0 + lambda
        goto 784
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=-321
        ik2=2212
        ic=325
c       k0 + cascade0 to k- + p
        goto 784
        endif
        ik1=-311
        ik2=2112
        ic=327
c       k0 + cascade0 to k0- + n
784     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10.
        goto 10
        endif

c       k0- + cascad0-
        if(kl.eq.-311 .and. kl1.eq.-3322)then
        if(isinel(304).eq.0)then
        si1=0.
        goto 785
        endif
        ik1=211
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3222),1)
        si1=s1724(ss,ilo1,0,the)*fac
785     if(isinel(305).eq.0)then
        si2=0.
        goto 786
        endif
        ik1=-211
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3112),1)
        si2=s1724(ss,ilo2,0,the)*fac
786     if(isinel(309).eq.0)then
        si3=0.
        goto 787
        endif
        ik1=111
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)*fac
787     if(isinel(310).eq.0)then
        si4=0.
        goto 788
        endif
        ik1=111
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &   0.5,0.5,0.,0.5,1.0,0.,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3122),1)
        si4=s1724(ss,ilo4,0,the)*fac
788     if(isinel(328).eq.0)then
        si5=0.
        goto 789
        endif
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac,
     &   0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
789     if(isinel(330).eq.0)then
        si6=0.
        goto 790
        endif
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac,
     &   0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
790     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   .and. ilo5.eq.0.and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=-3222
        ic=304
c       k-0 + cascade0- to pion+ + sigma-ba
        goto 791
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=-3112
        ic=305
c       k0- + cascade0- to pion- + sigma-ba
        goto 791
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=-3212
        ic=309
c       k0- + cascade0- to pion0 + sigma0-
        goto 791
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=111
        ik2=-3122
        ic=310
c       k0- + cascade0- to pion0 + lambda-
        goto 791
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=321
        ik2=-2212
        ic=328
c       k0- + cascade0- to k+ + p-
        goto 791
        endif
        ik1=311
        ik2=-2112
        ic=330
c       k0- + cascade0- to k0 + n-
791     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10.
        goto 10
        endif

c       k+ + omiga-
        if(kl.eq.321 .and. kl1.eq.3334)then
        if(isinel(331).eq.0)then
        si1=0.
        goto 792
        endif
        ik1=211
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
792     if(isinel(336).eq.0)then
        si2=0.
        goto 793
        endif
        ik1=111
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
793     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3312
        ic=331
c       k+ + omiga- to pion+ + cascade-
        goto 794
        endif
        ik1=111
        ik2=3322
        ic=336
c       k+ + omiga- to pion0 + cascade0
794     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k- + omiga-ba
        if(kl.eq.-321 .and. kl1.eq.-3334)then
        if(isinel(333).eq.0)then
        si1=0.
        goto 795
        endif
        ik1=-211
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
795     if(isinel(338).eq.0)then
        si2=0.
        goto 796
        endif
        ik1=111
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
796     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3312
        ic=333
c       k- + omiga-ba to pion- + cascade-ba
        goto 797
        endif
        ik1=111
        ik2=-3322
        ic=338
c       k- + omiga-ba to pion0 + cascade0-
797     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0 + omiga-
        if(kl.eq.311 .and. kl1.eq.3334)then
        if(isinel(332).eq.0)then
        si1=0.
        goto 798
        endif
        ik1=111
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
798     if(isinel(335).eq.0)then
        si2=0.
        goto 799
        endif
        ik1=-211
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(-211),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
799     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=3312
        ic=332
c       k0 + omiga- to pion0 + cascade-
        goto 800
        endif
        ik1=-211
        ik2=3322
        ic=335
c       k0 + omiga- to pion- + cascade0
800     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0- + omiga-ba
        if(kl.eq.-311 .and. kl1.eq.-3334)then
        if(isinel(334).eq.0)then
        si1=0.
        goto 801
        endif
        ik1=111
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(111),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
801     if(isinel(337).eq.0)then
        si2=0.
        goto 802
        endif
        ik1=211
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(pycomp(211),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
802     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=-3312
        ic=334
c       k0- + omiga-ba to pion0 + cascade-ba
        goto 803
        endif
        ik1=211
        ik2=-3322
        ic=337
c       k0- + omiga-ba to pion+ + cascade0-
803     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       pion+ + cascade0 to k0- + sigma+
        if(kl.eq.211 .and. kl1.eq.3322)then
        if(isinel(314).eq.0)goto 13
        the=pmas(pycomp(-311),1)+pmas(pycomp(3222),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=-311
        lc(icp,4)=3222
        lc(icp,5)=314
        ik3=-311
        ik4=3222
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion+ + cascade0- 
        if(kl.eq.211 .and. kl1.eq.-3322)then
        if(isinel(319).eq.0)then
        si1=0.
        goto 804
        endif
        ik1=321
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3212),1)
        si1=10.*s1724(ss,ilo1,0,the)*fac
804     if(isinel(320).eq.0)then
        si2=0.
        goto 805
        endif
        ik1=321
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3122),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
805     if(isinel(322).eq.0)then
        si3=0.
        goto 806
        endif
        ik1=311
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3112),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
806    if(isinel(137).eq.0)then
        si4=0.
        goto 1806
        endif
        the=pmas(pycomp(-311),1)+pmas(pycomp(-3334),1)
        si4=s1724(ss,ilo4,0,the)
1806    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &   ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=321
        ik2=-3212
        ic=319
c       pion+ + cascade0- to k+ + sigma0-
        goto 807
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=321
        ik2=-3122
        ic=320
c       pion+ + cascade0- to k+ + lambda-
        goto 807
        endif
       if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=311
        ik2=-3112
        ic=322
c       pion+ + cascade0- to k0 + sigma-ba
        goto 807
        endif
        ik1=-311
        ik2=-3334
        ic=137
cpion+ + cascade0- to k0- + omiga-ba
807     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + cascade0 
        if(kl.eq.-211 .and. kl1.eq.3322)then
        if(isinel(312).eq.0)then
        si1=0.
        goto 808
        endif
        ik1=-321
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3212),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
808     if(isinel(313).eq.0)then
        si2=0.
        goto 809
        endif
        ik1=-321
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3122),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
809     if(isinel(315).eq.0)then
        si3=0.
        goto 810
        endif
        ik1=-311
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3112),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
810     if(isinel(135).eq.0)then
        si4=0.
        goto 1810
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3334),1)
        si4=s1724(ss,ilo4,0,the)

1810    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &   ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-321
        ik2=3212
        ic=312
c       pion- + cascade0 to k- + sigma0
        goto 811
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-321
        ik2=3122
        ic=313
c       pion- + cascade0 to k- + lambda
        goto 811
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-311
        ik2=3112
        ic=315
c       pion- + cascade0 to k0- + sigma-
        goto 811
        endif
        ik1=311
        ik2=3334
        ic=135
c       pion- + cascade0 to k0 + omiga-
811     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + cascade0- to k0 + sigma+ba
        if(kl.eq.-211 .and. kl1.eq.-3322)then
        if(isinel(321).eq.0)goto 13
        the=pmas(pycomp(311),1)+pmas(pycomp(-3222),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=311
        lc(icp,4)=-3222
        lc(icp,5)=321
        ik3=311
        ik4=-3222
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion0 + cascade0 
        if(kl.eq.111 .and. kl1.eq.3322)then
        if(isinel(311).eq.0)then
        si1=0.
        goto 812
        endif
        ik1=-321
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3222),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
812     if(isinel(316).eq.0)then
        si2=0.
        goto 813
        endif
        ik1=-311
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3212),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
813     if(isinel(317).eq.0)then
        si3=0.
        goto 814
        endif
        ik1=-311
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3122),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
814     if(isinel(136).eq.0)then
        si4=0.
        goto 1814
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(2214),1)
        si4=s1724(ss,ilo4,0,the)
        
1814    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &   .and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-321
        ik2=3222
        ic=311
c       pion0 + cascade0 to k- + sigma+
        goto 815
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-311
        ik2=3212
        ic=316
c       pion0 + cascade0 to k0- + sigma0
        goto 815
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-311
        ik2=3122
        ic=317
c       pion0 + cascade0 to k0- + lambda
        goto 815
        endif
        ik1=321
        ik2=3334
        ic=136
c       pion0 + cascade0 to k+ + omiga-
815     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion0 + cascade0-
        if(kl.eq.111 .and. kl1.eq.-3322)then
        if(isinel(318).eq.0)then
        si1=0.
        goto 816
        endif
        ik1=321
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3222),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
816     if(isinel(323).eq.0)then
        si2=0.
        goto 817
        endif
        ik1=311
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3212),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
817     if(isinel(324).eq.0)then
        si3=0.
        goto 818
        endif
        ik1=311
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &   1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3122),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
818     if(isinel(138).eq.0)then
        si4=0.
        goto 1818
        endif
        the=pmas(pycomp(-321),1)+pmas(pycomp(-3334),1)
        si4=s1724(ss,ilo4,0,the)
        
1818    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &   ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=321
        ik2=-3222
        ic=318
c       pion0 + cascade0- to k+ + sigma+ba
        goto 819
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=311
        ik2=-3212
        ic=323
c       pion0 + cascade0- to k0 + sigma0-
        goto 819
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=311
        ik2=-3122
        ic=324
c       pion0 + cascade0- to k0 + lambda-
        goto 819
        endif
        ik1=-321
        ik2=-3334
        ic=138
c       pion0 + cascade0- to k- + omiga-ba
819     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion+ + omiga- to k0- + cascade0
        if(kl.eq.211 .and. kl1.eq.3334)then
        if(isinel(342).eq.0)goto 13
        the=pmas(pycomp(-311),1)+pmas(pycomp(3322),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=-311
        lc(icp,4)=3322
        lc(icp,5)=342
        ik3=-311
        ik4=3322
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion+ + omiga-ba to k+ + cascade-ba
        if(kl.eq.211 .and. kl1.eq.-3334)then
        if(isinel(343).eq.0)goto 13
        the=pmas(pycomp(321),1)+pmas(pycomp(-3312),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=321
        lc(icp,4)=-3312
        lc(icp,5)=343
        ik3=321
        ik4=-3312
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion- + omiga- to k- + cascade-
        if(kl.eq.-211 .and. kl1.eq.3334)then
        if(isinel(339).eq.0)goto 13
        the=pmas(pycomp(-321),1)+pmas(pycomp(3312),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=3312
        lc(icp,5)=339
        ik3=-321
        ik4=3312
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion- + omiga-ba to k0 + cascade0-
        if(kl.eq.-211 .and. kl1.eq.-3334)then
        if(isinel(346).eq.0)goto 13
        the=pmas(pycomp(311),1)+pmas(pycomp(-3322),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=311
        lc(icp,4)=-3322
        lc(icp,5)=346
        ik3=311
        ik4=-3322
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &   1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion0 + omiga- 
        if(kl.eq.111 .and. kl1.eq.3334)then
        if(isinel(340).eq.0)then
        si1=0.
        goto 820
        endif
        ik1=-311
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-311),1)+pmas(pycomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
820     if(isinel(341).eq.0)then
        si2=0.
        goto 821
        endif
        ik1=-321
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(-321),1)+pmas(pycomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
821     if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-311
        ik2=3312
        ic=340
c       pion0 + omiga- to k0- + cascade-
        goto 822
        endif
        ik1=-321
        ik2=3322
        ic=341
c       pion0 + omiga- to k- + cascade0
822     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si13/cspin
        goto 10
        endif

c       pion0 + omiga-ba 
        if(kl.eq.111 .and. kl1.eq.-3334)then
        if(isinel(344).eq.0)then
        si1=0.
        goto 823
        endif
        ik1=311
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(311),1)+pmas(pycomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
823     if(isinel(345).eq.0)then
        si2=0.
        goto 824
        endif
        ik1=321
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(pycomp(321),1)+pmas(pycomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
824     if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=311
        ik2=-3312
        ic=344
c       pion0 + omiga-ba to k0 + cascade-ba
        goto 825
        endif
        ik1=321
        ik2=-3322
        ic=345
c       pion0 + omiga-ba to k+ + cascade0-
825     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion+ + delta+ to pi+ + p
        if((kl.eq.211 .and. kl1.eq.2214).or.
     &   (kl1.eq.211 .and. kl.eq.2214))then
        if(isinel(351).eq.0)goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=351
        ww=sdelta(ss,ilo1,1,0.d0)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion0 + delta++ to pi+ + p
        if((kl.eq.111 .and. kl1.eq.2224).or.
     &   (kl1.eq.111 .and. kl.eq.2224))then
        if(isinel(352).eq.0)goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=352
        ww=sdelta(ss,ilo1,1,0.d0)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion+ + delta-
        if((kl.eq.211 .and. kl1.eq.1114).or.
     &   (kl1.eq.211 .and. kl.eq.1114))then
        if(isinel(347).eq.0)then
        si1=0.
        goto 1000
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion+ + delta- to pi- + p

1000    if(isinel(348).eq.0)then
        si2=0.
        goto 1002
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion+ + delta- to pi0 + n

1002    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=2212
        ic=347
c       cross section of pion+ + delta- to pi- + p
        goto 1004
        endif
        ik1=111
        ik2=2112
        ic=348
c       cross section of pion+ + delta- to pi0 + n
1004    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion+ + delta0
        if((kl.eq.211 .and. kl1.eq.2114).or.
     &   (kl1.eq.211 .and. kl.eq.2114))then
        if(isinel(349).eq.0)then
        si1=0.
        goto 1006
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion+ + delta0 to pi+ + n

1006    if(isinel(350).eq.0)then
        si2=0.
        goto 1008
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion+ + delta0 to pi0 + p

1008    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=2112
        ic=349
c       cross section of pion+ + delta0 to pi+ + n
        goto 1010
        endif
        ik1=111
        ik2=2212
        ic=350
c       cross section of pion+ + delta0 to pi0 + p
1010    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif

c       pion0 + delta+
        if((kl.eq.111 .and. kl1.eq.2214).or.
     &   (kl1.eq.111 .and. kl.eq.2214))then
        if(isinel(353).eq.0)then
        si1=0.
        goto 1012
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion0 + delta+ to pi0 + p

1012    if(isinel(354).eq.0)then
        si2=0.
        goto 1014
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion0 + delta+ to pi+ + n

1014    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2212
        ic=353
c       cross section of pion0 + delta+ to pi0 + p
        goto 1016
        endif
        ik1=211
        ik2=2112
        ic=354
c       cross section of pion0 + delta+ to pi+ + n
1016    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c       pion0 + delta0
        if((kl.eq.111 .and. kl1.eq.2114).or.
     &   (kl1.eq.111 .and. kl.eq.2114))then
        if(isinel(355).eq.0)then
        si1=0.
        goto 1018
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion0 + delta0 to pi0 + n

1018    if(isinel(356).eq.0)then
        si2=0.
        goto 1020
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion0 + delta0 to pi- + p

1020    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2112
        ic=355
c       cross section of pion0 + delta0 to pi0 + n
        goto 1022
        endif
        ik1=-211
        ik2=2212
        ic=356
c       cross section of pion0 + delta0 to pi- + p
1022    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion0 + delta- to pi- + n
        if((kl.eq.111 .and. kl1.eq.1114).or.
     &   (kl1.eq.111 .and. kl.eq.1114))then
        if(isinel(357).eq.0)goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=357
        ww=sdelta(ss,ilo1,1,0.d0)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion- + delta++
        if((kl.eq.-211 .and. kl1.eq.2224).or.
     &   (kl1.eq.-211 .and. kl.eq.2224))then
        if(isinel(358).eq.0)then
        si1=0.
        goto 1024
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion- + delta++ to pi0 + p

1024    if(isinel(359).eq.0)then
        si2=0.
        goto 1026
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion- + delta++ to pi+ + n

1026    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2212
        ic=358
c       cross section of pion- + delta++ to pi0 + p
        goto 1028
        endif
        ik1=211
        ik2=2112
        ic=359
c       cross section of pion- + delta++ to pi+ + n
1028    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion- + delta+
        if((kl.eq.-211 .and. kl1.eq.2214).or.
     &   (kl1.eq.-211 .and. kl.eq.2214))then
        if(isinel(360).eq.0)then
        si1=0.
        goto 1030
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion- + delta+ to pi- + p

1030    if(isinel(361).eq.0)then
        si2=0.
        goto 1032
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of pion- + delta+ to pi0 + n

1032    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=2212
        ic=360
c       cross section of pion- + delta+ to pi- + p
        goto 1034
        endif
        ik1=111
        ik2=2112
        ic=361
c       cross section of pion- + delta+ to pi0 + n
1034    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion- + delta0 to pi- + n
        if((kl.eq.-211 .and. kl1.eq.2114).or.
     &   (kl1.eq.-211 .and. kl.eq.2114))then
        if(isinel(362).eq.0)goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=362
        ww=sdelta(ss,ilo1,1,0.d0)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       prho0 + n
        if((kl.eq.113 .and. kl1.eq.2112).or.
     &   (kl1.eq.113 .and. kl.eq.2112))then
        if(isinel(363).eq.0)then
        si1=0.
        goto 1036
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=srho(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of rho0 + n to pi- + p

1036    if(isinel(364).eq.0)then
        si2=0.
        goto 1038
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=srho(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of rho0 + n to pi0 + n

1038    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=2212
        ic=363
c       cross section of rho0 + n to pi- + p
        goto 1040
        endif
        ik1=111
        ik2=2112
        ic=364
c       cross section of rho0 + n to pi0 + n
1040    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       prho+ + n
        if((kl.eq.213 .and. kl1.eq.2112).or.
     &   (kl1.eq.213 .and. kl.eq.2112))then
        if(isinel(366).eq.0)then
        si1=0.
        goto 1042
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=srho(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of rho+ + n to pi0 + p

1042    if(isinel(367).eq.0)then
        si2=0.
        goto 1044
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=srho(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of rho+ + n to pi+ + n

1044    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2212
        ic=366
c       cross section of rho+ + n to pi0 + p
        goto 1046
        endif
        ik1=211
        ik2=2112
        ic=367
c       cross section of rho+ + n to pi+ + n
1046    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       rho0 + p
        if((kl.eq.113 .and. kl1.eq.2212).or.
     &   (kl1.eq.113 .and. kl.eq.2212))then
        if(isinel(368).eq.0)then
        si1=0.
        goto 1048
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=srho(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of rho0 + p to pi0 + p

1048    if(isinel(369).eq.0)then
        si2=0.
        goto 1050
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=srho(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of rho0 + p to pi+ + n

1050    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2212
        ic=368
c       cross section of rho0 + p to pi0 + p
        goto 1052
        endif
        ik1=211
        ik2=2112
        ic=369
c       cross section of rho0 + p to pi+ + n
1052    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       rho- + p
        if((kl.eq.-213 .and. kl1.eq.2212).or.
     &   (kl1.eq.-213 .and. kl.eq.2212))then
        if(isinel(370).eq.0)then
        si1=0.
        goto 1054
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)

        if(ilo1.eq.0)then
        si1=0.
        else
        si1=srho(ss,ilo1,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of rho- + p to pi0 + n

1054    if(isinel(371).eq.0)then
        si2=0.
        goto 1056
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=srho(ss,ilo2,1,0.d0)*WEIGH(19)*fac
        endif
c       cross section of rho- + p to pi- + p

1056    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2112
        ic=370
c       cross section of rho- + p to pi0 + n
        goto 1058
        endif
        ik1=-211
        ik2=2212
        ic=371
c       cross section of rho- + p to pi- + p
1058    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       rho+ + p to pi+ + p
        if((kl.eq.213 .and. kl1.eq.2212).or.
     &   (kl1.eq.213 .and. kl.eq.2212))then
        if(isinel(372).eq.0)goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=372
        ww=srho(ss,ilo1,1,0.d0)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       rho- + n to pi- + n
        if((kl.eq.-213 .and. kl1.eq.2112).or.
     &   (kl1.eq.-213 .and. kl.eq.2112))then
        if(isinel(365).eq.0)goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=365
        ww=srho(ss,ilo1,1,0.d0)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta++ + n to p + p
        if((kl.eq.2224 .and. kl1.eq.2112).or.
     &   (kl1.eq.2224 .and. kl.eq.2112))then
        if(isinel(373).eq.0)goto 13
        ik3=2212
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2212
        lc(icp,5)=373
        ww=snn(ss,ilo1,1,0.d0)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta+ + n to p + n
        if((kl.eq.2214 .and. kl1.eq.2112).or.
     &   (kl1.eq.2214 .and. kl.eq.2112))then
        if(isinel(374).eq.0)goto 13
        ik3=2212
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,1.0)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2112
        lc(icp,5)=374
        ww=snn(ss,ilo1,1,0.d0)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta+ + p to p + p
        if((kl.eq.2214 .and. kl1.eq.2212).or.
     &   (kl1.eq.2214 .and. kl.eq.2212))then
        if(isinel(375).eq.0)goto 13
        ik3=2212
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2212
        lc(icp,5)=375
        ww=snn(ss,ilo1,1,0.d0)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta0 + p to p + n
        if((kl.eq.2114 .and. kl1.eq.2212).or.
     &   (kl1.eq.2114 .and. kl.eq.2212))then
        if(isinel(376).eq.0)goto 13
        ik3=2212
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,1.0)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2112
        lc(icp,5)=376
        ww=snn(ss,ilo1,1,0.d0)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta0 + n to n + n
        if((kl.eq.2114 .and. kl1.eq.2112).or.
     &   (kl1.eq.2114 .and. kl.eq.2112))then
        if(isinel(377).eq.0)goto 13
        ik3=2112
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2112
        lc(icp,4)=2112
        lc(icp,5)=377
        ww=snn(ss,ilo1,1,0.d0)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta- + p to n + n
        if((kl.eq.1114 .and. kl1.eq.2212).or.
     &   (kl1.eq.1114 .and. kl.eq.2212))then
        if(isinel(378).eq.0)goto 13
        ik3=2112
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &   1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2112
        lc(icp,4)=2112
        lc(icp,5)=378
        ww=snn(ss,ilo1,1,0.d0)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


c       lambda- p annihilation
        if((kl.eq.-3122.and.kl1.eq.2212)
     c   .or.(kl1.eq.-3122.and.kl.eq.2212))then
        if(isinel(593).eq.0)goto 13
        lc(icp,3)=323
        lc(icp,4)=223
        lc(icp,5)=593
          tw(icp)=para13/5./csnn
        goto 10
        endif

c       lambda- n annihilation
        if((kl.eq.-3122.and.kl1.eq.2112)
     c   .or.(kl1.eq.-3122.and.kl.eq.2112))then
        if(isinel(594).eq.0)goto 13
        tw(icp)=para13/5./csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=313
        lc(icp,4)=223
        lc(icp,5)=594
        goto 10
        endif

c       sigma0- p annihilation
        if((kl.eq.-3212.and.kl1.eq.2212)
     c   .or.(kl1.eq.-3212.and.kl.eq.2212))then
        if(isinel(595).eq.0)goto 13
        tw(icp)=para13/5./csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=323
        lc(icp,4)=223
        lc(icp,5)=595
        goto 10
        endif

c       sigma0- n annihilation
        if((kl.eq.-3212.and.kl1.eq.2112)
     c   .or.(kl1.eq.-3212.and.kl.eq.2112))then
        if(isinel(596).eq.0)goto 13
        tw(icp)=para13/5./csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=313
        lc(icp,4)=223
        lc(icp,5)=596
        goto 10
        endif

c       p- p annihilation
        if((kl.eq.-2212.and.kl1.eq.2212)
     c   .or.(kl1.eq.-2212.and.kl.eq.2212))then
        if(isinel(597).eq.0)goto 13
        tw(icp)=para13/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=113
        lc(icp,4)=223
        lc(icp,5)=597
        goto 10
        endif

c       p- n annihilation
        if((kl.eq.-2212.and.kl1.eq.2112)
     c   .or.(kl1.eq.-2212.and.kl.eq.2112))then
        if(isinel(598).eq.0)goto 13
        tw(icp)=para13/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=-213
        lc(icp,4)=223
        lc(icp,5)=598
        goto 10
        endif

c       n- p annihilation
        if((kl.eq.-2112.and.kl1.eq.2212)
     c   .or.(kl1.eq.-2112.and.kl.eq.2212))then
        if(isinel(599).eq.0)goto 13
        tw(icp)=para13/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=213
        lc(icp,4)=223
        lc(icp,5)=599
        goto 10
        endif

c       n- n annihilation
        if((kl.eq.-2112.and.kl1.eq.2112)
     c   .or.(kl1.eq.-2112.and.kl.eq.2112))then
        if(isinel(600).eq.0)goto 13
        tw(icp)=para13/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=113
        lc(icp,4)=223
        lc(icp,5)=600
        goto 10
        endif

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       J/Psi interacts with baryon or meson
        if(kl.eq.443 .and. kl1.eq.2112)then
        if(iabsb.eq.0)goto 13
        call jpsin(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,1)
        if(ioo.eq.0)goto 13
        goto 10
        endif

        if(kl.eq.443 .and. kl1.eq.2212)then
        if(iabsb.eq.0)goto 13
        call jpsip(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,1)
        if(ioo.eq.0)goto 13
        goto 10
        endif

        if((kl.eq.443 .and. kl1.eq.211)
     &   .or. (kl.eq.211 .and. kl1.eq.443))then
        if(iabsm.eq.0)goto 13
        call jpsip1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,1)
        if(ioo.eq.0)goto 13
        goto 10
        endif

        if((kl.eq.443 .and. kl1.eq.111)
     &   .or. (kl.eq.111 .and. kl1.eq.443))then
        if(iabsm.eq.0)goto 13
        call jpsip0(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,1)
        if(ioo.eq.0)goto 13
        goto 10
        endif

        if((kl.eq.443 .and. kl1.eq.-211)
     &   .or. (kl.eq.-211 .and. kl1.eq.443))then
        if(iabsm.eq.0)goto 13
        call jpsip2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,1)
        if(ioo.eq.0)goto 13
        goto 10
        endif

        if((kl.eq.443 .and. kl1.eq.213)
     &   .or. (kl.eq.213 .and. kl1.eq.443))then
        if(iabsm.eq.0)goto 13
        call jpsir1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,1)
        if(ioo.eq.0)goto 13
        goto 10
        endif

        if((kl.eq.443 .and. kl1.eq.113)
     &   .or. (kl.eq.113 .and. kl1.eq.443))then
        if(iabsm.eq.0)goto 13
        call jpsir0(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,1)
        if(ioo.eq.0)goto 13
        goto 10
        endif

        if((kl.eq.443 .and. kl1.eq.-213)
     &   .or. (kl.eq.-213 .and. kl1.eq.443))then
        if(iabsm.eq.0)goto 13
        call jpsir2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,1)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       Psi' interacts with baryon or meson
        if(kl.eq.100443 .and. kl1.eq.2112)then
        if(iabsb.eq.0)goto 13
        call jpsin(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,2)
c       argument '2' here refers to the Psi' induced reaction 
c       corresponding onse
        if(ioo.eq.0)goto 13
        icp5=lc(icp,5)
        if(icp5.le.186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
c       14 and 192 are the distance (in order number) of inelastic reaction 
c        between J/Psi induced reaction and Psi, induced corresponding 
c        reaction 
        goto 10
        endif

        if(kl.eq.100443 .and. kl1.eq.2212)then
        if(iabsb.eq.0)goto 13
        call jpsip(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,2)
        if(ioo.eq.0)goto 13
        icp5=lc(icp,5)
        if(icp5.le.186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl.eq.100443 .and. kl1.eq.211)
     &   .or. (kl.eq.211 .and. kl1.eq.100443))then
        if(iabsm.eq.0)goto 13
        call jpsip1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,2)
        if(ioo.eq.0)goto 13
        icp5=lc(icp,5)
        if(icp5.le.186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl.eq.100443 .and. kl1.eq.111)
     &   .or. (kl.eq.111 .and. kl1.eq.100443))then
        if(iabsm.eq.0)goto 13
        call jpsip0(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,2)
        if(ioo.eq.0)goto 13
        icp5=lc(icp,5)
        if(icp5.le.186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl.eq.100443 .and. kl1.eq.-211)
     &   .or. (kl.eq.-211 .and. kl1.eq.100443))then
        if(iabsm.eq.0)goto 13
        call jpsip2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,2)
        if(ioo.eq.0)goto 13
        icp5=lc(icp,5)
        if(icp5.le.186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl.eq.100443 .and. kl1.eq.213)
     &   .or. (kl.eq.213 .and. kl1.eq.100443))then
        if(iabsm.eq.0)goto 13
        call jpsir1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,2)
        if(ioo.eq.0)goto 13
        icp5=lc(icp,5)
        if(icp5.le.186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl.eq.100443 .and. kl1.eq.113)
     &   .or. (kl.eq.113 .and. kl1.eq.100443))then
        if(iabsm.eq.0)goto 13
        call jpsir0(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,2)
        if(ioo.eq.0)goto 13
        icp5=lc(icp,5)
        if(icp5.le.186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        if((kl.eq.100443 .and. kl1.eq.-213)
     &   .or. (kl.eq.-213 .and. kl1.eq.100443))then
        if(iabsm.eq.0)goto 13
        call jpsir2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,2)
        if(ioo.eq.0)goto 13
        icp5=lc(icp,5)
        if(icp5.le.186)then
        lc(icp,5)=icp5+14
        else
        lc(icp,5)=icp5+192
        endif
        goto 10
        endif

        write(9,*)'nothing is found in hadcas prod, kl, kl1=', kl, kl1
13      do m=3,5
        lc(icp,m)=0
        enddo
        tw(icp)=0.
        goto 999
10      if(lc(icp,5).lt.1593)then
        nchargef=PYCHGE(lc(icp,3))+PYCHGE(lc(icp,4))
        if(nchargei/3.ne.nchargef/3)then
        write(9,*)'initial,final,lc(icp,5)=',
     c   nchargei/3,nchargef/3,lc(icp,5)
        endif
        endif
999     return
        end



c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine spipi(lc3,lc4,ss,ilo)
c       a part of prod()
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        am1=pmas(pycomp(lc3),1)
        am2=pmas(pycomp(lc4),1)
        ilo=1
        if(ss.lt.(am1+am2))ilo=0
        return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine srev(kl,kl1,lc3,lc4,ss,ilo,fac,xii1,xii2,xsi1,xsi2,
     &   xif1,xif2,xsf1,xsf2,pauli)
c       a part of prob()
c181002
c       xii1,xii2: the isospin of of two particles in initial state
c       xsi1,xsi2: the spin of two particles in initial state
c       xif1,xif2: the isospin of of two particles in final state
c       xsf1,xsf2: the spin of two particles in final state
c       paul1: =1 if two particles are different in final state
c              =0.5 if identical two particles in final state
c181002
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        real*4 xii1,xii2,xsi1,xsi2,xif1,xif2,xsf1,xsf2,pauli
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        am1=pmas(pycomp(kl),1)
        am2=pmas(pycomp(kl1),1)
        am3=pmas(pycomp(lc3),1)
        am4=pmas(pycomp(lc4),1)
        am12=am1*am1
        am22=am2*am2
        am32=am3*am3
        am42=am4*am4
        ss2=ss*ss
        ss4=ss2*ss2
        ilo=1

        if(ss.lt.(am3+am4))then
        ilo=0
        fac=0.
        else
        pfp=ss4-2.*ss2*(am32+am42)+(am32-am42)*(am32-am42)
        pfp=pfp/4./ss2   ! 311002
        if(pfp.lt.0.)pfp=1.e-10
c311002 pfp=sqrt(pfp)/2./ss
        pip=ss4-2.*ss2*(am12+am22)+(am12-am22)*(am12-am22)
        pip=pip/4./ss2   ! 311002
        if(pip.lt.0.)pip=1.e-10
c311002 pip=sqrt(pip)/2./ss
c       phase=pauli*(2*xif1+1)*(2*xif2+1)*(2*xsf1+1)*(2*xsf2+1)/
c     &   ((2*xii1+1)*(2*xii2+1)*(2*xsi1+1)*(2*xsi2+1))
        phase=pauli*(2*xsf1+1)*(2*xsf2+1)/((2*xsi1+1)*(2*xsi2+1))
        fac=phase*pfp/pip
        endif
        return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine updpli(l,l1,icp,ss,pii,pjj,lc,tc,tw,winel,time,iia)
c       update particle list after inelastic collision &
c        truncate collision list correspondingly
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp   ! 250208 yan
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
c       note the name of the arrays in 'sa1' in this subroutine
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa19_h/coor(3)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,ecsspn,ecsspm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension pii(4),pjj(4),pp(4)
        integer winel
        kf1=k(l,2)
        kf2=k(l1,2)
c       do i=1,3
c       v(L,i)=(v(L,i)+v(L1,i))*0.5
c       v(L1,i)=v(L,i)
c       enddo
cc       n1=nwhole-n
        if(iia.ge.1593)goto 1300
        goto 1400
c       update the collision list for annihilation case
1300    ic=l
        jc=l1
        j=0
        do i=1,nctl
        i1=lc(i,1)
        j1=lc(i,2)
        if(i1.eq.ic .or. i1.eq.jc)goto 4000
        if(j1.eq.ic .or. j1.eq.jc)goto 4000
        if((tc(i)-time).le.ddt) goto 4000
c       through away the pairs with tc<= time
        j=j+1
        tc(j)=tc(i)
        tw(j)=tw(i)
        do m=1,5
        lc(j,m)=lc(i,m)
        enddo
4000    continue
        enddo
        do i=j+1,nctl
        tc(i)=0.0
        tw(i)=0.0
        do m=1,5
        lc(i,m)=0
        enddo
        enddo
        nctl=j
        goto 1200
1400    ik1=lc(icp,3)
        ik2=lc(icp,4)
c       put the scattered (produced) particles into particle list
c    & turncate collision list correspondingly.
        ll=l
        ll1=l1
        kf=ik1
        do i=1,4
        pp(i)=pii(i)
        enddo
        do 500 i=1,2
c220601
        if(n.ge.kszj)then
        write(9,*)'in updpli, n > kszj'
        stop 'updpli'
        endif
c220601
cc      nnn=0

        do 600 j=1,kfmax
cc       if(kf.ne.kfaco(j))goto 900
        if(kf.ne.kfaco(j))goto 600
        jj=numb(j)+1
c       update the particle list.
cc1000  do m=n+n1,jj,-1
1000    do m=n,jj,-1
        mm=m+1
        k(mm,2)=k(m,2)
        k(mm,1)=1
        k(mm,3)=k(m,3)
        do m1=1,5
        p(mm,m1)=p(m,m1)
        v(mm,m1)=v(m,m1)
        enddo
        ishp(mm)=ishp(m)
        tau(mm)=tau(m)
        enddo
        if(ll.ge.jj)ll=ll+1
        if(ll1.ge.jj)ll1=ll1+1
c       update the values of lc(m,1-2) 
        do m=1,nctl
        lc1=lc(m,1)
        if(lc1.ge.jj)lc(m,1)=lc1+1
        lc2=lc(m,2)
        if(lc2.ge.jj)lc(m,2)=lc2+1
        enddo
c       ll instead of jj all over the collision list
        do m=1,nctl
        lc1=lc(m,1)
        lc2=lc(m,2)
        if(lc1.eq.ll)lc(m,1)=jj
        if(lc2.eq.ll)lc(m,2)=jj
        enddo
c       give proper values to particle jj.
        k(jj,2)=kf
        k(jj,1)=1
        k(jj,3)=0
        do m=1,4
        p(jj,m)=pp(m)
        v(jj,m)=v(ll,m)
        enddo
        p(jj,4)=pp(4)
        p(jj,5)=pmas(pycomp(kf),1)
        ishp(jj)=ishp(ll)
c       add photon
        if(k(jj,2).eq.22 .or. kf.eq.44 .or. kf.eq.55 .or. kf.eq.66)then
        tau(jj)=time
        goto 301
        endif
c       tau(jj)=t0*p(jj,4)/p(jj,5)
        tau(jj)=time+t0*taup*p(jj,4)/p(jj,5)
c       give zero formation time to particles produced from
c        rescatttering 
301     lc(icp,i)=jj
cc      if(nnn.eq.1)goto 300
        do m=j,kfmax
        numb(m)=numb(m)+1
        enddo
        goto 2000
cc300   goto 200
cc900   if(j.lt.kfmax)goto 600
c       scattered particle with new flavor.
cc      nnn=1
cc      jj=numb(kfmax)+1
cc      kfmax=kfmax+1
cc      numb(kfmax)=jj
cc      kfaco(kfmax)=kf
cc      goto 1000
600     continue
c030706
c       give proper values to particle not considered in 'kfaco'
        jj=n+1
        k(jj,2)=kf
        k(jj,1)=1
        k(jj,3)=0
        do m=1,4
        p(jj,m)=pp(m)
        v(jj,m)=v(ll,m)
        enddo
        p(jj,5)=pmas(pycomp(kf),1)
        ishp(jj)=ishp(ll)
c       tau(jj)=t0*p(jj,4)/p(jj,5)
        tau(jj)=time+t0*taup*p(jj,4)/p(jj,5)
c       give zero formation time to particles produced from
c        rescatttering 
        lc(icp,i)=jj
c030706
cc200   continue
2000    n=n+1
        if(i.eq.2)goto 500
        ll2=ll
        ll=ll1
        ll1=ll2
        kf=ik2
        do j=1,4
        pp(j)=pjj(j)
        enddo
500     continue
        l=ll1
        l1=ll
1100    continue
c       remove the scattering particles from particle list &
c        truncate the collision list correspondingly.
1200    kf=kf1
        ll=l
        do 700 i=1,2
c011210
        if(ll.eq.n)then   !
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
cc       do j=ll+1,n+n1
        do j=ll+1,n
        jj=j-1
        k(jj,2)=k(j,2)
        k(jj,1)=1
        k(jj,3)=k(j,3)
        do m=1,5
        p(jj,m)=p(j,m)
        v(jj,m)=v(j,m)
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
        n=n-1
        if(l1.gt.ll)l1=l1-1
        if(i.eq.2)goto 700
        ll=l1
        kf=kf2
700     continue
1600    continue   ! 010706
        return   ! 010706 270407
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine coelas_h(ic,jc,eij,pi,pj)
c       perform elastic scattering
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
c       note the name of the arraies in 'sa1_h' in this subroutine
        dimension pi(4),pj(4)
        iic=k(ic,2)
        jjc=k(jc,2)
c100111 d=3.65*(eij-pmas(pycomp(iic),1)-pmas(pycomp(jjc),1))
c100111 if(d.lt.1.e-10)return
        pt=0.5   ! changed from 0.2 to 0.5 on 100111
        a=min(10.3,1./(1.12*pt)/(1.12*pt))
c100111 d6=d**6
c100111 b=d6*a/(1.+d6)
        b=a   ! 100111
c100111 if(b.lt.1.e-20)then
c100111 b=1.e-20
c100111 endif
        pm2=pi(1)**2+pi(2)**2+pi(3)**2
        pm=sqrt(pm2)
        t0=-4.*pm2
c100111 if(abs(t0).lt.1.e-20)then
c100111 cctas=1.
c100111 maximum of cos(\theta_s)=1
c100111 goto 100
c100111 endif
        cc=pyr(1)
c100111 if(abs(b*t0).lt.0.0001)then
c       abt=exp(b*t0),b*t0<0    100111
c100111 abt=1.
c       elseif(b*t0.lt.-50.)then
c       abt=0.
c100111 else
        abt=dexp(dmax1(-7.0D2,dble(b*t0)))   ! mathematical skill 100111
c       'dble': intrinsic function, to increase the precision from 'real 4'
c        to 'real 8' 100111
c100111 endif
        tt1=dlog(cc+(1.-cc)*abt)
c100111 if(abs(tt1).lt.1.e-30 .and. b.le.1.e-20)then
c100111 cctas=1.
c100111 goto 100
c100111 endif
        tt=tt1/b
c100111 if(abs(tt).lt.1.e-20)then
c100111 cctas=1.
c100111 goto 100
c100111 endif
        cctas=1.-tt*2./t0
        if(dabs(cctas).gt.1.d0)then
        cctas=dsign(1.d0,cctas)
        endif
100     continue
        sctas=sqrt(1.-cctas**2)
        fis=2.*3.1416*pyr(1)
        cfis=cos(fis)
        sfis=sin(fis)
        call rotate_h(cctas,sctas,cfis,sfis,pm,pi,pj)
        return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine rotate_h(cctas,sctas,cfis,sfis,pp3,pi,pj)
c       perform the rotation after elastic scattering
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        dimension pi(4),pj(4)
c       fi1=atan2(pi(2),pi(1))
c       cta1=atan2(sqrt(pi(1)**2+pi(2)**2),pi(3))
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

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine updple_h(ic,jc,b,pi,pj,time)
c       update the particle list after elastic scattering 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
c       note the name of the arrays in 'sa2' in this subroutine
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa19_h/coor(3)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,ecsspn,ecsspm
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

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine updatl_h(ic,jc,time,lc,tc,tw,winel,iii)
c       update the collision time list
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
c        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,ecsspn,ecsspm
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 141104
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/syspar_h/pio
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        integer winel

        dddt=adj1(11)   !141104
c       loop over old colliding pairs
        j=0
        if(nctl.eq.0)goto 370
        do i=1,nctl
        i1=lc(i,1)
        j1=lc(i,2)
c       ia=(i1-ic)*(j1-jc)*(i1-jc)*(j1-ic)
c       if(ia.eq.0) goto 400
        if(i1.eq.ic .or. i1.eq.jc)goto 400
        if(j1.eq.ic .or. j1.eq.jc)goto 400
        if(abs(tc(i)-time).le.dddt) goto 400   ! 200204
c       through away the pairs which have tc<= time
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

370     nctl=j+1
c       loop over particle list

        m2=numb(2)
        m4=numb(4)
        m7=numb(7)
        m9=numb(9)
        m17=numb(17)
        m19=numb(19)
        m25=numb(25)
        m29=numb(29)
        m32=numb(32)
        m34=numb(34)
c       m34=numb(kfmax-11)
c       subtract 11, since do not consider the rescattering of x0c, etc   
        j1=ic
        do ik=1,2

        if(j1.le.m19 
     c   .or. (j1.gt.m25 .and. j1.le.m34))goto 301
c       consider only the reinteraction among nucleon, pion, kaon,
c        sigma, lambda, delta, rho and psi
        goto 300

301     do i=1,m34
        if(j1.eq.ic .and. i.eq.jc)goto 600
        if(j1.eq.jc .and. i.eq.ic)goto 600
c       forbiden scattered particles colliding with each other immediately

        if(nctl.gt.nsize)then
        write(9,*)'2 nsa,nctl,dddt=',nsa,nctl,dddt
        write(9,*)'size of array "nsize" needs to be extended'
        write(9,*)'error is serious,stop running'
        stop 30000
        endif

        if(i.le.m19 
     c   .or. (i.gt.m25 .and. i.le.m34))goto 602
c       consider only the reinteraction among nucleon, pion, kaon,
c        sigma, lambda, delta, rho and psi
        goto 600

csa601     if(j1.le.m2 .and. i.le.m2.and.j1.lt.i)goto 600
c       no double counting for NN
csa        if((j1.gt.m2 .and. j1.le.m4) .and. i.gt.m25) goto 600
csa        if((i.gt.m2 .and. i.le.m4) .and. j1.gt.m25) goto 600
c       delta, rho, psi etc. do not scatter with Nbar
csa        if((j1.gt.m7.and.i.gt.m25).or.((j1.gt.m4 .and. j1.le.m7) .and. 
csa     c   (i.gt.m29 .and. i.le.m32)))goto 600
csa        if((i.gt.m7.and.j1.gt.m25).or.((i.gt.m4 .and. i.le.m7) .and. 
csa     c   (j1.gt.m29 .and. j1.le.m32)))goto 600
c       delta, rho etc. do not collid with the strange particles and
c       rho do not collid with pi
csa        if((j1.gt.m25 .and. j1.le.32) .and. (i.gt.m25 .and. i.le.32))
csa     c   goto 600
c       delta and rho do not scatter with each other

c       if((j1.gt.m32 .and. j1.le.m34) .and. ((i.gt.m7 .and. i.le.m9)
c     c   .or. (i.gt.m17 .and. i.le.m19)))goto 600
c       if((i.gt.m32 .and. i.le.m34) .and. ((j1.gt.m7 .and. j1.le.m9)
c     c   .or. (j1.gt.m17 .and. j1.le.m19)))goto 600
c       psi do not collid with k

602     i1=i
        iflag=0
        call rsfilt_h(j1,i1,iflag)
        if(iflag.eq.0)goto 100
        tc(nctl)=0.0
        call tcolij_h(j1,i1,time,nctl,lc,tc,tw)
c       if(tc(nctl).gt.1.0e-7) nctl=nctl+1   ! 141104
c141104
        tci=tc(nctl)
        if(tci.gt.1.0e-7)then
        if(nctl.gt.1)then   ! 031204
        do j1=1,nctl-1
        tcj=tc(j1)
        if(abs(tcj-tci).lt.dddt)goto 600
        enddo
        endif   ! 031204
        nctl=nctl+1
        endif
c141104
100     continue
600     enddo
300     if(ik.eq.2)goto 500
        j1=jc
500     enddo
700     if(tc(nctl).le.1.e-7) nctl=nctl-1
        return
        end



C*******************************************************************
        subroutine rsfilt_h(l,l1,iflag)
c       play the role of first range filter and guarantee the collision list 
c        is composed according to the entrance channels of considered 
c        inelastic reactions
c        subroutine intdis plays the role of second range filter
c       collision pairs not interested can not filter through both of rsfilt 
c        and intdis
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa19_h/coor(3)
        m2=numb(2)
        m4=numb(4)
        kl=k(l,2)
        kl1=k(l1,2)
        if(l.eq.l1) goto 10
        if(ishp(l).eq.0.or.ishp(l1).eq.0) goto 10

c       constraints on the direct reactions
        if(kl.eq.211 .and. (kl1.eq.-211 .or. kl1.eq.111 .or.
     &   abs(kl1).eq.2212 .or. abs(kl1).eq.2112 .or. kl1.eq.
     &   3112 .or. kl1.eq.-3122 .or. kl1.eq.-3222 .or. kl1
     &   .eq.-3212 .or. kl1.eq.3212 .or. kl1.eq.3122 .or. kl1.eq.3312
     &   .or. kl1.eq.-3322))goto 11
        if(kl1.eq.211 .and. (kl.eq.-211 .or. kl.eq.111 .or.
     &   abs(kl).eq.2212 .or. abs(kl).eq.2112 .or. kl.eq.
     &   3112 .or. kl.eq.-3122 .or. kl.eq.-3222 .or. kl
     &   .eq.-3212 .or. kl.eq.3212 .or. kl.eq.3122 .or. kl.eq.3312
     &   .or. kl.eq.-3322))goto 11
        if(kl.eq.-211 .and. (kl1.eq.111 .or.
     &   abs(kl1).eq.2212 .or. abs(kl1).eq.2112 .or. kl1.eq.
     &   -3112 .or. kl1.eq.3122 .or. kl1.eq.3222 .or. kl1
     &   .eq.3212 .or. kl1.eq.-3212 .or. kl1.eq.-3122 .or. kl1.eq.
     &   -3312 .or. kl1.eq.3322))goto 11
        if(kl1.eq.-211 .and. (kl.eq.111 .or.
     &   abs(kl).eq.2212 .or. abs(kl).eq.2112 .or. kl.eq.
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
     &   .or. kl1.eq.-3212 .or.kl1.eq.-3312 .or.kl1.eq.-3322))goto 11
        if(kl1.eq.321 .and. (kl.eq.-2212 .or. kl.eq.-2112 .or. 
     &   kl.eq.-3122 .or. kl.eq.-3222 .or. kl.eq.-3112
     &   .or. kl.eq.-3212 .or.kl.eq.-3312 .or.kl.eq.-3322))goto 11
        if(kl.eq.-321 .and. (kl1.eq.2212 .or. kl1.eq.2112 .or. 
     &   kl1.eq.3122 .or. kl1.eq.3222 .or. kl1.eq.3112
     &   .or. kl1.eq.3212 .or.kl1.eq.3312 .or. kl1.eq.3322))goto 11
        if(kl1.eq.-321 .and. (kl.eq.2212 .or. kl.eq.2112 .or. 
     &   kl.eq.3122 .or. kl.eq.3222 .or. kl.eq.3112
     &   .or. kl.eq.3212 .or.kl.eq.3312 .or. kl.eq.3322))goto 11
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
        if((kl.eq.443 .or. kl.eq.100443) .and. (kl1.eq.2212 .or. kl1.eq.
     &   2112 .or. kl1.eq.211 .or. kl1.eq.111 .or. kl1.eq.-211 
     &   .or. kl1.eq.213 .or. kl1.eq.113 .or. kl1.eq.-213))goto 11 
        if((kl1.eq.443 .or. kl1.eq.100443) .and. (kl.eq.2212 .or. kl.eq.
     &   2112 .or. kl.eq.211 .or. kl.eq.111 .or. kl.eq.-211 
     &   .or. kl.eq.213 .or. kl.eq.113 .or. kl.eq.-213))goto 11

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        if(kl.eq.2112.and.(kl1.eq.2112.or.kl1.eq.2212))goto 11
        if(kl.eq.2212.and.(kl1.eq.2112.or.kl1.eq.2212))goto 11
c NN scattering is chosen
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c       constraints on the annihilation reactions
        if(kl.eq.-3212 .and. (kl1.eq.2212 .or. kl1.eq.2112))goto 11
        if(kl.eq.-3122 .and. (kl1.eq.2212.or. kl1.eq.2112))goto 11
        if(kl1.eq.-3212 .and. (kl.eq.2212 .or. kl.eq.2112))goto 11
        if(kl1.eq.-3122 .and. (kl.eq.2212.or. kl.eq.2112))goto 11
        if(kl.eq.-2212 .and. (kl1.eq.2212 .or. kl1.eq.2112))goto 11
        if(kl1.eq.-2212 .and. (kl.eq.2212 .or. kl.eq.2112))goto 11
        if(kl.eq.-2112 .and. (kl1.eq.2212.or. kl1.eq.2112))goto 11
        if(kl1.eq.-2112 .and. (kl.eq.2212.or. kl.eq.2112))goto 11

c       constraints on the reverse reactions
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
     &   kl.eq.2214.or.abs(kl).eq.213.or.kl.eq.113))goto 11
        if(kl1.eq.2112.and.(kl.eq.2224.or.kl.eq.2114.or.
     &   kl.eq.2214.or.abs(kl).eq.213.or.kl.eq.113))goto 11
        if(kl.eq.2212.and.(kl1.eq.1114.or.kl1.eq.2114.or.
     &   kl1.eq.2214.or.abs(kl1).eq.213.or.kl1.eq.113))goto 11
        if(kl.eq.2112.and.(kl1.eq.2224.or.kl1.eq.2114.or.
     &   kl1.eq.2214.or.abs(kl1).eq.213.or.kl1.eq.113))goto 11
        goto 10

11      iflag=1
10      continue
        return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine tcolij_h(l,l1,time,icp,lc,tc,tw)
c       calculate the collision time & fill up lc(i,1-2),tc(i). 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/sa19_h/coor(3)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,ecsspn,ecsspm
        common/sa24/adj1(40),nnstop,non24,zstop   ! 031204
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension dr(3),db(3),pi(4),pj(4),vi(3),vj(3)
        dimension ri(4),rj(4),rfi(4),rfj(4),b(3)
        drmax=adj1(28)   ! 031204
        pel=p(l,4)
        pel1=p(l1,4)
        if(pel.lt.1.e-10)pel=1.e-10   ! 031204
        if(pel1.lt.1.e-10)pel1=1.e-10   ! 031204
        pi(4)=pel
        pj(4)=pel1
        do i=1,3
        pi(i)=p(l,i)
        pj(i)=p(l1,i)
        b(i)=(pi(i)+pj(i))/(pi(4)+pj(4))
        enddo
        ilo=0
c       perform Lorentz transf. to CMS frame
        call lorntz(ilo,b,pi,pj)
        bta=dsqrt(b(1)**2+b(2)**2+b(3)**2)
c       if boost is too violent,put particles on mass shell by hand.
        if(bta.gt.0.99999d+0)then
        bmi=pmas(pycomp(k(l,2)),1)
        bmj=pmas(pycomp(k(l1,2)),1)
        pi(4)=sqrt(bmi**2+pi(1)**2+pi(2)**2+pi(3)**2)
        pj(4)=sqrt(bmj**2+pj(1)**2+pj(2)**2+pj(3)**2)
        endif
        ss=pi(4)+pj(4)
c092200
        ecut0=0.02    ! P.Yang
        ecut=ss-bmi-bmj
        if(ecut.le.ecut0) return
c092200
c       if ss < threshold collision may not happen
c311002
c311002 for isothermal direct reactions
c311002 pipi->kk
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &   (abs(k(l1,2)).eq.211.or.k(l1,2).eq.111)).and.ss.le.
     &   2.*pmas(pycomp(321),1))goto 10
c311002 piN->ky
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &   (abs(k(l1,2)).eq.2112.or.abs(k(l1,2)).eq.2212)).and.ss.le.
     &   (pmas(pycomp(321),1)+1.2))goto 10   ! 311002
        if(((abs(k(l1,2)).eq.211.or.k(l1,2).eq.111).and.
     &   (abs(k(l,2)).eq.2112.or.abs(k(l,2)).eq.2212)).and.ss.le.
     &   (pmas(pycomp(321),1)+1.2))goto 10   !311002
c       piy->k(cascade)
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &   (abs(k(l1,2)).eq.3122.or.abs(k(l1,2)).eq.3212.or.k(l1,2)
     &   .eq.3112.or.k(l1,2).eq.3222)).and.ss.le.
     &   (0.5+1.322))goto 10
        if(((abs(k(l1,2)).eq.211.or.k(l1,2).eq.111).and.
     &   (abs(k(l,2)).eq.3122.or.abs(k(l,2)).eq.3212.or.k(l,2)
     &   .eq.3112.or.k(l,2).eq.3222)).and.ss.le.
     &   (0.5+1.322))goto 10
c       pi(cascade)->k(Omega)
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &   (abs(k(l1,2)).eq.3312.or.abs(k(l1,2)).eq.3322)).and.ss.le.
     &   (0.5+1.7))goto 10
        if(((abs(k(l1,2)).eq.211.or.k(l1,2).eq.111).and.
     &   (abs(k(l,2)).eq.3312.or.abs(k(l,2)).eq.3322)).and.ss.le.
     &   (0.5+1.7))goto 10
c       kN->k(cascade)
        if(((abs(k(l,2)).eq.321.or.abs(k(l,2)).eq.311).and.
     &   (abs(k(l1,2)).eq.2112.or.abs(k(l1,2)).eq.2212)).and.ss.le.
     &   (0.5+1.322))goto 10
        if(((abs(k(l1,2)).eq.321.or.abs(k(l1,2)).eq.311).and.
     &   (abs(k(l,2)).eq.2112.or.abs(k(l,2)).eq.2212)).and.ss.le.
     &   (0.5+1.322))goto 10
c       piN->pi(Delta)
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &   (abs(k(l1,2)).eq.2112.or.abs(k(l1,2)).eq.2212)).and.ss.le.
     &   (0.14+1.234))goto 10
        if(((abs(k(l1,2)).eq.211.or.k(l1,2).eq.111).and.
     &   (abs(k(l,2)).eq.2112.or.abs(k(l,2)).eq.2212)).and.ss.le.
     &   (0.14+1.234))goto 10
c       piN->(rho)N
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &   (abs(k(l1,2)).eq.2112.or.abs(k(l1,2)).eq.2212)).and.ss.le.
     &   1.71)goto 10
        if(((abs(k(l1,2)).eq.211.or.k(l1,2).eq.111).and.
     &   (abs(k(l,2)).eq.2112.or.abs(k(l,2)).eq.2212)).and.ss.le.
     &   1.71)goto 10
c       NN->N(Delta)
        if(((abs(k(l,2)).eq.2112.or.abs(k(l,2)).eq.2212).and.
     &   (abs(k(l1,2)).eq.2112.or.abs(k(l1,2)).eq.2212)).and.ss.le.
     &   2.172)goto 10
c       k(cascade)->pi(Omega)
        if(((abs(k(l,2)).eq.321.or.abs(k(l,2)).eq.311).and.
     &   (abs(k(l1,2)).eq.3312.or.abs(k(l1,2)).eq.3322)).and.ss.le.
     &   1.84)goto 10
        if(((abs(k(l1,2)).eq.321.or.abs(k(l1,2)).eq.311).and.
     &   (abs(k(l,2)).eq.3312.or.abs(k(l,2)).eq.3322)).and.ss.le.
     &   1.84)goto 10
c       for isothermal reverse reactions
c       piy->kN
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &   (abs(k(l1,2)).eq.3122.or.abs(k(l1,2)).eq.3212.or.k(l1,2)
     &   .eq.3112.or.k(l1,2).eq.3222)).and.ss.le.
     &   1.44)goto 10
        if(((abs(k(l1,2)).eq.211.or.k(l1,2).eq.111).and.
     &   (abs(k(l,2)).eq.3122.or.abs(k(l,2)).eq.3212.or.k(l,2)
     &   .eq.3112.or.k(l,2).eq.3222)).and.ss.le.
     &   1.44)goto 10
c       pi(cascade)->ky
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &   (abs(k(l1,2)).eq.3312.or.abs(k(l1,2)).eq.3322)).and.ss.le.
     &   1.7)goto 10
        if(((abs(k(l1,2)).eq.211.or.k(l1,2).eq.111).and.
     &   (abs(k(l,2)).eq.3312.or.abs(k(l,2)).eq.3322)).and.ss.le.
     &   1.7)goto 10
c311002

        do i=1,4
        ri(i)=v(l,i)
        rj(i)=v(l1,i)
        enddo
cc      ri(4)=time
cc      rj(4)=time
c       Lorentz transf. to CMS frame 
        call lorntz(ilo,b,ri,rj)
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
c       gamli=p(l,4)/p(l,5)
c       gamlj=p(l1,4)/p(l1,5)
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
c       dot=-1
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
c       for collision occurs,time must one step ahead
cTai
        do ik=1,3
        dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))+tcol*db(ik)
        rtai=rtai+dr(ik)*dr(ik)
        enddo
c       gamai=pi(4)/pmas(pycomp(k(l,2)),1)
c       gamaj=pj(4)/pmas(pycomp(k(l1,2)),1)

c TAIAN

c       when collision happens,particles should already be produced
c       we give a zero formation time for particles produced from
c        rescatttering
ctai
        endif
        sg=rtai

        dmin=sqrt(sg)
20      call intdis_h(l,l1,ss,rsig)
c       'intdis': calculate the interaction distance between 
c        particles l & l1.
        if(dmin.gt.rsig)goto 10
c       distance between the two particles should be smaller than rsig
        do ik=1,3
        ri(ik)=ri(ik)+vi(ik)*(tcol-ri(4))
        rj(ik)=rj(ik)+vj(ik)*(tcol-rj(4))
        enddo
c       move along Newton trajectory in CMS
        ri(4)=tcol
        rj(4)=tcol
        ilo=1
c       transform back to Lab.
        call lorntz(ilo,b,ri,rj)
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
c       collision happens in the future
c031204
        do i=1,3
        ri(i)=v(l,i)+p(l,i)*(tcol-time)/pel-coor(i)
        rj(i)=v(l1,i)+p(l1,i)*(tcol-time)/pel1-coor(i)
        enddo
        rri=sqrt(ri(1)*ri(1)+ri(2)*ri(2)+ri(3)*ri(3))
        rrj=sqrt(rj(1)*rj(1)+rj(2)*rj(2)+rj(3)*rj(3))
c       particles must be inside the largest region considered
        if(rri.gt.drmax)goto 10
        if(rrj.gt.drmax)goto 10
c031204
        if(tcol.le.drmax)goto 18   ! 031204
        return   ! 031204
18      tc(icp)=tcol
        lc(icp,1)=l
        lc(icp,2)=l1
10      return
        end


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine intdis_h(l,l1,ss,rsig)
c       calculate the interaction distance between particles l 
c        and l1.
c       it plays also the role of second range filter
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        common/sa1_h/n,non1,k(kszj,5),p(kszj,5),v(kszj,5)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,ecsspn,ecsspm
        rsig=0.
        kl=k(l,2)
        kl1=k(l1,2)
        if(abs(kl).eq.2212 .or. abs(kl).eq.2112)idpl=1

        if(abs(kl).eq.443 .or. abs(kl).eq.100443)idpl=2

        if(abs(kl).eq.211 .or. kl.eq.111)idpl=3
        if(abs(kl).eq.321 .or. abs(kl).eq.311)idpl=4
        if(abs(kl).eq.3212 .or. abs(kl).eq.3112 .or. abs(kl).eq.3222
     c   .or. abs(kl).eq.3122 .or. abs(kl).eq.3312 .or. abs(kl).eq.
     c   3322 .or. abs(kl).eq.3334)idpl=5
        if(abs(kl).eq.213 .or. kl.eq.113)idpl=6
        if(kl.eq.1114 .or. kl.eq.2114.or.kl.eq.2214 .or. kl.eq.2224)
     &   idpl=7

        if(abs(kl1).eq.2212 .or. abs(kl1).eq.2112)idpl1=1

        if(abs(kl1).eq.443 .or. abs(kl1).eq.100443)idpl1=2   ! 98/03/24

        if(abs(kl1).eq.211 .or. kl1.eq.111)idpl1=3
        if(abs(kl1).eq.321 .or. abs(kl1).eq.311)idpl1=4
        if(abs(kl1).eq.3212 .or. abs(kl1).eq.3112 .or. abs(kl1)
     c   .eq.3222 .or. abs(kl1).eq.3122 .or. abs(kl1).eq.3312
     c   .or. abs(kl1).eq.3322 .or. abs(kl1).eq.3334)idpl1=5
        if(abs(kl1).eq.213 .or. kl1.eq.113)idpl1=6
        if(kl1.eq.1114 .or. kl1.eq.2114.or.kl1.eq.2214 .or. kl1.eq.2224)
     c   idpl1=7

        if(idpl.eq.1 .and. idpl1.eq.1)rsig=ecsnn

        if(idpl.eq.2 .and. idpl1.eq.1)then
        rsig=ecspsn
        if(kl.eq.100443)rsig=ecsspn
        endif
        if(idpl.eq.1 .and. idpl1.eq.2)then
        rsig=ecspsn
        if(kl1.eq.100443)rsig=ecsspn
        endif
        if(idpl.eq.2 .and. idpl1.eq.3)then
        rsig=ecspsm
        if(kl.eq.100443)rsig=ecsspm
        endif
        if(idpl.eq.3 .and. idpl1.eq.2)then
        rsig=ecspsm
        if(kl1.eq.100443)rsig=ecsspm
        endif

        if(idpl.eq.2 .and. idpl1.eq.6)then
        rsig=ecspsm
        if(kl1.eq.113)rsig=ecspsm*1.414
c       at the case of rho0, cross section enlarges a factor 2 to 
c        consider the effect of omega
        if(kl.eq.100443)then
        rsig=ecsspm
        if(kl1.eq.113)rsig=ecsspm*1.414
        endif
        endif
        if(idpl.eq.6 .and. idpl1.eq.2)then
        rsig=ecspsm
        if(kl.eq.113)rsig=ecspsm*1.414
        if(kl1.eq.100443)then
        rsig=ecsspm
        if(kl.eq.113)rsig=ecsspm*1.414
        endif
        endif

        if(idpl.eq.3 .and. idpl1.eq.3)rsig=edipi
        if(idpl.eq.1 .and. idpl1.eq.3)rsig=epin
        if(idpl.eq.3 .and. idpl1.eq.1)rsig=epin
        if(idpl.eq.3 .and. idpl1.eq.5)rsig=epin
        if(idpl.eq.5 .and. idpl1.eq.3)rsig=epin
c       assume the total cross section of (pion)y,((pion)cascade) and 
c        ((pion)omiga) = (pion)n
        if(idpl.eq.4 .and. (idpl1.eq.1 .or. idpl1.eq.5))rsig=ekn
        if((idpl.eq.1 .or. idpl.eq.5) .and. idpl1.eq.4)rsig=ekn
        if(idpl.eq.4 .and. idpl1.eq.4)rsig=edipi

c       assume the total cross section of ky (k cascade) and (k omiga)
c        = kn
c       assume the total cross section of kk=(pion)(pion)
        if(idpl.eq.1 .and. idpl1.eq.6)rsig=epin
        if(idpl.eq.1 .and. idpl1.eq.7)rsig=ecsnn
        if(idpl.eq.3 .and. idpl1.eq.7)rsig=epin
        if(idpl1.eq.1 .and. idpl.eq.6)rsig=epin
        if(idpl1.eq.1 .and. idpl.eq.7)rsig=ecsnn
        if(idpl1.eq.3 .and. idpl.eq.7)rsig=epin
        if(idpl.eq.1 .and. idpl1.eq.5)rsig=ecsnn
        if(idpl.eq.5 .and. idpl1.eq.1)rsig=ecsnn
c       assume the total cross section of (rho)n=(pion)n
c       assume the total cross section of n(delta)=nn
c       assume the total cross section of (pion)(delta)=(pion)n
c       assume the total cross section of ny, n(cascade), and n(omiga)=nn
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        function s0715(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
c       pion- + p to k+ + sigma-,channel 4
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.9) goto 10
        IF(kjp20.EQ.0)then
        si=0.25*(1.-0.75*(ss-1.691))
        else
        si=vjp20
        endif
        goto 100
10      IF(kjp20.EQ.0)then
        si=309.1*exp(max(-40.,-3.77*ss))
        else
        si=vjp20
        endif
100     continue
        s0715=si
        return
        end

c******************************************************************
        function s07122(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
c       pion- + p to k0 + lambda,channel 5
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.684) goto 10
        IF(kjp20.EQ.0)THEN
        si=0.9/0.091*(ss-1.613)
        ELSE
        si=vjp20
        ENDIF
        goto 100
10      if(ss.ge.2.1) goto 20
        IF(kjp20.EQ.0)THEN
        si=436.3*exp(max(-40.,-4.154*ss))
        ELSE
        si=vjp20
        ENDIF
        goto 100
20      IF(kjp20.EQ.0)THEN
        si=0.314*exp(max(-40.,-0.301*ss))
        ELSE
        si=vjp20
        ENDIF
100     s07122=si
        return
        end

c*******************************************************************
        function s07123(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
c       pion- + p to k0 + sigma0,channel 6
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.722) goto 10
        IF(kjp20.EQ.0)THEN
        si=10.6*(ss-1.689)
        ELSE
        si=vjp20
        ENDIF
        goto 100
10      if(ss.ge.3.) goto 20
        IF(kjp20.EQ.0)THEN
        si=13.7*exp(max(-40.,-1.92*ss))
        ELSE
        si=vjp20
        ENDIF
        goto 100
20      IF(kjp20.EQ.0)THEN
        si=0.188*exp(max(-40.,-0.611*ss))
        ELSE
        si=vjp20
        ENDIF
100     s07123=si
        return
        end

c***********************************************************************
        function s1724(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
c       pion+ + n to k+ + sigma0,channel 8
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        IF(kjp20.EQ.0)THEN
        si=0.25*(s0715(ss,ilo,ii,the)+s07123(ss,ilo,ii,the)+
     c   s1713(ss,ilo,ii,the))
        ELSE
        si=vjp20
        ENDIF
100     s1724=si
        return
        end

c*******************************************************************
        function s1727(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
c       pion+ + n to k+ + lambda,channel 9
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.684) goto 10
        IF(kjp20.EQ.0)THEN
        si=0.9*(ss-1.613)/0.091
        ELSE
        si=vjp20
        ENDIF
        goto 100
10      if(ss.ge.2.1) goto 20
        IF(kjp20.EQ.0)THEN
        si=436.3*exp(max(-40.,-4.154*ss))
        ELSE
        si=vjp20
        ENDIF
        goto 100
20      IF(kjp20.EQ.0)THEN
        si=0.314*exp(max(-40.,-0.301*ss))
        ELSE
        si=vjp20
        ENDIF
c*SA
c100    s1727=si*0.25
100     s1727=si
        return
        end

c******************************************************************
        function s1713(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
c       pion+ + p to k+ + sigma+,channel 7
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.934) goto 10
        IF(kjp20.EQ.0)THEN
        si=0.7*(ss-1.683)/0.218
        ELSE
        si=vjp20
        ENDIF
        goto 100
10      if(ss.ge.3.) goto 20
        IF(kjp20.EQ.0)THEN
        si=60.26*exp(max(-40.,-2.31*ss))
        ELSE
        si=vjp20
        ENDIF
        goto 100
20      IF(kjp20.EQ.0)THEN
        si=0.36*exp(max(-40.,-0.605*ss))
        ELSE
        si=vjp20
        ENDIF
100     s1713=si
        return
        end

c**********************************************************************
        function s2325(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
c       pion0 + n to k+ + sigma-,channel 12
        ii=i
        s2325=s1724(ss,ilo,ii,the)
        return
        end

c**********************************************************************
        function s2314(ss,ilo,i,the)
c       pion0 + p to k+ + sigma0,channel 10
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        ii=i
        s2314=s1724(ss,ilo,ii,the)
        return
        end

c**********************************************************************
        function s2317(ss,ilo,i,the)
c       pion0 + p to k+ + lambda,channel 11
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        ii=i
        s2317=s1727(ss,ilo,ii,the)
        return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine ppdelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deal with pp -> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(174).eq.0)then
        sigma1=0.
        goto 522
        endif
        the=pmas(pycomp(2224),1)+pmas(pycomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(46)
c       cross section of p + p to delta++ + n

522     if(isinel(173).eq.0)then
        sigma2=0.
        goto 523
        endif
        the=pmas(pycomp(2214),1)+pmas(pycomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(45)
c       cross section of p+p to delta+ +p
523     if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6) goto 13
        ik1=2224
        ik2=2112
        ic=174
c       p+p to delta++ +  n
        sigm12=sigma1+sigma2
        if(pyr(1).gt.sigma1/sigm12)then
        ik1=2214
        ik2=2212
        ic=173
c       cross section of p+p to delta+ +p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1
13      return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pndelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deal with pn -> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(175).eq.0)then
        sigma1=0.
        goto 524
        endif
        the=pmas(pycomp(2214),1)+pmas(pycomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(47)
c       cross section of p + n to delta+ + n

524     if(isinel(176).eq.0)then
        sigma2=0.
        goto 525
        endif
        the=pmas(pycomp(2114),1)+pmas(pycomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(48)
c       cross section of p+n to delta0 +p
525     if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)goto 13
        ik1=2214
        ik2=2112
        ic=175
c       p+n to delta+ +  n
        sigm12=sigma1+sigma2
        if(pyr(1).gt.sigma1/sigm12)then
        ik1=2114
        ik2=2212
        ic=176
c       cross section of p+n to delta0 +p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1
13      return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine nndelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deal with nn-> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(177).eq.0)then
        sigma1=0.
        goto 527
        endif
        the=pmas(pycomp(2114),1)+pmas(pycomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(49)
c       cross section of n + n to delta0 + n

527     if(isinel(178).eq.0)then
        sigma2=0.
        goto 526

        endif
        the=pmas(pycomp(1114),1)+pmas(pycomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(50)
c       cross section of n+n to delta- + p
526     if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)goto 13
        ik1=2114
        ik2=2112
        ic=177
c       n+n to delta0 + n
        sigm12=sigma1+sigma2
        if(pyr(1).gt.sigma1/sigm12)then
        ik1=1114
        ik2=2212
        ic=178
c       cross section of n+n to delta- + p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pip1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deal with pion- + p -> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(11).eq.0)then
        sigma1=0.
        goto 401
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3112),1)
        sigma1=s0715(ss,ilo1,0,the)*WEIGH(15)
c       cross section of pion- + p to k+ + sigma-
c       cross section is here in the unit of mb
401     if(isinel(12).eq.0)then
        sigma2=0.
        goto 402
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3122),1)
        sigma2=s07122(ss,ilo2,0,the)*WEIGH(7)
c       cross section of pion- + p to k0 + lambda
402     if(isinel(13).eq.0)then
        sigma3=0.
        goto 403
        endif
        the=pmas(pycomp(311),1)+pmas(pycomp(3212),1)

        sigma3=s07123(ss,ilo3,0,the)*WEIGH(16)
c       cross section of pion- + p to k0 + sigma0
403     if(isinel(147).eq.0)then
        sigma4=0.
        goto 203
        endif
c       the: threshold energy of a reaction
        the=pmas(pycomp(1114),1)+pmas(pycomp(211),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(57)
c       cross section of pion- + p to delta- + pi+
203     if(isinel(148).eq.0)then
        sigma5=0.
        goto 204
        endif
        the=pmas(pycomp(113),1)+pmas(pycomp(2112),1)
        sigma5=srho(ss,ilo5,0,the)*WEIGH(63)
c       cross section of pion- + p to rho0 + n
204     if(isinel(149).eq.0)then
        sigma6=0.
        goto 205
        endif
        the=pmas(pycomp(-213),1)+pmas(pycomp(2212),1)
        sigma6=srho(ss,ilo6,0,the)*WEIGH(64)
c       cross section of pion- + p to rho- + p
205     if(isinel(150).eq.0)then
        sigma7=0.
        goto 310
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(2214),1)
        sigma7=sdelta(ss,ilo7,0,the)*WEIGH(58)
c       cross section of pion- + p to delta+ + pion-
310     if(isinel(151).eq.0)then
        sigma8=0.
        goto 805
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(2114),1)
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(58)
c       cross section of pion- + p to delta0 + pion0

805     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &   .and.ilo4.eq.0 .and. ilo5.eq.0 .and. ilo6.eq.0
     &   .and.ilo7.eq.0.and.ilo8.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &         .and.sigma4.lt.1.e-6 .and. sigma5.lt.1.e-6 .and.
     &          sigma6.lt.1.e-6.and. sigma7.lt.1.e-6
     &         .and.sigma8.lt.1.e-6)goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3112
c       3112 is the flavor code of sigma-
        ic=11
c       pion- + p to k+ + sigma-
        goto 416
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=311
        ik2=3122
c       3122 is the Kf code of lambda
        ic=12
c       pion- + p to k0 + lambda
        goto 416
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=311
        ik2=3212
c       3212 is the KF code of sigma0
        ic=13
c       pion- + p to k0 + sigma0
        goto 416
        endif
        if(rlus.gt.s3 .and. rlus.le.s4)then
        ik1=211
        ik2=1114
c       1114 is the Kf code of delta-
        ic=147
c       pion- + p to pi+ + delta-
        goto 416
        endif
        if(rlus.gt.s4 .and. rlus.le.s5)then
        ik1=113
        ik2=2112
c       113 is the Kf code of rho0
        ic=148
c       pion- + p to rho0 + n
        goto 416
        endif
        if(rlus.gt.s5 .and. rlus.le.s6)then
        ik1=-213
        ik2=2212
c       -213 is the Kf code of rho-
        ic=149
c       pion- + p to rho- + p
        goto 416
        endif
        if(rlus.gt.s6 .and. rlus.le.s7)then
        ik1=-211
        ik2=2214
c       2214 is the Kf code of delta+
        ic=150
c       pion- + p to delta++pion-
        goto 416
        endif
        ik1=111
        ik2=2114
c       2114 is the Kf code of delta0
        ic=151
c       pion- + p to delta0+pion0
        goto 416
416     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pin3(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deat with pion- + n -> ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(152).eq.0)then
        sigma1=0.
        goto 222
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(1114),1)
        sigma1=sdelta(ss,ilo1,0,the)*WEIGH(61)
c       cross section of pion- + n to delta- + pi0

222     if(isinel(153).eq.0)then
        sigma2=0.
        goto 223
        endif
        the=pmas(pycomp(-213),1)+pmas(pycomp(2112),1)
        sigma2=srho(ss,ilo2,0,the)*WEIGH(65)
c       cross section of pion- + n to rho- + n
223     if(isinel(154).eq.0)then
        sigma3=0.
        goto 228
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(2114),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(62)
c       cross section of pion- + n to delta0 + pion-
228     if(isinel(14).eq.0)then
        sigma4=0.
        goto 818
        endif
        ilo4=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3112),1)

        sigma4=s1724(ss,ilo4,0,the)*WEIGH(24)

c       cross section of pion- + n to k0 + sigma-
818     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0
     &   .and.ilo4.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6.and.
     &   sigma3.lt.1.e-6.and. sigma4.lt.1.e-6)goto 13
        sigm12=sigma1+sigma2
        sigm13=sigm12+sigma3
        sigm14=sigm13+sigma4
        s1=sigma1/sigm14
        s2=sigm12/sigm14
        s3=sigm13/sigm14
        rlu1=pyr(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=1114
        ic=152
c       pion- + n to delta-  +  pi0
        elseif(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-213
        ik2=2112
        ic=153
c       cross section of pion- + n to rho- + n
        elseif(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-211
        ik2=2114
c       cross section of pion- + n to delta0 + pion-
        ic=154
        else
        ik1=311
        ik2=3112
        ic=14
c       cross section of pion- + n to k0 + sigma-
        endif
810     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm14/cspin/10.
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pin2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deal with pion0 + n to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(172).eq.0)then
        sigma1=0.
        goto 216
        endif
        ilo1=1
        the=pmas(pycomp(-213),1)+pmas(pycomp(2212),1)
        sigma1=srho(ss,ilo1,0,the)*WEIGH(72)
c       cross section of pion0 + n to rho- + p
216     if(isinel(168).eq.0)then
        sigma2=0.
        goto 217
        endif
c       the--threshold energy of a reaction
        the=pmas(pycomp(-211),1)+pmas(pycomp(2214),1)
        sigma2=sdelta(ss,ilo2,0,the)*WEIGH(59)
c       cross section of pion0 + n to delta+ + pi-
217     if(isinel(171).eq.0)then
        sigma3=0.
        goto 218
        endif
        the=pmas(pycomp(113),1)+pmas(pycomp(2112),1)
        sigma3=srho(ss,ilo3,0,the)*WEIGH(70)
c       cross section of pion0 + n to rho0 + n
218     if(isinel(169).eq.0)then
        sigma4=0.
        goto 807
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(1114),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(60)
c       cross section of pion0 + n to delta- + pion+
807     if(isinel(19).eq.0)then
        sigma5=0.
        goto 811
        endif
        ilo5=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3122),1)

        sigma5=s1724(ss,ilo5,0,the)*WEIGH(10)

c       cross section of pion0 + n to k0 + lambda is assumed
c       to be cross section of pion0 + p to k+ + lambda
811     if(isinel(18).eq.0)then
        sigma6=0.
        goto 816
        endif
        ilo6=1
        the=pmas(pycomp(321),1)+pmas(pycomp(3112),1)

        sigma6=s2325(ss,ilo6,0,the)*WEIGH(22)

c       cross section of pion0 + n to k+ + sigma-
816     if(isinel(20).eq.0)then
        sigma7=0.
        goto 817
        endif
        ilo7=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3212),1)

        sigma7=s1724(ss,ilo7,0,the)*WEIGH(23)

c       cross section of pion0 + n to k0 + sigma0

817     if(isinel(170).eq.0)then
        sigma8=0.
        goto 827
        endif
        the=pmas(pycomp(111),1)+pmas(pycomp(2114),1)
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(60)
c       cross section of pion0 + n to delta0 + pion0
827     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &   ilo4.eq.0.and. ilo5.eq.0
     &   .and. ilo6.eq.0 .and. ilo7.eq.0 .and. ilo8.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &          .and. sigma4.lt.1.e-6.and. sigma5.lt.1.e-6
     &          .and. sigma6.lt.1.e-6.and. sigma7.lt.1.e-6
     &          .and. sigma8.lt.1.e-6)goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=-213
        ik2=2212
c       -213 is the flavor code of rho-
        ic=172
c       pion0 + n to p + rho-
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-211
        ik2=2214
c       2214 is the flavor code of delta+
        ic=168
c       pion0 + n to delta+ + pi-
        elseif(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=113
        ik2=2112
c       113 is the Kf code of rho
        ic=171
c       pion0 + n to rho  +  n
        elseif(rlus.gt.s3 .and. rlus.le.s4)then
        ik1=211
        ik2=1114
c       1114 is the Kf code of delta-
        ic=169
        elseif(rlus.gt.s4 .and. rlus.le.s5)then
        ik1=311
        ik2=3122
c       3122 is the Kf code of lambda
c       pion0 + n to k0 + lambda
        ic=19
        elseif(rlus.gt.s5 .and. rlus.le.s6)then
        ik1=321
        ik2=3112
c       3112 is the Kf code of sigma-
        ic=18
        elseif(rlus.gt.s6 .and. rlus.le.s7)then
c       pion0 + n to k+ + sigma-
        ik1=311
        ik2=3212
c       3212 is the Kf code of sigma0
        ic=20
c       pion0 + n to k0 + sigma0
        else
        ik1=111
        ik2=2114
c       2114 is the Kf code of delta0
        ic=170
c       pion0 + n to pi0 + delta0
        endif
219     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pin1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deal with pion+ + n to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(8).eq.0)then
        sigma1=0.
        goto 306
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3212),1)
        sigma1=s1724(ss,ilo1,0,the)*WEIGH(18)
c       cross section of pion+ + n to k+ + sigma0
306     if(isinel(9).eq.0)then
        sigma2=0.
        goto 207
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3122),1)
        sigma2=s1727(ss,ilo2,0,the)*WEIGH(8)
c       cross section of pion+ + n to k+ + lambda
207     if(isinel(158).eq.0)then
        sigma3=0.
        goto 208
        endif
c       the: threshold energy of a reaction
        the=pmas(pycomp(-211),1)+pmas(pycomp(2224),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(55)
c       cross section of pion+ + n to delta++ + pi-
208     if(isinel(161).eq.0)then
        sigma4=0.
        goto 209
        endif
        the=pmas(pycomp(113),1)+pmas(pycomp(2212),1)
        sigma4=srho(ss,ilo4,0,the)*WEIGH(66)
c       cross section of pion+ + n to rho0 + p
209     if(isinel(162).eq.0)then
        sigma5=0.
        goto 210
        endif
        the=pmas(pycomp(213),1)+pmas(pycomp(2112),1)
        sigma5=srho(ss,ilo5,0,the)*WEIGH(67)
c       cross section of pion+ + n to rho+ + n
210     if(isinel(159).eq.0)then
        sigma6=0.
        goto 806
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(2114),1)
        sigma6=sdelta(ss,ilo6,0,the)*WEIGH(56)
c       cross section of pion+ + n to delta0 + pion+
806     if(isinel(10).eq.0)then
        sigma7=0.
        goto 814
        endif
        ilo7=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3222),1)

        sigma7=s1724(ss,ilo7,0,the)*WEIGH(19)

c       cross section of pion+ + n to k0 + sigma+

814     if(isinel(160).eq.0)then
        sigma8=0.
        goto 818
        endif
        ilo8=1
        the=pmas(pycomp(111),1)+pmas(pycomp(2214),1)

        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(19)

c       cross section of pion+ + n to pi0 + delta+

818     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &   .and.ilo4.eq.0 .and. ilo5.eq.0.and.ilo6.eq.0
     &   .and.ilo7.eq.0.and.ilo8.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &         .and.sigma4.lt.1.e-6 .and. sigma5.lt.1.e-6.and. 
     &          sigma6.lt.1.e-6
     &         .and.sigma7.lt.1.e-6.and.sigma8.lt.1.e-6)goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3212
c       3212 is the flavor code of sigma0
        ic=8
c       pion+ + n to k+ + sigma0
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=321
        ik2=3122
c       3122 is the flavor code of lambda
        ic=9
c       pion+ + n to k+ + lambda
        elseif(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=-211
        ik2=2224
c       2224 is the Kf code of delta++
        ic=158
c       pion+ + n to pi- + delta++
        elseif(rlus.gt.s3 .and. rlus.le.s4)then
        ik1=113
        ik2=2212
c       113 is the Kf code of rho0
        ic=161
c       pion+ + n to rho0 + p
        elseif(rlus.gt.s4 .and. rlus.le.s5)then
        ik1=213
        ik2=2112
c       213 is the Kf code of rho+
        ic=162
c       pion+ + n to rho+ + n
        elseif(rlus.gt.s5 .and. rlus.le.s6)then
        ik1=211
        ik2=2114
c       2114 is the Kf code of delta0
        ic=159
        elseif(rlus.gt.s6 .and. rlus.le.s7)then
        ik1=311
        ik2=3222
c       3222 is the Kf code of sigma+
        ic=10
c       pion+ + n to k0 + sigma+
        else
        ik1=111
        ik2=2214
c       2214 is the Kf code of delta+
        ic=160
c       pion+ + n to pi0 + delta+
        endif
211     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1
13      return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pip2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deal with pion+ + p to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(7).eq.0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3222),1)
        sigma1=s1713(ss,ilo1,0,the)*WEIGH(17)
c       cross section of pion+ + p to k+ + sigma+
212     if(isinel(155).eq.0)then
        sigma2=0.
        goto 213
        endif
c       the: threshold energy of a reaction
        the=pmas(pycomp(111),1)+pmas(pycomp(2224),1)
        sigma2=sdelta(ss,ilo2,0,the)*WEIGH(51)
c       cross section of pion+ + p to delta++ + pi0
213     if(isinel(157).eq.0)then
        sigma3=0.
        goto 214
        endif
        the=pmas(pycomp(213),1)+pmas(pycomp(2212),1)
        sigma3=srho(ss,ilo3,0,the)*WEIGH(68)
c       cross section of pion+ + p to rho+ + p
214     if(isinel(156).eq.0)then
        sigma4=0.
        goto 804
        endif
        the=pmas(pycomp(211),1)+pmas(pycomp(2214),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(52)
c       cross section of pion+ + p to delta+ + pion+
804     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &   .and. ilo4.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &   .and. sigma4.lt.1.e-6 )goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        s1=sigma1/sigma14
        s2=sigma12/sigma14
        s3=sigma13/sigma14
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3222
c       3222 is the flavor code of sigma+
        ic=7
c       pion+ + p to k+ + sigma+
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=111
        ik2=2224
c       2224 is the flavor code of delta++
        ic=155
c       pion+ + p to delta++  +  pi0
        elseif(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=213
        ik2=2212
c       213 is the Kf code of rho+
        ic=157
c       pion+ + p to rho+ + p
        else
        ik1=211
        ik2=2214
c       2214 is the Kf code of delta+
        ic=156
c       pion+ + p to delta+ + pion+
        endif
215     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma14/cspin/10.
        ioo=1
13      return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        function sdelta(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
c       pion- + p to pi0 + delta0
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.6941) goto 10
        IF(kjp20.EQ.0)THEN
        si=-61.127+42.9365*ss
        ELSE
        si=vjp21
        ENDIF
        goto 100
10      IF(kjp20.EQ.0)THEN
        sit=-0.0186959*ss**3+0.310359*ss**2-0.755106*ss+0.565481
        si=1.0/sit
        ELSE
        si=vjp21
        ENDIF
100     sdelta=si
        return
        end

c*******************************************************************
        function srho(ss,ilo,i,the)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
c       pion- + p to pi0 + delta0
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.8837) goto 10
        IF(kjp20.EQ.0)THEN
        si=-23.3607+13.9936*ss
        ELSE
        si=vjp22
        ENDIF
        goto 100
10      IF(kjp20.EQ.0)THEN
        sit=0.331583*ss**3-1.86123*ss**2+3.81364*ss-2.50068
        si=1.0/sit
        ELSE
        si=vjp22
        ENDIF
100     srho=si
        return
        end

c*********************************************************************
        function snn(ss,ilo,i,the)
c       n + n to n + delta
c       parameterized x-section is not correct for n-n
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        IF(kjp20.EQ.0)THEN
        si=20*(ss-2.015)**2/(0.015+(ss-2.015)**2)
        ELSE
        si=vjp23
        ENDIF
100     snn=si
        return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pip3(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'prod' to deal with pion0 + p to ...
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(15).eq.0)then
        sigma1=0.
        goto 308
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3212),1)
        sigma1=s2314(ss,ilo1,0,the)*WEIGH(20)
c       cross section of pion0 + p to k+ + sigma0
308     if(isinel(16).eq.0)then
        sigma2=0.
        goto 239
        endif
        the=pmas(pycomp(321),1)+pmas(pycomp(3122),1)
        sigma2=s2317(ss,ilo2,0,the)*WEIGH(9)
c       cross section of pion0 + p to k+ + lambda
239     if(isinel(163).eq.0)then
        sigma3=0.
        goto 220
        endif
c       the: threshold energy of a reaction
        the=pmas(pycomp(211),1)+pmas(pycomp(2114),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(53)
c       cross section of pion0 + p to delta0 + pi+
220     if(isinel(165).eq.0)then
        sigma4=0.
        goto 240
        endif
        the=pmas(pycomp(213),1)+pmas(pycomp(2112),1)
        sigma4=srho(ss,ilo4,0,the)*WEIGH(69)
c       cross section of pion0 + p to rho+ + n
240     if(isinel(164).eq.0)then
        sigma5=0.
        goto 808
        endif
        the=pmas(pycomp(-211),1)+pmas(pycomp(2224),1)
        sigma5=sdelta(ss,ilo5,0,the)*WEIGH(54)
c       cross section of pion0 + p to delta++ + pion-
808     if(isinel(17).eq.0)then
        sigma6=0.
        goto 815
        endif
        ilo6=1
        the=pmas(pycomp(311),1)+pmas(pycomp(3222),1)

        sigma6=s1724(ss,ilo6,0,the)*WEIGH(21)

c       cross section of pion0 + P to k0 + sigma+

815     if(isinel(166).eq.0)then
        sigma7=0.
        goto 835
        endif
        ilo7=1
        the=pmas(pycomp(113),1)+pmas(pycomp(2212),1)
        sigma7=srho(ss,ilo7,0,the)*WEIGH(71)
c       cross section of pion0 + P to rho0 + p
835     if(isinel(167).eq.0)then
        sigma8=0.
        goto 845
        endif
        ilo8=1
        the=pmas(pycomp(111),1)+pmas(pycomp(2214),1)
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(54)
c       cross section of pion0 + P to pi0 + delta+
845     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &   .and.ilo4.eq.0.and.ilo5.eq.0.and.ilo6.eq.0
     &   .and.ilo7.eq.0.and.ilo8.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &   .and.sigma4.lt.1.e-6.and.sigma5.lt.1.e-6.and.sigma6.lt.1.e-6
     &   .and.sigma7.lt.1.e-6.and.sigma8.lt.1.e-6 )goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3212
        ic=15
c       3212 is the flavor code of sigma0
c       pion0 + p to k+ + sigma0
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=321
        ik2=3122
        ic=16
c       3122 is the flavor code of lambda
c       pion0 + p to k+ + lambda
        elseif(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=211
        ik2=2114
c       2214 is the Kf code of delta+
        ic=163
c       pion0 + p to pi+ + delta0
        elseif(rlus.gt.s3 .and. rlus.le.s4)then
        ik1=213
        ik2=2112
c       213 is the Kf code of rho+
        ic=165
c       pion0 + p to rho+ + n
        elseif(rlus.gt.s4 .and. rlus.le.s5)then
        ik1=-211
        ik2=2224
c       2224 is the Kf code of delta++
        ic=164
        elseif(rlus.gt.s5 .and. rlus.le.s6)then
        ik1=311
        ik2=3222
c       3222 is the Kf code of sigma+
c       pion0 + p to k0 + sigma+
        ic=17
        elseif(rlus.gt.s6 .and. rlus.le.s7)then
        ik1=113
        ik2=2212
c       113 is the Kf code of rho0
        ic=166
c       pion0 + p to rho0 + p
        else
        ik1=111
        ik2=2214
c       2214 is the Kf code of delta+
c       pion0 + p to pi0 + delta+
        ic=167
        endif
221     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine jpsin(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,ii)
c       a part of 'prod' to deal with J/Psi (psi') + n to ... 
c       ii = 1 for J/Psi, ii = 2 for Psi'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(179).eq.0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(4122),1)+pmas(pycomp(411),1)
c       the--threshold energy of a reaction (squart root of s)
        if(ss.le.the)then
        sigma1=0.
        goto 212
        endif
        sigma1=cspsn
c       cross section of J/Psi + n (fm^2) 
        if(ii.eq.2)sigma1=csspn
c       cross section of Psi' + n (fm^2) 
212     if(isinel(180).eq.0)then
        sigma2=0.
        goto 213
        endif
        the=pmas(pycomp(4212),1)+pmas(pycomp(411),1)
        if(ss.le.the)then
        sigma2=0.
        goto 213
        endif
        sigma2=cspsn
        if(ii.eq.2)sigma2=csspn
213     if(isinel(181).eq.0)then
        sigma3=0.
        goto 214
        endif
        the=pmas(pycomp(4112),1)+pmas(pycomp(421),1)
        if(ss.le.the)then
        sigma3=0.
        goto 214
        endif
        sigma3=cspsn
        if(ii.eq.2)sigma3=csspn
214     if(sigma1.lt.1.e-16 .and. sigma2.lt.1.e-16 .and. 
     &    sigma3.lt.1.e-16)goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        s1=sigma1/sigma13
        s2=sigma12/sigma13
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=4122
        ik2=-411
        ic=179
c       J/Psi (Psi') + n to lamdac + Dba 
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=4212
        ik2=-411
        ic=180
c       J/Psi (Psi') + n to sigmac+  + Dba
        else
        ik1=4112
        ik2=-421
        ic=181
c       J/Psi (Psi') + n to sigmac0 + D0ba
        endif
215     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma13*rcsit/cspsn
        if(ii.eq.2)tw(icp)=sigma13*rcsit/csspn
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine jpsip(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,ii)
c       a part of 'prod' to deal with J/Psi (Psi') + p to ... 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(182).eq.0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(4122),1)+pmas(pycomp(421),1)
c       the: threshold energy of a reaction (squart root of s)
        if(ss.le.the)then
        sigma1=0.
        goto 212
        endif
        sigma1=cspsn
        if(ii.eq.2)sigma1=csspn
212     if(isinel(183).eq.0)then
        sigma2=0.
        goto 213
        endif
        the=pmas(pycomp(4212),1)+pmas(pycomp(421),1)
        if(ss.le.the)then
        sigma2=0.
        goto 213
        endif
        sigma2=cspsn
        if(ii.eq.2)sigma2=csspn
213     if(isinel(184).eq.0)then
        sigma3=0.
        goto 214
        endif
        the=pmas(pycomp(4222),1)+pmas(pycomp(411),1)
        if(ss.le.the)then
        sigma3=0.
        goto 214
        endif
        sigma3=cspsn
        if(ii.eq.2)sigma3=csspn
214     if(sigma1.lt.1.e-16 .and. sigma2.lt.1.e-16 .and. 
     &    sigma3.lt.1.e-16)goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        s1=sigma1/sigma13
        s2=sigma12/sigma13
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=4122
        ik2=-421
        ic=182
c       J/Psi (psi') + p to lamdac + D0ba 
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=4212
        ik2=-421
        ic=183
c       J/Psi (psi') + p to sigmac+  + D0ba
        else
        ik1=4222
        ik2=-411
        ic=184
c       J/Psi (psi') + p to sigmac++ + Dba
        endif
215     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma13*rcsit/cspsn
        if(ii.eq.2)tw(icp)=sigma13*rcsit/csspn
        ioo=1
13      return
        end


c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine jpsip1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,ii)
c       a part of 'prod' to deal with J/Psi (Psi') + pion+ to ... 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(185).eq.0)goto 13
        the=pmas(pycomp(411),1)+pmas(pycomp(423),1)
        if(ss.le.the)goto 13
c98/03/23   ww=cspsm
c       J/Psi + pion+ to D + D*0ba
        lc(icp,3)=411
        lc(icp,4)=-423
        lc(icp,5)=185
        tw(icp)=rcsit
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine jpsip0(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,ii)
c       a part of 'prod' to deal with J/Psi (Psi') + pion0 to ... 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(186).eq.0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(421),1)+pmas(pycomp(423),1)
c       the: threshold energy of a reaction (squart root of s)
        if(ss.le.the)then
        sigma1=0.
        goto 212
        endif
        sigma1=cspsm   
c       cross section of J/Psi + m (fm^2) 
        if(ii.eq.2)sigma1=csspm
c       cross section of Psi' + m (fm^2) 
212     if(isinel(187).eq.0)then
        sigma2=0.
        goto 213
        endif
        the=pmas(pycomp(411),1)+pmas(pycomp(413),1) 
        if(ss.le.the)then
        sigma2=0.
        goto 213
        endif
        sigma2=cspsm
        if(ii.eq.2)sigma2=csspm
213     if(sigma1.lt.1.e-16 .and. sigma2.lt.1.e-16)goto 13
        sigma12=sigma1+sigma2
        s1=sigma1/sigma12
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=421
        ik2=-423
        ic=186
c       J/Psi (Psi') + pion0 to D0 + D*0ba 
        else
        ik1=411
        ik2=-413
        ic=187
c       J/Psi (Psi') + pion0 to D + D*ba
        endif
215     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma12*rcsit/cspsm
        if(ii.eq.2)tw(icp)=sigma12*rcsit/csspm
        ioo=1
13      return
        end


c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine jpsip2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,ii)
c       a part of 'prod' to deal with J/Psi (Psi') + pion- to ... 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(188).eq.0)goto 13
        the=pmas(pycomp(421),1)+pmas(pycomp(413),1)
        if(ss.le.the)goto 13
c98/03/23   ww=cspsm
c       J/Psi (Psi') + pion- to D0 + D*ba
        lc(icp,3)=421
        lc(icp,4)=-413
        lc(icp,5)=188
        tw(icp)=rcsit
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine jpsir1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,ii)
c       a part of 'prod' to deal with J/Psi (Psi') + rho+ to ... 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(189).eq.0)goto 13
        the=pmas(pycomp(411),1)+pmas(pycomp(421),1)
        if(ss.le.the)goto 13
c98/03/23   ww=cspsm
c       J/Psi (Psi,) + rho+ to D + D0ba
        lc(icp,3)=411
        lc(icp,4)=-421
        lc(icp,5)=189
        tw(icp)=rcsit
        ioo=1
13      return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine jpsir0(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,ii)
c       a part of 'prod' to deal with J/Psi (Psi') + rho0 to ... 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(190).eq.0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(pycomp(421),1)*2.
c       the--threshold energy of a reaction (squart root of s)
        if(ss.le.the)then
        sigma1=0.
        goto 212
        endif
        sigma1=cspsm
        if(ii.eq.2)sigma1=csspm
212     if(isinel(191).eq.0)then
        sigma2=0.
        goto 213
        endif
        the=pmas(pycomp(411),1)*2.
        if(ss.le.the)then
        sigma2=0.
        goto 213
        endif
        sigma2=cspsm
        if(ii.eq.2)sigma2=csspm
213     if(sigma1.lt.1.e-16 .and. sigma2.lt.1.e-16)goto 13
        sigma12=sigma1+sigma2
        s1=sigma1/sigma12
        rlus=pyr(1)
        if(rlus.le.s1)then
        ik1=421
        ik2=-421
        ic=190
c       J/Psi (Psi') + rho0 to D0 + D0ba 
        else
        ik1=411
        ik2=-411
        ic=191
c       J/Psi (Psi') + rho0 to D + Dba
        endif
215     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=2.*sigma12*rcsit/cspsm
        if(ii.eq.2)tw(icp)=2.*sigma12*rcsit/csspm
c       at the case of rho0, cross section enlarges a factor 2 to
c        consider the effect of omega
        ioo=1
13      return
        end


c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine jpsir2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo,ii)
c       a part of 'prod' to deal with J/Psi (Psi') + rho- to ... 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(NSIZE=750000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYCIDAT2/KFMAXT,NONT2,PARAM(20),WEIGH(600)
        SAVE /PYCIDAT2/
        common/sa10_h/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,csspn,csspm
        common/count_h/isinel(600)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        if(isinel(192).eq.0)goto 13
        the=pmas(pycomp(421),1)+pmas(pycomp(411),1)
        if(ss.le.the)goto 13
c98/03/23   ww=cspsm
c       J/Psi + rho- to D0 + Dba
        lc(icp,3)=421
        lc(icp,4)=-411
        lc(icp,5)=192
        tw(icp)=rcsit
        ioo=1
13      return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine decay(time,deltt,lc,tc,tw,iii,iiii)
c       deal with particle decay in transport processes
c       iii: number of loop within a event
c       iiii: number of event
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa17/nde,non17,kde(10,2),pde(10,5),vde(10,4)
        common/sa19_h/coor(3)
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension pp(5),rr(4)
c       idec: stord the line number (in the particle list) of 
c        particles after decay
c       iwide: stord the line number (in the particle list) of decaying
c        particles
c       the messages of decayed particles are
c        stored in the varibles and arraies in 'sa17'
c       pp, rr: momentum, position of decaying particle
        nde=0
        do i=1,10
        do j=1,2
        kde(i,j)=0
        enddo
        do j=1,5
        pde(i,j)=0.
        enddo
        do j=1,4
        vde(i,j)=0.
        enddo
        enddo
        do i=1,2000
        iwide(i)=0
        enddo
        ii=0
        do i=1,nsa
        kf=ksa(i,2)
        if(abs(kf).eq.213 .or. kf.eq.113)then
        ii=ii+1
        iwide(ii)=i
        endif
        enddo
        if(ii.eq.0)return
        do 100 i=1,ii
c       deal with decay of decaying particles one by one
        jj=iwide(i)
        kf=ksa(jj,2)
        do i1=1,5
        pp(i1)=psa(jj,i1)
        enddo
        do i1=1,3
        rr(i1)=vsa(jj,i1)
        enddo
        rr(4)=time
        ee=pp(4)
        if(ee.lt.1.e-15)ee=1.e-15
        amass=pp(5)
c010600 compo=amass*0.151*deltt/0.197/ee
c010600
        tauj=tau(jj)
        tdel=deltt
        if((time-deltt).lt.tauj)tdel=time-tauj
        compo=amass*0.151*tdel/0.197/ee
c010600
        if(compo.lt.1.e-15)compo=1.e-15
        prob=1.-exp(-compo)
c       0.151: full width in [GeV] of decaying particle (rho)
c       mass, full width, energy: in [GeV], time in [GeV^-1] originally
c       0.197: conversion factor from [GeV^-1] to [fm/c]
        if(pyr(1).le.prob)then
c       perform a particle decay
        call decpr(pp,rr,kf)
c       it gives momentum and position to the particles after decay 
        call upddep(jj,lc,tc,tw,time,iii)
c       update particle list after a particle decay 
        call upddet(jj,lc,tc,tw,time,iii)
c       update collision time list after a particle decay

c       if(nctl.eq.0)return

        endif
100     continue
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine decpr(pp,rr,kf)
c       give momentum and position to decayed particles
c       pp, rr: momentum, position of decaying particle
c       kf: flavour code of decaying particle
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa17/nde,non17,kde(10,2),pde(10,5),vde(10,4)
        dimension pp(5),rr(4),pi(4),pj(4),rd(3),b(3)
        amas=pp(5)
        amas2=amas*amas
        nde=2
c       assume two body decay

        if(kf.eq.213)then
c       rho+ to pion+ + pion0
        am1=pmas(pycomp(211),1)
        am2=pmas(pycomp(111),1)
        kde(1,2)=211
        kde(2,2)=111
        endif
        if(kf.eq.-213)then
c       rho- to pion- + pion0
        am1=pmas(pycomp(-211),1)
        am2=pmas(pycomp(111),1)
        kde(1,2)=-211
        kde(2,2)=111
        endif
        if(kf.eq.113)then
c       rho0 to pion0 + pion0
        am1=pmas(pycomp(111),1)
        am2=am1
        kde(1,2)=111
        kde(2,2)=111
        endif

        agum=(amas2-(am1+am2)**2)*(amas2-(am1-am2)**2)
        if(agum.le.0.)then
        write(9,*)'kf,kf1,kf2=',kf,kde(1,2),kde(2,2)
        write(9,*)'amas,am1,am2=',amas,am1,am2
        agum=1.e-15
        endif

        pab=sqrt(agum)/2./amas
        cita=2.*pyr(1)-1.
        if(cita.gt.0.9999)cita=0.9999
        sita=sqrt(1.-cita*cita)
        fi=2.*3.1416*pyr(1)
        px=pab*sita*cos(fi)
        py=pab*sita*sin(fi)
        pz=pab*cita
        pi(1)=px
        pi(2)=py
        pi(3)=pz
        pi(4)=sqrt(pab*pab+am1*am1)
        pj(1)=-px
        pj(2)=-py
        pj(3)=-pz
        pj(4)=sqrt(pab*pab+am2*am2)
        ee=pp(4)
        do i=1,3
        b(i)=pp(i)/ee
        enddo
        ilo=1
        call lorntz(ilo,b,pi,pj)
        do i=1,4
        pde(1,i)=pi(i)
        pde(2,i)=pj(i)
        enddo
        pde(1,5)=am1
        pde(2,5)=am2

        rrp=0.8
c       rrp: the radius (in fm) of decaying particle 
        do i=1,nde
        cita=2.*pyr(1)-1.
        if(cita.gt.0.9999)cita=0.9999
        sita=sqrt(1.-cita*cita)
        fi=2.*3.1416*pyr(1)
        rd(1)=rrp*sita*cos(fi)
        rd(2)=rrp*sita*sin(fi)
        rd(3)=rrp*cita
c       arrange the particle i after decay on the surface
c        of sphere with radius rrp and centered at the position of parent
        do j=1,3
        vde(i,j)=rd(j)+rr(j)
        enddo
        vde(i,4)=rr(4)
        enddo
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine upddep(ij,lc,tc,tw,time,iii)
c       update particle list after a particle decay and
c        truncate collision list correspondingly.
c       ij: line number (in particle list 'sa1_h') of decaying particle
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
        common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa8_h/tau(kszj),ishp(kszj)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa17/n,non17,k(10,2),p(10,5),v(10,4)
        common/sa19_h/coor(3)
        common/sa20_h/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,ecsspn,ecsspm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        dimension lc(nsize,5),tc(nsize),tw(nsize)
c       idec: stord the order number (in the particle list) of 
c        particles after decay
c        dimension peo(4)
        do m=1,2000
        idec(m)=0
        enddo
c       put decayed particles (in 'sa17') into particle list 'sa1_h'
c        and turncate collision list correspondingly.
        ll=ij
        kd=ksa(ll,2)
        do 500 i=1,n
        kf=k(i,2)
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
c       give proper values to particle jj.
        ksa(jj,2)=kf
        ksa(jj,1)=1
        ksa(jj,3)=0
        do m=1,5
        psa(jj,m)=p(i,m)
        vsa(jj,m)=v(i,m)
        enddo
        ishp(jj)=1
        tau(jj)=time
c       the values of 'ishp' and 'tau' of decayed particles 
c        are given here, its formation time are assumed to be zero 
        do m=j,kfmax
        numb(m)=numb(m)+1
        enddo
        idec(i)=jj
        do m=1,2000
        ipym=idec(m)
        iwi=iwide(m)
        if(ipym.gt.jj)idec(m)=ipym+1
        if(iwi.gt.jj)iwide(m)=iwi+1
        enddo
c       update the values of lc(m,1-2) with value.ge.jj
        do m=1,nctl
        lc1=lc(m,1)
        if(lc1.ge.jj)lc(m,1)=lc1+1
        lc2=lc(m,2)
        if(lc2.ge.jj)lc(m,2)=lc2+1
        enddo
1000    goto 200
600     continue
200     continue
        nsa=nsa+1
500     continue

c       through away the colli. pairs with partner of ll 
c        (decaying particle)
        j=0
        do i=1,nctl
        i1=lc(i,1)
        j1=lc(i,2)
        if(i1.eq.ll)goto 401
        if(j1.eq.ll)goto 401
        if((tc(i)-time).le.ddt) goto 401
c       through away the pairs with tc<= time
        j=j+1
        tc(j)=tc(i)
        tw(j)=tw(i)
        do m=1,5
        lc(j,m)=lc(i,m)
        enddo
401     continue
        enddo
        do i=j+1,nctl+1
        tc(i)=0.0
        tw(i)=0.0
        do m=1,5
        lc(i,m)=0
        enddo
        enddo
        nctl=j

c       remove decaying particle (ll) from particle list and 
c        truncate the collision list correspondingly.
        kf=ksa(ll,2)
cc      kf=kf1
cc      ll=l
        if(ll.eq.nsa)then
        kfd=numb(kfmax)-numb(kfmax-1)

c       following statements are added by sa on 07/March/97
        if(kfd.eq.0)then
        numbm=numb(kfmax)
        do i1=1,kfmax
        if(numb(i1).ne.numbm)goto 3000
        i2=i1
        goto 3001
3000    enddo
3001    do i1=i2,kfmax
        numb(i1)=numb(i1)-1
        enddo
        goto 400
        endif
        numb(kfmax)=numb(kfmax)-1
400     goto 100
        endif
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
        do m=1,2000
        ipym=idec(m)
        iwi=iwide(m)
        if(ipym.gt.ll)idec(m)=ipym-1
        if(iwi.gt.ll)iwide(m)=iwi-1
        enddo
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine upddet(ij,lc,tc,tw,time,iii)
c       update collision list after a particle decay
c       ij: oredr number (in particle list) of decaying particle
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000,NSIZE=750000)
c       common/sa1_h/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa9_h/kfmax,kfaco(100),numb(100),non9,disbe(100,100)
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa17/n,non17,k(10,2),p(10,5),v(10,4)
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 141104
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        dimension lc(nsize,5),tc(nsize),tw(nsize)
c       idec: store line number (in the particle list) of decayed particles 
        dddt=adj1(11)   ! 141104

        nctl=nctl+1
c       loop over particle list

        m2=numb(2)
        m4=numb(4)
        m7=numb(7)
        m9=numb(9)
        m17=numb(17)
        m19=numb(19)
        m25=numb(25)
        m29=numb(29)
        m32=numb(32)
        m34=numb(34)
c        m34=numb(kfmax-11)
c       subtract 11, since we do not consider the rescattering of x0c, etc

        j21=idec(1)
        j22=idec(2)
        do j11=1,n
        j1=idec(j11)
        if(j1.le.m19
     c   .or. (j1.gt.m25 .and. j1.le.m34))goto 301
c       consider only the reinteraction among nucleon, pion, kaon,
c        sigma, lambda, delta, rho and psi
        goto 300

301     do i=1,m34

        if(j11.eq.1 .and. i.eq.j22)goto 600
        if(j11.eq.2 .and. i.eq.j21)goto 600
c       decayed particles could not collide with each other immediately
        if(nctl.gt.nsize)then
        write(9,*)'3 nsa,nctl,dddt=',nsa,nctl,dddt
        write(9,*)'size of array "nsize" needs to be extended'
        write(9,*)'error is serious,stop running'
        stop 30000
        endif

        if(i.le.m19 
     c   .or. (i.gt.m25 .and. i.le.m34))goto 602
c       consider only the reinteraction among nucleon, pion, kaon,
c        sigma, lambda, delta, rho and psi
        goto 600

c       if((j1.le.m2 .or. (j1.gt.m4 .and. j1.le.m7)) .and.
c     c  (i.gt.m32 .and. i.le.m34))goto 602
c       if((i.le.m2 .or. (i.gt.m4 .and. i.le.m7)) .and.
c     c  (j1.gt.m32 .and. j1.le.m34))goto 602
c       goto 600
c       consider only the reinteraction between (nucleon, pion) and psi

csa601     if(j1.le.m2 .and. i.le.m2.and.j1.lt.i)goto 600
c       no double counting for NN
csa        if((j1.gt.m2 .and. j1.le.m4) .and. i.gt.m25) goto 600
csa     if((i.gt.m2 .and. i.le.m4) .and. j1.gt.m25) goto 600
c       delta, rho, psi etc. do not scatter with Nbar
csa        if((j1.gt.m7.and.i.gt.m25).or.((j1.gt.m4 .and. j1.le.m7) .and.
csa     c   (i.gt.m29 .and. i.le.m32)))goto 600
csa        if((i.gt.m7.and.j1.gt.m25).or.((i.gt.m4 .and. i.le.m7) .and.
csa     c   (j1.gt.m29 .and. j1.le.m32)))goto 600
c       delta, rho, psi etc. do not collid with the strange particles and
c       rho do not collid with pi
csa        if((j1.gt.m25 .and. j1.le.32) .and. (i.gt.m25 .and. i.le.32))
csa     c   goto 600
c       delta and rho do not scatter with each other

c       if((j1.gt.m32 .and. j1.le.m34) .and. ((i.gt.m7 .and. i.le.m9)
c     c  .or. (i.gt.m17 .and. i.le.m19)))goto 600
c       if((i.gt.m32 .and. i.le.m34) .and. ((j1.gt.m7 .and. j1.le.m9)
c     c  .or. (j1.gt.m17 .and. j1.le.m19)))goto 600
c       psi do not collid with k

602     i1=i
        iflag=0
        call rsfilt_h(j1,i1,iflag)
        if(iflag.eq.0)goto 100
        tc(nctl)=0.0
        call tcolij_h(i1,j1,time,nctl,lc,tc,tw)
c       if(tc(nctl).gt.1.0e-7) nctl=nctl+1   ! 141104
c141104
        tci=tc(nctl)
        if(tci.gt.1.0e-7)then
        do j1=1,nctl-1
        tcj=tc(j1)
        if(abs(tcj-tci).lt.dddt)goto 600
        enddo
        nctl=nctl+1
        endif
c141104
100     continue
600     enddo
300     enddo

700     if(tc(nctl).le.1.e-7) nctl=nctl-1
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
     &   3212,3112,3222,-3212,-3112,-3222,3122,-3122,311,
     &   321,3312,-3312,3322,-3322,3334,-3334,1114,2114,2214,2224,
     &   213,-213,113,443,100443,553,100553,10441,20443,445,411,-411,
     &   421,-421,4122,4112,4212,4222,223,323,313,413,-413,423,-423,
     &   431,-431,511,-511,521,-521,531,-531,541,-541,
     &   513,-513,523,-523,533,-533,543,-543,28*0/   ! 250420
        DATA DISDET/0.5,0.5,0.5,0.5,46*0.,0.5,0.5,0.5,0.5,46*0./
        DATA ISINELT/384*1,208*0,8*1/  ! with delta and rho
        ! DATA KFMAXT/72/   !Lei2023060 52 -> 72
        DATA KFMAXT/34/   !Lei2023060 to psi', for time-saving
        DATA PARAM/40.,25.,21.,10.,2.0,0.85,1.0,0.02,0.1,4.0,0.16,0.04,
     &        6.0,3.0,12.,6.,4*0/   ! 060813
        DATA WEIGH/600*1.0/
        data kjp20,vjp20,vjp21,vjp22,vjp23/1,0.3,4.0,1.5,8.0/

        END
C******************************************************************
C...........Main switches and parameters...........................
C \item[KFACOT] flavor order of considered particles
C \item[DISDET] allowable minimum distance between two
C       particles,=0.5 between two necleons,=0 otherwise
C \item[ISINELT] switch for i-th inelastic channel
C       =0 closed,=1,opened
C \item[KFMAXT](D=72) KFMAXT kinds of particles are involved in rescattering
C PARAM(1)(D=40.0mb) total cross-section of nucleon-nucleon
C PARAM(2)(D=25.0mb) total cross-section of pi-nucleon
C PARAM(3)(D=21.0mb) total cross-section of K-nucleon
C PARAM(4)(D=10.0mb) total cross-section of pi-pi
C PARAM(5)(D=2.0mb) cross-section of pi+pi --> K K
C PARAM(6)(D=0.85) ratio of inelastic cross-section to total cross-section
C PARAM(7)(D=1.0fm) formation time at rest-frame of particle
C PARAM(8)(D=0.02fm) time accuracy used in hadron cascade
C PARAM(9)(D=0.1) accuracy of four-momentum conservation
C PARAM(10)(D=4.0) size of effective rescattering region is product of 
C PARAM(10) and radius of target, origin is set on (0,0,0)
C PARAM(11)(D=0.16fm^-3) nucleon density of nucleus
C PARAM(12)(D=0.04 GeV^2/c^2) The <Pt^2> for the Gaussian distribution of 
C       spectator, not used anymore
C PARAM(13)(D=6.0mb) total cross-section of J/psi + n
C PARAM(14)(D=3.0mb) total cross-section of J/psi + meson
C PARAM(15)(D=12.0mb) total cross-section of psi' + n
C PARAM(16)(D=6.0mb) total cross-section of psi' + meson
C       kjp20 = 0 : energy dependent cross section
C             = 1 : constant cross section
C       vjp20 : constant cross section of strangeness production
C       vjp21 : cross section of pion + p to pion + Delta
C       vjp22 : cross section of pion + p to rho + p
C       vjp23 : cross section of n + n to n + Delta
Cccccccccccccccccccccccccccccccccc  end  ccccccccccccccccccccccccccc
