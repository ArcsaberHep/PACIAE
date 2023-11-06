        subroutine parcas(time_par,jjj,iijk,win,nap,rnt,rnp)
c       deals with parton cascade (partonic rescattering)
c       writen by Ben-Hao Sa at 19/11/2002
c       input messages are in 'pyjets'   ! 051122 
c       working block is 'pyjets'   ! 051122
c       output messages are in 'pyjets'   ! 051122
c160110 iiii: number of run
c       jjj: jjj-th loop interaction in a event   ! 180520
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MSCA=20000,MCLIS=280000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)  ! 240503
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)   ! 300623 Lei
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp   ! 120505
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 080520
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/scatt/pi(4),pj(4),ic,jc,n0   ! 051122
        common/work7/reac(9),crose(9)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
c       icol : current total number of collision pairs in collision time list
c       lc   : line number (in [1,n]) of colliding pair, eg.
c       lc(1,100): line number of first particle of 100-th colliding pair
c       tc   : the collision time of colliding pair
c       ic,jc: line number of colliding particles
c       n0: the n before current collision
c       pi,pj: four momentum of colliding particles
c       taup(i) : formation time of particle i.   ! 300623 Lei
c       ishp(i)=1 if i-th particle inside the simulated volume   ! 300623 Lei
c              =0 if i-th particle outside the simulated volume   ! 300623 Lei
        dimension rpo(kszj,4),ppo(kszj,4)    ! 051122
c051122 rpo,ppo: parton four coordinate before Neuton motion, four momentum 
c051122  before energy loss
c151203 iijk=0   ! 120603
        tl0=adj1(24)
c       tl0: cut off virtuality of time-like branching, i. e. Mu0**2
        if(ABS(adj1(1)).le.1D-15)return   ! 290505  300623 Lei i.e. adj1(1)=0
        time=time_par   ! 280910
        adj112=adj1(12)
        adj136=adj1(36)   ! 120505
        adj137=adj1(37)   ! 120505
        call reset_eve
c241104
c300623 ithroq_p=0   ! 300623 Lei
c300623 ithrob_p=0
c300623 ich_p=0
c300623 do i=1,4
c300623 throe_p(i)=0.
c300623 enddo
        do i1=1,n   ! 051122
        v(i1,4)=0.
        taup(i1)=0.   ! 300623 Lei
        enddo
c241104
c300623 Lei
        dpmax=adj1(27)
        drmax=PARAM(10)*dmax1(rnt,rnp)   ! 300623 Lei
        do i1=1,n
        if(k(i1,2).eq.88)then   ! Excludes junctions.
        ishp(i1)=0
        goto 200
        endif
        ppp1=p(i1,1)
        ppp2=p(i1,2)
        ppp3=p(i1,3)
        ppp4=p(i1,4)
        rrr1=v(i1,1)
        rrr2=v(i1,2)
        rrr3=v(i1,3)
        pppm=ppp1*ppp1+ppp2*ppp2+ppp3*ppp3
        if(pppm.lt.1D-28)pppm=1D-28
        if(pppm.gt.1D28)then
        ishp(i1)=0
        goto 200
        endif
        pppm=SQRT(pppm)
        rrrm=rrr1*rrr1+rrr2*rrr2+rrr3*rrr3
        if(rrrm.lt.1D-28)rrrm=1D-28
        if(rrrm.gt.1D28)then
        ishp(i1)=0
        goto 200
        endif
        rrrm=SQRT(rrrm)
        if((pppm.le.dpmax.and.ppp4.le.dpmax).and.rrrm.le.drmax)then
        ishp(i1)=1
        else
        ishp(i1)=0
        endif
200     enddo
c300623 Lei

c300623 calculate the position for the center of mass of the   ! 300623 Lei
c        non-freeze-out system. The distance of a particle, when checking
c        is it freezing out or not, is measured with respect to this center
        call copl_p(time)   ! 300623 Lei

c       step 1
c       create the parton-parton (initial) collision time list 
        call ctlcre_par(iijk)   ! 290803 
        if(iijk.eq.2)return   ! initial collis. list is empty 151203
c290803
c151203 if(iijk.eq.1)then
c151203 iiii=iiii-1
c151203 return
c151203 endif
c290803

c160110 the loop over parton-parton collisions within an event
        jjj=0
        icolo=icol   ! 120603
c       statistic of the number of loops in parton cascade within an event
24      jjj=jjj+1   ! 160110
c---------------------------------------------------------------------
        if(jjj.gt.100*icolo)then
        write(9,*)'infinite loop may have happened in'
        write(9,*)'parcas iiii,jjj,icolo=',iiii,jjj,icolo
        iiii=iiii-1
        iijk=1
        return
        endif
        n0=n   ! 051122

        if(jjj.gt.1) call copl_p(time)   ! 300623 Lei

c       step 2
c       find out the binary collision (icp) with the minimum colli. time
c       the collision time list is empty if icp=0
        call find_par(icp,tcp)
        if(icp.eq.0)goto 25 ! colli. list, empty
        ic=lc(1,icp)
        jc=lc(2,icp)
c131104
c080520 parton rescattering is assumed no longer than 10000 fm/c
        if(tcp.le.10000.)goto 27
        do i1=icp+1,icol
        lc(1,i1-1)=lc(1,i1)
        lc(2,i1-1)=lc(2,i1)
        tc(1,i1-1)=tc(1,i1)
        tc(2,i1-1)=tc(2,i1)
        enddo
        icol=icol-1
        jjj=jjj-1
        goto 24
27      continue
        time=tcp

c       step 3
c       perform classical Neuton motion for ic & jc partons
        do i=1,n   ! 051122
        do j=1,3
        rpo(i,j)=v(i,j)
c300623 vpij=p(i,j)/p(i,4)   ! 300623 Lei
c300623 v(i,j)=v(i,j)+vpij*(tcp-v(i,4))   ! 300623 Lei
        enddo
        rpo(i,4)=v(i,4)
c300623 v(i,4)=tcp   ! 300623 Lei
        enddo
        call his_p(time,rnt,rnp,istop)   ! 300623 Lei
        if(istop.eq.1) goto 25   ! 300623 Lei
c300623 istop=1 means all particles have get out of considered volume

c120603
c       consider parton energy loss phenomenologically
c120505
        if(adj136.eq.1)then
        bmax=rnt+rnp
        bmax2=bmax*bmax
        bp2=bp*bp
c       energy loss per unit distance
        dell=adj137*(nap/197.)**0.75*(1.-bp2/bmax2)
     c   *(win*win/40000)**0.25
        call eloss(dell,rpo,ppo)
        endif 
c120505
c       step 4
c       performs parton-parton collision & updates particle list 
c       if lmn=4,6,& 7 updates,'pyjets','sbe',diquark list,string 
c        list, & lc(1-2,m) either 
        kkk=0   ! 120603
        iway=0   ! 120505
        call collis(ic,jc,kf3,kf4,tcp,jjj,kkk,iway,icnew,jcnew
     c   ,lmn,time)   ! 120603 120505 160110
c120603
c120603

c       step 5

c060620 update collision list
        if(lmn.eq.4.or.lmn.eq.6.or.lmn.eq.7)goto 26 ! 070720
c       in above case 'update collision list' is executed in "collis'
        call update_ctl(ic,jc,kf3,kf4,iparres,lmn,iway,tcp)

c       'new' method
c1      if(adj136.eq.0)call update_ctl(tcp,iway)   
c1      if(adj136.eq.1 .and.iway.eq.1)call update_ctl(tcp,iway)
c       if collision does not happen remove the collision pairs (one of
c        company is colliding particle) from collision time list only
c1      if(adj136.eq.1 .and.iway.eq.0)then
c       create collision time list completely
c1      icol=0
c1      do i=1,mclis
c1      lc(1,i)=0
c1      lc(2,i)=0
c1      tc(1,i)=0.
c1      tc(2,i)=0.
c1      enddo
c1      call ctlcre_para(iijk,time)
c1      endif   
c120505
c       goto 25   ! it is actived temporally
26      continue   ! 120603
        goto 24   ! the loop over collisions within an event
25      continue
        time_par=time
c       time_par: is the time lasted in parton cascade hereafter
c250803
        do i=1,9
        reaci=reac(i)
        if(reaci.gt.0.)crose(i)=crose(i)/reaci
        enddo
c250803
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine reset_eve
c       initiate the collision time list
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MCLIS=280000)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc   ! 300623 Lei
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/scatt/pi(4),pj(4),ic,jc,n0   ! 051122
        common/work7/reac(9),crose(9)
c       reac and crose: the arraies to account for the number and
c        the value of cross section for 2->2 partonic processes
        ic=0
        jc=0
        if(iii.eq.1)then   ! 300623 Lei
        lc=0
        tc=0.
        endif   ! 300623 Lei
        icol=0
        reac=0.
        crose=0.
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ctlcre_par(iijk)
c       create the initial collision list 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MCLIS=280000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/collist/lc(2,mclis),tc(2,mclis),icol
        dddt=adj1(19)   ! 161104
        icol=1
        time=0.   ! 111599
c       every process (eg. parton casecade) starts from time equal to 0
        dminf=100.   ! 111599
        ijk=0   ! 010601
        do 100 i=2,n   ! upper diagonal 151203 051122
c       n00: total # of partons in projectile nucleus
c080603
        kfi=iabs(k(i,2))
c051122 if(kfi.gt.3 .and. kfi.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
        if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520 051122
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
c080603
        do 200 j=1,i-1   ! upper diagonal 151203 
c080603
        kfj=iabs(k(j,2))
c051122 if(kfj.gt.3 .and. kfj.ne.21)goto 200   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
        if(kfj.gt.6 .and. kfj.ne.21)goto 200   ! 080520 051122 300623 Lei kfi -> kfj, 100 -> 200
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
c080603
        if(icol.gt.mclis) then
        write(9,*)'icol over limit n,icol=',n,icol   ! sa
        stop 7777
        endif

        iflag=0   ! 300623 Lei
        call rsfilt_p(j,i,iflag)   ! 300623 Lei
        if(iflag.eq.0) goto 200   ! 300623 Lei

        tc(1,icol)=0.0
        tc(2,icol)=0.0
        kji=0   ! 240503
        call coij_p(i,j,time,icol,dminf,iff,jf,kji)
        if(kji.eq.1)then
        ijk=ijk+1   ! 010601
        goto 200
        endif
        tc1=tc(1,icol)
        tc2=tc(2,icol)
        tcicol=tc1
        if(tc1.gt.tc2)tcicol=tc2
c080304 if(tcicol.gt.0.0) icol=icol+1
c080304
        if(tcicol.gt.0.0)then   ! 121104
        tci=tcicol
        do j1=1,icol-1
        tcj=tc(1,j1)
        tck=tc(2,j1)
        if(tcj.gt.tck)tcj=tck
c161104 if(ddt.eq.0.)dddt=ddt+0.03
c161104 if(ddt.ne.0.)dddt=ddt*300
        if(dabs(tcj-tci).lt.dddt)goto 200
        enddo
        icol=icol+1
        endif     ! 121104
c080304
200     enddo
100     enddo
c300623 if(tcicol.eq.0.) icol=icol-1   ! 300623 Lei
        icol=icol-1   ! 300623 Lei
        if(icol.eq.0)then
c290803
        iijk=2   ! 1 151203
        return
c290803
c       at least one collision should occur, which has the smallest
c        'least approaching distance', that is guaranteed by the variable
c        'dminf'
c290803 icol=1
c290803 lc(1,icol)=iff
c290803 lc(2,icol)=jf
c290803 tc(1,icol)=0.02
c290803 tc(2,icol)=0.02
        endif
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine find_par(icp,tcp)
c       find out the binary collision (icp) with minimum colli. time (tcp)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MCLIS=280000)
        common/scatt/pi(4),pj(4),ic,jc,n0   ! 051122
        common/collist/lc(2,mclis),tc(2,mclis),icol
        icp=0
        tcp=10000.
        do i=1,icol
        ia=(iabs(lc(1,i)-ic)+iabs(lc(2,i)-jc))*
     c   (iabs(lc(1,i)-jc)+iabs(lc(2,i)-ic)) !it play role after first colli.
        if(ia.eq.0) goto 241
        tc1=tc(1,i)
        tc2=tc(2,i)
        tci=tc1
        if(tc1.gt.tc2)tci=tc2
c       tci=amax1(tc(1,i),tc(2,i))  ! alternative choice.
        if(tci.le.0.) goto 241
        if(tcp.lt.tci)  goto 241
        icp=i
        tcp=tci
241     continue
        enddo
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine update_ctl(ik1,ik2,kf3,kf4,iparres,lmn,iway,time)   ! 120505
c020512 update collision time list for both of w/o & w/ inelastic
c020512  parton-parton scattering 
c       ik1,ik2: line number of the colliding pair in parton list
c       kf3 and kf4: kf code of the collided pair 
c230520 if iway=1 the collision does not happen
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MCLIS=280000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/scatt/pi(4),pj(4),ic,jc,n0   ! 051122
        common/collist/lc(2,mclis),tc(2,mclis),icol
c-----------------------------------------------------------------------
        dddt=adj1(19)   ! 161104
        j=0
c       loop over old colliding pairs
        if(icol.eq.0)goto 370
        do i=1,icol
        i1=lc(1,i)
        j1=lc(2,i)
        iii=(i1-ic)*(j1-ic)
        jjj=(i1-jc)*(j1-jc)
        if(i1.eq.j1) goto 400   ! 061123 Lei Avoid the particle collideing with itself.
c       throw away the pairs composed of ic and/or jc
        if(iii.eq.0.or.jjj.eq.0) goto 400
380     continue
        tc1=tc(1,i)
        tc2=tc(2,i)
        tci=tc1
        if(tc1.gt.tc2)tci=tc2
c080104 if(tci.le.time) goto 400
c161104 if(ddt.eq.0.)dddt=ddt+0.03
c161104 if(ddt.ne.0.)dddt=ddt*300
c       throw away the pairs with tc-time<=dddt (time accuracy)
        if((tci-time).le.dddt) goto 400   ! 080104
c       proceeds for survivor
        j=j+1
        tc(1,j)=tc(1,i)
        tc(2,j)=tc(2,i)
        lc(1,j)=lc(1,i)
        lc(2,j)=lc(2,i)
400     continue
        enddo
        icol=j
        if(iway.eq.1)return   ! 120505
c-----------------------------------------------------------------
c       loop over ic,jc (new) and old partons (i.e. construct colli. pair 
c       by partons, one of which is ic or jc and another one is in parton list)
370     icol=j+1
        do 100 i=1,n0   ! 051122
        i1=i
        if(i1.eq.ic) goto 100
        if(i1.eq.jc) goto 100
c080603
        kfi=iabs(k(i,2))   ! 051122
c051122 if(kfi.gt.3 .and. kfi.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
        if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520 051122
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
c080603
        do 200 k1=1,2
        if(k1.eq.1)j1=ic
        if(k1.eq.2)j1=jc

        iflag=0   ! 300623 Lei
        call rsfilt_p(j1,i1,iflag)   ! 300623 Lei
        if(iflag.eq.0) goto 200   ! 300623 Lei

        tc(1,icol)=0.0
        tc(2,icol)=0.0
        call tcolij_par(i1,j1,time,icol)
        tc1=tc(1,icol)
        tc2=tc(2,icol)
        tcicol=tc1
        if(tc1.gt.tc2)tcicol=tc2
c020603
        if(tcicol.le.0.0)then
c       write(9,*)'i1,j1,tcicol=',i1,j1,tcicol   ! sa
        endif
c020603
c141104 if(tcicol.gt.0.0) icol=icol+1
        if(tcicol.gt.0.0)then   ! 141104
        tci=tcicol
        do j1=1,icol-1
        tcj=tc(1,j1)
        tck=tc(2,j1)
        if(tcj.gt.tck)tcj=tck
        if(dabs(tcj-tci).lt.dddt)goto 200
        enddo   ! 141104
        icol=icol+1
        endif
200     enddo
100     enddo
450     continue
c061123 if(tcicol.le.0.0) icol=icol-1   ! 061123 Lei
        icol=icol-1   ! 061123 Lei
        n0=n   ! 051122
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine update_ctlm(time,iway)   ! 120505 160110
c       a part of updating collision time list (throw away old collission 
c        pairs only) for inela. parton-parton scattering 7  ! 230520
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MCLIS=280000)
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/scatt/pi(4),pj(4),ic,jc,n0   ! 051122
        common/collist/lc(2,mclis),tc(2,mclis),icol
c-----------------------------------------------------------------------
        dddt=adj1(19)   ! 161104
        j=0
c       loop over old colliding pairs
        do i=1,icol
        i1=lc(1,i)
        j1=lc(2,i)
        iii=(i1-ic)*(j1-ic)
        jjj=(i1-jc)*(j1-jc)
        if(i1.eq.j1) goto 400   ! 061123 Lei Avoid the particle collideing with itself.
c       throw away the pairs composed of ic and/or jc
        if(iii.eq.0.or.jjj.eq.0) goto 400
380     continue
        tc1=tc(1,i)
        tc2=tc(2,i)
        tci=tc1
        if(tc1.gt.tc2)tci=tc2
c080104 if(tci.le.time) goto 400
c161104 if(ddt.eq.0.)dddt=ddt+0.03
c161104 if(ddt.ne.0.)dddt=ddt*300
c       throw away the pairs with tc-time<=dddt (time accuracy)
        if((tci-time).le.dddt) goto 400   ! 080104
        j=j+1
        tc(1,j)=tc(1,i)
        tc(2,j)=tc(2,i)
        lc(1,j)=lc(1,i)
        lc(2,j)=lc(2,i)
400     continue
        enddo
        icol=j
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine update_ctln(time,iway)   ! 120505 160110
c       update the collision time list (a part) for inela. parton-parton 
c        scattering 7
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MCLIS=280000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/scatt/pi(4),pj(4),ic,jc,n0   ! 051122
        common/collist/lc(2,mclis),tc(2,mclis),icol
c-----------------------------------------------------------------------
        dddt=adj1(19)   ! 161104
c       loop over ic (jc) and old 'pyjets' (i.e. construct colli. pair 
c        composed of partons one of which is ic (jc) and another one 
c        in old 'pyjets')
        icol=icol+1
        do 100 i=1,n0   ! 051122
        i1=i
        if(i1.eq.ic) goto 100
        if(i1.eq.jc) goto 100
c080603
        kfi=iabs(k(i,2))
c051122 if(kfi.gt.3 .and. kfi.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
        if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520 051122
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
c080603
        do 200 k1=1,2
        if(k1.eq.1)j1=ic
        if(k1.eq.2)j1=jc

        iflag=0   ! 300623 Lei
        call rsfilt_p(j1,i1,iflag)   ! 300623 Lei
        if(iflag.eq.0) goto 200   ! 300623 Lei

        tc(1,icol)=0.0
        tc(2,icol)=0.0
        call tcolij_par(i1,j1,time,icol)
        tc1=tc(1,icol)
        tc2=tc(2,icol)
        tcicol=tc1
        if(tc1.gt.tc2)tcicol=tc2
c020603
        if(tcicol.le.0.0)then
c       write(9,*)'i1,j1,tcicol=',i1,j1,tcicol   ! sa
        endif
c020603
c141104 if(tcicol.gt.0.0) icol=icol+1
        if(tcicol.gt.0.0)then   ! 141104
        tci=tcicol
        do j1=1,icol-1
        tcj=tc(1,j1)
        tck=tc(2,j1)
        if(tcj.gt.tck)tcj=tck
        if(abs(tcj-tci).lt.dddt)goto 200
        enddo   ! 141104
        icol=icol+1
        endif
200     enddo
100     enddo
450     continue
c061123 if(tcicol.le.0.0) icol=icol-1   ! 061123 Lei
        icol=icol-1   ! 061123 Lei
        n0=n   ! 051122
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coij_p(i,j,time,icp,dminf,iff,jf,kji)
c       calculate the collision time for construction
c        of initial collision time list
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MCLIS=280000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/syspar_p/rsig1,pio,tcut
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/scatt/pi(4),pj(4),ic,jc,n0   ! 051122
        dimension px(4),py(4),pij(4)
        dimension dr(3),db(3),vi(3),vj(3),pic(4),pjc(4),pxc(4),pyc(4)
        double precision b(3)
        pio=3.1416
        pi(4)=p(i,4)
        pj(4)=p(j,4)
        if(pi(4).lt.1.e-10)pi(4)=1.e-10   ! 031204
        if(pj(4).lt.1.e-10)pj(4)=1.e-10   ! 031204
        pij(4)=pi(4)+pj(4)
        pic(4)=pi(4)
        pjc(4)=pj(4)
        do k1=1,3
        pi(k1)=p(i,k1)
        pj(k1)=p(j,k1)
        pij(k1)=pi(k1)+pj(k1)
        pic(k1)=pi(k1)
        pjc(k1)=pj(k1)
        b(k1)=(pi(k1)+pj(k1))/pij(4)
        enddo
        rmi=p(i,5)
        rmj=p(j,5)
        eiej2=dot(pij,pij)
c230520 invariant mass
c       insert! energy cut
        if(eiej2.lt.0.) then
c240503
        kji=1   ! stop 2222
        return
c240503
        endif
        do n1=1,4
        px(n1)=v(i,n1)
        py(n1)=v(j,n1)
        pxc(n1)=px(n1)
        pyc(n1)=py(n1)
        enddo
        ilo=0
        kf1=k(i,2)
        kf2=k(j,2)
c160110
        ikf1=iabs(kf1)
        ikf2=iabs(kf2)
        if((ikf1.le.6.or.ikf1.eq.21).and.(ikf2.le.6.or.ikf2.eq.21))then
c       d,u,s,c,b,t quarks, their anti quarks, and gluon only
c       calculate the total cross section and decide the type of reaction
c        (when ilo=0) or sample the t value as well (when ilo=1)
        call fsig(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)! 230520
        else
        sig=0.
        endif
c160110
        if(sig.le.0.)return   ! 120603 250803
        if(ilo.eq.-2)then   ! 111999
        return   ! added by Sa on 24/06/96
        endif   ! 111999
        rsig1=dsqrt(sig/pio)

        call lorntz(0,b,pic,pjc)
        call lorntz(0,b,pxc,pyc)
        rb=0.
        bb=0.
        rr=0.
        rtai=0.
        do k1=1,3   ! 051122
        vi(k1)=pic(k1)/pic(4)
        vj(k1)=pjc(k1)/pjc(4)
        enddo
        do k1=1,3
        dr(k1)=pxc(k1)-pyc(k1)-(vi(k1)*pxc(4)-vj(k1)*pyc(4))
        db(k1)=vi(k1)-vj(k1)
        rb=rb+dr(k1)*db(k1)
        bb=db(k1)**2+bb
        rr=rr+dr(k1)*dr(k1)
        enddo
        if(bb.le.1.e-10)then
        return
        endif   ! sa
        tcol=0.-rb/bb
        do ik=1,3
        dr(ik)=dr(ik)+tcol*db(ik)
        rtai=rtai+dr(ik)*dr(ik)
        enddo
        sg=rtai
        dmin=dsqrt(sg)
        if(dmin.lt.dminf)then
        dminf=dmin
        iff=i
        jf=j
        endif
        if(dmin.gt.rsig1)then
        return
        endif
        do ik=1,3
        pxc(ik)=pxc(ik)+vi(ik)*(tcol-pxc(4))
        pyc(ik)=pyc(ik)+vj(ik)*(tcol-pyc(4))
        enddo
c       move along Newton trajectory in CMS
        pxc(4)=tcol
        pyc(4)=tcol
        call lorntz(1,b,pxc,pyc)
c       transform back to Lab.
c241104
        if(pxc(4).le.10000.and.pyc(4).le.10000.)goto 100
        return
c241104
100     lc(1,icp)=i
        lc(2,icp)=j
        tc(1,icp)=pxc(4)
        tc(2,icp)=pyc(4)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tcolij_par(i,j,time,icp)
c       calculate the collision time of i and j
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MCLIS=280000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        common/syspar_p/rsig1,pio,tcut
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/scatt/pi(4),pj(4),ic,jc,n0
        dimension px(4),py(4),dx(4),pij(4),pic(4),pjc(4),
     c   pxc(4),pyc(4)
        dimension dr(3),db(3),vi(3),vj(3),rfi(4),rfj(4)
        double precision b(3)
c300623 dimension taup(kszj)   ! 051122   ! 300623 Lei
c051122
c300623 taup=0.   ! 300623 Lei
c051122
        pi(4)=p(i,4)
        pj(4)=p(j,4)
        if(pi(4).lt.1.e-10)pi(4)=1.e-10   ! 031204
        if(pj(4).lt.1.e-10)pj(4)=1.e-10   ! 031204
        pij(4)=pi(4)+pj(4)
        pic(4)=pi(4)
        pjc(4)=pj(4)
        do k1=1,3
        pi(k1)=p(i,k1)
        pj(k1)=p(j,k1)
        pij(k1)=pi(k1)+pj(k1)
        pic(k1)=pi(k1)
        pjc(k1)=pj(k1)
        b(k1)=(pi(k1)+pj(k1))/(pi(4)+pj(4))
        enddo
        rmi=p(i,5)
        rmj=p(j,5)
        eiej2=dot(pij,pij)
c230520 squared invariant mass
c       insert! energy cut
        if(eiej2.lt.0.) then
c072200 stop 1111
        return   ! 072200
        endif
c081000 ecut=dsqrt(eiej2)-rmi-rmj
c081000 ecut0=0.02    ! P.Yang
c081000 if(ecut.le.ecut0) return
c Note! energy cut can be taken into account HERE to decide coll. pairs
c etc.  According to Y. Pang's opinions, May,1994, CCAST.(05/24/95)
        do n1=1,4
        px(n1)=v(i,n1)
        py(n1)=v(j,n1)
        pxc(n1)=px(n1)
        pyc(n1)=py(n1)
        enddo
        ilo=0
        kf1=k(i,2)
        kf2=k(j,2)
c160110
        ikf1=iabs(kf1)
        ikf2=iabs(kf2)
        if((ikf1.le.6.or.ikf1.eq.21).and.(ikf2.le.6.or.ikf2.eq.21))then
c       d,u,s,c,b,t quarks, their anti quarks, and gluon only
        call fsig(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)   ! 230520
        else
        sig=0.
        endif
c160110
cc      if(il.eq.2)return   ! put 'cc' by Sa on 24/06/96
        if(ilo.eq.-2)return   ! added by Sa on 24/06/96
        if(sig.le.0.)return   ! 120603 250803
        rsig1=dsqrt(sig/pio)

        do i1=1,3
        rfi(i1)=px(i1)+(taup(i)-time)*pi(i1)/pi(4)
        rfj(i1)=py(i1)+(taup(j)-time)*pj(i1)/pj(4)
        enddo
        rfi(4)=taup(i)
        rfj(4)=taup(j)
c       spatial coordinates of colliding particles at formation 
c        moment in Lab. frame
        call lorntz(0,b,rfi,rfj)
        ctaui=rfi(4)
        ctauj=rfj(4)
        tcol=ctaui
        if(ctaui.lt.ctauj)tcol=ctauj
c       for back to back collision the collision time is equal to 
c        the lager one of their formation times
        call lorntz(0,b,pic,pjc)
        call lorntz(0,b,pxc,pyc)
        rb=0.
        bb=0.
        rr=0.
        rtai=0.
        kflag=0
        do ik=1,3
        vi(ik)=pic(ik)/pic(4)
        vj(ik)=pjc(ik)/pjc(4)
        enddo
        do ik=1,3
        db(ik)=vi(ik)-vj(ik)
        dr(ik)=pxc(ik)-pyc(ik)-(vi(ik)*pxc(4)-vj(ik)*pyc(4))+db(ik)*tcol
        rtai=rtai+dr(ik)*dr(ik)
        enddo
        dott=0.
        do ik=1,3
        dott=dr(ik)*pic(ik)+dott
        enddo

        if(dott.ge.0.)then
        kflag=1
        if(tcol.le.pxc(4))return
        if(tcol.le.pyc(4))return
c       for the back to back collisions (collision happens in future)
        else
        rtai=0.
        do ik=1,3
        dr(ik)=pxc(ik)-pyc(ik)-(vi(ik)*pxc(4)-vj(ik)*pyc(4))
        rb=rb+dr(ik)*db(ik)
        bb=bb+db(ik)*db(ik)
        rr=rr+dr(ik)*dr(ik)
        enddo
        if(bb .le. 1.e-10)return
        tcol=0.-rb/bb
        if(tcol.le.pxc(4))return
        if(tcol.le.pyc(4))return
c020603 collision happens in future
        if(tcol-ctaui .le. 0.)return
        if(tcol-ctauj .le. 0.)return
c       collision must be after formation
        do ik=1,3
        dr(ik)=dr(ik)+db(ik)*tcol
        rtai=rtai+dr(ik)*dr(ik)
        enddo
c       for the face to face collisions

        endif

        sg=rtai
        dmin=dsqrt(sg)
        if(dmin.gt.rsig1) return
        do ik=1,3
        pxc(ik)=pxc(ik)+vi(ik)*(tcol-pxc(4))
        pyc(ik)=pyc(ik)+vj(ik)*(tcol-pyc(4))
        enddo
c       move along Newton trajectory in CMS
        pxc(4)=tcol
        pyc(4)=tcol
        call lorntz(1,b,pxc,pyc)
c       transform back to Lab.
        tcol1=pxc(4)
        tcol2=pyc(4)
        if(kflag.eq.0)then
        if(tcol1-taup(i).lt.0.)return
        if(tcol2-taup(j).lt.0.)return
        else
        if(tcol1-taup(i).lt.-1.E-4)return
        if(tcol2-taup(j).lt.-1.E-4)return
        endif
        if(tcol1.le.time)return
        if(tcol2.le.time)return
c       collision happens in the future
c241104
        if(tcol1.le.10000..and.tcol2.le.10000.)goto 100
        return
c241104
100     tc(1,icp)=tcol1
        tc(2,icp)=tcol2
        lc(1,icp)=i
        lc(2,icp)=j
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function dot(v1,v2)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        dimension v1(4),v2(4)
        dot=v1(4)*v2(4)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ppes(pei,il,ih,peo)
c       sum of momentum and energy
c       pei: two dimension array for input momentum and energy
c       il and ih: lower and higher limits of sum
c       peo: one dimension array of output momentum,energy & sqrt(s)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MPLIS=80000)   ! 030223 Yong & She
        dimension pei(mplis,5),peo(5)
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
c       write(9,*)'peo(5)=',peo(5)   ! w
        peo5=dabs(peo5)
        endif
        peo(5)=dsqrt(peo5)
        write(9,*)'sum=     px       py      pz      e    sqrt(s)' ! sa 
        write(9,100)(peo(i),i=1,5)   ! sa
100     format(4x,5(1x,f11.3))
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function ichge(kf)
c       calculate the charge (in unit of 1/3) of parton   ! 180520 
c       kf: parton flaver
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        if(kf.eq.21 .or. kf.eq.-21)then   ! 180520
        ichge=0   ! gluon charge
        elseif(kf.gt.0)then
        if(kf.eq.1)ichge=-1   ! d
        if(kf.eq.2)ichge=2   ! u
        if(kf.eq.3)ichge=-1   ! s
        if(kf.eq.4)ichge=2   ! c
        if(kf.eq.5)ichge=-1   ! b
        if(kf.eq.6)ichge=2   ! t
        elseif(kf.lt.0)then
        if(kf.eq.-1)ichge=1   ! dbar
        if(kf.eq.-2)ichge=-2   ! ubar
        if(kf.eq.-3)ichge=1   ! sbar
        if(kf.eq.-4)ichge=-2   ! cbar
        if(kf.eq.-5)ichge=1   ! bbar
        if(kf.eq.-6)ichge=-2   ! tbar
        else
        endif
        return
        end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fsig(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)!230520
c       calculate the total cross section and decide the type of reaction
c        (when ilo=0) or sample the t value as well (when ilo=1)
cc      tsmp: t value sampled.
cc      lmn: order number happened process
cc      sig: the total cross section of parton kf1 collides with kf2
cc      kf3 and kf4: kf code of the colliding pair after collision
cc      eiej2: squared invariant mass of colliding pair
c180520 jjj: jjj-th loop within a event
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)   ! 080520 

c       leadding order pQCD differential cross section of 2->2 processes
        external fs11_0,fs11_1,fs11_2,fs12_0,fs12_1,fsqq,fsgg_0,fsgg_1,
     c   fsqg

        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 250823 Lei
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc   ! 220803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 240803 181003
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 230520
        common/sa38/ i_mass, idummy, prob_ratio_q(6), am(6), amqq(6)   ! 290823 Lei
        common/syspar_p/rsig1,pio,tcut
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
c250803 common/work7/reac(9),crose(9)
        dimension ssig(0:3),eee(0:1000),dd(0:1000),eee1(0:1000),
     c   eee2(0:1000),eee3(0:1000), amqk(6)   ! 061123 Lei added amqk(6)
        logical lgg,lqg,lqq,lqaq,lqqp,lqaqp

c       classify the incident channel
        kf0=iabs(kf1)
        ikf2=iabs(kf2)   ! 090603
c090603 lqq=(kf1.eq.kf2).and.(kf2.ne.21)
        lqq=(kf1.eq.kf2).and.(kf0.ne.21)   ! 090603 ! qq (or q(-)q(-))->
c       lqq: two incoming partons are quarks (or antiquarks) with same 
c        flavor
        lgg=(kf1.eq.21).and.(kf2.eq.21)                    ! gg->
c       lgg: two gluons
c090603 lqg=((kf1.eq.21).and.(kf2.ne.21)).or.
c090603     c      ((kf2.eq.21).and.(kf1.ne.21))           ! qg (or q(-)g)->
        lqg=((kf0.eq.21).and.(kf2.ne.21)).or.
     c      ((ikf2.eq.21).and.(kf1.ne.21))           ! qg (or q(-)g)->
c       lqg: one quark (or antiquark) and one gluon
        lqaq=(kf1.eq.-kf2).and.(kf0.ne.21)          ! qq(-) (or q(-)q)->
c       lqaq: quark (antiquark) and its corresponding antiquark (quark)
c090603 lqqp=(kf0.ne.iabs(kf2)).and.(kf1.ne.21).and.(kf2.ne.21)
        lqqp=(kf0.ne.ikf2).and.(kf0.ne.21).and.(ikf2.ne.21)   ! 090603
c       lqqp: two quarks,or two antiquarks,or one quark and one antiquark
c        with different flavor.

c       calculate the total cross section for each incedent channel  
        sig=0.
        lmn=0   ! 250803
        do i=0,3
        ssig(i)=0.
        enddo
        idw=adj1(4)   ! 240803 changed from 1000
        adj120=adj1(20)   ! 230405

c061123 Lei
c250420
c       dmass=pymass(1)*2.   ! umass=pymass(2)=dmass
c       umass=pymass(2)*2.
c       smass=pymass(3)*2.
c       cmass=pymass(4)*2.
c       bmass=pymass(5)*2.
c       tmass=pymass(6)*2.
c250420
        amqk = 0D0
c       Kinematical mass
        if( i_mass.eq.1 )then
            do i=1,6,1
                am(i) = PMAS(i,1)   ! 1-6
            end do
c       Current algebra mass
        elseif( i_mass.eq.2 )then
            do i=1,6,1
                am(i) = PARF(90 + i) ! 91-96
            end do
c       Constituent mass
        elseif( i_mass.eq.3 )then
            do i=1,6,1
                am(i) = PARF(100 + i) ! 101-106
            end do
        end if
        dmass = amqk(1)*2.   ! umass <= dmass
        umass = amqk(2)*2.
        smass = amqk(3)*2.
        cmass = amqk(4)*2.
        bmass = amqk(5)*2.
        tmass = amqk(6)*2.
c061123 Lei

        xs=eiej2
c       xs: squared invariant mass of colliding pair   ! 230520
c080520
        if(eiej2.le.0.)eiej2=1.e-16
        eiej=dsqrt(eiej2)   ! invariant mass of colliding pair
c080520
        sm=0.
c       the sum of squared mass of partons in initial and final states, it 
c        is set to zero (zero mass appoxi. in kinematics treatment)
        nf=3
c112399 as=12.*pio/(33-2*nf)    ! alpha_s
c121099 as=12.*pio/(33-2*nf)/xs    ! alpha_s, 112399
c050603 as=.3    ! alpha_s, 121099, .47 (Bin)
        as=adj1(2)   ! 240803
c280202
c       uxt=0.   ! tmax
c       dxt=-xs    ! tmin
c       cf. my note book p.15
c280803 above uxt and dxt are used in limited and regularized cross section
c280202
c050603 tcut=.4   ! 120999, 5.86
csa's   tcut: the cut for avoidding divergence
csong's tcut: the cut value of t, a process is considered to be 'hard' 
csong's  when t>tcut, 'soft' when t<tcut.
        tcut=adj1(3)   ! 240803
c280803
        uxt=-tcut
cc      uxt: the upper limit (max.) of t
c       note,the integral variable t is negative 
c120699 dxt=(sm-xs)+tcut
        dxt=-xs-tcut   ! 291102 120505
c120505 dxt=-xs+tcut  ! 030903
cc      dxt: the lower limit (min.) of t
c       above uxt and dxt are used in LO pQCD cross section
c280803
        if(uxt.le.dxt) then
        ilo=-2
c       the process can not happen
        return
        endif

        if(lqaq.or.lgg) then   ! 1

        if(kf0.eq.21) then  ! 2 gg->
        call integ(fsgg_0,uxt,dxt,idw,eee,sum)  ! ->gg process 9
c       fsgg_0: differential cross section gg->gg
c160902
        sum1=sum
        do i1=0,idw
        eee1(i1)=eee(i1)
        enddo
c230520
        iparreso=iparres
        if(iparres.eq.1 .and. eiej.lt.umass)then   ! 061123 Lei dmass -> umass
        iparres=0
c       treates as elastic (iparres=0), process 9 happens
        call integ(fsgg_1,uxt,dxt,idw,eee,sum)   ! process 7 doesn't happen
c       add (common/sa33/smadel,ecce,secce,parecc,iparres) today in 'fsig'
        iparres=iparreso
        else
        call integ(fsgg_1,uxt,dxt,idw,eee,sum)  ! ->qq(-) process 7 (inelastic)
        endif
c230520
        sum2=sum
        do i1=0,idw
        eee2(i1)=eee(i1)
        enddo
        sums=sum1+sum2
        sum3=sum1/sums
        sele=pyr(1)
        if(sele.le.sum3)then   ! 3 gg->gg process 9
        do i1=0,idw
        eee(i1)=eee1(i1)
        enddo
        kf3=kf1
        kf4=kf2
c250803        reac(9)=reac(9)+1.
c250803        crose(9)=crose(9)+sum1
        lmn=9   ! 250803
        else   ! 3 gg->q*q(-) process 7 (inelastic)
        do i1=0,idw
        eee(i1)=eee2(i1)
        enddo
c       ->qq(-) branch
        call break_f(eiej,kf7,amq)   ! which is in coales_30.f 161022
        kf3=kf7   ! 161022
        kf4=-kf3   ! 250420
c250803        reac(7)=reac(7)+1.
c250803        crose(7)=crose(7)+sum2
c180520
        lmn=7   ! 250803
        endif   ! 3
        sum=sums
        sig=sums
c160902
c160902 sig=sum
c160902 kf3=kf1
c160902 kf4=kf2
c160902 reac(9)=reac(9)+1.
c160902 crose(9)=crose(9)+sum

        else       ! 2 qq(-) (or q(-)q)->
        call integ(fs11_2,uxt,dxt,idw,eee,sum) ! ->qq(-) (or q(-)q) process 5
c       fs11_2: differential cross section qq(-) (or q(-)q)->qq(-) 
c        (or q(-)q)
c160902
        sum1=sum
        do i1=0,idw
        eee1(i1)=eee(i1)
        enddo
        call integ(fsqq,uxt,dxt,idw,eee,sum) ! ->gg process 6 (inelastic)
        sum2=sum
        do i1=0,idw
        eee2(i1)=eee(i1)
        enddo
c300623 Lei
        iparreso=iparres
        if(iparres.eq.1 .and. eiej.lt.umass)then   ! 061123 Lei dmass -> umass
        iparres=0
c       treates as elastic (iparres=0), process 5 or 6 happens
        call integ(fs11_1,uxt,dxt,idw,eee,sum) ! ->q'q'(-) (or q'(-)q') process 4 (inelastic) doesn't happen.
c       add (common/sa33/smadel,ecce,secce,parecc,iparres) today in 'fsig'
        iparres=iparreso
        else
        call integ(fs11_1,uxt,dxt,idw,eee,sum) ! ->q'q'(-) (or q'(-)q') process 4 (inelastic)
        end if
c300623 Lei
        sum3=sum
        do i1=0,idw
        eee3(i1)=eee(i1)
        enddo
        sums=sum1+sum2+sum3
        sum4=sum1/sums
        sum5=(sum1+sum2)/sums
        sele=pyr(1)
        if(sele.le.sum4)then   ! 3 ->qq(-) (or q(-)q) process 5
c       sum=sum1
c       sig=sum1
        do i1=0,idw
        eee(i1)=eee1(i1)
        enddo
        kf3=kf1
        kf4=kf2
c250803        reac(5)=reac(5)+1.
c250803        crose(5)=crose(5)+sum1
        lmn=5   ! 250803
        elseif(sele.gt.sum4 .and. sele.le.sum5)then ! 3 ->gg process 6 (inelastic)
c       sum=sum1
c       sig=sum1
        do i1=0,idw
        eee(i1)=eee2(i1)
        enddo
        kf3=21
        kf4=21
c250803        reac(6)=reac(6)+1.
c250803        crose(6)=crose(6)+sum2
        lmn=6   ! 250803
c090603 elseif(sele.gt.sum5)then
        else   ! 3 ->q'q'(-) (or q'(-)q'), process 4 (inelastic)
        do i1=0,idw
        eee(i1)=eee3(i1)
        enddo
        aa=pyr(1)

c250420 branch

c300623 Lei
        do while(.true.)
            call break_f(eiej,kf7,amq)   ! which is in coales_30.f 161022
            if(kf7.ne.kf0) exit
        end do
        kf3=kf7
c300623 Lei

c250420
1002    kf4=-kf3   ! 250420
c250803        reac(4)=reac(4)+1.
c250803        crose(4)=crose(4)+sum3
        lmn=4   ! 250803
c090603 else
        endif   ! 3
c160902
        sum=sums   ! sum1 050605 Tan
        sig=sums   ! sum1 050605 Tan
        endif   ! 2

c090603 endif   ! 1

        elseif(lqq)then   ! 1 090603
        call integ(fs11_0,uxt,dxt,idw,eee,sum)   !  process 2
c       fs11_0 : differential x section of qq->qq (or q(-)q(-)->q(-)q(-))
        sig=sum
        kf3=kf1
        kf4=kf2
c250803 reac(2)=reac(2)+1.
c250803 crose(2)=crose(2)+sum
        lmn=2   ! 250803
c090603 endif

        elseif(lqqp)then   ! 1 090603
        call integ(fs12_0,uxt,dxt,idw,eee,sum)   ! process 1 and 3
c       fs_12 : differential x section of qq'->qq',differential x section 
c        of qq(-)'->qq(-)' (or q(-)q'->q(-)q')=differential x section of 
c        qq'->qq',differential x section of q(-)q(-)'->q(-)q(-)'=
c        differential x section of qq'->qq'
        if(kf1*kf2.gt.0)then   ! process 1
        sig=sum
        kf3=kf1
        kf4=kf2
c250803 reac(1)=reac(1)+1.
c250803 crose(1)=crose(1)+sum
        lmn=1   ! 250803
        endif
        if(kf1*kf2.lt.0)then   ! process 3
        sig=sum
        kf3=kf1
        kf4=kf2
c250803 reac(3)=reac(3)+1.
c250803 crose(3)=crose(3)+sum
        lmn=3   ! 250803
        endif
c090603 endif

        elseif(lqg)then   ! 1 090603 
        call integ(fsqg,uxt,dxt,idw,eee,sum)   ! process 8
c       fsqg : differential x section of qg->qg (or q(-)g->q(-)g)
        sig=sum
        kf3=kf1
        kf4=kf2
c250803 reac(8)=reac(8)+1.
c250803 crose(8)=crose(8)+sum
        lmn=8   ! 250803

        endif   ! 1

        sig=sig*0.04   ! 121099
c121099 0.04 is the transformation factor from (GeV)^-2 to (fm)^2
c280803
        if(adj120.eq.0 .or. adj120.eq.1)then   ! 230405
        do i=0,idw
        if(eee(i).lt.0.)then
        sig=0.
        return
        endif
        enddo
        endif   ! 230405
c280803

c       sample the t value
        if(ilo.eq.1) then
c230405
        if(adj120.eq.2 .or. adj120.eq.3)then
c       tsmp=-pyr(1)*xs
        seta=pyr(1)*pio
        tsmp=xs*(dcos(seta)-1)/2+tcut   ! 280407
        return
        endif
c230405
        call eee_dd(idw,eee,dd)
        call samp_integ(idw,uxt,dxt,tt,dd)
        tsmp=tt+tcut   ! 200407
c       tsmp: the t value sampled.
        endif
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine collis(ik1,ik2,kf3,kf4,tcp,jjj,kkk,iway,
     c   icnew,jcnew,lmn,time)   ! 120505 160110 
c       performs parton-parton collision & updates particle list
c       if lmn=4,6,& 7, updates 'pyjets',
c        'sbe',diquark list, string list, & lc(1-2,m) too
c       ik1,ik2: line number of the colliding pair in parton list
c180520 kf3 and kf4: kf code of the colliding pair after collision
c       tcp: collision time
c       jjj: jjj-th loop within a event
c       if kkk=1 throw away current collision
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MPLIS=80000,MSCA=20000)
        PARAMETER (KSZJ=80000,MCLIS=280000)   ! 160110  300623 Lei MCLIS
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 160110
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)   ! 160110
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5) ! 280620
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp   ! 180705
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2   ! 090820
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 160110
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 050620
        common/work7/reac(9),crose(9)
        common/ctllist_p/nreac(9),nrel
        common/show/vip(mplis),xap(mplis)
        common/collist/lc(2,MCLIS),tc(2,MCLIS),icol   ! 300623 Lei
c160110 ifcom(i): line number (in 'sbe') of first component of i-th diquark
c       nreac(i): statistics of the # of successful i-th collision
c       nrel: statistics of the # of collision blocked
cc      nsca: number of particles after collision
cc      pip(1-3,*): momentun of particle after collision
cc      pip(4,*): energy of particle after collision
cc      pip(5,*): virtuality of particle after collision
cc      pip(6,*): x value of particle after collision
cc      kpip(*): flavor code of particle after collision
cc      vip: virtuality of parton
cc      xap: momentum fraction of parton
c180705 ppsa(5): cumulate the charge losed in collision processes in an event
c       definition of ppsa(1), ..., and ppsa(4) are given in 'eloss'

        dimension pi(4),pj(4),pii(4),pjj(4),pi_o(4),pj_o(4),ps(4)
        dimension pss(mplis,5),px(4),py(4),pij(4)   ! 230520
        double precision b(3)
        adj12=adj1(12)
        kf1=k(ik1,2)
        kf2=k(ik2,2)
        stili=0.
c       stili: statistics of the times calling 'ti_li1'
c180520 kf1 and kf2: kf code of the colliding pair before interaction

c       statistics of the total charge of partons
        ichabe=0
        do i=1,n
        ik=k(i,2)
        ichabe=ichabe+ichge(ik)
        enddo

        nsca=0
        pi(4)=p(ik1,4)
        pj(4)=p(ik2,4)
        pi_o(4)=pi(4)
        pj_o(4)=pj(4)
        pij4=pi(4)+pj(4)
        do i=1,3
        pi(i)=p(ik1,i)
        pj(i)=p(ik2,i)
        pi_o(i)=pi(i)
        pj_o(i)=pj(i)
        b(i)=(pi(i)+pj(i))/pij4
        enddo
c       transfer to the current cms of colliding pair
        call lorntz(0,b,pi,pj)
c230520
        pij(1)=pi(1)+pj(1)
        pij(2)=pi(2)+pj(2)
        pij(3)=pi(3)+pj(3)
        pij(4)=pi(4)+pj(4)
        eiej2=dot(pij,pij)
c       squared invariant mass of colliding pair
c230520
c       call lorntz(1,b,pi,pj)
c120505
        if(eiej2.lt.0.)then
        iway=1
c       iway=1 means collision does not happen
        return
        endif
c120505
c       calculate the total cross section and decide the type of reaction
c        (when ilo=0) or sample the t value as well (when ilo=1)
        ilo=1
c160110
        ikf1=iabs(kf1)
        ikf2=iabs(kf2)
        if((ikf1.le.6.or.ikf1.eq.21).and.(ikf2.le.6.or.ikf2.eq.21))then
c       d,u,s,c,b,t quarks, their anti quarks, and gluon only
        call fsig(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)   ! 230520
        else
        sig=0.
        endif
c160110
c250803 lmn: order number of the process happened
c120603
        if(sig.le.0.)then
        jjj=jjj-1
        kkk=1   ! throw away that collision
        return
        endif
c120603
c       am1=amass(kf1)
c       am2=amass(kf2)
c       am3=amass(kf3)
c       am4=amass(kf4)
        am1=0.
        am2=0.
        am3=0.
        am4=0.
c       in consistent with sm in "fsig"
c230520 kinematics of 2->2 process,no matter what is the final state, is
c        treated in zero mass approximation from momentum point of view
        paa=dsqrt(pi(1)**2+pi(2)**2+pi(3)**2)
        xs=eiej2
        xt=tsmp
c       xu=am1**2+am2**2+am3**2+am4**2-xs-xt
        xu=-xs-xt   ! zero rest mass approximation
c       square of momentum transfer (Q**2)
        qc2=xt*xu/xs   ! pt**2, as the same in PYTHIA (Q^2_{hard})
c151018 qc2=4.*qc2   ! as the same in PYTHIA
c       qc2=min1(-xt,-xu)

c       calculate the direction cosines of momentum of one colliding
c        particle after scattering relative to the momentum
c        of corresponding particle before scattering
        fi=2.*pio*pyr(1)
c       relation between t and cos(ceta) in zero rest mass approximation
c        (i. e. treated as ela. scattering)
        if(paa.lt.1.e-10)paa=1.e-10   ! 031204
        cctas=1.+xt/2./paa**2
        if(dabs(cctas).gt.1.d0) then
        cctas=cctas/dabs(cctas)
        endif
        sctas=dsqrt(1.d0-cctas**2)
        cfis=dcos(fi)
        sfis=dsin(fi)
        pc1=paa
        call rotate(cctas,sctas,cfis,sfis,pc1,pi,pj)
c       pi and pj, as input in 'rotate', are four momentum of colliding 
c        pair (before scattering and rotation), as output are four 
c        momentum of that pair after scattering and rotation
c       pc1: momentum modulus of pi or pj, both are equal in their cms 
c        before and after scattering
c191103 pi(4)=dsqrt(pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)+am3*am3)
c191103 pj(4)=dsqrt(pj(1)*pj(1)+pj(2)*pj(2)+pj(3)*pj(3)+am4*am4)

        do i=1,4
        pii(i)=pi(i)   ! initial four momentum of parton-parton scattering
        pjj(i)=pj(i)
        enddo
c       boost back
        call lorntz(1,b,pi,pj)

c080520 proceeds for w/o time-like branching

c       final state of two parton scattering (nsca=2)
        nsca=2
        kpip(1)=kf3
        kpip(2)=kf4
        do i=1,4
        pip(i,1)=pi(i)
        pip(i,2)=pj(i)
        enddo

c       charge conservation
        ichbe=ichge(kf1)+ichge(kf2)
        ichaf=ichge(kf3)+ichge(kf4)
        if(ichbe.ne.ichaf)then
        write(9,*)'w/o time-like iiii,jjj,kf1,kf2=',
     c   iiii,jjj,kf1,kf2  ! sa 160110
        write(9,*)'nsca,kf3,kf4,ichbe,ichaf=',nsca,kf3,kf4,ichbe,ichaf! sa
        ppsa(5)=ppsa(5)-ichbe+ichaf   ! 180705 270407
        endif

c230520 updates particle list  (note: line numbers of the colliding 
c        pair after interaction are the same as ones before interaction)

c230520 updates 'pyjets', i.e. feedback final state of parton scattering
c        to ik1 & ik2
        l=ik1
        if(pip(4,1).lt.1.e-10)pip(4,1)=1.e-10   ! 031204
        do k1=1,3
        p(l,k1)=pip(k1,1)
        enddo
        p(l,4)=pip(4,1)
        k(l,2)=kpip(1)
        p(l,5)=amass(kpip(1))
        v(l,5)=0.

        l=ik2
        if(pip(4,2).lt.1.e-10)pip(4,2)=1.e-10   ! 031204
        do k1=1,3
        p(l,k1)=pip(k1,2)
        enddo
        p(l,4)=pip(4,2)
        k(l,2)=kpip(2)
        p(l,5)=amass(kpip(2))
        v(l,5)=0.
        reac(lmn)=reac(lmn)+1.
        crose(lmn)=crose(lmn)+sig
        nreac(lmn)=nreac(lmn)+1   ! 071103
        if(iparres.eq.0)return   ! 060520
        if(iparres.eq.1.and.(lmn.ne.4.and.lmn.ne.6.and.lmn.ne.7))return!080520
c       if(adj12.eq.1)return   ! 220820

c       proceeds for lmn=4 or =6 or =7 and hadronization by string 
c        fragmentation

        if(lmn.eq.4.or.lmn.eq.6.or.lmn.eq.7)then   !! 1

c070720 update collision list
        call update_ctlm(time,iway)   ! 070720

        if(i_time_shower.eq.1)goto 1004 ! w/ time-like branching
1008    continue
c240820 copies new string composed of final state parton pair (ik1 & ik2) 
c        to the end of 'pyjets'   ! 051122
        if(ik2.gt.ik1)then   ! 230520
        k(ik1,1)=2
        k(ik2,1)=1
c       copy ik1 to n+1
        call coend(ik1)
        icnew=n
c       copy ik2 to n+1
        call coend(ik2)
        jcnew=n
        endif
        if(ik2.lt.ik1)then   ! 230520
        k(ik1,1)=1
        k(ik2,1)=2
c       copy ik2 to n+1
        call coend(ik2)
        jcnew=n
c       copy ik1 to n+1
        call coend(ik1)
        icnew=n
        endif
        nstr1=nstr1+1
        nstr1a(nstr1)=n-1
        nstr1v(nstr1)=n
c080520 endif

c       treats scattering parton pair (scattered parton pair has been
c        treated as a new string above)

c300623 Lei-------------
c       Is ik1 a component of the normal/pure gluon string?
        call adjst(ik1,ik1str,ik1sa,ik1sv,i_g1)
        n_inStr1 = 0
        if(ik1str.ne.0) n_inStr1 = ik1sv - ik1sa
c       n_inStr1: number of entries (junction maybe) in the string to which ik1 belongs.
c       ik1str: oredr number of string to which ik1 belongs, equal 0 otherwise
c       ik1sa: line number of first component of above string, equal 0 otherwise
c       ik1sv: line number of last component of above string, equal 0 otherwise
c       i_g1: =1, ikl belongs to a pure gluon string, =0, ...normal string
        call adjst(ik2,ik2str,ik2sa,ik2sv,i_g2)
        n_inStr2 = 0
        if(ik2str.ne.0) n_inStr2 = ik2sv - ik2sa
        i_g = i_g1 + i_g2   ! 0: both ik1 and ik2 are in normal string.
                            ! 1: one of ik1 & ik2 is in gluon string.
                            ! 2: both ik1 and ik2 are in gluon string.
c       Is ik1 a component of a diquark?
        if(kf1.ne.21)then
            call adjdi(ik1,idi1,iway1)
            ifcom1 = 0
            npt1   = 0
            if(iway1.eq.2) npt1 = npt(idi1)
            if(iway1.eq.1) ifcom1 = ifcom(idi1)
c           npt1: line # (in 'sbe') of partner of idi1-th diquark
c           ifcom1: line # (in 'sbe') of first component of idi1-th diquark
        end if
        if(kf2.ne.21)then
            call adjdi(ik2,idi2,iway2)
            ifcom2 = 0
            npt2   = 0
            if(iway2.eq.2) npt2 = npt(idi2)
            if(iway2.eq.1) ifcom2 = ifcom(idi2)
        end if

        if(lmn.eq.7)then   !! 2 gg->qqbar
c       throws away scattering parton pair (ik1 & ik2)
        write(3,*)"****************************************************"
        write(3,*)"imn=7, parmov"
c300623 Lei
        if(i_g1.eq.0)then
c300623 If ik1 belongs to a normal string.
            call parmov(ik1,0,0,lmn)
        elseif(i_g1.eq.1)then
c300623 If ik1 belongs to a gluon string.   ! 300623 Lei
c        Copys (and removes) the rest gluon (not ik1) to n+1, n+2 ...
            do i_gStr=ik1sa,ik1sv,1
                if(i_gStr.ne.ik1)then   ! ikl has been moved
                    k(i_gStr,1) = 3   ! NOTE HERE!
                    call coend(i_gStr)
                    call parmov(i_g1,0,0,lmn)
                end if
            end do
        end if
        if(i_g2.eq.0)then
            call parmov(ik2,0,0,lmn)
        elseif(i_g2.eq.1)then
            do i_gStr=ik2sa,ik2sv,1
                if(i_gStr.ne.ik2)then   ! ikl has been moved
                    k(i_gStr,1) = 3   ! NOTE HERE!
                    call coend(i_gStr)
                    call parmov(i_gStr,0,0,lmn)
                end if
            end do
        end if
c300623 Lei

        if( INT(adj12).ne.0 .OR. (INT(adj12).eq.0 .AND. i_g.eq.0) )then   ! 300623 Lei
c300623 If coal or sfm with ik1 & ik2 being in normal string.   ! 300623 Lei
            if(ik1.gt.ik2)then
                call parmov(ik1,ik1,ik2,lmn)
                call parmov(ik2,ik1,ik2,lmn)
            endif
            if(ik1.lt.ik2)then
                call parmov(ik2,ik1,ik2,lmn)
                call parmov(ik1,ik1,ik2,lmn)
            endif
            n00=n-2   ! 300623 Lei New partons are above n-2 w/ or w/o branching.
            call updpli_p(n00,time)   ! 300623 Lei
            return
        elseif(INT(adj12).eq.0 .AND. i_g.ne.0)then
c300623 If sfm with ik1 &/or ik2 being in gluon string.   ! 300623 Lei
            if(i_g1.eq.1)then
c               Copys (and removes) the rest gluon (not ik1) to n+1, n+2 ...
                do i_gStr=ik1sa,ik1sv,1
                    if(i_gStr.ne.ik1)then   ! ikl has been moved
                        k(i_gStr,1) = 3   ! NOTE HERE!
                        call coend(i_gStr)
                        call parmov(i_g1,i_g1,i_g1,lmn)
                    end if
                end do
            end if
        end if
        endif   !! 2
c300623 For the case the string the ik1 and/or ik2 belongs to is pure gluon string. 
c        i.e. "g-g-g" string.
c300623 Lei------------

        if(lmn.eq.4.or.lmn.eq.6)then   !! 3, 4:q1q1bar->q2q2bar, 6:qqbar->gg

        if(adj12.ne.0)then   !! 4, fragments with coalescence 300623 Lei -> .ne.0
c       throws away scattering parton pair (ik1 & ik2)
        if(ik1.gt.ik2)then
        call parmov(ik1,ik1,ik2,lmn)
        call parmov(ik2,ik1,ik2,lmn)
        endif
        if(ik1.lt.ik2)then
        call parmov(ik2,ik1,ik2,lmn)
        call parmov(ik1,ik1,ik2,lmn)
        endif
        n00=n-2   ! 300623 Lei New partons are above n-2 w/ or w/o branching.
        call updpli_p(n00,time)   ! 300623 Lei
        return
        elseif(adj12.eq.0)then   !! 4
c       does ik1 (ik2) is a component of string
        call adjst(ik1,ik1str,ik1sa,ik1sv,i_g1)   ! 020620  300623 Lei
c       ik1str: oredr number of string to which ik1 belongs,equal 0 otherwise
c       ik1sa: line number of first component of above string,equal 0 otherwise
c       ik1sv: line number of last component of above string,equal 0 otherwise
        call adjst(ik2,ik2str,ik2sa,ik2sv,i_g2)   ! 020620  300623 Lei

c       does ik1 is a component of diquark ?
        call adjdi(ik1,idi1,iway1)
        ifcom1=0
        npt1=0
        if(iway1.eq.2)npt1=npt(idi1)   ! from 1->2 on 270620
c       npt1: line # (in 'sbe') of partner of idi1-th diquark
        if(iway1.eq.1)ifcom1=ifcom(idi1)   ! from 2->1 on 270620
c       ifcom1: line # (in 'sbe') of first component of idi1-th diquark

c       does ik2 is a component of diquark ?
        call adjdi(ik2,idi2,iway2)
        ifcom2=0
        npt2=0
        if(iway2.eq.2)npt2=npt(idi2)   ! from 1->2 on 270620
        if(iway2.eq.1)ifcom2=ifcom(idi2)   ! from 2->1 on 270620
c       ifcom1: line # of first component of idi1-th diquark

c       removes diquarks (to which ik1 & ik2 belong) from diquark list. It 
c        not effects to other lists

        if((idi1.ne.0 .and. idi2.ne.0) .and. idi1.eq.idi2)
     c   call diqmov(idi1)
        if((idi1.ne.0 .and. idi2.ne.0) .and. idi1.gt.idi2)then
        call diqmov(idi1)
        call diqmov(idi2)
        endif
        if((idi1.ne.0 .and. idi2.ne.0) .and. idi1.lt.idi2)then
        call diqmov(idi2)
        call diqmov(idi1)
        endif
        if(idi1.ne.0 .and. idi2.eq.0)call diqmov(idi1)
        if(idi1.eq.0 .and. idi2.ne.0)call diqmov(idi2)

c280620
c       moves conponents of string out & updates lists

        if(ik1str.eq.0 .and. ik2str.eq.0)then

        write(3,*)"****************************************************"
        write(3,*)"ik1str=0, ik2str=0, parmov"
c       moves ik1 & ik2 out as well as updates lists
        if(ik1.gt.ik2)then
        call parmov(ik1,ik1,ik2,lmn)   ! 070720
        call parmov(ik2,ik1,ik2,lmn)   ! 070720
        endif
        if(ik1.lt.ik2)then
        call parmov(ik2,ik1,ik2,lmn)   ! 070720
        call parmov(ik1,ik1,ik2,lmn)   ! 070720
        endif
        n00=n-2   ! 300623 Lei New partons are above n-2 w/ or w/o branching.
        call updpli_p(n00,time)   ! 300623 Lei
        return

        endif

        if(ik1str.eq.0 .and. ik2str.ne.0)then

        write(3,*)"****************************************************"
        write(3,*)"ik1str=0, ik2str /= 0, parmov"
        if(ik1.lt.ik2)then
        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)   ! 070720
        call parmov(ik1,ik1,ik2,lmn)   ! 070720
        endif
        if(ik1.gt.ik2)then
        call parmov(ik1,ik1,ik2,lmn)   ! 070720
        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)   ! 070720
        endif
        n00=n-(ik2sv-ik2sa+1)-1   ! 300623 Lei New partons are above n00 w/ or w/o branching.
        call updpli_p(n00,time)   ! 300623 Lei
        return

        endif

        if(ik1str.ne.0 .and. ik2str.eq.0)then

        write(3,*)"****************************************************"
        write(3,*)"ik1str /= 0, ik2str=0, parmov"
        if(ik1.lt.ik2)then
        call parmov(ik2,ik1,ik2,lmn)   ! 070720
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        else
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        call parmov(ik2,ik1,ik2,lmn)   ! 070720
        endif
        n00=n-(ik1sv-ik1sa+1)-1   ! 300623 Lei New partons are above n00 w/ or w/o branching.
        call updpli_p(n00,time)   ! 300623 Lei
        return

        endif

c       proceeds for both of ik1str & ik2str are not equal to zero
        if(ik1str.eq.ik2str)then
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        n00=n-(ik1sv-ik1sa+1)   ! 300623 Lei New partons are above n00 w/ or w/o branching.
        call updpli_p(n00,time)   ! 300623 Lei
        endif
        if(ik1str.gt.ik2str)then
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)   ! 070720
        n00=n-(ik1sv-ik1sa+1)-(ik2sv-ik2sa+1)   ! 300623 Lei New partons are above n00 w/ or w/o branching.
        call updpli_p(n00,time)   ! 300623 Lei
        endif
        if(ik1str.lt.ik2str)then
        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)   ! 070720
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        n00=n-(ik1sv-ik1sa+1)-(ik2sv-ik2sa+1)   ! 300623 Lei New partons are above n00 w/ or w/o branching.
        call updpli_p(n00,time)   ! 300623 Lei
        endif
        return
        endif  !! 4


        endif   !! 3

c280620

1004    continue
        write(3,*) "lmn, time-like branching"   !Lei
c0500   initial state of time-like branching
        nsca=2
        kpip(1)=kf3
        kpip(2)=kf4
        do i=1,4
        pip(i,1)=pi(i)
        pip(i,2)=pj(i)
        enddo
        pip(5,1)=qc2
        pip(5,2)=qc2
        pip(6,1)=1.
        pip(6,2)=1.
c0500   finished-------------------------------

c0700   performs time_like branching-------------------------------
        ichtb=0
        do i1=1,nsca
        kff=kpip(i1)
        ichtb=ichtb+ichge(kff)
        enddo

        nsca0=0
300     continue
        kmsca=nsca
        do i=nsca0+1,kmsca   ! do loop
        if(pip(5,i).gt.4*tl0)then   ! if loop
c       phase space for branching vanishes at a virtuality of 4*tl0
c        cf. Eq.25 in B. R. Webber, Ann. Rev. Nucl. Part. Sci., 36(86)253
        ea=pip(4,i)
        qa=pip(5,i)
        do j=1,3
        pa(j)=pip(j,i)
        enddo
        kfk=kpip(i)   !!
        xa=pip(6,i)
c       time-like branching for the final state partons & the induced
c        partons formed in time like branchings
        call ti_li1(stili,i)   ! 080520
        stili=stili+1
c220820 stili: statistics of the times calling 'ti_li1'
c-------update of i------------
c240820 kpip(i)=kfk   !!
c       do j=1,3
c       pip(j,i)=pa(j)
c       enddo
c       pip(4,i)=ea
c       pip(5,i)=qa
c240820 pip(6,i)=xa
c-------finished----------------
        endif   ! if loop
        enddo   ! do loop
        nsca0=kmsca
c220820 if(msca0.ne.nsca) goto 300

        icht=0
        do i1=1,nsca
        kff=kpip(i1)
        icht=icht+ichge(kff)
        enddo
        if(ichtb.ne.icht)then
        write(9,*)'af. ti_li1 iiii,nsca,ichtb,icht',iiii,nsca,ichtb,icht   ! sa
        ppsa(5)=ppsa(5)-ichtb+icht   ! 180705
        endif
c0700   performs time_like branching finished---------------------------

c240820
c0800   'nsca' to 'pyjets'
c       two-body scattering: kf1+kf2->kf3+kf4

        if(nsca.eq.2)then
c       updates 'pyjets', i.e. feedback final state of two parton scattering
c        to ik1 & ik2
        l=ik1
        if(pip(4,1).lt.1.e-10)pip(4,1)=1.e-10   ! 031204
        do k1=1,3
        p(l,k1)=pip(k1,1)
        enddo
        p(l,4)=pip(4,1)
        k(l,2)=kpip(1)
        p(l,5)=amass(kpip(1))
        v(l,5)=0.
        taup(l)=0.   ! 300623 Lei

        l=ik2
        if(pip(4,2).lt.1.e-10)pip(4,2)=1.e-10   ! 031204
        do k1=1,3
        p(l,k1)=pip(k1,2)
        enddo
        p(l,4)=pip(4,2)
        k(l,2)=kpip(2)
        p(l,5)=amass(kpip(2))
        v(l,5)=0.
        taup(l)=0.   ! 300623 Lei
        goto 1008
        endif

c       branching: g->gg, g->qqbar, & q-> qg (qbar->qbarg)
c       if considering only g->gg & no g->qqbar),then it is possible to 
c        construct chains of kf3->kf3+kf5, kf4->kf4+kf6
c       if nsca=3, makes 'nsca' in order of kf3 (or kf4),kf5,kf6
c       if nsca=4, makes 'nsca' in order of kf3,kf4,kf5,kf6

        if(nsca.eq.3)then   ! 1
        do i1=1,nsca
        if(kpip(i1).eq.kf3.or.kpip(i1).eq.kf4)then
c       moves particle i1 to the first position 'nsca'
        kpip1=kpip(1)
        pip1=pip(1,1)
        pip2=pip(2,1)
        pip3=pip(3,1)
        pip4=pip(4,1)
        pip5=pip(5,1)
        pip6=pip(6,1)
        kpip(1)=kpip(i1)
        do i2=1,6
        pip(i2,1)=pip(i2,i1)
        enddo
        kpip(i1)=kpip1
        pip(1,i1)=pip1
        pip(2,i1)=pip2
        pip(3,i1)=pip3
        pip(4,i1)=pip4
        pip(5,i1)=pip5
        pip(6,i1)=pip6   ! 051122 changes from 5 to 6
        endif
        enddo
        endif   ! 1

        if(nsca.eq.4)then   ! 2
        do i1=1,nsca
        if(kpip(i1).eq.kf3)then
c       moves particle i1 to the first position 'nsca'
        kpip1=kpip(1)
        pip1=pip(1,1)
        pip2=pip(2,1)
        pip3=pip(3,1)
        pip4=pip(4,1)
        pip5=pip(5,1)
        pip6=pip(6,1)
        kpip(1)=kpip(i1)
        do i2=1,6
        pip(i2,1)=pip(i2,i1)
        enddo
        kpip(i1)=kpip1
        pip(1,i1)=pip1
        pip(2,i1)=pip2
        pip(3,i1)=pip3
        pip(4,i1)=pip4
        pip(5,i1)=pip5
        pip(6,i1)=pip6   ! 051122 changes from 5 to 6
        endif
        enddo
        do i3=2,nsca
        if(kpip(i3).eq.kf4)then
c       moves particle i3 to the second position in 'nsca'
        kpip1=kpip(2)
        pip1=pip(1,2)
        pip2=pip(2,2)
        pip3=pip(3,2)
        pip4=pip(4,2)
        pip5=pip(5,2)
        pip6=pip(6,2)
        kpip(2)=kpip(i3)
        do i2=1,6
        pip(i2,2)=pip(i2,i3)
        enddo
        kpip(i3)=kpip1
        pip(1,i3)=pip1
        pip(2,i3)=pip2
        pip(3,i3)=pip3
        pip(4,i3)=pip4
        pip(5,i3)=pip5
        pip(6,i3)=pip6   ! 051122 changes from 5 to 6
        endif
        enddo
        endif   ! 2

        do i=1,3
        px(i)=v(ik1,i)   ! three position of ik1
        py(i)=v(ik2,i)   ! three position of ik2
        enddo
c       induced parton is distributed randomly in between colliding pair
        n0=n   ! 300623 Lei
        do i=3,nsca
        rl=pyr(1)
        n0=n0+1   ! 300623 Lei
        ! n=n+1   ! 300623 Lei Do not change n here, see the following
        do k1=1,3
        v(n0,k1)=px(k1)*rl+py(k1)*(1.-rl)   ! 300623 Lei i -> n0
        enddo
        v(n0,4)=tcp   ! 300623 Lei i -> n0
        v(n0,5)=0.    ! 300623 Lei
        taup(n0)=0.   ! 300623 Lei
        ishp(n0)=1   ! 300623 Lei
        enddo

c--------- keeps four momentum conservation-------
        do i2=1,4
        ps(i2)=pi_o(i2)+pj_o(i2)
c       pi_o (pj_o) four momentum of colliding pair (ik1 & ik2)
        enddo 
        do i1=1,nsca
        do i2=1,4
        pss(i1,i2)=pip(i2,i1)
        enddo
        pss(i1,5)=p(i1,5)
        enddo
c       adjust four momentum conservation by iteration,no more than
c        5000 iterations
        call cconse(pss,ps,1,nsca,1)
c       total momentum of partons after scattering (including induced partons
c        in time-like branche) is normalized to the total momentum of
c        scattering pair
        do i1=1,nsca
        do i2=1,4
        pip(i2,i1)=pss(i1,i2)
        enddo
        enddo
c---------- keep four momentum conservation finished -------------

c       a part (first & second) of 'nsca' to 'pyjets', i.e. 
c        feedback first & second partons in 'nsca' to ik1 & ik2
c       updates 'pyjets'
        l=ik1
        if(pip(4,1).lt.1.e-10)pip(4,1)=1.e-10   ! 031204
        do k1=1,3
        p(l,k1)=pip(k1,1)
        enddo
        p(l,4)=pip(4,1)
        k(l,2)=kpip(1)
        p(l,5)=amass(kpip(1))
        v(l,5)=0.
        taup(l)=0.   ! 300623 Lei
c       updates 'pyjets'
        l=ik2
        if(pip(4,2).lt.1.e-10)pip(4,2)=1.e-10   ! 031204
        do k1=1,3
        p(l,k1)=pip(k1,2)
        enddo
        p(l,4)=pip(4,2)
        k(l,2)=kpip(2)
        p(l,5)=amass(kpip(2))
        v(l,5)=0.
        taup(l)=0.   ! 300623 Lei
c0800   finished----------------------------------------------

c220820 deals with third & fourth partons in 'nsca'
        if(nsca.eq.3)then
        if(INT(adj12).eq.0.)then   ! 300623 Lei
c       moves third induced partons to 'trs' for sfm
        ntrs=ntrs+1
        ktrs(ntrs,1)=3
        ktrs(ntrs,2)=kpip(3)
        do j=1,4
        ptrs(ntrs,j)=pip(j,3)
        vtrs(ntrs,j)=v(n0,j)   ! 300623 Lei 0 -> v(n0,j)
        enddo
        ptrs(ntrs,5)=0.
        vtrs(ntrs,5)=0.
        elseif(INT(adj12).ne.0)then   ! 300623 Lei
c300623 Add this parton (gluon) into collision list for coal.   ! 300623 Lei
        n=n+1
        k(n,1)=2
        k(n,2)=kpip(3)
        k(n,3)=0
        k(n,4)=0
        k(n,5)=0
        do j=1,4
        p(n,j)=pip(j,3)
        enddo
        p(n,5)=0.
        endif   ! 300623 Lei
c300623
        goto 1008
        endif

        if(nsca.eq.4)then
c       composes third & fourth partons in 'nsca' as a string & moves it 
c        to the end of 'pyjets'
c       updates 'pyjets
        n=n+1
        l=n
        if(pip(4,3).lt.1.e-10)pip(4,3)=1.e-10   ! 031204
        do k1=1,3
        p(l,k1)=pip(k1,3)
        enddo
        p(l,4)=pip(4,3)
        k(l,1)=2   ! A
        k(l,2)=kpip(3)
        p(l,5)=amass(kpip(3))
        v(l,5)=0.
c       updates 'pyjets'
        n=n+1
        l=n
        if(pip(4,4).lt.1.e-10)pip(4,4)=1.e-10   ! 031204
        do k1=1,3
        p(l,k1)=pip(k1,4)
        enddo
        p(l,4)=pip(4,4)
        k(l,1)=1   ! V
        k(l,2)=kpip(4)
        p(l,5)=amass(kpip(4))
        v(l,5)=0.
        goto 1008
        endif
c220820

        endif   !! 1

        return   ! 160110
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ti_li1(stili,ij)   ! 080520
c       perform time-like branching (along main chain until tl0=Mn0**2
c        for branching parton)
c       stili: statistics the times calling 'ti_li1'
c       ij: order # of spliting parton in parton list
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000,MPLIS=80000)   ! 080520

        parameter (pio=3.1416)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)   ! 250823 Lei
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa38/ i_mass, idummy, prob_ratio_q(6), am(6), amqq(6)   ! 290823 Lei
        dimension pa1(4),pa2(4),p00(4)
        il=1
cc      il=1 means the branching will go on, il=0 means 'stop'.
        do i=1,4
        pa1(i)=0.
        pa2(i)=0.
        p00(i)=0.
        enddo
        sampq=0.   ! 211004
10      continue
        sampq=sampq+1.   ! 211004

c       varieties in time like branchs are all the squared mass (m**2) and 
c        m**2=E**2-|p|**2>0
c------------------------------------------------------------------
c       branching:   A  -> A1   +   A2
c       Vituality:  qqb    qqa     qqc
c       ID No.   :  kfk    kf1     kf2   !!
c        x       :  xa x'  xa1 x   xa2
c       momentum :  pa     pa1     pa2
c
c       equation : qqb=qqa/z+qqc/(1-z)+pt2/(z*(1-z))  (*)
c-------------------------------------------------------------------
c       above Eq. is deduced from relation of mass and energy (Einstein 
c        relation) and four momentum conservation, cf. note made in 
c        Webber's paper

c-------flavor of A1,A2 will be decided in the ensueing module.
        q2max=qa
        kf0=kfk   !!
        call suda_ti(ilo,kf0,ii,q2max,q2,z)
c       time_like Sudakov factor, input: kf0 & q2max, output: others
c       ilo: =0 fail in sampling q2 and z
c            =1 success in sampling q2 and z
        if(ilo.eq.0) then
        qa=tl0   ! forced
        return
        endif

        qqb=q2
        if(stili.gt.2.)qqb=qa
        zab=z
c       'zab' is the sampled ratio of x/x' in this branching.
        if(zab.eq.0.)return   ! 031204
        xa1=xa*zab
        xa2=xa-xa1

        if(kfk.ne.21) then   !! initial state is quark
c       branch q->qg
        kf1=kfk
        kf2=21
        if(pyr(1).ge.0.5) then
c       branch q->gq
        kf1=21
        kf2=kfk
        endif
        else   !! initial state is gluon
        if(ii.eq.2) then   !!!
c       branch: g->gg
        kf1=21
        kf2=21
        else   !!! ii=3
c220820 considering g->gg only
        goto 200   ! 220820
c       branch: g->qq(-)
c250420
c       ea=pip(4,ij): energy of spliting particle ij
c       pip(1-3,ij): three momentum of spliting particle ij
        eaa=ea*ea-pip(1,ij)*pip(1,ij)-pip(2,ij)*pip(2,ij)-
     c   pip(3,ij)*pip(3,ij)   ! 080520
        eaa=sqrt(eaa)   ! invariant mass of initial state gluon

c061123 Lei
c       if(eaa.lt.(PYMASS(1)*2.)) goto 200   ! 300623 Lei
        amu = 0D0
c       Kinematical mass
        if( i_mass.eq.1 )then
            amu = PMAS( 2, 1 )
c       Current algebra mass
        elseif( i_mass.eq.2 )then
            amu = PARF( 90 + 2 )
c       Constituent mass
        elseif( i_mass.eq.3 )then
            amu = PARF( 100 + 2 )
        end if
        if( eaa.lt.amu*2. ) goto 200
c061123 Lei

        call break_f(eaa,kff,amq)   ! which is in coales_30.f 161022
        kf1=kff   ! 161022
        kf2=-kf1   ! 250420
        if(pyr(1).ge.0.5) then
        kf1=kf2
        kf2=-kf2
        endif
200     endif   !!!
        endif   !!
c----------finish

100     continue
        q2max=zab*(qqb-tl0/(1.-zab))
c       the max. value of A1's virtuality decided by the eqution(*) with
c        qqc=tl0 & pt2=0 .
        call suda_ti(ilo,kf1,ii,q2max,q2,z)
        if(ilo.eq.0) then
        qqa=tl0
        il=0
c       qqa=tl0 means the forward evolution has finished, thus set il=0.
        else
        qqa=q2
        endif
        q2max=(qqb-qqa/zab)*(1.-zab)
c       the max. value of A2's virtuality decided by the eqution(*) with pt2=0
        call suda_ti(ilo,kf2,ii,q2max,q2,z)
        if(ilo.eq.0) then
        qqc=tl0
        else
        qqc=q2
        endif
        pt2=(qqb-qqa/zab-qqc/(1.-zab))*(1.-zab)*zab
c       pt2: the square of transverse momentum decided by the equation(*)

c0130   the following block gives the momentum of A1,A2 i.e. pa1,pa2
c        according to the definition of z & pt.
        pa02=pa(1)**2+pa(2)**2+pa(3)**2
        paz=dsqrt(pa02)   ! third momentum of A 
c       note Eq. (*) with assumption that z axis is set on momentum of A,
c        i. e. z axis of the frame, in branching, is set on the pa
c        direction
        pz1=zab*paz   ! third momentum of A1
        pz2=(1-zab)*paz   ! third momentum of A2
        ppa1=dsqrt(pz1**2+pt2)   ! momentum modulus of A1
        ppa2=dsqrt(pz2**2+pt2)   ! momentum modulus of A2
        if(ppa1.lt.1.e-10)ppa1=1.e-10   ! 031204
        if(ppa2.lt.1.e-10)ppa2=1.e-10   ! 031204
        ppt=dsqrt(pt2)
        sctas1=ppt/ppa1
        cctas1=pz1/ppa1
        fi=2.*pio*pyr(1)
        cfis=dcos(fi)
        sfis=dsin(fi)
c       direction cosines of A1 relative to momentum of A
c       note : transverse momentum of A1 & A2 must be in same plan in
c        order to keep transverse momentum conservation
        sctas2=-ppt/ppa2
        cctas2=pz2/ppa2
        do i=1,3
        pa1(i)=pa(i)
        pa2(i)=pa(i)
        enddo
        call rotate(cctas1,sctas1,cfis,sfis,ppa1,pa1,p00)
c       originally, momentum of A1 is relative to the
c        direction of momentum of A
c       after rotation, momentum of A1 is relative to the cms where A
c        belongs to
        call rotate(cctas2,sctas2,cfis,sfis,ppa2,pa2,p00)
c0130   finished---------------------------------------

c0140   new induced parton (A2) should be added to particle list.
        nsca=nsca+1
        kpip(nsca)=kf2
        pip(4,nsca)=amass(kf2)**2
        do i=1,3
        pip(4,nsca)=pip(4,nsca)+pa2(i)**2
        pip(i,nsca)=pa2(i)
        enddo
        pip(4,nsca)=dsqrt(pip(4,nsca))
c       pip(5,nsca)=qa2
        pip(5,nsca)=qqc
        pip(6,nsca)=xa2
c0140   finished--------------------------------------------

c------- update five very important values of KFK,QA,XA,PA,EA (time like)
c        & then forward evolution along main chain further
        kfk=kf1   !!
        qa=qqa
        do i=1,3
        pa(i)=pa1(i)
        enddo
        ea=dsqrt(pa(1)**2+pa(2)**2+pa(3)**2+amass(kfk)**2)
        xa=xa1

c080520 if(sampq.gt.10.)return   ! 211004
        if(sampq.ge.1.)return   ! once branching only 080520
        if(il.eq.1) goto 10
c       il=0 means the evolution 'stop' while il=1 'continue'
c       the time like branching of A2 parton will consider later, i.e.
c        in 'collis' after statement numbered 300
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine suda_ti(ilo,kf0,ii,q2max,q2,z)
c       time-like Sudakov factor
c       perform time-like forward branching,one step only
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/sa24/adj1(40),nnstop,non24,zstop
c----------------------------------------------------------------------
c       process    : A  ->   A1  +  A2
c       virtuality : qqb     qqa    qqc
c
c       equation: qqb=qqa/z+qqc/(1-z)+pt2/(1-z)z
c----------------------------------------------------------------------
        yint1(x)=-dlog(1.d0-x)          !q->qg cf. Sa' note p. 44
        yint2(x)=dlog(1.d0/(1.d0-x)-1.d0)   !g->gg
c       yint1,2: the integral functions of part of the splitting functions
        adj25=adj1(25)
c       adj25: lamda_QCD
        adj25=adj25*adj25
        ilo=0
        if(q2max.le.4*tl0) return   ! cf. W' paper p.264 Eq.25
        ikf=iabs(kf0)
        ii=0
        tmax=dlog(q2max/adj25)
        tmin=dlog(4.*tl0/adj25)   ! cf. W' paper p.264
c       tmax,tmin: bounds of t=dlog(Q**2/lamda) in the Sudakov form factor
        zmin=tl0/q2max     ! cf. W' paper p.263
        zmax=1-tl0/q2max   ! cf. W' paper p.263
c       approximated values of the bounds of z
        if(zmax.le.zmin) return
        if(ikf.ne.21) then   ! q->qg
        ciqg=4./3.*((2.*yint1(zmax)-zmax-0.5*zmax**2)-
     *              (2.*yint1(zmin)-zmin-0.5*zmin**2)) ! cf. Sa' note p. 43
c       ciqg: the integration of splitting function for q->qg
        cisum=ciqg
        ii=1       ! q->qg
        else       ! g->
        cigg=(6.*yint2(zmax)-12.*zmax+3.*zmax**2-2.*zmax**3)-
     *       (6.*yint2(zmin)-12.*zmin+3.*zmin**2-2.*zmin**3)
c       cigg: the integration of splitting function for g->gg
        ciqq=1./2.*((zmax-zmax**2+2./3.*zmax**3)-
     *              (zmin-zmin**2+2./3.*zmin**3))
c       ciqq: the integration of splitting function for g->qq(-)
        cisum=cigg+ciqq
        aaa=pyr(1)
        if(aaa.le.(cigg/cisum)) then
        ii=2       ! g->gg
        cisum=cigg
        else
        ii=3       ! g->qq(-)
        cisum=ciqq
        endif
        endif
        ce=9./2./cisum
        sampz=0.   ! 211004
100     continue
        aaa=pyr(1)
c       tt=tmax*aaa**ce   ! variable alpha_s
        tt=tmax+6.2832/as/cisum*dlog(aaa)   ! constant alpha_s
c       tt: the t value sampled from time-like Sudakov factor
        if(tt.le.tmin) return
        q2=adj25*dexp(tt)
c       adj25: (\Lambda_s)**2
        rzmax=0.5d0+0.5*dsqrt(1.-4.*tl0/q2)
c       cf. Eq.24 in B. R. Webber, Ann. Rev. Nucl. Part. Sci. 36(1986)253
c        that journal is simplified as W' elsewhere
        rzmin=1.-rzmax
c       rzmax,rzmin: exact value of zmax & zmin corresponding to the t 
c        sampled.
200     continue
c0170   sample z.
        aaa=pyr(1)
        bbb=pyr(1)
        if(ii.eq.1) then   ! q->qg
        zsp=aaa*yint1(zmax)+(1.-aaa)*yint1(zmin)
        zsp=1.d0-dexp(-zsp)
        ratio=(1.+zsp*zsp)/2.   ! accepting probability
c220702 if(ratio.le.bbb) goto 200
        if(bbb.gt.ratio) goto 200 ! the same as previous statement
        endif
        if(ii.eq.2) then   ! g->gg
        zsp=aaa*yint2(zmax)+(1.-aaa)*yint2(zmin)
        zsp=1.d0-1.d0/(1.d0+dexp(zsp))
        ratio=(1.d0-zsp*(1.-zsp))**2
c220702 if(ratio.le.bbb) goto 200
        if(bbb.gt.ratio) goto 200 ! the same as previous statement
        endif
        if(ii.eq.3) then   ! g->qq(-)
        zsp=aaa*zmax+(1.-aaa)*zmin
        ratio=1.-2.*zsp+2.*zsp**2
c220702 if(ratio.le.bbb) goto 200
        if(bbb.gt.ratio) goto 200 ! the same as previous statement
        endif
c0170   sample z finished--------------------------------------------------
        if(zsp.gt.rzmax.or.zsp.lt.rzmin) then
c       if the z sampled fall out of the exact region, reject the sampled t,
c       let tmax equal to the t sampled, and go back to resample.
        sampz=sampz+1.   ! 211004
        if(sampz.gt.1000.)then   ! 211004
        ilo=0   ! 211004
        return   ! forced stopping time-like branching 211004
        endif   ! 211004
        tmax=tt
        goto 100
        endif
        z=zsp
        ilo=1
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cconse(pp,ps,npl,np,nstep)
c       adjust four momentum conservation by iteration,no more than
c        5000 iterations
c       pp : four momentum of particles
c       ps : above four momenta should be conserved to ps
c       npl : line # of the first particle
c       np : line # of last particle
c       nstep : interval of the step
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MPLIS=80000)
        dimension pp(mplis,5),ps(4),ff(mplis),pxyz(3),arp(3)
        dep=1.e-5
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
c082200        if(dabs(1.-fr) .le. dep)goto 200
        do i=npl,np,nstep
        amas=pp(i,5)
        amas2=amas*amas
        ppm=pp(i,4)
        ppf=ppm/fr
        ff(i)=dsqrt(dabs(ppf*ppf-amas2)/(ppm*ppm-amas2))
        do j=1,3
        ppp=ff(i)*pp(i,j)
        pp(i,j)=ppp
        pxyz(j)=pxyz(j)+ppp
        enddo
        enddo
        do i=1,3
        arp(i)=dabs(1.-pxyz(i)/ps(i))
        enddo
        if(dabs(1.-fr).le.dep .and. arp(1).le.dep .and. arp(2).le.dep
     c   .and. arp(3).le.dep) goto 200
        do i=1,3
        pxyz(i)=pxyz(i)-ps(i)
        pxyz(i)=pxyz(i)/(dfloat(np-npl)/dfloat(nstep)+1)
        enddo
        do i=npl,np,nstep
        do j=1,3
        pp(i,j)=pp(i,j)-pxyz(j)
        enddo
        pp5=pp(i,5)
        pp52=pp5*pp5
        pp(i,4)=dsqrt(pp52+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
        enddo
        jj=jj+1
        if(jj.lt.5000)goto 100
200     return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fap_qq(zz)
c       split function for g-->q*q(-)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        fap_qq=(zz**2+(1.-zz)**2)/2.
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fap_gg(zz)
c       split function for g-->g*g
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        fap_gg=6.*(1.-zz*(1.-zz))**2/(zz*(1.-zz))
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fap_qg(zz)
c       split function for q-->q*g
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        fap_qg=4.*(1.+zz**2)/3./(1.-zz)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_ud(xa)
c       the struction function of valence u,d quark
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        external fun_beta
        ata1=0.419+0.004*dta-0.007*dta*dta
        ata2=3.46+0.724*dta-0.066*dta*dta
        rud=4.4-4.86*dta+1.33*dta*dta
        anud=3./(fun_beta(ata1,ata2+1)*(1+rud*ata1/(ata1+ata2+1.)))
        aaa=anud*xa**ata1*(1.-xa)**ata2*(1.+rud*xa)/xa
        fun_ud=aaa
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_d(xa)
c       the struction function of d quark
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        external fun_beta
        ata3=0.763-0.237*dta+0.026*dta*dta
        ata4=4.+0.627*dta-0.019*dta*dta
        rd=-0.421*dta+0.033*dta*dta
        andd=1./(fun_beta(ata3,ata4+1)*(1+rd*ata3/(ata3+ata4+1.)))
        aaa=andd*xa**ata3*(1.-xa)**ata4*(1.+rd*xa)/xa
        fun_d=aaa
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_u(xa)
c       the struction function of u quark
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        fun_u=fun_ud(xa)-fun_d(xa)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_g(xa)
c       the struction function of gluon
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        aag=1.56-1.71*dta+0.638*dta*dta
        ag=-0.949*dta+0.325*dta*dta
        bg=6.+1.44*dta-1.05*dta*dta
        a1g=9.-7.19*dta+0.255*dta*dta
        a2g=-16.5*dta+10.9*dta*dta
        a3g=15.3*dta-10.1*dta*dta
        fun_g=aag*xa**ag*(1.-xa)**bg*
     c   (1.+a1g*xa+a2g*xa*xa+a3g*xa*xa*xa)/xa
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_s(xa)
c       the struction function of sea quark
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        aas=1.265-1.132*dta+0.293*dta*dta
        ass=-0.372*dta-0.029*dta*dta
        bs=8.05+1.59*dta-0.153*dta*dta
        a1s=6.31*dta-0.273*dta*dta
        a2s=-10.5*dta-3.17*dta*dta
        a3s=14.7*dta+9.8*dta*dta
        fun_s=aas*xa**ass*(1.-xa)**bs*
     c   (1.+a1s*xa+a2s*xa*xa+a3s*xa*xa*xa)/xa/6.
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_gal(xa)
c       function of logarithm Gammma function 
c       real*8 coef(0:6)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        dimension coef(0:6)
        data coef/76.18009173,-86.50532033,24.01409822,-1.231739516,
     c   0.00120858003,-5.36382e-6,2.50662827565/
        xx=xa-1.0
        temp=xx+5.5
        temp=(xx+0.5)*dlog(temp)-temp
        ser=1.0
        do i=0,5
        xx=xx+1.
        ser=ser+coef(i)/xx
        enddo
        co=coef(6)*ser
        fun_gal=temp+dlog(co)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_beta(xa,ya)
c       function of Beta function
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        b=fun_gal(xa)+fun_gal(ya)-fun_gal(xa+ya)
        fun_beta=dexp(b)
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs12_0(xt)
c       differential cross section of q1*q2 -> q1*q2, process 1
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c       xs: s
c       sm: ma**2+mb**2+mc**2+md**2
c       xu: u
c       xt: t
c120699 ccc=dlog(xt*xu/xs/adj1(25).adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or.adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fs12_0=fsabcd*(xs**2+xu**2)*4./xt**2/9.*adj11   ! 210803
c240803 LO pQCD
        return
100     fs12_0=pioa*8./9./xt**2*adj11   !240803
c       limited and regularized (B. Zhang)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs11_0(xt)
c       differential cross section of q1*q1 -> q1*q1, process 2
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fs11_0=fsabcd*(4.*((xs**2+xu**2)/xt**2+(xs**2+xt**2)/xu**2)/9.-
     c   8.*xs**2/27./xu/xt)*adj11   ! 210803
c240803 LO pQCD
        return
100     fs11_0=pioa*8./9./xt**2*adj11   ! 240803
c       limited and regularized 
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs12_1(xt)
c       differential cross section of q1*q2(-) -> q1*q2(-), process 3
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fs12_1=fsabcd*(4.*(xs**2+xu**2)/xt**2/9.)*adj11   ! 210803
c240803 LO pQCD
        return
100     fs12_1=pioa*8./9./xt**2*adj11   ! 240803
c       limited and regularized
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs11_1(xt)
c       differential cross section of q1*q1(-) -> q2*q2(-), process 4
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2   ! 080820
c240412
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 131212
        if(iparres.eq.0 .or. (iparres.eq.1.and.i_inel_proc.eq.7))then
        fs11_1=0.   ! 160110
        return   ! 160110
        endif
c240412
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fs11_1=fsabcd*(4.*(xt**2+xu**2)/xs**2/9.)*adj11   ! 210803
c240803 LO pQCD
        return
100     fs11_1=0.   ! pioa*4./9.*(1.+xt*xt/xs*xs)*adj11   ! 240803
c       limited and regularized
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs11_2(xt)
c       differential cross section of q1*q1(-) -> q1*q1(-), process 5
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fs11_2=fsabcd*(4.*((xs**2+xu**2)/xt**2+(xu**2+xt**2)/xs**2)/9.-
     c   8.*xu**2/27./xs/xt)*adj11   ! 210803
c240803 LO pQCD
        return
100     fs11_2=pioa*8./9./xt**2*adj11   !
c       limited and regularized
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fsqq(xt)
c       differential cross section of q*q(-) -> g*g, process 6
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2   ! 080820
c240412
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 131212
        if(iparres.eq.0 .or. (iparres.eq.1.and.i_inel_proc.eq.7))then
        fsqq=0.   ! 160110
        return   ! 160110
        endif
c240412
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fsqq=fsabcd*(32.*(xu**2+xt**2)/xu/xt/27.-
     c   8.*(xu**2+xt**2)/xs**2/3.)*adj11   ! 210803
c240803 LO pQCD
        return
100     fsqq=-pioa*32./27./xs/xt*adj11   ! 240803
c       limited and regularized
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fsgg_1(xt)
c       differential cross section of g*g ->q*q(-), process 7
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
c240412
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 131212
        if(iparres.eq.0)then
        fsgg_1=0.   ! 160110
        return   ! 160110
        endif
c240412
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fsgg_1=fsabcd*((xu**2+xt**2)/xu/xt/6.-3.*(xu**2+xt**2)/xs**2/8.)
     c   *adj11   ! 210803
c240803 LO pQCD
        return
100     fsgg_1=-pioa*1./3./xs/xt*adj11   ! 240803
c       limited and regularized
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fsqg(xt)
c       differential cross section of q*g -> q*g, process 8
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fsqg=fsabcd*((xu**2+xs**2)/xt**2-4.*(xu**2+xs**2)/xu/xs/9.)
     c   *adj11   ! 210803
c240803 LO pQCD
        return
100     fsqg=pioa*2./xt**2*adj11   ! 240803 200407
c       limited and regularized
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fsgg_0(xt)
c       differential cross section of g*g -> g*g, process 9
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c   ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
        adj20=adj1(20)
c210803
c020903 xt=xt-tcut
c071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
c120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699 ccc=ccc**2
c120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
        fsabcd=pioa/xs**2   ! 120699
        fsgg_0=fsabcd*(9.*(3.-xu*xt/xs**2-xu*xs/xt**2-xs*xt/xu**2)/2.)
     c   *adj11   ! 210803
c240803 LO pQCD
        return
100     fsgg_0=pioa*9./2./xt**2*adj11   ! 240803
c       limited and regularized
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function amass(kf)
c       mass of partons
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        k1=iabs(kf)
        if(k1.eq.21) r=1.e-15
        if(k1.eq.2) r=0.333   ! 080520
        if(k1.eq.1) r=0.333   ! 080520
        if(k1.eq.3) r=0.5   ! 080520
        if(k1.eq.4) r=1.5   ! 080520
        if(k1.eq.5) r=4.8   ! 080520
        if(k1.eq.6) r=174   ! 080520
        amass=r
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine integ(fun,ux,dx,idw,eee,sum)
c       calculate distribution (eee) and integration (sum) of the 'fun',which
c081018  is a y=t-tcut distribution
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        dimension eee(0:idw)
        external fun
        dw=idw
        del=(ux-dx)/dw
c       calculate distribution function of external function 'fun'
        sum=0.
        do i=0,idw
        eee(i)=fun(dx+del*i)
        sum=sum+eee(i)
        enddo
        sum=sum-0.5*(eee(0)+eee(idw))
        sum=sum*del
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine eee_dd(idw,eee,dd)
c       calculate the relative integral distribution function 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        dimension dd(0:idw),eee(0:idw)
        dd(0)=0.
        do i=1,idw
        dd(i)=dd(i-1)+eee(i)+eee(i-1)
c280803 dd(i)=dd(i-1)+dabs(eee(i))+dabs(eee(i-1))
        enddo
        do i=1,idw
        dd(i)=dd(i)/dd(idw)
        enddo
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine samp_integ(idw,ux,dx,xf,dd)
c       samples a value from relative integral distribution function 'dd',
c        which is steaming from a corresponding differential y=t-tcut 
c        distribution function
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        dimension dd(0:idw)
        dw=idw
        del=(ux-dx)/dw
        ran=pyr(1)
        ir=0
        do i=1,idw
        if(dd(i).lt.ran) goto 10
        ir=i
        goto 20
10      continue
        enddo
20      continue
        xf=pyr(1)
        xf=dx+del*(ir-1)+xf*del
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine eloss(dell,rpo,ppo)
c       consider parton energy loss phenomenologically
c       dell: energy loss per unit distance
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MPLIS=80000,MSCA=20000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension rpo(mplis,4),ppo(mplis,4)
c       rpo,ppo: parton four coordinate before Newton motion, four momentum 
c        before energy loss
        dimension rr(3),pi(3),pia(3),deltp(3)
c       ppsa(1): cumulate the px losed in an event
c       ppsa(2): cumulate the py losed in an event
c       ppsa(3): cumulate the pz losed in an event
c       ppsa(4): cumulate the e losed in an event
        adj138=adj1(38)
        no=n
        do i=1,no   ! 1
        pei=p(i,4)   ! energy before loss
        ppo(i,4)=pei
        do j=1,3
        rr(j)=v(i,j)-rpo(i,j)
        pi(j)=p(i,j)   ! three momentum before loss
        ppo(i,j)=pi(j)
        enddo
        pt=pi(1)*pi(1)+pi(2)*pi(2)
        if(pt.gt.1.e20)pt=1.e20
        if(pt.lt.1.e-20)pt=1.e-20
        pt=dsqrt(pt)
        if(pt.le.adj138)goto 100
c       below ptmin='adj138' Gev/c jet can not loss energy
        rm=rr(1)*rr(1)+rr(2)*rr(2)+rr(3)*rr(3)
        if(rm.gt.1.e20)rm=1.e20
        if(rm.lt.1.e-20)rm=1.e-20
        rm=dsqrt(rm)
        delte=rm*dell
        if(delte.ge.pei)goto 100
c       'delte' should be less than energy of particle
        if(pei.lt.1.e-20)pei=1.e-20
        srtf=1.-delte/pei
        p(i,4)=pei-delte   ! energy after loss
        do j=1,3
        pia(j)=srtf*pi(j)   ! three momentum after loss
c       massless and collinear approximations
        deltp(j)=pi(j)-pia(j)   ! three momentum loss
        p(i,j)=pia(j)   ! three momentum after loss
        ppsa(j)=ppsa(j)+pi(j)-pia(j)   ! cumulate three momentum loss
        enddo
        ppsa(4)=ppsa(4)+delte   ! cumulate energy loss
c       tansfer losed energy to a gluon (added to the end of particle list)
        call eto1g(i,delte,deltp)
c       tansfer losed energy to two gluons (added to the end of particle list)
c       call eto2g(i,delte,deltp)
100     continue
        enddo   ! 1
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ctlcre_para(iijk,time)   ! 120505
c       creates the collision time list after energy loss
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MCLIS=280000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/collist/lc(2,mclis),tc(2,mclis),icol
        dddt=adj1(19)   ! 161104
c       adj136=adj1(36)
        icol=1
        do 100 i=2,n   ! upper diagonal
        kfi=iabs(k(i,2))
c051122 if(kfi.gt.3 .and. kfi.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
        if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
        do 200 j=1,i-1   ! upper diagonal
        kfj=iabs(k(j,2))
c051122 if(kfj.gt.3 .and. kfj.ne.21)goto 200   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
        if(kfj.gt.6 .and. kfj.ne.21)goto 200   ! 080520 051122 300623 Lei kfi -> kfj, 100 -> 200
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
        if(icol.gt.mclis) then
        write(9,*)'icol over limit n,icol=',n,icol   ! sa
        stop 77777
        endif

        iflag=0   ! 300623 Lei
        call rsfilt_p(j,i,iflag)   ! 300623 Lei
        if(iflag.eq.0) goto 200   ! 300623 Lei

        tc(1,icol)=0.0
        tc(2,icol)=0.0
        call tcolij_par(i,j,time,icol)
        tc1=tc(1,icol)
        tc2=tc(2,icol)
        tcicol=tc1
        if(tc1.gt.tc2)tcicol=tc2
        if(tcicol.le.0.0)then
        endif
        if(tcicol.gt.0.0)then   ! 141104
        tci=tcicol
        do j1=1,icol-1
        tcj=tc(1,j1)
        tck=tc(2,j1)
        if(tcj.gt.tck)tcj=tck
        if(dabs(tcj-tci).lt.dddt)goto 200
        enddo   ! 141104
        icol=icol+1
        endif
200     continue
100     continue
c300623 if(tcicol.le.0.0) icol=icol-1   ! 300623 Lei
        icol=icol-1   ! 300623 Lei
        n0=n
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine eto1g(ii,delte,deltp)
c       transfer the losed energy to a gluon (added to the end of 
c        particle list)
c       ii: the line number (in 'pyjets') of parton which lossing energy
c       delte: energy losed
c       deltp(3): three momentum losed
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        dimension deltp(3)
        nl=n+1
        if(nl.gt.kszj)then
        write(9,*)'in parcas, kszj needs to be enlarged nl=',nl   ! sa
        stop 9999
        endif
        pp4=delte
        p(n+1,4)=pp4
        do i=1,3
        ppi=deltp(i)
        p(n+1,i)=ppi
        enddo
        p(n+1,5)=0.
        k(n+1,2)=21
        v(n+1,4)=v(ii,4)
        taup(n+1)=0.   ! 300623 Lei
        ishp(n+1)=1    ! 300623 Lei
        do i=1,3
        v(n+1,i)=pyr(1)*v(ii,i)
        enddo
        n=n+1
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine eto2g(ii,delte,deltp)
c       transfer the losed energy to two gluons (added to the end of 
c        particle list)
c       ii: the line number (in 'pyjets') of parton which losed energy
c       delte: energy losed
c       deltp(3): three momentum losed
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        dimension deltp(3)
        nl=n+1
        if(nl.gt.kszj)then
        write(9,*)'in parcas, kszj needs to be enlarged nl=',nl   ! sa
        stop 9999
        endif
        nl=n+2
        if(nl.gt.kszj)then
        write(9,*)'in parcas, kszj needs to be enlarged nl=',nl   ! sa
        stop 99999
        endif
        pp4=pyr(1)*delte
        if(pp4.lt.1.e-20)pp4=1.e-20
        p(n+1,4)=pp4
        pp42=delte-pp4
        if(pp42.lt.1.e-20)pp42=1.e-20
        p(n+2,4)=pp42
        do i=1,3
        ppi=pyr(1)*deltp(i)
        p(n+1,i)=ppi
        ppi2=deltp(i)-ppi
        p(n+2,i)=ppi2
        enddo
        p(n+1,5)=0.
        p(n+2,5)=0.
        k(n+1,2)=21
        k(n+2,2)=21
        k(n+1,1)=1   ! 051122
        k(n+2,1)=1   ! 051122
        v(n+1,4)=v(ii,4)
        v(n+2,4)=v(ii,4)
        taup(n+1)=0.   ! 300623 Lei
        taup(n+2)=0.   ! 300623 Lei
        ishp(n+1)=1    ! 300623 Lei
        ishp(n+2)=1    ! 300623 Lei
        do i=1,3
        rpi=pyr(1)*v(ii,i)
        v(n+1,i)=rpi
        v(n+2,i)=v(ii,i)-rpi
        enddo
        n=n+2
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coend(ii)   ! 160110
c       copy parton ii to the end of 'pyjets'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        n=n+1
        do i3=1,5
        k(n,i3)=k(ii,i3)
        p(n,i3)=p(ii,i3)
        v(n,i3)=v(ii,i3)
        enddo
        taup(n)=taup(ii)   ! 300623 Lei
        ishp(n)=ishp(ii)   ! 300623 Lei
c        ndiq(n)=ndiq(ii)
c        do j1=1,idi
c        kk=ifcom(j1)
c        jj=npt(j1)
c        if(kk.eq.ii)ifcom(j1)=n
c        if(jj.eq.ii)npt(j1)=n
c        enddo
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine adjst(ik1,ik1str,ik1sa,ik1sv,i_g)   ! 160110
c300623 Added i_g for pure the gluon string "g-g-g".   ! 300623 Lei
c       finds order number etc. of ik1 string    ! 270620
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620
        i_g=0    ! 300623 Lei if ik1 does not belong to string
        ik1str=0   ! if ik1 not belong to string
        ik1sa=0   ! if ik1 not belong to string
        ik1sv=0   ! if ik1 not belong to string
        if(ik1.le.nbe)then   ! nbe is giving before call break
        do i1=1,nstr0   ! nstr0 is giving after call break
        i1a=nstr1a(i1)   ! i1a: line number of 'A' of i1-th string
        i1v=nstr1v(i1)   ! i1v: line number of 'V' of i1-th string
        kfa=kbe(i1a,2)   ! 300623 Lei kfa: KF code of 'A'
        kfv=kbe(i1v,2)   ! 300623 Lei kfa: KF code of 'V'
        if(kfa.eq.21 .AND. kfv.eq.21) i_g=1   ! 300623 Lei ik1 belongs to string
        if(ik1.ge.i1a .and. ik1.lt.i1v+1)then
        ik1str=i1   ! order number of string to which ik1 belongs
        ik1sa=i1a   ! line number of 'A' of above string
        ik1sv=i1v   ! line number of 'V' of above string
        goto 100
        endif
        enddo
        endif
        if(ik1.gt.nbe .and. nstr1.gt.nstr0)then
        do i1=nstr0+1,nstr1
        i1a=nstr1a(i1)
        i1v=nstr1v(i1)
        kfa=kbe(i1a,2)   ! 300623 Lei kfa: KF code of 'A'
        kfv=kbe(i1v,2)   ! 300623 Lei kfa: KF code of 'V'
        if(kfa.eq.21 .AND. kfv.eq.21) i_g=1   ! 300623 Lei ik1 belongs to string
        if(ik1.ge.i1a .and. ik1.lt.i1v+1)then
        ik1str=i1
        ik1sa=i1a
        ik1sv=i1v
        goto 100
        endif
        enddo
        endif
100     return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine adjdi(ii,idi1,idway)   ! 160110
c       does ii is component of a diquark
c       idway=0: if ii is not first component or it's partner of a diquark
c       idway=1: if ii is first component of a diquark
c       idway=2: if ii is partner of a diquark
c       idi1: order # of diquark to which ii belong
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
c       ifcom(i): line number (in 'sbe') of first component of i-th diquark
        idway=0 
        idi1=0
        do i1=1,idi
        i2=ifcom(i1)
        if(ii.eq.i2)then
c       ii is first component of idi1-th diquark
        idi1=i1
        idway=1
        goto 100
        endif
        enddo

        do i1=1,idi
        ii1=npt(i1)   ! line number of partner of i1-th diquark
        if(ii.eq.ii1)then   ! ii is partner of i1-th diquark
c       ii is partner of idi1-th diquark
        idi1=i1   ! ii is partner of idi1-th diquark
        idway=2
        goto 100
        endif
        enddo

100     return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine diqmov(idi1)   ! 160110
c       removes diquark idi1 out of list
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        if(idi.gt.0)then
        do i1=1,nbe
        if(ndiq(i1).eq.idi1)ndiq(i1)=0
c270620 ndiq(j): = 0 if j is quark (antiquark)
c                = order # of diquark if j is diquark (anti-diquark)
c270620       j: line number in 'sbe'
        enddo
        if(idi1.eq.idi)idi=idi-1
        if(idi1.lt.idi)then
c270620 updates diquark list
        do i1=idi1+1,idi
        i2=i1-1
        ifcom(i2)=ifcom(i1)
c270620 ifcom(idi): line number of first component of idi-th diquark
        npt(i2)=npt(i1)
        enddo
        idi=idi-1
        endif
        endif
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine parmov(ii,ik1,ik2,lmn)   ! 280620 070720
c       if 'ii' is not equal to ik1 (ik2) and lmn=4 or 6, moves parton 'ii' 
c        to 'trs' as well as updates lists
c       otherwise, throws away parton ii (ii is equal to ik1 (ik2)), 
c        i.e. updates lists
c       updates 'pyjets' one step downward from ii+1 to n
c       updates 'sbe' one step downward from ii+1 to nbe if ii.le.nbe
c       updates diquark list
c       updates string list
c       updates collision time list
c       ik1 & ik2 are the line number (in parton list) of colliding pair
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)   ! 051122
        PARAMETER (MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        common/sa24/adj1(40),nnstop,non24,zstop   ! 070720
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5)

        adj12=adj1(12)   ! 070720
c       moves parton ii (not equal to ik1 or ik2) to 'trs'
        if(adj12.eq.0 .and. (lmn.eq.4 .or. lmn.eq.6) .and. 
     c   (ii.ne.ik1 .or. ii.ne.ik2))then
        ntrs=ntrs+1
        ktrs(ntrs,1)=3
        ktrs(ntrs,2)=k(ii,2)
        do j=1,4
        ptrs(ntrs,j)=p(ii,j)
        vtrs(ntrs,j)=v(ii,j)
        enddo
        ptrs(ntrs,5)=p(ii,5)
        vtrs(ntrs,5)=0.
        endif

c       throws away parton ii, i.e. updates lists
c       updates diquark list
        do j1=1,idi
        jj=ifcom(j1)
        kk=npt(j1)
        if(jj.ge.ii)ifcom(j1)=jj-1
        if(kk.ge.ii)npt(j1)=kk-1
        enddo
        do i2=ii+1,n
        i3=i2-1
        ndiq(i3)=ndiq(i2)
        enddo

c       updates particle list 'pyjets'
        if(ii.eq.n)n=n-1
        if(ii.lt.n)then
        do i1=ii+1,n
        i2=i1-1
        do i3=1,5
        k(i2,i3)=k(i1,i3)
        p(i2,i3)=p(i1,i3)
        v(i2,i3)=v(i1,i3)
        enddo
        ishp(i2)=ishp(i1)   ! 300623 Lei
        taup(i2)=taup(i1)   ! 300623 Lei
        enddo
        n=n-1
        endif
c       if(ik1.gt.ii)ik1=ik1-1   ! 280620
c       if(ik2.gt.ii)ik2=ik2-1   ! 280620

c       updates particle list 'sbe'
        if(ii.eq.nbe)nbe=nbe-1
        if(ii.lt.nbe)then
        do i1=ii+1,nbe
        i2=i1-1
        do i3=1,5
        kbe(i2,i3)=kbe(i1,i3)
        pbe(i2,i3)=pbe(i1,i3)
        vbe(i2,i3)=vbe(i1,i3)
        enddo
        enddo
        nbe=nbe-1
        endif

c       updates string list
        do i1=1,nstr1
        nstraa=nstr1a(i1)
        if(nstraa.ge.ii)nstr1a(i1)=nstraa-1
        nstrvv=nstr1v(i1)
        if(nstrvv.ge.ii)nstr1v(i1)=nstrvv-1
        enddo

c       updates the values of lc(1-2,m) if which is .ge. ii
        do m=1,icol
        lc1=lc(1,m)
        if(lc1.ge.ii)lc(1,m)=lc1-1
        write(2,*)"parmov, m, lc1=", m, lc(1,m)   !Lei
        lc2=lc(2,m)
        if(lc2.ge.ii)lc(2,m)=lc2-1
        write(2,*)"parmov, m, lc2=", m, lc(2,m)   !Lei
        enddo

        write(3,*)"Af. parmov ------------------------------------"
        write(3,*)"icol, lmn=",icol,lmn   !Lei
        do iak=1,icol,1   !Lei
        write(3,*)"i, lc1, lc2, tc1, tc2=",
     &  iak, lc(1,iak), lc(2,iak), tc(1,iak), tc(2,iak)    !Lei
        end do

c300623 Lei
c       Removing colliding pairs with lc=0. This case happens when ii 
c        is the 1-st one in PYJETS, especially when removing string 
c        after inelastic collisions.
        jcol = 0
c       Loops over colliding pairs.
        do i=1,icol,1
            iic = lc(1,i)
            jjc = lc(2,i)
c           Throw away the pairs with lc=0.
            if(iic.ne.0 .AND. jjc.ne.0)then
                jcol = jcol + 1
                tc(1,jcol) = tc(1,i)
                tc(2,jcol) = tc(2,i)
                lc(1,jcol) = lc(1,i)
                lc(2,jcol) = lc(2,i)
            end if
        end do
        icol = jcol
c300623 Lei

        write(3,*)"Af. parmov remove ----------------------------------"
        write(3,*)"icol, lmn=",icol,lmn   !Lei
        do iak=1,icol,1   !Lei
        write(3,*)"i, lc1, lc2, tc1, tc2=",
     &  iak, lc(1,iak), lc(2,iak), tc(1,iak), tc(2,iak)    !Lei
        end do

        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine strmov(istr1,istr1a,istr1v,ik1,ik2,lmn)  ! 280620 070720
c       moves all conponents of 'istr1'-th string out
c       istr1: order number of string (in string list) to be moved
c       istr1a: line number (in parton list) of first component of above string 
c       istr1v: line number (in parton list) of last component of above string
c       ik1 & ik2 are line # of colliding pair in parton list 'pyjets'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)   ! 051122
        PARAMETER (MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620
c       moves components of istr1-th string out of parton list
        write(3,*)"****************************************************"
        write(3,*)"in strmov, parmov, istr1, istr1v, istr1a =",
     &               istr1, istr1v, istr1a   !Lei
        do ii=istr1v,istr1a,-1
        call parmov(ii,ik1,ik2,lmn)   ! 070720
        enddo
c       moves string istr1-th out of string list
        if(istr1.eq.nstr1)nstr1=nstr1-1
        if(istr1.lt.nstr1)then
        jj=istr1v-istr1a+1
        do ii=istr1+1,nstr1
c300623 if(ii.ge.istr1)then
c300623 nstr1v(ii-jj)=nstr1v(ii)   ! 300623 Lei
c300623 nstr1a(ii-jj)=nstr1a(ii)   ! 300623 Lei
        nstr1v(ii-1)=nstr1v(ii)   ! 300623 Lei
        nstr1a(ii-1)=nstr1a(ii)   ! 300623 Lei
c300623 endif
        enddo
        nstr1=nstr1-1
        endif
        return
        end




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine copl_p(tt)
c300623 Lei
c       calculate position of center of mass of the non-freeze-out system 
c       distance of a particle from this cms is used to checke whether
c        it freezes out or not 
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER(KSZJ=80000)
        common/PYJETS/N,NON1,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        do ii=1,3
        coor(ii)=0.
        enddo
        samass=0.
        do 110 ii=1,n
        if(ishp(ii).eq.0) goto 110
        kf=k(ii,2)
        if(kf.eq.88) goto 110   ! Exludes junctions.
        amass=pmas(pycomp(kf),1)
        samass=samass+amass
        do 100 jj=1,3
        coor(jj)=coor(jj)+amass*v(ii,jj)
100     continue
110     continue
        do ii=1,3
        coor(ii)=coor(ii)/max(0.33,samass)
        enddo
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine his_p(t1,rnt,rnp,istop)
c300623 Lei
c       classical Newton motion in AA/NN CMS/Lab. system
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)   ! 300623 Lei
        common/PYJETS/nsa,non1,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)   ! 300623 Lei
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei
        common/sa24/adj1(40),nnstop,non24,zstop
        common/collist/lc(2,mclis),tc(2,mclis),icol
        istop=1   ! 231104
        in=0   ! 231104
        r0=PARAM(10)*dmax1(rnt,rnp)   ! 300623 Lei
        do 200 i=1,nsa
c       if(t1.le.taup(i))goto 200
c       do move particles which have not produced
        if(ishp(i).eq.1) goto 10
        if(ksa(i,2).eq.88) goto 200   ! 300623 Lei
        in=in+1
        goto 200
10      aa=0.
        pp4=psa(i,4)
        do j=1,3
        vp=psa(i,j)/pp4
        vsa(i,j)=vsa(i,j)+vp*(t1-vsa(i,4))
        aa=aa+(vsa(i,j)-coor(j))**2
        enddo
c100505 vsa(i,4)=t1
        aa=sqrt(aa)
        if(aa.lt.r0) goto 300
c       if freeze-out already, deduct the distance between the last and 
c        current collisions
        do j=1,3
        vp=psa(i,j)/pp4
        vsa(i,j)=vsa(i,j)-vp*(t1-vsa(i,4))
        enddo
        ishp(i)=0
        do il=1,icol
        if(lc(1,il).eq.i) tc(1,il)=0.
        if(lc(2,il).eq.i) tc(2,il)=0.
        enddo
        goto 200
300     vsa(i,4)=t1   ! 100505
200     continue
        if(in.eq.nsa) return
        istop=0
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rsfilt_p(l,l1,iflag)
c300623 Lei
c       subroutine rsfilt_p plays the role of first range filter 
c       subroutine intdis plays the role of second range filter
c       collision pairs not interested can not filter through both of rsfilt 
c        and intdis
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)   ! 300623 Lei

        kl=k(l,2)
        kl1=k(l1,2)
        klab=iabs(kl)
        kl1ab=iabs(kl1)
        if(l.eq.l1) goto 10
        if(ishp(l).eq.0.or.ishp(l1).eq.0) goto 10

        if( ((klab.ge.1 .AND. klab.le.6) .OR. klab.eq.21) .AND.
     &      ((kl1ab.ge.1 .AND. kl1ab.le.6) .OR. kl1ab.eq.21) ) goto 11

        goto 10

11      iflag=1
10      continue
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updpli_p(n00,time)
c300623 Lei
c       Update collision time list for new partons when removing string 
c        and diqaruk after inelastic parton collisions.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (KSZJ=80000,MCLIS=280000)   ! 051122
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 051122
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/scatt/pi(4),pj(4),ic,jc,n0   ! 051122
        common/collist/lc(2,mclis),tc(2,mclis),icol
c-----------------------------------------------------------------------
        dddt = adj1(19)   ! 161104
        ! j = 0
c       Loop over new ii (greater than n00) and old jj partons (smaller than n00)
        icol = icol + 1
        do 100 ii=n00+1,n,1
            i1 = ii
            kfi = iabs(k(i,2))
            if(kfi.gt.6 .and. kfi.ne.21) goto 100
c           Consider d,u,s,c,b,t, their antiquarks, and gluon.
            do 200 jj=1,n00,1
                j1 = jj
                if(i1.eq.i1) goto 200

                iflag=0
                call rsfilt_p(j1,i1,iflag)
                if(iflag.eq.0) goto 200

                tc(1,icol) = 0.
                tc(2,icol) = 0.
                call tcolij_par(i1,j1,time,icol)
                tc1 = tc(1,icol)
                tc2 = tc(2,icol)
                tcicol = tc1
                if(tc1.gt.tc2) tcicol = tc2
                if(tcicol.gt.0.)then
                    tci = tcicol
                    do j2=1,icol-1
                        tcj = tc(1,j2)
                        tck = tc(2,j2)
                        if(tcj.gt.tck) tcj = tck
                        if(dabs(tcj-tci).lt.dddt) goto 200
                    enddo
                    icol = icol + 1
                endif
200         enddo
100     enddo
450     continue
        icol = icol - 1
        n0 = n
        return
        end



cccccccccccccccccccccccccccccccc end ccccccccccccccccccccccccccccc
