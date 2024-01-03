	subroutine parcas(time_par,jjj,iijk,win,nap,rnt,rnp) 
c       deals with parton cascade (partonic rescattering)
c	input messages are in 'parlist' ('pyjets' to 'parlist' in main.f)
c	working block is 'parlist'
c180520	output messages are in 'parlist' (in 'pyjets' either)  
c	writen by Ben-Hao Sa 19/11/2002 
c160110 iiii: number of run
c       jjj: jjj-th loop interaction in a event   ! 180520
	parameter (mplis=80000,msca=20000)
c181018	parameter(kszj=80000)   ! 160110
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
c181018	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 160110
	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)  ! 240503
	common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
	common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp   ! 120505
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 080520
      common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
c       iprl: current total number of partons in parton list
c       rp  : space and time coordinates of particle
c       pp  : momentum of particle
c       tp  : current time of particle
c       taup: formation time of particle
c       idp : flavor of particle
c       rmp : rest mass of particle
c       ep  : total energy of particle
c       vp  : velocity of particle

	parameter (mclis=280000)
	common/collist/lc(2,mclis),tc(2,mclis),icol

c       icol : current total number of collision pairs in collision time list
c       lc   : line number (in [1,iprl]) of colliding pair, eg.
c	 lc(1,100): line number of first particle of 100-th colliding pair
c       tc   : the collision time of colliding pair

	common/scatt/pi(4),pj(4),ic,jc,iprl0
c       ic,jc: line number of colliding particles
c       iprl0: the iprl before current collision
c	pi,pj: four momentum of colliding particles 
	common/work7/reac(9),crose(9)
      common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
      common/sa24/adj1(40),nnstop,non24,zstop
	common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)   ! 201104
c       ithroq_p : total # of quarks thrown away
c       ithrob_p : total # of antiquarks thrown away
c       throe_p : total momentum and energy of the partons thrown away
c       ichh_p : total charge of the partons thrown away
	dimension b(3),ppr(mplis,5),peo(5),bpp(20)
	dimension rpo(4,mplis),ppo(4,mplis)   ! 120603 120505
c151203	iijk=0   ! 120603
	tl0=adj1(24)
c       tl0: cut off virtuality of time-like branching, i. e. Mu0**2
	if(adj1(1).eq.0.)return   ! 290505
        time=time_par   ! 280910
c241104
	dpmax=adj1(27)
	drmax=adj1(28)
c241104
	adj112=adj1(12)
	adj136=adj1(36)   ! 120505
	adj137=adj1(37)   ! 120505
      call reset_eve
c201104
        ithroq_p=0
        ithrob_p=0
        ich_p=0
        do i=1,4
        throe_p(i)=0.
        enddo
	do i1=1,iprl
	rp(4,i1)=0.
	taup(i1)=0.
	enddo
c241104
c	throw away parton if its modular of three momentum > dpmax or energy 
c        > dpmax or modular of spatial coordinate > drmax
	i11=1   ! 0 Tan 050605
c200	do 300 i1=i11+1,iprl   ! 1 050605 Tan
200	do 300 i1=i11,iprl   ! 1 050605 Tan
	ppp1=pp(1,i1)
	ppp2=pp(2,i1)
	ppp3=pp(3,i1)
	ppp4=pp(4,i1)
	rrr1=rp(1,i1)
	rrr2=rp(2,i1)
	rrr3=rp(3,i1)
	pppm=ppp1*ppp1+ppp2*ppp2+ppp3*ppp3
	if(pppm.gt.1.e20)pppm=1.e20
	if(pppm.lt.1.e-20)pppm=1.e-20
	pppm=dsqrt(pppm)
	rrrm=rrr1*rrr1+rrr2*rrr2+rrr3*rrr3
	if(rrrm.gt.1.e20)rrrm=1.e20
        if(rrrm.lt.1.e-20)rrrm=1.e-20
	rrrm=dsqrt(rrrm)
	if(pppm.le.dpmax.and.ppp4.le.dpmax.and.rrrm.le.drmax)goto 300
	if(adj112.ne.0.)then   ! 140705
	kf=idp(i1)
	if(kf.gt.0.or.kf.eq.21)ithroq_p=ithroq_p+1
	if(kf.lt.0)ithrob_p=ithrob_p+1
	if(kf.ne.21)ich_p=ich_p+pychge(kf)
	do k1=1,4
	throe_p(k1)=throe_p(k1)+pp(k1,i1)
	enddo
	if(i1.eq.iprl)then
	iprl=iprl-1
	goto 500
	endif
	do k1=i1+1,iprl
	do k2=1,3
	rp(k2,k1-1)=rp(k2,k1)
	pp(k2,k1-1)=pp(k2,k1)
	vp(k2,k1-1)=vp(k2,k1)
	enddo
	rp(4,k1-1)=rp(4,k1)
	pp(4,k1-1)=pp(4,k1)
	idp(k1-1)=idp(k1)
	taup(k1-1)=taup(k1)
	rmp(k1-1)=rmp(k1)
	enddo
	iprl=iprl-1
	if(iprl.eq.0)return
	i11=i1
	goto 200   ! 140705
	endif   !
	if(adj112.eq.0.)then   !!
c	in case of string fragmentation should not throw away that parton
        if(pppm.gt.dpmax.and.ppp4.gt.dpmax)then   ! 311007
	cita=2*pyr(1)-1.
        fi=2.*3.1416*pyr(1)
        sita=dsqrt(1.-cita**2)
	pp1i=dpmax*sita*dcos(fi)
	pp2i=dpmax*sita*dsin(fi)
	pp3i=dpmax*cita
        pp(1,i1)=pp1i
        pp(2,i1)=pp2i
        pp(3,i1)=pp3i
	pp4i=rmp(i1)*rmp(i1)+pp1i*pp1i+pp2i*pp2i+pp3i*pp3i
	if(pp4i.gt.1.e20)pp4i=1.e20
	if(pp4i.lt.1.e-20)pp4i=1.e-20
	pp(4,i1)=dsqrt(pp4i)
        endif   ! 311007
        if(rrrm.gt.drmax)then   ! 311007
	cita=2*pyr(1)-1.
        fi=2.*3.1416*pyr(1)
        sita=dsqrt(1.-cita**2)
	rp(1,i1)=drmax*sita*dcos(fi)
	rp(2,i1)=drmax*sita*dsin(fi)
	rp(3,i1)=drmax**cita  
        endif   ! 311007
	goto 300   ! 050605
	endif   !!
c140705	goto 200
300	enddo   ! 1
500	if(ithroq_p.ne.0.or.ithrob_p.ne.0)then
c140705	write(9,*)'iiii,ithroq_p,ithrob_p=',iiii,ithroq_p,ithrob_p   ! sa
c140705	write(9,*)'throe_p=',throe_p   ! sa
	endif
c241104
c201104
c	step 1
c       create the parton-parton (initial) collision time list 
	call ctlcre_par(iijk)   ! 290803 
	if(iijk.eq.2)return   ! initial collis. list is empty 151203
c290803
c151203	if(iijk.eq.1)then
c151203	iiii=iiii-1
c151203	return
c151203	endif
c290803

c160110 the loop over parton-parton collisions within an event
	jjj=0   
	icolo=icol   ! 120603
c	statistic of the number of loops in parton cascade within an event
24      jjj=jjj+1   ! 160110   
c---------------------------------------------------------------------
	if(jjj.gt.100*icolo)then
	write(9,*)'infinite loop may have happened in'
        write(9,*)'parcas iiii,jjj,icolo=',iiii,jjj,icolo
	iiii=iiii-1
	iijk=1
	return
	endif
	iprl0=iprl

c	step 2
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
27	continue
	time=tcp

c	step 3
c	perform classical Newton motion for ic & jc partons
c120603	do j=1,3
c120603	rp(j,ic)=rp(j,ic)+vp(j,ic)*(tc(1,icp)-rp(4,ic))
c120603	rp(j,jc)=rp(j,jc)+vp(j,jc)*(tc(2,icp)-rp(4,jc))
c120603	enddo
c120603	rp(4,ic)=tc(1,icp)
c120603	rp(4,jc)=tc(2,icp)
c       perform classical Newton motion
	do i=1,iprl
	do j=1,3
	rpo(j,i)=rp(j,i)
	rp(j,i)=rp(j,i)+vp(j,i)*(tcp-rp(4,i))
	enddo
	rpo(4,i)=rp(4,i)
	rp(4,i)=tcp
	enddo
c120603
c	consider parton energy loss phenomenologically
c120505
	if(adj136.eq.1)then
	bmax=rnt+rnp
	bmax2=bmax*bmax
	bp2=bp*bp
c	energy loss per unit distance
	dell=adj137*(nap/197.)**0.75*(1.-bp2/bmax2)
     c	 *(win*win/40000)**0.25
c	write(9,*)'adj136,adj137,win=',adj136,adj137,win   ! sa
c	write(9,*)'nap,rnt,rnp,dell=',nap,rnt,rnp,dell   ! sa
	call eloss(dell,rpo,ppo)  
	endif 
c120505
c	step 4
c       performs parton-parton collision & updates particle list 
c       if lmn=4,6,& 7 updates 'parlist','pyjets','sbe',diquark list,string 
c        list, & lc(1-2,m) either 
	kkk=0   ! 120603
	iway=0   ! 120505
	call collis(ic,jc,kf3,kf4,tcp,jjj,kkk,iway,icnew,jcnew   
     c   ,lmn,time)   ! 120603 120505 160110
c120603
c120603
503     format(5(1x,e12.3))
 
c	step 5

c060620 update collision list
        if(lmn.eq.4.or.lmn.eq.6.or.lmn.eq.7)goto 26 ! 070720
c       in above case 'update collision list' is executed in "collis'    
        call update_ctl(ic,jc,kf3,kf4,iparres,lmn,iway,tcp)   

c	'new' method
c1	if(adj136.eq.0)call update_ctl(tcp,iway)   
c1	if(adj136.eq.1 .and.iway.eq.1)call update_ctl(tcp,iway)
c	if collision does not happen remove the collision pairs (one of
c	 company is colliding particle) from collision time list only
c1      if(adj136.eq.1 .and.iway.eq.0)then
c	create collision time list completely
c1	icol=0
c1	do i=1,mclis
c1	lc(1,i)=0
c1	lc(2,i)=0
c1	tc(1,i)=0.
c1	tc(2,i)=0.
c1	enddo
c1	call ctlcre_para(iijk,time)
c1	endif   
c120505
c	goto 25   ! it is actived temporally
26	continue   ! 120603
	goto 24   ! the loop over collisions within an event
25      continue
	time_par=time
c       time_par: is the time lasted in parton cascade hereafter
c       write(9,*)'af parton cascade time=',time   ! sa
c250803
	do i=1,9
	reaci=reac(i)
	if(reaci.gt.0.)crose(i)=crose(i)/reaci
	enddo
c250803
c	throw away parton if its modular of three momentum > dpmax or energy 
c	 > dpmax or modular of spatial coordinate > drmax
	i11=1   
201	do 301 i1=i11,iprl   ! 1 
	ppp1=pp(1,i1)
	ppp2=pp(2,i1)
	ppp3=pp(3,i1)
	ppp4=pp(4,i1)
	rrr1=rp(1,i1)
	rrr2=rp(2,i1)
	rrr3=rp(3,i1)
	pppm=ppp1*ppp1+ppp2*ppp2+ppp3*ppp3
	if(pppm.gt.1.e20)pppm=1.e20
	if(pppm.lt.1.e-20)pppm=1.e-20
	pppm=dsqrt(pppm)
	rrrm=rrr1*rrr1+rrr2*rrr2+rrr3*rrr3
	if(rrrm.gt.1.e20)rrrm=1.e20
        if(rrrm.lt.1.e-20)rrrm=1.e-20
	rrrm=dsqrt(rrrm)
	if(pppm.le.dpmax.and.ppp4.le.dpmax.and.rrrm.le.drmax)goto 301
	if(adj112.ne.0.)then   !
	kf=idp(i1)
	if(kf.gt.0.or.kf.eq.21)ithroq_p=ithroq_p+1
	if(kf.lt.0)ithrob_p=ithrob_p+1
	if(kf.ne.21)ich_p=ich_p+pychge(kf)
	do k1=1,4
	throe_p(k1)=throe_p(k1)+pp(k1,i1)
	enddo
	if(i1.eq.iprl)then
	iprl=iprl-1
	goto 501
	endif
	do k1=i1+1,iprl
	do k2=1,3
	rp(k2,k1-1)=rp(k2,k1)
	pp(k2,k1-1)=pp(k2,k1)
	vp(k2,k1-1)=vp(k2,k1)
	enddo
	rp(4,k1-1)=rp(4,k1)
	pp(4,k1-1)=pp(4,k1)
	idp(k1-1)=idp(k1)
	taup(k1-1)=taup(k1)
	rmp(k1-1)=rmp(k1)
	enddo
	iprl=iprl-1
	if(iprl.eq.0)return
	i11=i1
	goto 201   ! 140705
	endif   !
	if(adj112.eq.0.)then   !!
c	in case of string fragmentation should not throw away that parton
        if(pppm.gt.dpmax.and.ppp4.gt.dpmax)then   ! 311007
	cita=2*pyr(1)-1.
        fi=2.*3.1416*pyr(1)
        sita=dsqrt(1.-cita**2)
	pp1i=dpmax*sita*dcos(fi)
	pp2i=dpmax*sita*dsin(fi)
	pp3i=dpmax*cita
        pp(1,i1)=pp1i
        pp(2,i1)=pp2i
        pp(3,i1)=pp3i
	pp4i=rmp(i1)*rmp(i1)+pp1i*pp1i+pp2i*pp2i+pp3i*pp3i
	if(pp4i.gt.1.e20)pp4i=1.e20
	if(pp4i.lt.1.e-20)pp4i=1.e-20
	pp(4,i1)=dsqrt(pp4i)
        endif   ! 311007
        if(rrrm.gt.drmax)then   ! 311007
	cita=2*pyr(1)-1.
        fi=2.*3.1416*pyr(1)
        sita=dsqrt(1.-cita**2)
	rp(1,i1)=drmax*sita*dcos(fi)
	rp(2,i1)=drmax*sita*dsin(fi)
	rp(3,i1)=drmax**cita  
        endif   ! 311007
	goto 301   ! 050605
	endif   !!
c140705	goto 201
301	enddo   ! 1
501	if(ithroq_p.ne.0.or.ithrob_p.ne.0)then
c140705	write(9,*)'iiii,ithroq_p,ithrob_p=',iiii,ithroq_p,ithrob_p   ! sa
c140705	write(9,*)'throe_p=',throe_p   ! sa
	endif
c140705
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine reset_eve
c       initiate the collision time list
	parameter (mplis=80000,mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	common/collist/lc(2,mclis),tc(2,mclis),icol
	common/scatt/pi(4),pj(4),ic,jc,iprl0
	common/work7/reac(9),crose(9)
c       reac and crose: the arraies to account for the number and
c        the value of cross section for 2->2 partonic processes
	ic=0
	jc=0
	do i=1,mclis
	lc(1,i)=0
	lc(2,i)=0
	tc(1,i)=0.
	tc(2,i)=0.
	enddo
	icol=0
	do i=1,9
        reac(i)=0.
        crose(i)=0.
        enddo
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ctlcre_par(iijk)      
c       create the initial collision list 
	parameter (mplis=80000,mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/collist/lc(2,mclis),tc(2,mclis),icol
	dddt=adj1(19)   ! 161104
	icol=1
	time=0.   ! 111599
c	every process (eg. parton casecade) starts from time equal to 0
	dminf=100.   ! 111599
	ijk=0   ! 010601
	do 100 i=2,iprl   ! upper diagonal 151203   
c	iprl00: total # of partons in projectile nucleus
c080603
	kfi=iabs(idp(i))
	if(kfi.gt.3 .and. kfi.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
c       if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
c080603
	do 200 j=1,i-1   ! upper diagonal 151203 
c080603
        kfj=iabs(idp(j))
	if(kfj.gt.3 .and. kfj.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
c       if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
c080603
	if(icol.gt.mclis) then
	write(9,*)'icol over limit iprl,icol=',iprl,icol   ! sa
	stop 7777
	endif
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
c	write(9,*)'in ctlcre, i,j,tc1,tc2,tcicol',i,j,tc1,tc2,tcicol! sa
c080304	if(tcicol.gt.0.0) icol=icol+1
c080304
	if(tcicol.gt.0.0)then   ! 121104
        tci=tcicol
        do j1=1,icol-1
        tcj=tc(1,j1)
        tck=tc(2,j1)
        if(tcj.gt.tck)tcj=tck
c161104	if(ddt.eq.0.)dddt=ddt+0.03
c161104	if(ddt.ne.0.)dddt=ddt*300
        if(dabs(tcj-tci).lt.dddt)goto 200
        enddo
        icol=icol+1
	endif     ! 121104
c080304   
200	enddo
100	enddo
	if(tcicol.eq.0.) icol=icol-1
	if(icol.eq.0)then    
c290803
	iijk=2   ! 1 151203
	return
c290803
c       at least one collision should occur, which has the smallest
c        'least approaching distance', that is guaranteed by the variable
c        'dminf'
c290803	icol=1
c290803 lc(1,icol)=iff
c290803 lc(2,icol)=jf
c290803	tc(1,icol)=0.02
c290803	tc(2,icol)=0.02
	endif
	do i1=icol+1,mclis
	lc(1,i1)=0
	lc(2,i1)=0
	tc(1,i1)=0.
	tc(2,i1)=0.
	enddo
c	write(9,*)'in ctlcre,ijk=',ijk   ! sa
c	write(9,*)'in ctlcre,icol=',icol   ! sa
c	write(9,*)'i=',(lc(1,i1),i1=1,icol)   ! sa
c	write(9,*)'j=',(lc(2,i1),i1=1,icol)   ! sa
c	write(9,*)'ti=',(tc(1,i1),i1=1,icol)   ! sa
c	write(9,*)'tj=',(tc(2,i1),i1=1,icol)   ! sa
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find_par(icp,tcp)
c       find out the binary collision (icp) with minimum colli. time (tcp)
	parameter (mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	common/scatt/pi(4),pj(4),ic,jc,iprl0
	common/collist/lc(2,mclis),tc(2,mclis),icol
c	write(9,*)'in find, icol=',icol ! sa
c        do i1=1,icol   ! sa
c        write(9,*)lc(1,i1),lc(2,i1),tc(1,i1),tc(2,i1)   ! sa
c        enddo   ! sa
	icp=0
	tcp=10000.
	do i=1,icol
	ia=(iabs(lc(1,i)-ic)+iabs(lc(2,i)-jc))*
     c  (iabs(lc(1,i)-jc)+iabs(lc(2,i)-ic)) !it play role after first colli.
	if(ia.eq.0) goto 241
	tc1=tc(1,i)
	tc2=tc(2,i)
	tci=tc1
	if(tc1.gt.tc2)tci=tc2
c	write(9,*)'tci=',tci   ! sa
c	tci=amax1(tc(1,i),tc(2,i))  ! alternative choice.
	if(tci.le.0.) goto 241
	if(tcp.lt.tci)  goto 241
	icp=i
	tcp=tci
241     continue
	enddo
c	write(9,*)'in find, icp,tcp=',icp,tcp   ! sa
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update_ctl(ik1,ik2,kf3,kf4,iparres,lmn,iway,time)   ! 120505
c020512 update collision time list for both of w/o & w/ inelastic
c020512  parton-parton scattering 
c       ik1,ik2: line number of the colliding pair in parton list
c       kf3 and kf4: kf code of the collided pair 
c230520 if iway=1 the collision does not happen
	parameter (mplis=80000,mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/scatt/pi(4),pj(4),ic,jc,iprl0
	common/collist/lc(2,mclis),tc(2,mclis),icol
c-----------------------------------------------------------------------
c	write(9,*)'in upda ic,jc,iprl,icol=',ic,jc,iprl,icol   ! sa
	dddt=adj1(19)   ! 161104
	j=0
c       loop over old colliding pairs
	if(icol.eq.0)goto 370
	do i=1,icol
	i1=lc(1,i)
	j1=lc(2,i)
	iii=(i1-ic)*(j1-ic)
	jjj=(i1-jc)*(j1-jc)
c       throw away the pairs composed of ic and/or jc
	if(iii.eq.0.or.jjj.eq.0) goto 400
380     continue
	tc1=tc(1,i)
	tc2=tc(2,i)
	tci=tc1
	if(tc1.gt.tc2)tci=tc2
c080104	if(tci.le.time) goto 400
c161104	if(ddt.eq.0.)dddt=ddt+0.03
c161104	if(ddt.ne.0.)dddt=ddt*300
c       throw away the pairs with tc-time<=dddt (time accuracy)
	if((tci-time).le.dddt) goto 400   ! 080104
c       proceeds for serviver
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
370	icol=j+1
	do 100 i=1,iprl0
	i1=i
	if(i1.eq.ic) goto 100
	if(i1.eq.jc) goto 100
c080603
	kfi=iabs(idp(i))
	if(kfi.gt.3 .and. kfi.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
c       if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
c080603
	do 200 k=1,2
	if(k.eq.1)j1=ic
	if(k.eq.2)j1=jc
	tc(1,icol)=0.0
	tc(2,icol)=0.0
c	write(9,*)'be. call tcolij i1,j1=',i1,j1   ! sa
	call tcolij_par(i1,j1,time,icol)
c	write(9,*)'af. call tcolij i1,j1=',i1,j1   ! sa
	tc1=tc(1,icol)
	tc2=tc(2,icol)
	tcicol=tc1
	if(tc1.gt.tc2)tcicol=tc2
c020603
	if(tcicol.le.0.0)then
c	write(9,*)'i1,j1,tcicol=',i1,j1,tcicol   ! sa
	endif
c020603	
c141104	if(tcicol.gt.0.0) icol=icol+1
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
200	enddo
100	enddo
450     continue
	if(tcicol.le.0.0) icol=icol-1
	do i=icol+1,mclis
	lc(1,i)=0
	lc(2,i)=0
	tc(1,i)=0.
	tc(2,i)=0.
	enddo
	iprl0=iprl
c020603
c260620 write(9,*)'af. update_ctl iprl0,icol,ic,jc=',iprl0,icol,ic,jc   ! sa
c	do i1=1,icol   ! sa
c	write(9,*)lc(1,i1),lc(2,i1),tc(1,i1),tc(2,i1)   ! sa
c	enddo   ! sa
c020603
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update_ctlm(time,iway)   ! 120505 160110
c       a part of updating collision time list (throw away old collission 
c        pairs only) for inela. parton-parton scattering 7  ! 230520
	parameter (mplis=80000,mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/scatt/pi(4),pj(4),ic,jc,iprl0
	common/collist/lc(2,mclis),tc(2,mclis),icol
c-----------------------------------------------------------------------
c	write(9,*)'get in upda ic,jc,iprl,icol=',ic,jc,iprl,icol   ! sa
	dddt=adj1(19)   ! 161104
	j=0
c       loop over old colliding pairs
	do i=1,icol
	i1=lc(1,i)
	j1=lc(2,i)
	iii=(i1-ic)*(j1-ic)
	jjj=(i1-jc)*(j1-jc)
c       throw away the pairs composed of ic and/or jc
	if(iii.eq.0.or.jjj.eq.0) goto 400
380     continue
	tc1=tc(1,i)
	tc2=tc(2,i)
	tci=tc1
	if(tc1.gt.tc2)tci=tc2
c080104	if(tci.le.time) goto 400
c161104	if(ddt.eq.0.)dddt=ddt+0.03
c161104	if(ddt.ne.0.)dddt=ddt*300
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
	parameter (mplis=80000,mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/scatt/pi(4),pj(4),ic,jc,iprl0
	common/collist/lc(2,mclis),tc(2,mclis),icol
c-----------------------------------------------------------------------
        dddt=adj1(19)   ! 161104
c       loop over ic (jc) and old 'parlist' (i.e. construct colli. pair 
c        composed of partons one of which is ic (jc) and another one 
c        in old 'parlist')
        icol=icol+1
	do 100 i=1,iprl0   ! 
	i1=i
	if(i1.eq.ic) goto 100
	if(i1.eq.jc) goto 100
c080603
	kfi=iabs(idp(i))
	if(kfi.gt.3 .and. kfi.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
c       if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
c080603
	do 200 k1=1,2
	if(k1.eq.1)j1=ic
	if(k1.eq.2)j1=jc
	tc(1,icol)=0.0
	tc(2,icol)=0.0
c	write(9,*)'be. call tcolij i1,j1=',i1,j1   ! sa
	call tcolij_par(i1,j1,time,icol)
c	write(9,*)'af. call tcolij i1,j1=',i1,j1   ! sa
	tc1=tc(1,icol)
	tc2=tc(2,icol)
	tcicol=tc1
	if(tc1.gt.tc2)tcicol=tc2
c020603
	if(tcicol.le.0.0)then
c	write(9,*)'i1,j1,tcicol=',i1,j1,tcicol   ! sa
	endif
c020603	
c141104	if(tcicol.gt.0.0) icol=icol+1
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
200	enddo
100	enddo
450     continue
	if(tcicol.le.0.0) icol=icol-1
	do i=icol+1,mclis
	lc(1,i)=0
	lc(2,i)=0
	tc(1,i)=0.
	tc(2,i)=0.
	enddo
	iprl0=iprl
c020603
c	write(9,*)'in update_ctl af. loop iprl0,icol,ic,jc=',
c     c	 iprl0,icol,ic,jc   ! sa
c	do i1=1,icol   ! sa
c	write(9,*)lc(1,i1),lc(2,i1),tc(1,i1),tc(2,i1)   ! sa
c	enddo   ! sa
c020603
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coij_p(i,j,time,icp,dminf,iff,jf,kji)
c       calculate the collision time for construction
c	 of initial collision time list
	parameter (mplis=80000,mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/syspar_p/rsig1,pio,tcut
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	common/collist/lc(2,mclis),tc(2,mclis),icol
	common/scatt/pi(4),pj(4),ic,jc,iprl0
	dimension px(4),py(4),pij(4)
	dimension dr(3),db(3),vi(3),vj(3),pic(4),pjc(4),pxc(4),pyc(4)
	double precision b(3)
c	write(9,*)i,j,idp(i),idp(j)   ! sa
	pio=3.1416
	pi(4)=pp(4,i)
	pj(4)=pp(4,j)
	if(pi(4).lt.1.e-10)pi(4)=1.e-10   ! 031204
	if(pj(4).lt.1.e-10)pj(4)=1.e-10   ! 031204
	pij(4)=pi(4)+pj(4)
	pic(4)=pi(4)
        pjc(4)=pj(4)
	do k1=1,3
	pi(k1)=pp(k1,i)
	pj(k1)=pp(k1,j)
	pij(k1)=pi(k1)+pj(k1)
	pic(k1)=pi(k1)
        pjc(k1)=pj(k1)
	b(k1)=(pi(k1)+pj(k1))/pij(4)
	enddo
	rmi=rmp(i)
	rmj=rmp(j)
	eiej2=dot(pij,pij)
c230520 invariant mass
c	insert! energy cut
	if(eiej2.lt.0.) then
c	write(9,*)'pi=',pi    
c	write(9,*)'pj=',pj    
c	write(9,*)'pij=',pij    
c	write(9,*)'rmi,rmj=',rmi,rmj   
c	write(9,*)i,j,idp(i),idp(j)
c240503   
	kji=1   ! stop 2222
	return
c240503    
	endif
	do n1=1,4
	px(n1)=rp(n1,i)
	py(n1)=rp(n1,j)
	pxc(n1)=px(n1)
        pyc(n1)=py(n1)
	enddo
	ilo=0
	kf1=idp(i)
	kf2=idp(j)
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
c	write(9,*)'kf1,kf2,sig,eiej2=',kf1,kf2,sig,eiej2   ! sa
	if(sig.le.0.)return   ! 120603 250803
	if(ilo.eq.-2)then   ! 111999
c	write(9,*)'return ctlcre because ilo=-2, i,j=',i,j   ! sa
	return   ! added by Sa on 24/06/96
	endif   ! 111999
	rsig1=dsqrt(sig/pio)

	call lorntz(0,b,pic,pjc)
	call lorntz(0,b,pxc,pyc)
c	if(i.eq.3 .and. j.eq.9)then   ! sa
c	write(9,*)'after Lorentz'   ! sa
c	write(9,'3(f8.4)')(pic(i1),i1=1,3)   ! sa
c	write(9,'3(f8.4)')(pjc(i1),i1=1,3)   ! sa
c	endif   ! sa
	rb=0.
        bb=0.
        rr=0.
        rtai=0.
	do k=1,3
        vi(k)=pic(k)/pic(4)
        vj(k)=pjc(k)/pjc(4)
        enddo
        do k=1,3
        dr(k)=pxc(k)-pyc(k)-(vi(k)*pxc(4)-vj(k)*pyc(4))
        db(k)=vi(k)-vj(k)
        rb=rb+dr(k)*db(k)
        bb=db(k)**2+bb
        rr=rr+dr(k)*dr(k)
        enddo
c	if(i.eq.3 .and. j.eq.9)then   ! sa
c	write(9,*)'vi,vj,bb='   ! sa
c	write(9,'3(f8.4)')(vi(i1),i1=1,3)   ! sa
c	write(9,'3(f8.4)')(vj(i1),i1=1,3)   ! sa
c	write(9,'3(f8.4)')(db(i1),i1=1,3)   ! sa
c	endif   ! sa
        if(bb.le.1.e-10)then
c	write(9,*)'return ctlcre because bb too small, i,j',i,j    ! sa
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
c	write(9,*)'return ctlcre,i,j,dmin,rsig1=',i,j,dmin,rsig1   ! sa
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
100	lc(1,icp)=i
	lc(2,icp)=j
	tc(1,icp)=pxc(4)
	tc(2,icp)=pyc(4)
c	write(9,*)'in coij,iff,jf=',iff,jf   ! sa
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine tcolij_par(i,j,time,icp)
c       calculate the collision time of i and j
	parameter (mplis=80000,mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/syspar_p/rsig1,pio,tcut
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	common/collist/lc(2,mclis),tc(2,mclis),icol
	common/scatt/pi(4),pj(4),ic,jc,iprl0
	dimension px(4),py(4),dx(4),pij(4),pic(4),pjc(4),
     c	 pxc(4),pyc(4)
	dimension dr(3),db(3),vi(3),vj(3),rfi(4),rfj(4)
	double precision b(3)
	pi(4)=pp(4,i)
	pj(4)=pp(4,j)
	if(pi(4).lt.1.e-10)pi(4)=1.e-10   ! 031204
	if(pj(4).lt.1.e-10)pj(4)=1.e-10   ! 031204
	pij(4)=pi(4)+pj(4)
	pic(4)=pi(4)   
        pjc(4)=pj(4)   
	do k1=1,3
	pi(k1)=pp(k1,i)
	pj(k1)=pp(k1,j)
	pij(k1)=pi(k1)+pj(k1)
	pic(k1)=pi(k1)   
        pjc(k1)=pj(k1)   
        b(k1)=(pi(k1)+pj(k1))/(pi(4)+pj(4))   
	enddo
c	write(9,*)'in tcolij, i,j,idp(i),idp(j)=',i,j,idp(i),idp(j) ! sa
c	write(9,*)'pi=',(pi(i1),i1=1,3)   ! sa
c	write(9,*)'pj=',(pj(i1),i1=1,3)   ! sa
	rmi=rmp(i)
	rmj=rmp(j)
	eiej2=dot(pij,pij)
c230520 squared invariant mass
c	insert! energy cut
	if(eiej2.lt.0.) then
c	write(9,*)'eiej2 less than 0., give up that colli.,i,j=',i,j! sa
c	write(9,*)'pi=',pi
c	write(9,*)'pj=',pj
c	write(9,*)'pij=',pij
c	write(9,*)'rmi,rmj=',rmi,rmj
c	write(9,*)i,j,idp(i),idp(j)
c072200	stop 1111
	return   ! 072200
	endif
c081000	ecut=dsqrt(eiej2)-rmi-rmj
c081000	ecut0=0.02    ! P.Yang
c081000	if(ecut.le.ecut0) return
c Note! energy cut can be taken into account HERE to decide coll. pairs
c etc.  According to Y. Pang's opinions, May,1994, CCAST.(05/24/95) 
	do n=1,4
	px(n)=rp(n,i)
	py(n)=rp(n,j)
	pxc(n)=px(n)   
        pyc(n)=py(n)   
	enddo
	ilo=0
	kf1=idp(i)
	kf2=idp(j)
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
c	write(9,*)'kf1,kf2,sig,eiej2=',kf1,kf2,sig,eiej2
cc	if(il.eq.2)return   ! put 'cc' by Sa on 24/06/96
	if(ilo.eq.-2)return   ! added by Sa on 24/06/96
	if(sig.le.0.)return   ! 120603 250803 
	rsig1=dsqrt(sig/pio)

        do i1=1,3
        rfi(i1)=px(i1)+(taup(i)-time)*pi(i1)/pi(4)
        rfj(i1)=py(i1)+(taup(j)-time)*pj(i1)/pj(4)
	enddo
        rfi(4)=taup(i)
        rfj(4)=taup(j)
c	spatial coordinates of colliding particles at formation 
c	 moment in Lab. frame
        call lorntz(0,b,rfi,rfj)
        ctaui=rfi(4)
        ctauj=rfj(4)
        tcol=ctaui
        if(ctaui.lt.ctauj)tcol=ctauj
c	for back to back collision the collision time is equal to 
c	 the lager one of their formation times
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
c	write(9,*)'back to back,icp,tcol',icp,tcol   ! sa
	kflag=1
	if(tcol.le.pxc(4))return
	if(tcol.le.pyc(4))return
c	for the back to back collisions (collision happens in future)
	else   
c	write(9,*)'face to face,icp',icp ! sa
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
c020603	collision happens in future
	if(tcol-ctaui .le. 0.)return
	if(tcol-ctauj .le. 0.)return
c	collision must be after formation
	do ik=1,3
	dr(ik)=dr(ik)+db(ik)*tcol
	rtai=rtai+dr(ik)*dr(ik)
        enddo
c	for the face to face collisions

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
100	tc(1,icp)=tcol1
	tc(2,icp)=tcol2
	lc(1,icp)=i
	lc(2,icp)=j
c	write(9,*)'i,j,tcol1,tcol2=',i,j,tcol1,tcol2   ! sa
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
        parameter (mplis=80000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
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
c        write(9,*)'peo(5)=',peo(5)   ! w
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
c	calculate the charge (in unit of 1/3) of parton   ! 180520 
c	kf: parton flaver
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
cc	tsmp: t value sampled.
cc	lmn: order number of happening process 
cc	sig: the total cross section of parton kf1 colliding with kf2
cc	kf3 and kf4: kf code of the colliding pair after interaction
cc      eiej2: squared invariant mass of colliding pair
c180520 jjj: jjj-th loop within a event
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
c	leadding order pQCD differential cross section of 2->2 processes
	external fs11_0,fs11_1,fs11_2,fs12_0,fs12_1,fsqq,fsgg_0,fsgg_1,
     c	fsqg

	parameter (msca=20000,mplis=80000)   ! 080520

	common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc   ! 220803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 240803 181003
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 230520
        common/syspar_p/rsig1,pio,tcut
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)   ! 080520
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
c250803	common/work7/reac(9),crose(9)
	dimension ssig(0:3),eee(0:1000),dd(0:1000),eee1(0:1000),
     c	eee2(0:1000),eee3(0:1000)
	logical lgg,lqg,lqq,lqaq,lqqp,lqaqp

c       classify the incident channel
c	write(*,*)'kf1,kf2=',kf1,kf2
	kf0=iabs(kf1)
	ikf2=iabs(kf2)   ! 090603
c090603	lqq=(kf1.eq.kf2).and.(kf2.ne.21)
	lqq=(kf1.eq.kf2).and.(kf0.ne.21)   ! 090603 ! qq (or q(-)q(-))->
c       lqq: two incoming partons are quarks (or antiquarks) with same 
c	 flavor
	lgg=(kf1.eq.21).and.(kf2.eq.21)                    ! gg->
c       lgg: two gluons
c090603	lqg=((kf1.eq.21).and.(kf2.ne.21)).or.
c090603     c      ((kf2.eq.21).and.(kf1.ne.21))           ! qg (or q(-)g)->
	lqg=((kf0.eq.21).and.(kf2.ne.21)).or.
     c      ((ikf2.eq.21).and.(kf1.ne.21))           ! qg (or q(-)g)->
c       lqg: one quark (or antiquark) and one gluon
	lqaq=(kf1.eq.-kf2).and.(kf0.ne.21)          ! qq(-) (or q(-)q)->
c       lqaq: quark (antiquark) and its corresponding antiquark (quark)
c090603	lqqp=(kf0.ne.iabs(kf2)).and.(kf1.ne.21).and.(kf2.ne.21)
	lqqp=(kf0.ne.ikf2).and.(kf0.ne.21).and.(ikf2.ne.21)   ! 090603
c       lqqp: two quarks,or two antiquarks,or one quark and one antiquark
c        with different flavor. 

c	calculate the total cross section for each incedent channel  
	sig=0.
	lmn=0   ! 250803
	do i=0,3
	ssig(i)=0.
	enddo
        idw=adj1(4)   ! 240803 changed from 1000
	adj120=adj1(20)   ! 230405
c250420
        dmass=pymass(1)*2.   ! umass=pymass(2)=dmass
        umass=pymass(2)*2.
        smass=pymass(3)*2.
        cmass=pymass(4)*2.
        bmass=pymass(5)*2.
        tmass=pymass(6)*2.
c250420        
	xs=eiej2
c	xs: squared invariant mass of colliding pair   ! 230520
c080520
        if(eiej2.le.0.)eiej2=1.e-16
        eiej=dsqrt(eiej2)   ! invariant mass of colliding pair
c260620 if(ilo.eq.1 .and. jjj.eq.1)then
c       write(9,*)'in fsig gg-> ilo=',ilo
c       write(9,1006)dmass,umass,smass,cmass,bmass,tmass
c       endif
c200620 1006    format(6(1x,f12.2))
c260620 if(ilo.eq.1)write(9,*)'jjj,eiej2,eiej=',jjj,eiej2,eiej
c080520        
	sm=0.
c	the sum of squared mass of partons in initial and final states, it 
c	 is set to zero (zero mass appoxi. in kinematics treatment)
	nf=3
c112399	as=12.*pio/(33-2*nf)    ! alpha_s
c121099	as=12.*pio/(33-2*nf)/xs    ! alpha_s, 112399
c050603	as=.3    ! alpha_s, 121099, .47 (Bin)
	as=adj1(2)   ! 240803
c280202
c       uxt=0.   ! tmax
c       dxt=-xs    ! tmin
c       cf. my note book p.15
c280803	above uxt and dxt are used in limited and regularized cross section
c280202
c050603	tcut=.4   ! 120999, 5.86
csa's   tcut: the cut for avoidding divergence
csong's	tcut: the cut value of t, a process is considered to be 'hard' 
csong's	 when t>tcut, 'soft' when t<tcut.
	tcut=adj1(3)   ! 240803
c	write(9,*)'as,tcut=',as,tcut   ! sa
c280803
	uxt=-tcut
cc	uxt: the upper limit (max.) of t
c	note,the integral variable t is negative 
c120699	dxt=(sm-xs)+tcut
	dxt=-xs-tcut   ! 291102 120505
c120505	dxt=-xs+tcut  ! 030903
cc	dxt: the lower limit (min.) of t
c	above uxt and dxt are used in LO pQCD cross section
c280803
c	write(9,*)'as,dxt,uxt=',as,dxt,uxt   ! sa
	if(uxt.le.dxt) then
	ilo=-2
c	the process can not happen
	return
	endif
	if(lqaq.or.lgg) then   ! 1

	if(kf0.eq.21) then  ! 2	gg->
	call integ(fsgg_0,uxt,dxt,idw,eee,sum)  ! ->gg
c	fsgg_0: differential cross section gg->gg
c	write(9,*)'9,idw,eee=',idw,eee(i)   ! sa
c	write(9,*)'ux,dx,sum=',uxt,dxt,sum   ! sa
c160902
        sum1=sum
        do i1=0,idw
        eee1(i1)=eee(i1)
        enddo
c230520
        iparreso=iparres
        if(iparres.eq.1 .and. eiej.lt.dmass)then 
        iparres=0   
c       treates as elastic (iparres=0)
        call integ(fsgg_1,uxt,dxt,idw,eee,sum) 
c       add (common/sa33/smadel,ecce,secce,parecc,iparres) today in 'fsig'   
        iparres=iparreso
        else
	call integ(fsgg_1,uxt,dxt,idw,eee,sum)  ! ->qq(-)
        endif
c230520        
	sum2=sum
	do i1=0,idw
	eee2(i1)=eee(i1)
	enddo
	sums=sum1+sum2
	sum3=sum1/sums
	sele=pyr(1)
	if(sele.le.sum3)then   ! 3 gg->gg
        do i1=0,idw
        eee(i1)=eee1(i1)
        enddo
        kf3=kf1
        kf4=kf2
c250803        reac(9)=reac(9)+1.
c250803        crose(9)=crose(9)+sum1
	lmn=9   ! 250803
        else   ! 3 gg->q*q(-) 
	do i1=0,idw
	eee(i1)=eee2(i1)
	enddo
c       ->qq(-) branch
        aa=pyr(1)
c250420
	if(eiej.lt.smass)then ! d,u   
        if(aa.le.0.5)then
        kf3=1                                        ! ->dd(-)
        else
	kf3=2                                        ! ->uu(-)
        endif    
        goto 1001 
        endif

        if(eiej.ge.smass .and. eiej.lt.cmass)then ! d,u,s
        if(aa.le.0.3333)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.3333 .and. aa.le.0.6666)then
        kf3=2                                        ! ->uu(-)
        else
        kf3=3                                        ! ->ss(-)
        endif
        goto 1001
        endif

        if(eiej.ge.cmass .and. eiej.lt.bmass)then ! d,u,s,c
        if(aa.le.0.25)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.25 .and. aa.le.0.5)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.5 .and. aa.le.0.75)then
        kf3=3                                        ! ->ss(-)
        else
        kf3=4                                        ! ->cc(-)     
        endif
        goto 1001
        endif

        if(eiej.ge.bmass .and. eiej.lt.tmass)then ! d,u,s,c,b
        if(aa.le.0.2)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.2 .and. aa.le.0.4)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.4 .and. aa.le.0.6)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.6 .and. aa.le.0.8)then
        kf3=4                                        ! ->cc(-)
        else
        kf3=5                                        ! ->bb(-)   
        endif
        goto 1001
        endif

        if(eiej.ge.tmass)then ! d,u,s,c,b,t
        if(aa.le.0.1666)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.1666 .and. aa.le.0.3333)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.3333 .and. aa.le.0.4998)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.4998 .and. aa.le.0.6664)then
        kf3=4                                        ! ->cc(-)
        elseif(aa.gt.0.6664 .and. aa.le.0.833)then
        kf3=5                                        ! ->bb(-)
        else
        kf3=6                                        ! ->tt(-) 
        endif
        endif
c250420
1001	kf4=-kf3   ! 250420
c250803        reac(7)=reac(7)+1.
c250803        crose(7)=crose(7)+sum2
c180520
	lmn=7   ! 250803
c230520
c260620 if(ilo.eq.1)then
c       write(9,*)'in fsig jjj,lmn kf1-4=',jjj,lmn,kf1,kf2,kf3,kf4 
c260620 endif        
c230520
        endif   ! 3
        sum=sums 
        sig=sums 
c160902
c160902	sig=sum
c160902	kf3=kf1
c160902	kf4=kf2
c160902	reac(9)=reac(9)+1.
c160902	crose(9)=crose(9)+sum

	else       ! 2 qq(-) (or q(-)q)-> 
	call integ(fs11_2,uxt,dxt,idw,eee,sum) ! ->qq(-) (or q(-)q)
c	fs11_2: differential cross section qq(-) (or q(-)q)->qq(-) 
c	 (or q(-)q)
c	write(9,*)'5,idw,eee=',idw,eee(i)   ! sa
c	write(9,*)'ux,dx,sum=',uxt,dxt,sum   ! sa
c160902
	sum1=sum
	do i1=0,idw
	eee1(i1)=eee(i1)
	enddo
	call integ(fsqq,uxt,dxt,idw,eee,sum) ! ->gg
c	write(9,*)'6,idw,eee=',idw,eee(i)   ! sa
c	write(9,*)'ux,dx,sum=',uxt,dxt,sum   ! sa
	sum2=sum
	do i1=0,idw
	eee2(i1)=eee(i1)
	enddo
	call integ(fs11_1,uxt,dxt,idw,eee,sum) ! ->q'q'(-) (or q'(-)q')
	sum3=sum
	do i1=0,idw
	eee3(i1)=eee(i1)
	enddo
	sums=sum1+sum2+sum3
	sum4=sum1/sums
	sum5=(sum1+sum2)/sums
	sele=pyr(1)
	if(sele.le.sum4)then   ! 3 ->qq(-) (or q(-)q)
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
	elseif(sele.gt.sum4 .and. sele.le.sum5)then ! 3 ->gg   
c	sum=sum1
c	sig=sum1
	do i1=0,idw
	eee(i1)=eee2(i1)
	enddo
	kf3=21
	kf4=21
c250803        reac(6)=reac(6)+1.
c250803        crose(6)=crose(6)+sum2
	lmn=6   ! 250803
c090603	elseif(sele.gt.sum5)then
	else   ! 3 ->q'q'(-) (or q'(-)q'), process 4 (inelastic) 
	do i1=0,idw
	eee(i1)=eee3(i1)
	enddo
	aa=pyr(1)

c250420 branch
	if(eiej.lt.smass)then ! d,u

        if(kf0.eq.1)then   ! initial state is dd(-) (or d(-)d)        
        kf3=2                                        ! ->uu(-)
        elseif(kf0.eq.2)then   ! initial state is uu(-) (or u(-)u)
	kf3=1                                        ! ->dd(-)
        else
        endif   
        goto 1002

        endif


        if(eiej.ge.smass .and. eiej.lt.cmass)then ! u,d,s

        if(kf0.eq.1)then   ! initial state is dd(-) (or d(-)d)                
        if(aa.le.0.5)then
        kf3=2                                        ! ->uu(-)
        else
        kf3=3                                        ! ->ss(-)
        endif
        goto 1002
        endif

        if(kf0.eq.2)then   ! initial state is uu(-) (or u(-)u)                
        if(aa.le.0.5)then
        kf3=1                                        ! ->dd(-)
        else
        kf3=3                                        ! ->ss(-)
        endif
        goto 1002
        endif

        if(kf0.eq.3)then   ! initial state is ss(-) (or s(-)s)                
        if(aa.le.0.5)then
        kf3=1                                        ! ->dd(-)
        else
        kf3=2                                        ! ->uu(-)
        endif
        goto 1002
        endif

        endif


        if(eiej.ge.cmass .and. eiej.lt.bmass)then ! d,u,s,c

        if(kf0.eq.1)then   ! initial state is dd(-) (or d(-)d)
        if(aa.le.0.3333)then                                        
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.3333 .and. aa.le.0.6666)then
        kf3=3                                        ! ->ss(-)
        else
        kf3=4                                        ! ->cc(-)     
        endif
        goto 1002
        endif

        if(kf0.eq.2)then   ! initial state is uu(-) (or u(-)u)
        if(aa.le.0.3333)then                                        
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.3333 .and. aa.le.0.6666)then
        kf3=3                                        ! ->ss(-)
        else
        kf3=4                                        ! ->cc(-)     
        endif
        goto 1002   
        endif

        if(kf0.eq.3)then   ! initial state is ss(-) (or s(-)s)
        if(aa.le.0.3333)then                                        
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.3333 .and. aa.le.0.6666)then
        kf3=2                                        ! ->uu(-)
        else
        kf3=4                                        ! ->cc(-)     
        endif
        goto 1002
        endif

        if(kf0.eq.4)then   ! initial state is cc(-) (or c(-)c)
        if(aa.le.0.3333)then                                        
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.3333 .and. aa.le.0.6666)then
        kf3=2                                        ! ->uu(-)
        else
        kf3=3                                        ! ->ss(-)     
        endif
        goto 1002
        endif

        endif


        if(eiej.ge.bmass .and. eiej.lt.tmass)then ! d,u,s,c,b

        if(kf0.eq.1)then   ! initial state is dd(-) (or d(-)d)                
        if(aa.le.0.25)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.25 .and. aa.le.0.5)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.5 .and. aa.le.0.75)then
        kf3=4                                        ! ->cc(-)
        else
        kf3=5                                        ! ->bb(-)   
        endif
        goto 1002
        endif

        if(kf0.eq.2)then   ! initial state is uu(-) (or u(-)u)                
        if(aa.le.0.25)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.25 .and. aa.le.0.5)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.5 .and. aa.le.0.75)then
        kf3=4                                        ! ->cc(-)
        else
        kf3=5                                        ! ->bb(-)   
        endif
        goto 1002
        endif

        if(kf0.eq.3)then   ! initial state is ss(-) (or s(-)s)                
        if(aa.le.0.25)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.25 .and. aa.le.0.5)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.5 .and. aa.le.0.75)then
        kf3=4                                        ! ->cc(-)
        else
        kf3=5                                        ! ->bb(-)   
        endif
        goto 1002
        endif

        if(kf0.eq.4)then   ! initial state is cc(-) (or c(-)c)                
        if(aa.le.0.25)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.25 .and. aa.le.0.5)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.5 .and. aa.le.0.75)then
        kf3=3                                        ! ->ss(-)
        else
        kf3=5                                        ! ->bb(-)   
        endif
        goto 1002
        endif

        if(kf0.eq.5)then   ! initial state is bb(-) (or b(-)b)                
        if(aa.le.0.25)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.25 .and. aa.le.0.5)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.5 .and. aa.le.0.75)then
        kf3=3                                        ! ->ss(-)
        else
        kf3=4                                        ! ->cc(-)   
        endif
        goto 1002
        endif

        endif


        if(eiej.ge.tmass)then ! d,u,s,c,b,t

        if(kf0.eq.1)then   ! initial state is dd(-) (or d(-)d)                
        if(aa.le.0.2)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.2 .and. aa.le.0.4)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.4 .and. aa.le.0.6)then
        kf3=4                                        ! ->cc(-)
        elseif(aa.gt.0.6 .and. aa.le.0.8)then
        kf3=5                                        ! ->bb(-)
        else
        kf3=6                                        ! ->tt(-) 
        endif
        goto 1002
        endif

        if(kf0.eq.2)then   ! initial state is uu(-) (or u(-)u)                
        if(aa.le.0.2)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.2 .and. aa.le.0.4)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.4 .and. aa.le.0.6)then
        kf3=4                                        ! ->cc(-)
        elseif(aa.gt.0.6 .and. aa.le.0.8)then
        kf3=5                                        ! ->bb(-)
        else
        kf3=6                                        ! ->tt(-) 
        endif
        goto 1002
        endif

        if(kf0.eq.3)then   ! initial state is ss(-) (or s(-)s)                
        if(aa.le.0.2)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.2 .and. aa.le.0.4)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.4 .and. aa.le.0.6)then
        kf3=4                                        ! ->cc(-)
        elseif(aa.gt.0.6 .and. aa.le.0.8)then
        kf3=5                                        ! ->bb(-)
        else
        kf3=6                                        ! ->tt(-) 
        endif
        goto 1002
        endif

        if(kf0.eq.4)then   ! initial state is cc(-) (or c(-)c)                
        if(aa.le.0.2)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.2 .and. aa.le.0.4)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.4 .and. aa.le.0.6)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.6 .and. aa.le.0.8)then
        kf3=5                                        ! ->bb(-)
        else
        kf3=6                                        ! ->tt(-) 
        endif
        goto 1002
        endif

        if(kf0.eq.5)then   ! initial state is bb(-) (or b(-)b)                
        if(aa.le.0.2)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.2 .and. aa.le.0.4)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.4 .and. aa.le.0.6)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.6 .and. aa.le.0.8)then
        kf3=4                                        ! ->cc(-)
        else
        kf3=6                                        ! ->tt(-) 
        endif
        goto 1002
        endif

        if(kf0.eq.6)then   ! initial state is tt(-) (or t(-)t)                
        if(aa.le.0.2)then
        kf3=1                                        ! ->dd(-)
        elseif(aa.gt.0.2 .and. aa.le.0.4)then
        kf3=2                                        ! ->uu(-)
        elseif(aa.gt.0.4 .and. aa.le.0.6)then
        kf3=3                                        ! ->ss(-)
        elseif(aa.gt.0.6 .and. aa.le.0.8)then
        kf3=4                                        ! ->cc(-)
        else
        kf3=5                                        ! ->bb(-) 
        endif
        endif

        endif


c250420
1002	kf4=-kf3   ! 250420
c250803        reac(4)=reac(4)+1.
c250803        crose(4)=crose(4)+sum3
	lmn=4   ! 250803
c090603	else
	endif   ! 3
c160902
	sum=sums   ! sum1 050605 Tan
	sig=sums   ! sum1 050605 Tan
	endif   ! 2

c090603	endif   ! 1

	elseif(lqq)then   ! 1 090603
	call integ(fs11_0,uxt,dxt,idw,eee,sum)
c	fs11_0 : differential x section of qq->qq (or q(-)q(-)->q(-)q(-))
c	write(9,*)'2,idw,eee=',idw,eee(i)   ! sa
c	write(9,*)'ux,dx,sum=',uxt,dxt,sum   ! sa
	sig=sum
	kf3=kf1
	kf4=kf2
c250803	reac(2)=reac(2)+1.
c250803	crose(2)=crose(2)+sum
	lmn=2   ! 250803
c090603	endif

	elseif(lqqp)then   ! 1 090603
	call integ(fs12_0,uxt,dxt,idw,eee,sum)
c	fs_12 : differential x section of qq'->qq',differential x section 
c	 of qq(-)'->qq(-)' (or q(-)q'->q(-)q')=differential x section of 
c	 qq'->qq',differential x section of q(-)q(-)'->q(-)q(-)'=
c	 differential x section of qq'->qq'
	if(kf1*kf2.gt.0)then
c	write(9,*)'1,idw,eee=',idw,eee(i)   ! sa
c	write(9,*)'ux,dx,sum=',uxt,dxt,sum   ! sa
	sig=sum
	kf3=kf1
	kf4=kf2
c250803	reac(1)=reac(1)+1.
c250803	crose(1)=crose(1)+sum
	lmn=1   ! 250803
	endif
        if(kf1*kf2.lt.0)then
c	write(9,*)'3,idw,eee=',idw,eee(i)   ! sa
c	write(9,*)'ux,dx,sum=',uxt,dxt,sum   ! sa
	sig=sum
	kf3=kf1
	kf4=kf2
c250803	reac(3)=reac(3)+1.
c250803	crose(3)=crose(3)+sum
	lmn=3   ! 250803
	endif
c090603	endif

	elseif(lqg)then   ! 1 090603 
	call integ(fsqg,uxt,dxt,idw,eee,sum)
c	fsqg : differential x section of qg->qg (or q(-)g->q(-)g)
c	write(9,*)'8,idw,eee=',idw,eee(i)   ! sa
	sig=sum
c220803
c	if(iiii.eq.2 .and. sum.lt.0.)then
c	write(9,*)'idw,ux,dx,sum=',idw,uxt,dxt,sum ! sa
c	write(9,*)'as,tcut,sm,xs=',as,tcut,sm,xs ! sa
c	write(9,900)(eee(i),i=0,idw) 
c	endif
900	format(7(1x,f10.4)/)
c220803
	kf3=kf1
	kf4=kf2
c250803	reac(8)=reac(8)+1.
c250803	crose(8)=crose(8)+sum
	lmn=8   ! 250803
	else   ! 1 090603

	endif   ! 1

	sig=sig*0.04   ! 121099
c121099	0.04 is the transformation factor from (GeV)^-2 to (fm)^2
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
c	tsmp=-pyr(1)*xs
	seta=pyr(1)*pio
	tsmp=xs*(dcos(seta)-1)/2+tcut   ! 280407
	return
	endif   
c230405
	call eee_dd(idw,eee,dd)
	call samp_integ(idw,uxt,dxt,tt,dd)
	tsmp=tt+tcut   ! 200407
c	tsmp: the t value sampled.
	endif
	return
	end

	

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine collis(ik1,ik2,kf3,kf4,tcp,jjj,kkk,iway,
     c   icnew,jcnew,lmn,time)   ! 120505 160110 
c	performs parton-parton collision & updates particle list
c       if lmn=4,6,& 7, updates 'parlist','pyjets',
c        'sbe',diquark list,string list, & lc(1-2,m) either
c	ik1,ik2: line number of the colliding pair in parton list
c180520 kf3 and kf4: kf code of the colliding pair after interaction 
c	tcp: collision time
c	jjj: jjj-th loop within a event
c	if kkk=1 throw away current collision
	parameter (mplis=80000,msca=20000)
        parameter (kszj=80000)   ! 160110
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 160110
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
	common/syspar_p/rsig1,pio,tcut
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)   ! 160110
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5) ! 280620
	common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp   ! 180705
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   ! 090820
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 160110
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620   
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 050620
	common/work7/reac(9),crose(9)    
	common/ctllist_p/nreac(9),nrel
	common/show/vip(mplis),xap(mplis)
c160110 ifcom(i): line number (in 'sbe') of first component of i-th diquark
c	nreac(i): statistics of the # of successful i-th collision
c	nrel: statistics of the # of collision blocked
cc	nsca: number of particles after collision  
cc	pip(1-3,*)	: momentun of particle after collision
cc	pip(4,*)	: energy of particle after collision
cc	pip(5,*)	: virtuality of particle after collision
cc	pip(6,*)	: x value of particle after collision
cc	kpip(*)		: flavor code of particle after collision
cc      vip: virtuality of parton
cc      xap: momentum fraction of parton  
c180705	ppsa(5): cumulate the charge losed in collision processes in an event
c	definition of ppsa(1), ..., and ppsa(4) are given in 'eloss'
 
	dimension pi(4),pj(4),pii(4),pjj(4),pi_o(4),pj_o(4),ps(4)
	dimension pss(mplis,5),px(4),py(4),pij(4)   ! 230520
	double precision b(3)
	adj12=adj1(12)
c	write(9,*)'in collis adj12=',adj12   !sa
	kf1=idp(ik1)
	kf2=idp(ik2)
	stili=0.
c	stili: statistics of the times calling 'ti_li1'
c180520 kf1 and kf2: kf code of the colliding pair before interaction
c	write(9,*)'colli iiii,jjj,iprl=',iiii,jjj,iprl ! sa
c	write(9,*)'ik1,ik2,kf1,kf2=',ik1,ik2,kf1,kf2   ! sa
c	do i=1,iprl   ! sa
c	write(9,509)(pp(m,i),m=1,4),rmp(i)  ! sa
c	enddo   ! sa
509     format(5(1x,f10.4))

c       statistics of the total charge of partons
        ichabe=0.
        do i=1,iprl
        ik=idp(i)
        ichabe=ichabe+ichge(ik)
        enddo

	nsca=0
	pi(4)=pp(4,ik1)
	pj(4)=pp(4,ik2)
        pi_o(4)=pi(4)   
        pj_o(4)=pj(4)   
	pij4=pi(4)+pj(4)
	do i=1,3
	pi(i)=pp(i,ik1)
	pj(i)=pp(i,ik2)
        pi_o(i)=pi(i)   
        pj_o(i)=pj(i)   
	b(i)=(pi(i)+pj(i))/pij4
	enddo
c	if(jjj.eq.1)then
c	write(9,*)'before lorentz, b=',(b(i),i=1,3)   ! sa
c	write(9,510)(pi(i),i=1,4)   ! sa
c	write(9,510)(pj(i),i=1,4)   ! sa	
c	write(9,510)((pi(i)+pj(i)),i=1,4)   ! sa
c	endif
510     format(4(1x,f10.4))
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
c	if(jjj.eq.1)then 
c	write(9,*)'after lorentz'   ! sa
c	write(9,510)(pi(i),i=1,4)   ! sa
c	write(9,510)(pj(i),i=1,4)   ! sa
c	write(9,510)((pi(i)+pj(i)),i=1,4)   ! sa
c	write(9,*)'pij4,eiej2=',pij4,eiej2   ! sa
c	endif
c	call lorntz(1,b,pi,pj)
c	write(9,*)'boost back immediately, b=',(b(i),i=1,3)   ! sa
c	write(9,510)(pi(i),i=1,4)   ! sa
c	write(9,510)(pj(i),i=1,4)   ! sa
c120505
	if(eiej2.lt.0.)then
	iway=1
c	iway=1 means collision does not happen
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
c180520
c260620 write(9,*)'af. call fsig kf12,lmu=',kf1,kf2,lmn
c       write(9,1005)(nreac(i1),i1=1,9)
c260620 1005    format(9(1x,i5))   
c180520        
        else
        sig=0.
        endif
c160110   
c250803	lmn: order number of the process happened
c120603
	if(sig.le.0.)then   
	jjj=jjj-1   
	kkk=1   ! throw away that collision
	return
	endif
c120603
c	am1=amass(kf1)
c	am2=amass(kf2)
c	am3=amass(kf3)
c	am4=amass(kf4)
	am1=0.   
	am2=0.   
	am3=0.   
	am4=0.   
c	in consistent with sm in "fsig"
c230520	kinematics of 2->2 process,no matter what is the final state, is
c        treated in zero mass approximation from momentum point of view 
	paa=dsqrt(pi(1)**2+pi(2)**2+pi(3)**2)
c	write(9,*)'paa,xs,xt,xu=',paa,xs,xt,xu   ! sa
	xs=eiej2
	xt=tsmp
c	xu=am1**2+am2**2+am3**2+am4**2-xs-xt
	xu=-xs-xt   ! zero rest mass approximation
c	square of momentum transfer (Q**2)
	qc2=xt*xu/xs   ! pt**2, as the same in PYTHIA (Q^2_{hard})
c151018	qc2=4.*qc2   ! as the same in PYTHIA
c	qc2=min1(-xt,-xu)

c	calculate the direction cosines of momentum of one colliding
c        particle after scattering relative to the momentum
c        of corresponding particle before scattering
	fi=2.*pio*pyr(1)
c	relation between t and cos(ceta) in zero rest mass approximation
c	 (i. e. treated as ela. scattering)
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
c	 pair (before scattering and rotation), as output are four 
c	 momentum of that pair after scattering and rotation
c       pc1: momentum modulus of pi or pj, both are equal in their cms 
c	 before and after scattering
c191103	pi(4)=dsqrt(pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)+am3*am3)
c191103	pj(4)=dsqrt(pj(1)*pj(1)+pj(2)*pj(2)+pj(3)*pj(3)+am4*am4)
c	if(jjj.eq.1)then
c	write(9,*)'after retation'   ! sa
c	write(9,510)(pi(i),i=1,4)   ! sa
c	write(9,510)(pj(i),i=1,4)   ! sa
c	write(9,510)((pi(i)+pj(i)),i=1,4)   ! sa
c	endif

        do i=1,4
        pii(i)=pi(i)   ! initial four momentum of parton-parton scattering
        pjj(i)=pj(i)
        enddo
c	boost back
	call lorntz(1,b,pi,pj)
c	if(jjj.eq.1)then
c	write(9,*)'after boost back'   ! sa
c	write(9,510)(pi(i),i=1,4)   ! sa
c	write(9,510)(pj(i),i=1,4)   ! sa
c	write(9,510)((pi(i)+pj(i)),i=1,4)   ! sa
c	endif

c080520 proceeds for w/o time-like branching 

c       final state of two parton scattering (nsca=2)        
	nsca=2
	kpip(1)=kf3
	kpip(2)=kf4
	do i=1,4
	pip(i,1)=pi(i)
	pip(i,2)=pj(i)
	enddo
c	write(9,*)'pip='   ! sa
c	write(9,510)(pip(i,1),i=1,4)   ! sa
c	write(9,510)(pip(i,2),i=1,4)   ! sa

c       charge conservation
        ichbe=ichge(kf1)+ichge(kf2)
        ichaf=ichge(kf3)+ichge(kf4)
        if(ichbe.ne.ichaf)then
        write(9,*)'w/o time-like iiii,jjj,kf1,kf2=',
     c   iiii,jjj,kf1,kf2  ! sa 160110
        write(9,*)'nsca,kf3,kf4,ichbe,ichaf=',nsca,kf3,kf4,ichbe,ichaf! sa
        ppsa(5)=ppsa(5)-ichbe+ichaf   ! 180705 270407
        endif
c260620 write(9,*)'be. update particle list ik12,kf12=',ik1,ik2,kf1,kf2

c230520	updates particle list  (note: line numbers of the colliding 
c        pair after interaction are the same as ones before interaction)

c230520 updates 'parlist', i.e. feedback final state of two parton scattering
c        to ik1 & ik2
	l=ik1
	if(pip(4,1).lt.1.e-10)pip(4,1)=1.e-10   ! 031204
	do k1=1,3
	pp(k1,l)=pip(k1,1)
	vp(k1,l)=pip(k1,1)/pip(4,1)
	enddo
	pp(4,l)=pip(4,1)
	idp(l)=kpip(1)
	taup(l)=0.
	rmp(l)=amass(kpip(1))
c230520 updates 'pyjets'
        k(l,2)=idp(l)
        do j1=1,4
        v(l,j1)=rp(j1,l)
        p(l,j1)=pp(j1,l)
        enddo
        p(l,5)=rmp(l)
        v(l,5)=0.

	l=ik2
	if(pip(4,2).lt.1.e-10)pip(4,2)=1.e-10   ! 031204
	do k1=1,3
	pp(k1,l)=pip(k1,2)
	vp(k1,l)=pip(k1,2)/pip(4,2)
	enddo
	pp(4,l)=pip(4,2)
	idp(l)=kpip(2)
	taup(l)=0.
	rmp(l)=amass(kpip(2))
c230520 updates 'pyjets'
        k(l,2)=idp(l)
        do j1=1,4
        v(l,j1)=rp(j1,l)
        p(l,j1)=pp(j1,l)
        enddo
        p(l,5)=rmp(l)
        v(l,5)=0.
c	write(9,*)'in collis, iiii,jjj,lmn=',iiii,jjj,lmn   ! sa
	reac(lmn)=reac(lmn)+1.
	crose(lmn)=crose(lmn)+sig
	nreac(lmn)=nreac(lmn)+1   ! 071103
c180520
c260620 write(9,*)'af. update particle list jjj,lmn=',jjj,lmn
c       write(9,1005)(nreac(i1),i1=1,9)  
c       if(lmn.eq.7)then
c       write(9,*)'kf1-4=',kf1,kf2,kf3,kf4
c       write(9,*)'ik12,kpip12=',ik1,ik2,kpip(1),kpip(2)
c       write(9,*)'k(1-2)=',k(ik1,2),k(ik2,2)
c260620 endif   
c180520
	if(iparres.eq.0)return   ! 060520
        if(iparres.eq.1.and.(lmn.ne.4.and.lmn.ne.6.and.lmn.ne.7))return!080520
c       if(adj12.eq.1)return   ! 220820

c       proceeds for lmn=4 or =6 or =7 and hadronization by string 
c        fragmentation

        if(lmn.eq.4.or.lmn.eq.6.or.lmn.eq.7)then   !! 1

c070720 update collision list
        call update_ctlm(time,iway)   ! 070720

        if(mstj1_2.eq.1)goto 1004 ! w/ time-like branching
1008    continue
c240820 copies new string composed of final state parton pair (ik1 & ik2) 
c        to the end of 'parlist' ('pyjets') 
        if(ik2.gt.ik1)then   ! 230520
        k(ik1,1)=2
        k(ik2,1)=1
c       copy ik1 to iprl+1 (n+1)
        call coend(ik1)
        icnew=iprl
c       copy ik2 to iprl+1 (n+1)
        call coend(ik2)
        jcnew=iprl
        endif
        if(ik2.lt.ik1)then   ! 230520
        k(ik1,1)=1
        k(ik2,1)=2
c       copy ik2 to iprl+1 (n+1)
        call coend(ik2)
        jcnew=iprl
c       copy ik1 to iprl+1 (n+1)
        call coend(ik1)
        icnew=iprl
        endif
        nstr1=nstr1+1
        nstr1a(nstr1)=iprl-1
        nstr1v(nstr1)=iprl
c260620 write(9,*)'af. new string n-1,n,kf=',n-1,n,k(n-1,2),k(n,2)
c260620 write(9,*)'ijcnew,idp=',icnew,jcnew,idp(icnew),idp(jcnew)
c080520 endif   

c       treats scattering parton pair (scattered parton pair has been
c        treated as a new string above)
c        write(9,*)'iiii,iii,jjj,lmn,nbe,ik12,kf1,kf2,kf3,kf4=',
c     c   iiii,iii,jjj,lmn,nbe,ik1,ik2,kf1,kf2,kf3,kf4

        if(lmn.eq.7)then   !! 2 gg->qqbar
c       throws away scattering parton pair (ik1 & ik2) 
        if(ik1.gt.ik2)then
        call parmov(ik1,ik1,ik2,lmn)    
        call parmov(ik2,ik1,ik2,lmn)
        endif
        if(ik1.lt.ik2)then
        call parmov(ik2,ik1,ik2,lmn)
        call parmov(ik1,ik1,ik2,lmn)
        endif
        return
        endif   !! 2

        if(lmn.eq.4.or.lmn.eq.6)then   !! 3, 4:q1q1bar->q2q2bar, 6:qqbar->gg          
        if(adj12.eq.1)then   !! 4 fragments with coalescence
c       throws away scattering parton pair (ik1 & ik2) 
        if(ik1.gt.ik2)then
        call parmov(ik1,ik1,ik2,lmn)
        call parmov(ik2,ik1,ik2,lmn)
        endif
        if(ik1.lt.ik2)then
        call parmov(ik2,ik1,ik2,lmn)
        call parmov(ik1,ik1,ik2,lmn)
        endif
        return
        elseif(adj12.eq.0)then   !! 4        
c       does ik1 (ik2) is a component of string
        call adjst(ik1,ik1str,ik1sa,ik1sv)   ! 020620
c       ik1str: oredr number of string to which ik1 belongs,equal 0 otherwise
c       ik1sa: line number of first component of above string,equal 0 otherwise 
c       ik1sv: line number of last component of above string,equal 0 otherwise
        call adjst(ik2,ik2str,ik2sa,ik2sv)   ! 020620
c        write(9,*)'ik1str,ik1sa,ik1sv,ik2str,ik2sa,ik2sv=',
c     c   ik1str,ik1sa,ik1sv,ik2str,ik2sa,ik2sv

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
c        write(9,*)'idi1,iway1,ifcom1,npt1=',idi1,iway1,ifcom1,npt1
c        write(9,*)'idi2,iway2,ifcom2,npt2=',idi2,iway2,ifcom2,npt2

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

c       moves ik1 & ik2 out as well as updates lists                 
        if(ik1.gt.ik2)then
        call parmov(ik1,ik1,ik2,lmn)   ! 070720
        call parmov(ik2,ik1,ik2,lmn)   ! 070720
        endif
        if(ik1.lt.ik2)then
        call parmov(ik2,ik1,ik2,lmn)   ! 070720
        call parmov(ik1,ik1,ik2,lmn)   ! 070720
        endif
        return

        endif

        if(ik1str.eq.0 .and. ik2str.ne.0)then

        if(ik1.lt.ik2)then
        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)   ! 070720
        call parmov(ik1,ik1,ik2,lmn)   ! 070720
        endif
        if(ik1.gt.ik2)then    
        call parmov(ik1,ik1,ik2,lmn)   ! 070720
        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)   ! 070720
        endif    
        return

        endif

        if(ik1str.ne.0 .and. ik2str.eq.0)then

        if(ik1.lt.ik2)then
        call parmov(ik2,ik1,ik2,lmn)   ! 070720
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        else        
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        call parmov(ik2,ik1,ik2,lmn)   ! 070720
        endif    
        return

        endif

c       proceeds for both of ik1str & ik2str are not equal to zero        
        if(ik1str.eq.ik2str)then
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        endif
        if(ik1str.gt.ik2str)then
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)   ! 070720
        endif        
        if(ik1str.lt.ik2str)then
        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)   ! 070720
        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)   ! 070720
        endif    
        return  
        else   !! 4
        endif  !! 4  


        endif   !! 3
        
c280620   

1004    continue        
c0500   initial state of time-like branching        
        nsca=2
        kpip(1)=kf3
        kpip(2)=kf4
        do i=1,4
        pip(i,1)=pi(i)   
        pip(i,2)=pj(i)
        enddo
c       write(17,*)'pip='   ! sa
c       write(17,510)(pip(i,1),i=1,4)   ! sa
c       write(17,510)(pip(i,2),i=1,4)   ! sa
        pip(5,1)=qc2
        pip(5,2)=qc2
        pip(6,1)=1.  
        pip(6,2)=1.  
c0500	finished-------------------------------

c0700   performs time_like branching-------------------------------
	ichtb=0
	do i1=1,nsca
	kff=kpip(i1)
	ichtb=ichtb+ichge(kff)
	enddo

	nsca0=0
300	continue
	kmsca=nsca
c        write(9,*)'be. ti_li1 iiii,nsca,kmsca=',iiii,nsca,kmsca   ! bh
c       write(9,*)'pip(5,i)=',(pip(5,i),i=nsca0+1,kmsca)   ! sa
	do i=nsca0+1,kmsca   ! do loop
c	write(9,*)'time_like nsca,i=',nsca,i   ! sa
	if(pip(5,i).gt.4*tl0)then   ! if loop
c	phase space for branching vanishes at a virtuality of 4*tl0
c	 cf. Eq.25 in B. R. Webber, Ann. Rev. Nucl. Part. Sci., 36(86)253  
	ea=pip(4,i)
	qa=pip(5,i)
	do j=1,3
	pa(j)=pip(j,i)
	enddo
	kfk=kpip(i)   !!
	xa=pip(6,i)
c	write(*,*)'nsca,nsca0,kmsca,i,kf=',nsca,nsca0,kmsca,i,kfk   !!
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
c        write(9,*)'af. ti_li1 iiii,kmsca,nsca=',iiii,kmsca,nsca   ! bh
c        write(9,*)'ik1,ik2,iprl=',ik1,ik2,iprl   ! bh

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
c0800   'nsca' to 'parlist' ('pyjets')
c       two-body scattering: kf1+kf2->kf3+kf4

        if(nsca.eq.2)then
c        write(9,*)'nsca=2 ik1,ik2,iprl,n,nsca=',ik1,ik2,iprl,n,nsca   ! bh
c       updates 'parlist', i.e. feedback final state of two parton scattering
c        to ik1 & ik2
	l=ik1
	if(pip(4,1).lt.1.e-10)pip(4,1)=1.e-10   ! 031204
	do k1=1,3
	pp(k1,l)=pip(k1,1)
	vp(k1,l)=pip(k1,1)/pip(4,1)
	enddo
	pp(4,l)=pip(4,1)
	idp(l)=kpip(1)
	taup(l)=0.
	rmp(l)=amass(kpip(1))
c230520 updates 'pyjets'
        k(l,2)=idp(l)
        do j1=1,4
        v(l,j1)=rp(j1,l)
        p(l,j1)=pp(j1,l)
        enddo
        p(l,5)=rmp(l)
        v(l,5)=0.

	l=ik2
	if(pip(4,2).lt.1.e-10)pip(4,2)=1.e-10   ! 031204
	do k1=1,3
	pp(k1,l)=pip(k1,2)
	vp(k1,l)=pip(k1,2)/pip(4,2)
	enddo
	pp(4,l)=pip(4,2)
	idp(l)=kpip(2)
	taup(l)=0.
	rmp(l)=amass(kpip(2))
c230520 updates 'pyjets'
        k(l,2)=idp(l)
        do j1=1,4
        v(l,j1)=rp(j1,l)
        p(l,j1)=pp(j1,l)
        enddo
        p(l,5)=rmp(l)
        v(l,5)=0.
c        write(9,*)'af. nsca=2,ik1,ik2,iprl,n,nsca=',ik1,ik2,iprl,n,nsca   ! bh
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
        pip(5,i1)=pip6        
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
        pip(5,i1)=pip6        
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
        pip(5,i3)=pip6
        endif
        enddo

        endif   ! 2

        do i=1,3
        px(i)=rp(i,ik1)   ! three position of ik1
        py(i)=rp(i,ik2)   ! three position of ik2
        enddo
c	induced parton is distributed randumly in between colliding pair
        do i=3,nsca
	rl=pyr(1)   
	do k1=1,3   
	rp(k1,i)=px(k1)*rl+py(k1)*(1.-rl)   
	enddo
	rp(4,i)=tcp
        enddo   

c--------- keeps four momentum conservation-------
	do i2=1,4
	ps(i2)=pi_o(i2)+pj_o(i2)   
c       pi_o (pj_o) four momentum of colliding pair (ik1 & ik2)
	enddo 
c        write(9,*)'ps=',(ps(i2),i2=1,4)   ! sa
	do i1=1,nsca
        do i2=1,4
        pss(i1,i2)=pip(i2,i1)
        enddo
        pss(i1,5)=rmp(i1)
c       write(9,509)(pss(i1,m),m=1,4),rmp(l)  ! sa
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
c	write(9,*)'l=',l   ! sa
c	write(9,509)(pp(m,l),m=1,4),rmp(l)  ! sa
	enddo
c---------- keep four momentum conservation finished -------------

c       a part (first & second) of 'nsca' to 'parlist' ('pyjets'), i.e. 
c        feedback first & second partons in 'nsca' to ik1 & ik2
c       updates 'parlist'
	l=ik1
	if(pip(4,1).lt.1.e-10)pip(4,1)=1.e-10   ! 031204
	do k1=1,3
	pp(k1,l)=pip(k1,1)
	vp(k1,l)=pip(k1,1)/pip(4,1)
	enddo
	pp(4,l)=pip(4,1)
	idp(l)=kpip(1)
	taup(l)=0.
	rmp(l)=amass(kpip(1))
c       updates 'pyjets'
        k(l,2)=idp(l)
        do j1=1,4
        v(l,j1)=rp(j1,l)
        p(l,j1)=pp(j1,l)
        enddo
        p(l,5)=rmp(l)
        v(l,5)=0.
c       updates 'parlist'
	l=ik2
	if(pip(4,2).lt.1.e-10)pip(4,2)=1.e-10   ! 031204
	do k1=1,3
	pp(k1,l)=pip(k1,2)
	vp(k1,l)=pip(k1,2)/pip(4,2)
	enddo
	pp(4,l)=pip(4,2)
	idp(l)=kpip(2)
	taup(l)=0.
	rmp(l)=amass(kpip(2))
c       updates 'pyjets'
        k(l,2)=idp(l)
        do j1=1,4
        v(l,j1)=rp(j1,l)
        p(l,j1)=pp(j1,l)
        enddo
        p(l,5)=rmp(l)
        v(l,5)=0.
c0800	finished----------------------------------------------

c220820 deals with third & fourth partons in 'nsca'
        if(nsca.eq.3)then
c        write(9,*)'nsca=3 ik1,ik2,iprl,n,nsca=',ik1,ik2,iprl,n,nsca   ! bh
c       moves third induced partons to 'trs'
        ntrs=ntrs+1
        ktrs(ntrs,1)=3
        ktrs(ntrs,2)=kpip(3)
        do j=1,4
        ptrs(ntrs,j)=pip(j,3)
        vtrs(ntrs,j)=0.
        enddo
        ptrs(ntrs,5)=0.
        vtrs(ntrs,5)=0.
c        write(9,*)'af. nsca=3 ik1,ik2,iprl,n,nsca=',ik1,ik2,iprl,n,nsca   ! bh
        goto 1008
        endif

        if(nsca.eq.4)then
c        write(9,*)'nsca=4 ik1,ik2,iprl,n,nsca=',ik1,ik2,iprl,n,nsca   ! bh     
c       composes third & fourth partons in 'nsca' as a string & moves it 
c        to the end of 'parlist' ('pyjets')
c       updates 'parlist'
        iprl=iprl+1
        l=iprl
        if(pip(4,3).lt.1.e-10)pip(4,3)=1.e-10   ! 031204
        do k1=1,3
        pp(k1,l)=pip(k1,3)
        vp(k1,l)=pip(k1,3)/pip(4,3)
        enddo
        pp(4,l)=pip(4,3)
        idp(l)=kpip(3)
        taup(l)=0.
        rmp(l)=amass(kpip(3))
c       updates 'pyjets'
        n=n+1
        k(n,1)=2   ! A
        k(n,2)=idp(l)
        do j1=1,4
        v(n,j1)=rp(j1,l)
        p(n,j1)=pp(j1,l)
        enddo
        p(n,5)=rmp(l)
        v(n,5)=0.
c       updates 'parlist'
        iprl=iprl+1
        l=iprl
        if(pip(4,4).lt.1.e-10)pip(4,4)=1.e-10   ! 031204
        do k1=1,3
        pp(k1,l)=pip(k1,4)
        vp(k1,l)=pip(k1,4)/pip(4,4)
        enddo
        pp(4,l)=pip(4,4)
        idp(l)=kpip(4)
        taup(l)=0.
        rmp(l)=amass(kpip(4))
c       updates 'pyjets'
        n=n+1
        k(n,1)=1   ! V
        k(n,2)=idp(l)
        do j1=1,4
        v(n,j1)=rp(j1,l)
        p(n,j1)=pp(j1,l)
        enddo
        p(n,5)=rmp(l)
        v(n,5)=0.
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
c	stili: statistics the times calling 'ti_li1'
c       ij: order # of spliting parton in parton list
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP

	parameter (pio=3.1416)
	parameter (msca=20000,mplis=80000)   ! 080520
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)   ! 080520
	dimension pa1(4),pa2(4),p00(4)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
	common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
c250420
        dmass=pymass(1)   ! umass=pymass(2)=dmass
        umass=pymass(2)
        smass=pymass(3)
        cmass=pymass(4)
        bmass=pymass(5)
        tmass=pymass(6)
c250420        
	il=1
cc	il=1 means the branching will go on, il=0 means 'stop'.
	do i=1,4
	pa1(i)=0.
	pa2(i)=0.
	p00(i)=0.
	enddo
	sampq=0.   ! 211004
10	continue
	sampq=sampq+1.   ! 211004
c	write(9,*)'in ti_li1 sampq=',sampq   ! sa

c	varieties in time like branchs are all the squared mass (m**2) and  
c	 m**2=E**2-|p|**2>0
c------------------------------------------------------------------
c	branching:   A  -> A1   +   A2
c	Vituality:  qqb    qqa     qqc
c	ID No.   :  kfk    kf1     kf2   !!
c	 x       :  xa x'  xa1 x   xa2
c	momentum :  pa     pa1     pa2
c
c	equation : qqb=qqa/z+qqc/(1-z)+pt2/(z*(1-z))  (*)
c-------------------------------------------------------------------
c	above Eq. is deduced from relation of mass and energy (Einstein 
c	 relation) and four momentum conservation, cf. note made in 
c	 Webber's paper

c-------flavor of A1,A2 will be decided in the ensueing module.
	q2max=qa
	kf0=kfk   !!
	call suda_ti(ilo,kf0,ii,q2max,q2,z)
c       time_like Sudakov factor, input: kf0 & q2max, output: others
c	ilo: =0 fail in sampling q2 and z
c	     =1 success in sampling q2 and z							 
	if(ilo.eq.0) then
	qa=tl0   ! forced 
	return
	endif

	qqb=q2
	if(stili.gt.2.)qqb=qa
	zab=z
c	'zab' is the sampled ratio of x/x' in this branching.
	if(zab.eq.0.)return   ! 031204
	xa1=xa*zab
	xa2=xa-xa1
c	write(9,*)xa,zab,xa1,xa2

	if(kfk.ne.21) then   !! initial state is quark
c	branch q->qg
	kf1=kfk   
	kf2=21
	if(pyr(1).ge.0.5) then
c	branch q->gq
	kf1=21
	kf2=kfk   
	endif
	else   !! initial state is gluon	
	if(ii.eq.2) then   !!! 
c	branch: g->gg
	kf1=21
	kf2=21
	else   !!! ii=3
c220820 considering g->gg only
        goto 200   ! 220820
c	branch: g->qq(-)
	aa=pyr(1)
c250420
c       ea=pip(4,ij): energy of spliting particle ij
c       pip(1-3,ij): three momentum of spliting particle ij
        eaa=ea*ea-pip(1,ij)*pip(1,ij)-pip(2,ij)*pip(2,ij)-
     c   pip(3,ij)*pip(3,ij)   ! 080520
        eaa=sqrt(eaa)   ! invariant mass of initial state gluon

        if(eaa.lt.smass)then      ! ->d,u
        if(aa.le.0.5)then
        kf1=1                                      ! g->dd(-)
        else
        kf1=2                                      ! g->uu(-)
        endif
        goto 1003
        endif
       
        if(eaa.ge.smass .and. eaa.lt.cmass)then      ! ->d,u,s
        if(aa.le.0.3333)then
        kf1=1                                      ! g->dd(-)
        elseif(aa.gt.0.3333 .and. aa.lt.0.6666)then
        kf1=2                                      ! g->uu(-)
        else
        kf1=3                                      ! g->ss(-)
        endif
        goto 1003
        endif
        
        if(eaa.ge.cmass .and. eaa.lt.bmass)then        ! d,u,s,c
        if(aa.le.0.25)then
        kf1=1                                        ! ->dd(-)
        elseif(aa.gt.0.25 .and. aa.le.0.5)then
        kf1=2                                        ! ->uu(-)
        elseif(aa.gt.0.5 .and. aa.le.0.75)then
        kf1=3                                        ! ->ss(-)
        else
        kf1=4                                        ! ->cc(-)
        endif
        goto 1003
        endif

        if(eaa.ge.bmass .and. eaa.lt.tmass)then       ! d,u,s,c,b
        if(aa.le.0.2)then
        kf1=1                                        ! ->dd(-)
        elseif(aa.gt.0.2 .and. aa.le.0.4)then
        kf1=2                                        ! ->uu(-)
        elseif(aa.gt.0.4 .and. aa.le.0.6)then
        kf1=3                                        ! ->ss(-)
        elseif(aa.gt.0.6 .and. aa.le.0.8)then
        kf1=4                                        ! ->cc(-)
        else
        kf1=5                                        ! ->bb(-)   
        endif
        goto 1003   
        endif

        if(eaa.ge.tmass)then                          ! d,u,s,c,b,t
        if(aa.le.0.1666)then
        kf1=1                                        ! ->dd(-)
        elseif(aa.gt.0.1666 .and. aa.le.0.3333)then
        kf1=2                                        ! ->uu(-)
        elseif(aa.gt.0.3333 .and. aa.le.0.4998)then
        kf1=3                                        ! ->ss(-)
        elseif(aa.gt.0.4998 .and. aa.le.0.6664)then
        kf1=4                                        ! ->cc(-)
        elseif(aa.gt.0.6664 .and. aa.le.0.833)then
        kf1=5                                        ! ->bb(-)
        else
        kf1=6                                        ! ->tt(-) 
        endif
        endif

c250420
1003	kf2=-kf1   ! 250420
	if(pyr(1).ge.0.5) then
	kf1=kf2
	kf2=-kf2
	endif
200     endif   !!!
	endif   !!
c----------finish

100	continue
	q2max=zab*(qqb-tl0/(1.-zab))
c	the max. value of A1's virtuality decided by the eqution(*) with
c	 qqc=tl0 & pt2=0 .
	call suda_ti(ilo,kf1,ii,q2max,q2,z)
	if(ilo.eq.0) then
	qqa=tl0
	il=0
c	qqa=tl0 means the forward evolution has finished, thus set il=0.
	else
	qqa=q2
	endif
	q2max=(qqb-qqa/zab)*(1.-zab)
c	the max. value of A2's virtuality decided by the eqution(*) with pt2=0 
        call suda_ti(ilo,kf2,ii,q2max,q2,z)
	if(ilo.eq.0) then
	qqc=tl0
	else
	qqc=q2   
	endif
	pt2=(qqb-qqa/zab-qqc/(1.-zab))*(1.-zab)*zab
c	pt2: the square of transverse momentum decided by the equation(*)

c0130	the following block gives the momentum of A1,A2 i.e. pa1,pa2
c	 according to the definition of z & pt.
	pa02=pa(1)**2+pa(2)**2+pa(3)**2
        paz=dsqrt(pa02)   ! third momentum of A 
c	note Eq. (*) with assumption that z axis is set on momentum of A,
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
c0130	finished---------------------------------------

c0140   new induced parton (A2) should be added to particle list.
	nsca=nsca+1
	kpip(nsca)=kf2
	pip(4,nsca)=amass(kf2)**2
	do i=1,3
	pip(4,nsca)=pip(4,nsca)+pa2(i)**2
	pip(i,nsca)=pa2(i)
	enddo
	pip(4,nsca)=dsqrt(pip(4,nsca))
c	pip(5,nsca)=qa2
	pip(5,nsca)=qqc
	pip(6,nsca)=xa2
c0140	finished--------------------------------------------

c------- update five very important values of KFK,QA,XA,PA,EA (time like)
c        & then forward evolution along main chain further
	kfk=kf1   !!
	qa=qqa
	do i=1,3
	pa(i)=pa1(i)
	enddo
	ea=dsqrt(pa(1)**2+pa(2)**2+pa(3)**2+amass(kfk)**2)
	xa=xa1

c080520	if(sampq.gt.10.)return   ! 211004
        if(sampq.ge.1.)return   ! once branching only 080520
	if(il.eq.1) goto 10
c       il=0 means the evolution 'stop' while il=1 'continue'
c       the time like branching of A2 parton will consider later, i.e.
c        in 'collis' after statement numbered 300
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine suda_ti(ilo,kf0,ii,q2max,q2,z)
c	time-like Sudakov factor
c	perform time-like forward branching,one step only
	parameter (msca=20000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/sa24/adj1(40),nnstop,non24,zstop
c----------------------------------------------------------------------
c	process    : A  ->   A1  +  A2
c	virtuality : qqb     qqa    qqc
c
c	equation: qqb=qqa/z+qqc/(1-z)+pt2/(1-z)z
c----------------------------------------------------------------------
	yint1(x)=-dlog(1.d0-x)          !q->qg cf. Sa' note p. 44
	yint2(x)=dlog(1.d0/(1.d0-x)-1.d0)   !g->gg
c	yint1,2: the integral functions of part of the splitting functions	
	adj25=adj1(25)
c	adj25: lamda_QCD
	adj25=adj25*adj25
c	write(9,*)'adj25=',adj25   ! sa
	ilo=0
	if(q2max.le.4*tl0) return   ! cf. W' paper p.264 Eq.25    
	ikf=iabs(kf0)
	ii=0
	tmax=dlog(q2max/adj25)
	tmin=dlog(4.*tl0/adj25)   ! cf. W' paper p.264
c	tmax,tmin: bounds of t=dlog(Q**2/lamda) in the Sudakov form factor
        zmin=tl0/q2max     ! cf. W' paper p.263
	zmax=1-tl0/q2max   ! cf. W' paper p.263
c	approximated values of the bounds of z
	if(zmax.le.zmin) return
	if(ikf.ne.21) then   ! q->qg
	ciqg=4./3.*((2.*yint1(zmax)-zmax-0.5*zmax**2)-
     *             (2.*yint1(zmin)-zmin-0.5*zmin**2)) ! cf. Sa' note p. 43
c	ciqg: the integration of splitting function for q->qg
        cisum=ciqg
	ii=1       ! q->qg
	else       ! g->
	cigg=(6.*yint2(zmax)-12.*zmax+3.*zmax**2-2.*zmax**3)-
     *	     (6.*yint2(zmin)-12.*zmin+3.*zmin**2-2.*zmin**3)     
c	cigg: the integration of splitting function for g->gg
	ciqq=1./2.*((zmax-zmax**2+2./3.*zmax**3)-
     *	            (zmin-zmin**2+2./3.*zmin**3))
c	ciqq: the integration of splitting function for g->qq(-)
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
100	continue
	aaa=pyr(1)
c	tt=tmax*aaa**ce   ! variable alpha_s
	tt=tmax+6.2832/as/cisum*dlog(aaa)   ! constant alpha_s
c	tt: the t value sampled from time-like Sudakov factor 
	if(tt.le.tmin) return
	q2=adj25*dexp(tt)
c	adj25: (\Lambda_s)**2
	rzmax=0.5d0+0.5*dsqrt(1.-4.*tl0/q2)
c       cf. Eq.24 in B. R. Webber, Ann. Rev. Nucl. Part. Sci. 36(1986)253
c        that journal is simplified as W' elsewhere 
	rzmin=1.-rzmax
c	rzmax,rzmin: exact value of zmax & zmin corresponding to the t 
c	 sampled.
200	continue
c0170	sample z.
	aaa=pyr(1)
	bbb=pyr(1)
        if(ii.eq.1) then   ! q->qg
	zsp=aaa*yint1(zmax)+(1.-aaa)*yint1(zmin)   
	zsp=1.d0-dexp(-zsp)
	ratio=(1.+zsp*zsp)/2.   ! accepting probability
c220702	if(ratio.le.bbb) goto 200
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
c0170	sample z finished--------------------------------------------------
	if(zsp.gt.rzmax.or.zsp.lt.rzmin) then
c	if the z sampled fall out of the exact region, reject the sampled t,
c	let tmax equal to the t sampled, and go back to resample.
	sampz=sampz+1.   ! 211004
c	write(9,*)'in suda sampz=',sampz   ! sa
	if(sampz.gt.1000.)then   ! 211004
	ilo=0   ! 211004
	return   ! forced stopping time-like branching 211004
	endif   ! 211004
	tmax=tt
	goto 100
	endif
	z=zsp
	ilo=1
c	write(*,*)ii
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
        parameter (mplis=80000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
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
cc200   write(17,*)'jj=',jj   ! w
200     return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fap_qq(zz)
c	split function for g-->q*q(-)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	fap_qq=(zz**2+(1.-zz)**2)/2.
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fap_gg(zz)
c	split function for g-->g*g
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	fap_gg=6.*(1.-zz*(1.-zz))**2/(zz*(1.-zz))
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fap_qg(zz)
c	split function for q-->q*g
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	fap_qg=4.*(1.+zz**2)/3./(1.-zz)
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fun_ud(xa)
c	the struction function of valence u,d quark
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
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
c	the struction function of d quark
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
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
c	the struction function of u quark
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	fun_u=fun_ud(xa)-fun_d(xa)
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fun_g(xa)
c	the struction function of gluon
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
	aag=1.56-1.71*dta+0.638*dta*dta
	ag=-0.949*dta+0.325*dta*dta
	bg=6.+1.44*dta-1.05*dta*dta
	a1g=9.-7.19*dta+0.255*dta*dta
	a2g=-16.5*dta+10.9*dta*dta
	a3g=15.3*dta-10.1*dta*dta
	fun_g=aag*xa**ag*(1.-xa)**bg*
     c	(1.+a1g*xa+a2g*xa*xa+a3g*xa*xa*xa)/xa
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fun_s(xa)
c	the struction function of sea quark
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
	aas=1.265-1.132*dta+0.293*dta*dta
	ass=-0.372*dta-0.029*dta*dta
	bs=8.05+1.59*dta-0.153*dta*dta
	a1s=6.31*dta-0.273*dta*dta
	a2s=-10.5*dta-3.17*dta*dta
	a3s=14.7*dta+9.8*dta*dta
	fun_s=aas*xa**ass*(1.-xa)**bs*
     c	(1.+a1s*xa+a2s*xa*xa+a3s*xa*xa*xa)/xa/6.
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fun_gal(xa)
c	function of logarithm Gammma function 
c	real*8 coef(0:6)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	dimension coef(0:6)
	data coef/76.18009173,-86.50532033,24.01409822,-1.231739516,
     c	0.00120858003,-5.36382e-6,2.50662827565/
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
c	function of Beta function
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	b=fun_gal(xa)+fun_gal(ya)-fun_gal(xa+ya)	
	fun_beta=dexp(b)
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fs12_0(xt)
c	differential cross section of q1*q2 -> q1*q2, process 1
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
	common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
	adj11=adj1(1) 
	adj20=adj1(20)  
c	write(9,*)'adj11=',adj11   ! sa
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407
	xu=sm-(xt-tcut)-xs   ! 071018  
c	xs: s
c	sm: ma**2+mb**2+mc**2+md**2
c	xu: u
c	xt: t
c120699	ccc=dlog(xt*xu/xs/adj1(25).adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or.adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
        fs12_0=fsabcd*(xs**2+xu**2)*4./xt**2/9.*adj11   ! 210803
c240803	LO pQCD
	return
100	fs12_0=pioa*8./9./xt**2*adj11   !240803
c	limited and regularized (B. Zhang)
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fs11_0(xt)
c	differential cross section of q1*q1 -> q1*q1, process 2
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
	adj20=adj1(20)  
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407   
	xu=sm-(xt-tcut)-xs   ! 071018  
c120699	ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
	fs11_0=fsabcd*(4.*((xs**2+xu**2)/xt**2+(xs**2+xt**2)/xu**2)/9.-
     c	8.*xs**2/27./xu/xt)*adj11   ! 210803
c240803	LO pQCD
	return
100	fs11_0=pioa*8./9./xt**2*adj11   ! 240803 
c	limited and regularized 
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fs12_1(xt)
c	differential cross section of q1*q2(-) -> q1*q2(-), process 3
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
	adj20=adj1(20)  
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407   
	xu=sm-(xt-tcut)-xs   ! 071018  
c120699	ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
	fs12_1=fsabcd*(4.*(xs**2+xu**2)/xt**2/9.)*adj11   ! 210803
c240803	LO pQCD
	return
100	fs12_1=pioa*8./9./xt**2*adj11   ! 240803
c	limited and regularized 
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fs11_1(xt)
c	differential cross section of q1*q1(-) -> q2*q2(-), process 4
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   ! 080820
c240412
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 131212
        if(iparres.eq.0 .or. (iparres.eq.1.and.mstj1_1.eq.7))then
        fs11_1=0.   ! 160110
        return   ! 160110
        endif
c240412
        adj11=adj1(1)
	adj20=adj1(20)  
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407   
	xu=sm-(xt-tcut)-xs   ! 071018  
c120699	ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
	fs11_1=fsabcd*(4.*(xt**2+xu**2)/xs**2/9.)*adj11   ! 210803
c240803	LO pQCD
	return
100	fs11_1=0.   ! pioa*4./9.*(1.+xt*xt/xs*xs)*adj11   ! 240803
c	limited and regularized 
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fs11_2(xt)
c	differential cross section of q1*q1(-) -> q1*q1(-), process 5
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
	adj20=adj1(20)  
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407   
	xu=sm-(xt-tcut)-xs   ! 071018  
c120699	ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
	fs11_2=fsabcd*(4.*((xs**2+xu**2)/xt**2+(xu**2+xt**2)/xs**2)/9.-
     c	8.*xu**2/27./xs/xt)*adj11   ! 210803
c240803	LO pQCD
	return
100	fs11_2=pioa*8./9./xt**2*adj11   !
c	limited and regularized 
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fsqq(xt)
c	differential cross section of q*q(-) -> g*g, process 6
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   ! 080820
c240412
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 131212
        if(iparres.eq.0 .or. (iparres.eq.1.and.mstj1_1.eq.7))then
        fsqq=0.   ! 160110
        return   ! 160110
        endif
c240412
        adj11=adj1(1)
	adj20=adj1(20)  
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407   
	xu=sm-(xt-tcut)-xs   ! 071018  
c120699	ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
	fsqq=fsabcd*(32.*(xu**2+xt**2)/xu/xt/27.-
     c	8.*(xu**2+xt**2)/xs**2/3.)*adj11   ! 210803
c240803	LO pQCD
	return
100	fsqq=-pioa*32./27./xs/xt*adj11   ! 240803
c	limited and regularized 
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fsgg_1(xt)
c	differential cross section of g*g ->q*q(-), process 7
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
c240412
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 131212
c       write(9,*)'iparres=',iparres
        if(iparres.eq.0)then
        fsgg_1=0.   ! 160110
        return   ! 160110
        endif
c240412
        adj11=adj1(1)
	adj20=adj1(20)  
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407   
	xu=sm-(xt-tcut)-xs   ! 071018  
c120699	ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
	fsgg_1=fsabcd*((xu**2+xt**2)/xu/xt/6.-3.*(xu**2+xt**2)/xs**2/8.)
     c	 *adj11   ! 210803
c240803	LO pQCD
	return
100	fsgg_1=-pioa*1./3./xs/xt*adj11   ! 240803
c	limited and regularized 
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fsqg(xt)
c	differential cross section of q*g -> q*g, process 8
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
	adj20=adj1(20)  
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407   
	xu=sm-(xt-tcut)-xs   ! 071018  
c120699	ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
	fsqg=fsabcd*((xu**2+xs**2)/xt**2-4.*(xu**2+xs**2)/xu/xs/9.)
     c	 *adj11   ! 210803
c240803	LO pQCD
	return
100	fsqg=pioa*2./xt**2*adj11   ! 240803 200407
c	limited and regularized 
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fsgg_0(xt)
c	differential cross section of g*g -> g*g, process 9
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (msca=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa,
     c  ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
c210803
        common/sa24/adj1(40),nnstop,non24,zstop   ! 181003
        adj11=adj1(1)
	adj20=adj1(20)  
c210803
c020903	xt=xt-tcut
c071018	xu=sm-xt-xs-2.*tcut   ! 200407   
	xu=sm-(xt-tcut)-xs   ! 071018  
c120699	ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
c120699	ccc=ccc**2
c120699	fsabcd=pio*as**2/xs**2/ccc
	pioa=pio*as**2   ! 240803
	if(adj20.eq.1 .or. adj20.eq.3)goto 100   ! 230405
	fsabcd=pioa/xs**2   ! 120699
	fsgg_0=fsabcd*(9.*(3.-xu*xt/xs**2-xu*xs/xt**2-xs*xt/xu**2)/2.)
     c	 *adj11   ! 210803
c240803	LO pQCD
	return
100	fsgg_0=pioa*9./2./xt**2*adj11   ! 240803
c	limited and regularized 
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function amass(kf)
c	mass of partons
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
c	calculate distribution function of external function 'fun'
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
c	calculate the relative integral distribution function 
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	dimension dd(0:idw),eee(0:idw)
	dd(0)=0.
	do i=1,idw
	dd(i)=dd(i-1)+eee(i)+eee(i-1)
c280803	dd(i)=dd(i-1)+dabs(eee(i))+dabs(eee(i-1))
	enddo
	do i=1,idw
	dd(i)=dd(i)/dd(idw)
	enddo
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine samp_integ(idw,ux,dx,xf,dd)
c	samples a value from relative integral distribution function 'dd',
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
10	continue
	enddo
20	continue
	xf=pyr(1)
	xf=dx+del*(ir-1)+xf*del
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine eloss(dell,rpo,ppo)
c	consider parton energy loss phenomenologically
c	dell: energy loss per unit distance 
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (mplis=80000,msca=20000)
	common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
	common/sa24/adj1(40),nnstop,non24,zstop
	common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)   
	dimension rpo(4,mplis),ppo(4,mplis)   
c	rpo,ppo: parton four coordinate before Newton motion, four momentum 
c	 before energy loss
	dimension rr(3),pi(3),pia(3),deltp(3)
c	ppsa(1): cumulate the px losed in an event
c	ppsa(2): cumulate the py losed in an event
c	ppsa(3): cumulate the pz losed in an event
c	ppsa(4): cumulate the e losed in an event
	adj138=adj1(38)
c	write(9,*)'adj138=',adj1(38)   ! sa
	iprlo=iprl
	do i=1,iprlo   ! 1
	pei=pp(4,i)   ! energy before loss
	ppo(4,i)=pei
	do j=1,3
	rr(j)=rp(j,i)-rpo(j,i)
	pi(j)=pp(j,i)   ! three momentum before loss
	ppo(j,i)=pi(j)
	enddo
	pt=pi(1)*pi(1)+pi(2)*pi(2)
	if(pt.gt.1.e20)pt=1.e20
	if(pt.lt.1.e-20)pt=1.e-20
	pt=dsqrt(pt)
	if(pt.le.adj138)goto 100
c	below ptmin='adj138' Gev/c jet can not loss energy
	rm=rr(1)*rr(1)+rr(2)*rr(2)+rr(3)*rr(3)
	if(rm.gt.1.e20)rm=1.e20
	if(rm.lt.1.e-20)rm=1.e-20
	rm=dsqrt(rm)
	delte=rm*dell
	if(delte.ge.pei)goto 100
c	'delte' should be less than energy of particle
	if(pei.lt.1.e-20)pei=1.e-20
	srtf=1.-delte/pei
	pp(4,i)=pei-delte   ! energy after loss
	do j=1,3
	pia(j)=srtf*pi(j)   ! three momentum after loss
c	massless and collinear approximations
	deltp(j)=pi(j)-pia(j)   ! three momentum loss
	pp(j,i)=pia(j)   ! three momentum after loss
	ppsa(j)=ppsa(j)+pi(j)-pia(j)   ! cumulate three momentum loss
	enddo
	ppsa(4)=ppsa(4)+delte   ! cumulate energy loss
c	tansfer losed energy to a gluon (added to the end of particle list) 
	call eto1g(i,delte,deltp)
c	tansfer losed energy to two gluons (added to the end of particle list) 
c	call eto2g(i,delte,deltp)
100	continue
	enddo   ! 1
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ctlcre_para(iijk,time)   ! 120505
c       creates the collision time list after energy loss
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        parameter (mplis=80000,mclis=280000)
	common/sa24/adj1(40),nnstop,non24,zstop   ! 210803 181003 161104
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/collist/lc(2,mclis),tc(2,mclis),icol
	dddt=adj1(19)   ! 161104
c	adj136=adj1(36)
        icol=1
	do 100 i=2,iprl   ! upper diagonal
	kfi=iabs(idp(i))
	if(kfi.gt.3 .and. kfi.ne.21)goto 100   ! 120620
c       consider d,u,s their antiquarks, and gluon only   ! 120620
c       if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
	do 200 j=1,i-1   ! upper diagonal
	kfj=iabs(idp(j))
	if(kfj.gt.3 .and. kfj.ne.21)goto 100   ! 120620
c       consider d,u,s, their antiquarks, and gluon only   ! 120620
c       if(kfi.gt.6 .and. kfi.ne.21)goto 100   ! 080520
c       consider d,u,s,c,b,t, their antiquarks, and gluon only   ! 080520
	if(icol.gt.mclis) then
        write(9,*)'icol over limit iprl,icol=',iprl,icol   ! sa
        stop 77777
        endif
	tc(1,icol)=0.0
        tc(2,icol)=0.0
        call tcolij_par(i,j,time,icol)
        tc1=tc(1,icol)
        tc2=tc(2,icol)
        tcicol=tc1
        if(tc1.gt.tc2)tcicol=tc2
	if(tcicol.le.0.0)then
c       write(9,*)'i1,j1,tcicol=',i1,j1,tcicol   ! sa
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
200	continue
100	continue
	if(tcicol.le.0.0) icol=icol-1
        do i=icol+1,mclis
        lc(1,i)=0
        lc(2,i)=0
        tc(1,i)=0.
        tc(2,i)=0.
        enddo
	iprl0=iprl
	return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine eto1g(ii,delte,deltp)
c	transfer the losed energy to a gluon (added to the end of 
c	 particle list)
c	ii: the line number (in 'parlist') of parton which losed energy 
c	delte: energy losed
c	deltp(3): three momentum losed
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (mplis=80000)
	common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	dimension deltp(3)
	iprll=iprl+1
	if(iprll.gt.mplis)then
	write(9,*)'in parcas, mplis needs to be enlarged iprll=',iprll   ! sa
	stop 9999   
	endif
	pp4=delte
	pp(4,iprl+1)=pp4
	do i=1,3
	ppi=deltp(i)
	pp(i,iprl+1)=ppi
	vp(i,iprl+1)=ppi/pp4
	enddo
	taup(iprl+1)=0.
	rmp(iprl+1)=0.
	idp(iprl+1)=21
	rp(4,iprl+1)=rp(4,ii)
	do i=1,3
        rp(i,iprl+1)=pyr(1)*rp(i,ii)
	enddo
	iprl=iprl+1
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine eto2g(ii,delte,deltp)
c	transfer the losed energy to two gluons (added to the end of 
c	 particle list)
c	ii: the line number (in 'parlist') of parton which losed energy 
c	delte: energy losed
c	deltp(3): three momentum losed
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
	parameter (mplis=80000)
	common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
	dimension deltp(3)
	iprll=iprl+1
	if(iprll.gt.mplis)then
	write(9,*)'in parcas, mplis needs to be enlarged iprll=',iprll   ! sa
	stop 9999
	endif
	iprll=iprl+2
	if(iprll.gt.mplis)then
        write(9,*)'in parcas, mplis needs to be enlarged iprll=',iprll   ! sa
        stop 99999
        endif
	pp4=pyr(1)*delte
	if(pp4.lt.1.e-20)pp4=1.e-20
	pp(4,iprl+1)=pp4
	pp42=delte-pp4
	if(pp42.lt.1.e-20)pp42=1.e-20
	pp(4,iprl+2)=pp42
	do i=1,3
	ppi=pyr(1)*deltp(i)
	pp(i,iprl+1)=ppi
	ppi2=deltp(i)-ppi
	pp(i,iprl+2)=ppi2
	vp(i,iprl+1)=ppi/pp4
	vp(i,iprl+2)=ppi2/pp42
	enddo
	taup(iprl+1)=0.
	taup(iprl+2)=0.
	rmp(iprl+1)=0.
        rmp(iprl+2)=0.
	idp(iprl+1)=21
        idp(iprl+2)=21
	rp(4,iprl+1)=rp(4,ii)
        rp(4,iprl+2)=rp(4,ii)
	do i=1,3
        rpi=pyr(1)*rp(i,ii)
        rp(i,iprl+1)=rpi
        rp(i,iprl+2)=rp(i,ii)-rpi
	enddo
	iprl=iprl+2
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coend(ii)   ! 160110
c       copy parton ii to the end of 'parlist'
c       copy parton ii to the end of 'pyjets'
        parameter (mplis=80000,kszj=80000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   
        common/parlist/rp(4,mplis),pp(4,mplis), 
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis) 
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   
        iprl=iprl+1
        do i1=1,3
        rp(i1,iprl)=rp(i1,ii)
        pp(i1,iprl)=pp(i1,ii)
        vp(i1,iprl)=vp(i1,ii)
        enddo
        rp(4,iprl)=rp(4,ii)
        pp(4,iprl)=pp(4,ii)
        taup(iprl)=taup(ii)
        rmp(iprl)=rmp(ii)
        idp(iprl)=idp(ii)
        n=n+1
        do i3=1,5
        k(n,i3)=k(ii,i3)
        p(n,i3)=p(ii,i3)
        v(n,i3)=v(ii,i3)
        enddo
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
        subroutine adjst(ik1,ik1str,ik1sa,ik1sv)   ! 160110
c       finds order number etc. of ik1 string    ! 270620
        parameter(kszj=80000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620
        ik1str=0   ! if ik1 not belong to string    
        ik1sa=0   ! if ik1 not belong to string
        ik1sv=0   ! if ik1 not belong to string
        if(ik1.le.nbe)then   ! nbe is giving before call break 
        do i1=1,nstr0   ! nstr0 is giving after call break 
        i1a=nstr1a(i1)   ! i1a: line number of 'A' of i1-th string
        i1v=nstr1v(i1)   ! i1v: line number of 'V' of i1-th string
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
        parameter (kszj=80000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
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
        parameter (kszj=80000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
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
c       updates 'parlist' one step downward from ii+1 to iprl 
c       updates 'pyjets' one step downward from ii+1 to n
c       updates 'sbe' one step downward from ii+1 to nbe if ii.le.nbe
c       updates diquark list
c       updates string list
c       updates collision time list
c       ik1 & ik2 are the line number (in parton list) of colliding pair
        parameter (mplis=80000,kszj=80000)
        parameter (mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/collist/lc(2,mclis),tc(2,mclis),icol
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
        ktrs(ntrs,2)=idp(ii)
        do j=1,4
        ptrs(ntrs,j)=pp(j,ii)
        vtrs(ntrs,j)=rp(j,ii)
        enddo
        ptrs(ntrs,5)=rmp(ii)
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
        do i2=ii+1,iprl
        i3=i2-1
        ndiq(i3)=ndiq(i2)
        enddo  

c       updates particle list 'parlist' 
        if(ii.eq.iprl)iprl=iprl-1
        if(ii.lt.iprl)then
        do i2=ii+1,iprl
        i3=i2-1
        do i1=1,3
        rp(i1,i3)=rp(i1,i2)
        pp(i1,i3)=pp(i1,i2)
        vp(i1,i3)=vp(i1,i2)
        enddo
        rp(4,i3)=rp(4,i2)
        pp(4,i3)=pp(4,i2)
        taup(i3)=taup(i2)
        rmp(i3)=rmp(i2)
        idp(i3)=idp(i2)
        enddo
        iprl=iprl-1

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
        enddo
        n=n-1
        endif
c       if(ik1.gt.ii)ik1=ik1-1   ! 280620
c       if(ik2.gt.ii)ik2=ik2-1   ! 280620
        endif

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
        lc2=lc(2,m)
        if(lc2.ge.ii)lc(2,m)=lc2-1
        enddo

        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine strmov(istr1,istr1a,istr1v,ik1,ik2,lmn)  ! 280620 070720
c       moves all conponents of 'istr1'-th string out 
c       istr1: order number of string (in string list) to be moved
c       istr1a: line number (in parton list) of first component of above string 
c       istr1v: line number (in parton list) of last component of above string
c       ik1 & ik2 are line # of collidin pair in parton list 'parlist' 
        parameter (mplis=80000,kszj=80000)
        parameter (mclis=280000)
	IMPLICIT DOUBLE PRECISION(A-H, O-Z)
	IMPLICIT INTEGER(I-N)
	INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/parlist/rp(4,mplis),pp(4,mplis),
     c  taup(mplis),rmp(mplis),vp(3,mplis),iprl,idp(mplis)
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620
c       moves components of istr1-th string out of parton list  
        do ii=istr1v,istr1a,-1
        call parmov(ii,ik1,ik2,lmn)   ! 070720
        enddo
c       moves string istr1-th out of string list
        if(istr1.eq.nstr1)nstr1=nstr1-1
        if(istr1.lt.nstr1)then
        jj=istr1v-istr1a+1
        do ii=istr1+1,nstr1
        if(ii.ge.istr1)then
        nstr1v(ii-jj)=nstr1v(ii)
        nstr1a(ii-jj)=nstr1a(ii)
        endif
        enddo
        nstr1=nstr1-1
        endif
        return
        end



cccccccccccccccccccccccccccccccc end ccccccccccccccccccccccccccccc
