C***********************************************************************
	subroutine ctrsc(x,n,p,xbar,sdv,mvcode)
C Centers and scales the data matrix so that the observed data in every
C column have mean zero and variance one. If a column has zero variance
C or less than 2 observations then the data are centered (set equal
C to zero)
	integer n,p,count
        double precision x(n,p),mvcode
	double precision xbar(p),sdv(p),sum1,sum2
	do 10 j=1,p
           sum1=0
           sum2=0
           count=0
           do 5 i=1,n
              if(x(i,j).ne.mvcode) then
                 count=count+1
                 sum1=sum1+x(i,j)
                 sum2=sum2+x(i,j)**2
              endif
5	   continue
           if(count.gt.0) then
              xbar(j)=sum1/count
              sdv(j)=sqrt((sum2-(sum1**2)/count)/count)
              do 7 i=1,n
                 if(x(i,j).ne.mvcode) x(i,j)=x(i,j)-xbar(j)
7             continue
              if(sdv(j).gt.0.d0) then
                 do 9 i=1,n
                    if(x(i,j).ne.mvcode) x(i,j)=x(i,j)/sdv(j)
9                continue
               else
                 sdv(j)=1.d0
              endif
           else
              sdv(j)=1.d0
           endif
10	continue
	return
	end
C***********************************************************************
        subroutine mkpsi(p,psi)
C Generates a symmetric matrix of integers indicating the linear
C position in packed storage of the matrix elements
        integer p,psi(0:p,0:p),posn
        posn=0
        do 10 j=0,p
           posn=posn+1
           psi(j,j)=posn
           do 5 k=j+1,p
              posn=posn+1
              psi(j,k)=posn
              psi(k,j)=posn
5          continue
10      continue
        return
        end
C***********************************************************************
        subroutine swpm(q,psi,npsi,sigma,ncells,mu,p,pivot,
     /     submat,dir,what)
C Performs sweep on parameters of mixed normal-categorical model.
C Sweeps on pivot position. Sweeps only the (1:submat,1:submat)
C submatrix. If dir=1, performs ordinary sweep. If dir=-1, performs
C reverse sweep. Skips over any structural zero cells having
C a value in p of -999.0. If what=1, does full sweep. If what=0,
C sweeps only sigma. If what=2, sweeps only sigma and mu.
        integer q,psi(q,q),npsi,ncells,pivot,submat,dir
        integer what
        double precision sigma(npsi),mu(q,ncells),p(ncells),a,b,c
        a=sigma(psi(pivot,pivot))
        sigma(psi(pivot,pivot))=-1.d0/a
        do 10 j=1,submat
           if(j.ne.pivot)sigma(psi(j,pivot))=sigma(psi(j,pivot))/a*dir
10      continue
        if(what.ge.1)then
        do 11 j=1,ncells
           if(p(j).ne.-999.0d0) mu(pivot,j)=mu(pivot,j)/a*dir
11      continue
        endif
        do 30 i=1,submat
           if(i.ne.pivot)then
              b=sigma(psi(i,pivot))
              do 20 j=i,submat
                 if(j.ne.pivot)then
                    c=sigma(psi(j,pivot))
                    sigma(psi(i,j))=sigma(psi(i,j))-a*b*c
                 endif
20            continue
              if(what.ge.1)then
              do 21 j=1,ncells
                 if(p(j).ne.-999.0d0)then
                    c=mu(pivot,j)
                    mu(i,j)=mu(i,j)-a*b*c
                 endif
21            continue
              endif
           endif
30      continue
        if(what.eq.1)then
        do 40 j=1,ncells
           if(p(j).ne.-999.0d0)then
              c=mu(pivot,j)
              p(j)=p(j)-a*c*c
           endif
40      continue
        endif
        return
        end
C***********************************************************************
        subroutine initm(q,npsi,t1,ncells,t2,t3)
C Initializes workspace
        integer q,npsi,ncells
        double precision t1(npsi),t2(q,ncells),t3(ncells)
        do 1 i=1,npsi
           t1(i)=0
1       continue
        do 3 j=1,ncells
           t3(j)=0
           do 2 i=1,q
              t2(i,j)=0
2          continue
3       continue
        return
        end
C***********************************************************************
        subroutine tobsm(q,psi,npsi,t1,ncells,t2,t3,npattz,
     /    rz,mdpzgrp,npattw,p,rw,mdpwgrp,ngrp,mobs,mobsst,
     /    nmobs,n,z,ocw,ocz)
C Tabulates the known part of the sufficient stats for all missingness
C patterns.
        integer q,p,psi(q,q),npsi,ncells,npattz,rz(npattz,q)
        integer mdpzgrp(npattz),npattw,rw(npattw,p)
        integer mdpwgrp(npattw),ngrp
        integer mobs(ngrp),mobsst(ngrp),nmobs(ngrp),n
        integer ocz(q),ocw(p),nocz,nocw,pattz,pattw,grpno,a,b
        double precision z(n,q)
        double precision t1(npsi),t2(q,ncells),t3(ncells)
        call initm(q,npsi,t1,ncells,t2,t3)
        pattw=0
        grpno=0
        do 80 pattz=1,npattz
           call gtoc(q,npattz,rz,pattz,ocz,nocz,q)
           do 70 a=1,mdpzgrp(pattz)
              pattw=pattw+1
              call gtoc(p,npattw,rw,pattw,ocw,nocw,p)
              do 60 b=1,mdpwgrp(pattw)
                 grpno=grpno+1
                 if(nocw.eq.p) t3(mobs(grpno))=t3(mobs(grpno))+
     /                nmobs(grpno)
                    do 50 i=mobsst(grpno),(mobsst(grpno)+
     /                   nmobs(grpno)-1)
                       do 20 j=1,nocz
                          if(nocw.eq.p) t2(ocz(j),mobs(grpno))=
     /                         t2(ocz(j),mobs(grpno))+z(i,ocz(j))
                          do 10 k=j,nocz
                             t1(psi(ocz(j),ocz(k)))=
     /                           t1(psi(ocz(j),ocz(k)))+
     /                           z(i,ocz(j))*z(i,ocz(k))
10                        continue
20                     continue
50                  continue
60            continue
70         continue
80      continue
        return
        end
C************************************************************************
        subroutine gtmc(p,npatt,r,patt,mc,nmc,last)
C Finds the column numbers of the missing variables, and stores them
C in the first nmc elements of mc. Does not go beyond column=last.
        integer p,npatt,r(npatt,p),patt,mc(p),nmc,last
        nmc=0
        do 10 j=1,last
           if(r(patt,j).eq.0)then
              nmc=nmc+1
              mc(nmc)=j
           endif
10      continue
        return
        end
C************************************************************************
        subroutine gtoc(p,npatt,r,patt,oc,noc,last)
C Finds the column numbers of the observed variables, and stores them
C in the first noc elements of oc. Does not go beyond column=last.
        integer p,npatt,r(npatt,p),patt,oc(p),noc,last
        noc=0
        do 10 j=1,last
           if(r(patt,j).eq.1)then
              noc=noc+1
              oc(noc)=j
           endif
10      continue
        return
        end
C***********************************************************************
        subroutine stvlm(q,psi,npsi,t1,ncells,t2)
C Creates starting values mu=0 and sigma=identity.
        integer q,psi(q,q),npsi,ncells
        double precision t1(npsi),t2(q,ncells)
        do 1 i=1,npsi
           t1(i)=0
1       continue
        do 2 j=1,q
           t1(psi(j,j))=1
2       continue
        do 4 j=1,ncells
           do 3 i=1,q
              t2(i,j)=0
3          continue
4       continue
        return
        end
C***********************************************************************
        subroutine seteqm(q,npsi,ncells,sigma1,mu1,pi1,
     /       sigma2,mu2,pi2)
C sets theta2 equal to theta1
        integer q,npsi,ncells
        double precision sigma1(npsi),mu1(q,ncells),pi1(ncells)
        double precision sigma2(npsi),mu2(q,ncells),pi2(ncells)
        do 1 i=1,npsi
           sigma2(i)=sigma1(i)
1       continue
        do 3 j=1,ncells
           pi2(j)=pi1(j)
           do 2 i=1,q
              mu2(i,j)=mu1(i,j)
2          continue
3       continue
        return
        end
C***********************************************************************
        subroutine swpobsm(q,psi,npsi,ncells,sigma,mu,pii,npattz,rz,
     /       pattz,logdet,what)
C sweeps theta to condition on the observed variables in z. Updates
C determinant in det. The initial value of det should be 1.
C If what=1, does full parameter sweep. If what=0, sweeps
C only sigma.
        integer q,psi(q,q),npsi,ncells,npattz,rz(npattz,q),pattz
        integer what
        double precision sigma(npsi),mu(q,ncells),pii(ncells),logdet
        do 10 j=1,q
           if((rz(pattz,j).eq.1).and.(sigma(psi(j,j)).gt.0))then
              logdet=logdet+log(sigma(psi(j,j)))
              call swpm(q,psi,npsi,sigma,ncells,mu,pii,j,q,1,what)
           elseif((rz(pattz,j).eq.0).and.(sigma(psi(j,j)).lt.0))then
              call swpm(q,psi,npsi,sigma,ncells,mu,pii,j,q,-1,what)
              logdet=logdet-dlog(sigma(psi(j,j)))
           endif
10      continue
        return
        end
C***********************************************************************
        subroutine tk2log(ncells,pii)
C takes 2*log of cell probs, and sets zero cells to -999.0.
        integer ncells
        double precision pii(ncells)
        do 1 i=1,ncells
           if(pii(i).gt.0)then
              pii(i)=2.d0*log(pii(i))
           elseif(pii(i).eq.0)then
              pii(i)=-999.0d0
           endif
1       continue
        return
        end
C***********************************************************************
        subroutine estepm(q,psi,npsi,ncells,sigma,mu,pii,kn1,kn2,kn3,
     /       t1,t2,t3,npattz,rz,mcz,ocz,mdpzgrp,npattw,p,rw,mcw,
     /       mdpwgrp,ngrp,mobs,mobsst,nmobs,n,z,d,jmp,c,theta)
        integer q,psi(q,q),npsi,ncells,npattz,pattz,rz(npattz,q)
        integer mcz(q),ocz(q),nmcz,nocz,pattw,grpno,a,b,n,ngrp
        integer mdpzgrp(npattz),npattw,p,rw(npattw,p),mcw(p)
        integer nmcw,mdpwgrp(npattw),dmis,c(p),d(p),jmp(p)
        integer mobs(ngrp),mobsst(ngrp),nmobs(ngrp)
        double precision z(n,q)
        double precision sigma(npsi),mu(q,ncells),pii(ncells)
        double precision kn1(npsi),kn2(q,ncells),kn3(ncells)
        double precision t1(npsi),t2(q,ncells),t3(ncells)
        double precision theta(ncells),logdet
        logdet=0.0d0
        pattw=0
        grpno=0
        call seteqm(q,npsi,ncells,kn1,kn2,kn3,t1,t2,t3)
        call tk2log(ncells,pii)
        do 200 pattz=1,npattz
           call swpobsm(q,psi,npsi,ncells,sigma,mu,pii,npattz,rz,
     /          pattz,logdet,1)
           call gtmc(q,npattz,rz,pattz,mcz,nmcz,q)
           call gtoc(q,npattz,rz,pattz,ocz,nocz,q)
           do 180 a=1,mdpzgrp(pattz)
              pattw=pattw+1
              call gtmc(p,npattw,rw,pattw,mcw,nmcw,p)
              call gtdmis(p,d,mcw,nmcw,dmis)
              do 170 b=1,mdpwgrp(pattw)
                 grpno=grpno+1
                 do 150 i=mobsst(grpno),(mobsst(grpno)+nmobs(grpno)-1)
                    call gtprob(q,ncells,mu,pii,n,z,i,p,mcw,nmcw,c,
     /                   d,jmp,dmis,mobs(grpno),ocz,nocz,theta)
                    if(nmcw.eq.0)then
                       call addstat1(q,psi,npsi,ncells,sigma,mu,
     /                      theta,t1,t2,t3,n,z,i,p,mcw,nmcw,c,d,jmp,
     /                      dmis,mobs(grpno),ocz,nocz,mcz,nmcz)
                    else
                       call addstat2(q,psi,npsi,ncells,sigma,mu,
     /                      theta,t1,t2,t3,n,z,i,p,mcw,nmcw,c,d,jmp,
     /                      dmis,mobs(grpno),ocz,nocz,mcz,nmcz)
                    endif
150              continue
170           continue
180        continue
200     continue
        return
        end
C***********************************************************************
        subroutine mstepm(q,psi,npsi,ncells,t1,t2,t3,n,prior)
C converts t1, t2, t3 to ML estimates
        integer q,psi(q,q),npsi,ncells,m,n
        double precision t1(npsi),t2(q,ncells),t3(ncells),sum
	double precision prior(ncells)
        do 30 j=1,q
           do 20 k=j,q
              sum=0.0d0
              do 10 m=1,ncells
                 if(t3(m).ne.0) sum=sum+t2(j,m)*t2(k,m)/t3(m)
10            continue
              t1(psi(j,k))=(t1(psi(j,k))-sum)/dble(n)
20         continue
30      continue
	sum=0d0
        do 50 m=1,ncells
	   if(prior(m).ne.-999.0d0) sum=sum+t3(m)+prior(m)-1d0
           if(t3(m).ne.0)then
              do 40 j=1,q
                 t2(j,m)=t2(j,m)/t3(m)
40            continue
           endif
50      continue
	do 60 m=1,ncells
	   if(prior(m).ne.-999.0d0) t3(m)=(t3(m)+prior(m)-1d0)/sum
 60	continue
        return
        end
C***********************************************************************
        subroutine istepm(q,psi,npsi,ncells,sigma,mu,pii,kn1,kn2,kn3,
     /       t1,t2,t3,npattz,rz,mcz,ocz,mdpzgrp,npattw,p,rw,mcw,
     /       mdpwgrp,ngrp,mobs,mobsst,nmobs,n,z,d,jmp,c,theta,chf,w,zz)
        integer q,psi(q,q),npsi,ncells,npattz,pattz,rz(npattz,q)
        integer mcz(q),ocz(q),nmcz,nocz,pattw,grpno,a,b,n,ngrp
        integer mdpzgrp(npattz),npattw,p,rw(npattw,p),mcw(p)
        integer nmcw,mdpwgrp(npattw),dmis,c(p),d(p),jmp(p)
        integer mobs(ngrp),mobsst(ngrp),nmobs(ngrp),w(n,p)
        double precision z(n,q),junk
        double precision sigma(npsi),mu(q,ncells),pii(ncells)
        double precision t1(npsi),t2(q,ncells),t3(ncells)
        double precision kn1(npsi),kn2(q,ncells),kn3(ncells)
        double precision theta(ncells),logdet,chf(npsi),zz(q)
        junk=dble(gauss())
        logdet=0.0d0
        pattw=0
        grpno=0
        call tk2log(ncells,pii)
        call seteqm(q,npsi,ncells,kn1,kn2,kn3,t1,t2,t3)
        do 200 pattz=1,npattz
           call swpobsm(q,psi,npsi,ncells,sigma,mu,pii,npattz,rz,
     /          pattz,logdet,1)
           call gtmc(q,npattz,rz,pattz,mcz,nmcz,q)
           call gtoc(q,npattz,rz,pattz,ocz,nocz,q)
           call sigexm(npsi,sigma,chf,q,psi,mcz,nmcz)
           call cholsm(npsi,chf,q,psi,mcz,nmcz)
           do 180 a=1,mdpzgrp(pattz)
              pattw=pattw+1
              call gtmc(p,npattw,rw,pattw,mcw,nmcw,p)
              call gtdmis(p,d,mcw,nmcw,dmis)
              do 170 b=1,mdpwgrp(pattw)
                 grpno=grpno+1
                 do 150 i=mobsst(grpno),(mobsst(grpno)+nmobs(grpno)-1)
                    call gtprob(q,ncells,mu,pii,n,z,i,p,mcw,nmcw,c,
     /                   d,jmp,dmis,mobs(grpno),ocz,nocz,theta)
                    call istepim(q,psi,npsi,ncells,sigma,mu,
     /                   theta,t1,t2,t3,n,z,i,p,mcw,nmcw,c,d,jmp,dmis,
     /                   mobs(grpno),ocz,nocz,mcz,nmcz,chf,zz,w)
150              continue
170           continue
180        continue
200     continue
        return
        end
C***********************************************************************
        subroutine istepim(q,psi,npsi,ncells,sigma,mu,theta,t1,t2,
     /       t3,n,z,i,p,mcw,nmcw,c,d,jmp,dmis,mobs,ocz,nocz,mcz,nmcz,
     /       chf,zz,w)
C Draws missing data for unit i and increments sufficient statistics.
        integer q,psi(q,q),npsi,ncells,n,i,p,mcw(p),nmcw,d(p)
        integer jmp(p),dmis,mobs,ocz(q),nocz,mcz(q),nmcz,mmis,m,a
        integer c(p),w(n,p)
        double precision sigma(npsi),mu(q,ncells),theta(ncells)
        double precision t1(npsi),t2(q,ncells),t3(ncells),sum
        double precision chf(npsi),u,sum1,sum2,zz(q)
        double precision z(n,q)
        call initc(p,c,mcw,nmcw)
        mmis=0
        u=dble(rangen(0))
        sum1=0.0d0
        do 200 a=1,dmis
           if(a.ne.1)then
              call advc(p,c,d,mcw,nmcw)
              call gtmmis(p,c,mcw,nmcw,jmp,mmis)
           endif
           m=mobs+mmis
           if(theta(m).ne.-999.0d0)then
              sum1=sum1+theta(m)
              if((sum1.ge.u).or.(a.eq.dmis))then
                 if(nmcw.gt.0) t3(m)=t3(m)+1.0d0
                 do 5 j=1,nmcw
                    w(i,mcw(j))=c(mcw(j))
5                continue
                 do 20 j=1,nmcz
                    sum=mu(mcz(j),m)
                    do 10 k=1,nocz
                       sum=sum+sigma(psi(mcz(j),ocz(k)))*z(i,ocz(k))
10                  continue
                    z(i,mcz(j))=sum
20               continue
                 do 100 j=1,nmcz
                    zz(mcz(j))=dble(gauss())
                    sum2=0.0d0
                    do 80 k=1,j
                       sum2=sum2+zz(mcz(k))*chf(psi(mcz(j),mcz(k)))
80                  continue
                    z(i,mcz(j))=z(i,mcz(j))+sum2
100              continue
                 if(nmcw.eq.0)then
                    do 105 j=1,nmcz
                       t2(mcz(j),m)=t2(mcz(j),m)+z(i,mcz(j))
 105                continue
                 else
                    do 108 j=1,q
                       t2(j,m)=t2(j,m)+z(i,j)
 108                continue
                 endif
                 do 130 j=1,nmcz
                    do 110 k=1,nocz
                       t1(psi(mcz(j),ocz(k)))=t1(psi(mcz(j),ocz(k)))
     /                   + z(i,mcz(j))*z(i,ocz(k))
 110                continue
                    do 120 k=1,j
                       t1(psi(mcz(j),mcz(k)))=t1(psi(mcz(j),mcz(k)))
     /                   + z(i,mcz(j))*z(i,mcz(k))
 120                continue
 130             continue
              goto 210
              endif
           endif
200     continue
210     continue
        return
        end
C***********************************************************************
        subroutine pstepm(q,psi,npsi,ncells,t1,t2,t3,n,p,prior,chf,
     /     mx,zz,mcz,err)
C Converts t1, t2, t3 to draws from the complete-data posterior.
C If the observed count for any non-structural zero cell is zero, the
C posterior is not proper and err is set to 1.
        integer q,psi(q,q),npsi,ncells,p,mcz(q),n,m
        double precision err,junk
        double precision t1(npsi),t2(q,ncells),t3(ncells),df
        double precision prior(ncells),zz(q),sum,mx(q,q),chf(npsi)
        junk=dble(gauss())
        err=0
        df=dble(n)
        do 1 m=1,ncells
           if(prior(m).ne.-999.0d0)then
              if(t3(m).le.0d0)then
                 err=1
                 goto 200
              else
                 df=df-1d0
              endif
           endif
 1      continue
C calculate least-squares quantities
        do 30 j=1,q
           do 20 k=j,q
              sum=0d0
              do 10 m=1,ncells
                 if(prior(m).ne.-999.0d0)then
                    sum=sum+t2(j,m)*t2(k,m)/t3(m)
                 endif
 10           continue
              t1(psi(j,k))=t1(psi(j,k))-sum
 20        continue
 30     continue
        do 50 m=1,ncells
           if(prior(m).ne.-999.0d0)then
              do 40 j=1,q
                 t2(j,m)=t2(j,m)/t3(m)
 40           continue
           endif
 50     continue
C draw mu and sigma
        do 55 j=1,q
           mcz(j)=j
 55     continue
        call cholsm(npsi,t1,q,psi,mcz,q)
        call bfacm(npsi,chf,q,psi,df)
        call invtrm(npsi,chf,q,psi)
        call mmnm(npsi,chf,t1,q,psi,mx)
        do 80 m=1,ncells
           if(prior(m).ne.-999.0d0)then
              do 58 k=1,q
                 zz(k)=dble(gauss())
 58           continue
              do 70 i=1,q
                 sum=0d0
                 do 60 k=1,q
                    sum=sum+mx(k,i)*zz(k)
 60              continue
                 t2(i,m)=t2(i,m)+sum/dsqrt(t3(m))
 70           continue
           endif
 80     continue
        do 100 i=1,q
           do 90 j=i,q
              sum=0d0
              do 85 k=1,q
                 sum=sum+mx(k,i)*mx(k,j)
 85           continue
              t1(psi(i,j))=sum
 90        continue
 100    continue
C draw cell probs
        sum=0d0
        do 105 m=1,ncells
           if(prior(m).eq.-999.0d0)then
              t3(m)=0d0
           else
              df=df-1d0
              t3(m)=dble(gamm(sngl(t3(m)+prior(m))))
              sum=sum+t3(m)
           endif
 105      continue
        do 108 m=1,ncells
           t3(m)=t3(m)/sum
 108    continue
 200    continue
        return
        end
C***********************************************************************
        subroutine mstepcm(q,psi,npsi,ncells,t1,t2,t3,sigma,mu,n,
     /       r,design,wk,mcr,psir,npsir,wkr,wkd,beta)
        integer q,psi(q,q),npsi,ncells,m,n,r,mcr(r),npsir
        integer psir(r,r),mi,h
        double precision t1(npsi),t2(q,ncells),t3(ncells),sum
        double precision sigma(npsi),mu(q,ncells),beta(r,q)
        double precision design(ncells,r),wk(npsir),wkr(r),wkd(ncells)
C calculate least-squares quantities
        do 30 j=1,r
           do 20 k=j,r
              sum=0d0
              do 10 m=1,ncells
                 sum=sum+design(m,j)*design(m,k)*t3(m)
10            continue
              wk(psir(j,k))=sum
20         continue
30      continue
        call invsym(r,psir,npsir,wk,mcr)
        do 100 h=1,r
           do 60 mi=1,ncells
              sum=0d0
              do 50 j=1,r
                 sum=sum+wk(psir(h,j))*design(mi,j)
50            continue
              wkd(mi)=sum
60         continue
           do 90 k=1,q
              sum=0d0
              do 80 mi=1,ncells
                 sum=sum+wkd(mi)*t2(k,mi)
80            continue
              beta(h,k)=sum
90         continue
100     continue
        do 150 j=1,q
           do 110 h=1,r
              sum=0d0
              do 105 m=1,ncells
                 sum=sum+t2(j,m)*design(m,h)
105           continue
              wkr(h)=sum
110        continue
           do 140 k=j,q
              sum=0d0
              do 120 h=1,r
                 sum=sum+wkr(h)*beta(h,k)
120           continue
              sigma(psi(j,k))=(t1(psi(j,k))-sum)/dble(n)
140        continue
150     continue
C calculate mu from beta
        do 300 m=1,ncells
           do 240 j=1,q
              sum=0d0
              do 235 i=1,r
                 sum=sum+design(m,i)*beta(i,j)
235           continue
              mu(j,m)=sum
240         continue
300     continue
        return
        end
C***********************************************************************
        subroutine mmnm(d,l,u,p,psi,m)
C Multiplies lower triangular matrix l by upper-triangular matrix u,
C both in packed storage, and puts result into m which is unpacked
        integer d,p,psi(p,p)
        double precision l(d),u(d),sum,m(p,p)
        do 10 i=1,p
	   do 5 j=1,p
	     sum=0
             do 2 k=1,min(i,j)
	        sum=sum+l(psi(i,k))*u(psi(k,j))
2            continue
             m(i,j)=sum
5          continue
10      continue
        return
        end
C***********************************************************************
	subroutine cholsm(d,theta,p,psi,mc,nmc)
	integer d,p,psi(p,p),mc(p),nmc
	double precision theta(d),tmp
	do 40 i=1,nmc
	  tmp=0.0d0
	  do 10 k=1,i-1
	    tmp=tmp+theta(psi(mc(k),mc(i)))**2
10	  continue
	  theta(psi(mc(i),mc(i)))=sqrt(theta(psi(mc(i),mc(i)))-tmp)
	  do 30 j=i+1,nmc
	    tmp=0.d0
	    do 20 k=1,i-1
	      tmp=tmp+theta(psi(mc(k),mc(i)))*theta(psi(mc(k),mc(j)))
20	    continue
	    theta(psi(mc(i),mc(j)))=(theta(psi(mc(i),mc(j)))-tmp)
     /             /theta(psi(mc(i),mc(i)))
30	  continue
40	continue
	end
C***********************************************************************
        subroutine invtrm(npsi,t,q,psi)
C Inverts triangular matrix in packed storage
        integer npsi,q,psi(q,q)
        double precision t(npsi),sum
        t(psi(1,1))=1.d0/t(psi(1,1))
        do 10 k=2,q
           t(psi(k,k))=1.0d0/t(psi(k,k))
           do 5 j=1,k-1
              sum=0
              do 3 i=j,k-1
                 sum=sum+t(psi(i,j))*t(psi(i,k))
3             continue
              t(psi(j,k))=-sum*t(psi(k,k))
5          continue
10      continue
        return
        end
C***********************************************************************
        subroutine invsym(q,psi,npsi,mat,mc)
C inverts symmetric matrix in packed storage
        integer q,psi(q,q),npsi,mc(q)
        double precision mat(npsi),sum
        do 1 j=1,q
           mc(j)=j
1       continue
        call cholsm(npsi,mat,q,psi,mc,q)
        call invtrm(npsi,mat,q,psi)
        do 4 j=1,q
           do 3 k=1,j
              sum=0d0
              do 2 i=j,q
                 sum=sum+mat(psi(i,j))*mat(psi(i,k))
2             continue
              mat(psi(j,k))=sum
3          continue
4       continue
        return
        end
C***********************************************************************
        subroutine bfacm(npsi,b,q,psi,m)
C draws triangular square-root of a Wishart(m,I) using Bartlett
C decomposition, putting result into packed storage
        integer npsi,q,psi(q,q)
        double precision b(npsi),m
	do 10 j=1,q
	   b(psi(j,j))=dble(sqrt(2.*gamm((sngl(m)-float(j)+1.)/2.)))
10      continue
        do 30 j=1,q-1
	   do 20 k=j+1,q
             b(psi(j,k))=dble(gauss())
20         continue
30      continue
        return
        end
C***********************************************************************
        subroutine pstepcm(q,psi,npsi,ncells,t1,t2,t3,sigma,mu,n,
     /       r,design,wk,mcr,psir,npsir,wkr,wkd,mcz,chf,beta,mx)
        integer q,psi(q,q),npsi,ncells,m,r,mcr(r),npsir
        integer psir(r,r),mi,mcz(q),h
        double precision junk
        double precision t1(npsi),t2(q,ncells),t3(ncells),sum
        double precision sigma(npsi),mu(q,ncells)
        double precision design(ncells,r),wk(npsir),wkr(r),wkd(ncells)
        double precision chf(npsi),beta(r,q),zz,mx(q,q),df
        junk=dble(gauss())
C calculate least-squares quantities
        do 30 j=1,r
           do 20 k=j,r
              sum=0d0
              do 10 m=1,ncells
                 sum=sum+design(m,j)*design(m,k)*t3(m)
 10           continue
              wk(psir(j,k))=sum
 20        continue
 30     continue
        call invsym(r,psir,npsir,wk,mcr)
        do 100 h=1,r
           do 60 mi=1,ncells
              sum=0d0
              do 50 j=1,r
                 sum=sum+wk(psir(h,j))*design(mi,j)
 50           continue
              wkd(mi)=sum
 60        continue
           do 90 k=1,q
              sum=0d0
              do 80 mi=1,ncells
                 sum=sum+wkd(mi)*t2(k,mi)
 80           continue
              beta(h,k)=sum
 90        continue
 100    continue
        do 150 j=1,q
           do 110 h=1,r
              sum=0d0
              do 105 m=1,ncells
                 sum=sum+t2(j,m)*design(m,h)
 105          continue
              wkr(h)=sum
 110       continue
           do 140 k=j,q
              sum=0d0
              do 120 h=1,r
                 sum=sum+wkr(h)*beta(h,k)
 120          continue
              t1(psi(j,k))=t1(psi(j,k))-sum
 140       continue
 150    continue
C draw sigma
        do 155 j=1,q
           mcz(j)=j
 155    continue
        df=dble(n-r)
        call cholsm(npsi,t1,q,psi,mcz,q)
        call bfacm(npsi,chf,q,psi,df)
        call invtrm(npsi,chf,q,psi)
        call mmnm(npsi,chf,t1,q,psi,mx)
        do 200 i=1,q
           do 190 j=i,q
              sum=0d0
              do 185 k=1,q
                 sum=sum+mx(k,i)*mx(k,j)
 185          continue
              sigma(psi(i,j))=sum
 190       continue
 200    continue
C draw beta
        do 251 j=1,npsi
           chf(j)=sigma(j)
 251    continue
        call cholsm(npsi,chf,q,psi,mcz,q)
        call cholsm(npsir,wk,r,psir,mcr,r)
        do 300 j=1,q
           do 260 k=1,r
              wkr(k)=0d0
 260       continue
           do 270 k=1,r
              zz=dble(gauss())
              do 265 h=k,r
                 wkr(h)=wkr(h)+zz*wk(psir(h,k))
 265          continue
 270       continue
           do 290 k=j,q
              do 280 h=1,r
                 beta(h,k)=beta(h,k)+chf(psi(j,k))*wkr(h)
 280          continue
 290       continue
 300    continue
C calculate mu from beta
        do 400 m=1,ncells
           do 340 j=1,q
              sum=0d0
              do 335 i=1,r
                 sum=sum+design(m,i)*beta(i,j)
 335          continue
              mu(j,m)=sum
 340       continue
 400    continue
        return
        end
C************************************************************************
        subroutine initc(p,c,mc,nmc)
        integer p,c(p),mc(p),nmc
        do 1 j=1,nmc
           c(mc(j))=1
1       continue
        return
        end
C************************************************************************
        subroutine advc(p,c,d,mc,nmc)
C Advances c to next value
        integer p,c(p),d(p),mc(p),nmc
        do 1 j=1,nmc
           if(c(mc(j)).lt.d(mc(j)))then
              c(mc(j))=c(mc(j))+1
              goto 2
           else
              c(mc(j))=1
           endif
1       continue
2       continue
        return
        end
C************************************************************************
        subroutine gtmmis(p,c,mc,nmc,jmp,mmis)
C Calculates mmis from c
        integer p,c(p),mc(p),nmc,jmp(p),mmis
        mmis=0
        do 1 j=1,nmc
           mmis=mmis+(c(mc(j))-1)*jmp(mc(j))
1       continue
        return
        end
C************************************************************************
        subroutine gtdmis(p,d,mc,nmc,dmis)
        integer p,d(p),mc(p),nmc,dmis
        dmis=1
        do 1 j=1,nmc
           dmis=dmis*d(mc(j))
1       continue
        return
        end
C************************************************************************
        subroutine gtntab(ncon,con,ntab)
C find number of marginal tables to be fit
        integer ncon,con(ncon),ntab,flag,posn
        ntab=0
        flag=0
        do 10 posn=1,ncon
           if((con(posn).ne.0).and.(flag.eq.0))flag=1
           if((con(posn).eq.0).and.(flag.eq.1))then
              ntab=ntab+1
              flag=0
           endif
           if((flag.eq.1).and.(posn.eq.ncon))ntab=ntab+1
10      continue
        return
        end
C************************************************************************
        subroutine gtmarg(ncon,con,posn,p,marg,nmarg)
C extract the next set of margins to be fit, store them in the first
C nmarg elements of marg. For first set, posn should be 0.
        integer ncon,con(ncon),p,marg(p),nmarg,posn
1       continue
           posn=posn+1
           if(con(posn).eq.0)goto 1
        nmarg=0
11      continue
          if(con(posn).eq.0)goto 12
          nmarg=nmarg+1
          marg(nmarg)=con(posn)
          if(posn.eq.ncon)goto 12
          posn=posn+1
          goto 11
12      continue
        return
        end
C************************************************************************
        subroutine gtrest(p,marg,nmarg,rest,nrest)
        integer p,marg(p),nmarg,rest(p),nrest,flag
        nrest=0
        do 20 j=1,p
           flag=0
           do 10 k=1,nmarg
              if(j.eq.marg(k))then
                 flag=1
                 goto 15
              endif
10         continue
15         continue
           if(flag.eq.0)then
              nrest=nrest+1
              rest(nrest)=j
           endif
20      continue
        return
        end
C************************************************************************
        subroutine ipf(ncells,table,fit,ncon,con,p,d,jmp,c,marg,rest,
     /     eps)
C Performs one step of iterative proportional fitting, raking fit
C to the margins of table. Margins to be fit are specified in con.
        integer ncells,ncon,con(ncon),p,d(p),jmp(p),c(p),marg(p)
        integer rest(p),nmarg,nrest,m,mmarg,mrest,posn,ntab,tabno
        integer out,in,dmarg,drest
        double precision sumt,sumf,table(ncells),fit(ncells),eps
        call gtntab(ncon,con,ntab)
        posn=0
        do 100 tabno=1,ntab
           call gtmarg(ncon,con,posn,p,marg,nmarg)
           call gtrest(p,marg,nmarg,rest,nrest)
           call gtdmis(p,d,marg,nmarg,dmarg)
           drest=ncells/dmarg
           call initc(p,c,marg,nmarg)
           mmarg=1
           do 90 out=1,dmarg
              if(out.ne.1)then
                 call advc(p,c,d,marg,nmarg)
                 call gtmmis(p,c,marg,nmarg,jmp,mmarg)
                 mmarg=mmarg+1
              endif
              call sum2c(p,c,rest,nrest,d,jmp,mmarg,drest,ncells,
     /             table,sumt,fit,sumf)
              call initc(p,c,rest,nrest)
              if(sumf.ne.0)then
                 mrest=0
                 do 80 in=1,drest
                    if(in.ne.1)then
                       call advc(p,c,d,rest,nrest)
                       call gtmmis(p,c,rest,nrest,jmp,mrest)
                    endif
                    m=mmarg+mrest
                    if(fit(m).ge.eps) then
                       fit(m)=fit(m)*(sumt/sumf)
                    else
                       fit(m)=0d0
                    endif
80               continue
              endif
90         continue
100     continue
        return
        end
C************************************************************************
        subroutine bipf(ncells,table,theta,prior,ncon,con,p,d,jmp,c,
     /     marg,rest,err)
C Performs one cycle of Bayesian ipf. Cell counts are in table
C and starting value in theta.  Prior hyperparameters are in prior,
C with structural zeros denoted by -999. Replaces theta with an
C updated value.
        integer ncells,ncon,con(ncon),p,d(p),jmp(p),c(p),marg(p)
        integer rest(p),nmarg,nrest,m,mmarg,mrest,posn,ntab,tabno
        integer out,in,dmarg,drest,err,zflag
        double precision sumt,sumf,table(ncells),theta(ncells),g
        double precision prior(ncells),sum3
        call gtntab(ncon,con,ntab)
        err=0
        posn=0
        do 100 tabno=1,ntab
           sum3=0d0
           call gtmarg(ncon,con,posn,p,marg,nmarg)
           call gtrest(p,marg,nmarg,rest,nrest)
           call gtdmis(p,d,marg,nmarg,dmarg)
           drest=ncells/dmarg
           call initc(p,c,marg,nmarg)
           mmarg=1
           do 90 out=1,dmarg
              if(out.ne.1)then
                 call advc(p,c,d,marg,nmarg)
                 call gtmmis(p,c,marg,nmarg,jmp,mmarg)
                 mmarg=mmarg+1
              endif
              zflag=0
              call sum3c(p,c,rest,nrest,d,jmp,mmarg,drest,ncells,
     /             table,sumt,theta,sumf,prior,zflag)
              call initc(p,c,rest,nrest)
              if(sumt.le.0)then
                 err=1
                 goto 150
              endif
              if(zflag.eq.1)then
                 g=dble(gamm(sngl(sumt)))+1.0d-20
                 sum3=sum3+g
              endif
              mrest=0
              do 80 in=1,drest
                 if(in.ne.1)then
                    call advc(p,c,d,rest,nrest)
                    call gtmmis(p,c,rest,nrest,jmp,mrest)
                 endif
                 m=mmarg+mrest
                 theta(m)=theta(m)*g/sumf
 80           continue
 90        continue
           do 95 m=1,ncells
              theta(m)=theta(m)/sum3
 95        continue
100     continue
150     continue
        return
        end
C************************************************************************
        subroutine sum2c(p,c,mc,nmc,d,jmp,mobs,dmis,
     /      ncells,table1,sum1,table2,sum2)
        integer p,c(p),mc(p),nmc,d(p),jmp(p),mobs
        integer m,mmis,dmis,ncells
        double precision sum1,sum2,table1(ncells),table2(ncells)
        call initc(p,c,mc,nmc)
        sum1=0
        sum2=0
        mmis=0
        do 1 i=1,dmis
           if(i.ne.1)then
              call advc(p,c,d,mc,nmc)
              call gtmmis(p,c,mc,nmc,jmp,mmis)
           endif
           m=mobs+mmis
           sum1=sum1+table1(m)
           sum2=sum2+table2(m)
1       continue
        return
        end
C************************************************************************
        subroutine sum3c(p,c,mc,nmc,d,jmp,mobs,dmis,
     /      ncells,table1,sum1,table2,sum2,prior,zflag)
        integer p,c(p),mc(p),nmc,d(p),jmp(p),mobs
        integer m,mmis,dmis,ncells,zflag
        double precision sum1,sum2,table1(ncells),table2(ncells),
     /    prior(ncells)
        call initc(p,c,mc,nmc)
        sum1=0d0
        sum2=0d0
        mmis=0
        do 1 i=1,dmis
           if(i.ne.1)then
              call advc(p,c,d,mc,nmc)
              call gtmmis(p,c,mc,nmc,jmp,mmis)
           endif
           m=mobs+mmis
           sum2=sum2+table2(m)
           if(prior(m).ne.dble(-999)) then
              sum1=sum1+table1(m)+prior(m)
              zflag=1
           endif
1       continue
        return
        end
C************************************************************************
        subroutine addstat1(q,psi,npsi,ncells,sigma,mu,theta,t1,t2,
     /       t3,n,z,i,p,mcw,nmcw,c,d,jmp,dmis,mobs,ocz,nocz,mcz,nmcz)
        integer q,psi(q,q),npsi,ncells,n,i,p,mcw(p),nmcw,d(p)
        integer jmp(p),dmis,mobs,ocz(q),nocz,mcz(q),nmcz,c(p)
        double precision sigma(npsi),mu(q,ncells),theta(ncells)
        double precision t1(npsi),t2(q,ncells),t3(ncells),sum
        double precision z(n,q)
        do 20 j=1,nmcz
           sum=mu(mcz(j),mobs)
           do 10 k=1,nocz
              sum=sum+sigma(psi(mcz(j),ocz(k)))*z(i,ocz(k))
10         continue
           z(i,mcz(j))=sum
           t2(mcz(j),mobs)=t2(mcz(j),mobs)+sum
20      continue
        do 100 j=1,nmcz
           do 70 k=1,nocz
              t1(psi(mcz(j),ocz(k)))=t1(psi(mcz(j),ocz(k)))+
     /             z(i,mcz(j))*z(i,ocz(k))
70         continue
           do 80 k=j,nmcz
              t1(psi(mcz(j),mcz(k)))=t1(psi(mcz(j),mcz(k)))+
     /             z(i,mcz(j))*z(i,mcz(k))+sigma(psi(mcz(j),mcz(k)))
80         continue
100     continue
        return
        end
C***********************************************************************
        subroutine addstat2(q,psi,npsi,ncells,sigma,mu,theta,t1,t2,
     /       t3,n,z,i,p,mcw,nmcw,c,d,jmp,dmis,mobs,ocz,nocz,mcz,nmcz)
        integer q,psi(q,q),npsi,ncells,n,i,p,mcw(p),nmcw,c(p),d(p)
        integer jmp(p),dmis,mobs,ocz(q),nocz,mcz(q),nmcz,mmis,m,a
        double precision sigma(npsi),mu(q,ncells),theta(ncells)
        double precision t1(npsi),t2(q,ncells),t3(ncells),sum
        double precision z(n,q)
        call initc(p,c,mcw,nmcw)
        mmis=0
        do 200 a=1,dmis
           if(a.ne.1)then
              call advc(p,c,d,mcw,nmcw)
              call gtmmis(p,c,mcw,nmcw,jmp,mmis)
           endif
           m=mobs+mmis
           if(theta(m).ne.-999.0d0)then
              t3(m)=t3(m)+theta(m)
              do 20 j=1,nmcz
                 sum=mu(mcz(j),m)
                 do 10 k=1,nocz
                    sum=sum+sigma(psi(mcz(j),ocz(k)))*z(i,ocz(k))
10               continue
                 z(i,mcz(j))=sum
                 t2(mcz(j),m)=t2(mcz(j),m)+sum*theta(m)
20            continue
              do 30 j=1,nocz
                 t2(ocz(j),m)=t2(ocz(j),m)+z(i,ocz(j))*theta(m)
30            continue
              do 100 j=1,nmcz
                 do 70 k=1,nocz
                    t1(psi(mcz(j),ocz(k)))=t1(psi(mcz(j),ocz(k)))+
     /                   theta(m)*z(i,mcz(j))*z(i,ocz(k))
70               continue
                 do 80 k=j,nmcz
                    t1(psi(mcz(j),mcz(k)))=t1(psi(mcz(j),mcz(k)))+
     /                   theta(m)*z(i,mcz(j))*z(i,mcz(k))+
     /                   theta(m)*sigma(psi(mcz(j),mcz(k)))
80               continue
100           continue
           endif
200     continue
        return
        end
C***********************************************************************
        subroutine gtprob(q,ncells,mu,pii,n,z,i,p,mcw,nmcw,c,d,jmp,
     /       dmis,mobs,ocz,nocz,theta)
C Calculates cell probs corresponding to observation i in z, storing
C the result in theta. For structural zeros, theta is set to -999.0.
        integer i,q,ncells,n,p,mcw(p),nmcw,d(p),jmp(p),dmis,mobs
        integer mmis,m,a,ocz(q),nocz,c(p)
        double precision mu(q,ncells),pii(ncells),sum,theta(ncells)
        double precision z(n,q)
        call initc(p,c,mcw,nmcw)
        sum=0d0
        mmis=0
        do 30 a=1,dmis
           if(a.ne.1)then
              call advc(p,c,d,mcw,nmcw)
              call gtmmis(p,c,mcw,nmcw,jmp,mmis)
           endif
           m=mobs+mmis
           theta(m)=pii(m)
           if(theta(m).ne.-999.0d0)then
              theta(m)=0.5d0*theta(m)
              do 20 j=1,nocz
                 theta(m)=theta(m)+mu(ocz(j),m)*z(i,ocz(j))
20            continue
              theta(m)=exp(theta(m))
              sum=sum+theta(m)
           endif
30      continue
        call initc(p,c,mcw,nmcw)
        mmis=0
        do 40 a=1,dmis
           if(a.ne.1)then
              call advc(p,c,d,mcw,nmcw)
              call gtmmis(p,c,mcw,nmcw,jmp,mmis)
           endif
           m=mobs+mmis
           if(theta(m).ne.-999.0d0) theta(m)=theta(m)/sum
40      continue
        return
        end
C************************************************************************
        subroutine lobsm(q,psi,npsi,ncells,sigma,mu,pii,npattz,rz,mcz,
     /       ocz,mdpzgrp,npattw,p,rw,mcw,mdpwgrp,ngrp,mobs,mobsst,nmobs,
     /       n,z,d,jmp,c,ll)
        integer q,psi(q,q),npsi,ncells,npattz,pattz,rz(npattz,q)
        integer mcz(q),ocz(q),nmcz,nocz,pattw,grpno,a,b,n,ngrp
        integer mdpzgrp(npattz),npattw,p,rw(npattw,p),mcw(p)
        integer nmcw,mdpwgrp(npattw),dmis,c(p),d(p),jmp(p)
        integer mobs(ngrp),mobsst(ngrp),nmobs(ngrp)
        double precision z(n,q)
        double precision sigma(npsi),mu(q,ncells),pii(ncells)
        double precision ll,l2,l3,logdet
        l2=0d0
        l3=0d0
        logdet=0d0
        pattw=0
        grpno=0
        do 200 pattz=1,npattz
           call swpobsm(q,psi,npsi,ncells,sigma,mu,pii,npattz,rz,
     /          pattz,logdet,0)
           call gtmc(q,npattz,rz,pattz,mcz,nmcz,q)
           call gtoc(q,npattz,rz,pattz,ocz,nocz,q)
           do 180 a=1,mdpzgrp(pattz)
              pattw=pattw+1
              call gtmc(p,npattw,rw,pattw,mcw,nmcw,p)
              call gtdmis(p,d,mcw,nmcw,dmis)
              do 170 b=1,mdpwgrp(pattw)
                 grpno=grpno+1
                 do 150 i=mobsst(grpno),(mobsst(grpno)+nmobs(grpno)-1)
                    call qdfrm(q,psi,npsi,ncells,sigma,mu,pii,
     /                   n,z,i,p,mcw,nmcw,c,d,jmp,dmis,
     /                   mobs(grpno),ocz,nocz,mcz,nmcz,l3)
150              continue
                 l2=l2+logdet*dble(nmobs(grpno))
170           continue
180        continue
200     continue
        ll=-(l2/2d0)+l3
        return
        end
C************************************************************************
        subroutine qdfrm(q,psi,npsi,ncells,sigma,mu,pii,
     /       n,z,i,p,mcw,nmcw,c,d,jmp,dmis,mobs,ocz,nocz,mcz,nmcz,l3)
        integer q,psi(q,q),npsi,ncells,n,i,p,mcw(p),nmcw,d(p)
        integer jmp(p),dmis,mobs,ocz(q),nocz,mcz(q),nmcz,mmis,m,a,c(p)
        double precision sigma(npsi),mu(q,ncells),pii(ncells)
        double precision sum,l3,qf,lsum
        double precision z(n,q)
        lsum=0d0
        call initc(p,c,mcw,nmcw)
        mmis=0
        do 200 a=1,dmis
           if(a.ne.1)then
              call advc(p,c,d,mcw,nmcw)
              call gtmmis(p,c,mcw,nmcw,jmp,mmis)
           endif
           m=mobs+mmis
           if(pii(m).gt.0)then
              qf=0d0
              do 20 j=1,nocz
                 sum=0d0
                 do 10 k=1,nocz
                    sum=sum+sigma(psi(ocz(j),ocz(k)))*
     /                   (z(i,ocz(k))-mu(ocz(k),m))
10               continue
                 qf=qf+sum*(z(i,ocz(j))-mu(ocz(j),m))
20            continue
              lsum=lsum+pii(m)*exp(qf/2d0)
           endif
200     continue
        l3=l3+dlog(lsum)
        return
        end
C***********************************************************************
        subroutine sigexm(d,theta,extr,p,psi,mc,nmc)
C Extracts submatrix of theta corresponding to the columns of mc
        integer d,p,psi(p,p),mc(p),nmc
        double precision theta(d),extr(d)
        do 2 j=1,nmc
           do 1 k=j,nmc
              extr(psi(mc(j),mc(k)))=theta(psi(mc(j),mc(k)))
1          continue
2       continue
        return
        end
C************************************************************************
        real function rangen(init)
        integer a,p,ix,b15,b16,xhi,xalo,leftflo,fhi,k,init
        save ix
        data a/16807/,b15/32768/,b16/65536/,p/2147483647/
        if(init.ne.0) ix=init
        if(ix.eq.0) call rexit('rngseed has not been called')
        xhi=ix/b16
        xalo=(ix-xhi*b16)*a
        leftflo=xalo/b16
        fhi=xhi*a+leftflo
        k=fhi/b15
        ix=(((xalo-leftflo*b16)-p)+(fhi-k*b15)*b16)+k
        if (ix.lt.0)ix=ix+p
        rangen=float(ix)*4.656612875E-10
        return
        end
C***********************************************************************
        subroutine rngs(seed)
C initializes rangen with seed
        integer seed
        tmp=rangen(seed)
        return
        end
C***********************************************************************
        real function chisq(df)
C Generates a chisquare deviate with df degrees of freedom
        real df
        chisq=2.*gamm(df/2.)
        return
        end
C***********************************************************************
        real function gamm(a)
C Generates a random gamma(a) deviate. If a>=1, uses the method of
C Fishman (1976); if 0<a<1, the method of Ahrens (1974)
        real a,u,y,q,e,b,p,u1
        data e/2.718282/
        if(a.ge.1)then
1          continue
           u=rangen(0)
           y=-log(rangen(0))
           q=(y/exp(y-1))**(a-1)
           if(u.le.q)then
              gamm=a*y
           else
              goto 1
           endif
        else
2          continue
           u=rangen(0)
           b=(e+a)/e
           p=b*u
           if(p.gt.1) goto 4
           x=p**(1/a)
           u1=rangen(0)
           if(u1.gt.(e**(-x)))then
              goto 2
           else
              gamm=x
              goto 10
           endif
4          continue
           x=-log((b-p)/a)
           u1=rangen(0)
           if(u1.gt.(x**(a-1)))then
              goto 2
           else
              gamm=x
              goto 10
           endif
        endif
10      continue
        return
        end
C***********************************************************************
        real function gauss()
	integer alt
	real next
        save alt,next
	data pi/3.141593/
        if((alt.ne.0).and.(alt.ne.1)) alt=0
	if(alt.eq.0)then
	  u1=rangen(0)
	  u2=rangen(0)
	  gauss=sqrt(-2*log(u1))*cos(2*pi*u2)
	  next=sqrt(-2*log(u1))*sin(2*pi*u2)
	  alt=1
	else
	  gauss=next
	  alt=0
	endif
	return
	end
C***********************************************************************
