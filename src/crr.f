c Copyright (C) 2000 Robert Gray
c distributed under the terms of the GNU public license
      subroutine covt(ncov,x,n,ncov2,x2,tf,ndf,b,wk,xbt)
      double precision x(n,ncov),x2(n,ncov2),tf(ndf,ncov2),wk
      double precision b(ncov+ncov2),xbt(ncov+ncov2)
      integer ncov,ncov2,ndf,i,n
      wk=0
      if (ncov.gt.0) then
         do 10 i=1,ncov
            xbt(i)=x(1,i)
            wk=wk+x(1,i)*b(i)
 10      continue 
      endif
      if (ncov2.gt.0) then
         do 11 i=1,ncov2
            xbt(i+ncov)=x2(1,i)*tf(1,i)
            wk=wk+xbt(i+ncov)*b(ncov+i)
 11      continue 
      endif
      return
      end

      subroutine crrfsv(t2,ici,n,x,ncov,np,x2,ncov2,tf,ndf,wt,ncg,icg,b,
     $     lik,s,v,xb,xbt,vt) 
c all data sorted with t2 in ascending order
c t2, event/censoring time for each subject
c ici, ici(i)=1 if t2(i)  is a type 1 failure time (& =2 for other failures)
c            & =0 for censored
c x(n,ncov) ph covariates, x2(n,ncov2) covs multiplied by functions of time
c tf(i,j) is the value of the time function which multiplies the jth col of
c      x2, at the ith distinct type 1 failure time (ascending order); 
c      ndf is the number of distinct type 1 failures=# rows in tf
c wt are km est of censoring dists at times in t2- for each group formed by 
c      distinct values of icg.
c censoring groups in icg must be coded 1,2,...,ncg
c output in lik,s,v
c xb,xbt,vt are temporary storage
      double precision t2(n),b(np),s(np),v(np,np),lik,tf(ndf,ncov2)
      double precision xb(np),xb1,vt(np,np),wt(ncg,n),cft,twf
      double precision x(n,ncov),x2(n,ncov2),wk,xb1o,twt,xbt(np)
      integer ici(n),icg(n),n,np,i,j,iuc,ncov,ncov2,ncg,ndf,ldf,k,itmp
      lik=0.d0
      do 1 i=1,np
         s(i)=0
         do 2 j=i,np
            v(i,j)=0
 2       continue 
 1    continue 
c
c to see how this works, use the simple nested loop approach
c
      iuc=n
      ldf=ndf+1
c find next failure time
 98   itmp=iuc
      do 10 i=iuc,1,-1
         itmp=i
         if (ici(i).eq.1) then
            cft=t2(i)
            go to 11
         endif
 10   continue
c no more failures
      return
 11   iuc=itmp
      ldf=ldf-1
      twf=0
      do 13 i=iuc,1,-1
         if (t2(i).lt.cft) go to 14
         itmp=i
         if (ici(i).eq.1) then
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twf=twf+1
c using a minimization routine, so neg of objective fcn, scores, etc
            lik=lik-wk
            do 19 j=1,np
               s(j)=s(j)-xbt(j)
 19         continue 
         endif
 13   continue 
c calculate sums over risk set
 14   iuc=itmp
      xb1=0
      xb1o=xb1
      do 3 i=1,np
         xb(i)=0
         do 4 j=i,np
            vt(i,j)=0
 4       continue 
 3    continue 
      do 15 i=1,n
         if (t2(i).lt.cft) then
            if (ici(i).le.1) go to 15
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twt=exp(wk)*wt(icg(i),iuc)/wt(icg(i),i)
         else
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twt=exp(wk)
         endif
         xb1=xb1+twt
         do 17 j=1,np
            xb(j)=xb(j)+twt*xbt(j)
            xbt(j)=xbt(j)-xb(j)/xb1
 17      continue 
         if (xb1o.gt.0) then
            twt=xb1*twt/xb1o
            do 51 k=1,np
               do 52 j=k,np
                  vt(k,j)=vt(k,j)+twt*xbt(k)*xbt(j)
 52            continue 
 51         continue 
         endif
c     all dsyr('U',np,xb1*twt/xb1o,xbt,1,vt,np)
         xb1o=xb1
 15   continue 
      lik=lik+twf*log(xb1)
      twt=twf/xb1
      do 61 i=1,np
         s(i)=s(i)+twt*xb(i)
         do 62 j=i,np
            v(i,j)=v(i,j)+twt*vt(i,j)
            v(j,i)=v(i,j)
 62      continue 
 61   continue 
c      call daxpy(np,-twf/xb1,xb,1,s,1)
c      call daxpy(np*np,twf/xb1,vt,1,v,1)
      iuc=iuc-1
      if (iuc.gt.0) go to 98
      return
      end

      subroutine crrf(t2,ici,n,x,ncov,np,x2,ncov2,tf,ndf,wt,ncg,icg,b,
     $     lik,xbt) 
c all data sorted with t2 in ascending order
c t2, stop time for each subject
c ici, ici(i)=1 if t2(i)  is a type 1 failure time (& =2 for other failures)
c x(nrx,np) covariates
c wt are km est of censoring dists at times in t2-
c output in lik,s,v
c icrs,wk,xb,xbt,vt are temporary storage
      double precision t2(n),b(np),lik,tf(ndf,ncov2)
      double precision xb1,wt(ncg,n),cft,twf,xbt(np)
      double precision x(n,ncov),x2(n,ncov2),wk,twt
      integer ici(n),icg(n),n,np,i,iuc,ncov,ncov2,ncg,ndf,ldf,itmp
      lik=0.d0
      iuc=n
      ldf=ndf+1
c find next failure time
 98   itmp=iuc
      do 10 i=iuc,1,-1
         itmp=i
         if (ici(i).eq.1) then
            cft=t2(i)
            go to 11
         endif
 10   continue
c no more failures
      return
 11   iuc=itmp
      ldf=ldf-1
      twf=0
      do 13 i=iuc,1,-1
         if (t2(i).lt.cft) go to 14
         itmp=i
         if (ici(i).eq.1) then
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twf=twf+1
            lik=lik-wk
         endif
 13   continue 
c calculate sums over risk set
 14   iuc=itmp
      xb1=0
      do 15 i=1,n
         if (t2(i).lt.cft) then
            if (ici(i).le.1) go to 15
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twt=exp(wk)*wt(icg(i),iuc)/wt(icg(i),i)
         else
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twt=exp(wk)
         endif
         xb1=xb1+twt
 15   continue 
      lik=lik+twf*log(xb1)
      iuc=iuc-1
      if (iuc.gt.0) go to 98
      return
      end

      subroutine crrvv(t2,ici,n,x,ncov,np,x2,ncov2,tf,ndf,wt,ncg,icg,b,
     $     v,v2,vt,xb,xbt,qu,st,ss,ss2,icrsk) 
c all data sorted with t2 in ascending order
c t2, event/censoring time for each subject
c ici, ici(i)=1 if t2(i)  is a type 1 failure time (& =2 for other failures)
c            & =0 for censored
c x(n,ncov) ph covariates, x2(n,ncov2) covs multiplied by functions of time
c tf(i,j) is the value of the time function which multiplies the jth col of
c      x2, at the ith distinct type 1 failure time (ascending order); 
c      ndf is the number of distinct type 1 failures=# rows in tf
c wt are km est of censoring dists at times in t2- for each group formed by 
c      distinct values of icg.
c censoring groups in icg must be coded 1,2,...,ncg
c output in v,v2
c xb,xbt,vt,ss,st,qu,icrsk are temporary storage
      double precision t2(n),b(np),v(np,np),tf(ndf,ncov2)
      double precision xb(n,0:np),vt(np,np),wt(ncg,n),cft
      double precision x(n,ncov),x2(n,ncov2),wk,twt,xbt(np),ss2(np,ncg)
      double precision st(np,2),v2(np,np),ss(np,1),ss0,cft2,qu(np)
      integer ici(n),icg(n),n,np,i,j,ncov,ncov2,ncg,ndf,ldf,k
      integer lc,icrsk(ncg),j1,j2,ldf2
      do 3 i=1,ncg
         icrsk(i)=0
 3    continue 
      do 1 i=1,np
         do 123 j=1,ncg
            ss2(i,j)=0
 123     continue 
         do 2 j=i,np
            v(i,j)=0
            v2(i,j)=0
 2       continue 
 1    continue
      do 5 i=1,n
         icrsk(icg(i))=icrsk(icg(i))+1
         do 4 j=0,np
            xb(i,j)=0
 4       continue 
 5    continue 
      ldf=0
      cft=min(-1.d0,t2(1)*(1-1.d-5))
      do 6 i=1,n
         if (ici(i).ne.1) go to 6
         if (t2(i).gt.cft) then
            cft=t2(i)
            ldf=ldf+1
         endif
         do 7 j=1,n
            if (t2(j).lt.t2(i)) then
               if (ici(j).le.1) go to 7
               call covt(ncov,x(j,1),n,ncov2,x2(j,1),tf(ldf,1),ndf,b,wk,
     $              xbt)
               twt=exp(wk)*wt(icg(j),i)/wt(icg(j),j)
            else
               call covt(ncov,x(j,1),n,ncov2,x2(j,1),tf(ldf,1),ndf,b,wk,
     $              xbt)
               twt=exp(wk)
            endif
            xb(i,0)=xb(i,0)+twt
            do 8 k=1,np
               xb(i,k)=xb(i,k)+twt*xbt(k)
 8          continue 
 7       continue 
 6    continue 
c
      lc=1
      ldf2=0
      cft2=min(-1.d0,t2(1)*(1-1.d-5))
      do 10 i=1,n
         do 11 k=1,np
            st(k,1)=0
 11      continue 
         ldf=0
         cft=min(-1.d0,t2(1)*(1-1.d-5))
         do 15 j=1,n
            if (ici(j).ne.1) go to 15
            if (t2(j).gt.cft) then
               cft=t2(j)
               ldf=ldf+1
            endif
            if (t2(j).le.t2(i)) then
               call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $              xbt)
               twt=exp(wk)
            else if (t2(i).lt.t2(j).and.ici(i).gt.1) then
               call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $              xbt)
               twt=exp(wk)*wt(icg(i),j)/wt(icg(i),i)
            else
               go to 15
            endif
c d lambda hat portion of eta_i
            do 16 k=1,np
               st(k,1)=st(k,1)-(xbt(k)-xb(j,k)/xb(j,0))*twt/xb(j,0)
 16         continue 
 15      continue 
         if (ici(i).eq.1) then
            if (t2(i).gt.cft2) then
               cft2=t2(i)
               ldf2=ldf2+1
            endif
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf2,1),ndf,b,wk,
     $           xbt)
c d N_i portion of eta_i
            do 17 k=1,np
               st(k,1)=st(k,1)+xbt(k)-xb(i,k)/xb(i,0)
 17         continue 
c second derivatives
            do 201 j1=1,np
               do 202 j2=j1,np
                  vt(j1,j2)=0
 202           continue 
 201        continue 
            do 19 j=1,n
               if (t2(j).lt.t2(i)) then
                  if (ici(j).le.1) go to 19
                  call covt(ncov,x(j,1),n,ncov2,x2(j,1),tf(ldf2,1),ndf,
     $              b,wk,xbt)
                  twt=exp(wk)*wt(icg(j),i)/wt(icg(j),j)
               else
                  call covt(ncov,x(j,1),n,ncov2,x2(j,1),tf(ldf2,1),ndf,
     $              b,wk,xbt)
                  twt=exp(wk)
               endif
               do 18 k=1,np
                  xbt(k)=xbt(k)-xb(i,k)/xb(i,0)
 18            continue 
               do 51 k=1,np
                  do 52 j2=k,np
                     vt(k,j2)=vt(k,j2)+twt*xbt(k)*xbt(j2)
 52               continue 
 51            continue 
 19         continue 
            do 251 j1=1,np
               do 252 j2=j1,np
                  v(j1,j2)=v(j1,j2)+vt(j1,j2)/xb(i,0)
 252           continue 
 251        continue 
         endif
c         call dblepr('st1',3,st(1,1),np)
         if (i.eq.1.or.t2(i).gt.t2(i-1)) then
            do 40 j=i,n
               if (t2(j).gt.t2(i)) go to 39
               if (ici(j).eq.0) go to 38
 40         continue 
c calculate qu (q(u))
 38         ldf=ldf2
            cft=cft2
            do 281 k=1,np
               qu(k)=0
 281        continue 
c j1 indexes inner integral in q (t2(i) is lower limit of integration)
c dN_j2 part is 0, because integrand includes I(s>t2(j2)), where s is var
c of integration, so the sum over j1 is calculating the d lambda_1 hat part
            do 41 j1=lc,n
               if (ici(j1).ne.1) go to 41
               if (t2(j1).gt.cft) then
                  cft=t2(j1)
                  ldf=ldf+1
               endif
               ss0=0
               do 255 k=1,np
                  ss(k,1)=0
 255           continue 
c j2 indexes outer sum in q.  because of I(s>t2(j2)) and 
c  I(t2(j2)>=s)+I(t2(j2)<s,eps_j2=2), only type 2 failures contribute
               do 541 j2=1,n
                  if (t2(j2).ge.t2(i)) go to 542
                  if (ici(j2).le.1) go to 541
                  call covt(ncov,x(j2,1),n,ncov2,x2(j2,1),
     $                 tf(ldf,1),ndf,b,wk,xbt)
                  twt=exp(wk)*wt(icg(j2),j1)/wt(icg(j2),j2)
                  ss0=ss0+twt
                  do 256 k=1,np
                     ss(k,1)=ss(k,1)+xbt(k)*twt
 256              continue 
 541           continue 
c qu is q(t2(i))
 542           do 42 k=1,np
                  qu(k)=qu(k)+(ss(k,1)-xb(j1,k)*ss0/xb(j1,0))/
     $                 xb(j1,0)
 42            continue
 41         continue
c            call dblepr('qu',2,qu(1),1)
c update ss2 (note the sign, this is already negative)--ss2 is the 
c d lambda c hat portion of psi_i.  This is the addition
c to the d lambda c hat portion at the ith censoring time.
            do 43 j=i,n
               if (t2(j).gt.t2(i)) go to 39
               if (ici(j).eq.0) then
                  do 282 k=1,np
                     ss2(k,icg(j))=ss2(k,icg(j))-qu(k)/icrsk(icg(j))**2
 282              continue 
               endif
 43         continue 
         endif
c st(,2) is psi
 39      do 243 k=1,np
            st(k,2)=ss2(k,icg(i))
 243     continue 
c dNc portion of psi
         if (ici(i).eq.0) then
            do 45 k=1,np
               st(k,2)=st(k,2)+qu(k)/icrsk(icg(i))
 45         continue 
         endif
         do 271 j1=np,1,-1
            st(j1,1)=st(j1,1)+st(j1,2)
            do 272 j2=j1,np
               v2(j1,j2)=v2(j1,j2)+st(j1,1)*st(j2,1)
 272        continue 
 271     continue 
c         call dblepr('st1',3,st(1,1),np)
         if (i.lt.n) then
            if (t2(i+1).gt.t2(i)) then
               do 249 j=lc,i
                  icrsk(icg(j))=icrsk(icg(j))-1
 249           continue 
               lc=i+1
            endif
         endif
 10   continue
c      call intpr('icrsk',5,icrsk,2)
      do 333 j1=1,np-1
         do 334 j2=j1+1,np
            v(j2,j1)=v(j1,j2)
            v2(j2,j1)=v2(j1,j2)
 334     continue 
 333  continue 
      return
      end

      subroutine crrsr(t2,ici,n,x,ncov,np,x2,ncov2,tf,ndf,wt,ncg,icg,b,
     $     res,xb,xbt)
c compute contribution to scores at each type 1 failure time
c all data sorted with t2 in ascending order
c t2, event/censoring time for each subject
c ici, ici(i)=1 if t2(i)  is a type 1 failure time (& =2 for other failures)
c            & =0 for censored
c x(n,ncov) ph covariates, x2(n,ncov2) covs multiplied by functions of time
c tf(i,j) is the value of the time function which multiplies the jth col of
c      x2, at the ith distinct type 1 failure time (ascending order); 
c      ndf is the number of distinct type 1 failures=# rows in tf
c wt are km est of censoring dists at times in t2- for each group formed by 
c      distinct values of icg.
c score residuals in res
c xb,xbt are temporary storage
      double precision t2(n),b(np),tf(ndf,ncov2)
      double precision xb(np),wt(ncg,n),cft,res(np,ndf)
      double precision x(n,ncov),x2(n,ncov2),wk,twt,twf,xbt(np),xb1
      integer ici(n),icg(n),n,np,i,j,ncov,ncov2,ncg,ndf,ldf,iuc,itmp
      do 1 i=1,ndf
         do 2 j=1,np
            res(j,i)=0
 2       continue 
 1    continue 
      iuc=n
      ldf=ndf+1
c find next failure time
 98   itmp=iuc
      do 10 i=iuc,1,-1
         itmp=i
         if (ici(i).eq.1) then
            cft=t2(i)
            go to 11
         endif
 10   continue
c no more failures
      return
 11   iuc=itmp
      ldf=ldf-1
      twf=0
      do 13 i=iuc,1,-1
         if (t2(i).lt.cft) go to 14
         itmp=i
         if (ici(i).eq.1) then
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twf=twf+1
            do 241 j=1,np
               res(j,ldf)=res(j,ldf)+xbt(j)
 241        continue 
         endif
 13   continue 
c calculate sums over risk set
 14   iuc=itmp
      xb1=0
      do 3 i=1,np
         xb(i)=0
 3    continue 
      do 15 i=1,n
         if (t2(i).lt.cft) then
            if (ici(i).le.1) go to 15
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twt=exp(wk)*wt(icg(i),iuc)/wt(icg(i),i)
         else
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twt=exp(wk)
         endif
         xb1=xb1+twt
         do 17 j=1,np
            xb(j)=xb(j)+twt*xbt(j)
 17      continue 
 15   continue 
      twt=-twf/xb1
      do 61 i=1,np
         res(i,ldf)=res(i,ldf)+twt*xb(i)
 61   continue 
      iuc=iuc-1
c      call intpr('iuc',3,iuc,1)
      if (iuc.gt.0) go to 98
      return
      end

      subroutine crrfit(t2,ici,n,x,ncov,np,x2,ncov2,tf,ndf,wt,ncg,icg,b,
     $     res,xbt)
c computes jumps in estimate of cumulative underlying hazard
c at each type 1 failure time
c all data sorted with t2 in ascending order
c t2, event/censoring time for each subject
c ici, ici(i)=1 if t2(i)  is a type 1 failure time (& =2 for other failures)
c            & =0 for censored
c x(n,ncov) ph covariates, x2(n,ncov2) covs multiplied by functions of time
c tf(i,j) is the value of the time function which multiplies the jth col of
c      x2, at the ith distinct type 1 failure time (ascending order); 
c      ndf is the number of distinct type 1 failures=# rows in tf
c wt are km est of censoring dists at times in t2- for each group formed by 
c      distinct values of icg.
c res(i) is the jump in the estimated underlying cum pseudo haz for type 1
c      failures at the ith distinct (ascending) type 1 failure time
c xbt is temporary storage
      double precision t2(n),b(np),tf(ndf,ncov2)
      double precision wt(ncg,n),cft,res(ndf)
      double precision x(n,ncov),x2(n,ncov2),wk,twt,twf,xbt(np),xb1
      integer ici(n),icg(n),n,np,i,ncov,ncov2,ncg,ndf,ldf,iuc,itmp
      do 1 i=1,ndf
         res(i)=0
 1    continue 
      iuc=1
      ldf=0
c find next failure time
 98   itmp=iuc
      do 10 i=iuc,n
         itmp=i
         if (ici(i).eq.1) then
            cft=t2(i)
            go to 11
         endif
 10   continue
c no more failures
      return
 11   iuc=itmp
      ldf=ldf+1
      twf=0
      do 13 i=iuc,n
         if (t2(i).gt.cft) go to 14
         itmp=i
         if (ici(i).eq.1) twf=twf+1
 13   continue 
c calculate sums over risk set
 14   iuc=itmp
      xb1=0
      do 15 i=1,n
         if (t2(i).lt.cft) then
            if (ici(i).le.1) go to 15
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twt=exp(wk)*wt(icg(i),iuc)/wt(icg(i),i)
         else
            call covt(ncov,x(i,1),n,ncov2,x2(i,1),tf(ldf,1),ndf,b,wk,
     $           xbt)
            twt=exp(wk)
         endif
         xb1=xb1+twt
 15   continue 
      res(ldf)=res(ldf)+twf/xb1
      iuc=iuc+1
      if (iuc.le.n) go to 98
      return
      end
