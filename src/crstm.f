c Copyright (C) 2000 Robert Gray
c distributed under the terms of the GNU public license
      subroutine crstm(y,m,ig,ist,no,rho,nst,ng,s,vs,ys,
     &ms,igs,v,st,vt,wk,iwk)
c
c subroutine to calculate score and variance matrix for comparing
c cumulative incidence curves for a specific cause among groups.
c  test statistic given by s' inv(vs) s, dist approx chi-square(ng-1)
c
c  everything starting with i-n is integer, all others double precision
c
c  On input:
c    y is the failure times (sorted in increasing order)
c    m is coded 0 if censored, 1 if failed from the cause of interest
c       2 if failed from some other cause.
c    ig denotes group membership, must be coded 1,2,...,ng (ie
c       consecutive integers from 1 to ng), where ng is the number of groups
c    ist denotes strata membership, must be coded 1,2,...,nst, where nst
c       is the number of strata (code all 1's if only 1 strata)
c    no is the # of observations (length of y,m,ig,ist)
c    rho is the power used in the weight function in the test statistic
c    nst and ng are the # strata and # groups
c  Length of ys, ms, and igs must be at least as long as the size of the
c    largest strata
c  Length of s and st must be at least ng-1
c  Length of v and vt must be at least ng*(ng-1)/2
c  Length of vs must be at least (ng-1)^2
c  Length of wk must be at least ng*(4+3*ng)
c  Length of iwk must be at least 4*ng
c  
c  On output:
c    s gives the scores for the first ng-1 groups, and
c    vs the estimated variance covariance matrix of these scores
c
      implicit double precision (a-h,o-z)
      dimension y(no),m(no),ig(no),ist(no),ys(1),ms(1)
      dimension wk(ng*(4+3*ng)),iwk(4*ng)
      dimension igs(1),s(ng-1),v(1),st(ng-1),vt(1),vs(ng-1,ng-1)
      ng1=ng-1
      ng2=ng*ng1/2
      l=0
      do 12 i=1,ng1
      s(i)=0
      do 14 j=1,i
      l=l+1
      v(l)=0
  14  continue
  12  continue
      do 20 ks=1,nst
      n=0
      do 21 i=1,no
      if (ist(i).ne.ks) go to 21
      n=n+1
      ys(n)=y(i)
      ms(n)=m(i)
      igs(n)=ig(i)
  21  continue
      ng3=4*ng+1
      ng4=ng*ng
      call crst(ys(1),ms(1),igs(1),n,ng,rho,st,vt,ng1,ng2,wk(1),
     &wk(ng+1),wk(2*ng+1),wk(3*ng+1),wk(ng3),wk(ng3+ng4),wk(ng3+2*ng4),
     &wk(ng3+2*ng4+ng),iwk(1),iwk(ng+1))
      l=0
      do 23 i=1,ng1
      s(i)=s(i)+st(i)
      do 24 j=1,i
      l=l+1
      v(l)=v(l)+vt(l)
  24  continue
  23  continue
  20  continue
      l=0
      do 31 i=1,ng1
      do 332 j=1,i
      l=l+1
      vs(i,j)=v(l)
      vs(j,i)=vs(i,j)
 332  continue
  31  continue
 996  return
      end

      subroutine crst(y,m,ig,n,ng,rho,s,v,ng1,nv,f1m,f1,skmm,
     &skm,c,a,v3,v2,rs,d)
      implicit double precision (a-h,o-z)
      dimension y(n),s(ng1),f1m(ng),f1(ng),skmm(ng),skm(ng)
      dimension c(ng,ng),a(ng,ng),v(nv),v3(ng)
      dimension v2(ng1,ng)
      integer m(n),ig(n),rs(ng),d(0:2,ng)
c rs(j) will be the risk set size in group j at the current failure time
c  (initially the sample size in each group)
      do 15 i=1,ng
  15  rs(i)=0
      do 16 i=1,n
      j=ig(i)
  16  rs(j)=rs(j)+1
      l=0
      do 11 i=1,ng1
      s(i)=0
      do 12 j=1,i
      l=l+1
  12  v(l)=0
  11  continue
      do 14 i=1,ng
      f1m(i)=0
      f1(i)=0
      skmm(i)=1
      skm(i)=1
      v3(i)=0
      do 9 j=1,ng1
   9  v2(j,i)=0
      do 8 j=1,ng
   8  c(i,j)=0
  14  continue
      fm=0
      f=0
c begin looping over unique times:
      ll=1
      lu=ll
  50  lu=lu+1
      if (lu.gt.n) go to 55
      if (y(lu).gt.y(ll)) go to 55
      go to 50
  55  lu=lu-1
      nd1=0
      nd2=0
c  d will contain the # in each group censored, failed from
c  cause 1, and failing from cause 2, at this time
      do 56 i=1,ng
      d(0,i)=0
      d(1,i)=0
  56  d(2,i)=0
      do 57 i=ll,lu
      j=ig(i)
      k=m(i)
  57  d(k,j)=d(k,j)+1
      do 58 i=1,ng
      nd1=nd1+d(1,i)
      nd2=nd2+d(2,i)
  58  continue
      if (nd1.eq.0.and.nd2.eq.0) go to 90
      tr=0
      tq=0
      do 60 i=1,ng
      if (rs(i).le.0) go to 60
      td=d(1,i)+d(2,i)
c skmm is left continuous, and skm right continuous, km est.
      skm(i)=skmm(i)*(rs(i)-td)/rs(i)
c f1m is left continuous, and f1 right continuous, cuminc est.
      f1(i)=f1m(i)+(skmm(i)*d(1,i))/rs(i)
c in notation of the paper, tr is \sum_r\hat{h}_r, and tq is \sum_r R_r
      tr=tr+rs(i)/skmm(i)
      tq=tq+rs(i)*(1-f1m(i))/skmm(i)
  60  continue
      f=fm+nd1/tr
      fb=(1-fm)**rho
      do 66 i=1,ng
      do 166 j=i,ng
  166 a(i,j)=0
      if (rs(i).le.0) go to 66
      t1=rs(i)/skmm(i)
      a(i,i)=fb*t1*(1-t1/tr)
      c(i,i)=c(i,i)+a(i,i)*nd1/(tr*(1-fm))
      k=i+1
      if (k.gt.ng) go to 66
      do 67 j=k,ng
      if (rs(j).le.0) go to 67
      a(i,j)=-fb*t1*rs(j)/(skmm(j)*tr)
      c(i,j)=c(i,j)+a(i,j)*nd1/(tr*(1-fm))
  67  continue
  66  continue
      do 68 i=2,ng
      k=i-1
      do 69 j=1,k
      a(i,j)=a(j,i)
  69  c(i,j)=c(j,i)
  68  continue
      do 74 i=1,ng1
      if (rs(i).le.0) go to 74
      s(i)=s(i)+fb*(d(1,i)-nd1*rs(i)*(1-f1m(i))/(skmm(i)*tq))
  74  continue
      if (nd1.le.0) go to 77
      do 72 k=1,ng
      if (rs(k).le.0) go to 72
      t4=1
      if (skm(k).gt.0) t4=1-(1-f)/skm(k)
      t5=1
      if (nd1.gt.1) t5=1-(nd1-1)/(tr*skmm(k)-1)
      t3=t5*skmm(k)*nd1/(tr*rs(k))
      v3(k)=v3(k)+t4*t4*t3
      do 70 i=1,ng1
      t1=a(i,k)-t4*c(i,k)
      v2(i,k)=v2(i,k)+t1*t4*t3
      do 71 j=1,i
      l=i*(i-1)/2+j
      t2=a(j,k)-t4*c(j,k)
      v(l)=v(l)+t1*t2*t3
  71  continue
  70  continue
  72  continue
  77  if (nd2.eq.0) go to 90
      do 82 k=1,ng
      if (skm(k).le.0.or.d(2,k).le.0) go to 82
      t4=(1-f)/skm(k)
      t5=1
c following line changed 3-24-04 - had been performed as integer
      if (d(2,k).gt.1) t5=1-(d(2,k)-1.d0)/(rs(k)-1.d0)
      t3=t5*((skmm(k)**2)*d(2,k))/(rs(k)**2)
      v3(k)=v3(k)+t4*t4*t3
      do 80 i=1,ng1
      t1=t4*c(i,k)
      v2(i,k)=v2(i,k)-t1*t4*t3
      do 81 j=1,i
      l=i*(i-1)/2+j
      t2=t4*c(j,k)
      v(l)=v(l)+t1*t2*t3
  81  continue
  80  continue
  82  continue
  90  if (lu.ge.n) go to 30
      do 91 i=ll,lu
      j=ig(i)
  91  rs(j)=rs(j)-1
      fm=f
      do 92 i=1,ng
      f1m(i)=f1(i)
  92  skmm(i)=skm(i)
      ll=lu+1
      lu=ll
      go to 50
  30  l=0
      do 36 i=1,ng1
      do 37 j=1,i
      l=l+1
      do 38 k=1,ng
      v(l)=v(l)+c(i,k)*c(j,k)*v3(k)
      v(l)=v(l)+c(i,k)*v2(j,k)
      v(l)=v(l)+c(j,k)*v2(i,k)
  38  continue
  37  continue
  36  continue
      return
      end
