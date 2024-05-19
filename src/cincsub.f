c Copyright (C) 2000 Robert Gray
c distributed under the terms of the GNU public license
      subroutine cinc(y,ic,icc,n,x,f,v)
c  calculates estimates of cumulative incidence functions and their
c   variances
c
c  i-n are integer, all others double precision
c
C  DATA ARE ASSUMED TO BE SORTED (increasing order of y).  
c  On input:
c    y is the failure times (assumed >= 0)
c    ic is 1 if the case has failed (from any cause), 0 otherwise
c    icc is 1 if the case has failed from the specific cause for which
c      the estimate is being calculated, 0 otherwise
c    n is the length of y (and ic and icc)
c
c    x,f,v need to have length at least 2*nfc+2
c      where nfc is the # of unique failure times from the cause of
c      interest
c  On output:
c    x gives the times, f the estimates, and v the variances
c      This program was designed to produce estimates suitable for
c      drawing as a step function.  x(0) is set to 0, and f(0) is also.
c      X(1) and x(2) are the value of the smallest time where icc=1
c      f(1) is 0, and f(2) is the value f jumps to at x(1)=x(2).
c      in this way all the corners of the step function are given in the
c      vector, and all the plotting routine has to do to make a plot is
c      to connect up the points.
C     THE LAST ELEMENT OF X AND F CONTNUE OUT TO THE largest follow-up time
      implicit double precision (a-h,o-z)
      dimension x(*),f(*),ic(n),icc(n),v(*),y(n)
      fk=1
      nf=0
      v1=0
      v2=0
      v3=0
      x(1)=0
      f(1)=0
      v(1)=0
      lcnt=1
      l=1
      ll=1
      rs=n
      ty=y(1)
  10  l=l+1
      if (l.gt.n) go to 60
      if (y(l).eq.ty) go to 10
  60  l=l-1
      nd1=0
      nd2=0
      do 15 i=ll,l
      nd1=nd1+icc(i)
      nd2=nd2+ic(i)-icc(i)
 15   continue
      nd=nd1+nd2
      if (nd.eq.0) go to 40
      fkn=fk*(rs-nd)/rs
      if (nd1.gt.0) then
      lcnt=lcnt+2
      f(lcnt-1)=f(lcnt-2)
      f(lcnt)=f(lcnt-1)+fk*nd1/rs
      end if
      if (nd2.le.0.or.fkn.le.0) go to 30
      t5=1
      if (nd2.gt.1) t5=1-(nd2-1.)/(rs-1.)
      t6=fk*fk*t5*nd2/(rs*rs)
      t3=1./fkn
      t4=f(lcnt)/fkn
      v1=v1+t4*t4*t6
      v2=v2+t3*t4*t6
      v3=v3+t3*t3*t6
  30  if (nd1.le.0) go to 35
      t5=1
      if (nd1.gt.1) t5=1-(nd1-1.)/(rs-1.)
      t6=fk*fk*t5*nd1/(rs*rs)
      t3=0
      if (fkn.gt.0) t3=1/fkn
      t4=1+t3*f(lcnt)
      v1=v1+t4*t4*t6
      v2=v2+t3*t4*t6
      v3=v3+t3*t3*t6
c  changed 3-18-91: If largest follow-up time is a failure of type 1,
c  only v1 should be incremented (as above), but to get the right
c variance still need to incorporate v2 and v3 in v(lcnt)
c      t2=0
c      if (fkn.gt.0) t2=f(lcnt)
      t2=f(lcnt)
      x(lcnt-1)=y(l)
      x(lcnt)=y(l)
      v(lcnt-1)=v(lcnt-2)
      v(lcnt)=v1+t2*t2*v3-2*t2*v2
  35  fk=fkn
      nf=nf+nd1
  40  rs=n-l
      l=l+1
      if (l.gt.n) go to 50
      ll=l
      ty=y(l)
      go to 10
  50  lcnt=lcnt+1
      x(lcnt)=y(n)
      f(lcnt)=f(lcnt-1)
      v(lcnt)=v(lcnt-1)
      return
      end
