c Copyright (C) 2000 Robert Gray
c   Distributed under the terms of the GNU public license
      subroutine tpoi(x,n,ind,tp5,ntp)
      double precision x(n),tp5(ntp)
      integer n,ntp,ind(ntp),l,k,i
      l=ntp
      do 10 i=ntp,1,-1
         if (x(n).ge.tp5(i)) go to 9
         ind(l)=0
         l=l-1
 10   continue 
 9    if (l.le.0) return
      if (x(n).eq.tp5(l)) then
         ind(l)=n
         l=l-1
      endif
c assuming unique values sorted in ascending order
      k=n-1
 11   if (l.le.0) return
      do 20 i=k,1,-1
         if (x(k).le.tp5(l)) then
            ind(l)=k+1
            l=l-1
            go to 11
         else
            k=k-1
         endif
 20   continue 
      do 30 i=1,l
         ind(l)=0
 30   continue 
      return
      end
