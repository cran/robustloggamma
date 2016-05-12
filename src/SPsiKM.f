      SUBROUTINE SPsiKM(X, NX, TU, DPT, NT, DK, s)
      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)
      dimension x(nx), dnum(nx), den(nx), s(nx)
      dimension tu(nt), dpt(nt), y(nt)
      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      do 10 j=1,nx
        dnum(j) = dzero
        den(j) = dzero
 10   continue
      do 20 i=1,nt
        z = tu(i)/dk
        if (dabs(tu(i)).le.dk) then
          y(i) = (6.0d00*z-12.0d00*z**3.0d00+6.0d00*z**5.0d00)/dk
        else 
          y(i) = dzero
        endif
        do 30 j=1,nx
          if (tu(i).gt.x(j)) then
            dnum(j) = dnum(j) + y(i)*dpt(i)
          else
            den(j) = den(j) + dpt(i)
          endif
 30     continue
 20   continue
      do 40 j=1,nx
        s(j) = dnum(j)/(duno - den(j))        
 40   continue
      return
      end
