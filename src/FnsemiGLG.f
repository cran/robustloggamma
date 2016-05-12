      subroutine FnsemGLG(X, Y, NDELTA, N, NC, DMU, DSIGMA, DLAMBDA,
     + dp)

      implicit double precision (a-h,o,q-z)
      implicit integer (i,j,p,n)

      dimension y(n), ndelta(n), yc(nc), dpc(nc)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      external ploggamm
      dn = n
      dfirst = dzero
      ic = 0
      do 10 i=1,n
        if (y(i).le.x) then 
          if (ndelta(i).eq.1) then
            dfirst = dfirst + duno
          else 
            ic = ic + 1
            yc(ic) = (y(i)-dmu)/dsigma
          endif
        endif
 10   continue
      dsecond = dzero
      if (nc.gt.0) then
        xx = (x-dmu)/dsigma
        call ploggamm(xx, 1, dzero, duno, dlambda, dpx)
        call ploggamm(yc, nc, dzero, duno, dlambda, dpc)
        do 20 i=1,nc
          den = duno - dpc(i)
          if (xx.ge.yc(i)) then
            if (den.eq.dzero) then
              dsecond = dsecond + duno
            else
              dsecond = dsecond + (dpx - dpc(i))/den
            endif
          endif
 20     continue
      endif
      dp = (dfirst+dsecond)/dn      
      return
      end
