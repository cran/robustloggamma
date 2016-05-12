      SUBROUTINE QnsemGLG(X, XMIN, XMAX, DP, Y, NDELTA, N, 
     &    DMU, DSIGMA, DLAMBDA, TOL, xout)
      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)
      logical neverdown
      dimension y(n), ndelta(n)
      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      external FnsemGLG

      itermax = ceiling(-(dlog(tol)/dlog(2.0d00)))+100
      i = 0
      step = x
      xold = dzero
      diff = duno + tol 
      dinterva = xmax - xmin
      neverdown = .TRUE.
      ddp = dzero
 300  if (i.le.itermax.and.diff.gt.tol) then
        xold = x
        xout = x*dinterva+xmin
        nc = 0
        do 10 j=1,n
          if (y(j).le.xout.and.ndelta(j).eq.0) then
            nc = nc + 1
          endif
 10     continue
        call FnsemGLG(xout, y, ndelta, n, nc, 
     &                dmu, dsigma, dlambda, ddp)
        if (ddp.lt.dp) then
          if (neverdown) then
            x = min(ddue*x, duno)
            step = x
          else
            x = x+step/(ddue**i)
            i = i + 1
          endif
        else
          neverdown = .FALSE.
          x = x-step/(ddue**i)
          i = i + 1
        endif
        diff = dabs(x - xold)
        goto 300
      endif
      xout = x*dinterva+xmin 
      return
      end
