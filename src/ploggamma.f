      subroutine ploggamm(DQ, N, DMU, DSIGMA, DLAMBDA,
     + dp)

      implicit double precision (a-h,o,q-z)
      implicit integer (i,j,p,n)

      dimension dq(n),dp(n)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dqzero=0.0001)

      external dpgamma
      external dpnorm5

      do 10 i=1,n   
       ddq = (dq(i)-dmu)/dsigma
       if (dabs(dlambda) > dzero) then
         dalpha = duno/dlambda**ddue
         dx = dalpha*dexp(dlambda*ddq)
         call dpgamma(dx, dalpha, duno, 1, 0, dp(i))
       else
         call dpnorm5(ddq,dzero,duno,1,0,dp(i))
       endif
       if (dlambda.lt.-dqzero) then
         dp(i) = duno - dp(i)
       endif 
 10   continue
      return
      end
