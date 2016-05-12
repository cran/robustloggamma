# Score functions and Jacobians for the case with censored observations
# =====================================================================

Uscore.full <- function(y,w,mu,sigma,lambda) {
# Score functions vector / error model / complete observations components
  z  <- (y-mu)/sigma
  z1 <- csiLG(z,lambda)/sigma
  z2 <- (csiLG(z,lambda)*z+1)/sigma
  z3 <- psiLG(z,lambda)
  u1 <- mean(w*z1)
  u2 <- mean(w*z2)
  u3 <- mean(w*z3)
matrix(c(u1,u2,u3),nrow=3,byrow=TRUE)}

Uscore.cens <- function(y,w,mu,sigma,lambda) {
# Score functions vector /error model / censored observations components
  z  <- (y-mu)/sigma
  z1 <- csiLG(z,lambda)/sigma
  Fl.dot <- apply(as.matrix(z),1,d.FL,lambda=lambda)
  fl <- dloggamma(z,lambda=lambda)
  Sl <- 1-ploggamma(z,lambda=lambda)  
  r1 <- fl/Sl
  if (length(r1[is.na(r1) | abs(r1)==Inf])>0) r1[is.na(r1)| abs(r1)==Inf] = LL1(z,lambda)[is.na(r1)| abs(r1)==Inf]
  r2 <- Fl.dot/Sl
  if (length(r2[is.na(r2) | abs(r2)==Inf])>0) r2[is.na(r2)| abs(r2)==Inf] = LL2(z,lambda)[is.na(r2)| abs(r2)==Inf]
  s1 <- -r1/sigma
  s2 <- s1*z
  s3 <- r2  
  u1 <- mean(w*s1)
  u2 <- mean(w*s2)
  u3 <- mean(w*s3)
matrix(c(u1,u2,u3),nrow=3,byrow=TRUE)}

###########################

Uscore.full.r <- function(y,X,w,mu,beta,sigma,lambda) {
# Score functions vector / regression / complete observations components
  IX <- cbind(1,X)
  m <- length(y)
  p <- ncol(X)
  z  <- (y-X%*%beta-mu)/sigma
  z1 <- csiLG(z,lambda)/sigma
  z2 <- (csiLG(z,lambda)*z+1)/sigma
  z3 <- psiLG(z,lambda)
  u1 <- t(w*z1)%*%IX/m
  u2 <- mean(w*z2)
  u3 <- mean(w*z3)
matrix(c(u1,u2,u3),nrow=(3+p),byrow=TRUE)}

Uscore.cens.r <- function(y,X,w,mu,beta,sigma,lambda) {
# Score functions vector / regression / censored observations components
  IX <- cbind(1,X)
  m <- length(y)
  p <- ncol(X)
  z  <- (y-X%*%beta-mu)/sigma
  z1 <- csiLG(z,lambda)/sigma
  Fl.dot <- apply(as.matrix(z),1,d.FL,lambda=lambda)
  fl <- dloggamma(z,lambda=lambda)
  Sl <- 1-ploggamma(z,lambda=lambda)  
  r1 <- fl/Sl
  if (length(r1[is.na(r1) | abs(r1)==Inf])>0) r1[is.na(r1)| abs(r1)==Inf] = LL1(z,lambda)[is.na(r1)| abs(r1)==Inf]
  r2 <- Fl.dot/Sl
  if (length(r2[is.na(r2) | abs(r2)==Inf])>0) r2[is.na(r2)| abs(r2)==Inf] = LL2(z,lambda)[is.na(r2)| abs(r2)==Inf]
  s1 <- -r1/sigma
  s2 <- s1*z
  s3 <- r2  
  u1 <-  t(w*s1)%*%IX/m
  u2 <- mean(w*s2)
  u3 <- mean(w*s3)
matrix(c(u1,u2,u3),nrow=(3+p),byrow=TRUE)}

##########################################

Jacobian.full <- function(y,w,mu,sigma,lambda){
# Empirical Jacobian / error model / complete observations
z   =  (y-mu)/sigma; sigma2=sigma^2
z2  =  (csiLG(z,lambda)*z+1)/sigma
z11 = -csiLG.prime(z,lambda)/sigma2
z12 = -(csiLG.prime(z,lambda)*z + csiLG(z,lambda))/sigma2
z13 =  csiLG.dot(z,lambda)/sigma
z21 =  z12
z22 =  z12*z-z2/sigma
z23 =  z13*z
z31 = -psiLG.prime(z,lambda)/sigma
z32 =  z31*z
z33 =  psiLG.dot(z,lambda)
z11 =  mean(w*z11)
z12 =  mean(w*z12)
z13 =  mean(w*z13)
z21 =  mean(w*z21) 
z22 =  mean(w*z22)  
z23 =  mean(w*z23) 
z31 =  mean(w*z31) 
z32 =  mean(w*z32)
z33 =  mean(w*z33)
matrix(c(z11,z12,z13,z21,z22,z23,z31,z32,z33),nrow=3,byrow=TRUE)}

Jacobian.cens <- function(y,w,mu,sigma,lambda){
# Empirical Jacobian / error model / censored observations
z        = (y-mu)/sigma
fl       = dloggamma(z,lambda=lambda)
Sl       = 1-ploggamma(z,lambda=lambda)
fl.prime = p.loggamma(z,lambda=lambda)
Fl.dot   = apply(as.matrix(z),1,d.FL,lambda=lambda)
Fl.dot2  = apply(as.matrix(z),1,d2.FL,lambda=lambda)
fl.dot   = d.fL(z,lambda=lambda)
r1 = fl/Sl
r2 = Fl.dot/Sl
r3 = fl.prime/Sl
r4 = fl.dot/Sl
r5 = Fl.dot2/Sl
if (length(r1[is.na(r1) | abs(r1)==Inf])>0) r1[is.na(r1)| abs(r1)==Inf] = LL1(z,lambda)[is.na(r1)| abs(r1)==Inf]
if (length(r2[is.na(r2) | abs(r2)==Inf])>0) r2[is.na(r2)| abs(r2)==Inf] = LL2(z,lambda)[is.na(r2)| abs(r2)==Inf]
if (length(r3[is.na(r3) | abs(r3)==Inf])>0) r3[is.na(r3)| abs(r3)==Inf] = LL3(z,lambda)[is.na(r3)| abs(r3)==Inf]
if (length(r4[is.na(r4) | abs(r4)==Inf])>0) r4[is.na(r4)| abs(r4)==Inf] = LL4(z,lambda)[is.na(r4)| abs(r4)==Inf]
if (length(r5[is.na(r5) | abs(r5)==Inf])>0) r5[is.na(r5)| abs(r5)==Inf] = LL5(z,lambda)[is.na(r5)| abs(r5)==Inf]
a   = r3+r1^2
b   = r4+r1*r2
d   = r5+r2^2
s1  = -r1/sigma      
s2  = s1*z
s3  = r2
s11 = a/sigma^2 ;  s12 = (r1+a*z)/sigma^2 ;  s13 = -b/sigma
s21 = s12       ;  s22 = s12*z-s2/sigma   ;  s23 = s13*z
s31 = s13       ;  s32 = s23              ;  s33 = d
s11 = mean(w*s11)
s12 = mean(w*s12)
s13 = mean(w*s13)
s21 = mean(w*s21)
s22 = mean(w*s22)
s23 = mean(w*s23)
s31 = mean(w*s31)
s32 = mean(w*s32)
s33 = mean(w*s33)
matrix(c(s11,s12,s13,s21,s22,s23,s31,s32,s33),nrow=3,byrow=TRUE)}

##########################################

Exp.Jacobian.full <- function(mu,sigma,lambda,nmod=1000,x,y2) {
# Expected Jacobian / error model / complete observations
yy   =  qloggamma(ppoints(nmod), mu=mu, sigma=sigma, lambda=lambda)
Js0 <- splinefun(x,y2, method= "monoH.FC")
Js  <- function(x){ js0 = Js0(x); pmax(0, pmin(1,js0))}
Gy  =  1-Js(yy) 
z   =  (yy-mu)/sigma; sigma2=sigma^2
z2  =  (csiLG(z,lambda)*z+1)/sigma
z11 = -csiLG.prime(z,lambda)/sigma2
z12 = -(csiLG.prime(z,lambda)*z + csiLG(z,lambda))/sigma2
z13 =  csiLG.dot(z,lambda)/sigma
z21 =  z12
z22 =  z12*z-z2/sigma
z23 =  z13*z
z31 = -psiLG.prime(z,lambda)/sigma
z32 =  z31*z
z33 =  psiLG.dot(z,lambda)
z11 =  mean(z11*Gy)
z12 =  mean(z12*Gy)
z13 =  mean(z13*Gy)
z21 =  mean(z21*Gy) 
z22 =  mean(z22*Gy)  
z23 =  mean(z23*Gy) 
z31 =  mean(z31*Gy) 
z32 =  mean(z32*Gy)
z33 =  mean(z33*Gy)
matrix(c(z11,z12,z13,z21,z22,z23,z31,z32,z33),nrow=3,byrow=TRUE)}

Exp.Jacobian.cens <- function(mu,sigma,lambda,xmin,xmax,qemp,nmod=1000,x=x,y2=y2) {
# Expected Jacobian / error model / censored observations
yy   = J.Quant0(nmod=nmod,xmin,xmax,x=x,y2=y2)
Sy  = 1-ploggamma(yy, mu=mu, sigma=sigma, lambda=lambda)
z        = (yy-mu)/sigma
fl       = dloggamma(z,lambda=lambda)
Sl       = 1-ploggamma(z,lambda=lambda)
fl.prime = p.loggamma(z,lambda=lambda)
Fl.dot   = apply(as.matrix(z),1,d.FL, lambda=lambda)
Fl.dot2  = apply(as.matrix(z),1,d2.FL,lambda=lambda)
fl.dot   = d.fL(z,lambda=lambda)
r1 = fl/Sl
r2 = Fl.dot/Sl
r3 = fl.prime/Sl
r4 = fl.dot/Sl
r5 = Fl.dot2/Sl
if (length(r1[is.na(r1) | abs(r1)==Inf])>0) r1[is.na(r1)| abs(r1)==Inf] = LL1(z,lambda)[is.na(r1)| abs(r1)==Inf]
if (length(r2[is.na(r2) | abs(r2)==Inf])>0) r2[is.na(r2)| abs(r2)==Inf] = LL2(z,lambda)[is.na(r2)| abs(r2)==Inf]
if (length(r3[is.na(r3) | abs(r3)==Inf])>0) r3[is.na(r3)| abs(r3)==Inf] = LL3(z,lambda)[is.na(r3)| abs(r3)==Inf]
if (length(r4[is.na(r4) | abs(r4)==Inf])>0) r4[is.na(r4)| abs(r4)==Inf] = LL4(z,lambda)[is.na(r4)| abs(r4)==Inf]
if (length(r5[is.na(r5) | abs(r5)==Inf])>0) r5[is.na(r5)| abs(r5)==Inf] = LL5(z,lambda)[is.na(r5)| abs(r5)==Inf]
a   = r3+r1^2
b   = r4+r1*r2
d   = r5+r2^2
s1  = -r1/sigma      
s2  = s1*z
s3  = r2
s11 = a/sigma^2 ;  s12 = (r1+a*z)/sigma^2 ;  s13 = -b/sigma
s21 = s12       ;  s22 = s12*z-s2/sigma   ;  s23 = s13*z
s31 = s13       ;  s32 = s23              ;  s33 = d
Semp  = density(qemp, kernel="gaussian", bw=0.3, cut=3,n=512)
semp  = approxfun(x=Semp$x, y=Semp$y, rule=2)
qmodl = qloggamma(ppoints(nmod),lambda=lambda)
Smod  = density(qmodl, kernel="gaussian", bw=0.3, cut=3,n=512)
smod  = approxfun(x=Smod$x, y=Smod$y, rule=2)
pe2    = semp(z)
pm2    = smod(z)
Delta = pe2/pm2-1
w2     = pesi(x=Delta, raf="NED",tau=0.5)
w2[w2<0.04] = 0
s11 = mean(s11*Sy*w2)
s12 = mean(s12*Sy*w2)
s13 = mean(s13*Sy*w2)
s21 = mean(s21*Sy*w2)
s22 = mean(s22*Sy*w2)
s23 = mean(s23*Sy*w2)
s31 = mean(s31*Sy*w2)
s32 = mean(s32*Sy*w2)
s33 = mean(s33*Sy*w2)
matrix(c(s11,s12,s13,s21,s22,s23,s31,s32,s33),nrow=3,byrow=TRUE)}

##########################################

Jacobian.full.r <- function(y,X,w,mu,beta,sigma,lambda){
# Empirical Jacobian / regression / complete observations
  IX <- cbind(1,X)
  m <- length(y)
  p <- ncol(X)
  z  <- (y-X%*%beta-mu)/sigma  ; sigma2=sigma^2
z2  =  (csiLG(z,lambda)*z+1)/sigma
z11 = -csiLG.prime(z,lambda)/sigma2
z12 = -(csiLG.prime(z,lambda)*z + csiLG(z,lambda))/sigma2
z13 =  csiLG.dot(z,lambda)/sigma
z21 =  z12
z22 =  z12*z-z2/sigma
z23 =  z13*z
z31 = -psiLG.prime(z,lambda)/sigma
z32 =  z31*z
z33 =  psiLG.dot(z,lambda)
A = t(IX)%*%(IX*rep(w*z11,p+1))/m
B = matrix(c(mean(w*z22),mean(w*z23),mean(w*z32),mean(w*z33)),nrow=2,byrow=T)
C = rbind(t(w*z21)%*%IX/m,t(w*z31)%*%IX/m)
D = t(rbind(t(w*z12)%*%IX/m,t(w*z13)%*%IX/m))
rbind(cbind(A,D),cbind(C,B)) }

Jacobian.cens.r <- function(y,X,w,mu,beta,sigma,lambda){
# Empirical Jacobian / regression / censored observations
  IX <- cbind(1,X)
  m <- length(y)
  p <- ncol(X)
  z  <- (y-X%*%beta-mu)/sigma  ; sigma2=sigma^2
fl       = dloggamma(z,lambda=lambda)
Sl       = 1-ploggamma(z,lambda=lambda)
fl.prime = p.loggamma(z,lambda=lambda)
Fl.dot   = apply(as.matrix(z),1,d.FL,lambda=lambda)
Fl.dot2  = apply(as.matrix(z),1,d2.FL,lambda=lambda)
fl.dot   = d.fL(z,lambda=lambda)
r1 = fl/Sl
r2 = Fl.dot/Sl
r3 = fl.prime/Sl
r4 = fl.dot/Sl
r5 = Fl.dot2/Sl
if (length(r1[is.na(r1) | abs(r1)==Inf])>0) r1[is.na(r1)| abs(r1)==Inf] = LL1(z,lambda)[is.na(r1)| abs(r1)==Inf]
if (length(r2[is.na(r2) | abs(r2)==Inf])>0) r2[is.na(r2)| abs(r2)==Inf] = LL2(z,lambda)[is.na(r2)| abs(r2)==Inf]
if (length(r3[is.na(r3) | abs(r3)==Inf])>0) r3[is.na(r3)| abs(r3)==Inf] = LL3(z,lambda)[is.na(r3)| abs(r3)==Inf]
if (length(r4[is.na(r4) | abs(r4)==Inf])>0) r4[is.na(r4)| abs(r4)==Inf] = LL4(z,lambda)[is.na(r4)| abs(r4)==Inf]
if (length(r5[is.na(r5) | abs(r5)==Inf])>0) r5[is.na(r5)| abs(r5)==Inf] = LL5(z,lambda)[is.na(r5)| abs(r5)==Inf]
a   = r3+r1^2
b   = r4+r1*r2
d   = r5+r2^2
s1  = -r1/sigma      
s2  = s1*z
s3  = r2
s11 = a/sigma^2 ;  s12 = (r1+a*z)/sigma^2 ;  s13 = -b/sigma
s21 = s12       ;  s22 = s12*z-s2/sigma   ;  s23 = s13*z
s31 = s13       ;  s32 = s23              ;  s33 = d
A = t(IX)%*%(IX*rep(w*s11,p+1))/m
B = matrix(c(mean(w*s22),mean(w*s23),mean(w*s32),mean(w*s33)),nrow=2,byrow=T)
C = rbind(t(w*s21)%*%IX/m,t(w*s31)%*%IX/m)
D = t(rbind(t(w*s12)%*%IX/m,t(w*s13)%*%IX/m))
rbind(cbind(A,D),cbind(C,B)) }

Exp.Jacobian.full.reg <- function(yo,Xo,wo,mu,sigma,lambda,nmod=1000,x,v2) {
#
# Expected Jacobian, non-censored part
# note that yo, Xo, and wo are the non-censored components of y, X, and w
#
IX <- cbind(1,Xo)
m  <- length(yo)
p  <- ncol(Xo)
Js0 <- splinefun(x, v2, method= "monoH.FC")
Js  <- function(x){ js0 = Js0(x); pmax(0, pmin(1,js0))}
u   =  qloggamma(ppoints(nmod), mu=0, sigma=1, lambda=lambda)
Gu  =  1 - Js(mu+sigma*u) 
sigma2 = sigma^2
Au     = -csiLG.prime(u,lambda)/sigma2
z11    =  mean(Au*Gu)
Bu     = -(csiLG.prime(u,lambda)*u + csiLG(u,lambda))/sigma2
z12    =  mean(Bu*Gu)
Cu     =  csiLG.dot(u,lambda)/sigma
z13    =  mean(Cu*Gu)
z21    =  z12
Du     = -(csiLG.prime(u,lambda)*u^2 + 2*csiLG(u,lambda)*u + 1 )/sigma2
z22    =  mean(Du*Gu)
z23    =  mean(Cu*u*Gu)
z31    =  z13
z32    =  z23
z33    =  mean(psiLG.dot(u,lambda)*Gu)
#
A = t(IX)%*%(IX*rep(wo*z11,p+1))/m
B = matrix(c(mean(wo*z22),mean(wo*z23),mean(wo*z32),mean(wo*z33)),nrow=2,byrow=TRUE)
C = rbind(t(wo*z21)%*%IX/m,t(wo*z31)%*%IX/m)
D = t(rbind(t(wo*z12)%*%IX/m,t(wo*z13)%*%IX/m))
EJf = rbind(cbind(A,D),cbind(C,B)) 
EJf}

Exp.Jacobian.cens.reg <- function(yc,Xc,wc,mu,sigma,lambda,vmin,vmax,nmod=1000,x,v2) {
#
# Expected Jacobian, censored part
# Note that yc, Xc, and wc are the censored components of y, X, and w
#
IX <- cbind(1,Xc)
m  <- length(yc)
p  <- ncol(Xc)

y2 <- v2
q <- J.Quant0(nmod=nmod,vmin,vmax,tol=1e-4,x=x,y2=y2)
u <- (q-mu)/sigma
Su <- 1-ploggamma(u, mu=0, sigma=1, lambda=lambda)
#
sigma2 <- sigma^2
fl <- dloggamma(u,lambda=lambda)
Sl <- 1-ploggamma(u,lambda=lambda)
fl.prime <- p.loggamma(u,lambda=lambda)
Fl.dot <- apply(as.matrix(u),1,d.FL,lambda=lambda)
Fl.dot2 <- apply(as.matrix(u),1,d2.FL,lambda=lambda)
fl.dot <- d.fL(u,lambda=lambda)
r1 <- fl/Sl
r2 <- Fl.dot/Sl
r3 <- fl.prime/Sl
r4 <- fl.dot/Sl
r5 <- Fl.dot2/Sl
if (length(r1[is.na(r1) | abs(r1)==Inf])>0) r1[is.na(r1)| abs(r1)==Inf] <- LL1(u,lambda)[is.na(r1)| abs(r1)==Inf]
if (length(r2[is.na(r2) | abs(r2)==Inf])>0) r2[is.na(r2)| abs(r2)==Inf] <- LL2(u,lambda)[is.na(r2)| abs(r2)==Inf]
if (length(r3[is.na(r3) | abs(r3)==Inf])>0) r3[is.na(r3)| abs(r3)==Inf] <- LL3(u,lambda)[is.na(r3)| abs(r3)==Inf]
if (length(r4[is.na(r4) | abs(r4)==Inf])>0) r4[is.na(r4)| abs(r4)==Inf] <- LL4(u,lambda)[is.na(r4)| abs(r4)==Inf]
if (length(r5[is.na(r5) | abs(r5)==Inf])>0) r5[is.na(r5)| abs(r5)==Inf] <- LL5(u,lambda)[is.na(r5)| abs(r5)==Inf]
a   <- r3+r1^2
b   <- r4+r1*r2
d   <- r5+r2^2
#
s11  <- mean(a*Su)/sigma2
Au   <- (r1+a*u)/sigma^2
s12  <- mean(Au*Su)
s13  <- mean(-b*Su)/sigma
s21  <- s12 
Bu   <- (r1+a*u+r1*u)/sigma2
s22  <- mean(Bu*Su)
Cu   <- -b*u/sigma
s23  <- mean(Cu*Su)
s31  <- s13 
s32  <- s23 
s33  <- mean(d*Su)
#
A   <- t(IX)%*%(IX*rep(wc*s11,p+1))/m
B   <- matrix(c(mean(wc*s22),mean(wc*s23),mean(wc*s32),mean(wc*s33)),nrow=2,byrow=T)
C   <- rbind(t(wc*s21)%*%IX/m,t(wc*s31)%*%IX/m)
D   <- t(rbind(t(wc*s12)%*%IX/m,t(wc*s13)%*%IX/m))
EJc <- rbind(cbind(A,D),cbind(C,B)) 
EJc}
