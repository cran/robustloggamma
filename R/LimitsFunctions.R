# Approximate limits of ratios
# ============================

# 1. R = f/S
# ----------

R1 <- function(x,lambda) {dloggamma(x,lambda=lambda)/(1-ploggamma(x,lambda=lambda))}

L1 <- function(x,lambda){zero = 0.001
res  = x
if (abs(lambda) >  zero) res = (exp(lambda*x)-1)/lambda
res}

LL1 <- function(x,lambda) {
pu = 0.999
xu = qloggamma(pu,lambda=lambda)+1
corr = R1(xu,lambda)-L1(xu,lambda)
res=L1(x,lambda)+corr
res}

# 2. R = F.dot/S
# --------------

R2 <- function(x,lambda){
num = apply(as.matrix(x),1,d.FL,lambda=lambda)
res = num/(1-ploggamma(x,lambda=lambda))
res}

L2 <- function(x,lambda){psiLG(x,lambda)}

LL2 <- function(x,lambda) {
pu = 0.999
xu = qloggamma(pu,lambda=lambda)+1
corr = R2(xu,lambda)-L2(xu,lambda)
res=L2(x,lambda)+corr
res}

# 3. R = f.prime/S
# ----------------

R3 <- function(x,lambda){p.loggamma(x,lambda=lambda)/(1-ploggamma(x,lambda=lambda))}

L3 <- function(x,lambda){zero=0.001
res= -x^2 + 1
if (abs(lambda) > zero ) res= -(1-exp(lambda*x))^2/lambda^2 + exp(lambda*x)
res}

LL3 <- function(x,lambda) {
pu = 0.999
xu = qloggamma(pu,lambda=lambda)+1
corr = R3(xu,lambda)-L3(xu,lambda)
res=L3(x,lambda)+corr
res}

# 4. R = d.fL/S
# -------------

R4 <- function(x,lambda){ d.fL(x,lambda)/ (1-ploggamma(x,lambda=lambda)) }

L4 <- function(u,lambda){zero=0.001
res= -u^4/6+u^2/2
if (abs(lambda)>zero) {
exuu = exp(lambda*u)
T1 = psiLG(u,lambda)*(1-exuu)/lambda
T2 = lambda^(-2)*(exuu*(1-lambda*u)-1)
res = T1-T2}
res}

LL4 <- function(x,lambda) {
pu = 0.999
xu = qloggamma(pu,lambda=lambda)+1
corr = R4(xu,lambda)-L4(xu,lambda)
res=L4(x,lambda)+corr
res}

# 5. R = d2.FL/S
# --------------

R5 <- function(x,lambda){
num =  apply(as.matrix(x),1,d2.FL,lambda=lambda)
num/(1-ploggamma(x,lambda=lambda))}

L5 <- function(x,lambda){
psi  = psiLG(x,lambda)
dpsi = apply(as.matrix(x),1,psiLG.dot,lambda=lambda)
res=dpsi-psi^2
res}

LL5 <- function(x,lambda) {
pu = 0.999
xu = qloggamma(pu,lambda=lambda)+1
corr = R5(xu,lambda)-L5(xu,lambda)
res=L5(x,lambda)+corr
res}
