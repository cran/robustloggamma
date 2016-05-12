# Biweight functions

rhoBW <- function(x,k){
k2  <- k*k; k4 <- k2*k2; k6 <- k2*k4
x2  <- x*x; x4 <- x2*x2; x6 <- x4*x2
(3*x2/k2-3*x4/k4+x6/k6)*(abs(x)<k)+(abs(x)>=k)}

psiBW <- function(x,k){
(6/k)*(x/k)*(1-(x/k)^2)^2*(abs(x)<k)}

pspBW <- function(x,k){
k2  <- k*k; k4 <- k2*k2; k6 <- k2*k4
x2  <- x*x; x4 <- x2*x2
(6/k2-36*x2/k4+30*x4/k6)*(abs(x)<k)}

rhoBWdnorm <- function(x,k) {rhoBW(x,k)*dnorm(x)}
