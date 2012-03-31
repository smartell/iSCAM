#R-script for halibut growth estimations.
require(ggplot2); require(hacks)
source("../Read.ADMB.R")
A=read.rep("../../DATA/Halibut_2sex_develop.rep")
iage=1:30

theta <- c(linf=155, k=0.1, sig=1/5)

vonb <- function(theta)
{
	with(as.list(theta), {
		ii   <- which(lobs!=0)
		lhat <- linf*(1-exp(-k*age))
		nu   <- lobs[ii]-lhat[ii]
		f    <- (-1)*sum(dnorm(nu, 0, 1/sig, log=TRUE))
		return(list(f=f, lhat=lhat))
	})
	
}
fn	<- function(theta) vonb(theta)$f
solver <- function(fn, theta, hess=TRUE)
{
	fit <- optim(theta, fn, method="BFGS", hessian=hess)
	H	<- fit$hessian
	V	<- solve(H)
	std <- sqrt(diag(V))
	R	<- V/(std%o%std)
	
	fit$V=V; fit$std=std; fit$R=R
	return(fit)
}


par(mfrow=c(2, 2), las=1, mar=c(4.6,  4.1,  2.1,  2.1))
# Females
l_f  <- data.frame(t(rbind(iage, A$lt_obs[1:16, ])))
w_f  <- data.frame(t(rbind(iage, A$wt_obs[1:16, ])))
ml_f <- melt(l_f, id.var=1)
age  <- ml_f[, 1]
lobs <- ml_f[, 3]
ii = which(lobs!=0)
plot(age[ii], lobs[ii], xlim=c(0, 30), ylim=c(0, 160), pch=20, col=colr(2, 0.5), xlab="Age", ylab="Length (cm)")
fit_f <- solver(fn, theta)
age = iage
lhat_f = vonb(fit_f$par)$lhat
lines(age, lhat_f)
par<-fit_f$par; par[1]=fit_f$par[1]+0.1*fit_f$par[1]
lh = vonb(par)$lhat
lines(age,lh)
par<-fit_f$par; par[1]=fit_f$par[1]-0.1*fit_f$par[1]
lh = vonb(par)$lhat
lines(age,lh)
gletter(1)

# Males
l_m  <- data.frame(t(rbind(iage, A$lt_obs[33:45, ])))
ml_m <- melt(l_m, id.var=1)
age  <- ml_m[, 1]
lobs <- ml_m[, 3]
ii = which(lobs!=0)
plot(age[ii], lobs[ii], xlim=c(0, 30), ylim=c(0, 160), pch=20, col=colr(4, 0.5), xlab="Age", ylab="Length (cm)")
fit_m <- solver(fn, theta)
age=iage
lhat_m = vonb(fit_m$par)$lhat
lines(age, lhat_m)
par<-fit_m$par; par[1]=fit_m$par[1]+0.1*fit_m$par[1]
lh = vonb(par)$lhat
lines(age,lh)
par<-fit_m$par; par[1]=fit_m$par[1]-0.1*fit_m$par[1]
lh = vonb(par)$lhat
lines(age,lh)
gletter(2)



x=seq(-3, 3, length=100)
px = 2*plogis(x,0, 0.75)-1
linf_f=fit_f$par[1]
linf_m=fit_m$par[1]
plot(x,linf_f +0.1*px*linf_f,type="l", xlab="Relative cohort density", ylab="Asymptotic length (cm)", ylim=c(90, 170) )
lines(x,linf_f -0.1*px*linf_f, lty=2);segments(-3, linf_f, 3, linf_f, lty=3)
gletter(3)
plot(x,linf_m +0.1*px*linf_m,type="l", xlab="Relative cohort density", ylab="Asymptotic length (cm)", ylim=c(90, 170) )
lines(x,linf_m -0.1*px*linf_m, lty=2);segments(-3, linf_m, 3, linf_m, lty=3)
gletter(4)
dev.copy2pdf(file="../../FIGURES/figLengthAtAgeFit.pdf")


