#R-script for halibut growth estimations.
require(ggplot2); require(hacks)
source("../Read.ADMB.R")

#Careful if you uncomment this, will use the allometric relationship for weight at age
#A=read.rep("../../DATA/Halibut_2sex_develop.rep")
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


par(mfcol=c(2, 2), las=1, mar=c(4.6,  4.1,  2.1,  2.1))
# Females
l_f  <- data.frame(t(rbind(iage, A$lt_obs[1:16, ])))
w_f  <- data.frame(t(rbind(iage, A$wt_obs[1:16, ])))
ml_f <- melt(l_f, id.var=1)

fage <- ml_f[, 1]
lobs <- ml_f[, 3]

ii = which(lobs!=0)
plot(fage[ii], lobs[ii], xlim=c(0, 30), ylim=c(0, 160), pch=20, col=colr(2, 0.5), xlab="Age", ylab="Length (cm)")
age=fage
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

wt  =as.vector(A$wt_obs[1:16, ])
lt  =as.vector(A$lt_obs[1:16, ])
ii  = which(wt!=0)
xx = (lt[ii]); yy = (wt[ii])
plot(xx, yy, xlab="Length (in)", ylab="Weight (lb)", xlim=c(0, max(xx)), ylim=c(0, max(yy)))
#x2 <- xx^2
#x3 <- xx^3
#gf <- glm(yy~xx+x2+x3)
#points(xx,gf$fitted.values,col=2,pch=19)
xx=seq(0, max(xx))
yy=(9.321e-6*xx^3.16)
points(xx, yy, pch=19, col=2)
gletter(2)
print("female")
print(coef(gf))
gff <- gf

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
gletter(3)

# Males use colums 8:30,  otherwise the males are shrinking
wt  =as.vector(A$wt_obs[33:45, ])
lt  =as.vector(A$lt_obs[33:45, ])
ii  = which(wt!=0)
xx = subset(lt, wt!=0)
yy = subset(wt, wt!=0)


xx = (xx); yy = (yy)
plot(xx, yy, xlab="log Length (in)", ylab="Weight (lb)", xlim=c(0, max(xx)), ylim=c(0, max(yy)))
#x2 <- xx^2
#x3 <- xx^3
#gf <- glm(yy~xx+x2+x3)
#points(xx,gf$fitted.values,col=4,pch=19)
#y2<-predict(gff,data.frame(xx=xx,x2=x2,x3=x3))
xx=seq(0, max(xx))
yy=(9.321e-6*xx^3.16)
points(xx,yy,col=4,pch=19)
gletter(4)
print("male")
print(coef(gf))

dev.copy2pdf(file="../../FIGURES/figLengthAtAgeFit.pdf")


par(mfcol=c(1, 2))
x=seq(-3, 3, length=100)
px = 2*plogis(x,0, 0.75)-1
linf_f=fit_f$par[1]
linf_m=fit_m$par[1]
plot(x,linf_f +0.1*px*linf_f,type="l", xlab="Relative cohort density", ylab="Asymptotic length (cm)", ylim=c(90, 170) )
lines(x,linf_f -0.1*px*linf_f, lty=2);segments(-3, linf_f, 3, linf_f, lty=3)
gletter(1)
plot(x,linf_m +0.1*px*linf_m,type="l", xlab="Relative cohort density", ylab="Asymptotic length (cm)", ylim=c(90, 170) )
lines(x,linf_m -0.1*px*linf_m, lty=2);segments(-3, linf_m, 3, linf_m, lty=3)
gletter(2)
dev.copy2pdf(file="../../FIGURES/figLinf.pdf", height=4, width=8)


