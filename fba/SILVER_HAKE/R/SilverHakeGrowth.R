# R-script for estimating growth parameters for silver hake

df	<- read.table("SilverHakeGrowthData.txt", header=TRUE)
df	<- na.omit(df)
df$length <- df$length*100

# parameters
linf	<- 80
k		<- 0.4
to		<- -0.2
cv		<- 0.1
theta	<- c(linf=linf, k=k, to=to, cv=cv)

# routines
.vonb	<- function(theta, df)
{
	with(as.list(theta), {
		lhat     <- linf*(1.-exp(-k*(df$age-to)))
		sig      <- cv * lhat
		epsilon  <- df$length - lhat
		nll      <- (-1.0)*sum(dnorm(epsilon, 0, sig, log=TRUE))
		
		return(list(nll=nll, lhat=lhat, epsilon=epsilon))
	})
}

fn		<- function(theta, df) return(.vonb(theta, df)$nll)

fit		<- optim(theta, fn, method="BFGS", hessian=TRUE, df=na.omit(df))

F <- .vonb(fit$par, df)
plot(df$age,df$length)
points(na.omit(df$age),F$lhat,col=2)
lines(sort(na.omit(df$age)),sort(F$lhat),col=2)
fit$par
