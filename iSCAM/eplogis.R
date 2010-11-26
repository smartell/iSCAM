#Rcod for the exponential logistic selectivity model.
linf = 100
k = 0.25
m = 1.*k
A = 15
age= 1:A
la=linf*(1-exp(-k*age))
pa=exp(-m)^(age-1)*exp(rnorm(A, 0, 0.6))
pa[A]=pa[A]/(1-exp(-m))
pa=pa/sum(pa)

lbar=sum(pa*la)						#mean length of population
sigl=sum(la*pa*(1-pa))				#variance in length of population
abar = -log(-(lbar-linf)/linf)/k	#age of mean length in the population.

x=0:linf
plot(x, plogis(x, lbar, 0.2*lbar))
plot(age,plogis(la, lbar, 0.2*lbar))

#b=lbar
#a=sigl
#g=0.1
#x=0:linf
#sx=(1/(1-g))*((1-g)/g)^g*exp(1/a*g*(b-x))/(1+exp((b-x)/a))
#plot(x, sx)