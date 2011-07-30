## Makefile for running iscam
## Targets: default - run iscam with $(args)

args = -nox
datf = SOG.dat
niter = 10000
nsave = 100
mcmult = 2.0


.PHONY: default all mcmc mceval clean

default:
	./iscam -ind $(datf) $(args)

all: mcmc mceval

mcmc: 
	./iscam -ind $(datf) -mcmc $(niter) -mcsave $(nsave) -nosdmcmc -nox -mcmult $(mcmult)
	./iscam -ind $(datf) -mceval


mceval: 
	-./iscam -ind $(datf) -mceval



clean: 
	-rm -f iscam.* admodel.* variance eigv.rpt fmin.log