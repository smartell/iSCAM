## Makefile for running iscam
## Targets: default - run iscam with $(args)

args = -nox
datf = AREA2W.dat
niter = 10000
nsave = 100

.PHONY: default all mcmc mceval clean

default:
	./iscam -ind $(datf) $(args)

all: mcmc mceval

mcmc: 
	./iscam -ind $(datf) -mcmc $(niter) -mcsave $(nsave) -nosdmcmc -nox
	./iscam -ind $(datf) -mceval


mceval: 
	-./iscam -ind $(datf) -mceval



clean: 
	-rm -f iscam.* admodel.* variance eigv.rpt fmin.log