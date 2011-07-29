## Makefile for running iscam
## Targets: default - run iscam with $(args)

args = -nox
datf = WCVI.dat
niter = 10000
nsave = 100

.PHONY: default all mcmc mceval clean

default:
	./iscam -ind $(datf) $(args)

all: default mcmc mceval clean

mcmc: 
	./iscam -ind $(datf) -mcmc $(niter) -mcsave $(nsave) -nosdmcmc -nox



mceval: iscam.psv
	-./iscam -ind $(datf) -mceval



clean: 
	-rm -f iscam.* admodel.* variance eigv.rpt fmin.log