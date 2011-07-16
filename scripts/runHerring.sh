#!/bin/bash
runmle="off"
runmcmc="off"
runmceval="off"

niter=500000
nsave=500
nscale=100000

if [ $runmle = "on" ]; then
	./iscam -ind Examples/qciHerring/qci.dat -nox

	./iscam -ind Examples/prdHerring/prd.dat -nox

	./iscam -ind Examples/ccHerring/cc.dat -nox

	./iscam -ind Examples/sogHerring/sog.dat -nox

	./iscam -ind Examples/wcviHerring/wcvi.dat -nox	
	
	./iscam -ind Examples/Area2wHerring/Area2w.dat -nox	
	
	./iscam -ind Examples/Area27Herring/Area27.dat -nox	
fi

if [ $runmcmc = "on" ]; then
    ./iscam -ind Examples/qciHerring/qci.dat -mcmc $niter -mcsave $nsave -mcscale $nscale -nosdmcmc
    ./iscam -ind Examples/qciHerring/qci.dat -mceval
   
    ./iscam -ind Examples/prdHerring/prd.dat -mcmc $niter -mcsave $nsave -mcscale $nscale -nosdmcmc
    ./iscam -ind Examples/prdHerring/prd.dat -mceval
   
    ./iscam -ind Examples/ccHerring/cc.dat -mcmc $niter -mcsave $nsave -mcscale $nscale -nosdmcmc
    ./iscam -ind Examples/ccHerring/cc.dat -mceval

	./iscam -ind Examples/sogHerring/sog.dat -mcmc $niter -mcsave $nsave -mcscale $nscale -mcmult 2.0 -nosdmcmc
	./iscam -ind Examples/sogHerring/sog.dat -mceval

	./iscam -ind Examples/wcviHerring/wcvi.dat -mcmc $niter -mcsave $nsave -mcscale $nscale -nosdmcmc
	./iscam -ind Examples/wcviHerring/wcvi.dat -mceval
	
	./iscam -ind Examples/Area2wHerring/Area2w.dat -mcmc $niter -mcsave $nsave -mcscale $nscale -nosdmcmc
	./iscam -ind Examples/Area2wHerring/Area2w.dat -mceval
	
	./iscam -ind Examples/Area27Herring/Area27.dat -mcmc $niter -mcsave $nsave -mcscale $nscale -nosdmcmc
	./iscam -ind Examples/Area27Herring/Area27.dat -mceval
fi

if [ $runmceval = "on" ]; then
	cp iscam.psv tmp_iscam.psv
	
	cp examples/qciHerring/qciHerring2010HCAM.psv iscam.psv
	./iscam -ind Examples/qciHerring/qci.dat -mceval

	cp examples/prdHerring/prdHerring2010HCAM.psv iscam.psv
	./iscam -ind Examples/prdHerring/prd.dat -mceval

	cp examples/ccHerring/ccHerring2010HCAM.psv iscam.psv
	./iscam -ind Examples/ccHerring/cc.dat -mceval

	cp examples/sogHerring/sogHerring2010HCAM.psv iscam.psv
	./iscam -ind Examples/sogHerring/sog.dat -mceval

	cp examples/wcviHerring/wcviHerring2010HCAM.psv iscam.psv
	./iscam -ind Examples/wcviHerring/wcvi.dat -mceval
	
	cp tmp_iscam.psv iscam.psv
fi