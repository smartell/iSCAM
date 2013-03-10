# CAPAM.mak 
# Makefile for the simualtions for CAPAM selectivity workshop.
## Makefile for running iscam
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis
EXEC=iscam
prefix=../../../../../dist
DAT=RUN.dat
CTL=PHake2010
ARG=
MCFLAG=-mcmc 10000 -mcsave 100 -nosdmcmc
NR=4
NOSIM = 8

.PHONY = all run mcmc mceval retro clean data

ifdef DEBUG
  DIST=$(prefix)/debug/iscam
else
  DIST=$(prefix)/release/iscam
endif

all: $(EXEC) $(EXEC).rep

$(EXEC): $(DIST)
	cp $(DIST) $@

$(EXEC).rep: $(DIST) $(CTL).ctl
	./$(EXEC) -ind $(DAT) $(ARG)

run:  $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG)

mcmc: $(EXEC) $(EXEC).psv
	./$(EXEC) -ind $(DAT) $(ARG) -mceval

$(EXEC).psv: $(CTL).ctl
	./$(EXEC) -ind $(DAT) $(MCFLAG) $(ARG)

mceval: $(EXEC)
	cp $(CTL).psv $(EXEC).psv
	./$(EXEC) -ind $(DAT) $(ARG) -mceval

retro: $(EXEC) $(EXEC).ret1

$(EXEC).ret1:
	@echo $(RUNRETRO) | R --vanilla --slave

RUNRETRO = 'args = paste("-retro",c(1:$(NR),0)); \
            sapply(args,\
            function(a){ cmd=paste("./$(EXEC)","-ind $(DAT) $(ARG)",a);\
                        system(cmd)})'

clean: 
	-rm -rf 0* iscam.* admodel.* variance eigv.rpt fmin.log $(EXEC) variance


simdirs := $(shell echo 'cat(formatC(1:$(NOSIM), digits=3, flag="0"))' | R --slave)
datadone := $(foreach dir,$(simdirs),$(dir)/datadone)
simfiles := $(foreach dir,$(simdirs),$(dir)/simdone)

$(datadone):
	mkdir $(@D)
	cd $(@D); touch datadone

data: $(datadone)


$(simfiles): data 
	cp  ./PHake2010.[cdp]* ./CAPAM.mak ./RUN.dat $(@D)
	cd  $(@D); make clean --file=CAPAM.mak; \
	make mcmc ARG="-nox -sim $(@D) -ainp PHake2010.pin" --file=CAPAM.mak; \
	#cp  ../../../../dist/release/iscam $(@D)
	

est: $(simfiles)





