## Makefile for running iscam
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis
EXEC=iscam
prefix=../../../../dist
DAT=RUN.dat
CTL=PacificHake14
ARG=
MCFLAG=-mcmc 10000 -mcsave 100 -nosdmcmc
NR=4

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
	./$(EXEC) -ind $(DAT) -mceval

$(EXEC).psv: $(CTL).ctl
	./$(EXEC) -ind $(DAT) $(MCFLAG) $(ARG)

mceval: $(EXEC)
	cp $(CTL).psv $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval

retro: $(EXEC) $(EXEC).ret1

$(EXEC).ret1:
	@echo $(RUNRETRO) | R --vanilla --slave

RUNRETRO = 'args = paste("-retro",c(1:$(NR),0)); \
            sapply(args,\
            function(a){ cmd=paste("./$(EXEC)","-ind $(DAT)",a);\
                        system(cmd)})'

clean: 
	-rm -f iscam.* admodel.* variance eigv.rpt fmin.log $(EXEC) variance