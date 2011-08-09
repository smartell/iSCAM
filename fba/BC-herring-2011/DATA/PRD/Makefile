## Makefile for running iscam
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis
EXEC=iscam
prefix=../../../../dist
DAT=-ind PRD.dat
ARG=
MCFLAG=-mcmc 10000 -mcsave 100 -nosdmcmc
NR=4

ifdef DEBUG
  DIST=$(prefix)/debug/iscam
else
  DIST=$(prefix)/release/iscam
endif

all: $(EXEC) $(EXEC).rep

$(EXEC):
	cp $(DIST) $@

$(EXEC).rep:
	./$(EXEC) $(DAT) $(ARG)

run:  $(EXEC)
	./$(EXEC) $(DAT) $(ARG)

mcmc: $(EXEC) $(EXEC).psv
	./$(EXEC) $(DAT) -mceval

$(EXEC).psv:
	./$(EXEC) $(DAT) $(MCFLAG) $(ARG)

retro: $(EXEC) $(EXEC).ret1

$(EXEC).ret1:
	@echo $(RUNRETRO) | R --vanilla --slave

RUNRETRO = 'args = paste("-retro",c(1:$(NR),0)); \
            sapply(args,\
            function(a){ cmd=paste("./$(EXEC)","$(DAT)",a);\
                        system(cmd)})'

clean: 
	-rm -f iscam.* admodel.* variance eigv.rpt fmin.log $(EXEC) variance