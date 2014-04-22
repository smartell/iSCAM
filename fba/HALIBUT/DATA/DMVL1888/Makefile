## Makefile for running iscam
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis
EXEC=iscam
prefix=../../../../dist
DAT=RUN.dat
CTL=Halibut2014
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

# |——————————————————————————————————————————————————————————————————————————————————————|
# | Targets for Running MSE
# |——————————————————————————————————————————————————————————————————————————————————————|
# | use gsed -i 's/STRING_TO_REPLACE/STRING_TO_REPLACE_IT/g' Makefile
NOSIM = 8
REPDIR = reps
simdirs := $(shell echo 'cat(formatC(1:$(NOSIM), digits=3, flag="0"))' | R --vanilla --slave)
.datadone:= $(foreach dir,$(simdirs),$(dir)/datadone)
runsims := $(foreach dir,$(simdirs),$(dir)/runsims)

$(.datadone): $(DAT) 
	mkdir -p $(REPDIR)/$(@D);
	cp Makefile $(REPDIR)/$(@D)  
	cd $(REPDIR)/$(@D); gsed -i 's/..\/..\/..\/..\/dist/..\/..\/..\/..\/..\/..\/dist/g' Makefile
	cd $(REPDIR)/$(@D); gsed -i 's/CTL=$(CTL)/CTL=..\/..\/$(CTL)/g' Makefile
	cd $(REPDIR)/$(@D); echo ../../$(CTL).dat > RUN.dat;
	cd $(REPDIR)/$(@D); echo ../../$(CTL).ctl >>RUN.dat;
	cd $(REPDIR)/$(@D); echo ../../$(CTL).pfc >>RUN.dat;
	cd $(REPDIR)/$(@D); echo ../../$(CTL).mpc >>RUN.dat;
	cd $(REPDIR)/$(@D); touch .datadone;

data: $(EXEC) $(.datadone)

$(runsims): 
	cd $(REPDIR)/$(@D); make ARG='-est -ainp ../../$(CTL).par -maxfn 0 -mse $(@D)'

mse: data $(runsims)
	touch simdone

cleanmse:
	rm -r 0* *.Rdata $(REPDIR)

