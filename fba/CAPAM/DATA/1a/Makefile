## Makefile for running iscam
## Author: Steve Martell :)
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis

# |------------------------------------------------------------------------------------|
# | MACROS
# |------------------------------------------------------------------------------------|
# |
EXEC   = iscam
prefix =../../../../dist
DAT    = RUN.dat
CTL    = PHake2010
ARG    =
MCFLAG = -mcmc 20000 -mcsave 200 -nosdmcmc -nox
NR     = 4
NOSIM  = 16   


.PHONY = all run mcmc mceval retro clean data est
# |------------------------------------------------------------------------------------|
# | DEBUG FLAG
# |------------------------------------------------------------------------------------|
# |
ifdef DEBUG
  DIST=$(prefix)/debug/iscam
else
  DIST=$(prefix)/release/iscam
endif


# |------------------------------------------------------------------------------------|
# | COPY EXEC AND RUN MODEL
# |------------------------------------------------------------------------------------|
# |
all: $(EXEC) $(EXEC).par

$(EXEC): $(DIST)
	cp $(DIST) $@

$(EXEC).par: $(DIST) $(CTL).ctl
	./$(EXEC) -ind $(DAT) $(ARG)

run:  $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG)

# |------------------------------------------------------------------------------------|
# | MCMC and MCEVAL
# |------------------------------------------------------------------------------------|
# |
mcmc: $(EXEC) $(CTL).ctl $(EXEC).psv 
	./$(EXEC) -ind $(DAT) $(ARG) -mceval

$(EXEC).psv: $(CTL).ctl
	./$(EXEC) -ind $(DAT) $(MCFLAG) $(ARG)

mceval: $(EXEC)
	cp $(CTL).psv $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval

# |------------------------------------------------------------------------------------|
# | RETROSPECTIVE
# |------------------------------------------------------------------------------------|
# |
retro: $(EXEC) $(EXEC).ret1

$(EXEC).ret1:
	@echo $(RUNRETRO) | R --vanilla --slave

RUNRETRO = 'args = paste("-retro",c(1:$(NR),0),"-nox"); \
            sapply(args,\
            function(a){ cmd=paste("./$(EXEC)","-ind $(DAT) $(ARG)",a);\
                        system(cmd)})'


# |------------------------------------------------------------------------------------|
# | SIMULATIONS TO BE RUN IN PARALLEL IN NUMERIC DIRECTORIES
# |------------------------------------------------------------------------------------|
# | NOSIM determines the number of simulations.

simdirs := $(shell echo 'cat(formatC(1:$(NOSIM), digits=3, flag="0"))' | R --slave)
datadone:= $(foreach dir,$(simdirs),$(dir)/datadone)
simfiles:= $(foreach dir,$(simdirs),$(dir)/simdone)

$(datadone):
	mkdir $(@D)
	cp  ./PHake2010.[cdp]*[!v] ./CAPAM.mak ./RUN.dat $(@D)
	cd $(@D); make $(EXEC) --file=CAPAM.mak
	cd $(@D); touch datadone;

data: $(datadone)

$(simfiles):  
	cd  $(@D); make retro ARG="-nox -sim $(@D) -ainp PHake2010.pin" --file=CAPAM.mak
	cd  $(@D); make mcmc  ARG="-nox -sim $(@D) -ainp PHake2010.pin" --file=CAPAM.mak
	cd  $(@D); touch simdone
	
est: $(simfiles)

COLLECTALL = 	'dn<-dir(pattern="^[[:digit:]]"); \
				sims <- lapply(dn,function(d){require(Riscam);setwd(d);\
					A<-read.rep("PHake2010.rep");\
					B<-read.rep("iscam.sim");setwd("..");\
					bstatus <- log2(A$$sbt[33]/A$$bmsy)-log2(B$$sbt[33]/B$$bmsy) ;\
					fstatus <- log2(A$$ft[1,33]/A$$fmsy)-log2(B$$ft[1,33]/B$$fmsy);\
					c(A$$ENpar,A$$DIC,A$$fmsy,A$$msy,A$$bmsy,B$$fmsy,B$$msy,B$$bmsy,fstatus,bstatus)}); \
			   	save(sims,file="allSims.Rdata")'

allSims.Rdata:
	echo $(COLLECTALL) | R --vanilla --slave

collect: allSims.Rdata

COLLECTRETRO =  'dn<-dir(pattern="^[[:digit:]]"); \
				runs <- lapply(dn,function(d){setwd(d);\
					source("../../retroStat.R");\
					setwd("..");\
					bias <- c(mean=mean(bias),abs.mean=mean(abs(bias)),bias)});\
				save(runs,file="retroSims.Rdata")'

retroSims.Rdata:
	echo $(COLLECTRETRO) | R --vanilla --slave

retrostat: retroSims.Rdata

clean: 
	-rm -rf iscam.* admodel.* variance eigv.rpt fmin.log $(EXEC) variance 
	

cleanall:
	-rm -rf 0* allSims.Rdata allSims.R retroSims.Rdata
