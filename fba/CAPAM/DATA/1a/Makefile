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
NOSIM  = 8   

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
# | REMOVE TEMPORARY ADMB FILES
# |------------------------------------------------------------------------------------|
# |
dust:
	rm -f *.log *.rpt *.htp admodel.* variance *.bar *.mcm

# |------------------------------------------------------------------------------------|
# | SIMULATIONS TO BE RUN IN PARALLEL IN NUMERIC DIRECTORIES
# |------------------------------------------------------------------------------------|
# | NOSIM determines the number of simulations.
# | simdirs	: is the list of simulation directories to copy material to.
# | datadone: is a loop for looping over directories

simdirs := $(shell echo 'cat(formatC(1:$(NOSIM), digits=3, flag="0"))' | R --vanilla --slave)
datadone:= $(foreach dir,$(simdirs),$(dir)/datadone)
runsims := $(foreach dir,$(simdirs),$(dir)/runsims)
# |------------------------------------------------------------------------------------|
# | BUILD DIRECTORIES AND COPY FILES INTO EACH DIRECTORY (target = data)
# |------------------------------------------------------------------------------------|
# |
$(datadone): 
	mkdir $(@D);
	cp  ./PHake2010.[cdp]*[!v] ./CAPAM.mak ./RUN.dat ./SIM.dat ./Sim1_PHake2010.ctl ./Sim1_PHake2010.pin $(@D)
	cd $(@D); make $(EXEC) --file=CAPAM.mak
	cd $(@D); make run DAT=SIM.dat ARG="-sim $(@D) -ainp Sim1_PHake2010.pin -noest" --file=CAPAM.mak
	cd $(@D); touch datadone

data: $(datadone)

cleansims:
	rm -r 0* *.Rdata

# |------------------------------------------------------------------------------------|
# | RUN SIMULATIONS IN EACH DIRECTORY (target = sims)
# |------------------------------------------------------------------------------------|
# |
$(runsims):
	cd $(@D); make retro --file=CAPAM.mak; make --file=CAPAM.mak; make dust --file=CAPAM.mak

sims: data $(runsims)
	touch allDONE

# |------------------------------------------------------------------------------------|
# | COLLECT SUMMARY STATISTICS FROM SIMULATION RUNS (target = collect)
# |------------------------------------------------------------------------------------|
# | COLLECTALL is an R-script to open report files and save output to allSims.Rdata
COLLECTALL ='dn<-dir(pattern="^[[:digit:]]"); \
				sims <- lapply(dn,function(d){require(Riscam);setwd(d);\
				A<-read.rep("PHake2010.rep");\
				B<-read.rep("iscam.sim");setwd("..");\
				bstatus <- log2(A$$sbt[33]/A$$bmsy)-log2(B$$sbt[33]/B$$bmsy) ;\
				fstatus <- log2(A$$ft[1,33]/A$$fmsy)-log2(B$$ft[1,33]/B$$fmsy);\
				nu      <- subset(A$$A_nu,A$$A_nu[,2]==1)[,-c(1,2)];\
				effN    <- sum(nu^2)/A$$age_tau2[1];\
				c(Fmsy=A$$fmsy,MSY=A$$msy,Bmsy=A$$bmsy,\
				  hat.Fmsy=B$$fmsy,hat.MSY=B$$msy,hat.Bmsy=B$$bmsy,\
				  Bo = A$$bo, hat.Bo = B$$bo,\
				  Fstatus.err=fstatus,Bstatus.err=bstatus,\
				  EffectiveN=effN,\
				  Age.tau2=A$$age_tau2[1],\
				  sigma=A$$sig,\
				  tau=A$$tau)}); \
		   		save(sims,file="allSims.Rdata")'

COLLECTRETRO =	'dn<-dir(pattern="^[[:digit:]]"); \
				runs <- lapply(dn,function(d){setwd(d);\
				source("../../retroStat.R");\
				setwd("..");\
				bias <- c(mean=mean(bias),abs.mean=mean(abs(bias)),bias)});\
				save(runs,file="retroSims.Rdata")'

SPAWNBIO = 'dn<-dir(pattern="^[[:digit:]]"); \
			spbio <- lapply(dn,function(d){require(Riscam);setwd(d);\
			A<-read.rep("PHake2010.rep");\
			sbt <- A$$sbt[1:33];setwd("..");return(sbt)});\
			save(spbio,file="spawnbio.Rdata")'

allSims.Rdata: allDONE
	@echo $(COLLECTALL) | R --vanilla --slave

retroSims.Rdata: allDONE
	@echo $(COLLECTRETRO) | R --vanilla --slave

spawnbio.Rdata: allDONE
	@echo $(SPAWNBIO) | R --vanilla	--slave

collect: allSims.Rdata retroSims.Rdata spawnbio.Rdata



