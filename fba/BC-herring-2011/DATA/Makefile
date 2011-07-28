## Makefile for running the iscam program in each data directory
## Author: Steve Martell
## Usage: make <option>
## 
## Dependencies:
## Lines with the ":" in them are the dependency lines.
## They determine when the target is to be rebuilt.
## Example::  data : sim.R 
## if sim.R has been modified then data will be carried out.

ARG="-nox"
MCMC="'-nox -mcmc 5000 -mcsave 500 -mcscale 1000 -mcmult 2.0'"
#MCEVAL="'-mceval'"

SUBDIRS = AREA27 AREA2W CC PRD QCI SOG WCVI
RELEASE= "../../../dist/release/iscam"
.PHONY = est
.PHONY: est $(SUBDIRS)

all: release est clear


datadirs := $(shell echo 'd=dir()[file.info(dir())$$isdir]' | R --slave)
dirloop := $(foreach DIR,$(SUBDIRS),$(DIR)/cleaner)
clrloop := $(foreach DIR,$(SUBDIRS),$(DIR)/clr)
mcmcloop:= $(foreach DIR,$(SUBDIRS),$(DIR))


COPYISCAM = 'datadirs=dir()[file.info(dir())$$isdir]; \
             sapply(datadirs, \
               function(d){\
                 file.copy($(RELEASE),paste(d,"/iscam",sep=""))\
               }\
             )'


#RUNEST = $(foreach DIR,$(DIROBJS),cp $(RELEASE) )

release: 
	echo $(COPYISCAM) | R --vanilla --slave

est: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ ARG=$(ARG)

mcmc: $(mcmcloop) 
$(mcmcloop):
	$(MAKE) -C $@ ARG=$(MCMC); 





clear: $(clrloop)

$(clrloop):
	cd $(@D); rm -f admodel.* eigv.rpt fmin.log variance iscam.*


clean: $(dirloop) 

$(dirloop):
	cd $(@D); make clean;



