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

SUBDIRS = AREA27 AREA2W CC PRD QCI SOG WCVI
RELEASE= "../../../dist/release/iscam"
.PHONY = est
.PHONY: est $(SUBDIRS)

all: release est clear


datadirs := $(shell echo 'd=dir()[file.info(dir())$$isdir]' | R --slave)
dirloop := $(foreach DIR,$(SUBDIRS),$(DIR)/cleaner)
estloop := $(foreach DIR,$(SUBDIRS),$(DIR)/est)
clrloop := $(foreach DIR,$(SUBDIRS),$(DIR)/clr)


COPYISCAM = 'datadirs=dir()[file.info(dir())$$isdir]; \
             sapply(datadirs, \
               function(d){\
                 file.copy($(RELEASE),paste(d,"/iscam",sep=""))\
               }\
             )'


#RUNEST = $(foreach DIR,$(DIROBJS),cp $(RELEASE) )

release: 
	echo $(COPYISCAM) | R --vanilla --slave

test: $(estloop)

$(estloop):
	cd $(@D); make ARG=$(ARG);

est: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ ARG=$(ARG)

clear: $(clrloop)

$(clrloop):
	cd $(@D); rm -f admodel.* eigv.rpt fmin.log variance iscam.*


clean: $(dirloop) 

$(dirloop):
	cd $(@D); make clean;



