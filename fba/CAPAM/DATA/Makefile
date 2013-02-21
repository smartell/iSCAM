## Makefile for running models
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
SUBDIR = 1a 1b 1c 2a 2b 2c 3a 3b 3c
#SUBDIR = 3a 3b 3c
##SUBDIR = SOKeg
TARGET = 
SEED  = 999 
JOBS  = 8

.PHONY: default $(SUBDIR) mcmc


## Targets
default: $(SUBDIR)
$(SUBDIR):
	cd $@ && $(MAKE) $(TARGET)



.PHONY: clean
clean_files := $(foreach dir,$(SUBDIR),$(dir)/clean)

clean: $(clean_files)
$(clean_files):
	cd $(@D) && $(MAKE) clean

.PHONY: sim
sim_dirs := $(foreach dir,$(SUBDIR),$(dir)/sim)
run_dirs := $(foreach dir,$(SUBDIR),$(dir)/run)

sim: $(sim_dirs)
$(sim_dirs):
	cd $(@D); make clean; make data; make est -j$(JOBS); make collect #; make ARG="-sim $(SEED) -nox"

run: $(run_dirs)
$(run_dirs):
	cd $(@D); make clean; make ARG="-sim 999 -ainp PHake2010.pin -nox" mcmc
