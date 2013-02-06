## Makefile for running models
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
SUBDIR = 1a 1b
##SUBDIR = SOKeg
TARGET = 
SEED  = 999 

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

sim: $(sim_dirs)
$(sim_dirs):
	cd $(@D); make clean; make data; make est -j8; make collect #; make ARG="-sim $(SEED) -nox"
