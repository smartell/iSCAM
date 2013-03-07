## Makefile for running models
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
SUBDIR = 1a 1b 1c 2a 2b 2c 3a 3b 3c 
# SUBDIR = 1a 1b 1c
# SUBDIR = 2a 2b 2c
# SUBDIR = 3a 3b 3c

TARGET = 
SEED  = 991 
JOBS  = 8

.PHONY: default $(SUBDIR) mcmc


## Targets
default: $(SUBDIR)
$(SUBDIR):
	cd $@ && $(MAKE) $(TARGET)



.PHONY: clean
clean_files := $(foreach dir,$(SUBDIR),$(dir)/clean)
clean_dirs  := $(foreach dir,$(SUBDIR),$(dir)/cleanall)

clean: $(clean_files)
$(clean_files):
	cd $(@D) && $(MAKE) clean

cleanall: $(clean_dirs)

$(clean_dirs):
	cd $(@D) && $(MAKE) cleanall

.PHONY: sim
sim_dirs := $(foreach dir,$(SUBDIR),$(dir)/sim)
run_dirs := $(foreach dir,$(SUBDIR),$(dir)/run)



$(sim_dirs):
	cd $(@D);make data -j$(JOBS); make est -j$(JOBS); make collect; make retrostat; touch allDONE
	
sim: $(sim_dirs)

	#make clean; make data;
run: $(run_dirs)
$(run_dirs):
	cd $(@D); make ARG="-sim $(SEED) -ainp PHake2010.pin -nox";
