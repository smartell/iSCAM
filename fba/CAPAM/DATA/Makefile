## Makefile for running models
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
# SUBDIR = 1a 1b 1c 1d 2a 2b 2c 2d 3a 3b 3c 3d
SUBDIR = 1d 2d 3d
# SUBDIR = 1a 1b 1c 1d
# SUBDIR = 2a 2b 2c 2d
# SUBDIR = 3a 3b 3c 3d

TARGET = 
# SEED  = 98723 
# SEED  = 99123  
# SEED  = 12345   #Upward retrospective bias
SEED = 111411     #Downward retrospective bias
JOBS  = 8

.PHONY: default $(SUBDIR) mcmc


## Targets
default: $(SUBDIR)
$(SUBDIR):
	cd $@ && $(MAKE) $(TARGET)



.PHONY: clean
clean_files := $(foreach dir,$(SUBDIR),$(dir)/clean)
clean_dirs  := $(foreach dir,$(SUBDIR),$(dir)/cleanall)
collect_files := $(foreach dir,$(SUBDIR),$(dir)/collect)
run_dirs := $(foreach dir,$(SUBDIR),$dir/run)

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
	cd $(@D); make DAT=SIM.dat ARG="-sim $(SEED) -ainp Sim1_PHake2010.pin -noest" run;
	cd $(@D); cp Sim1_PHake2010.sim PHake2010.sim;
	# cd $(@D); make retro ARG="-nox -ainp PHake2010.pin"; touch tmp1;
	# cd $(@D); make mcmc  ARG="-nox -ainp PHake2010.pin"; touch tmp2;
	cd $(@D); make retro; touch tmp1;
	cd $(@D); make mcmc; touch tmp2;
	


.PHONY: collect
collect: $(collect_files)
$(collect_files):
	cd $(@D); rm allSims.Rdata
	cd $(@D); make collect
	cd $(@D); rm spawnbio.Rdata; make collectssb

