## Makefile for running examples
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
SUBDIR = DEMO/DATA             \
         ECODETECTIVE/DATA     \
         Herring4VWX/DATA      \
         LingcodSRA/DATA       \
         PacificHake/DATA      \
         BSPollock/DATA        \
         BCHerring/DATA

TARGET = make


## Targets
runall := $(foreach dir,$(SUBDIR),$(dir)/runall)
default: $(runall)	
$(runall):
	cd $(@D) && $(MAKE)




clean_files := $(foreach dir,$(SUBDIR),$(dir)/clean)

clean: $(clean_files)
$(clean_files):
	cd $(@D) && $(MAKE) clean

