## Makefile for running examples
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
SUBDIR = Hake PacificHake
TARGET = 

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

