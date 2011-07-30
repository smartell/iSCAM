## Makefile for running examples
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
subdirs = Hake PacificHake
makeopt = default

.PHONY: default $(subdirs)

## Targets
default: $(subdirs)
$(subdirs):
	@echo $@
	cd $@ && $(MAKE) all



.PHONY: clean
clean_files := $(foreach dir,$(subdirs),$(dir)/clean)

clean: $(clean_files)
$(clean_files):
	cd $(@D) && $(MAKE) clean

