## Makefile for running examples
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
SUBDIRS = CC PRD QCI SOG WCVI AREA27 AREA2W

.PHONY: default $(SUBDIRS)

## Targets
default: $(SUBDIRS)

$(SUBDIRS):
	@echo $@
	$(MAKE) mcmc -C $@



.PHONY: clean
clean_files := $(foreach dir,$(SUBDIRS),$(dir)/clean)

clean: $(clean_files)
$(clean_files):
	cd $(@D) && $(MAKE) clean

