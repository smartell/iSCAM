## Makefile for running examples
## Author: steven martell <martell.steve@gmail.com>
## targets:

## Macros
SUBDIRS = CC PRD QCI SOG WCVI AREA27 AREA2W
EXECUTABLE = ../../../dist/release/iscam

.PHONY: default $(SUBDIRS)

## Targets
default: $(SUBDIRS)

$(SUBDIRS):
	@echo $@
	$(MAKE) mcmc -C $@


dist: $(EXECUTABLE)
	cp $(EXECUTABLE) ./CC/
	cp $(EXECUTABLE) ./PRD/
	cp $(EXECUTABLE) ./QCI/
	cp $(EXECUTABLE) ./SOG/
	cp $(EXECUTABLE) ./WCVI/
	cp $(EXECUTABLE) ./AREA27/
	cp $(EXECUTABLE) ./AREA2W/

.PHONY: clean
clean_files := $(foreach dir,$(SUBDIRS),$(dir)/clean)

clean: $(clean_files)
$(clean_files):
	cd $(@D) && $(MAKE) clean

