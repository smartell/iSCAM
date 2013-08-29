## Make file for building iSCAM
## Author: Steven Martell
SRCS = msy.cpp OpMod.cpp Selex.cpp baranov.cpp

.PHONY: default opt clean

ifndef DISK
  DISK=../../dist
endif

default:
	admb -s -g iscam $(SRCS)
	cp iscam ${DISK}/debug

opt:
	admb -O iscam $(SRCS)
	cp iscam ${DISK}/release

clean:
	@rm -rvf ${DISK}/debug
	@rm -rvf ${DISK}/release
