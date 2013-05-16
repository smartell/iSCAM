## Make file for building iSCAM
## Author: Steven Martell

.PHONY: default opt clean

ifndef DISK
  DISK=../../dist
endif

default:
	admb -s -g iscam
	cp iscam ${DISK}/debug

opt:
	admb iscam
	cp iscam ${DISK}/release

clean:
	@rm -rvf ${DISK}/debug
	@rm -rvf ${DISK}/release
