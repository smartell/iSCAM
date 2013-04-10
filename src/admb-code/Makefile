## Make file for building iSCAM
## Author: Steven Martell
ADMB=~/admb-trunk/admb
.PHONY: default opt clean

ifndef DISK
  DISK=../../dist
endif

default:
	~/admb-trunk/admb -s -g iscam
	cp iscam ${DISK}/debug

opt:
	~/admb-trunk/admb iscam
	cp iscam ${DISK}/release

clean:
	@rm -rvf ${DISK}/debug
	@rm -rvf ${DISK}/release
