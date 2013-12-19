## Make file for building iSCAM
## Author: Steven Martell
## ADMB=~/admb-trunk/admb
.PHONY: default opt clean

ifndef DISK
  DISK=../../dist
endif

default:
	admb -s iscam logistic_normal.cpp
	cp iscam ${DISK}/debug

opt:
	admb -f iscam logistic_normal.cpp
	cp iscam ${DISK}/release

clean:
	@rm -rvf ${DISK}/debug
	@rm -rvf ${DISK}/release
