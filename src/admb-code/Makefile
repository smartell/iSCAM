#Make file for building iSCAM

.PHONY: default optimize clean

ifndef DISK
  DISK=../../dist/
endif

default:
	mkdir -p ${DISK}/safe
	admb -s iscam

optimize:
	mkdir -p ${DISK}/opt
	admb iscam
	