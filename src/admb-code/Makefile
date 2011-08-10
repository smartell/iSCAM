#Make file for building iSCAM

.PHONY: default opt clean

ifndef DISK
  DISK=../../dist
endif

default:
	mkdir -p ${DISK}/safe
	admb -s iscam
	cp iscam ${DISK}/safe

opt:
	mkdir -p ${DISK}/opt
	admb iscam
	cp iscam ${DISK}/opt

clean:
	@rm -rvf ${DISK}/safe
	@rm -rvf ${DISK}/opt
	@rm -rvf ${DISK}
