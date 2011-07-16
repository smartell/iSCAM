## Makefile for building distribtion folder for iscam
.PHONY: default clean


ifndef DISK
  DISK=dist/
endif

default:
	mkdir -p ${DISK}/debug
	mkdir -p ${DISK}/release
	
