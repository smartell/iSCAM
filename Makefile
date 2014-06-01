## Makefile for building iscam
## RULES:
## 		all

.PHONY: dist clean

ifndef DISK
  DISK=dist
endif

ifeq ($(mode),release)
	opt = TRUE
	dir = ${DISK}/release/
else
	mode = debug
	opt  = FALSE
	dir = ${DISK}/debug/
endif

all:
	mkdir -p $(dir)
	make  --directory=src/admb-code OPT=$(opt) -j
	cp    ./src/admb-code/iscam $(dir)


dist:
	mkdir -p ${DISK}/debug
	mkdir -p ${DISK}/R
	mkdir -p ${DISK}/release
	make     --directory=src/admb-code clean
	make     --directory=src/admb-code OPT=TRUE
	cp    ./src/admb-code/iscam ${DISK}/release/
	make     --directory=src/admb-code clean
	make     --directory=src/admb-code 
	cp    ./src/admb-code/iscam ${DISK}/debug/
	cp -r ./src/R/ ${DISK}/R/
	# cp -r ./src/r-code/ ${DISK}/R/

clean:
	#make --directory=src/admb-code --file=linux.mak clean
	make --directory=src/admb-code clean
	rm -r dist