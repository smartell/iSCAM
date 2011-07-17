## Makefile for building distribtion folder for iscam
.PHONY: dist clean


ifndef DISK
  DISK=dist/
endif


dist:
	mkdir -p ${DISK}/debug
	mkdir -p ${DISK}/R
	mkdir -p ${DISK}/release
	make --directory=src/admb-code --file=linux.mak 
	make --directory=src/admb-code --file=linux.mak opt
	cp ./src/r-code/iSCAM.R ${DISK}/R
	cp ./src/r-code/iSCAMViewTracker.txt ${DISK}/R
	cp ./src/r-code/iSCAMWin.txt ${DISK}/R
	cp ./src/r-code/read.admb.R ${DISK}/R
	cp ./src/r-code/iScamLogo.gif ${DISK}/R

clean:
	make --directory=src/admb-code --file=linux.mak clean