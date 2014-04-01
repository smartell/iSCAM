## Makefile for building distribtion folder for iscam
## TODO add verify target to run example models.
.PHONY: dist clean


ifndef DISK
  DISK=dist
endif


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