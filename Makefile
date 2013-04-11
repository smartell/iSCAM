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
	make     --directory=src/admb-code 
	make opt --directory=src/admb-code 
	cp -r ./src/r-code/ ${DISK}/R/

clean:
	make --directory=src/admb-code --file=linux.mak clean
	rm -r dist