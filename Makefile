## Makefile for building iscam
## RULES:
## 		all
##      clean

.PHONY: all clean

all:
	@make debug release --directory=./src/admb-code -j

clean:
	make --directory=./src/admb-code clean