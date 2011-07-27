## Author: Steve Martell
## Run the ADMB program
## Rules : est

ARG=-nox
DAT= Area2W.dat
REP= Area2WHerring2011

.PHONY = all est 

all: est 

est: $(REP).rep

$(REP).rep: $(REP).ctl $(DAT)
	chmod +x iscam;
	$(shell echo ./iscam -ind $(DAT) $(ARG))



clean:
	rm -f $(REP).rep $(REP).cor