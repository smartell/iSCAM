## Author: Steve Martell
## Run the ADMB program
## Rules : est

ARG=-nox
DAT= QCI.dat
REP= QCIHerring2011q1

.PHONY = all est 

all: est 

est: $(REP).rep

$(REP).rep: $(REP).ctl $(DAT)
	chmod +x iscam;
	$(shell echo ./iscam -ind $(DAT) $(ARG))



clean:
	rm -f $(REP).rep $(REP).cor