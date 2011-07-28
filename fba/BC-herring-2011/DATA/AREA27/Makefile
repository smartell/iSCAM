## Author: Steve Martell
## Run the ADMB program
## Rules : est

ARG=-nox
DAT= Area27.dat
REP= Area27Herring2011

.PHONY = all est 

all: est $(REP).mcmc

est: $(REP).rep 

$(REP).rep: $(REP).ctl $(DAT) 
	chmod +x iscam;
	$(shell echo ./iscam -ind $(DAT) $(ARG))

$(REP).mcmc: $(REP).ctl $(REP).psv
	$(shell echo ./iscam -ind $(DAT) -mceval)

clean:
	rm -f $(REP).rep $(REP).cor $(REP).mcmc