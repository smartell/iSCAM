## Author: Steve Martell
## Run the ADMB program
## Rules : est

ARG=-nox
DAT= WCVI.dat
REP= WCVIHerring2011q1

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