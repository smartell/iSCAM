# An R-script for reading the wobblesq.rep file and creating objects for halibut simulation

source("../read.admb.R")
fn	<- "../../MISC/wobblesq.rep"
A	<- read.rep(fn)


# Initial recruitment and init_log_rec_devs
age =2:30
for(h in 1:2)
{
	m = A$M[1, h]
	mj = m*(age-1)
	N = A$InitN[, h]
	ddotR = mean(log(N)+mj)
	wj = log(N)+mj - ddotR
	
	cat(ddotR, "\n")
	cat(wj, "\n")
}
