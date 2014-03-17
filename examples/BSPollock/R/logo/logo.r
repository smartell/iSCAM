#R-code for the BSCAM&D logo
#BSCAM&D: Bayesian Statistical Catch Age Model & Diagnostics
graphics.off()
quartz("SCAM logo", width=6, height=3)

logo <- function()
{
	set.seed(42008)
	n=5000; f=20; j=200; sd=0.4
	px = (1:n)/(2*n)
	ht = -.05 + (.05+.05)* sin((1:n-1)/(n-1)*12*pi)
	i=rnorm(n, 0.3, sd*px)
	s=rnorm(n, 0.2, sd*px)
	c=rnorm(n, 0.1, sd*px)
	a=rnorm(n, 0.0, sd*px)
	m=rnorm(n, -.1, sd*px)

	shift=0.1*n
	x = cbind(1:n, 1:n+shift, 1:n+2*shift, 1:n+3*shift, 1:n+4*shift)
	y = cbind(i, s, c, a, m)

	par(mfcol=c(1, 1), 
		xaxt="n", yaxt="n", mar=rep(0.25, length=4), oma=rep(0, length=4))

	matplot(x, y+ht, lty=1, type="l",col=1:4,  ylim=c(-0.75, 0.75), bty="n", ylab="", xlab="", xlim=c(-10, 0.5*n))

	id <- c(1:3, 5)
	text(x[1, id]-0.01*n, y[1, id]-.12, c("i", "S", "C", "M"), srt=0, col=1:3, cex=2, font=2)
	text(x[1, 4]-0.01*n, y[1, 4]-.12, c("A"), srt=180, col=4, cex=2, font=2)
}

# Encapsulate postscript file for conversion to png
logo()
dev.copy2eps(file="iscamLogo.eps", width=6, height=3)
dev.off()

#png(file="iScamLogo.png", bg="transparent")
#png(file="iScamLogo.png", bg="transparent",width=100, height=50, pointsize=6)
# png(file="iscamLogo.png", bg="transparent", width=440, height=220, res=220)
# par(mfcol=c(1, 1), 
# 	xaxt="n", yaxt="n", mar=rep(0.25, length=4), oma=rep(0, length=4))
# logo()
# 
# dev.off()
system("convert iscamLogo.eps -strip iscamLogo.png")
system("convert iscamLogo.eps -resize 300x150 -strip iscamLogoSmall.png")
system("convert iscamLogo.eps -resize 300x150 -strip iscamLogo.gif")

