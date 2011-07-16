## Prior distributions for Iscam UG
require(hacks)
par(las=1,cex.lab=1.5)
## Beta prior
x = seq(0,1,by=0.001)
y = dbeta((x-0.2)/0.8,1.01,1.01)
y2 = dbeta((x-0.2)/0.8,20,10)
y3 = dbeta((x-0.2)/0.8,2.0,3.0)

plot(x,y2,type="n",xlim=c(0.2,1.0),
xlab="Steepness",ylab="Prior density")

xx = c(x,rev(x))
yy = c(y,rep(0,length(x)))
polygon(xx,yy,col=colr(2,0.15))
text(0.95,0.8,"p1=1.01\np2=1.01")
yy = c(y2,rep(0,length(x)))
polygon(xx,yy,col=colr(3,0.15))
text(0.75,3,"p1=20\np2=10")
yy = c(y3,rep(0,length(x)))
polygon(xx,yy,col=colr(4,0.15))
text(0.45,1.25,"p1=2.0\np2=3.0")
dev.copy2pdf(file="betaPriors.pdf")

