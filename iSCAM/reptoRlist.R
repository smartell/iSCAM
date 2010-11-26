reptoRlist = function(fn)
{
	# Function to read a report file with file name fn
	# Example:             
	#	fn <-("~/admb/simple.rep")
	#	A <- reptoRlist(fn)
	#
	# Then the 'A' object contains a list structure
	# with all the elemements in the report file.
	# In the REPORT_SECTION of the AMDB template use 
	# the following format to output objects:
	#  	report<<"object \n"<<object<<endl;
	#
	# The part in quotations becomes the list name.
	# Created By Steven Martell
	options(warn=-1)
ifile=scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
idx=sapply(as.double(ifile),is.na)
vnam=ifile[idx] #list names
nv=length(vnam) #number of objects
A=list()
ir=0
for(i in 1:nv)
	{
		ir=match(vnam[i],ifile)
		if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
		dum=NA
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
		if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))

		if(is.numeric(dum))#Logical test to ensure dealing with numbers
		{
		A[[vnam[i]]]=dum
		}
	}
	options(warn=0)
return(A)
} 
 
read.fit<-function(file){ 
  # Function to read a basic AD Model Builder fit. 
  # Use for instance by: 
  #   simple.fit <- read.fit('c:/admb/examples/simple') 
  # 
  # Then the object 'simple.fit' is a list containing sub
  # 'names', 'est', 'std', 'cor', and 'cov' for all model
  # parameters and sdreport quantities.  
  #  
  ret<-list() 
  parfile<-as.numeric(scan(paste(file,'.par', sep=''),   
   what='', n=16, quiet=TRUE)[c(6,11,16)]) 
  ret$nopar<-as.integer(parfile[1]) 
  ret$nlogl<-parfile[2] 
  ret$maxgrad<-parfile[3] 
  file<-paste(file,'.cor', sep='') 
  lin<-readLines(file) 
  ret$npar<-length(lin)-2 
  ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
  sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
  ret$names<-unlist(lapply(sublin,function(x)x[2])) 
  ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
  ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
  ret$cor<-matrix(NA, ret$npar, ret$npar) 
  corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
  ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
  ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
  ret$cov<-ret$cor*(ret$std%o%ret$std) 
  return(ret) 
}

read.psv<-function(fn="SCAM.psv", nsamples=100000)
{
     filen <- file(fn, "rb")
	 nopar <- readBin(filen, what = integer(), n = 1)
	 mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
	 mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
	 close(filen)
	return(mcmc)
}


bubble.plot<-function(x,y,z,scale=1,xlab="Year",ylab="Age",log.scale=F,iyr=1:dim(z)[1],sample.size=F, leg=F,...)
{
	zo=z
	if(log.scale) zo=log(z+1)
	n=dim(z)[1]; ny=dim(z)[2]
	xo=outer(x,rep(1,length=length(y)))
	yo=t(outer(y,rep(1,length=length(x))))
	#scale = 2
	zo=apply(zo,2,"/",max(abs(zo)))*scale#*length(y)*scale;
	zo = z/max(abs(z))*scale
	nt=rowSums(z)
	matplot(xo,yo,type="n",xlab=xlab,ylab=ylab,...)
	#abline(v=pretty(x),lty=3,col="lightgrey")
	for(i in iyr){
		iclr=rep("honeydew",length=ny); iclr[zo[i,]<=0]="salmon"
		points(xo[i,1:ny],yo[i,1:ny],cex=abs(zo[i,]),pch=16,col=iclr)#, col="#FF000001")
		points(xo[i,1:ny],yo[i,1:ny],cex=abs(zo[i,]),pch=1,col="black")
		if(sample.size)text(x[i],ny+0.5,paste(nt[i]),cex=0.45)
		}
	pt = pretty(abs(z), n=2);
	#pt = seq(0,8,by=2)
	#pt=c(2,5,10)
	pts=pt/max(abs(z))*scale
	txt = paste(pt, "units")
	if(leg)legend("topleft", txt, pch=1, pt.cex=pts, inset=0.05,col="black",bty="n", cex=0.75)
}


require(MASS)
require( KernSmooth)
fried.egg=function(xx,yy,...)
{
	bw=25
	bwx=diff(extendrange(xx))/bw; bwy=diff(extendrange(yy))/bw
	#bwx=(max(xx)-min(xx))/bw
	#bwy=(max(yy)-min(yy))/bw
	est <- bkde2D(cbind(xx,yy),bandwidth=c(bwx,bwy),gridsize=c(81, 81))
	est$fhat=est$fhat/max(est$fhat)
	#plot(xx,yy,pch=".",col="dark grey",xlab=NA,ylab=NA,type="n")
	#text(max(xx),max(yy),labels="D",adj=c(1,1))
	lvs=c(0.05,0.25,0.75,0.95)
	maxct=max(lvs)
	nlvs=length(lvs)
	thelines=contourLines(est$x1,est$x2,est$fhat,levels=lvs)
	polygon(thelines[[nlvs-3]]$x,thelines[[nlvs-3]]$y,col="khaki",border="khaki",lwd=1)
	polygon(thelines[[nlvs-2]]$x,thelines[[nlvs-2]]$y,col="snow",border="snow1",lwd=2)
	polygon(thelines[[nlvs-1]]$x,thelines[[nlvs-1]]$y,col="yellow",border="yellow2",lwd=3)
	polygon(thelines[[nlvs]]$x,thelines[[nlvs]]$y,col="lightyellow",border="yellow",lwd=1)
	#contour(est$x1,est$x2,est$fhat,drawlabels=T,add=T,levels=lvs,lty=1,lwd=1,labcex= 0.7)
	#Add salt and pepper
	#xi=sample(1:length(xx),300)
	#points(xx[xi],yy[xi],pch=".",col=grey(0:10/10))
}