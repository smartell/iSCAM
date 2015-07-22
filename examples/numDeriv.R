#==================================================================================================================
#testing derivatives

### single fisheroes example

msycalc<-function(F){

	Ro<-1
	kappa<-5
	m<-0.2
	fe<-F

	ages<-1:70
	Winf<-3.165797
	k<-0.2728831
	to<- -3.2193539
	wa<-Winf*(1-exp(-k*(ages-to)))^3
	#wa<-c(1.0121,1.3858,1.725,2.0174,2.2609,2.4589,2.6171)
	
	f50<-2.0273497
	fsd<-0.4020049
	fa<-(1/(1+exp(-(ages-f50)/fsd)))*wa
	#fa<-c(0.07205217,0.48300621,0.91831884,0.99266382,0.99938078,0.99995933,1.00000000)
	
	v50<-2.9485810
	vsd<-0.5660189
	
	va<-(1/(1+exp(-(ages-v50)/vsd)))
	#va<-c(0.0041512,0.084884,0.60249,0.8,0.9,1,1)

	za <- m + fe * va
	sa <- exp(-za)
	oa <- 1-sa

	la<-NULL
	lx<-NULL
	la[1]<-1
	lx[1]<-1


	for(i in 2:(length(ages)))
	{
		la[i]<- la[i-1]* sa[i-1]
		lx[i]<- lx[i-1]* exp(-m)
	}
	LA<-la[length(ages)]
	la[length(ages)]<-la[length(ages)]/(oa[length(ages)])
	lx[length(ages)]<-lx[length(ages)]/(1-exp(-m))
		
	phiE<-sum(fa*lx)
	phie<-sum(fa*la)
	
	Re<-Ro*(kappa-phiE/phie)/(kappa-1)
	
	qa<- va*wa*oa/za

	phiq<-sum(la*qa)
	
	Ye<-fe*Re*phiq	

	#calculate derivatives

	#recruitment
	dla.dfe<-NULL

	dla.dfe[1]<-0

	for(i in 2:(length(ages)))
	{
		dla.dfe[i] <- dla.dfe[i-1]*sa[i-1] - la[i-1]*va[i-1]*sa[i-1]		
	}
	dLA.dfe<-dla.dfe[length(ages)]
	dla.dfe[length(ages)] <- dLA.dfe/oa[length(ages)] - LA*va[length(ages)]*sa[length(ages)]/(oa[length(ages)])^2


	dphie.dfe<- sum(fa*dla.dfe)

	dRe.dfe<- (dphie.dfe*Ro*phiE)/(phie^2*(kappa-1))

	#incidence functions

	dphiq.dfe<- sum(dla.dfe*qa+la*va^2*sa*wa/za-la*va^2*oa*wa/za^2)
	
	dYe.dfe<-Re*phiq+fe*Re*dphiq.dfe+fe*phiq*dRe.dfe

	#calculate second derivatives

	#recruitment

	ddla.ddfe<-NULL

	ddla.ddfe[1]<-0

	for(i in 2:(length(ages)))
	{
		ddla.ddfe[i] <- (ddla.ddfe[i-1] + la[i-1]*va[i-1]^2 - 2*dla.dfe[i-1]*va[i-1])*sa[i-1]		
	}
	ddLA.ddfe <- ddla.ddfe[length(ages)] 
	ddla.ddfe[length(ages)] <-  LA * (va[length(ages)])^2 * sa[length(ages)]    / ( oa[length(ages)]^2 )
							- 2 * dLA.dfe * va[length(ages)] * sa[length(ages)] / ( oa[length(ages)]^2 )
							+ ddLA.ddfe / oa[length(ages)] 
							+  2 * LA * (va[length(ages)])^2 * (sa[length(ages)])^2 /(oa[length(ages)]^3)
							
	ddphie.ddfe<- sum(fa*ddla.ddfe)

	ddRe.ddfe<- (phie* ddphie.ddfe - 2*dphie.dfe^2)*Ro*phiE/(phie^3*(kappa-1))

	#incidence functions

	ddphiq.ddfe<- sum(-la*va^3*sa*wa/za+2*dla.dfe*va^2*sa*wa/za-2*la*va^3*sa*wa/za^2
		+ddla.ddfe*qa-2*dla.dfe*va*qa/za+2*la*va^2*qa/za^2)


	ddYe.ddfe<- fe*phiq*ddRe.ddfe+2*fe*dphiq.dfe*dRe.dfe+2*phiq*dRe.dfe+fe*Re*ddphiq.ddfe+2*Re*dphiq.dfe

	return(list(F=fe,Ye=Ye, dYe.dfe=dYe.dfe,ddYe.ddfe=ddYe.ddfe, la=la, dla.dfe=dla.dfe, ddla.ddfe=ddla.ddfe))

}
msycalc(.3)

fe<-seq(0,2,0.0001)
calcYE<-NULL
calcdYE<-NULL
calcddYE<-NULL



for( a in 1:length(fe))
{
	result<-msycalc(fe[a])
	
	calcYE[a]<-result$Ye
	calcdYE[a]<-result$dYe.dfe
	calcddYE[a]<-result$ddYe.ddfe
	
}


par(mfrow=c(3,1))
plot(fe,calcYE,lwd=3,type="l")
plot(fe,calcdYE,lwd=3,type="l")
plot(fe,calcddYE,lwd=3,type="l")

which.max(calcYE)
which.min(abs(calcdYE))

calcYE[which.min(abs(calcdYE))]
fe[which.min(abs(calcdYE))]


Fstar<-0.1

#newton method
for(ct in 1:100)
{
	Fstar<- Fstar -  msycalc(Fstar)$dYe.dfe/msycalc(Fstar)$ddYe.ddfe
	print(Fstar)
}




numdYE<-NULL
calcdYE<-NULL
numddYE<-NULL
calcddYE<-NULL
calcdla<-NULL
numdla<-NULL
calcddla<-NULL
numddla<-NULL

fe<-seq(0,2,0.001)
h<-0.00001


for( a in 1:length(fe))
{
	result<-msycalc(fe[a])
	
	calcdYE[a]<-result$dYe.dfe
	
	numdYE[a]<-((msycalc(fe[a]+h)$Ye-msycalc(fe[a])$Ye)/h+ (msycalc(fe[a])$Ye-msycalc(fe[a]-h)$Ye)/h)/2

	calcddYE[a]<-result$ddYe.ddfe

	numddYE[a]<-(msycalc(fe[a]+2*h)$Ye - 2*msycalc(fe[a]+h)$Ye + msycalc(fe[a])$Ye)/h^2

	calcdla[a]<-sum(result$dla.dfe)

	numdla[a]<-sum(((msycalc(fe[a]+h)$la - msycalc(fe[a])$la)/h + (msycalc(fe[a])$la - msycalc(fe[a]-h)$la)/h)/2)

	calcddla[a]<-sum(result$ddla.ddfe)

	numddla[a]<-sum((msycalc(fe[a]+ h)$la - 2*msycalc(fe[a])$la + msycalc(fe[a]-h)$la)/h^2)
}

par(mfrow=c(2,2))
plot(fe,calcdYE, type="l", lwd=2)
lines(fe,numdYE, lwd=2, col="mediumorchid3" )
legend("topright", c("derivative", "numerical"),  col = c("black","mediumorchid3") , lwd =2, pch = NULL, bty = "n")
plot(fe,calcddYE, type="l", lwd=2)
lines(fe,numddYE,  lwd=2, col="mediumorchid3" )
plot(fe,calcdla, type="l", lwd=2)
lines(fe,numdla, lwd=2, col="mediumorchid3" )
plot(fe,calcddla, type="l", lwd=2)
lines(fe,numddla,  lwd=2, col="mediumorchid3" )


##wa<-c(1.0121,1.3858,1.725,2.0174,2.2609,2.4589,2.6171,2.6171,2.6171,2.6171)
##
##fx<-function(param){
##	Winf<-param[1]
##	k<-param[2]
##	to<-param[3]
##	wt<-Winf*(1-exp(-k*(1:7-to)))^3
##	wa<-c(1.0121,1.3858,1.725,2.0174,2.2609,2.4589,2.6171)
##	return(sum((wt-wa)^2))}
##
##fx(c(5,0.1,-1))
##
##(result<-optim(c(10,0.4,0), fx))
##
##Winf<-result$par[1]
##k<-result$par[2]
##to<-result$par[3]
##wt<-Winf*(1-exp(-k*(1:7-to)))^3
##
##plot(1:7,c(1.0121,1.3858,1.725,2.0174,2.2609,2.4589,2.6171))
##lines(1:7,wt, lwd=2)
##
##
##fx1<-function(param){
##	f50<-param[1]
##	fsd<-param[2]
##	ft<-(1/(1+exp(-(1:7-f50)/fsd)))
##	fa<-c(0.07205217,0.48300621,0.91831884,0.99266382,0.99938078,0.99995933,1.00000000)
##	return(sum((ft-fa)^2))}
##
##fx1(c(2,1))
##
##(result<-optim(c(2,1), fx1))
##
##fx2<-function(param){
##	v50<-param[1]
##	vsd<-param[2]
##	vt<-(1/(1+exp(-(1:7-v50)/vsd)))
##	va<-c(0.0041512,0.084884,0.60249,0.8,0.9,1,1)
##	return(sum((vt-va)^2))}
##
##fx2(c(2,1))
##
##(result<-optim(c(2.5,.5), fx2))
##v50<-result$par[1]
##vsd<-result$par[2]
##vt<-(1/(1+exp(-(1:7-v50)/vsd)))
##plot(1:7,c(0.0041512,0.084884,0.60249,0.8,0.9,1,1))
##lines(1:7,vt, lwd=2)



#==================================================================================================================
#testing derivatives

### two or more fisheries example

msycalc2<-function(F){

	Ro<-1
	kappa<-10
	m<-0.2
	fe<-F

	ages<-1:30
	Winf<-3.165797
	k<-0.2728831
	to<- -.7
	wa<-Winf*(1-exp(-k*(ages-to)))^3
	
	
	f50<-2.0273497
	fsd<-0.4020049
	fa<-(1/(1+exp(-(ages-f50)/fsd)))*wa
	
	
	v501<-2.9485810
	vsd1<-0.5660189
	va1<-(1/(1+exp(-(ages-v50)/vsd)))
	
	a2<-1
	b2<-3
	gamma2<-0.5
	va2<- (1/(1-gamma2))*(((1-gamma2)/gamma2)^gamma2)*(exp(a2*gamma2*(b2-ages))/(1+exp(a2*(b2-ages))))


	va<-cbind(va1,va2)

	#matplot(va)	

	za <- m
	for(n in 1:(ncol(fe))){
		za <- za + fe[,n] * va[,n]		
	}
	
	sa <- exp(-za)
	oa <- 1-sa

	la<-NULL
	lx<-NULL
	la[1]<-1
	lx[1]<-1

	for(i in 2:(length(ages)))
	{
		la[i]<- la[i-1]* sa[i-1]
		lx[i]<- lx[i-1]* exp(-m)
	}

	LA<-la[length(ages)]
	la[length(ages)]<-la[length(ages)]/(oa[length(ages)])
	lx[length(ages)]<-lx[length(ages)]/(1-exp(-m))
		
	phiE<-sum(fa*lx)
	phie<-sum(fa*la)
	
	Re<-Ro*(kappa-phiE/phie)/(kappa-1)
	
	qa<-matrix(NA,ncol=ncol(fe),nrow=length(ages))
	phiq<-matrix(NA,ncol=ncol(fe))
	Ye<-matrix(NA,ncol=ncol(fe))

	for(n in 1:(ncol(fe))){
		qa[,n]<- va[,n]*wa*oa/za
		phiq[,n]<-sum(la*qa[,n])
		Ye[,n]<-fe[,n]*Re*phiq[,n]	
	}

	#calculate derivatives

	#recruitment


	dla.dfe<- matrix(NA,ncol=ncol(fe),nrow=length(ages))
	dla.dfe[1,]<-0
	for(n in 1:(ncol(fe))){
		for(i in 2:(length(ages)))
		{
			dla.dfe[i,n] <- dla.dfe[i-1,n]*sa[i-1] - la[i-1]*va[i-1,n]*sa[i-1]		
		}
	}
	dLA.dfe<-dla.dfe[length(ages),]
	dla.dfe[length(ages),] <- dLA.dfe/oa[length(ages)] - LA*va[length(ages),]*sa[length(ages)]/(oa[length(ages)])^2
	
	

	dphie.dfe<- matrix(NA,ncol=ncol(fe),nrow=1)
	dRe.dfe<- matrix(NA,ncol=ncol(fe),nrow=1)


	for(n in 1:ncol(fe)){
		dphie.dfe[,n]<- sum(fa*dla.dfe[,n])
		dRe.dfe[,n]<- (dphie.dfe[,n]*Ro*phiE)/((kappa-1)*phie^2)
	}
		#incidence functions
	
	dqa.dfe<- list(matrix(NA,ncol=ncol(fe),nrow=length(ages)),matrix(NA,ncol=ncol(fe),nrow=length(ages)))
	dphiq.dfe<- matrix(NA,ncol=ncol(fe),nrow=ncol(fe))
	dYe.dfe<- NULL

	for(n in 1:ncol(fe)){
		for(nn in 1:ncol(fe)){

			if(n==nn){
				dqa.dfe[[n]][,nn]<- (va[,n]^2*wa)/za * (sa-oa/za)
				dphiq.dfe[n,nn]<- sum(qa[,n]*dla.dfe[,n]+dqa.dfe[[n]][,nn]*la)

			} else{
				dqa.dfe[[n]][,nn]<- (va[,n]*va[,nn]*wa)/za * (sa-oa/za)
				dphiq.dfe[n,nn]<- sum(qa[,n]*dla.dfe[,nn]+dqa.dfe[[n]][,nn]*la)
			}			
		
		}
		dYe.dfe[n]<-Re*phiq[,n]+fe[,n]*Re*dphiq.dfe[n,n]+fe[,n]*phiq[,n]*dRe.dfe[,n]
	}
	
	
	#calculate second derivatives

	#recruitment

	ddla.ddfe<-matrix(NA,ncol=ncol(fe),nrow=length(ages))
	ddphie.ddfe<-matrix(NA,ncol=ncol(fe),nrow=1)
	ddRe.ddfe<-matrix(NA,ncol=ncol(fe),nrow=1)

	ddla.ddfe[1,]<-0
	
	for(n in 1:(ncol(fe))){
		for(i in 2:(length(ages)))
		{
			ddla.ddfe[i,n] <- (ddla.ddfe[i-1,n] + la[i-1]*va[i-1,n]^2 - 2*dla.dfe[i-1,n]*va[i-1,n])*sa[i-1]		
		}
	}
	
	ddLA.ddfe <- ddla.ddfe[length(ages),] 

	ddla.ddfe[length(ages),] <-  (LA * (va[length(ages),])^2 * sa[length(ages)]    / ( oa[length(ages)]^2 )
								- 2 * dLA.dfe * va[length(ages),] * sa[length(ages)] / ( oa[length(ages)]^2 )
								+ ddLA.ddfe / oa[length(ages)] 
								+  2 * LA * (va[length(ages),])^2 * (sa[length(ages)])^2 /(oa[length(ages)]^3))
	

	for(n in 1:(ncol(fe))){	
	ddphie.ddfe[,n]<- sum(fa*ddla.ddfe[,n])

	ddRe.ddfe[,n]<- (phie* ddphie.ddfe[,n] - 2*dphie.dfe[,n]^2)*Ro*phiE/(phie^3*(kappa-1))
	}					
	

	

	#incidence functions
	ddqa.ddfe<- list(matrix(NA,ncol=ncol(fe),nrow=length(ages)),matrix(NA,ncol=ncol(fe),nrow=length(ages)))

	ddphiq.ddfe<-matrix(NA,ncol=ncol(fe),nrow=ncol(fe))

	for(n in 1:(ncol(fe))){
		for(nn in 1:(ncol(fe))){

			if(n==nn)
			{
				ddqa.ddfe[[n]][,nn]<- (va[,n]^3*sa*wa/za)*(-sa-2*sa/za+2*oa/za^2)

				
			}else{ 
				ddqa.ddfe[[n]][,nn]<- (va[,n]*va[,nn]^2*sa*wa/za)*(-sa-2*sa/za+2*oa/za^2)

			}
				ddphiq.ddfe[n,nn]<- sum(la*ddqa.ddfe[[n]][,nn]+2*dla.dfe[,n]*dqa.dfe[[n]][,nn]+qa[,n]*ddla.ddfe[,n])
		}
	}

	ddYe.ddfe<-matrix(NA,ncol=ncol(fe),nrow=ncol(fe))
	for(n in 1:(ncol(fe))){
		for(nn in 1:(ncol(fe))){
			if(n==nn)
			{
				ddYe.ddfe[n,nn]<- fe[,n]*phiq[,n]*ddRe.ddfe[,n]+2*fe[,n]*dphiq.dfe[n,nn]*dRe.dfe[,n]+2*phiq[,n]*dRe.dfe[,n]+fe[,n]*Re*ddphiq.ddfe[n,nn]+2*Re*dphiq.dfe[n,nn]
			}else{
				ddYe.ddfe[n,nn]<- fe[,n]*phiq[,n]*ddRe.ddfe[,nn]+2*fe[,n]*dphiq.dfe[n,nn]*dRe.dfe[,nn]+fe[,n]*Re*ddphiq.ddfe[n,nn]
			}
		}
	}
	return(list(F=fe,Ye=Ye, dYe.dfe=dYe.dfe,ddYe.ddfe=ddYe.ddfe, la=la, dla.dfe=dla.dfe, ddla.ddfe=ddla.ddfe))

}

F<-matrix(c(0.15,0.2),ncol=2)
msycalc2(F)


fe<-seq(0,1,0.01)

F<-expand.grid(seq(0,1,0.01),seq(0,1,0.01))

calcYE<-matrix(NA,ncol=ncol(F),nrow=nrow(F))
calcdYE<-matrix(NA,ncol=ncol(F),nrow=nrow(F))
calcddYE<-matrix(NA,ncol=ncol(F),nrow=nrow(F))



for( a in 1:nrow(F))
{
	result<-msycalc2(F[a,])
	
	calcYE[a,]<-result$Ye
	calcdYE[a,]<-result$dYe.dfe
	calcddYE[a,]<-diag(result$ddYe.ddfe)
	
}


par(mfrow=c(3,1))
plot(F[,1],calcYE[,1],lwd=3,type="l")
lines(F[,2],calcYE[,2],lwd=3,type="l",col="skyblue")
plot(F[,1],calcdYE[,1],lwd=3,type="l")
lines(F[,2],calcdYE[,2],lwd=3,type="l",col="skyblue")
plot(F[,1],calcddYE[,1],lwd=3,type="l")
lines(F[,2],calcddYE[,2],lwd=3,type="l",col="skyblue")



(cbind(F[,1],calcYE[,1]))[order(F[,1])]
plot((cbind(F[,1][order(F[,1])],calcYE[,1][order(F[,1])])),lwd=3,type="l")


#need a 3d plot in here

?order
