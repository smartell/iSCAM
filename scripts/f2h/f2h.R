# R-script for computing steepness (h) from Fmsy (f)
# Author: Steven Martell
# Date:   March 16,  2012
require(PBSmodelling)
require(ggplot2)

# CONSTANTS
  A    <- 55
  age  <- seq(1, A)
  fe   <- seq(0, 1.0, by=0.0025)

# PARAMETERS
theta <- c(
	msy  = 1, 
	fmsy = 0.215, 
	m    = 0.15, 
	linf = 130, 
	k    = 0.11,
	a    = 6.92e-6, 
	b    = 3.24, 
	a50  = 10.91, 
	k50  = 1.406, 
	l50  = 97.132, 
	g50  = 1/0.1666
	)

# PROCEDURES

# .fmsy2h calculates the fmsy to steepness transition
# Author: Steve Martell
.fmsy2h	<-
function(theta)
{
	
	with(as.list(theta), {
		fe	<- fmsy
		
		# Age-schedule information
		la	<- (linf*(1-exp(-k*age)))
		wa	<- a*la^b
		fa	<- plogis(age, a50, k50)*wa
		va	<- plogis(la, l50, g50)
		lx	<- exp(-m)^(age-1); lx[A]=lx[A]/(1-exp(-m))
		za	<- m+fe*va
		sa	<- (1-exp(-za))
		qa	<- va*sa/za
		t2	<- wa*va^2/za
		t3	<- exp(-za)-sa/za
		
		
		# Survivorship under fished conditions (fe>0)
		lz		<- rep(1, length=A)
		dlz.df	<- 0.
		dphie.df<- 0.
		dphiq.df<- t2[1]*t3[1]
		
		for(i in 2:A)
		{
			lz[i]	<- lz[i-1]*exp(-za[i-1])
			dlz.df	<- dlz.df *exp(-za[i-1]) - lz[i-1]*va[i-1]*exp(-za[i-1])
			if(i==A)
			{
				lz[i]	<- lz[i]/sa[i]
				dlz.df	<- dlz.df/sa[i] - lz[i-1]*exp(-za[i-1])*va[i]*exp(-za[i])/sa[i]^2
			}
			dphie.df	<- dphie.df + fa[i]*dlz.df
			dphiq.df	<- dphiq.df + wa[i]*qa[i]*dlz.df + lz[i]*t2[i]*t3[i]
		}
		
		# per recruit incidence functions
		phiE	<- sum(lx*fa)
		phie	<- sum(lz*fa)
		phiq	<- sum(lz*wa*qa)
		
		# steepness and unfished reference points (Beverton Holt)
		reck	<- phiE/phie - (fe*phiq*phiE/phie^2*dphie.df) / (phiq+fe*dphiq.df)
		re		<- msy/(fe*phiq)
		ro		<- re*(reck-1)/(reck-phiE/phie)
		bo		<- ro*phiE
		dre.df	<- ro/(reck-1)*phiE/phie^2*dphie.df
		h		<- reck/(4+reck)
		
		# named vector argument to return
		df		<- data.frame(cbind(age, lx, lz, fa, wa, va))
		theta	<- c(theta, h=h, reck=reck, ro=ro, bo=bo, dre.df=dre.df, df=df)
		return(theta)
	})
}

# .equilibrium returns steady state values for a given fishing rate fe
# Author. Steven Martell
.equilibrium	<-
function(theta, fe=0)
{
	fn	<- function(fe)
	{
	with(as.list(theta), {
		# Age-schedule information
		la	<- (linf*(1-exp(-k*age)))
		wa	<- a*la^b
		fa	<- plogis(age, a50, k50)*wa
		va	<- plogis(la, l50, g50)
		lx	<- exp(-m)^(age-1); lx[A]=lx[A]/(1-exp(-m))
		za	<- m+fe*va
		sa	<- (1-exp(-za))
		qa	<- va*sa/za
	
		# Survivorship under fished conditions (fe>0)
		lz		<- rep(1, length=A)
		for(i in 2:A)
		{
			lz[i]	<- lz[i-1]*exp(-za[i-1])
			if(i==A)
			{
				lz[i]	<- lz[i]/sa[i]
			}
		}
	
		# per recruit incidence functions
		phiE	<- sum(lx*fa)
		phie	<- sum(lz*fa)
		phiq	<- sum(lz*wa*qa)
	
		# Equilibrium values
		re	<- max(0, ro*(reck-phiE/phie)/(reck-1))
		ye	<- fe*re*phiq
		be	<- re*phie
		spr	<- phie/phiE
		
		# Return list
		phi	<- c(fe=fe, re=re, ye=ye, be=be, spr=spr, dep=be/bo, ypr=fe*phiq)
		return(phi)
	})
	}
	df	<- data.frame(t(as.matrix(sapply(fe,fn))))
	return(df)
}

# MAIN FUNCTION CALLS
phi	<- .fmsy2h(theta)
DF	<- .equilibrium(phi, fe)

dfe	<- data.frame(cbind(Age=phi$df.age,
						lx=phi$df.lx, 
						"Survivorship (F)"=phi$df.lz, 
						"Weight at age"=phi$df.wa, 
						"Fecundity at age"=phi$df.fa, 
						"Selectivity at age"=phi$df.va))
m.dfe <- melt(dfe,id.var=1)
nDF <- DF
names(nDF) <- c("fe", "Recruits", "Yield", "Biomass", "SPR", "Depletion", "YPR")
m.DF  <- melt(nDF, id.var=1); 

# Relationship between Fmsy & steepness
tmp <- function(theta) .fmsy2h(theta)[1:15]
xx	<- matrix(theta, ncol=length(theta), nrow=100, byrow=TRUE)
xx[,2] <- seq(0.01, 1.0, length=100)
colnames(xx) <- names(theta)
fh	<- matrix(unlist(apply(xx, 1, tmp)), ncol=15, nrow=100, byrow=TRUE)
colnames(fh) <- names(phi[1:15])
fh	<- as.data.frame(fh)

pfh <- ggplot(fh, aes(h, fmsy))+geom_line()+labs(x="Steepness (h)", y="F(MSY)") + BaseThemeX90(22)
pfh 
dev.copy2pdf(file="figFMSYSteepness.pdf")

s_fh <- subset(fh,select=c(fmsy,h,reck,ro,bo))
names(s_fh) <- c("Fmsy", "Steepness", "Recruitment Compensation", "Unfished Recruits", "Unfished Biomass")
m_fh <- melt(s_fh, id.var=2)
pffh <- ggplot(m_fh, aes(Steepness, value)) + labs(x="Steepness (h)",y="") + geom_line()
pffh + facet_wrap(~variable, scales="free") +BaseThemeX90(22)
dev.copy2pdf(file="figSteepVsOther.pdf")

# Prior density for steepness
require(MASS)
nn           = 2000
xx           = mvrnorm(nn,theta,0.05*0.05*diag(theta))
fx           = matrix(unlist(apply(xx, 1, tmp)), ncol=15, nrow=nn, byrow=TRUE)
colnames(fx) = names(phi[1:15])
fx           = data.frame(fx)
m_fx         = melt(fx)
s_fx         = subset(m_fx,subset=variable==c("fmsy","m","linf","k"))
gg           = ggplot(s_fx,aes(value)) +geom_density(adjust=1.5) +facet_wrap(~variable,scales="free") 

gg+BaseThemeX90(22); dev.copy2pdf(file="figThetaPrior.pdf")
qplot(h,data=fx,geom='density', adjust=2,xlim=c(0.2,1), xlab="Steepness", ylab="Prior density") +BaseThemeX90(22)
dev.copy2pdf(file="figSteepnessPrior.pdf")

# GRAPHICS
sr	<- ggplot(data=DF, aes(be, re))+labs(x="Spawning biomass", y="Recruitment")+opts(title="Stock Recruitment Relationship")
sr + geom_line(aes(size=ye/1.1, color=fe)) + BaseThemeX90(22); dev.copy2pdf(file="figStockRecruitII.pdf")

ye	<- ggplot(data=DF, aes(fe, ye))+labs(x="Fishing mortality", y="Yield")+opts(title="Equilibrium Yield vs. Fishing Mortality")
ye + geom_line(aes(size=dep, color=re)) + BaseThemeX90(22); dev.copy2pdf(file="figEquilYield.pdf")

lh	<- ggplot(subset(m.dfe,variable!="lx"),aes(age,value))+labs(x="Age (years)", y="Value")
lh +geom_line()+facet_wrap(~variable,scales="free") + BaseThemeX90(22)
dev.copy2pdf(file="figHalibutAgeSchedules.pdf")

q4	<- ggplot(subset(m.DF,subset=variable==c("Yield","YPR", "SPR", "Depletion")), aes(fe, value))+labs(x="Fishing mortality", y="")
q4 + geom_line()+facet_wrap(~variable, scales="free") + BaseThemeX90(22); dev.copy2pdf(file="figHalibutEquilYield.pdf")

pp	<- ggplot(data=DF, aes(dep, re/max(re))) +geom_line(colour="blue")+labs(x="Relative spawning biomass", y="Relative recruitment")
pp + geom_segment(aes(x=0.2,xend=0.2,y=0,yend=phi$h))+geom_segment(aes(x=0,xend=0.2,y=phi$h,yend=phi$h))+ BaseThemeX90(22)
dev.copy2pdf(file="figStockRecruitSteepness.pdf")
# GRAPHICS THEMES
require(hacks)
BaseThemeX90 <- function(base_size = 10) {

	structure(list(
		axis.line         = theme_blank(),
		axis.text.x       = theme_text(size = base_size * 0.8, lineheight = 0.9, colour = "grey80", hjust = 1, angle = 90),
		axis.text.y       = theme_text(size = base_size * 0.8, lineheight = 0.9, colour = "grey80", hjust = 1),
		axis.ticks        = theme_segment(colour = "grey50"),
		axis.title.x      = theme_text(size = base_size, colour = "grey50", face="bold"),
		axis.title.y      = theme_text(size = base_size, colour = "grey50", face="bold", angle = 90),
		axis.ticks.length = unit(0.15, "cm"),
		axis.ticks.margin = unit(0.1, "cm"),
		
		legend.background = theme_rect(colour= NA), 
		legend.key        = theme_rect(fill = "grey95", colour = "white"),
		legend.key.size   = unit(1.2, "lines"),
		legend.text       = theme_text(size = base_size * 0.7, colour = "grey50"),
		legend.title      = theme_text(size = base_size * 0.7, colour = "grey50", hjust = -1),
		legend.position   = "right",
		
		panel.background  = theme_rect(fill = "grey90", colour = NA), 
		panel.border      = theme_blank(), 
		panel.grid.major  = theme_line(colour = "white"),
		panel.grid.minor  = theme_line(colour = "grey95", size = 0.25),
		panel.margin      = unit(0.25, "lines"),
		
		strip.background  = theme_rect(fill = "grey80", colour = NA), 
		strip.label       = function(variable, value) value, 
		strip.text.x      = theme_text(size = base_size * 0.8),
		strip.text.y      = theme_text(size = base_size * 0.8, angle = -90),
		
		plot.background   = theme_rect(colour = NA),
		plot.title        = theme_text(size = base_size * 1.0, colour = "grey50"),
		plot.margin       = unit(c(1, 1, 0.5, 0.5), "lines")
), class = "options")

}


