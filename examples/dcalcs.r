#|-----------------------------------------------------------------------------------|#
#| NAME:dcalcs.r                                                                     |#
#| AUTHOR: Catarina Wor 										                     |#
#| PURPOSE: Keep track of functions, derivatives and second derivatives in MSY.hpp	 |#
#| DATE: JUL 2nd 2015	 															 |#
#|-----------------------------------------------------------------------------------|#


#Population Parameters
ro<-1
h<-0.75
m<-0.3
vbk<-0.2

age<-1:20
nage<-length(age)

#fisheries parameters

ak<-c(0.4,0.6)
fe<-c(0.222)
ah1<-4
ah2<-7

#newton method parameters
dh<-0.000001


#derived quantities

reck<-4*h/(1-h)
wa<-(1-exp(-vbk*age))^3
fa<-c(rep(0,sum(wa<=0.2)),wa[wa>0.2])
va1<-1/(1+exp(-(age-ah1)/(0.1*ah1)))
va2<-1/(1+exp(-(age-ah2)/(0.1*ah2)))

za<-m+fe*va1+fe*va2
sa<-exp(-za)
oa<-1-sa


lx<-exp(-m)^(age-1)
lx[nage]<-(exp(-m)^(nage-1))/(1-exp(-m)^nage)

lz[age[1]]<-1
for(i in 2:nage){
	lz[i]<-lz[i-1]*exp(-za[i-1])
}
lz[nage]<-lz[nage-1]*(exp(-za[nage-1]))/(1-exp(-za[nage]))



#expressions

ex_reck<- expression(4*h/(1-h))

ex_lx<-expression(exp(-m)^(age-1))
ex_lx_nage<-expression(exp(-m)^(age-1))


ex_za<-expression(m+fe*va1+fe*va2)


#=========================================================================================================================
# First derivative of lz

#age 1
ex_lz_sage<-expression(1)
D(ex_lz_sage,"fe")
#0

#age2
eq_lz2<-paste(ex_lz_sage,expression("*exp(-("),ex_za,expression("))"))
ex_lz2<- expression(1*exp(-( m + fe * va1 + fe * va2 )))
D(ex_lz2,"fe")
#r output
#dz2.df=-(exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1))

#final
#dz2.df=-(exp(-za1) * (va1a1 + va2a1))

#age3
ex_lz3<- expression( exp(-(m + fe *( va1a1 + va2a1))) *exp(-( m + fe * (va1a2 + va2a2))) )
D(ex_lz3,"fe")
#-(	exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + 
#    va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * 
#    (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + va2a2)))   )

#step1
#-(exp(-za1) * (exp(-za2) * (va1a2 + va2a2)) + exp(-za1) *  (va1a1 + va2a1) * exp(-za2))

#step2
#-(exp(-za1) * exp(-za2) * (va1a2 + va2a2) + exp(-za1) *  (va1a1 + va2a1) * exp(-za2))

#final
#dz3.df = -(lz3 * (va1a2 + va2a2) + (- dz2.df ) * exp(-za2))

#age4
ex_lz4<- expression( exp(-(m + fe *( va1a1 + va2a1))) *exp(-( m + fe * (va1a2 + va2a2))) * exp(-(m + fe *( va1a3 + va2a3))))
D(ex_lz4,"fe")
#r output
#-(exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
#    va2a2))) * (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + va2a3)) + 
#    (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + 
#        va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + 
#        va2a1))) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + 
#        va2a2)))) * exp(-(m + fe * (va1a3 + va2a3))))

#step1
#-(exp(-za1) * exp(-za2) * (exp(-za3) * (va1a3 + va2a3)) 
#	+ (exp(-za1) * (exp(-za2) * (va1a2 + va2a2)) 
#    + exp(-za1) * (va1a1 + va2a1) * exp(-za2)) * exp(-za3))


#step2
#dz4.df = -(exp(-za1) * exp(-za2) * exp(-za3) * (va1a3 + va2a3) 
#    + (exp(-za1) * exp(-za2) * (va1a2 + va2a2) + exp(-za1) * exp(-za2) * (va1a1 + va2a1) ) * exp(-za3)
#    )

#final
#dz4.df = -(lz4 * (va1a3 + va2a3) + (- dz3.df ) * exp(-za3))

#generalize for any age a
#dza.df = -(lza * sum(va-1) + (- dza-1.df ) * exp(-za-1))

#compare with MSY.hpp -  o parentesis não bate
#dza.df =  exp(-za-1) *  dza-1.df    -lza        * sum(va-1) 
#          sa(h)(j-1) * (dlz(k)(j-1) -lz(h)(j-1) *m_Va(h)(k)(j-1))


#Still need to calculate + group


