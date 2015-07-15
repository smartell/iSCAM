#|-----------------------------------------------------------------------------------|#
#| NAME:dcalcs.r                                                                     |#
#| AUTHOR: Catarina Wor 										                     |#
#| PURPOSE: Keep track of functions, derivatives and second derivatives in MSY.hpp	 |#
#| DATE: JUL 2nd 2015	 															 |#
#|-----------------------------------------------------------------------------------|#



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
ex_lz2<- expression(1*exp(-( m + fe * (va1a1  + va2a1 ))))
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
-(exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
    va2a2))) * (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + va2a3)) + 
    (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + 
        va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + 
        va2a1))) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + 
        va2a2)))) * exp(-(m + fe * (va1a3 + va2a3))))

#step1
-(exp(-za1) * exp(-za2) * (exp(-za3) * (va1a3 + va2a3)) 
	+ (exp(-za1) * (exp(-za2) * (va1a2 + va2a2)) 
    + exp(-za1) * (va1a1 + va2a1) * exp(-za2)) * exp(-za3))


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



#plus group - Assuming age4 is max
ex_lzplus<- expression((exp(-(m + fe *( va1a1 + va2a1))) *exp(-( m + fe * (va1a2 + va2a2))) * exp(-(m + fe *( va1a3 + va2a3))))/(1-exp(-(m + fe *( va1a4 + va2a4))))  )
D(ex_lzplus,"fe")

-((exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
    va2a2))) * (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + va2a3)) + 
    (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + 
        va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + 
        va2a1))) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + 
        va2a2)))) * exp(-(m + fe * (va1a3 + va2a3))))/(1 - exp(-(m + 
    fe * (va1a4 + va2a4)))) + (exp(-(m + fe * (va1a1 + va2a1))) * 
    exp(-(m + fe * (va1a2 + va2a2))) * exp(-(m + fe * (va1a3 + 
    va2a3)))) * (exp(-(m + fe * (va1a4 + va2a4))) * (va1a4 + 
    va2a4))/(1 - exp(-(m + fe * (va1a4 + va2a4))))^2)

#step2



-((exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3)) 
	+ (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) 
	+ exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3))/(1 - exp(-z4)) + (exp(-z1) * 
    exp(-z2) * exp(-z3)) * exp(-z4) * (va1a4 + va2a4)/(1 - exp(-z4))^2)


-(
	(exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3)) 
	+ (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) 
	+ exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3))
	/(1 - exp(-z4)) 
	+ exp(-z1) * exp(-z2) * exp(-z3) * exp(-z4) * (va1a4 + va2a4)
	/(1 - exp(-z4))^2)


-(
	(lz4* (va1a3 + va2a3) -dz3.df * exp(-z3))/(1 - exp(-z4)) 	
	+ exp(-z1) * exp(-z2) * exp(-z3) * exp(-z4) * (va1a4 + va2a4)
	/(1 - exp(-z4))^2)


-(
	-dlz4/(1 - exp(-z4)) 	
	+ exp(-z1) * exp(-z2) * exp(-z3) * exp(-z4) * (va1a4 + va2a4)
	/(1 - exp(-z4))^2)

-(
	-dlz4/(1 - exp(-z4)) 	
	+ lz4 * exp(-z4) * (va1a4 + va2a4)
	/(1 - exp(-z4))^2)


	dlz4/(1 - exp(-z4)) 	
	- lz4 * exp(-z4) * (va1a4 + va2a4)
	/(1 - exp(-z4))^2

#compare to msy.hpp
dlz4/(1 - exp(-z4)) 	
	- lz3 * exp(-z3)* exp(-z4) * (va1a4 + va2a4)
	/(1 - exp(-z4))^2

dlz(k)(j)/oa(h)(j) - 
lz(h)(j-1)*sa(h)(j-1)*m_Va(h)(k)(j)*sa(h)(j)
						             /square(oa(h)(j));



 #=========================================================================================================================
# Second derivative of lz

#age 2
dz2.df=expression(-(exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1)))
D(dz2.df,"fe")

#exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (va1a1 + va2a1)

#exp(-z1) * (va1a1 + va2a1) * (va1a1 + va2a1)


#age 3
#dz3.df=expression(-(exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + 
#    va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * 
#    (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + va2a2)))))
#
#D(dz3.df,"fe")
#
exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (exp(-(m + 
    fe * (va1a2 + va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * 
    (va1a1 + va2a1))) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-(m + 
    fe * (va1a2 + va2a2))) + (exp(-(m + fe * (va1a1 + va2a1))) * 
    (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + va2a2) * (va1a2 + 
        va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + 
    va2a1) * (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + va2a2)))
#
#exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + va2a2)) 
#+ exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + va2a2)))
#+ (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + va2a2) * (va1a2 + va2a2)) 
#	+ exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + va2a2)))
# 
#
#exp(-z1) * (va1a1 + va2a1) * exp(-z2) * (va1a2 + va2a2) 
#+ exp(-z1) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-z2)
#+ exp(-z1) * exp(-z2) * (va1a2 + va2a2) * (va1a2 + va2a2) 
#+ exp(-z1) * (va1a1 + va2a1) * exp(-z2) * (va1a2 + va2a2)
#
#
#
#exp(-z1) * (va1a1 + va2a1) * (exp(-z2) * (va1a2 + va2a2)) 
#+ exp(-z1) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-z2) 
#+ exp(-z1) * (exp(-z2) * (va1a2 + va2a2) * (va1a2 + va2a2)) 
#+ exp(-z1) * (va1a1 +  va2a1) * (exp(-z2) * (va1a2 + va2a2))
#
#
#d2lz2 * exp(-z2) 
#+ exp(-z1) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-z2) 
#+ exp(-z1) * exp(-z2) * (va1a2 + va2a2) * (va1a2 + va2a2) 
#+ exp(-z1) * (va1a1 +  va2a1) * exp(-z2) * (va1a2 + va2a2)
#
#
#
#exp(-z2) *(d2lz2
#+ exp(-z1) * (va1a1 + va2a1)^2 
#+ exp(-z1) * (va1a2 + va2a2)^2 
#+ exp(-z1) * (va1a1 + va2a1) * (va1a2 + va2a2))
#
#
#
#
#d2lz3=exp(-z2) *(d2lz2 + exp(-z1) * ((va1a1 + va2a1)^2 + (va1a2 + va2a2)^2 + (va1a1 + va2a1) * (va1a2 + va2a2)))
#
#d2lz3= exp(-z2) *d2lz2 
#	+ exp(-z2) * exp(-z1) * ((va1a1 + va2a1)^2 + (va1a2 + va2a2)^2 + (va1a1 + va2a1) * (va1a2 + va2a2))
#
#d2lz3= exp(-z2) *d2lz2 
#	+ lz3 * ((va1a1 + va2a1)^2 + (va1a2 + va2a2)^2 + (va1a1 + va2a1) * (va1a2 + va2a2))
#
#
#sa(h)(j-1)*(d2lz(k)(j-1)+lz(h)(j-1)*square(m_Va(h)(k)(j-1)));

#age 4

#dz4.df=expression(-(exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
#   va2a2))) * (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + va2a3)) + 
#   (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + 
#       va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + 
#        va2a1))) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + 
#       va2a2)))) * exp(-(m + fe * (va1a3 + va2a3)))))
#
#D(dz4.df,"fe")
#
#
#(exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + 
#    )    va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * 
#    (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + va2a2)))) * (exp(-(m + 
#    fe * (va1a3 + va2a3))) * (va1a3 + va2a3)) + (exp(-(m + fe * 
#    (va1a1 + va2a1))) * (va1a1 + va2a1) * (exp(-(m + fe * (va1a2 + 
#    va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * 
#    (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + 
#    va2a2))) + (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + 
#    fe * (va1a2 + va2a2))) * (va1a2 + va2a2) * (va1a2 + va2a2)) + 
#    exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (exp(-(m + 
#        fe * (va1a2 + va2a2))) * (va1a2 + va2a2)))) * exp(-(m + 
#    fe * (va1a3 + va2a3))) + (exp(-(m + fe * (va1a1 + va2a1))) * 
#    exp(-(m + fe * (va1a2 + va2a2))) * (exp(-(m + fe * (va1a3 + 
#    va2a3))) * (va1a3 + va2a3) * (va1a3 + va2a3)) + (exp(-(m + 
#    fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + va2a2))) * 
#    (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + 
#    va2a1) * exp(-(m + fe * (va1a2 + va2a2)))) * (exp(-(m + fe * 
#    (va1a3 + va2a3))) * (va1a3 + va2a3)))

##  (exp(-z1) * exp(-z2) * (va1a2 + va2a2) + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3) * (va1a3 + va2a3) 
##+ (exp(-z1) * (va1a1 + va2a1) * exp(-z2) * (va1a2 + va2a2) 
##+ exp(-z1) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-z2)
##+ (exp(-z1) * exp(-z2) * (va1a2 + va2a2) * (va1a2 + va2a2) 
##+ exp(-z1) * (va1a1 + va2a1) * (exp(-z2) * (va1a2 + va2a2)))) * exp(-z3) 
##+ exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3) * (va1a3 + va2a3)) 
##+ (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3) * (va1a3 + va2a3)

#
#
#  (exp(-z1) * exp(-z2) * (va1a2 + va2a2) + exp(-z1) * (va1a1 + va2a1) * exp(-z2))
#   * exp(-z3) * (va1a3 + va2a3) 
#+ d2lz3 * exp(-z3) 
#+ exp(-z1) * exp(-z2) * exp(-z3) * (va1a3 + va2a3) * (va1a3 + va2a3) 
#+ (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3) * (va1a3 + va2a3)
#
#
# exp(-z3) * (
# (exp(-z1) * exp(-z2) * (va1a2 + va2a2) + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * (va1a3 + va2a3) 
#+ d2lz3 
#+ exp(-z1) * exp(-z2)  * (va1a3 + va2a3) * (va1a3 + va2a3) 
#+ (exp(-z1) * exp(-z2) * (va1a2 + va2a2) + exp(-z1) * (va1a1 + va2a1) * exp(-z2))  * (va1a3 + va2a3))
#
#
#exp(-z3) * ( d2lz3 
#+ (exp(-z1) * exp(-z2) * (va1a2 + va2a2) 
#+  exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * (va1a3 + va2a3) 
#+  exp(-z1) * exp(-z2)  * (va1a3 + va2a3) * (va1a3 + va2a3) 
#+ (exp(-z1) * exp(-z2) * (va1a2 + va2a2) + exp(-z1) * (va1a1 + va2a1) * exp(-z2))  * (va1a3 + va2a3))
#
#
#exp(-z3) * ( d2lz3 + exp(-z1) * exp(-z2) *(
#+ ( (va1a2 + va2a2) +   (va1a1 + va2a1) ) * (va1a3 + va2a3) 
#+  (va1a3 + va2a3) * (va1a3 + va2a3) 
#+ (  (va1a2 + va2a2) +  (va1a1 + va2a1) )  * (va1a3 + va2a3)))
#
#
#exp(-z3) * d2lz3 + exp(-z1) * exp(-z2) * exp(-z3) *(
#+ ( (va1a2 + va2a2) +   (va1a1 + va2a1) ) * (va1a3 + va2a3) 
#+  (va1a3 + va2a3) * (va1a3 + va2a3) 
#+ (  (va1a2 + va2a2) +  (va1a1 + va2a1) )  * (va1a3 + va2a3)))
#
#
#exp(-z3) * d2lz3 + lz4 *(
#+ ( (va1a2 + va2a2) +   (va1a1 + va2a1) ) * (va1a3 + va2a3) 
#+  (va1a3 + va2a3) * (va1a3 + va2a3) 
#+ (  (va1a2 + va2a2) +  (va1a1 + va2a1) )  * (va1a3 + va2a3)))

#same as msy.cpp

#age plus group
#assume age 4 is plus group

ex_d2lzplus<- expression(-((exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
    va2a2))) * (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + va2a3)) + 
    (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + 
        va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + 
        va2a1))) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + 
        va2a2)))) * exp(-(m + fe * (va1a3 + va2a3))))/(1 - exp(-(m + 
    fe * (va1a4 + va2a4)))) + (exp(-(m + fe * (va1a1 + va2a1))) * 
    exp(-(m + fe * (va1a2 + va2a2))) * exp(-(m + fe * (va1a3 + 
    va2a3)))) * (exp(-(m + fe * (va1a4 + va2a4))) * (va1a4 + 
    va2a4))/(1 - exp(-(m + fe * (va1a4 + va2a4))))^2))

D(ex_d2lzplus,"fe")

#result
((exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
    va2a2))) * exp(-(m + fe * (va1a3 + va2a3)))) * (exp(-(m + 
    fe * (va1a4 + va2a4))) * (va1a4 + va2a4) * (va1a4 + va2a4)) + 
    (exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
        va2a2))) * (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + 
        va2a3)) + (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + 
        fe * (va1a2 + va2a2))) * (va1a2 + va2a2)) + exp(-(m + 
        fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * exp(-(m + 
        fe * (va1a2 + va2a2)))) * exp(-(m + fe * (va1a3 + va2a3)))) * 
        (exp(-(m + fe * (va1a4 + va2a4))) * (va1a4 + va2a4)))/(1 - 
    exp(-(m + fe * (va1a4 + va2a4))))^2 + (exp(-(m + fe * (va1a1 + 
    va2a1))) * exp(-(m + fe * (va1a2 + va2a2))) * exp(-(m + fe * 
    (va1a3 + va2a3)))) * (exp(-(m + fe * (va1a4 + va2a4))) * 
    (va1a4 + va2a4)) * (2 * (exp(-(m + fe * (va1a4 + va2a4))) * 
    (va1a4 + va2a4) * (1 - exp(-(m + fe * (va1a4 + va2a4))))))/((1 - 
    exp(-(m + fe * (va1a4 + va2a4))))^2)^2 + (((exp(-(m + fe * 
    (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + 
    va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * 
    exp(-(m + fe * (va1a2 + va2a2)))) * (exp(-(m + fe * (va1a3 + 
    va2a3))) * (va1a3 + va2a3)) + (exp(-(m + fe * (va1a1 + va2a1))) * 
    (va1a1 + va2a1) * (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + 
    va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * 
    (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + va2a2))) + (exp(-(m + 
    fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + va2a2))) * 
    (va1a2 + va2a2) * (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + 
    va2a1))) * (va1a1 + va2a1) * (exp(-(m + fe * (va1a2 + va2a2))) * 
    (va1a2 + va2a2)))) * exp(-(m + fe * (va1a3 + va2a3))) + (exp(-(m + 
    fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + va2a2))) * 
    (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + va2a3) * (va1a3 + 
        va2a3)) + (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + 
    fe * (va1a2 + va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * 
    (va1a1 + va2a1))) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + 
    va2a2)))) * (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + 
    va2a3))))/(1 - exp(-(m + fe * (va1a4 + va2a4)))) + (exp(-(m + 
    fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + va2a2))) * 
    (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + va2a3)) + (exp(-(m + 
    fe * (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + va2a2))) * 
    (va1a2 + va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + 
    va2a1) * exp(-(m + fe * (va1a2 + va2a2)))) * exp(-(m + fe * 
    (va1a3 + va2a3)))) * (exp(-(m + fe * (va1a4 + va2a4))) * 
    (va1a4 + va2a4))/(1 - exp(-(m + fe * (va1a4 + va2a4))))^2)

#replace f and m for z
(exp(-z1) * exp(-z2) * exp(-z3) * exp(-z4) * (va1a4 + va2a4) * (va1a4 + va2a4) 
+ (exp(-z1) * exp(-z2) * exp(-z3) * (va1a3 + va2a3) + (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) 
+ exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3)) * (exp(-z4) * (va1a4 + va2a4)))/(1 - exp(-z4))^2
+ (exp(-z1) * exp(-z2) * exp(-z3)) * (exp(-z4) * (va1a4 + va2a4)) * (2 * (exp(-z4) * (va1a4 + va2a4) * (1 - exp(-z4))))/((1 - exp(-z4))^2)^2 
+ (((exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * (exp(-z3) * (va1a3 + va2a3)) 
+ (exp(-z1) * (va1a1 + va2a1) * (exp(-z2) * (va1a2 + va2a2)) 
+ exp(-z1) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-z2) + (exp(-z1) * (exp(-z2) * (va1a2 + va2a2) * (va1a2 + va2a2)) 
+ exp(-z1) * (va1a1 + va2a1) * (exp(-z2) * (va1a2 + va2a2)))) * exp(-z3) 
+ (exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3) * (va1a3 + va2a3)) 
+ (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * (exp(-z3) * (va1a3 + va2a3))))/(1 - exp(-z4)) 
+ (exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3)) + (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) 
+ exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3)) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2)





#start to simplify- replace d2lz4
(exp(-z1) * exp(-z2) * exp(-z3) * exp(-z4) * (va1a4 + va2a4) * (va1a4 + va2a4) 
+ (exp(-z1) * exp(-z2) * exp(-z3) * (va1a3 + va2a3) + (exp(-z1) * exp(-z2) * (va1a2 + va2a2) 
+ exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3)) * (exp(-z4) * (va1a4 + va2a4)))/(1 - exp(-z4))^2
+ (exp(-z1) * exp(-z2) * exp(-z3)) * (exp(-z4) * (va1a4 + va2a4)) * (2 * (exp(-z4) * (va1a4 + va2a4) * (1 - exp(-z4))))/((1 - exp(-z4))^2)^2 
+ d2z4/(1 - exp(-z4)) 
+ (exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3)) + (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) 
+ exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3)) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2)


#replace lz4
(lz4 * exp(-z4) * (va1a4 + va2a4) * (va1a4 + va2a4) 
+ lz4 * (va1a3 + va2a3) * exp(-z4) * (va1a4 + va2a4)
+ lz4 * (va1a2 + va2a2) * exp(-z4) * (va1a4 + va2a4)
+ lz4 * (va1a1 + va2a1) * exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2
+ lz4 * exp(-z4) * (va1a4 + va2a4) * 2 * exp(-z4) * (va1a4 + va2a4) * (1 - exp(-z4))/((1 - exp(-z4))^2)^2 
+ d2z4/(1 - exp(-z4)) 
+ (lz4 * (va1a3 + va2a3) + lz4 * (va1a2 + va2a2) 
+ lz4 * (va1a1 + va2a1) ) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2)


(lz4 * exp(-z4) * (va1a4 + va2a4) * (va1a4 + va2a4))/(1 - exp(-z4))^2 
+ (lz4 * (va1a3 + va2a3) * exp(-z4) * (va1a4 + va2a4)
+ lz4 * (va1a2 + va2a2) * exp(-z4) * (va1a4 + va2a4)
+ lz4 * (va1a1 + va2a1) * exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2

+ lz4 * exp(-z4) * (va1a4 + va2a4) * 2 * exp(-z4) * (va1a4 + va2a4) * (1 - exp(-z4))/((1 - exp(-z4))^2)^2 

+ d2z4/(1 - exp(-z4)) 
+ (lz4 * (va1a3 + va2a3) * exp(-z4) * (va1a4 + va2a4)
+ lz4 * (va1a2 + va2a2) * exp(-z4) * (va1a4 + va2a4)
+ lz4 * (va1a1 + va2a1) * exp(-z4) * (va1a4 + va2a4) ) /(1 - exp(-z4))^2)


#compare with MSY.hpp -  not perfect but almost there
d2z4/(1 - exp(-z4)) 

+ 2*lz4  * exp(-z4) * (va1a4 + va2a4) *( va1a3 + va2a3
+  va1a2 + va2a2 +  va1a1 + va2a1 )/(1 - exp(-z4))^2

+ 2 *lz4 *  (va1a4 + va2a4)* (va1a4 + va2a4) * exp(-z4) * exp(-z4)  /(1 - exp(-z4))^3 

+ (lz4 * exp(-z4) * (va1a4 + va2a4) * (va1a4 + va2a4))/(1 - exp(-z4))^2 

 
#msy.hpp
d2lz(k)(j-1)*sa(h)(j-1)/oa(h)(j) 
+ 2*lz(h)(j-1)*V1*sa(h)(j-1)*V2*sa(h)(j)/oa2
+ 2*lz(h)(j-1)*sa(h)(j-1)*V2*V2*sa(h)(j)*sa(h)(j)/(oa(h)(j)*oa2)
+ lz(h)(j-1)*sa(h)(j-1)*V2*V2*sa(h)(j)/oa2;


#============================================================
# lw derivatives
# these derivatives assume that lw is calculated as follows:

#ex_lwa=lza*exp(-rho*( m + fe * (va1  + ... +van )))


#age1
ex_lw1<- expression(1*exp(-rho*( m + fe * (va1a1  + va2a1 ))))
D(ex_lw1,"fe")
#r output
#dw1df=-(exp(-rho * (m + fe * (va1a1 + va2a1))) * (rho * (va1a1 + va2a1)))

-(exp(-rho * z1) * rho * (va1a1 + va2a1))

#Compare with MSY.hpp
-psa(m_sage)*m_rho*m_Va(h)(k)(m_sage);  




#age2
ex_lw2<- expression( exp(-(m + fe *( va1a1 + va2a1))) *exp(-rho*( m + fe * (va1a2 + va2a2))) )
D(ex_lw2,"fe")
-(exp(-(m + fe * (va1a1 + va2a1))) * (exp(-rho * (m + fe * (va1a2 + 
    va2a2))) * (rho * (va1a2 + va2a2))) + exp(-(m + fe * (va1a1 + 
    va2a1))) * (va1a1 + va2a1) * exp(-rho * (m + fe * (va1a2 + 
    va2a2))))

-(exp(-z1) * (exp(-rho * z2) * (rho * (va1a2 + va2a2))) + exp(-z1) * (va1a1 + va2a1) * exp(-rho * z2))

-(exp(-z1) * (exp(-rho * z2) * rho * (va1a2 + va2a2)) + exp(-z1) * (va1a1 + va2a1) * exp(-rho * z2))

-(lz2 * exp(-rho * z2) * rho * (va1a2 + va2a2) + lz2 * (va1a1 + va2a1) * exp(-rho * z2))

-lz2 * exp(-rho * z2)*(rho * (va1a2 + va2a2) + (va1a1 + va2a1))


-lz(h)(j)*m_rho*m_Va(h)(k)(j)*psa(j); #almost the same, sum of previous selectivities is missing



#age3
ex_lw3<- expression( exp(-(m + fe *( va1a1 + va2a1)))* exp(-(m + fe *( va1a2 + va2a2))) *exp(-rho*( m + fe * (va1a3 + va2a3))) )
D(ex_lw3,"fe")
-(exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
    va2a2))) * (exp(-rho * (m + fe * (va1a3 + va2a3))) * (rho * 
    (va1a3 + va2a3))) + (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + 
    fe * (va1a2 + va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * 
    (va1a1 + va2a1))) * (va1a1 + va2a1) * exp(-(m + fe * (va1a2 + 
    va2a2)))) * exp(-rho * (m + fe * (va1a3 + va2a3))))

-(exp(-z1) * exp(-z2) * exp(-rho * z3) * rho * (va1a3 + va2a3) 
    + (exp(-z1) * exp(-z2) * (va1a2 + va2a2) 
    + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-rho * z3))


-(lz3 * exp(-rho * z3) * rho * (va1a3 + va2a3) 
    + (lz3 * (va1a2 + va2a2) + lz3 * (va1a1 + va2a1) ) * exp(-rho * z3))

-(lz3 * exp(-rho * z3) * rho * (va1a3 + va2a3) 
    + lz3 * (va1a2 + va2a2) * exp(-rho * z3) + lz3 * (va1a1 + va2a1) * exp(-rho * z3)  )

-lz3 * exp(-rho * z3) *(rho * (va1a3 + va2a3) + (va1a2 + va2a2) + (va1a1 + va2a1))

#almost again

-lz(h)(j)*m_rho*m_Va(h)(k)(j)*psa(j); 
#=================================================================================================

#plus group
#assuming  age 4 is plus group
ex_lwplus<-expression(exp(-(m + fe *( va1a1 + va2a1))) *exp(-( m + fe * (va1a2 + va2a2))) *exp(-( m + fe * (va1a3 + va2a3))) * exp(-rho*(m + fe *( va1a4 + va2a4)))/(1-exp(-(m + fe *( va1a4 + va2a4)))) )
D(ex_lwplus,"fe")

-((exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
    va2a2))) * exp(-(m + fe * (va1a3 + va2a3))) * (exp(-rho * 
    (m + fe * (va1a4 + va2a4))) * (rho * (va1a4 + va2a4))) + 
    (exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
        va2a2))) * (exp(-(m + fe * (va1a3 + va2a3))) * (va1a3 + 
        va2a3)) + (exp(-(m + fe * (va1a1 + va2a1))) * (exp(-(m + 
        fe * (va1a2 + va2a2))) * (va1a2 + va2a2)) + exp(-(m + 
        fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * exp(-(m + 
        fe * (va1a2 + va2a2)))) * exp(-(m + fe * (va1a3 + va2a3)))) * 
        exp(-rho * (m + fe * (va1a4 + va2a4))))/(1 - exp(-(m + 
    fe * (va1a4 + va2a4)))) + exp(-(m + fe * (va1a1 + va2a1))) * 
    exp(-(m + fe * (va1a2 + va2a2))) * exp(-(m + fe * (va1a3 + 
    va2a3))) * exp(-rho * (m + fe * (va1a4 + va2a4))) * (exp(-(m + 
    fe * (va1a4 + va2a4))) * (va1a4 + va2a4))/(1 - exp(-(m + 
    fe * (va1a4 + va2a4))))^2)


-((exp(-z1) * exp(-z2) * exp(-z3) * (exp(-rho * z4) * (rho * (va1a4 + va2a4))) 
    + (exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3)) 
    + (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) 
    + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3)) * exp(-rho * z4))/(1 - exp(-z4)) 
    + exp(-z1) * exp(-z2) * exp(-z3) * exp(-rho * z4) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2)

-((lz4 * exp(-rho * z4) * rho * (va1a4 + va2a4) 
    + (lz4 * (va1a3 + va2a3) 
    + (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) 
    + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3)) * exp(-rho * z4))/(1 - exp(-z4)) 
    + lz4 * exp(-rho * z4) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2)

-((lz4 * exp(-rho * z4) * rho * (va1a4 + va2a4) 
    + (lz4 * (va1a3 + va2a3) 
    + exp(-z1) * exp(-z2) * (va1a2 + va2a2) * exp(-z3) 
    + exp(-z1) * (va1a1 + va2a1) * exp(-z2) * exp(-z3) ) * exp(-rho * z4))/(1 - exp(-z4)) 
    + lz4 * exp(-rho * z4) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2)

-((lz4 * exp(-rho * z4) * rho * (va1a4 + va2a4) 
    + lz4 * exp(-rho * z4) *((va1a3 + va2a3) + (va1a2 + va2a2) + (va1a1 + va2a1)) )/(1 - exp(-z4)) 
    + lz4 * exp(-rho * z4) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2)

-(lz4 * exp(-rho * z4) * (rho * (va1a4 + va2a4) + ((va1a3 + va2a3) + (va1a2 + va2a2) + (va1a1 + va2a1)))/(1 - exp(-z4)) 
    + lz4 * exp(-rho * z4) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2)

-lz4 * exp(-rho * z4) * (rho * (va1a4 + va2a4) + ((va1a3 + va2a3) + (va1a2 + va2a2) + (va1a1 + va2a1)))/(1 - exp(-z4)) 
    - lz4 * exp(-rho * z4) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2

-lz4 * exp(-rho * z4) * (rho * (va1a4 + va2a4) + ((va1a3 + va2a3) + (va1a2 + va2a2) + (va1a1 + va2a1)))/(1 - exp(-z4)) 
    - lz4 * exp(-rho * z4) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2


-lz(h)(j-1)*sa(h)(j-1)*m_rho*m_Va(h)(k)(j)/oa(h)(j)
                                    - lz(h)(j-1)*psa(j)*m_Va(h)(k)(j)*sa(h)(j)
                                    /square(oa(h)(j));

#try again

#assuming  age 4 is plus group
ex_lwplus<-expression((exp(-(m + fe *( va1a1 + va2a1))) *exp(-( m + fe * (va1a2 + va2a2))) *exp(-( m + fe * (va1a3 + va2a3))))/(1-exp(-(m + fe *( va1a4 + va2a4)))) * exp(-rho*(m + fe *( va1a4 + va2a4))) )
D(ex_lwplus,"fe")

#output
-((exp(-(m + fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + 
    va2a2))) * exp(-(m + fe * (va1a3 + va2a3))))/(1 - exp(-(m + 
    fe * (va1a4 + va2a4)))) * (exp(-rho * (m + fe * (va1a4 + 
    va2a4))) * (rho * (va1a4 + va2a4))) + ((exp(-(m + fe * (va1a1 + 
    va2a1))) * exp(-(m + fe * (va1a2 + va2a2))) * (exp(-(m + 
    fe * (va1a3 + va2a3))) * (va1a3 + va2a3)) + (exp(-(m + fe * 
    (va1a1 + va2a1))) * (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + 
    va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * 
    exp(-(m + fe * (va1a2 + va2a2)))) * exp(-(m + fe * (va1a3 + 
    va2a3))))/(1 - exp(-(m + fe * (va1a4 + va2a4)))) + (exp(-(m + 
    fe * (va1a1 + va2a1))) * exp(-(m + fe * (va1a2 + va2a2))) * 
    exp(-(m + fe * (va1a3 + va2a3)))) * (exp(-(m + fe * (va1a4 + 
    va2a4))) * (va1a4 + va2a4))/(1 - exp(-(m + fe * (va1a4 + 
    va2a4))))^2) * exp(-rho * (m + fe * (va1a4 + va2a4))))

-((exp(-z1) * exp(-z2) * exp(-z3))/(1 - exp(-z4)) * (exp(-rho * z4) * (rho * (va1a4 + va2a4))) 
    + ((exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3)) 
    + (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) 
    + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3))/(1 - exp(-z4)) 
    + (exp(-z1) * exp(-z2) * exp(-z3)) * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2) * exp(-rho * z4))

-(lz4/(1 - exp(-z4)) * exp(-rho * z4) * rho * (va1a4 + va2a4) 
    + lz4 * ((va1a3 + va2a3) + (va1a2 + va2a2) +  (va1a1 + va2a1) )/(1 - exp(-z4)) 
    + exp(-rho * z4) *  lz4 * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2 )

-(lz4/(1 - exp(-z4)) * exp(-rho * z4) * rho * (va1a4 + va2a4) 
    + lz4 * ((va1a3 + va2a3) + (va1a2 + va2a2) +  (va1a1 + va2a1) )/(1 - exp(-z4)) 
    + exp(-rho * z4) *  lz4 * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2 )

-( exp(-rho * z4) * rho * (va1a4 + va2a4) * lz4/(1 - exp(-z4))
    + lz4 * ((va1a3 + va2a3) + (va1a2 + va2a2) +  (va1a1 + va2a1) )/(1 - exp(-z4)) 
    + exp(-rho * z4) *  lz4 * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2 )

-( lz4* ( exp(-rho * z4) * rho * (va1a4 + va2a4) +  ((va1a3 + va2a3) + (va1a2 + va2a2) + (va1a1 + va2a1)))/(1 - exp(-z4)) 
    + exp(-rho * z4) *  lz4 * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2 )

#compare to MSY.hpp
-( lz4* ( exp(-rho * z4) * rho * (va1a4 + va2a4) +  ((va1a3 + va2a3) + (va1a2 + va2a2) + (va1a1 + va2a1)))/(1 - exp(-z4)) 
    + exp(-rho * z4) *  lz4 * (exp(-z4) * (va1a4 + va2a4))/(1 - exp(-z4))^2 )



-lz(h)(j)*m_rho*m_Va(h)(k)(j)/oa(h)(j)
                                    - lz(h)(j-1)*psa(j)*m_Va(h)(k)(j)*sa(h)(j)
                                    /square(oa(h)(j));


#==============================================
#second derivative
dex_lw1<- expression(-(exp(-rho * (m + fe * (va1a1 + va2a1))) * (rho * (va1a1 + va2a1))))
D(dex_lw1,"fe")

#output
exp(-rho * (m + fe * (va1a1 + va2a1))) * (rho * (va1a1 + va2a1)) * 
    (rho * (va1a1 + va2a1))

exp(-rho * z1) * rho * (va1a1 + va2a1) * rho * (va1a1 + va2a1)

exp(-rho * z1) * (rho * (va1a1 + va2a1))^2

#compare with MSY.hpp
psa(m_sage)*square(m_rho)*square(m_Va(h)(k)(m_sage)); 

#age 2
dex_lw2<- expression( -(exp(-(m + fe * (va1a1 + va2a1))) * (exp(-rho * (m + fe * (va1a2 + 
    va2a2))) * (rho * (va1a2 + va2a2))) + exp(-(m + fe * (va1a1 + 
    va2a1))) * (va1a1 + va2a1) * exp(-rho * (m + fe * (va1a2 + 
    va2a2)))) )
D(dex_lw2,"fe")

#result
exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (exp(-rho * 
    (m + fe * (va1a2 + va2a2))) * (rho * (va1a2 + va2a2))) + 
    exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (va1a1 + 
        va2a1) * exp(-rho * (m + fe * (va1a2 + va2a2))) + (exp(-(m + 
    fe * (va1a1 + va2a1))) * (exp(-rho * (m + fe * (va1a2 + va2a2))) * 
    (rho * (va1a2 + va2a2)) * (rho * (va1a2 + va2a2))) + exp(-(m + 
    fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (exp(-rho * (m + 
    fe * (va1a2 + va2a2))) * (rho * (va1a2 + va2a2))))

exp(-z1) * (va1a1 + va2a1) * (exp(-rho * z2) * (rho * (va1a2 + va2a2))) 
+ exp(-z1) * (va1a1 + va2a1) * (va1a1 +  va2a1) * exp(-rho * (m + fe * (va1a2 + va2a2))) 
+ (exp(-z1) * (exp(-rho * z2) * (rho * (va1a2 + va2a2)) * (rho * (va1a2 + va2a2))) 
+ exp(-z1) * (va1a1 + va2a1) * exp(-rho * z2) * rho * (va1a2 + va2a2))


exp(-z1) * (va1a1 + va2a1) * exp(-rho * z2) * rho * (va1a2 + va2a2) 
+ exp(-z1) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-rho * z2) 
+ (exp(-z1) * exp(-rho * z2) * rho * (va1a2 + va2a2) * rho * (va1a2 + va2a2) 
+ exp(-z1) * (va1a1 + va2a1) * exp(-rho * z2) * rho * (va1a2 + va2a2))


