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
#exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + va2a1) * (exp(-(m + 
#    fe * (va1a2 + va2a2))) * (va1a2 + va2a2)) + exp(-(m + fe * 
#    (va1a1 + va2a1))) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-(m + 
#    fe * (va1a2 + va2a2))) + (exp(-(m + fe * (va1a1 + va2a1))) * 
#    (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + va2a2) * (va1a2 + 
#        va2a2)) + exp(-(m + fe * (va1a1 + va2a1))) * (va1a1 + 
#    va2a1) * (exp(-(m + fe * (va1a2 + va2a2))) * (va1a2 + va2a2)))
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
#
#  (exp(-z1) * exp(-z2) * (va1a2 + va2a2) + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3) * (va1a3 + va2a3) 
#+ (exp(-z1) * (va1a1 + va2a1) * exp(-z2) * (va1a2 + va2a2) 
#	+ exp(-z1) * (va1a1 + va2a1) * (va1a1 + va2a1) * exp(-z2)
#	+ (exp(-z1) * exp(-z2) * (va1a2 + va2a2) * (va1a2 + va2a2) 
#	+ exp(-z1) * (va1a1 + va2a1) * (exp(-z2) * (va1a2 + va2a2)))) * exp(-z3) 
#+ exp(-z1) * exp(-z2) * (exp(-z3) * (va1a3 + va2a3) * (va1a3 + va2a3)) 
#	+ (exp(-z1) * (exp(-z2) * (va1a2 + va2a2)) + exp(-z1) * (va1a1 + va2a1) * exp(-z2)) * exp(-z3) * (va1a3 + va2a3)
#
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


#=========================================================================================================================
# Second derivative of lz -- second derivative



