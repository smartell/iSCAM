# This is an example of finding the root of the Baranov Catch equation for 
# multiple fleets and a two sex population model, where the sex ratio of the
# observed catch is unknown.  Developed for the iSCAM operating model.
# 

# Steve Martell
# Sept. 10th, 2013.
# To be implemented in Baranov.h
set.seed(991)


# Population variables (Female, Male).
age  <- 1:30
sex  <- 1:2
M    <- c(0.2, 0.25)
surv <- function(x) exp(-x*(age-1))
N    <- sapply(M,surv)

# Fisheries
G    <- 3
gear <- 1:G
F    <- exp(rnorm(G,log(0.2),0.2))
ah   <- matrix(rpois(2*G,10),ncol=2)
gh   <- 2.5
sel  <- function(x) plogis(age,x,gh)
V    <- apply(ah,c(1,2),sel)

# True catch
# Cobs <- N*F*V*(1-exp(-M-F*V))/(M+F*V)
Fa   <- Ca <- V
Za   <- Sa <- Oa <- N
Cg   <- rep(0,length=G)
for(h in sex)
{
	Fa[,,h] <- t(F * t(V[,,h]))
	Za[,h]  <- M[h] + rowSums(Fa[,,h])
	Ca[,,h] <- N[,h]*Fa[,,h]*(1-exp(-Za[,h]))/(Za[,h]) 
	matplot(Ca[,,h])
}
Cobs <- apply(Ca,2,sum)

# ------------------------------------------ #
# Iterative soln for fishing mortality rate  #
# ------------------------------------------ #
# Start with Popes' approximation
TOL    <- 1.e-4
Ni     <- rowSums(sapply(sex,FUN=function(sex){(N[,sex]*exp(-0.5*M[sex]))%*%V[,,sex]}))
Fpope  <- Cobs / Ni
Fi     <- Fpope
Jacobi <- matrix(0,G,G) 

for (iter in 1:300)
{

	for(h in sex)
	{
		Fa[,,h] <- t(Fi * t(V[,,h]))
		Za[,h]  <- M[h] + rowSums(Fa[,,h])
		Sa[,h]  <- exp(-Za[,h])
		Oa[,h]	<- 1.0 - Sa[,h]
		Ca[,,h] <- N[,h]*Fa[,,h]*(1-exp(-Za[,h]))/(Za[,h])
		for(i in gear)
		{
			for(j in gear)
			{
				if(i==j)
				{
					k1   <- N[,h] * V[,i,h] / Za[,h]
					k2   <- k1 * V[,i,h]
					k3   <- k2 / Za[,h]
					dcdf <- -(k1%*%Oa[,h]) - Fi[i]*(k2%*%Sa[,h]) + Fi[i]*(k3%*%Oa[,h])

					Jacobi[i,j] <- Jacobi[i,j] + as.double(dcdf)
				}
				else
				{
					t1   <- Fi[i] * N[,h] * V[,i,h] / Za[,h]
					t2   <- t1 * V[,i,h]
					t3   <- t2 / Za[,h] 
					dcdf <- -(t2%*%Sa[,h]) + (t3%*%Oa[,h])

					Jacobi[j,i] <- Jacobi[j,i] + as.double(dcdf)
				}
			}
		}
	}
	Chat <- apply(Ca,2,sum)
	fx   <- Cobs - Chat
	invJ <- -solve(Jacobi)
	Fi   <- Fi + as.double( fx%*%invJ )
	# cat("fx = ", round(fx,4)," Sum = ",sum(abs(fx)),"\n")
	if(sum(abs(fx))<=TOL) break;
}



# ------------------------------------------ #
# | RESULTS                                | #
# ------------------------------------------ #
cat("| --------------------------------------------------------- |\n")
cat("| RESULTS * RESULTS * RESULTS * RESULTS * RESULTS * RESULTS |\n")
cat("| --------------------------------------------------------- |\n")
cat("| Observed catch vector               := ",round(Cobs,3),"|\n")
cat("| Predicted catch vector              := ",round(Chat,3),"|\n")
cat("| --------------------------------------------------------- |\n")
cat("| True fishing mortality              := ",round(F,3),"|\n")
cat("| Newton method fishing mortality     := ",round(Fi,3),"|\n")
cat("| Root of the function (should = 0)   := ",signif(sum(abs(fx)),3),"\n")
cat("| --------------------------------------------------------- |\n")
cat("| Fishing mortality from Popes approx := ",round(Fpope,3),"|\n")
cat("| --------------------------------------------------------- |\n")



