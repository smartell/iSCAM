# 
# SPATIAL EQUILIBRIUM MODEL
# 
nArea   <- 2
# Gravity model for movement kernel
res     <- 0.5
gw      <- runif(nArea); gw[nArea] <- 1 - sum(gw[-nArea])
G       <- matrix(log(gw),nArea,nArea,byrow=TRUE)
diag(G) <- diag(G) + res
G       <- exp(G)/rowSums(exp(G))


A   <- 5
nB  <- A * nArea
age <- 1:A
m   <- 0.85
lx  <- exp(-m)^(age-min(age)); lx[A] <- lx[A]/(1-exp(-m))
S   <- matrix(0,nB,nB)
aG  <- matrix(0,A,A)
r   <- rep(c(1,rep(0,length=A-1)),nArea)

# survival
for (i in 1:nB) 
{
  if(i %% A)
    S[i,i+1] = exp(-m)
  else
    S[i,i] = exp(-m)
}

# movement
# nA <- 2 * A
# G  <- matrix(0,nA,nA)
n  <- rep(lx,nArea)

# basis vectors
ej = diag(1,nArea*A)
