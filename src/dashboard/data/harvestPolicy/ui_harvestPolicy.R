# 
# SPATIAL EQUILIBRIUM MODEL
# 

A   <- 10
age <- 1:A
m   <- 0.35
lx  <- exp(-m)^(age-min(age)); lx[A] <- lx[A]/(1-exp(-m))
S   <- matrix(0,A,A)
r   <- c(1,rep(0,length=A-1))

# survival
for (i in age) 
{
  if(i < max(age))
    S[i,i+1] = exp(-m)
  else
    S[i,i] = exp(-m)
}

# movement
nA <- 2 * A
G  <- matrix(0,nA,nA)
