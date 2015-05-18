# 
# Load Libraries
# 
if(!require("shiny"))         install.packages("shiny")
if(!require("ggplot2"))       install.packages("ggplot2")
if(!require("reshape2"))      install.packages("reshape2")
if(!require("Rcpp"))          install.packages("Rcpp")
if(!require("markdown"))      install.packages("markdown")
if(!require("shinythemes"))   install.packages("shinythemes")
if(!require("dplyr"))         install.packages("dplyr")
if(!require("grid"))          install.packages("grid")



# Source R-scripts.
.REGAREA <- c("all","2A","2B","2C","3A","3B","4A","4B","4CDE")
.THEME      <- theme_bw(11)
.OVERLAY    <- FALSE
.UNITS      <- "(Mlb)"
.LIB        <- "data/"
.RFILES     <- list.files(.LIB,pattern="\\.[Rr]$",recursive=TRUE)
for(nm in .RFILES) source(file.path(.LIB, nm), echo=FALSE)


#
# Source Rcpp scripts
# 
sourceCpp("./data/src/halitosis.cpp")

# 
# LOAD estimated growth parameters from 2014 SAA data.
#
load("./data/equilibrium/halGrowthParams.Rdata")


# 
# LOAD MSE.Rdata 
# 
load("./data/mse/MSE.Rdata")


#
# Parameters names required for halitotis.cpp model
#
paramNames <- c("size_limit",
                "discard_mortality_rate",
                "selex_fishery",
                "selex_bycatch",
                "selex_bycatch_desc",
                "num_bycatch",
                "five",
                "ten",
                "twenty",
                "forty",
                "regArea",
                "linf_dev",
                "vonk_dev",
                "maternal_effect")