# 
# Load Libraries
# 

if(!require("shiny"))         install.packages("shiny")
if(!require("shinythemes"))   install.packages("shinythemes")
if(!require("ggplot2"))       install.packages("ggplot2")
if(!require("reshape2"))      install.packages("reshape2")


source("data/TMA.R",local=TRUE)

paramNames <- c("ng")



