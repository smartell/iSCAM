#|-----------------------------------------------------------------------------------|#
#| R4Milka.R                                                                         |#
#| AUTHOR: Catarina Wor (adapted from Steven Martell's R4iSCAM.r)                    |#
#|-----------------------------------------------------------------------------------|#
#| NOTES, if using sublime text 3.0
#|    use ⌘ + \ to set working directory to the directory of the current file.
#|    use ⌘ + enter to source a single line or highligted lines.
#|    use ⌘ + B to source entire file.

# |----------------------------------------------------------------------------------|
# | DEFINITIONS
# |----------------------------------------------------------------------------------|
# | .PWD       <- Global Parent Working Directory for R-scripts
# | .FIGUREDIR <- Directory for saving figures.
# | .RFILES    <- List of R functions to source from the lib directory.
# .PWD        needs to be changed when creating a new example
# .PWD        <- "/Users/catarinawor/Documents/iSCAM/examples/DEMO/R"
.PWD        <- "/Users/stevenmartell1/Documents/iSCAM-project/src/R"
.LIB        <- "../../../src/R/MSElib/"
setwd(.PWD)
.FIGUREDIR  <- "../FIGS/"
.RFILES     <- list.files(.LIB,pattern="\\.[Rr]$")
#.BOOLREADFN <- TRUE
.OVERLAY    <- FALSE
require(ggplot2)
.THEME      <- theme_bw(11)
.UNITS      <- "(mlb)"

# | Labels for gear sex area and group index.
.GEAR  = c("Directed","Wasteage","Bycatch","Sport","Personal","Setline Survey")
.SEX   = c("F & M","Female","Male")
.AREA  = c("Coast wide")
.GROUP = c("Pacific Halibut")


.MODELDIRS   <- "../DATA"
.MODELNAME   <- list.files(.MODELDIRS,pattern="\\.Rdata",full.name=TRUE)
for(i in 1:(length(.MODELNAME)))load(.MODELNAME[i])



M <- list(mse.data, rawmse.data)
 


for(nm in .RFILES) source(file.path(.LIB, nm), echo=FALSE)
.MSEplotMilkaSpawnBiomass( M )
.MSEplotDepletion( M )
.MSEplotCatch( M )
.MSEplotAAV( M )
.MSEplotSubLegal( M )
.MSEplotWastage( M )
.MSEplotEfficiency( M )
#.saveImages()




	

