# ——————————————————————————————————————————————————————————————————————————— #
# NOTES: The mse.data object is a list of data.frames that comes from the
#        saveMSEdataframe.R routine.  The following code disassbles this list
#        into several dataframes to be used in the Shiny application.
# DATA FRAMES:
#   - BIO.DF -> spawning biomass 
#   - CAT.DF -> catch related variables
#   - SUB.DF -> sublegal and wastage related variables
#   - MOT.DF -> data frame to be used with gvisMotionChart
# ——————————————————————————————————————————————————————————————————————————— #
BIO.DF <- mse.data$biomass.df
CAT.DF <- mse.data$catch.df
SUB.DF <- mse.data$sublegal.df
AAV.DF <- mse.data$AAV.df

MRG.DF <- merge(BIO.DF,CAT.DF,by=c("Scenario","Procedure","Year"))
MRG.DF <- merge(MRG.DF,AAV.DF,by=c("Scenario","Procedure","Year","gear","area","group"))
MSE.DF <- merge(MRG.DF,SUB.DF,by=c("Scenario","Procedure","Year","gear","area","sex","group"))
# Restricted data frame for gvisMotionChart for increased speed & less clutter.
hdr <- c("Scenario","Procedure","Year","t.Bt0.5","t.Dt0.5","ct50","AAV50")
MOT.DF <- MRG.DF[,which(names(MRG.DF) %in% hdr)]
MOT.DF$idvar <- paste0(MOT.DF$Scenario,MOT.DF$Procedure)
names(MOT.DF) <- c("Scenario","Procedure","Year","Median biomass","Depletion",
                   "Median catch","Median AAV","idvar")
