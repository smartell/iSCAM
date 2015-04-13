#tab_msyTable.R

.msyTable <- function(Scenario,objs)
{
	ss  <- Scenario %>% group_by(prefix) %>% filter(Yield==max(Yield))
	mss <- subset(melt(ss,id.vars=c("prefix","regarea")),variable %in% objs)
	css <- dcast(mss,prefix~variable)
	
	return(css)
}
