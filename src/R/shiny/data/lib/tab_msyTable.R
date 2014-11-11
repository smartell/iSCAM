#tab_msyTable.R

.msyTable <- function(Scenario,objs)
{
	ss  <- Scenario %>% group_by(prefix) %>% filter(Yield==max(Yield))
	mss <- subset(melt(ss,id.vars="prefix"),variable %in% objs)
	css <- dcast(mss,prefix~variable)
	
	return(css)
}
