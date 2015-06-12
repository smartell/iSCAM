# tab_meyTable.R

.meyTable <- function(Scenario,objs)
{
	ss  <- Scenario %>% group_by(prefix) %>% filter(Landed.Value==max(Landed.Value))
	mss <- subset(melt(ss,id.vars=c("prefix","regarea")),variable %in% objs)
	css <- dcast(mss,prefix~variable)
	
	return(css)
}