#tab_sprTable.R

.sprTable <- function(Scenario,objs)
{

	ss  <- Scenario %>% group_by(prefix) %>% filter(SPR<=0.46) %>% filter(SPR==max(SPR))
	mss <- subset(melt(ss,id.vars=c("prefix","regarea")),variable %in% objs)
	css <- dcast(mss,prefix~variable)
	
	return(css)
}