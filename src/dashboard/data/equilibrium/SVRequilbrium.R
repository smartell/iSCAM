#SVRequilbrium.R



getParams <- function(prefix,input) {
      print("Oh boy")
      # input[[paste0(prefix, "_recalc")]]

      params <- lapply(paramNames, function(p) {
        input[[paste0(prefix, "_", p)]]
      })
      names(params) <- paramNames
      params <- c(params,prefix=prefix)
      # print(params)
      params
    }

.updateDiscardMortality <- function(prefix,input)
{
  
}


.plotEquilFe <- function(Scenario,objs)
{
  SS  <- Scenario  %>% group_by(prefix) %>% filter(Yield == max(Yield))
  mSS <- subset(melt(SS,id.vars=1:7),variable %in% objs)
  sdf <- subset(melt(Scenario,id.vars=1:7),variable %in% objs)

  p <- ggplot(sdf,(aes(Fe,value,col=prefix))) + geom_line(size=1.25)
  
  p <- p + geom_segment(aes(x=Fe,xend=Fe,y=value,yend=0,col=prefix),data=mSS, arrow = arrow(length = unit(0.25,"cm")),alpha=0.4,size=1.25)
  p <- p + geom_segment(aes(x=Fe,xend=0,y=value,yend=value,col=prefix),data=mSS, arrow = arrow(length = unit(0.25,"cm")),alpha=0.4,size=1.25)
  p <- p + facet_wrap(~variable,scales="free")
  p <- p + labs(x="Fishing Intensity (Fe)",col="Procedure",y="")
  p <- p + theme_bw(18) + theme(legend.position="top") 
  print(p)
  

}


.plotFishSelex <- function(pars)
{
    # assume units are in inches
  n <- length(pars)
  x = seq(10,200,by=2.5)/2.54
  df <- NULL
  for( i in 1:n )
  {
    pref <- pars[[i]]$prefix
    bL50 <- pars[[i]]$selex_fishery[1]
    bL95 <- pars[[i]]$selex_fishery[2]
    # bR50 <- pars[[i]]$selex_bycatch_desc[1]
    # bR95 <- pars[[i]]$selex_bycatch_desc[2]

    sel  <- plogis95_cpp(x,bL50,bL95)#*plogis95_cpp(x,bR95,bR50)
    df   <- rbind(df,data.frame("prefix"=pref,"len"=x,"sel"=sel))
  }
  
  p <- ggplot(df,aes(len,sel,col=prefix)) + geom_line(size=1.1)
  p <- p + ylim(c(0,1))
  p <- p + labs(x="Length (in.)",y="Selectivity",col="Procedure")
  print(p + theme_bw())

}


.plotBycatchSelex <- function(pars)
{
  # assume units are in inches
  n <- length(pars)
  x = seq(10,200,by=2.5)/2.54
  df <- NULL
  for( i in 1:n )
  {
    pref <- pars[[i]]$prefix
    bL50 <- pars[[i]]$selex_bycatch[1]
    bL95 <- pars[[i]]$selex_bycatch[2]
    bR50 <- pars[[i]]$selex_bycatch_desc[1]
    bR95 <- pars[[i]]$selex_bycatch_desc[2]

    sel  <- plogis95_cpp(x,bL50,bL95)*plogis95_cpp(x,bR95,bR50)
    df   <- rbind(df,data.frame("prefix"=pref,"len"=x,"sel"=sel))
  }
  
  p <- ggplot(df,aes(len,sel,col=prefix)) + geom_line(size=1.1)
  p <- p + ylim(c(0,1))
  p <- p + labs(x="Length (in.)",y="Selectivity",col="Procedure")
  print(p + theme_bw())
}