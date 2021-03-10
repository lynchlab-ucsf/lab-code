plot.alphaDiversity <- function(varname,alphadivs,grpvar,grps,mydata,stars=TRUE) {
empty.frame <- NULL
for(i in alphadivs) {
   for(j in grps) {
      if(j=="Overall") { #Allows you to use "Overall" as a group, and include it as a category in your plot
         loop.data <- mydata
      } else { # For the rest of the groups
      loop.data <- mydata[mydata[,grpvar] == j,]
      }
      lev.tot <- length(levels(as.factor(loop.data[,varname])))
      my.lm <- lmer(loop.data[,i] ~ factor(loop.data[,varname]) + (1|studyid), data=loop.data)
         est <- NULL
         se <- NULL
         ci.upper <- NULL
         ci.lower <- NULL
         level <- NULL
         p.val <- NULL
      for(k in 1:lev.tot) {
#We want to grab the intercept's estimate first so that I can compare the rest of the estimates to that number
         if(k==1) {
         est[k] <- summary(my.lm)$coef[k,1]
         se[k] <- summary(my.lm)$coef[k,2]
         } else {
         est[k] <- est[1] + summary(my.lm)$coef[k,1]
         se[k] <- summary(my.lm)$coef[k,2]
         p.val[k] <- summary(my.lm)$coef[k,5]
         }
         ci.lower[k] <- est[k]-1.96*se[k]
         ci.upper[k] <- est[k]+1.96*se[k]
         level[k] <- levels(as.factor(loop.data[,varname]))[k]
      }
      aline <- cbind(i,j,level,est,se,ci.lower,ci.upper,p.val)
      empty.frame <- rbind(aline,empty.frame)
   }
}
my.plot.df <- data.frame(empty.frame)
my.plot.df$est <- as.numeric(as.character(my.plot.df$est))
my.plot.df$se <- as.numeric(as.character(my.plot.df$se))
my.plot.df$ci.lower <- as.numeric(as.character(my.plot.df$ci.lower))
my.plot.df$ci.upper <- as.numeric(as.character(my.plot.df$ci.upper))
my.plot.df$stars[as.numeric(as.character(my.plot.df$p.val))<0.05] <- "*"
my.plot.df$stars[as.numeric(as.character(my.plot.df$p.val))<0.01] <- "**"
my.plot.df$stars[as.numeric(as.character(my.plot.df$p.val))<0.005] <- "***"

#Relevel grps and alphadivs in the specified order (according to function code)
my.plot.df$j <- factor(my.plot.df$j,levels=grps)
my.plot.df$i <- factor(my.plot.df$i,levels=alphadivs)
#So that values don't overlap:
pd <- position_dodge(0.5)

#Define non-variable labels for each of your alpha diversity measures (ie, instead of showing "chao1", you label it as "Richness")
labels <- c(chao1 = "Richness", equitability="Evenness",PD_whole_tree = "Phylogenetic Diversity")
p1 <- ggplot(my.plot.df,aes(x=j,y=est,group=level)) + geom_point(position=pd)+ geom_errorbar(aes(ymin=ci.lower,ymax=ci.upper),width=0.1,position=pd) + facet_wrap(~i,nrow=1,scales="free",labeller=labeller(i = labels)) + aes(color=as.factor(level)) + xlab("Treatment Group") + ylab("Diversity Score") + theme(legend.title=element_blank())
if(stars == TRUE) {
pf <- p1 + geom_text(data=my.plot.df,aes(label=stars,x=as.numeric(j),y=Inf,vjust=1.5),na.rm=TRUE,col="black",size=6)
}
if(stars == FALSE) {
pf <- p1 + geom_text(data=my.plot.df,aes(label=p.val,x=as.numeric(j),y=Inf,vjust=1.5),na.rm=TRUE,col="black",size=6)
}

pf
}
#Example:
#plot.alphaDiversity("exac",c("chao1","equitability","PD_whole_tree"),"group",c("Overall","Placebo","ICS","Xolair"),mapping,stars=FALSE)
