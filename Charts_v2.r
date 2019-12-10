##############################################################################
# This reads in simulation results and plots them
###########################################

#Libraries
library(ggplot2)
library(hrbrthemes)
library(ggpubr)
library(Hmisc)

# Functions
#-------------
rename <- function(df){
  names(df) <- c('incorrect','joint','correlation','marginal','matching_variable','z_vars','corr_scen','yz_HD')
  df$MatchMeth <- 'HotDeck'
  df[grepl("DPMPM",rownames(df)),9] <- 'DPMPM'
 return(df)
}

#Read data
#--------------

z2 <-  readRDS(file = 'StatMatch_v2_Checks_z2.rds')
z2 <- rename(z2)

z8 <-  readRDS(file = 'StatMatch__v2_Checks_z8.rds')
z8 <- rename(z8)


to_plot <- rbind(z2,z8)
to_plot$z_vars <- as.factor(to_plot$z_vars)

to_plot <- to_plot[!(to_plot$corr_scen=="low" | to_plot$corr_scen=='high'),]


#Graphs
#------------

measures <- c('incorrect', 'joint', 'correlation', 'marginal','yz_HD')
ylabels <- c('proportion of incorrect','Hellinger Distance','squared error of correlation difference','Hellinger Distance', 'Hellinger Distance')
# New facet label names for supp variable
supp.labs <- c("2 Z variables", "8 Z variables")
names(supp.labs) <- c("2", "8")

for(m in 1:length(measures)){
  ggplot(to_plot, aes_string(x='corr_scen', y=measures[m], color='MatchMeth', shape='MatchMeth',label='matching_variable')) + #[to_plot$z_vars==8,]
    stat_summary(fun.y = mean, fun.ymax = mean, fun.ymin = mean,#fun.data="median_hilow",# fun.args = list(mult=1),
                     geom="crossbar", width = 0.6,
                     position = position_dodge(0.5)) +
    geom_point(size=4,position = position_dodge(0.5),alpha = 1/3.5) +
    theme(legend.position="top") + xlab('correlation scenario') +
    labs(color='Method' ,shape='Method') +  ylab(ylabels[m]) +
    facet_wrap(~z_vars, labeller = labeller(z_vars = supp.labs)) +
    #geom_text(hjust=0, vjust=0,position = position_dodge(0.5)) + #ggtitle("8 Z variables") +
    scale_colour_manual(values=c("#f7180c","#0e1247"),aesthetics = "colour")
  ggsave(paste0(measures[m],'_Z.png'))
}


#Separate plots for HD and DP
plot_DP <- to_plot[grepl("DPMPM",rownames(to_plot)),]
plot_HD <- to_plot[grepl("HotDeck",to_plot$MatchMeth),]

for(m in 1:length(measures)){

  ggplot(to_plot, aes_string(y=measures[m], x='corr_scen', fill='z_vars')) +
    geom_boxplot(color="black") +
    theme(legend.position="top") +
    facet_wrap(~MatchMeth) +
    scale_fill_manual(values=c("#E69F00","#999999"))
  ggsave(paste0(measures[m],'.png'))
}


library(stargazer)
stargazer(plot_DP[(plot_DP$z_vars==2 & plot_DP$corr_scen=='highhigh'),])
stargazer(plot_DP[(plot_DP$z_vars==8 & plot_DP$corr_scen=='highhigh'),])
