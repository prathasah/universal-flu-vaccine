library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(lattice)
library(reshape)
library(RColorBrewer)   # for brewer.pal(...)
library(cowplot)
library(Hmisc)
library(viridis)
################################################################

df1_1 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.csv", header=T)
df1_2 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.2.csv", header=T)
df1_3 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.4.csv", header=T)
df1_4 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.6.csv", header=T)
df1_5 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.8.csv", header=T)
df1_6 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_1.csv", header=T)

df2_1 <- read.csv("../data/optimization_results/optimization_results_hospitalizations_1000iter_0.csv", header=T)
df2_2 <- read.csv("../data/optimization_results/optimization_results_hospitalizations_1000iter_0.2.csv", header=T)
df2_3 <- read.csv("../data/optimization_results/optimization_results_hospitalizations_1000iter_0.4.csv", header=T)
df2_4 <- read.csv("../data/optimization_results/optimization_results_hospitalizations_1000iter_0.6.csv", header=T)
df2_5 <- read.csv("../data/optimization_results/optimization_results_hospitalizations_1000iter_0.8.csv", header=T)
df2_6 <- read.csv("../data/optimization_results/optimization_results_hospitalizations_1000iter_1.csv", header=T)

df3_1 <- read.csv("../data/optimization_results/optimization_results_deaths_1000iter_0.csv", header=T)
df3_2 <- read.csv("../data/optimization_results/optimization_results_deaths_1000iter_0.2.csv", header=T)
df3_3 <- read.csv("../data/optimization_results/optimization_results_deaths_1000iter_0.4.csv", header=T)
df3_4 <- read.csv("../data/optimization_results/optimization_results_deaths_1000iter_0.6.csv", header=T)
df3_5 <- read.csv("../data/optimization_results/optimization_results_deaths_1000iter_0.8.csv", header=T)
df3_6 <- read.csv("../data/optimization_results/optimization_results_deaths_1000iter_1.csv", header=T)

df4_1 <- read.csv("../data/optimization_results/optimization_results_DALY_1000iter_0.csv", header=T)
df4_2 <- read.csv("../data/optimization_results/optimization_results_DALY_1000iter_0.2.csv", header=T)
df4_3 <- read.csv("../data/optimization_results/optimization_results_DALY_1000iter_0.4.csv", header=T)
df4_4 <- read.csv("../data/optimization_results/optimization_results_DALY_1000iter_0.6.csv", header=T)
df4_5 <- read.csv("../data/optimization_results/optimization_results_DALY_1000iter_0.8.csv", header=T)
df4_6 <- read.csv("../data/optimization_results/optimization_results_DALY_1000iter_1.csv", header=T)

df1 <- rbind(df1_1,df1_2, df1_3,df1_4,df1_5, df1_6)
df2 <- rbind(df2_1,df2_2, df2_3,df2_4,df2_5, df2_6)
df3 <- rbind(df3_1,df3_2, df3_3,df3_4,df3_5, df3_6)
df4 <- rbind(df4_1,df4_2, df4_3,df4_4,df4_5, df4_6)

dt1 <- read.csv("../data/random_vaccination_results/Incidence_19Feb2018_v1.csv")
dt2 <- read.csv("../data/random_vaccination_results/Hospitalizations_19Feb2018_v1.csv")
dt3 <- read.csv("../data/random_vaccination_results/Mortality_19Feb2018_v1.csv")
dt4 <- read.csv("../data/random_vaccination_results/DALY_19Feb2018_v1.csv")

di1 <- read.csv("../data/typical_vaccination_results/Incidence_19Feb2018_v1.csv")
di2 <- read.csv("../data/typical_vaccination_results/Hospitalizations_19Feb2018_v1.csv")
di3 <- read.csv("../data/typical_vaccination_results/Mortality_19Feb2018_v1.csv")
di4 <- read.csv("../data/typical_vaccination_results/DALY_19Feb2018_v1.csv")

############################################
#3subset data
str(dt1)
optimization_columns <- c("iter1", "scenario", "objective_1", "objective_2",
                          "objective_3", "objective_4", "objective_5", "objective_6", "objective_7", "objective_8", "objective_9",
                          "objective_10", "objective_11", 
                          "objective_12", "objective_13", "objective_14", "objective_15", "objective_16", "objective_17")
df1 <- df1[,optimization_columns]
df2 <- df2[,optimization_columns]
df3 <- df3[,optimization_columns]
df4 <- df4[,optimization_columns]

colnames(df1)[3:length(colnames(df1))] <- c( 'age0', 'age1', 'age5', 'age10', 'age15', 'age20', 'age25' , 'age30', 'age35', 'age40',
                      'age45', 'age50', 'age55', 'age60', 'age65', 'age70', 'age75')
colnames(df2)[3:length(colnames(df2))] <- c( 'age0', 'age1', 'age5', 'age10', 'age15', 'age20', 'age25' , 'age30', 'age35', 'age40',
                                             'age45', 'age50', 'age55', 'age60', 'age65', 'age70', 'age75')
colnames(df3)[3:length(colnames(df3))] <- c( 'age0', 'age1', 'age5', 'age10', 'age15', 'age20', 'age25' , 'age30', 'age35', 'age40',
                                             'age45', 'age50', 'age55', 'age60', 'age65', 'age70', 'age75')
colnames(df4)[3:length(colnames(df4))] <- c( 'age0', 'age1', 'age5', 'age10', 'age15', 'age20', 'age25' , 'age30', 'age35', 'age40',
                                             'age45', 'age50', 'age55', 'age60', 'age65', 'age70', 'age75')


vax_columns <- c('iter', 'scenario', 'age0', 'age1', 'age5', 'age10', 'age15', 'age20', 'age25' , 'age30', 'age35', 'age40',
                 'age45', 'age50', 'age55', 'age60', 'age65', 'age70', 'age75')

dt1 <- dt1[,vax_columns]
dt2 <- dt2[,vax_columns]
dt3 <- dt3[,vax_columns]
dt4 <- dt4[,vax_columns]

di1 <- di1[,vax_columns]
di2 <- di2[,vax_columns]
di3 <- di3[,vax_columns]
di4 <- di4[,vax_columns]
##########################################
## compute mean over multiple iterations
df11 <- data.frame(aggregate(.~ df1$scenario, df1[,3:19], mean))
colnames(df11)[1] <- "scenario"
df22 <- data.frame(aggregate(.~ df2$scenario, df2[,3:19], mean))
colnames(df22)[1] <- "scenario"
df33 <- data.frame(aggregate(.~ df3$scenario, df3[,3:19], mean))
colnames(df33)[1] <- "scenario"
df44 <- data.frame(aggregate(.~ df4$scenario, df4[,3:19], mean))
colnames(df44)[1] <- "scenario"

dt11 <- data.frame(aggregate(.~ dt1$scenario, dt1[,3:19], mean))
colnames(dt11)[1] <- "scenario"
dt22 <- data.frame(aggregate(.~ dt2$scenario, dt2[,3:19], mean))
colnames(dt22)[1] <- "scenario"
dt33 <- data.frame(aggregate(.~ dt3$scenario, dt3[,3:19], mean))
colnames(dt33)[1] <- "scenario"
dt44 <- data.frame(aggregate(.~ dt4$scenario, dt4[,3:19], mean))
colnames(dt44)[1] <- "scenario"

di11 <- data.frame(aggregate(.~ di1$scenario, di1[,3:19], mean))
colnames(di11)[1] <- "scenario"
di22 <- data.frame(aggregate(.~ di2$scenario, di2[,3:19], mean))
colnames(di22)[1] <- "scenario"
di33 <- data.frame(aggregate(.~ di3$scenario, di3[,3:19], mean))
colnames(di33)[1] <- "scenario"
di44 <- data.frame(aggregate(.~ di4$scenario, di4[,3:19], mean))
colnames(di44)[1] <- "scenario"
##########################################
##convert wide to long
library(tidyr)
age_list <- c( 'age0', 'age1', 'age5', 'age10', 'age15', 'age20', 'age25' , 'age30', 'age35', 'age40',
               'age45', 'age50', 'age55', 'age60', 'age65', 'age70', 'age75')

df11<-gather(df11, age_group, value, age0:age75, factor_key=TRUE)
df22<-gather(df22, age_group, value, age0:age75, factor_key=TRUE)
df33<-gather(df33, age_group, value, age0:age75, factor_key=TRUE)
df44<-gather(df44, age_group, value, age0:age75, factor_key=TRUE)

dt11<-gather(dt11, age_group, value, age0:age75, factor_key=TRUE)
dt22<-gather(dt22, age_group, value, age0:age75, factor_key=TRUE)
dt33<-gather(dt33, age_group, value, age0:age75, factor_key=TRUE)
dt44<-gather(dt44, age_group, value, age0:age75, factor_key=TRUE)

di11<-gather(di11, age_group, value, age0:age75, factor_key=TRUE)
di22<-gather(di22, age_group, value, age0:age75, factor_key=TRUE)
di33<-gather(di33, age_group, value, age0:age75, factor_key=TRUE)
di44<-gather(di44, age_group, value, age0:age75, factor_key=TRUE)

######################################
df11$type <- "Optimized"
df22$type <- "Optimized"
df33$type <- "Optimized"
df44$type <- "Optimized"

dt11$type <- "Random"
dt22$type <- "Random"
dt33$type <- "Random"
dt44$type <- "Random"

di11$type <- "Typical"
di22$type <- "Typical"
di33$type <- "Typical"
di44$type <- "Typical"

df11$objective <- "Incidence"
dt11$objective <- "Incidence"
di11$objective <- "Incidence"

df22$objective <- "Hospitalizations"
dt22$objective <- "Hospitalizations"
di22$objective <- "Hospitalizations"

df33$objective <- "Deaths"
dt33$objective <- "Deaths"
di33$objective <- "Deaths"

df44$objective <- "DALYs"
dt44$objective <- "DALYs"
di44$objective <- "DALYs"

####################################
##combine dataframe
d1 <- rbind(df11,dt11,di11)
d2 <- rbind(df22,dt22,di22)
d3 <- rbind(df33,dt33,di33)
d4 <- rbind(df44,dt44,di44)

d1$value <- d1$value/1e6
d1$age_group <- gsub('age', '', d1$age_group)
d1$age_group <- paste(d1$age_group, "-", sep="")
d1$age_group <- factor(d1$age_group, levels = c("0-"	,"1-"	,"5-"	,"10-"	,"15-"	,"20-"	,"25-"	,"30-",
                                                "35-"	,"40-"	,"45-"	,"50-"	,"55-"	,"60-"	,"65-"	,"70-"	,"75-"))
d1$type <- factor(d1$type, levels = c("Random",  "Typical", "Optimized"))
d1$scenario <- d1$scenario*100

d2$value <- d2$value/1e3
d2$age_group <- gsub('age', '', d2$age_group)
d2$age_group <- paste(d2$age_group, "-", sep="")
d2$age_group <- factor(d2$age_group, levels = c("0-"	,"1-"	,"5-"	,"10-"	,"15-"	,"20-"	,"25-"	,"30-",
                                                "35-"	,"40-"	,"45-"	,"50-"	,"55-"	,"60-"	,"65-"	,"70-"	,"75-"))
d2$type <- factor(d2$type, levels = c("Random",  "Typical", "Optimized"))
d2$scenario <- d2$scenario*100


d3$value <- d3$value/1e3
d3$age_group <- gsub('age', '', d3$age_group)
d3$age_group <- paste(d3$age_group, "-", sep="")
d3$age_group <- factor(d3$age_group, levels = c("0-"	,"1-"	,"5-"	,"10-"	,"15-"	,"20-"	,"25-"	,"30-",
                                                "35-"	,"40-"	,"45-"	,"50-"	,"55-"	,"60-"	,"65-"	,"70-"	,"75-"))
d3$type <- factor(d3$type, levels = c("Random",  "Typical", "Optimized"))
d3$scenario <- d3$scenario*100



d4$value <- d4$value/1e6
d4$age_group <- gsub('age', '', d4$age_group)
d4$age_group <- paste(d4$age_group, "-", sep="")
d4$age_group <- factor(d4$age_group, levels = c("0-"	,"1-"	,"5-"	,"10-"	,"15-"	,"20-"	,"25-"	,"30-",
                                                "35-"	,"40-"	,"45-"	,"50-"	,"55-"	,"60-"	,"65-"	,"70-"	,"75-"))
d4$type <- factor(d4$type, levels = c("Random",  "Typical", "Optimized"))
d4$scenario <- d4$scenario*100
#####################################

p1 <- ggplot(d1, aes(x=age_group,y=scenario, fill = value))+
  geom_tile() + ylab("Scenario") + xlab("Age group")+ 
  facet_grid(.~type)+
  scale_fill_viridis(option = "viridis", limits =c(0,15))+
  scale_y_continuous(breaks=c(0, 20, 40,  60, 80,100))+
  theme(legend.title=element_blank())+ggtitle("Incidence (millions)")+
  theme(plot.title = element_text(size=17))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.y = element_text(size=18))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.x = element_text(size=15, angle=45))+
  guides(fill= guide_colorbar(barheight=10))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())

p2 <- ggplot(d2, aes(x=age_group,y=scenario, fill = value))+
  geom_tile() + ylab("Scenario") + xlab("Age group")+ 
  facet_grid(.~type)+
  scale_y_continuous(breaks=c(0, 20, 40,  60, 80,100))+
  scale_fill_gradientn(colours=rev(brewer.pal(9, "RdPu")),limits = c(0,120))+
  theme(legend.title=element_blank())+ggtitle("Hospitalizations (thousands)")+
  theme(plot.title = element_text(size=17))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.y = element_text(size=18))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.x = element_text(size=15, angle=45))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  guides(fill= guide_colorbar(barheight=10))

p3 <- ggplot(d3, aes(x=age_group,y=scenario, fill = value))+
  geom_tile() + ylab("High efficacy vaccine distribution (%)") + xlab("Age group")+ 
  facet_grid(.~type)+
  scale_fill_viridis(option = "inferno", limits = c(0,15))+
  scale_y_continuous(breaks=c(0, 20, 40,  60, 80,100))+
  theme(legend.title=element_blank())+ggtitle("Deaths (thousands)")+
  theme(plot.title = element_text(size=17))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.y = element_text(size=18))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.x = element_text(size=15, angle=45))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  guides(fill= guide_colorbar(barheight=10))

p4 <- ggplot(d4, aes(x=age_group,y=scenario, fill = value))+
  geom_tile() + ylab("High efficacy vaccine distribution (%)") + xlab("Age group")+ 
  facet_grid(.~type)+
  scale_fill_viridis(option = "cividis", limits = c(0,0.8))+
  scale_y_continuous(breaks=c(0, 20, 40,  60, 80,100))+
  theme(legend.title=element_blank())+ggtitle("DALY (millions)")+
  theme(plot.title = element_text(size=17))+
  theme(axis.title.y=element_blank())+
  theme(axis.text.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=15, angle=45))+
  theme(axis.title.x = element_text(size=18))+
  guides(fill= guide_colorbar(barheight=10))

p <- grid.arrange(p1,p2,p3,p4, ncol=1, heights = c(1,1,1,1.2))
g <- arrangeGrob(p,  left=grid::textGrob("High efficacy vaccine distribution (%)", rot = 90, vjust = 1,gp = gpar(fontsize = 18)))
ggsave("../plots/agewise_impact.png", plot=g, width = 12, height = 14)
#ggsave("../plots/agewise_impact.pdf", plot=p, width = 9, height = 12)
