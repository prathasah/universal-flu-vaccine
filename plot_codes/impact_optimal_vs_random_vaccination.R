library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(lattice)
library(reshape)
library(RColorBrewer)   # for brewer.pal(...)
library(cowplot)
library(Hmisc)
library(tidyr)
library('stringr')
library(dplyr)
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
df1 <-df1[,c("iter1", "scenario", "objective")]
df2 <-df2[,c("iter1", "scenario", "objective")]
df3 <-df3[,c("iter1", "scenario", "objective")]
df4 <-df4[,c("iter1", "scenario", "objective")]

df11 <- data.frame(df1 %>% group_by(scenario) %>% summarise(Incidence = mean(objective),  se = sd(objective)/sqrt(100)))
df22 <- data.frame(df2 %>% group_by(scenario) %>% summarise(Hospitalization = mean(objective),  se = sd(objective/sqrt(100))))
df33 <- data.frame(df3 %>% group_by(scenario) %>% summarise(Mortality = mean(objective),  se = sd(objective/sqrt(100))))
df44 <- data.frame(df4 %>% group_by(scenario) %>% summarise(DALY = mean(objective),  se = sd(objective/sqrt(100))))

df11$type <- "Optimized"
df22$type <- "Optimized"
df33$type <- "Optimized"
df44$type <- "Optimized"
#########################################
str(dt4)
dt1 <- dt1[,c("iter", "scenario", "total_infections")]
dt2 <- dt2[,c("iter", "scenario", "total_hospitalizations")]
dt3 <- dt3[,c("iter", "scenario", "total_mortality")]
dt4 <- dt4[,c("iter", "scenario", "total_DALY")]

dt11 <- data.frame(dt1 %>% group_by(scenario) %>% summarise(Incidence = mean(total_infections),  se = sd(total_infections/sqrt(1000))))
dt22 <- data.frame(dt2 %>% group_by(scenario) %>% summarise(Hospitalization = mean(total_hospitalizations),  se = sd(total_hospitalizations/sqrt(1000))))
dt33 <- data.frame(dt3 %>% group_by(scenario) %>% summarise(Mortality = mean(total_mortality),  se = sd(total_mortality/sqrt(1000))))
dt44 <- data.frame(dt4 %>% group_by(scenario) %>% summarise(DALY = mean(total_DALY),  se = sd(total_DALY/sqrt(1000))))

dt11$type <- "Random"
dt22$type <- "Random"
dt33$type <- "Random"
dt44$type <- "Random"
#########################################
str(dt4)
di1 <- di1[,c("iter", "scenario", "total_infections")]
di2 <- di2[,c("iter", "scenario", "total_hospitalizations")]
di3 <- di3[,c("iter", "scenario", "total_mortality")]
di4 <- di4[,c("iter", "scenario", "total_DALY")]

di11 <- data.frame(di1 %>% group_by(scenario) %>% summarise(Incidence = mean(total_infections),  se = sd(total_infections/sqrt(1000))))
di22 <- data.frame(di2 %>% group_by(scenario) %>% summarise(Hospitalization = mean(total_hospitalizations),  se = sd(total_hospitalizations/sqrt(1000))))
di33 <- data.frame(di3 %>% group_by(scenario) %>% summarise(Mortality = mean(total_mortality),  se = sd(total_mortality/sqrt(1000))))
di44 <- data.frame(di4 %>% group_by(scenario) %>% summarise(DALY = mean(total_DALY),  se = sd(total_DALY/sqrt(1000))))

di11$type <- "Typical"
di22$type <- "Typical"
di33$type <- "Typical"
di44$type <- "Typical"
###################################
d1 <- rbind(df11, dt11, di11)
d2 <- rbind(df22, dt22, di22)
d3 <- rbind(df33, dt33, di33)
d4 <- rbind(df44, dt44, di44)

us_population_2016 <- 323127513

d1$scenario <- d1$scenario *100
d1$Incidence <- d1$Incidence/ 1e6
d1$se <- d1$se/1e6
d1$type <- factor(d1$type, levels = c("Random", "Typical",  "Optimized"))

d2$scenario <- d2$scenario *100
d2$Hospitalization <- d2$Hospitalization / 1e3
d2$se <- d2$se/1e3
d2$type <- factor(d2$type, levels = c("Random", "Typical", "Optimized"))

d3$scenario <- d3$scenario *100
d3$Mortality <- d3$Mortality / 1e3
d3$se <- d3$se/1e3
d3$type <- factor(d3$type, levels = c("Random", "Typical", "Optimized"))

d4$scenario <- d4$scenario *100
d4$DALY <- d4$DALY / 1e6
d4$se <- d4$se/1e6
d4$type <- factor(d4$type, levels = c("Random",  "Typical", "Optimized"))
###########################

p1 <- ggplot(d1, aes(x=scenario, y=Incidence, fill = type))+
  geom_bar(position=position_dodge(),stat="identity",  color = "black", width = 10)+
  geom_errorbar(aes(ymin = Incidence -se, 
                    ymax = Incidence+se), width=0.5, position=position_dodge(9))+
  scale_x_continuous(breaks=c(0, 20, 40,  60, 80,100))+
  scale_fill_manual(values=c("#66c2a5", "#8da0cb", "#fc8d62"))+
  ggtitle("A")+ 
  theme(plot.title = element_text(hjust=-0.15))+
  theme(legend.position = c(0.65,0.9), legend.title = element_blank())+
  theme(axis.title.x=element_blank())+
  ylab("Infections (millions)")
p1

p2 <- ggplot(d2, aes(x=scenario, y=Hospitalization, fill = type))+
  geom_bar(position=position_dodge(),stat="identity",  color = "black", width = 10)+
  geom_errorbar(aes(ymin = Hospitalization -se, 
                    ymax = Hospitalization+se), width=0.5, position=position_dodge(9))+
  scale_x_continuous(breaks=c(0, 20, 40,  60, 80,100))+
  scale_fill_manual(values=c("#66c2a5", "#8da0cb", "#fc8d62"))+
  ggtitle("B")+ 
  theme(plot.title = element_text(hjust=-0.15))+
  guides(fill=FALSE)+
  theme(axis.title.x=element_blank())+
  ylab("Hospitalizations (thousands)")
p2


p3 <- ggplot(d3, aes(x=scenario, y=Mortality, fill = type))+
  geom_bar(position=position_dodge(),stat="identity",  color = "black", width = 10)+
  geom_errorbar(aes(ymin = Mortality -se, 
                    ymax = Mortality+se), width=0.5, position=position_dodge(9))+
  scale_x_continuous(breaks=c(0, 20, 40,  60, 80,100))+
  scale_fill_manual(values=c("#66c2a5", "#8da0cb", "#fc8d62"))+
  ggtitle("C")+ 
  theme(plot.title = element_text(hjust=-0.15))+
  guides(fill=FALSE)+
  theme(axis.title.x=element_blank())+
  ylab("Deaths (thousands)")
p3

p4 <- ggplot(d4, aes(x=scenario, y=DALY, fill = type))+
  geom_bar(position=position_dodge(),stat="identity",  color = "black", width = 10)+
  geom_errorbar(aes(ymin = DALY -se, 
                    ymax = DALY+se), width=0.5, position=position_dodge(9))+
  scale_x_continuous(breaks=c(0, 20, 40,  60, 80,100))+
  scale_fill_manual(values=c("#66c2a5", "#8da0cb", "#fc8d62"))+
  ggtitle("D")+ 
  theme(plot.title = element_text(hjust=-0.15))+
  guides(fill=FALSE)+
  theme(axis.title.x=element_blank())+
  ylab("DALYs (millions)")
p4

p <- grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)
g <- arrangeGrob(p,  bottom=grid::textGrob("High efficacy vaccine distribution (%)",gp = gpar(fontsize = 14)))
ggsave(filename="../plots/optimization_impact_population_level.png", plot=g, width = 7, height = 7)
#ggsave(filename="../plots/optimization_impact_population_level.pdf", plot=g, width = 7, height = 7)


