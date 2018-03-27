library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(lattice)
library(reshape)
library(RColorBrewer)   # for brewer.pal(...)
library(cowplot)
library(Hmisc)
################################################################
demography <- c(3970145,3995008,3992154,3982074,3987656,4032515,4029655,4029991,4159114,4178524,4144019,
                 4131222,4139558,4109703,4093731,4196991,4265224,4205001,4219303,4243480,4286221,
                 4386854,4480904,4552952,4674097,4758352,4747254,4559206,4451507,4374565,4392155,
                 4423807,4283076,4345786,4341535,4293125,4375562,4103498,4022119,3979601,3862150,
                 3979596,3854553,3912307,4087645,4319616,4371961,4142964,4055074,4058008,4131293,4375892,
                 4452636,4445605,4433630,4484565,4518758,4362807,4332953,4281025,4123968,4083962,3911429,
                 3757382,3606295,3490890,3403647,3295266,3251936,3378344,2487211,2445650,2371252,2413647,
                 2092487,1900211,1784266,1664566,1578915,1439937,1358260,1284298,1135109,1079082,1008890,
                 938333,876190,766479,693882,610338,520906,449986,372625,300000,239313,186408,135797,94311,
                 68972,44895,81896)

######################################################
dt1_1 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.csv", header=T)
dt1_2 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.2.csv", header=T)
dt1_3 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.4.csv", header=T)
dt1_4 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.6.csv", header=T)
dt1_5 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_0.8.csv", header=T)
dt1_6 <- read.csv("../data/optimization_results/optimization_results_infections_1000iter_1.csv", header=T)

dt <- rbind(dt1_1,dt1_2, dt1_3,dt1_4,dt1_5, dt1_6)

str(dt)
#dt <- rbind(dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8)

#[0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
# group numbers 0: 0- , 1:1-, 2:5-, 3:10-, 4:15,5:20, 6:25,7:30, 8:35,9:40,10:45,11:50,12:55, 13:60,14:65,15:70: 16:75
#grouped_ages = {"0-": [1,2], "5-":[3,4,5], "20-":[6,7], "30-":[8,9], "40-":[10,11], "50-":[12,13,14], "65+":[15,16,17]}

#ages 0 to 4
df1 <- dt[,c( "scenario", "efficacy1", "efficacy2" ,  "max_doses1", "max_doses2" , "iter1", 
              "tot_low_efficacy_vac1", "prop_low_efficacy_vac1", "tot_low_efficacy_vac2", "prop_low_efficacy_vac2",
              "tot_high_efficacy_vac1", "prop_high_efficacy_vac1", "tot_high_efficacy_vac2", "prop_high_efficacy_vac2")]

df1$low_doses <- df1$tot_low_efficacy_vac1 + df1$tot_low_efficacy_vac2
df1$high_doses <- df1$tot_high_efficacy_vac1 + df1$tot_high_efficacy_vac2
df1$low_proportion_vaccinated <- df1$low_doses/ sum(demography[1:5])
df1$high_proportion_vaccinated <- df1$high_doses/ sum(demography[1:5])
df1 <- df1[,c("scenario", "iter1", "low_doses", "high_doses" ,"low_proportion_vaccinated", "high_proportion_vaccinated")]
df1$age_group <- "0-"
df1$age_group_xpos <- 2
df1$age_group_xtick <- 2.5

#ages 5 to 19
df2 <- dt[,c( "scenario", "efficacy1", "efficacy2" ,  "max_doses1", "max_doses2" , "iter1", 
              "tot_low_efficacy_vac3", "prop_low_efficacy_vac3", "tot_low_efficacy_vac4", "prop_low_efficacy_vac4", "tot_low_efficacy_vac5", "prop_low_efficacy_vac5",
              "tot_high_efficacy_vac3", "prop_high_efficacy_vac3", "tot_high_efficacy_vac4", "prop_high_efficacy_vac4", "tot_high_efficacy_vac5", "prop_high_efficacy_vac5")]

df2$low_doses <- df2$tot_low_efficacy_vac3 + df2$tot_low_efficacy_vac4 + df2$tot_low_efficacy_vac5
df2$high_doses <- df2$tot_high_efficacy_vac3 + df2$tot_high_efficacy_vac4 + df2$tot_high_efficacy_vac5
df2$low_proportion_vaccinated <- df2$low_doses/ sum(demography[6:20])
df2$high_proportion_vaccinated <- df2$high_doses/ sum(demography[6:20])
df2 <- df2[,c("scenario", "iter1", "low_doses", "high_doses" ,"low_proportion_vaccinated", "high_proportion_vaccinated")]
df2$age_group <- "5-"
df2$age_group_xpos <- 3
df2$age_group_xtick <- 4

#ages 20-29
df3 <- dt[,c( "scenario", "efficacy1", "efficacy2" ,  "max_doses1", "max_doses2" , "iter1", 
              "tot_low_efficacy_vac6", "prop_low_efficacy_vac6", "tot_low_efficacy_vac7", "prop_low_efficacy_vac7",
              "tot_high_efficacy_vac6", "prop_high_efficacy_vac6", "tot_high_efficacy_vac7", "prop_high_efficacy_vac7")]

df3$low_doses <- df3$tot_low_efficacy_vac6 + df3$tot_low_efficacy_vac7
df3$high_doses <- df3$tot_high_efficacy_vac6 + df3$tot_high_efficacy_vac7
df3$low_proportion_vaccinated <- df3$low_doses/ sum(demography[21:30])
df3$high_proportion_vaccinated <- df3$high_doses/ sum(demography[21:30])
df3 <- df3[,c("scenario", "iter1", "low_doses", "high_doses" ,"low_proportion_vaccinated", "high_proportion_vaccinated")]
df3$age_group <- "20-"
df3$age_group_xpos <- 6
df3$age_group_xtick <- 6.5

#ages 30-39
df4 <- dt[,c( "scenario", "efficacy1", "efficacy2" ,  "max_doses1", "max_doses2" , "iter1", 
              "tot_low_efficacy_vac8", "prop_low_efficacy_vac8", "tot_low_efficacy_vac9", "prop_low_efficacy_vac9",
              "tot_high_efficacy_vac8", "prop_high_efficacy_vac8", "tot_high_efficacy_vac9", "prop_high_efficacy_vac9")]

df4$low_doses <- df4$tot_low_efficacy_vac8 + df4$tot_low_efficacy_vac9
df4$high_doses <- df4$tot_high_efficacy_vac8 + df4$tot_high_efficacy_vac9
df4$low_proportion_vaccinated <- df4$low_doses/ sum(demography[31:40])
df4$high_proportion_vaccinated <- df4$high_doses/ sum(demography[31:40])
df4 <- df4[,c("scenario", "iter1", "low_doses", "high_doses" ,"low_proportion_vaccinated", "high_proportion_vaccinated")]
df4$age_group <- "30-"
df4$age_group_xpos <- 8
df4$age_group_xtick <- 8.5

#ages 40-49
df5 <- dt[,c( "scenario", "efficacy1", "efficacy2" ,  "max_doses1", "max_doses2" , "iter1", 
              "tot_low_efficacy_vac10", "prop_low_efficacy_vac10", "tot_low_efficacy_vac11", "prop_low_efficacy_vac11",
              "tot_high_efficacy_vac10", "prop_high_efficacy_vac10", "tot_high_efficacy_vac11", "prop_high_efficacy_vac11")]

df5$low_doses <- df5$tot_low_efficacy_vac10 + df5$tot_low_efficacy_vac11
df5$high_doses <- df5$tot_high_efficacy_vac10 + df5$tot_high_efficacy_vac11
df5$low_proportion_vaccinated <- df5$low_doses/ sum(demography[41:50])
df5$high_proportion_vaccinated <- df5$high_doses/ sum(demography[41:50])
df5 <- df5[,c("scenario", "iter1", "low_doses", "high_doses" ,"low_proportion_vaccinated", "high_proportion_vaccinated")]
df5$age_group <- "40-"
df5$age_group_xpos <- 10
df5$age_group_xtick <- 10.5

#ages 50-64
df6 <- dt[,c( "scenario", "efficacy1", "efficacy2" ,  "max_doses1", "max_doses2" , "iter1", 
              "tot_low_efficacy_vac12", "prop_low_efficacy_vac12", "tot_low_efficacy_vac13", "prop_low_efficacy_vac13", "tot_low_efficacy_vac14", "prop_low_efficacy_vac14",
              "tot_high_efficacy_vac12", "prop_high_efficacy_vac12", "tot_high_efficacy_vac13", "prop_high_efficacy_vac13", "tot_high_efficacy_vac14", "prop_high_efficacy_vac14")]

df6$low_doses <- df6$tot_low_efficacy_vac12 + df6$tot_low_efficacy_vac13 + df6$tot_low_efficacy_vac14
df6$high_doses <- df6$tot_high_efficacy_vac12 + df6$tot_high_efficacy_vac13 + df6$tot_high_efficacy_vac14
df6$low_proportion_vaccinated <- df6$low_doses/ sum(demography[51:65])
df6$high_proportion_vaccinated <- df6$high_doses/ sum(demography[51:65])
df6 <- df6[,c("scenario", "iter1", "low_doses", "high_doses" ,"low_proportion_vaccinated", "high_proportion_vaccinated")]
df6$age_group <- "50-"
df6$age_group_xpos <- 12
df6$age_group_xtick <- 13


#ages 65+
df7 <- dt[,c( "scenario", "efficacy1", "efficacy2" ,  "max_doses1", "max_doses2" , "iter1", 
              "tot_low_efficacy_vac15", "prop_low_efficacy_vac15", "tot_low_efficacy_vac16", "prop_low_efficacy_vac16", "tot_low_efficacy_vac17", "prop_low_efficacy_vac17",
              "tot_high_efficacy_vac15", "prop_high_efficacy_vac15", "tot_high_efficacy_vac16", "prop_high_efficacy_vac16", "tot_high_efficacy_vac17", "prop_high_efficacy_vac17")]

df7$low_doses <- df7$tot_low_efficacy_vac15 + df7$tot_low_efficacy_vac16 + df7$tot_low_efficacy_vac17
df7$high_doses <- df7$tot_high_efficacy_vac15 + df7$tot_high_efficacy_vac16 + df7$tot_high_efficacy_vac17
df7$low_proportion_vaccinated <- df7$low_doses/ sum(demography[66:length(demography)])
df7$high_proportion_vaccinated <- df7$high_doses/ sum(demography[66:length(demography)])
df7 <- df7[,c("scenario", "iter1", "low_doses", "high_doses" ,"low_proportion_vaccinated", "high_proportion_vaccinated")]
df7$age_group <- "65-"
df7$age_group_xpos <- 15
df7$age_group_xtick <- 16


df <- rbind(df1,df2,df3,df4,df5,df6,df7)
str(df)
dt1 <- df[,c("scenario", "iter1", "age_group" , "age_group_xpos", "age_group_xtick",  "low_doses", "low_proportion_vaccinated")]
dt2 <- df[,c("scenario", "iter1",  "age_group" , "age_group_xpos", "age_group_xtick", "high_doses", "high_proportion_vaccinated")]
dt1$type <- "low_efficacy"
dt2$type <- "high_efficacy"
str(dt1) ####PS rename low and high to combine into 1 dataframe
colnames(dt1)[6:7] <- c("doses", "proportion_vaccinated")
colnames(dt2)[6:7] <- c("doses", "proportion_vaccinated")

dt <- rbind(dt1,dt2)
dt$age_group <- factor(dt$age_group, levels = c("0-", "5-", "20-", "30-", "40-", "50-", "65-"))
dt$doses <- dt$doses/1000000
dt$percentage_vaccinated <-dt$proportion_vaccinated * 100

summary(factor(dt$scenario))
dt1<- aggregate(x= dt$doses, list(dt$scenario, dt$type, dt$age_group, dt$age_group_xpos, dt$age_group_xtick), mean)
dt2<- aggregate(x= dt$percentage_vaccinated, list(dt$scenario, dt$type, dt$age_group, dt$age_group_xpos, dt$age_group_xtick), mean)



str(dt1)
colnames(dt1)[1] <- "scenario"
colnames(dt1)[2] <- "type"
colnames(dt1)[3] <- "age_group"
colnames(dt1)[4] <- "age_group_xpos"
colnames(dt1)[5] <- "age_group_xtick"
colnames(dt1)[6] <- "doses"
write.csv(dt1, "../data/optimization_results/cleaned_optimization_infection_response_doses_100iter.csv")


colnames(dt2)[1] <- "scenario"
colnames(dt2)[2] <- "type"
colnames(dt2)[3] <- "age_group"
colnames(dt2)[4] <- "age_group_xpos"
colnames(dt2)[5] <- "age_group_xtick"
colnames(dt2)[6] <- "percentage_vaccinated"


write.csv(dt2, "../data/optimization_results/cleaned_optimization_infection_response_coverage_100iter.csv")
##################################################################################
