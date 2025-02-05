library(ggplot2)
library(rstatix)
library(tidyverse)
library(ggpubr)

# Load in data - change path here to Priming.csv

data <- read.csv("Priming.csv", header = T)
data$Sub <- as.factor(data$Sub)
data$Rep <- as.factor(data$Rep)
data$Metal <- as.factor(data$Metal)

level_order1 <- c('control', 'Zn', 'Cd') 

ggplot(data = data,
       aes(y=averageGrowth, x=Sub, 
           colour = factor(Metal, level = level_order1), group = factor(Metal, level = level_order1), 
           ymin = averageGrowth-stdev, 
           ymax = averageGrowth+stdev)) +
  geom_point(stat="summary", size = 3) +
  geom_line(stat="summary") +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", color = "black", width = 0.01) +
  ggtitle("Growth during experiment") +
  theme(plot.title = element_text(face = "bold.italic"), 
        axis.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1.25, 'cm'),
        legend.text = element_text(size=12)) +
  scale_color_manual(values=c('#00BFC4', '#00BA38', '#C77CFF'))

## graph for each replicate split up
ggplot(data = data,
       aes(y=averageGrowth, x=Sub, 
           group = interaction(Metal, Rep), color = Metal, shape = Metal)
) +
  geom_point(stat="identity", size = 3) +
  geom_line(stat="identity") +
  geom_errorbar(
    aes(
      ymin = averageGrowth-stdev, 
      ymax = averageGrowth+stdev,
      width = 0.2)) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", color = "black", width = 0.01) +
  ggtitle("Growth during experiment") +
  theme(plot.title = element_text(face = "bold.italic"), 
        axis.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1.25, 'cm'),
        legend.text = element_text(size=12)) +
  scale_color_manual(values=c('#00BFC4', '#00BA38', '#C77CFF'))


## compare week 1 growth to week 10 growth

versusGrowth <- data[data$Sub=="1" | data$Sub=="10",]
versusGrowth <- mutate(versusGrowth, St = case_when(Sub == 1 ~ "St.0",
                              Sub == 10 & Metal == "control" ~ "St.H2O",
                              Sub == 10 & Metal == "Zn" ~ "St.Zn",
                              Sub == 10 & Metal == "Cd" ~ "St.Cd"))
versusGrowth$St <- as.factor(versusGrowth$St)

ggplot(versusGrowth, aes(group = St, fill=St, y=averageGrowth,
                                                 x=factor(Metal, level = level_order1))) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width=0.9)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width=0.9),  
               geom = "errorbar", color = "black", width = 0.05) +
  ggtitle("Growth - Week 1 vs. Week 10") +
  ylim(0, 1.75) +
  theme(plot.title = element_text(face = "bold.italic"), 
        axis.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(1.25, 'cm'),
        legend.text = element_text(size=12)) +
  scale_fill_manual(breaks = c("St.0","St.H2O","St.Zn","St.Cd"),
                     values=c('#F8766D','#00BFC4','#00BA38','#C77CFF'))

## STATS
dataforstats <- data[data$Sub=="1" | data$Sub=="10",]
#Summary stats
dataforstats %>%
  group_by(interaction(Metal, Sub)) %>%
  get_summary_stats(averageGrowth, type = "mean_sd")
#Check for outliers
dataforstats %>% 
  group_by(interaction(Metal, Sub)) %>%
  identify_outliers(averageGrowth)
#Check for normality
# Build the linear model
model  <- lm(averageGrowth ~ interaction(Metal, Sub), data = dataforstats)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
# p value NS = normality
shapiro_test(residuals(model))
#Levene’s test to check the homogeneity of variances
# p values NS = homogeneity
dataforstats %>% levene_test(averageGrowth ~ interaction(Metal, Sub))
#ANOVA
res.aov <- dataforstats %>% anova_test(averageGrowth ~ Metal * Sub)
res.aov
# Pairwise comparisons
pwc <- dataforstats %>% tukey_hsd(averageGrowth ~ Metal * Sub)
pwc
pwc %>% print(n=470)

# what about t-tests for the before and after comparisons??
dataforttest <- data[data$Sub=="1" & data$Metal=="Cd" | data$Sub=="10" & data$Metal=="Cd",]

# Compute t-test
res <- t.test(averageGrowth ~ Sub, data = dataforttest)
res



## Input data from crosschecks - change path here to PrimingCrosschecks.csv

xdata <- read.csv("PrimingCrosschecks.csv", header = T)
crossdata <- xdata[!xdata$Isolate=="week1",]
crossdata$Isolate <- as.factor(crossdata$Isolate)
crossdata$Rep <- as.factor(crossdata$Rep)
crossdata$Metal <- as.factor(crossdata$Metal)

level_order <- c('week1', 'control', 'zinc', 'cadmium') 

# this plot shows the growth of each isolate (either isoalte cultured on 0, zn or Cd for 10 wks)
# showing the growth of these isolates on each metal
# includes week one as control control
ggplot(xdata, aes(fill=factor(Metal, level = level_order1), y=averageGrowth, x=factor(Isolate, level = level_order))) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width=0.9)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width=0.9),  
               geom = "errorbar", color = "black", width = 0.05) +
  ggtitle("Tolerance crosschecks") +
  theme(plot.title = element_text(face = "bold.italic"), 
        axis.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(1.25, 'cm'),
        legend.text = element_text(size=12))

# plot for each metal (plate)

ggplot(xdata[xdata$Metal=="control",], aes(fill=factor(Isolate, level = level_order), y=averageGrowth, 
                                           x=factor(Isolate, level = level_order))) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width=1)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width=1),  
               geom = "errorbar", color = "black", width = 0.05) +
  ggtitle("Growth of subcultured isolates on control Fries Media") +
  ylim(0,1.75) +
  theme(plot.title = element_text(face = "bold.italic"), 
        axis.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(1.25, 'cm'),
        legend.text = element_text(size=12)) +
  scale_fill_manual(values=c('#F8766D','#00BFC4','#00BA38', '#C77CFF'))

ggplot(xdata[xdata$Metal=="Zn",], aes(fill=factor(Isolate, level = level_order), y=averageGrowth, 
                                      x=factor(Isolate, level = level_order))) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width=1)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width=1),  
               geom = "errorbar", color = "black", width = 0.05) +
  ggtitle("Growth of subcultured isolates on Zinc") +
  ylim(0,1.5) +
  theme(plot.title = element_text(face = "bold.italic"), 
        axis.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(1.25, 'cm'),
        legend.text = element_text(size=12)) +
  scale_fill_manual(values=c('#F8766D','#00BFC4','#00BA38', '#C77CFF'))

ggplot(xdata[xdata$Metal=="Cd",], aes(fill=factor(Isolate, level = level_order), y=averageGrowth, 
                                      x=factor(Isolate, level = level_order))) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width=1)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width=1),  
               geom = "errorbar", color = "black", width = 0.05) +
  ggtitle("Growth of subcultured isolates on Cadmium") +
  ylim(0,1) +
  theme(plot.title = element_text(face = "bold.italic"), 
        axis.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(1.25, 'cm'),
        legend.text = element_text(size=12)) +
  scale_fill_manual(values=c('#F8766D','#00BFC4','#00BA38', '#C77CFF'))


## STATS
dataforstats <- xdata[xdata$Metal=="Cd",]
#Summary stats
dataforstats %>%
  group_by(Isolate) %>%
  get_summary_stats(averageGrowth, type = "mean_sd")
#Check for outliers
dataforstats %>% 
  group_by(Isolate) %>%
  identify_outliers(averageGrowth)
#Check for normality
# Build the linear model
model  <- lm(averageGrowth ~ Isolate, data = dataforstats)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
# p value NS = normality
shapiro_test(residuals(model))
#Levene’s test to check the homogeneity of variances
# p values NS = homogeneity
dataforstats %>% levene_test(averageGrowth ~ Isolate)
#ANOVA
res.aov <- dataforstats %>% anova_test(averageGrowth ~ Isolate)
res.aov
# Pairwise comparisons
pwc <- dataforstats %>% tukey_hsd(averageGrowth ~ Isolate)
pwc
pwc %>% print(n=470)

# Compute t-test
res <- t.test(averageGrowth ~ Isolate, 
              data = xdata[xdata$Metal=="Cd" & xdata$Isolate=="week1" | xdata$Isolate=="cadmium" & xdata$Metal=="Cd",])
res

xdata[xdata$Metal=="Cd" & xdata$Isolate=="week1" | xdata$Isolate=="cadmium" & xdata$Metal=="Cd",]

## get letter groups for comparison
mod_means_contr <- emmeans::emmeans(object = model,
                                    pairwise ~ "Isolate",
                                    adjust = "tukey")

mod_means <- multcomp::cld(object = mod_means_contr$emmeans,
                             Letters = letters)
mod_means_contr
mod_means
