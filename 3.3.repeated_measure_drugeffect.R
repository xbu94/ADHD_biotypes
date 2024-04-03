library(readr)
library(tidyr)
library(ggplot2)
library(car)
library(dplyr)
library(HH)
setwd("D:/Helab/ADHDsubtyping/Stats/subtyping_nm_gender/subgroup_analysis_nmgender/drugeffect")

## Before-after treatment response by repeated measure ANOVA
# full model: 3 factor ANOVA
mod_aov = aov(score~drug*time*subgroups+Error(id/time),data=T_long)
summary(mod_aov) 

# within subtype1
mod1 = aov(score~drug*time+Error(id/time),data=subset(T_long,subgroups == 'subtype1'))
summary(mod1)
interaction2wt(score ~ drug*time, data = subset(T_long,subgroups == 'subtype1'))
tiff('./fig/Drugeffets_total_subtype1.tif',units="cm",res=200,height=10,width=12)
with(subset(T_long,subgroups == 'subtype1'),interaction.plot(time,drug,score,
                                                              type='b',col=c('slategray4','darkgrey'),pch=c(16,18),lwd = 3,
                                                              legend = F,trace.label = 'drug',
                                                             ylim = range(15,35),
                                                              main='drug effets on total for subtype1'))
dev.off()

mod1.1 = aov(score~time+Error(id/time),data=subset(T_long,subgroups == 'subtype1'|drug == 'MPH'))
summary(mod1.1)
mod1.2 = aov(score~time+Error(id/time),data=subset(T_long,subgroups == 'subtype1'|drug == 'ATX'))
summary(mod1.2)
summary(aov(score~time,data=subset(T_long,subgroups == 'subtype1'|drug == 'MPH')))
summary(aov(score~time,data=subset(T_long,subgroups == 'subtype1'|drug == 'ATX')))


# within sub-type2
mod2 = aov(score~drug*time+Error(id/time),data=subset(T_long,subgroups == 'subtype2'))
summary(mod2)
interaction2wt(score ~ time*drug, data = subset(T_long,subgroups == 'subtype2'))
tiff('./fig/Drugeffets_total_subtype2.tif',units="cm",res=200,height=10,width=12)
with(subset(T_long,subgroups == 'subtype2'),interaction.plot(time,drug,score,
                                                    type='b',col=c('slategray4','darkgrey'),pch=c(16,18),lwd = 3,
                                                    legend = F,trace.label = 'drug',
                                                    ylim = range(15,35),
                                                   main='drug effets on symptoms for subtype2'))
dev.off()
mod2.1 = aov(score~time+Error(id/time),data=subset(T_long,subgroups == 'subtype2'|drug == 'MPH'))
summary(mod2.1)
mod2.2 = aov(score~time+Error(id/time),data=subset(T_long,subgroups == 'subtype2'|drug == 'ATX'))
summary(mod2.2)

