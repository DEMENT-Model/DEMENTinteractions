library(ggplot2)
library(cowplot)
library(tidyverse)
## Looking at interdependence of microbes
setwd("/Users/brittnibertolet/OneDrive - UC Irvine/Postdoc stuff/Projects/the untouched DEMENTpy-master/working_DEMENT/")
setwd("output/output-20231019/Data")

#### Read in carbon data ####

# Read in the carbon data frame 
carbon=read.csv("carbonDF.csv", stringsAsFactors = F)
# Get only year 3 
carbon=carbon[carbon$time>730 & carbon$time<1096,]

#### Calculate percent change in substrate degraded #### 
# Get substrate beginning of year 3: time = 731
carbon730=carbon[carbon$time==731,]
carbon730$substrate730=carbon730$substrate
# Get % substrate degraded at the end of year 3: time = 1095
carbon1095=carbon[carbon$time==1095,]
carbon1095$substrate1095=carbon1095$substrate

# Merge the two
temp=merge(carbon1095[,c("file", "substrate1095")],carbon730[,c("file", "substrate730")], by="file")
temp$percSubDegraded=(1-(temp$substrate1095/temp$substrate730))*100

# Merge back with 1096 
carbon1095$substrate1095=NULL
carbon1095=merge(carbon1095, temp, by="file")


### Calculate interdependence at 1095
allTaxa=unique(carbon1095$file[!grepl("exclude",carbon1095$file)])
carbon1095all=carbon1095[carbon1095$file%in%allTaxa,]
carbon1095$taxaEffect=NA
strings=gsub("_all.pickle", "", allTaxa)
# Percent change 
for(i in 1:length(strings)){
  carbon1095$taxaEffect[grepl(strings[i],carbon1095$file)]=
    (carbon1095all$percSubDegraded[grepl(strings[i] ,carbon1095all$file)]-
       carbon1095$percSubDegraded[grepl(strings[i],carbon1095$file)])/
    carbon1095all$percSubDegraded[grepl(strings[i] ,carbon1095all$file)]*100
}

# Get rid of effect of their own exclusion
length(carbon1095$taxaEffect[carbon1095$taxaEffect==0])
carbon1095=carbon1095[carbon1095$taxaEffect!=0,]

# Get community and population trait data
traits=read.csv("traitsDF.csv", stringsAsFactors = F)
# Get Taxa ID
traits$taxaID=gsub("_all.pickle", "", traits$file)
traits$taxaID=paste(traits$taxaID, traits$Taxa, sep="_")
# Get Taxa ID for carbon 1095 data
temp=as.data.frame(str_split_fixed(carbon1095$file, '_', 5)[,1:3])
carbon1095$taxaID=paste(temp$V1, temp$V2, temp$V3, sep="_")
temp=str_split_fixed(carbon1095$file, '_', 5)[,5]
temp=paste0("Tax",as.numeric(str_split_fixed(temp, '\\.', 2)[,1])+1)
carbon1095$taxaID=paste(carbon1095$taxaID, temp, sep="_")


# Calculation carbon allocation
# Look at allocation to osmolyte
# Drought tolerance is the relative allocation compared to other taxa (0-1)
traits$Osmo_Alloc=traits$Osmo_Consti_Cost+traits$Osmo_Induci_Cost
# Look at allocation to enzymes
traits$Enz_Alloc=(traits$Enz_Consti_Cost+traits$Enz_Induci_Cost)*traits$Enz_Gene
# Look at allocation to uptake
traits$Uptake_Alloc=(traits$Uptake_Cost)*traits$Uptake_Gene
# Calculate ratio
traits$ratioEU=traits$Enz_Alloc/traits$Uptake_Alloc

#Get average weighted traits of the community 
comAvg=data.frame()
seeds=unique(carbon1095$microbe)
for(i in 1:length(seeds)){
  temp=traits[carbon1095$microbe==seeds[i],]
  temp=temp[grepl("grass", temp$file),]
  temp3=data.frame(microbe=seeds[i],
                   avgEnz_Gene=weighted.mean(temp$Enz_Gene, temp$Avg_Abund),
                   avgEnz_Alloc=weighted.mean(temp$Enz_Alloc, temp$Avg_Abund),
                   avgDrought_Tol=weighted.mean(temp$Drought_tolerance, temp$Avg_Abund),
                   avgOsmo_Alloc=weighted.mean(temp$Osmo_Alloc, temp$Avg_Abund),
                   avgUptake_Alloc=weighted.mean(temp$Uptake_Alloc, temp$Avg_Abund),
                   avgRation=weighted.mean(temp$ratioEU, temp$Avg_Abund)
                   
  )
  comAvg=rbind(comAvg, temp3)
}

# Merge trait data and carbon data
carbon1095=merge(carbon1095, traits, by=c("taxaID"))
# Merge with carbon data
carbon1095=merge(carbon1095, comAvg, by="microbe")


# Plot 3A - Ecosystem-scale population impact
# Fix climate factor order
carbon1095$climate=factor(carbon1095$climate, labels=c("Ambient", "Drought", "Moist"))

carbon1095$climate=factor(carbon1095$climate, levels=c("Moist", "Ambient", "Drought"))
# Add a factor for faceting to indicate ecosystem-scale 
carbon1095$eco="Ecosystem-scale impact"

# Plot 2A
carbon1095$avgEnz_Alloc2=paste(round(carbon1095$avgEnz_Alloc, digits=5))
plot2A=ggplot(carbon1095, aes(x=climate, y=taxaEffect))+
  geom_hline(yintercept = 0, color="black")+
  #geom_jitter(aes(color=avgEnz_Alloc2),alpha=0.6, size=1,width = 0.3)+
  geom_point(aes(color=avgEnz_Alloc2), position=position_dodge(width = .7), alpha=0.4, size=1)+
  #geom_violin()+
  #geom_boxplot(outlier.shape = NA, fill="NA")+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+theme(panel.grid = element_blank(),
                   strip.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.text.y = element_text(size = 8, angle=90, hjust = 0.5),
                   legend.title = element_text(size = 10), 
                   legend.text = element_text(size = 8),
                   legend.key.height = unit(0.4, 'cm'),
                   plot.margin = margin(0.5,0.3,0.1,0.3, "cm"))+
  ylab(expression(Delta*" substrate degradation (%)")) + xlab("Climate Scenario")+
  facet_grid(~eco)+
  labs(color="Community-average\nenzyme production")+
  NULL


# Are there differences across climates, treating taxa id as a random effect?
carbon1095$taxaID2=str_split_fixed(carbon1095$taxaID, '_', 2)[,2]
unique(carbon1095$taxaID2)

library(rstatix)
# Two-way mixed factors ANOVA
res.aov <- anova_test(
  data = carbon1095, dv = taxaEffect, wid = taxaID2,
  between = microbe, within = climate)
get_anova_table(res.aov)

# Post-hoc tests 
res.aov <- anova_test(
  data = carbon1095[carbon1095$climate%in%c("Ambient", "Drought"),], dv = taxaEffect, wid = taxaID2, within = climate)
get_anova_table(res.aov)

# Post-hoc tests 
res.aov <- anova_test(
  data = carbon1095[carbon1095$climate%in%c("Ambient", "Moist"),], dv = taxaEffect, wid = taxaID2, within = climate)
get_anova_table(res.aov)

# Post-hoc tests 
res.aov <- anova_test(
  data = carbon1095[carbon1095$climate%in%c("Drought", "Moist"),], dv = taxaEffect, wid = taxaID2, within = climate)
get_anova_table(res.aov)

# Look at simple correltations across the climate scenarios 
drought=carbon1095[carbon1095$climate=="Drought",]
drought=drought[order(drought$taxaID),]

ambient=carbon1095[carbon1095$climate=="Ambient",]
ambient=ambient[order(ambient$taxaID),]

moist=carbon1095[carbon1095$climate=="Moist",]
moist=moist[order(moist$taxaID),]

summary(lm(drought$taxaEffect~ambient$taxaEffect))
summary(lm(drought$taxaEffect~moist$taxaEffect))


#### Calculate pairwise interactions ####
# Read in microbial abundance data
microbe=read.csv("microbeDF.csv", stringsAsFactors = F)

# Get exclusion and full community seperately
microbe_full=microbe[grepl("all", microbe$file), ]
unique(microbe_full$file)
microbe_ex=microbe[!grepl("all", microbe$file), ]
pops=unique(microbe_ex$file)

# Get climate treatments 
climates=c("grass", "drought", "moist")

biomassChange=c()
for(i in 1:length(pops)){
  # Get the file of a single population being excluded
  temp=microbe_ex[microbe_ex$file==pops[i],]
  # Get average biomass of all other taxa in third year
  temp=temp[temp$X>730,]
  # Get rid of unnecessary columns
  temp$X=NULL
  temp$file=NULL
  # Get column averages 
  temp=colMeans(temp)
  
  # Now get full simulation
  fullSim=str_split_fixed(pops[i], '_', 4)[,3]
  cli=str_split_fixed(pops[i], '_', 4)[,1]
  tempFull=microbe_full[grepl(fullSim, microbe_full$file),]
  tempFull=tempFull[grepl(cli, tempFull$file),]
  # Get average biomass of all other taxa in third year
  tempFull=tempFull[tempFull$X>730,]
  # Get rid of unnecessary columns
  tempFull$X=NULL
  tempFull$file=NULL
  # Get column averages 
  tempFull=colMeans(tempFull)
  
  # Get effect of taxa presence on biomass of neighboring populations
  # Full communnity - community when excluded
  tempOut=data.frame(population=pops[i], taxa=names(temp), climate=cli, biomassImpact=tempFull-temp)
  
  # Concatonate to output
  biomassChange=rbind(biomassChange, tempOut)
}


# Get rid of NAs
biomassChange=biomassChange[complete.cases(biomassChange),]

# Clean up environment to open up memory
rm(carbon, carbon1095all, tempOut, cli, allTaxa, climates, fullSim, pops, temp, tempFull)
rm(microbe, microbe_ex, microbe_full)

# Plot 2B - Community-scale population impact
# Fix climate factor order
biomassChange$climate=factor(carbon1095$climate, labels=c("Drought", "Ambient", "Moist"))
biomassChange$climate=factor(carbon1095$climate, levels=c( "Moist", "Ambient", "Drought"))
# Add a factor for faceting to indicate ecosystem-scale 
biomassChange$com="Community-scale impact"

# Get microbe ID column
biomassChange$microbe=str_split_fixed(biomassChange$population, '_', 5)[,3]
temp=carbon1095[,c("microbe", "avgEnz_Alloc2")]
temp=temp[!duplicated(temp),]
biomassChange=left_join(biomassChange, temp, by="microbe")

plot2B=ggplot(biomassChange, aes(x=climate, y=biomassImpact))+
  geom_hline(yintercept = 0, color="grey", linetype="dashed")+
  geom_point(aes(color=avgEnz_Alloc2), position=position_dodge(width = .7), alpha=0.4, size=1)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+theme(panel.grid = element_blank(),
                   strip.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.text.y = element_text(size = 8, angle=90, hjust = 0.5),
                   legend.title = element_text(size = 10), 
                   legend.text = element_text(size = 8),
                   legend.key.height = unit(0.4, 'cm'),
                   #legend.position = c(0.2, 0.8), # c(0,0) bottom left, c(1,1) top-right.
                   legend.background = element_rect(fill = "transparent", colour = NA),
                   plot.margin = margin(0.5,0.3,0.1,0.3, "cm"))+
  ylab(expression(Delta*" neighbor biomass (mg C)")) + xlab("Climate Scenario")+
  facet_grid(~com)
plot2B


left=plot_grid(plot3A+guides(color="none"), plot2B+guides(color="none"), align="hv", axis="lr", labels=c("A", "B"))
plot_grid(left, get_legend(plot3A), rel_widths = c(1,0.3))
#ggsave("~/Desktop/DEMENT_Figures/Figures_20231115/Fig2.png", height=2.8, width=6.5)


# Two-way mixed factors ANOVA
biomassChange$population2=paste(str_split_fixed(biomassChange$population, '_', 2)[,2], biomassChange$taxa, sep="_")

res.aov <- anova_test(
  data = biomassChange, dv = biomassImpact, wid = population2,
  between = microbe, within = climate)
get_anova_table(res.aov)

res.aov <- anova_test(
  data = carbon1095, dv = taxaEffect, wid = taxaID2,
  between = microbe, within = climate)
get_anova_table(res.aov)



#### Relationship between traits ####
### Look at ambient data 
ambient=carbon1095[carbon1095$climate=="Ambient",]
drought=carbon1095[carbon1095$climate=="Drought",]

plot3=ggplot(ambient, aes(x=Enz_Alloc, y=taxaEffect))+
  geom_point(alpha=0.8, size=2, aes(color=avgEnz_Alloc2))+
  #geom_point(size=3, alpha=0.5, aes(color=avgEnz_Alloc))+
  stat_smooth(method="lm", se=F, color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.key.height = unit(0.4, 'cm'),
        axis.text.x = element_text(size = 8, angle=90, vjust=0.5),
        plot.title = element_text(hjust = 0.5))+
  xlab("Population Relative Enzyme Production")+
  ylab(expression(Delta*" substrate degradation (%)")) +
  guides(color="none", alpha="none", size="none")+
  scale_alpha_continuous(limits=c(0, 12000))+
  #annotate("text", x=5, y=3, label="R2 = 0.47")+
  facet_wrap(~avgEnz_Alloc2, nrow = 2)+
  NULL
plot3
#ggsave("~/Desktop/DEMENT_Figures/Figures_20231115/Fig3.png", height=4, width=6.8)

# Plot supplemental figure of all climate scenarios
plot3=ggplot(carbon1095, aes(x=Enz_Alloc, y=taxaEffect))+
  geom_point(alpha=0.8, size=2, aes(color=avgEnz_Alloc2))+
  #geom_point(size=3, alpha=0.5, aes(color=avgEnz_Alloc))+
  stat_smooth(method="lm", se=F, color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        #legend.key.height = unit9(0.4, 'cm'),
        axis.text.x = element_text(size = 8, angle=90, vjust=0.5),
        plot.title = element_text(hjust = 0.5))+
  xlab("Population Relative Enzyme Production")+
  ylab(expression(Delta*" substrate degradation (%)")) +
  guides(color="none", alpha="none", size="none")+
  scale_alpha_continuous(limits=c(0, 12000))+
  #annotate("text", x=5, y=3, label="R2 = 0.47")+
  facet_grid(climate~round(avgEnz_Alloc, digits=5))+
  NULL
plot3
#ggsave("~/Desktop/DEMENT_Figures/Figures_20231115/FigS5.png", height=5, width=9)

# Ambient linear regressions ####
summary(lm(taxaEffect~Enz_Alloc, data=ambient))
AIC(lm(taxaEffect~Enz_Alloc, data=ambient))

summary(lm(taxaEffect~Uptake_Alloc, data=ambient))
AIC(lm(taxaEffect~Uptake_Alloc, data=ambient))

summary(lm(taxaEffect~Osmo_Alloc, data=ambient))
AIC(lm(taxaEffect~Osmo_Alloc, data=ambient))

summary(lm(taxaEffect~Enz_Alloc*Uptake_Alloc, data=ambient))
AIC(lm(taxaEffect~Enz_Alloc*Uptake_Alloc, data=ambient))

summary(lm(taxaEffect~Enz_Alloc*Osmo_Alloc, data=ambient))
AIC(lm(taxaEffect~Enz_Alloc*Osmo_Alloc, data=ambient))

summary(lm(taxaEffect~Uptake_Alloc*Osmo_Alloc, data=ambient))
AIC(lm(taxaEffect~Uptake_Alloc*Osmo_Alloc, data=ambient))

summary(lm(taxaEffect~avgEnz_Alloc, data=ambient))
AIC(lm(taxaEffect~avgEnz_Alloc, data=ambient))

summary(lm(taxaEffect~Enz_Alloc*avgEnz_Alloc, data=ambient))
AIC(lm(taxaEffect~Enz_Alloc*avgEnz_Alloc, data=ambient))

# Drought linear regressions ####
drought=carbon1095[carbon1095$climate=="Drought",]

summary(lm(taxaEffect~Enz_Alloc, data=drought))
AIC(lm(taxaEffect~Enz_Alloc, data=drought))

summary(lm(taxaEffect~Uptake_Alloc, data=drought))
AIC(lm(taxaEffect~Uptake_Alloc, data=drought))

summary(lm(taxaEffect~Osmo_Alloc, data=drought))
AIC(lm(taxaEffect~Osmo_Alloc, data=drought))

summary(lm(taxaEffect~Enz_Alloc*Uptake_Alloc, data=drought))
AIC(lm(taxaEffect~Enz_Alloc*Uptake_Alloc, data=drought))

summary(lm(taxaEffect~Enz_Alloc*Osmo_Alloc, data=drought))
AIC(lm(taxaEffect~Enz_Alloc*Osmo_Alloc, data=drought))

summary(lm(taxaEffect~Uptake_Alloc*Osmo_Alloc, data=drought))
AIC(lm(taxaEffect~Uptake_Alloc*Osmo_Alloc, data=drought))

summary(lm(taxaEffect~avgEnz_Alloc, data=drought))
AIC(lm(taxaEffect~avgEnz_Alloc, data=drought))

summary(lm(taxaEffect~Enz_Alloc*avgEnz_Alloc, data=drought))
AIC(lm(taxaEffect~Enz_Alloc*avgEnz_Alloc, data=drought))

# Moist linear regressions ####
moist=carbon1095[carbon1095$climate=="Moist",]

summary(lm(taxaEffect~Enz_Alloc, data=moist))
AIC(lm(taxaEffect~Enz_Alloc, data=moist))

summary(lm(taxaEffect~Uptake_Alloc, data=moist))
AIC(lm(taxaEffect~Uptake_Alloc, data=moist))

summary(lm(taxaEffect~Osmo_Alloc, data=moist))
AIC(lm(taxaEffect~Osmo_Alloc, data=moist))

summary(lm(taxaEffect~Enz_Alloc*Uptake_Alloc, data=moist))
AIC(lm(taxaEffect~Enz_Alloc*Uptake_Alloc, data=moist))

summary(lm(taxaEffect~Enz_Alloc*Osmo_Alloc, data=moist))
AIC(lm(taxaEffect~Enz_Alloc*Osmo_Alloc, data=moist))

summary(lm(taxaEffect~Uptake_Alloc*Osmo_Alloc, data=moist))
AIC(lm(taxaEffect~Uptake_Alloc*Osmo_Alloc, data=moist))

summary(lm(taxaEffect~avgEnz_Alloc, data=moist))
AIC(lm(taxaEffect~avgEnz_Alloc, data=moist))

summary(lm(taxaEffect~Enz_Alloc*avgEnz_Alloc, data=moist))
AIC(lm(taxaEffect~Enz_Alloc*avgEnz_Alloc, data=moist))


#### Look at relationships between drought tolerance and taxa effect ####

# Plot Figure 4
plot4=ggplot(drought, aes(x=Osmo_Alloc, y=taxaEffect))+
  geom_point(alpha=0.8, aes(color=ratioEU*100))+
  geom_hline(yintercept = 0, color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.4, 'cm'))+ 
  ylab(expression(Delta*" substrate degradation (%)")) +
  xlab("Drought tolerance")+
  #facet_grid(~avgEnz_Alloc2)+
  scale_color_continuous(breaks = c(0.1, 0.4, 0.8), low="tan1", high="blue4", name=expression(over("Enzyme", "Monomer")))
plot4
#ggsave("~/Desktop/DEMENT_Figures/Figures_20231115/Fig4.png", height=2.8, width=3.2)

# Plot supplemental figure across all climate scenario
plot4=ggplot(carbon1095, aes(x=Osmo_Alloc, y=taxaEffect))+
  #geom_point(data=ambient, alpha=0.8, aes(color=ratioEU*100))+
  #geom_point(data=drought, alpha=0.8, aes(color=ratioEU*100))+
  geom_point(alpha=0.8, aes(color=ratioEU*100))+
  geom_hline(yintercept = 0, color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.4, 'cm'))+ 
  ylab(expression(Delta*" substrate degradation (%)")) +
  xlab("Drought tolerance")+
  facet_grid(~climate)+
  ylim(-20, 20)+
  scale_color_continuous(breaks = c(0.1, 0.4, 0.8), low="tan1", high="blue4", name=expression(over("Enzyme", "Monomer")))
#plot4
#ggsave("~/Desktop/DEMENT_Figures/AGU/Fig1d.png", height=3.5, width=7)




#### Look at both ecosystem- and commmunity scale impacts on same figure ####
# Get the presence impact from carbon1095
temp=carbon1095[,c("file.x","climate", "taxaEffect", "ratioEU", "Osmo_Alloc", "avgEnz_Alloc")]
colnames(temp)=c("population", "climate", "taxaEffect", "ratioEnzymeUptake", "Osmo_Alloc", "avgEnz_Alloc")

# Merge
temp=merge(temp, biomassChange, by=c("population"))
biomassChange=temp

# Plot Figure 5A
plot5A=ggplot(biomassChange, aes(x=taxaEffect, y=biomassImpact))+
  geom_point(alpha=0.5, aes(color=ratioEnzymeUptake*100))+
  #geom_point(alpha=0, aes(color=ratioEnzymeUptake*100))+
  geom_hline(yintercept = 0, color="black")+
  geom_vline(xintercept = 0, color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.direction="vertical",
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.4, 'cm'),
        #legend.position = c(0.08,0.77), 
        #legend.background = element_rect(fill = "transparent")
  )+
  xlab(expression(Delta*" substrate degradation (%)")) +
  ylab(expression(Delta*" neighbor biomass (mg C)"))+
  xlim(-20,20)+ylim(-2216,2216)+
  facet_grid(~climate.x)+
  scale_color_continuous(breaks = c(0.1, 0.4, 0.8), low="tan1", high="blue4", name=expression(over("Enzyme", "Monomer")))+
  NULL
plot5A

# Plot 5B conceptual figure
plot5B=ggplot(biomassChange, aes(x=c(-10,10), y=c(-10,10)))+
  geom_hline(yintercept = 0, color="grey50")+
  geom_vline(xintercept = 0, color="grey50")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0.5,4,0.5,3, "cm"))+
  xlab(expression(Delta*" ecosystem function")) +
  ylab(expression(Delta*" neighbor fitness"))+
  annotate("text", x=-7.5, y=-9, label="Exploitation", size = 10/.pt)+
  annotate("text", x=7.5, y=9, label="Facilitation", size = 10/.pt)+
  annotate("text", x=7.5, y=-9, label="Competition", size = 10/.pt)+
  xlim(-10,10)+ylim(-10,10)+
  NULL
plot5B

plot5=plot_grid(plot5A, plot5B, nrow=2, #rel_widths = c(1, 0.35), 
                #align="hv", axis="l", 
          labels = c("A", "B"))

#ggsave("~/Desktop/DEMENT_Figures/Figures_20231115/Fig5.png", plot = plot5, height=5, width=6.8)


