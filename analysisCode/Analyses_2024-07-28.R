library(ggplot2)
library(cowplot)
library(tidyverse)

#setwd("/Users/brittnibertolet/OneDrive - UC Irvine/GitHub/DEMENTinteractions/simulationDataFinal/")
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
temp$subDegraded=(temp$substrate730-temp$substrate1095)

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

#Percent change 
for(i in 1:length(strings)){
  carbon1095$taxaEffect[grepl(strings[i],carbon1095$file)]=
    (carbon1095all$subDegraded[grepl(strings[i] ,carbon1095all$file)]-
       carbon1095$subDegraded[grepl(strings[i],carbon1095$file)])
}

# Get rid of effect of their own exclusion
length(carbon1095$taxaEffect[carbon1095$taxaEffect==0])
carbon1095=carbon1095[carbon1095$taxaEffect!=0,]

# Get community and population trait data
traits=read.csv("traitOutput.csv", stringsAsFactors = F)
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

# Visualize trait data
a=ggplot(traits[grepl("moist", traits$file),], aes(x=Enz_Alloc))+geom_histogram(fill="grey80", bins=20)+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Enzyme allocation")+ylab("Count")
range(traits$Enz_Alloc)
b=ggplot(traits[grepl("moist", traits$file),], aes(x=Uptake_Alloc))+geom_histogram(fill="grey80", bins=20)+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Uptake allocation")+ylab("Count")
range(traits$Uptake_Alloc)
c=ggplot(traits[grepl("moist", traits$file),], aes(x=Osmo_Alloc))+geom_histogram(fill="grey80", bins=20)+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Osmolyte allocation")+ylab("Count")
range(traits$Osmo_Alloc)
d=ggplot(traits[grepl("moist", traits$file),], 
         aes(x=Enz_Alloc, y=Uptake_Alloc))+geom_point(alpha=0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Enzyme allocation")+ylab("Uptake allocation")
e=ggplot(traits[grepl("moist", traits$file),], 
         aes(x=Enz_Alloc, y=Osmo_Alloc))+geom_point(alpha=0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Enzyme allocation")+ylab("Osmolyte allocation")
f=ggplot(traits[grepl("moist", traits$file),], 
         aes(x=Uptake_Alloc, y=Osmo_Alloc))+geom_point(alpha=0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Uptake allocation")+ylab("Osmolyte allocation")

plot_grid(a, b, c, d, e, f, nrow=2, align="hv")
#ggsave("~/Desktop/DEMENT_Figures/Figures_20240715/FigS2-revision.png", height=6, width=8)

#Get average weighted traits of the community 
comAvg=data.frame()
seeds=unique(carbon1095$microbe)
for(i in 1:length(seeds)){
  #temp=traits[carbon1095$microbe==seeds[i],]
  temp=traits[grepl(seeds[i], traits$file),]
  
  temp=temp[grepl("grass", temp$file),]
  temp3=data.frame(microbe=seeds[i],
                   avgEnz_Gene=weighted.mean(temp$Enz_Gene, temp$avgAbund),
                   avgEnz_Alloc=weighted.mean(temp$Enz_Alloc, temp$avgAbund),
                   avgDrought_Tol=weighted.mean(temp$Drought_tolerance, temp$avgAbund),
                   avgOsmo_Alloc=weighted.mean(temp$Osmo_Alloc, temp$avgAbund),
                   avgUptake_Alloc=weighted.mean(temp$Uptake_Alloc, temp$avgAbund),
                   avgRation=weighted.mean(temp$ratioEU, temp$avgAbund)
                   
  )
  comAvg=rbind(comAvg, temp3)
}

#Traits
traits$taxaID=gsub("_all.pickle", "", as.character(traits$file))
traits$taxaID=paste(traits$taxaID, traits$X, sep="_")

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
carbon1095$avgEnz_Alloc2=paste(round(carbon1095$avgEnz_Alloc, digits=8))
length(unique(carbon1095$avgEnz_Alloc))
length(unique(carbon1095$microbe))

library(viridis)

# Plot change in substrate degraded when taxa is excluded
plot2A=ggplot(carbon1095, aes(x=climate, y=taxaEffect/10000))+
  geom_hline(yintercept = 0, color="black")+
  #geom_point(aes(color=avgEnz_Alloc2), position=position_dodge(width = .7), alpha=0.4, size=1)+
  geom_point(aes(color=avgEnz_Alloc, fill=avgEnz_Alloc2), position=position_dodge(width = .7), size=1, alpha=0.5)+
  geom_boxplot(outlier.shape = NA, alpha=0.7)+
  theme_bw()+theme(panel.grid = element_blank(),
                   strip.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.text.y = element_text(size = 8, angle=90, hjust = 0.5),
                   legend.title = element_text(size = 10), 
                   legend.text = element_text(size = 8),
                   legend.key.height = unit(0.4, 'cm'),
                   plot.margin = margin(0.5,0.3,0.1,0.3, "cm"))+
  ylab(expression(paste(Delta*" substrate degraded (mg C ", cm^-3,")")))+
  xlab("Climate Scenario")+
  labs(color="Community-average\nenzyme production")+
  guides(fill="none")+
  scale_color_viridis()
plot2A

# Plot total substrate degraded
plot2Ar=ggplot(carbon1095, aes(x=climate, y=subDegraded/10000))+
  #geom_point(aes(color=avgEnz_Alloc2), position=position_dodge(width = .7), alpha=0.4, size=1)+
  geom_point(aes(color=avgEnz_Alloc, fill=avgEnz_Alloc2), position=position_dodge(width = .7), size=1, alpha=0.5)+
  geom_boxplot(outlier.shape = NA, alpha=0.7)+
  #geom_violin()+
  scale_colour_viridis()+
  guides(fill="none")+
  #geom_boxplot(outlier.shape = NA, fill="NA")+
  theme_bw()+theme(panel.grid = element_blank(),
                   strip.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.text.y = element_text(size = 8, angle=90, hjust = 0.5),
                   legend.title = element_text(size = 10), 
                   legend.text = element_text(size = 8),
                   legend.key.height = unit(0.4, 'cm'),
                   plot.margin = margin(0.5,0.3,0.1,0.3, "cm"))+
  ylab(expression(paste("Total substrate degraded (mg C ", cm^-3,")")))+
  xlab("Climate Scenario")+
  labs(color="Community-average\nenzyme production")+
  NULL


mean(carbon1095$subDegraded[carbon1095$climate=="Moist"])/10000
sd(carbon1095$subDegraded[carbon1095$climate=="Moist"])/10000

mean(carbon1095$subDegraded[carbon1095$climate=="Ambient"])/10000
sd(carbon1095$subDegraded[carbon1095$climate=="Ambient"])/10000

mean(carbon1095$subDegraded[carbon1095$climate=="Drought"])/10000
sd(carbon1095$subDegraded[carbon1095$climate=="Drought"])/10000

# Look at taxon impacts under different scenarios 
subdegraded=carbon1095[, c("microbe", "grid", "Taxa","climate", "subDegraded")]
subdegraded=spread(subdegraded, key="climate",value="subDegraded")

plota=ggplot(subdegraded, aes(x=Ambient/10000, y=Drought/10000))+
  geom_point(alpha=0.3)+
  stat_smooth(method="lm")+
  theme_bw()+theme(panel.grid = element_blank(),
                   strip.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.text.y = element_text(size = 8, angle=90, hjust = 0.5),
                   legend.title = element_text(size = 10), 
                   legend.text = element_text(size = 8),
                   legend.key.height = unit(0.4, 'cm'),
                   plot.margin = margin(0.5,0.3,0.1,0.3, "cm"))+
  xlab(expression(paste("Ambient Condition - Substrate degraded (mg C ", cm^-3,")")))+
  ylab(expression(paste("Drought Condition - Substrate degraded (mg C ", cm^-3,")")))+
  annotate(geom = "text", x=240, y=240, label="R^2 == 0.85", parse=T)

summary(lm(subdegraded$Drought~subdegraded$Ambient))  
  

plotb=ggplot(subdegraded, aes(x=Moist/10000, y=Drought/10000))+
  geom_point(alpha=0.3)+
  stat_smooth(method="lm")+
  theme_bw()+theme(panel.grid = element_blank(),
                   strip.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.text.y = element_text(size = 8, angle=90, hjust = 0.5),
                   legend.title = element_text(size = 10), 
                   legend.text = element_text(size = 8),
                   legend.key.height = unit(0.4, 'cm'),
                   plot.margin = margin(0.5,0.3,0.1,0.3, "cm"))+
  xlab(expression(paste("Moist Condition - Substrate degraded (mg C ", cm^-3,")")))+
  ylab(expression(paste("Drought Condition - Substrate degraded (mg C ", cm^-3,")")))+
  annotate(geom = "text", x=283, y=240, label="R^2 == 0.40", parse=T)
plotb
summary(lm(subdegraded$Drought~subdegraded$Moist))  


plotc=ggplot(subdegraded, aes(x=Ambient/10000, y=Moist/10000))+
  geom_point(alpha=0.3)+
  stat_smooth(method="lm")+
  theme_bw()+theme(panel.grid = element_blank(),
                   strip.text = element_text(size = 10),
                   axis.title = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.text.y = element_text(size = 8, angle=90, hjust = 0.5),
                   legend.title = element_text(size = 10), 
                   legend.text = element_text(size = 8),
                   legend.key.height = unit(0.4, 'cm'),
                   plot.margin = margin(0.5,0.3,0.1,0.3, "cm"))+
  xlab(expression(paste("Ambient Condition - Substrate degraded (mg C ", cm^-3,")")))+
  ylab(expression(paste("Moist Condition - Substrate degraded (mg C ", cm^-3,")")))+
  annotate(geom = "text", x=240, y=287, label="R^2 == 0.60", parse=T)
plotc
summary(lm(subdegraded$Moist~subdegraded$Ambient))

plot_grid(plota, plotc, plotb, nrow=1, labels=c("a", "b", "c"))
#ggsave("~/Desktop/DEMENT_Figures/Figures_20240715/FigS4-revision.png", height=4, width=11)

# Are there differences across climates, treating taxa id as a random effect?
carbon1095$taxaID2=str_split_fixed(carbon1095$taxaID, '_', 2)[,2]
unique(carbon1095$taxaID2)

library(rstatix)
# Two-way mixed factors ANOVA
res.aov <- anova_test(
  data = carbon1095, dv = taxaEffect, wid = taxaID2,
  between = microbe, within = climate)
get_anova_table(res.aov)


res.aov <- anova_test(
  data = carbon1095, dv = subDegraded, wid = taxaID2,
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
  tempOut=data.frame(population=pops[i], taxa=names(temp), climate=cli, biomassImpact=(tempFull-temp))
  
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
temp=carbon1095[,c("microbe", "avgEnz_Alloc", "avgEnz_Alloc2")]
temp=temp[!duplicated(temp),]
biomassChange=left_join(biomassChange, temp, by="microbe")
biomassChange=biomassChange[!biomassChange$microbe%in%c("microbe432", "microbe410"),]

plot2B=ggplot(biomassChange, aes(x=climate, y=biomassImpact/10000))+
  geom_hline(yintercept = 0, color="grey", linetype="dashed")+
  #geom_point(aes(color=avgEnz_Alloc2), position=position_dodge(width = .7), alpha=0.4, size=1)+
  geom_point(aes(color=avgEnz_Alloc, fill=avgEnz_Alloc2), position=position_dodge(width = .7), size=1, alpha=0.5)+
  #geom_violin()+
  scale_colour_viridis()+
  guides(fill="none")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
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
  ylab(expression(paste(Delta*" associate biomass (mg C ", cm^-3,")")))+
  xlab("Climate Scenario")+
  guides(fill="none")
plot2B


top=plot_grid(plot2Ar+guides(color="none"), plot2A+guides(color="none"), align="hv", axis="lr", labels=c("A", "B" ), nrow=1)
bottom=plot_grid(plot2B+guides(color="none")+theme(plot.margin = margin(0.5,0.5,0.5,3, "cm")),get_legend(plot2A),  
                 nrow=1, labels=c("C", ""), rel_widths = c(1,0.3))

plot_grid(top, bottom, nrow=2)
ggsave("~/Desktop/DEMENT_Figures/Figures_20240715/Fig2.png", height=6.5, width=6.5)

# Calculating how many interactions are positive versus negative for each community
intDir=c()
microbes=unique(biomassChange$microbe)
for(i in 1:length(microbes)){
  temp=biomassChange[biomassChange$microbe==microbes[i],]
  perNeg=nrow(temp[temp$biomassImpact<0,])/nrow(temp)
  perPos=nrow(temp[temp$biomassImpact>0,])/nrow(temp)
  tempOut=data.frame(microbe=microbes[i], perNeg=perNeg, perPos=perPos)
  # Concatonate to output
  intDir=rbind(intDir, tempOut)
  
}
mean(intDir$perNeg)

sd(intDir$perNeg)


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
# Plot supplemental figure of all climate scenarios

# Get example microbial communities to plot
# Community with low community average enzyme production and large impact on drought
#"microbe1510" "microbe3365"
# Community with high community average enzyme production and low impact on drought
#"microbe440" "microbe26152"
extremes=c("microbe3365", "microbe26152")

extremes=carbon1095[carbon1095$microbe%in%extremes,]
extremes$avgEnz_Alloc3=paste0("Cenz = ",extremes$avgEnz_Alloc2)
labels=as_labeller(c("Cenz = 0.00198294" = "Low Cenz\n(Cenz = 0.0019)", 
                     "Cenz = 0.0032246" = "High Cenz\n(Cenz = 0.0032)",
                     "Moist"="Moist", "Ambient"="Ambient", "Drought"="Drought"))
plot3=ggplot(extremes, aes(x=Enz_Alloc, y=taxaEffect/10000))+
  geom_point(aes(color=climate), size=2)+
  scale_color_manual(values=c("#00BFC4", "#7CAE00", "#F8766D"))+
  #geom_point(size=3, alpha=0.5, aes(color=avgEnz_Alloc))+
  stat_smooth(method="lm", color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        #legend.key.height = unit9(0.4, 'cm'),
        axis.text.x = element_text(size = 8, angle=90, vjust=0.5),
        plot.title = element_text(hjust = 0.5))+
  xlab("Taxon-specific relative enzyme production")+
  ylab(expression(paste(Delta*" substrate degraded (mg C ", cm^-3,")")))+
  guides(color="none", alpha="none", size="none")+
  scale_alpha_continuous(limits=c(0, 12000))+
  #annotate("text", x=5, y=3, label="R2 = 0.47")+
  facet_grid(climate~avgEnz_Alloc3, labeller = labels)+guides(shape="none", linetype="none")+
  NULL
plot3
ggsave("~/Desktop/DEMENT_Figures/Figures_20240509/FigS5.png", height=5, width=9)


# Get slope of each relationship
reps=unique(carbon1095$microbe)
cli=unique(carbon1095$climate)
slopes=c()
for(i in 1:length(reps)){
  for(j in 1:length(cli)){
    temp=carbon1095[carbon1095$microbe==reps[i],]
    temp=temp[temp$climate==cli[j],]
    out=summary(lm(temp$taxaEffect/10000~temp$Enz_Alloc))
    slope=out$coefficients[2,1]
    slopes=rbind(slopes, data.frame(microbe=reps[i], climate=cli[j], slope=slope))
  }
}

temp=merge(carbon1095, slopes, by=c("microbe", "climate"))
temp=temp[!duplicated(temp$slope),]
plot3B=ggplot(temp, aes(x=avgEnz_Alloc, y=slope))+
  geom_point(aes(color=climate), size=2)+
  scale_color_manual(values=c("#00BFC4", "#7CAE00", "#F8766D"))+
  stat_smooth(method="lm", aes(color=climate))+
  #scale_linetype_manual(values = c("Moist" = "dotdash", "Ambient" = "solid", "Drought"="dashed")) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.key.height = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5))+
  xlab("Community-average relative enzyme production")+
  #ylab(expression(paste("Effect of taxon-specific enzyme production on ",Delta*" substrate degraded (mg C ", cm^-3,")")))+
  ylab("Slope")+
  labs(color="Climate")
plot3B
  
summary(lm(slope~avgEnz_Alloc*climate, data=temp))

bottom=plot_grid(plot3B+theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))+guides(color="none", shape="none", linetype="none"), 
                 get_legend(plot3B), nrow=1, rel_widths = c(1,0.3))
bottom
plot_grid(plot3, bottom, nrow=2, labels=c("A", "B"))

ggsave("~/Desktop/DEMENT_Figures/Figures_20240715/Fig3.png", height=6.5, width=6)

# Get the slopes to put labels on the figures
slopeEx=temp[temp$microbe%in%extremes$microbe,]

# Ambient linear regressions ####
ambient=carbon1095[carbon1095$climate=="Ambient",]

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
plot4=ggplot(drought, aes(x=Osmo_Alloc, y=taxaEffect/10000))+
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
  ylab(expression(paste(Delta*" substrate degraded (mg C ", cm^-3,")")))+
  xlab("Drought tolerance")+
  #facet_grid(~avgEnz_Alloc2)+
  scale_color_continuous(breaks = c(0.1, 0.4, 0.8), low="tan1", high="blue4", name=expression(over("Enzyme", "Monomer")))
plot4
ggsave("~/Desktop/DEMENT_Figures/Figures_20240715/Fig4.png", height=2.8, width=3.2)

# Plot supplemental figure across all climate scenario
plot4=ggplot(carbon1095, aes(x=Osmo_Alloc, y=taxaEffect/10000))+
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
  ylab(expression(paste(Delta*" substrate degraded (mg C ", cm^-3,")")))+
  xlab("Drought tolerance")+
  facet_grid(~climate)+
  #ylim(-20, 20)+
  scale_color_continuous(breaks = c(0.1, 0.4, 0.8), low="tan1", high="blue4", name=expression(over("Enzyme", "Monomer")))
plot4
ggsave("~/Desktop/DEMENT_Figures/Figures_20240715/FigS6.png", height=3.5, width=7)




#### Look at both ecosystem- and commmunity scale impacts on same figure ####
# Get the presence impact from carbon1095
temp=carbon1095[,c("file.x","climate", "taxaEffect", "ratioEU", "Osmo_Alloc", "avgEnz_Alloc")]
colnames(temp)=c("population", "climate", "taxaEffect", "ratioEnzymeUptake", "Osmo_Alloc", "avgEnz_Alloc")

# Merge
temp=merge(temp, biomassChange, by=c("population"))
biomassChange=temp

# Plot Figure 5A
plot5A=ggplot(biomassChange, aes(x=taxaEffect/10000, y=biomassImpact/10000))+
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
  xlab(expression(paste(Delta*" substrate degraded (mg C ", cm^-3,")")))+
  ylab(expression(paste(Delta*" biomass (mg C ", cm^-3,")")))+
  xlim(-70,70)+ylim(-0.31,0.31)+
  facet_grid(~climate.x)+
  scale_color_continuous(breaks = c(0.1, 0.4, 0.8), low="tan1", high="blue4", name=expression(over("Enzyme", "Monomer")))+
  NULL
plot5A

range(biomassChange$taxaEffect/10000)
range(biomassChange$biomassImpact/10000)

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
  ylab(expression(Delta*" associate fitness"))+
  #annotate("text", x=-7.5, y=-9, label="Exploitation", size = 10/.pt)+
  #annotate("text", x=7.5, y=9, label="Facilitation", size = 10/.pt)+
  #annotate("text", x=7.5, y=-9, label="Competition", size = 10/.pt)+
  xlim(-10,10)+ylim(-10,10)+
  NULL
plot5B

plot5=plot_grid(plot5A, plot5B, nrow=2, #rel_widths = c(1, 0.35), 
                #align="hv", axis="l", 
          labels = c("A", "B"))

ggsave("~/Desktop/DEMENT_Figures/Figures_20240715/Fig5.png", plot = plot5, height=5, width=6.8)



