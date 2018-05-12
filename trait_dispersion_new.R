library(maps)
library(ggplot2)
library(ape)
library(caper)
library(gridExtra)
library(dplyr)
library(vegan)
library(FD)

setwd("C:/Users/User/Dropbox/UBC/AVILES LAB/CommComp Project/Analysis/")

data<-read.csv("data.csv", stringsAsFactors=TRUE) 
species<-read.csv("species.csv", row.names = 1) 
#Arcsine-transform sociality index
species$Sociality.index <- asin(species$Sociality.index)
species$Sociality.index.scaled <- scale(species$Sociality.index)
comm<-read.csv("communitydata.csv", row.names = 1) 
env<-read.csv("sites2016.csv", row.names = 1) 

# Remove spcies found alone:
#species <- species[-c(which(row.names(species)=="vierae"), which(row.names(species)=="analyticus")),]

# Phylogenetic datasets
species.phylo<-read.csv("species_phylo2.csv", row.names = 1)
phylo <- read.tree("anelosimus_tree.tre") 
phylo<-drop.tip(phylo,setdiff(as.character(phylo$tip.label),species.phylo$Species)) #remove species not in our data

#### PGLS ANALYSES ####

#Phylogenetic tree of species with prey size data
#phylo$tip.label <- paste("A.", phylo$tip.label) #For label in fig.2
ts_tree <- chronopl(phylo, lambda=1)
anel.tree <- plot.phylo(ts_tree, label.offset = 0.05)



#Switch sociality to continuous index (scaled and arcsin-transformed)
species.phylo$Sociality.index <- asin(species.phylo$Sociality.index)
species.phylo$Sociality.index.scaled <- scale(species.phylo$Sociality.index)

'
tiff(filename = "C:/Users/User/Dropbox/UBC/AVILES LAB/CommComp Project/Graphs and figures/FigS1_phylo.tiff", 
     res = 300, width = 7, height = 6, units = "in")
plot.phylo(ts_tree, label.offset = 0.05)
dev.off()

pdf(file = "C:/Users/User/Dropbox/UBC/AVILES LAB/CommComp Project/Graphs and figures/phylo1.pdf", 
     width = 5, height = 4)
plot.phylo(ts_tree, label.offset = 0.05)
dev.off()

#write.csv(species.phylo[,c(1:4)], "traits_table2.csv")
'

#Combined phylo and trait data set
combined <- comparative.data(phy = ts_tree, data = species.phylo, names.col = Species)

### Sociality PGLS (all species) ###
(soc.pgls <- pgls(formula =  mean.prey.size ~ Sociality.index.scaled, data = combined, lambda = "ML"))

#Null PGLS
(summary(pgls(formula =  Sociality.index.scaled ~ 1, data = combined, lambda = "ML")))

anova(soc.pgls)
(summ.pgls.soc <- summary(soc.pgls))

#lambda of model
soc.pgls$param[2] 

#AIC of model
soc.pgls$aic


### Body size PGLS (only sociality 1-2) ###
species.phylo.sol<-species.phylo[species.phylo$Sociality<3,]
phylo.sol<-drop.tip(phylo,setdiff(as.character(phylo$tip.label),species.phylo.sol$Species)) #remove species not in our data
ts_tree.sol <- chronopl(phylo.sol, lambda=1)
combined.sol <- comparative.data(phy = ts_tree.sol, data = species.phylo.sol, names.col = Species)

(bs.pgls.sol <- pgls(formula =  mean.prey.size ~ Body.size.scaled, data = combined.sol, lambda = "ML"))

anova(bs.pgls.sol)
(summ.bs.pgls.sol <- summary(bs.pgls.sol))

#lambda of model
bs.pgls.sol$param[2] 


### Body size PGLS (only sociality 3-4) ###
species.phylo.soc<-species.phylo[species.phylo$Sociality>=3,]
phylo.soc<-drop.tip(phylo,setdiff(as.character(phylo$tip.label),species.phylo.soc$Species)) #remove species not in our data
ts_tree.soc <- chronopl(phylo.soc, lambda=1)
combined.soc <- comparative.data(phy = ts_tree.soc, data = species.phylo.soc, names.col = Species)

(bs.pgls.soc <- pgls(formula =  mean.prey.size ~ Body.size.scaled, data = combined.soc, lambda = "ML"))

anova(bs.pgls.soc)

#lambda of model
bs.pgls.soc$param[2] 

(summ.bs.pgls.soc <- summary(bs.pgls.soc))


### Body size PGLS (OVERALL, controlled for sociality level) ###
  #Intermediates removed because not enough data points (2)

species.phylo.all<-species.phylo[-which(species.phylo$Sociality==3),]
phylo.all<-drop.tip(phylo,setdiff(as.character(phylo$tip.label),species.phylo.all$Species)) #remove species not in our data
ts_tree.all <- chronopl(phylo.all, lambda=1)
combined.all <- comparative.data(phy = ts_tree.all, data = species.phylo.all, names.col = Species)

(bs.pgls <- pgls(formula =  mean.prey.size ~ sociality.cat+Body.size.scaled, data = combined.all, lambda = "ML"))

anova(bs.pgls)
(summ.bs.pgls <- summary(bs.pgls))

#lambda of model
bs.pgls$param[2] 
bs.pgls$aic


### Compare LM and PGLS ###

#Sociality
fit1 <- lm(mean.prey.size ~ Sociality.index.scaled, data = species)
anova(fit1)
summary(fit1)

  # Delta AIC
abs(soc.pgls$aic-AIC(fit1))

#Body size
fit2 <- lm(mean.prey.size ~ sociality.cat+Body.size.scaled, data = species[-which(species$Sociality==3),])
anova(fit2)
summary(fit2)

  #Delta AIC
abs(bs.pgls$aic-AIC(fit2))



#### TRAITS CONTRIBUTION TO PREY CAPTURE ####

#data is scaled:
#function scale() calculates the mean and standard deviation of the entire vector, 
#then "scale" each element by those values by subtracting the mean and dividing by the sd

#Replace sociality by continuous index
range(species$Sociality.index.scaled[!is.na(species$mean.prey.size)])
plot(species$mean.prey.size ~ species$Sociality.index.scaled)

#Data set for plotting with ggplot2 
  #Add Species column and remove species without prey size data
species.gg <- species[which(!is.na(species$mean.prey.size)),]
species.gg$species <- row.names(species.gg)
species.gg$species <- paste("A.", species.gg$species)

### SOCIALITY ###
fit1 <- lm(mean.prey.size ~ Sociality.index.scaled, data = species)
anova(fit1)
(summ.fit1<-summary(fit1))

plot1 <- ggplot(species.gg, aes(x=Sociality.index.scaled,y=mean.prey.size))+
  geom_point(aes(shape=factor(species)),size=3)+
  geom_abline(intercept = summ.pgls.soc$coefficients[1,1], 
              slope = summ.pgls.soc$coefficients[2,1], size = 1)+
  scale_shape_manual(name  ="Anelosimus species",
                     values=c(0,16,17,1,18,7,9,2,5))+
  geom_smooth(method = "lm", se = FALSE, fullrange=FALSE, 
              linetype = 2, size = 0.5, color = "grey40")+
  xlab("Sociality index (0.0-1.0,\narcsine-transformed and scaled)")+
  ylab("Mean prey size (mm)")+
  ylim(0,16)+
  geom_text(x = min(species$Sociality.index.scaled)+0.1, y = 15, label = "a", size = 12)+
  #geom_text(x = 1.2, y = 7, label = paste("slope = ", round(summ.pgls.soc$coefficients[2,1], 2), sep = "")) +
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=12))

plot1

### BODY SIZE ###

# Effect of BS controlling for sociality (interaction not significant)
fit.all <- lm(mean.prey.size ~ sociality.cat+Body.size.scaled, data = species[-which(species$Sociality==3),])
anova(fit.all)
(summ.fit.all <- summary(fit.all))

# By sociality category
fit.sol <- lm(mean.prey.size ~ Body.size.scaled, data = species[species$sociality.cat=="solitary",])
anova(fit.sol)
summary(fit.sol)

fit.soc <- lm(mean.prey.size ~ Body.size.scaled, data = species[species$Sociality==4,])
anova(fit.soc)
summary(fit.soc)

plot2 <- ggplot(species.gg, aes(Body.size.scaled,mean.prey.size))+ 
  geom_point(aes(shape=factor(species)), size=3)+
  stat_smooth(data = subset(species.gg, Sociality<=2), se = FALSE, fullrange=FALSE, 
              method="lm", linetype = 2, size = 0.5, color = "grey40")+
  stat_smooth(data = subset(species.gg, Sociality==4), se = FALSE, fullrange=FALSE, 
              method="lm", linetype = 2, size = 0.5, color = "grey40")+
  geom_abline(intercept = summ.bs.pgls$coefficients[1,1]+(summ.bs.pgls$coefficients[2,1]/2), 
              slope = summ.bs.pgls$coefficients[3,1], size = 1)+
  scale_shape_manual(name  ="Species",
                     values=c(0,16,17,1,18,7,9,2,5))+
  xlab("Body size (mm, scaled)\n ")+
  xlim(-1,1.5)+
  ylab("Mean prey size (mm)")+
  ylim(0,16)+
  geom_text(x = -0.9, y = 15, label = "b", size = 12)+
  #geom_text(x = 1, y = 8.5, label = paste("slope = ", round(summ.bs.pgls$coefficients[3,1], 2), sep = "")) +
  theme_classic()+
  theme(axis.title=element_text(size=12),
        axis.title.y=element_blank())+
  theme(legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(face = "italic"))

plot2


# Plot both graphs in the same figure
grid.arrange(plot1, plot2, ncol=2, widths = c(475,625))

'
tiff(filename = "C:/Users/User/Dropbox/UBC/AVILES LAB/CommComp Project/Graphs and figures/Fig2_ab_index.tiff", 
     res = 300, pointsize = 1, width = 7, height = 3.1, units = "in")
grid.arrange(plot1, plot2, ncol=2, widths = c(2.9,4.1))
dev.off()

pdf(file = "C:/Users/User/Dropbox/UBC/AVILES LAB/CommComp Project/Graphs and figures/Fig2_ab.pdf", 
     width = 7, height = 2.9)
grid.arrange(plot1, plot2, ncol=2, widths = c(2.9,4.1))
dev.off()
'


#### HABITAT FILTER ####

#Remove species not occuring at any site ("analyticus" and "vierae")?
#species<- species[-which(rownames(species) %in% setdiff(rownames(species), colnames(comm))),]

#Species list to draw from in null model
filt <- list()
filt[[1]]<-row.names(species)[which(between(species$Sociality, 1, 2))] 
filt[[2]]<-row.names(species)[which(between(species$Sociality, 1, 3))] 
filt[[3]]<-row.names(species)[which(between(species$Sociality, 1, 4))] 
filt[[4]]<-row.names(species)[which(between(species$Sociality, 4, 4))]

for(i in 1:length(levels(as.factor(env$Community.Type)))){
  names(filt)[[i]] <- levels(as.factor(env$Community.Type))[i]
}

#Filter of each site
site.filter <- env$Community.Type


#### FUNCTIONAL DISPERSION ####
# in gowdis() each variable (column) is first standardized 
# by subtracting the minimum value and dividing each entry 
# by the range of the corresponding variable; 
# consequently the rescaled variable has range [0,1], exactly.
# Using Podani's method for ordinal variables. 
# ord.method specifies the method for ordinal variable.
# "podani" refers to Eqs. 2a-b of Podani (1999),
# "metric" refers to his Eq. 3 (see 'details'); both options convert ordinal variables to ranks.
# "classic" simply treats ordinal variables as continuous variables.

source("null_model_FD.R")

#Sociality level as ordinal?
#species$Sociality <- ordered(species$Sociality)

### Trait dispersion (FD) ###
#Using PGLS coefficient estimates to calculate trait weights 

# RATIO when solitary/subsocial and social species combined were considered 
(ratio <- c(as.numeric(coef(soc.pgls)[2]/coef(bs.pgls)[3]),1))

null.FD <- null_model_FD(comm.matrix = comm, 
                   species.data = species[c(9,3)], # Dataframe of sociality and BS (needs to be in this order)
                   env.filter = site.filter,
                   filter.list = filt, 
                   trait.weights = ratio, 
                   n.iter = 1000) 

null.FD$mean.null
sd(null.FD$disp.null)
null.FD$mean.obs
sd(null.FD$disp.obs)
null.FD$pvalue

# RATIO when only the solitary-subsocial species
(ratio2 <- c(as.numeric(coef(soc.pgls)[2]/coef(bs.pgls.sol)[2]),1))

null.FD2 <- null_model_FD(comm.matrix = comm, 
                          species.data = species[c(1,3)], # Dataframe of sociality and BS (needs to be in this order)
                          env.filter = site.filter,
                          filter.list = filt, 
                          trait.weights = ratio2,
                          n.iter = 1000) 

null.FD2$mean.null
sd(null.FD2$disp.null)
null.FD2$mean.obs
sd(null.FD2$disp.obs)
null.FD2$pvalue

#Only sociality
null.FD.soc <- null_model_FD(comm.matrix = comm, 
                         species.data = species[c(9)], 
                         env.filter = site.filter,
                         filter.list = filt, 
                         trait.weights = c(1),
                         n.iter = 10) 

null.FD.soc$mean.null
sd(null.FD.soc$disp.null)
null.FD.soc$mean.obs
sd(null.FD.soc$disp.obs)
null.FD.soc$pvalue

#Only Body Size
null.FD.bs <- null_model_FD(comm.matrix = comm, 
                         species.data = species[c(3)], 
                         env.filter = site.filter,
                         filter.list = filt, 
                         trait.weights = c(1),
                         n.iter = 10) 

null.FD.bs$mean.null
sd(null.FD.bs$disp.null)
null.FD.bs$mean.obs
sd(null.FD.bs$disp.obs)
null.FD.bs$pvalue



#### COMMUNITY COMPOSITION OF EACH SITE ####
species<-read.csv("species.csv", row.names = 1) 
commcomp<-read.csv("commcompPFF_NEW1.csv", row.names = 1, stringsAsFactors=TRUE) 
commcomp$Which.species <- paste("A. ", commcomp$Which.species, sep = "")

commcomp.gg<-commcomp %>%
  mutate(Site.label=paste(commcomp$Site.abbrev, ",", commcomp$Country.abbrev, sep="")) %>%
  arrange(plot.order)
#OR
commcomp.gg<-commcomp %>%
  mutate(Site.label=paste(commcomp$Site.abbrev2, ", ", commcomp$Country.abbrev, sep="")) %>%
  arrange(plot.order)

commcomp.gg$Site.label <- factor(commcomp.gg$Site.label, 
                                 levels = unique(commcomp.gg$Site.label))

shapes = c(0,4,35,2,11,5,6,7,8,9,10,3,12,13,14,15,16,36,17,38,18,20,1)

# Sizing requirements:
#Recommended max height: 54 picas / 9" / 22.5 cm.
#Use one of the following widths:
#   1 column wide (20.5 picas / 3.42" / 8.7 cm)
#   1.5 columns wide (27 picas / 4.5" / 11.4 cm)
#   2 columns wide (42.125 picas / 7" / 17.8 cm)

tiff(filename = "C:/Users/User/Dropbox/UBC/AVILES LAB/CommComp Project/Graphs and figures/Fig3_communities.tiff", 
     res = 300, width = 7, height = 6, units = "in")

g<-ggplot(commcomp.gg, aes(Sociality.index,Body.size))
g+geom_point(aes(shape=factor(Which.species)), size=2)+
  guides(shape = guide_legend(ncol = 1)) +
  scale_shape_manual(name  ="Species", values=shapes) +
  facet_wrap(~Site.label, ncol=8)+
  scale_x_continuous(lim = c(0,1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  xlab("Sociality index (0.0-1.0)")+
  ylab("Body size (mm)")+
  theme_bw()+ 
  theme(strip.text = element_text(size=6),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.spacing.x=unit(0.1, "lines"),
        panel.spacing.y=unit(0, "lines"),
        legend.text = element_text(face="italic", size=8),
        legend.key.height=unit(0.75,"line"),
        axis.title=element_text(size=10),
        axis.text=element_text(size=8))

dev.off()


# other dimensions
tiff(filename = "C:/Users/User/Dropbox/UBC/AVILES LAB/CommComp Project/Graphs and figures/Fig3_communities_small.tiff", 
     res = 300, width = 7, height = 4, units = "in")

g<-ggplot(commcomp.gg, aes(Sociality.index,Body.size))
g+geom_point(aes(shape=factor(Which.species)), size=1.5)+
  guides(shape = guide_legend(ncol = 1)) +
  scale_shape_manual(name="Species", values=shapes) +
  facet_wrap(~Site.label, ncol=10)+
  scale_x_continuous(lim = c(0,1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  xlab("Sociality index (0.0-1.0)")+
  ylab("Body size (mm)")+
  theme_bw()+ 
  theme(strip.text = element_text(size=6),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.spacing.x=unit(0.1, "lines"),
        panel.spacing.y=unit(0, "lines"),
        legend.text = element_text(face="italic", size=6),
        legend.title = element_text(size=10),
        legend.key.height=unit(0.5,"line"),
        axis.title=element_text(size=10),
        axis.text=element_text(size=6))

dev.off()






