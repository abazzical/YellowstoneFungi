library(vegan)
library("FactoMineR")
library("factoextra")
library(viridis)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(dplyr)
library(ggord)
library(compositions)
library(wesanderson)


# RDA ---------------------------------------------------------------------


# read env table

soils <- read.csv("Yellowstone_Soils_temp.csv")
row.names(soils) <- soils$SampleID
transects = soils$Transect
cores = soils$Core
soilpca <- prcomp(soils[ ,3:16], scale. = TRUE) #SE
# PC1 ph 
vecs <- eigen(cov(soilpca))$vectors
## NOTE: Z <- X %*% V # principal components are the matrix * eigenvectors on each column
PCs <- soilpca %*% vecs

plotcoords=data.frame(soilpca$x)
envloading1=data.frame(soilpca$rotation)

envloading <- envloading1*10

# plot PCA for paper color code based on the core number in transect
pcadata <- soils[c(17), ]
soils$Core <- factor(soils$Core,levels = c("5", "4", "3", "2", "1"))
# prettier plot
PCAsoil<-ggplot()+
  geom_segment(aes(x=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0) ,y=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0), xend= c(envloading[1,1],envloading[2,1],envloading[3,1],envloading[4,1],envloading[5,1],envloading[6,1],envloading[7,1],envloading[8,1],envloading[9,1],envloading[10,1],envloading[11,1],envloading[12,1],envloading[13,1],envloading[14,1]), yend=c(envloading[1,2],envloading[2,2],envloading[3,2],envloading[4,2],envloading[5,2],envloading[6,2],envloading[7,2],envloading[8,2],envloading[9,2],envloading[10,2],envloading[11,2],envloading[12,2],envloading[13,2],envloading[14,2])), color="#999999", arrow=arrow(angle = 20, length = unit(0.1,"cm"), ends = "last", type = "open"), size = 1)+
  scale_color_viridis(discrete = soils$Core)+
  #scale_color_manual(values = wes_palette("BottleRocket2", 100, type="continuous"))+
  geom_point(aes(x=plotcoords$PC1, y=plotcoords$PC2, color = soils$Core), size = 4)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("PC 1 -- 34% of variance")+
  ylab("PC 2 -- 24% of variance")+
  labs(colour = "Core in Transect")+
  # geom_text(aes(x=plotcoords$PC1, y=plotcoords$PC2, label = pcadata$Sample, color = pcadata$Transect), nudge_x =0.25)+
  
  annotate("text", 
           x=c(envloading[1,1]-0.07,
               envloading[2,1],
               envloading[3,1]+0.1,
               envloading[4,1],
               envloading[5,1],
               envloading[6,1]+0.1,
               envloading[7,1]+0.25,
               envloading[8,1],
               envloading[9,1]+0.25,
               envloading[10,1],
               envloading[11,1]+0.2,
               envloading[12,1]+0.25,
               envloading[13,1],
               envloading[14,1]-0.2), 
           y=c(envloading[1,2]+0.1,
               envloading[2,2]+0.2,
               envloading[3,2]+0.2,
               envloading[4,2]+0.2,
               envloading[5,2]+0.2,
               envloading[6,2]+0.35,
               envloading[7,2],
               envloading[8,2]+0.2,
               envloading[9,2],
               envloading[10,2]-0.2,
               envloading[11,2],
               envloading[12,2],
               envloading[13,2]+0.2,
               envloading[14,2]+0.2), 
           label =  c("Moisture Content","pH", "PO4","Total Nitrogen","Total Carbon","Cd", "Co", "Cu", "Fe", "Mn", "Ni", "Pb", "Zn", "Temperature"), color = "#999999")
ggsave("PCoAsoilcore_cont.pdf", height=4, width=6.5, device="pdf") # save a PDF 3 inches by 4 inches


PC1cont = fviz_contrib(soilpca, choice = "var", axes = 1)
ggsave("PC1cont.pdf", height=4, width=6.5, device="pdf")
PC2cont = fviz_contrib(soilpca, choice = "var", axes = 2)
ggsave("PC2cont.pdf", height=4, width=6.5, device="pdf")

# Make a figure with both the PCA and the variable contributions for the first 2 PCs using cowplot
contr = plot_grid(PC1cont, PC2cont, labels = c('B',"C"))
Fig1<-plot_grid(PCAsoil, contr, labels = c('A'), nrow = 2, label_size = 16)
ggsave("Fig1.pdf", height=8, width=7, device="pdf") #


# A red reference dashed line is also shown on the barplot. 
# This reference line corresponds to the expected value if the contribution where uniform


# Prepare PCA ordination to use in RDA
# use the first two PCs of the environnmental variables 
plotcoords=data.frame(soilpca$x)
envi_rabbit = plotcoords[,1:2]
envi_df<-as.data.frame(envi_rabbit)
# envi_df<-as.matrix(envi_df[-c(2,17)])
# row.names(envi_df) = envi_df[ , 1]
# envi_df = envi_df[ , -1]

# set transects and blocks to test in anova and color in the plot
envi_df$Transect = soils$Transect
envi_df$Core = soils$Core




# transpose OTU table and env table 
otu_rabbit<-read.table("outs-table.txt", sep = "\t", header = TRUE)
rabbit_df<-as.data.frame(otu_rabbit)
otuID <- rabbit_df$OTU.ID
rabbit_transpose <- as.data.frame(t(as.matrix(rabbit_df[,-1])))
colnames(rabbit_transpose) <- otuID
rabbit_transpose<-rabbit_transpose[-c(70), ] # remove taxonomy at the bottom row BUT maybe should keep it for the plotting????

# remove not overlapping samples to have two matrices (OTU table and environmental variable/first 2 PCs) with the same number of rows
rabbit_trim<-rabbit_transpose[-c(1,4,14,66,67), ]
# env_trim<-as.data.frame(env_df[-c(24,25,36,37,38,39,43,44,48,49,50,51,54,55,58,59,61,64), ])
envi_trim<-as.data.frame(envi_df[-c(50,54), ])
# check order of rows is the same
all.equal(rownames(rabbit_trim),rownames(envi_df))
all.equal(rownames(rabbit_trim),rownames(envi_trim))

#make all the values in the dataframes numeric check with sapply(rabbit, class), sapply(env_rabbit, class)
# env_rabbit <- mutate_all(env_trim, function(x) as.numeric(as.character(x)))
envi_rabbit <- mutate_all(envi_trim, function(x) as.numeric(as.character(x)))
rabbit <- mutate_all(rabbit_trim, function(x) as.numeric(as.character(x)))

rownames(rabbit) = rownames(rabbit_trim)
rownames(envi_rabbit) = rownames(envi_trim)


rabbit_trim
envi_trim

#### #### #### #### #### #### #### 
#### before the RDA I would use a centered log ratio (CLR) transform on otu data. 
# Usually it makes the data more normal/less skewed.
#### #### #### #### #### #### #### 
# get compositions package (can also do the transformation manually)
# The CLR transform basically does this:
# 1) if there are zeros in the dataset, use a heuristic to remove them so the log transform is possible. 
        # I add a small value (1). I wonder if it distorts the data a bit, 
        # but the stats journal articles do this adding 1 thing.
# 2) get the relative abundance ("closure") of non-zero OTU table
# 3) log2 transform the relative abundance table ("closed" data)

rabbit_transform <- as.data.frame(clr((rabbit+1)/rowSums(rabbit+1)))

# the CLR transformation makes the sequencing data work in euclidean space/statistics.

#compute RDA
rda_tree = rda(rabbit_transform ~ PC1 + PC2, data=envi_rabbit) # with the CLR transform 
# rda_tree = rda(rabbit ~ . , data=envi_trim)


# initial plot RDA
plot(rda_tree, type='n', scaling=1)
orditorp(rda_tree, display = 'sp', cex=1, scaling=0.5, col='navy')
orditorp(rda_tree, display = 'sites', cex=1, scaling=0.5, col='deeppink')
#lc seems better
text(rda_tree, display='cn', cex=1, col='black')
text(rda_tree, display='wa', cex=1, col='deeppink')
plot(rda_tree)

#### 
## yeah, this looks fine - what do you think is the best way to show what soil variables load on each PC?

#R2 values and adjusted R2 
RsquareAdj(rda_tree)

# transect metadata to test the blocks of transects
# rownames of sites S

#RDA testing anova 1000 permutations
anova(rda_tree, permutations=1000, strata = envi_trim$Transect) # strata = by transect
anova(rda_tree, by='margin', permutations=1000, strata = envi_trim$Transect) ##### ####  hey! its significant!

RsquareAdj(rda_tree) #####  its like 5% - that's not unusual and its good.
summary(eigenvals(rda_tree, model = "constrained"))
signif.full <- anova.cca(rda_tree, parallel=getOption("mc.cores")) # default is permutation=999
signif.full # #### yay! it's significant!

screeplot(rda_tree)

# contribution of species
spcont = sort(round(100*scores(rda_tree, display = "sp", scaling = 0)[,1]^2, 3), decreasing = TRUE)
unsortsp = round(100*scores(rda_tree, display = "sp", scaling = 0)[,1]^2, 3)
write.table(unsortsp, "spunsortcont.txt", sep = "\t")
write.table(x = spcont,"spcont.txt", sep = "\t")

# plot nice RDA

# axes 1 & 2

# scores=data.frame(scores(rda_tree)$sites)
scores = data.frame(rda_tree$CCA$wa)
sites = row.names(rda_tree)
env_arrows <- as.data.frame(rda_tree$CCA$biplot) # these are the principal component arrow coordinates (environmental PCs)
env_arrows$RDA0 <- c(0,0) # need to add origins to make arrows in the plot

envi_trim$Core <- factor(envi_trim$Core,levels = c("5", "4", "3", "2", "1"))
rdafig = ggplot()+
  geom_point(aes(x = scores[,1], y = scores[,2], color = envi_trim$Core), size = 3)+
  scale_color_viridis(discrete = TRUE)+
  xlab("RDA1--18.1 % variation explained")+ # usually have this info on the plot
  ylab("RDA2--7.6 % variation explained")+
  # labs(colour = "env_arrows")+
  theme_bw() + 
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  labs(colour = "Core in Transect")+
  geom_segment(aes(x = env_arrows$RDA0, y = env_arrows$RDA0, xend = env_arrows$RDA1, yend = env_arrows$RDA2), 
               color="darkgray", 
               arrow=arrow(angle = 20, length = unit(0.15,"cm"), 
                           ends = "last", type = "open"), size = 0.6) +
  annotate("text", 
         x=c(env_arrows[1,1]-0.05,
             env_arrows[2,1]+0.1),
         y=c(env_arrows[1,2]+0.04,
             env_arrows[2,2]+0.05),
         label =  c("PC1","PC2"), color = "#999999")
ggsave("Fig2_RDA.pdf", height=6, width=7, device="pdf")
#annotate("text", x=c(biplot[1,1],biplot[2,1],biplot[3,1],biplot[4,1],biplot[5,1],biplot[6,1],biplot[7,1]), y=c(biplot[1,2],biplot[2,2],biplot[3,2],biplot[4,2],biplot[5,2],biplot[6,2],biplot[7,2]), label =  c("KS","SD","EP","RP","Poor","Sensitive","Tolerant"), color = "darkgray", size = 4)+
# annotate above was an alternative way to label ends of arrows. its more work but you can manually adjust each label in case you don't have photoshop;).

# theme(axis.text = element_text(size = 10))




# EM RDA

otu_EM<-read.table("EM-table.txt", sep = "\t", header = TRUE)

EM_df<-as.data.frame(otu_EM)
otuID <- EM_df$OTUID
EM_transpose <- as.data.frame(t(as.matrix(EM_df[,-1])))
colnames(EM_transpose) <- otuID
#EM_transpose<-EM_transpose[-c(52), ] # remove taxonomy at the bottom row
EM <- mutate_all(EM_transpose, function(x) as.numeric(as.character(x)))
rownames(EM) = rownames(EM_transpose)

# remove not overlapping samples to have two matrices with the same number of rows
EM_trim<-EM_transpose[-c(1,4,14), ]
EMenvi_trim<-as.data.frame(envi_df[-c(24,25,36,37,38,39,43,44,48,49,50,51,54,55,58,59,61,64), ])
EMenvi_trim<-as.data.frame(EMenvi_trim[-c(36,40,41,42,43,44,45,46,47,48), ]) #remove even more empty sites
# check order of rows is the same
all.equal(rownames(EM_trim),rownames(EMenvi_trim))
#make all the values in the dataframes numeric check with sapply(rabbit, class), sapply(env_rabbit, class)
# env_rabbit <- mutate_all(env_trim, function(x) as.numeric(as.character(x)))
EMenvi_rabbit <- mutate_all(EMenvi_trim, function(x) as.numeric(as.character(x)))
EMrabbit <- mutate_all(EM_trim, function(x) as.numeric(as.character(x)))

rownames(EMrabbit) = rownames(EM_trim)
rownames(EMenvi_rabbit) = rownames(EMenvi_trim)


# the CLR transformation makes the sequencing data work in euclidean space/statistics.
EMrabbit_transform <- as.data.frame(clr((EMrabbit+1)/rowSums(EMrabbit+1)))

#compute RDA
EMrda_tree = rda(EMrabbit_transform ~ PC1 + PC2, data=EMenvi_rabbit) # with the CLR transform 
# rda_tree = rda(rabbit ~ . , data=envi_trim)
#RDA testing anova 1000 permutations
anova(EMrda_tree, permutations=1000, strata = EMenvi_trim$Transect) # strata = by transect
anova(EMrda_tree, by='margin', permutations=1000, strata = EMenvi_trim$Transect) ##### ####  hey! its significant!

RsquareAdj(EMrda_tree) #####  its like 5% - that's not unusual and its good.
summary(eigenvals(EMrda_tree, model = "constrained"))
signif.full <- anova.cca(EMrda_tree, parallel=getOption("mc.cores")) # default is permutation=999
signif.full # #### yay! it's significant!

screeplot(EMrda_tree)

EMscores = data.frame(EMrda_tree$CCA$wa)
EMsites = row.names(EMrda_tree)
EMenv_arrows <- as.data.frame(EMrda_tree$CCA$biplot) # these are the principal component arrow coordinates (environmental PCs)
EMenv_arrows$RDA0 <- c(0,0) # need to add origins to make arrows in the plot

EMenvi_trim$Core <- factor(EMenvi_trim$Core,levels = c("5", "4", "3", "2", "1"))
EMrdafig = ggplot()+
  geom_point(aes(x = EMscores[,1], y = EMscores[,2], color = EMenvi_trim$Core), size = 3)+
  scale_color_viridis(discrete = TRUE)+
  xlab("RDA1--3.4 % variation explained")+ # usually have this info on the plot
  ylab("RDA2--1.7 % variation explained")+
  # labs(colour = "env_arrows")+
  theme_bw() + 
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  labs(colour = "Core in Transect")+
  geom_segment(aes(x = EMenv_arrows$RDA0, y = EMenv_arrows$RDA0, xend = EMenv_arrows$RDA1, yend = EMenv_arrows$RDA2), 
               color="darkgray", 
               arrow=arrow(angle = 20, length = unit(0.15,"cm"), 
                           ends = "last", type = "open"), size = 0.6) +
  annotate("text", 
           x=c(EMenv_arrows[1,1]-0.05,
               EMenv_arrows[2,1]+0.1),
           y=c(EMenv_arrows[1,2]+0.04,
               EMenv_arrows[2,2]+0.05),
           label =  c("PC1","PC2"), color = "#999999")
ggsave("Fig3_EMRDA.pdf", height=6, width=7, device="pdf")

# SAP RDA
otu_SAP<-read.table("SAP-table.txt", sep = "\t", header = TRUE)

SAP_df<-as.data.frame(otu_SAP)
otuID <- SAP_df$OTUID
SAP_transpose <- as.data.frame(t(as.matrix(SAP_df[,-1])))
colnames(SAP_transpose) <- otuID
#EM_transpose<-EM_transpose[-c(52), ] # remove taxonomy at the bottom row
SAP <- mutate_all(SAP_transpose, function(x) as.numeric(as.character(x)))
rownames(SAP) = rownames(SAP_transpose)

# remove not overlapping samples to have two matrices with the same number of rows
SAP_trim<-SAP_transpose[-c(1,4,14), ]
SAPenvi_trim<-as.data.frame(envi_df[-c(24,25,36,37,38,39,43,44,48,49,50,51,54,55,58,59,61,64), ])
SAPenvi_trim<-as.data.frame(SAPenvi_trim[-c(36,40,41,42,43,44,45,46,47,48), ]) #remove even more empty sites
# check order of rows is the same
all.equal(rownames(SAP_trim),rownames(SAPenvi_trim))
#make all the values in the dataframes numeric check with sapply(rabbit, class), sapply(env_rabbit, class)
# env_rabbit <- mutate_all(env_trim, function(x) as.numeric(as.character(x)))
SAPenvi_rabbit <- mutate_all(SAPenvi_trim, function(x) as.numeric(as.character(x)))
SAPrabbit <- mutate_all(SAP_trim, function(x) as.numeric(as.character(x)))

rownames(SAPrabbit) = rownames(SAP_trim)
rownames(SAPenvi_rabbit) = rownames(SAPenvi_trim)


# the CLR transformation makes the sequencing data work in euclidean space/statistics.
SAPrabbit_transform <- as.data.frame(clr((SAPrabbit+1)/rowSums(SAPrabbit+1)))

#compute RDA
SAPrda_tree = rda(SAPrabbit_transform ~ PC1 + PC2, data=SAPenvi_rabbit) # with the CLR transform 
# rda_tree = rda(rabbit ~ . , data=envi_trim)
#RDA testing anova 1000 permutations
anova(SAPrda_tree, permutations=1000, strata = SAPenvi_trim$Transect) # strata = by transect
anova(SAPrda_tree, by='margin', permutations=1000, strata = SAPenvi_trim$Transect) ##### ####  hey! its significant!

RsquareAdj(SAPrda_tree) #####  its like 5% - that's not unusual and its good.
summary(eigenvals(SAPrda_tree, model = "constrained"))
signif.full <- anova.cca(SAPrda_tree, parallel=getOption("mc.cores")) # default is permutation=999
signif.full # #### yay! it's significant!

SAPscores = data.frame(SAPrda_tree$CCA$wa)
SAPsites = row.names(SAPrda_tree)
SAPenv_arrows <- as.data.frame(SAPrda_tree$CCA$biplot) # these are the principal component arrow coordinates (environmental PCs)
SAPenv_arrows$RDA0 <- c(0,0) # need to add origins to make arrows in the plot

SAPenvi_trim$Core <- factor(SAPenvi_trim$Core,levels = c("5", "4", "3", "2", "1"))
SAPrdafig = ggplot()+
  geom_point(aes(x = SAPscores[,1], y = SAPscores[,2], color = SAPenvi_trim$Core), size = 3)+
  scale_color_viridis(discrete = TRUE)+
  xlab("RDA1--6.2 % variation explained")+ # usually have this info on the plot
  ylab("RDA2--2.5 % variation explained")+
  # labs(colour = "env_arrows")+
  theme_bw() + 
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  labs(colour = "Core in Transect")+
  geom_segment(aes(x = SAPenv_arrows$RDA0, y = SAPenv_arrows$RDA0, xend = SAPenv_arrows$RDA1, yend = SAPenv_arrows$RDA2), 
               color="darkgray", 
               arrow=arrow(angle = 20, length = unit(0.15,"cm"), 
                           ends = "last", type = "open"), size = 0.6) +
  annotate("text", 
           x=c(SAPenv_arrows[1,1]-0.05,
               SAPenv_arrows[2,1]+0.1),
           y=c(SAPenv_arrows[1,2]+0.04,
               SAPenv_arrows[2,2]+0.05),
           label =  c("PC1","PC2"), color = "#999999")
ggsave("Fig3_SAPRDA.pdf", height=6, width=7, device="pdf")


# plot both EM and SAP RDAs
Fig3<-plot_grid(SAPrdafig+ theme(legend.position="none"), EMrdafig, labels = c('A', 'B'), label_size = 16, rel_widths = c(.75, 1))
ggsave("Fig3_RDAs.pdf", height=5, width=10, device="pdf") #





# load rarefied OTU tables to report OTu numbers and estimate Shannon Index of richness

# transpose OTU table and env table 
otu_rabbit_rare<-read.table("otu-table.txt", sep = "\t", header = TRUE)
rabbit_df_rare<-as.data.frame(otu_rabbit_rare)
otuID <- rabbit_df_rare$OTUID
rabbit_transpose_rare <- as.data.frame(t(as.matrix(rabbit_df_rare[,-1])))
colnames(rabbit_transpose_rare) <- otuID
rabbit_transpose_rare<-rabbit_transpose_rare[-c(52), ] # remove taxonomy at the bottom row BUT maybe should keep it for the plotting????

# remove not overlapping samples to have two matrices (OTU table and environmental variable/first 2 PCs) with the same number of rows
rabbit_trim_rare<-rabbit_transpose_rare[-c(1,4,14), ]
# env_trim<-as.data.frame(env_df[-c(24,25,36,37,38,39,43,44,48,49,50,51,54,55,58,59,61,64), ])
envi_trim<-as.data.frame(envi_df[-c(24,25,36,37,38,39,43,44,48,49,50,51,54,55,58,59,61,64), ])
# check order of rows is the same
all.equal(rownames(rabbit_trim_rare),rownames(envi_trim))
all.equal(rownames(rabbit_trim_rare),rownames(envi_trim))

#make all the values in the dataframes numeric check with sapply(rabbit, class), sapply(env_rabbit, class)
# env_rabbit <- mutate_all(env_trim, function(x) as.numeric(as.character(x)))
envi_rabbit <- mutate_all(envi_trim, function(x) as.numeric(as.character(x)))
rabbit_rare <- mutate_all(rabbit_trim_rare, function(x) as.numeric(as.character(x)))



# richness Shannon Index and OTu numbers
allshannon = as.data.frame(diversity(rabbit_rare, index="shannon"))
allshannon$Transect = envi_trim$Transect
allshannon$Core = envi_trim$Core
allshannon$PC1 = envi_trim$PC1
allshannon$PC2 = envi_trim$PC2
allshannon$allno = rowSums(rabbit_rare != 0)
EMshannon = as.data.frame(diversity(EMrabbit, index="shannon"))
EMshannon$EMno = rowSums(EMrabbit != 0)
SAPshannon = as.data.frame(diversity(SAPrabbit, index="shannon"))
SAPshannon$SAPno = rowSums(SAPrabbit != 0)
redshannon = as.data.frame(allshannon[-c(36,40,41,42,43,44,45,46,47,48), ]) #reduce to the same sites as the other groups
redshannon$Core = SAPenvi_trim$Core
redshannon$Transect = SAPenvi_trim$Transect
shannon = cbind(redshannon, SAPshannon, EMshannon)
colnames(shannon)
names(shannon)[names(shannon) == "diversity(rabbit_rare, index = \"shannon\")"] <- "allOTUs"
names(shannon)[names(shannon) == "diversity(EMrabbit, index = \"shannon\")"] <- "EMOTUs"
names(shannon)[names(shannon) == "diversity(SAPrabbit, index = \"shannon\")"] <- "SAPOTUs"
names(shannon)[names(shannon) == "diversity(rabbit_rare, index = \"shannon\")"] <- "shannon"

shannon$Core <- factor(shannon$Core,levels = c("1", "2", "3", "4", "5"))

coreEM <- ggplot(shannon, aes(factor(Core), EMOTUs))
test<- coreEM + 
  geom_violin() + 
  geom_jitter(height = 0, width = 0.1) +
  ylim(0,4) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="grey") + 
  theme_bw() +
  geom_text(aes(x=3,y=4, label = "mean = 0.71")) +
  geom_text(aes(x=3,y=3.75, label = "median = 0.77")) +
  labs(
    x = "",
    y = "Ectomycorrhizal"
  ) 

noEM = ggplot(shannon, aes(factor(Core), EMno))
numEM = noEM +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1) +
  ylim(0,15) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="grey") + 
  theme_bw() +
  geom_text(aes(x=3,y=15, label = "mean = 4.5")) +
  geom_text(aes(x=3,y=14, label = "median = 4")) +
  labs(
    x = "",
    y = "Ectomycorrhizal"
  ) 


coreSAP <- ggplot(shannon, aes(factor(Core), SAPOTUs))
sapro<- coreSAP + 
  geom_violin() + 
  geom_jitter(height = 0, width = 0.1) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="grey") + 
  ylim(0,4) +
  theme_bw() +
  geom_text(aes(x=3,y=4, label = "mean = 1.27")) +
  geom_text(aes(x=3,y=3.75, label = "median = 1.36")) +
  labs(
    x = "",
    y = "Saprotrophic"
  ) 
noSAP = ggplot(shannon, aes(factor(Core), SAPno))
numSAP = noSAP +
  geom_violin() +
  ylim(0,25) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="grey") + 
  geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  geom_text(aes(x=3,y=25, label = "mean = 9.32")) +
  geom_text(aes(x=3,y=23, label = "median = 8.5")) +
  labs(
    x = "Core Number of Each Transect",
    y = "Saprotrophic"
  ) 

coreOTUs <- ggplot(shannon, aes(factor(Core), allOTUs))
OTUs<- coreOTUs + 
  geom_violin() + 
  geom_jitter(height = 0, width = 0.1) +
  ylim(0,4) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="grey") + 
  theme_bw() +
  geom_text(aes(x=3,y=4, label = "mean = 1.93")) +
  geom_text(aes(x=3,y=3.75, label = "median = 1.84")) +
  labs(
    x = "",
    y = "Shannon Index per Core"
  ) 
noall = ggplot(shannon, aes(factor(Core), allno))
numALL = noall +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1) +
  ylim(0,80) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="grey") + 
  theme_bw() +
  geom_text(aes(x=3,y=80, label = "mean = 39.66")) +
  geom_text(aes(x=3,y=75, label = "median = 37")) +
  labs(
    x = "",
    y = "OTUs per Core"
  ) 

plot_grid(OTUs, sapro, test, numALL, numSAP, numEM,  ncol = 3, labels = c("A", "B", "C", "D", "E", "F"))

# save figure 4
ggsave("Fig4_diversity.pdf", height=5.5, width=10, device="pdf")

mean(shannon$allOTUs)
mean(shannon$allno)
mean(shannon$EMOTUs)
mean(shannon$EMno)
mean(shannon$SAPOTUs)
mean(shannon$SAPno)

median(shannon$allOTUs)
median(shannon$allno)
median(shannon$EMOTUs)
median(shannon$EMno)
median(shannon$SAPOTUs)
median(shannon$SAPno)

# All OTUs GLM
# combine the envi_rabbit with the shannon numbers
shannon = as.data.frame(shannon)
shannonglm = glm(allOTUs ~ PC2* PC1* Core, family=gaussian, data=shannon)
aov.sha=anova(shannonglm)

numglm = glm(allno ~ Core, family=gaussian, data=shannon)
summary(numglm)
shannonlm = lm(allOTUs ~ Core, data=shannon)
numlm = lm(allno ~ Core, data=shannon)
summary(shannonlm)
summary(numlm)

EMshannonlm = lm(EMOTUs ~ Core, data=shannon)
EMnumlm = lm(EMno ~ Core, data=shannon)
summary(EMshannonlm)
summary(EMnumlm)
SAPshannonlm = lm(SAPOTUs ~ Core, data=shannon)
SAPnumlm = lm(SAPno ~ Core, data=shannon)
summary(SAPshannonlm)
summary(SAPnumlm)


