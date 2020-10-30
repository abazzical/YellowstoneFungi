## Rabbit Creek

#setwd("~/Documents/Gannon_Transfer_of_Work/demux/RDA")
library(viridis)
library(cowplot)
library(tidyverse)
library(qiime2R)
library(wesanderson)
names(wes_palettes)
library(vegan)
library(ggrepel)
library(dplyr)
library(ggord)
library(VennDiagram)


pal <- wes_palette("Zissou1", 100, type = "continuous")
br2  = c("#FAD510", "#CB2314", "#273046")
setwd("~/Documents/Gannon_Transfer_of_Work/demux/RDA")
metadata<-read.table("~/Documents/Gannon_Transfer_of_Work/demux/metadataRabbitCreek20184PCA.txt", header = TRUE)
metadata$pH <- factor(metadata$pH,levels = c("<6", "6-8", ">8"))
rownames(metadata)<-rownames(metadata$SampleID)

## metadata only for the cores used in the NMDS
SAP_meta<-read_csv("SAPmetadata.csv")
SAP_meta$ph.bins <- factor(SAP_meta$`pH bins`,levels = c("<6", "6-8", ">8"))
rownames(SAP_meta)<-rownames(SAP_meta$SampleID)

## metadata only for the cores used in the NMDS
EM_meta<-read_csv("EMmetadata.csv")
EM_meta$ph.bins <- factor(EM_meta$`pH bins`,levels = c("<6", "6-8", ">8"))
rownames(EM_meta)<-rownames(SAP_meta$SampleID)

## Venn diagram
thing<-read.csv("~/Documents/Gannon_Transfer_of_Work/demux/OTUnumbers/vennOTUsall.csv")
thingEM<-read.csv("~/Documents/Gannon_Transfer_of_Work/demux/OTUnumbers/vennOTUsEM.csv")
thingSAP<-read.csv("~/Documents/Gannon_Transfer_of_Work/demux/OTUnumbers/vennOTUsSAP.csv")



#NMDS in R with vegan package
#otu table columns are OTUs and rows are soil cores (community)
#use the metaMDS command to run the NMDS

# NMDS --------------------------------------------------------------------
##NMDS from filtered OTU table (whole)
otu_rabbit<-read.table("otu-table.txt", sep = "\t", header = TRUE)

rabbit_df<-as.data.frame(otu_rabbit)
otuID <- rabbit_df$OTUID
rabbit_transpose <- as.data.frame(t(as.matrix(rabbit_df[,-1])))
colnames(rabbit_transpose) <- otuID
rabbit_transpose<-rabbit_transpose[-c(52), ] # remove taxonomy at the bottom row
rabbit <- mutate_all(rabbit_transpose, function(x) as.numeric(as.character(x)))
rownames(rabbit) = rownames(rabbit_transpose)

rabbit_NMDS=metaMDS(rabbit,k=2,trymax=100)
stressplot(rabbit_NMDS)
plot(rabbit_NMDS)
ordiplot(rabbit_NMDS,type="n")
orditorp(rabbit_NMDS,display="species",col="red",air=0.01)
orditorp(rabbit_NMDS,display="sites",cex=1.25,air=0.01)

#NMDS w/o samples from surface water
rabbit_wo_14<-rabbit[-c(47,48),]
rabbit_wo_14_NMDS=metaMDS(rabbit_wo_14,k=2,trymax=100)
stressplot(rabbit_NMDS)
plot(rabbit_wo_14_NMDS)
ordiplot(rabbit_wo_14_NMDS,type="n")
orditorp(rabbit_wo_14_NMDS,display="species",col="red",air=0.01)
orditorp(rabbit_wo_14_NMDS,display="sites",cex=1.25,air=0.01)


## presence/absence NMDS
abs_rabbit<-rabbit
abs_rabbit[abs_rabbit >0]<-1
abs_rabbit_NMDS=metaMDS(abs_rabbit,k=2,trymax=100)
stressplot(abs_rabbit_NMDS)
plot(abs_rabbit_NMDS)
ordiplot(abs_rabbit_NMDS,type="n")
orditorp(abs_rabbit_NMDS,display="species",col="red",air=0.01)
orditorp(abs_rabbit_NMDS,display="sites",cex=1.25,air=0.01)
# plot shows no difference with the other matrix

## metadata only for the cores used in the NMDS
NMDS_meta<-read_csv("NMDSmetadata.csv")
NMDS_meta$ph.bins <- factor(NMDS_meta$`pH bins`,levels = c("<6", "6-8", ">8"))
rownames(NMDS_meta)<-rownames(NMDS_meta$SampleID)


# make a nice looking plot of this NMDS!!
# color with binned pHs? with gradient??
data.scores <- as.data.frame(scores(rabbit_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores


# pH S1 3 clusters  ----------------------------------------------------------------


##ph bins
data.scores$ph.bins <- NMDS_meta$ph.bins  #  add the grp variable created earlier
head(data.scores)  #look at the data
species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data
data.scores$pH <- NMDS_meta$pH  #  add the grp variable created earlier
head(data.scores)  #look at the data

# add polygons to the graph based on ph

grp.a <- data.scores[data.scores$ph.bins == ">8", ][chull(data.scores[data.scores$ph.bins == 
                                                                        ">8", c("NMDS1", "NMDS2")]), ]  # hull values for >8
grp.b <- data.scores[data.scores$ph.bins == "6-8", ][chull(data.scores[data.scores$ph.bins == 
                                                                         "6-8", c("NMDS1", "NMDS2")]), ]  # hull values for 6-8
grp.c <- data.scores[data.scores$ph.bins == "<6", ][chull(data.scores[data.scores$ph.bins == 
                                                                         "<6", c("NMDS1", "NMDS2")]), ]  # hull values for 6-8
hull.data <- rbind(grp.a, grp.b, grp.c)  #combine 
hull.data

#because the ph above 8 is the only bin that matches the community-based cluters, this graph looks pretty stupid


NMDS3clust<- ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=pH), size=5) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis()+
  #scale_color_manual(values = wes_palette("BottleRocket2"))
  
  #coord_equal() +
  theme_bw()
ggsave("NMDS3clust.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches




# ph gradient -------------------------------------------------------------

data.scores$ph <- NMDS_meta$pH  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=ph)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()

# temperature gradient ----------------------------------------------------

data.scores$temp <- NMDS_meta$Temperature  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=temp)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Moisture gradient -------------------------------------------------------

data.scores$moist <- NMDS_meta$Moisture  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=moist)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Phosphorus --------------------------------------------------------------

data.scores$PO4 <- NMDS_meta$PO4  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=PO4)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Total Nitrogen ----------------------------------------------------------

data.scores$TN <- NMDS_meta$TN  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=TN)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Total Carbon ------------------------------------------------------------

data.scores$TC <- NMDS_meta$TC  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=TC)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Cadmium -----------------------------------------------------------------

data.scores$Cd <- NMDS_meta$Cd  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Cd)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Cobalt ------------------------------------------------------------------

data.scores$Co <- NMDS_meta$Co  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Co)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Copper ------------------------------------------------------------------

data.scores$Cu <- NMDS_meta$Cu  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Cu)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Iron --------------------------------------------------------------------

data.scores$Fe <- NMDS_meta$Fe  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Fe)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()



# Manganese ---------------------------------------------------------------

data.scores$Mn <- NMDS_meta$Mn  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Mn)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Nickel ------------------------------------------------------------------

data.scores$Ni <- NMDS_meta$Ni  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Ni)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Lead --------------------------------------------------------------------

data.scores$Pb <- NMDS_meta$Pb  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Pb)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Zinc --------------------------------------------------------------------

data.scores$Zn <- NMDS_meta$Zn  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=Zn)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()
# Indicator species 3 clusters -------------------------------------------------------

# Indicator species analysis pH bins
#install.packages("indicspecies")
library(indicspecies)
pc = read.table("trans_rabbit.txt", header= TRUE)
abund = pc[,4:ncol(pc)]
ph = pc$ph
clust = pc$clust
inv_ph = multipatt(abund, ph, func = "r.g", control = how(nperm=9999))
inv_clust = multipatt(abund, clust, func = "r.g", control = how(nperm=9999))
summary(inv_ph)
summary(inv_clust)


# Violin plots of env. var. with 3 clusters using NMDS_meta ---------------------------------

pH <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$pH))
pH<- pH + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "pH"
  ) 
pH

tempe <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Temperature))
tempe<- tempe + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Temperature (˚C)"
  ) 
tempe
TN <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$TN))
TN<- TN + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Total Nitrogen"
  ) 
TN
TC <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$TC))
TC<- TC + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Total Carbon"
  ) 
TC
moist <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Moisture))
moist<- moist + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Moisture"
  ) 
moist

PO4 <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$PO4))
PO4<- PO4 + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Phosphorus"
  ) 
PO4

# metal measurement units (mg/kg soil)

Cd <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Cd))
Cd<- Cd + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Cadmium (mg/kg soil)"
  ) 
Cd

Co <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Co))
Co<- Co + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Cobalt (mg/kg soil)"
  ) 
Co

Cu <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Cu))
Cu<- Cu + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Copper (mg/kg soil)"
  ) 
Cu

Fe <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Fe))
Fe<- Fe + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Iron (mg/kg soil)"
  ) 
Fe

Mn <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Mn))
Mn<- Mn + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Manganese (mg/kg soil)"
  ) 
Mn

Ni <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Ni))
Ni<- Ni + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Nickel (mg/kg soil)"
  ) 
Ni

Pb <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Pb))
Pb<- Pb + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Lead (mg/kg soil)"
  ) 
Pb

Zn <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$Zn))
Zn<- Zn + geom_violin() + geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Zinc (mg/kg soil)"
  ) 
Zn

prow<-plot_grid(pH, tempe, TN, TC, PO4, Cd, Co, Cu, Fe, Mn, Ni, Pb, Zn, 
                labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M'), 
                ncol = 5,
                label_size = 16)
prow
ggsave("plotViolin.pdf", height=10, width=12, device="pdf") # save a PDF 3 inches by 4 inches



numOTUs <- ggplot(NMDS_meta, aes(factor(NMDS_meta$cluster), NMDS_meta$`# OTUs`))
numOTUs<- numOTUs + geom_violin(trim=FALSE) + geom_jitter(height = 0, width = 0.1) +
  geom_boxplot(width = 0.05) +
  theme_bw() +
  labs(
    x = "Clusters",
    y = "Number of OTUs per core"
  ) 
numOTUs
ggsave("plotNumOTUs.pdf", height=4, width=7, device="pdf") # save a PDF 3 inches by 4 inches


# Violin plots of numbers of OTUs per core --------------------------------
OTUnums<-read_csv("OTUnumsViolin.csv")
allOTUs<-ggplot(OTUnums, aes(x=clusts, y=All)) + 
  scale_fill_manual(values = c("#273046","#FAD510"))+
  geom_violin(trim = FALSE, aes(fill = factor(clusts))) + 
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="white") + 
  theme_bw()+ 
  xlab("")+
  ylim(0,80)+
  ylab("Number of OTUs in a soil core")
allOTUs

EMOTUs<-ggplot(OTUnums, aes(x=clusts, y=EM)) + 
  scale_fill_manual(values = c("#273046","#FAD510"))+
  geom_violin(trim = FALSE, aes(fill = factor(clusts))) + 
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="white") + 
  theme_bw()+ 
  xlab("")+
  ylim(0,12)+
  ylab("Number of EM OTUs in a soil core")
EMOTUs

SAPOTUs<-ggplot(OTUnums, aes(x=clusts, y=Saprotrophic)) + 
  scale_fill_manual(values = c("#273046","#FAD510"))+
  geom_violin(trim = FALSE, aes(fill = factor(clusts))) + 
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange",width=0.1, col="white") + 
  theme_bw()+ 
  xlab("")+
  ylim(0,25)+
  ylab("Number of saprotrophic OTUs in a soil core")
SAPOTUs

legend_b <- get_legend(
  allOTUs + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text( size = 16))+
    theme(legend.title = element_text( size = 16))
)

sprow<-plot_grid(allOTUs+ theme(legend.position="none"), EMOTUs+ theme(legend.position="none"), SAPOTUs+ theme(legend.position="none"),labels = c('A', 'B','C'), ncol = 3, label_size = 16)

plot_grid(sprow)
ggsave("plotsOTUsnum.pdf", height=4, width=6, device="pdf") # save a PDF 3 inches by 4 inches

group1 <- subset(OTUnums, clusts=="pH >8")
group2 <- subset(OTUnums, clusts=="pH <8")
t.test(group1$All, group2$All)
t.test(group1$EM, group2$EM)
t.test(group1$Saprotrophic, group2$Saprotrophic)

#Bonferroni correction of 3 tests = looking for p-values <0.01667

# Exclude Sites with unIDed taxa clustering -------------------------------

# the right-most group is determined by taxa that are not identified by UNITE
# NMDS excluding sites 5,6,7,8,9, 60, 63, 64, 67,69

#NMDS w/o samples clustering because of taxa with no IDs
rabbit_wo_ID<-rabbit[-c(39,43,44,45,46,47,48,49,50,51),]
rabbit_wo_ID_NMDS=metaMDS(rabbit_wo_ID,k=2,trymax=100)
stressplot(rabbit_wo_ID_NMDS)
plot(rabbit_wo_ID_NMDS)
ordiplot(rabbit_wo_ID_NMDS,type="n")
#orditorp(rabbit_wo_ID_NMDS,display="species",col="red",air=0.01)
orditorp(rabbit_wo_ID_NMDS,display="sites",cex=1.25,air=0.01)

# make a nice looking plot of this NMDS!!
# color with binned pHs? with gradient??
data.scores_wo_ID_ <- as.data.frame(scores(rabbit_wo_ID_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores_wo_ID_$site <- rownames(data.scores_wo_ID_)  # create a column of site names, from the rownames of data.scores


# pH gradient 2 clusters ----------------------------------------------------------------


##ph bins
data.scores_wo_ID_$ph.bins <- SAP_meta$ph.bins  #  add the grp variable created earlier
head(data.scores_wo_ID_)  #look at the data
species.scores_wo_ID_ <- as.data.frame(scores(rabbit_wo_ID_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
data.scores_wo_ID_$clusts <- SAP_meta$clusts  #  add the grp variable created earlier
head(data.scores_wo_ID_)  #look at the data


clust2ALL<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores_wo_ID_,aes(x=NMDS1,y=NMDS2,colour=pH),size=4) + # add the point markers
  #geom_text(data=data.scores_wo_ID_,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis()+
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()
ggsave("NMDS2clust.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

vennALL<-venn.diagram(
  x = list(
    thing %>% filter(cluster=="Group 1") %>% select(OTUID) %>% unlist() , 
    thing %>% filter(cluster=="Group 2") %>% select(OTUID) %>% unlist() 
  ),
  category.names = c("pH >8" , "pH <8"),
  filename = 'vennClust.png',
  output = TRUE ,
  imagetype = "png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 500,
  compression = "lzw",
  lwd = 1,
  col=c("#FAD510", "#273046"),
  fill = c(alpha("#FAD510",0.3), alpha("#273046",0.3)),
  cex = 0.75,
  fontfamily = "sans",
  cat.cex = 0.75,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#FAD510", "#273046"),
  #rotation = 1
)

FIG2<-plot_grid(clust2ALL,grobTree(vennALL), 
                labels = "AUTO", 
                ncol = 2,
                label_size = 16)

ggsave("Fig2_allOTUs.pdf", height=10, width=12, device="pdf") # save a PDF 3 inches by 4 inches


# Indicator species 2 clusters -------------------------------------------------------

# Indicator species analysis pH bins
pc = read.table("trans_rabbit_2clust.txt", header= TRUE)
abund = pc[,4:ncol(pc)]
ph = pc$ph
clust = pc$clust
inv_ph = multipatt(abund, ph, func = "r.g", control = how(nperm=9999))
inv_clust = multipatt(abund, clust, func = "r.g", control = how(nperm=9999))
summary(inv_ph)
summary(inv_clust)


# MANOVA of 2 clusters ----------------------------------------------------


clust2man<-manova(cbind(data.scores_wo_ID_$NMDS1,data.scores_wo_ID_$NMDS2) ~ as.factor(data.scores_wo_ID_$clusts), data=data.scores_wo_ID_, subset=as.factor(data.scores_wo_ID_$clusts) %in% c("a","b"))

summary(clust2man)

# Guilds NMDS -------------------------------------------------------------

##NMDS from filtered OTU table (EM fungi)
otu_EM<-read.table("EM-table.txt", sep = "\t", header = TRUE)

EM_df<-as.data.frame(otu_EM)
otuID <- EM_df$OTUID
EM_transpose <- as.data.frame(t(as.matrix(EM_df[,-1])))
colnames(EM_transpose) <- otuID
#EM_transpose<-EM_transpose[-c(52), ] # remove taxonomy at the bottom row
EM <- mutate_all(EM_transpose, function(x) as.numeric(as.character(x)))
rownames(EM) = rownames(EM_transpose)

EM_NMDS=metaMDS(EM,k=2,trymax=1000)
stressplot(EM_NMDS)
plot(EM_NMDS)
ordiplot(EM_NMDS,type="n")
orditorp(EM_NMDS,display="species",col="red",air=0.01)
orditorp(EM_NMDS,display="sites",cex=1.25,air=0.01)

EM_transpose_woS42<-EM_transpose[-c(35), ] # remove Site 42 which only has ceratobasidium as the only taxon and uniquely has that taxon
EM <- mutate_all(EM_transpose_woS42, function(x) as.numeric(as.character(x)))
rownames(EM) = rownames(EM_transpose_woS42)

EM_NMDS_woS42=metaMDS(EM,k=2,trymax=100)
stressplot(EM_NMDS_woS42)
plot(EM_NMDS_woS42)
ordiplot(EM_NMDS_woS42,type="n")
orditorp(EM_NMDS_woS42,display="species",col="red",air=0.01)
orditorp(EM_NMDS_woS42,display="sites",cex=1.25,air=0.01)

##NMDS from filtered OTU table (SAP fungi)
otu_SAP<-read.table("SAP-table.txt", sep = "\t", header = TRUE)

SAP_df<-as.data.frame(otu_SAP)
otuID <- SAP_df$OTUID
SAP_transpose <- as.data.frame(t(as.matrix(SAP_df[,-1])))
colnames(SAP_transpose) <- otuID
#EM_transpose<-EM_transpose[-c(52), ] # remove taxonomy at the bottom row
SAP <- mutate_all(SAP_transpose, function(x) as.numeric(as.character(x)))
rownames(SAP) = rownames(SAP_transpose)

SAP_NMDS=metaMDS(SAP,k=2,trymax=100)
stressplot(SAP_NMDS)
plot(SAP_NMDS)
ordiplot(SAP_NMDS,type="n")
orditorp(SAP_NMDS,display="species",col="red",air=0.01)
orditorp(SAP_NMDS,display="sites",cex=1.25,air=0.01)


# Guild EM graphs w metadata -------------------------------------------------
## metadata only for the cores used in the NMDS
EM_meta<-read_csv("EMmetadata.csv")
EM_meta$ph.bins <- factor(EM_meta$`pH bins`,levels = c("<6", "6-8", ">8"))
rownames(EM_meta)<-rownames(EM_meta$SampleID)

# make a nice looking plot of this NMDS!!
# color with binned pHs? with gradient??
data.scoresEM <- as.data.frame(scores(EM_NMDS_woS42))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scoresEM$site <- rownames(data.scoresEM)  # create a column of site names, from the rownames of data.scores


# pH bins  ----------------------------------------------------------------


##ph bins
data.scoresEM$ph.bins <- EM_meta$ph.bins  #  add the grp variable created earlier
head(data.scoresEM)  #look at the data
species.scoresEM <- as.data.frame(scores(EM, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scoresEM$species <- rownames(species.scoresEM)  # create a column of species, from the rownames of species.scores
head(species.scoresEM)  #look at the data

# add polygons to the graph based on ph

grp.aEM <- data.scoresEM[data.scoresEM$ph.bins == ">8", ][chull(data.scoresEM[data.scoresEM$ph.bins == 
                                                                        ">8", c("NMDS1", "NMDS2")]), ]  # hull values for >8
grp.bEM <- data.scoresEM[data.scoresEM$ph.bins == "6-8", ][chull(data.scoresEM[data.scoresEM$ph.bins == 
                                                                         "6-8", c("NMDS1", "NMDS2")]), ]  # hull values for 6-8
grp.cEM <- data.scoresEM[data.scoresEM$ph.bins == "<6", ][chull(data.scoresEM[data.scoresEM$ph.bins == 
                                                                        "<6", c("NMDS1", "NMDS2")]), ]  # hull values for 6-8
hull.dataEM <- rbind(grp.aEM, grp.bEM, grp.cEM)  #combine 
hull.dataEM

#because the ph above 8 is the only bin that matches the community-based cluters, this graph looks pretty stupid


ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,shape=ph.bins,colour=ph.bins), size=5) + # add the point markers
  geom_text(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()

ggplot() + 
  geom_polygon(data=hull.dataEM,aes(x=NMDS1,y=NMDS2,fill=ph.bins,group=ph.bins),alpha=0.30) + # add the convex hulls
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,shape=ph.bins,colour=ph.bins),size=4) + # add the point markers
  scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


# ph gradient and Venn EM-------------------------------------------------------------

data.scoresEM$ph <- EM_meta$pH  #  add the grp variable created earlier
head(data.scoresEM)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=ph),size=4) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()
ggsave("NMDS2clustEM.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

venn.diagram(
  x = list(
    thingEM %>% filter(cluster=="Group 1") %>% select(OTUID) %>% unlist() , 
    thingEM %>% filter(cluster=="Group 2") %>% select(OTUID) %>% unlist() 
  ),
  category.names = c("pH >8" , "pH <8"),
  filename = 'vennClustEM.png',
  output = TRUE ,
  imagetype = "png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 500,
  compression = "lzw",
  lwd = 1,
  col=c("#FAD510", "#273046"),
  fill = c(alpha("#FAD510",0.3), alpha("#273046",0.3)),
  cex = 0.75,
  fontfamily = "sans",
  cat.cex = 0.75,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#FAD510", "#273046"),
  #rotation = 1
)
data.scoresEM$clusts <- EM_meta$clusts  #  add the grp variable created earlier
head(data.scoresEM)  #look at the data

clust2EM<-manova(cbind(data.scoresEM$NMDS1,data.scoresEM$NMDS2) ~ as.factor(data.scoresEM$clusts), data=data.scoresEM, subset=as.factor(data.scoresEM$clusts) %in% c("a","b"))

summary(clust2EM)


# temperature gradient EM----------------------------------------------------

data.scoresEM$temp <- EM_meta$Temperature  #  add the grp variable created earlier
head(data.scoresEM)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=temp)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Moisture gradient EM-------------------------------------------------------

data.scoresEM$moist <- EM_meta$Moisture  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=moist)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Phosphorus EM--------------------------------------------------------------

data.scoresEM$PO4 <- EM_meta$PO4  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=PO4)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Total Nitrogen EM----------------------------------------------------------

data.scoresEM$TN <- EM_meta$TN  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=TN)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Total Carbon EM------------------------------------------------------------

data.scoresEM$TC <- EM_meta$TC  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=TC)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Cadmium EM-----------------------------------------------------------------

data.scoresEM$Cd <- EM_meta$Cd  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=Cd)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Cobalt EM------------------------------------------------------------------

data.scoresEM$Co <- EM_meta$Co  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=Co)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Copper EM------------------------------------------------------------------

data.scoresEM$Cu <- EM_meta$Cu  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=Cu)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Iron EM--------------------------------------------------------------------

data.scoresEM$Fe <- EM_meta$Fe  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=Fe)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()



# Manganese EM---------------------------------------------------------------

data.scoresEM$Mn <- EM_meta$Mn  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=Mn)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Nickel EM------------------------------------------------------------------

data.scoresEM$Ni <- EM_meta$Ni  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=Ni)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Lead EM--------------------------------------------------------------------

data.scoresEM$Pb <- EM_meta$Pb  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=Pb)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Zinc EM --------------------------------------------------------------------

data.scoresEM$Zn <- EM_meta$Zn  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresEM,aes(x=NMDS1,y=NMDS2,colour=Zn)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Guild SAP metadata ------------------------------------------------------



# make a nice looking plot of this NMDS!!
# color with binned pHs? with gradient??
data.scoresSAP <- as.data.frame(scores(SAP_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scoresSAP$site <- rownames(data.scoresSAP)  # create a column of site names, from the rownames of data.scores


# pH bins  ----------------------------------------------------------------


##ph bins
data.scoresSAP$ph.bins <- SAP_meta$ph.bins  #  add the grp variable created earlier
head(data.scoresSAP)  #look at the data
species.scoresSAP <- as.data.frame(scores(SAP, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scoresSAP$species <- rownames(species.scoresSAP)  # create a column of species, from the rownames of species.scores
head(species.scoresSAP)  #look at the data

# add polygons to the graph based on ph

grp.aSAP <- data.scoresSAP[data.scoresSAP$ph.bins == ">8", ][chull(data.scoresSAP[data.scoresSAP$ph.bins == 
                                                                                    ">8", c("NMDS1", "NMDS2")]), ]  # hull values for >8
grp.bSAP <- data.scoresSAP[data.scoresSAP$ph.bins == "6-8", ][chull(data.scoresSAP[data.scoresSAP$ph.bins == 
                                                                                     "6-8", c("NMDS1", "NMDS2")]), ]  # hull values for 6-8
grp.cSAP <- data.scoresSAP[data.scoresSAP$ph.bins == "<6", ][chull(data.scoresSAP[data.scoresSAP$ph.bins == 
                                                                                    "<6", c("NMDS1", "NMDS2")]), ]  # hull values for 6-8
hull.dataSAP <- rbind(grp.aSAP, grp.bSAP, grp.cSAP)  #combine 
hull.dataSAP

#because the ph above 8 is the only bin that matches the community-based cluters, this graph looks pretty stupid


ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,shape=ph.bins,colour=ph.bins)) + # add the point markers
  geom_text(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()

br2 =wes_palette("BottleRocket2", 18, type="continuous")
ggplot() + 
  geom_polygon(data=hull.dataSAP,aes(x=NMDS1,y=NMDS2,group=ph.bins),alpha=0.20) + # add the convex hulls
  
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=ph.bins),size=4) + # add the point markers
  scale_color_manual(values = wes_palette("BottleRocket2")) +
  coord_equal() +
  theme_bw() + 
  #labs(colour = "pH", )+
  theme(axis.text.x = element_text(size=18),  # size of the axis text
        axis.text.y = element_text(size=18), #
        #axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # size x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


# ph gradient and Venn SAP-------------------------------------------------------------

data.scoresSAP$ph <- SAP_meta$pH  #  add the grp variable created earlier
head(data.scoresSAP)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=ph), size=4) + # add the point markers
  #geom_text(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()
ggsave("NMDS2clustSAP.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

vennSAP<-venn.diagram(
  x = list(
    thingSAP %>% filter(cluster=="Group 1") %>% select(OTUID) %>% unlist() , 
    thingSAP %>% filter(cluster=="Group 2") %>% select(OTUID) %>% unlist() 
  ),
  category.names = c("pH >8" , "pH <8"),
  filename = 'vennClustSAP.png',
  output = TRUE ,
  imagetype = "png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 500,
  compression = "lzw",
  lwd = 1,
  col=c("#FAD510", "#273046"),
  fill = c(alpha("#FAD510",0.3), alpha("#273046",0.3)),
  cex = 0.75,
  fontfamily = "sans",
  cat.cex = 0.75,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#FAD510", "#273046"),
  #rotation = 1
)

data.scoresSAP$clusts <- SAP_meta$clusts  #  add the grp variable created earlier
head(data.scoresSAP)  #look at the data

clust2SAP<-manova(cbind(data.scoresSAP$NMDS1,data.scoresSAP$NMDS2) ~ as.factor(data.scoresSAP$clusts), data=data.scoresSAP, subset=as.factor(data.scoresSAP$clusts) %in% c("a","b"))

summary(clust2SAP)


# temperature gradient SAP----------------------------------------------------

data.scoresSAP$temp <- SAP_meta$Temperature  #  add the grp variable created earlier
head(data.scoresSAP)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=temp)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Moisture gradient SAP-------------------------------------------------------

data.scoresSAP$moist <- SAP_meta$Moisture  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=moist)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Phosphorus SAP--------------------------------------------------------------

data.scoresSAP$PO4 <- SAP_meta$PO4  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=PO4)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Total Nitrogen SAP----------------------------------------------------------

data.scoresSAP$TN <- SAP_meta$TN  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=TN)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Total Carbon SAP------------------------------------------------------------

data.scoresSAP$TC <- SAP_meta$TC  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=TC)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Cadmium SAP-----------------------------------------------------------------

data.scoresSAP$Cd <- SAP_meta$Cd  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=Cd)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Cobalt SAP------------------------------------------------------------------

data.scoresSAP$Co <- SAP_meta$Co  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=Co)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Copper SAP------------------------------------------------------------------

data.scoresSAP$Cu <- SAP_meta$Cu  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=Cu)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Iron SAP--------------------------------------------------------------------

data.scoresSAP$Fe <- SAP_meta$Fe  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=Fe)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()



# Manganese SAP---------------------------------------------------------------

data.scoresSAP$Mn <- SAP_meta$Mn  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=Mn)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Nickel SAP------------------------------------------------------------------

data.scoresSAP$Ni <- SAP_meta$Ni  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=Ni)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Lead SAP--------------------------------------------------------------------

data.scoresSAP$Pb <- SAP_meta$Pb  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=Pb)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()


# Zinc SAP --------------------------------------------------------------------

data.scoresSAP$Zn <- SAP_meta$Zn  #  add the grp variable created earlier
head(data.scores)  #look at the data
#species.scores <- as.data.frame(scores(rabbit_NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scoresSAP,aes(x=NMDS1,y=NMDS2,colour=Zn)) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site)) +  # add the site labels
  scale_color_viridis() +
  #scale_color_manual(values = wes_palette("BottleRocket2")) +
  #coord_equal() +
  theme_bw()






# Soil chemistry ----------------------------------------------------------


## SOIL chemistry PCoA
#missing values for samples 1.1, 3.2, 5.1, 10.5 - removed for PCA
soils <- read.csv("~/Documents/Gannon_Transfer_of_Work/demux/Yellowstone_Soils_temp.csv")
row.names(soils) <- soils$SampleID
#pcadata <- scale(soils[c(17)])

#colnames(soils[-c(1,12,21,50),c(2:14,17)])

#soilpca <- prcomp(soils[-c(1,12,21,50), c(2:14,17)], scale. = TRUE) #SE
soilpca <- prcomp(soils[ ,3:16], scale. = TRUE) #SE
biplot(soilpca, xlab = "PC1: 34% of variance", ylab = "PC2: 24% of variance", col = c("navy", "deeppink3"), pc.biplot = FALSE)
summary(soilpca)

#vecs <- eigen(cov(pcadata))$vectors
## NOTE: Z <- X %*% V # principal components are the matrix * eigenvectors on each column
#PCs <- pcadata %*% vecs

plotcoords=data.frame(soilpca$x)
envloading1=data.frame(soilpca$rotation)

envloading <- envloading1*10

pcadata <- soils[c(17), ]
soils$pH <- factor(soils$pH,levels = c("<6", "6-8", ">8"))
# prettier plot
PCAsoil<-ggplot()+
  geom_segment(aes(x=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0) ,y=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0), xend= c(envloading[1,1],envloading[2,1],envloading[3,1],envloading[4,1],envloading[5,1],envloading[6,1],envloading[7,1],envloading[8,1],envloading[9,1],envloading[10,1],envloading[11,1],envloading[12,1],envloading[13,1],envloading[14,1]), yend=c(envloading[1,2],envloading[2,2],envloading[3,2],envloading[4,2],envloading[5,2],envloading[6,2],envloading[7,2],envloading[8,2],envloading[9,2],envloading[10,2],envloading[11,2],envloading[12,2],envloading[13,2],envloading[14,2])), color="#999999", arrow=arrow(angle = 20, length = unit(0.1,"cm"), ends = "last", type = "open"), size = 1)+
  scale_color_viridis()+
  #scale_color_manual(values = wes_palette("BottleRocket2", 100, type="continuous"))+
  geom_point(aes(x=plotcoords$PC1, y=plotcoords$PC2, color = soils$pH_), size = 4)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("PC 1 -- 34% of variance")+
  ylab("PC 2 -- 24% of variance")+
  labs(colour = "pH")+
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
ggsave("PCoAsoilpH_cont.pdf", height=4, width=6.5, device="pdf") # save a PDF 3 inches by 4 inches



# In vitro ph and temperature culture growth ------------------------------


phYNP <- read.csv("pH.csv")
phYNP$abs_mean <- apply(X = phYNP[ ,7:10], MARGIN = 1, FUN = mean) # the command apply will average the columns with the radius measurements
#R is usually more reliable than excel to make measurements, excel will do funny things sometimes!

#YNP1 - Agaricus, your culture!
YNP1 <- phYNP[which(phYNP$isolate == "YNP1"), ] # only consider the isolates YNP1 from the spreadsheet
YNP1p4 <- YNP1[which(YNP1$pH == "4"), ] #of the YNP1, select the ones grown on pH 4
avg_YNP1p4 <- aggregate(YNP1p4$abs_mean, by = list(YNP1p4$day_measured), FUN = mean) #get an average of the replicates to be used to plot the line
YNP1p7 <- YNP1[which(YNP1$pH == "7"), ]
avg_YNP1p7 <- aggregate(YNP1p7$abs_mean, by = list(YNP1p7$day_measured), FUN = mean)
YNP1p9 <- YNP1[which(YNP1$pH == "9"), ]
avg_YNP1p9 <- aggregate(YNP1p9$abs_mean, by = list(YNP1p9$day_measured), FUN = mean)
pH<-as.character(phYNP$pH) # can't remember if this is actually useful or needed
phYNP1plot<-ggplot(YNP1) +
  scale_color_manual(values = wes_palette("FantasticFox1")) + #you can choose whatever colors you prefer, I like Wes Anderson palettes, and Bottle Rocket has nice distinctive colors - I like the Zissou palette, but they are difficult to see for color blind folk
  geom_point(aes(day_measured, abs_mean, color = as.factor(pH)), size = 4) + # plots points
  geom_line(aes(Group.1, x), avg_YNP1p4)+ #plots lines
  geom_line(aes(Group.1, x), avg_YNP1p7)+
  geom_line(aes(Group.1, x), avg_YNP1p9)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("Day Measured")+
  ylab("Fungus Radius (mm)")+
  labs(colour = "pH")
#To save a graph only with your species, use something like this:
# + ggsave("plotAgaricusPHseamus.pdf", height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#YNP6 - Pisolithus, the dog turd fungus
YNP6 <- phYNP[which(phYNP$isolate == "YNP6"), ]
YNP6p4 <- YNP6[which(YNP6$pH == "4"), ]
avg_YNP6p4 <- aggregate(YNP6p4$abs_mean, by = list(YNP6p4$day_measured), FUN = mean)
YNP6p7 <- YNP6[which(YNP6$pH == "7"), ]
avg_YNP6p7 <- aggregate(YNP6p7$abs_mean, by = list(YNP6p7$day_measured), FUN = mean)
YNP6p9 <- YNP6[which(YNP6$pH == "9"), ]
avg_YNP6p9 <- aggregate(YNP6p9$abs_mean, by = list(YNP6p9$day_measured), FUN = mean)
pH<-as.character(phYNP$pH)
phYNP6plot<-ggplot(YNP6) +
  scale_color_manual(values = wes_palette("FantasticFox1")) +
  geom_point(aes(day_measured, abs_mean, color = as.factor(pH)), size = 4) + 
  geom_line(aes(Group.1, x), avg_YNP6p4)+
  geom_line(aes(Group.1, x), avg_YNP6p7)+
  geom_line(aes(Group.1, x), avg_YNP6p9)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("Day Measured")+
  ylab("Fungus Radius (mm)")+
  labs(colour = "pH")

#YNP56 - a soil mold
YNP56 <- phYNP[which(phYNP$isolate == "YNP56"), ]
YNP56p4 <- YNP56[which(YNP56$pH == "4"), ]
avg_YNP56p4 <- aggregate(YNP56p4$abs_mean, by = list(YNP56p4$day_measured), FUN = mean)
YNP56p7 <- YNP56[which(YNP56$pH == "7"), ]
avg_YNP56p7 <- aggregate(YNP56p7$abs_mean, by = list(YNP56p7$day_measured), FUN = mean)
YNP56p9 <- YNP56[which(YNP56$pH == "9"), ]
avg_YNP56p9 <- aggregate(YNP56p9$abs_mean, by = list(YNP56p9$day_measured), FUN = mean)
pH<-as.character(phYNP$pH)
phYNP56plot<-ggplot(YNP56) +
  scale_color_manual(values = wes_palette("FantasticFox1")) +
  geom_point(aes(day_measured, abs_mean, color = as.factor(pH)), size = 4) + 
  geom_line(aes(Group.1, x), avg_YNP56p4)+
  geom_line(aes(Group.1, x), avg_YNP56p7)+
  geom_line(aes(Group.1, x), avg_YNP56p9)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("Day Measured")+
  ylab("Fungus Radius (mm)")+
  labs(colour = "pH")

#YNP67 - another soil mold
YNP67 <- phYNP[which(phYNP$isolate == "YNP67"), ]
YNP67p4 <- YNP67[which(YNP67$pH == "4"), ]
avg_YNP67p4 <- aggregate(YNP67p4$abs_mean, by = list(YNP67p4$day_measured), FUN = mean)
YNP67p7 <- YNP67[which(YNP67$pH == "7"), ]
avg_YNP67p7 <- aggregate(YNP67p7$abs_mean, by = list(YNP67p7$day_measured), FUN = mean)
YNP67p9 <- YNP67[which(YNP67$pH == "9"), ]
avg_YNP67p9 <- aggregate(YNP67p9$abs_mean, by = list(YNP67p9$day_measured), FUN = mean)
pH<-as.character(phYNP$pH)
phYNP67plot<-ggplot(YNP67) +
  scale_color_manual(values = wes_palette("FantasticFox1")) +
  geom_point(aes(day_measured, abs_mean, color = as.factor(pH)),size=4) + 
  geom_line(aes(Group.1, x), avg_YNP67p4)+
  geom_line(aes(Group.1, x), avg_YNP67p7)+
  geom_line(aes(Group.1, x), avg_YNP67p9)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("Day Measured")+
  ylab("Fungus Radius (mm)")+
  labs(colour = "pH")

legend_b <- get_legend(
  phYNP6plot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text( size = 16))+
    theme(legend.title = element_text( size = 16))
)
prowPH<-plot_grid(phYNP1plot+ theme(legend.position="none"), phYNP6plot+ theme(legend.position="none"), phYNP56plot+ theme(legend.position="none"), phYNP67plot + theme(legend.position="none"), labels = c('A', 'B', 'C', 'D'), label_size = 16)

plot_grid(prowPH, legend_b, ncol = 1, rel_heights = c(1, .1))
ggsave("plotsCultPH.pdf", height=7, width=6, device="pdf") # save a PDF 3 inches by 4 inches



# 

#===================================================================================
# TEMPERATURE

tempYNP <- read.csv("temp.csv")

#YNP1
YNP1 <- tempYNP[which(tempYNP$isolate == "YNP1"), ]
YNP1t28 <- YNP1[which(YNP1$temp == "30"), ]
avg_YNP1t28 <- aggregate(YNP1t28$area_measured, by = list(YNP1t28$day_measured), FUN = mean)
YNP1t35 <- YNP1[which(YNP1$temp == "35"), ]
avg_YNP1t35 <- aggregate(YNP1t35$area_measured, by = list(YNP1t35$day_measured), FUN = mean)
YNP1t40 <- YNP1[which(YNP1$temp == "40"), ]
avg_YNP1t40 <- aggregate(YNP1t40$area_measured, by = list(YNP1t40$day_measured), FUN = mean)
temp<-as.character(temp)
YNP1plot<-ggplot(YNP1) +
  scale_color_manual(values = wes_palette("GrandBudapest1")) +
  geom_point(aes(day_measured, area_measured, color = as.factor(temp)), size = 4) + 
  geom_line(aes(Group.1, x), avg_YNP1t28)+
  geom_line(aes(Group.1, x), avg_YNP1t35)+
  geom_line(aes(Group.1, x), avg_YNP1t40)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("Day Measured")+
  ylab("Fungus Area (mm)")+
  labs(colour = "Temperature")

#YNP6
YNP6 <- tempYNP[which(tempYNP$isolate == "YNP6"), ]
ynp6t28 <- YNP6[which(YNP6$temp == "30"), ]
avg_ynp6t28 <- aggregate(ynp6t28$area_measured, by = list(ynp6t28$day_measured), FUN = mean)
ynp6t35 <- YNP6[which(YNP6$temp == "35"), ]
avg_ynp6t35 <- aggregate(ynp6t35$area_measured, by = list(ynp6t35$day_measured), FUN = mean)
ynp6t40 <- YNP6[which(YNP6$temp == "40"), ]
avg_ynp6t40 <- aggregate(ynp6t40$area_measured, by = list(ynp6t40$day_measured), FUN = mean)
temp<-as.character(temp)
YNP6plot<-ggplot(YNP6) +
  scale_color_manual(values = wes_palette("GrandBudapest1")) +
  geom_point(aes(day_measured, area_measured, color = as.factor(temp)), size = 4) + 
  geom_line(aes(Group.1, x), avg_ynp6t28)+
  geom_line(aes(Group.1, x), avg_ynp6t35)+
  geom_line(aes(Group.1, x), avg_ynp6t40)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("Day Measured")+
  ylab("Fungus Area (mm)")+
  labs(colour = "Temperature")

#YNP56
YNP56 <- tempYNP[which(tempYNP$isolate == "YNP56"), ]
ynp56t28 <- YNP56[which(YNP56$temp == "30"), ]
avg_ynp56t28 <- aggregate(ynp56t28$area_measured, by = list(ynp56t28$day_measured), FUN = mean)
ynp56t35 <- YNP56[which(YNP56$temp == "35"), ]
avg_ynp56t35 <- aggregate(ynp56t35$area_measured, by = list(ynp56t35$day_measured), FUN = mean)
ynp56t40 <- YNP56[which(YNP56$temp == "40"), ]
avg_ynp56t40 <- aggregate(ynp56t40$area_measured, by = list(ynp56t40$day_measured), FUN = mean)
temp<-as.character(temp)
YNP56plot<-ggplot(YNP56) +
  scale_color_manual(values = wes_palette("GrandBudapest1")) +
  geom_point(aes(day_measured, area_measured, color = as.factor(temp)), size = 4) + 
  geom_line(aes(Group.1, x), avg_ynp56t28)+
  geom_line(aes(Group.1, x), avg_ynp56t35)+
  geom_line(aes(Group.1, x), avg_ynp56t40)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("Day Measured")+
  ylab("Fungus Area (mm)")+
  labs(colour = "Temperature")

#YNP67
YNP67 <- tempYNP[which(tempYNP$isolate == "YNP67"), ]
ynp67t28 <- YNP67[which(YNP67$temp == "30"), ]
avg_ynp67t28 <- aggregate(ynp67t28$area_measured, by = list(ynp67t28$day_measured), FUN = mean)
ynp67t35 <- YNP67[which(YNP67$temp == "35"), ]
avg_ynp67t35 <- aggregate(ynp67t35$area_measured, by = list(ynp67t35$day_measured), FUN = mean)
ynp67t40 <- YNP67[which(YNP67$temp == "40"), ]
avg_ynp67t40 <- aggregate(ynp67t40$area_measured, by = list(ynp67t40$day_measured), FUN = mean)
temp<-as.character(temp)
YNP67plot<-ggplot(YNP67) +
  scale_color_manual(values = wes_palette("GrandBudapest1")) +
  geom_point(aes(day_measured, area_measured, color = as.factor(temp)), size = 4) + 
  geom_line(aes(Group.1, x), avg_ynp67t28)+
  geom_line(aes(Group.1, x), avg_ynp67t35)+
  geom_line(aes(Group.1, x), avg_ynp67t40)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("Day Measured")+
  ylab("Fungus Area (mm)")+
  labs(colour = "Temperature")

legend_b <- get_legend(
  YNP6plot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text( size = 16))+
    theme(legend.title = element_text( size = 16))
)
prowTemp<-plot_grid(YNP1plot+ theme(legend.position="none"), YNP6plot+ theme(legend.position="none"), YNP56plot+ theme(legend.position="none"), YNP67plot+ theme(legend.position="none") + theme(legend.position="none"), labels = c('A', 'B', 'C', 'D'), label_size = 16)

plot_grid(prowTemp, legend_b, ncol = 1, rel_heights = c(1, .1))
ggsave("plotsTempYNP.pdf", height=7, width=6, device="pdf") #

