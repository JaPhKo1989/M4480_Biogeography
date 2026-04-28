
library(vegan)
library(terra)
library(rnaturalearth)
library(robust)
library(qvalue)
library(ggplot2)



setwd("~/Schreibtisch/RDA-landscape-genomics-main./RDA-landscape-genomics-main")

Genotypes <- read.table("./Data/Pine_AllNatural_GCandTotemIndivs_GWAS_SNPs_June8th2019.txt", header = T)

## Function to transform the alleles into counts
genotypes2geno <- function(vec) {
  vec <- as.character(vec)
  lev <- levels(as.factor(vec))[-which(levels(as.factor(vec))=="00")]
  if(length(lev)==3){
    vec[vec==lev[1]] <- 0
    vec[vec==lev[2]] <- 1
    vec[vec==lev[3]] <- 2
    vec[vec=="00"] <- NA
  }
  else if(length(lev)==2){
    vec[!vec%in%c("AC","AG","AT","CG","CT","GT","00")] <- 0
    vec[vec%in%c("AC","AG","AT","CG","CT","GT")] <- 1
    vec[vec=="00"] <- NA
  }
  return(as.integer(vec))
}

## Running the function on the lodgepole pine dataset 
names_ind <- row.names(Genotypes)
Genotypes <- as.data.frame(apply(Genotypes, 2, genotypes2geno))
row.names(Genotypes) <- names_ind

## Loading samples metadata 
InfoInd <- read.table("./Data/Pine_TotemField_AllNaturalIndivResiduals_Jan20th2017.csv", sep = ",", header = T)
Genotypes <- Genotypes[match(InfoInd$Internal_ID, row.names(Genotypes), nomatch = 0),]

## Estimating population allele frequencies
AllFreq <- aggregate(Genotypes, by = list(InfoInd$ProvSeedlotCode), function(x) mean(x, na.rm = T)/2)
row.names(AllFreq) <- as.character(AllFreq$Group.1)

## Delete loci with missing data in > 12/281 populations
na_pop <- apply(AllFreq[,-1], 2, function(x) sum(is.na(x)))
AllFreq <- AllFreq[,(which(na_pop<12)+1)]

## impute mean

for (i in 1:ncol(AllFreq))
{
  AllFreq[which(is.na(AllFreq[,i])),i] <- median(AllFreq[-which(is.na(AllFreq[,i])),i], na.rm=TRUE)
}

## Filtering on MAF
freq_mean <- colMeans(AllFreq)
AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)]

## Ordering loci based on their scaffold
AllFreq <- AllFreq[,order(colnames(AllFreq))]

################################################################################

path_to_climate_data = "~/Schreibtisch/RDA-landscape-genomics-main./RDA-landscape-genomics-main/Data/ClimateNA/output_tiles/significant"


## Loading the climatic rasters
files <- list.files(path_to_climate_data, "\\.tif$", full.names = TRUE)
ras <- rast(files) # stack rasters
names(ras) <- tools::file_path_sans_ext(basename(files)) #name stack layers

## store rasters for each time point in substack
ras_6190 <- ras[[grep("6190", names(ras))]]
names(ras_6190) <- unlist(strsplit(names(ras_6190), split = "_6190"))

ras_2050 <- ras[[grep("2050_85", names(ras))]]
names(ras_2050) <- unlist(strsplit(names(ras_2050), split = "_2050_85"))

ras_2080 <- ras[[grep("2080_85", names(ras))]]
names(ras_2080) <- unlist(strsplit(names(ras_2080), split = "_2080_85"))

## check
names(ras_6190)
names(ras_2050)
names(ras_2080)

#################################################################################
## Dieser Block ist ein klassischer Trick, um sicherzustellen, dass alle Raster denselben NA-Muster haben, 
## also dass an jeder Zelle entweder für alle Layer Werte vorhanden sind oder NA, damit spätere Analysen sauber laufen.

remove_NAs_stack <- function(rast_stack){
  nom <- names(rast_stack)
  
  # 1. Summe über alle Layer
  test1 <- app(rast_stack, fun = sum, na.rm = FALSE)
  
  # 2. Maske: 1 für gültige Zellen, NA für fehlende
  test1[!is.na(test1)] <- 1
  
  # 3. Multiplikation Layerweise
  test2 <- rast_stack * test1
  
  # 4. Namen wiederherstellen
  names(test2) <- nom
  
  return(test2)
}

ras_6190 <- remove_NAs_stack(ras_6190)
ras_2050 <- remove_NAs_stack(ras_2050)
ras_2080 <- remove_NAs_stack(ras_2080)

# Climate Variables Plot
plot(ras_6190[["Eref"]],
     col = hcl.colors(100, "YlOrRd"),
     main = "DD18 (Westkanada)")

################################################################################

## Source populations coordinates
Coordinates <- read.table("./Data/PlSeedlots.csv", sep = ",", header = T, row.names = 1)
Coordinates <- Coordinates[match(row.names(AllFreq), Coordinates$id2, nomatch = 0),]
colnames(Coordinates) <- c("Population", "Latitude", "Longitude", "Elevation")


# Show Coordinates
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(data = Coordinates[, c("Latitude", "Longitude", "Elevation")],
             aes(x = Longitude, y = Latitude),
             color = "red",
             size = 2) +
  coord_sf() +
  theme_minimal()

# Zoom in
ggplot(world) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(data = Coordinates[, c("Latitude", "Longitude", "Elevation")],
             aes(x = Longitude, y = Latitude),
             color = "red", size = 2) +
  coord_sf(xlim = c(-140, -110),
           ylim = c(45, 65)) +
  theme_minimal()



## The coordinates were then used to extract climatic data for each pixel/site from the 27 climatic rasters.

# Koordinaten extrahieren: Longitude = Spalte 3, Latitude = Spalte 2 

Env <- data.frame(extract(ras_6190, Coordinates[, c("Longitude", "Latitude")]))
Env$ID <- NULL

# Populationsnamen als Zeilen setzen
rownames(Env) <- Coordinates$Population

## Standardization of the variables
Env <- scale(Env, center=TRUE, scale=TRUE) 

## Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')

## Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- c(Coordinates$Population)

################################################################################

## Load Phenotypic data

## Estimating mean trait value per population
traits <- aggregate(InfoInd[,19:24], by = list(InfoInd$ProvSeedlotCode), function(x) mean(x, na.rm = T))
traits <- as.data.frame(scale(traits[match(Coordinates$Population, traits$Group.1, nomatch = 0),-1]))
colnames(traits) <- c("Height","GthRate","ShootDryMass","Gth5pct","Gth95pct","ColdInjury")

################################################################################

##Population structure analysis

## Loading the intergenic SNPs dataset
Neutral <- read.table("./Data/Pine_AllNatural_GCandTotemIndivs_ControlSNPs_June8th2019.txt", header = T)
names_ind_neutral <- row.names(Neutral)

## Formatting genotypes
Neutral <- as.data.frame(apply(Neutral, 2, genotypes2geno))
row.names(Neutral) <- names_ind_neutral

## Sorting genetic data
Neutral <- Neutral[match(InfoInd$Internal_ID, row.names(Neutral), nomatch = 0),]

## Estimating allele frequencies for each source population
AllFreq_neutral <- aggregate(Neutral, by = list(InfoInd$ProvSeedlotCode), function(x) mean(x, na.rm = T)/2)

## Imputation of missing population frequencies by the median across the complete sampling
for (i in 2:ncol(AllFreq_neutral))
{
  AllFreq_neutral[which(is.na(AllFreq_neutral[,i])),i] <- median(AllFreq_neutral[-which(is.na(AllFreq_neutral[,i])),i], na.rm=TRUE)
}
AllFreq_neutral[,-1] <- AllFreq_neutral[,-1][,-which(is.na(colMeans(AllFreq_neutral[,-1])))]

## Running a PCA on neutral genetic markers
pca <- rda(AllFreq_neutral[,-1], scale=T) # PCA in vegan uses the rda() call without any predictors
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")

## Neutral population structure table
PCs <- scores(pca, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = AllFreq_neutral[,1], PCs)
colnames(PopStruct) <- c("Population", "PC1", "PC2", "PC3")

## Table gathering all variables
Variables <- data.frame(Coordinates, PopStruct[,-1], Env, traits)


# load boundaries of all countries as polygone objects in admin
# Parameters:
# scale = "medium" → mittlere Detailgenauigkeit der Ländergrenzen.
# returnclass = "sf" → die Daten werden als sf-Objekt (Simple Features) zurückgegeben, ein moderner Standard für Geodaten in R.

admin_sf <- ne_countries(scale = "medium", returnclass = "sf")
# In terra-Objekt umwandeln
admin <- vect(admin_sf)

# Lade Artenverbreitung
range <- vect("./Data/pinucon/pinucon_o.shp")

################################################################################
###Variance partitioning: disentangling the drivers of genetic variation 

## Full model
pRDAfull <- rda(AllFreq ~ PC1 + PC2 + PC3 + Longitude + Latitude + MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS,  Variables)
RsquareAdj(pRDAfull)

anova(pRDAfull)

## Pure climate model
pRDAclim <- rda(AllFreq ~ MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + Condition(Longitude + Latitude + PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAclim)

anova(pRDAclim)

## Pure neutral population structure model  
pRDAstruct <- rda(AllFreq ~ PC1 + PC2 + PC3 + Condition(Longitude + Latitude + MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS),  Variables)
RsquareAdj(pRDAstruct)

anova(pRDAstruct)

##Pure geography model 
pRDAgeog <- rda(AllFreq ~ Longitude + Latitude + Condition(MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAgeog)

anova(pRDAgeog)


library(corrplot)
cor_mat <- cor(Variables[,c("MAR", "EMT", "MWMT", "CMD", "Tave_wt", "DD_18", "MAP", "Eref", "PAS" )], use = "pairwise.complete.obs")
corrplot(cor_mat,
         method = "color",
         type="upper",
         col = colorRampPalette(c("red", "white", "blue"))(200),
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         number.cex = 0.7)

Variables2 <- as.data.frame(Variables)

Variables2["Elevation"] <- NULL
library(dplyr)
rownames(Variables2) <- Variables2$Population
Variables2 <- Variables2 %>% select("PC1", "PC2", "PC3", "Longitude", "Latitude", "MAR", "EMT", "MWMT", "CMD", "Tave_wt", "DD_18", "MAP", "Eref", "PAS")

ccr<-cor(Variables2)
corrplot::corrplot(ccr)
diag(ccr)<-0
thr<-range(abs(ccr))[2]
while(thr>0.6){
  ccr<-cor(ev1)
  diag(ccr)<-0
  thr<-range(abs(ccr))[2]
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}

lyrs<-colnames(ev1)
corrplot::corrplot(ccr)
Variables<-Variables[,lyrs]
#bio8, bio9, bio18

evars<-lyrs
eq<-paste0(evars,collapse = "+")
RDAfull<-rda(as.formula(paste0("AllF ~",eq)),Variables)
(vf<-vif.cca(RDAfull))

#keep only the vars with VIF < 5 # 10 für ein relaxteren threshold
while(max(vf)>5){
  evars<-names(vf)
  evars<-evars[-which.max(vf)]
  eq<-paste0(evars,collapse = "+")
  RDAfull<-rda(as.formula(paste0("AllF ~",eq)),Variables)
  (vf<-vif.cca(RDAfull))
}

################################################################################

rownames(Variables2) <- Variables2$Population
Variables2$Population <- NULL
Variables2 <- as.data.frame(scale(Variables2))


## Genotype-Environment Associations: identifying loci under selection

RDA_env <- rda(AllFreq ~ MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + Condition(PC1 + PC2 + PC3),  Variables2)

screeplot(RDA_env, main="Eigenvalues of constrained axes")

## Function rdadapt
source("./src/rdadapt.R")

## Running the function with K = 2
rdadapt_env<-rdadapt(RDA_env, 2)


## P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rdadapt_env$p.values)

## Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

## Top hit outlier per contig
outliers <- outliers[order(outliers$contig, outliers$p.value),]

## List of outlier names
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

## Manhattan plot
Outliers <- rep("Neutral", length(colnames(AllFreq)))
Outliers[colnames(AllFreq)%in%outliers$Loci] <- "All outliers"
Outliers[colnames(AllFreq)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(AllFreq)), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outliers)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))


## Not accounting for population structure

## Running a simple RDA model
RDA_env_unconstrained <- rda(AllFreq ~ MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS,  Variables)

## Running the rdadapt function
rdadapt_env_unconstrained <- rdadapt(RDA_env_unconstrained, 2)

## Setting the p-value threshold 
thres_env <- 0.01/length(rdadapt_env_unconstrained$p.values)

## Identifying the outliers for the simple RDA
outliers_unconstrained <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env_unconstrained$p.values<thres_env)], p.value = rdadapt_env_unconstrained$p.values[which(rdadapt_env_unconstrained$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env_unconstrained$p.values<thres_env)], split = "_"), function(x) x[1])))
outliers_unconstrained <- outliers_unconstrained[order(outliers_unconstrained$contig, outliers_unconstrained$p.value),]
outliers_rdadapt_env_unconstrained <- as.character(outliers_unconstrained$Loci[!duplicated(outliers_unconstrained$contig)])

## For all the outliers
list_outliers_RDA_all <- list(RDA_constrained = as.character(outliers$Loci), RDA_unconstrained = as.character(outliers_unconstrained$Loci))

## Only for the top hit locus per contig
list_outliers_RDA_top <- list(RDA_constrained = outliers_rdadapt_env, RDA_unconstrained = outliers_rdadapt_env_unconstrained)

# Outliers found by both methods
common_outliers_RDA_top <- Reduce(intersect, list_outliers_RDA_top)



################################################################################

## Adaptive landscape: projecting adaptive gradient(s) across space 

## Adaptively enriched RDA
RDA_outliers <- rda(AllFreq[,common_outliers_RDA_top] ~ MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS,  Variables)

## RDA biplot
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))


## Function to predict the adaptive index across the landscape
source("./src/adaptive_index_terra.R")

## Running the function for all the climatic pixels of lodgepole pine distribution range
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = ras_6190, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)


