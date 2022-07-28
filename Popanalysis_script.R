## workspace set up

library(vcfR)
library(adegenet)
library(poppr)

library(ggplot2)
library(ggrepel)

library(fields)
library(mapplots)
library(LEA)

library(SNPRelate)
library(dartR)
library(reshape2)
library(ggpubr)

library(svglite)

## Load some LEA functions
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

## Keep the default margins before we change them
default_par_mar = par("mar")
default_par_oma = par("oma")


########################################################################################
## ========   Basic required input files  ========
# 1. vcf file containing filtered SNP
# 2. Population information file containing in csv format

# e.g.
# name        location         lat      long      Altitude grp_plains grp_long       patch    grp_river region
# 1 AndyLT-1-1 Andy Lane Trail 37.45592 -80.01412      453 Moutains Central  Central AndyLT.1 Catawba Creek center
# 2 AndyLT-1-2 Andy Lane Trail 37.45592 -80.01412      453 Moutains Central  Central AndyLT.1 Catawba Creek center
# 3 AndyLT-1-3 Andy Lane Trail 37.45592 -80.01412      453 Moutains Central  Central AndyLT.1 Catawba Creek center
# 4 AndyLT-1-4 Andy Lane Trail 37.45592 -80.01412      453 Moutains Central  Central AndyLT.1 Catawba Creek center
# 5 AndyLT-1-5 Andy Lane Trail 37.45592 -80.01412      453 Moutains Central  Central AndyLT.1 Catawba Creek center
# 6 AndyLT-2-1 Andy Lane Trail 37.45561 -80.01260      441 Moutains Central  Central AndyLT.2 Catawba Creek center


########################################################################################
## Set up the working directory
setwd("/home/mydir/")

#########################################################################################
## LOAD THE SAMPLES FOR THE POPULATION STRUCTURE ANALYSIS
#########################################################################################

## 1.1- Load the samples
## Package (vcfR)

AstriPanel = read.vcfR("filteredSNP.vcf")


# convert the vcf file to a genlight object 
GL_AstriPanel = vcfR2genlight(AstriPanel)
## Provide the ploidy of the samples
ploidy(GL_AstriPanel) = 2

# 1.2- Load the population information
pop.data = read.csv("pop.data.csv", header = TRUE )
# confirm that sample names in popdata and vcf match
all(colnames(AstriPanel@gt)[-1] == pop.data$name)
# assign population for analysis
pop(GL_AstriPanel) = pop.data$river_basin

# ============ Clonal correction =========
## Package (poppr)
## If clones are identified they can be removed to avoid potential bias
## To separate by strata, add a column with strata separated by "_"
# e.g. Zone_River_Region_Patch

# Convert the strata column to a dataframe
popstrata <- as.data.frame(pop.data$strata) 
Astri_sc <- as.snpclone(GL_AstriPanel, parallel = F, ) # Convert genlight object to SNPclone format and assign strata to the SNPclone object
Astri_sc
strata(Astri_sc) <- popstrata
## now its possible to subset the data by strata
splitStrata(Astri_sc) <- ~Patch/River/Zone

## the sampling by strata
Astristrata <- strata(Astri_sc) %>%
  group_by(Patch, River, Zone) %>%
  summarize(Count = n())  

nameStrata(Astri_sc) # The order of our variables
Astristrata$Count    # The variable used for the block size

treemap(dtf = Astristrata, index = nameStrata(Astri_sc), vSize = "Count",
        type = "categorical", vColor = "Patch", title = "A. triloba population stratifications",
        align.labels = label_position, fontsize.labels = label_size)

itreemap(dtf = Astristrata, index = nameStrata(sc), vSize = "Count",
         type = "categorical", vColor = "Region")

Astricc <- clonecorrect(Astri_sc, strata = ~Patch/Zone , keep = 1:2)

Astricc
Astri_sc

setPop(Astri_sc) <- ~Patch/Zone

cc <- locus_table(x=Astricc, info=F, )
scc <- locus_table(GL_Astri_sc, info=F)
#####################################################################################
# Package (LEA)
## 1.3- Create the Structure file
# FASTRUCTURE preforms a Bayesian analysis assuming Hardy–Weinberg equilibrium and 
# linkage equilibrium between loci within population.This is not suitable for clonal population 
# but may still provide a good indication of population structure

## Export the GL object to a Structure file for further analysis (may take some time)
gl2faststructure(GL_AstriPanel, outfile = "AstriDiversityPanel.Variants.fs.str", outpath = ".")

#########################################################################################
## STRUCTURE ANALYSIS WITH FASTSTRUCTURE
#########################################################################################
struct2geno=struct2geno("./AstriDiversityPanel.Variants.fs.str" , FORMAT = 2,  extra.row = 0, extra.col = 0)
## 2.2- Run the Structure analysis with the desired number of clusters (K) e.g., from 1 to 30
obj_str = snmf("genotype.geno", K = 1:30, ploidy = 2, entropy = T, CPU = 5, project = "new")

## 2.3- Check for the best K
plot(obj_str, col = "blue4", cex = 1.8, pch = 19 ,main="<Plot title>", font.main = 1)
## add abline at notable mimimas
abline(v=c(24,29), col="red")

## the cross entropy score for a given K value
cross.entropy(obj_str, K = c(1)) 
cross.entropy(obj_str, K = c(10))
cross.entropy(obj_str, K = c(15))
cross.entropy(obj_str, K = c(24))
cross.entropy(obj_str, K = c(29))

#############  Admixture ##################################
## LEA admixture is run from the FASTRUCTURE object, cross validation using ADMIXTURE software can
# be run using the guide at https://speciationgenomics.github.io/ADMIXTURE/

par(mar=c(6,3,1,1))
qmatrix = Q(obj_str, K = 4)
row.names(qmatrix)=pop.data$name
#qmatrix <- qmatrix[order(qmatrix[,2],qmatrix[,1] ),]
barplot(t(as.matrix(qmatrix)), las =2 , horiz=F, col = Duet , las=2, cex.names = 0.5,
        main = "<Plot title>")

#########################################################################################
## Neighbor-joining tree
#########################################################################################
## 3.1- Generate an NJ tree
library(ape)
#package(poppr) ## generate tree
#package (ape) ## plot tree

pop(GL_AstriPanel) = pop.data$river_basin

Astri_tree = aboot(GL_AstriPanel , tree = "nj", distance = bitwise.dist(), sample = 100, showtree = TRUE, cutoff = 50, quiet = T)

## if more colours are needed we can use rainbow
rainbow(n = nPop(GL_AstriPanel))
cols = rainbow(n = nPop(GL_AstriPanel))

par(oma=c(3,2,1,1), mar=c(0,0,0,0))
par(bg="white")
plot.phylo(Astri_tree, cex = 0.8, font = 4, adj = 0, tip.color = cols[pop(GL_AstriPanel)] )
nodelabels(Astri_tree_tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend('bottomright', legend = levels(pop(GL_AstriPanel)), fill = cols, border = FALSE, bty = "n", cex = 0.7)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

# Write out the tree in netwick format to use with other software like Figtree or splitstree
ape::write.tree(Astri_tree, file = "mydir/treename.newick")


#########################################################################################
## PCA ANALYSIS
#########################################################################################
## package (adegenet)
## 4.1- Run PCA analysis
pop(GL_AstriPanel) = pop.data$river_basin

Astri.pca = glPca(GL_AstriPanel, nf=3)
str(Astri.pca)
## 4.2-Plot the PCA

## Set up new margins
par(oma=c(3,2,1,1), mar=c(1,1,1,1), font=1 )
## Plot the barplot for the eigenvalues
barplot(100*Astri.pca$eig/sum(Astri.pca$eig), col = rep(c("red","grey"),c(3,1000)), main="<Eigenvalues plot title>", xlab="Eigenvalues", ylab="Percent of variance explained")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)


Astri.pca_scores = as.data.frame(Astri.pca$scores)
Astri.pca_scores$pop = pop(GL_AstriPanel)

options(ggrepel.max.overlaps = Inf)
# choose which PCs to compare
p <- ggplot(Astri.pca_scores, aes(x=PC1, y=PC2, colour = pop.data$river_basin))
p <- p + geom_point(aes(shape=pop.data$grp_river), size = 11, alpha=0.8 ) +
          scale_shape_manual(values=c(15, 16, 17, 18, 21, 20 , 12, 8, 13))
## To show the labels we will use the geom_text_repel function
p <- p + geom_text_repel(aes(label= pop.data$name), show.legend = FALSE)
p <- p + scale_color_manual(values = <mycolorpalette>)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + theme_classic()
## to find the PCA percentage use the eignvalue - eigenvalue of first divided by the sum of all eigenvalues multiplied by 100
# sum(Astri.pca$eig[1]/sum(Astri.pca$eig))*100
#p <- p + xlab("PC 1 (5.9%)")
p <- p + xlab("PC 2 (5.0%)")
p <- p + ylab("PC 3 (4.5%)")
#p + geom_label_repel(aes(label=rownames(Astri.pca_scores)), show.legend = FALSE)
p <- p + ggtitle("<Plot title>")
p <- p + theme(
  plot.title = element_text(color="Black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12, face="bold"),
  axis.title.y = element_text(color="black", size=12, face="bold"),
  legend.title = element_blank(),
  legend.position="none")
p


#########################################################################################
## DAPC ANALYSIS
#########################################################################################
## Tutorial https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

# Package (reshape2)
# Package (adegenet)

par(mar=c(2,2,3,1))
pop(GL_AstriPanel) = pop.data$river_basin
dapc_optimize <- adegenet::dapc(GL_AstriPanel, n.da = 100)
## before checking for clusters optim.a.score to see how many pcas to include
temp <- optim.a.score(dapc_optimize)
## use find.clusters to identify clusters with lowest associate Bayesian Information Criterion (BIC)

names(grp)
head(grp$Kstat, 30) ## BIC values of k
grp$stat           ## selected number of clusters and associated values
head(grp$grp)      ## Group membership
grp$size           ## Group sizes

table(pop(GL_AstriPanel), grp$grp)

## to check the representation of each group with the suggested number of pcs (x) use:
temp <- summary(dapc(GL_AstriPanel, n.da=100, n.pca=6))$assign.per.pop*150
par(mar=c(3,12,1,1))
barplot(temp, xlab= "% of reassignment to group",
        horiz=TRUE, las=1, main = "<Plot title>")

## Run the DAPC analysis with n components and groups
pop(GL_AstriPanel) = pop.data$river_basin
Astri_DAPC <- adegenet::dapc(GL_AstriPanel, n.pca = 6, n.da = 100)

## Plot components
par(oma=c(1,1,1,1), mar=c(5,1,1,1), font=2 )


scatter(Astri_DAPC, col = lmute9, cex = 2, legend = FALSE,
        clabel = TRUE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "bottomleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 0.5)

par(mar=c(5,1,1,1))
compoplot(Astri_DAPC, col = lmute9, space = 0.09, posi = 'top', show.lab = TRUE, legend=TRUE,
          cex.names=0.5, cleg=0.6, size.leg=8 )



## Generate the barplot with populations assigned to each sector
dapc.results <- as.data.frame(Astri_DAPC$posterior)
dapc.results$pop <- pop(GL_AstriPanel)
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Individual","River Basin","Posterior membership probability")

p <- ggplot(dapc.results, aes(x=Individual, y=`Posterior membership probability`, fill=`River Basin`))
p <- p + geom_bar(stat='identity')
p <- p + scale_fill_manual(values = lmute9)
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p <- p + reorder(dapc.results$`Posterior membership probability`)
p




## The DAPC can be used also to infer the number of groups and then plot them
## in a similar fashion than the Structure. To do it:

## Check the number of groups and their variation to find the optimal K

default_par_mar = par("mar")
default_par_oma = par("oma")

## Use as maximum K=30 and create a matrix with the Kstat results of find clusters 
maxK = 30
myMat = matrix(nrow=5, ncol=maxK)
colnames(myMat) = 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp = find.clusters(GL_AstriPanel, n.pca = 8, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

## Generate a data.frame for the plot
my_df = melt(myMat)
colnames(my_df)[1:3] = c("Group", "K", "BIC")
my_df$K = as.factor(my_df$K)
head(my_df)

## Plot the BIC values associated with find.clusters
p1 = ggplot(my_df, aes(x = K, y = BIC))
p1 = p1 + geom_boxplot()
p1 = p1 + theme_bw()
p1 = p1 + xlab("Number of groups (K)")
p1


## Create the accession group assignments based in the Linear Discriminants (LD) 

## We will assay 4 groups
pop(GL_AstriPanel) = pop.data$river_basin

my_k = c(4,5)

grp_l = vector(mode = "list", length = length(my_k))
dapc_l = vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] = find.clusters(GL_AstriPanel, n.pca = 6, n.clust = my_k[i])
  dapc_l[[i]] = dapc(GL_AstriPanel, pop = grp_l[[i]]$grp, n.pca = 6, n.da = my_k[i])
}

my_df = as.data.frame(dapc_l[[1 ]]$ind.coord)
my_df$Group = dapc_l[[ 1 ]]$grp
#my_df = as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
#my_df$Group = dapc_l[[ length(dapc_l) ]]$grp
head(my_df)



p2 = ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 = p2 + geom_point(size = 8, shape = 16)
p2 = p2 + theme_bw()
p2 = p2 + scale_color_manual(values=c(lmute9))
p2 = p2 + scale_fill_manual(values=c(paste(mega_col, "66", sep = "")))
p2

##  Create the accession group assignments based in the posterior probability
tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop(GL_AstriPanel)
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop(GL_AstriPanel)
  
  my_df <- rbind(my_df, tmp)
}


grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "free",
                      labeller = labeller(K = grp.labs) )
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + xlab(element_blank ())
p3 <- p3 + theme(legend.position='none')
p3 <- p3 + scale_fill_manual(values=c(lmute9))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3


###############################################################################################################
## Minimum spanning networks (MSN) poppR
## MSN clusters multilocus genotypes (MLG) by genetic distances between them
# https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html

# Package (igraph)
# Package (poppr)

pop(GL_AstriPanel) = pop.data$river_basin

rubi.dist <- bitwise.dist(GL_AstriPanel)
rubi.msn <- poppr.msn(GL_AstriPanel, rubi.dist, showplot = FALSE, include.ties = T)

node.size <- rep(2, times = nInd(GL_AstriPanel))
names(node.size) <- indNames(GL_AstriPanel)
vertex.attributes(rubi.msn$graph)$size <- node.size


set.seed(99)
plot_poppr_msn(GL_AstriPanel, rubi.msn ,  palette = mega_col , gadj = 90, 
               wscale = TRUE , layfun = layout_with_graphopt)


#########################################################################################
## BASIC STATS FOR GROUPS (based in the groups produced by DAPC)
#########################################################################################

########################################################################################
## BASIC STATS USING PDartR and POPGENOME
#########################################################################################
## Once the different individuals have assigned to different populations, last step will be to redo the population 
## assignment and calculate some population genetic parameters.
pop(GL_AstriPanel) = pop.data$river_basin
dartR.fst <- dartR::gl.fst.pop(GL_AstriPanel, nboots = 1000, percent = 95, nclusters = 3)
head(dartR.fst)
# $Fsts
# center       west     south east
# center         NA         NA        NA   NA
# west   0.04539081         NA        NA   NA
# south  0.10982481 0.12210849        NA   NA
# east   0.04460933 0.07036082 0.1423037   NA

He <- dartR::gl.He(gl = GL_AstriPanel)

my_dfK = my_df[my_df$K == 4,]
Kgroups = my_dfK[my_dfK$Posterior > 0.8,]
Kgroups = my_dfK
length(Kgroups$Isolate)

## Asign group from the my_df object or from the pop.data datatable
pop(GL_AstriPanel) = pop.data$river_basin

## Generate the stats
GL_Astri_BStats = gl.basic.stats(GL_AstriPanel)
GL_Astri_BStats
He_dartR <- GL_Astri_BStats$Ho

write.table(He_dartR, "He_dartr.txt", sep = "\t")
He_df = read.csv("./He_dartr.csv", header = TRUE )


ggplot(He_df, aes(x=River.Basin, y=He)) + 
  geom_boxplot()+
  #geom_jitter(size = 2, alpha=0.4, position = position_jitter(height=0.02)) + 
  stat_smooth(method = "loess", colour = "blue", size = 1.5)

## Calculate the mean Fst per group [,1], [,2]....
mean(GL_Astri_BStats$Ho[,1])
mean(GL_Astri_BStats$Ho[,2])
mean(GL_Astri_BStats$Ho[,3])
mean(GL_Astri_BStats$Ho[,4])

mean(GL_Astri_BStats$Hs[,1])
mean(GL_Astri_BStats$Hs[,2])
mean(GL_Astri_BStats$Hs[,3])
mean(GL_Astri_BStats$Hs[,4])
### It is necessary to break the VCF into sequences (for readData). To do it it is necessary to use bcftools
### in Linux:
### 1- Compress the file: "bgzip -c myfile.vcf > myfile.vcf.gz"
### 2- Index the compressed file "tabix -p vcf myfile.vcf.gz"
### 3- Create a script to divide the VCF file per sequences
###    
### 4- Executate the script   

library(PopGenome)
## 6.1- Load the data
AstriPanel_PG = readData("untitled folder/", format = "VCF")


## Check the total number of SNPs
AstriPanel_PG_sumdata = as.data.frame(get.sum.data(AstriPanel_PG))
sum(AstriPanel_PG_sumdata$n.biallelic.sites)
## It gives 20,000 SNPs less and I am not sure why. Nevertheless, the stats should be similar

## Add the populations
Kgroup1list = pop.data$name[pop.data$river_basin == "New River"]
Kgroup2list = pop.data$name[pop.data$river_basin == "Holsten River"]
Kgroup3list = pop.data$name[pop.data$river_basin == "James River"]
Kgroup4list = pop.data$name[pop.data$river_basin == "York River"]
AstriPanel_PG = set.populations(AstriPanel_PG, list( Kgroup1list, Kgroup2list, Kgroup3list, Kgroup4list), 
                                diploid = TRUE)


#Kgroup1list = pop.data$name[pop.data$grp_jr == "VA"]
#Kgroup2list = pop.data$name[pop.data$grp_jr == "JR"]
#AstriPanel_PG = set.populations(AstriPanel_PG, list( Kgroup1list, Kgroup2list), 
#                                diploid = TRUE)

## 6.2- Estimate the different parameters
PG_Astri_concatenated = concatenate.regions(AstriPanel_PG)
PG_Astri_concatenated = neutrality.stats(PG_Astri_concatenated, FAST=TRUE)
get.neutrality(PG_Astri_concatenated)
PG_Astri_concatenated = diversity.stats(PG_Astri_concatenated)
get.diversity(PG_Astri_concatenated)
PG_Astri_concatenated = F_ST.stats(PG_Astri_concatenated)
get.F_ST(PG_Astri_concatenated)
PG_Astri_concatenated = detail.stats(PG_Astri_concatenated, site.spectrum=TRUE, site.FST=TRUE) 
PG_Astri_results = get.detail(PG_Astri_concatenated) 

## Get the segregating sites
PG_Astri_concatenated@n.segregating.sites

## Get nucleotide diversity, as Pi, for each population
PG_Astri_concatenated@Pi

## Get the Tajima D for each of the populations
PG_Astri_concatenated@Tajima.D

## Get the Waterson Theta for each of the populations
PG_Astri_concatenated@theta_Watterson

## Rm11DP3Q1000
## $overall
## Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
## 0.3364  0.2760  0.2844  0.0084  0.2927  0.0168  0.0295  0.0573 -0.2189  0.0232 



# --------------------------------------#
#  mean(GL_OleaPanelCp_BStats$Ho[,1])   #
#  [1] 0.3626456                        #
#  mean(GL_OleaPanelCp_BStats$Ho[,2])   #
#  [1] 0.3100806                        #
# --------------------------------------#
# --------------------------------------#
#  mean(GL_OleaPanelCp_BStats$Ho[,1])   #
#  [1] 0.2942231                        #
#  mean(GL_OleaPanelCp_BStats$Ho[,2])   #
#  [1] 0.257731                         #
# --------------------------------------#



## Now it will extract some basic stats per population group
## +-------+-------------+----+-----------+-------------+
## | Group |  Group_Name | N  |    Ho     |    Hs       |
## +-------+-------------+----+-----------+-------------+
## |   1   |  OldDomest  | 17 | 0.09302979 | NaN        |
## |   2   |  NewDomest  | 14 | 0.0871276  | 0.07737128 |
## |   3   |  WildTypes  | 27 | 0.07663825 | 0.07217912 |
## +-------+-------------+----+-----------+-------------+

## +-------+-------------+----+-----------------+-----------------+-----------------+-----------+
## | Group |  Group_Name | N  | Segr. Sites (S) | Nucl. div. (Pi) | Theta Watterson | TajimaD   | 
## +-------+-------------+----+-----------------+-----------------+-----------------+-----------+
## |   1   |  OldDomest  | 17 | 29636           | 7557.859        | 7248.095        | 0.1645776 |
## |   2   |  NewDomest  | 14 | 18813           | 5328.153        | 4834.436        | 0.4045994 | 
## |   3   |  WildTypes  | 27 | 27404           | 5496.317        | 6013.722        | -0.312761 |
## +-------+-------------+----+-----------------+-----------------+-----------------+-----------+
##
## Intertretation:
## 1- Segregating Sites (S): Estimation of the mutation rate assuming no selection (mostly used to esti,ate other parameters).
## 2- Nucleotide Diversity (Pi): Degree of polymorphism in a population. It is measure of the genetic variation, so the 
## higher the Pi, the more genetic variation there is contained in a population. It is similar to expected heterozygosity.
## 3- Theta Watterson: It is also an estimator of the genetic diversity. 
## 4- Tajima D: It is an estimator of selection (https://en.wikipedia.org/wiki/Tajima%27s_D#Interpreting_Tajima's_D)

## So, according these values, the group 1, in which we have the monumental Vouves olive tree (5,000 year old) as well 
## as Santander (1,200 years old, variety Farga) and some old Italian varieties presents the higher estimators for 
## nucleotide diversity. Then the other groups present similar values, with slightly higher values in the group of the
## wild types. Nevertheless this group present a negative Tajima D representing a recent expansion after a genetic
## bottleneck.

###################################
## https://www-users.york.ac.uk/~dj757/popgenomics/workshop5.html#plotting_tajima%E2%80%99s_d

pi.all <- read.table("Pi_minDP5_prunned.windowed.pi",header=T)
head(pi.all)
str(pi.all)
summary(pi.all)
taj.all <- read.table("TajimaD_minDP5_prunned.Tajima.D",header=T)
nrow(taj.all)
ncol(taj.all)
str(taj.all)





###################################
## Heterozygosity
###################################
angsd=read.csv("het.csv", header = TRUE)

ggplot(data = angsd, mapping = aes(x = name, y = He , color = river_basin)) +
  geom_point(aes(shape = grp_river) , size = 6) +
  scale_shape_manual(values=c(15, 16, 17, 18, 21, 20 , 12, 8, 13))+
  scale_color_manual(values = lmute9) +
  labs(x = "Individual sample name", y = "ANGSD estimated Heterozygosity") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12, face="bold"))
theme_bw()




ggplot(data = pi, mapping = aes(x = pi$Patch, y = pi$Nucleotide_Diversity)) +
  geom_point( size = 6) +
  #scale_shape_manual(values=c(15, 16, 17, 18, 21, 20 , 12, 8, 13))+
  scale_color_manual(values = lmute9) +
  labs(x = "Individual sample name", y = "ANGSD estimated Heterozygosity") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12, face="bold"))
theme_bw()




hist(taj.all$TajimaD,br=20)

library(ggstatsplot)
#Patil, I. (2021). Visualizations with statistical details: The
#'ggstatsplot' approach. Journal of Open Source Software, 6(61), 3167,
#doi:10.21105/joss.03167
library(palmerpenguins)
data("penguins", package = "palmerpenguins")

penguins <- drop_na(penguins)
plt <- ggbetweenstats(
  data = penguins,
  x = species,
  y = bill_length_mm
)


pi=read_csv("Pi_diversity.csv")
pi2=as.data.frame(pi)

#https://r-graph-gallery.com/web-violinplot-with-ggstatsplot.html
plt <- ggbetweenstats(
  data = pi,
  x = patch,
  y = pi)
plt


#https://paulvanderlaken.com/2018/08/26/ggstatsplot-creating-graphics-including-statistical-details/


ggstatsplot::ggscatterstats(
  data = pi, 
  x = pi, 
  y = mlg,
  title = "Observed heterozygosity in Virginia populations by longitude ",
  messages = FALSE
)


# ===================== Isolation by distance =============


library(adegenet)
library(genepop)
########################################################################
## Load file
AstriPanel = read.vcfR("Astri.minDP5.LDpruned_renamedSorted.vcf")
GL_AstriPanel = vcfR2genlight(AstriPanel)
ploidy(GL_AstriPanel) = 2

## Pop data
pop.data = read.csv("pop.data.csv", header = TRUE )
all(colnames(AstriPanel@gt)[-1] == pop.data$name)
pop(GL_AstriPanel) = pop.data$name

## Coordinates
coor <- read.table("latlong.txt" , header =T , sep = "\t")
colnames(coor)<-c("x","y")
GL_AstriPanel@other$latlong <- coor



Astri.gen <- gl2gi(GL_AstriPanel)
Astri.pop <- genind2genpop(Astri.gen)


Dgen <- dist.genpop(Astri.pop,method=2)
Dgeo <- dist(Astri.pop$other$latlong)
ibd <- mantel.randtest(Dgen,Dgeo)
ibd
par(mar=c(3,2,1,1))
plot(ibd)
# Cline or distant patches?
plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo), col="red",lty=2)

library(MASS)
dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")

