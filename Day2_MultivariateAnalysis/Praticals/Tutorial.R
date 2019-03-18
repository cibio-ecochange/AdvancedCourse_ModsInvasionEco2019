##############################################
##MULTIVARIATE ANALYSIS AND NON-NATIVE TREES## 
######PRACTICALS - PART I ####################
##############################################

################
##INTRODUCTION##
################

## This document contains a set of exercises aiming to test some multivariate analysis used in ecology
## These exercises were adapted from Vegan to meet the goals of the course
## To conduct the exercises, "vegan" (URL https://cran.r-project.org, https://github.com/vegandevs/vegan) and "MASS" packages are required
## The set of exercises here proposed, consist of four different parts: unconstrained ordination analysis, constrained ordination modelling, tests of dissimilarity, cluster analysis.

###################
##BEFORE STARTING##
###################
# Set your working directory (mind that you should add the files for this tutorial in the same page as you decide as working directory)
# If you're using R studio: Session -> Set Working Directory -> Choose Directory

setwd("./Day2_MultivariateAnalysis/Praticals") # or you can add this code


# Download package "vegan" (it might need also packages "permute" and "lattice")
# If you're using R studio: Tools -> Install Packages -> search for packages
install.packages("vegan") 
install.packages("MASS")

# Then, recall the packages to start working
library(vegan)
library(MASS)

############################################
##PART I - EXPLORATORY ORDINATION ANALYSIS##
############################################

########################################
##Example 1: Unconstrained ordination ##
########################################
# In this exercice, species data is not constrained by environmental data; 
# It's usefull as an exploratory tool for an indirect analysis, when: 
# (1) you have tree species data only
# (2) you're not sure if your tree data relates with an environmental dataset

trees<-read.table("trees.txt", sep=",", header=TRUE) # calling the dataset containing the abundance (%) of several tree species
View(trees) # view the dataset
str(trees) # check the properties of the dataset

trees.dis <- vegdist(trees) # Calculate distance
help(vegdist) #To check for other dissimilarity metrics. The default is Bray-Curtis

#dist<-as.data.frame(as.matrix(trees.dis)) #To save the distance matrix convert first to a data frame
#write.csv(dist, file = "treesdist.csv") #Then, save it

trees.nmds <- isoMDS(trees.dis) # perform the NMDS
help("isoMDS") #To check the details of NMDS calculation


#The results of isoMDS is a list of points for the configuration and the stress. 
#Stress S is a goodness of fit measure, as a function of and non-linear transformations of observed dissimilarities and ordination distances

stressplot(trees.nmds, trees.dis) # Shepard plot showing ordination distances against community dissimilarities

ordiplot(trees.nmds, type = "t") # Visualising the results of NMDS, considering dissimilarities among sites only
scores(trees.nmds) # Check score values for the different sites

trees.mds <- metaMDS(trees) # Calculate NMDS based on dissimiliarities among species
trees.mds # Check details on the calculation for NMDS
plot(trees.mds, type = "t") # Visualise the NMDS results
scores(trees.mds) # get species scores on NMDS

#Fit environmental data in NMDS
soil<-read.table("soil.txt", sep=",", header=TRUE) # calling the dataset containing several soil measures
View(soil) #Data set containg environmental data
efnmds <- envfit(trees.mds, soil, permu = 999) #Fitting data under a permutation test with 999 iteractions
efnmds #Check the direction of axes, correlation coefficient, and significance

plot(trees.mds, display = "sites") #NMDS plot
plot(efnmds, p.max = 0.05) #NMDS plot with fitted environmental data (p>0.05)

#You may now check the variation of environmental data in the ordination space
tmp <- with(soil, ordisurf(trees.mds, Al, add = TRUE))
with(soil, ordisurf(trees.mds, pH, add = TRUE, col = "green4"))

#You can try re-analysing this data, based on different dissimilarity measures, using: 
help(vegdist)

#Check which of the dissimilarity indices best separate communities along gradients using rank correlation 
rankindex(scale(trees), soil, c("euc","man","bray","jac","kul"))

#It is often useful to transform or standardize data, specially if there is a difference between smallest non-zero abundance and largest abundance
#Data can be standardised sites to equal sum of squares using a vector norm, or use the Hellinger distance (based on square roots of sites standardized to unit total)
dis1 <- vegdist(decostand(trees, "norm"), "bray")
dis2 <- vegdist(decostand(trees, "hell"), "bray")

#Function vegdist in vegan contains distance measures to be used as input.
##You cand define your own distance metric:
##help(designdist)

##You may ask how much did you win by using metaMDS (different start solutions) against isoMDS?
vare.mds0 <- isoMDS(dis1, trace = 0)
pro <- procrustes(trees.mds, trees.nmds)
pro
plot(pro) # Visualise differences between the different points
plot(pro, kind = 2) #Visualise the sum of squared arrows in the Procrustes plot


#If we only want to consider certain types of dissimilarities and make a linear mapping: eigenvalue methods
#Depending on your data, you may need to use different ordination methods (linear or unimodal):
  #PCA (linear method, based on Euclidean distances), CA (unimodal method, based on chi-quares), MDS (or principal coordinates, linear method, based on either of the two)
#To best decide on a linear or unimodal (weighted method), you may check first the lenght of the most explanatory gradient in your data (or the eigenvalue of the first ordination axis)
#To check this value, you should first conduct a DCA

trees.dca <- decorana(trees)
trees.dca

#If the lenghts of gradient of the first axis is higher than 4 SD, indicate high turnover levels and is more suitable for a weighted/unimodel distribution. 
#Values less than 3SD indicate a linear distribution, more suitable for linear analysis 
#Values betweeen 3-4SD may be suitable for both linear and weighted/unimodal analyses

#PCA
trees.pca <- rda(trees) #Run PCA with tree dataset
trees.pca #Check inertia = variance explained by the 23 gradients used; Inertia is the sum of all species variances explained by the gradients (eigenvalues)
plot(trees.pca) #See original plot of the PCA
scores(trees.pca)

biplot(trees.pca, scaling = -1) #PCA biplot with arrows for species and dots for sites: we can see that even with the scaling, the plot seems crowded

trees.pca2 <- rda(trees, scale = TRUE) #We may prefer to standarde all species to unit variance, or using correlation coefficients instead of covariances to give a more balanced ordination
trees.pca2 #Now, note that inertia is correlation (and the correlation with itself is one, providing the 44 species as the total inertia)
plot(trees.pca2, scaling = 3)

#Fit environmental data in PCA
efPCA <- envfit(trees.pca, soil, permutations = 999)
efPCA
plot(trees.pca2, display = "sites", type = "p")
plot(efPCA)

##CA
trees.ca <- cca(trees)
trees.ca #Since CA is based on Chi-squared distance, the inertia is the Chi-squared statistic of a data matrix standardized to unit total

chisq.test(trees/sum(trees)) #P-values which not be noticed here, but chech that the reported X-squared, similar to the CA inertia
plot(trees.ca, scaling = 1) #Scaling=1 to display the site scores as weighted averages of species scores

#Fit environmental data in CA (continuous data)
efCA <- envfit(trees.ca, soil, permutations = 999)
efCA
plot(trees.ca, display = "sites")
plot(efCA)
plot(efCA, p.max = 0.05) # only Al and Humdepth showed significant results for p < 0.05

#Fit environmental data in CA (categorical data = factors)
soil2<-read.table("soil2.txt", sep=",", header=TRUE)
trees.CA2 <- cca(trees)
efCA2 <- envfit(trees.CA2, soil2, permutations = 999)
efCA2

plot(trees.CA2, display = "sites")
plot(efCA2) #Not a great way to visualise data

plot(trees.CA2, display = "sites", type = "p")
with(soil2, ordiellipse(trees.CA2, Manag, conf = 0.95, col = "green", label = TRUE)) #draws ellipses for class standard deviations, standard er- rors or confidence areas
with(soil2, ordispider(trees.CA2, Manag, col = "pink")) #combines items to their (weighted) class centroid
with(soil2, ordihull(trees.CA2, Manag, col="red", lty=2)) #draws an enclosing convex hull for the items in a class

#CA may fail with very long gradients, sites appearing at extremes, and with the presence of rare species
#In this case, one may visualise the DCA
plot(trees.dca)
plot(trees.dca, display="sites")

######################################
##Example 2: Constrained ordination ##
######################################
# In this exercice species data is constrained by environmental data; 
# It's usefull as an exploratory tool for an direct analysis, when: 
# (1) you have both tree species and environmental data
# (2) you want to test for relations between tree species and environment

#CCA
trees.cca <- cca(trees ~ Al + Humdepth, soil) #constructing the multivariate model with the two significant environmental variables Al and Humdepth
trees.cca 

#The output is similar as in unconstrained ordination, but total inertia is decomposed into constrained and unconstrained components.
#There were 2 constraints, and the rank of constrained component is 2 
#The rank of unconstrained component is 20, when it was 23 in the previous analysis (e.g. CA)
#In some cases, the ranks may be lower than the number of constraints, since some of the constraints are dependent on each other, and an informative message is printed with the result
plot(trees.cca)

trees.cca2 <- cca(trees ~ Invasion, soil2)
plot(trees.cca2)
trees.cca2

#Test significance of all predictors in the model
anova(trees.cca2)

#Test significance of each predictor in the model
trees.cca3 <- cca(trees ~ pH + Invasion, soil2)
anova(trees.cca3)
anova(trees.cca3, by = "term", step=200)

#Building models
mod1 <- cca(trees ~ ., soil2)
mod1 #constrained model with all variables

mod0 <- cca(trees ~ 1, soil2)
mod <- step(mod0, scope = formula(mod1), test = "perm")
mod #constrained model using stepwise selection

plot(procrustes(cca(trees), mod)) #Comparing with an unconstrained ordination

#Conditioning models
trees.cca <- cca(trees ~ Manure + Condition(pH), soil2) #constrained model that evaluates the effect of management, after eliminating the natural effect of pH
plot(trees.cca)
trees.cca #The total inertia is decomposed into three components: inertia explained by conditions, inertia explained by constraints and the remaining unconstrained inertia.

anova(trees.cca, perm.max = 2000) #permutation for manure in the conditional model
anova(cca(trees ~ Manure, soil2)) #permutation for manure alone

######################################
##Example 3: Tests of dissimilarity ##
######################################

#MANOVA - partitions dissimilarities for the sources of variation, and uses permutation tests to inspect the significances of those partitions
betad <- betadiver(trees, "z")
adonis(betad ~ Manure, soil2, perm=200) #for differences among factors
adonis(betad ~ pH*Manure, soil2, perm = 200) #for interactions between variables

#Mantel test - it is the correlation between dissimilarity entries
pc <- prcomp(trees, scale = TRUE)
pc<- scores(pc, display = "sites", choices = 1:4)
edis <- vegdist(pc, method = "euclid")
vare.dis <- vegdist(wisconsin(sqrt(trees)))
mantel(vare.dis, edis)

plot(vare.dis, edis) #Plot two dissimilarity entities against each other


################################
##Example 4: Cluster analysis ##
################################

#There are numerous specific packages for clustering analysis, and derivates
#Most popular methods are single linkage (or nearest neighbour), complete linkage (furthest neighbour), and various brands of average linkage methods

dis <- vegdist(trees) #single linkage similarity among species
clus <- hclust(dis, "single")
plot(clus)

cluc <- hclust(dis, "complete") #complete linkage
plot(cluc)

clua <- hclust(dis, "average") #average linkage
plot(clua)

range(dis) #dissimilarity ranges

#Cophenetic correlation measures the similarity between original dissimilarities and dissimilarities estimated from the tree
cor(dis, cophenetic(clus))
cor(dis, cophenetic(cluc))
cor(dis, cophenetic(clua)) #This performs better!

#Dendograms

#Based on CA that structures table optimally into a diagonal structure. We can use its first axis to reorder the tree:
wa <- scores(trees, display = "sites", choices = 1)
den <- as.dendrogram(clua) #Tree 1
oden <- reorder(den, wa, mean) #Tree 2

plot(den)
plot(oden)
#Similar dendograms with distinct number of members

###########################


