
## Models in Invasion Ecology Advanced Course 2019
##
## R Installation script ----------------------------------------------------------------------


# Install SegOptim package and dependencies ---------------------------------------------------
				 
if(!("devtools" %in% installed.packages()[,1])){ 
	
	install.packages("devtools") 
}

library(devtools) 

install_bitbucket("joao_goncalves/segoptim")


## Check if packages can be loaded -----------------------------------------------------------

library(raster) 

library(randomForest) 

library(SegOptim)


## Check if help info is available -----------------------------------------------------------
## (this info will appear in the 'Help' tab in RStudio)

?segmentation_OTB_LSMS 

# or 

?calibrateClassifier


## Install other necessary packages ----------------------------------------------------------

install.packages(c("biomod2","vegan"), dependencies=TRUE)

## Check if packages can be loaded

library(biomod2) 

library(vegan) 

