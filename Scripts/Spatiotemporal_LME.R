#installing and loading package
install.packages("fslmer")
devtools::install_github("Deep-MI/fslmer", build_vignettes=TRUE)
library(fslmer)

#set wd
setwd("/p01-hdd/dsb/asever/For_Aleksandar/Both_TimePoints/")

##### lh thickness analysis #####
# read lh.thickness file
lh.thickness <- lme_openmgh("lh.thickness_sm10.mgh")

# read template surface
lh.sphere <- lme_readsurf("fsaverage/surf/lh.sphere")

# read cortical label file
lh.cortex <- lme_readlabel("fsaverage/label/lh.cortex.label")[,1]

# read qdec table
dat <- read.csv("long.qdec.2tps.manually.dat", sep = "\t")

# order qdec table based on fsid.base and years
dat <- dat[order(dat$fsid.base, dat$years), ]

# create vector of number of observations per subject
ni <- matrix(unname(table(dat$fsid.base)), ncol=1)

# create nSubjects-times-nVertices matrix for thickness data
Y <- t(drop(lh.thickness$x))

# create (1-based) vector of cortical label indices (to exclude non-cortical vertices)
maskvtx <- sort(lh.cortex)+1

# create model matrix
X <- model.matrix(~years*group, dat)

# create contrasts
C <- matrix(c(0, 0, 0, 0, 1, 0), nrow=1)

# determine the type of mixed-effects model (assuming the intercept is column 1 
# in the design matrix and the time variable is column 2).
#Zcols <- 1 # random-intercept 
Zcols <- c(1, 2) # random-slope, random-intercept

#mass univariate analysis
stats <- lme_mass_fit_vw(X, Zcols, Y, ni, numcore=6)
F_C <- lme_mass_F(stats, C)
# obtain initial estimates
FitInit <- lme_mass_fit_init(X=X, Zcols=Zcols, Y=Y, ni=ni, maskvtx=maskvtx, numcore=8)

# run algorithm to identify spatially homogeneous regions
RgGrow <- lme_mass_RgGrow(lh.sphere, FitInit$Re0, FitInit$Theta0, maskvtx=maskvtx, nst=2, prc=95)

# fit model
FitRgw <- lme_mass_fit_Rgw(X, Zcols, Y, ni, FitInit$Theta0, RgGrow$Regions, lh.sphere, prs=4)


##### rh thickness analysis #####

# read rh.thickness file
rh.thickness <- lme_openmgh("rh.thickness_sm10.mgh")

# read template surface
rh.sphere <- lme_readsurf("fsaverage/surf/rh.sphere")

# read cortical label file
rh.cortex <- lme_readlabel("fsaverage/label/rh.cortex.label")[,1]

# read qdec table
dat <- read.csv("long.qdec.2tps.manually.dat", sep = "\t")

# order qdec table based on fsid.base and years
dat <- dat[order(dat$fsid.base, dat$years), ]

# create vector of number of observations per subject
ni <- matrix(unname(table(dat$fsid.base)), ncol=1)

# create nSubjects-times-nVertices matrix for thickness data
Y <- t(drop(rh.thickness$x))

# create (1-based) vector of cortical label indices (to exclude non-cortical vertices)
maskvtx <- sort(rh.cortex)+1

# create model matrix
X <- model.matrix(~years*group, dat)

# create contrasts
C <- matrix(c(0, 0, 0, 0, 1, 1), nrow=1)

# determine the type of mixed-effects model (assuming the intercept is column 1 
# in the design matrix and the time variable is column 2).
#Zcols <- 1 # random-intercept 
Zcols <- c(1, 2) # random-slope, random-intercept

stats <- lme_mass_fit_vw(X, Zcols, Y, ni, numcore=6)

# obtain initial estimates
FitInit_rh.thickness <- lme_mass_fit_init(X=X, Zcols=Zcols, Y=Y, ni=ni, maskvtx=maskvtx, numcore=4)

# run algorithm to identify spatially homogeneous regions
RgGrow_rh.thickness <- lme_mass_RgGrow(rh.sphere, FitInit_rh.thickness$Re0, FitInit_rh.thickness$Theta0, maskvtx=maskvtx, nst=2, prc=95)

# fit model
FitRgw_rh.thickness <- lme_mass_fit_Rgw(X, Zcols, Y, ni, FitInit_rh.thickness$Theta0, RgGrow_rh.thickness$Regions, rh.sphere, prs=4)


##### lh area analysis #####

# read lh.area file
lh.area <- lme_openmgh("lh.area_sm10.mgh")

# read template surface
lh.sphere <- lme_readsurf("fsaverage/surf/lh.sphere")

# read cortical label file
lh.cortex <- lme_readlabel("fsaverage/label/lh.cortex.label")[,1]

# read qdec table
dat <- read.csv("long.qdec.2tps.manually.dat", sep = "\t")

# order qdec table based on fsid.base and years
dat <- dat[order(dat$fsid.base, dat$years), ]

# create vector of number of observations per subject
ni <- matrix(unname(table(dat$fsid.base)), ncol=1)

# create nSubjects-times-nVertices matrix for area data
Y <- t(drop(lh.area$x))

# create (1-based) vector of cortical label indices (to exclude non-cortical vertices)
maskvtx <- sort(lh.cortex)+1

# create model matrix
X <- model.matrix(~years*group, dat)

# create contrasts
C <- matrix(c(0, 0, 0, 0, 1, 1), nrow=1)

# determine the type of mixed-effects model (assuming the intercept is column 1 
# in the design matrix and the time variable is column 2).
#Zcols <- 1 # random-intercept 
Zcols <- c(1, 2) # random-slope, random-intercept

# obtain initial estimates
FitInit_lh.area <- lme_mass_fit_init(X=X, Zcols=Zcols, Y=Y, ni=ni, maskvtx=maskvtx, numcore=4)

# run algorithm to identify spatially homogeneous regions
RgGrow_lh.area <- lme_mass_RgGrow(lh.sphere, FitInit_lh.area$Re0, FitInit_lh.area$Theta0, maskvtx=maskvtx, nst=2, prc=95)

# fit model
FitRgw_lh.area <- lme_mass_fit_Rgw(X, Zcols, Y, ni, FitInit_lh.area$Theta0, RgGrow_lh.area$Regions, lh.sphere, prs=4)

##### rh area analysis #####

# read rh.area file
rh.area <- lme_openmgh("rh.area_sm10.mgh")

# read template surface
rh.sphere <- lme_readsurf("fsaverage/surf/rh.sphere")

# read cortical label file
rh.cortex <- lme_readlabel("fsaverage/label/rh.cortex.label")[,1]

# read qdec table
dat <- read.csv("long.qdec.2tps.manually.dat", sep = "\t")

# order qdec table based on fsid.base and years
dat <- dat[order(dat$fsid.base, dat$years), ]

# create vector of number of observations per subject
ni <- matrix(unname(table(dat$fsid.base)), ncol=1)

# create nSubjects-times-nVertices matrix for area data
Y <- t(drop(rh.area$x))

# create (1-based) vector of cortical label indices (to exclude non-cortical vertices)
maskvtx <- sort(rh.cortex)+1

# create model matrix
X <- model.matrix(~years*group, dat)

# create contrasts
C <- matrix(c(0, 0, 0, 0, 1, 1), nrow=1)

# determine the type of mixed-effects model (assuming the intercept is column 1 
# in the design matrix and the time variable is column 2).
#Zcols <- 1 # random-intercept 
Zcols <- c(1, 2) # random-slope, random-intercept

# obtain initial estimates
FitInit_rh.area <- lme_mass_fit_init(X=X, Zcols=Zcols, Y=Y, ni=ni, maskvtx=maskvtx, numcore=4)

# run algorithm to identify spatially homogeneous regions
RgGrow_rh.area <- lme_mass_RgGrow(rh.sphere, FitInit_rh.area$Re0, FitInit_rh.area$Theta0, maskvtx=maskvtx, nst=2, prc=95)

# fit model
FitRgw_rh.area <- lme_mass_fit_Rgw(X, Zcols, Y, ni, FitInit_rh.area$Theta0, RgGrow_rh.area$Regions, rh.sphere, prs=4)





