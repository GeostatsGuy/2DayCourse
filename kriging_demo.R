# Kriging Tutorial in R for Engineers and Geoscientists 

# Michael Pyrcz, University of Texas at Austin, Twitter @GeostatsGuy



# This will be used in my Introduction to Geostatistics undergraduate class 

# It is assumed that students have no previous R experience.  

# This utilizes the gstat library by Edzer Pedesma, appreciation to Dr. Pedesma for assistance.

# Load the required libraries, you may have to first go to "Tools/Install Packages..." to install these first

library(gstat)                                 # geostatistical methods by Edzer Pebesma
library(sp)                                    # spatial points addition to regular data frames
library(plyr)                                  # splitting, applying and combining data by Hadley Wickham 

# Specify the grid parameters (same as GSLIB / GEO-DAS parameterization in 2D)

nx = 400
ny = 400
xmin = 5.0
ymin = 5.0
xsize = 10.0
ysize = 10.0

# Declare functions

# This function completes standard normal transform on a data vector
nscore <- function(x) {                        # written by Ashton Shortridge, May/June, 2008
  # Takes a vector of values x and calculates their normal scores. Returns 
  # a list with the scores and an ordered table of original values and
  # scores, which is useful as a back-transform table. See backtr().
  nscore <- qqnorm(x, plot.it = FALSE)$x  # normal score 
  trn.table <- data.frame(x=sort(x),nscore=sort(nscore))
  return (list(nscore=nscore, trn.table=trn.table))
}

# This function builds a spatial points dataframe with the locations for estimation / simulation 
addcoord <- function(nx,xmin,xsize,ny,ymin,ysize) { # Michael Pyrcz, March, 2018                      
  # makes a 2D dataframe with coordinates based on GSLIB specification
  coords = matrix(nrow = nx*ny,ncol=2)
  ixy = 1
  for(iy in 1:nx) {
    for(ix in 1:ny) {
      coords[ixy,1] = xmin + (ix-1)*xsize  
      coords[ixy,2] = ymin + (iy-1)*ysize 
      ixy = ixy + 1
    }
  }
  coords.df = data.frame(coords)
  colnames(coords.df) <- c("X","Y")
  coordinates(coords.df) =~X+Y
  return (coords.df)
}  

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/PGE337")

# Read the data table from a comma delimited file - data on GitHub/GeostatsGuy/GeoDataSets
mydata = read.csv("2D_MV_200Wells.csv")        # read in comma delimited data file

# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(mydata)                                   # show the first several rows of a data table in the console

# The columns are variables with variable names at the top and the rows are samples

# Convert the dataframe to a spatial points dataframe
class(mydata)                                  # confirms that it is a dataframe
coordinates(mydata) = ~X+Y                     # indicate the X, Y spatial coordinates
summary(mydata)                                # confirms that it is now a spatial points dataframe
head(coordinates(mydata))                      # check the first several coordinates

# Normal scores transform of the porosity data to assist with variogram calculation
npor.trn = nscore(mydata$porosity)             # normal scores transform
mydata[["NPorosity"]]<-npor.trn$nscore         # append the normal scores transform into the spatial data table
head(mydata)

# Now let's check the porosity data distribution
par(mfrow=c(2,2))                              # set up a 2x2 matrix of plots 
hist(mydata$porosity,main="Porosity (%)",xlab="Porosity (%)",nclass = 15) # histogram

# ecdf makes a cdf object and plot command plots it
plot(ecdf(mydata$porosity),main="Porosity",xlab="Porosity (%",ylab="Cumulative Probability") # CDF

# Now let's check the normal score transform of the porosity data distribution
hist(mydata$NPorosity,main="N[Porosity (%)]",xlab="N[Porosity (%)]",nclass = 15) # histogram
# ecdf makes a cdf object and plot command plots it
plot(ecdf(mydata$NPorosity),main="N[Porosity]",xlab="N[Porosity (%)]",ylab="Cumulative Probability") #CDF

# We can plot the porosity data location map
spplot(mydata, "porosity", do.log = TRUE,      # location map of porosity data
       key.space=list(x=1.05,y=0.97,corner=c(0,1)),
       scales=list(draw=T),xlab = "X (m)", ylab = "Y (m)",main ="Porosity (%)")

# We can also plot a simple bubble plot of the porosity data
bubble(mydata, "NPorosity", fill = FALSE, maxsize = 2, main ="Porosity (%)",identify = FALSE,xlab = "X (m)", ylab = "Y (m)")

# Anisotropic variogram models from last tutorial, note anisotropy is parameterized as c(azimuth,dip,plunge,hratio,vratio) in #3D and c(azimuth,hratio) in 2D.
par(mfrow=c(1,1))

por.vm.ani <- vgm(psill = 0.6, "Exp", 800, anis = c(035, 0.5),nugget=0.4)
por.vm.ani                                     # check the variogram model parameters 

# Anisotropic variogram plots
name = c("035","125")                          # make name matrix
color = c("blue","red")                        # make color matrix

por.vg.035 = variogram(NPorosity~1,mydata,cutoff = 3000,width =500,alpha = 35.0,tol.hor=22.5) # 035 directional 
por.vg.125 = variogram(NPorosity~1,mydata,cutoff = 3000,width =500,alpha = 125.0,tol.hor=22.5) # 125 directional

plot(por.vg.035$dist,por.vg.035$gamma,main="Porosity Anisotropic Variogram",xlab="  Lag Distance (m) ",ylab=" Semivariogram ", col=color[1],ylim=c(0,1.2))
points(por.vg.125$dist,por.vg.125$gamma,col=color[2])
abline(h = 1.0)

unit_vector = c(sin(35*pi/180),cos(35*pi/180),0) # unit vector for 035 azimuth
vm.ani.035 <- variogramLine(por.vm.ani,maxdist=3000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 035 variogram model
lines(vm.ani.035$dist,vm.ani.035$gamma,col=color[1]) # include variogram model 

unit_vector = c(sin(55*pi/180),-1*cos(35*pi/180),0) # unit vector for 125 azimuth
vm.ani.125 <- variogramLine(por.vm.ani,maxdist=3000,min=0.0001,n=100,dir=unit_vector,covariance=FALSE) # calculate 125 variogram model
lines(vm.ani.125$dist,vm.ani.125$gamma,col=color[2]) # include variogram model

legend(2000,.8,name, cex=0.8, col=color,pch=c(21,21,21),lty=c(1,1,1)) # add legend

# Make the model 2D grid

coords <- addcoord(nx,xmin,xsize,ny,ymin,ysize) # make a dataframe with all the estimation locations
summary(coords)                                # check the coordinates

# First attempt inverse distance on this regular grid
# We will try 3 different powers to demonstrate the influence on the interpolation.  Higher powers result 
# in more locally specific models, because distance has a greater influence.

cuts = c(.05,.07,.09,.11,.13,.15,.17,.19,.21,.23)
cuts.var = c(0.05,.1,.15,.20,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95)
cuts.var

porosity.idw.2.0 = idw(porosity~1, idp = 2.0,mydata, coords) # inverse distance power 2.0
class(porosity.idw.2.0)                        # check the inverse distance object
spplot(porosity.idw.2.0["var1.pred"],main = "Porosity Inverse Distance (p=2.0)", key.space = "right",cuts = cuts,xlab = "X (m)", ylab = "Y (m)")

porosity.idw.1.0 = idw(porosity~1, idp = 1.0,mydata, coords) # inverse distance power 1.0
spplot(porosity.idw.1.0["var1.pred"],main = "Porosity Inverse Distance (p=1.0)", key.space = "right",cuts = cuts,xlab = "X (m)", ylab = "Y (m)")

porosity.idw.3.0 = idw(porosity~1, idp = 3.0,mydata, coords) # inverse distance power 3.0
spplot(porosity.idw.3.0["var1.pred"],main = "Porosity Inverse Distance (p=3.0)", key.space = "right",cuts = cuts,xlab = "X (m)", ylab = "Y (m)")

# Let's try some ordinary kriging with unlimited search

porosity.kriged = krige(porosity~1, mydata, coords, model = por.vm.ani,maxdist = Inf,nmin = 0,omax=Inf) # ordinary kriging
summary(porosity.kriged)

# Visualize the estimates
spplot(porosity.kriged["var1.pred"],main = "Porosity Ordinary Kriging", key.space = "right",cuts = cuts,xlab = "X (m)", ylab = "Y (m)")

# Visualize the estimation variance, given sill = 1.0, the maximum estimation variance is 1.0.
spplot(porosity.kriged["var1.var"],main = "Porosity Ordinary Kriging Variance", key.space = "right",cuts = cuts.var,xlab = "X (m)", ylab = "Y (m)")

# Set up limited search parameters repeat the kriging
maxdist = 800                                  # maximum distance to look for data
nmin = 3                                       # minimum number of data for an estimate
omax = 1                                       # maximum number of data per octant    

porosity.kriged.sp = krige(porosity~1, mydata, coords, model = por.vm.ani,maxdist = maxdist,nmin = nmin,omax=omax) # ordianry kriging
spplot(porosity.kriged.sp["var1.pred"],main = "Porosity Ordinary Kriging / Limited Search", key.space = "right",cuts = cuts,xlab = "X (m)", ylab = "Y (m)")

# Observed the search artifacts (lines and bands) due to too restrictive search 

# Let's remove the nugget effect and krige again once again with unlimited search
# Anisotropic variogram with nugget effect removed
por.vm.ani.nonugget <- vgm(psill = 1.0, "Exp", 800, anis = c(035, 0.5),nugget=0.0)
por.vm.ani.nonugget  

porosity.kriged.nonugget = krige(porosity~1, mydata, coords, model = por.vm.ani.nonugget,maxdist = Inf,nmin = 0,omax=Inf) # ordianry kriging
spplot(porosity.kriged.nonugget["var1.pred"],main = "Porosity Ordinary Kriging No Nugget", key.space = "right",cuts = cuts,xlab = "X (m)", ylab = "Y (m)")
spplot(porosity.kriged.nonugget["var1.var"],main = "Porosity Ordinary Kriging Variance No Nugget", key.space = "right",cuts = cuts.var,xlab = "X (m)", ylab = "Y (m)")

# On your own try changing the variogram parameters and observe the results.  Also consider kriging with a trend.

# Hope this was helpful,

# Michael


