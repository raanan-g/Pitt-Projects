# R-Script Project 2 Spatial Data Analytics:

# Packages
library(plyr)
library(dplyr)
library(data.table)
library(rgdal)
library(spatstat)

# Read shapefile
OGdata <- readOGR(getwd(),"OilGasLocationPA")
plot(OGdata)
#####################################
###### Quadrat Counts Method ########
#####################################

# Get the boundaries of the data set 
# to calculate the area of the study region
minx <- min(OGdata$LATITUDE)
maxx <- max(OGdata$LATITUDE)
miny <- min(OGdata$LONGITUDE)
maxy <- max(OGdata$LONGITUDE)
# Area = a
a <- ((maxy-miny)*(maxx-minx))

# Quadrat area = {a/X, a(sqrt(2))/X}, where X = number of observations and 
# sqrt(2) represents regularization
X <- length(OGdata$LATITUDE)
aq1 <- a/X
aq2 <- (a/X)*sqrt(2)

#########################################
###### Regular sampling approach ########
#########################################

# Generate point pattern from OGdata and  
coords <- ppp(OGdata$LATITUDE,OGdata$LONGITUDE,c(minx,maxx),c(miny,maxy))
qcounts <- quadratcount(coords,10,10)
# For second quadrat size
qcounts <- quadratcount(coords,20,20)

# Plot contingency table of qcounts with original points
plot(coords, pch=16,cex=0.5,main="Oil-Gas Locations Quadrat Counts")
plot(qcounts, add=TRUE, entries=NULL)

# Plot intensity to visualize qcounts 
plot(intensity(qcounts, image=TRUE))

# Generate table of k, X, k - u, (k-u)^2, X(k-u)^2
qcframe <- data.frame(qcounts)
qctable <- data.frame(table(qcframe$Freq, exclude=NULL))
mu <- mean(qctable$Freq)
Xsum <- sum(qctable$Freq)

##########################################
######## Random sampling approach ########
##########################################

# Set initial boundaries

# Number of randomly generated quadrats
totalnumber=100

plot(OGdata)

x<- list()

i<-0
# Generate 100 randomly placed quadrats
while(i < 100) {
  # Generate random coordinates within study area
  rX<-runif(1,minx,maxx)
  rY<-runif(1,miny,maxy)
  
  # Make sure the whole square is within the boundaries
  if(rX+((maxx-minx)/10)<maxx && rY+((maxx-minx)/10) < maxy) {
    
    # Generate the square 
    rect(rX, rY, rX+((maxx-minx)/10), rY+((maxx-minx)/10), border="red")
    
    # Represent square as w (window)
    w<-owin(c(rX,rX+((maxx-minx)/10)),c(rY, rY+((maxx-minx)/10)))
    
    # Return list (x) of quadrat counts
    ok<-inside.owin(coords$x, coords$y, w)
    m <-  length(which(ok, arr.ind = TRUE))
    x[[paste0("val1", i)]] <- m
  
    i <- i+1
  }
}

# Convert list x to dataframe to generate table of statistics
df <- ldply(x, data.frame)
setnames(df, old=c(".id","X..i.."), new=c("id", "num"))
w2 = table(df$num)
class(w2)
t2 = as.data.frame(w2)

# Add columns with desired statistics
t2$Var1 <-as.numeric(as.character(t2$Var1))
t2$mu <- with(t2, Var1 - 156)
t2$mu_squared <-with(t2, mu^2)
t2$x_to_musq <- with(t2, mu_squared*Freq)

# Print out table with statistics
t2










