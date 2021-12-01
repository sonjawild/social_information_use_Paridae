
# ## EXPERIMENTAL EVIDENCE THAT TITS (PARIDAE SP.) USE SOCIAL INFORMATION WHEN LOCATING AND CHOOSING NEST LINING MATERIAL-----------------------------------------------

# 1) Load necessary libraries -------------------------------------------------------

library(asnipe)
library(sp)
library(geosphere)
library(Imap)
library(data.table)
library(ggplot2)
library(gridExtra)


# 2) Create foraging network ------------------------------------------------

# 2.1. Read gmm data (association data) --------------------------------------------------

load("gmm.spring.RData")
# the object contains three slots:
gmm.spring$gbi # group by individual matrix
gmm.spring$metadata
gmm.spring$B

# 2.2. Calculate association network using SRI --------------------------------------

foraging_network <-
  get_network(
    association_data = gmm.spring$gbi,
    data_format = "GBI",
    association_index = "SRI"
  )

dim(foraging_network)


# 3) Create distance matrix ----------------------------------------------------


# 3.1. Read GPS data ------------------------------------------------------

GPS_data <- read.csv("coordinates_boxes+dispensers_ALL.csv")

# 3.2. Load functions to calculate distances in m -------------------------

##  GeoDistanceInMetresMatrix() function that generates the distance matrix
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}


# 3.3. calculate distances between nest boxes and dispensers --------------

head(GPS_data)
coordinates <- data.frame(name = GPS_data$Box,
                          lat  = GPS_data$Lat,
                          lon  = GPS_data$Long)

distance_matrix<-round(GeoDistanceInMetresMatrix(coordinates))
head(distance_matrix)


# 3.4. Extract which boxes are within 200m of dispensers ------------------
neighbour_matrix <- distance_matrix
# m contains the boxes within 200m of the dispensers
m <- neighbour_matrix[(length(rownames(neighbour_matrix))-4):length(rownames(neighbour_matrix)),1:(length(rownames(neighbour_matrix))-5)]<=200

boxes.to.include <- unique(colnames(m)[col(m)[which(m)]])
length(boxes.to.include)
# 175 boxes are within a 200m radius of the dispensers

# 3.5. Calculate the neighbour matrix -------------------------------------

# we set all values above 50m to 0
# and then use inverted distances
# so that boxes closer together have higher values

neighbour_matrix <- distance_matrix

neighbour_matrix[neighbour_matrix>=50] <- 0

neighbour_matrix <- 1/neighbour_matrix

# set all infinity values to 0
neighbour_matrix[neighbour_matrix==Inf] <- 0
hist(neighbour_matrix[neighbour_matrix!=0])


# 4) Load data from wool dispensers -----------------------------------------

load("dispenser.data.RDA")

dispenser.data.comb <- rbind(dispenser.data$dispenser.data.1, 
                             dispenser.data$dispenser.data.2,
                             dispenser.data$dispenser.data.3,
                             dispenser.data$dispenser.data.4,
                             dispenser.data$dispenser.data.5)

# How many different birds have visited (this includes all species and also potentially males)
length(unique(dispenser.data.comb$PIT))
# 39


# # 5) Load individual-level variables (ILVs) ---------------------------------------------------------

load("ILVs.combined.RDA")

# 6) NBDA - social information to use to find lining material -----------------------------------------------------------------
# load NBDa via devtools
#install.packages("devtools")
#install.packages("rtools40")
#install.packages("ade4")
library(devtools)
# install_github("whoppitt/NBDA")
library(NBDA)

# 6.1. prepare matrices -------------------------------------------------

# extract females that were breeding in our study area and were part of the foraging network
IDs.to.include.in.NBDA <- intersect(ILVs.combined$PIT_f, rownames(foraging_network))
length(IDs.to.include.in.NBDA)
# 46 tagged females were seen at the foraging network that were also breeding in our study area

# extract how many of those females have visited the dispensers
ILVs.combined.sub <- subset(ILVs.combined, ILVs.combined$PIT_f %in% IDs.to.include.in.NBDA)
# remove three boxes (as they were double occupations by the same females - c("D04", "R06", "G33")) 
ILVs.combined.sub <- subset(ILVs.combined.sub, !(ILVs.combined.sub$Box %in% c("D04", "R06", "G33")))
length(subset(ILVs.combined.sub$PIT_f, ILVs.combined.sub$D1.visited==1 | ILVs.combined.sub$D2.visited==1| ILVs.combined.sub$D3.visited==1| ILVs.combined.sub$D4.visited==1| ILVs.combined.sub$D5.visited==1))
# [1] 21

# extract how many were breeding in nest boxes
length(subset(ILVs.combined.sub$PIT_f, ILVs.combined.sub$Box!="no_box"))
# [1] 36 breeding in the boxes, and 10 presumably in natural tree cavities

# extract how many visits these 21 females have made to the dispenser
dispenser.data.comb.sub <- subset(dispenser.data.comb, dispenser.data.comb$PIT %in% IDs.to.include.in.NBDA)
count <- 1
visits <- NULL
for(i in IDs.to.include.in.NBDA){
  sub <- subset(dispenser.data.comb.sub, dispenser.data.comb.sub$PIT==i)
  num.visits <- length(sub[,1])
  visits[count] <- num.visits
  count <- count+1
}
# extract the average and max number of visits of the females that have been recorded at the dispenser
mean(visits[visits!=0])
max(visits)

# we now get the two networks into the right shape - i.e. subsetting them to the 46 females we want to include
# first the foraging network
foraging.network.NBDA <- foraging_network[rownames(foraging_network) %in% IDs.to.include.in.NBDA, colnames(foraging_network) %in% IDs.to.include.in.NBDA]
foraging.network.NBDA <- foraging.network.NBDA[order(rownames(foraging.network.NBDA)), order(colnames(foraging.network.NBDA))]
dim(foraging.network.NBDA)
hist(foraging.network.NBDA)

# subset the neighbour network to these boxes that these females breed in
boxes.to.include.in.NBDA <- unique(subset(ILVs.combined$Box, ILVs.combined$PIT_f %in% IDs.to.include.in.NBDA))
boxes.to.include.in.NBDA <- boxes.to.include.in.NBDA[boxes.to.include.in.NBDA!="no_box"]
new.names.all <- NULL
neighbour_matrix.NBDA <- neighbour_matrix[rownames(neighbour_matrix) %in% boxes.to.include.in.NBDA, colnames(neighbour_matrix) %in% boxes.to.include.in.NBDA]
for(i in rownames(neighbour_matrix.NBDA)){
  new.name <- subset(ILVs.combined$PIT_f, ILVs.combined$Box==i)
  new.names.all[which(rownames(neighbour_matrix.NBDA)==i)] <- new.name
}
rownames(neighbour_matrix.NBDA) <- new.names.all
colnames(neighbour_matrix.NBDA) <- new.names.all
dim(neighbour_matrix.NBDA)
# this matrix only contains 36 (+3 double ups of females starting to build nests in two boxes) individuals 
# (excludes the ones that weren't breeding in nest boxes)
# we add those individuals and assign 0 in their neighbour connections, as we do not know who their neighbours are

neighbour_matrix.NBDA.new <- matrix(0, ncol=length(rownames(foraging.network.NBDA)), nrow=length(rownames(foraging.network.NBDA)), dimnames = list(rownames(foraging.network.NBDA),rownames(foraging.network.NBDA)))
# now fill in the real distance values for those available
cols <- colnames(neighbour_matrix.NBDA.new)[colnames(neighbour_matrix.NBDA.new) %in% colnames(neighbour_matrix.NBDA)]
rows <- rownames(neighbour_matrix.NBDA.new)[rownames(neighbour_matrix.NBDA.new) %in% rownames(neighbour_matrix.NBDA)]

neighbour_matrix.NBDA.new[rows, cols] <- neighbour_matrix.NBDA[rows,cols]
dim(neighbour_matrix.NBDA.new)
hist(neighbour_matrix.NBDA.new)
# to ensure they are on a similar scale, we multiply the values in the neighbour matrix *10
neighbour_matrix.NBDA.new <- neighbour_matrix.NBDA.new*10
hist(neighbour_matrix.NBDA.new)

# 6.2. Check for correlation between foraging and neighbour matrix ----------------------------------------------

library(vegan)
 mantel(neighbour_matrix.NBDA.new, foraging.network.NBDA, permutations = 9999)
#  Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = neighbour_matrix.NBDA.new, ydis = foraging.network.NBDA,      permutations = 9999) 
# 
# Mantel statistic r: 0.04815 
#       Significance: 0.0826 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0429 0.0610 0.0783 0.1010 
# Permutation: free
# Number of permutations: 9999

# non-significant correlation between foraging and neighbour network, which means we can include them both in the NBDA analysis at the same time
 
# 6.3. Prepare NBDA data ------------------------------------------------------------------

 # this function automatically gets the data into shape for NBDA analysis
 # it subsets the data to those individuals breeding around each dispenser area (within 200m)
 # and extracts ILVs for these individuals and gets them into the right format

 
prepare.NBDA.data <- function(dispenser.data, include.all, ILVs.include){
  dispenser.data <- dispenser.data[order(dispenser.data$date.time),] # ensure it is sorted according to date
  location <- unique(dispenser.data$Location)
  
  # We extract which females were breeding in boxes within 200m of the respective dispenser
  IDs.included <- subset(ILVs.combined$PIT_f, ILVs.combined[,location]<=200 | ILVs.combined[,paste(location, ".visited", sep="")]==1 ) 
  IDs.included <- subset(IDs.included, IDs.included %in% IDs.to.include.in.NBDA)
  # we remove boxes D04, R06, G33 (females breeding in two boxes - we retained the ones closer to the dispenser)
  ILVs.sub.disp <- subset(ILVs.combined, ILVs.combined$PIT_f %in% IDs.included & !(ILVs.combined$Box %in% c("D04", "R06", "G33")))

  # subset the dispenser data to only those females
  dispenser.data <- subset(dispenser.data, dispenser.data$PIT %in%  IDs.included)
  
  # order data in ascending order (according to PIT tag)
  ILVs.sub.disp <- ILVs.sub.disp[order(ILVs.sub.disp$PIT_f),] # order ascending according to PIT tag
  
  # subset the two networks to those IDs
  forage.net <- foraging.network.NBDA[rownames(foraging.network.NBDA) %in% IDs.included, colnames(foraging.network.NBDA) %in% IDs.included]
  neighbour.net <- neighbour_matrix.NBDA.new[rownames(neighbour_matrix.NBDA.new) %in% IDs.included, colnames(neighbour_matrix.NBDA.new) %in% IDs.included]
  # create an array with the two matrices
  assMatrix.nbda <- array(data=c(forage.net, neighbour.net), dim=c(nrow(forage.net), ncol(forage.net), 2))
  
  # create objects in the global environment for each ILV
  species.nbda <- as.matrix(ILVs.sub.disp$Species) 
  species.nbda[species.nbda!="GRETI"] <- -0.5 # we assign -0.5 for non great tits
  species.nbda[species.nbda=="GRETI"] <- 0.5 # we assign 0.5 for great tits
  species.nbda <- as.matrix(as.numeric(species.nbda)) # ILVs need to be defined as matrices
  
  age.nbda <- as.matrix(ILVs.sub.disp$Age) 
  age.nbda[age.nbda=="first.year"] <- -0.5 # -0.5 for juveniles
  age.nbda[age.nbda=="adult"] <- 0.5 # 0.5 for adult birds
  age.nbda <- as.matrix(as.numeric(age.nbda))
  
  distance.nbda <- as.matrix(as.numeric(ILVs.sub.disp[,location])) # we extract the distance of the box to the respective dispenser
  distance.nbda[is.na(distance.nbda)] <- mean(distance.nbda[!is.na(distance.nbda)]) # for those not breeding in boxes, add the average distance to the respective dispenser
  # we standardize the distance for better model fitting - using the standard deviation and mean
  distance.nbda <- as.matrix(as.numeric(scale(distance.nbda)))
  
  # we need to assign these ILVs as objects to the global environment
  assign(paste("species", location, sep="_"), species.nbda, envir = .GlobalEnv)
  assign(paste("age", location, sep="_"), age.nbda, envir = .GlobalEnv)
  assign(paste("distance", location, sep="_"), distance.nbda, envir = .GlobalEnv)
  

  ILVs <- paste(ILVs.include, location, sep="_")
  assign(paste("ILVs", location, sep="_"), ILVs, envir = .GlobalEnv)
  
  # extract the order of finding the dispenser from the dispenser data
  t.first <- dispenser.data[match(unique(dispenser.data$PIT), dispenser.data$PIT),]
  order <- NULL
  time <- NULL
  num.visits <- NULL
  for( i in t.first$PIT){
    order[which(t.first$PIT==i)] <- which(rownames(forage.net)==i)
    time[which(t.first$PIT==i)] <- as.POSIXct(as.character(subset(t.first$date.time, t.first$PIT==i)), format="%y%m%d%H%M%S", origin="1970-01-01")-as.POSIXct("21032612300000", format="%y%m%d%H%M%S") # difference in days
    # 26.03.21 12:30 CEST was the start date of the experiment
  }

  object <- NULL
  object$forage.net <- forage.net
  object$neighbour.net <- neighbour.net
  object$ILVs.full <- ILVs.sub.disp
  object$assMatrix <- assMatrix.nbda

  object$OAc <- order
  object$TAc <- time
  

  # finally, we create the NBDA data object, which is returned from running the function
  object2 <- nbdaData(label=location,                        
           assMatrix=object$assMatrix,          
           asoc_ilv=get(paste("ILVs", location, sep="_")),            
           int_ilv=get(paste("ILVs", location, sep="_")),            
           multi_ilv="ILVabsent",        
           orderAcq=object$OAc,          
           timeAcq=object$TAc,           
           endTime=41
  )

  return(object2)
}


# 6.4. Including all females ----------------------------------------------

 
# we run the code to prepare NBDA data objects including all females (those breeding in boxes + the ones who have visited the dispenser) 
# Note that the function will not work for D4 - as there was only one tagged learner - we exclude it from the NBDA
nbdaData_D1.all <- prepare.NBDA.data(dispenser.data = dispenser.data$dispenser.data.1, include.all = TRUE, ILVs.include = c("age", "species", "distance"))
nbdaData_D2.all <- prepare.NBDA.data(dispenser.data = dispenser.data$dispenser.data.2, include.all = TRUE, ILVs.include = c("age", "species", "distance"))
nbdaData_D3.all <- prepare.NBDA.data(dispenser.data = dispenser.data$dispenser.data.3, include.all = TRUE, ILVs.include = c("age", "species", "distance"))
# nbdaData_D4.all <- prepare.NBDA.data(dispenser.data = dispenser.data$dispenser.data.4, include.all = TRUE, ILVs.include = c("age", "species", "distance"))
nbdaData_D5.all <- prepare.NBDA.data(dispenser.data = dispenser.data$dispenser.data.5, include.all = TRUE, ILVs.include = c("age", "species", "distance"))

# extract the number of birds in each dispenser area and the number of learners
# D1
dim(nbdaData_D1.all@assMatrix[,,2,1])
length(nbdaData_D1.all@orderAcq)
# 9 birds, 5 learners
# D2
dim(nbdaData_D2.all@assMatrix[,,2,1])
length(nbdaData_D2.all@orderAcq)
# 26 birds, 7 learners
#D3
dim(nbdaData_D3.all@assMatrix[,,2,1])
length(nbdaData_D3.all@orderAcq)
# 19 birds, 5 learners

# skipping D4 as only one learner

#D5
dim(nbdaData_D5.all@assMatrix[,,2,1])
length(nbdaData_D5.all@orderAcq)
# 5 birds, 4 learners


# 6.5. Create constraints Vector Matrix ---------------------------------

# here, we provide a function to generate the constraints vector matrix
# which defines the parameter combinations of all NBDA models that are to be run

create.constraints.Vect.Matrix <- function(NBDA_data_object, num_networks, num_ILVs){
  suppressWarnings(
    if(NBDA_data_object@asoc_ilv=="ILVabsent"){
      num.ILV.asoc <- 0
    } else {num.ILV.asoc <- length(NBDA_data_object@asoc_ilv)})
  
  suppressWarnings(
    if(NBDA_data_object@int_ilv=="ILVabsent"){
      num.ILV.int<- 0
    } else {num.ILV.int<- length(NBDA_data_object@int_ilv)})
  
  suppressWarnings(
    if(NBDA_data_object@multi_ilv=="ILVabsent"){
      num.ILV.multi <- 0
    } else {num.ILV.multi <- length(NBDA_data_object@multi_ilv)})
  
  vector <- seq(1:(num_networks+num.ILV.asoc+num.ILV.int+num.ILV.multi))
  
  count <- 0 # create an object 'count', which starts on 0
  
  constraintsVect <- matrix(nrow = 10000000, ncol=(num_networks+num.ILV.asoc+num.ILV.int+num.ILV.multi)) # create a matrix to save the combination of parameters in
  constraintsVect[1,] <- vector # the first row gets filled with a sequence from 1:8 (all parameters will be estimated, none are set to 0)
  
  for (i in 1:(length(vector)-1)){ # a loop for each number of parameters to be estimated
    array <- combn(vector, i, FUN = NULL, simplify = TRUE) # for each number of paramters to be estiamted (e.g. 2) create all possible combinations of numbers between 1:12 (e.g. 2&8, 1&5 etc)
    
    for (j in 1:length(array[1,])){ # for each of those combinations
      vector2 <- seq(1:((num_networks+(num.ILV.asoc+num.ILV.int+num.ILV.multi))-i)) # create a second vector with 11-i free spaces
      position <- array[,j] # for each created combination
      count <- count+1 # add +1 to the count
      
      for (k in position){ # at each possible position
        vector2 <- append(vector2, 0, after=k-1) # add a 0 (e.g. 1 0 2 3 ...; 1 2 0 3 4 5 ...; 1 2 3 0 4 5 ....)
      }
      constraintsVect[count+1,] <- vector2 # and save the resulting order in a matrix
    }
  }
  
  
  constraintsVect <- na.omit(constraintsVect) # remove all NAs from the matrix
  
  # extract which columns are networks
  col.networks <- c(1:num_networks)
  
  col.names <- NULL
  
  if(num.ILV.asoc!=0){
    col.names <- rep("asoc", num.ILV.asoc)
  }
  
  if(num.ILV.int!=0){
    col.names <- c(col.names, rep("int", num.ILV.int))
  }
  
  if(num.ILV.multi!=0){
    col.names <- c(col.names, rep("multi", num.ILV.multi))
  }
  
  colnames(constraintsVect) <- c(rep("network", num_networks), col.names)
  
  constraintsVect <- as.matrix(as.data.frame(constraintsVect))
  
  # extract the models containing any social network
  
  social.models <- rep(NA, length(constraintsVect[,1]))
  
  for (k in 1:length(constraintsVect[,1])){
    sum <- sum(constraintsVect[k,1:num_networks])
    if(sum!=0){
      social.models[k] <- k
    }
  }
  social.models <- as.vector(na.omit(social.models))
  
  social.models.matrix <- constraintsVect[social.models,]
  
  # if multiplicative models are fit, we need to adjust the matrix
  # if the multiplicative slots are filled, it automatically fits the parameter for asoc and social (just constrained to be the same)
  # meaning that we can remove it from the asoc and int slot
  
  if(num.ILV.multi!=0){
    social.models.retain <- rep(NA, length(social.model.matrix[,1]))
    multi.models <- rep(NA, length(social.models.matrix[,1]))
    for (k in 1:length(social.models.matrix[,1])){
      sum <- sum(social.models.matrix[k,which(colnames(social.models.matrix)=="multi")])
      sum2 <- sum(social.models.matrix[k, c(which(colnames(social.models.matrix)=="asoc"),which(colnames(social.models.matrix)=="int"))])
      if(sum!=0 & sum2==0){ # if multi models are fit and int and asoc are set to 0
        multi.models[k] <- k # then retain the model
      } else if (sum==0){
        social.models.retain[k] <- k
      }
    }
    
    multi.models <- as.vector(na.omit(multi.models))
    social.models.retain <- as.vector(na.omit(social.models.retain))
    
    models.to.retain <- c(multi.models, social.models.retain)
    
    # these models are retained
    retain.matrix.soc <- social.models.matrix[models.to.retain,]
    
    social.models.matrix <- retain.matrix.soc
  }
  
  # extract the models containing no social network
  
  asocial.models <- rep(NA, length(constraintsVect[,1]))
  
  for (k in 1:length(constraintsVect[,1])){
    sum <- sum(constraintsVect[k,1:num_networks])
    if(sum==0){
      asocial.models[k] <- k
    }
  }
  asocial.models <- as.vector(na.omit(asocial.models))
  
  asocial.models.matrix <- constraintsVect[asocial.models,]
  
  cols.asoc <- which(colnames(constraintsVect)=="asoc")
  
  asocial.retain <- rep(NA, length(asocial.models))
  for (k in 1:length(asocial.models)){
    sum <- sum(asocial.models.matrix[k,which(colnames(constraintsVect)!="asoc")])
    if(sum==0){
      asocial.retain[k] <- k
    }
  }
  
  
  asocial.retain <- as.vector(na.omit(asocial.retain))
  
  asocial.models.to.retain <- asocial.models.matrix[asocial.retain, ]
  asocial.models.to.retain.matrix <- as.matrix(asocial.models.to.retain)
  constraintsVectMatrix <- rbind(social.models.matrix,asocial.models.to.retain)
  
  # add the Null model (without social learning, and no ILVs)
  constraintsVectMatrix <- rbind(constraintsVectMatrix, rep(0, length(constraintsVectMatrix[1,])))
  
  row.names(constraintsVectMatrix) <- NULL
  return(constraintsVectMatrix)
}

# we can run it on any NBDA data object
# it does not matter which NBDA data object we choose - the matrix is the same for all
# we simply specify the number of networks used (2) and the number of ILVs (3)
constraintsVectMatrix <- create.constraints.Vect.Matrix(NBDA_data_object = nbdaData_D2.all, num_networks = 2, num_ILVs = 3)

colnames(constraintsVectMatrix) <- c("foraging network", "neighbour network", "asoc_species",
                                     "asoc_age", "asoc_distance", "soc_species", "soc_age",
                                     "soc_distance")
# we have a look at the ouput
head(constraintsVectMatrix)

# each row represents a model, each column a parameters
# the first two columns refer to the two networks
# columns 3-5 to the ILVs influencing asocial learning
# columns 6-8 to the ILVs influencing social learning
# if a parameter is set to 0, it is not estimated in that model
# if it is a number >0, then the parameter is estimated 
# we could in theory constrain model parameters to be the same 
# by setting equal numbers (e.g. foraging netowrk=1, neighbour network =1)
# but here, we want to estimate all parameters independently (hence, consecutive numbers for each additional parameter)



# 6.6. Run TADA on all females  --------------------------------------------------------------------

# we run TADA with multiple diffusions
TADA.finding.all <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1.all,
      nbdaData_D2.all,
      nbdaData_D3.all,
      nbdaData_D5.all),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
)

# we have a look at the resulting AICc table
# each row corresponds to a model
# Akaike weights are given in the last column
# we can see that the top model is model 187 with an Akaike weight of 0.122
print(TADA.finding.all@printTable)

# we can extract network support via summed Akaike weights
networksSupport(TADA.finding.all)

#       support        numberOfModels
# 0:0 0.13807078              8
# 0:1 0.06135748             64
# 1:0 0.67768236             64
# 1:2 0.12288938             64

# most evidential support for transmission along the foraging network (0.68), 
# followed by asocial models (0.14)
# followed by transmission through both the neighbour and the foraging network (0.12)
# transmission through the neighbour network alone getslittle support (0.06)

variableSupport(TADA.finding.all)

#            s1        s2      ASOC:age_D1   ASOC:species_D1 A SOC:distance_D1 SOCIAL:age_D1   SOCIAL:species_D1  SOCIAL:distance_D1
# support 0.8005717 0.1842469   0.1706535       0.1849043         0.672019      0.272989         0.1887481          0.1766985


# support for distance influencing the asocial learning rate (weight >0.5)
# none of the other ILVs influencing social or asocial learning rate (weight < 0.5)


# extracting effect sizes: model averaged estimates
mle <- modelAverageEstimates(TADA.finding.all , averageType = "median")
mle
# s1                 s2              ASOCIALage_D1     ASOCIALspecies_D1 ASOCIALdistance_D1    SOCIALage_D1   SOCIALspecies_D1  SOCIALdistance_D1 
# 3.2331317          0.0000000          0.0000000          0.0000000         -0.6600527          0.0000000          0.0000000          0.0000000 

# we can see that the value for the distance influencing the asocial learning rate is negative
# which means that as the distance increases, the asocial learning rate decreases

# 6.7. Extract effect sizes ------------------------------------------------

# we extract effect sizes conditional on the best model
# the best model is model 187 (top model in AIC table)
constraintsVectMatrix[187,]
# foraging network neighbour network      asoc_species          asoc_age     asoc_distance       soc_species           soc_age      soc_distance 
# 1                 0                 0                 0                 2                 0                 0                 0
# it contains the foraging network and distance influencing the asocial learning rate

# we create constrained NBDA Data Objects for that specific model
bestModelData1 <- constrainedNBDAdata(nbdadata=nbdaData_D1.all,constraintsVect =constraintsVectMatrix[187,])
bestModelData2 <- constrainedNBDAdata(nbdadata=nbdaData_D2.all,constraintsVect =constraintsVectMatrix[187,])
bestModelData3 <- constrainedNBDAdata(nbdadata=nbdaData_D3.all,constraintsVect =constraintsVectMatrix[187,])
bestModelData5 <- constrainedNBDAdata(nbdadata=nbdaData_D5.all,constraintsVect =constraintsVectMatrix[187,])


# and run TADA on the best model
model.best.social <-
  tadaFit(
    list(
      bestModelData1,
      bestModelData2,
      bestModelData3,
    #  bestModelData4,
      bestModelData5
    )
  )

cbind.data.frame(model.best.social@varNames, model.best.social@outputPar)

# model.best.social@varNames model.best.social@outputPar
# 1            Scale (1/rate):                 160.6239013
# 2    1 Social transmission 1                   3.3275671
# 3     2 Asocial: distance_D1                  -0.6600527


# extract the % of events occurring through social learning
prop.solve.social.byevent <-
  oadaPropSolveByST.byevent(
    nbdadata = list(
      bestModelData1,
      bestModelData2,
      bestModelData3,
  #   bestModelData4,
      bestModelData5
    ),
    model = model.best.social
  )
prop.solve.social.byevent 
# this gives an estimate of the likelihood of each event occurring through social learning

# this extracts the overall percentage that have learned socially
prop.solve.social <-
  oadaPropSolveByST(
    nbdadata = list(
      bestModelData1,
      bestModelData2,
      bestModelData3,
 #     bestModelData4,
      bestModelData5
    ),
    model = model.best.social
  )
prop.solve.social # P=0.4191

# this means that 41.9% of birds have found the dispensers through social learning 
# the remaining 58.1% have done so through asocial learning

# extract profile likelihood. which=1 extracts the first parameter 
# (in this case s for the foraging network)
plotProfLik(which=1,model=model.best.social,range=c(0,20), resolution=10) 
# we check where the profile likelihood crosses the dotted line to get the
# range for the lower and upper interval - set the ranges accordingly
CIs <- profLikCI(which=1,model=model.best.social, lowerRange = c(0,2), upperRange = c(10,20)) # extract confidence intervals
CIs
# Lower CI   Upper CI 
# 0.1522006 15.5464630

# we also want to extract what this means in %

#To get the estimates for the lower bound 
# we have to compute the corresponding values of the other parameters for that model
# if s is constrained to the value of the lower bound 
bestModelDataS1LowerBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D1.all,
  constraintsVect = constraintsVectMatrix[187, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D2.all,
  constraintsVect = constraintsVectMatrix[187, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D3.all,
  constraintsVect = constraintsVectMatrix[187, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D5.all,
  constraintsVect = constraintsVectMatrix[187, ],
  offset = c(CIs[1] , rep(0, 7))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then the value of s at the lower bound is added to s as an offset
bestModelS1LowerBound <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.D1,
      bestModelDataS1LowerBound.D2,
      bestModelDataS1LowerBound.D3,
#      bestModelDataS1LowerBound.D4,
      bestModelDataS1LowerBound.D5
    ) ,
    type = "asocial"
  )
bestModelS1LowerBound@outputPar
# [1] 99.1662098 -0.4286967

#Now we plug this into the prop solve function to get %
prop.solve.social.lower <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound,
    nbdadata = list(
      bestModelDataS1LowerBound.D1,
      bestModelDataS1LowerBound.D2,
      bestModelDataS1LowerBound.D3,
      bestModelDataS1LowerBound.D4,
      bestModelDataS1LowerBound.D5
    )
  )
prop.solve.social.lower
# P(S offset) 
# 0.05132
# lower bound for % of birds having learned the dial task through social learning is 5.1%

# We repeat it for the upper bound
bestModelDataS1upperBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D1.all,
  constraintsVect = constraintsVectMatrix[187, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D2.all,
  constraintsVect = constraintsVectMatrix[187, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D3.all,
  constraintsVect = constraintsVectMatrix[187, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D5.all,
  constraintsVect = constraintsVectMatrix[187, ],
  offset = c(CIs[2] , rep(0, 7))
)


# we again the fit the 'asocial' model with the offset to constrain s to the value of the upper bound
bestModelS1upperBound <-
  tadaFit(
    list(
      bestModelDataS1upperBound.D1,
      bestModelDataS1upperBound.D2,
      bestModelDataS1upperBound.D3,
      #      bestModelDataS1upperBound.D4,
      bestModelDataS1upperBound.D5
    ) ,
    type = "asocial"
  )
bestModelS1upperBound@outputPar
# [1] 407.109918  -1.188658
#Now plug into the prop solve function
prop.solve.social.upper <-
  oadaPropSolveByST(
    model = bestModelS1upperBound,
    nbdadata = list(
      bestModelDataS1upperBound.D1,
      bestModelDataS1upperBound.D2,
      bestModelDataS1upperBound.D3,
      bestModelDataS1upperBound.D4,
      bestModelDataS1upperBound.D5
    )
  )
prop.solve.social.upper
# P(S offset) 
# 0.63167 
# upper bound for % of birds having learned about the location of the dispensers through social learning is 63.2%

# Finally, we extract the effect size for the distance influencing asocial learning
# again, we base this off the best performing model
# all we have to do is set which=2 (for the second parameter in the model)
model.best.social@varNames
# [1] "Scale (1/rate):"         "1 Social transmission 1" "2 Asocial: distance_D1" 
plotProfLik(which=2,model=model.best.social,range=c(-2,1), resolution=10) 
# we check where the profile likelihood crosses the dotted line to get the
# range for the lower and upper interval - set the ranges accordingly
CIs <- profLikCI(which=2,model=model.best.social, lowerRange = c(-2,-1), upperRange = c(-0.5,0.5)) # extract confidence intervals
CIs

# Lower CI    Upper CI 
# -1.44397740 -0.06182192

# all we need to do now is back-transform the estimates

exp(c(mle[5], CIs[1], CIs[2]))
# the unit here is standard deviation - as the distance was standardized for each dispenser area separately
# the standard deviations are slightly different in each area - we therefore retain the unit as 'standard deviation'

# ASOCIALdistance_D1           Lower CI           Upper CI 
# 0.5168241                   0.2359873          0.9400503 


# 7) Wool choice ----------------------------------------------------------
# here, we investigate whether birds showed a colour preference for the colour initially provided (that the demonstrators were restricted to)

wool.choice <- read.delim("Wool choice.txt", sep="\t")

wool.choice.learners <- subset(wool.choice, wool.choice$Demos!="yes")
wool.choice.learners
# this df shows the picked colour (first_colour), as well as the initially provided colour (initial.col)
# it only contains birds that had a choice (excluding demonstrators)
# it also shows whether birds have incorporate a second colour and which dispenser area they belong to

# We conduct a Fisher's exact test to see whether there is non-random clustering within each dispenser are
# we expect a preference for the first introduced colour if birds use social information

fisher <- fisher.test(wool.choice.learners$first_color, wool.choice.learners$Initial_col_provided, alternative = "greater")
fisher


# Fisher's Exact Test for Count Data

# data:  wool.choice.learners$first_color and wool.choice.learners$Initial_col_provided
# p-value = 0.02498
# alternative hypothesis: greater

# to investigate influence of age and species, we subset it to know individuals
wool.choice.learners.glm <- subset(wool.choice.learners, wool.choice.learners$PIT!=0)
glm.wool.choice <-
  glm(
    wool.choice.learners.glm$Matched ~ wool.choice.learners.glm$Age + wool.choice.learners.glm$Species,
    family = binomial(link = "logit")
  )
summary(glm.wool.choice)

Call:
  glm(formula = wool.choice.learners.glm$Matched ~ wool.choice.learners.glm$Age + 
        wool.choice.learners.glm$Species, family = binomial(link = "logit"))
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.5829  -1.0108   0.8203   0.8203   1.3537  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                             -17.566   3956.180  -0.004    0.996
# wool.choice.learners.glm$Agefirst.year   -1.322      1.238  -1.067    0.286
# wool.choice.learners.glm$SpeciesGRETI    18.482   3956.180   0.005    0.996
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 17.945  on 12  degrees of freedom
# Residual deviance: 15.106  on 10  degrees of freedom
# AIC: 21.106
# 
# Number of Fisher Scoring iterations: 16


# 8) Visualization --------------------------------------------------------

# create a network with polgygons around dispenser areas

library(igraph)

# for each PIT tag, extract whether it was seen at a network feeder
col.plot <- NULL

for(i in rownames(foraging.network.NBDA)){
  box <- subset(ILVs.combined$Box, ILVs.combined$PIT_f==i)
  if(length(box)>1){
    box <- subset(box, box %in% wool.choice$Box)
  } 
  if(length(box)==0){
    col <- "none"
    box <- "no"
  }
    if(box=="no_box"){
    wool.used <- "unknown"
  } else {
    wool.used <- subset(wool.choice$first_color, wool.choice$Box==box)  
    
  }
  if(length(wool.used)==0){
    wool.used <- "none"
  }
    col.plot[which(rownames(foraging.network.NBDA)==i)] <- wool.used
}

# reassign the correct hex codes for each colour
col.plot[col.plot=="unknown"] <- "#9a9999" # grey
col.plot[col.plot=="none"] <- "#faf7f7" # white
col.plot[col.plot=="Pi"] <- "#fb0edc" # pink
col.plot[col.plot=="O"] <- "#fdb633" # orange 
col.plot[col.plot=="Pu"] <- "#b87be3" # purple
col.plot[col.plot=="B"] <- "#33CCFF" # blue



g.net <- graph_from_adjacency_matrix(foraging.network.NBDA, mode = "undirected",
                                     weighted = TRUE, diag = TRUE)

E(g.net)$width <- E(g.net)$weight

V(g.net)$colour <- col.plot

# # assign to which disepnser area(s) each bird belongs
# # for those not breeding in boxes, we assign the dispenser they visited
# list.D <- NULL
# for(i in c("D1", "D2", "D3","D4", "D5")){
#   # We extract which females were breeding in boxes within 200m of the respective dispenser or have visited it
#   IDs.included <- subset(ILVs.combined$PIT_f, ILVs.combined[,i]<=200 | ILVs.combined[,paste(i, ".visited", sep="")]==1 ) 
#   IDs.included <- subset(IDs.included, IDs.included %in% IDs.to.include.in.NBDA) # cut down to IDs in the foraging network
#   D <- which(rownames(foraging.network.NBDA) %in% IDs.included)
#   list.D[[i]] <- D
# }
# 
# list.D
 

# set transparent polgygons
# col.adj <- grDevices::adjustcolor(c("#a08f00","#62dab9", "#b738bd", "#7a3f63", "#c31910"), alpha=0.15)


png( "network.png", units="in", width=12, height=4, res=400)

set.seed(4)
igraph::plot.igraph( g.net,
      vertex.size = 8,
      edge.curved = 0.25,
      edge.color =  "#8c8989",
      vertex.color = V(g.net)$colour,
      vertex.label = NA,
      vertex.frame.colour = "black",
      edge.width = E(g.net)$width*4,
      frame = FALSE,
      layout=layout.graphopt(g.net),
      asp = 1,
 #     mark.groups = list.D,
  #    mark.border =NA,
  #  mark.col=col.adj
      
      
)



dev.off()


