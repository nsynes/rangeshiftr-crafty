library(rgdal)
library(raster)
library(RangeshiftR)
library(sf)

wd <- 'C:/dev/ForestResearch/OPM'
dirData <- file.path(wd, 'for-RangeshiftR')

rasterizeRes <- 2 # resolution at which to rasterise the habitat polygons
habitatRes <- 200 # Habitat resolution for rangeshifter

# Need to create these folders for rangeshiftr to work
dirRsftr <- file.path(wd, 'test-example')
dirRsftrInput <- file.path(dirRsftr,"Inputs")
dirRsftrOuput <- file.path(dirRsftr,"Outputs")
dirRsftrOutputMaps <- file.path(dirRsftr,"Output_Maps")
dir.create(dirRsftr)
dir.create(dirRsftrInput)
dir.create(dirRsftrOuput)
dir.create(dirRsftrOutputMaps)


######################################
######################################
# GET HABITAT DATA AS A RASTER
######################################
######################################

shpHabitat <- st_read(dirData, layer = "OSMM_gspaceAOI")

# Create a unique number of each 'priFunc' category - we'll use this to define habitat types.
shpHabitat$priFuncN <- as.numeric(shpHabitat$priFunc)

# Slow process tp rasterise at fine resolution, so check if file already exists
if ( !file.exists(file.path(dirRsftrInput, sprintf('Habitat-%sm.tif', rasterizeRes)))) {
  rstHabitat <- rasterize(x=shpHabitat,
                          y=raster(extent(shpHabitat), res=rasterizeRes),
                          field='priFuncN',
                          fun=min, # where there are overlapping polygons, use lowest value (is there a better way to choose?)
                          update=TRUE)
  writeRaster(rstHabitat, file.path(dirRsftrInput, sprintf('Habitat-%sm.tif', rasterizeRes)))
}

# Aggregate the raster pixels - good to start from fine res raster, but don't need to model at such fine resolution.
rstHabitat <- raster(file.path(dirRsftrInput, sprintf('Habitat-%sm.tif', rasterizeRes)))
rstHabitat <- aggregate(rstHabitat, fact=habitatRes/rasterizeRes, fun=min)

# Export as ascii file for RangeShifter.
writeRaster(rstHabitat, file.path(dirRsftrInput, sprintf('Habitat-%sm.asc', habitatRes)), format="ascii", overwrite=TRUE)


# Create a dataframe of habitat codes and names, then add a quality column. This defines the quality (in terms of carrying capacity)
# for the species, and should range from 0-100%.
# The values I have used are just made up to create an example.
dfHabitatTypes <- as.data.frame(shpHabitat[c('priFuncN','priFunc')])
dfHabitatTypes$geometry <- NULL
dfHabitatTypes <- dfHabitatTypes[!duplicated(dfHabitatTypes), ]
dfHabitatTypes$quality <- 50
dfHabitatTypes[dfHabitatTypes$priFunc %in% c('Natural', 'Private Garden', 'Public Park Or Garden', 'Cemetery', 'Allotments Or Community Growing Spaces'),]$quality <- 100
dfHabitatTypes[dfHabitatTypes$priFunc %in% c('Amenity - Residential Or Business', 'Land Use Changing'),]$quality <- 10

# Reclassify the habitat type raster using the quality values. We run RangeShiftR with the habitat quality raster.
rstHabitatQuality <- reclassify(rstHabitat, dfHabitatTypes[c('priFuncN','quality')])
writeRaster(rstHabitatQuality, file.path(dirRsftrInput, sprintf('HabitatQuality-%sm.asc', habitatRes)), format="ascii", overwrite=TRUE)

######################################
######################################





######################################
######################################
# GET SPECIES LOCATIONS AS SHAPEFILE
######################################
######################################

csvSpecies <- read.csv(file=file.path(dirData, 'opm_sites.csv'), header=TRUE, sep=",")

# Subset to only sites that were infested in 2012
csvSpecies <- subset(csvSpecies, Status_2012 == 'Infested')

# Convert to shapefile
csvSpecies <- csvSpecies[!is.na(csvSpecies$Easting), ]
csvSpecies <- csvSpecies[!is.na(csvSpecies$Northing), ]
shpInitialIndividuals <- csvSpecies
csvSpecies <- NULL
coordinates(shpInitialIndividuals) <- ~Easting+Northing
projection(shpInitialIndividuals) <- crs(shpHabitat)
shapefile(shpInitialIndividuals, file.path(dirData, "opm_sites.shp"), overwrite=TRUE)

# Crop to species locations only within the study landscape
shpInitialIndividuals <- crop(shpInitialIndividuals, extent(shpHabitat))
shpInitialIndividuals$n <- 1

# Need to rasterise then extract the species locations to get the xy, row/col (not spatial) indices for rangeshifter.
dfInitialIndividuals <- extract(rasterize(shpInitialIndividuals, rstHabitat, field='n', background=0), shpInitialIndividuals, cellnumbers=T, df=TRUE)

# RangeShiftR requires a specific format for the individuals file, so add the required columns here,
# and convert 'cells' value to x/y, row/col values.
# As an example we are just initialising each cell with 100 individuals - we may need to adjust this later.
dfInitialIndividuals$Year <- 0
dfInitialIndividuals$Species <- 0
dfInitialIndividuals$X <- dfInitialIndividuals$cells %% ncol(rstHabitat)
dfInitialIndividuals$Y <- nrow(rstHabitat) - (floor(dfInitialIndividuals$cells / ncol(rstHabitat)))
dfInitialIndividuals$Ninds <- 100
dfInitialIndividuals <- dfInitialIndividuals[ , !(names(dfInitialIndividuals) %in% c('ID', 'cells', 'layer'))]

write.table(dfInitialIndividuals, file.path(dirRsftrInput, 'initial_inds.txt'), row.names = F, quote = F, sep = '\t')

######################################




######################################
######################################
# LOAD LONDON BOROUGH SHAPEFILE TO USE AS PROXY FOR CRAFTY MANAGEMENT STRATEGIES
######################################
######################################

# Load the boroughs shapefile and make sure it has the same coordinate system as the habitat data.
shpBoroughs <- st_read(dirData, layer = "case_studies")
shpBoroughs <- st_transform(shpBoroughs, crs(shpHabitat))

######################################
######################################








######################################
######################################
# RANGESHIFTER PARAMETER SETUP
######################################
######################################

# We need to simulate at least 2 years of RangeShifR for each CRAFTY iteration.
rangeshiftrYears <- 2

sim <- Simulation(Simulation = 1,
                  Years = rangeshiftrYears,
                  Replicates = 1,
                  OutIntPop = 1, # interval for output of population data
                  #OutIntInd = 1, # interval for output of individual data
                  ReturnPopRaster=TRUE) # We need RangeShiftR to return a raster with the population so we can use it for CRAFTY

land <- ImportedLandscape(LandscapeFile=sprintf('HabitatQuality-%sm.asc', habitatRes),
                          Resolution=habitatRes,
                          HabitatQuality=TRUE,
                          K=100) # carrying capacity (individuals per hectare) when habitat at 100% quality

# We have often used ReproductionType=0 for invertebrate species in RangeShifter.
# It doesn't mean the species is asexual, just that we are modelling only the females - assumming that the species mates
# upon emergence in the natal patch, then fertilised females disperse.
# Rmax (intrinsic growth rate) - we don't have a good basis for this value yet.
demo <- Demography(Rmax = 10,
                   ReproductionType = 0) # 0 = asexual / only female; 1 = simple sexual; 2 = sexual model with explicit mating system

# Basing dispersal kernel distance on the alpha value given in Cowley et al 2015. Would be good to test sensitivty around this value
# and to find other sources for dispersal distance estimates.
disp <-  Dispersal(Emigration = Emigration(EmigProb = 0.2), 
                   Transfer   = DispersalKernel(Distances = 800),
                   Settlement = Settlement() )

# InitType=2 tells RangeShiftR that we will provide a txt file of individuals to initialise from.
init <- Initialise(InitType=2,
                   InitIndsFile='initial_inds.txt')

# Setup the simulation with the parameters defined above.
s <- RSsim(simul = sim, land = land, demog = demo, dispersal = disp, init = init)
validateRSparams(s)

######################################
######################################





######################################
######################################
# MAIN MODEL LOOP
######################################
######################################

# Create empty data frame and raster stack to store the output data
dfRangeShiftrData <- data.frame()
outRasterStack <- stack()

# we run RangeShiftR and CRAFTY once per iteration.
for (iteration in 1:20) {
  
  # Set up RangeShiftR for current iteration
  sim <- Simulation(Simulation = iteration,
                    Years = rangeshiftrYears,
                    Replicates = 1,
                    #OutIntPop = 1, # interval for output of population data
                    #OutIntInd = 1, # interval for output of individual data
                    ReturnPopRaster=TRUE)
  s <- RSsim(simul = sim, land = land, demog = demo, dispersal = disp, init = init)
  
  # Run RangeShiftR. Use result to store our output population raster.
  result <- RunRS(s, sprintf('%s/',dirRsftr))
  
  # COMMENTED OUT - we don't need the population file for the current setup, but may be useful to check results later.
  #dfPop <- read.table(file.path(dirRsftrOuput, sprintf('Batch1_Sim%s_Land1_Pop.txt', iteration)), header=TRUE)
  #dfPop <- subset(dfPop, Year == rangeshiftrYears-1)
  
  # Store RangeShiftR's population raster in our output stack.
  # Set the coordinate reference system for RangshiftR's output population raster.
  crs(result) <- crs(rstHabitat)
  extent(result) <- extent(rstHabitat)
  outRasterStack <- addLayer(outRasterStack, result[[rangeshiftrYears]])
  
  # Store RangeShiftR's population data in our output data frame.
  dfRange <- readRange(s, sprintf('%s/',dirRsftr))
  dfRange$iteration <- iteration
  dfRangeShiftrData <- rbind(dfRangeShiftrData, dfRange[1,])
  
  # Extract the population raster to a shapefile of the individuals
  shpIndividuals <- rasterToPoints(result[[rangeshiftrYears]], fun=function(x){x > 0}, spatial=TRUE) %>%st_as_sf()
  shpIndividuals <- st_transform(shpIndividuals, crs(shpHabitat))
  shpIndividuals$id <- 1:nrow(shpIndividuals)
  
  
  
  #############
  # RUN CRAFTY HERE...
  # Could pass the current population data to CRAFTY as updated capitals for each borough, or whichever management
  # structure we end up using.
  #############
  
  
  
  #############
  # Demonstration of how an action might remove OPM in certain locations:
  # This is basically simulating how you might use management actions in different boroughs to change the
  # habitat quality in the next iteration of RangeShiftR.
  #############
  
  # Randomly select a borough.
  randomBorough <- sample(shpBoroughs$objectid, 1)
  
  # COMMENTED OUT - could use the management action to remove individuals from the simulation:
  #shpIndividuals <- subset(st_intersection(shpIndividuals, shpBoroughs), objectid != randomBorough)
  
  # We can use the manamgement action to adjust the habitat quality of specific locations. Here
  # we are just changing habitat quality based on the randomly selected borough, so the selected
  # borough has all habitat qualities reduced by half for the next iteration.
  selectedBorough <- mask(rstHabitatQuality, subset(shpBoroughs, objectid == randomBorough))
  otherBoroughs <- mask(rstHabitatQuality, subset(shpBoroughs, objectid != randomBorough))
  writeRaster(merge(selectedBorough * 0.5, otherBoroughs), file.path(dirRsftrInput, sprintf('HabitatQuality-%sm.asc', habitatRes)), format="ascii", overwrite=TRUE)
 
  #############
  #############
  
  
  #############
  # If we have manually removed any individuals/populations, we need to create a new individuals table to pass
  # to the next iteration of RangeShiftR.
  #############
  
  dfNewIndsTable <- extract(rasterize(shpIndividuals, rstHabitat, field=sprintf('rep0_year%s', rangeshiftrYears-1)), shpIndividuals, cellnumbers=T, df=TRUE)
  
  dfNewIndsTable$Year <- 0
  dfNewIndsTable$Species <- 0
  dfNewIndsTable$X <- dfNewIndsTable$cells %% ncol(rstHabitat)
  dfNewIndsTable$Y <- nrow(rstHabitat) - (floor(dfNewIndsTable$cells / ncol(rstHabitat)))
  dfNewIndsTable$Ninds <- dfNewIndsTable$layer
  dfNewIndsTable <- dfNewIndsTable[ , !(names(dfNewIndsTable) %in% c('ID', 'cells', 'layer'))]
  dfNewIndsTable <- dfNewIndsTable[!is.na(dfNewIndsTable$Ninds),]
  
  write.table(dfNewIndsTable, file.path(dirRsftrInput, sprintf('inds%s.txt', iteration)),row.names = F, quote = F, sep = '\t')
  
  # Update the initialisation parameters for RangeShiftR to point to the new individuals file.
  init <- Initialise(InitType=2, InitIndsFile=sprintf('inds%s.txt', iteration))
  #############
  #############
  
}

######################################
######################################



######################################
# Export the final raster stack and view the output dataframe
######################################
writeRaster(outRasterStack, file.path(dirRsftrOutputMaps, 'final-result.tif'), overwrite=TRUE)

plot(outRasterStack)

View(dfRangeShiftrData)


######################################
######################################


