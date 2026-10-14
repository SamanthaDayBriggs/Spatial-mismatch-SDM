library(raster); library(sp); library(megaSDM); library(sdm); library(geodata); 
library(ENMTools);library(dismo); library(earth); library(terra);
library(gbm); library(kernlab); library(rJava); library(parallel)

setwd("~/Desktop/Data")

########################################################################
############################ SDMs ######################################
########################################################################

# Get recent historical climate data and elevation from WorldClim
bio_cur <- worldclim_global(var='bio', res=0.5, download=TRUE, 
                            path="~/Desktop/Data/climate")
elev <- worldclim_global(var='elev', res=0.5, download=TRUE, 
                         path="~/Desktop/Data/climate")
# Combine climate and elevation
bio_curr <- c(bio_cur,elev)
# Simplify variable names
names(bio_curr) = c(paste0("bio", 1:19), "elev")


#Repeat data collection, add elevation, and renaming for other scenarios

fut_2050_245 <- geodata::cmip6_world("GISS-E2-1-G", "245", var="bioc", 
                    "2041-2060", res=0.5, path="~/Desktop/Data/climate")
fut_2050_585 <- geodata::cmip6_world("GISS-E2-1-G", "585", var="bioc", 
                    "2041-2060", res=0.5, path="~/Desktop/Data/climate")
fut_2070_245 <- geodata::cmip6_world("GISS-E2-1-G", "245", var="bioc", 
                    "2061-2080", res=0.5, path="~/Desktop/Data/climate")
fut_2070_585 <- geodata::cmip6_world("GISS-E2-1-G", "585", var="bioc", 
                    "2061-2080", res=0.5, path="~/Desktop/Data/climate")
fut_2090_245 <- geodata::cmip6_world("GISS-E2-1-G", "245", var="bioc", 
                    "2081-2100", res=0.5, path="~/Desktop/Data/climate")
fut_2090_585 <- geodata::cmip6_world("GISS-E2-1-G", "585", var="bioc", 
                    "2081-2100", res=0.5, path="~/Desktop/Data/climate")

biofut_2050_245 <- c(fut_2050_245, elev)
names(biofut_2050_245) = c(paste0("bio", 1:19), "elev")
biofut_2050_585 <- c(fut_2050_585, elev)
names(biofut_2050_585) = c(paste0("bio", 1:19), "elev")
biofut_2070_245 <- c(fut_2070_245, elev)
names(biofut_2070_245) = c(paste0("bio", 1:19), "elev")
biofut_2070_585 <- c(fut_2070_585, elev)
names(biofut_2070_585) = c(paste0("bio", 1:19), "elev")
biofut_2090_245 <- c(fut_2090_245, elev)
names(biofut_2090_245) = c(paste0("bio", 1:19), "elev")
biofut_2090_585 <- c(fut_2090_585, elev)
names(biofut_2090_585) = c(paste0("bio", 1:19), "elev")

# Confirm that rasters are stacked
bio_cur2 <- stack(bio_curr)
bio_fut.50.585s <- stack(biofut_2050_585)
bio_fut.50.245s <- stack(biofut_2050_245)
bio_fut.70.585s <- stack(biofut_2070_585)
bio_fut.70.245s <- stack(biofut_2070_245)
bio_fut.90.585s <- stack(biofut_2090_585)
bio_fut.90.245s <- stack(biofut_2090_245)

# Climate data is currently world wide, crop down to the extent of North America
# spatial extent of the study
max.lat <- 75
min.lat <- 10
max.lon <- -50
min.lon <- -175
study <- extent(min.lon, max.lon, min.lat, max.lat)

# crop climate data to the spatial extent of the study
bio_cur2 <- crop(bio_cur2, study)
bio_fut.50.245 <- crop(bio_fut.50.245s, study)
bio_fut.50.585 <- crop(bio_fut.50.585s, study)
bio_fut.70.245 <- crop(bio_fut.70.245s, study)
bio_fut.70.585 <- crop(bio_fut.70.585s, study)
bio_fut.90.245 <- crop(bio_fut.90.245s, study)
bio_fut.90.585 <- crop(bio_fut.90.585s, study)


##### Species specific from here on, modify for each species #################
### For Iris missouriensis ###################################################

# Read in CSV downloaded from GBIF.org
occ_Imiss <- read.csv("~/Data/NorthAmerica_occurrences/Iris_missouriensis.csv")

# Generate buffers
oc_Imiss<-occ_Imiss [c("decimalLongitude", "decimalLatitude")]
coordinates(oc_Imiss) <- ~decimalLongitude+decimalLatitude
buf_Imiss<-raster::buffer(oc_Imiss, width=100000)
### Generate background points
BackgroundPoints(spplist=c('Imiss'), envdata=bio_cur2, output='data', 
                 nbg= nrow(occ_Imiss), 
                 method='Varela', spatial_weights = 0.5, 
                 buffers=buf_Imiss, ncores=15)
bg.points.Imiss <- read.csv('~/Data/data/Imiss_background.csv')
# Only keep the lat/longs
bg.points.Imiss <- as.data.frame(bg.points.Imiss)[,c('x','y')]
colnames(bg.points.Imiss) <- c('longitude', 'latitude')

# Combine into one dataframe with binary presence/absence
oc_Imiss <- occ_Imiss[c("decimalLongitude", "decimalLatitude")]
oc_Imiss$presence <- 1
bg.points.Imiss$presence <- 0
colnames(oc_Imiss) <- colnames(bg.points.Imiss)
all_points_Imiss <- rbind(oc_Imiss, bg.points.Imiss)
# Convert to a SpatialPointsDataFrame
species_Imiss <- SpatialPointsDataFrame(coords = all_points_Imiss
                                        [, c("longitude", "latitude")],
                                       data = all_points_Imiss,
                                       proj4string = CRS
                                       ("+proj=longlat +datum=WGS84"))
bio_cur2<-stack(bio_curr)
e <- as(extent(bio_cur2), 'SpatialPolygons')
crs(e) <- crs(species_Imiss)

# create the sdmData for the sdm function
dat_Imiss<- sdmData(train = species_Imiss, predictors = bio_cur2, 
                    parallelSettings=list(ncore=15, method='parallel'))

#Fit the models
percent <- 30
rep <- 3

# Begin with the full model with all variables
formula_all<- presence ~ bio1+bio2+bio3+bio4+bio5+bio6+bio7+bio8+bio9+bio10+
  bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19+elev
Mod.all.Imiss<-sdm(formula = formula_all, data = dat_Imiss, methods = 
                     c("glm", "maxlike", "mars","maxent"), 
                  replication = 'sub', test.percent = percent, n=rep,
                  parallelSettings=list(ncore=15, method='parallel'))

# Check the importance of each variable
viImiss <- getVarImp(Mod.all.Imiss, id=1:rep)
plot(viImiss, main= "viImiss")
# Save the VI plot
Destination="~/Data/Output/ImissVI.pdf"
pdf(file=Destination)
plot(viImiss, main = "Bioclim all")
dev.off()

# Generate a correlation matrix
raster.cor.matrix(bio_curr, method = "pearson")

# From the full model with generated VI plot and correlation matrix, 
#determine appropriate variables 
formula_Imiss <- presence ~ bio1 + bio2 + bio7 + elev
mod.Imiss <- sdm(formula = formula_Imiss, data = dat_Imiss, methods = 
                   c("glm", "maxlike", "mars","maxent"), 
              replication = 'sub', test.percent = percent, n=rep,
              parallelSettings=list(ncore=15, method='parallel'))

# Evaluate the fit of the model
modeval.Imiss <- getEvaluation(mod.Imiss)
write.csv(modeval.Imiss,file="~/Data/Output/modeval.Imiss.csv", 
          row.names = TRUE)

# Create an ensemble prediction based on current climate
mod.Imiss.ensemble <- ensemble(mod.Imiss,
                              newdata = bio_cur2,
                              setting = list(method ='weighted', stat = 'TSS'),
                              parallelSettings=list(ncore=15, method='parallel'))

# Create ensemble models for future climates and time scenarios
mod.Imiss.2090.585 <- ensemble(mod.Imiss,
                              newdata = bio_fut.90.585s,
                              setting = list(method ='weighted', stat = 'TSS'),
                              parallelSettings=list(ncore=15, method='parallel'))
mod.Imiss.2090.245 <- ensemble(mod.Imiss,
                              newdata = bio_fut.90.245s,
                              setting = list(method ='weighted', stat = 'TSS'),
                              parallelSettings=list(ncore=15, method='parallel'))
mod.Imiss.2070.585 <- ensemble(mod.Imiss,
                              newdata = bio_fut.70.585s,
                              setting = list(method ='weighted', stat = 'TSS'),
                              parallelSettings=list(ncore=15, method='parallel'))
mod.Imiss.2070.245 <- ensemble(mod.Imiss,
                              newdata = bio_fut.70.245s,
                              setting = list(method ='weighted', stat = 'TSS'),
                              parallelSettings=list(ncore=15, method='parallel'))
mod.Imiss.2050.585 <- ensemble(mod.Imiss,
                              newdata = bio_fut.50.585s,
                              setting = list(method ='weighted', stat = 'TSS'),
                              parallelSettings=list(ncore=15, method='parallel'))
mod.Imiss.2050.245 <- ensemble(mod.Imiss,
                              newdata = bio_fut.50.245s,
                              setting = list(method ='weighted', stat = 'TSS'),
                              parallelSettings=list(ncore=15, method='parallel'))


# Save ensembles as pdfs and rasters

Destination="~/Data/Output/ImissRplot.pdf"
pdf(file=Destination)
plot(mod.Imiss.ensemble, main = "Bioclim Contemporary Ensemble")
plot(mod.Imiss.2090.585, main = "SSP585 2081-2100 Ensemble")
plot(mod.Imiss.2090.245, main = "SSP245 2081-2100 Ensemble")
plot(mod.Imiss.2070.585, main = "SSP585 2061-2080 Ensemble")
plot(mod.Imiss.2070.245, main = "SSP245 2061-2080 Ensemble")
plot(mod.Imiss.2050.585, main = "SSP585 2041-2060 Ensemble")
plot(mod.Imiss.2050.245, main = "SSP245 2041-2060 Ensemble")
dev.off()

writeRaster(mod.Imiss.2050.245, filename="~/Data/Output/Imiss50245.tif", 
            overwrite=TRUE)
writeRaster(mod.Imiss.2050.585, filename="~/Data/Output/Imiss50585.tif", 
            overwrite=TRUE)
writeRaster(mod.Imiss.2070.245, filename="~/Data/Output/Imiss70245.tif", 
            overwrite=TRUE)
writeRaster(mod.Imiss.2070.585, filename="~/Data/Output/Imiss70585.tif", 
            overwrite=TRUE)
writeRaster(mod.Imiss.2090.245, filename="~/Data/Output/Imiss90245.tif", 
            overwrite=TRUE)
writeRaster(mod.Imiss.2090.585, filename="~/Data/Output/Imiss90585.tif", 
            overwrite=TRUE)
writeRaster(mod.Imiss.ensemble, filename="~/Data/Output/Imisscontemp.tif", 
            overwrite=TRUE)

# Generate Binary models to overlay

# Pull occurrences and absences
presence_predictions <- extract(mod.Imiss.ensemble, oc_Imiss
                                [, c("longitude", "latitude")])
background_predictions <- extract(mod.Imiss.ensemble, bg.points.Imiss
                                  [, c("longitude", "latitude")])
Imissbin <- mod.Imiss.ensemble

# Calculate the threshold for occurrence
ev <- evaluate(p = presence_predictions, a = background_predictions)
threshold <- threshold(ev, stat = "spec_sens")
print(threshold)

# Use the calculated threshold to binarize the models
Imissbin[] <- ifelse(Imissbin[] >= threshold, 1, 0) 
Imissbin90245<-mod.Imiss.2090.245
Imissbin90245[] <- ifelse(Imissbin90245[] >= threshold, 1, 0)
Imissbin90585<-mod.Imiss.2090.585
Imissbin90585[] <- ifelse(Imissbin90585[] >= threshold, 1, 0)
Imissbin70245<-mod.Imiss.2070.245
Imissbin70245[] <- ifelse(Imissbin70245[] >= threshold, 1, 0)
Imissbin70585<-mod.Imiss.2070.585
Imissbin70585[] <- ifelse(Imissbin70585[] >= threshold, 1, 0)
Imissbin50245<-mod.Imiss.2050.245
Imissbin50245[] <- ifelse(Imissbin50245[] >= threshold, 1, 0)
Imissbin50585<-mod.Imiss.2050.585
Imissbin50585[] <- ifelse(Imissbin50585[] >= threshold, 1, 0)

#Save the binary rasters
writeRaster(Imissbin50245, filename="~/Data/Output/Imissbin50245.tif", 
            overwrite=TRUE)
writeRaster(Imissbin50585, filename="~/Data/Output/Imissbin50585.tif", 
            overwrite=TRUE)
writeRaster(Imissbin70245, filename="~/Data/Output/Imissbin70245.tif", 
            overwrite=TRUE)
writeRaster(Imissbin70585, filename="~/Data/Output/Imissbin70585.tif", 
            overwrite=TRUE)
writeRaster(Imissbin90245, filename="~/Data/Output/Imissbin90245.tif", 
            overwrite=TRUE)
writeRaster(Imissbin90585, filename="~/Data/Output/Imissbin90585.tif", 
            overwrite=TRUE)
writeRaster(Imissbin, filename="~/Data/Output/Imissbin.tif", 
            overwrite=TRUE)

#############################################################################
##### Overlay binary rasters ################################################
#############################################################################

######### Recent historical raster overlay only #############################
# Repeat for each climate/time scenario

setwd("~/Data/Output")
raster_folder <- "~/Data/Output"

# First, identify the "base" raster that everything will be layered onto
Imiss <- rast("~/Data/Output/Imissbin.tif") 

# List all binary raster files in the folder
binary_raster_files <- list.files(path=raster_folder, pattern="bin.tif$", 
                                full.names=TRUE)

# Create a vector to store overlap areas
overlap_areas <- numeric(length(binary_raster_files))

# Create a loop that layers all of the raster files onto the "base" file
for (i in seq_along(binary_raster_files)) {
  # Load the binary rasters
  binary_raster <- rast(binary_raster_files[i])
  # Calculate the overlap of the binary SDMs when both are present (=1)
  overlap <- binary_raster == 1 & Imiss == 1
  # Calculate the area of overlap in km 
  cell_areas <- cellSize(binary_raster, unit= "km")
  overlap_area <- sum(cell_areas[values(overlap)==1], na.rm=TRUE)
  # Store the overlap areas in the vector created previously
  overlap_areas[i] <- overlap_area
}

# Generate a data frame for the results and save them
resultscont <- data.frame(Raster_File=basename(binary_raster_files),
  Overlap_con_Area_km2=overlap_areas)

write.csv(resultscont, file= "~/Data/Output/Overlap/overlap_areas_con.csv", 
          row.names=FALSE)

########################################################################
################### Elevation range shifts #############################
########################################################################

# Set the extent for regions of interest
# Here we used a northern, mid (Colorado), and southern range extent
# Northern extent 
extent_no <- c(-110, -100, 44.5, 49.5)
# Colorado extent
extent_co <- c(-110, -100, 36.5, 41.5)
# Southern extent
extent_so <- c(-110, -100, 28.5, 33.5)

#Northern extent elevational shifts (change out extents for different regions)

# Load in Iris recent historical and end of century binaries along with elevation
Iris_con_raster <- rast("~/Desktop/Data/Datanew/Imissbin.tif")
Iris90245_raster <- rast("~/Desktop/Data/Datanew/Imissbin90245.tif")
Iris90585_raster <- rast("~/Desktop/Data/Datanew/Imissbin90585.tif")
elev_raster <- rast("~/Desktop/Data/climate/wc2.1_30s/wc2_elev.tif")

# Crop the rasters to the spatial extent of interest
Iris_con_no <- crop(Iris_con_raster, extent_no)
Iris_90245_no <- crop(Iris90245_raster, extent_no)
Iris_90585_no <- crop(Iris90585_raster, extent_no)
elev_no <- crop(elev_raster, extent_no)

# Mask elevation so that elevation is 0 where the binary raster is 0
elev_at_Irisbinno <- mask(elev_no , Iris_con_no , maskvalues=0)
elev_at_Irisbin90245no <- mask(elev_no , Iris_90245_no , maskvalues=0)
elev_at_Irisbin90585no <- mask(elev_no , Iris_90585_no , maskvalues=0)

# Get elevational range of occurrence points
Irisbinno_elev_range <- global(elev_at_Irisbinno , range, na.rm=TRUE)
Irisbin90245no_elev_range <- global(elev_at_Irisbin90245no , range, na.rm=TRUE)
Irisbin90585no_elev_range <- global(elev_at_Irisbin90585no , range, na.rm=TRUE)

############# Repeat for every species ##################################


# Combine data from recent historical elevational ranges
contno_list<- list(Irisbinno_elev_range,Bappbinno_elev_range, Bflabinno_elev_range, Bgribinno_elev_range, Bhunbinno_elev_range,
                   Bmelbinno_elev_range, Bnevbinno_elev_range, Binsbinno_elev_range, Hthebinno_elev_range, Mtribinno_elev_range, 
                   Prutbinno_elev_range, Splabinno_elev_range)

# Combine data from end of century SSP2-4.5 elevational ranges
bin90245no_list<- list(Irisbin90245no_elev_range, Bappbin90245no_elev_range, Bflabin90245no_elev_range, Bgribin90245no_elev_range, Bhunbin90245no_elev_range,
                       Bmelbin90245no_elev_range, Bnevbin90245no_elev_range, Binsbin90245no_elev_range, Hthebin90245no_elev_range, Mtribin90245no_elev_range, 
                       Prutbin90245no_elev_range, Splabin90245no_elev_range)

# Combine data from end of century SSP5-8.5 elevational ranges
bin90585no_list<- list(Irisbin90585no_elev_range, Bappbin90585no_elev_range, Bflabin90585no_elev_range, Bgribin90585no_elev_range, Bhunbin90585no_elev_range,
                       Bmelbin90585no_elev_range, Bnevbin90585no_elev_range, Binsbin90585no_elev_range, Hthebin90585no_elev_range, Mtribin90585no_elev_range, 
                       Prutbin90585no_elev_range, Splabin90585no_elev_range)

# Combine all to create a single dataframe and write the csv
elevrangesno <- do.call(rbind, c(contno_list, bin90245no_list, bin90585no_list))
elevrangesno <- cbind(Species = rownames(elevrangesno), elevrangesno)
write.csv(elevrangesno, "elevationno_ranges_by_scenario.csv", row.names = TRUE)

# Create elevational plot
elevrangesno$Species <- factor(elevrangesno$Species, levels = rev(c("Iris missouriensis", "Bombus appositus", "Bombus flavidus", "Bombus griseocollis", "Bombus huntii", "Bombus insularis", "Bombus melanopygus", "Bombus nevadensis", "Hemaris thetis", "Melissodes tristis", "Papilio rutulus", "Selasphorus platycercus")))    
north_elev<- ggplot()+geom_segment(data=elevrangesno, aes(x=s90585_min, xend=s90585_max, y=Species, yend=Species), col="darkorange2", size=1)+ 
  geom_segment(data=elevrangesno, aes(x=s90245_min, xend=s90245_max, y = as.numeric(factor(Species)) + 0.1, 
                                      yend = as.numeric(factor(Species)) + 0.1), col="deepskyblue2", size=1)+
  geom_segment(data=elevrangesno, aes(x=cont_min, xend=cont_max, y = as.numeric(factor(Species)) + 0.2, 
                                      yend = as.numeric(factor(Species)) + 0.2), col="black", size=1)+
  scale_x_continuous(name= "Elevation (m)", limits= c(-200,6000))+ theme_classic() + theme(
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14, color= "black", face="italic"),   
    axis.title.y = element_text(size = 16), 
    axis.title.x = element_text(size = 16)  
  )


