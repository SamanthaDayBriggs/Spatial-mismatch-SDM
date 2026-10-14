library(vegan); library(car); library(stats); library(dplyr); library(tidyr); 
library(visreg); library(bipartite)  

# Set working directory
setwd("~/Desktop/poll diversity")

################# Import data #########################

# Pollinator community data by elevation/year
newdata <- read.csv("poll_obs_year.csv")
newdata$ID <- paste(newdata$Elevation, newdata$Year, sep = "_")
newdata <- column_to_rownames(newdata, var = "ID")
newdata <- newdata[ , !(names(newdata) %in% c("Elevation", "Year"))]

# Create matrix
newdatamatrix <- as.matrix(newdata)

################### Rarefaction ###########################

# Minimum sample size 
min_ind <- min(rowSums(newdatamatrix))

# Rarefy community matrix to equal sampling depth
newdatarare <- rrarefy(newdatamatrix, min_ind)

# Write the csv
write.csv(newdatarare, file = "visitrare.csv")

############# Shannon diversity ##########################

# Diversity (Shannon)
visitdiv <- diversity(visitrare, index = "shannon")

# Format output into dataframe
visitdiv <- data.frame(Site = names(visitdiv), 
                       Diversity = visitdiv, 
                       row.names = NULL)

# Separate Elevation and Year back into columns
visitdiv <- visitdiv %>%
  tidyr::separate(Site, into = c("Elevation", "Year"), sep = "_")

# Write the csv
write.csv(visitdiv, file = "visitdiv.csv", row.names = FALSE)

################ Species richness #############################

# Richness
visitrich <- specnumber(visitrare)

# Format richness output
visitrich <- data.frame(Site = names(visitrich), 
                       Richness = visitrich, 
                       row.names = NULL)

# Separate Elevation and Year back into columns
visitrich <- visitrich %>%
  tidyr::separate(Site, into = c("Elevation", "Year"), sep = "_")

# Write the csv
write.csv(visitrich, file = "visitrich.csv", row.names = FALSE)


################ Statistical analyses ######################

# Read saved diversity/richness csv
visitdiv <- read.csv("visitdiv.csv")
visitrich <- read.csv("visitrich.csv")
visitdiv$Elevation <- as.numeric(visitdiv$Elevation)
visitrich$Elevation <- as.numeric(visitrich$Elevation)

# GLMs for rarefied diversity/richness
gammarich <- glm(Richness ~ Elevation + I(Elevation^2) + Year, 
                 family = Gamma(link = "log"), data = visitrich)

gammadiv <- glm(Diversity ~ Elevation + I(Elevation^2) + Year, 
                family = Gamma(link = "log"), data = visitdiv)

# Model stats
Anova(gammarich, type = "III")
Anova(gammadiv, type = "III")


# Model visualizations
visreg(gammadiv, "Elevation", scale = "response", gg = TRUE, partial = TRUE) + theme_bw()
visreg(gammarich, "Elevation", scale = "response", gg = TRUE, partial = TRUE) + theme_bw()


################ Beta diversity #####################################

# Community data by site
visitdata <- read.csv("poll_obs_all.csv", row.names = 1)

# Create matrix
visitmatrix <- as.matrix(visitdata)

# Minimum sample size
min_indbeta <- min(rowSums(visitmatrix))

# Rarefy community matrix to equal sampling depth
visitrare <- rrarefy(visitmatrix, min_indbeta)

# Beta diversity
betadiv <- betadiver(visitrare, method = 15)
betadiv

################ Generating the network ###############################

# Community data by site
visitdata <- read.csv("poll_obs_all.csv", row.names = 1)

# Plot the web
plotweb(visitdata, spacing = "auto", srt = 90, text_size = 0.5)
