
# generate the files to run the country and TDWG extract functions

library(sf)
library(dplyr)
library(leaflet)
library(rnaturalearth)
library(rredlist)
library(purrr)

#############################################
### Natural Earth - Country boundary data ###
#############################################
# grab the Natural Earth data - https://www.naturalearthdata.com/downloads/10m-cultural-vectors/
get_ne() 

# unzip the ne countries shapefile
unzip(zipfile = "data-raw/ne_10m_admin_0_countries.zip",
      exdir = "data-raw")

# now convert to sf
ne_10 <- st_read("data-raw/ne_10m_admin_0_countries.shp")

# save as RDS to data folder
saveRDS(ne_10, file <- "data/countries_ne_10.rds")

################################################
### WGSRPD - botanical country boundary data ###
################################################
# grab the TDWG data - https://github.com/tdwg/wgsrpd/
tdwgl3 <- get_tdwgl3()

# now convert to sf
l3 <- st_read("data-raw/level3.shp")

# save as RDS to data folder
saveRDS(l3, file <- "data/wgsrpd_l3.rds")

###############################################
### IUCN Red List - country names and codes ###
###############################################
# grab the IUCN Red List country list
# note this uses the IUCN demo token to access the API 
# generate your own token here: http://apiv3.iucnredlist.org/api/v3/token
rl <- rl_countries(key = "9bb4facb6d23f48efbf424bb05c0c1ef1cf6f468393bc745d42179ac4aca5fee")
rl <- rl$results

# save as RDS to data-raw folder
saveRDS(rl, file <- "data/rl.rds")


