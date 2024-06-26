## Download ices areas shapefile
library(rgdal)
icesAreasURL <- "https://gis.ices.dk/shapefiles/ICES_areas.zip"
download.file(icesAreasURL, destfile = "ICES_areas.zip")
## Unzip
shapefiles <- normalizePath(unzip("ICES_areas.zip"))
## Read shapefile
icesShapefile <- shapefiles[grepl(".shp$", shapefiles)]
shape <- readOGR(icesShapefile)
## Change projection to lon, lat
crs <- CRS("+proj=longlat +datum=WGS84")
shape <- spTransform(shape, crs)
proj4string(shape) <- CRS("")
## Save ices areas
saveRDS(shape, "icesAreas.Rds")