library(icesTAF)
library(icesDatras)
library(maps); library(mapdata)
taf.library(DATRAS)
taf.library(mgcv)
taf.library(surveyIndex)

## Species specific parameters:
lastYear <- 2023 ## as.integer(substr(Sys.Date(), 1, 4)) - 1 ## find last year
yearsibts <- (1983:lastYear)
yearsbits <- 1991:lastYear
yearsbts <- 1987:lastYear
genus="Psetta"
bfamily="maxima"
mc.cores=parallel::detectCores() - 2
outFolder <- path.expand(".")
icesAreas <- c("3.a.20", "3.a.21", "4.a", "4.b", "4.c", "22", "23", "24")
icesShapefile <- readRDS("../icesAreas/icesAreas.Rds")
latrange <-  c(53,59)
lonmin <- 5

## 00. Helper functions ----

## Turbot used to have 2 scientific names Psetta maxima and Scophtalmus maximus
## Use Aphia (turbot has 127149) to save all entries as P. maxima
addPsetta <- function(x) {
  if(! "Psetta maxima" %in% levels(x[["HL"]]$Species)) {
    levels(x[["HL"]]$Species) <- c(levels(x[["HL"]]$Species), "Psetta maxima")
  }
  if(! "Psetta maxima" %in% levels(x[["CA"]]$Species)) {
    levels(x[["CA"]]$Species) <- c(levels(x[["CA"]]$Species), "Psetta maxima")
  }
  x[["HL"]]$Species[x[["HL"]]$Valid_Aphia == 127149] <- "Psetta maxima"
  x[["CA"]]$Species[x[["CA"]]$Valid_Aphia == 127149] <- "Psetta maxima"
  x[["HL"]]$Species[(x[["HL"]]$SpecCode == 127149 & x[["HL"]]$SpecCodeType == "W")] <- "Psetta maxima"
  x[["CA"]]$Species[(x[["CA"]]$SpecCode == 127149 & x[["CA"]]$SpecCodeType == "W")] <- "Psetta maxima"
  subset(x,Species=="Psetta maxima")
}


## 01. Read in data ----

bits <- readRDS("../BITS_1991-2023.Rds")
bts <- readRDS("../BTS_1987-2023.Rds")
ibts <- readRDS("../NS-IBTS_1983-2023.Rds")
tng <- readRDS("../TNG_2004-2023.Rds")
tor <- readRDS("../TOR_2008-2023.Rds")


## 02. Subset ----
d <- c(bits, bts, ibts, tng, tor)
d <- addPsetta(d)

#dorig <- d
d = subset(d, HaulVal=="V")
d <- addSpatialData(d, icesShapefile)
d <- subset(d, Area_27 %in% c("4.a","4.b","4.c","3.a.20","3.a.21","3.c.22","3.b.23","3.d.24"))
geartab <- sort(table(d$Gear), TRUE)
table(d$Gear, d$Year)

d <- subset(d,
            Gear %in% names(geartab[geartab > 120]),
            lat >= latrange[1]  & lat <= latrange[2],
            lon > lonmin)
## d <- subset(d, !(d$lon < 8.5 & d$lat < 56))

d <- addSpectrum(d)
d <- addWeightByHaul(d)

d$ShipG <- factor(paste(d$Ship,d$Gear,sep=":"))

table(d$Year, d$ShipG)

# DATRAS:::plot.DATRASraw(d, col = adjustcolor(as.integer(d$Survey), 0.1),
#                         plot.squares = FALSE)
# maps::map("mapdata::worldHires", add = TRUE, fill = TRUE, col = "#12345622")
# legend("topright", legend = levels(d$Survey), col= 1:5, pch = 20)

## 03. Add missing depths ----
md <- gam(Depth ~ s(lon, lat, bs="ds"), data = d[["HH"]])
##plot(md)
d[["HH"]][is.na(d[['HH']]$Depth), "Depth"] <- predict(md, d[["HH"]][is.na(d[['HH']]$Depth),c("lon", "lat")] , type = "response")


## 04. Save all data
saveRDS(d, "tur.27.3a_surveyData.Rds")




## XX. Code that uses the ICES DATRAS web services ----
# There are unfortunatelly issues, probasbly have to do with dummy data returned by the WS

## From DATRAS package - filters out dummy data
# getDatrasExchange <- function (survey, years, quarters, strict = TRUE) {
#   ## Exclude dummy data - contacted ICES DATRAS about that - dummy data do not appear in the exchange zip files
#   ## Exclude hauls with country being missing or DUM, or ship being unknown (code: AA36)
#   filterDummy <- function(x) {subset(x, !is.na(Country) & !startsWith(Country, "DUM") & !startsWith(Ship, "AA36"))}
#   ca <- getDATRAS("CA", survey = survey, years = years, quarters = quarters)
#   if (identical(ca, FALSE)) {
#     stop()
#   }
#   hh <- getDATRAS("HH", survey = survey, years = years, quarters = quarters)
#   hl <- getDATRAS("HL", survey = survey, years = years, quarters = quarters)
#   d <- list(CA = filterDummy(ca), HH = filterDummy(hh), HL = filterDummy(hl))
#   cat("Classes of the variables\n")
#   print(lapply(d, sapply, class))
#   for (i in 1:3) d[[i]] <- renameDATRAS(d[[i]])
#   d <- minus9toNA(d)
#   if (is.null(d$CA$StatRec))
#     d$CA$StatRec <- d$CA$AreaCode
#   d <- DATRAS:::addExtraVariables(d)
#   d <- DATRAS:::fixMissingHaulIds(d, strict = strict)
#   for (i in 1:3) {
#     for (k in 1:ncol(d[[i]])) {
#       if (is.character(d[[i]][, k]))
#         d[[i]][, k] <- factor(d[[i]][, k])
#     }
#   }
#   class(d) <- "DATRASraw"
#   d
# }

# get_DATRAS <- function(surveys, years, genus, bfamily, common, quarter,
#                        icesAreas, icesShapefile,
#                        path, exchangePath) {
#   fn <- paste0(path, "/", common, "_", paste(surveys, collapse="+"), "_Q",
#                paste(quarter, collapse = "+"), "_", paste(range(years), collapse = "-"),".RData")
#
#   if (! file.exists(fn)) {
#     dds <- list()
#     for (survey in surveys) {
#       dds[[survey]] <- getDatrasExchange(survey, years, quarters = quarter, strict = FALSE)
#       dds[[survey]] <- addSpatialData(dds[[survey]], icesShapefile)
#     }
#
#     if(length(dds) > 1) {
#       dAll <- Reduce(c, dds)
#     } else {
#       dAll <- dds[[1]]
#     }
#     dAll <- addPsetta(dAll)
#     dAll <- subset(dAll, Year %in% years,
#                    HaulVal == "V", StdSpecRecCode == 1, Quarter %in% quarter)
#
#     ## impute missing depths
#     if (any(is.na(dAll[[2]]$Depth))) {
#       dmodel <- gam(log(Depth) ~ s(lon, lat, k = 200), data = dAll[[2]])
#       sel <- subset(dAll, is.na(Depth))
#       sel$Depth <- 0; ## Guard against NA-error
#       dAll$Depth[is.na(dAll$Depth)] <- exp(predict(dmodel, newdata = sel[[2]]))
#       sel <- dmodel <- NULL; gc(verbose = FALSE)
#     }
#
#     d <- subset(dAll, Month %in% quarter2months(quarter),
#                 Area_27 %in% icesAreas)
#     save(d, file = fn)
#   }
#   var <- load(fn)
#   if (length(var) > 1) warning("More variables found in ", fn, ". Only the first is returned.")
#   assign("res", get(var))
#   res
# }
#
# quarter2months <- function(quarter) as.vector(sapply(quarter, function(q) ((q - 1) * 3 + 1) : (q * 3)))
#
#
# system.time(ibts <- get_DATRAS(surveys = c("NS-IBTS"), years = yearsibts,
#                                genus = genus, bfamily = bfamily, common = "turbot",
#                                quarter = 1:4, icesAreas = icesAreas, icesShapefile = icesShapefile,
#                                path = ".", exchangePath = "."))
# system.time(bits <- get_DATRAS(surveys = c("BITS"), years = yearsbits,
#                                genus = genus, bfamily = bfamily, common = "turbot",
#                                quarter = 1:4, icesAreas = icesAreas, icesShapefile = icesShapefile,
#                                path = ".", exchangePath = "."))
# system.time(bts <- get_DATRAS(surveys = c("BTS"), years = yearsbts,
#                               genus = genus, bfamily = bfamily, common = "turbot",
#                               quarter = 1:4, icesAreas = icesAreas, icesShapefile = icesShapefile,
#                               path = ".", exchangePath = "."))
#
# dAll <- c(ibts, bits, bts)
# dAll <- addSpectrum(dAll)
# dAll <- addWeightByHaul(dAll)
#
# save(dAll, file = "turbot_NS-IBTS+BITS+BITS_1983-2021.RData")