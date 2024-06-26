## Run analysis, write model results

## Before: tur.27.3a_catch.csv, tur.27.3a_index.csv
## After: tur.27.3a_fit.Rds


library(icesTAF)
taf.library(TMB)
taf.library(spict, messages = TRUE)

mkdir("model")

icesTAF::msg("Model: read in data")

## Load catch and index data
index <- read.taf("data/tur_27_3a_index.csv")
catch <- read.taf("data/tur_27_3a_catch.csv")

## What is the last year in the data?
lastYear <- max(catch$Year)

## Management period and evaluation time
maninterval <- seq(lastYear + 2, lastYear + 3)
maneval <- lastYear + 3

## Priors
logbkfrac <- log(0.5)
priors <- list(
  logn = c(0,0,0), ## No prior for n (fixed to Schaefer below)
  logbkfrac = c(logbkfrac, 0.5, 1)
)

baseinp <- inp <- list(
  timeC = catch$Year,
  obsC  = catch$Total,
  timeI = index$Year,
  obsI  = index$Index,
  stdevfacI = index$sdlogI / mean(index$sdlogI),
  stdevfacC = c(rep(3, length(1975:2001)), rep(1, length(2002:lastYear))),
  priors = priors,
  maninterval = maninterval,
  maneval = maneval,
  ini = list(logn = log(2)), phases = list(logn = -1), ## Schaefer model
  optimiser.control = list(iter.max = 1e4, eval.max = 1e4)
)

inp <- check.inp(inp)

icesTAF::msg("Model: model fit")
fit <- fit.spict(inp)
fit <- calc.osa.resid(fit)
fit <- calc.process.resid(fit)


### ###############################33

#fitold <- readRDS("~/git-repos/TAF/2022_tur.27.3a_assessment/model/tur.27.3a_fit.Rds")

#plotspict.compare(list(fitold = fitold, fit = fit), plot.unc = 2)

# par(mfrow = c(1,1))
# plot(fit$inp$timeI[[1]], fit$inp$obsI[[1]], type = "b", col = adjustcolor(1, 0.8), pch = 20)
# lines(fitold$inp$timeI[[1]], fitold$inp$obsI[[1]], type = "b", col = adjustcolor(2, 0.8), pch = 20)
#
# indapr5 <- read.csv("data/tur_27_3a_index_Apr5.csv")
# indapr16 <- read.csv("data/tur_27_3a_index_Apr16.csv")
# indold <- read.csv("~/git-repos/TAF/2022_tur.27.3a_assessment/data/tur_27_3a_index.csv")
# lwd <- 2
# plot(fit$inp$timeI[[1]], fit$inp$obsI[[1]], type = "b", col = adjustcolor(2, 0.8), pch = 20, lwd = lwd)
# lines(indold$Year, indold$Index, type = "b", col = adjustcolor(1, 0.8), pch = 20, lwd = lwd)
# lines(indapr5$Year, indapr5$Index, type = "b", col = adjustcolor(4, 0.8), pch = 20, lwd = lwd)
# lines(indapr16$Year, indapr16$Index, type = "b", col = adjustcolor(5, 0.8), pch = 20, lwd = lwd)
# legend("top", legend = c("Last year", "Current", "Apr5", "Apr16"), col = c(1, 2, 4, 5), pch = 20, lty = 1, lwd = lwd)

### #################################

fitExtraFoptions <- fit

icesTAF::msg("Model: management scenarios")
fit <- add.man.scenario(fit, "F=Fmsy_C_fractile", fractiles = list(catch = 0.35), breakpointB = c(1/2))
fit <- add.man.scenario(fit, "F=Fmsy", breakpointB = c(1/2))
fit <- add.man.scenario(fit, "F=Fsq", ffac = 1)
fit <- add.man.scenario(fit, "F=0", ffac = 0)
fit <- add.man.scenario(fit, "F=Fmsy_All_fractiles", fractiles = list(catch = 0.35, bbmsy = 0.35, ffmsy = 0.35), breakpointB = c(1/2))

icesTAF::msg("Model: retro")
fit <- retro(fit)

icesTAF::msg("Model: saving the results")
saveRDS(fit, "model/tur.27.3a_fit.Rds")

icesTAF::msg("Model: extra Foption management scenarios")
## F options from 0.01 to upper 95% bound of Fmsy estimate
for (fopt in seq(0.01, get.par("Fmsy", fit)[3], 0.01 )) {
  fitExtraFoptions <- add.man.scenario(fitExtraFoptions, paste0("F=", fopt), fractiles = list(catch = 0.35), fabs = fopt)
}

saveRDS(fitExtraFoptions, "model/tur.27.3a_fit_extraFoptions.Rds")


icesTAF::msg("Model: correct retro")
correctRetro <- fit
retroindex <- readRDS("bootstrap/data/retroindex.Rds")
retroinps <- lapply(retroindex, function(indx) {
  retroinp <- baseinp
  retroinp$timeI <- indx$Year
  retroinp$obsI <- indx$Index
  retroinp$stdevfacI <- indx$sdlogI / mean(indx$sdlogI)
  keep <- retroinp$timeC %in% c(1975:1982, indx$Year)
  retroinp$timeC <- retroinp$timeC[keep]
  retroinp$obsC<- retroinp$obsC[keep]
  lstyr <- max(indx$Year)
  retroinp$stdevfacC <- c(rep(3, length(1975:2001)), rep(1, length(2002:lstyr)))
  retroinp$maneval <- NULL
  retroinp$maninterval <- NULL
  check.inp(retroinp)
})
inpretro <- c(list(baseinp), retroinps)
correctRetro$retro <- lapply(inpretro, fit.spict)
saveRDS(correctRetro, file = "model/tur.27.3a_fit_correctRetro.Rds")

# icesTAF::msg("Model: extra run using sdI and sdB priors")
# from1975 <- catch$Year >= 1975
# inp <- list(
#   timeC = catch$Year[from1975],
#   obsC  = catch$Total[from1975],
#   timeI = index$Year,
#   obsI  = index$Index,
#   stdevfacI = index$sdlogI,
#   stdevfacC = c(rep(3, length(1975:2001)), rep(1, length(2002:lastYear))),
#   priors = list(
#     logn = c(0,0,0),
#     logalpha = c(0,0,0),
#     logsdi = c(log(0.5), 0.5, 1),
#     logsdb = c(log(0.5), 0.5, 1),
#     logbkfrac = c(logbkfrac, 0.5, 1)
#   ),
#   maninterval = maninterval,
#   maneval = maneval,
#   ini = list(logn = log(2)),
#   phases = list(logn = -1),
#   optimiser.control = list(iter.max = 1e4, eval.max = 1e4)
# )
#
# inp3c <- check.inp(inp)
# fit3c <- fit.spict(inp3c)
# fit3c <- calc.osa.resid(fit3c)
# fit3c <- retro(fit3c)
# fit3c <- add.man.scenario(fit3c, "F=Fmsy_C_fractile", fractiles = list(catch = 0.35), breakpointB = c(1/2))
# saveRDS(fit3c, "model/tur.27.3a_fit_2021_priors_for_sdb_sdi.Rds")

