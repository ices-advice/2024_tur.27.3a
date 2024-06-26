## Survey index calculations
library(icesTAF)
taf.library(DATRAS)
taf.library(mgcv)
taf.library(marmap) ##########
taf.library(surveyIndex)
library(scam)

set.seed(123)

doRetroLo <- F
mc.cores <- 1 # floor(parallel::detectCores() / 2)

d <- readRDS("../surveyData/tur.27.3a_surveyData.Rds")

d <- addSpectrum(d,by = 1)
d <- addWeightByHaul(d,to1min=FALSE)

grid <- getGrid(d, nLon=60)

load("../lfq_no_raising2012_2018.RData")
lfq_all_years <- Reduce(function(x, y) {
  lengths <- sort(unique(c(x$length, y$length)))
  freq <- numeric(length(lengths))
  for (i in seq_along(lengths)) {
    freq[i] <- sum(c(x$freq[x$length == lengths[i]], y$freq[y$length == lengths[i]]), na.rm = TRUE)
  }
  data.frame(length = lengths, freq = freq)
}, lfqs)

lfq_all_years$length = lfq_all_years$length/10
lfq_all_years$length = factor(lfq_all_years$length,levels=as.character(15:76))
d = addSpectrum(d,cm.breaks=15:77)
ds = data.frame( length=15:76,com = as.vector(xtabs(freq~length,lfq_all_years)), surv=colSums(d$N))
dsf = ds[ds$com>0 & ds$surv>0,]

m2 <- scam( log(com / surv )  ~ s(length,bs="mpi"), data=dsf )

# pdf("ESBcorr.pdf")
# par(mfrow=c(3,1),mar=c(4,4,2,2))
# plot(ds$length,ds$com/ds$surv,xlab="Length",ylab="Commercial / Survey")
# title("Selectivity ratio")
mp2 <- exp(predict(m2, newdata = ds))
# lines(ds$length,mp2,lty=2)
corr <- c(rep(0,15),as.numeric(mp2/max(mp2)),rep(1,23))

# plot(corr,xlab="Length",ylab="weighting")
#
# plot(ds$length,ds$com/sum(ds$com),type="l",xlab="Length",ylab="Proportion",lwd=2.5,ylim=c(0,0.08))
# lines(ds$length,ds$surv/sum(ds$surv),lty=2,lwd=2.5)
#
# tmp = ds$surv *mp2
# lines(ds$length,tmp/sum(tmp) ,lty=3,lwd=2.5,col=3)
# legend("topright",legend=c("Survey","Commercial", "Survey re-weighted"),lty=c(2,1,3),col=c(1,1,3),lwd=2.5)
#
# dev.off()
d = addSpectrum(d,cm.breaks=0:100)


# par(mfrow=c(2,1))
# plot(colSums(d$N),type="h")
d$N = t(t(d$N)*corr)
# plot(colSums(d$N),type="h")

d<-addWeightByHaul(d,to1min=FALSE)

d$ctime <- as.numeric(as.character(d$Year)) + round(d$timeOfYear,1)

d$Nage <- matrix(d$HaulWgt, nrow=nrow(d$N),ncol=1)
colnames(d$Nage) <- 1
ages <- 1

knots <- list(timeOfYear=seq(0,1,length=6))
tt <- median(d$timeOfYear[d$Quarter=="1"])

## table(d$ShipG, d$Year)

k.ctime <- 15 ## Increase the number of knots to a higher value than benchmark, 19 last year

##floor(length(unique(d$Year)) / 3) + 2 ## start from the benchmark value of 15 (benchmark one) and add one more knot every 3 years

fm6 <- paste("Gear +",
             "te(lon, lat, bs=c('ds','ds'), k = c(20, 15), m = c(1, 0)) +",
             "te(timeOfYear, lon, lat, bs = c('cc','ds','ds'), k=c(6, 6, 5), m=c(1,0)) +",
             "te(ctime,lon,lat,bs=c('ds','ds','ds'), k=c(", k.ctime, ",8,6), m=c(1,0)) + ", # k of ctime was 15
             "s(Depth, bs='ds', k=5, m=c(1,0)) +",
             "s(ShipG, bs='re', by=dum) +",
             "offset(log(HaulDur))")

system.time(
  m0 <- getSurveyIdx(d,ages,myids=grid[[3]],cutOff=0.01,fam="Tweedie",
                     modelP=fm6,knotsP=knots,nBoot=0,predfix=list(timeOfYear=tt),
                     control= gam.control(trace = TRUE, maxit = 15, nthreads = mc.cores) ) )

save(m0,file="m0.RData")

sink("modelSummary.txt")
summary(m0$pModels[[1]])
sink()

AIC(m0$pModels[[1]])
# qts <- 0:3/4+1/8
# mqs <- list()
# for(i in 1)
#   mqs[[i]] <- redoSurveyIndex(x = d,model = m0,myids=grid[[3]],predfix=list(timeOfYear=qts[i],ShipG = "77AR:GOV"), mc.cores = mc.cores)

getBathyGrid<-function(d,minDepth=15, maxDepth=150, resolution=2,maxDist=Inf, keep=TRUE, shapefile=NULL, select=NULL){
  if (!requireNamespace("marmap", quietly = TRUE)) {
    stop("Package \"marmap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  bathy <- getNOAA.bathy(lon1=min(d$lon),lon2=max(d$lon),lat1=min(d$lat),lat2=max(d$lat),resolution=resolution, keep=keep)
  xyzbathy <- as.data.frame( as.xyz(bathy) )
  colnames(xyzbathy) <- c("lon","lat","Depth")
  xyzbathy$Depth <- -xyzbathy$Depth
  xyzbathy <- subset(xyzbathy, Depth>minDepth & Depth<maxDepth)
  if(is.finite(maxDist)){
    if (!requireNamespace("RANN", quietly = TRUE)) {
      stop("Package \"RANN\" needed for this function to work when 'maxDist' is used. Please install it.",
           call. = FALSE)
    }
    nearest <- RANN::nn2(d[[2]][,c("lon","lat")],xyzbathy[,1:2],k=1)
    xyzbathy <- xyzbathy[ nearest$nn.dists < maxDist ,]
  }
  xyzbathy <- as.data.frame(xyzbathy)
  if (is.character(shapefile)) {
    if (file.exists(shapefile)) {
      shape <- readRDS(shapefile)
    } else {
      stop("shapefile not found")
    }
    tmp <- xyzbathy
    sp::coordinates(tmp) <- ~lon + lat
    xtra <- sp::over(tmp, shape)
    if (!is.null(select))
      xtra <- xtra[select]
    xyzbathy <- cbind(xyzbathy, xtra)
  }
  
  xyzbathy
}

bgrid <- getBathyGrid(d,minDepth=5,maxDist=0.2,resolution=3,shapefile="../icesAreas/icesAreas.Rds",select="Area_27") 

IIIagrid = subset(bgrid,Area_27 %in% c("3.a.20","3.a.21"))
save(IIIagrid, file="IIIagrid.RData")

qts <- 0:3/4+1/8
m3a=list()
for(qq in 1){
  cat("doing quarter ",qq,"\n")
  m3a[[qq]] = redoSurveyIndex(d,m0,myids=NULL,predD=IIIagrid,predfix=list(timeOfYear=qts[qq],ShipG = "77AR:GOV"),nBoot=1001,mc.cores=mc.cores)
}

resQ1 <- data.frame(Index=m3a[[1]]$idx[,1]/mean(m3a[[1]]$idx[,1]),
                    Year=rownames(m3a[[1]]$idx),
                    Quarter=1,
                    lo=as.vector(m3a[[1]]$lo/mean(m3a[[1]]$idx[,1])),
                    up=as.vector(m3a[[1]]$up/mean(m3a[[1]]$idx[,1])),
                    sdlogI=as.vector((log(m3a[[1]]$up)-log(m3a[[1]]$lo))/4 ))

rownames(resQ1)<-NULL
write.csv(resQ1,file="indexQ1.csv", row.names = FALSE)

save(m3a, file="m3a.RData")

if (doRetroLo) {
  system.time({
    retro3a = retro.surveyIdx(m3a[[1]],d = d, grid = NULL,predD=IIIagrid,npeels=5,control=list(trace=TRUE,maxit=10), mc.cores = mc.cores)
    retroindex <- lapply(retro3a, function(x) data.frame(Index = x$idx[,1] / mean(x$idx[,1]),
                                                         Year = as.numeric(rownames(x$idx)),
                                                         lo = as.vector(x$lo/mean(x$idx[,1])),
                                                         up = as.vector(x$up/mean(x$idx[,1])),
                                                         sdlogI = as.vector((log(x$up)-log(x$lo))/4)))
    saveRDS(retroindex, file = "retroindex.Rds")
    saveRDS(retroindex, file = "../../initial/data/retroindex.Rds")
    
    save(retro3a, file = "retro.RData")
    d$Survey2 = as.character(d$Survey)
    d$Survey2[ d$Survey2 %in% c("TNG","TOR") ] <- "TN/TOR"
    d$Survey2=factor(d$Survey2)
    lo3a = leaveout.surveyIdx(m3a[[1]],d,NULL,predD=IIIagrid,fac=d$Survey2,control=list(trace=TRUE,maxit=10), mc.cores = mc.cores)
    save(lo3a, file = "lo.RData")
  })
}

## Save session info
sink("sessionInfo.txt");sessionInfo();sink()