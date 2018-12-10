## peperomiaGlobal.R
# Produces an annotated map of the Pacific

## DIRECTORIES ============================
main.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"
fig.dir <- file.path(main.dir, "figures")
## PACKAGES ============================
library(sp); library(rgdal)
library(maps)
library(maptools) # to import maps
library(ggplot2) # plotting
library(rgeos)
library(gpclib)

## CREATE MAP ============================
#https://stackoverflow.com/questions/10620862/use-different-center-than-the-prime-meridian-in-plotting-a-world-map
makeSLmap <- function() {
  llCRS <- CRS("+proj=longlat +ellps=WGS84")
  wrld <- map("world", interior = FALSE, plot=FALSE,
              xlim = c(-179.9, 179.9), ylim = c(-89.9, 89.9))
  wrld_p <- pruneMap(wrld, xlim = c(-179.9, 179.9))
  map2SpatialLines(wrld_p, proj4string = llCRS)
}

## Clip SpatialLines neatly along the antipodal meridian
sliceAtAntipodes <- function(SLmap, lon_0) {
  ## Preliminaries
  long_180 <- (lon_0 %% 360) - 180
  llCRS  <- CRS("+proj=longlat +ellps=WGS84")  ## CRS of 'maps' objects
  eqcCRS <- CRS("+proj=eqc")
  ## Reproject the map into Equidistant Cylindrical/Plate Caree projection 
  SLmap <- spTransform(SLmap, eqcCRS)
  ## Make a narrow SpatialPolygon along the meridian opposite lon_0
  L  <- Lines(Line(cbind(long_180, c(-89.9, 89.9))), ID="cutter")
  SL <- SpatialLines(list(L), proj4string = llCRS)
  SP <- gBuffer(spTransform(SL, eqcCRS), 10, byid = TRUE)
  ## Use it to clip any SpatialLines segments that it crosses
  ii <- which(gIntersects(SLmap, SP, byid=TRUE))
  # Replace offending lines with split versions
  # (but skip when there are no intersections (as, e.g., when lon_0 = 0))
  if(length(ii)) { 
    SPii <- gDifference(SLmap[ii], SP, byid=TRUE)
    SLmap <- rbind(SLmap[-ii], SPii)  
  }
  return(SLmap)
}

## re-center, and clean up remaining streaks
recenterAndClean <- function(SLmap, lon_0) {
  llCRS <- CRS("+proj=longlat +ellps=WGS84")  ## map package's CRS
  #newCRS <- CRS(paste("+proj=eqc +lon_0=", lon_0, sep=""))
  newCRS <- CRS(paste("+proj=moll +lon_0=", lon_0, " +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs", sep = ""))
  
  ## Recenter 
  SLmap <- spTransform(SLmap, newCRS)
  ## identify remaining 'scratch-lines' by searching for lines that
  ## cross 2 of 3 lines of longitude, spaced 120 degrees apart
  # v1 <-spTransform(readWKT("LINESTRING(-62 -89, -62 89)", p4s=llCRS), newCRS)
  # v2 <-spTransform(readWKT("LINESTRING(58 -89, 58 89)",   p4s=llCRS), newCRS)
  # v3 <-spTransform(readWKT("LINESTRING(178 -89, 178 89)", p4s=llCRS), newCRS)
  # ii <- which((gIntersects(v1, SLmap, byid=TRUE) +
  #                gIntersects(v2, SLmap, byid=TRUE) +
  #                gIntersects(v3, SLmap, byid=TRUE)) >= 2)
  # SLmap[-ii]
}

## Put it all together:
Recenter <- function(lon_0 = -100, grid=FALSE, ...) {                        
  SLmap <- makeSLmap()
  SLmap2 <- sliceAtAntipodes(SLmap, lon_0)
  recenterAndClean(SLmap2, lon_0)
}

## Try it out
pacMap <- Recenter(150)
# test <- SpatialLines2PolySet(pacMap)
# test2 <- PolySet2SpatialPolygons(test, close_polys = TRUE)
# test <- SpatialLinesDataFrame(pacMap, data=)
pdf(width = 10, height = 5, file = file.path(fig.dir, "pacificMap.pdf"))
plot(pacMap, col="grey40", lwd = 1) ## Centered on International Date Line
dev.off()