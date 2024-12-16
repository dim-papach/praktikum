library(INLA)
library(IDPmisc)
require(rasterVis)
library(viridis)
library(latex2exp)
library(fields) #;library(imager);
library(lattice);require(latticeExtra) #;require(INLAutils)
require(classInt);require(reshape2)
library(FITSio)
library(MASS)


setwd("~/Documents/praktikum/inla_project/compare")
source(paste("./inla_fct.R",sep=''))
save_dir <- "./r-inla/"

gal_center <- "galactic_center"
fitsname <- gal_center

threshold <- seq(0, 40, by = 10)  # Creating threshold array

for (i in threshold) {
  
  fitsfile <- sprintf("%s_masked_%s.fits", fitsname, i)
  outfile <- sprintf("%s%s_filled_%s",save_dir, fitsname, i)
  
  fits <- readFITS(file=fitsfile)
  img <- fits$imDat[,]
  
  subset_img <- img[1:40, 1:40]  # This will extract the top-left 20x20 portion of the img array
  # Create a masked array from the data, masking NaN values
  # Your further operations with mask_array or masked_data here
  # Get dimensions of the img array
  dims <- dim(img)

# Create x and y arrays using matrix indexing
  x <- matrix(rep(1:dims[1], dims[2]), nrow = dims[1], ncol = dims[2])
  y <- matrix(rep(1:dims[2], each = dims[1]), nrow = dims[1], ncol = dims[2])


##valid data
  valid <- which(!(is.na(img)))
  xsize <- dims[2]
  ysize <- dims[1]
  xfin <- xsize
  yfin <- ysize


  logimg = log10(img)

## INLA
#imginla <- nonparametric_inla(x[valid],y[valid],logimg[valid], xsize=xsize,ysize=ysize)

  m_weights <- rep(1, nrow(cbind(x, y)))
  imginla <- stationary_inla(x[valid],y[valid],logimg[valid],zoom=1,xsize=xsize,shape='none',xfin=xfin,yfin=yfin,ysize=ysize,cutoff = 0.9, tolerance=1e-6,restart=0L)

  save(x,y,logimg,imginla,file=paste(outfile,'Rdata',sep='.'))
  writeFITSim(imginla$out,file=paste(outfile,'fits',sep='.'))
  
  ##PLOT
  imginla[is.nan(imginla$out[valid])] <- 0
  logimg[is.nan(logimg[valid])] <- 0
  name <-"dim"
  
  Van_Gogh <- c("#263C8B","#4E74A6", "#BDBF78", "#BFA524", "#2E231F") 
  colpal<-colorRampPalette(Van_Gogh)(100)
  cutColor <- 100
  zoom <- 1
  
  ## PLOT
  fj5_img <- classIntervals(c(logimg[valid],imginla$out[valid]), n = cutColor, style= "fisher")
  
  inimg <- levelplot(imginla$image,xlim=c(0,ysize), ylim=c(0,xsize),col.regions=colpal,main = list(name,side=1,line=0.5,cex=1.4),ylab=list('y [pixels]',cex=1.5),cut=cutColor,at=fj5_img$brks,xlab=list('x [pixels]',cex=1.3))
  
  outimg <- levelplot(imginla$out,xlim=c(0,ysize), ylim=c(0,xsize),
                      col.regions=colpal,xlab='x', ylab='y',cut=cutColor,at=fj5_img$brks)
  
  fj5_erimg <- classIntervals(imginla$outsd[valid], n = cutColor, style= "fisher")
  
  outsdimg <- levelplot(imginla$outsd,xlim=c(0,ysize), ylim=c(0,xsize),col.regions=colpal,xlab='x', ylab='y',cut=cutColor,at=fj5_erimg$brks)
  
  save(x,y,logimg,imginla,file=paste(outfile,'Rdata',sep='.'))
  writeFITSim(imginla$out,file=paste(outfile,'fits',sep='.'))
  writeFITSim(imginla$outsd,file=paste(eroutfile,'fits',sep='.'))
  
  png(filename=paste(outfile,'png',sep='.'))
  c(inimg,outimg,outsdimg,layout=c(3,1))
  dev.off()
  
}
  

