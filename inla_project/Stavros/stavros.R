# Library Loading

library(INLA)
library(IDPmisc)
require(rasterVis)
library(viridis)
library(latex2exp)
library(fields) #;library(imager);
library(lattice);require(latticeExtra) #;require(INLAutils)
require(classInt);require(reshape2)
library(FITSio)
library(rmarkdown)



# File and Directory Setup


home <- ("/home/dp/")
source(paste("./inla_fct.R",sep=''))
name <- "stavros"


dir <- paste("./first_results/",sep='')
file <- paste(dir,"HI_6563s",sep='')
outfile <- paste(file,"inla",sep='-')
eroutfile <- paste(file,"erinla",sep='-')

#read file
fits <- readFITS(file=paste(file,"fits",sep='.'))
img <- fits$imDat[,]


## Starry night or how to set a new color pallete

Van_Gogh <- c("#263C8B","#4E74A6", "#BDBF78", "#BFA524", "#2E231F") 
colpal<-colorRampPalette(Van_Gogh)(100)
cutColor <- 100
zoom <- 1


### The *cutoff* variable

cutoff <- 0.9 #bit lower than minimum distance


subset_img <- img[1:20, 1:20]  # This will extract the top-left 20x20 portion of the img array


## x and y
dims <- dim(img)
x <- array(0,dim=dim)
y <- array(0,dim=dim)
for (i in 1:dim[1]){
    for (j in 1:dim[2]){
        x[i,j] <- i
        y[i,j] <- j
    }
}


##valid data
valid <- which(!(img < 27 | is.na(img)))
xsize <- dims[2]
ysize <- dims[1]
xfin <- xsize
yfin <- ysize

##normalize data
logimg = log10(img)

summary(logimg)
## INLA
imginla <- stationary_inla(x[valid],y[valid],logimg[valid],zoom=zoom,xsize=xsize,shape='none',xfin=xfin,yfin=yfin,ysize=ysize,cutoff=cutoff,tolerance=1e-6,restart=0)

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

