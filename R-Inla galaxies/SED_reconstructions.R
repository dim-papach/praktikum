rm(list = ls())
options(digits = 16)
require(compiler)
require(beepr)
require(rgl)
require(lattice)
require(stringr)
require(INLA)
require(sp)
require(fields)
require(spam)
require(ggplot2)
require(viridis)
require(yaml)
require(colorspace)
require(spatstat)
require(data.table)
require(gtools)
require(FITSio)
require(IDPmisc)
library(reticulate)
corTesteA <- colorspace::diverge_hsv
set.seed(134775813)
#weighted.median com w = inla
library(gridExtra)
library(grid)
library(plotrix)

library(RColorBrewer)



################################################################### 
######################## USER INTERACTION #########################
###################################################################
 
# number of PCA/NMF components to be used for the reconstruction
num.comp<-10

fit<-'halo16_ph3e08_Inc=0.0_Azi=165.34_total.fits'
fit_ref<-'Au16_Inc=0.0_Azi=165.34_total.fits'


#original fits file

input_fits <- readFITS(fit)
dims =  dim(input_fits$imDat)
features <- 1:dims[3]#
original_pre <- input_fits$imDat[1201:1800,1201:1800, features]
rm(input_fits)
# 
dim_cube<-dim(original_pre)
input_dataset_size <- dim_cube[1]*dim_cube[2]
# 
original_post <- array_reshape(original_pre, dim = c(input_dataset_size,dim_cube[3]))
original_post[which(original_post < 0)] <- 0

rm(original_pre)


###########################3
#original reference

input_fits <- readFITS(fit_ref)
dims <-  dim(input_fits$imDat)    
original_ref_pre <- input_fits$imDat[1201:1800,1201:1800, features]
rm(input_fits)
original_ref_post <- array_reshape(original_ref_pre, dim = c(input_dataset_size,dim_cube[3]))
original_ref_post[which(original_ref_post < 0)] <- 0

rm(original_ref_pre)

#sed
sed<-'halo16_ph3e08_Inc=0.0_Azi=165.34_sed.dat'
dat<-read.table(sed)
inputScale <- dat$V1


# pure INLA

latent_size<-dim_cube[3]

input_fits <- readFITS( 'pure_inla_Mean.fits')
pure_inla_pre <- input_fits$imDat[1:dim_cube[1],1:dim_cube[2],1:latent_size]
pure_inla_post <- array_reshape(pure_inla_pre, dim = c(input_dataset_size,dim_cube[3]))


#PCA+INLA

latent_size<-num.comp
input_fits <- readFITS('pca+inla_Mean.fits')


pca_pre <- input_fits$imDat[1:dim_cube[1],1:dim_cube[2],1:dim_cube[3]]
pca_post <- array_reshape(pca_pre, dim = c(input_dataset_size,dim_cube[3]))
rm(input_fits)



#  NMF+INLA

input_fits <- readFITS('nmf+inla_Mean.fits')


nmf_pre <- input_fits$imDat[1:dim_cube[1],1:dim_cube[2],1:dim_cube[3]]


nmf_post <- array_reshape(nmf_pre, dim = c(input_dataset_size,dim_cube[3]))

rm(input_fits)




xSize <- dim_cube[1]
ySize <- dim_cube[2]
X <- rep(1:xSize, each = ySize)
Y <- rep(1:ySize, times = xSize)


#### figure_5

xs<-c(321,363,461, 535, 256,578)
ys<-c(347,453,294,427, 497,522)

    pdf(file='spaxel_face_on.pdf',width = 10, height = 6) #x=178 y=68
    par(mfrow=c(2,3))
    par(mar=c(4,4.5,1,0.5))
    
    
    for (i in 1:6) {
      
      xs_temp<-xs[i]
      ys_temp<-ys[i]
    
      
      if ((i==1) || (i==4)){
        ylab_name = "log10(Flux density) (MJy/sr)"
      } else {ylab_name = ""}
      
      if ((i==4) || (i==5) || (i==6)){
        xlab_name =  expression(paste("log10(Wavelength (", mu, "m))"))

      } else {xlab_name = ""}
      
      
    
    a<-which(X == xs_temp)
    b<-which(Y[a] == ys_temp)
    
    temp <- log10(c(original_ref_post[a[b],], original_post[a[b],],pure_inla_post[a[b],],pca_post[a[b],], nmf_post[a[b],] ))                    #
    temp[which(is.infinite(temp))] <- NaN                #
    min_int <- min(temp, na.rm = TRUE)                   #
    max_int <- max(temp, na.rm = TRUE)                   #
    #
    matplot(log10(inputScale), cbind(log10(original_ref_post[a[b],]),log10(original_post[a[b],]) ,log10(pure_inla_post[a[b],]),
                                     log10(pca_post[a[b],]),log10(nmf_post[a[b],])), col = 1:5,type = "l",lwd =2,
            xlab = xlab_name, ylab = ylab_name,
            ylim = c(min_int, max_int) , cex.lab=1.4, cex.axis=1.4, cex.main=1.4)                 #
    #axis(3, at = log10(inputScale), labels = round(inputScale))#
    if (i == 1 ){
    legend('bottomright', median(temp, na.rm = TRUE), legend = c("HPN ref.", "LPN input", "pure INLA", "PCA+INLA", "NMF+INLA"),
           col = 1:5, lty = 1:5, cex=1)
    }
    mtext(paste0("Spaxel ", i,": (",xs_temp, ", ",ys_temp, ")"), side=3, 
           col = "blue" , line = -1.5, cex=0.9, at=par("usr")[1]+0.4*diff(par("usr")[1:2]))    #

    
    }
    
    dev.off()
    
  