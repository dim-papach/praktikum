rm(list = ls())
options(digits = 16)
require(ggplot2)
require(gtools)
require(FITSio)
require(keras)


require(jjb)
require(classInt)
library(reticulate)
library ( pls )
library(ggfortify)
library( lattice )
library(gridExtra)
library(grid)

library(cowplot)
################################################################### 
######################## USER INTERACTION #########################
###################################################################
###################################################################

angle<-c('face','0.0')
#angle<-c('edge','90')

num.comp<-10
n_sample<-0.25
fit<-'halo16_ph3e08_Inc=0.0_Azi=165.34_total.fits'
fit_ref<-'Au16_Inc=0.0_Azi=165.34_total.fits'

input_fits <- readFITS(fit)

dims =  dim(input_fits$imDat)
features <- 1:dims[3]#

original_pre <- input_fits$imDat[1201:1800,1201:1800, features]
#dims =  dim(original_pre)

original_pre[which(original_pre < 0)] <- 0

dim_cube<-dim(original_pre)
input_dataset_size <- dim_cube[1]*dim_cube[2]

original_post <- array_reshape(original_pre, dim = c(input_dataset_size,dim_cube[3]))


rm(input_fits)


###########################
#original 1e8 fits file 

input_fits <- readFITS(fit_ref)

dims <-  dim(input_fits$imDat)    

original_ref_pre <- input_fits$imDat[1201:1800,1201:1800, features]


rm(input_fits)
original_ref_post <- array_reshape(original_ref_pre, dim = c(input_dataset_size,dim_cube[3]))

original_ref_post[which(original_ref_post < 0)] <- 0



#sed
sed<-'halo16_ph3e08_Inc=0.0_Azi=165.34_sed.dat'
dat<-read.table(sed)
inputScale <- dat$V1



# pure INLA

latent_size<-dim_cube[3]


input_fits <- readFITS( 'pure_inla_Mean.fits')

pure_inla_pre <- input_fits$imDat[1:dim_cube[1],1:dim_cube[2],1:dim_cube[3]]

pure_inla_post <- array_reshape(pure_inla_pre, dim = c(input_dataset_size,dim_cube[3]))

#true zeroes zeroes at all lambda for reference 

true_zero_ref<-apply(original_ref_post, 1, sum )

true_zero_pixels_ref <- which(true_zero_ref == 0)

pure_inla_post[true_zero_pixels_ref,] <- 0
pure_inla_pre <- array_reshape(pure_inla_post, dim = dim_cube)


 original_ref_pre_r <- original_ref_pre                       #
 original_ref_pre_r[which(original_ref_pre_r == 0)] <- 1

########## color maps

# sampled original input
set.seed(134775813)
spatial_size<-dim_cube[1]*dim_cube[2]
N <- n_sample * spatial_size
sampled_ind<-c(1:spatial_size)
input_Net_indeces <- sample(sampled_ind, N)#

#######################################
input_Net <- array(NaN, dim = c(spatial_size,dim_cube[3])) 


input_Net[input_Net_indeces, ] <- original_post[input_Net_indeces, ]


input_Net<-  array_reshape(input_Net,
                           dim = dim_cube)

rm(original_post, original_ref_post)

lambda<-c(1:50) #index

#for (i in 1:1) {
for (i in 1:length(lambda)) {
  png(file=paste('color_map',lambda[i],'_part1.png',sep=''))
  
  res_map<- abs(original_ref_pre[,,lambda[i]] - pure_inla_pre[,,lambda[i]]) / original_ref_pre_r[,,lambda[i]] * 100
  res_map[which(is.infinite(res_map))] <- NA  
  
  res_median<-res_map
  res_median[which(res_median == 0)] <- NA
  res_median<-median(res_median, na.rm = TRUE)
  
  res_map<-log10(res_map)
  
  
  
  res_original_map<- abs(original_ref_pre[,,lambda[i]] - original_pre[,,lambda[i]]) / original_ref_pre_r[,,lambda[i]] * 100
  res_original_map[which(is.infinite(res_original_map))] <- NA
  
  res_median_original<-res_original_map
  res_median_original[which(res_median_original == 0)] <- NA
  res_median_original<-median(res_median_original, na.rm = TRUE)
  
  res_original_map<-log10(res_original_map)

  
  res_original<- abs(original_ref_pre[,,lambda[i]] - original_pre[,,lambda[i]]) / original_ref_pre_r[,,lambda[i]] * 100
  
  res_median_original<-res_original
  res_median_original[which(res_median_original == 0)] <- NA
  res_median_original<-median(res_median_original, na.rm = TRUE)
  
  
  res_original<-log10(res_original)
  
  
  #inla residuals
  
  res_map_granice<-res_map   
  res_map_granice[which(is.infinite(res_map_granice))] <- NA 
  
  res_map_granice[which(res_map < ( mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)
  res_map_granice[which(res_map > ( mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)
  
  
  granice<-seq(min(res_map_granice, na.rm = TRUE), max(res_map_granice, na.rm = TRUE),length.out = 100)
  
  
  ######
  res_map_granice_original<-res_original_map    
  
  res_map_granice_original[which(res_original_map < ( mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)
  res_map_granice_original[which(res_original_map > ( mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)
  
  # 
  # 
  ref_granice_pre<-log10(original_ref_pre[,,lambda[i]])
  ref_granice_pre[which(is.infinite(ref_granice_pre))] <- NA 
  
  ref_granice_pre[which(ref_granice_pre < ( mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)
  ref_granice_pre[which(ref_granice_pre > ( mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)
  
  granice_ref<-seq(min(ref_granice_pre, na.rm = TRUE), max(ref_granice_pre, na.rm = TRUE),length.out = 100)
  
  #####################
  
  ref_granice<-log10(input_Net[,,lambda[i]])
  ref_granice[which(is.infinite(ref_granice))] <- NA 
  
  ref_granice[which(ref_granice < ( mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)
  ref_granice[which(ref_granice > ( mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)
  
  ####################################
  ref_granice_inla<-log10(pure_inla_pre[,,lambda[i]])
  ref_granice_inla[which(is.infinite(ref_granice_inla))] <- NA 
  
  ref_granice_inla[which(ref_granice_inla < ( mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)
  ref_granice_inla[which(ref_granice_inla > ( mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)
  #
  
  
  Van_Gogh <- c("#263C8B","#4E74A6", "#BDBF78", "#BFA524", "#2E231F")
  colpal<-colorRampPalette(Van_Gogh)(100)
  cutColor <- 100
  

  my.padding <- list(layout.heights= list(top.padding=0.0, bottom.padding=0.0),
       layout.widths= list(left.padding=0.0, right.padding=0.0))
  

  
  p1<-levelplot(ref_granice,col.regions=colpal,
                cut=cutColor,at=granice_ref, main=list('LPN input'), par.settings = my.padding, 
                scales=list(draw=FALSE),
                xlab='', ylab=''
                #xlab=expression(bold(paste(Delta, 't= 2d 18h 30m'),))
  ) 
  
  
  p2<-levelplot(ref_granice_inla,col.regions=colpal,
                cut=cutColor,at=granice_ref, main=list(paste('pure INLA')), par.settings = my.padding, 
                scales=list(draw=FALSE),
                xlab='', ylab=''
                #xlab=expression(bold(paste(Delta, 't= 5h 40m'),))
  ) 
  
  p3<-levelplot(res_map_granice_original,col.regions=colpal,
                cut=cutColor,at=granice, 
                main=list(paste('Res. [%]',
                                'median=',  trunc(res_median_original)),cex=0.8)
                # main=list(paste('Res. [%] (HPN ref. to LPN input)',
                #     'median=',  format(round(res_median_original, 2), nsmall = 2) ),cex=0.8)
                , par.settings = my.padding, 
                scales=list(draw=FALSE),
                xlab='', ylab='')
  
  
  p4<-levelplot(res_map_granice,col.regions=colpal,
                cut=cutColor,at=granice, 
                main=list(paste('Res. [%] ',
                   'median=', trunc(res_median)),cex=0.8),
                # main=list(paste('Res. [%] (INLA to HPN ref.)', 
                #    'median=', format(round(res_median, 2), nsmall = 2)),cex=0.8),
                
                par.settings = my.padding, scales=list(draw=FALSE),
                xlab='', ylab='')
  
  
  
  grid.arrange(grobs = list(p1,p2, p3,p4), ncol=2, common.legend = TRUE, legend="right")
  
  
  
  dev.off()
  
}




